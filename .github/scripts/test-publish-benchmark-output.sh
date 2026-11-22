#!/bin/bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
fixture=$(mktemp -d)
trap 'rm -rf "$fixture"' EXIT
cd "$fixture"
export HOME="$fixture/home" TMPDIR="$fixture/temp with spaces"
mkdir -p "$HOME" "$TMPDIR" markdown
cp "$script_dir/publish_benchmark_output.sh" .
ssh-keygen -q -t ed25519 -N '' -f host_key
ssh-keygen -q -t ed25519 -N '' -f deploy_key

# Feed each SSH connection directly to sshd over stdio, entirely offline.
cat > sshd_config <<EOF_CONFIG
HostKey "$fixture/host_key"
StrictModes no
AuthorizedKeysFile "$fixture/deploy_key.pub"
PasswordAuthentication no
KbdInteractiveAuthentication no
EOF_CONFIG
sshd=$(command -v sshd)
export TEST_PROXY="ProxyCommand=\"$sshd\" -i -f \"$fixture/sshd_config\""
export TEST_OUTPUT_REPO="$fixture/output.git"
cat > test-publish-ssh <<'EOF_TRANSPORT'
#!/bin/bash
set -eu
git_host=$1
# Emulate Git for Windows: a POSIX program path in GIT_SSH_COMMAND is mapped
# onto the install root and respawned unquoted, splitting a spaced root.
eval "set -- ${GIT_SSH_COMMAND:-ssh}"
if [[ -n "${FAKE_GIT_ROOT:-}" && "$1" == /* && "$1" != "$FAKE_GIT_ROOT"/* && -f "$FAKE_GIT_ROOT$1" ]]; then
    read -ra split_prog <<<"$FAKE_GIT_ROOT$1"
    set -- "${split_prog[@]}" "${@:2}"
fi
exec "$@" -o "$TEST_PROXY" "$(id -un)@github.com" "$git_host '$TEST_OUTPUT_REPO'"
EOF_TRANSPORT
chmod +x test-publish-ssh
# A Git "install" under a spaced root, laid out like Git for Windows.
export FAKE_GIT_ROOT="$fixture/Program Files/Git"
mkdir -p "$FAKE_GIT_ROOT/cmd" "$FAKE_GIT_ROOT/usr/bin"
printf '#!/bin/sh\nexec "%s" "$@"\n' "$(command -v git)" > "$FAKE_GIT_ROOT/cmd/git"
chmod +x "$FAKE_GIT_ROOT/cmd/git"
ln -s "$(command -v ssh)" "$FAKE_GIT_ROOT/usr/bin/ssh"
export PATH="$FAKE_GIT_ROOT/cmd:$fixture:$PATH"
git init -q -b main seed
git -C seed -c user.name=Test -c user.email=test@example.org commit -q --allow-empty -m seed
git clone -q --bare seed output.git
export GIT_ALLOW_PROTOCOL=ext PUBLISH_OUTPUT_REPO='ext::test-publish-ssh %S'
printf 'benchmark result\n' > markdown/page.md
export DEPLOY_KEY GITHUB_REPOSITORY=owner/benchmarks GITHUB_SHA=fixture
DEPLOY_KEY=$(cat deploy_key)
printf 'github.com %s\n' "$(cat host_key.pub)" > github_known_hosts
test -s github_known_hosts
grep -q '^github.com ssh-ed25519' github_known_hosts

./publish_benchmark_output.sh
git --git-dir=output.git show HEAD:markdown/page.md | cmp markdown/page.md -
echo 'Published the expected artifact over SSH with an empty runner home.'

test -s "$script_dir/github_known_hosts"
grep -q '^github.com ' "$script_dir/github_known_hosts"
cp "$script_dir/github_known_hosts" github_known_hosts
if ./publish_benchmark_output.sh > rejected.log 2>&1; then
    echo 'Unexpectedly accepted a different SSH host key.' >&2
    exit 1
fi
grep -F 'Host key verification failed.' rejected.log
test -z "$(ls -A "$TMPDIR")"
echo 'Rejected the untrusted host and cleaned up temporary credentials.'
