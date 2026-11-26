#!/bin/bash
set -euo pipefail

sshd=$(command -v sshd || command -v /usr/sbin/sshd || true)
if [[ -z "$sshd" ]]; then
    echo 'Skipping publish output test: sshd is unavailable.'
    exit 0
fi
if [[ "$sshd" != /* ]]; then
    sshd="$PWD/$sshd"
fi

script_dir=$(cd "$(dirname "$0")" && pwd)
fixture=$(mktemp -d)
trap '{ [[ -f "$fixture/sshd.pid" ]] && kill "$(cat "$fixture/sshd.pid")" || true; } 2>/dev/null; rm -rf "$fixture"' EXIT
cd "$fixture"
export HOME="$fixture/home" TMPDIR="$fixture/temp with spaces"
mkdir -p "$HOME" "$TMPDIR" markdown

cp "$script_dir/publish_benchmark_output.sh" .
ssh-keygen -q -t ed25519 -N '' -f host_key
ssh-keygen -q -t ed25519 -N '' -f deploy_key

git init -q -b main seed
git -C seed -c user.name=Test -c user.email=test@example.org commit -q --allow-empty -m seed
git clone -q --bare seed output.git
printf 'benchmark result\n' > markdown/page.md
export GITHUB_REPOSITORY=owner/benchmarks GITHUB_SHA=fixture
DEPLOY_KEY=$(cat deploy_key)
export DEPLOY_KEY

# cygpath shim for the emulation leg: -m/-w give drive-letter forms, -u and a
# missing flag pass the POSIX argument through, like the real tool.
mkdir -p bin
cat > bin/cygpath <<'EOF'
#!/bin/sh
mode=
case "${1:-}" in
    -*) mode=$1; shift ;;
esac
case "$mode" in
    -m) printf '%s\n' "$FAKE_GIT_NATIVE${1#"$FAKE_GIT_ROOT"}" ;;
    -w) printf '%s\n' "$FAKE_GIT_NATIVE${1#"$FAKE_GIT_ROOT"}" | tr '/' '\' ;;
    *)  printf '%s\n' "$1" ;;
esac
EOF
chmod +x bin/cygpath

# --- Leg 1: emulated Git-for-Windows install under a spaced root. ----------
export FAKE_GIT_ROOT="$fixture/Program Files/Git"
export FAKE_GIT_NATIVE='C:/Program Files/Git'
mkdir -p "$FAKE_GIT_ROOT/cmd" "$FAKE_GIT_ROOT/usr/bin"
printf '#!/bin/sh\nexec "%s" "$@"\n' "$(command -v git)" > "$FAKE_GIT_ROOT/cmd/git"
printf '#!/bin/sh\nexec "%s" "$@"\n' "$(command -v ssh)" > "$FAKE_GIT_ROOT/usr/bin/ssh.exe"
chmod +x "$FAKE_GIT_ROOT/cmd/git" "$FAKE_GIT_ROOT/usr/bin/ssh.exe"

cat > sshd_config <<EOF_CONFIG
HostKey "$fixture/host_key"
StrictModes no
AuthorizedKeysFile "$fixture/deploy_key.pub"
PasswordAuthentication no
KbdInteractiveAuthentication no
EOF_CONFIG
export TEST_PROXY="ProxyCommand=\"$sshd\" -i -f \"$fixture/sshd_config\""
export TEST_OUTPUT_REPO="$fixture/output.git"
cat > test-publish-ssh <<'EOF_TRANSPORT'
#!/bin/bash
set -eu
service=$1
# Emulate native Git parsing GIT_SSH_COMMAND: drive-letter paths run as-is,
# POSIX paths are resolved under the install root the way MSYS maps them.
eval "set -- ${GIT_SSH_COMMAND:-ssh}"
case "$1" in
    [A-Za-z]:/*) set -- "$FAKE_GIT_ROOT${1#"$FAKE_GIT_NATIVE"}" "${@:2}" ;;
    /*)          set -- "$FAKE_GIT_ROOT$1" "${@:2}" ;;
esac
exec "$@" -o "$TEST_PROXY" "$(id -un)@github.com" "$service '$TEST_OUTPUT_REPO'"
EOF_TRANSPORT
chmod +x test-publish-ssh

emulate_publish() {
    PATH="$FAKE_GIT_ROOT/cmd:$FAKE_GIT_ROOT/usr/bin:$fixture/bin:$fixture:$PATH" \
        GIT_ALLOW_PROTOCOL=ext PUBLISH_OUTPUT_REPO='ext::test-publish-ssh %S' \
        ./publish_benchmark_output.sh
}
printf 'github.com %s\n' "$(cat host_key.pub)" > github_known_hosts
test -s github_known_hosts
grep -q '^github.com ssh-ed25519' github_known_hosts

emulate_publish
git --git-dir=output.git show HEAD:markdown/page.md | cmp markdown/page.md -
echo 'Published the expected artifact over SSH with an empty runner home.'

test -s "$script_dir/github_known_hosts"
grep -q '^github.com ' "$script_dir/github_known_hosts"
pinned=$(ssh-keygen -lf "$script_dir/github_known_hosts" | awk '$3=="github.com"{print $2}')
[[ "$pinned" == 'SHA256:+DiY3wvvV6TuJJhbpZisF/zLDA0zPMSvHdkr4UvCOqU' ]]
cp "$script_dir/github_known_hosts" github_known_hosts
if emulate_publish > rejected.log 2>&1; then
    echo 'Unexpectedly accepted a different SSH host key.' >&2
    exit 1
fi
grep -F 'Host key verification failed.' rejected.log
test -z "$(ls -A "$TMPDIR")"
echo 'Rejected the untrusted host and cleaned up temporary credentials.'

# --- Leg 2: real Git for Windows under Wine ---------------------------------
if [[ -z "${GIT_FOR_WINDOWS:-}" && -z "${WINE:-}" ]]; then
    echo 'Set GIT_FOR_WINDOWS to run the native Git Bash publish test under Wine.'
    exit 0
fi
: "${GIT_FOR_WINDOWS:?Set GIT_FOR_WINDOWS to a Git for Windows installation}"
wine=${WINE:-$(command -v wine || command -v wine64 || true)}
if ! wine=$(command -v "$wine"); then
    echo 'Wine is required for the configured native Git Bash test.' >&2
    exit 1
fi
for exe in bin/bash.exe cmd/git.exe usr/bin/{ssh,cygpath,mktemp}.exe; do
    [[ -x "$GIT_FOR_WINDOWS/$exe" ]] || { echo "Missing executable: $GIT_FOR_WINDOWS/$exe" >&2; exit 1; }
done

wport=$((30000 + RANDOM % 20000))
cat > sshd_listen <<EOC
Port $wport
ListenAddress 127.0.0.1
HostKey "$fixture/host_key"
AuthorizedKeysFile "$fixture/deploy_key.pub"
PasswordAuthentication no
KbdInteractiveAuthentication no
PermitRootLogin yes
StrictModes no
UsePAM no
PidFile "$fixture/sshd.pid"
EOC
"$sshd" -f "$fixture/sshd_listen"
git clone -q --bare seed output2.git
printf '[127.0.0.1]:%s %s\n' "$wport" "$(cat host_key.pub)" > github_known_hosts
PUBLISH_OUTPUT_REPO="ssh://$(id -un)@127.0.0.1:$wport$fixture/output2.git"
export PUBLISH_OUTPUT_REPO
ln -s "$GIT_FOR_WINDOWS" "$fixture/Program Files/Native Git"
cat > git-bash-env <<'EOF'
unset BASH_ENV
export HOME="$PWD/home" TMPDIR="$PWD/temp with spaces"
ssh_path=$(cygpath -m /usr/bin/ssh.exe)
if [[ "$ssh_path" != [A-Za-z]:/*" "*/usr/bin/ssh.exe ]]; then
    echo "Unexpected native SSH path: $ssh_path" >&2
    exit 1
fi
printf 'Git Bash uses %s\n' "$ssh_path"
EOF
(
    export WINEPREFIX="$fixture/wine-prefix" WINEDEBUG=-all BASH_ENV=./git-bash-env
    wineserver=$(command -v "$(dirname "$wine")/wineserver" || command -v wineserver)
    trap '"$wineserver" -w' EXIT
    # bin/bash.exe initializes the bundled MSYS environment. Git 2.56.1 runs with Wine Staging 11.22.
    "$wine" "$fixture/Program Files/Native Git/bin/bash.exe" --noprofile --norc ./publish_benchmark_output.sh
)
git --git-dir=output2.git show HEAD:markdown/page.md | cmp markdown/page.md -
echo 'Published the expected artifact through Git Bash under Wine (byte-for-byte match).'
