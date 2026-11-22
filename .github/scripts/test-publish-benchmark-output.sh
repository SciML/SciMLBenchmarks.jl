#!/bin/bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
fixture=$(mktemp -d)
trap '{ [[ -f "$fixture/sshd.pid" ]] && kill "$(cat "$fixture/sshd.pid")" || true; } 2>/dev/null; rm -rf "$fixture"' EXIT
cd "$fixture"
export HOME="$fixture/home" TMPDIR="$fixture/temp with spaces"
mkdir -p "$HOME" "$TMPDIR" markdown

cp "$script_dir/publish_benchmark_output.sh" .
ssh-keygen -q -t ed25519 -N '' -f host_key
ssh-keygen -q -t ed25519 -N '' -f deploy_key
sshd=$(command -v sshd || command -v /usr/sbin/sshd || true)
[[ -n "$sshd" ]] || { echo 'sshd is required for this test.' >&2; exit 1; }

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

export PATH="$FAKE_GIT_ROOT/cmd:$FAKE_GIT_ROOT/usr/bin:$fixture/bin:$fixture:$PATH"
export GIT_ALLOW_PROTOCOL=ext PUBLISH_OUTPUT_REPO='ext::test-publish-ssh %S'
printf 'github.com %s\n' "$(cat host_key.pub)" > github_known_hosts
test -s github_known_hosts
grep -q '^github.com ssh-ed25519' github_known_hosts

./publish_benchmark_output.sh
git --git-dir=output.git show HEAD:markdown/page.md | cmp markdown/page.md -
echo 'Published the expected artifact over SSH with an empty runner home.'

test -s "$script_dir/github_known_hosts"
grep -q '^github.com ' "$script_dir/github_known_hosts"
pinned=$(ssh-keygen -lf "$script_dir/github_known_hosts" | awk '$3=="github.com"{print $2}')
[[ "$pinned" == 'SHA256:+DiY3wvvV6TuJJhbpZisF/zLDA0zPMSvHdkr4UvCOqU' ]]
cp "$script_dir/github_known_hosts" github_known_hosts
if ./publish_benchmark_output.sh > rejected.log 2>&1; then
    echo 'Unexpectedly accepted a different SSH host key.' >&2
    exit 1
fi
grep -F 'Host key verification failed.' rejected.log
test -z "$(ls -A "$TMPDIR")"
echo 'Rejected the untrusted host and cleaned up temporary credentials.'
unset GIT_ALLOW_PROTOCOL

# --- Leg 2: real Git for Windows under Wine ---------------------------------
wine=${WINE:-$(command -v wine || command -v wine64 || true)}
if [[ -n "$wine" && -n "${GIT_FOR_WINDOWS:-}" && -x "$GIT_FOR_WINDOWS/cmd/git.exe" && -x "$GIT_FOR_WINDOWS/usr/bin/ssh.exe" && -x "$GIT_FOR_WINDOWS/usr/bin/cygpath.exe" ]]; then
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
    if "$sshd" -f sshd_listen; then
        export FAKE_GIT_ROOT="$fixture/Program Files/Git-for-Windows"
        export FAKE_GIT_NATIVE="Z:$FAKE_GIT_ROOT"
        export GIT_FOR_WINDOWS
        export TEST_WINE="$wine" TEST_GIT_EXE="$GIT_FOR_WINDOWS/cmd/git.exe"
        export WINEPREFIX="$fixture/wine-prefix" WINEDEBUG=-all
        mkdir -p "$FAKE_GIT_ROOT/cmd" bin-wine
        ln -s "$GIT_FOR_WINDOWS/usr" "$FAKE_GIT_ROOT/usr"
        cat > bin-wine/cygpath <<'EOF'
#!/bin/bash
# Bridge to the real cygpath.exe: /z is the host root in the msys namespace,
# so translate in and out and let the real binary apply the requested mode.
mode=
case "${1:-}" in
    -*) mode=$1; shift ;;
esac
out=$("$TEST_WINE" "$GIT_FOR_WINDOWS/usr/bin/cygpath.exe" ${mode:+$mode} "/z$1")
rc=$?
out=${out%$'\r'}
[[ "$out" == /z/* ]] && out=${out#/z}
printf '%s\n' "$out"
exit "$rc"
EOF
        chmod +x bin-wine/cygpath
        cat > "$FAKE_GIT_ROOT/cmd/git" <<'EOF'
#!/bin/bash
set -e
# Bridge namespaces for the Windows binary: cd for -C, and rewrite the -F
# config (and the absolute paths inside it) into Z:-drive form. The ssh
# program token is left exactly as the publish script produced it.
if [[ "${1:-}" == "-C" ]]; then cd "$2"; shift 2; fi
git_args=("$@")
eval "set -- ${GIT_SSH_COMMAND:-}"
prog=$1
cfg=
while (($#)); do
    [[ "$1" == "-F" ]] && { cfg=${2:-}; break; }
    shift
done
if [[ -n "$cfg" && "$cfg" == /* ]]; then
    sed 's|"/|"Z:/|g' "$cfg" > "$cfg.win"
    printf -v GIT_SSH_COMMAND '%q -F %q' "$prog" "Z:$cfg.win"
    export GIT_SSH_COMMAND
fi
exec "$TEST_WINE" "$TEST_GIT_EXE" "${git_args[@]}"
EOF
        chmod +x "$FAKE_GIT_ROOT/cmd/git"
        git clone -q --bare seed output2.git
        printf '[127.0.0.1]:%s %s\n' "$wport" "$(cat host_key.pub)" > github_known_hosts
        export PATH="$FAKE_GIT_ROOT/cmd:$fixture/bin-wine:$fixture/bin:$fixture:$PATH"
        export PUBLISH_OUTPUT_REPO="ssh://$(id -un)@127.0.0.1:$wport$fixture/output2.git"
        ./publish_benchmark_output.sh
        git --git-dir=output2.git show HEAD:markdown/page.md | cmp markdown/page.md -
        echo 'Published the expected artifact through Git for Windows under Wine.'
    else
        echo 'Could not start a local sshd listener; skipped the Git-for-Windows leg.' >&2
    fi
else
    echo 'Git for Windows under Wine unavailable; skipped that leg (set GIT_FOR_WINDOWS to a Git-for-Windows install and install Wine to run it).'
fi
