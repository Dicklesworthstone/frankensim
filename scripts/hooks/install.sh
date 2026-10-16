#!/bin/sh
# Install the FrankenSim pre-commit guard.
#
#   scripts/hooks/install.sh            install (preserves any existing hook)
#   scripts/hooks/install.sh --status   report what is installed
#
# Installation is explicit rather than automatic. `.git/hooks` is not tracked,
# it is the developer's own machine state, and silently writing an executable
# into someone's git plumbing is not something a repository should do on its
# own.
#
# An existing pre-commit hook is never clobbered: it is moved to
# `pre-commit.local`, which this guard chains before running its own checks. A
# non-zero exit from the chained hook still aborts the commit, so Agent Mail's
# reservation guard keeps working exactly as before.

set -eu

root=$(git rev-parse --show-toplevel)
hooks="$root/.git/hooks"
source_hook="$root/scripts/hooks/pre-commit"
target="$hooks/pre-commit"

if [ "${1:-}" = "--status" ]; then
    if [ ! -e "$target" ]; then
        echo "pre-commit: NOT installed"
    elif cmp -s "$source_hook" "$target"; then
        echo "pre-commit: installed and current"
    else
        echo "pre-commit: installed but DIFFERS from scripts/hooks/pre-commit"
    fi
    [ -e "$hooks/pre-commit.local" ] && echo "pre-commit.local: present (chained first)"
    exit 0
fi

mkdir -p "$hooks"

if [ -e "$target" ] && ! cmp -s "$source_hook" "$target"; then
    if [ -e "$hooks/pre-commit.local" ]; then
        echo "refusing to install: both pre-commit and pre-commit.local already exist." >&2
        echo "Resolve by hand so no one's hook is lost." >&2
        exit 1
    fi
    mv "$target" "$hooks/pre-commit.local"
    echo "preserved the existing hook as pre-commit.local (it will be chained first)"
fi

cp "$source_hook" "$target"
chmod +x "$target"
echo "installed $target"
echo "warns by default; set FRANKENSIM_HOOK_STRICT=1 to make warnings refusals"
