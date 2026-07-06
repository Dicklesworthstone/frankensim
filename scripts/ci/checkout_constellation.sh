#!/usr/bin/env bash
# Check out the Franken constellation siblings at the exact pins recorded in
# constellation.lock, laid out adjacent to the frankensim checkout — the same
# relative layout the workspace's path dependencies and `xtask
# check-constellation` expect:
#
#   <parent>/frankensim        (this repo)
#   <parent>/asupersync        @ pinned git_head
#   <parent>/frankensqlite     @ pinned git_head
#   ...
#
# Usage: scripts/ci/checkout_constellation.sh [parent-dir]
#   parent-dir defaults to the parent of this repository's root.
#
# Each sibling is fetched shallowly at its pinned SHA (depth 1), so CI cost
# stays proportional to tree size, not history. Existing checkouts are reused
# when already at the pin; a mismatched existing checkout is a hard error
# (never silently rebuilt against drifted sources — Decalogue P9).
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
parent="${1:-$(dirname "$repo_root")}"
lock="$repo_root/constellation.lock"

if [[ ! -f "$lock" ]]; then
    echo "FATAL: constellation.lock not found at $lock" >&2
    exit 1
fi

# lock JSON -> "lib head" lines (python3 is present on all CI runners).
entries="$(python3 - "$lock" <<'PY'
import json, sys
lock = json.load(open(sys.argv[1]))
for lib in lock["libraries"]:
    print(lib["lib"], lib["git_head"])
PY
)"

status=0
while read -r lib head; do
    dir="$parent/$lib"
    if [[ -d "$dir/.git" ]]; then
        have="$(git -C "$dir" rev-parse HEAD)"
        if [[ "$have" == "$head" ]]; then
            echo "{\"constellation\":\"$lib\",\"verdict\":\"reused\",\"head\":\"$head\"}"
            continue
        fi
        echo "{\"constellation\":\"$lib\",\"verdict\":\"drift\",\"want\":\"$head\",\"have\":\"$have\"}"
        echo "FATAL: $dir exists at $have but the lock pins $head" >&2
        status=1
        continue
    fi
    mkdir -p "$dir"
    git -C "$dir" init --quiet
    git -C "$dir" remote add origin "https://github.com/Dicklesworthstone/$lib.git"
    git -C "$dir" fetch --quiet --depth 1 origin "$head"
    git -C "$dir" checkout --quiet --detach FETCH_HEAD
    echo "{\"constellation\":\"$lib\",\"verdict\":\"cloned\",\"head\":\"$head\"}"
done <<<"$entries"

exit "$status"
