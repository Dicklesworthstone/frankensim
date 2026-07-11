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
# Usage:
#   scripts/ci/checkout_constellation.sh [parent-dir]
#   scripts/ci/checkout_constellation.sh --verify-only [parent-dir]
#   scripts/ci/checkout_constellation.sh --snapshot [parent-dir]
#
# parent-dir defaults to the parent of this repository's root. Snapshot mode
# prints one SHA-256 identity covering this checkout's HEAD + tracked/untracked
# changes, constellation.lock, and every clean pinned sibling tree.
#
# Each sibling is fetched shallowly at its pinned SHA (depth 1), so CI cost
# stays proportional to tree size, not history. Existing checkouts are reused
# when already at the pin; a mismatched existing checkout is a hard error
# (never silently rebuilt against drifted sources — Decalogue P9).
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mode="checkout"
case "${1:-}" in
    --verify-only|--snapshot)
        mode="${1#--}"
        shift
        ;;
esac
if [[ $# -gt 1 ]]; then
    echo "FATAL: usage: ${BASH_SOURCE[0]} [--verify-only|--snapshot] [parent-dir]" >&2
    exit 2
fi
parent="${1:-$(dirname "$repo_root")}"
lock="$repo_root/constellation.lock"

if [[ ! -f "$lock" ]]; then
    echo "FATAL: constellation.lock not found at $lock" >&2
    exit 1
fi

# Lock JSON -> validated tab-separated rows. Tabs/newlines are forbidden so the
# shell loop cannot reinterpret a repository identity.
if ! entries="$(python3 - "$lock" <<'PY'
import json
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
if path.stat().st_size > 1_048_576:
    raise SystemExit("constellation lock exceeds the 1 MiB parser bound")

def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result

document = json.loads(path.read_text(), object_pairs_hook=unique_object)
expected_top_keys = {"schema", "lock_hash", "note", "libraries"}
if set(document) != expected_top_keys:
    raise SystemExit("constellation lock has missing or unknown top-level fields")
if document.get("schema") != "frankensim-constellation-lock-v2":
    raise SystemExit("unsupported constellation lock schema")
note = (
    "lock_hash covers (lib, version, git_head) only — paths are per-machine; "
    "remote is transport for bootstrap-constellation (content identity is the git head)"
)
if document.get("note") != note:
    raise SystemExit("constellation lock identity note is not canonical")
lock_hash = document.get("lock_hash")
if not isinstance(lock_hash, str) or not re.fullmatch(r"[0-9a-f]{16}", lock_hash):
    raise SystemExit("constellation lock hash is not canonical lowercase hex")
rows = document.get("libraries")
if not isinstance(rows, list) or not rows:
    raise SystemExit("constellation lock has no libraries")
expected = {
    "asupersync", "franken_networkx", "franken_numpy", "frankenpandas",
    "frankenscipy", "frankensqlite", "frankentorch",
}
seen = set()
for row in rows:
    if not isinstance(row, dict) or set(row) != {"lib", "version", "git_head", "remote", "path"}:
        raise SystemExit("constellation lock row has missing or unknown fields")
    values = [row["lib"], row["version"], row["git_head"], row["remote"], row["path"]]
    if not all(isinstance(value, str) and value for value in values):
        raise SystemExit("constellation lock row has a missing/non-string field")
    lib, version, head, remote, _local_path = values
    if lib in seen:
        raise SystemExit(f"duplicate constellation library: {lib}")
    seen.add(lib)
    if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", head):
        raise SystemExit(f"invalid git head for {lib}")
    if any(any(ord(character) < 32 or ord(character) == 127 for character in value) for value in values):
        raise SystemExit(f"control character in lock row for {lib}")
    print("\t".join((lib, version, head, remote)))
if seen != expected:
    missing = sorted(expected - seen)
    extra = sorted(seen - expected)
    raise SystemExit(f"constellation library set mismatch: missing={missing}, extra={extra}")

identity = "".join(
    f'{row["lib"]}={row["version"]}@{row["git_head"]}\n'
    for row in sorted(rows, key=lambda candidate: candidate["lib"])
)
value = 0xcbf29ce484222325
for byte in identity.encode():
    value ^= byte
    value = (value * 0x100000001b3) & 0xffffffffffffffff
if f"{value:016x}" != lock_hash:
    raise SystemExit("constellation lock hash disagrees with its declared rows")
PY
)"; then
    echo "FATAL: could not parse $lock" >&2
    exit 1
fi
if [[ -z "$entries" ]]; then
    echo "FATAL: constellation lock produced no entries" >&2
    exit 1
fi

json_row() { # library verdict expected actual detail
    python3 - "$1" "$2" "$3" "$4" "$5" <<'PY'
import json
import sys

library, verdict, expected, actual, detail = sys.argv[1:]
print(json.dumps({
    "check": "constellation-checkout",
    "constellation": library,
    "verdict": verdict,
    "expected_head": expected,
    "actual_head": actual,
    "detail": detail,
}, separators=(",", ":")))
PY
}

working_tree_status() {
    local dir="$1" tracked untracked index_flags hidden_index
    tracked="$(git -C "$dir" -c core.fileMode=true -c core.excludesFile=/dev/null \
        status --porcelain --untracked-files=all)" || return 1
    untracked="$(git -C "$dir" -c core.excludesFile=/dev/null \
        ls-files --others --exclude-per-directory=.gitignore)" || return 1
    index_flags="$(git -C "$dir" ls-files -v)" || return 1
    hidden_index="$(printf '%s\n' "$index_flags" \
        | LC_ALL=C awk 'substr($0, 1, 1) == "S" || substr($0, 1, 1) ~ /[a-z]/ { print; exit }')"
    printf '%s%s%s' "$tracked" "$untracked" \
        "${hidden_index:+index flag hides worktree state: $hidden_index}"
}

verify_clean() { # directory expected-head
    local dir="$1" expected="$2" actual status confirmed
    if ! actual="$(git -C "$dir" rev-parse HEAD 2>/dev/null)"; then
        echo "FATAL: $dir is not a readable git checkout" >&2
        return 1
    fi
    if [[ "$actual" != "$expected" ]]; then
        echo "FATAL: $dir exists at $actual but the lock pins $expected" >&2
        return 1
    fi
    if ! status="$(working_tree_status "$dir")"; then
        echo "FATAL: could not inspect working-tree status for $dir" >&2
        return 1
    fi
    if [[ -n "$status" ]]; then
        echo "FATAL: $dir is dirty at pinned head $expected" >&2
        return 1
    fi
    if ! confirmed="$(git -C "$dir" rev-parse HEAD 2>/dev/null)" \
        || [[ "$confirmed" != "$expected" ]]; then
        echo "FATAL: $dir moved while its pinned state was being verified" >&2
        return 1
    fi
}

snapshot_identity() {
    python3 - "$repo_root" "$parent" "$lock" <<'PY'
import hashlib
import json
import os
import stat
import subprocess
import sys

root, parent, lock_path = sys.argv[1:]
digest = hashlib.sha256()

def add(label, data):
    if isinstance(label, str):
        label = label.encode()
    if isinstance(data, str):
        data = data.encode()
    digest.update(len(label).to_bytes(8, "big"))
    digest.update(label)
    digest.update(len(data).to_bytes(8, "big"))
    digest.update(data)

def git(repo, *args):
    completed = subprocess.run(
        ["git", "-C", repo, *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        detail = completed.stderr.decode(errors="replace").strip()
        raise RuntimeError(f"git {args!r} in {repo} failed: {detail}")
    return completed.stdout

def file_digest(path):
    hashed = hashlib.sha256()
    length = 0
    with open(path, "rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            length += len(chunk)
            hashed.update(chunk)
    return length, hashed.digest()

try:
    add("schema", "frankensim-ci-content-snapshot-v2")
    add("root-head", git(root, "rev-parse", "HEAD").strip())
    index = git(root, "ls-files", "--stage", "-z")
    add("root-index", index)
    listed = git(
        root,
        "-c",
        "core.excludesFile=/dev/null",
        "ls-files",
        "-z",
        "--cached",
        "--others",
        "--exclude-per-directory=.gitignore",
    )
    root_bytes = os.fsencode(root)
    for relative in sorted(set(filter(None, listed.split(b"\0")))):
        add("root-path", relative)
        absolute = os.path.join(root_bytes, relative)
        try:
            metadata = os.lstat(absolute)
        except FileNotFoundError:
            add("root-entry-kind", "missing")
            continue
        if stat.S_ISREG(metadata.st_mode):
            mode = "100755" if metadata.st_mode & stat.S_IXUSR else "100644"
            length, content = file_digest(absolute)
            add("root-entry-kind", "file")
            add("root-entry-mode", mode)
            add("root-entry-length", length.to_bytes(8, "big"))
            add("root-entry-sha256", content)
        elif stat.S_ISLNK(metadata.st_mode):
            add("root-entry-kind", "symlink")
            add("root-entry-mode", "120000")
            add("root-entry-target", os.fsencode(os.readlink(absolute)))
        elif stat.S_ISDIR(metadata.st_mode):
            add("root-entry-kind", "git-directory")
            add("root-entry-head", git(absolute, "rev-parse", "HEAD").strip())
            add(
                "root-entry-status",
                git(
                    absolute,
                    "-c",
                    "core.fileMode=true",
                    "-c",
                    "core.excludesFile=/dev/null",
                    "status",
                    "--porcelain",
                    "--untracked-files=all",
                ),
            )
        else:
            raise RuntimeError(
                f"unsupported working-tree entry type: {os.fsdecode(relative)!r}"
            )

    lock_bytes = open(lock_path, "rb").read()
    add("constellation-lock", lock_bytes)
    document = json.loads(lock_bytes)
    for row in sorted(document["libraries"], key=lambda candidate: candidate["lib"]):
        sibling = os.path.join(parent, row["lib"])
        actual = git(sibling, "rev-parse", "HEAD").decode().strip()
        if actual != row["git_head"]:
            raise RuntimeError(
                f"{sibling} is at {actual}, lock pins {row['git_head']}"
            )
        tracked_status = git(
            sibling,
            "-c",
            "core.fileMode=true",
            "-c",
            "core.excludesFile=/dev/null",
            "status",
            "--porcelain",
            "--untracked-files=all",
            "--ignore-submodules=none",
        )
        untracked = git(
            sibling,
            "-c",
            "core.excludesFile=/dev/null",
            "ls-files",
            "--others",
            "--exclude-per-directory=.gitignore",
        )
        index_flags = git(sibling, "ls-files", "-v")
        hidden_index = next(
            (
                line
                for line in index_flags.splitlines()
                if line[:1] == b"S" or line[:1].islower()
            ),
            None,
        )
        if tracked_status or untracked or hidden_index:
            raise RuntimeError(f"{sibling} is dirty at pinned head {actual}")
        confirmed = git(sibling, "rev-parse", "HEAD").decode().strip()
        if confirmed != actual:
            raise RuntimeError(
                f"{sibling} moved during snapshot: before={actual}, after={confirmed}"
            )
        add("sibling-lib", row["lib"])
        add("sibling-version", row["version"])
        add("sibling-expected-head", row["git_head"])
        add("sibling-actual-head", actual)
        add(
            "sibling-tree",
            git(sibling, "rev-parse", f"{row['git_head']}^{{tree}}").strip(),
        )
        add("sibling-remote", row["remote"])
except (OSError, RuntimeError, KeyError, TypeError, json.JSONDecodeError) as error:
    print(f"FATAL: cannot establish content snapshot: {error}", file=sys.stderr)
    raise SystemExit(1)

print(digest.hexdigest())
PY
}

if [[ "$mode" == "snapshot" ]]; then
    if snapshot_identity; then
        exit 0
    fi
    exit 1
fi

bootstrap_marker_key="frankensim.bootstrapIncomplete"

is_git_checkout() {
    local dir="$1" top dir_real top_real
    [[ -d "$dir" && ! -L "$dir" ]] || return 1
    top="$(git -C "$dir" rev-parse --show-toplevel 2>/dev/null)" || return 1
    dir_real="$(cd "$dir" && pwd -P)" || return 1
    top_real="$(cd "$top" && pwd -P)" || return 1
    [[ "$dir_real" == "$top_real" ]]
}

has_bootstrap_marker() {
    [[ "$(git -C "$1" config --local --get "$bootstrap_marker_key" 2>/dev/null || true)" == "true" ]]
}

materialize_checkout() { # directory expected-head remote
    local dir="$1" expected="$2" remote="$3" status origin
    if [[ "$remote" == "no-remote" ]]; then
        echo "FATAL: lock declares no remote for $dir" >&2
        return 1
    fi
    if ! status="$(working_tree_status "$dir")"; then
        echo "FATAL: could not inspect incomplete checkout $dir" >&2
        return 1
    fi
    if [[ -n "$status" ]]; then
        echo "FATAL: incomplete checkout $dir has worktree changes; refusing to overwrite it" >&2
        return 1
    fi
    git -C "$dir" config --local "$bootstrap_marker_key" true || return 1
    if origin="$(git -C "$dir" remote get-url origin 2>/dev/null)"; then
        if [[ "$origin" != "$remote" ]]; then
            echo "FATAL: incomplete checkout $dir origin is $origin, expected $remote" >&2
            return 1
        fi
    else
        git -C "$dir" remote add origin "$remote" || return 1
    fi
    git -C "$dir" fetch --quiet --depth 1 origin "$expected" || return 1
    git -C "$dir" checkout --quiet --detach "$expected" || return 1
    verify_clean "$dir" "$expected" || return 1
    if ! git -C "$dir" config --local --unset-all "$bootstrap_marker_key" >/dev/null 2>&1; then
        echo "FATAL: verified checkout $dir retained its incomplete marker" >&2
        return 1
    fi
}

status=0
while IFS=$'\t' read -r lib _version head remote; do
    [[ -n "$lib" ]] || continue
    dir="$parent/$lib"
    if is_git_checkout "$dir"; then
        have="$(git -C "$dir" rev-parse HEAD 2>/dev/null || true)"
        if [[ "$have" == "$head" ]] && verify_clean "$dir" "$head"; then
            if has_bootstrap_marker "$dir" \
                && ! git -C "$dir" config --local --unset-all "$bootstrap_marker_key" >/dev/null 2>&1; then
                json_row "$lib" "refused" "$head" "$have" "verified checkout retained incomplete marker"
                echo "FATAL: verified checkout $dir retained its incomplete marker" >&2
                status=1
                continue
            fi
            json_row "$lib" "verified" "$head" "$have" "pinned head and clean tree"
            continue
        fi
        if [[ "$mode" == "verify-only" ]]; then
            json_row "$lib" "refused" "$head" "${have:-<unborn>}" "existing checkout is incomplete, drifted, dirty, or unreadable"
            echo "FATAL: required constellation sibling $dir is not a complete clean pin" >&2
            status=1
            continue
        fi
        if [[ -n "$have" ]] && ! has_bootstrap_marker "$dir"; then
            json_row "$lib" "refused" "$head" "$have" "ordinary existing checkout is at the wrong head"
            echo "FATAL: $dir exists at $have but the lock pins $head" >&2
            status=1
            continue
        fi
        if [[ -z "$have" ]] && ! has_bootstrap_marker "$dir"; then
            origin="$(git -C "$dir" remote get-url origin 2>/dev/null || true)"
            if [[ "$origin" != "$remote" ]]; then
                json_row "$lib" "refused" "$head" "<unborn>" \
                    "unmarked unborn checkout does not have the exact locked origin"
                echo "FATAL: refusing to adopt unmarked unborn checkout $dir without locked origin $remote" >&2
                status=1
                continue
            fi
        fi
        if materialize_checkout "$dir" "$head" "$remote"; then
            json_row "$lib" "resumed" "$head" "$head" "resumed incomplete checkout and verified clean pin"
        else
            json_row "$lib" "refused" "$head" "${have:-<unborn>}" "incomplete checkout could not be resumed safely"
            status=1
        fi
        continue
    fi
    if [[ "$mode" == "verify-only" ]]; then
        json_row "$lib" "refused" "$head" "<missing>" "required sibling checkout is missing"
        echo "FATAL: required constellation sibling $dir is missing" >&2
        status=1
        continue
    fi
    if [[ -L "$dir" ]] \
        || [[ -e "$dir" && ! -d "$dir" ]] \
        || { [[ -d "$dir" ]] && [[ -n "$(find "$dir" ! -path "$dir" -print -quit 2>/dev/null)" ]]; }; then
        json_row "$lib" "refused" "$head" "<non-git>" "non-empty path is not a git checkout"
        echo "FATAL: refusing to initialize non-empty non-git directory $dir" >&2
        status=1
        continue
    fi
    mkdir -p "$dir"
    if ! git -C "$dir" init --quiet \
        || ! git -C "$dir" config --local "$bootstrap_marker_key" true \
        || ! materialize_checkout "$dir" "$head" "$remote"; then
        json_row "$lib" "refused" "$head" "<fetch-failed>" "clone/fetch/checkout verification failed"
        echo "FATAL: could not materialize $lib at $head from $remote" >&2
        status=1
        continue
    fi
    json_row "$lib" "cloned" "$head" "$head" "fetched pinned head and verified clean tree"
done <<<"$entries"

if [[ "$status" -ne 0 ]]; then
    exit 1
fi
exit 0
