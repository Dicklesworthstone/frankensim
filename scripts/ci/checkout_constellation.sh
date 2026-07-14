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

git_at() { # repository git-arguments...
    local repository="$1"
    shift
    env \
        -u GIT_DIR \
        -u GIT_WORK_TREE \
        -u GIT_INDEX_FILE \
        -u GIT_COMMON_DIR \
        -u GIT_OBJECT_DIRECTORY \
        -u GIT_ALTERNATE_OBJECT_DIRECTORIES \
        -u GIT_NAMESPACE \
        -u GIT_PREFIX \
        -u GIT_REPLACE_REF_BASE \
        -u GIT_SHALLOW_FILE \
        -u GIT_CEILING_DIRECTORIES \
        -u GIT_DISCOVERY_ACROSS_FILESYSTEM \
        -u GIT_CONFIG \
        -u GIT_CONFIG_COUNT \
        -u GIT_CONFIG_PARAMETERS \
        -u GIT_CONFIG_GLOBAL \
        -u GIT_CONFIG_SYSTEM \
        -u GIT_CONFIG_NOSYSTEM \
        -u GIT_ATTR_NOSYSTEM \
        -u GIT_ATTR_SOURCE \
        -u GIT_TEMPLATE_DIR \
        -u GIT_DEFAULT_HASH \
        -u GIT_DEFAULT_REF_FORMAT \
        -u GIT_REFERENCE_BACKEND \
        -u GIT_EXEC_PATH \
        -u GIT_EXTERNAL_DIFF \
        -u GIT_ASKPASS \
        -u SSH_ASKPASS \
        -u GIT_SSH \
        -u GIT_SSH_COMMAND \
        -u GIT_PROXY_COMMAND \
        -u GIT_ALLOW_PROTOCOL \
        -u GIT_PROTOCOL_FROM_USER \
        -u GIT_TERMINAL_PROMPT \
        -u GIT_QUARANTINE_PATH \
        -u GIT_OPTIONAL_LOCKS \
        -u GIT_REDIRECT_STDIN \
        -u GIT_REDIRECT_STDOUT \
        -u GIT_REDIRECT_STDERR \
        -u GIT_LITERAL_PATHSPECS \
        -u GIT_GLOB_PATHSPECS \
        -u GIT_NOGLOB_PATHSPECS \
        -u GIT_ICASE_PATHSPECS \
        -u GIT_TRACE \
        -u GIT_TRACE_CURL \
        -u GIT_TRACE_CURL_NO_DATA \
        -u GIT_TRACE_FSMONITOR \
        -u GIT_TRACE_PACK_ACCESS \
        -u GIT_TRACE_PACKET \
        -u GIT_TRACE_PACKFILE \
        -u GIT_TRACE_PERFORMANCE \
        -u GIT_TRACE_REFS \
        -u GIT_TRACE_SETUP \
        -u GIT_TRACE_SHALLOW \
        -u GIT_TRACE2 \
        -u GIT_TRACE2_EVENT \
        -u GIT_TRACE2_PERF \
        -u GIT_NO_LAZY_FETCH \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_NO_LAZY_FETCH=1 \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_ATTR_NOSYSTEM=1 \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_TERMINAL_PROMPT=0 \
        git \
            -c core.hooksPath=/dev/null \
            -c core.attributesFile=/dev/null \
            -c credential.helper= \
            -c protocol.allow=never \
            -c protocol.file.allow=always \
            -c protocol.https.allow=always \
            -c protocol.ssh.allow=always \
            -C "$repository" "$@"
}

if [[ ! -f "$lock" ]]; then
    echo "FATAL: constellation.lock not found at $lock" >&2
    exit 1
fi

# Lock JSON -> validated tab-separated rows. Tabs/newlines are forbidden so the
# shell loop cannot reinterpret a repository identity.
if ! validated_lock="$(python3 -I - "$lock" <<'PY'
import hashlib
import json
import pathlib
import re
import sys
import unicodedata

path = pathlib.Path(sys.argv[1])
with path.open("rb") as source:
    raw = source.read(1_048_577)
if len(raw) > 1_048_576:
    raise SystemExit("constellation lock exceeds the 1 MiB parser bound")

def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result

text = raw.decode("utf-8")
document = json.loads(text, object_pairs_hook=unique_object)
expected_top_keys = {
    "schema", "identity_domain", "identity_version", "lock_hash", "note", "libraries",
}
if set(document) != expected_top_keys:
    raise SystemExit("constellation lock has missing or unknown top-level fields")
if document.get("schema") != "frankensim-constellation-lock-v2":
    raise SystemExit("unsupported constellation lock schema")
identity_domain = document.get("identity_domain")
if not isinstance(identity_domain, str):
    raise SystemExit("constellation lock identity domain must be a string")
if identity_domain != "org.frankensim.xtask.constellation-lock.v1":
    raise SystemExit("unsupported constellation lock identity domain")
identity_version = document.get("identity_version")
# bool is an int subclass in Python; exact type admission keeps true from
# impersonating identity version 1.
if type(identity_version) is not int:
    raise SystemExit("constellation lock identity version must be an integer")
if identity_version != 1:
    raise SystemExit("unsupported constellation lock identity version")
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
entry_lines = []
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
    if any(
        any(unicodedata.category(character) in {"Cc", "Cs"} for character in value)
        for value in values
    ):
        raise SystemExit(f"control character in lock row for {lib}")
    entry_lines.append("\t".join((lib, version, head, remote)))
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

# The Rust consumers accept one exact byte grammar, not merely an equivalent
# JSON object. Re-render with the same key order, spacing, escaping, and
# caller-supplied library order so whitespace or object-key rearrangement
# cannot create a shell-only accepted lock.
def quoted(value):
    return json.dumps(value, ensure_ascii=False)

canonical_rows = []
for row in rows:
    canonical_rows.append(
        "    {\"lib\": " + quoted(row["lib"])
        + ", \"version\": " + quoted(row["version"])
        + ", \"git_head\": " + quoted(row["git_head"])
        + ", \"remote\": " + quoted(row["remote"])
        + ", \"path\": " + quoted(row["path"])
        + "}"
    )
canonical = (
    "{\n"
    + "  \"schema\": " + quoted(document["schema"]) + ",\n"
    + "  \"identity_domain\": " + quoted(identity_domain) + ",\n"
    + f"  \"identity_version\": {identity_version},\n"
    + "  \"lock_hash\": " + quoted(lock_hash) + ",\n"
    + "  \"note\": " + quoted(document["note"]) + ",\n"
    + "  \"libraries\": [\n"
    + ",\n".join(canonical_rows)
    + "\n  ]\n}\n"
)
if text != canonical:
    raise SystemExit("constellation lock is valid JSON but not canonical")
print(f"@lock-sha256\t{hashlib.sha256(raw).hexdigest()}")
print("\n".join(entry_lines))
PY
)"; then
    echo "FATAL: could not parse $lock" >&2
    exit 1
fi
lock_header="${validated_lock%%$'\n'*}"
entries="${validated_lock#*$'\n'}"
validated_lock_sha256="${lock_header#@lock-sha256$'\t'}"
if [[ "$lock_header" == "$validated_lock" ]] \
    || [[ "$lock_header" != @lock-sha256$'\t'* ]] \
    || [[ ! "$validated_lock_sha256" =~ ^[0-9a-f]{64}$ ]]; then
    echo "FATAL: constellation lock validator produced no canonical identity" >&2
    exit 1
fi
if [[ -z "$entries" ]]; then
    echo "FATAL: constellation lock produced no entries" >&2
    exit 1
fi

json_row() { # library verdict expected actual detail
    python3 -I - "$1" "$2" "$3" "$4" "$5" <<'PY'
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
    python3 -I - "$1" <<'PY'
import hashlib
import os
import stat
import subprocess
import sys

root = os.fsencode(sys.argv[1])
FILE_MODE_OBSERVABLE = os.name != "nt"

git_environment = os.environ.copy()
git_environment["GIT_OPTIONAL_LOCKS"] = "0"
git_environment["GIT_NO_REPLACE_OBJECTS"] = "1"
for inherited in (
    "GIT_DIR",
    "GIT_WORK_TREE",
    "GIT_INDEX_FILE",
    "GIT_COMMON_DIR",
    "GIT_OBJECT_DIRECTORY",
    "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    "GIT_NAMESPACE",
    "GIT_PREFIX",
    "GIT_REPLACE_REF_BASE",
    "GIT_SHALLOW_FILE",
    "GIT_CEILING_DIRECTORIES",
    "GIT_DISCOVERY_ACROSS_FILESYSTEM",
    "GIT_CONFIG",
    "GIT_CONFIG_COUNT",
    "GIT_CONFIG_PARAMETERS",
    "GIT_CONFIG_GLOBAL",
    "GIT_CONFIG_SYSTEM",
    "GIT_CONFIG_NOSYSTEM",
    "GIT_ATTR_NOSYSTEM",
    "GIT_ATTR_SOURCE",
    "GIT_TEMPLATE_DIR",
    "GIT_DEFAULT_HASH",
    "GIT_DEFAULT_REF_FORMAT",
    "GIT_REFERENCE_BACKEND",
    "GIT_EXEC_PATH",
    "GIT_EXTERNAL_DIFF",
    "GIT_ASKPASS",
    "SSH_ASKPASS",
    "GIT_SSH",
    "GIT_SSH_COMMAND",
    "GIT_PROXY_COMMAND",
    "GIT_ALLOW_PROTOCOL",
    "GIT_PROTOCOL_FROM_USER",
    "GIT_TERMINAL_PROMPT",
    "GIT_QUARANTINE_PATH",
    "GIT_REDIRECT_STDIN",
    "GIT_REDIRECT_STDOUT",
    "GIT_REDIRECT_STDERR",
    "GIT_LITERAL_PATHSPECS",
    "GIT_GLOB_PATHSPECS",
    "GIT_NOGLOB_PATHSPECS",
    "GIT_ICASE_PATHSPECS",
    "GIT_TRACE",
    "GIT_TRACE_CURL",
    "GIT_TRACE_CURL_NO_DATA",
    "GIT_TRACE_FSMONITOR",
    "GIT_TRACE_PACK_ACCESS",
    "GIT_TRACE_PACKET",
    "GIT_TRACE_PACKFILE",
    "GIT_TRACE_PERFORMANCE",
    "GIT_TRACE_REFS",
    "GIT_TRACE_SETUP",
    "GIT_TRACE_SHALLOW",
    "GIT_TRACE2",
    "GIT_TRACE2_EVENT",
    "GIT_TRACE2_PERF",
    "GIT_NO_LAZY_FETCH",
):
    git_environment.pop(inherited, None)
git_environment["GIT_CONFIG_NOSYSTEM"] = "1"
git_environment["GIT_ATTR_NOSYSTEM"] = "1"
git_environment["GIT_CONFIG_GLOBAL"] = os.devnull
git_environment["GIT_TERMINAL_PROMPT"] = "0"
git_environment["GIT_NO_LAZY_FETCH"] = "1"

def git_argv(repo, *args):
    return [
        "git",
        "-c", f"core.hooksPath={os.devnull}",
        "-c", f"core.attributesFile={os.devnull}",
        "-c", "credential.helper=",
        "-c", "protocol.allow=never",
        "-c", "protocol.file.allow=always",
        "-c", "protocol.https.allow=always",
        "-c", "protocol.ssh.allow=always",
        "-C", repo,
        *args,
    ]

def run_git(repo, *args):
    return subprocess.run(
        git_argv(repo, *args),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=git_environment,
    )

def git(repo, *args):
    completed = run_git(repo, *args)
    if completed.returncode != 0:
        detail = completed.stderr.decode(errors="backslashreplace").strip()
        raise RuntimeError(f"git {args!r} in {repo!r} failed: {detail}")
    return completed.stdout

def git_stdin(repo, data, *args):
    completed = subprocess.run(
        git_argv(repo, *args),
        input=data,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=git_environment,
    )
    if completed.returncode != 0:
        detail = completed.stderr.decode(errors="backslashreplace").strip()
        raise RuntimeError(f"git {args!r} in {repo!r} failed: {detail}")
    return completed.stdout

def git_line(repo, *args):
    output = git(repo, *args)
    if not output.endswith(b"\n"):
        raise RuntimeError(f"git {args!r} in {repo!r} returned no terminal newline")
    return output[:-1]

def optional_head(repo):
    completed = run_git(repo, "rev-parse", "--verify", "-q", "HEAD")
    if completed.returncode == 1 and not completed.stdout:
        return None
    if completed.returncode != 0:
        detail = completed.stderr.decode(errors="backslashreplace").strip()
        raise RuntimeError(f"cannot read HEAD in {repo!r}: {detail}")
    if not completed.stdout.endswith(b"\n"):
        raise RuntimeError(f"git rev-parse HEAD in {repo!r} returned no terminal newline")
    return completed.stdout[:-1]

def executable_config_key(key):
    key = key.lower()
    return (
        key.startswith((b"include.", b"includeif.", b"filter."))
        or key in {
            b"core.fsmonitor",
            b"core.hookspath",
            b"core.sshcommand",
            b"core.gitproxy",
            b"core.alternaterefscommand",
            b"core.askpass",
            b"credential.helper",
            b"extensions.refstorage",
            b"fetch.bundleuri",
            b"fetch.uriprotocols",
            b"gc.recentobjectshook",
            b"protocol.allow",
        }
        or (key.startswith(b"protocol.") and key.endswith(b".allow"))
        or (key.startswith(b"credential.") and key.endswith(b".helper"))
        or (
            key.startswith(b"remote.")
            and key.endswith((b".uploadpack", b".receivepack", b".vcs"))
        )
        or (key.startswith(b"submodule.") and key.endswith(b".update"))
        or (
            key.startswith(b"url.")
            and key.endswith((b".insteadof", b".pushinsteadof"))
        )
    )

def require_no_executable_config(repo):
    worktree_config = run_git(
        repo,
        "config",
        "--local",
        "--no-includes",
        "--type=bool",
        "--get-all",
        "extensions.worktreeConfig",
    )
    if (
        worktree_config.returncode == 1
        and not worktree_config.stdout
        and not worktree_config.stderr
    ):
        worktree_config_enabled = False
    elif worktree_config.returncode != 0:
        detail = worktree_config.stderr.decode(errors="backslashreplace").strip()
        raise RuntimeError(
            "cannot inspect extensions.worktreeConfig in "
            f"{repo!r}: {detail}"
        )
    elif worktree_config.stderr:
        detail = worktree_config.stderr.decode(errors="backslashreplace").strip()
        raise RuntimeError(
            "extensions.worktreeConfig inspection emitted an unexpected "
            f"diagnostic in {repo!r}: {detail}"
        )
    elif worktree_config.stdout == b"true\n":
        worktree_config_enabled = True
    elif worktree_config.stdout == b"false\n":
        worktree_config_enabled = False
    else:
        raise RuntimeError(
            "extensions.worktreeConfig must have exactly one canonical "
            "boolean value"
        )
    for config_scope in ("--local", "--worktree"):
        if config_scope == "--worktree" and not worktree_config_enabled:
            continue
        inventory = git(
            repo,
            "config",
            config_scope,
            "--null",
            "--name-only",
            "--no-includes",
            "--list",
        )
        if inventory and not inventory.endswith(b"\0"):
            raise RuntimeError(
                f"{config_scope} Git configuration inventory is not NUL-terminated"
            )
        for key in filter(None, inventory.split(b"\0")):
            if executable_config_key(key):
                raise RuntimeError(
                    f"{config_scope} Git configuration authority {key!r} may "
                    "execute code or redirect history; remove it before bootstrap "
                    "verification"
                )

def windows_reparse_point(metadata):
    if os.name != "nt":
        return False
    attributes = getattr(metadata, "st_file_attributes", None)
    if attributes is None:
        raise RuntimeError("Windows file metadata lacks reparse-point attributes")
    return bool(attributes & 0x400)

def git_marker(repo):
    marker = os.path.join(repo, b".git")
    try:
        metadata = os.lstat(marker)
    except FileNotFoundError:
        return None
    if windows_reparse_point(metadata):
        raise RuntimeError(f"unsupported reparse-point .git marker at {marker!r}")
    if stat.S_ISREG(metadata.st_mode):
        return b"file"
    if stat.S_ISDIR(metadata.st_mode):
        return b"directory"
    raise RuntimeError(f"unsupported .git marker type at {marker!r}")

def safe_repository_path(repo, require_marker=True):
    repo = os.fsencode(repo)
    try:
        metadata = os.lstat(repo)
    except FileNotFoundError as error:
        raise RuntimeError(f"repository root is absent at {repo!r}") from error
    if stat.S_ISLNK(metadata.st_mode) or windows_reparse_point(metadata):
        raise RuntimeError(f"repository root must not be a link or reparse point: {repo!r}")
    if not stat.S_ISDIR(metadata.st_mode):
        raise RuntimeError(f"repository root is not a directory: {repo!r}")
    repo = os.path.realpath(repo)
    if require_marker and git_marker(repo) is None:
        raise RuntimeError(f"repository has no .git marker at {repo!r}")
    return repo

def repository_boundary(repo, require_marker=True):
    repo = safe_repository_path(repo, require_marker)
    if git_line(repo, "rev-parse", "--is-inside-work-tree") != b"true":
        raise RuntimeError(f"repository is not a worktree at {repo!r}")
    prefix = git_line(repo, "rev-parse", "--show-prefix")
    if prefix:
        raise RuntimeError(
            f"repository path is not its exact worktree root: path={repo!r}, prefix={prefix!r}"
        )
    return repo

HFS_IGNORED_UTF8 = (
    b"\xe2\x80\x8c", b"\xe2\x80\x8d", b"\xe2\x80\x8e", b"\xe2\x80\x8f",
    b"\xe2\x80\xaa", b"\xe2\x80\xab", b"\xe2\x80\xac", b"\xe2\x80\xad",
    b"\xe2\x80\xae", b"\xe2\x81\xaa", b"\xe2\x81\xab", b"\xe2\x81\xac",
    b"\xe2\x81\xad", b"\xe2\x81\xae", b"\xe2\x81\xaf", b"\xef\xbb\xbf",
)

def windows_ntfs_fallback_alias(component, prefix):
    candidate = component.rstrip(b" .").lower()
    if len(candidate) != 8:
        return False
    saw_tilde = False
    for index, byte in enumerate(candidate):
        if saw_tilde:
            if not 0x30 <= byte <= 0x39:
                return False
        elif byte == 0x7E:
            if index + 1 >= len(candidate) or not 0x31 <= candidate[index + 1] <= 0x39:
                return False
            saw_tilde = True
        elif index >= 6 or byte != prefix[index]:
            return False
    return saw_tilde

def windows_gitmodules_alias(component):
    candidate = component.rstrip(b" .").lower()
    return (
        candidate == b".gitmodules"
        or (
            len(candidate) == 8
            and candidate[:7] == b"gitmod~"
            and candidate[7:8] in {b"1", b"2", b"3", b"4"}
        )
        or windows_ntfs_fallback_alias(candidate, b"gi7eba")
    )

def validate_git_worktree_path(relative, symlink=False):
    if not relative or b"\0" in relative:
        raise RuntimeError(f"git emitted an empty or NUL-bearing path: {relative!r}")
    components = relative.split(b"/")
    if relative.startswith(b"/") or relative.endswith(b"/") or any(
        component in {b"", b".", b".."} for component in components
    ):
        raise RuntimeError(f"git emitted an unsafe working-tree path: {relative!r}")
    for component in components:
        if component.lower() == b".git":
            raise RuntimeError(
                f"git emitted a path through repository metadata: {relative!r}"
            )
        if symlink and component.lower() == b".gitmodules":
            raise RuntimeError(
                f"git emitted a forbidden .gitmodules symlink: {relative!r}"
            )
        if any(sequence in component for sequence in HFS_IGNORED_UTF8):
            raise RuntimeError(
                f"git emitted an HFS-ambiguous Unicode path: {relative!r}"
            )
        if os.name == "nt":
            windows_alias = component.rstrip(b" .").lower()
            short_name = windows_alias[4:] if windows_alias.startswith(b"git~") else b""
            if windows_alias == b".git" or (
                1 <= len(short_name) <= 6 and short_name.isdigit()
            ):
                raise RuntimeError(
                    f"git emitted a Windows .git alias: {relative!r}"
                )
            if symlink and windows_gitmodules_alias(component):
                raise RuntimeError(
                    f"git emitted a Windows .gitmodules symlink alias: {relative!r}"
                )
    if os.name == "nt" and (b"\\" in relative or b":" in relative):
        raise RuntimeError(
            "git emitted a rooted, drive-qualified, backslash, or "
            f"alternate-stream path: {relative!r}"
        )
    return relative

def parse_index(index):
    parsed = {}
    for record in filter(None, index.split(b"\0")):
        if b"\t" not in record:
            raise RuntimeError("git ls-files --stage emitted a malformed record")
        metadata, relative = record.split(b"\t", 1)
        fields = metadata.split()
        if len(fields) != 3 or fields[2] not in {b"0", b"1", b"2", b"3"}:
            raise RuntimeError("git ls-files --stage emitted malformed metadata")
        mode, object_id, stage = fields
        validate_git_worktree_path(relative, symlink=mode == b"120000")
        parsed.setdefault(relative, []).append((stage, mode, object_id))
    return parsed

def parse_tagged_paths(output, command):
    parsed = {}
    for record in filter(None, output.split(b"\0")):
        if len(record) < 3 or record[1:2] != b" ":
            raise RuntimeError(f"{command} emitted a malformed tagged path")
        parsed[record[2:]] = record[:1]
    return parsed

def parse_status(output):
    records = []
    for record in filter(None, output.split(b"\0")):
        if len(record) < 4 or record[2:3] != b" ":
            raise RuntimeError("git status --porcelain=v1 emitted a malformed record")
        records.append((record[3:], record[:2]))
    return records

def parse_paths(output, command):
    paths = list(filter(None, output.split(b"\0")))
    if len(paths) != len(set(paths)):
        raise RuntimeError(f"{command} emitted duplicate paths")
    return paths

def safe_worktree_path(repo, relative):
    return os.path.join(repo, validate_git_worktree_path(relative))

MAX_HASH_BATCH_PATHS = 128
MAX_HASH_BATCH_BYTES = 65_536
MAX_GRAFT_BYTES = 1_048_576
MAX_REPOSITORY_DEPTH = 32

def metadata_identity(metadata):
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_mode,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )

primary_index_git_line = git_line
MAX_RAW_PRIMARY_INDEX_BYTES = 1_073_741_824

def primary_index_git_dir(repo):
    reported = primary_index_git_line(
        repo,
        "rev-parse",
        "--absolute-git-dir",
    )
    if not os.path.isabs(reported):
        raise RuntimeError(f"reported primary git dir is not absolute: {reported!r}")
    metadata = os.lstat(reported)
    if windows_reparse_point(metadata) or not stat.S_ISDIR(metadata.st_mode):
        raise RuntimeError(f"primary git dir is not an ordinary directory: {reported!r}")
    canonical = os.path.realpath(reported)
    if os.name != "nt":
        if reported != canonical:
            raise RuntimeError(
                f"reported primary git dir redirects from {reported!r} to {canonical!r}"
            )
        return reported, metadata_identity(metadata)
    canonical_metadata = os.lstat(canonical)
    if windows_reparse_point(canonical_metadata) \
            or not stat.S_ISDIR(canonical_metadata.st_mode) \
            or metadata_identity(canonical_metadata) != metadata_identity(metadata):
        raise RuntimeError(
            f"resolved primary git dir is not the reported ordinary directory: {canonical!r}"
        )
    return canonical, metadata_identity(canonical_metadata)

def primary_index_path(repo, git_dir):
    reported = primary_index_git_line(
        repo,
        "rev-parse",
        "--path-format=absolute",
        "--git-path",
        "index",
    )
    if not os.path.isabs(reported):
        raise RuntimeError(f"reported primary index path is not absolute: {reported!r}")
    if os.path.basename(reported) != b"index":
        raise RuntimeError(f"reported primary index basename is not 'index': {reported!r}")
    parent_authority = os.path.realpath(os.path.dirname(reported))
    if parent_authority != git_dir:
        raise RuntimeError(
            f"reported primary index parent {parent_authority!r} is not {git_dir!r}"
        )
    if os.name != "nt" and reported != os.path.join(git_dir, b"index"):
        raise RuntimeError(
            f"reported primary index is not a direct canonical child: {reported!r}"
        )
    return os.path.join(git_dir, b"index")

def raw_primary_index_state(repo):
    object_format = primary_index_git_line(
        repo,
        "rev-parse",
        "--show-object-format=storage",
    )
    if object_format == b"sha1":
        object_id_width = 20
        checksum = hashlib.sha1
    elif object_format == b"sha256":
        object_id_width = 32
        checksum = hashlib.sha256
    else:
        raise RuntimeError(
            f"unsupported primary-index object format: {object_format!r}"
        )
    git_dir, git_dir_identity = primary_index_git_dir(repo)
    index_path = primary_index_path(repo, git_dir)
    try:
        metadata = os.lstat(index_path)
    except FileNotFoundError:
        confirmed_object_format = primary_index_git_line(
            repo,
            "rev-parse",
            "--show-object-format=storage",
        )
        try:
            os.lstat(index_path)
        except FileNotFoundError:
            pass
        else:
            raise RuntimeError("primary index appeared while inspecting its absence")
        confirmed_git_dir, confirmed_git_dir_identity = primary_index_git_dir(repo)
        confirmed_index_path = primary_index_path(repo, confirmed_git_dir)
        if confirmed_object_format != object_format \
                or confirmed_git_dir != git_dir \
                or confirmed_git_dir_identity != git_dir_identity \
                or confirmed_index_path != index_path:
            raise RuntimeError("primary-index authority moved while inspecting its absence")
        return b"absent\0" + object_format, None
    if windows_reparse_point(metadata) or not stat.S_ISREG(metadata.st_mode):
        raise RuntimeError(f"primary index is not an ordinary file: {index_path!r}")
    if metadata.st_size > MAX_RAW_PRIMARY_INDEX_BYTES:
        raise RuntimeError("primary index exceeds the 1 GiB raw inspection bound")
    identity = metadata_identity(metadata)
    with open(index_path, "rb") as source:
        opened_metadata = os.fstat(source.fileno())
        if metadata_identity(opened_metadata) != identity:
            raise RuntimeError("primary index moved before raw inspection")
        raw_index = source.read(metadata.st_size)
        grew_during_read = source.read(1)
        confirmed_opened_metadata = os.fstat(source.fileno())
    confirmed_metadata = os.lstat(index_path)
    if windows_reparse_point(confirmed_metadata) \
            or metadata_identity(confirmed_opened_metadata) != identity \
            or metadata_identity(confirmed_metadata) != identity:
        raise RuntimeError("primary index moved during raw inspection")
    if len(raw_index) != metadata.st_size or grew_during_read:
        raise RuntimeError("primary index changed length during raw inspection")
    confirmed_object_format = primary_index_git_line(
        repo,
        "rev-parse",
        "--show-object-format=storage",
    )
    confirmed_git_dir, confirmed_git_dir_identity = primary_index_git_dir(repo)
    confirmed_index_path = primary_index_path(repo, confirmed_git_dir)
    if confirmed_object_format != object_format \
            or confirmed_git_dir != git_dir \
            or confirmed_git_dir_identity != git_dir_identity \
            or confirmed_index_path != index_path:
        raise RuntimeError("primary-index authority moved during raw inspection")

    if len(raw_index) < 12 + object_id_width:
        raise RuntimeError("primary index is truncated before its checksum")
    if raw_index[:4] != b"DIRC":
        raise RuntimeError("primary index has an invalid DIRC signature")
    version = int.from_bytes(raw_index[4:8], "big")
    if version == 4:
        raise RuntimeError("primary index version 4 path compression is unsupported")
    if version not in {2, 3}:
        raise RuntimeError(f"unsupported primary index version: {version}")
    entry_count = int.from_bytes(raw_index[8:12], "big")
    checksum_offset = len(raw_index) - object_id_width
    expected_checksum = raw_index[checksum_offset:]
    actual_checksum = checksum(raw_index[:checksum_offset]).digest()
    if actual_checksum != expected_checksum:
        raise RuntimeError(
            f"primary index {object_format.decode('ascii')} checksum is invalid"
        )

    cursor = 12
    entry_order = []
    fixed_entry_size = 40 + object_id_width + 2
    for entry_number in range(entry_count):
        entry_start = cursor
        fixed_entry_end = entry_start + fixed_entry_size
        if fixed_entry_end > checksum_offset:
            raise RuntimeError(
                f"primary index entry {entry_number} is truncated before its flags"
            )
        entry_mode = int.from_bytes(raw_index[entry_start + 24:entry_start + 28], "big")
        flags = int.from_bytes(raw_index[fixed_entry_end - 2:fixed_entry_end], "big")
        cursor = fixed_entry_end
        if flags & 0x4000:
            if version != 3:
                raise RuntimeError("primary index v2 entry has unsupported extended flags")
            if cursor + 2 > checksum_offset:
                raise RuntimeError(
                    f"primary index entry {entry_number} has truncated extended flags"
                )
            extended_flags = int.from_bytes(raw_index[cursor:cursor + 2], "big")
            if extended_flags & 0x9FFF:
                raise RuntimeError(
                    f"primary index entry {entry_number} has unsupported extended flags"
                )
            cursor += 2
        path_start = cursor
        encoded_name_length = flags & 0x0FFF
        if encoded_name_length < 0x0FFF:
            path_end = path_start + encoded_name_length
            if path_end >= checksum_offset:
                raise RuntimeError(
                    f"primary index entry {entry_number} is truncated in its path"
                )
            if b"\0" in raw_index[path_start:path_end]:
                raise RuntimeError(
                    f"primary index entry {entry_number} path contains an early NUL"
                )
        else:
            path_end = raw_index.find(b"\0", path_start, checksum_offset)
            if path_end < 0:
                raise RuntimeError(
                    f"primary index entry {entry_number} has no path terminator"
                )
            if path_end - path_start < 0x0FFF:
                raise RuntimeError(
                    f"primary index entry {entry_number} has an ambiguous long-path length"
                )
        stage = (flags >> 12) & 0x3
        path = raw_index[path_start:path_end]
        validate_git_worktree_path(
            path,
            symlink=entry_mode & 0o170000 == 0o120000,
        )
        entry_order.append((path, stage))
        padding_size = 8 - ((path_end - entry_start) % 8)
        entry_end = path_end + padding_size
        if entry_end > checksum_offset:
            raise RuntimeError(
                f"primary index entry {entry_number} has truncated path padding"
            )
        if raw_index[path_end:entry_end] != b"\0" * padding_size:
            raise RuntimeError(
                f"primary index entry {entry_number} has malformed path padding"
            )
        cursor = entry_end

    extensions = []
    while cursor < checksum_offset:
        if checksum_offset - cursor < 8:
            raise RuntimeError("primary index has a truncated extension header")
        signature = raw_index[cursor:cursor + 4]
        extension_size = int.from_bytes(raw_index[cursor + 4:cursor + 8], "big")
        extension_end = cursor + 8 + extension_size
        if extension_end > checksum_offset:
            raise RuntimeError(
                f"primary index extension {signature!r} is truncated"
            )
        if signature == b"FSMN":
            raise RuntimeError(
                "primary index contains forbidden FSMonitor FSMN extension"
            )
        if signature == b"link":
            raise RuntimeError("split-index link extension is unsupported")
        if signature == b"sdir":
            raise RuntimeError("sparse-index sdir extension is unsupported")
        if not 0x41 <= signature[0] <= 0x5A:
            raise RuntimeError(
                f"unsupported required primary-index extension: {signature!r}"
            )
        extensions.append(signature + extension_size.to_bytes(4, "big"))
        cursor = extension_end
    if cursor != checksum_offset:
        raise RuntimeError("primary index extension layout is ambiguous")
    previous_entry = None
    for path, stage in entry_order:
        current_entry = (path, stage)
        if not path:
            raise RuntimeError("primary index contains an empty entry path")
        if previous_entry is not None and current_entry == previous_entry:
            raise RuntimeError(
                f"primary index contains duplicate entry order: {current_entry!r}"
            )
        if previous_entry is not None and current_entry < previous_entry:
            raise RuntimeError(
                "primary index entries are not strictly sorted: "
                f"previous={previous_entry!r} current={current_entry!r}"
            )
        previous_entry = current_entry

    state = b"present\0" + object_format + b"\0" \
        + version.to_bytes(4, "big") \
        + entry_count.to_bytes(4, "big") \
        + b"".join(extensions)
    return state, entry_count

def bind_raw_primary_index(raw_state, index_inventory):
    state, raw_entry_count = raw_state
    if index_inventory:
        if not index_inventory.endswith(b"\0"):
            raise RuntimeError("git ls-files --stage inventory is not NUL-terminated")
        inventory_records = index_inventory[:-1].split(b"\0")
        if any(not record for record in inventory_records):
            raise RuntimeError("git ls-files --stage inventory contains an empty record")
        inventory_entry_count = len(inventory_records)
    else:
        inventory_entry_count = 0
    if raw_entry_count is None:
        if inventory_entry_count:
            raise RuntimeError(
                "primary index is absent but git ls-files reported index entries"
            )
    elif raw_entry_count != inventory_entry_count:
        raise RuntimeError(
            "primary-index entry count disagrees with git ls-files: "
            f"raw={raw_entry_count} inventory={inventory_entry_count}"
        )
    return state

def worktree_kind(metadata):
    if stat.S_ISREG(metadata.st_mode):
        return b"regular"
    if stat.S_ISLNK(metadata.st_mode):
        return b"symlink"
    if stat.S_ISDIR(metadata.st_mode):
        return b"directory"
    return b"special"

def index_materialization_snapshot(repo, paths):
    states = {}
    identity_owners = {}
    for path in paths:
        safe_worktree_path(repo, path)
        logical = b""
        materialized = repo
        components = path.split(b"/")
        for component_index, component in enumerate(components):
            parent = materialized
            logical = component if not logical else logical + b"/" + component
            materialized = os.path.join(materialized, component)
            try:
                metadata = os.lstat(materialized)
            except FileNotFoundError:
                break

            if os.name == "nt":
                try:
                    with os.scandir(parent) as entries:
                        exact_matches = sum(
                            os.fsencode(entry.name) == component for entry in entries
                        )
                except OSError as error:
                    raise RuntimeError(
                        f"cannot enumerate tracked index prefix parent {parent!r}: {error}"
                    ) from error
                if exact_matches != 1:
                    raise RuntimeError(
                        f"tracked index prefix {logical!r} is not materialized with "
                        "one exact directory-entry spelling"
                    )

            is_final = component_index + 1 == len(components)
            is_reparse_point = windows_reparse_point(metadata)
            if (is_reparse_point and not stat.S_ISLNK(metadata.st_mode)) or (
                not is_final
                and (
                    stat.S_ISLNK(metadata.st_mode)
                    or not stat.S_ISDIR(metadata.st_mode)
                )
            ):
                raise RuntimeError(
                    f"tracked index prefix {logical!r} is an ancestor link, directory "
                    "reparse point, or non-directory"
                )

            state = (worktree_kind(metadata), metadata_identity(metadata))
            if logical in states and states[logical] != state:
                raise RuntimeError(
                    f"tracked index prefix {logical!r} moved during materialization inspection"
                )
            states[logical] = state
            if os.name != "nt":
                identity = (metadata.st_dev, metadata.st_ino)
                owner = identity_owners.get(identity)
                if owner is not None and owner != logical:
                    raise RuntimeError(
                        f"tracked index prefixes {owner!r} and {logical!r} resolve to one "
                        "filesystem identity; case-folding, normalization, and hard-link "
                        "aliases are not admissible"
                    )
                identity_owners[identity] = logical
    return tuple(
        (logical, kind, identity)
        for logical, (kind, identity) in sorted(states.items())
    )

def replace_ref_inventory(repo):
    output = git(
        repo,
        "for-each-ref",
        "--format=%(refname)%00%(objectname)",
        "refs/replace/",
    )
    if output and not output.endswith(b"\n"):
        raise RuntimeError("git for-each-ref returned no terminal newline")
    records = []
    for record in output.splitlines():
        fields = record.split(b"\0")
        if len(fields) != 2 or not fields[0].startswith(b"refs/replace/"):
            raise RuntimeError("git for-each-ref emitted malformed replace authority")
        object_id = fields[1]
        if len(object_id) not in {40, 64} or any(
            byte not in b"0123456789abcdef" for byte in object_id
        ):
            raise RuntimeError("git for-each-ref emitted malformed replace object id")
        records.append((fields[0], object_id))
    if len(records) != len(set(records)):
        raise RuntimeError("git for-each-ref emitted duplicate replace authority")
    return tuple(sorted(records))

def graft_authority_inventory(repo):
    grafts = git_line(repo, "rev-parse", "--git-path", "info/grafts")
    if not os.path.isabs(grafts):
        grafts = os.path.abspath(os.path.join(repo, grafts))
    try:
        metadata = os.lstat(grafts)
    except FileNotFoundError:
        return (b"absent",)
    if not stat.S_ISREG(metadata.st_mode):
        return (b"unsupported", worktree_kind(metadata))
    identity = metadata_identity(metadata)
    with open(grafts, "rb") as source:
        opened = os.fstat(source.fileno())
        if metadata_identity(opened) != identity:
            raise RuntimeError(f"graft authority moved before reading: {grafts!r}")
        content = source.read(MAX_GRAFT_BYTES + 1)
    confirmed = os.lstat(grafts)
    if metadata_identity(confirmed) != identity:
        raise RuntimeError(f"graft authority moved while reading: {grafts!r}")
    if len(content) > MAX_GRAFT_BYTES:
        raise RuntimeError("graft authority exceeds the 1 MiB inspection bound")
    return (b"regular", len(content), hashlib.sha256(content).digest())

def parse_hash_lines(output, expected_count, command):
    if expected_count and not output.endswith(b"\n"):
        raise RuntimeError(f"{command} returned no terminal newline")
    hashes = output.splitlines()
    if len(hashes) != expected_count:
        raise RuntimeError(
            f"{command} returned {len(hashes)} hashes for {expected_count} paths"
        )
    for object_id in hashes:
        if not object_id or any(byte not in b"0123456789abcdef" for byte in object_id):
            raise RuntimeError(f"{command} returned a malformed object id")
    return hashes

def bounded_path_batches(paths):
    batch = []
    batch_bytes = 0
    for path in paths:
        path_bytes = len(path) + 1
        if path_bytes > MAX_HASH_BATCH_BYTES:
            raise RuntimeError(f"tracked path exceeds hash batch bound: {path!r}")
        if batch and (
            len(batch) >= MAX_HASH_BATCH_PATHS
            or batch_bytes + path_bytes > MAX_HASH_BATCH_BYTES
        ):
            yield batch
            batch = []
            batch_bytes = 0
        batch.append(path)
        batch_bytes += path_bytes
    if batch:
        yield batch

def hash_regular_paths(repo, paths):
    hashed = {}
    for batch in bounded_path_batches(paths):
        output = git(repo, "hash-object", "--no-filters", "--", *batch)
        object_ids = parse_hash_lines(output, len(batch), "git hash-object --no-filters")
        hashed.update(zip(batch, object_ids))
    return hashed

def hash_symlink_target(repo, target):
    output = git_stdin(repo, target, "hash-object", "--no-filters", "--stdin")
    return parse_hash_lines(output, 1, "git hash-object --no-filters --stdin")[0]

def raw_tracked_sources(repo, index):
    parsed_index = parse_index(index)
    materialization_before = index_materialization_snapshot(repo, sorted(parsed_index))
    observations = []
    pending_regular = []
    regular_metadata = {}
    for path in sorted(parsed_index):
        records = parsed_index[path]
        if len(records) != 1 or records[0][0] != b"0":
            observations.append((b"conflicted", path, tuple(sorted(records))))
            continue
        _stage, expected_mode, expected_object_id = records[0]
        if expected_mode == b"160000":
            observations.append(
                (b"gitlink", path, expected_mode, expected_object_id)
            )
            continue
        if expected_mode not in {b"100644", b"100755", b"120000"}:
            observations.append(
                (b"unsupported-index-mode", path, expected_mode, expected_object_id)
            )
            continue
        absolute = safe_worktree_path(repo, path)
        try:
            metadata = os.lstat(absolute)
        except FileNotFoundError:
            observations.append((b"missing", path, expected_mode, expected_object_id))
            continue
        identity = metadata_identity(metadata)
        if expected_mode in {b"100644", b"100755"}:
            if not stat.S_ISREG(metadata.st_mode):
                observations.append(
                    (
                        b"wrong-type",
                        path,
                        expected_mode,
                        expected_object_id,
                        worktree_kind(metadata),
                        identity,
                    )
                )
                continue
            pending_regular.append(path)
            regular_metadata[path] = (expected_mode, expected_object_id, identity)
            continue
        if not stat.S_ISLNK(metadata.st_mode):
            observations.append(
                (
                    b"wrong-type",
                    path,
                    expected_mode,
                    expected_object_id,
                    worktree_kind(metadata),
                    identity,
                )
            )
            continue
        target_before = os.fsencode(os.readlink(absolute))
        actual_object_id = hash_symlink_target(repo, target_before)
        confirmed_metadata = os.lstat(absolute)
        target_after = os.fsencode(os.readlink(absolute))
        if metadata_identity(confirmed_metadata) != identity or target_after != target_before:
            raise RuntimeError(f"tracked symlink moved while hashing: {path!r}")
        observations.append(
            (
                b"symlink",
                path,
                expected_mode,
                expected_object_id,
                b"120000",
                actual_object_id,
                identity,
            )
        )

    regular_hashes = hash_regular_paths(repo, pending_regular)
    for path in pending_regular:
        expected_mode, expected_object_id, identity = regular_metadata[path]
        confirmed_metadata = os.lstat(safe_worktree_path(repo, path))
        if metadata_identity(confirmed_metadata) != identity:
            raise RuntimeError(f"tracked file moved while hashing: {path!r}")
        actual_mode = (
            b"100755" if confirmed_metadata.st_mode & stat.S_IXUSR else b"100644"
        ) if FILE_MODE_OBSERVABLE else b"<unavailable>"
        observations.append(
            (
                b"regular",
                path,
                expected_mode,
                expected_object_id,
                actual_mode,
                regular_hashes[path],
                identity,
            )
        )
    materialization_after = index_materialization_snapshot(repo, sorted(parsed_index))
    if materialization_after != materialization_before:
        raise RuntimeError("tracked index path materialization moved while hashing")
    observations.append((b"index-materialization", b"", materialization_before))
    return tuple(sorted(observations, key=lambda observation: (observation[1], observation[0])))

def forced_visible_status(repo):
    return git(
        repo,
        "-c",
        "core.fileMode=true" if FILE_MODE_OBSERVABLE else "core.fileMode=false",
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "status",
        "--porcelain=v1",
        "-z",
        "--untracked-files=no",
        "--ignore-submodules=none",
        "--no-renames",
    )

def observe_repository(repo):
    repo = safe_repository_path(repo)
    require_no_executable_config(repo)
    boundary = repository_boundary(repo)
    require_no_executable_config(boundary)
    primary_index = raw_primary_index_state(boundary)
    replace_refs = replace_ref_inventory(boundary)
    graft_authority = graft_authority_inventory(boundary)
    index = git(
        boundary,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "ls-files",
        "--stage",
        "-z",
        "--",
    )
    primary_index = bind_raw_primary_index(primary_index, index)
    raw_sources = raw_tracked_sources(boundary, index)
    staged = git(
        boundary,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "diff",
        "--cached",
        "--name-only",
        "-z",
        "--no-renames",
        "--no-ext-diff",
        "--ignore-submodules=none",
        "--",
    )
    # Intentionally name only .gitignore as an exclusion source. This keeps
    # committed project ignore policy while bypassing .git/info/exclude and
    # core.excludesFile, including a caller's global ignore file.
    untracked = git(
        boundary,
        "-c",
        "core.excludesFile=/dev/null",
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "ls-files",
        "-z",
        "--others",
        "--exclude-per-directory=.gitignore",
        "--",
    )
    untracked_gitignores = git(
        boundary,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "ls-files",
        "-z",
        "--others",
        "--",
        ".gitignore",
        ":(glob)**/.gitignore",
    )
    flags_t = git(
        boundary,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "ls-files",
        "-t",
        "-z",
        "--",
    )
    flags_v = git(
        boundary,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "ls-files",
        "-v",
        "-z",
        "--",
    )
    return (
        boundary,
        optional_head(boundary),
        replace_refs,
        graft_authority,
        index,
        raw_sources,
        staged,
        untracked,
        untracked_gitignores,
        flags_t,
        flags_v,
        primary_index,
    )

def finding(scope, kind, path, detail, action):
    return (scope, kind, path, detail, action)

def raw_source_findings(scope, raw_sources):
    findings = []
    for observation in raw_sources:
        state = observation[0]
        path = observation[1]
        if state in {b"gitlink", b"index-materialization"}:
            continue
        if state in {b"regular", b"symlink"}:
            (
                _state,
                _path,
                expected_mode,
                expected_object_id,
                actual_mode,
                actual_object_id,
                _identity,
            ) = observation
            mode_matches = (
                state == b"regular" and not FILE_MODE_OBSERVABLE
            ) or actual_mode == expected_mode
            if mode_matches and actual_object_id == expected_object_id:
                continue
            findings.append(
                finding(
                    scope,
                    b"raw-tracked-source-mismatch",
                    path,
                    b"expected="
                    + expected_mode
                    + b":"
                    + expected_object_id
                    + b" actual="
                    + actual_mode
                    + b":"
                    + actual_object_id,
                    b"restore the exact raw indexed bytes and mode deliberately",
                )
            )
            continue
        if state == b"conflicted":
            findings.append(
                finding(
                    scope,
                    b"conflicted-index-entry",
                    path,
                    b"no unique stage-0 tracked source exists",
                    b"resolve the index conflict before verification",
                )
            )
            continue
        if state == b"missing":
            findings.append(
                finding(
                    scope,
                    b"missing-tracked-source",
                    path,
                    b"stage-0 source is absent from the worktree",
                    b"restore the indexed source before verification",
                )
            )
            continue
        if state == b"wrong-type":
            findings.append(
                finding(
                    scope,
                    b"tracked-source-type-mismatch",
                    path,
                    b"expected=" + observation[2] + b" actual=" + observation[4],
                    b"restore the indexed file type before verification",
                )
            )
            continue
        if state == b"unsupported-index-mode":
            findings.append(
                finding(
                    scope,
                    b"unsupported-index-mode",
                    path,
                    b"mode=" + observation[2],
                    b"use a canonical regular, symlink, or gitlink index mode",
                )
            )
            continue
        raise RuntimeError(f"unknown raw tracked-source state: {state!r}")
    return findings

def authority_findings(scope, replace_refs, graft_authority):
    findings = []
    for refname, object_id in replace_refs:
        findings.append(
            finding(
                scope,
                b"replace-ref-authority",
                refname,
                b"replacement-object=" + object_id,
                b"remove the replace ref deliberately before verification",
            )
        )
    if graft_authority[0] == b"regular" and graft_authority[1] > 0:
        findings.append(
            finding(
                scope,
                b"graft-authority",
                b".git/info/grafts",
                b"nonempty-byte-length=" + str(graft_authority[1]).encode("ascii"),
                b"empty or remove the graft authority deliberately before verification",
            )
        )
    elif graft_authority[0] == b"unsupported":
        findings.append(
            finding(
                scope,
                b"graft-authority-type",
                b".git/info/grafts",
                b"unsupported-kind=" + graft_authority[1],
                b"replace it with an empty regular file or no file before verification",
            )
        )
    return findings

def findings_for_observation(scope, observation, expected_head, status):
    (
        _boundary,
        head,
        replace_refs,
        graft_authority,
        index,
        raw_sources,
        staged,
        untracked,
        untracked_gitignores,
        flags_t,
        flags_v,
        _primary_index,
    ) = observation
    findings = authority_findings(scope, replace_refs, graft_authority)
    findings.extend(raw_source_findings(scope, raw_sources))
    if expected_head is not None and head != expected_head:
        findings.append(
            finding(
                scope,
                b"gitlink-head-mismatch",
                b".",
                b"expected=" + expected_head + b" actual=" + (head or b"<unborn>"),
                b"check out the parent-recorded gitlink deliberately",
            )
        )
    for path in parse_paths(staged, "git diff --cached --name-only"):
        findings.append(
            finding(
                scope,
                b"staged-index-change",
                path,
                b"index differs from repository HEAD",
                b"inspect and restore or commit the staged change deliberately",
            )
        )
    for path, code in parse_status(status):
        findings.append(
            finding(
                scope,
                b"tracked-or-index-change",
                path,
                b"porcelain=" + code,
                b"inspect and restore or commit the tracked change deliberately",
            )
        )
    visible_untracked = parse_paths(untracked, "git ls-files --others")
    for path in visible_untracked:
        findings.append(
            finding(
                scope,
                b"untracked-not-project-ignored",
                path,
                b"local and global excludes do not exempt constellation source",
                b"remove, relocate, project-ignore, or commit the path deliberately",
            )
        )
    visible_untracked_set = set(visible_untracked)
    for path in parse_paths(untracked_gitignores, "untracked .gitignore inventory"):
        if path in visible_untracked_set:
            continue
        findings.append(
            finding(
                scope,
                b"untracked-ignore-policy",
                path,
                b"only tracked .gitignore files may define project ignore semantics",
                b"remove or commit the .gitignore deliberately",
            )
        )
    tagged_t = parse_tagged_paths(flags_t, "git ls-files -t")
    tagged_v = parse_tagged_paths(flags_v, "git ls-files -v")
    if tagged_t.keys() != tagged_v.keys():
        raise RuntimeError(f"index flag inventories disagree in repository {scope!r}")
    for path in sorted(tagged_t):
        if tagged_t[path] == b"S":
            findings.append(
                finding(
                    scope,
                    b"skip-worktree",
                    path,
                    b"index flag can conceal worktree state",
                    b"clear with update-index --no-skip-worktree before verification",
                )
            )
        if tagged_v[path].islower():
            findings.append(
                finding(
                    scope,
                    b"assume-unchanged",
                    path,
                    b"index flag can conceal worktree state",
                    b"clear with update-index --no-assume-unchanged before verification",
                )
            )
    return findings, parse_index(index)

def inspect_repository(
    repo,
    scope,
    expected_head,
    active_roots,
    constellation_root,
    depth,
):
    if depth > MAX_REPOSITORY_DEPTH:
        raise RuntimeError(
            f"recursive repository depth exceeds {MAX_REPOSITORY_DEPTH}: {scope!r}"
        )
    before = observe_repository(repo)
    repo_root = before[0]
    try:
        common = os.path.commonpath((constellation_root, repo_root))
    except ValueError as error:
        raise RuntimeError(f"repository escaped the constellation root: {error}") from error
    if common != constellation_root:
        raise RuntimeError(
            f"repository escaped the constellation root: root={constellation_root!r}, repo={repo_root!r}"
        )
    if repo_root in active_roots:
        raise RuntimeError(f"recursive repository cycle at {repo_root!r}")
    active_roots.add(repo_root)
    try:
        parsed_index = parse_index(before[4])
        findings = []
        child_snapshots = []
        for relative in sorted(parsed_index):
            stage_zero = [
                (mode, object_id)
                for stage, mode, object_id in parsed_index[relative]
                if stage == b"0"
            ]
            if len(stage_zero) > 1:
                raise RuntimeError(f"multiple stage-0 index records for {relative!r}")
            if not stage_zero or stage_zero[0][0] != b"160000":
                continue
            expected_child_head = stage_zero[0][1]
            child = safe_worktree_path(repo_root, relative)
            child_scope = scope + b"/" + relative
            try:
                metadata = os.lstat(child)
            except FileNotFoundError:
                child_snapshots.append((child_scope, b"uninitialized-missing"))
                continue
            if not stat.S_ISDIR(metadata.st_mode):
                kind = b"symlink" if stat.S_ISLNK(metadata.st_mode) else b"non-directory"
                findings.append(
                    finding(
                        scope,
                        b"gitlink-worktree-obstruction",
                        relative,
                        b"kind=" + kind,
                        b"restore an empty uninitialized path or the exact initialized submodule",
                    )
                )
                child_snapshots.append((child_scope, b"obstructed-" + kind))
                continue
            marker = git_marker(child)
            if marker is None:
                with os.scandir(child) as entries:
                    is_empty = next(entries, None) is None
                if is_empty:
                    child_snapshots.append((child_scope, b"uninitialized-empty"))
                    continue
                findings.append(
                    finding(
                        scope,
                        b"gitlink-worktree-obstruction",
                        relative,
                        b"uninitialized gitlink directory is not empty",
                        b"relocate its bytes or initialize the exact recorded submodule deliberately",
                    )
                )
                child_snapshots.append((child_scope, b"uninitialized-nonempty"))
                continue
            require_no_executable_config(child)
            child_root = repository_boundary(child)
            try:
                child_common = os.path.commonpath((constellation_root, child_root))
            except ValueError as error:
                raise RuntimeError(f"submodule escaped the constellation root: {error}") from error
            if child_common != constellation_root:
                raise RuntimeError(
                    f"submodule escaped the constellation root: path={child!r}, root={child_root!r}"
                )
            child_snapshot, child_findings = inspect_repository(
                child_root,
                child_scope,
                expected_child_head,
                active_roots,
                constellation_root,
                depth + 1,
            )
            child_snapshots.append(child_snapshot)
            findings.extend(child_findings)
        # A forced-visible status may enter initialized submodules. Every
        # descendant has now passed executable-config admission, so it is safe
        # to override ignore policy at this parent boundary.
        status = forced_visible_status(repo_root)
        after = observe_repository(repo_root)
        confirmed_status = forced_visible_status(repo_root)
        if after != before or confirmed_status != status:
            raise RuntimeError(f"repository moved during cleanliness inspection: {repo_root!r}")
        repository_findings, _ = findings_for_observation(
            scope,
            before,
            expected_head,
            status,
        )
        findings.extend(repository_findings)
        snapshot = (scope, before, status, tuple(child_snapshots))
        return snapshot, sorted(findings)
    finally:
        active_roots.remove(repo_root)

def inspect_tree(path):
    path = safe_repository_path(path)
    require_no_executable_config(path)
    constellation_root = repository_boundary(path)
    return inspect_repository(
        constellation_root,
        b".",
        None,
        set(),
        constellation_root,
        0,
    )

def quoted(value):
    return repr(value)

try:
    first_snapshot, first_findings = inspect_tree(root)
    confirmed_snapshot, confirmed_findings = inspect_tree(root)
    if confirmed_snapshot != first_snapshot or confirmed_findings != first_findings:
        raise RuntimeError("repository tree moved between complete cleanliness observations")
    for scope, kind, path, detail, action in first_findings:
        print(
            "  repository={} kind={} path={} detail={} action={}".format(
                quoted(scope),
                kind.decode("ascii"),
                quoted(path),
                quoted(detail),
                quoted(action),
            )
        )
except (OSError, RuntimeError) as error:
    print(
        f"FATAL: recursive cleanliness inspection failed for {root!r}: {error}",
        file=sys.stderr,
    )
    raise SystemExit(1)
PY
}

verify_clean() { # directory expected-head
    local dir="$1" expected="$2" actual status confirmed
    if ! actual="$(git_at "$dir" rev-parse HEAD 2>/dev/null)"; then
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
        printf 'FATAL: %s is dirty at pinned head %s:\n%s\n' \
            "$dir" "$expected" "$status" >&2
        return 1
    fi
    if ! confirmed="$(git_at "$dir" rev-parse HEAD 2>/dev/null)" \
        || [[ "$confirmed" != "$expected" ]]; then
        echo "FATAL: $dir moved while its pinned state was being verified" >&2
        return 1
    fi
}

snapshot_identity() {
    python3 -I - "$repo_root" "$parent" "$lock" "$validated_lock_sha256" <<'PY'
import hashlib
import json
import os
import stat
import subprocess
import sys

root, parent, lock_path, validated_lock_sha256 = sys.argv[1:]
FILE_MODE_OBSERVABLE = os.name != "nt"
git_environment = os.environ.copy()
git_environment["GIT_OPTIONAL_LOCKS"] = "0"
git_environment["GIT_NO_REPLACE_OBJECTS"] = "1"
for inherited in (
    "GIT_DIR",
    "GIT_WORK_TREE",
    "GIT_INDEX_FILE",
    "GIT_COMMON_DIR",
    "GIT_OBJECT_DIRECTORY",
    "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    "GIT_NAMESPACE",
    "GIT_PREFIX",
    "GIT_REPLACE_REF_BASE",
    "GIT_SHALLOW_FILE",
    "GIT_CEILING_DIRECTORIES",
    "GIT_DISCOVERY_ACROSS_FILESYSTEM",
    "GIT_CONFIG",
    "GIT_CONFIG_COUNT",
    "GIT_CONFIG_PARAMETERS",
    "GIT_CONFIG_GLOBAL",
    "GIT_CONFIG_SYSTEM",
    "GIT_CONFIG_NOSYSTEM",
    "GIT_ATTR_NOSYSTEM",
    "GIT_ATTR_SOURCE",
    "GIT_TEMPLATE_DIR",
    "GIT_DEFAULT_HASH",
    "GIT_DEFAULT_REF_FORMAT",
    "GIT_REFERENCE_BACKEND",
    "GIT_EXEC_PATH",
    "GIT_EXTERNAL_DIFF",
    "GIT_ASKPASS",
    "SSH_ASKPASS",
    "GIT_SSH",
    "GIT_SSH_COMMAND",
    "GIT_PROXY_COMMAND",
    "GIT_ALLOW_PROTOCOL",
    "GIT_PROTOCOL_FROM_USER",
    "GIT_TERMINAL_PROMPT",
    "GIT_QUARANTINE_PATH",
    "GIT_REDIRECT_STDIN",
    "GIT_REDIRECT_STDOUT",
    "GIT_REDIRECT_STDERR",
    "GIT_LITERAL_PATHSPECS",
    "GIT_GLOB_PATHSPECS",
    "GIT_NOGLOB_PATHSPECS",
    "GIT_ICASE_PATHSPECS",
    "GIT_TRACE",
    "GIT_TRACE_CURL",
    "GIT_TRACE_CURL_NO_DATA",
    "GIT_TRACE_FSMONITOR",
    "GIT_TRACE_PACK_ACCESS",
    "GIT_TRACE_PACKET",
    "GIT_TRACE_PACKFILE",
    "GIT_TRACE_PERFORMANCE",
    "GIT_TRACE_REFS",
    "GIT_TRACE_SETUP",
    "GIT_TRACE_SHALLOW",
    "GIT_TRACE2",
    "GIT_TRACE2_EVENT",
    "GIT_TRACE2_PERF",
    "GIT_NO_LAZY_FETCH",
):
    git_environment.pop(inherited, None)
git_environment["GIT_CONFIG_NOSYSTEM"] = "1"
git_environment["GIT_ATTR_NOSYSTEM"] = "1"
git_environment["GIT_CONFIG_GLOBAL"] = os.devnull
git_environment["GIT_TERMINAL_PROMPT"] = "0"
git_environment["GIT_NO_LAZY_FETCH"] = "1"
try:
    with open(lock_path, "rb") as source:
        lock_bytes = source.read(1_048_577)
    if len(lock_bytes) > 1_048_576:
        raise RuntimeError("constellation lock exceeds the 1 MiB parser bound")
    if hashlib.sha256(lock_bytes).hexdigest() != validated_lock_sha256:
        raise RuntimeError("constellation lock moved after strict validation")
    document = json.loads(lock_bytes)
except (OSError, RuntimeError, TypeError, json.JSONDecodeError) as error:
    print(f"FATAL: cannot capture validated constellation lock: {error}", file=sys.stderr)
    raise SystemExit(1)

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

def run_git(repo, *args):
    return subprocess.run(
        [
            "git",
            "-c", f"core.hooksPath={os.devnull}",
            "-c", f"core.attributesFile={os.devnull}",
            "-c", "credential.helper=",
            "-c", "protocol.allow=never",
            "-c", "protocol.file.allow=always",
            "-c", "protocol.https.allow=always",
            "-c", "protocol.ssh.allow=always",
            "-C", repo,
            *args,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=git_environment,
    )

def git(repo, *args):
    completed = run_git(repo, *args)
    if completed.returncode != 0:
        detail = completed.stderr.decode(errors="replace").strip()
        raise RuntimeError(f"git {args!r} in {repo} failed: {detail}")
    return completed.stdout

def git_stdin(repo, data, *args):
    completed = subprocess.run(
        [
            "git",
            "-c", f"core.hooksPath={os.devnull}",
            "-c", f"core.attributesFile={os.devnull}",
            "-c", "credential.helper=",
            "-c", "protocol.allow=never",
            "-c", "protocol.file.allow=always",
            "-c", "protocol.https.allow=always",
            "-c", "protocol.ssh.allow=always",
            "-C", repo,
            *args,
        ],
        input=data,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=git_environment,
    )
    if completed.returncode != 0:
        detail = completed.stderr.decode(errors="replace").strip()
        raise RuntimeError(f"git {args!r} in {repo} failed: {detail}")
    return completed.stdout

def read_git_line(repo, *args):
    output = git(repo, *args)
    if not output.endswith(b"\n"):
        raise RuntimeError(f"git {args!r} in {repo!r} returned no terminal newline")
    return output[:-1]

def executable_config_key(key):
    key = key.lower()
    return (
        key.startswith((b"include.", b"includeif.", b"filter."))
        or key in {
            b"core.fsmonitor",
            b"core.hookspath",
            b"core.sshcommand",
            b"core.gitproxy",
            b"core.alternaterefscommand",
            b"core.askpass",
            b"credential.helper",
            b"extensions.refstorage",
            b"fetch.bundleuri",
            b"fetch.uriprotocols",
            b"gc.recentobjectshook",
            b"protocol.allow",
        }
        or (key.startswith(b"protocol.") and key.endswith(b".allow"))
        or (key.startswith(b"credential.") and key.endswith(b".helper"))
        or (
            key.startswith(b"remote.")
            and key.endswith((b".uploadpack", b".receivepack", b".vcs"))
        )
        or (key.startswith(b"submodule.") and key.endswith(b".update"))
        or (
            key.startswith(b"url.")
            and key.endswith((b".insteadof", b".pushinsteadof"))
        )
    )

def require_no_executable_config(repo, scope):
    worktree_config = run_git(
        repo,
        "config",
        "--local",
        "--no-includes",
        "--type=bool",
        "--get-all",
        "extensions.worktreeConfig",
    )
    if (
        worktree_config.returncode == 1
        and not worktree_config.stdout
        and not worktree_config.stderr
    ):
        worktree_config_enabled = False
    elif worktree_config.returncode != 0:
        detail = worktree_config.stderr.decode(errors="replace").strip()
        raise RuntimeError(
            "cannot inspect extensions.worktreeConfig in "
            f"{repo!r}: {detail}"
        )
    elif worktree_config.stderr:
        detail = worktree_config.stderr.decode(errors="replace").strip()
        raise RuntimeError(
            "extensions.worktreeConfig inspection emitted an unexpected "
            f"diagnostic in {repo!r}: {detail}"
        )
    elif worktree_config.stdout == b"true\n":
        worktree_config_enabled = True
    elif worktree_config.stdout == b"false\n":
        worktree_config_enabled = False
    else:
        raise RuntimeError(
            "extensions.worktreeConfig must have exactly one canonical "
            "boolean value"
        )
    for config_scope in ("--local", "--worktree"):
        if config_scope == "--worktree" and not worktree_config_enabled:
            continue
        inventory = git(
            repo,
            "config",
            config_scope,
            "--null",
            "--name-only",
            "--no-includes",
            "--list",
        )
        if inventory and not inventory.endswith(b"\0"):
            raise RuntimeError(
                f"{config_scope} Git configuration inventory is not NUL-terminated"
            )
        for key in filter(None, inventory.split(b"\0")):
            if executable_config_key(key):
                raise RuntimeError(
                    f"repository {scope!r} has {config_scope} Git configuration "
                    f"authority {key!r} that may execute code or redirect history"
                )

def windows_reparse_point(metadata):
    if os.name != "nt":
        return False
    attributes = getattr(metadata, "st_file_attributes", None)
    if attributes is None:
        raise RuntimeError("Windows file metadata lacks reparse-point attributes")
    return bool(attributes & 0x400)

def metadata_identity(metadata):
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_mode,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )

primary_index_git_line = read_git_line
MAX_RAW_PRIMARY_INDEX_BYTES = 1_073_741_824

def primary_index_git_dir(repo):
    reported = primary_index_git_line(
        repo,
        "rev-parse",
        "--absolute-git-dir",
    )
    if not os.path.isabs(reported):
        raise RuntimeError(f"reported primary git dir is not absolute: {reported!r}")
    metadata = os.lstat(reported)
    if windows_reparse_point(metadata) or not stat.S_ISDIR(metadata.st_mode):
        raise RuntimeError(f"primary git dir is not an ordinary directory: {reported!r}")
    canonical = os.path.realpath(reported)
    if os.name != "nt":
        if reported != canonical:
            raise RuntimeError(
                f"reported primary git dir redirects from {reported!r} to {canonical!r}"
            )
        return reported, metadata_identity(metadata)
    canonical_metadata = os.lstat(canonical)
    if windows_reparse_point(canonical_metadata) \
            or not stat.S_ISDIR(canonical_metadata.st_mode) \
            or metadata_identity(canonical_metadata) != metadata_identity(metadata):
        raise RuntimeError(
            f"resolved primary git dir is not the reported ordinary directory: {canonical!r}"
        )
    return canonical, metadata_identity(canonical_metadata)

def primary_index_path(repo, git_dir):
    reported = primary_index_git_line(
        repo,
        "rev-parse",
        "--path-format=absolute",
        "--git-path",
        "index",
    )
    if not os.path.isabs(reported):
        raise RuntimeError(f"reported primary index path is not absolute: {reported!r}")
    if os.path.basename(reported) != b"index":
        raise RuntimeError(f"reported primary index basename is not 'index': {reported!r}")
    parent_authority = os.path.realpath(os.path.dirname(reported))
    if parent_authority != git_dir:
        raise RuntimeError(
            f"reported primary index parent {parent_authority!r} is not {git_dir!r}"
        )
    if os.name != "nt" and reported != os.path.join(git_dir, b"index"):
        raise RuntimeError(
            f"reported primary index is not a direct canonical child: {reported!r}"
        )
    return os.path.join(git_dir, b"index")

def raw_primary_index_state(repo):
    object_format = primary_index_git_line(
        repo,
        "rev-parse",
        "--show-object-format=storage",
    )
    if object_format == b"sha1":
        object_id_width = 20
        checksum = hashlib.sha1
    elif object_format == b"sha256":
        object_id_width = 32
        checksum = hashlib.sha256
    else:
        raise RuntimeError(
            f"unsupported primary-index object format: {object_format!r}"
        )
    git_dir, git_dir_identity = primary_index_git_dir(repo)
    index_path = primary_index_path(repo, git_dir)
    try:
        metadata = os.lstat(index_path)
    except FileNotFoundError:
        confirmed_object_format = primary_index_git_line(
            repo,
            "rev-parse",
            "--show-object-format=storage",
        )
        try:
            os.lstat(index_path)
        except FileNotFoundError:
            pass
        else:
            raise RuntimeError("primary index appeared while inspecting its absence")
        confirmed_git_dir, confirmed_git_dir_identity = primary_index_git_dir(repo)
        confirmed_index_path = primary_index_path(repo, confirmed_git_dir)
        if confirmed_object_format != object_format \
                or confirmed_git_dir != git_dir \
                or confirmed_git_dir_identity != git_dir_identity \
                or confirmed_index_path != index_path:
            raise RuntimeError("primary-index authority moved while inspecting its absence")
        return b"absent\0" + object_format, None
    if windows_reparse_point(metadata) or not stat.S_ISREG(metadata.st_mode):
        raise RuntimeError(f"primary index is not an ordinary file: {index_path!r}")
    if metadata.st_size > MAX_RAW_PRIMARY_INDEX_BYTES:
        raise RuntimeError("primary index exceeds the 1 GiB raw inspection bound")
    identity = metadata_identity(metadata)
    with open(index_path, "rb") as source:
        opened_metadata = os.fstat(source.fileno())
        if metadata_identity(opened_metadata) != identity:
            raise RuntimeError("primary index moved before raw inspection")
        raw_index = source.read(metadata.st_size)
        grew_during_read = source.read(1)
        confirmed_opened_metadata = os.fstat(source.fileno())
    confirmed_metadata = os.lstat(index_path)
    if windows_reparse_point(confirmed_metadata) \
            or metadata_identity(confirmed_opened_metadata) != identity \
            or metadata_identity(confirmed_metadata) != identity:
        raise RuntimeError("primary index moved during raw inspection")
    if len(raw_index) != metadata.st_size or grew_during_read:
        raise RuntimeError("primary index changed length during raw inspection")
    confirmed_object_format = primary_index_git_line(
        repo,
        "rev-parse",
        "--show-object-format=storage",
    )
    confirmed_git_dir, confirmed_git_dir_identity = primary_index_git_dir(repo)
    confirmed_index_path = primary_index_path(repo, confirmed_git_dir)
    if confirmed_object_format != object_format \
            or confirmed_git_dir != git_dir \
            or confirmed_git_dir_identity != git_dir_identity \
            or confirmed_index_path != index_path:
        raise RuntimeError("primary-index authority moved during raw inspection")

    if len(raw_index) < 12 + object_id_width:
        raise RuntimeError("primary index is truncated before its checksum")
    if raw_index[:4] != b"DIRC":
        raise RuntimeError("primary index has an invalid DIRC signature")
    version = int.from_bytes(raw_index[4:8], "big")
    if version == 4:
        raise RuntimeError("primary index version 4 path compression is unsupported")
    if version not in {2, 3}:
        raise RuntimeError(f"unsupported primary index version: {version}")
    entry_count = int.from_bytes(raw_index[8:12], "big")
    checksum_offset = len(raw_index) - object_id_width
    expected_checksum = raw_index[checksum_offset:]
    actual_checksum = checksum(raw_index[:checksum_offset]).digest()
    if actual_checksum != expected_checksum:
        raise RuntimeError(
            f"primary index {object_format.decode('ascii')} checksum is invalid"
        )

    cursor = 12
    entry_order = []
    fixed_entry_size = 40 + object_id_width + 2
    for entry_number in range(entry_count):
        entry_start = cursor
        fixed_entry_end = entry_start + fixed_entry_size
        if fixed_entry_end > checksum_offset:
            raise RuntimeError(
                f"primary index entry {entry_number} is truncated before its flags"
            )
        entry_mode = int.from_bytes(raw_index[entry_start + 24:entry_start + 28], "big")
        flags = int.from_bytes(raw_index[fixed_entry_end - 2:fixed_entry_end], "big")
        cursor = fixed_entry_end
        if flags & 0x4000:
            if version != 3:
                raise RuntimeError("primary index v2 entry has unsupported extended flags")
            if cursor + 2 > checksum_offset:
                raise RuntimeError(
                    f"primary index entry {entry_number} has truncated extended flags"
                )
            extended_flags = int.from_bytes(raw_index[cursor:cursor + 2], "big")
            if extended_flags & 0x9FFF:
                raise RuntimeError(
                    f"primary index entry {entry_number} has unsupported extended flags"
                )
            cursor += 2
        path_start = cursor
        encoded_name_length = flags & 0x0FFF
        if encoded_name_length < 0x0FFF:
            path_end = path_start + encoded_name_length
            if path_end >= checksum_offset:
                raise RuntimeError(
                    f"primary index entry {entry_number} is truncated in its path"
                )
            if b"\0" in raw_index[path_start:path_end]:
                raise RuntimeError(
                    f"primary index entry {entry_number} path contains an early NUL"
                )
        else:
            path_end = raw_index.find(b"\0", path_start, checksum_offset)
            if path_end < 0:
                raise RuntimeError(
                    f"primary index entry {entry_number} has no path terminator"
                )
            if path_end - path_start < 0x0FFF:
                raise RuntimeError(
                    f"primary index entry {entry_number} has an ambiguous long-path length"
                )
        stage = (flags >> 12) & 0x3
        path = raw_index[path_start:path_end]
        validate_git_worktree_path(
            path,
            symlink=entry_mode & 0o170000 == 0o120000,
        )
        entry_order.append((path, stage))
        padding_size = 8 - ((path_end - entry_start) % 8)
        entry_end = path_end + padding_size
        if entry_end > checksum_offset:
            raise RuntimeError(
                f"primary index entry {entry_number} has truncated path padding"
            )
        if raw_index[path_end:entry_end] != b"\0" * padding_size:
            raise RuntimeError(
                f"primary index entry {entry_number} has malformed path padding"
            )
        cursor = entry_end

    extensions = []
    while cursor < checksum_offset:
        if checksum_offset - cursor < 8:
            raise RuntimeError("primary index has a truncated extension header")
        signature = raw_index[cursor:cursor + 4]
        extension_size = int.from_bytes(raw_index[cursor + 4:cursor + 8], "big")
        extension_end = cursor + 8 + extension_size
        if extension_end > checksum_offset:
            raise RuntimeError(
                f"primary index extension {signature!r} is truncated"
            )
        if signature == b"FSMN":
            raise RuntimeError(
                "primary index contains forbidden FSMonitor FSMN extension"
            )
        if signature == b"link":
            raise RuntimeError("split-index link extension is unsupported")
        if signature == b"sdir":
            raise RuntimeError("sparse-index sdir extension is unsupported")
        if not 0x41 <= signature[0] <= 0x5A:
            raise RuntimeError(
                f"unsupported required primary-index extension: {signature!r}"
            )
        extensions.append(signature + extension_size.to_bytes(4, "big"))
        cursor = extension_end
    if cursor != checksum_offset:
        raise RuntimeError("primary index extension layout is ambiguous")
    previous_entry = None
    for path, stage in entry_order:
        current_entry = (path, stage)
        if not path:
            raise RuntimeError("primary index contains an empty entry path")
        if previous_entry is not None and current_entry == previous_entry:
            raise RuntimeError(
                f"primary index contains duplicate entry order: {current_entry!r}"
            )
        if previous_entry is not None and current_entry < previous_entry:
            raise RuntimeError(
                "primary index entries are not strictly sorted: "
                f"previous={previous_entry!r} current={current_entry!r}"
            )
        previous_entry = current_entry

    state = b"present\0" + object_format + b"\0" \
        + version.to_bytes(4, "big") \
        + entry_count.to_bytes(4, "big") \
        + b"".join(extensions)
    return state, entry_count

def bind_raw_primary_index(raw_state, index_inventory):
    state, raw_entry_count = raw_state
    if index_inventory:
        if not index_inventory.endswith(b"\0"):
            raise RuntimeError("git ls-files --stage inventory is not NUL-terminated")
        inventory_records = index_inventory[:-1].split(b"\0")
        if any(not record for record in inventory_records):
            raise RuntimeError("git ls-files --stage inventory contains an empty record")
        inventory_entry_count = len(inventory_records)
    else:
        inventory_entry_count = 0
    if raw_entry_count is None:
        if inventory_entry_count:
            raise RuntimeError(
                "primary index is absent but git ls-files reported index entries"
            )
    elif raw_entry_count != inventory_entry_count:
        raise RuntimeError(
            "primary-index entry count disagrees with git ls-files: "
            f"raw={raw_entry_count} inventory={inventory_entry_count}"
        )
    return state

MAX_GRAFT_BYTES = 1_048_576
MAX_REPOSITORY_DEPTH = 32

def snapshot_worktree_kind(metadata):
    if stat.S_ISREG(metadata.st_mode):
        return b"regular"
    if stat.S_ISLNK(metadata.st_mode):
        return b"symlink"
    if stat.S_ISDIR(metadata.st_mode):
        return b"directory"
    return b"special"

def replace_ref_inventory(repo):
    output = git(
        repo,
        "for-each-ref",
        "--format=%(refname)%00%(objectname)",
        "refs/replace/",
    )
    if output and not output.endswith(b"\n"):
        raise RuntimeError("git for-each-ref returned no terminal newline")
    records = []
    for record in output.splitlines():
        fields = record.split(b"\0")
        if len(fields) != 2 or not fields[0].startswith(b"refs/replace/"):
            raise RuntimeError("git for-each-ref emitted malformed replace authority")
        object_id = fields[1]
        if len(object_id) not in {40, 64} or any(
            byte not in b"0123456789abcdef" for byte in object_id
        ):
            raise RuntimeError("git for-each-ref emitted malformed replace object id")
        records.append((fields[0], object_id))
    if len(records) != len(set(records)):
        raise RuntimeError("git for-each-ref emitted duplicate replace authority")
    return tuple(sorted(records))

def graft_authority_inventory(repo):
    grafts = read_git_line(repo, "rev-parse", "--git-path", "info/grafts")
    if not os.path.isabs(grafts):
        grafts = os.path.abspath(os.path.join(os.fsencode(repo), grafts))
    try:
        metadata = os.lstat(grafts)
    except FileNotFoundError:
        return (b"absent",)
    if not stat.S_ISREG(metadata.st_mode):
        return (b"unsupported", snapshot_worktree_kind(metadata))
    identity = metadata_identity(metadata)
    with open(grafts, "rb") as source:
        opened = os.fstat(source.fileno())
        if metadata_identity(opened) != identity:
            raise RuntimeError(f"graft authority moved before reading: {grafts!r}")
        content = source.read(MAX_GRAFT_BYTES + 1)
    confirmed = os.lstat(grafts)
    if metadata_identity(confirmed) != identity:
        raise RuntimeError(f"graft authority moved while reading: {grafts!r}")
    if len(content) > MAX_GRAFT_BYTES:
        raise RuntimeError("graft authority exceeds the 1 MiB inspection bound")
    return (b"regular", len(content), hashlib.sha256(content).digest())

def require_no_history_overrides(repo, scope):
    replace_refs = replace_ref_inventory(repo)
    graft_authority = graft_authority_inventory(repo)
    if replace_refs:
        raise RuntimeError(
            f"repository {scope!r} has replace-ref authority: {replace_refs!r}"
        )
    if graft_authority[0] == b"regular" and graft_authority[1] > 0:
        raise RuntimeError(
            f"repository {scope!r} has nonempty .git/info/grafts authority"
        )
    if graft_authority[0] == b"unsupported":
        raise RuntimeError(
            f"repository {scope!r} has unsupported .git/info/grafts kind: "
            f"{graft_authority[1]!r}"
        )
    return replace_refs, graft_authority

def file_digest(path, expected_metadata):
    hashed = hashlib.sha256()
    length = 0
    with open(path, "rb") as source:
        opened_metadata = os.fstat(source.fileno())
        if metadata_identity(opened_metadata) != metadata_identity(expected_metadata):
            raise RuntimeError(f"working-tree file moved before hashing: {path!r}")
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            length += len(chunk)
            hashed.update(chunk)
    confirmed_metadata = os.lstat(path)
    if metadata_identity(confirmed_metadata) != metadata_identity(expected_metadata):
        raise RuntimeError(f"working-tree file moved while hashing: {path!r}")
    if length != expected_metadata.st_size:
        raise RuntimeError(f"working-tree file changed length while hashing: {path!r}")
    return length, hashed.digest()

def symlink_target(path, expected_metadata):
    target = os.fsencode(os.readlink(path))
    confirmed_metadata = os.lstat(path)
    confirmed_target = os.fsencode(os.readlink(path))
    if metadata_identity(confirmed_metadata) != metadata_identity(expected_metadata) \
            or confirmed_target != target:
        raise RuntimeError(f"working-tree symlink moved while hashing: {path!r}")
    return target

HFS_IGNORED_UTF8 = (
    b"\xe2\x80\x8c", b"\xe2\x80\x8d", b"\xe2\x80\x8e", b"\xe2\x80\x8f",
    b"\xe2\x80\xaa", b"\xe2\x80\xab", b"\xe2\x80\xac", b"\xe2\x80\xad",
    b"\xe2\x80\xae", b"\xe2\x81\xaa", b"\xe2\x81\xab", b"\xe2\x81\xac",
    b"\xe2\x81\xad", b"\xe2\x81\xae", b"\xe2\x81\xaf", b"\xef\xbb\xbf",
)

def windows_ntfs_fallback_alias(component, prefix):
    candidate = component.rstrip(b" .").lower()
    if len(candidate) != 8:
        return False
    saw_tilde = False
    for index, byte in enumerate(candidate):
        if saw_tilde:
            if not 0x30 <= byte <= 0x39:
                return False
        elif byte == 0x7E:
            if index + 1 >= len(candidate) or not 0x31 <= candidate[index + 1] <= 0x39:
                return False
            saw_tilde = True
        elif index >= 6 or byte != prefix[index]:
            return False
    return saw_tilde

def windows_gitmodules_alias(component):
    candidate = component.rstrip(b" .").lower()
    return (
        candidate == b".gitmodules"
        or (
            len(candidate) == 8
            and candidate[:7] == b"gitmod~"
            and candidate[7:8] in {b"1", b"2", b"3", b"4"}
        )
        or windows_ntfs_fallback_alias(candidate, b"gi7eba")
    )

def validate_git_worktree_path(relative, symlink=False):
    if not relative or b"\0" in relative:
        raise RuntimeError(f"git emitted an empty or NUL-bearing path: {relative!r}")
    components = relative.split(b"/")
    if relative.startswith(b"/") or relative.endswith(b"/") or any(
        component in {b"", b".", b".."} for component in components
    ):
        raise RuntimeError(f"git emitted an unsafe working-tree path: {relative!r}")
    for component in components:
        if component.lower() == b".git":
            raise RuntimeError(
                f"git emitted a path through repository metadata: {relative!r}"
            )
        if symlink and component.lower() == b".gitmodules":
            raise RuntimeError(
                f"git emitted a forbidden .gitmodules symlink: {relative!r}"
            )
        if any(sequence in component for sequence in HFS_IGNORED_UTF8):
            raise RuntimeError(
                f"git emitted an HFS-ambiguous Unicode path: {relative!r}"
            )
        if os.name == "nt":
            windows_alias = component.rstrip(b" .").lower()
            short_name = windows_alias[4:] if windows_alias.startswith(b"git~") else b""
            if windows_alias == b".git" or (
                1 <= len(short_name) <= 6 and short_name.isdigit()
            ):
                raise RuntimeError(
                    f"git emitted a Windows .git alias: {relative!r}"
                )
            if symlink and windows_gitmodules_alias(component):
                raise RuntimeError(
                    f"git emitted a Windows .gitmodules symlink alias: {relative!r}"
                )
    if os.name == "nt" and (b"\\" in relative or b":" in relative):
        raise RuntimeError(
            "git emitted a rooted, drive-qualified, backslash, or "
            f"alternate-stream path: {relative!r}"
        )
    return relative

def parse_index(index):
    parsed = {}
    for record in filter(None, index.split(b"\0")):
        if b"\t" not in record:
            raise RuntimeError("git ls-files --stage emitted a malformed record")
        metadata, relative = record.split(b"\t", 1)
        fields = metadata.split()
        if len(fields) != 3 or fields[2] not in {b"0", b"1", b"2", b"3"}:
            raise RuntimeError("git ls-files --stage emitted malformed metadata")
        mode, _object_id, stage = fields
        validate_git_worktree_path(relative, symlink=mode == b"120000")
        parsed.setdefault(relative, []).append((stage, mode, record))
    return parsed

def encoded_index_records(records):
    return b"".join(record + b"\0" for _stage, _mode, record in sorted(records))

def git_listing(repo, *args):
    return git(
        repo,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        *args,
    )

def repository_listing(repo):
    return git_listing(
        repo,
        "-c",
        "core.excludesFile=/dev/null",
        "ls-files",
        "-z",
        "--cached",
        "--others",
        "--exclude-per-directory=.gitignore",
    )

def untracked_gitignore_inventory(repo):
    output = git_listing(
        repo,
        "ls-files",
        "-z",
        "--others",
        "--",
        ".gitignore",
        ":(glob)**/.gitignore",
    )
    paths = tuple(sorted(filter(None, output.split(b"\0"))))
    if len(paths) != len(set(paths)):
        raise RuntimeError("untracked .gitignore inventory emitted duplicate paths")
    return paths

def repository_top(repo):
    return os.path.realpath(read_git_line(repo, "rev-parse", "--show-toplevel"))

def has_git_marker(path):
    marker = os.path.join(os.fsencode(path), b".git")
    try:
        metadata = os.lstat(marker)
    except FileNotFoundError:
        return False
    if windows_reparse_point(metadata) or not (
        stat.S_ISREG(metadata.st_mode) or stat.S_ISDIR(metadata.st_mode)
    ):
        raise RuntimeError(f"unsupported .git marker type at {marker!r}")
    return True

def is_initialized_repository(path):
    path = os.fsencode(path)
    try:
        metadata = os.lstat(path)
    except FileNotFoundError:
        return False
    if stat.S_ISLNK(metadata.st_mode) or windows_reparse_point(metadata):
        raise RuntimeError(f"repository root must not be a link or reparse point: {path!r}")
    if not stat.S_ISDIR(metadata.st_mode):
        return False
    if not has_git_marker(path):
        return False
    require_no_executable_config(path, b"repository-preflight")
    if read_git_line(path, "rev-parse", "--is-inside-work-tree") != b"true":
        raise RuntimeError(f"nested repository is not a worktree: {path!r}")
    expected_top = os.path.realpath(os.fsencode(path))
    actual_top = repository_top(path)
    if actual_top != expected_top:
        raise RuntimeError(
            f"nested repository escaped its gitlink path: path={path!r}, top={actual_top!r}"
        )
    return True

def worktree_path(repo, relative):
    return os.path.join(os.fsencode(repo), validate_git_worktree_path(relative))

def capture_repository(
    repo,
    scope,
    parent_index_records,
    active_git_dirs,
    depth,
):
    if depth > MAX_REPOSITORY_DEPTH:
        raise RuntimeError(
            f"recursive repository depth exceeds {MAX_REPOSITORY_DEPTH}: {scope!r}"
        )
    repo_bytes = os.fsencode(repo)
    expected_top = os.path.realpath(repo_bytes)
    if not is_initialized_repository(repo_bytes):
        raise RuntimeError(f"repository is not initialized at {repo_bytes!r}")
    require_no_executable_config(repo_bytes, scope)
    git_dir, git_dir_identity = primary_index_git_dir(repo_bytes)
    if git_dir in active_git_dirs:
        raise RuntimeError(f"nested repository cycle through git dir {git_dir!r}")
    active_git_dirs.add(git_dir)
    try:
        primary_index_before = raw_primary_index_state(repo_bytes)
        replace_refs_before, graft_authority_before = require_no_history_overrides(
            repo_bytes,
            scope,
        )
        head_before = read_git_line(repo_bytes, "rev-parse", "HEAD")
        index_before = git_listing(repo_bytes, "ls-files", "--stage", "-z")
        primary_index_before = bind_raw_primary_index(
            primary_index_before,
            index_before,
        )
        index_flags_before = git_listing(repo_bytes, "ls-files", "-v", "-z")
        listed_before = repository_listing(repo_bytes)
        untracked_gitignores_before = untracked_gitignore_inventory(repo_bytes)
        if untracked_gitignores_before:
            raise RuntimeError(
                f"repository {scope!r} has untracked .gitignore authority: "
                f"{untracked_gitignores_before!r}"
            )
        parsed_index = parse_index(index_before)
        materialization_before = cleanliness_index_materialization_snapshot(
            repo_bytes, sorted(parsed_index)
        )

        frames = []

        def frame(label, data):
            if isinstance(label, str):
                label = label.encode()
            if isinstance(data, str):
                data = data.encode()
            frames.append((label, data))

        frame("repository-begin", scope)
        frame("repository-parent-index-records", parent_index_records)
        frame("repository-head", head_before)
        frame("repository-replace-refs", b"")
        frame("repository-graft-authority", graft_authority_before[0])
        frame("repository-index", index_before)
        frame("repository-index-flags", index_flags_before)
        frame("repository-primary-index-state", primary_index_before)
        frame("repository-untracked-gitignores", b"")

        paths = sorted(set(filter(None, listed_before.split(b"\0"))))
        for relative in paths:
            absolute = worktree_path(repo_bytes, relative)
            records = parsed_index.get(relative, [])
            records_bytes = encoded_index_records(records)
            is_gitlink = any(mode == b"160000" for _stage, mode, _record in records)
            frame("repository-entry-begin", relative)
            frame("repository-entry-index-records", records_bytes)
            try:
                metadata = os.lstat(absolute)
            except FileNotFoundError:
                if is_gitlink:
                    frame("repository-entry-kind", "gitlink")
                    frame("repository-gitlink-state", "uninitialized")
                    frame("repository-gitlink-worktree", "missing")
                else:
                    frame("repository-entry-kind", "missing")
                frame("repository-entry-end", relative)
                continue

            if windows_reparse_point(metadata) and not stat.S_ISLNK(metadata.st_mode):
                raise RuntimeError(
                    f"working-tree entry is a non-symbolic reparse point: {relative!r}"
                )

            if is_gitlink:
                frame("repository-entry-kind", "gitlink")
                if stat.S_ISDIR(metadata.st_mode):
                    if is_initialized_repository(absolute):
                        frame("repository-gitlink-state", "initialized")
                        child_scope = scope + b"/" + relative
                        frames.extend(
                            capture_repository(
                                absolute,
                                child_scope,
                                records_bytes,
                                active_git_dirs,
                                depth + 1,
                            )
                        )
                    else:
                        with os.scandir(absolute) as entries:
                            if next(entries, None) is not None:
                                raise RuntimeError(
                                    "uninitialized gitlink contains unframed bytes at "
                                    f"{absolute!r}"
                                )
                        frame("repository-gitlink-state", "uninitialized")
                        frame("repository-gitlink-worktree", "empty-directory")
                elif stat.S_ISREG(metadata.st_mode):
                    frame("repository-gitlink-state", "uninitialized")
                    frame("repository-gitlink-worktree", "obstructing-file")
                    mode = "100755" if metadata.st_mode & stat.S_IXUSR else "100644"
                    length, content = file_digest(absolute, metadata)
                    frame("repository-entry-mode", mode)
                    frame("repository-entry-length", length.to_bytes(8, "big"))
                    frame("repository-entry-sha256", content)
                elif stat.S_ISLNK(metadata.st_mode):
                    frame("repository-gitlink-state", "uninitialized")
                    frame("repository-gitlink-worktree", "obstructing-symlink")
                    frame("repository-entry-mode", "120000")
                    frame("repository-entry-target", symlink_target(absolute, metadata))
                else:
                    raise RuntimeError(
                        f"unsupported gitlink worktree entry type: {relative!r}"
                    )
                frame("repository-entry-end", relative)
                continue

            if stat.S_ISREG(metadata.st_mode):
                mode = "100755" if metadata.st_mode & stat.S_IXUSR else "100644"
                length, content = file_digest(absolute, metadata)
                frame("repository-entry-kind", "file")
                frame("repository-entry-mode", mode)
                frame("repository-entry-length", length.to_bytes(8, "big"))
                frame("repository-entry-sha256", content)
            elif stat.S_ISLNK(metadata.st_mode):
                frame("repository-entry-kind", "symlink")
                frame("repository-entry-mode", "120000")
                frame("repository-entry-target", symlink_target(absolute, metadata))
            elif stat.S_ISDIR(metadata.st_mode) and is_initialized_repository(absolute):
                frame("repository-entry-kind", "embedded-git-repository")
                frame("repository-gitlink-state", "initialized")
                child_scope = scope + b"/" + relative
                frames.extend(
                    capture_repository(
                        absolute,
                        child_scope,
                        b"",
                        active_git_dirs,
                        depth + 1,
                    )
                )
            else:
                raise RuntimeError(
                    f"unsupported working-tree entry type: {os.fsdecode(relative)!r}"
                )
            frame("repository-entry-end", relative)

        head_after = read_git_line(repo_bytes, "rev-parse", "HEAD")
        replace_refs_after, graft_authority_after = require_no_history_overrides(
            repo_bytes,
            scope,
        )
        index_after = git_listing(repo_bytes, "ls-files", "--stage", "-z")
        primary_index_after = raw_primary_index_state(repo_bytes)
        primary_index_after = bind_raw_primary_index(
            primary_index_after,
            index_after,
        )
        index_flags_after = git_listing(repo_bytes, "ls-files", "-v", "-z")
        listed_after = repository_listing(repo_bytes)
        untracked_gitignores_after = untracked_gitignore_inventory(repo_bytes)
        materialization_after = cleanliness_index_materialization_snapshot(
            repo_bytes, sorted(parsed_index)
        )
        confirmed_git_dir, confirmed_git_dir_identity = primary_index_git_dir(
            repo_bytes
        )
        confirmed_top = repository_top(repo_bytes)
        if (
            head_after != head_before
            or replace_refs_after != replace_refs_before
            or graft_authority_after != graft_authority_before
            or index_after != index_before
            or index_flags_after != index_flags_before
            or primary_index_after != primary_index_before
            or listed_after != listed_before
            or untracked_gitignores_after != untracked_gitignores_before
            or materialization_after != materialization_before
            or confirmed_git_dir != git_dir
            or confirmed_git_dir_identity != git_dir_identity
            or confirmed_top != expected_top
        ):
            raise RuntimeError(f"repository moved during snapshot: {repo_bytes!r}")
        frame("repository-end", scope)
        return frames
    finally:
        active_git_dirs.remove(git_dir)

def parse_cleanliness_status(output):
    parsed = []
    for record in filter(None, output.split(b"\0")):
        if len(record) < 4 or record[2:3] != b" ":
            raise RuntimeError("git status --porcelain=v1 emitted a malformed record")
        parsed.append((record[3:], record[:2]))
    return parsed

def parse_cleanliness_paths(output, command):
    parsed = list(filter(None, output.split(b"\0")))
    if len(parsed) != len(set(parsed)):
        raise RuntimeError(f"{command} emitted duplicate paths")
    return parsed

def parse_cleanliness_tags(output, command):
    parsed = {}
    for record in filter(None, output.split(b"\0")):
        if len(record) < 3 or record[1:2] != b" ":
            raise RuntimeError(f"{command} emitted a malformed tagged path")
        parsed[record[2:]] = record[:1]
    return parsed

MAX_HASH_BATCH_PATHS = 128
MAX_HASH_BATCH_BYTES = 65_536

def cleanliness_worktree_kind(metadata):
    if stat.S_ISREG(metadata.st_mode):
        return b"regular"
    if stat.S_ISLNK(metadata.st_mode):
        return b"symlink"
    if stat.S_ISDIR(metadata.st_mode):
        return b"directory"
    return b"special"

def cleanliness_index_materialization_snapshot(repo, paths):
    states = {}
    identity_owners = {}
    for path in paths:
        worktree_path(repo, path)
        logical = b""
        materialized = repo
        components = path.split(b"/")
        for component_index, component in enumerate(components):
            parent = materialized
            logical = component if not logical else logical + b"/" + component
            materialized = os.path.join(materialized, component)
            try:
                metadata = os.lstat(materialized)
            except FileNotFoundError:
                break

            if os.name == "nt":
                try:
                    with os.scandir(parent) as entries:
                        exact_matches = sum(
                            os.fsencode(entry.name) == component for entry in entries
                        )
                except OSError as error:
                    raise RuntimeError(
                        f"cannot enumerate tracked index prefix parent {parent!r}: {error}"
                    ) from error
                if exact_matches != 1:
                    raise RuntimeError(
                        f"tracked index prefix {logical!r} is not materialized with "
                        "one exact directory-entry spelling"
                    )

            is_final = component_index + 1 == len(components)
            is_reparse_point = windows_reparse_point(metadata)
            if (is_reparse_point and not stat.S_ISLNK(metadata.st_mode)) or (
                not is_final
                and (
                    stat.S_ISLNK(metadata.st_mode)
                    or not stat.S_ISDIR(metadata.st_mode)
                )
            ):
                raise RuntimeError(
                    f"tracked index prefix {logical!r} is an ancestor link, directory "
                    "reparse point, or non-directory"
                )

            state = (cleanliness_worktree_kind(metadata), metadata_identity(metadata))
            if logical in states and states[logical] != state:
                raise RuntimeError(
                    f"tracked index prefix {logical!r} moved during materialization inspection"
                )
            states[logical] = state
            if os.name != "nt":
                identity = (metadata.st_dev, metadata.st_ino)
                owner = identity_owners.get(identity)
                if owner is not None and owner != logical:
                    raise RuntimeError(
                        f"tracked index prefixes {owner!r} and {logical!r} resolve to one "
                        "filesystem identity; case-folding, normalization, and hard-link "
                        "aliases are not admissible"
                    )
                identity_owners[identity] = logical
    return tuple(
        (logical, kind, identity)
        for logical, (kind, identity) in sorted(states.items())
    )

def parse_cleanliness_hashes(output, expected_count, command):
    if expected_count and not output.endswith(b"\n"):
        raise RuntimeError(f"{command} returned no terminal newline")
    hashes = output.splitlines()
    if len(hashes) != expected_count:
        raise RuntimeError(
            f"{command} returned {len(hashes)} hashes for {expected_count} paths"
        )
    for object_id in hashes:
        if not object_id or any(byte not in b"0123456789abcdef" for byte in object_id):
            raise RuntimeError(f"{command} returned a malformed object id")
    return hashes

def bounded_cleanliness_batches(paths):
    batch = []
    batch_bytes = 0
    for path in paths:
        path_bytes = len(path) + 1
        if path_bytes > MAX_HASH_BATCH_BYTES:
            raise RuntimeError(f"tracked path exceeds hash batch bound: {path!r}")
        if batch and (
            len(batch) >= MAX_HASH_BATCH_PATHS
            or batch_bytes + path_bytes > MAX_HASH_BATCH_BYTES
        ):
            yield batch
            batch = []
            batch_bytes = 0
        batch.append(path)
        batch_bytes += path_bytes
    if batch:
        yield batch

def hash_cleanliness_regular_paths(repo, paths):
    hashed = {}
    for batch in bounded_cleanliness_batches(paths):
        output = git(repo, "hash-object", "--no-filters", "--", *batch)
        object_ids = parse_cleanliness_hashes(
            output,
            len(batch),
            "git hash-object --no-filters",
        )
        hashed.update(zip(batch, object_ids))
    return hashed

def hash_cleanliness_symlink_target(repo, target):
    output = git_stdin(repo, target, "hash-object", "--no-filters", "--stdin")
    return parse_cleanliness_hashes(
        output,
        1,
        "git hash-object --no-filters --stdin",
    )[0]

def cleanliness_index_record(record):
    metadata, _relative = record.split(b"\t", 1)
    fields = metadata.split()
    if len(fields) != 3:
        raise RuntimeError("git ls-files --stage emitted malformed metadata")
    return fields[1]

def raw_cleanliness_sources(repo, index):
    parsed_index = parse_index(index)
    materialization_before = cleanliness_index_materialization_snapshot(
        repo, sorted(parsed_index)
    )
    observations = []
    pending_regular = []
    regular_metadata = {}
    for path in sorted(parsed_index):
        records = parsed_index[path]
        normalized_records = tuple(
            sorted(
                (stage, mode, cleanliness_index_record(record))
                for stage, mode, record in records
            )
        )
        if len(normalized_records) != 1 or normalized_records[0][0] != b"0":
            observations.append((b"conflicted", path, normalized_records))
            continue
        _stage, expected_mode, expected_object_id = normalized_records[0]
        if expected_mode == b"160000":
            observations.append(
                (b"gitlink", path, expected_mode, expected_object_id)
            )
            continue
        if expected_mode not in {b"100644", b"100755", b"120000"}:
            observations.append(
                (b"unsupported-index-mode", path, expected_mode, expected_object_id)
            )
            continue
        absolute = worktree_path(repo, path)
        try:
            metadata = os.lstat(absolute)
        except FileNotFoundError:
            observations.append((b"missing", path, expected_mode, expected_object_id))
            continue
        identity = metadata_identity(metadata)
        if expected_mode in {b"100644", b"100755"}:
            if not stat.S_ISREG(metadata.st_mode):
                observations.append(
                    (
                        b"wrong-type",
                        path,
                        expected_mode,
                        expected_object_id,
                        cleanliness_worktree_kind(metadata),
                        identity,
                    )
                )
                continue
            pending_regular.append(path)
            regular_metadata[path] = (expected_mode, expected_object_id, identity)
            continue
        if not stat.S_ISLNK(metadata.st_mode):
            observations.append(
                (
                    b"wrong-type",
                    path,
                    expected_mode,
                    expected_object_id,
                    cleanliness_worktree_kind(metadata),
                    identity,
                )
            )
            continue
        target_before = os.fsencode(os.readlink(absolute))
        actual_object_id = hash_cleanliness_symlink_target(repo, target_before)
        confirmed_metadata = os.lstat(absolute)
        target_after = os.fsencode(os.readlink(absolute))
        if metadata_identity(confirmed_metadata) != identity or target_after != target_before:
            raise RuntimeError(f"tracked symlink moved while hashing: {path!r}")
        observations.append(
            (
                b"symlink",
                path,
                expected_mode,
                expected_object_id,
                b"120000",
                actual_object_id,
                identity,
            )
        )

    regular_hashes = hash_cleanliness_regular_paths(repo, pending_regular)
    for path in pending_regular:
        expected_mode, expected_object_id, identity = regular_metadata[path]
        confirmed_metadata = os.lstat(worktree_path(repo, path))
        if metadata_identity(confirmed_metadata) != identity:
            raise RuntimeError(f"tracked file moved while hashing: {path!r}")
        actual_mode = (
            b"100755" if confirmed_metadata.st_mode & stat.S_IXUSR else b"100644"
        ) if FILE_MODE_OBSERVABLE else b"<unavailable>"
        observations.append(
            (
                b"regular",
                path,
                expected_mode,
                expected_object_id,
                actual_mode,
                regular_hashes[path],
                identity,
            )
        )
    materialization_after = cleanliness_index_materialization_snapshot(
        repo, sorted(parsed_index)
    )
    if materialization_after != materialization_before:
        raise RuntimeError("tracked index path materialization moved while hashing")
    observations.append((b"index-materialization", b"", materialization_before))
    return tuple(sorted(observations, key=lambda observation: (observation[1], observation[0])))

def forced_visible_cleanliness_status(repo):
    return git(
        repo,
        "-c",
        "core.fileMode=true" if FILE_MODE_OBSERVABLE else "core.fileMode=false",
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "status",
        "--porcelain=v1",
        "-z",
        "--untracked-files=no",
        "--ignore-submodules=none",
        "--no-renames",
    )

def cleanliness_observation(repo, scope):
    repo_bytes = os.fsencode(repo)
    if not is_initialized_repository(repo_bytes):
        raise RuntimeError(f"repository is not initialized at {repo_bytes!r}")
    require_no_executable_config(repo_bytes, scope)
    top = repository_top(repo_bytes)
    primary_index = raw_primary_index_state(repo_bytes)
    head = read_git_line(repo_bytes, "rev-parse", "HEAD")
    replace_refs, graft_authority = require_no_history_overrides(repo_bytes, scope)
    index = git_listing(repo_bytes, "ls-files", "--stage", "-z", "--")
    primary_index = bind_raw_primary_index(primary_index, index)
    raw_sources = raw_cleanliness_sources(repo_bytes, index)
    staged = git(
        repo_bytes,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "diff",
        "--cached",
        "--name-only",
        "-z",
        "--no-renames",
        "--no-ext-diff",
        "--ignore-submodules=none",
        "--",
    )
    # The committed per-directory .gitignore policy remains authoritative.
    # Local .git/info/exclude and global core.excludesFile are deliberately
    # bypassed so they cannot conceal constellation source.
    untracked = git(
        repo_bytes,
        "-c",
        "core.excludesFile=/dev/null",
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "ls-files",
        "-z",
        "--others",
        "--exclude-per-directory=.gitignore",
        "--",
    )
    untracked_gitignores = git(
        repo_bytes,
        "-c",
        "core.fsmonitor=false",
        "-c",
        "core.untrackedCache=false",
        "ls-files",
        "-z",
        "--others",
        "--",
        ".gitignore",
        ":(glob)**/.gitignore",
    )
    flags_t = git_listing(repo_bytes, "ls-files", "-t", "-z", "--")
    flags_v = git_listing(repo_bytes, "ls-files", "-v", "-z", "--")
    return (
        top,
        head,
        replace_refs,
        graft_authority,
        index,
        raw_sources,
        staged,
        untracked,
        untracked_gitignores,
        flags_t,
        flags_v,
        primary_index,
    )

def cleanliness_finding(scope, kind, path, detail, action):
    return (scope, kind, path, detail, action)

def raw_cleanliness_findings(scope, raw_sources):
    findings = []
    for observation in raw_sources:
        state = observation[0]
        path = observation[1]
        if state in {b"gitlink", b"index-materialization"}:
            continue
        if state in {b"regular", b"symlink"}:
            (
                _state,
                _path,
                expected_mode,
                expected_object_id,
                actual_mode,
                actual_object_id,
                _identity,
            ) = observation
            mode_matches = (
                state == b"regular" and not FILE_MODE_OBSERVABLE
            ) or actual_mode == expected_mode
            if mode_matches and actual_object_id == expected_object_id:
                continue
            findings.append(
                cleanliness_finding(
                    scope,
                    b"raw-tracked-source-mismatch",
                    path,
                    b"expected="
                    + expected_mode
                    + b":"
                    + expected_object_id
                    + b" actual="
                    + actual_mode
                    + b":"
                    + actual_object_id,
                    b"restore the exact raw indexed bytes and mode deliberately",
                )
            )
            continue
        if state == b"conflicted":
            findings.append(
                cleanliness_finding(
                    scope,
                    b"conflicted-index-entry",
                    path,
                    b"no unique stage-0 tracked source exists",
                    b"resolve the index conflict before verification",
                )
            )
            continue
        if state == b"missing":
            findings.append(
                cleanliness_finding(
                    scope,
                    b"missing-tracked-source",
                    path,
                    b"stage-0 source is absent from the worktree",
                    b"restore the indexed source before verification",
                )
            )
            continue
        if state == b"wrong-type":
            findings.append(
                cleanliness_finding(
                    scope,
                    b"tracked-source-type-mismatch",
                    path,
                    b"expected=" + observation[2] + b" actual=" + observation[4],
                    b"restore the indexed file type before verification",
                )
            )
            continue
        if state == b"unsupported-index-mode":
            findings.append(
                cleanliness_finding(
                    scope,
                    b"unsupported-index-mode",
                    path,
                    b"mode=" + observation[2],
                    b"use a canonical regular, symlink, or gitlink index mode",
                )
            )
            continue
        raise RuntimeError(f"unknown raw tracked-source state: {state!r}")
    return findings

def observation_findings(scope, observation, expected_head, status):
    (
        _top,
        head,
        _replace_refs,
        _graft_authority,
        index,
        raw_sources,
        staged,
        untracked,
        untracked_gitignores,
        flags_t,
        flags_v,
        _primary_index,
    ) = observation
    findings = raw_cleanliness_findings(scope, raw_sources)
    if head != expected_head:
        findings.append(
            cleanliness_finding(
                scope,
                b"gitlink-head-mismatch",
                b".",
                b"expected=" + expected_head + b" actual=" + head,
                b"check out the parent-recorded gitlink deliberately",
            )
        )
    for path in parse_cleanliness_paths(staged, "git diff --cached --name-only"):
        findings.append(
            cleanliness_finding(
                scope,
                b"staged-index-change",
                path,
                b"index differs from repository HEAD",
                b"inspect and restore or commit the staged change deliberately",
            )
        )
    for path, code in parse_cleanliness_status(status):
        findings.append(
            cleanliness_finding(
                scope,
                b"tracked-or-index-change",
                path,
                b"porcelain=" + code,
                b"inspect and restore or commit the tracked change deliberately",
            )
        )
    visible_untracked = parse_cleanliness_paths(untracked, "git ls-files --others")
    for path in visible_untracked:
        findings.append(
            cleanliness_finding(
                scope,
                b"untracked-not-project-ignored",
                path,
                b"local and global excludes do not exempt constellation source",
                b"remove, relocate, project-ignore, or commit the path deliberately",
            )
        )
    visible_untracked_set = set(visible_untracked)
    for path in parse_cleanliness_paths(
        untracked_gitignores,
        "untracked .gitignore inventory",
    ):
        if path in visible_untracked_set:
            continue
        findings.append(
            cleanliness_finding(
                scope,
                b"untracked-ignore-policy",
                path,
                b"only tracked .gitignore files may define project ignore semantics",
                b"remove or commit the .gitignore deliberately",
            )
        )
    tagged_t = parse_cleanliness_tags(flags_t, "git ls-files -t")
    tagged_v = parse_cleanliness_tags(flags_v, "git ls-files -v")
    if tagged_t.keys() != tagged_v.keys():
        raise RuntimeError(f"index flag inventories disagree in repository {scope!r}")
    for path in sorted(tagged_t):
        if tagged_t[path] == b"S":
            findings.append(
                cleanliness_finding(
                    scope,
                    b"skip-worktree",
                    path,
                    b"index flag can conceal worktree state",
                    b"clear with update-index --no-skip-worktree before verification",
                )
            )
        if tagged_v[path].islower():
            findings.append(
                cleanliness_finding(
                    scope,
                    b"assume-unchanged",
                    path,
                    b"index flag can conceal worktree state",
                    b"clear with update-index --no-assume-unchanged before verification",
                )
            )
    return findings, parse_index(index)

def stage_zero_gitlink(records):
    stage_zero = []
    for stage, mode, record in records:
        if stage != b"0":
            continue
        metadata, _relative = record.split(b"\t", 1)
        fields = metadata.split()
        if len(fields) != 3:
            raise RuntimeError("git ls-files --stage emitted malformed metadata")
        stage_zero.append((mode, fields[1]))
    if len(stage_zero) > 1:
        raise RuntimeError("git ls-files --stage emitted multiple stage-0 records")
    if not stage_zero or stage_zero[0][0] != b"160000":
        return None
    return stage_zero[0][1]

def inspect_cleanliness(
    repo,
    scope,
    expected_head,
    active_roots,
    constellation_root,
    depth,
):
    if depth > MAX_REPOSITORY_DEPTH:
        raise RuntimeError(
            f"recursive repository depth exceeds {MAX_REPOSITORY_DEPTH}: {scope!r}"
        )
    before = cleanliness_observation(repo, scope)
    repo_root = before[0]
    try:
        common = os.path.commonpath((constellation_root, repo_root))
    except ValueError as error:
        raise RuntimeError(f"repository escaped the constellation root: {error}") from error
    if common != constellation_root:
        raise RuntimeError(
            f"repository escaped the constellation root: root={constellation_root!r}, repo={repo_root!r}"
        )
    if repo_root in active_roots:
        raise RuntimeError(f"recursive repository cycle at {repo_root!r}")
    active_roots.add(repo_root)
    try:
        parsed_index = parse_index(before[4])
        findings = []
        child_snapshots = []
        for relative in sorted(parsed_index):
            expected_child_head = stage_zero_gitlink(parsed_index[relative])
            if expected_child_head is None:
                continue
            child = worktree_path(repo_root, relative)
            child_scope = scope + b"/" + relative
            try:
                metadata = os.lstat(child)
            except FileNotFoundError:
                child_snapshots.append((child_scope, b"uninitialized-missing"))
                continue
            if not stat.S_ISDIR(metadata.st_mode):
                kind = b"symlink" if stat.S_ISLNK(metadata.st_mode) else b"non-directory"
                findings.append(
                    cleanliness_finding(
                        scope,
                        b"gitlink-worktree-obstruction",
                        relative,
                        b"kind=" + kind,
                        b"restore an empty uninitialized path or the exact initialized submodule",
                    )
                )
                child_snapshots.append((child_scope, b"obstructed-" + kind))
                continue
            if not has_git_marker(child):
                with os.scandir(child) as entries:
                    is_empty = next(entries, None) is None
                if is_empty:
                    child_snapshots.append((child_scope, b"uninitialized-empty"))
                    continue
                findings.append(
                    cleanliness_finding(
                        scope,
                        b"gitlink-worktree-obstruction",
                        relative,
                        b"uninitialized gitlink directory is not empty",
                        b"relocate its bytes or initialize the exact recorded submodule deliberately",
                    )
                )
                child_snapshots.append((child_scope, b"uninitialized-nonempty"))
                continue
            child_root = repository_top(child)
            try:
                child_common = os.path.commonpath((constellation_root, child_root))
            except ValueError as error:
                raise RuntimeError(f"submodule escaped the constellation root: {error}") from error
            if child_common != constellation_root:
                raise RuntimeError(
                    f"submodule escaped the constellation root: path={child!r}, root={child_root!r}"
                )
            child_snapshot, child_findings = inspect_cleanliness(
                child_root,
                child_scope,
                expected_child_head,
                active_roots,
                constellation_root,
                depth + 1,
            )
            child_snapshots.append(child_snapshot)
            findings.extend(child_findings)
        # Run forced-visible status only after recursive admission has ruled out
        # executable child configuration at every initialized descendant.
        status = forced_visible_cleanliness_status(repo_root)
        after = cleanliness_observation(repo_root, scope)
        confirmed_status = forced_visible_cleanliness_status(repo_root)
        if after != before or confirmed_status != status:
            raise RuntimeError(f"repository moved during cleanliness inspection: {repo_root!r}")
        repository_findings, _ = observation_findings(
            scope,
            before,
            expected_head,
            status,
        )
        findings.extend(repository_findings)
        return (scope, before, status, tuple(child_snapshots)), sorted(findings)
    finally:
        active_roots.remove(repo_root)

def inspect_cleanliness_tree(repo, expected_head, scope):
    repo_bytes = os.fsencode(repo)
    if not is_initialized_repository(repo_bytes):
        raise RuntimeError(f"repository is not initialized at its exact path: {repo_bytes!r}")
    constellation_root = repository_top(repo_bytes)
    return inspect_cleanliness(
        constellation_root,
        scope,
        expected_head,
        set(),
        constellation_root,
        0,
    )

def format_cleanliness_findings(findings):
    return "\n".join(
        "  repository={} kind={} path={} detail={} action={}".format(
            repr(scope),
            kind.decode("ascii"),
            repr(path),
            repr(detail),
            repr(action),
        )
        for scope, kind, path, detail, action in findings
    )

def read_validated_lock_bytes():
    with open(lock_path, "rb") as source:
        current = source.read(1_048_577)
    if len(current) > 1_048_576:
        raise RuntimeError("constellation lock exceeds the 1 MiB parser bound")
    if current != lock_bytes \
            or hashlib.sha256(current).hexdigest() != validated_lock_sha256:
        raise RuntimeError("constellation lock moved during content snapshot")
    return current

def capture_constellation():
    captured_lock = read_validated_lock_bytes()
    frames = list(capture_repository(root, b"root", b"", set(), 0))
    frames.append((b"constellation-lock", captured_lock))
    sibling_observations = []
    for row in sorted(document["libraries"], key=lambda candidate: candidate["lib"]):
        sibling = os.path.join(parent, row["lib"])
        expected_head = row["git_head"].encode("ascii")
        scope = b"sibling:" + os.fsencode(row["lib"])
        cleanliness, findings = inspect_cleanliness_tree(
            sibling,
            expected_head,
            scope,
        )
        actual = read_git_line(sibling, "rev-parse", "HEAD")
        tree = read_git_line(
            sibling,
            "rev-parse",
            f"{row['git_head']}^{{tree}}",
        )
        frames.extend(
            (
                (b"sibling-lib", row["lib"].encode()),
                (b"sibling-version", row["version"].encode()),
                (b"sibling-expected-head", expected_head),
                (b"sibling-actual-head", actual),
                (b"sibling-tree", tree),
                (b"sibling-remote", row["remote"].encode()),
            )
        )
        sibling_observations.append(
            (
                row["lib"],
                row["git_head"],
                cleanliness,
                tuple(findings),
                actual,
                tree,
            )
        )
    confirmed_lock = read_validated_lock_bytes()
    if confirmed_lock != captured_lock:
        raise RuntimeError("constellation lock moved during complete snapshot capture")
    return tuple(frames), tuple(sibling_observations), captured_lock

try:
    first_capture = capture_constellation()
    confirmed_capture = capture_constellation()
    if confirmed_capture != first_capture:
        raise RuntimeError(
            "constellation moved between complete root-plus-sibling snapshot observations"
        )
    frames, sibling_observations, _captured_lock = first_capture
    for library, expected_head, _cleanliness, findings, actual, _tree in sibling_observations:
        if findings:
            raise RuntimeError(
                f"constellation sibling {library!r} is not recursively clean at pinned "
                f"head {expected_head}:\n"
                + format_cleanliness_findings(findings)
            )
        if actual != expected_head.encode("ascii"):
            raise RuntimeError(
                f"constellation sibling {library!r} is at {actual!r}, lock pins "
                f"{expected_head}"
            )

    add("identity-domain", "org.frankensim.ci.content-snapshot.v3")
    add("identity-version", (3).to_bytes(4, "big"))
    add("schema", "frankensim-ci-content-snapshot-v3")
    for label, data in frames:
        add(label, data)
    read_validated_lock_bytes()
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

destination_is_ordinary_or_absent() {
    python3 -I -c '
import os
import stat
import sys

path = os.fsencode(sys.argv[1])
try:
    metadata = os.lstat(path)
except FileNotFoundError:
    raise SystemExit(0)
is_reparse_point = False
if os.name == "nt":
    attributes = getattr(metadata, "st_file_attributes", None)
    if attributes is None:
        print(f"FATAL: Windows metadata lacks reparse attributes for {path!r}", file=sys.stderr)
        raise SystemExit(1)
    is_reparse_point = bool(attributes & 0x400)
if stat.S_ISLNK(metadata.st_mode) or is_reparse_point or not stat.S_ISDIR(metadata.st_mode):
    print(f"FATAL: bootstrap destination is not an ordinary directory: {path!r}", file=sys.stderr)
    raise SystemExit(1)
' "$1"
}

is_git_checkout() {
    local dir="$1"
    [[ -d "$dir" && ! -L "$dir" ]] || return 1
    working_tree_status "$dir" >/dev/null 2>&1
}

has_bootstrap_marker() {
    [[ "$(git_at "$1" config --local --get "$bootstrap_marker_key" 2>/dev/null || true)" == "true" ]]
}

directory_is_empty() {
    local directory="$1" first_entry
    if ! first_entry="$(find "$directory" -mindepth 1 -maxdepth 1 -print -quit)"; then
        echo "FATAL: cannot enumerate bootstrap destination $directory" >&2
        return 1
    fi
    [[ -z "$first_entry" ]]
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
        printf 'FATAL: incomplete checkout %s has worktree changes; refusing to overwrite it:\n%s\n' \
            "$dir" "$status" >&2
        return 1
    fi
    git_at "$dir" config --local "$bootstrap_marker_key" true || return 1
    if origin="$(git_at "$dir" remote get-url origin 2>/dev/null)"; then
        if [[ "$origin" != "$remote" ]]; then
            echo "FATAL: incomplete checkout $dir origin is $origin, expected $remote" >&2
            return 1
        fi
    else
        git_at "$dir" remote add origin "$remote" || return 1
    fi
    git_at "$dir" fetch --no-auto-maintenance --no-recurse-submodules --quiet --depth 1 origin "$expected" || return 1
    git_at "$dir" checkout --no-recurse-submodules --quiet --no-overwrite-ignore --detach "$expected" || return 1
    verify_clean "$dir" "$expected" || return 1
    if ! git_at "$dir" config --local --unset-all "$bootstrap_marker_key" >/dev/null 2>&1; then
        echo "FATAL: verified checkout $dir retained its incomplete marker" >&2
        return 1
    fi
}

status=0
while IFS=$'\t' read -r lib _version head remote; do
    [[ -n "$lib" ]] || continue
    dir="$parent/$lib"
    if ! destination_is_ordinary_or_absent "$dir"; then
        json_row "$lib" "refused" "$head" "<non-git>" "destination is a link, reparse point, or non-directory"
        status=1
        continue
    fi
    if [[ -e "$dir/.git" || -L "$dir/.git" ]]; then
        if ! preflight="$(working_tree_status "$dir" 2>&1)"; then
            json_row "$lib" "refused" "$head" "<unreadable>" "repository preflight refused local authority or redirection"
            printf 'FATAL: repository preflight refused %s:\n%s\n' "$dir" "$preflight" >&2
            status=1
            continue
        fi
    fi
    if is_git_checkout "$dir"; then
        have="$(git_at "$dir" rev-parse HEAD 2>/dev/null || true)"
        if [[ "$have" == "$head" ]]; then
            if verify_clean "$dir" "$head"; then
                if has_bootstrap_marker "$dir" \
                    && ! git_at "$dir" config --local --unset-all "$bootstrap_marker_key" >/dev/null 2>&1; then
                    json_row "$lib" "refused" "$head" "$have" "verified checkout retained incomplete marker"
                    echo "FATAL: verified checkout $dir retained its incomplete marker" >&2
                    status=1
                    continue
                fi
                json_row "$lib" "verified" "$head" "$have" "pinned head and clean tree"
            else
                json_row "$lib" "refused" "$head" "$have" \
                    "checkout is at the pinned head but is dirty or unreadable"
                echo "FATAL: required constellation sibling $dir is pinned but not clean" >&2
                status=1
            fi
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
            origin="$(git_at "$dir" remote get-url origin 2>/dev/null || true)"
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
        || { [[ -d "$dir" ]] && ! directory_is_empty "$dir"; }; then
        json_row "$lib" "refused" "$head" "<non-git>" "non-empty or unreadable path is not a git checkout"
        echo "FATAL: refusing to initialize non-empty or unreadable non-git directory $dir" >&2
        status=1
        continue
    fi
    mkdir -p "$dir"
    if ! git_at "$dir" init --quiet --template= \
        || ! git_at "$dir" config --local "$bootstrap_marker_key" true \
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
