#!/usr/bin/env bash
# M8-W6a: fetch or verify the reference recordings named by references.json.
#   --verify (default)  verify present files; fetch absent ones, checking the remote bytes'
#                       sha256 and length BEFORE they are moved into place. Any mismatch refuses.
#   --pin               print the observed sha256 and bytes of every entry for review; writes nothing.
# Never overwrites and never deletes. Audio is never committed: FRANKENSIM_REFERENCE_DIR must be an
# existing directory, and a directory on the Mac Data volume (/System/Volumes/Data) is refused.
# One JSONL line per file on stdout: utc, event, id, url, sha256_expected, sha256_observed, bytes,
# dest_mount, elapsed_ms, host. FRANKENSIM_REFERENCES_MANIFEST and FRANKENSIM_TEST_MOUNT exist for
# scripts/tests/test_instrument_references.py only.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
mode="${1:---verify}"
case "$mode" in
--verify | --pin) ;;
*)
    echo "usage: fetch-references.sh [--verify|--pin]" >&2
    exit 2
    ;;
esac
dest="${FRANKENSIM_REFERENCE_DIR:?set FRANKENSIM_REFERENCE_DIR to an existing directory off the Mac Data volume}"
[ -d "$dest" ] || {
    echo "REFUSED: FRANKENSIM_REFERENCE_DIR is not an existing directory: $dest" >&2
    exit 1
}
mount="${FRANKENSIM_TEST_MOUNT:-$(df -P "$dest" | awk 'NR == 2 {print $6}')}"
if [ "$mount" = "/System/Volumes/Data" ]; then
    echo "REFUSED: $dest is on the Mac Data volume ($mount); reference audio lives off-volume" >&2
    exit 1
fi
manifest="${FRANKENSIM_REFERENCES_MANIFEST:-$here/references.json}"

exec python3 -I - "$mode" "$manifest" "$dest" "$mount" <<'PY'
import datetime
import hashlib
import json
import os
import socket
import sys
import time
import urllib.request

mode, manifest, dest, mount = sys.argv[1:]
USER_AGENT = "frankensim-agent (+https://github.com/Dicklesworthstone/frankensim)"
host = socket.gethostname()
refused = 0


def emit(event, entry, observed, size, start):
    print(json.dumps({
        "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "event": event, "id": entry["id"], "url": entry["url"],
        "sha256_expected": entry["sha256"], "sha256_observed": observed, "bytes": size,
        "dest_mount": mount, "elapsed_ms": round((time.monotonic() - start) * 1000.0, 3), "host": host,
    }), flush=True)


def digest_file(path):
    h, size = hashlib.sha256(), 0
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
            size += len(block)
    return h.hexdigest(), size


def download(url, part):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    h, size = hashlib.sha256(), 0
    with urllib.request.urlopen(request, timeout=120) as response, open(part, "wb") as out:
        for block in iter(lambda: response.read(1 << 20), b""):
            h.update(block)
            size += len(block)
            out.write(block)
    return h.hexdigest(), size


files = json.load(open(manifest))["files"]
for entry in files:
    start = time.monotonic()
    parts = entry["path"].split("/")
    if os.path.isabs(entry["path"]) or any(p in ("", ".", "..") for p in parts):
        print(f"REFUSED {entry['id']}: storage path escapes the destination: {entry['path']!r}", file=sys.stderr)
        emit("refused", entry, None, None, start)
        refused += 1
        continue
    path = os.path.join(dest, *parts)
    good = lambda d, n: d == entry["sha256"] and n == entry["bytes"]
    if os.path.exists(path):
        observed, size = digest_file(path)
        if mode == "--pin":
            emit("pinned", entry, observed, size, start)
        elif good(observed, size):
            emit("verified", entry, observed, size, start)
        else:
            emit("refused", entry, observed, size, start)
            refused += 1
        continue
    if mode == "--pin":
        observed, size = download(entry["url"], os.devnull)
        emit("pinned", entry, observed, size, start)
        continue
    os.makedirs(os.path.dirname(path), exist_ok=True)
    part = path + ".part"
    try:
        observed, size = download(entry["url"], part)
    except OSError as error:
        print(f"REFUSED fetch {entry['id']}: {error}", file=sys.stderr)
        emit("refused", entry, None, None, start)
        refused += 1
        continue
    if not good(observed, size):
        # The unverified bytes stay at .part for inspection; nothing is deleted.
        print(f"REFUSED remote {entry['id']}: {observed} {size} B kept at {part}", file=sys.stderr)
        emit("refused", entry, observed, size, start)
        refused += 1
        continue
    try:
        os.link(part, path)  # refuses, never overwrites, if the path appeared meanwhile
    except FileExistsError:
        print(f"REFUSED {entry['id']}: {path} appeared during the fetch", file=sys.stderr)
        emit("refused", entry, observed, size, start)
        refused += 1
        continue
    os.unlink(part)  # drops the temporary name only; the verified bytes stay at path
    retained, retained_size = digest_file(path)
    if not good(retained, retained_size):
        emit("refused", entry, retained, retained_size, start)
        refused += 1
        continue
    emit("fetched", entry, retained, retained_size, start)
sys.exit(1 if refused else 0)
PY
