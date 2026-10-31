#!/usr/bin/env bash
# M8-W6a: fetch or verify the reference recordings named by references.json.
#   --verify (default)  verify present files; fetch absent ones, checking the remote bytes'
#                       sha256 and length BEFORE they are moved into place. Any mismatch refuses.
#   --pin               print the observed sha256 and bytes of every entry for review; writes nothing.
#   --instrument NAME   select an instrument (repeatable, case-insensitive exact name).
#   --id ID             select a recording (repeatable); combined selectors intersect.
#   --list              print selected manifest entries without downloading or requiring a destination.
#   --manifest PATH     use another reference manifest with the same file fields.
#   --members           acquire checksum-pinned archive members. Tar sources stream without being
#                       retained; ZIP sources are cached for random access, without a full unzip.
#   --zip-ranges        with --members, read ZIP members by HTTP ranges instead of fetching the
#                       archive. Verifies each member; does NOT verify the complete archive hash.
# Never overwrites and never deletes. Audio is never committed: FRANKENSIM_REFERENCE_DIR must be an
# existing directory, and a directory on the Mac Data volume (/System/Volumes/Data) is refused.
# One JSONL line per file on stdout: utc, event, id, url, sha256_expected, sha256_observed, bytes,
# dest_mount, elapsed_ms, host. No selectors means the complete corpus (over 12 GB).
# FRANKENSIM_TEST_MOUNT exists for scripts/tests/test_instrument_references.py only.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
dest="${FRANKENSIM_REFERENCE_DIR:-}"
mount="${FRANKENSIM_TEST_MOUNT:-}"
if [ -d "$dest" ] && [ -z "$mount" ]; then
    mount="$(df -P "$dest" | awk 'NR == 2 {print $6}')"
fi
manifest="${FRANKENSIM_REFERENCES_MANIFEST:-$here/references.json}"

exec python3 -I - "$manifest" "$dest" "$mount" "$@" <<'PY'
import argparse
import datetime
import hashlib
import io
import json
import os
import socket
import stat
import sys
import tarfile
import time
import urllib.request
import zipfile

manifest, dest, mount = sys.argv[1:4]
parser = argparse.ArgumentParser(description="Fetch checksum-pinned real instrument recordings; no selectors downloads the complete corpus.")
modes = parser.add_mutually_exclusive_group()
modes.add_argument("--verify", action="store_true", help="verify present files and fetch absent ones (default)")
modes.add_argument("--pin", action="store_true", help="report remote or local digests without retaining downloads")
parser.add_argument("--instrument", action="append", default=[], help="exact instrument name, repeatable and case-insensitive")
parser.add_argument("--id", action="append", default=[], help="exact recording ID, repeatable")
parser.add_argument("--list", action="store_true", help="list selected manifest entries without downloading")
parser.add_argument("--manifest", default=manifest, help="reference manifest JSON")
parser.add_argument("--members", action="store_true", help="extract only pinned archive members; tar streams, ZIP caches the source for random access")
parser.add_argument("--zip-ranges", action="store_true", help="with --members, fetch ZIP ranges and verify member pins only, not the complete archive hash")
args = parser.parse_args(sys.argv[4:])
if args.zip_ranges and (not args.members or args.pin):
    parser.error("--zip-ranges requires --members and existing member pins (not --pin)")
mode = "--pin" if args.pin else "--verify"
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
    with urllib.request.urlopen(request, timeout=120) as response, open(part, "wb" if part == os.devnull else "xb") as out:
        for block in iter(lambda: response.read(1 << 20), b""):
            h.update(block)
            size += len(block)
            out.write(block)
    return h.hexdigest(), size


def storage_path(relative):
    parts = relative.split("/")
    if os.path.isabs(relative) or any(p in ("", ".", "..") for p in parts):
        raise ValueError(f"storage path escapes the destination: {relative!r}")
    path = os.path.join(dest, *parts)
    if os.path.commonpath((os.path.realpath(dest), os.path.realpath(path))) != os.path.realpath(dest):
        raise ValueError(f"storage path escapes through a symlink: {relative!r}")
    return path


class HashedReader:
    """Hash compressed source bytes, including padding after tar's end marker."""
    def __init__(self, source):
        self.source, self.digest, self.size = source, hashlib.sha256(), 0

    def read(self, count=-1):
        data = self.source.read(count)
        self.digest.update(data)
        self.size += len(data)
        return data


class RangedZip(io.RawIOBase):
    """Seek within a public ZIP, refusing servers that ignore or alter a range.

    A successful read binds only the requested bytes. Member SHA-256 checks,
    not ZIP metadata or the source's published MD5, authorize publication.
    """
    def __init__(self, entry, member_bytes):
        self.url, self.size = entry["url"], entry["bytes"]
        self.position, self.transferred = 0, 0
        self.budget = 32 * (1 << 20) + 2 * member_bytes

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        if whence not in (0, 1, 2):
            raise ValueError("invalid ZIP seek origin")
        position = (0, self.position, self.size)[whence] + offset
        if position < 0 or position > self.size:
            raise ValueError("ZIP seek outside declared source length")
        self.position = position
        return position

    def read(self, count=-1):
        count = min(self.size - self.position, self.size if count < 0 else count)
        if count == 0:
            return b""
        if count > 16 * (1 << 20) or self.transferred + count > self.budget:
            raise ValueError("ZIP range acquisition exceeds its bounded transfer budget")
        end = self.position + count - 1
        requested = f"bytes {self.position}-{end}/{self.size}"
        request = urllib.request.Request(self.url, headers={
            "User-Agent": USER_AGENT, "Accept-Encoding": "identity",
            "Range": f"bytes={self.position}-{end}",
        })
        with urllib.request.urlopen(request, timeout=120) as response:
            if response.status != 206 or response.headers.get("Content-Range") != requested:
                raise ValueError("source did not honor the exact ZIP byte range")
            data = response.read(count + 1)
        if len(data) != count:
            raise ValueError("ZIP byte range was truncated or oversized")
        self.position += count
        self.transferred += count
        return data


def fetch_members(entry, start):
    members = entry.get("members", [])
    if not members:
        raise ValueError("archive has no checksum-pinned members")
    targets, paths = {}, set()
    for member in members:
        name = member["archive_member"]
        path = storage_path(member["path"])
        if name in targets or path in paths:
            raise ValueError("duplicate archive member or destination")
        paths.add(path)
        targets[name] = dict(member, id=entry["id"] + ":" + name, url=entry["url"])
    pending = {}
    for name, member in targets.items():
        path = storage_path(member["path"])
        if os.path.exists(path):
            observed, size = digest_file(path)
            if mode != "--pin" and (observed != member["sha256"] or size != member["bytes"]):
                raise ValueError(f"existing member does not match its pin: {member['path']}")
            emit("pinned" if mode == "--pin" else "verified", member, observed, size, start)
        else:
            pending[name] = member
    if not pending:
        return

    archive_path = storage_path(entry["path"])
    staged, seen = [], set()

    def stage(name, member_size, regular, open_audio):
        member = pending[name]
        if name in seen or not regular:
            raise ValueError(f"duplicate or non-regular archive member: {name}")
        seen.add(name)
        if mode != "--pin" and member_size != member["bytes"]:
            raise ValueError(f"archive member length does not match its pin: {name}")
        path = storage_path(member["path"])
        part = path + ".part"
        if mode != "--pin":
            os.makedirs(os.path.dirname(path), exist_ok=True)
        digest, size = hashlib.sha256(), 0
        # Exclusive creation preserves incomplete earlier acquisitions for inspection.
        with open_audio() as audio, open(os.devnull if mode == "--pin" else part, "wb" if mode == "--pin" else "xb") as out:
            for block in iter(lambda: audio.read(1 << 20), b""):
                digest.update(block)
                size += len(block)
                out.write(block)
        observed = digest.hexdigest()
        if mode != "--pin" and (observed != member["sha256"] or size != member["bytes"]):
            raise ValueError(f"archive member does not match its pin; bytes kept at {part}")
        staged.append((member, path, part, observed, size))

    archive_format = entry.get("archive_format", "tar")
    ranged = args.zip_ranges and archive_format == "zip"
    if ranged:
        reader = RangedZip(entry, sum(member["bytes"] for member in pending.values()))
        with zipfile.ZipFile(reader) as archive:
            for item in archive.infolist():
                if item.filename in pending:
                    file_type = stat.S_IFMT(item.external_attr >> 16)
                    regular = not item.is_dir() and file_type in (0, stat.S_IFREG)
                    stage(item.filename, item.file_size, regular, lambda: archive.open(item))
        observed, source_size = None, reader.transferred
    elif not entry.get("sha256"):
        raise ValueError("source has no complete SHA-256 pin; use --members --zip-ranges for pinned ZIP members")
    elif archive_format == "zip":
        # ZIP's directory lives at the end. Keep the verified compressed source on the
        # reference volume; never unpack its unselected files or buffer it all in RAM.
        if os.path.exists(archive_path):
            observed, source_size = digest_file(archive_path)
        else:
            if mode == "--pin":
                raise ValueError("ZIP --pin --members requires a local source archive for random access")
            os.makedirs(os.path.dirname(archive_path), exist_ok=True)
            archive_part = archive_path + ".part"
            observed, source_size = download(entry["url"], archive_part)
            if observed != entry["sha256"] or source_size != entry["bytes"]:
                raise ValueError(f"source archive does not match its pin; bytes kept at {archive_part}")
            os.link(archive_part, archive_path)
            os.unlink(archive_part)  # only the temporary name; verified source remains
        if mode != "--pin" and (observed != entry["sha256"] or source_size != entry["bytes"]):
            raise ValueError("source archive does not match its pin")
        with zipfile.ZipFile(archive_path) as archive:
            for item in archive.infolist():
                if item.filename in pending:
                    file_type = stat.S_IFMT(item.external_attr >> 16)
                    regular = not item.is_dir() and file_type in (0, stat.S_IFREG)
                    stage(item.filename, item.file_size, regular, lambda: archive.open(item))
    elif archive_format == "tar":
        request = urllib.request.Request(entry["url"], headers={"User-Agent": USER_AGENT})
        source = open(archive_path, "rb") if os.path.exists(archive_path) else urllib.request.urlopen(request, timeout=120)
        with source:
            reader = HashedReader(source)
            with tarfile.open(fileobj=reader, mode="r|*") as archive:
                for item in archive:
                    if item.name in pending:
                        stage(item.name, item.size, item.isfile(), lambda: archive.extractfile(item))
            while reader.read(1 << 20):
                pass
        observed, source_size = reader.digest.hexdigest(), reader.size
    else:
        raise ValueError(f"unsupported archive format: {archive_format}")
    if seen != set(pending):
        raise ValueError("archive is missing members: " + ", ".join(sorted(set(pending) - seen)))
    if not ranged and mode != "--pin" and (observed != entry["sha256"] or source_size != entry["bytes"]):
        raise ValueError("source archive does not match its pin; unverified members remain at .part")
    emit("archive-ranges-read" if ranged else "pinned" if mode == "--pin" else "archive-verified", entry, observed, source_size, start)
    # Range mode verifies every member; other modes also verify the complete source.
    for member, path, part, observed, size in staged:
        if mode != "--pin":
            os.link(part, path)
            os.unlink(part)  # only the temporary name; verified bytes remain at path
        emit("pinned" if mode == "--pin" else "fetched", member, observed, size, start)


with open(args.manifest) as source:
    files = json.load(source)["files"]
instruments = {f.get("instrument", "").casefold() for f in files}
ids = {f["id"] for f in files}
unknown = [name for name in args.instrument if name.casefold() not in instruments]
unknown += [name for name in args.id if name not in ids]
if unknown:
    parser.error("unknown recording selector(s): " + ", ".join(unknown))
selected_instruments = {name.casefold() for name in args.instrument}
selected_ids = set(args.id)
files = [f for f in files
         if (not selected_instruments or f.get("instrument", "").casefold() in selected_instruments)
         and (not selected_ids or f["id"] in selected_ids)]
if not files:
    parser.error("selectors match no recordings")
if args.list:
    for entry in files:
        print(json.dumps(entry))
    sys.exit(0)
if not os.path.isdir(dest):
    sys.exit("REFUSED: set FRANKENSIM_REFERENCE_DIR to an existing directory off the Mac Data volume")
if mount == "/System/Volumes/Data":
    sys.exit(f"REFUSED: {dest} is on the Mac Data volume ({mount}); reference audio lives off-volume")
for entry in files:
    start = time.monotonic()
    if args.members and entry.get("kind") == "archive":
        try:
            fetch_members(entry, start)
        except (OSError, ValueError, tarfile.TarError, zipfile.BadZipFile, RuntimeError) as error:
            print(f"REFUSED archive {entry['id']}: {error}", file=sys.stderr)
            emit("refused", entry, None, None, start)
            refused += 1
        continue
    if not entry.get("sha256"):
        print(f"REFUSED {entry['id']}: source has no complete SHA-256 pin; use --members --zip-ranges", file=sys.stderr)
        emit("refused", entry, None, None, start)
        refused += 1
        continue
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
