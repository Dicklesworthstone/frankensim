"""M8-W6a: data/instruments/references.json and fetch-references.sh (Python stdlib only).

Run: python3 -I scripts/tests/test_instrument_references.py -v
Every test prints one IRPDATA line for scripts/e2e/instrument_data.sh. Artifact digests are sha256
(the stdlib has no BLAKE3). Tolerance everywhere: an exact sha256 match and an exact byte count.
"""
import contextlib
import functools
import hashlib
import http.server
import io
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
import threading
import time
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "data/instruments/references.json"
FETCH = ROOT / "data/instruments/fetch-references.sh"
LINKS = ROOT / "apps/instrument-lab/tools/instrgrade/uiowa-links.json"
LEDGER = ROOT / "data/instruments/ledger.tsv"
BEAD = "frankensim-irp-v5-instrument-realism-pexfp.217"
FIELDS = ("id", "link", "path", "page", "url", "sha256", "bytes", "kind", "instrument", "dynamic", "range",
          "chain", "licence", "ledger_source_id", "channels", "pair", "gradable")
SHA256 = re.compile(r"^[0-9a-f]{64}$")


@functools.lru_cache(maxsize=None)
def manifest():
    return json.loads(MANIFEST.read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def irpdata(test, rows, refusals, observed, start, source_ids=("uiowa-mis",)):
    artifacts = [{"path": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in (MANIFEST, FETCH)]
    print("IRPDATA " + json.dumps({
        "bead": BEAD, "test": "ReferenceTests." + test, "source_ids": list(source_ids),
        "artifacts": artifacts, "rows_checked": rows, "refusals_exercised": refusals, "ladder": None,
        "tolerance": "exact sha256 and byte count", "observed": observed, "verdict": "pass",
        "elapsed_ms": round((time.monotonic() - start) * 1000.0, 3),
    }), flush=True)


def ledger_rows():
    lines = LEDGER.read_text().splitlines()
    header = lines[1].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[2:] if line.strip()]


def fetch(directory, entries, mount="/fake-external", mode=None, selectors=()):
    """Run the real script against a temporary manifest; FRANKENSIM_TEST_MOUNT stands in for df."""
    test_manifest = Path(directory) / "manifest.json"
    test_manifest.write_text(json.dumps({"files": entries}))
    env = dict(os.environ, FRANKENSIM_REFERENCE_DIR=str(Path(directory) / "refs"),
               FRANKENSIM_REFERENCES_MANIFEST=str(test_manifest), FRANKENSIM_TEST_MOUNT=mount)
    (Path(directory) / "refs").mkdir(exist_ok=True)
    command = ["bash", str(FETCH)] + ([mode] if mode else []) + list(selectors)
    result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=120)
    events = [json.loads(line) for line in result.stdout.splitlines() if line.startswith("{")]
    return result, events


def entry(payload, url="http://127.0.0.1:9/unreachable"):
    return {"id": "test-entry", "link": "sub/dir/a file.aiff", "path": "sub/dir/a file.aiff", "url": url,
            "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}


def archive_entry(directory, compression=""):
    source = Path(directory) / "source.tar"
    selected = b"original audio bytes"
    with tarfile.open(source, "w" + compression) as archive:
        for name, payload in (("samples/selected.wav", selected), ("other.wav", b"unselected")):
            item = tarfile.TarInfo(name)
            item.size = len(payload)
            archive.addfile(item, io.BytesIO(payload))
    result = dict(entry(source.read_bytes(), source.as_uri()), kind="archive", path="container.tar")
    result["members"] = [{"archive_member": "samples/selected.wav", "path": "selected/note.wav",
                          "sha256": hashlib.sha256(selected).hexdigest(), "bytes": len(selected)}]
    return result, selected


def zip_entry(directory):
    source = Path(directory) / "source.zip"
    payload = b"original violin recording"
    with zipfile.ZipFile(source, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("notes/A4.wav", payload)
        archive.writestr("unselected.wav", b"not requested")
    result = dict(entry(source.read_bytes(), source.as_uri()), kind="archive",
                  archive_format="zip", path="container.zip")
    result["members"] = [{"archive_member": "notes/A4.wav", "path": "selected/A4.wav",
                          "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}]
    return result, payload


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


@contextlib.contextmanager
def range_server(payload, fault=None):
    class Ranges(Quiet):
        def do_GET(self):
            first, last = map(int, self.headers["Range"].removeprefix("bytes=").split("-"))
            body = payload[first:last + 1]
            self.send_response(200 if fault == "ignored" else 206)
            self.send_header("Content-Range", f"bytes {first + (fault == 'shifted')}-{last}/{len(payload)}")
            if fault == "truncated":
                body = body[:-1]
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Ranges)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        yield f"http://127.0.0.1:{server.server_port}/source.zip"
    finally:
        server.shutdown()
        server.server_close()


def tree(directory):
    return sorted((str(p.relative_to(directory)), p.read_bytes()) for p in Path(directory).rglob("*") if p.is_file())


class ReferenceTests(unittest.TestCase):
    def test_zip_ranges_verify_members_without_claiming_full_archive_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, payload = zip_entry(directory)
            archive["sha256"] = None  # Source publishes MD5; no complete SHA-256 was measured.
            with range_server((Path(directory) / "source.zip").read_bytes()) as url:
                archive["url"] = url
                result, events = fetch(directory, [archive], selectors=("--members", "--zip-ranges"))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([e["event"] for e in events], ["archive-ranges-read", "fetched"])
            self.assertIsNone(events[0]["sha256_observed"])
            self.assertEqual(events[1]["sha256_observed"], archive["members"][0]["sha256"])
            self.assertEqual(tree(Path(directory) / "refs"), [("selected/A4.wav", payload)])
            archive["url"] = "http://127.0.0.1:9/unreachable"
            result, events = fetch(directory, [archive], selectors=("--members", "--zip-ranges"))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([e["event"] for e in events], ["verified"])

    def test_zip_ranges_refuse_bad_ranges_and_member_pins(self):
        for fault in ("ignored", "shifted", "truncated", "member-hash"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                archive, _ = zip_entry(directory)
                if fault == "member-hash":
                    archive["members"][0]["sha256"] = "0" * 64
                with range_server((Path(directory) / "source.zip").read_bytes(), fault) as url:
                    archive["url"] = url
                    result, events = fetch(directory, [archive], selectors=("--members", "--zip-ranges"))
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual([e["event"] for e in events], ["refused"])
                self.assertFalse((Path(directory) / "refs/selected/A4.wav").exists())

    def test_unpinned_source_cannot_trigger_a_complete_archive_download(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, _ = zip_entry(directory)
            archive["sha256"] = None
            archive["url"] = "http://127.0.0.1:9/unreachable"
            for selectors in ((), ("--members",)):
                result, events = fetch(directory, [archive], selectors=selectors)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn("source has no complete SHA-256 pin", result.stderr)
                self.assertEqual([e["event"] for e in events], ["refused"])
                self.assertEqual(tree(Path(directory) / "refs"), [])

    def test_zip_caches_verified_source_and_extracts_only_selected_members(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, payload = zip_entry(directory)
            result, events = fetch(directory, [archive], selectors=("--members",))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([e["event"] for e in events], ["archive-verified", "fetched"])
            refs = Path(directory) / "refs"
            self.assertEqual(tree(refs), [("container.zip", (Path(directory) / "source.zip").read_bytes()),
                                          ("selected/A4.wav", payload)])
            archive["url"] = "http://127.0.0.1:9/unreachable"
            result, events = fetch(directory, [archive], selectors=("--members",))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([e["event"] for e in events], ["verified"])

    def test_zip_refuses_corrupt_source_or_member_before_publication(self):
        for bad_target in ("archive", "member"):
            with self.subTest(bad_target=bad_target), tempfile.TemporaryDirectory() as directory:
                archive, _ = zip_entry(directory)
                (archive if bad_target == "archive" else archive["members"][0])["sha256"] = "0" * 64
                result, events = fetch(directory, [archive], selectors=("--members",))
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual([e["event"] for e in events], ["refused"])
                self.assertFalse((Path(directory) / "refs/selected/A4.wav").exists())

    def test_zip_pin_requires_cached_source_and_does_not_write_members(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, _ = zip_entry(directory)
            result, _ = fetch(directory, [archive], mode="--pin", selectors=("--members",))
            self.assertEqual(result.returncode, 1)
            self.assertIn("requires a local source archive", result.stderr)
            refs = Path(directory) / "refs"
            self.assertEqual(tree(refs), [])
            (refs / "container.zip").write_bytes((Path(directory) / "source.zip").read_bytes())
            before = tree(refs)
            result, events = fetch(directory, [archive], mode="--pin", selectors=("--members",))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([e["event"] for e in events], ["pinned", "pinned"])
            self.assertEqual(tree(refs), before)

    def test_zip_preserves_existing_partial_download(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, _ = zip_entry(directory)
            refs = Path(directory) / "refs"
            refs.mkdir()
            (refs / "container.zip.part").write_bytes(b"earlier incomplete download")
            before = tree(refs)
            result, events = fetch(directory, [archive], selectors=("--members",))
            self.assertEqual(result.returncode, 1)
            self.assertEqual([e["event"] for e in events], ["refused"])
            self.assertEqual(tree(refs), before)

    def test_archive_members_stream_tar_xz_and_bzip_without_retaining_container(self):
        for compression in ("", ":xz", ":bz2"):
            with self.subTest(compression=compression), tempfile.TemporaryDirectory() as directory:
                archive, payload = archive_entry(directory, compression)
                result, events = fetch(directory, [archive], selectors=("--members",))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual([e["event"] for e in events], ["archive-verified", "fetched"])
                self.assertEqual(tree(Path(directory) / "refs"), [("selected/note.wav", payload)])
                # Pinned extracted audio is usable without another archive download.
                archive["url"] = "http://127.0.0.1:9/unreachable"
                result, events = fetch(directory, [archive], selectors=("--members",))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual([e["event"] for e in events], ["verified"])

    def test_archive_must_verify_before_any_member_is_published(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, payload = archive_entry(directory)
            archive["sha256"] = "0" * 64
            result, events = fetch(directory, [archive], selectors=("--members",))
            self.assertEqual(result.returncode, 1)
            self.assertEqual([e["event"] for e in events], ["refused"])
            self.assertIn("source archive does not match", result.stderr)
            self.assertEqual(tree(Path(directory) / "refs"), [("selected/note.wav.part", payload)])

    def test_archive_member_pin_and_missing_member_are_checked(self):
        for defect in ("sha256", "archive_member"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                archive, _ = archive_entry(directory)
                archive["members"][0][defect] = "0" * 64 if defect == "sha256" else "missing.wav"
                result, events = fetch(directory, [archive], selectors=("--members",))
                self.assertEqual(result.returncode, 1)
                self.assertEqual([e["event"] for e in events], ["refused"])
                self.assertFalse((Path(directory) / "refs/selected/note.wav").exists())

    def test_archive_member_destinations_must_not_escape_or_alias(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, _ = archive_entry(directory)
            archive["members"][0]["path"] = "../escape.wav"
            result, events = fetch(directory, [archive], selectors=("--members",))
            self.assertEqual(result.returncode, 1)
            self.assertIn("storage path escapes", result.stderr)
            self.assertEqual(tree(Path(directory) / "refs"), [])
            member = dict(archive["members"][0], path="note.wav")
            archive["members"] = [member, dict(member, archive_member="other.wav")]
            result, events = fetch(directory, [archive], selectors=("--members",))
            self.assertEqual(result.returncode, 1)
            self.assertIn("duplicate", result.stderr)
            self.assertEqual(tree(Path(directory) / "refs"), [])

    def test_archive_member_pin_mode_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, _ = archive_entry(directory, ":bz2")
            result, events = fetch(directory, [archive], mode="--pin", selectors=("--members",))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([e["event"] for e in events], ["pinned", "pinned"])
            self.assertEqual(events[0]["sha256_observed"], archive["sha256"])
            self.assertEqual(events[1]["sha256_observed"], archive["members"][0]["sha256"])
            self.assertEqual(tree(Path(directory) / "refs"), [])

    def test_list_selection_needs_no_audio_destination(self):
        env = dict(os.environ, FRANKENSIM_REFERENCE_DIR="/nonexistent-reference-directory")
        result = subprocess.run(["bash", str(FETCH), "--list", "--instrument", "violin"],
                                env=env, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        entries = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertGreater(len(entries), 0)
        self.assertEqual({entry["instrument"] for entry in entries}, {"Violin"})

    def test_id_selection_does_not_fetch_unselected_files(self):
        with tempfile.TemporaryDirectory() as directory:
            good = dict(entry(b"right bytes"), instrument="Violin")
            other = dict(entry(b"unavailable"), id="other", path="other.aiff", instrument="Cello")
            target = Path(directory) / "refs" / good["path"]
            target.parent.mkdir(parents=True)
            target.write_bytes(b"right bytes")
            result, events = fetch(directory, [good, other], selectors=("--id", good["id"]))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([event["id"] for event in events], [good["id"]])
            self.assertFalse((Path(directory) / "refs" / other["path"]).exists())

    def test_bad_or_disjoint_selectors_refuse_before_fetch(self):
        with tempfile.TemporaryDirectory() as directory:
            entries = [dict(entry(b"x"), instrument="Violin"),
                       dict(entry(b"y"), id="cello", path="cello.aiff", instrument="Cello")]
            for selectors in (("--instrument", "misspelt"), ("--id", "missing"),
                              ("--instrument", "Violin", "--id", "cello")):
                result, events = fetch(directory, entries, selectors=selectors)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(events, [])
                self.assertEqual(tree(Path(directory) / "refs"), [])

    def test_schema_and_required_fields(self):
        start = time.monotonic()
        m = manifest()
        self.assertEqual(m["schema"], "frankensim.instrument-references.v1")
        licences = {licence["id"] for licence in m["licences"]}
        chains = {chain["id"] for chain in m["chains"]}
        for f in m["files"]:
            missing = [k for k in FIELDS if k not in f]
            self.assertEqual(missing, [], f["id"])
            self.assertRegex(f["sha256"], SHA256, f["id"])
            self.assertIsInstance(f["bytes"], int, f["id"])
            self.assertGreater(f["bytes"], 0, f["id"])
            self.assertIn(f["licence"], licences)
            self.assertIn(f["chain"], chains)
            self.assertIn(f["kind"], ("chromatic-scale", "ambient-silence"))
            self.assertIn(f["channels"], ("mono", "stereo"))
            self.assertIsInstance(f["gradable"], bool)
            self.assertTrue(f["url"].startswith("https://theremin.music.uiowa.edu/"))
        for licence in m["licences"]:
            for k in ("text", "url", "retrieved_utc", "decision"):
                self.assertTrue(licence[k], k)
        irpdata("test_schema_and_required_fields", len(m["files"]), 0, len(m["files"]), start)

    def test_ids_and_urls_unique(self):
        start = time.monotonic()
        files = manifest()["files"]
        for key in ("id", "url", "link", "path"):
            values = [f[key] for f in files]
            self.assertEqual(len(values), len(set(values)), key)
        irpdata("test_ids_and_urls_unique", len(files), 0, len(files), start)

    def test_every_uiowa_link_is_manifested_or_excluded(self):
        start = time.monotonic()
        m = manifest()
        links = json.loads(LINKS.read_text())
        self.assertEqual(m["source_index"]["sha256"], digest(LINKS))
        self.assertEqual(m["source_index"]["links"], sum(len(v) for v in links.values()))
        covered = {f["link"]: f["page"] for f in m["files"]}
        covered.update({x["link"]: x["page"] for x in m["excluded"]})
        for x in m["excluded"]:
            self.assertTrue(x["reason"])
        total = 0
        for page, page_links in links.items():
            for link in page_links:
                total += 1
                self.assertIn(link, covered, link)
                self.assertEqual(covered[link], page, link)
        self.assertEqual(len(covered), m["source_index"]["unique_links"])
        irpdata("test_every_uiowa_link_is_manifested_or_excluded", total, 0, len(covered), start)

    def test_mono_stereo_pairing_is_explicit(self):
        start = time.monotonic()
        files = manifest()["files"]
        by_id = {f["id"]: f for f in files}
        unpaired = {}
        for f in files:
            if not f["page"].endswith("2012.html") and f["page"] != "MISguitar.html":
                continue
            if f["pair"] is None:
                self.assertTrue(f.get("unpaired"), f["id"])
                unpaired[f["page"]] = unpaired.get(f["page"], 0) + 1
            else:
                partner = by_id[f["pair"]]
                self.assertEqual(partner["pair"], f["id"])
                self.assertNotEqual(partner["channels"], f["channels"])
        # The pack's link-level imbalances (guitar 50 stereo/48 mono, viola2012 64/65) are
        # surfaced and explained: guitar's is duplicate hrefs only (its 90 files pair 45/45);
        # viola2012 keeps one genuinely unpaired file.
        links = json.loads(LINKS.read_text())
        pages = {p["page"]: p for p in manifest()["pages"]}
        counts = {}
        for page in ("MISguitar.html", "MISviola2012.html"):
            link_level = {c: sum(f".{c}." in link for link in links[page]) for c in ("mono", "stereo")}
            unique = {c: sum(f["page"] == page and f["channels"] == c for f in files) for c in ("mono", "stereo")}
            duplicates = {c: sum(f".{c}." in link for link in pages[page]["duplicate_links"]) for c in ("mono", "stereo")}
            self.assertNotEqual(link_level["mono"], link_level["stereo"], page)
            self.assertEqual(pages[page]["links"], len(links[page]))
            counts[page] = {"links": link_level, "files": unique, "duplicate_links": duplicates}
        self.assertEqual(counts["MISguitar.html"]["files"], {"mono": 45, "stereo": 45})
        self.assertNotIn("MISguitar.html", unpaired)
        self.assertGreater(unpaired.get("MISviola2012.html", 0), 0)
        irpdata("test_mono_stereo_pairing_is_explicit", len(files), 0, {"counts": counts, "unpaired": unpaired}, start)

    def test_stereo_not_gradable_until_derivation_documented(self):
        start = time.monotonic()
        files = manifest()["files"]
        stereo = [f for f in files if f["channels"] == "stereo"]
        self.assertTrue(stereo)
        for f in stereo:
            if f["chain"] == "uiowa-mis-qtc40-2011":
                self.assertFalse(f["gradable"], f["id"])
                self.assertIn("not documented", f["gradable_reason"])
            else:
                # Stereo-only sessions have no mono partner: the channel choice is declared.
                self.assertIn("undocumented", f["channel_policy"], f["id"])
        irpdata("test_stereo_not_gradable_until_derivation_documented", len(stereo), len(stereo), len(stereo), start)

    def test_chain_assignment_by_era(self):
        start = time.monotonic()
        m = manifest()
        chains = {c["id"]: c for c in m["chains"]}
        for c in chains.values():
            self.assertTrue(c["source_url"] and c["retrieved_utc"], c["id"])
            self.assertIsNone(c["low_cutoff_hz"])
            self.assertEqual(c["low_cutoff_status"], "undocumented")
            self.assertEqual(c["distance_status"] == "stated", c["distance_m"] is not None, c["id"])
        for page in m["pages"]:
            if "chain" not in page:
                continue
            era = chains[page["chain"]]["era"]
            if "2012" in page["page"] or (page["session_date"] or "") >= "2011":
                self.assertEqual(era, "2011+", page["page"])
            else:
                self.assertEqual(era, "pre-2011", page["page"])
            for f in m["files"]:
                if f["page"] == page["page"]:
                    self.assertEqual(f["chain"], page["chain"], f["id"])
        irpdata("test_chain_assignment_by_era", len(m["pages"]), 0, sorted(chains), start)

    def test_ambient_entry_or_absence_per_session(self):
        start = time.monotonic()
        m = manifest()
        ids = {f["id"]: f for f in m["files"]}
        present = 0
        sessions = [p for p in m["pages"] if "chain" in p]
        self.assertEqual(len(sessions), len(m["pages"]) - 1)
        for page in sessions:
            if page["ambient"] == "absent-on-page":
                self.assertFalse(any(f["kind"] == "ambient-silence" for f in m["files"] if f["page"] == page["page"]))
            else:
                self.assertEqual(ids[page["ambient"]]["kind"], "ambient-silence")
                self.assertEqual(ids[page["ambient"]]["page"], page["page"])
                present += 1
        self.assertGreaterEqual(present, 1, "non-vacuity: the Bb clarinet session has one")
        irpdata("test_ambient_entry_or_absence_per_session", len(sessions), 0, present, start)

    def test_ledger_rows_resolve_both_ways(self):
        start = time.monotonic()
        m = manifest()
        ids = {f["id"] for f in m["files"]}
        rows = {r["source_id"]: r for r in ledger_rows()}
        refs = [a.split(":", 1)[1] for r in rows.values() for a in r["artifacts"].split(";") if a.startswith("reference:")]
        self.assertTrue(refs, "non-vacuity: the ledger names reference ids")
        for ref in refs:
            self.assertIn(ref, ids, ref)
        for source in {f["ledger_source_id"] for f in m["files"]}:
            self.assertEqual(rows[source]["decision"], "admit", source)
            artifacts = rows[source]["artifacts"].split(";")
            for path in ("data/instruments/references.json", "data/instruments/fetch-references.sh"):
                self.assertIn("table:" + path, artifacts)
        irpdata("test_ledger_rows_resolve_both_ways", len(rows), 0, len(refs), start)

    def test_measured_content_defects_are_recorded(self):
        start = time.monotonic()
        defects = [f for f in manifest()["files"] if f.get("defects")]
        kinds = {d["kind"] for f in defects for d in f["defects"]}
        self.assertEqual(kinds, {"sample-rate-header", "note-absent"})
        for f in defects:
            if any(d["kind"] == "sample-rate-header" for d in f["defects"]):
                self.assertFalse(f["gradable"], f["id"])
        irpdata("test_measured_content_defects_are_recorded", len(defects), len(defects), sorted(f["id"] for f in defects), start)

    def test_fetch_refuses_data_volume_destination(self):
        start = time.monotonic()
        with tempfile.TemporaryDirectory() as directory:
            result, events = fetch(directory, [entry(b"x")], mount="/System/Volumes/Data")
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("Mac Data volume", result.stderr)
            self.assertEqual(events, [])
            self.assertEqual(tree(Path(directory) / "refs"), [])
        irpdata("test_fetch_refuses_data_volume_destination", 1, 1, result.returncode, start)

    def test_fetch_refuses_mismatched_existing_file_and_deletes_nothing(self):
        start = time.monotonic()
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "refs/sub/dir/a file.aiff"
            target.parent.mkdir(parents=True)
            target.write_bytes(b"wrong bytes")
            before = tree(Path(directory) / "refs")
            result, events = fetch(directory, [entry(b"right bytes")])
            self.assertEqual(result.returncode, 1)
            self.assertEqual([e["event"] for e in events], ["refused"])
            self.assertEqual(events[0]["sha256_observed"], hashlib.sha256(b"wrong bytes").hexdigest())
            self.assertEqual(tree(Path(directory) / "refs"), before)
        irpdata("test_fetch_refuses_mismatched_existing_file_and_deletes_nothing", 1, 1, events[0]["event"], start)

    def test_fetch_verify_mode_never_downloads_when_present(self):
        start = time.monotonic()
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "refs/sub/dir/a file.aiff"
            target.parent.mkdir(parents=True)
            target.write_bytes(b"right bytes")
            before = tree(Path(directory) / "refs")
            # The URL is unreachable: any download attempt would refuse.
            result, events = fetch(directory, [entry(b"right bytes")])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([e["event"] for e in events], ["verified"])
            self.assertEqual(tree(Path(directory) / "refs"), before)
        irpdata("test_fetch_verify_mode_never_downloads_when_present", 1, 0, events[0]["event"], start)

    def test_fetch_refuses_escaping_paths(self):
        start = time.monotonic()
        with tempfile.TemporaryDirectory() as directory:
            entries = [dict(entry(b"x"), id=f"escape-{i}", path=p) for i, p in enumerate(("../out.bin", "/abs.bin", "a//b.bin"))]
            result, events = fetch(directory, entries)
            self.assertEqual(result.returncode, 1)
            self.assertEqual([e["event"] for e in events], ["refused"] * 3)
            self.assertFalse((Path(directory) / "out.bin").exists())
            self.assertEqual(tree(Path(directory) / "refs"), [])
        irpdata("test_fetch_refuses_escaping_paths", 3, 3, [e["event"] for e in events], start)

    def test_fetch_checks_remote_bytes_before_retaining(self):
        start = time.monotonic()
        payload = b"served payload"
        with tempfile.TemporaryDirectory() as directory:
            served = Path(directory) / "served"
            served.mkdir()
            (served / "f.bin").write_bytes(payload)
            handler = functools.partial(Quiet, directory=str(served))
            server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
            threading.Thread(target=server.serve_forever, daemon=True).start()
            url = f"http://127.0.0.1:{server.server_port}/f.bin"
            try:
                good = dict(entry(payload, url), path="good.bin")
                bad = dict(entry(b"other bytes", url), id="bad-entry", path="bad.bin")
                result, events = fetch(directory, [good, bad])
                pinned, pins = fetch(directory, [bad], mode="--pin")
            finally:
                server.shutdown()
                server.server_close()
            refs = Path(directory) / "refs"
            self.assertEqual(result.returncode, 1)
            self.assertEqual([e["event"] for e in events], ["fetched", "refused"])
            self.assertEqual((refs / "good.bin").read_bytes(), payload)
            self.assertFalse((refs / "bad.bin").exists())
            self.assertEqual((refs / "bad.bin.part").read_bytes(), payload)
            self.assertEqual(pinned.returncode, 0)
            self.assertEqual(pins[0]["event"], "pinned")
            self.assertEqual(pins[0]["sha256_observed"], hashlib.sha256(payload).hexdigest())
            self.assertFalse((refs / "bad.bin").exists())
        irpdata("test_fetch_checks_remote_bytes_before_retaining", 2, 1, [e["event"] for e in events], start)


if __name__ == "__main__":
    unittest.main()
