#!/usr/bin/env python3
"""Make the public copy of the Instrument Explorer served at https://frankensim.org/instrument-tests/.

Usage:
  python3 apps/instrument-explorer/tools/publish.py [--data DIR] [--jobs N]
  npx wrangler deploy --config apps/instrument-explorer/cloudflare/wrangler.jsonc

Reads the app (this directory) and the data built by tools/build.py (default ~/.local/share/instrument-explorer/site) and
writes dist-cf/instrument-tests/ (ignored by git, like the audio itself):
  - every FLAC the data points at becomes AAC in an .m4a (Apple's encoder, 160 kbit/s); the spectrograms stay as built,
    from the lossless audio;
  - recordings whose licence does not allow redistribution (YouTube videos without a reuse licence) lose their audio and
    stay listed with their source link; licence texts say what each licence allows;
  - build-host paths are dropped, and the build refuses to finish if any text names the private product. Its names stay
    out of this public repository: one regular expression per line in ~/.config/instrument-explorer/private-names, and
    publishing refuses to start without that file.
Only files the data points at are copied, so a withheld recording never reaches the output. Existing outputs are reused
while they carry their source's mtime; nothing is deleted. Media paths keep build.py's ?v=<source version> query.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
DIST = APP / "dist-cf"
SITE = "instrument-tests"
DEFAULT_DATA = Path.home() / ".local" / "share" / "instrument-explorer" / "site"
URL = f"https://frankensim.org/{SITE}/"

UIOWA = "University of Iowa Musical Instrument Samples: free to use without restriction"
WITHHELD = ("This video is under YouTube's standard licence (all rights reserved), so it is not played here. "
            "It was used only to grade; follow the source link to hear it.")
HOST_PATH = re.compile(r"(?:/data/tmp|/home/ubuntu|/Users/[^/\s]+)/[^\s'\"`),;]*")
PRIVATE_NAMES = Path.home() / ".config" / "instrument-explorer" / "private-names"
DROP_KEYS = {"real_audio"}


def private_pattern() -> re.Pattern:
    """The private product's names, from PRIVATE_NAMES (kept out of the repository): one regular expression per line."""
    if not PRIVATE_NAMES.is_file():
        sys.exit(f"missing {PRIVATE_NAMES}: one regular expression per line naming the private product")
    lines = [s.strip() for s in PRIVATE_NAMES.read_text().splitlines() if s.strip() and not s.lstrip().startswith("#")]
    if not lines:
        sys.exit(f"{PRIVATE_NAMES} names nothing")
    return re.compile("|".join(f"(?:{s})" for s in lines), re.I)


PRIVATE = private_pattern()


def withheld(ref: dict) -> bool:
    return "YouTube standard" in (ref.get("licence") or "")


def fix_reference(r: dict) -> None:
    lic = r.get("licence") or ""
    if lic.startswith("University of Iowa"):
        r["licence"], r["licence_notes"] = UIOWA, None
    elif "youtube" in (r.get("licence_id") or "") or "YouTube" in lic:
        r["licence_notes"] = None
    if withheld(r):
        r["audio"], r["withheld"] = None, WITHHELD


def fix_pair(p: dict) -> None:
    real = p.get("real")
    if isinstance(real, dict) and (real.get("licence") or "").startswith("UIowa"):
        real["licence"] = UIOWA


def scrub(x, media: set[str], where: str):
    """Rewrite one JSON value for the public build: .flac -> .m4a, host paths out; record every media path kept."""
    if isinstance(x, dict):
        return {k: scrub(v, media, where) for k, v in x.items() if k not in DROP_KEYS}
    if isinstance(x, list):
        return [scrub(v, media, where) for v in x]
    if isinstance(x, str):
        if PRIVATE.search(x):
            sys.exit(f"refusing to publish: {where} names the private product: {x[:160]!r}")
        if x.startswith(("audio/", "img/")) and " " not in x:
            path, _, version = x.partition("?")  # build.py's ?v=<source version>, kept so a re-render is fetched again
            media.add(path)
            if path.endswith(".flac"):
                path = path[:-5] + ".m4a"
            return f"{path}?{version}" if version else path
        return HOST_PATH.sub("…", x)
    return x


def encode(src: Path, dst: Path) -> str | None:
    # The .m4a takes its FLAC's mtime and is remade whenever they differ (a replacement FLAC may be the older file)
    st = src.stat()
    if dst.exists() and dst.stat().st_mtime_ns == st.st_mtime_ns:
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.stem + ".part.m4a")
    base = ["ffmpeg", "-nostdin", "-hide_banner", "-loglevel", "error", "-y", "-i", str(src), "-map_metadata", "-1",
            "-c:a", "aac_at", "-b:a", "160k", "-movflags", "+faststart"]
    for extra in ([], ["-ar", "48000"]):  # Apple's AAC tops out at 48 kHz
        r = subprocess.run(base + extra + [str(tmp)], capture_output=True, text=True)
        if r.returncode == 0:
            os.utime(tmp, ns=(st.st_atime_ns, st.st_mtime_ns))
            os.replace(tmp, dst)
            return None
    return f"{src}: {r.stderr.strip()[-300:]}"


def copy_if_newer(src: Path, dst: Path) -> None:
    # copy2 keeps the source's mtime, so any difference means a different file (an older replacement included)
    if dst.exists() and dst.stat().st_mtime_ns == src.stat().st_mtime_ns and dst.stat().st_size == src.stat().st_size:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


HEADERS = f"""/{SITE}/data/audio/*
  Cache-Control: public, max-age=86400
/{SITE}/data/img/*
  Cache-Control: public, max-age=86400
/{SITE}/share/*
  Cache-Control: public, max-age=86400
/{SITE}/data/*.json
  Cache-Control: public, max-age=300
/{SITE}/js/*
  Cache-Control: public, max-age=300
/{SITE}/css/*
  Cache-Control: public, max-age=300
/{SITE}/
  Cache-Control: public, max-age=60
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    data, out = args.data, DIST / SITE
    out.mkdir(parents=True, exist_ok=True)

    # The app: unchanged files, plus the public URL in the page head.
    app_files = ["css/app.css", *sorted(str(p.relative_to(APP)) for p in (APP / "js").rglob("*.js"))]
    for rel in app_files:
        copy_if_newer(APP / rel, out / rel)
    page = (APP / "index.html").read_text()
    marker = '<link rel="preconnect"'
    if marker not in page:
        sys.exit(f"index.html has no {marker!r} to put the canonical URL before")
    # Share cards (share/og.png 1200x630, share/twitter.png 1200x600): absolute URLs, each with a content hash so a
    # new image is fetched again by the platforms that cache by URL.
    share = []
    for name, (w, hgt) in {"og.png": (1200, 630), "twitter.png": (1200, 600)}.items():
        src = APP / "share" / name
        if not src.is_file():
            sys.exit(f"missing {src}: render it from tools/share/og.html (see the README)")
        copy_if_newer(src, out / "share" / name)
        share.append((f"share/{name}", f"{URL}share/{name}?v={hashlib.sha256(src.read_bytes()).hexdigest()[:10]}", w, hgt))
    alt = ("A spectrogram of a simulated trumpet fanfare above one of the same notes played on a real trumpet, with "
           "their pitch and loudness traced beneath: Instruments computed from physics, played next to the real thing.")
    (_, og_url, ow, oh), (_, tw_url, _, _) = share
    tags = [f'<link rel="canonical" href="{URL}">', f'<meta property="og:url" content="{URL}">',
            '<meta property="og:type" content="website">', '<meta property="og:site_name" content="FrankenSim">',
            f'<meta property="og:image" content="{og_url}">', '<meta property="og:image:type" content="image/png">',
            f'<meta property="og:image:width" content="{ow}">', f'<meta property="og:image:height" content="{oh}">',
            f'<meta property="og:image:alt" content="{alt}">', '<meta name="twitter:card" content="summary_large_image">',
            '<meta name="twitter:title" content="Instrument tests · FrankenSim">',
            '<meta name="twitter:description" content="Physics-simulated instruments next to real recordings, with every note measured: pitch, timbre, noise, attack.">',
            f'<meta name="twitter:image" content="{tw_url}">', f'<meta name="twitter:image:alt" content="{alt}">']
    page = page.replace(marker, "\n  ".join(tags) + f"\n  {marker}", 1)
    (out / "index.html").write_text(page)

    # The data.
    media: set[str] = set()
    counts = {"withheld": 0, "references": 0}
    docs = sorted(p.relative_to(data) for p in [*data.glob("*.json"), *(data / "family").glob("*.json"),
                                                *(data / "phrase").glob("*.json")])
    for rel in docs:
        doc = json.loads((data / rel).read_text())
        for p in doc.get("phrases", []) if isinstance(doc, dict) else []:
            ref = p.get("ref")
            if ref and any("standard" in (s.get("licence") or "").lower() for s in ref.get("sources", [])):
                sys.exit(f"refusing to publish: {rel} phrase {p['id']} splices an all-rights-reserved recording")
            for s in (ref or {}).get("sources", []):
                if s.get("licence") == "UIowa":
                    s["licence"] = UIOWA
        if rel.parts[0] == "family":
            for r in doc.get("references", []):
                fix_reference(r)
                counts["references"] += 1
                counts["withheld"] += r["audio"] is None
            for p in doc.get("pairs", []):
                fix_pair(p)
        doc = scrub(doc, media, str(rel))
        dst = out / "data" / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(json.dumps(doc, ensure_ascii=False, separators=(",", ":")))

    missing = sorted(m for m in media if not (data / m).is_file())
    audio = sorted(m for m in media if m.endswith(".flac") and m not in missing)
    other = sorted(m for m in media if not m.endswith(".flac") and m not in missing)
    for m in other:
        copy_if_newer(data / m, out / "data" / m)
    failures = []
    with cf.ThreadPoolExecutor(args.jobs) as pool:
        for i, err in enumerate(pool.map(lambda m: encode(data / m, out / "data" / (m[:-5] + ".m4a")), audio), 1):
            if err:
                failures.append(err)
            if i % 200 == 0:
                print(f"  encoded {i}/{len(audio)}", flush=True)
    (DIST / "_headers").write_text(HEADERS)

    # Files an earlier build left that this one does not point at (a recording since withheld, a renamed phrase, a
    # half-written encode) stay on disk but are kept out of the upload by .assetsignore.
    wanted = {"index.html", *app_files, *(s[0] for s in share), *(f"data/{rel}" for rel in docs), *(f"data/{m}" for m in other),
              *(f"data/{m[:-5]}.m4a" for m in audio)}
    stale = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file() and str(p.relative_to(out)) not in wanted)
    escape = lambda s: re.sub(r"([\\*?\[\]!#])", r"\\\1", s)  # noqa: E731
    (DIST / ".assetsignore").write_text("".join(f"/{SITE}/{escape(s)}\n" for s in stale))

    size = sum((out / w).stat().st_size for w in wanted if (out / w).is_file())
    print(f"{out}: {len(wanted)} files to publish, {size / 1e6:.1f} MB; {len(audio)} clips as AAC, {len(other)} images, "
          f"{counts['withheld']} of {counts['references']} recordings withheld (licence), {len(docs)} JSON documents")
    if stale:
        print(f"{len(stale)} files from earlier builds are left out of the upload (dist-cf/.assetsignore); "
              f"first: {stale[:3]}")
    if missing or failures:
        if missing:
            print(f"{len(missing)} paths in the data have no file (first: {missing[:3]})")
        if failures:
            print(f"{len(failures)} encodes failed:\n  " + "\n  ".join(failures[:10]))
        sys.exit(1)


if __name__ == "__main__":
    main()
