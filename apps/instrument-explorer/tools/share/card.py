#!/usr/bin/env python3
"""Render the share cards (share/og.png 1200x630, share/twitter.png 1200x600) from one phrase's real analysis.

Usage:
  python3 apps/instrument-explorer/tools/share/card.py [--phrase trumpet__fanfare] [--title "Trumpet · fanfare"]

The card is og.html: the phrase's two spectrograms, its measured pitch and loudness lanes, and the site's counts read
from the last build. The page is served from a staging directory and photographed by headless Chrome.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE.parents[1]
HOME_DATA = Path.home() / ".local" / "share" / "instrument-explorer"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phrases", type=Path, help="a phrasekit output; default the newest phrases-<UTC>")
    ap.add_argument("--site", type=Path, default=HOME_DATA / "site", help="the built data, for the counts")
    ap.add_argument("--phrase", default="trumpet__fanfare")
    ap.add_argument("--title", default="Trumpet · fanfare")
    ap.add_argument("--playhead", type=float, default=0.62, help="where the playhead crosses, as a fraction of the phrase")
    args = ap.parse_args()
    ph = args.phrases or max(HOME_DATA.glob("phrases-*"), default=None)
    if ph is None:
        sys.exit("no phrases-<UTC> directory; pass --phrases")
    p = next((x for x in json.loads((ph / "phrases.json").read_text())["phrases"] if x["id"] == args.phrase), None)
    if p is None or not (p.get("ref") or {}).get("map"):
        sys.exit(f"{ph}: no phrase {args.phrase} with a real track")
    an = json.loads((ph / f"{args.phrase}.analysis.json").read_text())
    idx = json.loads((args.site / "index.json").read_text())
    n_phrases = len(json.loads((args.site / "phrases.json").read_text())["phrases"])

    stage = HOME_DATA / "share-card"
    stage.mkdir(parents=True, exist_ok=True)
    shutil.copy(HERE / "og.html", stage / "og.html")
    shutil.copy(ph / p["map"], stage / "sim.webp")
    shutil.copy(ph / p["ref"]["map"], stage / "ref.webp")
    fr = an["frames"]
    (stage / "data.json").write_text(json.dumps(dict(
        title=args.title, sub="simulation above, real notes below", sim="sim.webp", ref="ref.webp",
        playhead=args.playhead, hop=an["hop_s"], duration=p["duration_s"],
        stats=dict(i=f"{idx['totals']['families']:,}", r=f"{idx['totals']['references']:,}", p=f"{n_phrases:,}"),
        notes=[dict(at=n["at"], dur=n["dur"], target_hz=n.get("target_hz")) for n in p.get("notes") or []],
        frames={k: dict(f0=fr[k]["f0"], level=fr[k]["level"]) for k in ("sim", "ref")})))

    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(Quiet, directory=str(stage)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    port = server.server_address[1]
    (APP / "share").mkdir(exist_ok=True)
    try:
        for name, h in (("og.png", 630), ("twitter.png", 600)):
            out = APP / "share" / name
            subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars", "--force-device-scale-factor=1",
                            f"--window-size=1200,{h}", f"--screenshot={out}", "--virtual-time-budget=15000",
                            f"http://127.0.0.1:{port}/og.html?h={h}"], check=True, capture_output=True)
            print("wrote", out.relative_to(APP))
    finally:
        server.shutdown()


if __name__ == "__main__":
    main()
