#!/usr/bin/env python3
"""Serve the Instrument Explorer locally: the app from this directory's parent, the built data under /data/.

Usage:
  python3 apps/instrument-explorer/tools/serve.py [--port 8742] [--data DIR] [--open]

Binds 127.0.0.1 only; the recordings it serves are for local listening, not for publishing.
"""
from __future__ import annotations

import argparse
import mimetypes
import os
import re
import sys
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit

APP = Path(__file__).resolve().parents[1]
DEFAULT_DATA = Path.home() / ".local" / "share" / "instrument-explorer" / "site"

for ext, kind in ((".flac", "audio/flac"), (".webp", "image/webp"), (".js", "text/javascript"), (".json", "application/json"),
                  (".wav", "audio/wav"), (".svg", "image/svg+xml")):
    mimetypes.add_type(kind, ext)


class Handler(SimpleHTTPRequestHandler):
    data_dir: Path = DEFAULT_DATA

    def translate_path(self, path: str) -> str:
        rel = unquote(urlsplit(path).path)
        root, rel = (self.data_dir, rel[len("/data/"):]) if rel.startswith("/data/") else (APP, rel.lstrip("/"))
        target = (root / rel).resolve()
        if root.resolve() not in (target, *target.parents):
            return str(root / "__outside__")  # refuses ../ escapes with a 404
        return str(target)

    def end_headers(self):
        path = urlsplit(self.path).path
        if path.endswith((".json", ".js", ".css", ".html")) or path == "/":
            self.send_header("Cache-Control", "no-cache")
        else:
            self.send_header("Cache-Control", "max-age=3600")
        self.send_header("Accept-Ranges", "bytes")
        super().end_headers()

    def send_head(self):
        """Plain files with single byte ranges (audio seeking in any player); everything else as the stdlib does it."""
        rng = self.headers.get("Range")
        path = Path(self.translate_path(self.path))
        if not rng or not path.is_file():
            return super().send_head()
        m = re.fullmatch(r"bytes=(\d*)-(\d*)", rng.strip())
        size = path.stat().st_size
        if not m or (not m.group(1) and not m.group(2)):
            return super().send_head()
        if m.group(1):
            start = int(m.group(1))
            end = min(int(m.group(2)), size - 1) if m.group(2) else size - 1
        else:
            start, end = max(size - int(m.group(2)), 0), size - 1
        if start > end or start >= size:
            self.send_response(416)
            self.send_header("Content-Range", f"bytes */{size}")
            self.end_headers()
            return None
        f = open(path, "rb")
        f.seek(start)
        self.send_response(206)
        self.send_header("Content-Type", self.guess_type(str(path)))
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(end - start + 1))
        self.end_headers()
        self._range_left = end - start + 1
        return f

    def copyfile(self, source, outputfile):
        left = getattr(self, "_range_left", None)
        if left is None:
            return super().copyfile(source, outputfile)
        while left > 0:
            chunk = source.read(min(65536, left))
            if not chunk:
                break
            outputfile.write(chunk)
            left -= len(chunk)
        self._range_left = None

    def log_message(self, fmt, *args):
        if os.environ.get("IE_QUIET"):
            return
        super().log_message(fmt, *args)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--port", type=int, default=8742)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA, help=f"built site data (default {DEFAULT_DATA})")
    ap.add_argument("--open", action="store_true", help="open the explorer in the default browser")
    args = ap.parse_args()
    if not (args.data / "index.json").is_file():
        print(f"no built data in {args.data}: run tools/sync.sh and tools/build.py first", file=sys.stderr)
        return 1
    Handler.data_dir = args.data
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    url = f"http://localhost:{args.port}/"
    print(f"Instrument Explorer on {url} (app {APP}, data {args.data}); Ctrl-C stops it.", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
