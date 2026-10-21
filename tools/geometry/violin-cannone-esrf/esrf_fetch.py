"""esrf_fetch.py LIST OUTDIR : fetch the listed files of the public ESRF il Cannone data collection (DOI 10.15151/ESRF-DC-2413577471, CC BY 4.0) through the portal's own
unauthenticated route /ids/data/download?datafileIds=<id>. LIST lines are `<datafile id> <size> <path in the dataset>` (paths may contain spaces). Each file is written to
OUTDIR/<id>_<basename with spaces as underscores>; its size must equal the listed size or it is reported and kept as .partial. Refuses a file of ESRF_FETCH_MAX bytes (default 1.5 GB) or more. Existing files are skipped."""
import os
import sys
import urllib.request

UA = "FrankenSim-geometry-fetch/1.0 (open-source acoustics project; github.com/Dicklesworthstone/frankensim; one-off CC BY data collection fetch)"
MAX = int(os.environ.get("ESRF_FETCH_MAX", "1500000000"))  # bytes; larger files are refused unless the caller raises it for one named file


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **k):
        return None


op = urllib.request.build_opener(NoRedirect)
lst, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for line in open(lst, encoding="utf-8").read().splitlines():
    if not line.strip():
        continue
    did, size, path = line.split(" ", 2)
    size = int(size)
    name = "%s_%s" % (did, os.path.basename(path).replace(" ", "_"))
    dst = os.path.join(out, name)
    if size >= MAX:
        print("REFUSED (>= 1.5 GB)", name, size, flush=True)
        continue
    if os.path.exists(dst) and os.path.getsize(dst) == size:
        print("have", name, flush=True)
        continue
    try:
        op.open(urllib.request.Request("https://icatplus.esrf.fr/ids/data/download?datafileIds=%s" % did, headers={"User-Agent": UA}), timeout=60)
        loc = None
    except urllib.error.HTTPError as e:
        loc = e.headers.get("Location")
    if not loc:
        print("NO REDIRECT", name, flush=True)
        continue
    got = 0
    tmp = dst + ".partial"
    with urllib.request.urlopen(urllib.request.Request(loc, headers={"User-Agent": UA}), timeout=300) as r, open(tmp, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b)
            got += len(b)
    if got == size:
        os.replace(tmp, dst)
        print("ok", name, got, flush=True)
    else:
        print("SIZE MISMATCH", name, "got", got, "listed", size, "(kept as .partial)", flush=True)
