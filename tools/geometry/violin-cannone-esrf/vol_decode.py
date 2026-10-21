"""vol_decode.py ZIP K0 K1 OUT.npy NPROC : decode frontal slices K0..K1-1 of the reduced 8-bit volume through the JPEG draft mode at 1/8 scale (one pixel = 8 x 8 voxels = 0.28528 mm in the slice
plane, slice spacing 0.03566 mm) into one uint8 array (K1-K0, rows, cols). Read-only on the zip; the output is a scratch file on the analysis host, never committed."""
import io
import os
import sys
import time
import zipfile
from multiprocessing import Pool

import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
zp, k0, k1, out, nproc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
SHAPE = (2117, 744)


def work(span):
    a, b = span
    z = zipfile.ZipFile(zp)
    prefix = [n for n in z.namelist() if n.endswith("_0000.jpg")][0][:-8]
    mm = np.load(out, mmap_mode="r+")
    for k in range(a, b):
        im = Image.open(io.BytesIO(z.read("%s%04d.jpg" % (prefix, k))))
        im.draft("L", (im.width // 8, im.height // 8))
        arr = np.asarray(im.convert("L"))
        assert arr.shape == SHAPE, arr.shape
        mm[k - k0] = arr
    mm.flush()
    return b - a


if __name__ == "__main__":
    t0 = time.time()
    mm = np.lib.format.open_memmap(out, mode="w+", dtype=np.uint8, shape=(k1 - k0,) + SHAPE)
    del mm
    step = 8
    spans = [(a, min(a + step, k1)) for a in range(k0, k1, step)]
    done = 0
    with Pool(nproc) as p:
        for n in p.imap_unordered(work, spans):
            done += n
            if done % 160 == 0:
                print("decoded %d / %d slices, %.0f s" % (done, k1 - k0, time.time() - t0), flush=True)
    print("done %d slices in %.0f s -> %s" % (done, time.time() - t0, out), flush=True)
