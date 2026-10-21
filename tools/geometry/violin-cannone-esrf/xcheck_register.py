"""xcheck_register.py ZIP RIBS_GLB SLICE : registration of the frontal-slice grid of the reduced 8-bit volume to the mesh frame (x across, z along the body), in translation and the four axis flips.
The slice K (at 1/8 scale through the JPEG draft mode, 8 x 0.03566 mm per pixel) is thresholded (grey > 40), each row is filled between its first and last wood pixel (a violin outline is not
concave along a row), and the filled footprint is cross-correlated by FFT with the filled xz projection of the ribs mesh (vertices rasterised at the same pitch, rows filled the same way).
Prints, per flip, the best overlap (intersection over union) and the shift, and the mesh coordinates of the soundpost-like small round feature if requested. Read-only."""
import io
import os
import sys
import zipfile

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cannone_io as cio  # noqa: E402

Image.MAX_IMAGE_PIXELS = None
VOX = 0.03566
zp, glb, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
S = 8 * VOX

z = zipfile.ZipFile(zp)
prefix = [n for n in z.namelist() if n.endswith("_0000.jpg")][0][:-8]
im = Image.open(io.BytesIO(z.read("%s%04d.jpg" % (prefix, k))))
im.draft("L", (im.width // 8, im.height // 8))
a = np.asarray(im.convert("L"))
print("slice", k, "shape", a.shape, "pixel pitch mm", S)
wood = a > 40


def fill_rows(m):
    out = np.zeros_like(m)
    for r in range(m.shape[0]):
        c = np.nonzero(m[r])[0]
        if len(c):
            out[r, c[0]:c[-1] + 1] = True
    return out


# rows of the slice that belong to the body ring (the neck and bridge cuts lie above it): take the largest row run whose fill exceeds 100 mm
widths = np.array([ (np.nonzero(r)[0][-1] - np.nonzero(r)[0][0] + 1) * S if r.any() else 0 for r in wood])
body_rows = np.nonzero(widths > 60)[0]
# the ring is the lowest connected block of such rows
runs = np.split(body_rows, np.nonzero(np.diff(body_rows) > 3)[0] + 1)
blk = max(runs, key=lambda r: len(r) * 1.0 if widths[r].max() > 150 else 0)
r0, r1 = blk[0], blk[-1]
print("body rows %d..%d (%.1f mm tall), max width %.1f mm" % (r0, r1, (r1 - r0 + 1) * S, widths[blk].max()))
img_mask = np.zeros_like(wood)
img_mask[r0:r1 + 1] = fill_rows(wood[r0:r1 + 1])
# columns extent of the ring
cols = np.nonzero(img_mask.any(axis=0))[0]
c0, c1 = cols[0], cols[-1]
print("body cols %d..%d (%.1f mm)" % (c0, c1, (c1 - c0 + 1) * S))

v, t, info = cio.glb_triangles(glb)
xs, zs = v[:, 0], v[:, 2]
print("ribs mesh x %.2f..%.2f z %.2f..%.2f (extent %.2f x %.2f)" % (xs.min(), xs.max(), zs.min(), zs.max(), xs.max() - xs.min(), zs.max() - zs.min()))
H, W = a.shape
best = None
for fx in (1, -1):
    for fz in (1, -1):
        # mesh grid at pitch S, index i = (x - xmin)/S etc., orientation chosen so that image col grows with fx*x and image row with fz*z... fz=-1: row grows as z falls (neck, +z, at the top)
        gx = fx * xs
        gz = fz * zs
        ci = np.floor((gx - gx.min()) / S).astype(int)
        ri = np.floor((gz - gz.min()) / S).astype(int)
        m = np.zeros((ri.max() + 3, ci.max() + 3), dtype=bool)
        m[ri + 1, ci + 1] = True
        m = fill_rows(m)
        # place into the big canvas and cross-correlate via FFT
        Hc, Wc = H + m.shape[0], W + m.shape[1]
        nh, nw = 1 << (Hc - 1).bit_length(), 1 << (Wc - 1).bit_length()
        fa = np.fft.rfft2(img_mask.astype(np.float32), s=(nh, nw))
        fm = np.fft.rfft2(m.astype(np.float32), s=(nh, nw))
        cc = np.fft.irfft2(fa * np.conj(fm), s=(nh, nw))
        inter = cc.max()
        iy, ix = np.unravel_index(cc.argmax(), cc.shape)
        # shift (dy, dx) such that image index = mesh-grid index + shift
        dy = iy if iy < nh // 2 else iy - nh
        dx = ix if ix < nw // 2 else ix - nw
        union = img_mask.sum() + m.sum() - inter
        iou = inter / union
        print("fx %+d fz %+d: IoU %.4f, shift rows %d cols %d, mesh area %.0f mm2 image area %.0f mm2" % (fx, fz, iou, dy, dx, m.sum() * S * S, img_mask.sum() * S * S))
        if best is None or iou > best[0]:
            best = (iou, fx, fz, dy, dx, gx.min(), gz.min())
iou, fx, fz, dy, dx, gxmin, gzmin = best
print("BEST fx=%+d fz=%+d IoU=%.4f" % (fx, fz, iou))
# image pixel (row, col) centre -> mesh x = fx * (gxmin + (col - dx - 1 + 0.5) * S), z = fz * (gzmin + (row - dy - 1 + 0.5) * S)
# since mesh cell i covers gx in [gxmin + i S, gxmin + (i+1) S) and was placed at canvas index i+1, image index = i + 1 + shift
print("MAP x_mm = %+d * (%.5f + (col - %d + 0.5) * %.6f)" % (fx, gxmin, dx + 1, S))
print("MAP z_mm = %+d * (%.5f + (row - %d + 0.5) * %.6f)" % (fz, gzmin, dy + 1, S))
# self-consistency: translated mesh mask vs image mask extents
