"""xcheck_plates.py VOL.npy JOB [JOB ...] with JOB = MAPS.tsv:top|back:THR:Y0:OUT.tsv (Y0 may be "none") : cross-check of the plate maps (mesh) against the reduced 8-bit volume (vol_decode.py).
For every 3 mm block of the maps table with at least 30 of 36 valid cells and at least 8 mm from the outline, the pixel columns of the volume inside the block (about 10 x 10 at 0.285 mm) are read
along the slice index k (the vertical axis, 0.03566 mm per slice): the wood run is the run of grey > THR (runs separated by 6 slices or fewer are joined) whose outer end is nearest the expected outer surface
k = (outer y - Y0) / 0.03566 within 0.7 mm, or, with Y0 none, the first (top plate) or last (back plate) run.
Per block: median run start and end and length, the plane fit of the outer end across the block (slope in x and z), and the maps-table values, one row per block in OUT.tsv (scratch, not committed).
The volume is read once and the jobs run in forked workers that share it (six readers of the same compressed file at once are very slow)."""
import sys

import numpy as np

S = 8 * 0.03566
VOX = 0.03566
MAP_X0, MAP_DX = -101.49264, 13.5
MAP_Z0, MAP_DZ = 88.69451, 824.5
GAP = 6
# the window around the expected outer surface when Y0 is given: the offset Y0 is fitted per plate from a first pass without it
WIN = int(0.7 / 0.03566)


def scan(vol, mapsp, which, outp, THR, Y0):
    top = which == "top"
    rows = [l.rstrip("\n").split("\t") for l in open(mapsp) if not l.startswith("frankensim.")]
    hdr, data = rows[0], [[float(v) for v in r] for r in rows[1:]]
    ix = {h: i for i, h in enumerate(hdr)}
    blocks = [r for r in data if r[ix["valid_cells_of_36"]] >= 30 and r[ix["dist_to_outline_mm"]] >= 8.0]
    print("%s plate: %d map rows, %d used (valid >= 30 of 36 and >= 8 mm from the outline)" % (which, len(data), len(blocks)))


    def run_of(prof, k_exp=None):
        nz = np.nonzero(prof > THR)[0]
        if len(nz) == 0:
            return None
        if k_exp is not None:
            # all gap-joined runs; take the one whose outer end is nearest the expected outer surface, if within the window
            cut = np.nonzero(np.diff(nz) > GAP)[0]
            starts = np.r_[nz[0], nz[cut + 1]]
            ends = np.r_[nz[cut], nz[-1]]
            outer = starts if top else ends
            j = int(np.argmin(np.abs(outer - k_exp)))
            return (starts[j], ends[j]) if abs(outer[j] - k_exp) <= WIN else None
        if top:
            k0 = nz[0]
            brk = np.nonzero(np.diff(nz) > GAP)[0]
            k1 = nz[brk[0]] if len(brk) else nz[-1]
            return k0, k1
        k1 = nz[-1]
        brk = np.nonzero(np.diff(nz) > GAP)[0]
        k0 = nz[brk[-1] + 1] if len(brk) else nz[0]
        return k0, k1


    out = open(outp, "w")
    out.write("x_mm\tz_mm\tmap_outer_y\tmap_inner_y\tmap_thick_normal\tn_cols\tk_start_med\tk_end_med\tlen_med_slices\tlen_p10\tlen_p90\tgrad_x\tgrad_z\tfrac_ok\n")
    by_row = {}
    for b in blocks:
        by_row.setdefault(b[ix["z_mm"]], []).append(b)
    for zc in sorted(by_row):
        r_lo = int(np.ceil((MAP_Z0 - (zc + 1.5)) / S + MAP_DZ - 0.5))
        r_hi = int(np.floor((MAP_Z0 - (zc - 1.5)) / S + MAP_DZ - 0.5))
        slab = np.ascontiguousarray(vol[:, r_lo:r_hi + 1, :].transpose(1, 2, 0))  # (rows, cols, K)
        for b in by_row[zc]:
            xc = b[ix["x_mm"]]
            c_lo = int(np.ceil((xc - 1.5 - MAP_X0) / S + MAP_DX - 0.5))
            c_hi = int(np.floor((xc + 1.5 - MAP_X0) / S + MAP_DX - 0.5))
            pts = []
            for ri in range(r_hi - r_lo + 1):
                for ci in range(c_lo, c_hi + 1):
                    rr = run_of(slab[ri, ci], None if Y0 is None else (b[ix["outer_y_mm"]] - Y0) / 0.03566)
                    if rr is not None:
                        pts.append((ri + r_lo, ci, rr[0], rr[1]))
            if len(pts) < 20:
                continue
            p = np.array(pts, dtype=float)
            ln = p[:, 3] - p[:, 2] + 1
            ref = p[:, 2] if top else p[:, 3]
            # plane fit of the outer surface (top: run start; back: run end) over the block
            A = np.c_[np.ones(len(p)), (p[:, 1] - p[:, 1].mean()) * S, -(p[:, 0] - p[:, 0].mean()) * S]
            coef = np.linalg.lstsq(A, ref * VOX, rcond=None)[0]
            okf = np.mean(np.abs(ln - np.median(ln)) <= 0.5 / VOX)
            out.write("%.1f\t%.1f\t%.3f\t%.3f\t%.3f\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.4f\t%.4f\t%.3f\n" % (
                xc, zc, b[ix["outer_y_mm"]], b[ix["inner_y_mm"]], b[ix["thickness_normal_mm"]], len(p), np.median(p[:, 2]), np.median(p[:, 3]), np.median(ln),
                np.percentile(ln, 10), np.percentile(ln, 90), coef[1], coef[2], okf))
    out.close()
    print("wrote", outp)


VOL = None


def runjob(j):
    mapsp, which, thr, y0, outp = j
    scan(VOL, mapsp, which, outp, int(thr), None if y0 == "none" else float(y0))
    return outp


if __name__ == "__main__":
    # xcheck_plates2.py VOL.npy JOB [JOB ...] with JOB = MAPS.tsv:top|back:THR:Y0:OUT.tsv (Y0 may be "none"); the volume is read once and the jobs run in forked workers that share it
    from multiprocessing import get_context
    VOL = np.load(sys.argv[1])
    jobs = [j.split(":") for j in sys.argv[2:]]
    with get_context("fork").Pool(min(len(jobs), 6)) as p:
        for o in p.imap_unordered(runjob, jobs):
            print("finished", o, flush=True)
