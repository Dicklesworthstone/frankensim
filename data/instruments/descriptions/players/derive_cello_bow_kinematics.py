#!/usr/bin/env python3
"""Bow position, bow velocity and bow-bridge distance of a cello player, read off an open motion-capture file by script (not by hand).

    python3 derive_cello_bow_kinematics.py DataMocap_Player.csv cello-open-strings-bow-kinematics-pamies-vila-2024.tsv

Input: DataMocap_Player.csv of Zenodo 10.5281/zenodo.10696680 (Pamies-Vila, Mayer, Matusiak, Chatziioannou, "A method for the
reproduction of cello bow kinematics using a robotic arm and motion capture", Acta Acustica 2024; CC BY 4.0). One professional
cellist, about one minute of the cello part of Vivaldi's "Spring" on the open strings A3, D3, G2, OptiTrack, 240 frames/s. The file
(md5 7cc2513edfc033108700a3fab3f35dce on the Zenodo record, sha-256 2dcf299c...f29f1c as downloaded) is NOT committed; fetch it from
https://zenodo.org/api/records/10696680/files/DataMocap_Player.csv/content and this script refuses any file whose header differs.

Geometry (all from the file's own markers, metres, one frame at a time):
  bow axis    the line from MFrog to MTip (the two markers on the stick's ends; the hair is not marked);
  string      the line from MNut to MBridge (the nut and the top of the bridge);
  contact     the closest approach of the two lines. It gives
    bow_pos_m           distance from MFrog along the bow axis to that point,
    bow_vel_m_s         d(bow_pos_m)/dt by central differences (one-sided at the two ends), positive toward the tip; it is a velocity
                        relative to the string line, so the cello's own motion is already removed,
    bridge_over_string  (distance from that point to MBridge along the string line) / |MNut MBridge|, a relative bow-bridge distance
                        (the usual beta) up to the offset between the bridge marker and the string's end on the bridge,
    line_sep_m          the distance between the two lines at their closest approach,
    skew_deg            the angle between the two lines, 90 degrees when the bow is square to the string.
What this does NOT give: bow force; whether the hair is on the string. The markers are on the stick, so line_sep_m is a few centimetres
when the bow is on the string and larger when it is off, and the file has no contact signal. The nut-bridge marker distance is not the
open string length to a stated tolerance. Nothing here is a fit to a recording.
"""
from __future__ import annotations

import csv
import hashlib
import math
import sys
from pathlib import Path

EXPECTED_HEADER = [
    "Time_s", "MTip_x", "MTip_y", "MTip_z", "MFrog_x", "MFrog_y", "MFrog_z", "MBow_A_x", "MBow_A_y", "MBow_A_z", "MBow_B_x", "MBow_B_y", "MBow_B_z",
    "MBow_C_x", "MBow_C_y", "MBow_C_z", "MNut_x", "MNut_y", "MNut_z", "MBridge_x", "MBridge_y", "MBridge_z",
]
SCHEMA = "frankensim.player-cello-bow-kinematics.v1"
COLUMNS = ["frame", "bow_pos_m", "bow_vel_m_s", "bridge_over_string", "line_sep_m", "skew_deg"]
Vec = tuple[float, float, float]


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a: Vec) -> float:
    return math.sqrt(dot(a, a))


def unit(a: Vec) -> Vec:
    n = norm(a)
    return (a[0] / n, a[1] / n, a[2] / n)


def closest_approach(frog: Vec, tip: Vec, nut: Vec, bridge: Vec) -> tuple[float, float, float, float]:
    """(bow_pos_m, bridge_over_string, line_sep_m, skew_deg) for one frame."""
    u = unit(sub(tip, frog))
    string_len = norm(sub(bridge, nut))
    w = unit(sub(bridge, nut))
    r0 = sub(frog, nut)
    b = dot(u, w)
    d, e = dot(u, r0), dot(w, r0)
    den = 1.0 - b * b  # |u| = |w| = 1
    if den < 1e-12:
        raise ValueError("the bow axis is parallel to the string line in a frame; there is no closest approach")
    s = (b * e - d) / den  # along the bow from the frog
    r = (e - b * d) / den  # along the string from the nut
    p_bow = (frog[0] + s * u[0], frog[1] + s * u[1], frog[2] + s * u[2])
    p_str = (nut[0] + r * w[0], nut[1] + r * w[1], nut[2] + r * w[2])
    return s, (string_len - r) / string_len, norm(sub(p_bow, p_str)), math.degrees(math.acos(min(1.0, abs(b))))


def gradient(y: list[float], t: list[float]) -> list[float]:
    """Central differences on the actual time stamps, one-sided at the ends."""
    n = len(y)
    out = [0.0] * n
    out[0] = (y[1] - y[0]) / (t[1] - t[0])
    out[-1] = (y[-1] - y[-2]) / (t[-1] - t[-2])
    for i in range(1, n - 1):
        out[i] = (y[i + 1] - y[i - 1]) / (t[i + 1] - t[i - 1])
    return out


def derive(rows: list[list[float]]) -> list[list[float]]:
    t = [r[0] for r in rows]
    col = {n: i for i, n in enumerate(EXPECTED_HEADER)}
    pt = lambda r, m: (r[col[f"{m}_x"]], r[col[f"{m}_y"]], r[col[f"{m}_z"]])  # noqa: E731
    approach = [closest_approach(pt(r, "MFrog"), pt(r, "MTip"), pt(r, "MNut"), pt(r, "MBridge")) for r in rows]
    s = [a[0] for a in approach]
    v = gradient(s, t)
    return [[i, s[i], v[i], approach[i][1], approach[i][2], approach[i][3]] for i in range(len(rows))]


def read(path: Path) -> list[list[float]]:
    with path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        if header != EXPECTED_HEADER:
            raise SystemExit(f"{path}: the header is not the one this script was written for: {header}")
        rows = [[float(x) for x in line] for line in reader]
    if len(rows) < 3:
        raise SystemExit(f"{path}: too few frames")
    return rows


def write(derived: list[list[float]], path: Path) -> None:
    with path.open("w", newline="") as f:
        f.write(SCHEMA + "\n" + "\t".join(COLUMNS) + "\n")
        for fr, s, v, beta, sep, skew in derived:
            f.write(f"{int(fr)}\t{s:.5f}\t{v:.4f}\t{beta:.5f}\t{sep:.5f}\t{skew:.2f}\n")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    rows = read(src)
    write(derive(rows), dst)
    print(f"{len(rows)} frames from {src.name} (sha-256 {hashlib.sha256(src.read_bytes()).hexdigest()}) -> {dst}")


if __name__ == "__main__":
    main()
