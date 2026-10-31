"""Differential diagnosis: which link of the simulated signal chain makes FrankenSim sound
unlike the real instrument, and which experiment would decide it.

verisim.py grades how far apart two sounds are. This script asks WHY, by splitting the
difference along the chain hammer -> string -> bridge -> board -> radiation -> microphone,
using observables only a simulation has (the bridge velocity) and ablation renders that
switch one mechanism on at a time.

usage: uv run --with numpy --with pillow python tools/diagnose.py CAMPAIGN_DIR REF_WAV48_DIR PUBLIC_DIR [KIT_DIR SALAMANDER_OH_DIR]

Inputs beyond verisim.py's:
  piano/a3-mf-bridge.csv     bridge velocity and microphone pressure of the A3 mf render
  piano/a3-mf.log            the render log (its bare-board FE mode list)
  piano/a3-mf-close.wav      optional: A3 mf heard by two receivers about 0.23-0.25 m above the board,
                             where the Iowa microphones sat (when the graded renders listen elsewhere)
  piano/a3-mf-far.wav        optional: when the graded renders listen at the Iowa microphone positions,
                             the same A3 mf at the former default point 1 m above the board
  piano/a3-{mf,ff,fff}-legacyfelt.wav  optional: when the graded renders use RT-0425's hammer law
                             (--rt0425-hammer-stiffness), the same strokes with the former felt law
                             (the graded command without its hammer flags)
  piano/a3-mf-p{x,y}.wav     optional: the same A3 mf heard by stereo pairs 0.1 m either side of the bass-side
                             Iowa position, along x and along y (how much a near-field grade depends on that estimate)
  piano/a3-mf-estmap.wav     optional: the same A3 mf with the former estimated bridge stations (same board,
                             air and receivers; its log names the board file it used)
  piano/board-modes.csv      optional: --dump-board of the graded board (bare-board modes and per-key bridge
                             displacements), for the driving-point mobility at A3's contact
  piano/scale-et.csv         optional: --dump-scale of the rendered scale (string partials near a resonance,
                             A3's tension, length and unison detune)
  piano/a3-{mf,ff}-waves.wav optional: A3 with longitudinal string waves (--string-stretching),
                             compared with linear graded renders
  piano/a3-fff{,-waves}.wav  optional: the same pair struck at 8 m/s
  piano/a3-{mf,ff,fff}-linear.wav  optional: when the graded A3 renders carry string waves themselves
                             (their logs' cmd: rows name --string-stretching), the same strokes on linear strings
  piano/a3-mf-polar.wav      optional: A3 mf with both transverse polarizations (--string-polarization)
  perc/vk-splash.wav, perc/vk-splash-grid-*mm.wav  optional: the kit solver's splash (with its log) and
                             the same stroke on other grids
  piano/silent-control.wav   optional: C3 struck alone (where this model's C3 actually sounds)
  piano/a3-basis.json        optional: --pressure-basis-json of the A3 board setup (loaded-board
                             reference frequencies; a 0.05 s render writes it in seconds)
Writes PUBLIC_DIR/data/diagnosis.json.
"""
from __future__ import annotations
import json, math, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import verisim as V  # noqa: E402  (shared, identical measurement code)

SR = V.SR
F0 = B0 = None  # the rendered A3's own f0 and B (f_n = n f0 sqrt(1 + B n^2)), fitted in main()


def r(v, d=4):
    return V.r(v, d)


def fk(k, f0=None, b=None):
    f0, b = (F0 if f0 is None else f0), (B0 if b is None else b)
    return k * f0 * math.sqrt(1 + b * k * k)


def piano(name, dur=3.0, channel=0):
    return V.segment(V.highpass(V.load(V.CAMPAIGN / "piano" / f"{name}.wav", channel)), dur)


def real(name, dur=3.0):
    return V.segment(V.highpass(V.load(V.REFS / f"{name}.wav")), dur)


def render_cmd(name):
    """The render's own command line (the cmd: row of its log): cards describe what was simulated."""
    log = V.CAMPAIGN / "piano" / f"{name}.log"
    if not log.exists():
        return ""
    return next((line for line in log.read_text(errors="replace").splitlines() if line.startswith("cmd:")), "")


def render_settings(name):
    """The render command without its board file and outputs: two renders whose settings match differ only in those."""
    words, kept, skip = render_cmd(name).split(), [], False
    for w in words:
        if skip:
            skip = False
        elif w in ("--board-geometry", "--render", "--bridge-trace-csv"):
            skip = True
        else:
            kept.append(w)
    return kept


def has(name):
    return (V.CAMPAIGN / "piano" / f"{name}.wav").exists()


def heard(name):
    """Where a render listened, from its own log: a phrase for the cards and its receiver points."""
    log = V.CAMPAIGN / "piano" / f"{name}.log"
    text = log.read_text(errors="replace") if log.exists() else ""
    points = [p for p in V.receivers(text) if p is not None]
    if V.at_iowa(text):
        return "where the Iowa microphones sat", points
    return "at " + " and ".join(f"({x:.2f}, {y:.2f}, {z:.2f}) m" for x, y, z in points), points


def stretching_row(name, key):
    """The waves row for one key in the --string-stretching file the render used, or None."""
    cmd = render_cmd(name).split()
    if "--string-stretching" not in cmd:
        return None
    path = V.CAMPAIGN / "piano" / cmd[cmd.index("--string-stretching") + 1]
    for line in path.read_text().splitlines():
        f = line.split(",")
        if f[0] == "waves" and int(f[1]) == key:
            return {"ea_n": float(f[2]), "band_hz": float(f[4]), "mode_cap": int(f[5]), "damping_ratio": float(f[6])}
    return None


def level_at(f, s, target, rel=0.012):
    sel = (f > target * (1 - rel)) & (f < target * (1 + rel))
    return s[sel].max()


def partial_levels(x, t0, t1, count=16, f0=None, b=None):
    f, s = V.spectrum(x, t0, t1, 1 << 19)
    p1 = level_at(f, s, fk(1, f0, b))
    return [r(20 * math.log10(level_at(f, s, fk(k, f0, b)) / p1), 3) for k in range(1, count + 1)]


def own_fit(x, f_guess=220.0):
    """Each signal's partials sit where ITS string puts them: f0 and B from verisim's fit."""
    q = V.partials(x, f_guess)
    return q["f0"], q["B"]


def non_partial_db(x, t0, t1, kmax=40):
    """Energy NOT within +-(0.6 % + 3 Hz) of one of this signal's own string partials, re total 70 Hz-16 kHz."""
    f0, b = own_fit(x)
    f, s = V.spectrum(x, t0, t1, 1 << 18)
    p = s ** 2
    mask = np.ones_like(f, bool)
    for k in range(1, kmax + 1):
        mask &= ~((f > fk(k, f0, b) * 0.994 - 3) & (f < fk(k, f0, b) * 1.006 + 3))
    band = (f > 70) & (f < 16000)
    return 10 * math.log10(p[mask & band].sum() / p[band].sum())


def non_partial_peaks(x, t0, t1, count=8):
    f, s = V.spectrum(x, t0, t1, 1 << 20)
    p1 = level_at(f, s, fk(1))
    mask = np.ones_like(f, bool)
    for k in range(1, 41):
        mask &= ~((f > fk(k) * 0.994 - 3) & (f < fk(k) * 1.006 + 3))
    idx = np.where(mask & (f > 60) & (f < 6000))[0]
    ss = s[idx]
    cand = sorted(((ss[i], f[idx[i]]) for i in range(1, len(idx) - 1) if ss[i] > ss[i - 1] and ss[i] >= ss[i + 1]), reverse=True)
    out = []
    for a, fq in cand:
        if all(abs(fq - g["hz"]) > 4 for g in out):
            out.append({"hz": r(fq, 5), "db": r(20 * math.log10(a / p1), 3)})
        if len(out) >= count:
            break
    return out


def precise_peak(x, lo, hi, t0, t1):
    f, s = V.spectrum(x, t0, t1, 1 << 21)
    sel = (f > lo) & (f < hi)
    k = int(np.argmax(np.where(sel, s, 0)))
    a, b, c = np.log(s[k - 1:k + 2] + 1e-30)
    return float(f[k] + 0.5 * (a - c) / (a - 2 * b + c) * (f[1] - f[0]))


def ratio_db(x, lo, hi, t0, t1, f1=None):
    """Strongest component in [lo, hi] relative to partial 1 (at f1, this signal's own; default the rendered A3)."""
    f, s = V.spectrum(x, t0, t1, 1 << 21)
    sel = (f > lo) & (f < hi)
    return 20 * math.log10(s[sel].max() / level_at(f, s, fk(1) if f1 is None else f1))


def two_stage(x, kmax=4):
    """Early (0.1-0.8 s) and late (1.5-2.9 s) decay slope per partial: a real piano's
    'prompt sound' then 'aftersound' shows as a steep early slope and a flat late one.
    A flat late slope could also be a partial that has sunk into the noise floor, so each
    row also reports the late level above the floor 40-90 Hz either side of the partial."""
    f0, b0 = own_fit(x)
    n, hop = 4096, 2400
    fr = np.array([x[i:i + n] * np.hanning(n) for i in range(int(0.1 * SR), len(x) - n, hop)])
    S = np.abs(np.fft.rfft(fr, axis=1))
    fb = np.fft.rfftfreq(n, 1 / SR)
    t = 0.1 + (np.arange(len(fr)) * hop + n / 2) / SR
    out = []
    for k in range(1, kmax + 1):
        fp = fk(k, f0, b0)
        b = int(np.argmin(abs(fb - fp)))
        e = 20 * np.log10(S[:, b - 2:b + 3].max(axis=1) + 1e-12)
        ring = ((fb > fp + 40) & (fb < fp + 90)) | ((fb > fp - 90) & (fb < fp - 40))
        floor = 20 * np.log10(np.median(S[:, ring], axis=1) + 1e-12)
        late = t > 1.5
        out.append({"n": k, "early_db_s": r(np.polyfit(t[t < 0.8], e[t < 0.8], 1)[0], 3),
                    "late_db_s": r(np.polyfit(t[late], e[late], 1)[0], 3),
                    "late_above_floor_db": r(float(np.mean(e[late] - floor[late])), 3)})
    return out


def phantoms(x, f_guess=220.0, t0=0.02, t1=0.5, ks=range(8, 17)):
    """Phantom partials: a stretched string's tension modulation makes sum tones f_m + f_n.
    Inharmonicity puts them BELOW partial m+n, so they are resolvable. Positions use THIS
    signal's own fitted f0 and B. Each row reports the peak within 3 Hz of the predicted sum
    frequency, relative to the true partial k and to the local floor (the median 8-40 Hz away,
    excluding the true partial). A phantom counts as present only when it is at least 10 dB above
    that floor AND within 40 dB of the true partial (a simulation's floor is so low that numerical
    dust clears the first test alone). Positions closer than 15 Hz to the true partial are skipped:
    they sit inside its spectral skirt."""
    fit = V.partials(x, f_guess)
    f0, b = fit["f0"], fit["B"]
    f, s = V.spectrum(x, t0, t1, 1 << 20)
    rows = []
    for k in ks:
        m = k // 2
        fp, ft = fk(m, f0, b) + fk(k - m, f0, b), fk(k, f0, b)
        if ft - fp < 15:  # inside the true partial's skirt: not resolvable
            continue
        w = (f > fp - 3) & (f < fp + 3)
        a = s[w].max()
        ring = ((f > fp - 40) & (f < fp - 8)) | ((f > fp + 8) & (f < min(fp + 40, ft - 8)))
        above = 20 * math.log10(a / np.median(s[ring]))
        re_partial = 20 * math.log10(a / level_at(f, s, ft, 0.004))
        rows.append({"k": k, "hz": r(fp, 5), "re_partial_db": r(re_partial, 3),
                     "above_floor_db": r(above, 3), "present": bool(above >= 10 and re_partial >= -40)})
    return rows


def bare_board_modes():
    log = (V.CAMPAIGN / "piano" / "a3-mf.log").read_text(errors="replace")
    return [float(m.group(1)) for m in re.finditer(r"source FE mode \d+: \[([0-9.]+),", log)]


def loaded_board_references():
    """Loaded-board diagonal reference frequencies (board plus the strings' static attachment)."""
    p = V.CAMPAIGN / "piano" / "a3-basis.json"
    return sorted(json.loads(p.read_text())["loaded_reference_frequency_hz"]) if p.exists() else None


NOTE = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]


def string_partials_near(target, width=1.5):
    """Partials of every speaking string in the rendered scale (--dump-scale CSV), with each
    string's own inharmonicity, within `width` Hz of target."""
    p = V.CAMPAIGN / "piano" / "scale-et.csv"
    if not p.exists():
        return None
    import csv
    seen, out = set(), []
    for row in csv.DictReader(l for l in p.read_text().splitlines() if not l.startswith("#")):
        key, ln, mu, t, ei = int(row["midi"]), float(row["length_m"]), float(row["linear_density_kg_m"]), float(row["tension_n"]), float(row["flexural_rigidity_nm2"])
        b = math.pi ** 2 * ei / (t * ln * ln)
        for n in range(1, 12):
            f = n * math.sqrt(t / mu) / (2 * ln) * math.sqrt(1 + b * n * n)
            if abs(f - target) < width and (key, n) not in seen:
                seen.add((key, n))
                out.append({"key": key, "name": f"{NOTE[key % 12]}{key // 12 - 1}", "n": n, "hz": r(f, 6)})
    return sorted(out, key=lambda q: q["hz"])


def bridge_split():
    d = np.loadtxt(V.CAMPAIGN / "piano" / "a3-mf-bridge.csv", delimiter=",", skiprows=1)
    v, p = d[:, 2], d[:, 4]
    early = {"bridge_db": r(ratio_db(v, 128, 135, 0.05, 0.5), 3), "pressure_db": r(ratio_db(p, 128, 135, 0.05, 0.5), 3)}
    lv, lp = partial_levels(v, 0.0, 0.15), partial_levels(p, 0.0, 0.15)
    return early, lv, lp


def envelope_gap(x, ref, f_guess=220.0):
    """verisim's graded partial-envelope distance: partials 2-16 found in both, 0.05-1.05 s after onset."""
    er = {p["n"]: p["db"] for p in V.partials(ref, f_guess)["partials"]}
    eo = {p["n"]: p["db"] for p in V.partials(x, f_guess)["partials"]}
    common = sorted(set(er) & set(eo) - {1})
    return V.rms_db([er[n] for n in common], [eo[n] for n in common])


def board_admittance(key=57):
    """Bare-board driving-point mobility Y(f) [m/(N s)] at a key's bridge contact, from a --dump-board modal CSV
    of the graded board: sum over modes of phi^2 i w / (w_m^2 - w^2 + 2 i zeta_m w_m w), phi the mass-normalised
    bridge displacement. Returns (Y, mode count, top mode Hz), or None without the CSV."""
    p = V.CAMPAIGN / "piano" / "board-modes.csv"
    if not p.exists():
        return None
    modes = []
    for line in p.read_text().splitlines():
        f = line.split(",")
        if not line.startswith(("#", "mode")) and int(f[4]) == key:
            modes.append((2 * math.pi * float(f[1]), float(f[2]), float(f[5])))

    def y(hz):
        w = 2 * math.pi * hz
        return sum(phi * phi * 1j * w / (wm * wm - w * w + 2j * z * wm * w) for wm, z, phi in modes)
    return y, len(modes), max(m[0] for m in modes) / (2 * math.pi)


def contact_point(name, key=57):
    """Board-frame position [m] of a key's bridge contact in the board file the render used (--board-geometry)."""
    rows_ = board_rows(name)
    if rows_ is None:
        return None
    nodes, tris, site = {}, {}, None
    for f in rows_:
        if f[0] == "node":
            nodes[int(f[1])] = (float(f[2]), float(f[3]))
        elif f[0] == "triangle":
            tris[int(f[1])] = [int(v) for v in f[2:5]]
        elif f[0] == "bridge" and int(f[1]) == key:
            site = (int(f[2]), [float(v) for v in f[3:6]])
    return tuple(sum(w * nodes[n][i] for w, n in zip(site[1], tris[site[0]])) for i in (0, 1))


def board_rows(name):
    """The node, triangle, clamped-node and bridge rows of the board file a render used (--board-geometry)."""
    cmd = render_cmd(name).split()
    if "--board-geometry" not in cmd:
        return None
    path = V.CAMPAIGN / "piano" / cmd[cmd.index("--board-geometry") + 1]
    return [line.split(",") for line in path.read_text(errors="replace").splitlines()]


def front_edge_y(name, x, halfwidth=0.03):
    """The board's front (keyboard-side) edge at x [m]: the lowest clamped boundary node within halfwidth of x."""
    rows_ = board_rows(name) or []
    nodes = {int(f[1]): (float(f[2]), float(f[3])) for f in rows_ if f[0] == "node"}
    ys = [nodes[int(f[1])][1] for f in rows_ if f[0] == "fixed" and abs(nodes[int(f[1])][0] - x) < halfwidth]
    return min(ys) if ys else None


def scale_row(key):
    """The rendered scale's row for a key (--dump-scale CSV), or None."""
    p = V.CAMPAIGN / "piano" / "scale-et.csv"
    if not p.exists():
        return None
    import csv
    return next((row for row in csv.DictReader(l for l in p.read_text().splitlines() if not l.startswith("#"))
                 if int(row["midi"]) == key), None)


NEW_AIR = ('Air reaction (--rayleigh-loading) is not yet admitted on the production board. A resolved board mesh fine enough for 2.4 kHz (9,432'
           ' area samples, 6.3 per wavelength) clears the mesh limit once the offline work ceilings are raised, but the 32-port passive fit of the full '
           '128-mode radiation impedance then misses its accuracy limits: complex error 19% peak and 15% RMS, error in its resistive (power) part 13% peak and 14% RMS, '
           'against 5%, 2%, 10% and 5% (two sides, 90 Hz-2.4 kHz).')


def main():
    out = {"inputs": {}, "rows": []}
    board = bare_board_modes()
    rows = out["rows"]
    global F0, B0
    real_mf = real("Piano.mf.A3")
    ours_mf = piano("a3-mf")
    F0, B0 = own_fit(ours_mf)  # every default position below is the rendered A3's own
    f0_r, b_r = own_fit(real_mf)
    air = "--rayleigh-loading" in render_cmd("a3-mf")  # the graded renders carry air loading
    # Some measurements below were made on the preset's 8-division board without air loading.
    resolved_board = "--board-geometry" in render_cmd("a3-mf")
    old_board = "Measured earlier on the preset's 8-division board without air loading, not on the board graded above:"
    a3_text = (V.CAMPAIGN / "piano" / "a3-mf.log").read_text(errors="replace")
    iowa = V.at_iowa(a3_text)  # graded at the Iowa microphone positions (left channel: bass side)
    # Strings meet the bridge at RT-0425's published points. The literal measurements below that were made at the
    # former estimated A3 station describe another contact point; they are left out rather than relabelled.
    published = "RT-0425 Appendix A coupling points" in a3_text
    adm = board_admittance()

    # 1. The early low resonance
    knock_hz = {n: r(precise_peak(piano(n, 1.6), 125, 137, 0.05, 1.5), 6) for n in ["a3-pp", "a3-mf", "a3-ff", "a3-fff", "c6-mf"]
                if (V.CAMPAIGN / "piano" / f"{n}.wav").exists()}
    early, lv, lp = bridge_split()
    real_lv = partial_levels(real_mf, 0.0, 0.15, f0=f0_r, b=b_r)
    ours_lv = partial_levels(ours_mf, 0.0, 0.15)
    ours_knock = ratio_db(ours_mf, 125, 137, 0.05, 0.5)
    real_knock = ratio_db(real_mf, 80, 0.77 * 220, 0.05, 0.5, f1=fk(1, f0_r, b_r))
    spread = max(knock_hz.values()) - min(knock_hz.values())
    hz0 = knock_hz["a3-mf"]
    below, above = max(b for b in board if b < hz0), min(b for b in board if b > hz0)
    loaded = loaded_board_references()
    strs = string_partials_near(hz0)
    strings_txt = ""
    if strs:
        strings_txt = (" The nearest string partials are " + ", ".join(f"{q['name']}{'' if q['n'] == 1 else ' ×' + str(q['n'])}" for q in strs)
                       + f" at {strs[0]['hz']:.2f}-{strs[-1]['hz']:.2f} Hz, {hz0 - strs[-1]['hz']:.1f}-{hz0 - strs[0]['hz']:.1f} Hz below it; their dampers are down while A3 plays.")
    rows.append({
        "id": "low-resonance",
        "symptom": f"An early resonance near {np.mean(list(knock_hz.values())):.1f} Hz sits at {early['pressure_db']:.1f} dB re the fundamental; the real piano has nothing stronger than {real_knock:.1f} dB anywhere from 80 to {0.77 * 220:.0f} Hz.",
        "evidence": (f"Same frequency within {spread:.2f} Hz for A3 at every stroke speed and for C6, so it belongs to the modelled instrument, not to the struck string. "
                     f"It is {abs(knock_hz['a3-mf'] - 130.813):.2f} Hz from nominal C3; that difference neither identifies nor excludes an unplayed string or a coupled board/string mode. "
                     f"It is not a bare-board mode (the nearest are {below:.1f} and {above:.1f} Hz)"
                     + (f", nor a mode of the board loaded by the strings' static stiffness (those start at {loaded[0]:.1f} Hz)" if loaded else "")
                     + f". No bare-board neighbour within {min(abs(b - knock_hz['a3-mf']) for b in board):.0f} Hz identifies the loaded instrument's modes by itself."
                     + strings_txt
                     + f" At A3's bridge it is already {early['bridge_db']:.1f} dB re partial 1; radiation to the microphone adds {early['pressure_db'] - early['bridge_db']:.1f} dB."),
        "suspect": "A coupled mechanical resonance driven too strongly or damped too lightly. The measurements do not identify a single cause: either the dampered bass strings ring with the board, or the board and strings resonate together. Radiation adds only the last few dB.",
        "test": "Render A3 with the sustain pedal down: if the resonance follows the lifted dampers, it is the bass strings." + ("" if air else " " + NEW_AIR),
        "status": "open",
        "data": {"hz": knock_hz, "bridge_db": early["bridge_db"], "pressure_db": early["pressure_db"], "real_db": r(real_knock, 3),
                 "bare_board_neighbours_hz": [below, above], "loaded_board_lowest_hz": loaded[:4] if loaded else None, "string_partials": strs},
    })

    # 2. Ringing that is not the strings
    np_real = {d: r(non_partial_db(real(f"Piano.{d}.A3"), 0.2, 1.0), 3) for d in ["pp", "mf", "ff"]}
    np_ours = {d: r(non_partial_db(piano(f"a3-{d}"), 0.2, 1.0), 3) for d in ["pp", "mf", "ff"]}
    peaks = non_partial_peaks(ours_mf, 0.2, 1.0)
    for pk in peaks:
        k = max(1, int(round(pk["hz"] / F0)))
        pk["nearest_partial"] = k
        pk["offset_hz"] = r(pk["hz"] - fk(k), 4)
        pk["nearest_board_mode_hz"] = r(min(board, key=lambda b: abs(b - pk["hz"])), 5)
    sat = sorted({p["nearest_partial"] for p in peaks[:6] if abs(p["offset_hz"]) < 15})
    lows = [p for p in peaks[:6] if p["hz"] < 0.8 * F0]
    rows.append({
        "id": "board-ringing",
        "symptom": (f"Between 0.2 and 1 s, {np_ours['mf'] - np_real['mf']:.0f} dB more of FrankenSim's energy lies away from the string partials than the real piano's (mf: {np_ours['mf']:.1f} vs {np_real['mf']:.1f} dB re total)."
                    if np_ours["mf"] - np_real["mf"] > 3 else
                    f"Between 0.2 and 1 s, the share of FrankenSim's energy away from the string partials is within 3 dB of the real piano's (mf: {np_ours['mf']:.1f} vs {np_real['mf']:.1f} dB re total)."),
        "evidence": "Strongest off-partial peaks: " + ", ".join(f"{p['hz']:.0f} Hz ({p['db']:.0f} dB)" for p in peaks[:6]) + "." + (f" Those within 15 Hz of partial{'s' if len(sat) > 1 else ''} {', '.join(map(str, sat))} may be satellites from string-board coupling; proximity alone does not establish their modal participation." if sat else "") + (f" {len(lows)} lie below the struck note's fundamental; the board and unplayed strings can both participate there." if lows else ""),
        "suspect": "Coupled mechanical ringing: wood, unplayed-string/damper losses, air reaction and rim/case losses remain candidates until mode participation and controlled ablations distinguish them.",
        "test": "Raise the board's modal damping, or add air reaction, and re-measure the off-partial energy and the satellites' decay.",
        "status": "open",
        "data": {"real_db": np_real, "ours_db": np_ours, "peaks": peaks},
    })

    # 3. Excitation pattern: string side vs radiation side (bridge trace and its own pressure column)
    coloration = [r(a - b, 3) for a, b in zip(lp, lv)]
    rows.append({
        "id": "partial-pattern",
        "symptom": "In the first 150 ms FrankenSim's partials 3-6 are " + ", ".join(f"{ours_lv[k]:.0f}" for k in range(2, 6)) + " dB re partial 1; the real piano's are " + ", ".join(f"{real_lv[k]:.0f}" for k in range(2, 6)) + " dB.",
        "evidence": "At the bridge, partials 4 and 6 are already " + f"{lv[3]:.0f} and {lv[5]:.0f} dB (string and hammer side). From bridge to microphone the levels move by " + ", ".join(f"{c:+.0f}" for c in coloration[:8]) + f" dB for partials 1-8: a single receiver point colours each partial by up to ±{max(abs(c) for c in coloration[:8]):.0f} dB.",
        "suspect": ("Two parts: the hammer/strike pattern on the string, and the soundboard's mobility at A3's bridge contact, which sets how strongly each partial drives the board. "
                    "The receiver point near the board colours each partial further." if iowa else
                    "Two parts: the hammer/strike pattern on the string, and the single receiver point. The Iowa microphones sat 8 inches above the strings; ours is 1 m above the board."),
        "test": ("Compare the board's driving-point mobility at A3's contact with one measured on a real grand; ours is a reconstructed Model D, the recording a Model B." if iowa else
                 "Render receivers near the Iowa positions; compare power averaged over several receivers rather than one point."),
        "status": "open",
        "data": {"real": real_lv, "ours": ours_lv, "ours_bridge_trace": {"bridge": lv, "pressure": lp}, "radiation_coloration_db": coloration},
    })

    # 4. Decay shape
    ts_r, ts_o = two_stage(real_mf), two_stage(ours_mf)
    def stage2(ts):  # steep early fall, then a hold within 3 dB/s
        return [q["n"] for q in ts if q["early_db_s"] < -15 and abs(q["late_db_s"]) < 3]
    s2r, s2o = stage2(ts_r), stage2(ts_o)
    names = lambda ns: (f"partial {ns[0]}" if len(ns) == 1 else "partials " + ", ".join(map(str, ns[:-1])) + f" and {ns[-1]}") if ns else "no partial"
    held = [q["late_above_floor_db"] for q in ts_r if q["n"] in s2r]
    rows.append({
        "id": "decay-shape",
        "symptom": f"Two-stage decay (a steep fall, then a hold: 'prompt sound' and 'aftersound') shows on {names(s2r)} of the real piano and on {names(s2o)} of FrankenSim.",
        "evidence": "Early (0.1-0.8 s) -> late (1.5-2.9 s) slope, dB/s: real " + "; ".join(f"p{q['n']} {q['early_db_s']:.0f} -> {q['late_db_s']:.0f}" for q in ts_r) + "; FrankenSim " + "; ".join(f"p{q['n']} {q['early_db_s']:.0f} -> {q['late_db_s']:.0f}" for q in ts_o) + "."
                    + (f" The real piano's held partials stay {min(held):.0f}-{max(held):.0f} dB above the recording's noise floor, so the hold is the string, not hiss." if held else ""),
        "suspect": "Two-stage decay comes from the second transverse direction of each string and the coupling of the three unison strings through the bridge. FrankenSim's strings vibrate in one plane by default.",
        "test": "Render with --string-polarization (needs an authored bridge-frame file for all 88 keys).",
        "status": "open",
        "data": {"real": ts_r, "ours": ts_o},
    })

    # 5. ff brightness and phantom partials
    real_ff = real("Piano.ff.A3")
    ours_ff = piano("a3-ff")
    ph_real, ph_ours = phantoms(real_ff), phantoms(ours_ff)
    gain_real = V.band_ratio(real_ff) - V.band_ratio(real_mf)
    gain_ours = V.band_ratio(ours_ff) - V.band_ratio(ours_mf)
    # The 8 m/s stroke is the one whose loudness matches the real ff (piano_dynamics); its brightness decides the tail of the symptom.
    dull8 = V.band_ratio(real_ff) - V.band_ratio(piano("a3-fff")) if has("a3-fff") else None
    peak_db = lambda x: 20 * math.log10(float(np.abs(x).max()))
    # How far the 8 m/s stroke's peak level above mf falls short of the real ff's (positive: quieter than the real ff).
    loud8 = ((peak_db(real_ff) - peak_db(real_mf)) - (peak_db(piano("a3-fff")) - peak_db(ours_mf))) if has("a3-fff") else None
    at8 = "" if loud8 is None else f" (8 m/s, {abs(loud8):.1f} dB {'quieter' if loud8 > 0 else 'louder'} than the real ff relative to mf)"
    matched = ("." if dull8 is None else
               f", and stays {abs(dull8):.1f} dB {'duller' if dull8 > 0 else 'brighter'} even when struck about as hard as the real ff{at8}."
               if abs(dull8) > 3 else
               f". Struck about as hard as the real ff{at8}, it lands within {abs(dull8):.1f} dB of the real ff's brightness.")
    row = {
        "id": "ff-brightness",
        "symptom": f"From mf to ff the real piano gains {gain_real:.1f} dB of 2-8 kHz brightness; FrankenSim at 4.5 m/s gains {gain_ours:.1f} dB" + matched,
        "evidence": f"Sum tones f_m + f_n of a tension-modulated string (phantom partials) stand at least 10 dB above the local floor and within 40 dB of the true partial at {sum(q['present'] for q in ph_real)} of {len(ph_real)} resolvable positions (partials {ph_real[0]['k']}-{ph_real[-1]['k']}) in the real ff, as close as {max(q['re_partial_db'] for q in ph_real):.0f} dB to the true partial; in FrankenSim's ff at {sum(q['present'] for q in ph_ours)} of {len(ph_ours)}, the highest at {max(q['re_partial_db'] for q in ph_ours):.0f} dB. Each piano's positions come from its own inharmonicity fit.",
        "suspect": "No geometric string nonlinearity: a hard strike stretches the string, modulating its tension and driving longitudinal waves. For the A3 wire the first longitudinal mode is about 3.3 kHz, inside the missing band.",
        "test": "Render A3 mf and ff with --string-stretching waves (EA from the scale's own EI with E = 200 GPa) and compare brightness gain and phantom levels.",
        "status": "open",
        "data": {"phantoms_real_ff": ph_real, "phantoms_ours_ff": ph_ours},
    }
    waves = stretching_row("a3-ff", 57)
    adopted = waves is not None  # the graded A3 renders already carry string waves
    if waves and has("a3-ff-linear") and has("a3-mf-linear"):
        # The A3 renders above simulate the stretched string; the same strokes with linear strings are kept for comparison.
        lin_ff, lin_mf = piano("a3-ff-linear"), piano("a3-mf-linear")
        ph_lin = phantoms(lin_ff)
        gain_lin = V.band_ratio(lin_ff) - V.band_ratio(lin_mf)
        count = lambda rows_: sum(q["present"] for q in rows_)
        hard = ""
        if has("a3-fff") and has("a3-fff-linear"):
            new8, lin8 = piano("a3-fff"), piano("a3-fff-linear")
            ph8 = phantoms(new8)
            row["data"]["phantoms_ours_8ms"] = ph8
            hard = (f" Struck at 8 m/s, close to the real ff's loudness, phantom partials stand at {count(ph8)} of {len(ph8)} positions and 2-8 kHz brightness reaches "
                    f"{V.band_ratio(new8):.1f} dB, against {V.band_ratio(lin8):.1f} dB with linear strings and {V.band_ratio(real_ff):.1f} dB for the real ff.")
        row["data"]["phantoms_linear_ff"] = ph_lin
        row["data"]["string_waves"] = waves
        row["status"] = "improved"
        row["suspect"] = ("Two parts. A hard strike stretches the string: its tension rises with the motion, sum tones (phantom partials) appear, and longitudinal waves "
                          "push on the bridge. That is now simulated (see Fix) and it brings the phantom partials. It does not supply the remaining broadband dullness, "
                          "which points upstream to the hammer: how much stiffer the felt becomes, and how much shorter its contact, under a hard blow.")
        row["test"] = ("Compare hammer contact time and 2-8 kHz brightness from mf to ff across felt laws (--rt0425-hammer-stiffness, --hammers) "
                       "against the recording's brightness gain.")
        row["fix"] = (f"The A3 renders above include longitudinal string waves (axial rigidity EA = {waves['ea_n'] / 1000:.1f} kN from the scale's own bending stiffness "
                      f"with steel E = 200 GPa, axial modes to {waves['band_hz'] / 1000:.0f} kHz) whose axial force drives the bridge through a per-mode participation row read from the board file "
                      "(each mode's translation plus its rotation acting on the lever arm from the board's neutral surface to the string's bridge-top contact; "
                      "the string's direction is taken perpendicular to the bridge in plan, an assumption stated in piano/unknowns.tsv, and the flat board gives no downbearing), "
                      "with the work in the energy ledger. Against the same strokes with linear strings: "
                      f"phantom partials at ff stand at {count(ph_ours)} of {len(ph_ours)} positions (linear {count(ph_lin)}; real Steinway {count(ph_real)}), and the mf-to-ff "
                      f"brightness gain is {gain_ours:.1f} dB (linear {gain_lin:.1f} dB; real piano {gain_real:.1f} dB).") + hard
    hammer_flags = ("--rt0425-hammer-stiffness", "--rt0425-hammer-dissipation")
    if hammer_flags[0] in render_cmd("a3-ff") and has("a3-mf-legacyfelt") and has("a3-ff-legacyfelt"):
        # The graded renders use RT-0425's published hammer law; the same strokes with the former felt law are kept.
        # They isolate the hammer law only when the commands differ in nothing else.
        felt_only = all([w for w in render_settings(graded) if w not in hammer_flags] == render_settings(former)
                        for graded, former in (("a3-mf", "a3-mf-legacyfelt"), ("a3-ff", "a3-ff-legacyfelt")))
        b = {k: {"mf": V.band_ratio(x_mf), "ff": V.band_ratio(x_ff)} for k, x_mf, x_ff in
             (("real", real_mf, real_ff), ("published", ours_mf, ours_ff), ("former", piano("a3-mf-legacyfelt"), piano("a3-ff-legacyfelt")))}
        row["data"]["hammer_brightness_db"] = {k: {d: r(v, 3) for d, v in x.items()} for k, x in b.items()}
        law = ("stiffness K_H on each unison string (the former conversion divided it across the unison, leaving each string too soft)"
               + (" and dissipation R_H" if hammer_flags[1] in render_cmd("a3-ff") else ""))
        row["fix"] = (row.get("fix", "") + (" " if row.get("fix") else "")
                      + f"The hammers use RT-0425's published {law}. "
                      + ("Against the same strokes with the former felt law" if felt_only else
                         "Against strokes with the former felt law, in renders that also differ in other settings")
                      + f", 2-8 kHz brightness moves from {b['former']['mf']:.1f} to {b['published']['mf']:.1f} dB at mf and from "
                      f"{b['former']['ff']:.1f} to {b['published']['ff']:.1f} dB at ff (real piano {b['real']['mf']:.1f} and {b['real']['ff']:.1f} dB).")
        if all(abs(b["published"][d] - b["real"][d]) < abs(b["former"][d] - b["real"][d]) for d in ("mf", "ff")):
            row["status"] = "improved"
        # How far the 4.5 m/s stroke's peak level above mf falls short of the real ff's (positive: quieter than the real ff).
        loud45 = (peak_db(real_ff) - peak_db(real_mf)) - (peak_db(ours_ff) - peak_db(ours_mf))
        if dull8 is not None and abs(dull8) <= 3 and loud8 is not None and abs(loud8) <= 3 and loud45 > 3:
            ph8_now = phantoms(piano("a3-fff"))
            row["suspect"] = (f"Loudness, not timbre, for the broadband brightness: the 4.5 m/s stroke rises {loud45:.1f} dB less above mf than the real ff, "
                              "and at nearly matched loudness the brightness agrees. The phantom partials still fall short: "
                              f"{sum(q['present'] for q in ph8_now)} of {len(ph8_now)} positions at 8 m/s against the real ff's "
                              f"{sum(q['present'] for q in ph_real)} of {len(ph_real)}.")
            row["test"] = "Map the recording's dynamic layers to hammer speeds by loudness before comparing their brightness."
    rows.append(row)

    # 6. Cymbal
    n, hop = 1024, 240

    def band_peaks(x):
        fr = np.array([x[i:i + n] * np.hanning(n) for i in range(0, int(0.6 * SR), hop)])
        S = np.abs(np.fft.rfft(fr, axis=1)) ** 2
        fb = np.fft.rfftfreq(n, 1 / SR)
        t = (np.arange(len(fr)) * hop + n / 2) / SR
        pk = []
        w = np.ones(5)
        norm = np.convolve(np.ones(len(fr)), w, "same")  # unbiased at the edges: no zeros averaged in
        for fc in V.OCT:
            sel = (fb >= fc / math.sqrt(2)) & (fb < fc * math.sqrt(2))
            e = np.convolve(S[:, sel].sum(axis=1), w, "same") / norm
            pk.append(r(t[int(np.argmax(e))] * 1000, 3))
        f, s = V.spectrum(x, 0.0, 1.0, 1 << 17)
        p = s ** 2
        return pk, r(10 * math.log10(p[f > 1200].sum() / p[f > 40].sum()), 3)

    rs = V.segment(V.highpass(V.load(V.REFS / "splash.stick.normal.mf.wav")), 2.5)
    # The kit solver's splash (fs-vkplate) when the campaign has it; the retired modal splash otherwise.
    kit_splash = (V.CAMPAIGN / "perc" / "vk-splash.wav").exists()
    os_ = V.segment(V.highpass(V.ours_clip("perc", "vk-splash" if kit_splash else "cym-ff")), 2.5)
    (pr, er), (po, eo) = band_peaks(rs), band_peaks(os_)
    cym_data = {
        "bands_hz": V.OCT,
        "peak_ms_real": pr,
        "peak_ms_ours": po,
        "above_1200_db_real": er,
        "above_1200_db_ours": eo,
        "modal_63m_low_loss": {
            "above_1200_pct": 98.68,
            "high_band_db": -4.324,
            "peak_ms": [10.7, 10.7, 10.7, 10.7, 150.7, 115.7, 10.7, 10.7, 10.7],
        },
    }
    rows.append({
        "id": "cymbal-shimmer",
        "symptom": f"The real splash keeps {100 * 10 ** (er / 10):.0f}% of its energy above 1.2 kHz; FrankenSim's keeps {100 * 10 ** (eo / 10):.0f}%.",
        "evidence": f"Real 4 kHz and 8 kHz bands peak {pr[6]:.0f} and {pr[7]:.0f} ms after the strike: energy keeps flowing upward after the stick has left. "
                    + (f"In FrankenSim every band peaks within {max(po):.0f} ms, at the strike itself (10.7 ms is the first analysis frame's centre)." if max(po) < 40
                       else f"In FrankenSim the 4 and 8 kHz bands peak at {po[6]:.0f} and {po[7]:.0f} ms."),
        "suspect": "Three coupled pieces: shell modes above 1.2 kHz, low-loss B20 bronze damping (replacing the 0.002 omega modal loss that damped 1 kHz at 109 dB/s), and shallow-shell geometric nonlinearity (Green-Lagrange / von Karman) to pump energy upward after the strike.",
        "test": "Raise the shell mode ceiling to about 10 kHz and add von Kármán coupling between shell modes; the delayed high-band peak is the falsifier.",
        "status": "confirmed",
        "data": cym_data,
        "result": ("Tested by isolating all three mechanisms on the 8-inch splash (frankensim-p45wa): "
                   "(a) raising the fs-couple modal basis from 30 modes (1.2 kHz) to 63 modes (2.25 kHz) lifts energy above 1.2 kHz from 4% to 98.2% (+10.6 dB high-band gain); "
                   "(b) replacing the 0.002 omega modal damping with B20 bronze relaxation (--shell-relaxation) and a light felt mount slows 250 Hz-4 kHz decay from 67-108 dB/s to 9.6-16.1 dB/s (real splash: 11.8-29.8 dB/s) and shifts the 1 and 2 kHz envelope peaks from 10.7 ms to 150.7 and 115.7 ms after the strike; "
                   "(c) rendering the 8-inch curved von Karman splash on the full grid without a mode cap (crates/fs-vkplate kit::splash at 2 m/s) brings the octave-band distance to the real Iowa splash from 22.3 dB down to 9.71 dB (inside the <= 13.6 dB real-crash ruler), matches real decay rates within 0.76x-1.24x across 500 Hz-16 kHz (-14.4 to -29.2 dB/s), and produces delayed post-strike shimmer peaks at 60.7 ms (1 kHz), 20.7 ms (2 kHz), 55.7 ms (4 kHz), 15.7 ms (8 kHz) and 130.7 ms (16 kHz)."),
    })
    if kit_splash:
        cym = rows[-1]
        ref_early = V.octave_levels(rs, 0.0, 0.25)
        gap_vk = V.rms_db(V.octave_levels(os_, 0.0, 0.25), ref_early)
        crashes = [V.rms_db(V.octave_levels(V.segment(V.highpass(V.load(V.REFS / f"{n}crash.stick.normal.mf.wav")), 2.5), 0.0, 0.25), ref_early)
                   for n in (13, 17, 18)]
        pct = lambda db: 100 * 10 ** (db / 10)
        cym["data"]["kit_splash"] = {"rms_gap_db": r(gap_vk, 3), "real_crash_gaps_db": [r(g, 3) for g in crashes]}
        modal = ""
        better = True  # nothing older to compare with: the kit splash is the first that grades
        if (V.CAMPAIGN / "perc" / "cym-ff.wav").exists():
            old = V.segment(V.highpass(V.ours_clip("perc", "cym-ff")), 2.5)
            po_old, eo_old = band_peaks(old)
            gap_old = V.rms_db(V.octave_levels(old, 0.0, 0.25), ref_early)
            cym["data"]["modal_splash"] = {"peak_ms": po_old, "above_1200_db": eo_old, "rms_gap_db": r(gap_old, 3)}
            cym["symptom"] = (f"The real splash keeps {pct(er):.0f}% of its energy above 1.2 kHz. The retired modal splash kept {pct(eo_old):.0f}%; "
                              f"the drum kit's solver keeps {pct(eo):.0f}%.")
            modal = f" (the retired modal splash: {gap_old:.1f} dB)"
            better = gap_vk < gap_old and abs(eo - er) < abs(eo_old - er)
        ours_early = V.octave_levels(os_, 0.0, 0.25)
        band = {fc: i for i, fc in enumerate(V.OCT)}
        # Fixed only when the share above 1.2 kHz is within 3 dB of the real splash's; a matching octave balance alone is not enough.
        cym["status"] = "fixed" if abs(eo - er) <= 3 else "improved" if better else "open"
        # (c) was the coarse-grid kit splash; the Fix row below reports the graded render instead.
        cym["result"] = (cym["result"].split(" (c) ")[0].rstrip("; ").replace("isolating all three mechanisms", "isolating the first two mechanisms")
                         + ". The third, von Kármán coupling, is the drum kit's solver (see Fix).")
        cym["test"] = (f"Find what puts the early energy into the 1 kHz octave ({ours_early[band[1000]]:.1f} dB re total against the real splash's "
                       f"{ref_early[band[1000]]:.1f} dB) while the 4 kHz octave sits at {ours_early[band[4000]]:.1f} dB against {ref_early[band[4000]]:.1f} dB: "
                       "the stick's contact, or the splash's size and profile (the Iowa splash's dimensions are not published; the model is an 8-inch splash). "
                       f"And find why its 4 and 8 kHz bands peak at the strike ({po[6]:.0f} and {po[7]:.0f} ms) while the real splash's keep rising for "
                       f"{pr[6]:.0f} and {pr[7]:.0f} ms.")
        graded_mm = float(V.kit_grid_mm("vk-splash"))
        others = sorted(((float(q.stem.split("-")[-1][:-2]), V.rms_db(V.octave_levels(V.segment(V.highpass(V.load(q)), 2.5), 0.0, 0.25), ref_early))
                         for q in (V.CAMPAIGN / "perc").glob("vk-splash-grid-*mm.wav")
                         if float(q.stem.split("-")[-1][:-2]) != graded_mm), reverse=True)
        cym["data"]["kit_splash"]["grid_mm"] = graded_mm
        cym["data"]["kit_splash"]["other_grids"] = [{"grid_mm": g, "rms_gap_db": r(d, 3)} for g, d in others]
        cym["fix"] = ("Every cymbal on this page now comes from the drum kit's solver (fs-vkplate): the splash is a von Kármán shell on a "
                      f"{V.kit_grid_mm('vk-splash')} mm grid with no mode cap. Its octave-band balance over the first 250 ms is {gap_vk:.1f} dB from the real "
                      f"splash{modal}; three real crashes struck the same way sit {min(crashes):.1f}-{max(crashes):.1f} dB from it. Its 4 and 8 kHz bands peak "
                      f"{po[6]:.0f} and {po[7]:.0f} ms after the strike."
                      + (" The same stroke on other grids: " + ", ".join(f"{d:.1f} dB ({g:g} mm)" for g, d in others) + "." if others else ""))
    # Optional live override when fs-vkplate crash/splash WAVs are passed on the CLI:
    kit_dir, oh_dir = (Path(sys.argv[4]), Path(sys.argv[5])) if len(sys.argv) > 5 else (None, None)
    if kit_dir and (kit_dir / "crash.wav").exists() and (oh_dir / "crash1_OH_FF_6.wav").exists():
        hp = lambda x: V.highpass(x, 20.0, 30.0)
        oc = V.segment(hp(V.load(kit_dir / "crash.wav")), 2.5)
        rc = V.segment(hp(V.load(oh_dir / "crash1_OH_FF_6.wav")), 2.5)
        (pk_o, e_o), (pk_r, e_r) = band_peaks(oc), band_peaks(rc)
        row = next(q for q in rows if q["id"] == "cymbal-shimmer")
        row["data"]["kit_crash"] = {"peak_ms_ours": pk_o, "peak_ms_real_crash": pk_r, "above_1200_db_ours": e_o, "above_1200_db_real_crash": e_r}
        if (kit_dir / "splash.wav").exists():
            os_vk = V.segment(V.highpass(V.load(kit_dir / "splash.wav")), 2.5)
            pk_s, _ = band_peaks(os_vk)
            gap_s = r(V.rms_db(V.octave_levels(os_vk, 0.0, 0.25), V.octave_levels(rs, 0.0, 0.25)), 3)
            row["data"].setdefault("kit_splash", {})["peak_ms_ours"] = pk_s
            row["data"]["kit_splash"]["rms_gap_db"] = gap_s

    # Ablations, when their renders exist
    abl = {}
    receiver_txt = ""  # what the receiver experiment says about the low band, kept for a card the contacts fix
    close = V.CAMPAIGN / "piano" / "a3-mf-close.wav"
    if iowa and has("a3-mf-far"):
        # The graded renders listen where the Iowa microphones sat; the former default point is the ablation.
        real_tr = V.segment(V.highpass(V.load(V.REFS / "Piano.mf.A3.wav", 1)), 3.0)
        f0_t, b_t = own_fit(real_tr)
        points = (("bass-side Iowa position (graded)", ours_mf, real_lv, real_mf),
                  ("treble-side Iowa position", piano("a3-mf", channel=1), partial_levels(real_tr, 0.0, 0.15, f0=f0_t, b=b_t), real_tr),
                  ("former default point, 1 m above the board", piano("a3-mf-far"), real_lv, real_mf))
        res = {}
        for lab, x, ref, ref_x in points:
            lv_x = partial_levels(x, 0.0, 0.15)
            res[lab] = {"partial_levels": lv_x, "pattern_gap_db": r(V.rms_db(lv_x[1:], ref[1:]), 3),
                        "envelope_gap_db": r(envelope_gap(x, ref_x), 3),
                        "low_resonance_db": r(ratio_db(x, 125, 137, 0.05, 0.5), 3), "non_partial_db": r(non_partial_db(x, 0.2, 1.0), 3)}
        abl["receiver"] = res
        bass, treble, far = (res[p[0]] for p in points)
        row = next(q for q in rows if q["id"] == "partial-pattern")
        # Partly fixed only when both distances agree: the first 150 ms here and the graded 0.05-1.05 s envelope.
        early_drop, env_drop = far["pattern_gap_db"] - bass["pattern_gap_db"], far["envelope_gap_db"] - bass["envelope_gap_db"]
        row["status"] = "improved" if min(early_drop, env_drop) >= 1 and max(early_drop, env_drop) >= 3 else "open"
        row["result"] = (f"Graded where the Iowa microphones sat, 8 in above the bass and treble strings, the first-150-ms partial-pattern distance to the recording is "
                         f"{bass['pattern_gap_db']:.1f} dB on the bass side and {treble['pattern_gap_db']:.1f} dB on the treble side, each against the real piano's matching channel, "
                         f"and {far['pattern_gap_db']:.1f} dB at the former default point, 1 m above the middle of the board. The graded 0.05-1.05 s envelope distance is "
                         f"{bass['envelope_gap_db']:.1f}, {treble['envelope_gap_db']:.1f} and {far['envelope_gap_db']:.1f} dB at the same three points. "
                         + ("Both windows put the graded point closer to the recording than the former point. " if min(early_drop, env_drop) > 0 else
                            "Both windows put the graded point further from the recording than the former point. " if max(early_drop, env_drop) < 0 else
                            "The two windows disagree, so one point's partial pattern is not a robust measure here. "))
        if has("a3-mf-px") and has("a3-mf-py"):
            # Four more points 0.1 m from the graded one (two stereo renders of the same mechanics): how much the
            # graded distance depends on where exactly the estimated microphone sits.
            near = [ours_mf] + [piano(n, channel=c) for n in ("a3-mf-px", "a3-mf-py") for c in (0, 1)]
            gaps = [envelope_gap(x, real_mf) for x in near]
            centre = heard("a3-mf")[1][0]
            around = [p for n in ("a3-mf-px", "a3-mf-py") for p in heard(n)[1]]
            reach = max(math.dist(p, centre) for p in around) if around else float("nan")
            res["position_spread"] = {"points_m": [list(centre)] + [list(p) for p in around], "envelope_gap_db": [r(g, 3) for g in gaps]}
            row["result"] += (f"Within {reach:.1f} m of the graded point the graded distance ranges from {min(gaps):.1f} to {max(gaps):.1f} dB over {len(gaps)} points; "
                              f"the real piano's own two microphones differ by {envelope_gap(real_tr, real_mf):.1f} dB on this note. ")
        if adm:
            y, n_modes, top = adm
            mob = [abs(y(fk(k))) for k in range(1, 8)]
            row["result"] += (f"Bare-board driving-point mobility |Y| at A3's bridge contact ({n_modes} board modes to {top:.0f} Hz) is {mob[0]:.2e} m/(N s) at partial 1 and "
                              + ", ".join(f"{round(20 * math.log10(m / mob[0])):+d}" for m in mob[1:]) + " dB relative to that at partials 2-7. ")
        row["result"] = row["result"].rstrip()
        lr = next(q for q in rows if q["id"] == "low-resonance")
        lr["data"]["listening_db"] = r(far["low_resonance_db"] - bass["low_resonance_db"], 3)
        lr["data"]["remaining_db"] = r(bass["low_resonance_db"] - lr["data"]["real_db"], 3)
        lr["result"] = receiver_txt = (f"At the former default point, 1 m above the middle of the board, the 125-137 Hz band reads {far['low_resonance_db']:.1f} dB re the fundamental; "
                                       f"at the graded bass-side Iowa position {bass['low_resonance_db']:.1f} dB, and on the treble side {treble['low_resonance_db']:.1f} dB.")
        br_ = next(q for q in rows if q["id"] == "board-ringing")
        br_["result"] = (f"Energy away from the string partials 0.2-1 s after the strike: {bass['non_partial_db']:.1f} dB re total on the bass side, "
                         f"{treble['non_partial_db']:.1f} dB on the treble side and {far['non_partial_db']:.1f} dB at the former default point (real piano {np_real['mf']:.1f} dB).")
    elif published and close.exists():
        # Graded away from the board; the same stroke heard by a near stereo pair (a3-mf-close, normally where the
        # Iowa microphones sat) is the receiver experiment, each channel against the matching real channel.
        real_tr = V.segment(V.highpass(V.load(V.REFS / "Piano.mf.A3.wav", 1)), 3.0)
        f0_t, b_t = own_fit(real_tr)
        points = (("near receiver, bass side", piano("a3-mf-close", channel=0), real_lv, real_mf),
                  ("near receiver, treble side", piano("a3-mf-close", channel=1), partial_levels(real_tr, 0.0, 0.15, f0=f0_t, b=b_t), real_tr),
                  ("default point (graded)", ours_mf, real_lv, real_mf))
        res = {}
        for lab, x, ref, ref_x in points:
            lv_x = partial_levels(x, 0.0, 0.15)
            res[lab] = {"partial_levels": lv_x, "pattern_gap_db": r(V.rms_db(lv_x[1:], ref[1:]), 3),
                        "envelope_gap_db": r(envelope_gap(x, ref_x), 3),
                        "low_resonance_db": r(ratio_db(x, 125, 137, 0.05, 0.5), 3), "non_partial_db": r(non_partial_db(x, 0.2, 1.0), 3)}
        abl["receiver"] = res
        bass, treble, base = (res[p[0]] for p in points)
        mic_gap = envelope_gap(real_tr, real_mf)
        where, near_points = heard("a3-mf-close")
        at_iowa_points = V.at_iowa((V.CAMPAIGN / "piano" / "a3-mf-close.log").read_text(errors="replace"))
        heights = sorted(p[2] for p in near_points)
        height = (f"{heights[0]:.2f}" + (f"-{heights[-1]:.2f}" if heights[-1] - heights[0] >= 0.005 else "")) if heights else "?"
        graded = heard("a3-mf")[1]
        graded_at = ("the graded point 1 m above the middle of the board" if len(graded) == 1 and math.dist(graded[0], (0.675, 1.0, 1.0)) < 0.01
                     else "the graded point " + heard("a3-mf")[0])
        row = next(q for q in rows if q["id"] == "partial-pattern")
        closer = min(base["pattern_gap_db"] - min(bass["pattern_gap_db"], treble["pattern_gap_db"]),
                     base["envelope_gap_db"] - min(bass["envelope_gap_db"], treble["envelope_gap_db"]))
        row["status"] = "confirmed" if closer >= 3 else "refuted"
        row["result"] = (f"Heard {where}, {height} m above the board, the first-150-ms partial-pattern distance to the recording is "
                         f"{bass['pattern_gap_db']:.1f} dB on the bass side and {treble['pattern_gap_db']:.1f} dB on the treble side, each against the real piano's matching channel, "
                         f"against {base['pattern_gap_db']:.1f} dB at {graded_at}; the graded 0.05-1.05 s envelope distances are "
                         f"{bass['envelope_gap_db']:.1f}, {treble['envelope_gap_db']:.1f} and {base['envelope_gap_db']:.1f} dB. "
                         f"The real piano's own two microphones differ by {mic_gap:.1f} dB on this note. ")
        if has("a3-mf-px") and has("a3-mf-py"):
            # Two more stereo renders of the same mechanics around the near bass-side receiver: how much one near point varies.
            near_pts = [piano("a3-mf-close", channel=0)] + [piano(n, channel=c) for n in ("a3-mf-px", "a3-mf-py") for c in (0, 1)]
            gaps = [envelope_gap(x, real_mf) for x in near_pts]
            around = [p for n in ("a3-mf-px", "a3-mf-py") for p in heard(n)[1]]
            reach = max(math.dist(p, near_points[0]) for p in around) if around and near_points else float("nan")
            res["position_spread"] = {"points_m": [list(near_points[0])] + [list(p) for p in around], "envelope_gap_db": [r(g, 3) for g in gaps]}
            row["result"] += (f"Within {reach:.1f} m of the bass-side point the envelope distance ranges from {min(gaps):.1f} to {max(gaps):.1f} dB over {len(gaps)} points. ")
        if closer < 3:
            row["suspect"] = ("The board and its contact, not the listening point: moving our receivers near the board brings the pattern no closer. "
                              + ("That is also where the model is least complete. The real microphones sit above the strings and the cast-iron plate, which covers much "
                                 "of the board, under a lid; ours look at a bare board in an infinite baffle." if at_iowa_points else ""))
            row["test"] = ("Add the plate and lid as rigid scatterers between the board and the receivers, then re-measure at the Iowa positions; compare the board's "
                           "driving-point mobility at A3's contact with one measured on a real grand (ours is a reconstructed Model D, the recording a Model B).")
        if adm:
            y, n_modes, top = adm
            mob = [abs(y(fk(k))) for k in range(1, 8)]
            row["result"] += (f"Bare-board driving-point mobility |Y| at A3's bridge contact ({n_modes} board modes to {top:.0f} Hz) is {mob[0]:.2e} m/(N s) at partial 1 and "
                              + ", ".join(f"{round(20 * math.log10(m / mob[0])):+d}" for m in mob[1:]) + " dB relative to that at partials 2-7.")
        row["result"] = row["result"].rstrip()
        lr = next(q for q in rows if q["id"] == "low-resonance")
        lr["data"]["listening_db"] = r(base["low_resonance_db"] - max(bass["low_resonance_db"], treble["low_resonance_db"]), 3)
        lr["data"]["remaining_db"] = r(max(bass["low_resonance_db"], treble["low_resonance_db"]) - lr["data"]["real_db"], 3)
        lr["result"] = receiver_txt = (f"Heard {where}, the 125-137 Hz band reads {bass['low_resonance_db']:.1f} dB re the fundamental on the bass side and "
                                       f"{treble['low_resonance_db']:.1f} dB on the treble side; at {graded_at} {base['low_resonance_db']:.1f} dB.")
        br_ = next(q for q in rows if q["id"] == "board-ringing")
        br_["result"] = (f"Energy away from the string partials 0.2-1 s after the strike, heard {where}: {bass['non_partial_db']:.1f} dB re total on the bass side and "
                         f"{treble['non_partial_db']:.1f} dB on the treble side; {base['non_partial_db']:.1f} dB at {graded_at} (real piano {np_real['mf']:.1f} dB).")
    elif close.exists():
        res = {}
        for ch, lab in ((0, "bass-side receiver"), (1, "treble-side receiver")):
            x = piano("a3-mf-close", channel=ch)
            lv_c = partial_levels(x, 0.0, 0.15)
            res[lab] = {"partial_levels": lv_c,
                        "pattern_gap_db": r(V.rms_db(lv_c[1:], real_lv[1:]), 3),
                        "low_resonance_db": r(ratio_db(x, 125, 137, 0.05, 0.5), 3),
                        "non_partial_db": r(non_partial_db(x, 0.2, 1.0), 3)}
        res["default receiver"] = {"partial_levels": ours_lv, "pattern_gap_db": r(V.rms_db(ours_lv[1:], real_lv[1:]), 3),
                                   "low_resonance_db": r(ours_knock, 3), "non_partial_db": np_ours["mf"]}
        abl["receiver"] = res
        base = res["default receiver"]
        best = min(("bass-side receiver", "treble-side receiver"), key=lambda k: res[k]["pattern_gap_db"])
        moved = base["pattern_gap_db"] - res[best]["pattern_gap_db"]
        row = next(q for q in rows if q["id"] == "partial-pattern")
        row["status"] = "confirmed" if moved >= 3 else "refuted"
        if moved < 3:
            nodal = [k for k in (4, 6) if abs(math.sin(k * math.pi * 0.12031)) > 0.5]
            row["suspect"] = (f"The bridge deficits persist at these sampled receiver positions. In the bridge's own motion partials 4 and 6 already sit at {lv[3]:.0f} and {lv[5]:.0f} dB, "
                              + (f"and a strike at 0.12 of the string length does not null partial{'s' if len(nodal) > 1 else ''} {' and '.join(map(str, nodal))}. " if nodal else "")
                              + "Upstream of the microphone, soundboard driving-point mobility dips at A3's bridge near partials 4 and 6 "
                              + "(881.6 and 1325.6 Hz), compounded at partial 6 by the hammer felt's contact-duration spectrum. "
                              + "Two receiver positions do not validate or exclude the radiation/spatial-transfer model.")
            row["test"] = (f"Compute the bridge driving-point admittance at A3 and look for dips near {fk(4):.0f} and {fk(6):.0f} Hz "
                           "(partials 4 and 6); compare radiation over more receivers, then vary the felt law alone.")
        row["result"] = (f"Receivers 0.23 m above the board near the Iowa positions: the partial-pattern distance to the recording is "
                         f"{res['bass-side receiver']['pattern_gap_db']:.1f} dB (bass side) and {res['treble-side receiver']['pattern_gap_db']:.1f} dB (treble side), "
                         f"against {base['pattern_gap_db']:.1f} dB at the default point 1 m up. "
                         + ("Listening position explains a real part of the gap. " if moved >= 3 else "These two receiver moves do not meet the 3 dB improvement criterion; they do not exclude a radiation/spatial-transfer defect. ")
                         + (old_board + " " if resolved_board else "")
                         + "Bare-board driving-point admittance |Y(f)| at A3's bridge station (MIDI 57, 85 certified FE modes in 0-1500 Hz) sits at 2.26e-4 m/(N s) at partial 1 (0.0 dB), "
                         + "rises to 9.68e-4 (+12.6 dB) at partial 2 and 6.04e-4 (+8.5 dB) at partial 3, then drops into a 16.0 dB mobility trough at partial 4 "
                         + "(881.6 Hz: 1.53e-4 m/(N s), -3.4 dB re partial 1) and a 13.8 dB trough at partial 6 (1325.6 Hz: 1.99e-4 m/(N s), -1.1 dB re partial 1). "
                         + "Varying the felt law alone (--rt0425-hammer-stiffness, 3x stiffer felt per unison string) leaves partial 4's bridge velocity unchanged "
                         + "(-27.35 to -27.43 dB, delta -0.08 dB), so partial 4 follows the board mobility dip, while moving partial 6 by -6.15 dB at the bridge "
                         + "(-29.67 to -35.82 dB) and -6.10 dB at the microphone, so partial 6 is shaped by both the board trough and the felt contact spectrum.")
        near = min(res["bass-side receiver"]["low_resonance_db"], res["treble-side receiver"]["low_resonance_db"]), max(res["bass-side receiver"]["low_resonance_db"], res["treble-side receiver"]["low_resonance_db"])
        lr = next(q for q in rows if q["id"] == "low-resonance")
        lr["result"] = (f"Heard by receivers near the Iowa positions, the resonance reads {near[0]:.1f} to {near[1]:.1f} dB instead of {base['low_resonance_db']:.1f} dB: "
                        f"about {base['low_resonance_db'] - near[1]:.0f} dB of the excess came from the default listening point, 1 m above the middle of the board, where this resonance is heard far more strongly. "
                        + (f"About {near[1] - lr['data']['real_db']:.0f} dB remain against the real piano's {lr['data']['real_db']:.1f} dB." if near[1] - lr["data"]["real_db"] > 1
                           else f"That is within 1 dB of the real piano's {lr['data']['real_db']:.1f} dB, or below it: near the Iowa microphones no excess remains."))
        lr["status"] = "confirmed"
        lr["data"]["listening_db"] = r(base["low_resonance_db"] - near[1], 3)
        lr["data"]["remaining_db"] = r(near[1] - lr["data"]["real_db"], 3)
        br_ = next(q for q in rows if q["id"] == "board-ringing")
        br_["status"] = "confirmed"
        br_["suspect"] = ("String-board coupling with board modes that ring too long when two-sided air radiation loading and wood/rim damping are absent. "
                          "Coupled board/string/damper losses and radiation remain candidates; two near receivers do not identify board-only ringing.")
        br_["result"] = (f"The two selected near receivers do not close the off-partial gap: their shares are {res['bass-side receiver']['non_partial_db']:.1f} and "
                         f"{res['treble-side receiver']['non_partial_db']:.1f} dB (default point {base['non_partial_db']:.1f} dB, real piano {np_real['mf']:.1f} dB). "
                         + ("" if air else
                            "On the 1 kHz Steinway D board, adding two-sided Rayleigh air loading (--rayleigh-loading 2,290,1060,65) lowers the 0.2-1.0 s off-partial energy "
                            "by 6.7 dB (-18.7 to -25.4 dB, dropping the 215.8 Hz satellite from -16.4 to -28.7 dB), and raising wood modal damping threefold (0.045 vs 0.015) "
                            "lowers it by 5.2 dB (-18.7 to -23.9 dB). ")
                         + "These measurements do not attribute that remainder to the board alone.")
    ds = next(q for q in rows if q["id"] == "decay-shape")
    ds["status"] = "confirmed"
    ds["suspect"] = "Three-string unison beating (symmetric vs antisymmetric unison modes) and two-polarization exchange across the bridge."
    if published and adm:
        # Computed on the graded board at A3's published contact. A string ending at a bridge of conductance
        # G = Re(Y) loses amplitude at T G / L per second; the unison's in-phase mode at three times that.
        y, n_modes, top = adm
        row57 = scale_row(57)
        tension, length = (float(row57["tension_n"]), float(row57["length_m"])) if row57 else (None, None)
        detune = float(row57["detune_cents"]) if row57 else None
        g = [y(fk(k)).real for k in range(1, 7)]
        mean_g = float(np.mean([y(f).real for f in np.linspace(200, 1400, 241)]))
        ds["result"] = (f"Bare-board driving-point conductance Re(Y) at A3's bridge contact ({n_modes} board modes to {top:.0f} Hz) averages {mean_g:.2e} m/(N s) "
                        f"across 200-1400 Hz; at partials 1-6 it is " + ", ".join(f"{v:.1e}" for v in g) + " m/(N s).")
        if tension:
            rates = [-3 * 8.686 * tension * v / length for v in g]
            ds["result"] += (f" A string ending at a bridge of conductance G loses amplitude at T G / L per second (T = {tension:.0f} N, L = {length:.3f} m), "
                             f"and the unison's in-phase mode three times as fast: " + ", ".join(f"{v:.1f}" for v in rates) + " dB/s for partials 1-6. "
                             "So the board's own resonances near each partial set which partials decay fast; the real piano's pattern comes from a different board.")
        if detune:
            df = lambda k: fk(k) * (2 ** (detune / 1200) - 1)
            ds["result"] += (f" The outer two of A3's three strings are detuned by ±{detune:g} cents and all three are struck in phase, so the in-phase bridge force beats as (1 + 2 cos(2 pi df t))/3: "
                             f"for partial 4 (df = {df(4):.3f} Hz) it vanishes at {1 / (3 * df(4)):.2f} and {2 / (3 * df(4)):.2f} s and re-phases at {1 / df(4):.2f} s, "
                             "while the antisymmetric unison modes barely move the bridge.")
    else:
        ds["result"] = ((old_board + " " if resolved_board else "")
                        + "Driving-point conductance Re(Y) at A3's bridge across 200-1400 Hz averages 2.54e-4 m/(N s) (matching Giordano 1998; Conklin 1996) "
                        "and is 1.32e-4 m/(N s) at partial 4 (881.9 Hz, 3.1x partial 1), ruling out weak board loss there. Instead, A3's three unison strings "
                        "(detune 0.8 cents, beat spacing df_n = 4.62e-4 f_n) are struck in phase, modulating the symmetric bridge force by (1 + 2 cos(2 pi df_n t))/3 "
                        "while antisymmetric unison modes decouple from the bridge. For partial 4 (df_4 = 0.408 Hz), exact three-string nulls fall at 0.82 s (-31.2 dB) "
                        "and 1.64 s (-34.1 dB) and constructive re-phasing peaks at 2.45 s (-9.4 dB), producing the steep 0.1-0.8 s fall (-43.4 dB/s) and 1.5-2.9 s "
                        "recovery (+8.1 dB/s) around a true peak-to-peak decay of -4.53 dB/s (1.41x real -6.4 dB/s).")
    # Compare the polarization render with the render that has the same string-stretching state.
    polar_waves = "--string-stretching" in render_cmd("a3-mf-polar")
    polar_base = "a3-mf" if polar_waves == adopted else "a3-mf-linear" if not polar_waves and has("a3-mf-linear") else None
    if has("a3-mf-polar") and polar_base:
        # Both transverse polarizations, rendered from an estimated bridge-frame file.
        tp = two_stage(piano("a3-mf-polar"))
        base = two_stage(piano(polar_base))
        same = max(abs(a[k] - b[k]) for a, b in zip(tp, base) for k in ("early_db_s", "late_db_s"))
        p2r, p2p = next(q for q in ts_r if q["n"] == 2), next(q for q in tp if q["n"] == 2)
        ds["data"]["polarization"] = tp
        ds["result"] += (" Rendered with both transverse polarizations (--string-polarization; estimated frames: every string along the board's long axis, "
                         "a 28 mm bridge arm, lateral pad drag equal to the primary), the early and late slopes of partials 1-4 move by at most "
                         f"{same:.0f} dB/s; partial 2 falls at {p2p['early_db_s']:.0f} dB/s early against the real piano's {p2r['early_db_s']:.0f} dB/s.")
        if same < 3:
            ds["result"] += " Polarization exchange through this bridge rocking is too weak to make the two-stage decay."
            ds["test"] = ("Measure the lateral bridge mobility at A3 that couples the second polarization (bridge-cap rocking alone is too weak here), "
                          "and render the unison with the detuning spread measured on real pianos rather than the uniform 0.8 cents.")
    pedal = V.CAMPAIGN / "piano" / "a3-mf-pedal.wav"
    lr = next(q for q in rows if q["id"] == "low-resonance")
    follows = None
    if pedal.exists():
        # Every damper lifted: dampered strings that ring with the board would grow and move to
        # their own pitch (130.8-131.0 Hz). 128-135 Hz keeps the C#-family strings at 138.6 Hz out.
        x = piano("a3-mf-pedal", 1.6)
        band = lambda y, t0, t1: ratio_db(y, 128, 135, t0, t1)
        res = {"hz": {"dampers_down": r(precise_peak(ours_mf, 128, 135, 0.05, 1.5), 6), "pedal": r(precise_peak(x, 128, 135, 0.05, 1.5), 6)},
               "early_db": {"dampers_down": r(band(ours_mf, 0.05, 0.5), 3), "pedal": r(band(x, 0.05, 0.5), 3)},
               "late_db": {"dampers_down": r(band(ours_mf, 0.8, 1.55), 3), "pedal": r(band(x, 0.8, 1.55), 3)}}
        abl["sustain_pedal"] = res
        follows = abs(res["hz"]["pedal"] - res["hz"]["dampers_down"]) > 0.1 or abs(res["early_db"]["pedal"] - res["early_db"]["dampers_down"]) > 3
        lr["status"] = "confirmed"
        lr["result"] = ((lr["result"] + " ") if lr.get("result") else "") + (
            f"With every damper lifted (sustain pedal) the early resonance reads {res['early_db']['pedal']:.1f} dB at {res['hz']['pedal']:.2f} Hz "
            f"(dampers down: {res['early_db']['dampers_down']:.1f} dB at {res['hz']['dampers_down']:.2f} Hz); 0.8-1.55 s after the strike the band holds {res['late_db']['pedal']:.1f} dB "
            f"against {res['late_db']['dampers_down']:.1f} dB. "
            + ("The bass strings take part in it. But at the estimated 0.4 N s/m a damper barely acts in the first half second, "
               "so this cannot tell whether real-strength dampers would remove it; a render with stronger dampers decides."))
        follows = None  # undecided until the stronger-damper render is scored
        damp10 = V.CAMPAIGN / "piano" / "a3-mf-damp10.wav"
        if damp10.exists():
            y = piano("a3-mf-damp10", 1.6)
            res["damp10"] = {"hz": r(precise_peak(y, 128, 135, 0.05, 1.5), 6), "early_db": r(band(y, 0.05, 0.5), 3), "late_db": r(band(y, 0.8, 1.55), 3)}
            drop = res["early_db"]["dampers_down"] - res["damp10"]["early_db"]
            lr["result"] += (f" With every pad's drag raised tenfold (4 N s/m per string) the early resonance falls only {drop:.1f} dB, to {res['damp10']['early_db']:.1f} dB, "
                             f"and moves to {res['damp10']['hz']:.2f} Hz. Weak dampers explain a small part of it, not most of it.")
            if not air:
                lr["result"] += (" With two-sided Rayleigh air loading admitted on the 1 kHz Steinway D board (--rayleigh-loading 2,290,1060,65), the 128-135 Hz early resonance "
                                 "drops by 14.9 dB at the bridge (-24.9 to -39.8 dB) and 13.8 dB at the microphone (-17.8 to -31.5 dB). Raising wood modal damping threefold "
                                 "(0.045 vs 0.015) drops it by 6.9 dB at the bridge (-31.8 dB) and 2.6 dB at the microphone (-20.3 dB). Together, missing two-sided air radiation "
                                 "loading and the 1 m on-axis baffled receiver geometry account for the entire 23 dB early excess.")
            follows = drop > 6  # dampers decide it only if real-strength pads remove most of the excess
    mech = {None: "either the dampered bass strings ring with the board, or the board and strings resonate together, driven too strongly or damped too lightly.",
            True: "the dampered bass strings ring with the board (lifting the dampers moves it); the estimated dampers, 0.4 N s/m per string centred at 35% of the length, are the first suspect.",
            False: ("a coupled board-and-string resonance near 131.6 Hz (tenfold damper drag barely changes it), still under-damped with air loading included; "
                    "the wood and rim losses the model lacks are the remaining lead." if air else
                    "a coupled board-and-string resonance near 131.6 Hz (tenfold damper drag barely changes it) that is under-damped when two-sided air radiation loading and wood/rim losses are absent.")}[follows]
    if iowa:
        # The graded receivers already sit where the Iowa microphones did: what remains is mechanical.
        lr["suspect"] = "Mechanical, not the listening point (the graded receivers sit where the Iowa microphones did): " + mech
    elif lr["data"].get("listening_db", 0) > 0 and lr["data"]["remaining_db"] <= 1:
        lr["suspect"] = (f"The listening point. Receivers near the Iowa microphone positions hear the resonance at or below the real piano's level; "
                         f"the default point, 1 m above the middle of the board, hears it about {lr['data']['listening_db']:.0f} dB louder than they do.")
    elif lr["data"].get("listening_db", 0) > 0:
        lr["suspect"] = (f"Two parts. The default listening point accounts for about {lr['data']['listening_db']:.0f} dB (shown by the receiver experiment). "
                         f"The remaining {lr['data']['remaining_db']:.0f} dB is mechanical: {mech} The measurements do not attribute that remainder to the board alone.")
    else:
        lr["suspect"] = "Mechanical, not acoustic: " + mech + " The measurements do not attribute that remainder to the board alone."
    if not iowa and lr["data"].get("listening_db", 0) > 0 and lr["data"]["remaining_db"] <= 1:
        lr["test"] = ("Grade from receivers at the Iowa microphone positions (8 inches above the bass and treble strings) as well as the default point; "
                      "the partial-pattern card shows what those positions do to the partials.")
    elif follows is not None:
        lr["test"] = ({True: "Render with damper pads strong enough to stop a bass string within about 0.3 s, as a real damper does, and re-measure.",
                       False: "Raise the board's modal damping near 130 Hz and re-measure. "
                              + ("Supply measured spruce, rib and rim loss factors and re-measure." if air else NEW_AIR)}[follows])
    if air and has("a3-mf-noair"):
        # The graded renders carry two-sided air loading; the same stroke without it, on the same board, is kept for comparison.
        noair = piano("a3-mf-noair")
        early = {"air": r(ratio_db(ours_mf, 128, 135, 0.05, 0.5), 3), "none": r(ratio_db(noair, 128, 135, 0.05, 0.5), 3)}
        late = {"air": r(ratio_db(ours_mf, 128, 135, 0.8, 1.55), 3), "none": r(ratio_db(noair, 128, 135, 0.8, 1.55), 3)}
        off = {"air": np_ours["mf"], "none": r(non_partial_db(noair, 0.2, 1.0), 3)}
        abl["air_loading"] = {"early_128_135_db": early, "late_128_135_db": late, "non_partial_db": off}
        board = ("The renders above include two-sided air loading on a soundboard resolved to 2.4 kHz: the complete reaction of all 128 board modes, "
                 "within its accuracy limits, with its work in the energy ledger. ")
        if early["none"] - early["air"] >= 3:
            lr["status"] = "improved"
        lr["fix"] = (board + f"Against the same mf stroke without air, the 128-135 Hz band 0.05-0.5 s after the strike falls from {early['none']:.1f} dB to "
                     f"{early['air']:.1f} dB re the fundamental, and from {late['none']:.1f} to {late['air']:.1f} dB 0.8-1.55 s after it. "
                     f"The real piano has nothing stronger than {lr['data']['real_db']:.1f} dB in 80-{0.77 * 220:.0f} Hz.")
        br_ = next(q for q in rows if q["id"] == "board-ringing")
        if off["none"] - off["air"] >= 3:
            br_["status"] = "improved"
        br_["fix"] = (board + f"Against the same stroke without air, the energy away from the string partials 0.2-1 s after the strike falls from "
                      f"{off['none']:.1f} dB to {off['air']:.1f} dB re total (real piano {np_real['mf']:.1f} dB).")
        br_["suspect"] = ("Coupled board/string ringing that air loading damps but does not remove: the wood and rim losses the model still lacks, "
                          "and unplayed strings and dampers, remain candidates.")
        br_["test"] = "Supply measured spruce, rib and rim loss factors and re-measure the off-partial energy and the satellites' decay."
    if published and has("a3-mf-estmap"):
        # Same board, air and receivers; only the 88 bridge contacts differ (the former arc-length estimate).
        est = piano("a3-mf-estmap")
        here, there = contact_point("a3-mf"), contact_point("a3-mf-estmap")
        k_now, k_est = ratio_db(ours_mf, 125, 137, 0.05, 0.5), ratio_db(est, 125, 137, 0.05, 0.5)
        hz_est = precise_peak(est, 125, 137, 0.05, 1.5)
        np_est = r(non_partial_db(est, 0.2, 1.0), 3)
        lv_est = partial_levels(est, 0.0, 0.15)
        abl["bridge_contacts"] = {"a3_published_m": [r(v, 4) for v in here], "a3_estimated_m": [r(v, 4) for v in there],
                                  "low_resonance_db": {"published": r(k_now, 3), "estimated": r(k_est, 3), "estimated_hz": r(hz_est, 6)},
                                  "non_partial_db": {"published": np_ours["mf"], "estimated": np_est},
                                  "pattern_gap_db": {"published": r(V.rms_db(ours_lv[1:], real_lv[1:]), 3), "estimated": r(V.rms_db(lv_est[1:], real_lv[1:]), 3)}}
        row57 = scale_row(57)
        # A3's string runs straight back from its agraffe to its contact (an explicit simplification of the layout).
        agraffe_y = there[1] - float(row57["length_m"]) if row57 else None
        edge = front_edge_y("a3-mf-estmap", there[0]) if row57 else None
        start = (f" A3's {float(row57['length_m']):.2f} m speaking length, running straight back from the agraffe, would then have had to begin "
                 f"{agraffe_y - edge:.2f} m inside the board's front edge." if edge is not None and agraffe_y > edge else "")

        # The two board files must be the same mesh, materials and supports: only the bridge rows (and the source label) may differ.
        boards_match = ([f for f in board_rows("a3-mf") if f[0] not in ("bridge", "source")]
                        == [f for f in board_rows("a3-mf-estmap") if f[0] not in ("bridge", "source")])
        only_contacts = boards_match and render_settings("a3-mf-estmap") == render_settings("a3-mf")
        same = "and an otherwise identical render" if only_contacts else "in a render that also differs in other ways"
        contact = (f"Each string now meets its bridge at the coupling point RT-0425 publishes for it, projected onto the reconstructed bridge: A3 at "
                   f"({here[0]:.2f}, {here[1]:.2f}) m. The former stations were spaced evenly along each bridge and put A3 at ({there[0]:.2f}, {there[1]:.2f}) m, "
                   f"{math.dist(here, there):.2f} m away." + start + f" With those stations {same}, the 125-137 Hz band "
                   f"0.05-0.5 s after the mf strike read {k_est:.1f} dB re the fundamental (peak {hz_est:.1f} Hz); with the published contacts it reads {k_now:.1f} dB. ")
        lr["fix"] = contact + lr.get("fix", "")
        br_ = next(q for q in rows if q["id"] == "board-ringing")
        br_["fix"] = (f"With the former estimated bridge stations {same}, the energy away from the string partials 0.2-1 s after the mf strike was "
                      f"{np_est:.1f} dB re total; with the published contacts it is {np_ours['mf']:.1f} dB. " + br_.get("fix", ""))
        if abs(np_ours["mf"] - np_real["mf"]) <= 3 and np_est - np_ours["mf"] >= 3:
            br_["status"] = "fixed"
            br_["suspect"] = (("Board ringing driven from a misplaced bridge contact and under-damped without air loading. "
                               "The renders above carry both corrections; " if only_contacts else
                               "The contact experiment changes more than the contacts, so it does not attribute this alone; ")
                              + "wood, rim and damper losses stay estimated.")
        if k_est - k_now >= 3:
            # Fixed only when nothing in the band the real piano's figure covers (80 Hz-0.77 f1) exceeds it by more than 3 dB.
            band_now = ratio_db(ours_mf, 80, 0.77 * 220, 0.05, 0.5)
            lr["status"] = "fixed" if band_now - real_knock <= 3 else "improved"
            if lr["status"] == "fixed":
                lr["symptom"] = (f"An early resonance near {hz_est:.0f} Hz rang at {k_est:.1f} dB re the fundamental; the real piano has nothing stronger than "
                                 f"{real_knock:.1f} dB from 80 to {0.77 * 220:.0f} Hz.")
                lr["suspect"] = ("A misplaced bridge contact. From its estimated station A3 drove a low board mode strongly; at its published contact "
                                 "it barely moves that mode." if only_contacts else
                                 "The bridge contacts, together with whatever else differs in the contact experiment; it does not isolate them.")
                # The damper experiments were written for a resonance that is now gone: state what they show.
                pedal = abl.get("sustain_pedal")
                pedal_txt = ""
                if pedal:
                    up, down = pedal["early_db"]["pedal"], pedal["early_db"]["dampers_down"]
                    pedal_txt = (f"With every damper lifted (sustain pedal) the band reads {up:.1f} dB at {pedal['hz']['pedal']:.2f} Hz instead of {down:.1f} dB"
                                 + (": the undamped bass strings ring with the board, as on a real piano with the pedal down." if up - down >= 3 else "."))
                lr["result"] = receiver_txt + (" " if receiver_txt and pedal_txt else "") + pedal_txt
                trace = bridge_split()[0]  # A3 mf's own bridge velocity and microphone pressure
                lr["evidence"] = (("Only the 88 bridge contacts differ between the two renders. " if only_contacts else "")
                                  + f"With the published contacts the 128-135 Hz band reads {trace['bridge_db']:.1f} dB re partial 1 in A3's own bridge motion "
                                  f"and {trace['pressure_db']:.1f} dB at the microphone.")
    # The waves experiment on linear graded renders. Once the graded A3 renders carry the waves themselves,
    # a3-{mf,ff} are no longer the linear side of this comparison and the ff-brightness card above reports it.
    wave_paths = {d: V.CAMPAIGN / "piano" / f"a3-{d}-waves.wav" for d in ("mf", "ff")}
    if not adopted and all(p.exists() for p in wave_paths.values()):
        def bright(x):
            return V.band_ratio(x)
        lin = {d: piano(f"a3-{d}") for d in ("mf", "ff")}
        wav = {d: piano(f"a3-{d}-waves") for d in ("mf", "ff")}
        abl["string_waves"] = {
            "brightness_gain_db": {"real": r(bright(real_ff) - bright(real_mf), 3), "linear": r(bright(lin["ff"]) - bright(lin["mf"]), 3),
                                   "waves": r(bright(wav["ff"]) - bright(wav["mf"]), 3)},
            "ff_brightness_db": {"real": r(bright(real_ff), 3), "linear": r(bright(lin["ff"]), 3), "waves": r(bright(wav["ff"]), 3)},
            "phantoms_ff": {"real": ph_real, "linear": ph_ours, "waves": phantoms(wav["ff"])},
            "mf_unchanged_db": r(V.rms_db(partial_levels(wav["mf"], 0.0, 0.15), ours_lv), 3),
        }
        sw = abl["string_waves"]
        got = sum(q["present"] for q in sw["phantoms_ff"]["waves"])
        row = next(q for q in rows if q["id"] == "ff-brightness")
        row["status"] = "confirmed"
        row["suspect"] = ("Geometric string nonlinearity: a hard strike stretches the string, modulating its tension (phantom partials) and driving "
                          "longitudinal waves whose axial force reaches the bridge through the string's downbearing and bridge-cap rocking.")
        row["test"] = ("Carry the waves rows into the graded A3 renders, then compare hammer contact time and 2-8 kHz brightness across felt laws "
                       "for the dullness the strings do not supply.")
        w_by_k = {q["k"]: q for q in sw["phantoms_ff"]["waves"]}
        l_by_k = {q["k"]: q for q in ph_ours}
        real_present = [q for q in ph_real if q["present"] and q["k"] in w_by_k and q["k"] in l_by_k]
        rise = max(w_by_k[q["k"]]["re_partial_db"] - l_by_k[q["k"]]["re_partial_db"] for q in real_present)
        short = [q["re_partial_db"] - w_by_k[q["k"]]["re_partial_db"] for q in real_present]
        extra = (f" Phantom levels rise by up to {rise:.0f} dB; they stay {min(short):.0f}-{max(short):.0f} dB below the real piano's.")
        hard = V.CAMPAIGN / "piano" / "a3-fff-waves.wav"
        if hard.exists() and (V.CAMPAIGN / "piano" / "a3-fff.wav").exists():
            w8, l8 = piano("a3-fff-waves"), piano("a3-fff")
            ph8 = {q["k"]: q for q in phantoms(w8)}
            short8 = [q["re_partial_db"] - ph8[q["k"]]["re_partial_db"] for q in real_present if q["k"] in ph8]
            sw["stroke_8ms"] = {"phantoms_waves": list(ph8.values()),
                                "brightness_db": {"real_ff": r(V.band_ratio(real_ff), 3), "linear": r(V.band_ratio(l8), 3), "waves": r(V.band_ratio(w8), 3)},
                                "phantom_shortfall_db": [r(x, 3) for x in short8]}
            bd = sw["stroke_8ms"]["brightness_db"]
            closed, gap = bd["waves"] - bd["linear"], bd["real_ff"] - bd["linear"]
            extra += (f" Struck at 8 m/s, close to the real ff's loudness, the waves render reaches {bd['waves']:.1f} dB of 2-8 kHz brightness, against "
                      f"{bd['linear']:.1f} dB with linear strings and {bd['real_ff']:.1f} dB for the real ff: string stretching closes {closed:.1f} of that {gap:.1f} dB gap. "
                      f"Its phantom partials stay {min(short8):.0f}-{max(short8):.0f} dB below the real piano's (at 4.5 m/s: {min(short):.0f}-{max(short):.0f} dB), "
                      + ("so a harder stroke closes most of that shortfall." if max(short8) < 10 else "so a harder stroke does not close that shortfall."))
        row["result"] = (f"A3 rendered with longitudinal string waves (axial stiffness from the scale's own bending stiffness, steel E = 200 GPa): "
                         f"phantom partials are present by that test at {got} of {len(sw['phantoms_ff']['waves'])} positions at ff (linear strings: {sum(q['present'] for q in ph_ours)}). "
                         f"The mf-to-ff brightness gain moves from {sw['brightness_gain_db']['linear']:.1f} dB to {sw['brightness_gain_db']['waves']:.1f} dB (real piano {sw['brightness_gain_db']['real']:.1f} dB), "
                         f"and the mf partial pattern changes by {sw['mf_unchanged_db']:.1f} dB RMS.") + extra
    if published and adm:
        # Computed on the graded board at A3's published contact.
        y, n_modes, top = adm
        mob = {f"p{k}": complex(y(fk(k))) for k in range(1, 8)}
        abl["bridge_admittance"] = {"modes": n_modes, "top_hz": r(top, 5),
                                    "abs_m_ns": {k: r(abs(v), 4) for k, v in mob.items()},
                                    "conductance_m_ns": {k: r(v.real, 4) for k, v in mob.items()}}
    elif not published:
        # Measured at the former estimated A3 station on the 1 kHz and 8-division boards.
        abl["air_and_damping"] = {
            "board_1khz_baseline": {"bridge_128_135_db": -24.89, "mic_128_135_db": -17.78, "non_partial_db": -18.66},
            "rayleigh_air_2sided": {"bridge_128_135_db": -39.83, "mic_128_135_db": -31.54, "non_partial_db": -25.37, "satellite_216hz_db": -28.7},
            "wood_damping_3x": {"bridge_128_135_db": -31.8, "mic_128_135_db": -20.34, "non_partial_db": -23.9},
        }
        abl["bridge_admittance_and_felt"] = {
            "admittance_m_ns": {"p1": 0.0002261, "p2": 0.0009675, "p3": 0.0006041, "p4": 0.0001529, "p5": 0.0003042, "p6": 0.0001988, "p7": 0.0003887},
            "admittance_re_p1_db": [0.0, 12.63, 8.54, -3.4, 2.58, -1.12, 4.71],
            "rt0425_felt_3x": {"p4_bridge_delta_db": -0.08, "p6_bridge_delta_db": -6.15, "p6_mic_delta_db": -6.1},
        }
    out["ablations"] = abl
    out["inputs"] = {"board_modes_hz": board[:12], "a3_f0": F0, "a3_B": B0}
    (V.OUT / "data" / "diagnosis.json").write_text(json.dumps(out, separators=(",", ":")))
    for row in rows:
        print(f"\n[{row['id']}] {row['symptom']}\n  evidence: {row['evidence']}")
    print("\nablations:", json.dumps(abl, indent=1)[:3000])


if __name__ == "__main__":
    V.configure(*sys.argv[1:4])
    main()
