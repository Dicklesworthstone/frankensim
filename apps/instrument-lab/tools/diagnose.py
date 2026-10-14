"""Differential diagnosis: which link of the simulated signal chain makes FrankenSim sound
unlike the real instrument, and which experiment would decide it.

verisim.py grades how far apart two sounds are. This script asks WHY, by splitting the
difference along the chain hammer -> string -> bridge -> board -> radiation -> microphone,
using observables only a simulation has (the bridge velocity) and ablation renders that
switch one mechanism on at a time.

usage: uv run --with numpy --with pillow python tools/diagnose.py CAMPAIGN_DIR REF_WAV48_DIR PUBLIC_DIR

Inputs beyond verisim.py's:
  piano/a3-mf-bridge.csv     bridge velocity and microphone pressure of the A3 mf render
  piano/a3-mf.log            the render log (its bare-board FE mode list)
  piano/a3-mf-close.wav      optional: A3 mf heard by two receivers 0.23 m above the board,
                             near the Iowa microphone positions
  piano/a3-{mf,ff}-waves.wav optional: A3 with longitudinal string waves (--string-stretching)
Writes PUBLIC_DIR/data/diagnosis.json.
"""
from __future__ import annotations
import json, math, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import verisim as V  # noqa: E402  (shared, identical measurement code)

SR = V.SR
F0, B0 = 220.25, 2.45e-4  # A3 as rendered (verisim.py measures both)


def r(v, d=4):
    return V.r(v, d)


def fk(k, f0=F0, b=B0):
    return k * f0 * math.sqrt(1 + b * k * k)


def piano(name, dur=3.0, channel=0):
    return V.segment(V.highpass(V.load(V.CAMPAIGN / "piano" / f"{name}.wav", channel)), dur)


def real(name, dur=3.0):
    return V.segment(V.highpass(V.load(V.REFS / f"{name}.wav")), dur)


def level_at(f, s, target, rel=0.012):
    sel = (f > target * (1 - rel)) & (f < target * (1 + rel))
    return s[sel].max()


def partial_levels(x, t0, t1, count=16, f0=F0, b=B0):
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


def ratio_db(x, lo, hi, t0, t1):
    """Strongest component in [lo, hi] relative to partial 1."""
    f, s = V.spectrum(x, t0, t1, 1 << 21)
    sel = (f > lo) & (f < hi)
    return 20 * math.log10(s[sel].max() / level_at(f, s, fk(1)))


def two_stage(x, kmax=4):
    """Early (0.1-0.8 s) and late (1.5-2.9 s) decay slope per partial: a real piano's
    'prompt sound' then 'aftersound' shows as a steep early slope and a flat late one."""
    n, hop = 4096, 2400
    fr = np.array([x[i:i + n] * np.hanning(n) for i in range(int(0.1 * SR), len(x) - n, hop)])
    S = np.abs(np.fft.rfft(fr, axis=1))
    fb = np.fft.rfftfreq(n, 1 / SR)
    t = 0.1 + (np.arange(len(fr)) * hop + n / 2) / SR
    out = []
    for k in range(1, kmax + 1):
        b = int(np.argmin(abs(fb - fk(k))))
        e = 20 * np.log10(S[:, b - 2:b + 3].max(axis=1) + 1e-12)
        out.append({"n": k, "early_db_s": r(np.polyfit(t[t < 0.8], e[t < 0.8], 1)[0], 3),
                    "late_db_s": r(np.polyfit(t[t > 1.5], e[t > 1.5], 1)[0], 3)})
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


def bridge_split():
    d = np.loadtxt(V.CAMPAIGN / "piano" / "a3-mf-bridge.csv", delimiter=",", skiprows=1)
    v, p = d[:, 2], d[:, 4]
    early = {"bridge_db": r(ratio_db(v, 128, 135, 0.05, 0.5), 3), "pressure_db": r(ratio_db(p, 128, 135, 0.05, 0.5), 3)}
    lv, lp = partial_levels(v, 0.0, 0.15), partial_levels(p, 0.0, 0.15)
    return early, lv, lp


def main():
    out = {"inputs": {}, "rows": []}
    board = bare_board_modes()
    rows = out["rows"]

    # 1. The early low resonance
    knock_hz = {n: r(precise_peak(piano(n, 1.6), 125, 137, 0.05, 1.5), 6) for n in ["a3-pp", "a3-mf", "a3-ff", "a3-fff", "c6-mf"]
                if (V.CAMPAIGN / "piano" / f"{n}.wav").exists()}
    early, lv, lp = bridge_split()
    real_mf = real("Piano.mf.A3")
    ours_mf = piano("a3-mf")
    real_lv = partial_levels(real_mf, 0.0, 0.15, f0=own_fit(real_mf)[0], b=own_fit(real_mf)[1])
    ours_lv = partial_levels(ours_mf, 0.0, 0.15)
    ours_knock = ratio_db(ours_mf, 125, 137, 0.05, 0.5)
    real_knock = ratio_db(real_mf, 80, 0.77 * 220, 0.05, 0.5)
    spread = max(knock_hz.values()) - min(knock_hz.values())
    rows.append({
        "id": "low-resonance",
        "symptom": f"An early resonance near {np.mean(list(knock_hz.values())):.1f} Hz sits at {early['pressure_db']:.1f} dB re the fundamental; the real piano's strongest component in that band is at {real_knock:.1f} dB.",
        "evidence": f"Same frequency within {spread:.2f} Hz for A3 at every stroke speed and for C6, so it belongs to the modelled instrument, not a string. It is {abs(knock_hz['a3-mf'] - 130.813):.2f} Hz from C3, so not a sympathetic C string, and no bare-board mode lies within {min(abs(b - knock_hz['a3-mf']) for b in board):.0f} Hz. At A3's bridge it is already {early['bridge_db']:.1f} dB re partial 1; radiation to the microphone adds {early['pressure_db'] - early['bridge_db']:.1f} dB.",
        "suspect": "Mostly mechanical: a string-loaded board resonance driven too strongly or damped too lightly. The infinite-baffle radiation adds only the last few dB.",
        "test": "Render with air reaction on the board (--rayleigh-loading). Today the production board cannot admit it: 128 modes exceed the work ceiling, the mesh is too coarse above 1.1 kHz, and the 1 kHz board's passive fit misses its accuracy limits.",
        "status": "open",
        "data": {"hz": knock_hz, "bridge_db": early["bridge_db"], "pressure_db": early["pressure_db"], "real_db": r(real_knock, 3)},
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
        "symptom": f"Between 0.2 and 1 s, {np_ours['mf'] - np_real['mf']:.0f} dB more of FrankenSim's energy lies away from the string partials than the real piano's (mf: {np_ours['mf']:.1f} vs {np_real['mf']:.1f} dB re total).",
        "evidence": "Strongest off-partial peaks: " + ", ".join(f"{p['hz']:.0f} Hz ({p['db']:.0f} dB)" for p in peaks[:6]) + "." + (f" Those within 15 Hz of partial{'s' if len(sat) > 1 else ''} {', '.join(map(str, sat))} are satellites: a string partial split by coupling to a nearby board mode." if sat else "") + (f" {len(lows)} lie below the fundamental, where only the board resonates." if lows else ""),
        "suspect": "String-board coupling with board modes that ring too long: too little damping in the wood, no air reaction, or no losses at the rim and case.",
        "test": "Raise the board's modal damping, or add air reaction, and re-measure the off-partial energy and the satellites' decay.",
        "status": "open",
        "data": {"real_db": np_real, "ours_db": np_ours, "peaks": peaks},
    })

    # 3. Excitation pattern: string side vs radiation side (bridge trace and its own pressure column)
    coloration = [r(a - b, 3) for a, b in zip(lp, lv)]
    rows.append({
        "id": "partial-pattern",
        "symptom": "In the first 150 ms FrankenSim's partials 3-6 are " + ", ".join(f"{ours_lv[k]:.0f}" for k in range(2, 6)) + " dB re partial 1; the real piano's are " + ", ".join(f"{real_lv[k]:.0f}" for k in range(2, 6)) + " dB.",
        "evidence": "At the bridge, partials 4 and 6 are already " + f"{lv[3]:.0f} and {lv[5]:.0f} dB (string and hammer side). From bridge to microphone the levels move by " + ", ".join(f"{c:+.0f}" for c in coloration[:8]) + " dB for partials 1-8: a single receiver point colours each partial by up to ±10 dB.",
        "suspect": "Two parts: the hammer/strike pattern on the string, and the single receiver point. The Iowa microphones sat 8 inches above the strings; ours is 1 m above the board.",
        "test": "Render receivers near the Iowa positions; compare power averaged over several receivers rather than one point.",
        "status": "open",
        "data": {"real": real_lv, "ours": ours_lv, "ours_bridge_trace": {"bridge": lv, "pressure": lp}, "radiation_coloration_db": coloration},
    })

    # 4. Decay shape
    ts_r, ts_o = two_stage(real_mf), two_stage(ours_mf)
    def stage2(ts):  # steep early fall, then a hold within 3 dB/s
        return [q["n"] for q in ts if q["early_db_s"] < -15 and abs(q["late_db_s"]) < 3]
    s2r, s2o = stage2(ts_r), stage2(ts_o)
    names = lambda ns: (f"partial {ns[0]}" if len(ns) == 1 else "partials " + ", ".join(map(str, ns[:-1])) + f" and {ns[-1]}") if ns else "no partial"
    rows.append({
        "id": "decay-shape",
        "symptom": f"Two-stage decay (a steep fall, then a hold: 'prompt sound' and 'aftersound') shows on {names(s2r)} of the real piano and on {names(s2o)} of FrankenSim.",
        "evidence": "Early (0.1-0.8 s) -> late (1.5-2.9 s) slope, dB/s: real " + "; ".join(f"p{q['n']} {q['early_db_s']:.0f} -> {q['late_db_s']:.0f}" for q in ts_r) + "; FrankenSim " + "; ".join(f"p{q['n']} {q['early_db_s']:.0f} -> {q['late_db_s']:.0f}" for q in ts_o) + ".",
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
    row = {
        "id": "ff-brightness",
        "symptom": f"From mf to ff the real piano gains {gain_real:.1f} dB of 2-8 kHz brightness; FrankenSim gains {gain_ours:.1f} dB, and stays dull even when struck hard enough to match the real ff's loudness.",
        "evidence": f"Sum tones f_m + f_n of a tension-modulated string (phantom partials) stand at least 10 dB above the local floor at {sum(q['present'] for q in ph_real)} of {len(ph_real)} resolvable positions (partials {ph_real[0]['k']}-{ph_real[-1]['k']}) in the real ff, as close as {max(q['re_partial_db'] for q in ph_real):.0f} dB to the true partial; in FrankenSim's ff at {sum(q['present'] for q in ph_ours)} of {len(ph_ours)}, the highest at {max(q['re_partial_db'] for q in ph_ours):.0f} dB. Each piano's positions come from its own inharmonicity fit.",
        "suspect": "No geometric string nonlinearity: a hard strike stretches the string, modulating its tension and driving longitudinal waves. For the A3 wire the first longitudinal mode is about 3.3 kHz, inside the missing band.",
        "test": "Render A3 mf and ff with --string-stretching waves (EA from the scale's own EI with E = 200 GPa) and compare brightness gain and phantom levels.",
        "status": "open",
        "data": {"phantoms_real_ff": ph_real, "phantoms_ours_ff": ph_ours},
    }
    rows.append(row)

    # 6. Cymbal
    n, hop = 1024, 240

    def band_peaks(x):
        fr = np.array([x[i:i + n] * np.hanning(n) for i in range(0, int(0.6 * SR), hop)])
        S = np.abs(np.fft.rfft(fr, axis=1)) ** 2
        fb = np.fft.rfftfreq(n, 1 / SR)
        t = (np.arange(len(fr)) * hop + n / 2) / SR
        pk = []
        for fc in V.OCT:
            sel = (fb >= fc / math.sqrt(2)) & (fb < fc * math.sqrt(2))
            e = np.convolve(S[:, sel].sum(axis=1), np.ones(5) / 5, "same")
            pk.append(r(t[int(np.argmax(e))] * 1000, 3))
        f, s = V.spectrum(x, 0.0, 1.0, 1 << 17)
        p = s ** 2
        return pk, r(10 * math.log10(p[f > 1200].sum() / p[f > 40].sum()), 3)

    rs = V.segment(V.highpass(V.load(V.REFS / "splash.stick.normal.mf.wav")), 2.5)
    os_ = V.segment(V.highpass(V.ours_clip("perc", "cym-ff")), 2.5)
    (pr, er), (po, eo) = band_peaks(rs), band_peaks(os_)
    rows.append({
        "id": "cymbal-shimmer",
        "symptom": f"The real splash keeps {100 * 10 ** (er / 10):.0f}% of its energy above 1.2 kHz; FrankenSim's keeps {100 * 10 ** (eo / 10):.0f}%.",
        "evidence": f"Real 4 kHz and 8 kHz bands peak {pr[6]:.0f} and {pr[7]:.0f} ms after the strike: energy keeps flowing upward after the stick has left. In FrankenSim every band peaks at the strike ({po[6]:.0f} ms).",
        "suspect": "Two missing pieces: shell modes above 1.2 kHz are not in the model, and there is no geometric nonlinearity to pump energy from low modes into high ones (the cascade that makes a cymbal shimmer).",
        "test": "Raise the shell mode ceiling to about 10 kHz and add von Kármán coupling between shell modes; the delayed high-band peak is the falsifier.",
        "status": "open",
        "data": {"bands_hz": V.OCT, "peak_ms_real": pr, "peak_ms_ours": po, "above_1200_db_real": er, "above_1200_db_ours": eo},
    })

    # Ablations, when their renders exist
    abl = {}
    close = V.CAMPAIGN / "piano" / "a3-mf-close.wav"
    if close.exists():
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
            row["suspect"] = (f"Upstream of the microphone. In the bridge's own motion partials 4 and 6 already sit at {lv[3]:.0f} and {lv[5]:.0f} dB, "
                              + (f"and a strike at 0.12 of the string length does not null partial{'s' if len(nodal) > 1 else ''} {' and '.join(map(str, nodal))}. " if nodal else "")
                              + "The remaining candidates are the board's mobility at A3's bridge and the felt's contact spectrum.")
            row["test"] = "Compute the bridge driving-point admittance at A3 (grand_piano's admittance mode) and look for dips near 882 and 1323 Hz; then vary the felt law alone."
        row["result"] = (f"Receivers 0.23 m above the board near the Iowa positions: the partial-pattern distance to the recording is "
                         f"{res['bass-side receiver']['pattern_gap_db']:.1f} dB (bass side) and {res['treble-side receiver']['pattern_gap_db']:.1f} dB (treble side), "
                         f"against {base['pattern_gap_db']:.1f} dB at the default point 1 m up. "
                         + ("Listening position explains a real part of the gap." if moved >= 3 else "Moving the receiver does not close the gap, so the pattern is made upstream, in the hammer, string or board.")
                         )
        near = min(res["bass-side receiver"]["low_resonance_db"], res["treble-side receiver"]["low_resonance_db"]), max(res["bass-side receiver"]["low_resonance_db"], res["treble-side receiver"]["low_resonance_db"])
        lr = next(q for q in rows if q["id"] == "low-resonance")
        lr["result"] = (f"Heard by receivers near the Iowa positions, the resonance reads {near[0]:.1f} to {near[1]:.1f} dB instead of {base['low_resonance_db']:.1f} dB: "
                        f"about {base['low_resonance_db'] - near[1]:.0f} dB of the excess came from the default listening point, 1 m above the middle of the board, where that low mode is heard far more strongly. "
                        f"About {near[1] - lr['data']['real_db']:.0f} dB remain against the real piano's {lr['data']['real_db']:.1f} dB. The air-loading test for the remainder is still blocked.")
        lr["status"] = "testing"
        lr["suspect"] = (f"Two parts. The default listening point accounts for about {base['low_resonance_db'] - near[1]:.0f} dB (shown by the receiver experiment). "
                         f"The remaining {near[1] - lr['data']['real_db']:.0f} dB is in the board: a string-loaded resonance driven too strongly or damped too lightly, possibly for lack of air loading.")
        br_ = next(q for q in rows if q["id"] == "board-ringing")
        br_["result"] = (f"Not a listening-position effect: at receivers near the Iowa positions the off-partial share is {res['bass-side receiver']['non_partial_db']:.1f} and "
                         f"{res['treble-side receiver']['non_partial_db']:.1f} dB (default point {base['non_partial_db']:.1f} dB, real piano {np_real['mf']:.1f} dB). The board itself rings too much; damping and air loading remain the candidates.")
    waves = {d: V.CAMPAIGN / "piano" / f"a3-{d}-waves.wav" for d in ("mf", "ff")}
    if all(p.exists() for p in waves.values()):
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
        lift = sw["brightness_gain_db"]["waves"] - sw["brightness_gain_db"]["linear"]
        row = next(q for q in rows if q["id"] == "ff-brightness")
        row["status"] = "confirmed" if got >= len(sw["phantoms_ff"]["waves"]) // 2 and lift >= 2 else ("testing" if got or lift >= 1 else "refuted")
        real_present = [q for q in ph_real if q["present"]]
        w_by_k = {q["k"]: q for q in sw["phantoms_ff"]["waves"]}
        l_by_k = {q["k"]: q for q in ph_ours}
        rise = max(w_by_k[q["k"]]["re_partial_db"] - l_by_k[q["k"]]["re_partial_db"] for q in real_present)
        short = [q["re_partial_db"] - w_by_k[q["k"]]["re_partial_db"] for q in real_present]
        extra = (f" Phantom levels rise by up to {rise:.0f} dB but stay {min(short):.0f}-{max(short):.0f} dB below the real piano's. "
                          "So the mechanism is right and its strength is not. In this model the longitudinal waves have fixed ends, so their force never reaches the bridge "
                          "(STRING_STRETCHING.md: no longitudinal force port to the soundboard). The phantom-partial literature (Conklin 1999; Bank and Sujbert 2005) "
                          "points to that axial force on the bridge as the main way phantom partials reach the board.")
        if row["status"] == "testing":
            row["suspect"] = ("Geometric string nonlinearity, confirmed in direction: switching on longitudinal waves creates phantom partials and adds brightness at ff. "
                              "Missing in strength: the axial string force's path into the bridge.")
            row["test"] = ("Add a longitudinal bridge port (the axial string force acting through the string's angle at the bridge, with its work in the energy ledger) "
                           "and check that phantom levels approach the real piano's at partials 13-16 (bead frankensim-na4ju).")
        if row["status"] == "refuted":
            row["suspect"] = "Not the longitudinal waves alone: with them switched on the brightness and phantom partials barely move. The hammer felt's stiffening law and the board's high-frequency response remain."
        hard = V.CAMPAIGN / "piano" / "a3-fff-waves.wav"
        if hard.exists() and (V.CAMPAIGN / "piano" / "a3-fff.wav").exists():
            w8, l8 = piano("a3-fff-waves"), piano("a3-fff")
            ph8 = {q["k"]: q for q in phantoms(w8)}
            short8 = [q["re_partial_db"] - ph8[q["k"]]["re_partial_db"] for q in real_present]
            sw["stroke_8ms"] = {"phantoms_waves": list(ph8.values()),
                                "brightness_db": {"real_ff": r(V.band_ratio(real_ff), 3), "linear": r(V.band_ratio(l8), 3), "waves": r(V.band_ratio(w8), 3)},
                                "phantom_shortfall_db": [r(x, 3) for x in short8]}
            got8 = sum(q["present"] for q in ph8.values())
            bd = sw["stroke_8ms"]["brightness_db"]
            closed, gap = bd["waves"] - bd["linear"], bd["real_ff"] - bd["linear"]
            extra += (f" Struck at 8 m/s, close to the real ff's loudness, the waves render reaches {bd['waves']:.1f} dB of 2-8 kHz brightness, against "
                      f"{bd['linear']:.1f} dB with linear strings and {bd['real_ff']:.1f} dB for the real ff: string stretching closes {closed:.1f} of that {gap:.1f} dB gap. "
                      f"Its phantom partials stay {min(short8):.0f}-{max(short8):.0f} dB below the real piano's (at 4.5 m/s: {min(short):.0f}-{max(short):.0f} dB), "
                      + ("so a harder stroke closes most of that shortfall." if max(short8) < 10 else "so a harder stroke does not close that shortfall; the missing axial force path at the bridge remains the lead."))
        row["result"] = (f"A3 rendered with longitudinal string waves (axial stiffness from the scale's own bending stiffness, steel E = 200 GPa): "
                         f"phantom partials stand above the floor at {got} of {len(sw['phantoms_ff']['waves'])} positions at ff (linear strings: {sum(q['present'] for q in ph_ours)}). "
                         f"The mf-to-ff brightness gain moves from {sw['brightness_gain_db']['linear']:.1f} dB to {sw['brightness_gain_db']['waves']:.1f} dB (real piano {sw['brightness_gain_db']['real']:.1f} dB), "
                         f"and the mf partial pattern changes by {sw['mf_unchanged_db']:.1f} dB RMS.") + extra
    out["ablations"] = abl
    out["inputs"] = {"board_modes_hz": board[:12], "a3_f0": F0, "a3_B": B0}
    (V.OUT / "data" / "diagnosis.json").write_text(json.dumps(out, separators=(",", ":")))
    for row in rows:
        print(f"\n[{row['id']}] {row['symptom']}\n  evidence: {row['evidence']}")
    print("\nablations:", json.dumps(abl, indent=1)[:3000])


if __name__ == "__main__":
    V.configure(*sys.argv[1:4])
    main()
