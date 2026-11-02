"""identity: Q2 Alive, the registered Q2 statistic, with the retired definition and the candidates it was chosen from.

REGISTERED 2026-10-22: Q2 Alive (STATISTIC), the owner's choice (a) on bead pexfp.357, relayed by MagentaBadger (mail 5833).
"Moves like a real recording of the instrument", for SUSTAINED families only (bowed strings arco, winds, brass, voice):
  q2_alive(c, S, kind) = judge(c, S, LIFE_MEASURES): the note's mean spectral flux after the attack, as a level, lies within the
  spread of its instrument's real recordings S (below). A held-out note passes at margin >= 0 dB; a family passes Q2 when >= 80 %
  of its held-out graded notes pass and gate_run.py's controls hold in the same run (real recordings from a second library pass at
  >= 80 %; the five null baselines of freeze.null_signals all fail; +20 c moves no margin by 1 dB or more).
It is NOT an identity test: a real recording passes as another instrument most of the time. Which instrument it is goes to the
owner's blind test (Q4); nearest_set() is reported beside Q2 Alive, never gating. Struck and plucked notes (piano, bars, guitar,
pizzicato) get NO claim (q2_alive returns None): after the attack a real decay is as still as a null, and the flux of a short
pizzicato is not stable under +20 c (identity_eval.py). It is a final validation and must never be used to choose, tune or fit
physics (owner, 2026-10-21).

RETIRED, bead_q2(): m = d_o - d_t on the spectral distance spec_db (instrgrade.distance: third-octave envelope and harmonic
levels), d_t to the ONE reference recording the pair is graded against, d_o to the nearest note of another instrument of the same
kind within 2 semitones, any dynamic; pass at m >= BEAD_MARGIN_DB. It passed real recordings from a second library on 1 of 46 and
the piano on 0 of 3, and the frozen null on 28 of 56 (gate_run.py, 2026-10-22): the distance to one recording measures its room,
microphone, player and specimen. Kept so that table can be reproduced.

The candidates the owner chose from judge a note c of instrument I against the SET S of the library's recordings of I with the
same articulation and kind within SEMITONES of c's note: every dynamic, take and string (at least MIN_SET of them, otherwise no
claim). On a measure m:

  nearest_m(c) = min over e in S of gap_m(c, e)                                 how far c is from the nearest real recording
  ruler_m(S)   = max over e in S of min over e' in S, e' != e, of gap_m(e, e')    how far a real recording of S can sit from
                                                                                its own nearest neighbour in S
  margin_m     = ruler_m - nearest_m                                            judge(): c passes when every margin_m >= 0

  A  set membership    judge() on SET_MEASURES (spec_db, cepstral_env_db, flux_level_db): the note sits among its instrument's
                       recordings on its spectrum and on how much it moves. Not chosen: across specimens it passes 34.5 %.
  B  life in the set   judge() on LIFE_MEASURES (flux_level_db alone): the note moves as much as its instrument's recordings do.
                       A frame-to-frame difference of log spectra, so a fixed recording chain cancels out of it; it says nothing
                       about WHICH instrument (stationary noise in the map moves a spectrum too). CHOSEN: Q2 Alive, sustained only.
  C  life + identity   B, and nearest_set(): c's nearest recording on IDENTITY_MEASURE (the cepstral envelope with its best
                       straight line in log frequency taken out, so a smooth recording-chain tilt does not count) is one of I's.
                       Not chosen as a gate (36.2 % across specimens); its identity is reported beside Q2 Alive.

identity_eval.py measures each on library recordings only, held out across specimens (UIowa's two violins, violas, cellos and
basses), singers (VocalSet) or files, with the five nulls of freeze.null_signals and a +20 c copy of every note. On sustained
notes B passed 83.5 % of them, 1 of 2155 nulls, and moved 2 of 428 notes by 1.0-1.2 dB under +20 c.

flux_level_db = |20 log10| of the ratio of two notes' mean spectral flux (decompose.flux_curve over the sounding frames of
[FLUX_SKIP_S, FLUX_SPAN_S]): a level, so a note that does not move (a periodic render, a frozen recording) is far from every real
one whatever the units.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import freeze  # noqa: E402
from decompose import HOP, SR, Note  # noqa: E402
from instrgrade import rms_gap  # noqa: E402

IDENTITY_VERSION = "identity-1"  # the features' version (a cache of them stays valid while it does not change)
STATISTIC = "q2-alive-1"  # the registered Q2 statistic (pexfp.357)
BEAD_MARGIN_DB = 3.0  # the retired definition's pass margin
SEMITONES = 2  # a candidate's set: notes within this many semitones of its own
MIN_SET = 3  # fewer real recordings than this in the set: no claim
SET_MEASURES = ("spec_db", "cepstral_env_db", "flux_level_db")  # candidate A
LIFE_MEASURES = ("flux_level_db",)  # candidates B and C
IDENTITY_MEASURE = "cep_detilt_db"  # candidate C's identity, by the nearest set
# The mean flux is taken over the sounding frames of [FLUX_SKIP_S, FLUX_SPAN_S]: after the attack. An attack is a note's shape,
# which a render and the frozen null (freeze.ATTACK_S kept as recorded, handed over in as long again) share by construction, and its
# one or two steep frames weigh with where the 10 ms frame grid falls on the onset. The kept attack and its handover end 0.17 s
# into the aligned note; a frame centred at 0.25 s (4096 points, 85 ms) and its predecessor lie after it. This does NOT make the
# flux a measure for struck and plucked notes (identity_eval.py, 2026-10-22): after the attack a real bar's decay is as still as a
# null (UIowa xylophone sets reach 0.018 dB per frame, the nulls 0.046), and a short pizzicato has a handful of sounding frames
# left, whose mean a +20 c copy moved by up to 25 dB.
FLUX_SKIP_S, FLUX_SPAN_S = 0.25, 2.0
FLUX_MIN_DB = 1e-3  # a mean flux below this (a perfectly periodic tone reads 0) is held here before its level is taken
NHARM = 20  # harmonic levels compared, as instrgrade.distance does


def features(n: Note) -> dict:
    """What the measures need from one aligned note, as plain numbers (JSON-able, so a library's features can be cached)."""
    t = np.arange(n.frames) * HOP / SR
    use = n.sound[: n.frames] & (t >= FLUX_SKIP_S) & (t < FLUX_SPAN_S)
    flux = float(np.mean(n.flux_db[: n.frames][use])) if use.any() else float("nan")
    return dict(version=IDENTITY_VERSION, f0_hz=float(n.f0), kind=n.kind, duration_s=float(n.dur), level_db=float(n.level_db),
                bands_db=[float(v) for v in n.desc["bands_db"]], harm_db=[float(v) for v in n.desc["harm_db"][:NHARM]],
                env_rel_db=[float(v) for v in n.env_db - n.level_db], mean_flux_db=flux)


def flux_level_db(f):
    return 20 * math.log10(max(f["mean_flux_db"], FLUX_MIN_DB))


# Decay-domain liveness for struck and plucked notes (bead pexfp.360; candidates, NOT registered: the owner chooses). A real decay
# moves in ways the five nulls do not: its partials' levels wobble about their own smooth decay (two polarisations, two strings or a
# bar's mode doublet beating; the room), and their frequencies move (tension modulation, the same beating). freeze's frozen null keeps
# every partial's level smoothed over SMOOTH_S and holds its frequency; the four tones carry one smoothed envelope at fixed
# frequencies. Each partial is a spectral peak of the decay (so a bar's or a stiff string's partials are found as they are), heard
# through a Hann band (hann_band, bw = min(30 Hz, 0.35 f0)) that follows the partial's slow frequency (so a glide does not bend
# its level on the band's slope) and sampled at DECAY_RATE. Only frames 30 dB over the partial's own floor and within 40 dB of its
# own peak count, so the recording's noise is not read as movement.
DECAY_VERSION = "decay-1"
DECAY_MEASURES = ("decay_amp_db", "decay_freq_db", "decay_spread_db")  # D-amp, D-freq, D-dec; compared as levels
DECAY_SKIP_S = 2 * freeze.ATTACK_S  # past the frozen null's kept attack and its handover to the smoothed envelope
DECAY_SPAN_S = 2.0  # at most this long after DECAY_SKIP_S
DECAY_RATE = 200  # Hz, the partial envelopes' sample rate (bw <= 30 Hz)
DECAY_PEAKS = 12  # partials followed: the strongest spectral peaks within DECAY_RANGE_DB of the strongest
DECAY_RANGE_DB = 40.0
DECAY_DROP_DB = 40.0  # a partial counts while within this of its own peak ...
DECAY_SNR_DB = 30.0  # ... and this far over its own floor (decay_floor_db)
DECAY_FLOOR_NARROW = 3  # the floor's band is bw / this, so the partial itself leaks into it more than 50 dB down
DECAY_MIN_S = 0.15  # a partial with fewer counted seconds than this is left out
DECAY_PHASE_S = 0.05  # D-freq: the slope of the partial's phase over this long
DECAY_FLOORS = {"decay_amp_db": 1e-4, "decay_freq_db": 1e-4, "decay_spread_db": 1e-3}  # values held here before the level is taken


def hann_band(z, bw):
    """Low-pass of a complex signal by a Hann kernel 2/bw long (first null at bw Hz), by FFT, centred. Not freeze's brick-wall:
    a brick wall rings after a struck note's attack with tails falling only as 1/t (about +-0.6 dB of level 0.16 s after the onset
    at bw = 30 Hz), which is larger than the movement measured here; this kernel's response ends 1/bw after any event."""
    half = max(1, int(SR / bw))
    w = np.hanning(2 * half + 1)
    w /= w.sum()
    n = 1 << int(math.ceil(math.log2(len(z) + 2 * half + 1)))
    return np.fft.ifft(np.fft.fft(z, n) * np.fft.fft(w, n))[half:half + len(z)]


def local_line(L, use, half):
    """At every sample, the value of the Hann-weighted (half-width `half` samples) least-squares straight line through the samples
    of L where `use` holds. A straight decay in dB is reproduced exactly, at the window's ends too; the residual L - local_line
    is the movement faster than the window (freeze keeps everything slower in its frozen null). NaN where the fit is singular."""
    d = np.arange(-half, half + 1, dtype=np.float64)
    h = np.hanning(2 * half + 1)
    m = use.astype(np.float64)
    Lm = np.where(use, L, 0.0)
    # sums over the window centred on i of h(d) d^p m(i + d) [L(i + d)]: h is even, so convolving with d^p h flips the sign of odd p
    s0, s1, s2 = (np.convolve(m, h * d ** p, mode="same") * (-1) ** p for p in (0, 1, 2))
    t0, t1 = (np.convolve(Lm, h * d ** p, mode="same") * (-1) ** p for p in (0, 1))
    det = s0 * s2 - s1 * s1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(det > 1e-9 * np.maximum(s0 * s2, 1e-30), (t0 * s2 - t1 * s1) / det, np.nan)


def decay_peaks(y, f0_hz, s0):
    """The frequencies [Hz] of the strongest spectral peaks of the decay that starts at sample s0: at most DECAY_PEAKS, within
    DECAY_RANGE_DB of the strongest, above 0.8 f0 and below freeze.F_TOP, at least 0.5 f0 apart."""
    seg = y[s0:s0 + int(0.5 * SR)]
    n = 1 << int(math.ceil(math.log2(4 * len(seg))))
    X = 20 * np.log10(np.abs(np.fft.rfft(seg * np.hanning(len(seg)), n)) + 1e-30)
    fr = np.fft.rfftfreq(n, 1 / SR)
    i = np.arange(1, len(X) - 1)
    i = i[(X[i] > X[i - 1]) & (X[i] >= X[i + 1]) & (fr[i] > 0.8 * f0_hz) & (fr[i] < freeze.F_TOP)]
    picked = []
    for j in i[np.argsort(-X[i], kind="stable")]:
        if picked and X[j] < X[picked[0]] - DECAY_RANGE_DB or len(picked) == DECAY_PEAKS:
            break
        if all(abs(fr[j] - fr[p]) >= 0.5 * f0_hz for p in picked):
            picked.append(j)
    return [float(fr[j]) for j in picked]


def decay_floor_db(y, t, fk, f0_hz, bw, step, i0, i1):
    """Partial fk's floor in its own band's units: what the band hears midway to the neighbouring harmonics (fk +- f0/2, the
    quieter side), over the decay window [i0, i1) of the envelope samples. The noise, the room and the other partials' leakage
    there, not the partial's own tail (a floor read from the end of the cut is the partial itself when a note is still ringing).
    Heard through a band DECAY_FLOOR_NARROW times narrower, then raised by 10 log10 of that ratio (a noise floor's power scales
    with the band)."""
    lv = []
    for fo in (fk - 0.5 * f0_hz, fk + 0.5 * f0_hz):
        if fo > 0:
            b = np.abs(hann_band(y * np.exp(-2j * np.pi * fo * t), bw / DECAY_FLOOR_NARROW))[::step][i0:i1]
            lv.append(float(np.median(20 * np.log10(b + 1e-30))))
    return min(lv) + 10 * math.log10(DECAY_FLOOR_NARROW)


def decay_features(y, f0_hz):
    """The decay-domain measures of one aligned note `y` (the module comment above DECAY_VERSION): {decay_amp_db, decay_freq_db,
    decay_spread_db, partials, ...}. Each is the median over the partials counted (spread: the std of their fitted decay rates);
    NaN when fewer than one (spread: three) partials have DECAY_MIN_S of counted frames."""
    y = np.asarray(y, dtype=np.float64)
    out = dict(version=DECAY_VERSION, partials=0, decay_amp_db=float("nan"), decay_freq_db=float("nan"), decay_spread_db=float("nan"))
    s0 = freeze.onset_of(y) + int(DECAY_SKIP_S * SR)
    if s0 + int(DECAY_MIN_S * SR) >= len(y):
        return out
    step = SR // DECAY_RATE
    bw = min(30.0, 0.35 * f0_hz)
    t = np.arange(len(y)) / SR
    i0, i_max = -(-s0 // step), (s0 + int(DECAY_SPAN_S * SR)) // step
    half, lag = int(freeze.SMOOTH_S * DECAY_RATE), int(DECAY_PHASE_S * DECAY_RATE)
    amps, freqs, slopes = [], [], []
    for fk in decay_peaks(y, f0_hz, s0):
        floor = decay_floor_db(y, t, fk, f0_hz, bw, step, i0, i_max)
        theta = np.zeros(len(y))
        for _ in range(2):
            # pass 1 at fk; pass 2 along the partial's slow frequency (pass 1's phase slope through the same local straight line as the
            # levels), so a gliding partial stays at the band's centre and its level is not bent by the band's slope
            a = hann_band(y * np.exp(-1j * (2 * np.pi * fk * t + theta)), bw)[::step]
            L = 20 * np.log10(np.abs(a) + 1e-30)
            thr = max(float(np.max(L[i0:])) - DECAY_DROP_DB, floor + DECAY_SNR_DB)
            use = np.zeros(len(L), dtype=bool)
            use[i0:i_max] = L[i0:i_max] >= thr
            ph = np.unwrap(np.angle(a)) + theta[::step]
            if theta.any() or use.sum() < DECAY_MIN_S * DECAY_RATE:
                break
            dv = np.full(len(a), np.nan)
            dv[lag // 2:len(a) - lag + lag // 2] = (ph[lag:] - ph[:-lag]) * DECAY_RATE / (2 * np.pi * lag)  # Hz off fk, centred
            ok = use & np.isfinite(dv)
            sl = local_line(np.where(ok, dv, 0.0), ok, half)
            k = np.flatnonzero(ok & np.isfinite(sl))
            if len(k) < 2:
                break
            theta = 2 * np.pi * np.cumsum(np.interp(np.arange(len(y)) / step, k, sl[k])) / SR
        if use.sum() < DECAY_MIN_S * DECAY_RATE:
            continue
        resid = (L - local_line(L, use, half))[use]
        resid = resid[np.isfinite(resid)]
        if len(resid) < DECAY_MIN_S * DECAY_RATE:
            continue
        amps.append(float(np.sqrt(np.mean(resid ** 2))))
        dev = np.full(len(a), np.nan)
        dev[lag:] = (ph[lag:] - ph[:-lag]) * DECAY_RATE / (2 * np.pi * lag)  # Hz off fk, over the last DECAY_PHASE_S
        ok = use & np.isfinite(dev)
        ok[lag:] &= use[:-lag]  # both ends of the phase span counted
        if ok.sum() >= 2:
            c = 1200 * np.log2(np.maximum(fk + dev[ok], 1e-9) / fk)
            freqs.append(float(np.sqrt(np.mean((c - np.median(c)) ** 2))))
        tt = np.flatnonzero(use) / DECAY_RATE
        slopes.append(float(np.polyfit(tt, L[use], 1)[0]))
    out["partials"] = len(amps)
    if amps:
        out["decay_amp_db"] = float(np.median(amps))
    if freqs:
        out["decay_freq_db"] = float(np.median(freqs))
    if len(slopes) >= 3:
        out["decay_spread_db"] = float(np.std(slopes))
    return out


def decay_level_db(f, m):
    v = f.get(m, float("nan"))
    return 20 * math.log10(max(v, DECAY_FLOORS[m])) if np.isfinite(v) else float("nan")


def gap(m, a, b):
    """The gap between notes a and b on measure m [dB], from their features(). None when a note lacks what m needs."""
    if m == "spec_db":
        env = rms_gap(a["bands_db"], b["bands_db"])
        harm = rms_gap(a["harm_db"], b["harm_db"])
        return float(math.hypot(env, harm)) if np.isfinite(harm) else float(env * math.sqrt(2))
    if m in ("cepstral_env_db", "cep_detilt_db"):
        d = np.asarray(a["env_rel_db"]) - np.asarray(b["env_rel_db"])
        if m == "cep_detilt_db":  # the best straight line in log2(f) taken out; the envelope grid is 1/6 octave from 100 Hz
            x = np.arange(len(d)) / 6.0
            A = np.stack([np.ones_like(x), x], axis=1)
            d = d - A @ np.linalg.lstsq(A, d, rcond=None)[0]
        return float(np.sqrt(np.mean(d * d)))
    if m == "flux_level_db":
        if not (np.isfinite(a["mean_flux_db"]) and np.isfinite(b["mean_flux_db"])):
            return None
        return abs(flux_level_db(a) - flux_level_db(b))
    if m in DECAY_MEASURES:
        la, lb = decay_level_db(a, m), decay_level_db(b, m)
        return abs(la - lb) if np.isfinite(la) and np.isfinite(lb) else None
    raise KeyError(m)


def ruler(members, m):
    """(ruler_m, the label of the most isolated member) for a set of (label, features) pairs: the largest nearest-neighbour gap."""
    worst, who = None, None
    for i, (la, fa) in enumerate(members):
        nn = [g for j, (_, fb) in enumerate(members) if j != i and (g := gap(m, fa, fb)) is not None]
        if nn and (worst is None or min(nn) > worst):
            worst, who = min(nn), la
    return worst, who


def q2_alive(cand, members, kind, rulers=None):
    """Q2 Alive, the registered Q2 statistic (module docstring): judge() on LIFE_MEASURES against the set `members` of the
    instrument's real recordings. A note of kind "decaying" (struck or plucked) gets no claim: recognisable None and the reason."""
    if kind != "sustained":
        return dict(statistic=STATISTIC, version=IDENTITY_VERSION, recognisable=None, margin_db=None, measures={}, set_size=len(members),
                    reason="no claim: Q2 Alive makes no realism claim for struck and plucked notes (owner, 2026-10-22)")
    return dict(judge(cand, members, LIFE_MEASURES, rulers), statistic=STATISTIC)


def bead_q2(cand, ref, others, margin_db=BEAD_MARGIN_DB):
    """The RETIRED definition (module docstring): d_t = spec_db to the pair's one reference `ref` (features), d_o = spec_db to the
    nearest of `others` [(label, features)]: the notes of other instruments of the same kind within 2 semitones, any dynamic.
    Returns {d_t_db, d_o_db, nearest_other, margin_db = d_o - d_t, recognisable: margin >= margin_db} or None without others."""
    best = min(((gap("spec_db", cand, f), la) for la, f in others), default=None)
    if best is None or ref is None:
        return None
    d_t = gap("spec_db", cand, ref)
    return dict(definition="bead", d_t_db=d_t, d_o_db=best[0], nearest_other=best[1], margin_db=best[0] - d_t,
                recognisable=bool(best[0] - d_t >= margin_db))


def judge(cand, members, measures, rulers=None):
    """Judge one candidate's features against its set on `measures` (SET_MEASURES: candidate A; LIFE_MEASURES: B and C): a list of
    (label, features) of the real recordings (see the module docstring). Returns {recognisable: True | False | None, margin_db,
    measures: {m: {nearest_db, nearest, ruler_db, isolated, margin_db}}, set_size}; recognisable is None (no claim) when the set is
    smaller than MIN_SET. `rulers` ({m: ruler(members, m)}) saves recomputing a set's rulers when many candidates meet it."""
    out = dict(version=IDENTITY_VERSION, set_size=len(members), measures={})
    if len(members) < MIN_SET:
        return dict(out, recognisable=None, margin_db=None, reason=f"set of {len(members)} real recordings, fewer than {MIN_SET}")
    for m in measures:
        gs = [(g, la) for la, f in members if (g := gap(m, cand, f)) is not None]
        r, who = rulers[m] if rulers else ruler(members, m)
        if not gs or r is None:
            out["measures"][m] = dict(nearest_db=None, ruler_db=r, margin_db=None)
            continue
        g, la = min(gs)
        out["measures"][m] = dict(nearest_db=g, nearest=la, ruler_db=r, isolated=who, margin_db=r - g)
    margins = [v["margin_db"] for v in out["measures"].values()]
    if any(v is None for v in margins):
        return dict(out, recognisable=None, margin_db=None, reason="a measure has no value for the candidate or the set")
    return dict(out, recognisable=bool(min(margins) >= 0.0), margin_db=float(min(margins)))


def nearest_set(cand, own, others, m=IDENTITY_MEASURE):
    """Candidate C's identity: is c's nearest real recording on measure m one of its own instrument's? `own` is c's set
    [(label, features)], `others` {instrument: [(label, features)]} the other instruments' sets at c's pitch. Returns {measure,
    own_db, other_db, nearest_other, margin_db = other_db - own_db, own_nearest}, or None when either side has no value."""
    d_own = min((g for _, f in own if (g := gap(m, cand, f)) is not None), default=None)
    best = min(((g, inst) for inst, ms in others.items() for _, f in ms if (g := gap(m, cand, f)) is not None), default=None)
    if d_own is None or best is None:
        return None
    return dict(measure=m, own_db=d_own, other_db=best[0], nearest_other=best[1], margin_db=best[0] - d_own,
                own_nearest=bool(d_own < best[0]))


def in_set(entry, inst, art, kind, midi):
    """Whether library entry `entry` belongs to the set of a candidate of instrument `inst`, articulation `art`, kind `kind`
    at MIDI note `midi`."""
    return entry["inst"] == inst and entry["art"] == art and entry["kind"] == kind and abs(entry["note"] - midi) <= SEMITONES
