"""identity: the Q2 ("recognisable") statistic as the beads define it today, and the candidate redefinitions put to the owner.

UNREGISTERED. None of the candidates below is a gate until the owner chooses one on bead pexfp.357. Until then no render may be
graded with them: gate_run.py grades renders with the CURRENT definition only and reports the candidates' numbers on library
recordings alone (identity_eval.py). None of this may be used to choose, tune or fit physics (owner, 2026-10-21).

The current definition, bead_q2(): m = d_o - d_t on the spectral distance spec_db (instrgrade.distance: third-octave envelope and
harmonic levels), d_t to the ONE reference recording the pair is graded against, d_o to the nearest note of another instrument
of the same kind within 2 semitones, any dynamic; a note passes at m >= BEAD_MARGIN_DB, a family at >= 80 % of its notes. In
WhiteGull's audit (2026-10-22) it passed real recordings from a second library on 1 of 46 and the accepted piano on 0 of 2: the
distance to one recording measures its room, microphone, player and specimen. And it passes a recording frozen (freeze.py),
whose spectrum is the recording's.

The candidates judge a note c of instrument I against the SET S of the library's recordings of I with the same articulation and
kind within SEMITONES of c's note: every dynamic, take and string (at least MIN_SET of them, otherwise no claim). On a measure m:

  nearest_m(c) = min over e in S of gap_m(c, e)                                 how far c is from the nearest real recording
  ruler_m(S)   = max over e in S of min over e' in S, e' != e, of gap_m(e, e')    how far a real recording of S can sit from
                                                                                its own nearest neighbour in S
  margin_m     = ruler_m - nearest_m                                            judge(): c passes when every margin_m >= 0

  A  set membership    judge() on SET_MEASURES (spec_db, cepstral_env_db, flux_level_db): the note sits among its instrument's
                       recordings on its spectrum and on how much it moves.
  B  life in the set   judge() on LIFE_MEASURES (flux_level_db alone): the note moves as much as its instrument's recordings do.
                       A frame-to-frame difference of log spectra, so a fixed recording chain cancels out of it; it says nothing
                       about WHICH instrument (stationary noise moves a spectrum too).
  C  life + identity   B, and nearest_set(): c's nearest recording on IDENTITY_MEASURE (the cepstral envelope with its best
                       straight line in log frequency taken out, so a smooth recording-chain tilt does not count) is one of I's.

identity_eval.py measures each on library recordings only, held out across specimens (UIowa's two violins, violas, cellos and
basses), singers (VocalSet) or files, with the five nulls of freeze.null_signals and a +20 c copy of every note.

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
from decompose import HOP, SR, Note  # noqa: E402
from instrgrade import rms_gap  # noqa: E402

IDENTITY_VERSION = "identity-1"
BEAD_MARGIN_DB = 3.0  # the current definition's pass margin
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
    raise KeyError(m)


def ruler(members, m):
    """(ruler_m, the label of the most isolated member) for a set of (label, features) pairs: the largest nearest-neighbour gap."""
    worst, who = None, None
    for i, (la, fa) in enumerate(members):
        nn = [g for j, (_, fb) in enumerate(members) if j != i and (g := gap(m, fa, fb)) is not None]
        if nn and (worst is None or min(nn) > worst):
            worst, who = min(nn), la
    return worst, who


def bead_q2(cand, ref, others, margin_db=BEAD_MARGIN_DB):
    """The CURRENT definition (module docstring): d_t = spec_db to the pair's one reference `ref` (features), d_o = spec_db to the
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
