"""Join the phrase renders, draw their spectrograms and measure every note's pitch (run on trj: numpy + PIL).

usage: python3 phrasekit.py JOBS_DIR PLAN.json EXTRA.json OUT_DIR
  JOBS_DIR   the per-note and native renders (jobs2/)
  PLAN.json  phrases.py's plan
  EXTRA.json [{id, family, kind, title, audio, how, notes?: [[at, dur, note, wav?]]}]: clips that already exist
             (the piano and kit phrases of frankensim.org/music, the playable notes cut into a row)
OUT_DIR/<id>.wav (16-bit, peak -1 dBFS), OUT_DIR/<id>.webp (log-frequency spectrogram), OUT_DIR/phrases.json.
"""
import json
import math
import subprocess
import sys
import wave
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from pitchlib import f0_samples, midi_of, read_wav, resample as pitchlib_resample  # noqa: E402

JOBS = PLAN = EXTRA = OUT = None  # set by main(); analyze.py imports spectrogram() without them
SR = 48000
PX_PER_S, HEIGHT, RANGE_DB = 110, 380, 80.0
VIBRATO_SETTLED_S = 0.26  # a sung note's onset is over: the unmodulated vowels hold one pitch from here to their release
# magma, 9 anchors
MAGMA = np.array([[0, 0, 4], [28, 16, 68], [79, 18, 123], [129, 37, 129], [181, 54, 122], [229, 80, 100],
                  [251, 135, 97], [254, 194, 135], [252, 253, 191]], float)


def hz(note):
    m = midi_of(note)
    return None if m is None else 440.0 * 2 ** ((m - 69) / 12)


def resample(x, sr):
    return pitchlib_resample(x, sr, SR)


def decode(path):
    raw = subprocess.run(["ffmpeg", "-nostdin", "-v", "error", "-i", str(path), "-f", "f32le", "-ac", "1", "-ar", str(SR), "-"],
                         check=True, capture_output=True, stdin=subprocess.DEVNULL).stdout
    return np.frombuffer(raw, "<f4").astype(float)


def load(rel):
    x, sr = read_wav(JOBS / rel)
    return resample(x, sr), sr


def spectrogram(x, f_lo, f_hi):
    n, hop, nfft = 2048, SR // 200, 8192
    win = np.hanning(n)
    pad = np.concatenate([np.zeros(n // 2), x, np.zeros(n)])
    frames = 1 + (len(pad) - n) // hop
    idx = np.arange(n)[None, :] + hop * np.arange(frames)[:, None]
    mag = np.abs(np.fft.rfft(pad[idx] * win, nfft, axis=1))
    db = 20 * np.log10(mag + 1e-12)
    db -= db.max()
    freqs = np.fft.rfftfreq(nfft, 1 / SR)
    rows = f_lo * (f_hi / f_lo) ** (1 - np.arange(HEIGHT) / (HEIGHT - 1))
    img = np.stack([np.interp(rows, freqs, fr) for fr in db], axis=1)  # HEIGHT x frames
    width = max(600, min(2600, round(len(x) / SR * PX_PER_S)))
    cols = np.linspace(0, frames - 1, width)
    img = np.stack([np.interp(cols, np.arange(frames), r) for r in img])
    v = np.clip((img + RANGE_DB) / RANGE_DB, 0, 1) * (len(MAGMA) - 1)
    lo = np.floor(v).astype(int).clip(0, len(MAGMA) - 2)
    frac = (v - lo)[..., None]
    rgb = MAGMA[lo] * (1 - frac) + MAGMA[lo + 1] * frac
    return Image.fromarray(rgb.astype(np.uint8), "RGB")


def pitch(x, a, b, want):
    """Pitch of x[a:b] seconds; cents against the wanted note (None for chords or no target)."""
    seg = x[int(a * SR):int(b * SR)]
    if len(seg) < 0.05 * SR:
        return None
    f, clarity = f0_samples(seg, SR)
    if not f:
        return None
    # A fundamental must have odd harmonics. When f itself sits 20 dB under 2f and the levels at f, 3f, 5f stand 12 dB
    # or more below those at 2f, 4f, 6f, the period the autocorrelation chose is twice the true one (it happened on
    # bowed-bass notes whose 2f carries the energy): halve it, and say so.
    n = 1 << 17
    mag = np.abs(np.fft.rfft(seg[:n] * np.hanning(min(len(seg), n)), n))
    db = 20 * np.log10(mag / (mag.max() or 1) + 1e-12)
    freqs = np.fft.rfftfreq(n, 1 / SR)

    def level(fx):
        sel = (freqs >= fx * 0.98) & (freqs <= fx * 1.02)
        return float(db[sel].max()) if sel.any() else -120.0

    f_ac = f
    folded = 0
    while 6 * f < SR / 2 and folded < 3:
        odd = np.mean([level(k * f) for k in (1, 3, 5)])
        even = np.mean([level(k * f) for k in (2, 4, 6)])
        if odd >= even - 12 or level(f) >= level(2 * f) - 20:
            break
        f *= 2
        folded += 1

    # The converse: vibrato spreads the k-th harmonic k times as wide as the fundamental, so the autocorrelation at the
    # true period can sink under the candidates' bar while a strong upper harmonic stays sharp (a sung /a/ C3 at 131
    # Hz, its 7th harmonic 11 dB over the 1st, read 455 Hz and was folded up to that 7th harmonic). When every slot
    # j f/k (j = 1..k-1) carries a line within 25 dB of the level at f and 10 dB or more above the spectrum midway to
    # both neighbouring slots, f is the k-th harmonic of f/k: divide by the largest such k, and say so. A true pitch has
    # nothing between its own harmonics, so no k passes. After an octave fold only odd k are tried: an even one would
    # put back the period the fold found to have no odd harmonics.
    def floor(fx):
        sel = (freqs >= fx * 0.98) & (freqs <= fx * 1.02)
        return float(np.median(db[sel])) if sel.any() else -120.0

    divided, f_fold = 1, f
    ref = level(f)
    for k in range(8, 1, -1):
        g = f / k
        if (k % 2 or not folded) and g >= 30 and all(
                level(j * g) >= ref - 25 and level(j * g) - max(floor((j - 0.5) * g), floor((j + 0.5) * g)) >= 10
                for j in range(1, k)):
            divided = k
            break
    f /= divided
    t = hz(want) if want else None
    out = dict(f0_hz=round(f, 2), cents=None if not t else round(1200 * math.log2(f / t), 1), periodicity=round(clarity, 3))
    if divided > 1:
        moved = f", which the octave check moved to {f_fold:.1f} Hz" if folded else ""
        out["harmonic_check"] = (f"the autocorrelation read {f_ac:.1f} Hz{moved}: harmonic {divided} of {f:.1f} Hz, "
                                 f"whose every lower harmonic carries a line; corrected")
    elif folded:
        out["octave_check"] = f"the autocorrelation's period was {2 ** folded}x too long (no odd harmonics); corrected"
    return out


def whole_cycles(x, a, b, rate, m, want):
    """m, a pitch() row, read again as the mean frequency of x[a:b] seconds over the whole cycles of a vibrato at `rate` Hz.

    One autocorrelation across a swing does not read its centre: it favours the loud moments, and sung notes are louder
    in the low half of the swing (/e/ C4 read -75 cents with its swing centred near -5). So count the fundamental's
    cycles instead: its phase advance over a whole number of vibrato cycles weighs every moment alike, whatever the
    level does. m's f0_hz only has to be within a tenth of the note, to find the band. This is a mean of frequency, not
    of cents: at +-90 cents it stands about 1.2 cents above the mean of cents. On the ten unmodulated vowels it agrees
    with pitch() within 0.15 cents."""
    cycles = int((b - a) * rate)
    if cycles < 1:
        return m
    f, n = m["f0_hz"], int(round(cycles / rate * SR))
    i = int((a + b) / 2 * SR) - n // 2
    lo = max(0, i - SR // 10)  # a margin either side keeps the band filter's edges off the span
    seg = x[lo:i + n + SR // 10 + 1]
    off = np.abs(np.fft.fftfreq(len(seg), 1 / SR) - f) / f
    gain = np.where(off <= 0.2, 1.0, np.where(off >= 0.3, 0.0, 0.5 + 0.5 * np.cos(np.pi * (off - 0.2) / 0.1)))
    phase = np.unwrap(np.angle(np.fft.ifft(np.fft.fft(seg - seg.mean()) * gain)))
    f = float(phase[i - lo + n] - phase[i - lo]) / (2 * math.pi) * SR / n
    t = hz(want) if want else None
    # periodicity too: one span across a swing has no single period to repeat, and /e/ C4 read as no steady pitch for it.
    # Windows this short read lower on any tone: the unmodulated vowels give 0.81 to 0.90 here and 0.99 on the long span.
    w = SR // 25
    clear = np.median([f0_samples(x[j:j + w], SR)[1] for j in range(i, i + n - w, w)])
    return dict(m, f0_hz=round(f, 2), cents=None if not t else round(1200 * math.log2(f / t), 1),
                periodicity=round(float(clear), 3), whole_cycles=cycles)


def fade_out(y, at, secs):
    """Silence y from `at` on, over `secs` (a stopped string, a tongued note's end)."""
    i = int(round(at * SR))
    if i >= len(y):
        return
    f = int(secs * SR)
    seg = y[i:i + f]
    y[i:i + len(seg)] = seg * np.linspace(1, 0, f)[:len(seg)]
    y[i + f:] = 0


def place(out, x, at):
    i = int(round(at * SR))
    end = i + len(x)
    if end > len(out):
        out = np.concatenate([out, np.zeros(end - len(out))])
    out[i:end] += x
    return out


def finish(pid, fam, kind, title, how, x, notes, f_hi, extra=None):
    peak = np.max(np.abs(x)) or 1.0
    y = x / peak * 10 ** (-1 / 20)
    with wave.open(str(OUT / f"{pid}.wav"), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes((y * 32767).astype("<i2").tobytes())
    f_lo = 30.0 if fam in ("double-bass", "kit", "piano") else 50.0
    spectrogram(x, f_lo, f_hi).save(OUT / f"{pid}.webp", quality=88, method=6)
    doc = dict(id=pid, family=fam, kind=kind, title=title, how=how, audio=f"{pid}.wav", map=f"{pid}.webp",
               duration_s=round(len(x) / SR, 3), f_lo_hz=f_lo, f_hi_hz=f_hi, range_db=RANGE_DB, notes=notes)
    if extra:
        doc.update(extra)
    return doc


HOW = {
    "detache": "Each note is its own render: the open-string description with the string shortened by 2^(-n/12) "
               "(a rigid stop, no finger on the string) and the bow lifted at the end of the note. The renders are "
               "placed at their note times; a note on the same string stops the one before (40 ms fade), a note on "
               "another string lets it ring.",
    "tongued": "Each note is its own render of that fingering's committed description, only its duration changed; "
               "the renders are placed at their note times with a 30 ms fade where one note meets the next.",
    "cut": "Every note this family can play today, each from its own render (the committed description of that "
           "fingering, string or setting, unchanged), cut to its slot with a 40 ms fade. The pitch under each note "
           "is measured on the render.",
    "legato": "Each vowel is its own render (tract shape and fold settings are fixed within a render); consecutive "
              "vowels are joined with 80 ms crossfades, the later vowel entering after its own onset.",
}


def main():
    global JOBS, PLAN, EXTRA, OUT
    JOBS, PLAN, EXTRA, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
    OUT.mkdir(parents=True, exist_ok=True)
    docs = []
    for pid, p in json.loads(PLAN.read_text()).items():
        docs.append(render_planned(pid, p))
    for e in json.loads(EXTRA.read_text()) if EXTRA.exists() else []:
        docs.append(render_extra(e))
    (OUT / "phrases.json").write_text(json.dumps(dict(phrases=docs, pitch_method=(
        "Candidate periods from the autocorrelation of the note's steady part, checked against its spectrum: a period "
        "with no line at its own frequency whose lines all fall at multiples of k is k times too long and is moved; "
        "of the rest, the lowest with 70 % of the best harmonic contrast wins. A pitch f whose every lower slot j f/k "
        "(j < k) also carries a line, as when vibrato blurs the true period and an upper harmonic's own period wins, "
        "is the k-th harmonic of f/k and is moved down, and the note says so. A note whose autocorrelation peak "
        "stays under 0.5 has no steady pitch.")), indent=1))
    print(len(docs), "phrases")


def render_planned(pid, p):
    fam, kind = p["family"], p["kind"]
    sid = pid.replace("/", "__")
    f_hi = 12000.0
    if "native" in p:  # one lip render; notes from its schedule
        x, sr = load(p["native"])
        f_hi = min(f_hi, 0.45 * sr)
        notes, t = [], 0.25
        for note, d in p["seq"]:
            m = pitch(x, t + min(0.3, d * 0.45), t + d - 0.06, note)
            notes.append(dict(at=round(t, 3), dur=d, note=note, target_hz=round(hz(note), 2), **(m or {})))
            t += d
        how = ("One render: the player's mouth pressure and lip pre-stress change during it, so every note "
               "transition happens inside the simulation. " + ("Each note is tongued: the pressure drops to zero "
               "for 30 ms, then attacks at 14 kPa and settles to that note's sustain pressure."
               if kind != "slur" else "The pressure stays at 14 kPa; only the lip tightens and relaxes, "
               "so the lip and the bore find each note themselves."))
        return finish(sid, fam, kind, p["title"], how, x, notes, f_hi, dict(native=True))
    out = np.zeros(1)  # legato: one track, each vowel crossfaded into the last
    tracks = []  # otherwise one track per note, silenced when a later note takes its string (or, for winds, any note)
    notes = []
    for k, n in enumerate(p["notes"]):
        if n.get("refused"):  # the solver refused this note: its slot stays silent (its string is still stopped) and
            # the note row carries the refusal, so the page counts it as asked for and not played
            strings = set(n["string"].split("+")) if p["join"] == "detache" and n.get("string") else {None}
            for held, y in tracks:
                if held & strings:
                    fade_out(y, n["at"], 0.04 if p["join"] == "detache" else 0.03)
            row = dict(at=n["at"], dur=n["dur"], note=n["note"], refused=n["refused"])
            if "+" not in n["note"] and hz(n["note"]):
                row["target_hz"] = round(hz(n["note"]), 2)
            row.update({key: n[key] for key in ("string", "semitones") if key in n})
            notes.append(row)
            continue
        x, sr = load(n["wav"])
        if p["join"] == "cut":  # a long held render, cut to its slot with a 40 ms fade
            x = x[:int((n["dur"] + 0.04) * SR)].copy()
            x[-int(0.04 * SR):] *= np.linspace(1, 0, int(0.04 * SR))
        if p["join"] == "legato":
            if k:
                x = x[int(0.12 * SR):]
                i = int(round(n["at"] * SR))
                ramp = np.linspace(0, 1, int(0.08 * SR))
                x[:len(ramp)] *= ramp
                tail = out[i:i + len(ramp)]
                out[i:i + len(tail)] = tail * ramp[::-1][:len(tail)]
                out[i + len(ramp):] = 0
            out = place(out, x, n["at"])
        else:
            # A bowed note stops whatever its own string (or strings, for a double stop) was playing and lets the
            # other strings ring; a wind note ends the one before.
            strings = set(n["string"].split("+")) if p["join"] == "detache" and n.get("string") else {None}
            for held, y in tracks:
                if held & strings:
                    fade_out(y, n["at"], 0.04 if p["join"] == "detache" else 0.03)
            tracks.append((strings, place(np.zeros(1), x, n["at"])))
        is_chord = "+" in n["note"]
        want = None if is_chord or p["family"] == "voice" else n["note"]
        if p["family"] == "voice":
            want = n["note"].split()[-1]
        lead = 0.16 if p["join"] == "legato" else (0.12 if p["join"] == "detache" else 0.08)
        if p["join"] == "cut":  # the held render's steady middle, as the playable-notes check measured it
            full = load(n["wav"])[0]
            m = pitch(full, 0.25 * len(full) / SR, 0.75 * len(full) / SR, want)
        else:
            m = None if is_chord else pitch(x if p["join"] != "legato" else load(n["wav"])[0], lead, max(lead + 0.06, n["dur"] - 0.04), want)
            if m and p.get("vibrato_hz"):  # legato, so this reads the note's own render and its times are the render's
                m = whole_cycles(load(n["wav"])[0], VIBRATO_SETTLED_S, max(lead + 0.06, n["dur"] - 0.04), p["vibrato_hz"], m, want)
        row = dict(at=n["at"], dur=n["dur"], note=n["note"])
        if not is_chord and want and hz(want):
            row["target_hz"] = round(hz(want), 2)
        if n.get("string"):
            row["string"] = n["string"]
        if "semitones" in n:
            row["semitones"] = n["semitones"]
        if m:
            row.update(m)
        notes.append(row)
    if tracks:
        # at least as long as the last note's slot, which a refused last note would otherwise leave off the map
        slots = int(round(max(n["at"] + n["dur"] for n in p["notes"]) * SR))
        out = np.zeros(max(slots, *(len(y) for _, y in tracks)))
        for _, y in tracks:
            out[:len(y)] += y
    # a plan may say how its renders were made (another engine than the cards); otherwise the join's own text
    return finish(sid, fam, kind, p["title"], p.get("how") or HOW[p["join"]], out, notes, f_hi)


def render_extra(e):
    x = decode(e["audio"])
    notes = []
    for at, dur, note, *_ in e.get("notes", []):
        row = dict(at=at, dur=dur, note=note)
        if hz(note):
            row["target_hz"] = round(hz(note), 2)
            m = pitch(x, at + 0.15, at + dur - 0.05, note)
            if m:
                row.update(m)
        notes.append(row)
    doc = finish(e["id"], e["family"], e["kind"], e["title"], e["how"], x, notes, e.get("f_hi", 12000.0),
                 {k: v for k, v in e.items() if k in ("real_audio", "credit", "licence", "source")})
    if e.get("real_audio"):
        xr = decode(e["real_audio"])
        spectrogram(xr, doc["f_lo_hz"], doc["f_hi_hz"]).save(OUT / f"{e['id']}-real.webp", quality=88, method=6)
        peak = np.max(np.abs(xr)) or 1.0
        with wave.open(str(OUT / f"{e['id']}-real.wav"), "wb") as w:
            w.setnchannels(1)
            w.setsampwidth(2)
            w.setframerate(SR)
            w.writeframes((xr / peak * 10 ** (-1 / 20) * 32767).astype("<i2").tobytes())
        doc.update(real=dict(audio=f"{e['id']}-real.wav", map=f"{e['id']}-real.webp", duration_s=round(len(xr) / SR, 3)))
    return doc


if __name__ == "__main__":
    main()
