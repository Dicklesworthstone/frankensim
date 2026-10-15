// From a MIDI file's parts to one pianist's part.
//
// A MIDI file written for a General MIDI synthesizer (most game and pop files) gives every part its
// own instrument and balances them with channel volume (CC7) and expression (CC11). Played as it
// stands on one piano, every part lands at full force on the same strings: leads doubled on three
// instruments, delayed copies standing in for echo, chiptune arpeggios of 20 ms notes, basses in
// the piano's bottom octave, and velocities that are all 127 because the synthesizer's volume did
// the shaping. This stage reads each part's role and turns the file into something one piano can
// play: it keeps the file's own volume, leaves out copies, thins ornaments, lifts basses and gives
// flat parts dynamics by role. Every decision is reported per part.

export type Role = "percussion" | "bass" | "melody" | "harmony" | "silent" | "ornament" | "pad";

export interface SourceNote {
  t: number;
  off: number;
  key: number;
  /** written velocity 1..127 */
  vel: number;
  /** the channel's volume x expression at the note's start, 0..1 (1 when the file sends none) */
  level: number;
}
export interface SourcePart {
  track: number;
  channel: number;
  name: string;
  program: number | null;
  /** the role from the file's names and statistics (M1 rules) */
  baseRole: Role;
  notes: SourceNote[];
}
export interface PlayedNote { key: number; t: number; off: number; vel: number; part: number; prio: number }
export interface PartPlan {
  role: Role;
  /** what happened to the part, in words, for the page */
  actions: string[];
  notesIn: number;
  notesOut: number;
  leftOut: boolean;
}

/** Below this effective velocity (written x volume x expression) a note is silent. */
const SILENT = 6;
/** Copies: share of a part's notes that match a kept part (same onset or a constant delay). */
const COPY_SHARE = 0.7;
/** Ornaments (chiptune arpeggios, echo trails): notes this short, this many per second. */
const ORNAMENT_DUR = 0.05;
const ORNAMENT_RATE = 6;
/** An ornament keeps at most one note per this many seconds. */
const ORNAMENT_GAP = 0.125;
/** The shortest note a key plays (a shorter one lets the damper fall on the hammer's echo). */
const MIN_DUR = 0.06;
/** Dynamics for parts whose velocities are all the same (machine velocities), by role. */
const TARGET: Partial<Record<Role, number>> = { melody: 84, bass: 78, harmony: 62, ornament: 46, pad: 42 };
const PRIO: Partial<Record<Role, number>> = { melody: 5, bass: 4, harmony: 3, pad: 2, ornament: 1 };
const ROLE_WORD: Partial<Record<Role, string>> = { melody: "melody", bass: "bass", harmony: "accompaniment", pad: "pad", ornament: "ornament" };

const median = (v: number[]): number => { const s = [...v].sort((a, b) => a - b); return s.length ? s[s.length >> 1]! : 0; };
const quantile = (v: number[], q: number): number => { const s = [...v].sort((a, b) => a - b); return s.length ? s[Math.min(s.length - 1, Math.floor(q * s.length))]! : 0; };

export function orchestrate(parts: SourcePart[]): { notes: PlayedNote[]; plans: PartPlan[]; silenced: number } {
  let silenced = 0;
  const plans: PartPlan[] = parts.map((p) => ({ role: p.baseRole, actions: [], notesIn: p.notes.length, notesOut: 0, leftOut: false }));
  // 1. the file's own volume: a note faded out by CC7/CC11 is not played
  const live: SourceNote[][] = parts.map((p, i) => {
    const kept = p.notes.filter((n) => n.vel * n.level >= SILENT);
    const quiet = p.notes.length - kept.length;
    silenced += quiet;
    if (quiet) plans[i]!.actions.push(`${quiet} notes silenced by the file's own volume`);
    return kept.sort((a, b) => a.t - b.t);
  });

  // 2. roles: the M1 role, refined by what the notes do and the GM program the file asks for
  parts.forEach((p, i) => {
    const ns = live[i]!;
    const plan = plans[i]!;
    if (!ns.length) { plan.role = "silent"; return; }
    const durs = ns.map((n) => n.off - n.t);
    const span = Math.max(0.5, ns[ns.length - 1]!.t - ns[0]!.t);
    const rate = ns.length / span;
    const prog = p.program;
    const poly = polyphony(ns);
    if (median(durs) < ORNAMENT_DUR && rate > ORNAMENT_RATE) plan.role = "ornament";
    else if (prog !== null && prog >= 32 && prog <= 39 && poly <= 2 && median(ns.map((n) => n.key)) < 55) plan.role = "bass";
    else if (prog !== null && prog >= 88 && prog <= 103 && parts.length > 2) plan.role = "pad";
    else plan.role = p.baseRole === "silent" ? "harmony" : p.baseRole;
  });

  // 3. copies: a part that doubles a kept part (same onsets, same pitch classes, within two
  // octaves) or repeats it a constant delay later (an echo) is left out
  const order = parts.map((_, i) => i).filter((i) => plans[i]!.role !== "silent")
    .sort((a, b) => (PRIO[plans[b]!.role] ?? 0) - (PRIO[plans[a]!.role] ?? 0) || live[b]!.length - live[a]!.length);
  const kept: number[] = [];
  for (const i of order) {
    let copy: string | null = null;
    for (const j of kept) {
      const r = relation(live[i]!, live[j]!);
      const label = `${parts[j]!.track + 1}·${parts[j]!.channel + 1}`;
      if (r.same >= COPY_SHARE) { copy = `doubles part ${label}`; break; }
      if (r.echo >= COPY_SHARE) { copy = `an echo of part ${label}, ${Math.round(r.delay * 1000)} ms later`; break; }
    }
    if (copy) { plans[i]!.leftOut = true; plans[i]!.actions.push(`left out: ${copy}`); }
    else kept.push(i);
  }

  // 4. per part: register, ornaments, dynamics, shortest note
  const out: PlayedNote[] = [];
  for (const i of kept) {
    const plan = plans[i]!;
    let ns = live[i]!.map((n) => ({ ...n }));
    if (plan.role === "bass") {
      let lift = 0;
      // only a bass that lives in the bottom octave (median below G1) moves: a walking bass from
      // E1 up is a pianist's left hand, but A0-A#1 on the strings is mostly rumble
      while (median(ns.map((n) => n.key)) + lift < 31 && lift < 24) lift += 12;
      if (lift) { ns.forEach((n) => (n.key += lift)); plan.actions.push(`raised ${lift / 12} octave${lift > 12 ? "s" : ""} out of the piano's muddiest register`); }
    }
    if (plan.role === "ornament" || plan.role === "pad") {
      const high = ns.filter((n) => n.key > 108).length;
      if (high) { ns = ns.filter((n) => n.key <= 108); plan.actions.push(`${high} notes above the keyboard left out`); }
    }
    if (plan.role === "ornament") {
      const before = ns.length;
      const thin: SourceNote[] = [];
      for (const n of ns) if (!thin.length || n.t - thin[thin.length - 1]!.t >= ORNAMENT_GAP) thin.push(n);
      ns = thin;
      if (before > ns.length) plan.actions.push(`thinned from ${before} to ${ns.length} notes (at most ${Math.round(1 / ORNAMENT_GAP)} a second)`);
    }
    const written = ns.map((n) => n.vel);
    const flat = quantile(written, 0.9) - quantile(written, 0.1) <= 2;
    const eff = ns.map((n) => n.vel * n.level);
    const target = TARGET[plan.role];
    let vels: number[];
    if (flat && target) {
      const mid = Math.max(1, median(eff));
      vels = eff.map((e) => target * (e / mid));
      plan.actions.push(`every note written at velocity ${median(written)}: played as ${ROLE_WORD[plan.role] ?? plan.role}, around ${target}`);
    } else {
      vels = eff;
      if (ns.some((n) => n.level < 0.999)) plan.actions.push("dynamics from the file's own velocity, volume and expression");
    }
    let lengthened = 0;
    ns.forEach((n, k) => {
      const off = Math.max(n.off, n.t + MIN_DUR);
      if (off > n.off) lengthened++;
      out.push({ key: n.key, t: n.t, off, vel: Math.max(8, Math.min(118, Math.round(vels[k]!))), part: i, prio: PRIO[plan.role] ?? 3 });
    });
    if (lengthened > ns.length * 0.2) plan.actions.push(`${lengthened} notes shorter than ${MIN_DUR * 1000} ms held to ${MIN_DUR * 1000} ms`);
    plan.notesOut = ns.length;
    if (!plan.actions.length) plan.actions.push("played as written");
  }
  for (const [i, plan] of plans.entries()) if (plan.role === "silent" && !plan.actions.length && parts[i]!.notes.length) plan.actions.push("left out: silent");
  return { notes: out.sort((a, b) => a.t - b.t || a.key - b.key), plans, silenced };
}

/**
 * How much of part `a` is a copy of part `b`: the share of its notes that start within 15 ms of a
 * note of `b` with the same pitch class within two octaves (a doubling), and the largest share
 * that follows such a note by one constant delay of 60-600 ms (an echo).
 */
function relation(a: SourceNote[], b: SourceNote[]): { same: number; echo: number; delay: number } {
  if (!a.length || !b.length) return { same: 0, echo: 0, delay: 0 };
  let same = 0;
  const delays = new Map<number, number>();
  for (const n of a) {
    let lo = 0, hi = b.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (b[m]!.t < n.t - 0.6) lo = m + 1; else hi = m; }
    let matched = false;
    const seen = new Set<number>();
    for (let j = lo; j < b.length && b[j]!.t <= n.t + 0.015; j++) {
      const q = b[j]!, d = n.t - q.t, dk = Math.abs(n.key - q.key);
      if (dk % 12 !== 0 || dk > 24) continue;
      if (Math.abs(d) <= 0.015) matched = true;
      else if (d >= 0.06) { const bin = Math.round(d * 50); if (!seen.has(bin)) { seen.add(bin); delays.set(bin, (delays.get(bin) ?? 0) + 1); } }
    }
    if (matched) same++;
  }
  let best = 0, delay = 0;
  for (const [bin, c] of delays) if (c > best) { best = c; delay = bin / 50; }
  return { same: same / a.length, echo: best / a.length, delay };
}

function polyphony(ns: SourceNote[]): number {
  const e: [number, number][] = [];
  for (const n of ns) e.push([n.t, 1], [n.off, -1]);
  e.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  let c = 0, p = 0;
  for (const [, d] of e) { c += d; p = Math.max(p, c); }
  return p;
}
