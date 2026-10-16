// Chord names for what the roll is playing.
//
// Provenance: the template table and the ranking law are copied from jazz_chord_progression_editor_html
// (same owner, commit 0c3c7aa): src/theory/reverse-chord-templates.ts (values and order unchanged)
// and `resolveSonority` / `groupSonorities` in src/export/midi-import.ts (the 40 ms simultaneity
// window in ticks from the file's own tempo; candidate roots only among the pitch classes present;
// exact before omitted-fifth, root position before slash, fewer tones, table order). The symbol
// formatter at the bottom is new: that project spells symbols through its T1 grammar, which the
// player does not carry.

interface Template { formula: string; alt: string | null; offsets: readonly number[]; omissibleFifth: boolean }
const T = (formula: string, alt: string | null, offsets: number[], omissibleFifth: boolean): Template => ({ formula, alt, offsets, omissibleFifth });

/** M0 reverse-T1 templates, in ranking-tiebreak order. */
const TEMPLATES: readonly Template[] = [
  T("base-power", null, [0, 7], false),
  T("base-major", null, [0, 4, 7], false),
  T("base-minor", null, [0, 3, 7], false),
  T("base-diminished", null, [0, 3, 6], false),
  T("base-augmented", null, [0, 4, 8], false),
  T("base-sus2", null, [0, 2, 7], false),
  T("base-sus4", null, [0, 5, 7], false),
  T("sixth-major", null, [0, 4, 7, 9], true),
  T("sixth-minor", null, [0, 3, 7, 9], true),
  T("seventh-major", null, [0, 4, 7, 11], true),
  T("seventh-dominant", null, [0, 4, 7, 10], true),
  T("seventh-minor", null, [0, 3, 7, 10], true),
  T("seventh-minor-major", null, [0, 3, 7, 11], true),
  T("seventh-half-diminished", null, [0, 3, 6, 10], false),
  T("seventh-diminished", null, [0, 3, 6, 9], false),
  T("seventh-augmented-major", null, [0, 4, 8, 11], false),
  T("extension-suspended-dominant", null, [0, 5, 7, 10], true),
  T("extension-major", null, [0, 2, 4, 7, 11], true),
  T("extension-dominant", null, [0, 2, 4, 7, 10], true),
  T("extension-minor", null, [0, 2, 3, 7, 10], true),
  T("altered-dominant", "alt-b9-b5", [0, 1, 4, 6, 10], false),
  T("altered-dominant", "alt-b9-sharp5", [0, 1, 4, 8, 10], false),
  T("altered-dominant", "alt-sharp9-b5", [0, 3, 4, 6, 10], false),
  T("altered-dominant", "alt-sharp9-sharp5", [0, 3, 4, 8, 10], false),
];

const SUFFIX: Record<string, string> = {
  "base-power": "5", "base-major": "", "base-minor": "m", "base-diminished": "dim", "base-augmented": "+",
  "base-sus2": "sus2", "base-sus4": "sus4", "sixth-major": "6", "sixth-minor": "m6", "seventh-major": "maj7",
  "seventh-dominant": "7", "seventh-minor": "m7", "seventh-minor-major": "m(maj7)", "seventh-half-diminished": "m7♭5",
  "seventh-diminished": "dim7", "seventh-augmented-major": "maj7♯5", "extension-suspended-dominant": "7sus4",
  "extension-major": "maj9", "extension-dominant": "9", "extension-minor": "m9",
};
const ALT: Record<string, string> = { "alt-b9-b5": "7♭9♭5", "alt-b9-sharp5": "7♭9♯5", "alt-sharp9-b5": "7♯9♭5", "alt-sharp9-sharp5": "7♯9♯5" };
/** The flat-preferring jazz spelling M0 uses for names it cannot otherwise spell. */
const NAMES = ["C", "D♭", "D", "E♭", "E", "F", "G♭", "G", "A♭", "A", "B♭", "B"];

/** Name a sonority (its pitch classes and its lowest pitch class), or null when no template fits. */
export function nameSonority(pitchClasses: readonly number[], bass: number): string | null {
  const pcs = [...new Set(pitchClasses.map((p) => ((p % 12) + 12) % 12))].sort((a, b) => a - b);
  if (pcs.length < 3) return null; // two notes are an interval, not a chord: say nothing
  let best: Rank | null = null;
  for (const root of pcs) {
    const off = pcs.map((p) => (p - root + 12) % 12).sort((a, b) => a - b);
    for (let order = 0; order < TEMPLATES.length; order++) {
      const t = TEMPLATES[order]!;
      let exact = -1;
      if (same(off, t.offsets)) exact = 0;
      else if (t.omissibleFifth && same(off, t.offsets.filter((o) => o !== 7))) exact = 1;
      if (exact < 0) continue;
      const cand: Rank = { exact, root: root === bass ? 0 : 1, tones: t.offsets.length, order, rootPc: root, t };
      if (!best || cmp(cand, best) < 0) best = cand;
    }
  }
  if (!best) return null;
  const name = NAMES[best.rootPc]! + (best.t.alt ? ALT[best.t.alt]! : SUFFIX[best.t.formula] ?? "");
  return best.rootPc === bass ? name : `${name}/${NAMES[bass]}`;
}
interface Rank { exact: number; root: number; tones: number; order: number; rootPc: number; t: Template }
function cmp(a: Rank, b: Rank): number {
  return a.exact - b.exact || a.root - b.root || a.tones - b.tones || a.order - b.order || a.rootPc - b.rootPc;
}
function same(a: readonly number[], b: readonly number[]): boolean {
  return a.length === b.length && a.every((v, i) => v === b[i]);
}

export interface TimedNote { t: number; key: number; off: number }
/**
 * Chord labels along the music: notes whose onsets fall within 40 ms of an anchor form one
 * sonority (M0's window, anchored, never chained); the label is computed from every pitch sounding
 * at that moment (held notes included), and emitted when it changes.
 */
export function chordTimeline(notes: readonly TimedNote[]): { t: number; label: string }[] {
  const sorted = [...notes].sort((a, b) => a.t - b.t);
  const out: { t: number; label: string }[] = [];
  let i = 0;
  while (i < sorted.length) {
    const anchor = sorted[i]!.t;
    let j = i;
    while (j < sorted.length && sorted[j]!.t - anchor <= 0.04) j++;
    const at = anchor + 0.02;
    const sounding: number[] = [];
    for (let k = Math.max(0, i - 64); k < j; k++) {
      const n = sorted[k]!;
      if (n.t <= at && n.off > at) sounding.push(n.key);
    }
    if (sounding.length >= 3) {
      const bass = Math.min(...sounding);
      const label = nameSonority(sounding, bass % 12);
      if (label && out[out.length - 1]?.label !== label) out.push({ t: anchor, label });
    }
    i = j;
  }
  return out;
}
