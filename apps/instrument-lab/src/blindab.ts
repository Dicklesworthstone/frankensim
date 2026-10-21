// The owner's blind A/B on a pair of renders that differ in one thing: which of X and Y is closer to the real instrument, whether either
// sounds like it, and a note. Pure functions: the random order, the text a recorded answer adds to his notes, and reading it back, so what is
// stored (and so what "Copy all my notes" carries) can be tested without a page. Nothing here scores anything or says which clip is better.
export type Slot = "A" | "B"; // the lane's own names for its two renders (pair.json's a and b); the owner never sees them before he answers
export type Order = [x: Slot, y: Slot];
export type Closer = "X" | "Y" | "same";
export type Sounds = "none" | "X" | "Y" | "both";
export interface Answer { order: Order; closer: Closer; sounds: Sounds }
export interface Meta {
  id: string; set: string; placeholder: boolean;
  a: { label: string; commit: string }; b: { label: string; commit: string };
}

export const CLOSER: readonly Closer[] = ["X", "Y", "same"];
export const SOUNDS: readonly Sounds[] = ["none", "X", "Y", "both"];

/** The notes key of a pair. */
export const abKey = (id: string): string => `blind-ab/${id}`;

/** A fresh random order. `rand` is injectable so the two outcomes can be tested; the page passes a crypto-backed one. */
export function randomOrder(rand: () => number): Order {
  return rand() < 0.5 ? ["A", "B"] : ["B", "A"];
}

/** A random number in [0, 1) from the browser's crypto source, else Math.random. */
export function cryptoRandom(): number {
  try {
    const v = new Uint32Array(1);
    crypto.getRandomValues(v);
    return v[0]! / 2 ** 32;
  } catch { return Math.random(); }
}

/** Which of the lane's two renders plays under a letter. */
export const slotOf = (letter: "X" | "Y", order: Order): Slot => order[letter === "X" ? 0 : 1];

const label = (m: Meta, s: Slot): string => (s === "A" ? m.a.label : m.b.label);

/**
 * The note for a recorded answer. Everything above "note:" is generated and is what parseNote reads back; the owner's own words come
 * after "note:" and can say anything without being mistaken for an answer. The header names both renders, their commits and the order
 * he was shown, so the analysis can decode X and Y from the copied text alone.
 */
export function composeNote(m: Meta, a: Answer, comment: string): string {
  const who = (l: "X" | "Y"): string => `${slotOf(l, a.order)}: ${label(m, slotOf(l, a.order))}`;
  const closer = a.closer === "same" ? "same (equally close)" : `${a.closer} (${who(a.closer)})`;
  const sounds = a.sounds === "none" ? "none" : a.sounds === "both" ? `both (X = ${who("X")}; Y = ${who("Y")})` : `${a.sounds} (${who(a.sounds)})`;
  const lines = [
    `blind A/B pair: ${m.id} (${m.set})${m.placeholder ? " PLACEHOLDER, NOT A TEST" : ""}`,
    `A: ${m.a.label} (commit ${m.a.commit})`,
    `B: ${m.b.label} (commit ${m.b.commit})`,
    `shown as: X=${a.order[0]}, Y=${a.order[1]}`,
    `closer to the real instrument: ${closer}`,
    `sounds like the real instrument: ${sounds}`,
  ];
  return comment.trim() ? `${lines.join("\n")}\nnote: ${comment.trim()}` : lines.join("\n");
}

/** The answer a note holds and the words after "note:"; answer is null when the note is not a recorded answer. */
export function parseNote(note: string | undefined): { answer: Answer | null; comment: string } {
  const text = note ?? "";
  const cut = text.search(/\nnote:/);
  const head = cut === -1 ? text : text.slice(0, cut);
  const comment = cut === -1 ? "" : text.slice(cut + "\nnote:".length).replace(/^[ \t]/, "");
  const shown = /^shown as: X=([AB]), Y=([AB])[ \t]*$/m.exec(head);
  const closer = /^closer to the real instrument: (X|Y|same)\b/m.exec(head);
  const sounds = /^sounds like the real instrument: (none|X|Y|both)\b/m.exec(head);
  if (!shown || !closer || !sounds || shown[1] === shown[2]) return { answer: null, comment: text };
  return { answer: { order: [shown[1] as Slot, shown[2] as Slot], closer: closer[1] as Closer, sounds: sounds[1] as Sounds }, comment };
}
