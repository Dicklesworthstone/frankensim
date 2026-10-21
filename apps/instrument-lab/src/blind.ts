// The blind A/B on an #instruments card: which clip plays as "X", and the line a guess adds to the owner's own notes.
// Pure functions, so the order and the stored line can be tested without a page. Nothing here scores anything.
export type Clip = "real" | "ours";
export type Order = [x: Clip, y: Clip];

/** 32-bit FNV-1a of a string. */
function fnv1a(s: string): number {
  let h = 0x811c9dc5;
  for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 0x01000193); }
  return h >>> 0;
}

/**
 * Which clip is X and which is Y for a note. Seeded by the note's key and by how many times the owner has asked for another
 * shuffle of it, so a note opens in the same order every time until he presses "Shuffle again".
 */
export function blindOrder(key: string, shuffles: number): Order {
  let h = fnv1a(`${key}#${shuffles}`);
  h ^= h >>> 16; h = Math.imul(h, 0x85ebca6b); h ^= h >>> 13; // FNV's low bit barely moves for similar short keys: mix before taking one
  return (h & 1) === 0 ? ["real", "ours"] : ["ours", "real"];
}

/** The letter a clip plays under. */
export const letterOf = (clip: Clip, order: Order): "X" | "Y" => (order[0] === clip ? "X" : "Y");

/** The line a guess adds to the note's notes, so "Copy all my notes" carries it. `picked` is the slot he chose: 0 = X, 1 = Y. */
export const blindLine = (order: Order, picked: 0 | 1): string => `blind: picked ${picked === 0 ? "X" : "Y"}, real was ${letterOf("real", order)}`;
