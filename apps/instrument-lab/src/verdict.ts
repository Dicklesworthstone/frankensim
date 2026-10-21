// The owner's verdict on a reference candidate, kept inside the same local notes as everything else he writes: the first line of
// a candidate's note is "verdict: VALID" or "verdict: INVALID", and the rest is whatever he typed. Pure functions, so what is
// stored (and so what "Copy all my notes" carries) can be tested without a page. Nothing here decides anything for him.
export type Verdict = "VALID" | "INVALID";

const LINE = /^verdict: (VALID|INVALID)[ \t]*(?:\r?\n|$)/;

/** The verdict at the top of a note, or null while he has not decided. */
export function verdictOf(note: string | undefined): Verdict | null {
  const m = LINE.exec(note ?? "");
  return m ? (m[1] as Verdict) : null;
}

/** What he typed apart from the verdict line. */
export const commentOf = (note: string | undefined): string => (note ?? "").replace(LINE, "");

/** A note holding this verdict (null clears it) and this comment. An empty result is "" so the caller can drop the key. */
export function compose(verdict: Verdict | null, comment: string): string {
  const text = comment.replace(/^\s*\n/, "");
  if (verdict === null) return text;
  return text.trim() ? `verdict: ${verdict}\n${text}` : `verdict: ${verdict}`;
}

/** The same note with its verdict replaced, comment kept. */
export const withVerdict = (note: string | undefined, verdict: Verdict | null): string => compose(verdict, commentOf(note));

/** The same note with its comment replaced, verdict kept. */
export const withComment = (note: string | undefined, comment: string): string => compose(verdictOf(note), comment);

/** The notes key of a candidate. */
export const candidateKey = (id: string): string => `candidate/${id}`;

/** How many of these candidates he has marked each way. */
export function tally(ids: string[], notes: Record<string, string>): { valid: number; invalid: number; undecided: number } {
  let valid = 0, invalid = 0;
  for (const id of ids) {
    const v = verdictOf(notes[candidateKey(id)]);
    if (v === "VALID") valid++; else if (v === "INVALID") invalid++;
  }
  return { valid, invalid, undecided: ids.length - valid - invalid };
}
