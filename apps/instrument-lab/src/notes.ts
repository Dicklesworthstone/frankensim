// The owner's own notes, one store for the whole page: a card's note, a blind guess, a candidate's verdict all live in this one map
// (key -> text) and one localStorage entry, so "Copy all my notes" from any section carries all of them. Two copies of the map would
// overwrite each other on save, hence the single instance.
export const NOTES_KEY = "struck:instruments:notes";

function read(): Record<string, string> {
  try { return JSON.parse(localStorage.getItem(NOTES_KEY) ?? "{}") as Record<string, string>; } catch { return {}; }
}

export const notes: Record<string, string> = read();

export function saveNotes(): void {
  try { localStorage.setItem(NOTES_KEY, JSON.stringify(notes)); } catch { /* storage unavailable: notes stay on the page only */ }
}

/** Everything he has written, as one block of text for the clipboard; "" when there is nothing. */
export function allNotesText(): string {
  const entries = Object.entries(notes).filter(([, v]) => v.trim());
  return entries.length ? `# Instrument listening notes\n\n${entries.map(([k, v]) => `## ${k}\n${v.trim()}`).join("\n\n")}\n` : "";
}

export const noteCount = (): number => Object.values(notes).filter((v) => v.trim()).length;
