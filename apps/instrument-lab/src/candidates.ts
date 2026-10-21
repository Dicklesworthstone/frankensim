// The "Reference candidates" section: real recordings found for every instrument family (tools/build_candidates.py), each with the
// facts known about it and two buttons for the owner, VALID and INVALID. His marks live in the same local notes as everything else
// (verdict.ts), so "Copy all my notes" carries them. Nothing here is used for grading until he has marked it.
import "./instruments.css";
import { GroupPlayer } from "./audio";
import { allNotesText, notes, saveNotes } from "./notes";
import { candidateKey, commentOf, tally, verdictOf, withComment, withVerdict, type Verdict } from "./verdict";
import { asset, el, fmt, json, onVisible } from "./util";

interface SteadyNote { note: string; f0_hz: number; cents_from_et_a440: number; pitch_std_cents: number; level_dbfs: number; clip_t0: number; clip_t1: number }
interface Facts {
  codec?: string; channels_in_file?: number; peak_dbfs?: number; rms_dbfs?: number; bandwidth_99pct_hz?: number; noise_floor_dbfs_p5_of_50ms_frames?: number;
  level_range_db?: number; decay_after_loudest_db_per_s?: number | null;
  pitch?: { median_hz: number; nearest_note: string; cents_from_note: number; fraction_within_30_cents_of_median: number } | null;
  steady_notes_by_pitch_tracker?: SteadyNote[];
}
interface Candidate {
  id: string; instrument: string; notes: string[]; dynamic: string; articulation: string; source_type: string; url: string; title: string; channel_or_author: string;
  licence: string; licence_evidence: string; start_s: number; end_s: number; conditions: string; sample_rate_hz: number; reason: string; render_coverage: { graded: string[]; reference_only: string[]; other: string[]; text: string | null };
  creative_commons: boolean; grader_reference: boolean; audio: string; clip_duration_s: number; playback_s: number; playback_gain_db: number; facts: Facts | null;
}
interface Looked { source: string; query_or_url: string; result: string; file?: string }
interface Family {
  id: string; title: string; candidates: Candidate[]; omitted: string[]; looked: Looked[]; made_by: string[];
  none_found: { reason: string; looked_summary: string } | null;
  render_notes: { rendered_and_graded: string[]; reference_only: string[]; listening_only_or_refused: string[] } | null;
  refused: { id: string; problems: string[]; file?: string }[];
}
interface Index { built_utc: string; families: Family[] }

const dataUrl = (path: string): string => asset(`data/candidates/${path}`);
const clock = (s: number): string => `${Math.floor(s / 60)}:${(Math.floor(s) % 60).toString().padStart(2, "0")}`;
const plural = (n: number, one: string, many = `${one}s`): string => `${n} ${n === 1 ? one : many}`;
const db = (v: number | undefined): string => (v === undefined ? "—" : `${v.toFixed(1)} dBFS`);
const SOURCE_NAME: Record<string, string> = { "open-library": "open library", youtube: "YouTube", freesound: "Freesound", other: "other source" };

function facts(c: Candidate): HTMLElement {
  const f = c.facts ?? {};
  const row = (k: string, ...v: (string | Node)[]): HTMLElement => el("div", { class: "ir-fact" }, el("dt", {}, k), el("dd", {}, ...v));
  const lossy = c.source_type === "youtube" ? " (a YouTube stream: lossy, whatever the uploader recorded)" : "";
  const channels = f.channels_in_file === undefined ? "" : f.channels_in_file === 1 ? ", mono" : `, ${f.channels_in_file} channels`;
  const p = f.pitch;
  const steady = f.steady_notes_by_pitch_tracker ?? [];
  const cov = c.render_coverage;
  const covers = [cov.graded.length ? `graded: ${cov.graded.join(" ")}` : "", cov.reference_only.length ? `reference-only notes: ${cov.reference_only.join(" ")}` : "", cov.text && !cov.graded.length ? cov.text : ""].filter(Boolean).join("; ")
    || "none of the notes our renders are graded on";
  return el("dl", { class: "ir-facts small" },
    row("Source", el("a", { href: c.url, rel: "noopener" }, c.title), ` — ${c.channel_or_author} (${SOURCE_NAME[c.source_type] ?? c.source_type})`),
    row("Licence", el("span", { class: c.creative_commons ? "ir-cc" : "ir-std" }, c.creative_commons ? "Creative Commons" : "Not Creative Commons"), ` ${c.licence}. Read from: ${c.licence_evidence}.`),
    row("Where in the source", `${clock(c.start_s)}–${clock(c.end_s)} (${fmt(c.clip_duration_s)} s excerpt${c.playback_s < c.clip_duration_s ? `; this page plays only its first ${fmt(c.playback_s)} s` : ""})`),
    row("Clip", `${fmt(c.sample_rate_hz / 1000)} kHz${channels}${f.codec ? `, ${f.codec}` : ""}${lossy}`),
    row("Dynamic, articulation", `${c.dynamic}; ${c.articulation}`),
    row("Covers our render notes", covers),
    row("Level (measured)", `peak ${db(f.peak_dbfs)}, RMS ${db(f.rms_dbfs)}, quietest 5% of 50 ms frames ${db(f.noise_floor_dbfs_p5_of_50ms_frames)}${f.decay_after_loudest_db_per_s != null ? `, level falls ${fmt(-f.decay_after_loudest_db_per_s)} dB/s after the loudest part` : ""}; played here at ${c.playback_gain_db >= 0 ? "+" : ""}${c.playback_gain_db} dB so all candidates are equally loud`),
    row("Bandwidth", `99% of the energy is below ${f.bandwidth_99pct_hz === undefined ? "—" : fmt(f.bandwidth_99pct_hz)} Hz`),
    p ? row("Pitch (measured)", `median ${fmt(p.median_hz, 4)} Hz, nearest note ${p.nearest_note} (${p.cents_from_note >= 0 ? "+" : ""}${p.cents_from_note} cents from A440 equal temperament); ${Math.round(p.fraction_within_30_cents_of_median * 100)}% of the frames within 30 cents of the median`) : row("Pitch (measured)", "no steady pitch found"),
    steady.length ? row("Steady notes found", el("details", { class: "ir-notes-found" }, el("summary", {}, plural(steady.length, "note")),
      el("ul", {}, ...steady.map((n) => el("li", {}, `${n.note} at ${fmt(n.clip_t0)}–${fmt(n.clip_t1)} s: ${fmt(n.f0_hz, 4)} Hz (${n.cents_from_et_a440 >= 0 ? "+" : ""}${n.cents_from_et_a440} cents from A440 equal temperament, pitch wanders ${n.pitch_std_cents} cents), ${db(n.level_dbfs)}`))))) : null,
    row("Conditions", c.conditions));
}

function candidateCard(c: Candidate, refresh: () => void): HTMLElement {
  const key = candidateKey(c.id);
  const root = el("article", { class: "panel ir-cand", "data-id": c.id });
  const play = el("button", { class: "btn sm", type: "button" }, "▶ Play the excerpt");
  const player = new GroupPlayer({ clip: dataUrl(c.audio) }, "clip");
  player.onState = () => { play.textContent = player.state === "playing" ? "■ Stop" : "▶ Play the excerpt"; };
  play.addEventListener("click", () => { if (player.state === "playing") player.pause(); else void player.play(0).catch(() => undefined); });
  const valid = el("button", { class: "btn sm ir-valid", type: "button" }, "VALID");
  const invalid = el("button", { class: "btn sm ir-invalid", type: "button" }, "INVALID");
  const clear = el("button", { class: "btn ghost sm", type: "button" }, "Clear");
  const comment = el("textarea", { class: "ir-feedback", rows: "2", placeholder: "Why? (optional). What would you need to hear to call it valid? Kept in this browser; “Copy all my notes” sends it.", "aria-label": `Your note on ${c.id}` });
  const state = el("span", { class: "small muted ir-verdict-state", "aria-live": "polite" });
  const paint = (): void => {
    const v = verdictOf(notes[key]);
    valid.setAttribute("aria-pressed", String(v === "VALID"));
    invalid.setAttribute("aria-pressed", String(v === "INVALID"));
    root.dataset["verdict"] = v ?? "none";
    clear.hidden = v === null;
    state.textContent = v === null ? "Not marked yet." : `Marked ${v}. Saved in this browser's notes as “verdict: ${v}”.`;
  };
  const store = (next: string): void => {
    if (next.trim()) notes[key] = next; else delete notes[key];
    saveNotes();
  };
  const mark = (v: Verdict | null): void => { store(withVerdict(notes[key], v)); paint(); refresh(); };
  valid.addEventListener("click", () => mark(verdictOf(notes[key]) === "VALID" ? null : "VALID"));
  invalid.addEventListener("click", () => mark(verdictOf(notes[key]) === "INVALID" ? null : "INVALID"));
  clear.addEventListener("click", () => mark(null));
  comment.value = commentOf(notes[key]);
  comment.addEventListener("input", () => store(withComment(notes[key], comment.value)));
  const names = c.notes.length > 10 ? `${c.notes.slice(0, 10).join(" ")} +${c.notes.length - 10} more` : c.notes.join(" ");
  root.append(
    el("div", { class: "panel-head" }, el("h4", {}, c.instrument), el("div", { class: "tags" }, c.grader_reference ? el("span", { class: "tag offline", "data-grader-reference": "1" }, "our grader's reference") : null, el("span", { class: "tag plain" }, names), el("span", { class: `tag ${c.creative_commons ? "offline" : "plain"}` }, c.creative_commons ? "Creative Commons" : "not Creative Commons"))),
    el("p", { class: "small muted" }, c.reason),
    el("div", { class: "btn-row" }, play, valid, invalid, clear, state),
    facts(c),
    comment);
  paint();
  return root;
}

function looked(f: Family): HTMLElement | null {
  if (!f.looked.length) return null;
  return el("details", { class: "small muted ir-looked-box" }, el("summary", {}, `Where it was looked for this family (${f.looked.length}${f.made_by.length ? `; by ${f.made_by.join(", ")}` : ""})`),
    el("ul", { class: "ir-looked" }, ...f.looked.map((l) => el("li", {}, el("b", {}, `${l.source}: `), `${l.query_or_url} — ${l.result}`))));
}

export async function mountCandidates(root: HTMLElement): Promise<void> {
  const section = root.closest("section");
  let index: Index;
  try { index = await json<Index>(dataUrl("index.json")); } catch {
    if (section) section.hidden = true; // a build without candidate data has no such section
    return;
  }
  const families = index.families;
  const all = families.flatMap((f) => f.candidates.map((c) => c.id));
  const summary = el("p", { class: "ir-cand-summary", "aria-live": "polite" });
  const copied = el("span", { class: "small muted", "aria-live": "polite" });
  const copy = el("button", { class: "btn sm", type: "button" }, "Copy all my notes");
  const tiles = el("div", { class: "ir-tiles", role: "group", "aria-label": "Choose an instrument family" });
  const host = el("div", { class: "ir-cands" });
  const tileText = (f: Family): string => {
    if (!f.candidates.length) return f.none_found ? "none found" : f.refused.length ? "records refused" : "nothing recorded yet";
    const t = tally(f.candidates.map((c) => c.id), notes);
    return `${plural(f.candidates.length, "candidate")} · ${t.valid} valid · ${t.invalid} invalid`;
  };
  const refresh = (): void => {
    const t = tally(all, notes);
    summary.textContent = `${plural(all.length, "candidate")} in ${families.filter((f) => f.candidates.length).length} families. You have marked ${t.valid} VALID and ${t.invalid} INVALID; ${t.undecided} not marked yet.`;
    for (const b of tiles.querySelectorAll("button")) {
      const f = families.find((x) => x.id === b.dataset["family"]);
      const sm = b.querySelector("small");
      if (f && sm) sm.textContent = tileText(f);
    }
  };
  const show = (id: string): void => {
    for (const b of tiles.querySelectorAll("button")) b.setAttribute("aria-pressed", String(b.dataset["family"] === id));
    const f = families.find((x) => x.id === id)!;
    const parts: HTMLElement[] = [];
    const rn = f.render_notes;
    if (rn) parts.push(el("p", { class: "small muted" }, `Our renders of this family are graded on ${rn.rendered_and_graded.length ? rn.rendered_and_graded.join(" ") : "no note yet"}; the reference library has recordings of ${rn.reference_only.length ? rn.reference_only.join(" ") : "no other note"}. A candidate that covers the graded notes is the most useful.`));
    const graderRefs = f.candidates.filter((c) => c.grader_reference).length;
    if (graderRefs) parts.push(el("p", { class: "small muted", "data-grader-reference-line": "1" }, `${plural(graderRefs, "recording")} below ${graderRefs === 1 ? "is" : "are"} the reference our grader has used for a graded note of this family. Every grade so far was measured against ${graderRefs === 1 ? "it" : "them"}, and none has been confirmed by the owner yet.`));
    if (f.refused.length) parts.push(el("div", { class: "ir-refused small", role: "status" }, el("b", {}, `${plural(f.refused.length, "record")} refused, not shown: `), ...f.refused.map((r) => el("div", {}, `${r.id}${r.file ? ` (${r.file})` : ""}: ${r.problems.join("; ")}`))));
    if (f.candidates.length) parts.push(...f.candidates.map((c) => candidateCard(c, refresh)));
    else if (f.none_found) parts.push(el("div", { class: "ir-empty" }, el("h4", {}, `None found for ${f.title.toLowerCase()}`), el("p", {}, f.none_found.reason), f.none_found.looked_summary ? el("p", { class: "small muted" }, f.none_found.looked_summary) : null));
    else parts.push(el("div", { class: "ir-empty" }, el("h4", {}, `Nothing recorded yet for ${f.title.toLowerCase()}`), el("p", {}, "No search and no candidate has been recorded for this family. That is not the same as none found.")));
    if (f.omitted.length) parts.push(el("div", { class: "ir-refused small", role: "status" }, el("b", {}, `${plural(f.omitted.length, "more candidate")} found, not in this build: `),
      `to keep within the page host's file limit only the ${f.candidates.length} best ranked are shown (${f.omitted.join(", ")}). They are kept with the others and are not marked.`));
    const box = looked(f);
    if (box) parts.push(box);
    host.replaceChildren(...parts);
  };
  for (const f of families) {
    const kind = f.candidates.length ? "ungraded" : f.none_found ? "blocked" : "none";
    const b = el("button", { type: "button", class: `ir-tile ir-${kind}`, "data-family": f.id, "aria-pressed": "false" }, el("b", {}, f.title), el("small", {}, tileText(f)));
    b.addEventListener("click", () => show(f.id));
    tiles.append(b);
  }
  copy.addEventListener("click", () => {
    const text = allNotesText();
    if (!text) { copied.textContent = "Nothing written yet."; return; }
    void navigator.clipboard.writeText(text).then(() => { copied.textContent = `Copied ${plural(text.split("\n## ").length - 1, "note")}.`; }, () => { copied.textContent = "The browser would not copy; select the text of each box instead."; });
  });
  root.replaceChildren(summary, el("div", { class: "btn-row" }, copy, copied), tiles, host,
    el("p", { class: "small muted", style: "margin-top: 18px" }, `Candidates built ${index.built_utc.replace("T", " ").replace("Z", " UTC")}. Everything shown is measured from the clip or copied from the source's own page; none of it was listened to by whoever found it.`));
  root.style.minHeight = "0";
  refresh();
  show((families.find((f) => f.candidates.length) ?? families[0]!).id);
}

const root = document.getElementById("candidates-root");
if (root) onVisible(root, () => void mountCandidates(root), "800px");
