// The "Blind A/B" section: pairs of renders that differ in one thing (tools/build_blind_ab.py), each judged by the owner without knowing which
// is which. Two clips play as X and Y in a random order drawn in this browser, the real reference has its own button, and he answers
// (a) which is closer to the real instrument (X / Y / same) and (b) whether either sounds like it (neither / X / Y / both), plus a free note.
// Only after the answer is recorded does the page say what X and Y were and what differed. The answer is stored in the same local notes as
// everything else he writes (notes.ts), with the pair id, both labels, both commits and the order he was shown, so "Copy all my notes"
// carries it. Nothing here scores or interprets an answer.
import "./instruments.css";
import { GroupPlayer } from "./audio";
import { allNotesText, notes, saveNotes } from "./notes";
import { abKey, composeNote, cryptoRandom, parseNote, randomOrder, slotOf, type Answer, type Closer, type Meta, type Order, type Slot, type Sounds } from "./blindab";
import { asset, el, json, onVisible } from "./util";

interface Level { audio: string; gain_db: number; rms_dbfs: number; peak_dbfs: number; source_duration_s: number; playback_s: number }
interface Render extends Level { label: string; commit: string; recipe: string }
interface Pair {
  id: string; set: string; title: string; question: string; change: string; made_by: string; placeholder: boolean;
  level_match: "rms" | "none"; target_rms_dbfs: number | null;
  a: Render; b: Render; reference: Level & { label: string; source: string; licence?: string };
}
interface Index { built_utc: string; pairs: Pair[]; refused: { id: string; problems: string[] }[] }

const dataUrl = (path: string): string => asset(`data/blind-ab/${path}`);
const plural = (n: number, one: string, many = `${one}s`): string => `${n} ${n === 1 ? one : many}`;

// The order a pair was shown in: drawn at random the first time the pair is opened and kept, so it does not change under him while he
// listens or after a reload; once he has answered, the order is the one written in his note.
const ORDER_KEY = "struck:blindab:order";
const memory: Record<string, Order> = {};
const storedOrders = (): Record<string, string> => { try { return JSON.parse(localStorage.getItem(ORDER_KEY) ?? "{}") as Record<string, string>; } catch { return {}; } };
function keepOrder(id: string, order: Order): Order {
  memory[id] = order;
  try { const all = storedOrders(); all[id] = order.join(""); localStorage.setItem(ORDER_KEY, JSON.stringify(all)); } catch { /* storage unavailable: the order lives for this page view */ }
  return order;
}
function orderOf(id: string): Order {
  const recorded = parseNote(notes[abKey(id)]).answer;
  if (recorded) return recorded.order;
  const known = memory[id];
  if (known) return known;
  const s = storedOrders()[id];
  if (s === "AB" || s === "BA") return keepOrder(id, [s[0] as Slot, s[1] as Slot]);
  return keepOrder(id, randomOrder(cryptoRandom));
}

const CLOSER_TEXT: Record<Closer, string> = { X: "Sound X", Y: "Sound Y", same: "Same: equally close" };
const SOUNDS_TEXT: Record<Sounds, string> = { none: "Neither", X: "Only X", Y: "Only Y", both: "Both" };

function choice(group: string, legend: string, texts: Record<string, string>): { box: HTMLElement; value: () => string | null; lock: (v: string | null, locked: boolean) => void; onChange: (f: () => void) => void } {
  const inputs = Object.entries(texts).map(([v, t]) => {
    const input = el("input", { type: "radio", name: group, value: v });
    return { v, input, label: el("label", { class: "toggle" }, input, t) };
  });
  const box = el("fieldset", { class: "ir-ab-ask" }, el("legend", {}, legend), el("div", { class: "ir-ab-choices" }, ...inputs.map((i) => i.label)));
  return {
    box,
    value: () => inputs.find((i) => i.input.checked)?.v ?? null,
    lock: (v, locked) => { for (const i of inputs) { i.input.checked = i.v === v; i.input.disabled = locked; } },
    onChange: (f) => { for (const i of inputs) i.input.addEventListener("change", f); },
  };
}

function pairCard(p: Pair, refresh: () => void): HTMLElement {
  const key = abKey(p.id);
  const meta: Meta = { id: p.id, set: p.set, placeholder: p.placeholder, a: { label: p.a.label, commit: p.a.commit }, b: { label: p.b.label, commit: p.b.commit } };
  const render = (s: Slot): Render => (s === "A" ? p.a : p.b);
  const root = el("article", { class: "panel ir-ab", "data-id": p.id, "data-state": "open" });
  let order = orderOf(p.id);
  let player: GroupPlayer | null = null;

  const playX = el("button", { class: "btn sm", type: "button" });
  const playY = el("button", { class: "btn sm", type: "button" });
  const playRef = el("button", { class: "btn sm", type: "button" });
  const stop = el("button", { class: "btn ghost sm", type: "button" }, "Stop");
  const closer = choice(`closer-${p.id}`, "Which is closer to the real instrument?", CLOSER_TEXT);
  const sounds = choice(`sounds-${p.id}`, "Does either sound like the real instrument?", SOUNDS_TEXT);
  const comment = el("textarea", { class: "ir-feedback", rows: "3", placeholder: "Anything to add? What do you hear that is wrong, or right? Kept in this browser; “Copy all my notes” sends it.", "aria-label": `Your note on ${p.title}` });
  const record = el("button", { class: "btn sm ir-ab-record", type: "button" }, "Record my answer and reveal");
  const again = el("button", { class: "btn ghost sm", type: "button" }, "Clear my answer and test again");
  const state = el("span", { class: "small muted", "aria-live": "polite" });
  const reveal = el("div", { class: "small ir-ab-reveal", "aria-live": "polite" });
  reveal.hidden = true;

  const labelButtons = (): void => {
    const now = player?.state === "playing" ? player.current : null;
    playX.textContent = `${now === "x" ? "■" : "▶"} Sound X`;
    playY.textContent = `${now === "y" ? "■" : "▶"} Sound Y`;
    playRef.textContent = `${now === "ref" ? "■" : "▶"} Real: ${p.reference.label}`;
  };
  const makePlayer = (): void => {
    player?.pause();
    player = new GroupPlayer({ x: dataUrl(render(order[0]).audio), y: dataUrl(render(order[1]).audio), ref: dataUrl(p.reference.audio) }, "x");
    player.onState = labelButtons;
    labelButtons();
  };
  const playClip = (id: "x" | "y" | "ref"): void => {
    if (!player) return;
    if (player.state === "playing" && player.current === id) { player.pause(); return; }
    player.select(id);
    if (player.state !== "playing") void player.play(0).catch(() => undefined);
  };
  playX.addEventListener("click", () => playClip("x"));
  playY.addEventListener("click", () => playClip("y"));
  playRef.addEventListener("click", () => playClip("ref"));
  stop.addEventListener("click", () => player?.pause());

  const store = (text: string): void => { if (text.trim()) notes[key] = text; else delete notes[key]; saveNotes(); };
  const readyToRecord = (): boolean => closer.value() !== null && sounds.value() !== null;

  const showReveal = (a: Answer): void => {
    const x = render(slotOf("X", a.order)), y = render(slotOf("Y", a.order));
    const sha = (r: Render): string => r.commit.slice(0, 12);
    reveal.replaceChildren(
      el("p", {}, el("b", {}, "Recorded. "), `Closer to the real instrument: ${CLOSER_TEXT[a.closer]}. Sounds like it: ${SOUNDS_TEXT[a.sounds]}.`),
      el("p", {}, `X was ${x.label} (render ${sha(x)}); Y was ${y.label} (render ${sha(y)}).`),
      el("p", {}, el("b", {}, "What differed: "), p.change),
      el("details", { class: "small muted" }, el("summary", {}, "Recipes and levels"),
        el("p", {}, `X: ${x.recipe}`), el("p", {}, `Y: ${y.recipe}`),
        el("p", {}, `Reference: ${p.reference.label}, from ${p.reference.source}${p.reference.licence ? `. Licence: ${p.reference.licence}` : ""}.`),
        el("p", {}, p.level_match === "rms"
          ? `Level-matched: X, Y and the reference were all brought to ${p.target_rms_dbfs} dBFS RMS (gains: X ${render(slotOf("X", a.order)).gain_db} dB, Y ${render(slotOf("Y", a.order)).gain_db} dB, reference ${p.reference.gain_db} dB).`
          : "Not level-matched: the clips play as rendered.")),
      el("p", { class: "muted" }, `Saved in this browser's notes as “${key}”, with the pair id, both labels, both commits and the order you were shown; “Copy all my notes” sends it.`));
    reveal.hidden = false;
  };

  const paint = (): void => {
    const parsed = parseNote(notes[key]);
    const a = parsed.answer;
    root.dataset["state"] = a ? "answered" : "open";
    closer.lock(a?.closer ?? closer.value(), a !== null);
    sounds.lock(a?.sounds ?? sounds.value(), a !== null);
    record.hidden = a !== null;
    record.disabled = !readyToRecord();
    again.hidden = a === null;
    if (a) { comment.value = parsed.comment; showReveal(a); state.textContent = ""; }
    else { reveal.hidden = true; reveal.replaceChildren(); state.textContent = readyToRecord() ? "" : "Answer both questions to record."; }
  };

  const answerChanged = (): void => { record.disabled = !readyToRecord(); state.textContent = readyToRecord() ? "" : "Answer both questions to record."; };
  closer.onChange(answerChanged);
  sounds.onChange(answerChanged);
  record.addEventListener("click", () => {
    const c = closer.value() as Closer | null, s = sounds.value() as Sounds | null;
    if (c === null || s === null) return;
    store(composeNote(meta, { order, closer: c, sounds: s }, comment.value));
    paint();
    refresh();
  });
  again.addEventListener("click", () => {
    player?.pause();
    store("");
    order = keepOrder(p.id, randomOrder(cryptoRandom));
    comment.value = "";
    closer.lock(null, false);
    sounds.lock(null, false);
    makePlayer();
    paint();
    refresh();
  });
  comment.addEventListener("input", () => {
    const a = parseNote(notes[key]).answer;
    if (a) store(composeNote(meta, a, comment.value)); // before he records, the words wait in the box and go in with the answer
  });

  root.append(
    el("div", { class: "panel-head" }, el("h4", {}, p.title), el("div", { class: "tags" }, el("span", { class: "tag plain" }, p.set), p.placeholder ? el("span", { class: "tag offline" }, "placeholder") : null)),
    ...(p.placeholder ? [el("p", { class: "ir-ab-placeholder small" }, "Placeholder: two existing renders so this page can be tried before the first real pair arrives. It tests nothing; an answer here is stored with the words PLACEHOLDER, NOT A TEST.")] : []),
    el("p", {}, p.question),
    el("p", { class: "small muted" }, p.level_match === "rms"
      ? "The two sounds and the real reference play at the same loudness, so loudness is not a clue. Listen as often as you like, in any order."
      : "These clips play as rendered: their loudness has not been matched. Listen as often as you like, in any order."),
    el("div", { class: "btn-row" }, playX, playY, playRef, stop),
    closer.box, sounds.box, comment,
    el("div", { class: "btn-row" }, record, again, state),
    reveal);
  makePlayer();
  paint();
  return root;
}

export async function mountBlindAb(root: HTMLElement): Promise<void> {
  const section = root.closest("section");
  let index: Index;
  try { index = await json<Index>(dataUrl("index.json")); } catch {
    if (section) section.hidden = true; // a build without blind A/B pairs has no such section
    return;
  }
  if (!index.pairs.length && !index.refused.length) { if (section) section.hidden = true; return; }
  const summary = el("p", { class: "ir-cand-summary", "aria-live": "polite" });
  const copied = el("span", { class: "small muted", "aria-live": "polite" });
  const copy = el("button", { class: "btn sm", type: "button" }, "Copy all my notes");
  const refresh = (): void => {
    const answered = index.pairs.filter((p) => parseNote(notes[abKey(p.id)]).answer !== null).length;
    const real = index.pairs.filter((p) => !p.placeholder).length;
    summary.textContent = `${plural(index.pairs.length, "pair")}${real !== index.pairs.length ? ` (${real} real, ${index.pairs.length - real} placeholder)` : ""}. You have answered ${answered}; ${index.pairs.length - answered} still to do.`;
  };
  copy.addEventListener("click", () => {
    const text = allNotesText();
    if (!text) { copied.textContent = "Nothing written yet."; return; }
    void navigator.clipboard.writeText(text).then(() => { copied.textContent = `Copied ${plural(text.split("\n## ").length - 1, "note")}.`; }, () => { copied.textContent = "The browser would not copy; select the text of each box instead."; });
  });
  const refused = index.refused.length
    ? el("div", { class: "ir-refused small", role: "status" }, el("b", {}, `${plural(index.refused.length, "pair")} refused, not shown: `), ...index.refused.map((r) => el("div", {}, `${r.id}: ${r.problems.join("; ")}`)))
    : null;
  root.replaceChildren(summary, el("div", { class: "btn-row" }, copy, copied), refused ?? "", ...index.pairs.map((p) => pairCard(p, refresh)),
    el("p", { class: "small muted", style: "margin-top: 18px" }, `Pairs built ${index.built_utc.replace("T", " ").replace("Z", " UTC")}. Nothing on this page scores an answer.`));
  root.style.minHeight = "0";
  refresh();
}

const root = document.getElementById("blind-ab-root");
if (root) onVisible(root, () => void mountBlindAb(root), "800px");
