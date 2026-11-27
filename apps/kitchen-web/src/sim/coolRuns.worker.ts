/**
 * The cooling figure's runs, off the page's main thread (sim/coolDown.ts): the four ways for a batch, open or
 * covered. A full pan cools the same however many there are, so its run is kept, and so is each part-filled pan's by
 * what it holds; only the pot runs again for every new batch. Only the newest ask is kept: a slider's steps come
 * faster than a pot is run.
 */
import { PAGE_RUN, PAN_2, PAN_3, WAYS, coolKey, coolRest, coolWay, pansFor, restLitres, wayOf, type CoolRun, type WayRun } from "./coolDown";

export type CoolAsk = { litres: number; covered: boolean };
export type CoolDone = { key: string; runs: WayRun[] };

const scope = self as unknown as { onmessage: ((e: MessageEvent<CoolAsk>) => void) | null; postMessage(m: CoolDone): void };
const pans = new Map<string, WayRun>();
const rests = new Map<string, CoolRun | null>();
let ask: CoolAsk | null = null;
let busy = false;

function runs(a: CoolAsk): WayRun[] {
  return WAYS.map((w) => {
    if (w.key === "pot") return coolWay(w.key, a.litres, a.covered, PAGE_RUN);
    // The first pan's run: a full pan's, or the one pan's when the batch is less than a pan.
    const first = Math.min(a.litres, (w.key === "deep" ? PAN_3 : PAN_2).litres);
    const k = `${w.key}|${a.covered}`;
    const fk = `${k}|${first.toFixed(3)}`;
    let r = pans.get(fk);
    if (!r) pans.set(fk, (r = coolWay(w.key, first, a.covered, PAGE_RUN)));
    const rk = `${k}|${restLitres(w.key, a.litres).toFixed(3)}`;
    let rest = rests.get(rk);
    if (rest === undefined) rests.set(rk, (rest = coolRest(w.key, a.litres, a.covered, r.run.movedMin, PAGE_RUN)));
    return wayOf(w, r.run, pansFor(w.key, a.litres), rest);
  });
}

function pump(): void {
  if (busy) return;
  busy = true;
  // A turn of the event loop first, so a newer ask can take this one's place.
  setTimeout(() => {
    const a = ask;
    ask = null;
    if (a) scope.postMessage({ key: coolKey(a.litres, a.covered), runs: runs(a) });
    busy = false;
    if (ask) pump();
  }, 0);
}

scope.onmessage = (e) => {
  ask = e.data;
  pump();
};
