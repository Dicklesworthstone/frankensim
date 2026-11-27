/**
 * The cooling figure's runs, off the page's main thread (sim/coolDown.ts): the four ways for a batch, open or
 * covered. The pans cool the same however many there are, so their runs are kept and only the pot runs again for a
 * new batch. Only the newest ask is kept: a slider's steps come faster than a pot is run.
 */
import { DAY_MIN, WAYS, coolKey, coolWay, pansFor, type WayRun } from "./coolDown";

export type CoolAsk = { litres: number; covered: boolean };
export type CoolDone = { key: string; runs: WayRun[] };

const scope = self as unknown as { onmessage: ((e: MessageEvent<CoolAsk>) => void) | null; postMessage(m: CoolDone): void };
const pans = new Map<string, WayRun>();
let ask: CoolAsk | null = null;
let busy = false;

function runs(a: CoolAsk): WayRun[] {
  return WAYS.map((w) => {
    const o = { minMin: DAY_MIN, maxMin: 72 * 60 };
    if (w.key === "pot") return coolWay(w.key, a.litres, a.covered, o);
    const k = `${w.key}|${a.covered}`;
    let r = pans.get(k);
    if (!r) pans.set(k, (r = coolWay(w.key, a.litres, a.covered, o)));
    return { ...r, pans: pansFor(w.key, a.litres) };
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
