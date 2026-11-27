/**
 * The carriers figure's runs, off the page's main thread (sim/holdHot.ts): the four carriers for chili packed at one
 * temperature. Only two are ever asked for, so each is kept; the newest ask goes first.
 */
import { fToC } from "./coolDown";
import { holdAll, type HoldRun } from "./holdHot";

export type HoldAsk = { startF: number };
export type HoldDone = { startF: number; runs: HoldRun[] };

const scope = self as unknown as { onmessage: ((e: MessageEvent<HoldAsk>) => void) | null; postMessage(m: HoldDone): void };
const made = new Map<number, HoldRun[]>();
const asks: number[] = [];
let busy = false;

function pump(): void {
  if (busy) return;
  busy = true;
  setTimeout(() => {
    const f = asks.shift();
    if (f !== undefined) {
      let r = made.get(f);
      if (!r) made.set(f, (r = holdAll(fToC(f))));
      scope.postMessage({ startF: f, runs: r });
    }
    busy = false;
    if (asks.length) pump();
  }, 0);
}

scope.onmessage = (e) => {
  const f = e.data.startF;
  const i = asks.indexOf(f);
  if (i >= 0) asks.splice(i, 1);
  asks.unshift(f);
  pump();
};
