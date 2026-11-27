/**
 * The combi-trip figure's runs, off the page's main thread: for a menu and an order, the morning's plan and what each
 * trip on the grid costs it, for both combis (sim/combiTrip.ts). Only the newest ask is kept: a slider's steps come
 * faster than a map is made.
 */
import { planMorning, type Morning } from "./cateringMorning";
import { tripBase, tripKey, tripMap, type TripMap, type TripSettings } from "./combiTrip";

export type TripAsk = { s: TripSettings };
export type TripDone = { key: string; plan: Morning; map: TripMap };

const scope = self as unknown as { onmessage: ((e: MessageEvent<TripAsk>) => void) | null; postMessage(m: TripDone): void };
let ask: TripSettings | null = null;
let busy = false;

function pump(): void {
  if (busy) return;
  busy = true;
  // A turn of the event loop first, so a newer ask can take this one's place.
  setTimeout(() => {
    const s = ask;
    ask = null;
    if (s) {
      const plan = planMorning(tripBase(s));
      scope.postMessage({ key: tripKey(s), plan, map: tripMap(plan) });
    }
    busy = false;
    if (ask) pump();
  }, 0);
}

scope.onmessage = (e) => {
  ask = e.data.s;
  pump();
};
