/**
 * The two-kitchen figure's plans, off the page's main thread: for an order, the hub alone and the best split between
 * the hub and the kitchen across town (bestPair); for a menu and crew, the most meals each gets to lunch (pairLimit).
 * Only the newest ask of each kind is answered; plans are kept between asks, so a split seen before is not planned
 * again.
 */
import { bestPair, limitKey, pairLimit, planKey, type PairPlan, type PairSettings, type PlanCache } from "./twoKitchens";
import type { Crew } from "./cateringMorning";

export type PairAsk = { kind: "plan"; s: PairSettings } | { kind: "limit"; menu: number; hubCooks: Crew; farCooks: Crew };
export type PairDone =
  | { kind: "plan"; key: string; alone: PairPlan; both: PairPlan }
  | { kind: "limit"; key: string; alone: number; both: number; parts: number[] };

const scope = self as unknown as { onmessage: ((e: MessageEvent<PairAsk>) => void) | null; postMessage(m: PairDone): void };
const cache: PlanCache = new Map();
let plan: Extract<PairAsk, { kind: "plan" }> | null = null;
let limit: Extract<PairAsk, { kind: "limit" }> | null = null;
let busy = false;

function pump(): void {
  if (busy) return;
  busy = true;
  // A turn of the event loop first, so a newer ask can take this one's place.
  setTimeout(() => {
    // A big plan cache is dropped rather than kept growing (each plan holds the morning's every job).
    if (cache.size > 3000) cache.clear();
    if (plan) {
      const { s } = plan;
      plan = null;
      const alone = bestPair(s, true, cache);
      const both = bestPair(s, false, cache);
      scope.postMessage({ kind: "plan", key: planKey(s), alone, both });
    } else if (limit) {
      const { menu, hubCooks, farCooks } = limit;
      limit = null;
      const r = pairLimit({ menu, hubCooks, farCooks }, 4000, cache);
      scope.postMessage({ kind: "limit", key: limitKey(menu, hubCooks, farCooks), alone: r.alone, both: r.meals, parts: r.parts });
    }
    busy = false;
    if (plan || limit) pump();
  }, 0);
}

scope.onmessage = (e) => {
  if (e.data.kind === "plan") plan = e.data;
  else limit = e.data;
  pump();
};
