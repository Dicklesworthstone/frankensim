/**
 * The two-kitchen figure's plans, off the page's main thread: for an order, the hub alone and the best split between
 * the hub and the kitchen across town (bestPair); for a menu and crew, the most meals each gets to lunch (pairLimit).
 * Only the newest ask of each kind is answered, a plan before the rest of a limit's search (which goes a split at a
 * time); plans are kept between asks, so a split seen before is not planned again.
 */
import { bestPair, limitKey, pairLimitSteps, planKey, type PairLimit, type PairPlan, type PairSettings, type PlanCache } from "./twoKitchens";
import type { Crew } from "./cateringMorning";

export type PairAsk = { kind: "plan"; s: PairSettings } | { kind: "limit"; menu: number; hubCooks: Crew; farCooks: Crew };
export type PairDone =
  | { kind: "plan"; key: string; alone: PairPlan; both: PairPlan }
  | { kind: "limit"; key: string; alone: number; both: number; parts: number[] };

const scope = self as unknown as { onmessage: ((e: MessageEvent<PairAsk>) => void) | null; postMessage(m: PairDone): void };
const cache: PlanCache = new Map();
let plan: Extract<PairAsk, { kind: "plan" }> | null = null;
/** The newest limit asked for, searched a split at a time between plans. */
let limit: { key: string; steps: Generator<void, PairLimit, void> } | null = null;
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
      const job = limit;
      const r = job.steps.next();
      if (r.done) {
        limit = null;
        scope.postMessage({ kind: "limit", key: job.key, alone: r.value.alone, both: r.value.meals, parts: r.value.parts });
      }
    }
    busy = false;
    if (plan || limit) pump();
  }, 0);
}

scope.onmessage = (e) => {
  const a = e.data;
  if (a.kind === "plan") plan = a;
  else {
    const key = limitKey(a.menu, a.hubCooks, a.farCooks);
    if (limit?.key !== key) limit = { key, steps: pairLimitSteps({ menu: a.menu, hubCooks: a.hubCooks, farCooks: a.farCooks }, 4000, cache) };
  }
  pump();
};
