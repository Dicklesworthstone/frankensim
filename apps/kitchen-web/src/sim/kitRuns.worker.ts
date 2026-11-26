/**
 * The buy-one-thing figure's plans, off the page's main thread: for an order, the morning with the equipment and crew
 * bought and the morning as built (planMorning); for a menu, the most meals each choice of purchases gets onto the van
 * with every portion at its best (morningLimit), posted one at a time as they come. A plan is answered before the
 * next limit, and only the newest ask of each kind is kept.
 */
import { morningLimit, planMorning, type Morning } from "./cateringMorning";
import { BUYS, limitKey, morningOf, planKey, type Buy, type KitSettings } from "./buyOne";

export type KitAsk = { kind: "plan"; s: KitSettings } | { kind: "limits"; menu: number };
export type KitDone = { kind: "plan"; key: string; m: Morning; built: Morning } | { kind: "limit"; key: string; meals: number };

const scope = self as unknown as { onmessage: ((e: MessageEvent<KitAsk>) => void) | null; postMessage(m: KitDone): void };
let plan: KitSettings | null = null;
/** The newest menu's limits still to find. */
let limits: { menu: number; left: Buy[] } | null = null;
const found = new Set<string>();
let busy = false;

function pump(): void {
  if (busy) return;
  busy = true;
  // A turn of the event loop first, so a newer ask can take this one's place.
  setTimeout(() => {
    if (plan) {
      const s = plan;
      plan = null;
      const m = planMorning(morningOf(s.meals, s.menu, s));
      const built = s.combi || s.range || s.cook ? planMorning(morningOf(s.meals, s.menu, BUYS[0]!)) : m;
      scope.postMessage({ kind: "plan", key: planKey(s), m, built });
    } else if (limits) {
      const { menu } = limits;
      const b = limits.left.shift()!;
      const key = limitKey(menu, b);
      if (!limits.left.length) limits = null;
      if (!found.has(key)) {
        const { meals: _, ...s } = morningOf(0, menu, b);
        found.add(key);
        scope.postMessage({ kind: "limit", key, meals: morningLimit(s) });
      }
    }
    busy = false;
    if (plan || limits) pump();
  }, 0);
}

scope.onmessage = (e) => {
  const a = e.data;
  if (a.kind === "plan") plan = a.s;
  else if (limits?.menu !== a.menu) limits = { menu: a.menu, left: [...BUYS] };
  pump();
};
