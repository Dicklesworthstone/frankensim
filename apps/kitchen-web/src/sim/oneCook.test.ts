import { describe, expect, test } from "bun:test";
import { DISHES, menuDishes, planMorning, type Morning } from "./cateringMorning";
import { COOK_DISHES, cardOf, endOf, oneCook } from "./oneCook";

const ORDERS = [20, 60, 120, 160] as const;
const PAIRS: [number, number][] = [];
for (let i = 0; i < COOK_DISHES.length; i++) for (let j = 0; j < COOK_DISHES.length; j++) if (i !== j) PAIRS.push([i, j]);
const RUNS = ORDERS.flatMap((p) => PAIRS.map(([i, j]) => oneCook(COOK_DISHES[i]!, COOK_DISHES[j]!, p)));

/** The most of one cook's hands any moment of a plan asks for. */
const peakHands = (m: Morning) => {
  let peak = 0;
  for (const j of m.jobs) {
    let h = 0;
    for (const k of m.jobs) if (k.t0 <= j.t0 + 1e-9 && k.t1 > j.t0 + 1e-9) h += k.hands;
    peak = Math.max(peak, h);
  }
  return peak;
};

describe("one cook, two dishes", () => {
  test("every plan finishes, with one cook", () => {
    for (const r of RUNS) {
      for (const m of [r.together, r.apart]) {
        expect(m.finished).toBe(true);
        expect(m.jobs.every((j) => j.cook === 0)).toBe(true);
      }
    }
  });

  test("together is never longer than one after the other, nor shorter than either dish alone", () => {
    for (const r of RUNS) {
      expect(r.end).toBeLessThanOrEqual(r.endApart + 1e-6);
      expect(r.end).toBeGreaterThanOrEqual(Math.max(r.endA, r.endB) - 1e-6);
      expect(endOf(r.apart)).toBeCloseTo(r.endApart, 6);
    }
  });

  test("the same work either way: hands-on minutes, batches and pans", () => {
    for (const r of RUNS) {
      expect(r.together.handsMin).toBeCloseTo(r.apart.handsMin, 6);
      expect(r.together.batches.length).toBe(r.apart.batches.length);
      for (const slot of [0, 1]) {
        const p = (m: Morning) => m.batches.filter((b) => b.slot === slot).reduce((a, b) => a + b.portions, 0);
        expect(p(r.together)).toBe(r.portions);
        expect(p(r.apart)).toBe(r.portions);
      }
      expect(r.together.pans).toBe(r.apart.pans);
      expect(r.free).toBeCloseTo(r.end - r.together.handsMin, 9);
    }
  });

  test("the cook never needs more than one pair of hands", () => {
    for (const r of RUNS) {
      expect(peakHands(r.together)).toBeLessThanOrEqual(1 + 1e-9);
      expect(peakHands(r.apart)).toBeLessThanOrEqual(1 + 1e-9);
    }
  });

  test("one after the other: the second dish starts once the first is done", () => {
    for (const r of RUNS) {
      const first = r.apart.batches.filter((b) => b.slot === 0).reduce((a, b) => Math.max(a, b.packed), 0);
      expect(first).toBeCloseTo(r.endA, 6);
      for (const j of r.apart.jobs) if (j.slot === 1) expect(j.t0).toBeGreaterThanOrEqual(r.endA - 1e-6);
      for (const l of r.apart.loads) if (l.slot === 1) expect(l.t0).toBeGreaterThanOrEqual(r.endA - 1e-6);
      // Each job's batch is one of its own dish's.
      for (const j of r.apart.jobs) expect(r.apart.batches[j.batch]!.slot).toBe(j.slot);
    }
  });

  test("two dishes that need the cook's hands at the range the whole time they cook save nothing together", () => {
    for (const p of ORDERS)
      for (const [a, b] of [
        [DISHES.risotto, DISHES.piccata],
        [DISHES.risotto, DISHES.tacos],
      ] as const) {
        const r = oneCook(a, b, p);
        expect(r.end).toBeCloseTo(r.endApart, 6);
      }
  });

  test("a dish that roasts on its own takes the other's hands-on work inside it", () => {
    for (const p of ORDERS) {
      const r = oneCook(DISHES.thighs, DISHES.risotto, p);
      // As much as the shorter of the risotto's work and the 35-minute roast, less a few minutes at the oven door.
      expect(r.endApart - r.end).toBeGreaterThan(Math.min(r.endB, DISHES.thighs.stages[1].min) - 5);
    }
  });

  test("cards: the table needs the hands all the time, the range its share, a combi only to load and unload", () => {
    for (const d of COOK_DISHES)
      d.stages.forEach((s, k) => {
        const c = cardOf(d)[k]!;
        expect(c.verb).toBe(s.verb);
        if (s.at === "table") expect(c.hands).toBe(1);
        if (s.at === "range") expect(c.hands).toBe(s.hands ?? 1);
        if (s.at === "oven") {
          expect(c.hands).toBeGreaterThan(0);
          expect(c.hands).toBeLessThan(0.15);
        }
      });
  });

  test("a menu's dishes listed in the settings plan the same morning as the menu", () => {
    for (const menu of [0, 3, 7])
      for (const ahead of [false, true]) {
        const s = { meals: 700, cooks: 3, menu, ahead } as const;
        const plain = planMorning(s);
        const listed = planMorning({ ...s, dishes: menuDishes(menu) });
        expect(JSON.stringify({ ...listed, settings: null })).toBe(JSON.stringify({ ...plain, settings: null }));
      }
  });
});
