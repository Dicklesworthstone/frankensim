import { describe, expect, test } from "bun:test";
import { DOOR, MENUS, SERVE, madeIt, menuDishes, morningLimit, planMorning, portionsFor, type Crew } from "./cateringMorning";
import { FAR_DOOR, bestPair, pairLimit, planPair, splitPortions, splitsFor, travels, type PairSettings, type PlanCache } from "./twoKitchens";

const cache: PlanCache = new Map();
const ORDERS: PairSettings[] = [];
for (const meals of [400, 1000, 1600]) for (let menu = 0; menu < MENUS; menu++) for (const farCooks of [2, 3] as Crew[]) ORDERS.push({ meals, menu, hubCooks: 3, farCooks });

describe("twoKitchens", () => {
  test("a split hands every portion to one kitchen or the other, the odd one to the hub", () => {
    for (const meals of [50, 333, 1400, 2399])
      for (const parts of splitsFor(7)) {
        const { here, there } = splitPortions(meals, parts);
        portionsFor(meals).forEach((p, k) => {
          expect(here[k]! + there[k]!).toBe(p);
          expect(there[k]).toBe(Math.floor(p * parts[k]!));
          expect(here[k]!).toBeGreaterThanOrEqual(0);
        });
      }
  });

  test("only a dish that keeps past the far van's drive is ever split, and every split is tried once", () => {
    for (let menu = 0; menu < MENUS; menu++) {
      const dishes = menuDishes(menu);
      const splits = splitsFor(menu);
      expect(splits[0]!.every((f) => f === 0)).toBe(true);
      expect(splits.length).toBe(3 ** dishes.filter((d) => travels(d.holdMin)).length);
      expect(new Set(splits.map((s) => s.join())).size).toBe(splits.length);
      for (const s of splits) dishes.forEach((d, k) => !travels(d.holdMin) && expect(s[k]).toBe(0));
    }
    // The line: a dish must keep longer than the hour and a half from 10:00 to lunch.
    expect(travels(SERVE - FAR_DOOR)).toBe(false);
    expect(travels(SERVE - FAR_DOOR + 1)).toBe(true);
  });

  test("the hub alone is the morning figure's own plan", () => {
    for (const s of ORDERS.filter((o) => o.farCooks === 2)) {
      const p = bestPair(s, true, cache);
      const m = planMorning({ meals: s.meals, cooks: s.hubCooks, menu: s.menu, ahead: false });
      expect(p.far).toBeNull();
      // The same plan; its settings record only spells out the portions and van the morning takes by default.
      expect(JSON.stringify({ ...p.hub, settings: null })).toBe(JSON.stringify({ ...m, settings: null }));
      expect(p.hub.settings).toEqual({ ...m.settings, portions: portionsFor(s.meals), door: DOOR });
      expect(p.lunchLate).toBe(m.late);
      expect(p.pastBest).toBe(m.pastBest);
      expect(p.made).toBe(madeIt(m));
    }
  });

  test("each kitchen cooks and packs exactly its portions, and neither van leaves before its time", () => {
    for (const s of ORDERS) {
      const p = bestPair(s, false, cache);
      for (const [m, want, door] of [[p.hub, p.here, DOOR], ...(p.far ? [[p.far, p.there, FAR_DOOR] as const] : [])] as const) {
        const got = want.map(() => 0);
        for (const b of m.batches) got[b.slot]! += b.portions;
        expect(got).toEqual(want);
        expect(m.depart).toBeGreaterThanOrEqual(door);
        if (m.finished) expect(m.late).toBeCloseTo(m.depart - door, 9);
      }
      expect(p.far === null).toBe(p.there.every((x) => x === 0));
    }
  });

  test("lunch waits for the later van, and a portion is past its best by the time it is actually served", () => {
    for (const s of ORDERS) {
      const p = bestPair(s, false, cache);
      const lunch = SERVE + p.lunchLate;
      expect(p.lunchLate).toBe(Math.max(p.hub.late, p.far?.late ?? 0));
      const dishes = menuDishes(s.menu);
      const by = dishes.map(() => 0);
      for (const m of p.far ? [p.hub, p.far] : [p.hub]) for (const b of m.batches) if (b.cooked >= 0 && lunch - b.cooked > dishes[b.slot]!.holdMin + 1e-6) by[b.slot]! += b.portions;
      expect(p.pastBestBySlot).toEqual(by);
      expect(p.made).toBe(madeIt(p.hub) && (!p.far || madeIt(p.far)) && p.lunchLate <= 0.5 && p.pastBest === 0);
    }
  });

  test("the best split is at least as good as the hub alone, and as every other split", () => {
    for (const s of ORDERS.filter((o) => o.meals === 1000)) {
      const best = bestPair(s, false, cache);
      const alone = bestPair(s, true, cache);
      if (alone.made) expect(best.made).toBe(true);
      for (const parts of splitsFor(s.menu)) {
        const p = planPair(s, parts, cache);
        if (p.made) expect(best.made).toBe(true);
        if (p.made === best.made) expect(best.lunchLate).toBeLessThanOrEqual(p.lunchLate + 0.5);
      }
    }
  });

  test("the order on the page: the hub alone is late, both kitchens get lunch out on time at its best", () => {
    const s: PairSettings = { meals: 1400, menu: MENUS - 1, hubCooks: 3, farCooks: 2 };
    const alone = bestPair(s, true, cache);
    const both = bestPair(s, false, cache);
    expect(alone.made).toBe(false);
    expect(alone.lunchLate).toBeGreaterThan(15);
    expect(both.made).toBe(true);
    expect(both.far).not.toBeNull();
  });

  test("the most meals: never fewer than the hub alone; more only when every dish keeps", () => {
    const most = new Map<Crew, number>();
    for (const farCooks of [2, 3] as Crew[])
      for (let menu = 0; menu < MENUS; menu++) {
        const c: PlanCache = new Map();
        const r = pairLimit({ menu, hubCooks: 3, farCooks }, 4000, c);
        expect(r.alone).toBe(morningLimit({ cooks: 3, menu, ahead: false }));
        expect(r.meals).toBeGreaterThanOrEqual(r.alone);
        expect(r.meals % 10).toBe(0);
        if (r.meals > r.alone) expect(planPair({ meals: r.meals, menu, hubCooks: 3, farCooks }, r.parts, c).made).toBe(true);
        const allKeep = menuDishes(menu).every((d) => travels(d.holdMin));
        if (allKeep) expect(r.meals).toBeGreaterThan(r.alone);
        else expect(r.meals).toBe(r.alone);
        if (menu === MENUS - 1) most.set(farCooks, r.meals);
      }
    // A third cook across town carries more.
    expect(most.get(3)!).toBeGreaterThan(most.get(2)!);
  }, 60_000);

  test("deterministic, with or without the plan cache", () => {
    const s: PairSettings = { meals: 1800, menu: MENUS - 1, hubCooks: 3, farCooks: 2 };
    const a = JSON.stringify(bestPair(s, false, new Map()));
    expect(JSON.stringify(bestPair(s, false, null))).toBe(a);
    expect(JSON.stringify(bestPair(s, false, cache))).toBe(a);
  });
});
