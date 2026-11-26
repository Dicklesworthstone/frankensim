import { describe, expect, test } from "bun:test";
import { KIT, MENUS, RANGE_BURNERS, RANGE_PLACES, TABLE_SPOTS, kitOf, morningLimit, planMorning, type Morning } from "./cateringMorning";
import { BUYS, buyKey, morningOf } from "./buyOne";

const ORDERS = [200, 700, 1400];

/** Whether the jobs and burns of a morning keep to its kit at every moment. */
function keepsToKit(m: Morning): void {
  const kit = kitOf(m.settings);
  const ts = new Set<number>();
  for (const j of m.jobs) ts.add(j.t0 + 1e-4);
  for (const l of m.loads) ts.add(l.t0 + 1e-4);
  for (const t of ts) {
    const on = (a: number, b: number) => a <= t && b > t;
    const range = m.jobs.filter((j) => j.at === "range" && on(j.t0, j.t1));
    // Each half of each range has at most one cook, and a cook's pans are on the range of the half they stand at.
    for (let side = 0; side < kit.ranges * RANGE_PLACES; side++) expect(new Set(range.filter((j) => j.side === side).map((j) => j.cook)).size).toBeLessThanOrEqual(1);
    for (const j of range) expect(Math.floor(j.spot / RANGE_BURNERS)).toBe(Math.floor(j.side! / RANGE_PLACES));
    const burners = m.burns.filter((b) => on(b.t0, b.t1)).map((b) => b.burner);
    expect(new Set(burners).size).toBe(burners.length);
    for (const b of burners) expect(b).toBeLessThan(kit.ranges * RANGE_BURNERS);
    const ovens = m.loads.filter((l) => on(l.t0, l.t3)).map((l) => l.oven);
    expect(new Set(ovens).size).toBe(ovens.length);
    for (const o of ovens) expect(o).toBeLessThan(kit.ovens);
    expect(m.jobs.filter((j) => j.at === "table" && on(j.t0, j.t1)).length).toBeLessThanOrEqual(TABLE_SPOTS);
  }
}

describe("buyOne", () => {
  test("the kitchen as built is the morning figure's own, plan for plan", () => {
    for (const meals of ORDERS)
      for (let menu = 0; menu < MENUS; menu++) {
        const a = planMorning({ meals, menu, cooks: 3, ahead: false });
        const b = planMorning(morningOf(meals, menu, BUYS[0]!));
        expect(JSON.stringify({ ...b, settings: null })).toBe(JSON.stringify({ ...a, settings: null }));
        expect(kitOf(b.settings)).toEqual(KIT);
      }
  });

  test("every purchase keeps to its kitchen: a cook to a range half, pans on that range, one load to a combi", () => {
    for (const b of BUYS) for (const meals of ORDERS) for (const menu of [0, 6, 7]) keepsToKit(planMorning(morningOf(meals, menu, b)));
  });

  test("a cook waits for room only with free hands, while work is ready for a station that is full", () => {
    for (const b of BUYS)
      for (const menu of [0, 6, 7]) {
        const m = planMorning(morningOf(1000, menu, b));
        for (const w of m.waits) {
          expect(w.t1).toBeGreaterThan(w.t0);
          const mid = (w.t0 + w.t1) / 2;
          // Nobody waiting has a job at the time.
          for (const c of w.cooks) expect(m.jobs.some((j) => j.cook === c && j.t0 <= mid && j.t1 > mid)).toBe(false);
          // And the station waited for is full then.
          const kit = kitOf(m.settings);
          if (w.at === "oven") expect(m.loads.filter((l) => l.t0 <= mid && l.t3 > mid).length).toBe(kit.ovens);
          if (w.at === "table") expect(m.jobs.filter((j) => j.at === "table" && j.t0 <= mid && j.t1 > mid).length).toBe(TABLE_SPOTS);
          if (w.at === "range")
            for (let r = 0; r < kit.ranges; r++) {
              const halves = new Set(m.jobs.filter((j) => j.at === "range" && j.t0 <= mid && j.t1 > mid && Math.floor(j.side! / RANGE_PLACES) === r).map((j) => j.side)).size;
              const lit = m.burns.filter((x) => x.t0 <= mid && x.t1 > mid && Math.floor(x.burner / RANGE_BURNERS) === r).length;
              expect(halves === RANGE_PLACES || lit === RANGE_BURNERS).toBe(true);
            }
        }
      }
  });

  test("what each purchase buys: the combis on the oven menu, the range on the rest, a cook only with room for one", () => {
    const lim = (menu: number, k: string) => {
      const b = BUYS.find((x) => buyKey(x) === k)!;
      const { meals: _, ...s } = morningOf(0, menu, b);
      return morningLimit(s);
    };
    // Every dish in the combis: a third combi raises the most meals, a second range does not.
    const ovens = MENUS - 1;
    expect(lim(ovens, "100")).toBeGreaterThan(lim(ovens, "000"));
    expect(lim(ovens, "010")).toBe(lim(ovens, "000"));
    // Piccata, mac and braise: the range is full. A second range raises it, a combi does not, and a fourth cook adds
    // nothing alone but more with the range.
    const mixed = 6;
    const base = lim(mixed, "000");
    expect(lim(mixed, "010")).toBeGreaterThan(base);
    expect(lim(mixed, "100")).toBe(base);
    expect(lim(mixed, "001")).toBe(base);
    expect(lim(mixed, "011")).toBeGreaterThan(lim(mixed, "010"));
  }, 60_000);

  test("buying more never lowers the most meals", () => {
    for (const menu of [0, 3, 6, 7]) {
      const v = new Map(BUYS.map((b) => {
        const { meals: _, ...s } = morningOf(0, menu, b);
        return [buyKey(b), morningLimit(s)] as const;
      }));
      for (const a of BUYS)
        for (const b of BUYS)
          if ((!a.combi || b.combi) && (!a.range || b.range) && (!a.cook || b.cook)) expect(v.get(buyKey(b))!).toBeGreaterThanOrEqual(v.get(buyKey(a))!);
    }
  }, 60_000);
});
