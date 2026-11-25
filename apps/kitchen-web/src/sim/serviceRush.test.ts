import { describe, expect, test } from "bun:test";
import { KEEPERS, RUSH_MENU, SERVICE, modesFor, runRush, rushMenu, simulateService, type RushSettings } from "./serviceRush";

const toOrder: RushSettings = { modes: modesFor([]), day: "busy", sigma: 0 };
const held: RushSettings = { modes: modesFor([...KEEPERS]), day: "busy", sigma: 0 };

describe("one service", () => {
  test("every order is plated once, after it came in, and no food appears or vanishes", () => {
    for (const s of [toOrder, held, { ...held, pool: true, menu: "balanced" as const, sigma: 0.3 }])
      for (const seed of [1, 2, 3]) {
        const d = simulateService(s, seed);
        const menu = rushMenu(s.menu ?? "written");
        expect(d.orders.length).toBeGreaterThan(100);
        const served = menu.map(() => 0);
        for (const o of d.orders) {
          expect(Number.isFinite(o.plated)).toBe(true);
          expect(o.food).toBeGreaterThanOrEqual(o.arrive - 1e-9);
          expect(o.plated).toBeGreaterThanOrEqual(o.food + SERVICE.handMin - 1e-9);
          expect(o.plated).toBeLessThan(SERVICE.end);
          served[o.dish]!++;
        }
        // Portions cooked = portions served + portions thrown out, dish by dish.
        menu.forEach((_, i) => expect(d.cooked[i]).toBe(served[i]! + d.wasted[i]!));
      }
  });

  test("a cook is never busy more than the minute, and idle time is a share of the crew's", () => {
    const d = simulateService({ ...held, pool: true }, 5);
    for (const b of d.busy) for (const v of b) expect(v).toBeLessThanOrEqual(1 + 1e-6);
    expect(d.idle).toBeGreaterThan(0);
    expect(d.idle).toBeLessThan(1);
  });

  test("nothing is thrown out when nothing is held, and the same seed gives the same day", () => {
    const a = simulateService(toOrder, 9);
    expect(a.wasteTotal).toBe(0);
    const b = simulateService(toOrder, 9);
    expect(b.waits).toEqual(a.waits);
  });

  test("a held dish waits in the cabinet: its guests wait less than for one cooked to order", () => {
    const beef = RUSH_MENU.findIndex((d) => d.key === "beef");
    const p90 = (s: RushSettings) => {
      const d = simulateService(s, 4);
      const w = d.orders.filter((o) => o.dish === beef).map((o) => o.plated - o.arrive).sort((x, y) => x - y);
      return w[Math.floor(w.length * 0.9)]!;
    };
    expect(p90(held)).toBeLessThan(p90(toOrder) / 3);
  });
});

describe("the levers, across many days", () => {
  const order = runRush(toOrder, 80);
  const hold = runRush(held, 80);
  const pooled = runRush({ ...held, pool: true }, 80);

  test("holding what keeps moves work into the quiet hour before the doors open", () => {
    const early = (s: RushSettings) => simulateService(s, 3).busy.reduce((a, b) => a + b.slice(0, SERVICE.open).reduce((x, y) => x + y, 0), 0);
    expect(early(toOrder)).toBe(0);
    expect(early(held)).toBeGreaterThan(30);
  });

  test("on the busy day, holding and pooled plating cut the long waits; holding costs some food", () => {
    expect(hold.waitP90.p50).toBeLessThan(order.waitP90.p50);
    expect(pooled.waitP90.p50).toBeLessThan(order.waitP90.p50 / 2);
    expect(order.waste.p50).toBe(0);
    expect(hold.waste.p50).toBeGreaterThan(0);
    expect(hold.waste.p50).toBeLessThan(0.1 * 420);
  });

  test("moving the thighs to the combi takes load off the griddle in the rush", () => {
    const peak = (s: RushSettings) => {
      let sum = 0;
      for (let seed = 1; seed <= 20; seed++) {
        const b = simulateService(s, seed).busy[1]!;
        for (let m = 120; m < 180; m++) sum += b[m]!;
      }
      return sum;
    };
    expect(peak({ ...held, menu: "balanced" })).toBeLessThan(peak(held));
  });

  test("a less certain forecast throws out more held food", () => {
    const sure = runRush(held, 80);
    const unsure = runRush({ ...held, sigma: 0.3 }, 80);
    expect(unsure.waste.p90).toBeGreaterThan(sure.waste.p90);
  });
});
