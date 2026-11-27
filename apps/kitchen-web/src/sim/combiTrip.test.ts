import { describe, expect, test } from "bun:test";
import { OVENS, carryOn, menuDishes, planMorning, type Morning } from "./cateringMorning";
import { AT_N, FOR_N, TRIP_AT, TRIP_FOR, cellOf, free, tripBase, tripMap, tripMorning, type Trip } from "./combiTrip";

/** Menus and orders around their limits with three cooks, the all-combi menu (7) on both sides of its step. */
const PLANS: Morning[] = (
  [
    [7, 800],
    [7, 820],
    [7, 1200],
    [7, 1400],
    [6, 320],
    [3, 200],
    [5, 170],
    [0, 80],
  ] as const
).map(([menu, meals]) => planMorning(tripBase({ menu, meals })));
const TRIPS: Trip[] = [];
for (let oven = 0; oven < OVENS; oven++) for (let at = TRIP_AT.from; at <= TRIP_AT.to; at += 40) for (let dur = TRIP_FOR.from; dur <= TRIP_FOR.to; dur += 40) TRIPS.push({ oven, at, mins: dur });
const overlap = (a0: number, a1: number, b0: number, b1: number) => Math.max(0, Math.min(a1, b1) - Math.max(a0, b0));
/** A morning's record without its settings, for comparing two runs. */
const record = (m: Morning) => JSON.stringify({ ...m, settings: null });

describe("combiTrip", () => {
  test("with no combi out, carrying on is the plan itself, job for job", () => {
    for (const plan of PLANS) expect(record(carryOn(plan.settings, plan))).toBe(record(plan));
  });

  test("a combi that is out takes nothing in, and holds up only the program running in it", () => {
    for (const plan of PLANS) {
      const dishes = menuDishes(plan.settings.menu);
      for (const trip of TRIPS) {
        const m = tripMorning(plan, trip);
        const t1 = trip.at + trip.mins;
        for (const l of m.loads) {
          if (l.oven === trip.oven) expect(l.t0 >= trip.at - 1e-6 && l.t0 < t1 - 1e-6).toBe(false);
          // The program's own minutes, plus the minutes it stood stopped in the combi that was out.
          const st = dishes[l.slot]!.stages.find((x) => x.at === "oven" && x.verb === l.verb)!;
          const stopped = l.oven === trip.oven ? overlap(l.t1, l.t2, trip.at, t1) : 0;
          expect(l.t2 - l.t1 - stopped).toBeCloseTo(st.min, 6);
          expect(l.tu).toBeGreaterThanOrEqual(l.t2 - 1e-9);
        }
      }
    }
  });

  test("up to the trip the morning is the plan's, and every dish keeps the plan's hold on it", () => {
    for (const plan of PLANS)
      for (const trip of TRIPS) {
        const m = tripMorning(plan, trip);
        const before = (x: { t0: number }) => x.t0 < trip.at - 1e-6;
        expect(JSON.stringify(m.jobs.filter(before))).toBe(JSON.stringify(plan.jobs.filter(before)));
        const start = ({ oven, batches, slot, verb, pans, t0, t1 }: Morning["loads"][number]) => ({ oven, batches, slot, verb, pans, t0, t1 });
        expect(m.loads.filter(before).map(start)).toEqual(plan.loads.filter(before).map(start));
        expect(m.batches.map((b) => b.release)).toEqual(plan.batches.map((b) => b.release));
      }
  });

  test("each map cell is its trip's morning: the van's minutes late and the portions past their best, against the plan", () => {
    const plan = PLANS[0]!;
    const map = tripMap(plan);
    expect(map.late.length).toBe(OVENS);
    for (const trip of TRIPS) {
      const k = cellOf(trip.at, trip.mins);
      expect(k).toBeGreaterThanOrEqual(0);
      expect(k).toBeLessThan(AT_N * FOR_N);
      const m = tripMorning(plan, trip);
      expect(map.late[trip.oven]![k]!).toBeCloseTo(Math.max(0, m.depart - plan.depart), 4);
      expect(map.pastBest[trip.oven]![k]!).toBe(Math.max(0, m.pastBest - plan.pastBest));
    }
  });

  test("the all-combi menu: at 800 meals nine trips in ten cost nothing; at 820 the braise takes both combis and four in five cost the van", () => {
    const share = (plan: Morning) => {
      const map = tripMap(plan);
      let n = 0;
      for (let o = 0; o < OVENS; o++) for (let k = 0; k < AT_N * FOR_N; k++) if (free({ late: map.late[o]![k]!, pastBest: map.pastBest[o]![k]! })) n++;
      return n / (OVENS * AT_N * FOR_N);
    };
    const [p800, p820] = [PLANS[0]!, PLANS[1]!];
    expect(share(p800)).toBeGreaterThan(0.85);
    expect(share(p820)).toBeLessThan(0.2);
    // At 800 the braise is one load, in one combi; at 820 it is two, one in each, both in by 6:30.
    const braise = (p: Morning) => p.loads.filter((l) => l.verb === "Braise");
    expect(braise(p800).length).toBe(1);
    expect(new Set(braise(p820).map((l) => l.oven)).size).toBe(2);
    for (const l of braise(p820)) expect(l.t0).toBeLessThan(30);
  });
});
