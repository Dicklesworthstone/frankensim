import { describe, expect, test } from "bun:test";
import { COLUMNS, FREEZER_PLACES, ITEM, LEVELS, ITEMS, MENU, SHELVES, THIS_WEEK, WALKIN_PLACES, casesOf, layoutWeek, manyWeeks, runWeeks, thisWeekCounts, type Policy, type Spot, type Trucks } from "./walkInWeek";

const POLICIES: Policy[] = [];
// Every setting the figure offers.
for (const frozen of [false, true]) for (const trucks of [2, 3, 5] as Trucks[]) for (let k = 0; k <= 5; k++) POLICIES.push({ trucks, frozen, cushion: k * 0.05 });
const isSpot = (e: unknown): e is Spot => typeof e === "object" && e !== null;

describe("a week of meat and fish through the walk-in", () => {
  test("the chart's keeping times, the day it comes in counted as day 1, and the Food Code's cooking temperatures", () => {
    expect(ITEMS.map((i) => [i.key, i.freshDays, i.thawedDays])).toEqual([
      ["salmon", 3, 2],
      ["ribs", 5, 5],
      ["beef", 2, 2],
      ["thighs", 2, 2],
    ]);
    // Top to bottom on the shelves: the hottest cooked lowest.
    for (let i = 1; i < ITEMS.length; i++) expect(ITEMS[i]!.cookF).toBeGreaterThanOrEqual(ITEMS[i - 1]!.cookF);
    expect(WALKIN_PLACES).toBe(SHELVES * COLUMNS * LEVELS);
    expect(THIS_WEEK.length).toBe(MENU.length);
  });

  test("a fresh case is used no later than its last day, and what is left after it is thrown out that night", () => {
    for (const p of POLICIES) {
      const w = runWeeks(p, thisWeekCounts());
      for (const ev of w.events) {
        for (const l of ev.stock) {
          expect(l.lb).toBeGreaterThan(0);
          // Nothing on hand past its last day once the night's check is done.
          if (ev.kind === "check") expect(l.lastDay).toBeGreaterThan(ev.day + 7);
        }
        if (ev.kind === "prep") for (const u of ev.lots) expect(w.events.flatMap((e) => e.stock).find((l) => l.id === u.id)?.lastDay ?? Infinity).toBeGreaterThanOrEqual(ev.day + 7);
      }
    }
  });

  test("two fresh trucks a week leave Wednesday's thighs short: Monday's are past their day by Tuesday night", () => {
    const w = runWeeks({ trucks: 2, frozen: false, cushion: 0.1 }, thisWeekCounts());
    const wed = w.events.find((e) => e.kind === "prep" && e.day === 2)!;
    expect(wed.shortLunches).toBe(THIS_WEEK[2]);
    expect(w.shortBy.thighs).toBe(THIS_WEEK[2]);
    // A truck every weekday knows each day's final count by its order: nothing short, and the cushion changes nothing.
    const five = [0, 0.1, 0.25].map((c) => runWeeks({ trucks: 5, frozen: false, cushion: c }, thisWeekCounts()));
    for (const r of five) expect(r.shortLunches).toBe(0);
    expect(new Set(five.map((r) => r.thrownLb)).size).toBe(1);
  });

  test("fresh, what a truck brings for its own day's lunch is that day's final count, rounded up to whole cases", () => {
    const w = runWeeks({ trucks: 5, frozen: false, cushion: 0 }, thisWeekCounts());
    for (const m of MENU) {
      const truck = w.events.find((e) => e.kind === "truck" && e.day === m.day)!;
      const lb = truck.lots.filter((l) => l.item === m.item).reduce((s, l) => s + l.lb, 0);
      const it = ITEM[m.item];
      // Less what is still good from before.
      expect(lb).toBeLessThanOrEqual(Math.ceil((THIS_WEEK[m.day]! * it.lbPer) / it.caseLb) * it.caseLb);
      expect(lb % it.caseLb).toBe(0);
    }
  });

  test("frozen, the thaw goes by the plan a day ahead: with no cushion a count over the plan comes up short", () => {
    const none = runWeeks({ trucks: 3, frozen: true, cushion: 0 }, thisWeekCounts());
    const some = runWeeks({ trucks: 3, frozen: true, cushion: 0.1 }, thisWeekCounts());
    expect(none.shortLunches).toBeGreaterThan(0);
    expect(some.shortLunches).toBe(0);
    expect(some.thrownLb).toBeGreaterThan(none.thrownLb);
    for (const ev of some.events) if (ev.kind === "thaw") expect(ev.hour % 24).toBe(8);
  });

  test("frozen on two trucks fills the freezer's shelves; on five they stay half empty", () => {
    const two = manyWeeks({ trucks: 2, frozen: true, cushion: 0.1 }, 100);
    const five = manyWeeks({ trucks: 5, frozen: true, cushion: 0.1 }, 100);
    expect(two.freezerPeak.mid).toBeGreaterThan(0.85 * FREEZER_PLACES);
    expect(five.freezerPeak.p90).toBeLessThan(0.6 * FREEZER_PLACES);
    expect(five.anyFloor.freezer).toBe(0);
  });

  test("the many weeks are the same weeks for every policy, and more cushion never leaves more lunches short", () => {
    const a = manyWeeks({ trucks: 3, frozen: false, cushion: 0.1 }, 120);
    const b = manyWeeks({ trucks: 3, frozen: false, cushion: 0.1 }, 120);
    expect(a).toEqual(b);
    for (const frozen of [false, true]) {
      let last = Infinity;
      for (const c of [0, 0.05, 0.1, 0.15, 0.2, 0.25]) {
        const m = manyWeeks({ trucks: 3, frozen, cushion: c }, 120);
        expect(m.shortLunches.mid).toBeLessThanOrEqual(last);
        last = m.shortLunches.mid;
      }
    }
  });

  test("every case has a place: chicken never above beef or fish, the floor only once the shelves are full", () => {
    const weeks: number[][] = [thisWeekCounts()];
    // And a spread of weeks either side of the plan.
    for (const f of [0.8, 0.9, 1.1, 1.2]) weeks.push(thisWeekCounts().map((c, d) => (d >= 7 && d < 14 ? Math.round(c * f) : c)));
    for (const p of POLICIES)
      for (const counts of weeks) {
        const w = runWeeks(p, counts);
        const L = layoutWeek(w);
        // On the floor only what more than the shelves hold.
        if (L.floorMost.walkin) expect(w.walkinPeak).toBeGreaterThan(SHELVES * COLUMNS);
        if (L.floorMost.freezer) expect(w.freezerPeak).toBeGreaterThan(SHELVES * COLUMNS);
        for (const st of L.events) {
          // The cases on hand are the lots' cases, each in a place of its own.
          const want = st.ev.stock.reduce((s, l) => s + casesOf(l), 0);
          expect(st.after.length).toBe(want);
          const key = (box: string, floor: boolean, tier: number, col: number, level: number) => [box, floor, tier, col, level].join("|");
          const keys = st.after.map((c) => key(c.at.box, c.at.floor, c.at.tier, c.at.col, c.at.level));
          expect(new Set(keys).size).toBe(keys.length);
          // A case on top stands on one of its own lot.
          for (const c of st.after)
            if (c.at.level === 1) expect(st.after.find((x) => key(x.at.box, x.at.floor, x.at.tier, x.at.col, x.at.level) === key(c.at.box, c.at.floor, c.at.tier, c.at.col, 0))?.lot).toBe(c.lot);
          // A case put on the floor found the shelves nearly full.
          for (const m of st.moves)
            if (isSpot(m.to) && m.to.floor) {
              const box = m.to.box;
              expect(st.after.filter((c) => c.at.box === box && !c.at.floor).length).toBeGreaterThan(0.75 * WALKIN_PLACES);
            }
          for (const box of ["walkin", "freezer"] as const) {
            const here = st.after.filter((c) => c.at.box === box);
            const shelved = here.filter((c) => !c.at.floor);
            for (const a of shelved) for (const b of shelved) if (a.at.col === b.at.col && a.at.tier > b.at.tier) expect(ITEM[a.item].cookF).toBeLessThanOrEqual(ITEM[b.item].cookF);
            // Each case in the box its lot is in.
            for (const c of here) expect(st.ev.stock.find((l) => l.id === c.lot)!.where).toBe(box);
          }
          // Each move's ends: trucks from the dock, prep and the check to the kitchen, a thaw from the freezer to the walk-in.
          for (const m of st.moves) {
            // A case moved down a shelf in its column, to make room.
            if (isSpot(m.from) && isSpot(m.to) && m.from.box === m.to.box) {
              expect(st.ev.kind === "truck" || st.ev.kind === "thaw").toBe(true);
              expect([m.to.col, m.to.tier, m.to.floor]).toEqual([m.from.col, m.from.tier - 1, false]);
              continue;
            }
            if (st.ev.kind === "truck") expect(m.from).toBe("dock");
            else if (st.ev.kind === "thaw") expect(isSpot(m.from) && m.from.box === "freezer" && isSpot(m.to) && m.to.box === "walkin").toBe(true);
            else expect(m.to).toBe("kitchen");
          }
        }
      }
  });
});
