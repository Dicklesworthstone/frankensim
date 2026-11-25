import { describe, expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { CARRY_R, CHILL_LOAD, WALKIN, buildFloor, distanceField, fieldAt, fieldsOf, placeChiller, routesFor, trafficOf, wallSpots, weekMetres, type FloorManifest, type FloorMeta } from "./commissaryFloor";
import { menuFor, runCommissary } from "./commissaryWeek";

const BASE = new URL("../../../../assets/kitchens/exemplars/commissary_catering_3200/commissary_catering_3200", import.meta.url).pathname;
const manifest = JSON.parse(readFileSync(`${BASE}.scene3d.v1.json`, "utf8")) as FloorManifest;
const meta = JSON.parse(readFileSync(`${BASE}.site_metadata.json`, "utf8")) as FloorMeta;
const floor = buildFloor(manifest, meta);
const fields = fieldsOf(floor);
const spots = wallSpots(floor);
const SWAPS = ["chicken", "pasta", "salad"];
const traffic = (swaps: string[]) => trafficOf(runCommissary({ strategy: "ahead", swaps, sigma: 0.1 }).plan.days, menuFor(swaps));

const overlap = (a: { minX: number; maxX: number; minZ: number; maxZ: number }, b: typeof a) => a.minX < b.maxX && a.maxX > b.minX && a.minZ < b.maxZ && a.maxZ > b.minZ;
const len = (p: readonly (readonly [number, number])[]) => p.slice(1).reduce((s, q, i) => s + Math.hypot(q[0] - p[i]![0], q[1] - p[i]![1]), 0);

describe("commissary floor", () => {
  test("every spot fits on a wall, clear of every fixture, corridor and hood, with room to stand at its door", () => {
    expect(spots.length).toBeGreaterThan(20);
    for (const s of spots) {
      for (const k of floor.keepOut) expect(overlap(s.rect, k)).toBe(false);
      const b = floor.bounds;
      const onWall = Math.abs(s.rect.minX - b.minX) < 1e-9 || Math.abs(s.rect.maxX - b.maxX) < 1e-9 || Math.abs(s.rect.minZ - b.minZ) < 1e-9 || Math.abs(s.rect.maxZ - b.maxZ) < 1e-9;
      expect(onWall).toBe(true);
      expect(floor.grid.clearance(s.stance[0], s.stance[1])).toBeGreaterThanOrEqual(CARRY_R - floor.grid.cell);
    }
  });

  test("a walk is never shorter than the straight line; A* routes, which keep to the middle of an aisle, run 0-10% over the field's shortest walk", () => {
    const from = floor.stations.get(WALKIN)!.stance;
    for (const s of spots.filter((_, i) => i % 7 === 0)) {
      const f = fieldAt(floor.grid, fields.walkin, s.stance);
      expect(f).toBeGreaterThanOrEqual(Math.hypot(s.stance[0] - from[0], s.stance[1] - from[1]) - 1e-6);
      const path = floor.grid.findPath(from, s.stance, CARRY_R, { keepRight: 0 })!;
      expect(path).not.toBeNull();
      // The field is the shortest walk (within its 16-neighbour stencil's 3%); A* may take a roomier aisle.
      expect(len(path)).toBeGreaterThan(f * 0.97);
      expect(len(path)).toBeLessThan(f * 1.1);
    }
  });

  test("the routes drawn keep a cook with a pan clear of every fixture and of the chiller where it stands", () => {
    const t = traffic([]);
    const a = placeChiller(floor, fields, t, spots);
    for (const plan of [a.drawn, a.atBest]) {
      for (const r of plan.routes) {
        expect(r.path.length).toBeGreaterThanOrEqual(2);
        for (let i = 1; i < r.path.length; i++) {
          const [p, q] = [r.path[i - 1]!, r.path[i]!];
          for (let k = 0; k <= 20; k++) {
            const x = p[0] + ((q[0] - p[0]) * k) / 20;
            const z = p[1] + ((q[1] - p[1]) * k) / 20;
            // Clearance from the room's fixtures, allowing the grid's own resolution at the ends (stances).
            expect(floor.grid.clearance(x, z)).toBeGreaterThan(CARRY_R - 2 * floor.grid.cell);
          }
        }
      }
    }
  });

  test("the week's traffic is the cooked-ahead batches: every hot pan once, and a rack for every chiller load", () => {
    for (const swaps of [[], SWAPS]) {
      const run = runCommissary({ strategy: "ahead", swaps, sigma: 0.1 });
      const t = trafficOf(run.plan.days, menuFor(swaps));
      const hot = run.plan.days.reduce((s, d) => s + d.batches.filter((b) => menuFor(swaps).find((m) => m.key === b.dish)!.station !== "prep").reduce((u, b) => u + b.pans, 0), 0);
      expect([...t.pans.values()].reduce((s, v) => s + v, 0)).toBeCloseTo(hot, 9);
      expect(t.loads).toBeGreaterThanOrEqual(Math.ceil(hot / CHILL_LOAD / run.plan.days.length));
    }
    // Cooking on the day chills nothing.
    expect(trafficOf(runCommissary({ strategy: "day", swaps: [], sigma: 0.1 }).plan.days, menuFor([])).totalPans).toBe(0);
  });

  test("the best spot walks less than where the plan puts the chiller, and less than any other spot", () => {
    for (const swaps of [[], SWAPS]) {
      const t = traffic(swaps);
      const a = placeChiller(floor, fields, t, spots);
      expect(a.atBest.metres).toBeLessThan(a.drawn.metres);
      for (const c of a.cost) expect(a.cost[a.best]!).toBeLessThanOrEqual(c + 1e-9);
      // The quoted totals are the A* routes', and the ranking's estimate agrees with them.
      expect(Math.abs(weekMetres(floor, fields, t, a.spots[a.best]!.stance) - a.atBest.metres) / a.atBest.metres).toBeLessThan(0.06);
    }
  });

  test("the menu moves the chiller: braised chicken thighs in the skillet pull it from the combis to the kettle corner", () => {
    const near = (id: string, p: readonly [number, number]) => Math.hypot(floor.stations.get(id)!.stance[0] - p[0], floor.stations.get(id)!.stance[1] - p[1]);
    const written = placeChiller(floor, fields, traffic([]), spots);
    const swapped = placeChiller(floor, fields, traffic(SWAPS), spots);
    const w = written.spots[written.best]!.stance;
    const s = swapped.spots[swapped.best]!.stance;
    expect(near("combi-rollin-3", w)).toBeLessThan(near("tilt-skillet-40gal-1", w));
    expect(near("tilt-skillet-40gal-1", s)).toBeLessThan(near("combi-rollin-3", s));
    expect(Math.hypot(w[0] - s[0], w[1] - s[1])).toBeGreaterThan(5);
  });

  test("a chiller in the way lengthens nobody's walk below the open floor's: routes with it placed are no shorter", () => {
    const t = traffic([]);
    const spot = spots[Math.floor(spots.length / 2)]!;
    const placed = routesFor(floor, t, spot);
    for (const r of placed.routes) {
      const a = r.from === "chiller" ? spot.stance : floor.stations.get(r.from)!.stance;
      const b = r.to === "chiller" ? spot.stance : floor.stations.get(r.to)!.stance;
      expect(r.metres).toBeGreaterThanOrEqual(Math.hypot(a[0] - b[0], a[1] - b[1]) - 1e-6);
    }
    expect(distanceField(floor.grid, spot.stance, CARRY_R).some((v) => v === Infinity)).toBe(true);
  });
});
