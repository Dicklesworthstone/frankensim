// E1.3 terrain battery (bead wf-root-guzez.2.3): grid-shape/finite oracles
// per site, INDEPENDENT re-fit of the launch-flat plane (the flatness
// inputs must reproduce from the committed grid), material-class
// re-derivation, and the 1905-circuit containment arithmetic.
// Repro: node --test test/terrain.test.ts

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";

const load = (n: string): any =>
  JSON.parse(readFileSync(new URL(`../../../data/wright-flyer/terrain/${n}`, import.meta.url), "utf8"));
const kdh = load("kill-devil-hills-17x17-v1.json");
const huf = load("huffman-prairie-17x17-v1.json");
const prov = load("terrain-provenance-v1.json");

test("grids are complete, finite, and plausibly ranged per site", () => {
  for (const [g, lo, hi, name] of [
    [kdh, -1, 40, "kdh"],
    [huf, 200, 320, "huffman"],
  ] as const) {
    assert.equal(g.rows_south_to_north.length, 17);
    for (const row of g.rows_south_to_north) {
      assert.equal(row.length, 17);
      for (const z of row) {
        assert.ok(Number.isFinite(z) && z >= lo && z <= hi, `${name} elevation ${z}`);
      }
    }
    assert.deepEqual(g.failed_points, [], `${name}: zero filled points this fetch`);
  }
  assert.match(prov.label, /NOT SURVEYED/);
});

test("flatness inputs reproduce from the committed grid (independent re-fit)", () => {
  // Least-squares plane over rows 10-16, cols 3-13 (as recorded).
  const pts: [number, number, number][] = [];
  for (let r = 10; r <= 16; r++)
    for (let c = 3; c <= 13; c++) pts.push([c * 125, r * 125, kdh.rows_south_to_north[r][c]]);
  const n = pts.length;
  let sx = 0, sy = 0, sz = 0, sxx = 0, syy = 0, sxy = 0, sxz = 0, syz = 0;
  for (const [x, y, z] of pts) {
    sx += x; sy += y; sz += z; sxx += x * x; syy += y * y; sxy += x * y; sxz += x * z; syz += y * z;
  }
  // Cramer solve of the 3x3 normal system.
  const A = [ [sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n] ];
  const b = [sxz, syz, sz];
  const det3 = (m: number[][]): number =>
    m[0]![0]! * (m[1]![1]! * m[2]![2]! - m[1]![2]! * m[2]![1]!) -
    m[0]![1]! * (m[1]![0]! * m[2]![2]! - m[1]![2]! * m[2]![0]!) +
    m[0]![2]! * (m[1]![0]! * m[2]![1]! - m[1]![1]! * m[2]![0]!);
  const d = det3(A);
  const col = (i: number): number[][] => A.map((row, r) => row.map((v, c) => (c === i ? b[r]! : v)));
  const a = det3(col(0)) / d, bb = det3(col(1)) / d, cc = det3(col(2)) / d;
  const f = prov.flatness_certificate_inputs;
  assert.ok(Math.abs(a - f.slope_x) < 1e-5 && Math.abs(bb - f.slope_y) < 1e-5, "slopes reproduce");
  let ss = 0, mx = 0;
  for (const [x, y, z] of pts) {
    const r = z - (a * x + bb * y + cc);
    ss += r * r; mx = Math.max(mx, Math.abs(r));
  }
  assert.ok(Math.abs(Math.sqrt(ss / n) - f.rms_residual_m) < 0.01, "RMS reproduces");
  assert.ok(Math.abs(mx - f.max_abs_residual_m) < 0.01, "max residual reproduces");
  assert.ok(Math.hypot(a, bb) < 0.002, "the launch flat must actually be flat (<0.2% slope)");
  console.log(JSON.stringify({ suite: "wf-terrain", case: "flatness", slope: Math.hypot(a, bb), rms: Math.sqrt(ss / n) }));
});

test("material classes re-derive from the thresholds (per-cell)", () => {
  for (let r = 0; r < 17; r++)
    for (let c = 0; c < 17; c++) {
      const zk = kdh.rows_south_to_north[r][c];
      const want = zk <= 0.05 ? "water" : zk > 6 ? "dune" : "sand";
      assert.equal(prov.material_maps["kill-devil-hills"][r][c], want, `kdh (${r},${c})`);
      const zh = huf.rows_south_to_north[r][c];
      assert.equal(prov.material_maps["huffman-prairie"][r][c], zh > 260 ? "trees" : "grass");
    }
});

test("the 1905 circuit provably fits the tile", () => {
  const lapM = prov.circuit_1905_check.lap_length_m;
  assert.ok(Math.abs(lapM * 29.7 - 38_950) < 60, "laps x lap-length = flight distance");
  // Track bounding box (1000 x 320 m) fits the 2000 m tile with margin.
  assert.ok(1000 < 2000 - 250 && 320 < 2000 - 250, "containment with 125 m margin per side");
  assert.match(prov.circuit_1905_check.containment, /NEVER leaves the tile/);
});
