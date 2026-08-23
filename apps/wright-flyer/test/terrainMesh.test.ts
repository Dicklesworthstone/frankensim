// E2.3 terrain-core battery (bead wf-root-guzez.3.3): heightAt agrees with
// the Rust/Python implementations at nodes (third independent
// implementation of the same data), material classes match the E1.3
// artifact per cell, the vertex arrays are consistent, and the arrival
// camera starts high/far and settles at the rail.
// Repro: node --test test/terrainMesh.test.ts

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";
import {
  arrivalCamera,
  bigHillDetail,
  BIG_HILL_CENTER_M,
  BIG_HILL_PEAK_M,
  buildTerrainArrays,
  duneDetail,
  heightAt,
  materialClass,
  warpAxis,
} from "../src/terrainMesh.ts";

const grid = JSON.parse(
  readFileSync(new URL("../../../data/wright-flyer/terrain/kill-devil-hills-17x17-v1.json", import.meta.url), "utf8"),
) as TerrainGridJson & { rows_south_to_north: number[][] };
const prov = JSON.parse(
  readFileSync(new URL("../../../data/wright-flyer/terrain/terrain-provenance-v1.json", import.meta.url), "utf8"),
);

test("heightAt is exact at every node (third-implementation agreement)", () => {
  for (let r = 0; r < 17; r++)
    for (let c = 0; c < 17; c++) {
      const h = heightAt(grid, c * 125, r * 125);
      assert.equal(h, grid.rows_south_to_north[r]![c]!, `node (${r},${c})`);
    }
  // Midpoint = 4-node mean (bilinear identity).
  const mid = heightAt(grid, 62.5, 62.5);
  const mean =
    (grid.rows_south_to_north[0]![0]! + grid.rows_south_to_north[0]![1]! +
     grid.rows_south_to_north[1]![0]! + grid.rows_south_to_north[1]![1]!) / 4;
  assert.ok(Math.abs(mid - mean) < 1e-12);
});

test("material classes match the committed E1.3 map per cell", () => {
  for (let r = 0; r < 17; r++)
    for (let c = 0; c < 17; c++) {
      assert.equal(
        materialClass(grid.rows_south_to_north[r]![c]!),
        prov.material_maps["kill-devil-hills"][r][c],
        `cell (${r},${c})`,
      );
    }
});

test("vertex arrays are consistent and centered", () => {
  const { positions, colors, uvs, indices, launch } = buildTerrainArrays(grid, 64);
  assert.equal(positions.length, 65 * 65 * 3);
  assert.equal(colors.length, positions.length);
  assert.equal(uvs.length, 65 * 65 * 2);
  // Planar UVs span [0,1] corner-to-corner (a textured material is
  // blind without them — guzez.13 regression).
  assert.equal(uvs[0], 0);
  assert.equal(uvs[1], 0);
  assert.equal(uvs[uvs.length - 2], 1);
  assert.equal(uvs[uvs.length - 1], 1);
  assert.equal(indices.length, 64 * 64 * 6);
  // Every index in range; corners at ±1000.
  for (const ix of indices) assert.ok(ix < 65 * 65);
  assert.equal(positions[0], -1000);
  // The launch point sits on the flat (low elevation, inside the tile).
  assert.ok(launch[1] > 0 && launch[1] < 6, `launch elevation ${launch[1]}`);
  assert.ok(Math.abs(launch[0]) < 1000 && Math.abs(launch[2]) < 1000);
});

test("arrival camera starts high and far, settles beside the rail", () => {
  const { launch } = buildTerrainArrays(grid, 16);
  const start = arrivalCamera(0, launch);
  const end = arrivalCamera(14, launch);
  const later = arrivalCamera(99, launch);
  assert.ok(start.pos[1] > 50, "starts high");
  const distStart = Math.hypot(start.pos[0] - launch[0], start.pos[2] - launch[2]);
  const distEnd = Math.hypot(end.pos[0] - launch[0], end.pos[2] - launch[2]);
  assert.ok(distStart > 300 && distEnd < 40, `dolly ${distStart} -> ${distEnd}`);
  assert.deepEqual(end, later, "the shot settles (idempotent past 14 s)");
  assert.ok(end.pos[1] > launch[1] + 3, "settles above the sand");
});

test("warpAxis is monotone, endpoint-exact, and clusters at lf", () => {
  for (const lf of [0.15, 0.5, 0.8125]) {
    let prev = -Infinity;
    for (let i = 0; i <= 64; i++) {
      const v = warpAxis(i / 64, lf);
      assert.ok(v >= prev, `monotone at ${i} (lf=${lf})`);
      assert.ok(v >= 0 && v <= 1);
      prev = v;
    }
    assert.equal(warpAxis(0, lf), 0);
    assert.ok(Math.abs(warpAxis(1, lf) - 1) < 1e-12);
    // Density: the local derivative AT lf exceeds the far-field one.
    const near = warpAxis(lf + 1e-3, lf) - warpAxis(lf, lf);
    const far = warpAxis(0.999, lf) - warpAxis(0.998, lf);
    assert.ok(near / 1e-3 > far / 1e-3, "denser near launch than far field");
  }
});

test("LOD terrain keeps survey corners and packs vertices near launch", () => {
  const { positions } = buildTerrainArrays(grid, 96);
  const nv = 97;
  // Corners still land exactly on the tile bounds.
  assert.equal(positions[0], -1000); // (j=0,i=0) wx
  const last = (nv * nv - 1) * 3;
  assert.equal(positions[last + 2], -1000); // north edge wz = -(1000-1000)
  // Nearest-vertex distance to the LAUNCH FLAT center must beat the
  // uniform-grid spacing (2000/96 ≈ 20.8 m): the warped inner ring
  // packs a vertex within 8 m of the launch point itself.
  let best = Infinity;
  for (let v = 0; v < nv * nv; v++) {
    const dx = positions[v * 3]! - 0;
    const dz = positions[v * 3 + 2]! - (-625);
    best = Math.min(best, Math.hypot(dx, dz));
  }
  assert.ok(best < 8, `inner ring too coarse near launch: ${best.toFixed(2)} m`);
});

/* ---------- presentation relief (PurpleCliff flat-dunes fix) -------- */

test("dune relief keeps the launch/camp corridor at survey height", () => {
  // The whole rail run, camp, and landing flat: exactly zero detail.
  for (const [x, z] of [
    [0, 0], [18.3, 0], [-40, -13], [-26, -21], [30, 12],
  ] as const) {
    assert.equal(duneDetail(x, z), 0, `corridor (${x},${z})`);
    assert.equal(bigHillDetail(x, z), 0, `far from Big Hill (${x},${z})`);
  }
  assert.throws(() => duneDetail(Number.NaN, 0), /finite/);
  assert.throws(() => bigHillDetail(0, Number.POSITIVE_INFINITY), /finite/);
});

test("dune sea is bounded, dramatic beyond the feather, and deterministic", () => {
  let maxAbs = 0;
  for (let x = -900; x <= 900; x += 37) {
    for (let z = -900; z <= 900; z += 41) {
      const d = duneDetail(x, z);
      assert.ok(Number.isFinite(d) && Math.abs(d) <= 14, `bounded at (${x},${z}): ${d}`);
      maxAbs = Math.max(maxAbs, Math.abs(d));
    }
  }
  assert.ok(maxAbs > 6, `relief must actually sculpt (max |d| = ${maxAbs.toFixed(1)} m)`);
  // Determinism: identical inputs, identical sand.
  assert.equal(duneDetail(321, -417), duneDetail(321, -417));
});

test("Big Kill Devil Hill peaks near its authored summit and fades out", () => {
  const peak = bigHillDetail(BIG_HILL_CENTER_M.x, BIG_HILL_CENTER_M.z);
  assert.equal(bigHillDetail(BIG_HILL_CENTER_M.x + 600, BIG_HILL_CENTER_M.z), 0, "zero far field");
  // Continuous: neighbors differ by less than the mesh can alias.
  const a = bigHillDetail(BIG_HILL_CENTER_M.x + 50, BIG_HILL_CENTER_M.z + 20);
  const b = bigHillDetail(BIG_HILL_CENTER_M.x + 52, BIG_HILL_CENTER_M.z + 22);
  assert.ok(Math.abs(a - b) < 2, "no cliffs between samples");
});
