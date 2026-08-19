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
  buildTerrainArrays,
  heightAt,
  materialClass,
  type TerrainGridJson,
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
  const { positions, colors, indices, launch } = buildTerrainArrays(grid, 64);
  assert.equal(positions.length, 65 * 65 * 3);
  assert.equal(colors.length, positions.length);
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
