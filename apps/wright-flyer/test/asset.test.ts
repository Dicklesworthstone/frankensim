// E2.1 asset battery (bead wf-root-guzez.3.1): the in-repo model matches
// its provenance record byte-for-byte (size + sha256 recomputed), is a
// structurally valid draco GLB, and the provenance carries the license
// confirmations + the visual-only role boundary.
// Repro: node --test test/asset.test.ts

import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { test } from "node:test";

const load = (rel: string): Buffer =>
  readFileSync(new URL(`../../../data/wright-flyer/${rel}`, import.meta.url));

const provenance = JSON.parse(load("assets/flyer-model-provenance-v1.json").toString());
const glb = load("assets/flyer-smithsonian.glb");

test("model bytes match the provenance record exactly", () => {
  assert.equal(glb.length, provenance.asset.bytes, "byte length must match the record");
  const sha = createHash("sha256").update(glb).digest("hex");
  assert.equal(sha, provenance.asset.sha256, "sha256 must match the record");
  console.log(JSON.stringify({ suite: "wf-asset", case: "bytes", bytes: glb.length, sha256: sha }));
});

test("the file is a structurally valid glTF-2 binary with draco required", () => {
  assert.equal(glb.subarray(0, 4).toString("ascii"), "glTF", "GLB magic");
  assert.equal(glb.readUInt32LE(4), 2, "glTF version 2");
  assert.equal(glb.readUInt32LE(8), glb.length, "declared length equals file size");
  // First chunk is JSON; parse it and check the draco extension is declared.
  const jsonLen = glb.readUInt32LE(12);
  assert.equal(glb.readUInt32LE(16), 0x4e4f534a, "first chunk must be JSON");
  const gltf = JSON.parse(glb.subarray(20, 20 + jsonLen).toString("utf8"));
  assert.ok(
    gltf.extensionsRequired?.includes("KHR_draco_mesh_compression"),
    "draco must be declared required (the loader contract)",
  );
  assert.equal(gltf.meshes?.length, 1);
  console.log(
    JSON.stringify({ suite: "wf-asset", case: "structure", meshes: gltf.meshes.length }),
  );
});

test("STL companion matches its record and is a valid binary STL", () => {
  const stl = load(`assets/smithsonian-nasm-1903-flyer.cc0.stl`);
  assert.equal(stl.length, provenance.asset_stl.bytes, "STL byte length matches the record");
  // Binary STL: 80-byte header + u32 tri count + 50 bytes/tri.
  const tris = stl.readUInt32LE(80);
  assert.equal(84 + tris * 50, stl.length, "binary-STL structure must be self-consistent");
  assert.ok(tris > 100_000, `high-res source expected (${tris} tris)`);
  console.log(JSON.stringify({ suite: "wf-asset", case: "stl", tris }));
});

test("patent facsimile matches its provenance sha256", () => {
  const pdf = load(`patent/us-821393-wright-flyer.pdf`);
  const rec = JSON.parse(load(`patent/patent-us821393-v1.json`).toString());
  const sha = createHash("sha256").update(pdf).digest("hex");
  assert.equal(sha, rec.source_identity.facsimile_sha256, "patent facsimile immutable");
  const transcript = load(`patent/us-821393-transcript.md`).toString();
  assert.match(transcript, /twisted or warped in opposite directions/);
  assert.match(transcript, /rudder in conjunction with the movement/);
  assert.match(rec.claim_boundaries.forbidden[0], /QUANTITATIVE/);
  console.log(JSON.stringify({ suite: "wf-asset", case: "patent", sha256: sha.slice(0, 8) }));
});

test("provenance carries the license verdict, caveat, and role boundary", () => {
  assert.equal(provenance.license.verdict, "CC0 1.0 (Smithsonian Open Access)");
  assert.ok(provenance.license.confirmations.length >= 3, "three independent confirmations");
  assert.match(provenance.license.caveat, /Voyager/, "the boilerplate-header caveat is recorded");
  assert.match(provenance.role_boundary.visual_only, /NEVER as a source of section geometry/);
  assert.match(provenance.role_boundary.physics_binding, /FAVOR OF THE DOSSIERS/);
  assert.ok(provenance.alternatives_ranked.length >= 3, "fallbacks recorded per DONE-WHEN");
  assert.ok(provenance.no_claims.length > 40);
});
