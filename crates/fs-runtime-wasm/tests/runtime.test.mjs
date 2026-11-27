import test from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(here, "../../..");
const pkgDir = process.env.FS_RUNTIME_WASM_PKG
  ? path.resolve(process.env.FS_RUNTIME_WASM_PKG)
  : path.resolve(here, "../pkg");
const require = createRequire(import.meta.url);
const { WasmRuntime } = require(path.join(pkgDir, "fs_runtime_wasm.js"));

const f64Bits = (arr, idx) => new BigUint64Array(arr.buffer, arr.byteOffset, arr.length)[idx];
const f32Bits = (arr, idx) => new Uint32Array(arr.buffer, arr.byteOffset, arr.length)[idx];

test("forced_oscillator.fsim wasm matches native Rust bit-for-bit", () => {
  const fsim = fs.readFileSync(path.join(repoRoot, "examples/assemblies/forced_oscillator.fsim"));
  const steel = fs.readFileSync(
    path.join(repoRoot, "examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"),
  );
  const rt = WasmRuntime.load(fsim, steel, new Uint32Array([steel.length]), undefined);
  try {
    const inputs = JSON.parse(rt.inputs());
    assert.equal(inputs.length, 1);
    assert.equal(inputs[0].name, "load");
    assert.equal(inputs[0].unit, "N");
    assert.equal(inputs[0].default_value, 10);

    // Dynamic offset lookup from generated SnapshotDescriptor (zero hard-coded offsets).
    const desc = JSON.parse(rt.probes());
    assert.equal(desc.schema, "frankensim.runtime.snapshot_descriptor.v1");
    const pos = desc.probes.find((p) => p.name === "mass_pos");
    const vel = desc.probes.find((p) => p.name === "mass_vel");
    assert.ok(pos && vel && pos.offset !== vel.offset);

    for (let k = 0; k < 50; k++) {
      assert.equal(rt.step(1e-3), BigInt(k + 1));
    }
    const snap50 = new Float64Array(rt.snapshot_len());
    assert.equal(rt.read_probes(snap50), rt.snapshot_len());
    assert.equal(f64Bits(snap50, pos.offset), 0x3f78d8c356036c4cn);
    assert.equal(f64Bits(snap50, vel.offset), 0x3fcde0a297440181n);

    // Mid-run live load reversal + bit-for-bit parity + closed kernel ledger.
    rt.set_input_by_name("load", -5.0);
    for (let k = 0; k < 25; k++) {
      rt.step(1e-3);
    }
    const snap75 = new Float64Array(rt.snapshot_len());
    rt.read_probes(snap75);
    assert.equal(f64Bits(snap75, pos.offset), 0x3f85f471fa9d6701n);
    assert.equal(f64Bits(snap75, vel.offset), 0x3fc1a15bdf0f7a33n);

    // AudioWorklet render_f32 bit-for-bit parity.
    const block = new Float32Array(8);
    assert.equal(rt.render_f32(pos.index, 1e-3, 1, block), 8);
    assert.equal(f32Bits(block, 7), 0x3c3f9dbb);

    const ledger = JSON.parse(rt.ledger());
    assert.equal(ledger.schema, "frankensim.window_ledger.v1");
    assert.ok(Math.abs(ledger.totals.energy_residual_j) < 1e-12);

    // Refusal halts cleanly with exact silence and once-delivered refusal JSON.
    const bad = new Float32Array([1, 2, 3, 4]);
    assert.equal(rt.render_f32(pos.index, -1e-3, 1, bad), 0);
    assert.deepEqual([...bad], [0, 0, 0, 0]);
    assert.equal(rt.is_halted(), true);
    const ref1 = JSON.parse(rt.take_refusal());
    assert.equal(ref1.code, "invalid-dt");
    assert.equal(rt.take_refusal(), undefined);
  } finally {
    rt.free();
  }
});

test("bad pack framing refuses with canonical JSON", () => {
  const fsim = fs.readFileSync(path.join(repoRoot, "examples/assemblies/forced_oscillator.fsim"));
  const steel = fs.readFileSync(
    path.join(repoRoot, "examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"),
  );
  assert.throws(
    () => WasmRuntime.load(fsim, steel, new Uint32Array([steel.length - 1]), undefined),
    (err) => JSON.parse(String(err)).code === "invalid-pack-framing",
  );
  assert.throws(
    () => WasmRuntime.load(fsim, new Uint8Array(0), new Uint32Array(129), undefined),
    (err) => JSON.parse(String(err)).code === "invalid-pack-framing",
  );
});
