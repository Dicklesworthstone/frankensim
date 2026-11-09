import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import initWasm, { PocketWatchSim } from "../src/wasm-pkg/fs_pocket_watch.js";

const __dirname = dirname(fileURLToPath(import.meta.url));
const wasmBytes = readFileSync(join(__dirname, "../src/wasm-pkg/fs_pocket_watch_bg.wasm"));

test("WASM PocketWatchSim — chronometer timekeeping, winding, and all 5 multiphysics perturbations", async () => {
  await initWasm({ module_or_path: wasmBytes });

  // 1. Chronometer timekeeping at 293.15 K (18,000 vph = 5 beats/sec)
  const sim = new PocketWatchSim();
  const snap0 = sim.snapshot();
  const tWatch0 = snap0[1];
  for (let i = 0; i < 100; i++) {
    sim.step(0.05); // 5.0 seconds total
  }
  const snap5 = sim.snapshot();
  const beats5 = snap5[14];
  const dtWatch = snap5[1] - tWatch0;
  assert.ok(beats5 >= 24 && beats5 <= 26, `expected ~25 beats in 5s, got ${beats5}`);
  assert.ok(Math.abs(dtWatch - 5.0) <= 0.25, `expected watch hands to advance ~5.0s, got ${dtWatch}s`);
  assert.equal(sim.status_label(), "RUNNING — CHRONOMETER CERTIFIED");

  // 2. Virtual Microphone 48 kHz synthesis
  sim.set_microphone(0, -11, 10, 1.5);
  const audio = sim.render_audio_block(1024);
  assert.equal(audio.length, 1024);
  let peak = 0;
  for (const s of audio) peak = Math.max(peak, Math.abs(s));
  assert.ok(peak > 1e-4, `expected non-zero modal audio peak, got ${peak}`);
  assert.ok(Number.isFinite(sim.mic_spl_db()), "expected finite mic_spl_db");
  assert.ok(sim.mic_loudness_sone() >= 0.0, "expected non-negative Zwicker loudness (sone)");
  assert.ok(sim.mic_loudness_phon() >= 0.0, "expected non-negative loudness level (phon)");
  assert.ok(sim.mic_sharpness_acum() >= 0.0, "expected non-negative DIN sharpness (acum)");
  assert.ok(sim.rate_sensitivity_s_day_per_k() < -1.0, "expected negative thermoelastic rate sensitivity d(rate)/dT");
  assert.ok(sim.rate_interval_lo_s_day() <= sim.rate_interval_hi_s_day(), "expected valid interval rate bounds");
  assert.ok(sim.balance_aero_drag_s_inv() > 0.0, "expected positive Stokes-Stuart oscillating-cylinder aero drag");
  assert.ok(sim.cosc_e_value() > 0.0, "expected positive COSC e-process wealth");
  assert.ok(sim.gdt_allocated_cost() > 0.0, "expected positive Lagrange GD&T allocated cost");
  const auditSummary = sim.multiphysics_crate_audit_summary();
  assert.equal(auditSummary.length, 9, "expected 9-element 51-crate audit summary");

  // 3. Blowtorch heating -> thermal expansion timekeeping drift -> 18k gold case melting -> liquid gold fouls gears
  sim.set_blowtorch(1.0, 0, -25, 0);
  for (let i = 0; i < 160; i++) {
    sim.step(0.05);
  }
  const snapMelt = sim.snapshot();
  assert.ok(snapMelt[26] > 0.25, `expected molten gold fouling factor > 0.25, got ${snapMelt[26]}`);
  const caseShell = sim.case_shell_snapshot();
  assert.equal(caseShell.length, 512 * 6, "expected 512 3D Lagrangian shell nodes * 6 floats");
  const spatialGrid = sim.spatial_grid_snapshot();
  assert.equal(spatialGrid.length, 2304 * 7, "expected 48x48 (2304) PDE cells * 7 floats");
  const coils = sim.mainspring_coils_mm();
  assert.equal(coils.length, 18, "expected 18 discretized elastica mainspring coil radii");

  // 4. Pour Liquid Nitrogen (LN2, 77.36 K) -> stops melting and freezes gold in deformed shape
  sim.set_blowtorch(0.0, 0, -25, 0);
  sim.set_liquid_nitrogen(1.0, 0, -18, 2);
  for (let i = 0; i < 110; i++) {
    sim.step(0.05);
  }
  assert.equal(sim.status_label(), "LOCKED — ENCASED IN SOLIDIFIED GOLD");

  // 5. Reset & test Watch Oil flooding + ultrasonic clean
  sim.reset();
  sim.set_oil_pour(1.0, 0, -11, 2);
  for (let i = 0; i < 20; i++) {
    sim.step(0.05);
  }
  const snapOil = sim.snapshot();
  assert.ok(snapOil[21] > 1.02, `expected capillary coil bridging to raise hairspring stiffness ratio, got ${snapOil[21]}`);
  sim.set_oil_pour(0.0, 0, -11, 2);
  sim.clean_excess_oil();

  // 6. External NdFeB N42 Super-Magnet pins pallet fork and leaves remanence until degaussed
  sim.set_magnet(true, 1.0, 10, -11, 2);
  sim.step(0.05);
  assert.equal(sim.status_label(), "STALLED — PALLET FORK MAGNETICALLY PINNED");
  sim.set_magnet(false, 1.0, 60, -11, 2);
  sim.step(0.05);
  assert.ok(sim.snapshot()[25] > 5.0, "expected residual ferromagnetic remanence");
  sim.trigger_degauss();
  for (let i = 0; i < 40; i++) {
    sim.step(0.05);
  }
  assert.equal(sim.snapshot()[25], 0.0, "expected AC degauss coil to clear remanence");
});
