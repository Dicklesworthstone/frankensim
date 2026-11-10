// Behavioural tests of the wasm engine the page ships. Each scenario mirrors a guided experiment in
// src/main.ts and checks the outcome the physics produces (thresholds carry margin around the values
// the native calibration runs `presets_*` in engine-watch/src/lib.rs print).
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import initWasm, { PocketWatchSim } from "../src/wasm-pkg/fs_pocket_watch.js";

const __dirname = dirname(fileURLToPath(import.meta.url));
const wasmBytes = readFileSync(join(__dirname, "../src/wasm-pkg/fs_pocket_watch_bg.wasm"));
await initWasm({ module_or_path: wasmBytes });

const HEADER = 32;
const STRIDE = 12;
const node = (snap, i, field) => snap[HEADER + i * STRIDE + field];
const N_CASE_6H = 3;
const N_BALANCE = 10;
const T = 0;

const run = (sim, labS, dt = 0.1) => {
  for (let i = 0; i < Math.round(labS / dt); i++) sim.step(dt);
  return sim.snapshot();
};

const maxShellMelt = (sim) => {
  const shell = sim.case_shell_snapshot();
  let m = 0;
  for (let i = 5; i < shell.length; i += 6) m = Math.max(m, shell[i]);
  return m;
};

const withSim = (fn) => {
  const sim = new PocketWatchSim();
  try {
    return fn(sim);
  } finally {
    sim.free();
  }
};

test("a fresh watch keeps time at 18 000 vph and certifies", () =>
  withSim((sim) => {
    const t0 = sim.snapshot()[1];
    const s = run(sim, 5.0);
    assert.ok(s[14] >= 24 && s[14] <= 26, `expected ~25 beats in 5 s, got ${s[14]}`);
    assert.ok(Math.abs(s[1] - t0 - 5.0) <= 0.25, `hands advanced ${s[1] - t0} s in 5 s`);
    assert.ok(s[13] > 260 && s[13] < 320, `amplitude ${s[13]}°`);
    assert.equal(sim.status_label(), "RUNNING — CHRONOMETER CERTIFIED");
  }));

test("the buffers the page decodes have the documented layout", () =>
  withSim((sim) => {
    run(sim, 2.0);
    assert.equal(sim.snapshot().length, HEADER + 12 * STRIDE);
    assert.equal(sim.grid_n(), 64);
    assert.equal(sim.grid_half_span_mm(), 28.8);
    assert.equal(sim.spatial_grid_snapshot().length, 64 * 64 * 7);
    assert.equal(sim.case_shell_snapshot().length, 64 * 8 * 6);
    assert.equal(sim.mainspring_coils_mm().length, 18);
    const esc = sim.escapement_state();
    assert.equal(esc.length, 16);
    assert.ok(esc[13] > 30 && esc[13] < 60, `power reserve ${esc[13]} h`);
    const drained = sim.timegrapher_drain();
    assert.equal(drained.length % 3, 0);
    assert.ok(drained.length / 3 >= 8, `expected lock impacts from 2 s of running, got ${drained.length / 3}`);
    assert.equal(sim.timegrapher_drain().length, 0, "a drain empties the queue");
  }));

test("readouts follow the inputs at once, even while paused", () =>
  withSim((sim) => {
    sim.set_compensated_hairspring(false);
    const steel = sim.rate_sensitivity_s_day_per_k();
    assert.ok(steel < -8 && steel > -14, `carbon-steel hairspring d(rate)/dT = ${steel} s/d/K`);
    sim.set_compensated_hairspring(true);
    const comp = sim.rate_sensitivity_s_day_per_k();
    assert.ok(Math.abs(comp) < 0.5, `compensated hairspring d(rate)/dT = ${comp} s/d/K`);
    sim.set_ambient_pressure_atm(0.0);
    assert.equal(sim.balance_aero_drag_s_inv(), 0.0, "no aerodynamic drag in a vacuum");
    sim.set_ambient_pressure_atm(1.0);
    assert.ok(sim.balance_aero_drag_s_inv() > 0.0);
    assert.ok(sim.rate_interval_lo_s_day() <= sim.rate_interval_hi_s_day());
    assert.ok(sim.cosc_e_value() > 0.0);
  }));

test("same inputs replay to the same state digest", () => {
  const digest = (pour) =>
    withSim((sim) => {
      sim.set_blowtorch(0.5, 0, -22, 4);
      if (pour) sim.set_oil_pour(0.3, 0, -11, 2);
      run(sim, 3.0);
      return sim.state_digest_hex();
    });
  assert.equal(digest(false), digest(false));
  assert.notEqual(digest(false), digest(true));
});

test("the microphone hears the ticks without clipping, flame or not", () =>
  withSim((sim) => {
    sim.set_microphone(0, -11, 11, 1.4);
    const listen = (labS) => {
      let peak = 0;
      let rails = 0;
      for (let k = 0; k < Math.round((labS * 48000) / 1024); k++) {
        sim.step(1024 / 48000);
        for (const x of sim.render_audio_block(1024)) {
          peak = Math.max(peak, Math.abs(x));
          if (Math.abs(x) >= 0.975) rails++;
        }
      }
      return { peak, rails };
    };
    const quiet = listen(1.0);
    assert.ok(quiet.peak > 0.1 && quiet.peak < 0.9, `tick peak ${quiet.peak}`);
    assert.equal(quiet.rails, 0, "ticks never reach the rails at the default gain");
    assert.ok(Number.isFinite(sim.mic_spl_db()) && sim.mic_spl_db() > 30, `SPL ${sim.mic_spl_db()} dB`);
    assert.ok(sim.mic_loudness_sone() > 0 && sim.mic_loudness_phon() > 0 && sim.mic_sharpness_acum() > 0);
    sim.set_blowtorch(1.0, 0, -24, 2);
    const torch = listen(0.5);
    assert.ok(torch.peak < 0.95, `torch + ticks peak ${torch.peak}`);
  }));

test("a gentle flame warms the balance and the watch loses", () =>
  withSim((sim) => {
    sim.set_blowtorch(0.25, 0.0, -22.0, 6.0);
    const s = run(sim, 10.0);
    const tBal = node(s, N_BALANCE, T);
    assert.ok(tBal > 360 && tBal < 470, `balance at ${tBal} K`);
    assert.ok(s[16] < -50, `theoretical rate ${s[16]} s/d`);
    assert.ok(s[15] < -20, `measured rate ${s[15]} s/d`);
    assert.equal(sim.status_label(), "RUNNING — MULTIPHYSICS PERTURBED");
  }));

test("an oil flood stalls the balance; cleaning lets it restart", () =>
  withSim((sim) => {
    sim.set_oil_pour(0.85, 0, -11, 2);
    let s = run(sim, 8.0);
    assert.equal(sim.status_label(), "STALLED — VISCOUS OIL FLOOD OVERDAMPING");
    assert.ok(s[13] < 40, `flooded amplitude ${s[13]}°`);
    sim.set_oil_pour(0, 0, -11, 2);
    sim.clean_excess_oil();
    s = run(sim, 10.0);
    assert.match(sim.status_label(), /^RUNNING/);
    assert.ok(s[13] > 200, `amplitude after cleaning ${s[13]}°`);
  }));

test("a close magnet pins the lever and magnetises the spring until degaussed", () =>
  withSim((sim) => {
    sim.set_magnet(true, 1.0, 13, -11, 2.5);
    let s = run(sim, 3.0);
    assert.equal(sim.status_label(), "STALLED — PALLET FORK MAGNETICALLY PINNED");
    assert.ok(s[24] > 30, `field at the escapement ${s[24]} mT`);
    sim.set_magnet(false, 1.0, 60, -11, 2.5);
    s = run(sim, 1.0);
    assert.ok(s[24] < 1, `field after removal ${s[24]} mT`);
    assert.ok(s[25] > 50, `remanence after removal ${s[25]} mT`);
    sim.trigger_degauss();
    s = run(sim, 4.0);
    assert.ok(s[25] < 5, `remanence after degauss ${s[25]} mT`);
  }));

test("in a vacuum the amplitude climbs to the banking", () =>
  withSim((sim) => {
    sim.set_ambient_pressure_atm(0.0);
    const s = run(sim, 10.0);
    assert.ok(s[13] > 330, `vacuum amplitude ${s[13]}°`);
  }));

test("a 500 g shock under 4 g leaves the watch running", () =>
  withSim((sim) => {
    sim.set_gravity_vector(0, -4, 0);
    sim.trigger_mechanical_shock(500);
    const s = run(sim, 2.0);
    assert.match(sim.status_label(), /^RUNNING/);
    assert.ok(s[13] > 200, `amplitude ${s[13]}°`);
  }));

test("LN₂ pooled on a dial-up watch congeals the oil and seizes it", () =>
  withSim((sim) => {
    sim.set_gravity_vector(0, 0, -1);
    sim.set_liquid_nitrogen(1.0, 0, -11, 3);
    const s = run(sim, 55.0);
    assert.equal(sim.status_label(), "SEIZED — THERMAL / CRYOGENIC CLEARANCE LOSS");
    assert.ok(node(s, N_BALANCE, T) < 200, `balance at ${node(s, N_BALANCE, T)} K`);
    assert.ok(s[27] > 10, `LN₂ boiling power ${s[27]} W`);
  }));

test("a full flame melts the 6h band; LN₂ then quenches through its vapour film", () =>
  withSim((sim) => {
    sim.set_blowtorch(1.0, 0, -25, 1.5);
    let s = run(sim, 75.0);
    assert.ok(s[28] > 100, `torch power absorbed ${s[28]} W`);
    assert.ok(node(s, N_CASE_6H, T) > 1050, `6h case at ${node(s, N_CASE_6H, T)} K`);
    assert.ok(maxShellMelt(sim) > 0.5, `case-shell melt fraction ${maxShellMelt(sim)}`);
    assert.equal(sim.status_label(), "RUNNING — MULTIPHYSICS PERTURBED");
    sim.set_blowtorch(0, 0, -25, 1.5);
    sim.set_gravity_vector(0, 0, -1);
    sim.set_liquid_nitrogen(1.0, 0, -11, 3);
    s = run(sim, 20.0);
    assert.ok(s[27] > 100, `film-boiling power ${s[27]} W`);
    assert.ok(node(s, N_CASE_6H, T) < 800, `6h case after 20 s of LN₂ at ${node(s, N_CASE_6H, T)} K`);
    assert.equal(maxShellMelt(sim), 0, "the case has refrozen");
  }));
