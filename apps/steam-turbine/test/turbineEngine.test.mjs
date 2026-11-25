// Node integration & physical emergence test suite against the REAL compiled
// release WebAssembly bundle (`src/wasm-pkg/fs_steam_turbine_bg.wasm`).
// Asserts non-input emergent physics, L2 MORPH geometry certification,
// fs-system::Ledger conservation laws, BLAKE3 determinism, and 48 kHz audio integrity.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";

import initWasm, { SteamTurbineSim } from "../src/wasm-pkg/fs_steam_turbine.js";

const here = dirname(fileURLToPath(import.meta.url));
const wasmBytes = readFileSync(resolve(here, "../src/wasm-pkg/fs_steam_turbine_bg.wasm"));
await initWasm({ module_or_path: wasmBytes });

test("L2 MORPH geometry catalog exports 24 valid FrankenSim meshes and watertight certified rotor", () => {
  const sim = new SteamTurbineSim();
  assert.equal(sim.mesh_count(), 24, "expected 24 exported FrankenSim meshes");

  for (let id = 0; id < 24; id++) {
    const pos = sim.mesh_positions(id);
    const norm = sim.mesh_normals(id);
    const idx = sim.mesh_indices(id);
    assert.ok(pos.length >= 9 && pos.length % 3 === 0, `mesh ${id} positions invalid`);
    assert.equal(norm.length, pos.length, `mesh ${id} normals length mismatch`);
    assert.ok(idx.length >= 3 && idx.length % 3 === 0, `mesh ${id} indices invalid`);
    for (let i = 0; i < Math.min(pos.length, 300); i++) {
      assert.ok(Number.isFinite(pos[i]), `mesh ${id} non-finite vertex`);
      assert.ok(Number.isFinite(norm[i]), `mesh ${id} non-finite normal`);
    }
  }

  const receipt = sim.geometry_receipt_flat();
  assert.equal(receipt.length, 14);
  const totalMassKg = receipt[2];
  const polarInertia = receipt[4];
  const hausdorffM = receipt[9];
  const eulerChar = receipt[10];
  const watertight = receipt[11];
  const sceneTris = receipt[13];

  assert.ok(totalMassKg > 2500 && totalMassKg < 9000, `unexpected rotor mass ${totalMassKg}`);
  assert.ok(polarInertia > 100 && polarInertia < 900, `unexpected J_p ${polarInertia}`);
  assert.ok(hausdorffM > 0 && hausdorffM < 0.005, `Hausdorff bound too coarse: ${hausdorffM}`);
  assert.equal(eulerChar, 2, "forged rotor must have Euler characteristic 2");
  assert.equal(watertight, 1, "forged rotor must be certified watertight");
  assert.ok(sceneTris > 100_000, `expected >100k instanced scene triangles, got ${sceneTris}`);

  const rows = sim.blade_rows_flat();
  assert.equal(rows.length, 34 * 9, `expected 34 blade/vane rows x 9 floats, got ${rows.length}`);
  sim.free();
});

test("Nominal 50 MW steady state closes fs-system::Ledger First Law, Mass, and Second Law", () => {
  const sim = new SteamTurbineSim();
  for (let i = 0; i < 40; i++) {
    assert.equal(sim.step(0.05), true, `step ${i} failed: ${sim.last_error()}`);
  }
  const s = sim.snapshot();
  assert.equal(s.length, 80);

  const rpm = s[2];
  const elecMw = s[5];
  const elecEff = s[7];
  const chpEff = s[8];
  const extBar = s[15];
  const condKpa = s[19];
  const lpQuality = s[21];
  const energyResJ = s[71];
  const massResKg = s[72];
  const entropyGen = s[73];
  const boundaryMj = s[74];

  assert.ok(Math.abs(rpm - 3000) < 35, `nominal rpm out of band: ${rpm}`);
  assert.ok(elecMw > 25 && elecMw < 80, `nominal electric MW out of band: ${elecMw}`);
  assert.ok(elecEff > 25 && elecEff < 48, `electric efficiency out of band: ${elecEff}`);
  assert.ok(chpEff > elecEff, `CHP eff (${chpEff}) must exceed electric eff (${elecEff})`);
  assert.ok(extBar > 10.5 && extBar < 22, `12-bar extraction header pressure out of band: ${extBar}`);
  assert.ok(condKpa > 4 && condKpa < 20, `condenser pressure out of band: ${condKpa}`);
  assert.ok(lpQuality > 0.80 && lpQuality < 0.99, `LP exhaust dryness x out of band: ${lpQuality}`);

  const relEnergy = Math.abs(energyResJ) / (Math.abs(boundaryMj) * 1e6 + 1);
  assert.ok(relEnergy < 1e-5, `First-Law relative residual too large: ${relEnergy}`);
  assert.ok(Math.abs(massResKg) < 1e-5, `Mass residual too large: ${massResKg}`);
  assert.ok(entropyGen > 0, `Second-Law entropy generation must be positive: ${entropyGen}`);

  // High-extraction CHP swing must raise combined utilization noticeably
  sim.set_load_setpoint(1.02);
  sim.set_extraction_demand(0.90);
  for (let i = 0; i < 35; i++) assert.equal(sim.step(0.05), true);
  const sChp = sim.snapshot();
  assert.ok(sChp[8] > chpEff + 5.0, `high CHP extraction must boost utilization: ${chpEff} -> ${sChp[8]}`);
  sim.free();
});

test("Emergent transients: Load Rejection, Condenser Air Ingress, Wet Steam Baumann, and Oil Whip", () => {
  // 1. Load Rejection (11 kV breaker trip -> speed surge + governor catch)
  const simLr = new SteamTurbineSim();
  for (let i = 0; i < 15; i++) assert.equal(simLr.step(0.05), true);
  const rpm0 = simLr.snapshot()[2];
  simLr.set_grid_breaker(false);
  let peakRpm = rpm0;
  for (let i = 0; i < 70; i++) {
    assert.equal(simLr.step(0.05), true);
    peakRpm = Math.max(peakRpm, simLr.snapshot()[2]);
  }
  const sLr = simLr.snapshot();
  assert.ok(peakRpm > rpm0 + 25, `expected speed surge after breaker open: ${rpm0} -> ${peakRpm}`);
  assert.ok(sLr[34] < 0.45, `expected governor to close CV after load rejection, got ${sLr[34]}`);
  assert.ok(sLr[5] < 1.0, `expected near-zero MW with breaker open, got ${sLr[5]}`);
  simLr.free();

  // 2. Condenser Air Ingress (Colburn-Hougen HTC collapse -> vacuum loss -> MW drop)
  const simAir = new SteamTurbineSim();
  for (let i = 0; i < 20; i++) assert.equal(simAir.step(0.05), true);
  const baseCond = simAir.snapshot()[19];
  const baseHtc = simAir.snapshot()[29];
  const baseMw = simAir.snapshot()[5];
  simAir.set_condenser_air_fraction(0.20);
  for (let i = 0; i < 100; i++) assert.equal(simAir.step(0.05), true);
  const afterAir = simAir.snapshot();
  assert.ok(afterAir[29] < baseHtc * 0.5, `HTC must drop with air ingress: ${baseHtc} -> ${afterAir[29]}`);
  assert.ok(afterAir[19] > baseCond + 4.0, `vacuum pressure must rise: ${baseCond} -> ${afterAir[19]}`);
  assert.ok(afterAir[5] < baseMw - 1.5, `output MW must drop under backpressure: ${baseMw} -> ${afterAir[5]}`);
  simAir.free();

  // 3. Wet Steam & Baumann Penalty (dropping live superheat cuts dryness x and IP/LP efficiency)
  const simWet = new SteamTurbineSim();
  for (let i = 0; i < 15; i++) assert.equal(simWet.step(0.05), true);
  const dryX = simWet.snapshot()[21];
  const dryEta = simWet.snapshot()[23];
  simWet.set_live_steam_temp_c(360);
  simWet.set_live_steam_bar(125);
  for (let i = 0; i < 40; i++) assert.equal(simWet.step(0.05), true);
  const wetSnap = simWet.snapshot();
  assert.ok(wetSnap[21] < dryX - 0.04, `LP dryness x must drop with cold steam: ${dryX} -> ${wetSnap[21]}`);
  assert.ok(wetSnap[23] < dryEta - 1.5, `Baumann moisture penalty must cut IP/LP eta: ${dryEta} -> ${wetSnap[23]}`);
  simWet.free();

  // 4. 4-DOF Rotordynamics Oil Whip & Unbalance (eigenvalues destabilize & orbit grows)
  const simRot = new SteamTurbineSim();
  for (let i = 0; i < 25; i++) assert.equal(simRot.step(0.05), true);
  const baseRe = simRot.snapshot()[47];
  const baseVib = simRot.snapshot()[43];
  simRot.set_oil_whip_destabilize(2.2);
  simRot.set_lube_oil_pressure_frac(0.40);
  simRot.set_rotor_unbalance_g_mm(120_000);
  for (let i = 0; i < 45; i++) assert.equal(simRot.step(0.05), true);
  const badRot = simRot.snapshot();
  assert.ok(badRot[47] > baseRe + 2.0, `4-DOF max Re(lambda) must destabilize: ${baseRe} -> ${badRot[47]}`);
  assert.ok(badRot[43] > baseVib + 25.0, `vibration amplitude must grow: ${baseVib} -> ${badRot[43]}`);
  simRot.free();
});

test("All 9 guided experiments, ESV manual trip coastdown, and extreme slider corners stay finite and close the Ledger", () => {
  const cases = [
    { name: "nominal", steps: 30, dt: 0.05, setup: (s) => s.reset(false) },
    {
      name: "load_rejection",
      steps: 45,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_grid_breaker(false);
      },
    },
    {
      name: "air_ingress",
      steps: 60,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_condenser_air_fraction(0.18);
      },
    },
    {
      name: "oil_whip",
      steps: 40,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_oil_whip_destabilize(1.95);
        s.set_lube_oil_pressure_frac(0.42);
        s.set_rotor_unbalance_g_mm(75_000);
      },
    },
    {
      name: "cold_start",
      steps: 60,
      dt: 0.05,
      setup: (s) => {
        s.reset(true);
        s.set_load_setpoint(0.28);
        s.set_rotor_unbalance_g_mm(18_000);
      },
    },
    {
      name: "wet_steam",
      steps: 40,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_live_steam_bar(125);
        s.set_live_steam_temp_c(365);
        s.set_extraction_demand(0.05);
      },
    },
    {
      name: "labyrinth_rub",
      steps: 50,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_rotor_unbalance_g_mm(280_000);
        s.set_oil_whip_destabilize(1.8);
        s.set_lube_oil_pressure_frac(0.50);
      },
    },
    {
      name: "chp_swing",
      steps: 40,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_load_setpoint(1.02);
        s.set_extraction_demand(0.90);
      },
    },
    {
      name: "lube_failure",
      steps: 60,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_lube_oil_pressure_frac(0.24);
        s.set_cw_inlet_temp_c(35.0);
      },
    },
    {
      name: "manual_esv_trip",
      steps: 60,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_manual_trip(true);
        s.set_grid_breaker(false);
      },
    },
    {
      name: "max_165bar_565C_full_load",
      steps: 30,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_live_steam_bar(165);
        s.set_live_steam_temp_c(565);
        s.set_load_setpoint(1.10);
      },
    },
    {
      name: "min_40bar_340C_zero_load",
      steps: 30,
      dt: 0.05,
      setup: (s) => {
        s.reset(false);
        s.set_live_steam_bar(40);
        s.set_live_steam_temp_c(340);
        s.set_load_setpoint(0.0);
        s.set_extraction_demand(0.0);
      },
    },
  ];

  const sim = new SteamTurbineSim();
  for (const c of cases) {
    c.setup(sim);
    let maxRubKn = 0;
    for (let i = 0; i < c.steps; i++) {
      assert.equal(sim.step(c.dt), true, `${c.name} step ${i} failed: ${sim.last_error()}`);
      maxRubKn = Math.max(maxRubKn, sim.snapshot()[46]);
    }
    const s = sim.snapshot();
    for (let k = 0; k < s.length; k++) {
      assert.ok(Number.isFinite(s[k]), `${c.name} snapshot[${k}] is non-finite: ${s[k]}`);
    }
    const relE = Math.abs(s[71]) / (Math.abs(s[74]) * 1e6 + 1);
    assert.ok(relE < 1e-4, `${c.name} First-Law relative residual too large: ${relE}`);
    assert.ok(Math.abs(s[72]) < 1e-4, `${c.name} mass residual too large: ${s[72]}`);
    assert.ok(s[73] >= 0, `${c.name} Second-Law entropy must be non-negative: ${s[73]}`);

    if (c.name === "manual_esv_trip") {
      assert.ok(s[2] > 2000 && s[2] < 3150, `ESV trip coastdown rpm unexpected: ${s[2]}`);
      assert.ok(s[35] < 0.05, `ESV must be closed after trip: ${s[35]}`);
    }
    if (c.name === "cold_start") {
      assert.notEqual(Math.round(s[9]), 5, `cold_start must cross 1st critical without tripping vibration, got vib=${s[43]}`);
    }
    if (c.name === "labyrinth_rub") {
      assert.ok(maxRubKn > 1.0, `labyrinth_rub must produce real Stribeck rub contact force > 1 kN, got ${maxRubKn} kN`);
    }
  }
  sim.free();
});

test("BLAKE3 bit-exact determinism and 48 kHz acoustic synthesis produce valid finite audio", () => {
  const a = new SteamTurbineSim();
  const b = new SteamTurbineSim();
  for (let i = 0; i < 20; i++) {
    assert.equal(a.step(0.04), true);
    assert.equal(b.step(0.04), true);
  }
  const digestA = a.state_digest_hex();
  const digestB = b.state_digest_hex();
  assert.equal(digestA.length, 64);
  assert.equal(digestA, digestB, "identical deterministic trajectories must share identical BLAKE3 digest");

  const pcm = a.render_audio_block(2048);
  assert.equal(pcm.length, 2048);
  let maxAbs = 0;
  for (let i = 0; i < pcm.length; i++) {
    assert.ok(Number.isFinite(pcm[i]), `non-finite audio sample at ${i}`);
    maxAbs = Math.max(maxAbs, Math.abs(pcm[i]));
  }
  assert.ok(maxAbs > 0.005 && maxAbs <= 1.0, `audio peak out of bounds: ${maxAbs}`);

  const spl = a.snapshot()[76];
  assert.ok(spl > 55 && spl < 120, `SPL dB out of industrial band: ${spl}`);
  a.free();
  b.free();
});

test("v2 engineering APIs: UVs, cutaway cap indices, 17-stage meanline, 2D PDE, fs-plate FEM & fs-dfo optimizer", () => {
  const sim = new SteamTurbineSim();

  // 1. Verify UV attributes and cutaway cap index splits on all 24 FrankenSim meshes
  for (let id = 0; id < 24; id++) {
    const pos = sim.mesh_positions(id);
    const uvs = sim.mesh_uvs(id);
    const idx = sim.mesh_indices(id);
    const capStart = sim.mesh_cap_start(id);
    assert.equal(uvs.length, (pos.length / 3) * 2, `mesh ${id} UV count mismatch`);
    if (id >= 2 && id <= 16) {
      assert.ok(capStart > 0 && capStart < idx.length, `cutaway mesh ${id} must split exterior vs cap indices`);
    }
  }

  for (let i = 0; i < 15; i++) assert.equal(sim.step(0.05), true);

  // 2. 17-stage row-by-row meanline telemetry
  const stages = sim.stage_profile_flat();
  assert.equal(stages.length, 17 * 16, `expected 17*16 stage floats, got ${stages.length}`);
  const p0 = stages[3];
  const p16 = stages[16 * 16 + 3];
  assert.ok(p0 > p16 * 100, `pressure must expand monotonically across 17 stages: ${p0} -> ${p16} bar`);

  // 3. 2D Axisymmetric thermal/stress PDE continuum field + surface probe
  const dims = sim.continuum_grid_dims();
  assert.equal(dims[0], 64);
  assert.equal(dims[1], 24);
  const cont = sim.continuum_field_flat();
  assert.equal(cont.length, 64 * 24 * 8);
  for (let i = 0; i < cont.length; i += 17) {
    assert.ok(Number.isFinite(cont[i]), `continuum[${i}] non-finite`);
  }
  const probeHp = sim.probe_continuum(1.45, 0.22);
  const probeLp = sim.probe_continuum(4.4, 0.90);
  assert.equal(probeHp.length, 8);
  // probe_continuum returns [temp_c, sigma_vm_mpa, p_bar, mach, q_mag_kw_m2, wetness_pct, creep_log, region_id]
  assert.ok(probeHp[0] > 300 && probeHp[0] < 580, `HP rotor temp_c out of band: ${probeHp[0]} °C`);
  assert.ok(probeLp[0] > 20 && probeLp[0] < 130, `LP hood temp_c out of band: ${probeLp[0]} °C`);
  assert.ok(probeHp[0] > probeLp[0] + 150, `HP rotor must be >150 °C hotter than LP hood: ${probeHp[0]} vs ${probeLp[0]}`);
  assert.ok(probeHp[4] >= 0 && probeLp[4] >= 0, `heat flux magnitude |q| must be non-negative`);
  assert.ok(probeHp[5] >= 0 && probeHp[5] <= 5, `superheated HP wetness_pct must be in [0, 5]%: ${probeHp[5]}`);
  assert.ok(probeLp[5] >= 0 && probeLp[5] <= 25, `LP wetness_pct must be in [0, 25]%: ${probeLp[5]}`);

  // 4. fs-plate::shell FEM Exhaust Hood eigenmode certificate + fs-duct Levine-Schwinger resistance
  const fem = sim.fem_shell_certificate_flat();
  assert.equal(fem.length, 6);
  assert.equal(fem[5], 1, "exhaust hood FEM certificate must be certified via Sylvester inertia");
  assert.ok(fem[0] > 20 && fem[1] > fem[0] && fem[2] > fem[1], `invalid FEM frequencies: ${fem.slice(0, 3)}`);
  assert.ok(fem[4] > 0.02 && fem[4] <= 2.0, `invalid Levine-Schwinger resistance: ${fem[4]}`);

  // 5. Parametric hardware geometry rebuild (live LP last-stage blade height change!)
  sim.apply_design_params(0.45, 0.46, 0.62, 185, 0.42, 85);
  const paramsAfter = sim.design_params_flat();
  assert.ok(Math.abs(paramsAfter[2] - 0.62) < 1e-6, "LP blade height must update live");
  const rows = sim.blade_rows_flat();
  const lp16TipR = rows[33 * 9 + 5];
  const lp16HubR = rows[33 * 9 + 4];
  assert.ok(Math.abs(lp16TipR - lp16HubR - 0.62) < 0.015, `3D LP Stage 16 geometry must reflect 0.62 m blade height`);

  // 6. L4 ASCENT fs-dfo CMA-ES + Nelder-Mead design optimization MUST lift LIVE fs-system plant output too!
  sim.apply_design_params(1.65, 0.20, 0.43, 280, 0.08, 130);
  const liveBefore = sim.snapshot();
  const liveBeforeMw = liveBefore[5];
  const liveBeforeEff = liveBefore[7];
  const opt = sim.run_design_optimization(0, 120);
  assert.ok(opt.length >= 28, "optimization report too short");
  const baseEff = opt[6];
  const optEff = opt[16];
  assert.ok(optEff > baseEff + 0.5, `fs-dfo must improve net electric efficiency: ${baseEff.toFixed(2)}% -> ${optEff.toFixed(2)}%`);
  assert.ok(opt[26] >= 60, `expected >= 60 fs-dfo evaluations, got ${opt[26]}`);
  sim.apply_design_params(opt[0], opt[1], opt[2], opt[3], opt[4], opt[5]);
  for (let i = 0; i < 20; i++) assert.equal(sim.step(0.05), true);
  const liveAfter = sim.snapshot();
  assert.ok(
    liveAfter[5] > liveBeforeMw + 0.5,
    `applying fs-dfo optimum must raise LIVE fs-system MW: ${liveBeforeMw.toFixed(2)} -> ${liveAfter[5].toFixed(2)} MW`,
  );
  assert.ok(
    liveAfter[7] > liveBeforeEff + 0.3,
    `applying fs-dfo optimum must raise LIVE fs-system eff: ${liveBeforeEff.toFixed(2)}% -> ${liveAfter[7].toFixed(2)}%`,
  );
  sim.free();
});
