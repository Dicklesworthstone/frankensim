import { describe, expect, test } from "bun:test";
import {
  createInitialCarrierState,
  createInitialOvenState,
  createInitialStoveState,
  createInitialWalkinState,
  evaluateCarrierSurrogate,
  stepCarrierKernel,
  stepOvenKernel,
  stepStoveKernel,
  stepWalkinKernel,
} from "./physicsKernels";
import * as refKernels from "../../../kitchen-lab/public/live/physics-kernels.mjs";

describe("conservation physics kernels parity & invariants", () => {
  test("stove kernel deglaze drops pan temp, releases vapor, and matches Rust-parity reference <1e-12", () => {
    const cfg = {
      source: "induction" as const,
      pan: "cast_iron" as const,
      fat: "canola" as const,
      powerFrac: 0.72,
      steakMassG: 340,
    };
    const s = createInitialStoveState(cfg);
    const sRef = refKernels.createInitialStoveState(cfg);
    stepStoveKernel(s, cfg, 45);
    refKernels.stepStoveKernel(sRef, cfg, 45);

    const tBefore = s.pan_surface_c;
    s.deglaze_liquid_g += 120.0;
    sRef.deglaze_liquid_g += 120.0;
    stepStoveKernel(s, cfg, 15.0);
    refKernels.stepStoveKernel(sRef, cfg, 15.0);

    expect(s.pan_surface_c).toBeLessThan(tBefore);
    expect(s.vapor_released_g).toBeGreaterThan(20);
    expect(s.fond_dissolved).toBeGreaterThan(0.05);
    expect(Math.abs(s.pan_surface_c - sRef.pan_surface_c)).toBeLessThan(1e-12);
    expect(Math.abs(s.vapor_released_g - sRef.vapor_released_g)).toBeLessThan(1e-12);
    expect(Math.abs(s.maillard_progress - sRef.maillard_progress)).toBeLessThan(1e-12);
  });

  test("oven kernel steam injection increases oven spring and matches reference <1e-12", () => {
    const dryCfg = { cavityTempC: 235, steamRh: 0.05, fanSpeed: 0.5, hydration: 0.75 };
    const steamCfg = { cavityTempC: 235, steamRh: 0.85, fanSpeed: 0.5, hydration: 0.75 };
    const sDry = createInitialOvenState(dryCfg);
    const sSteam = createInitialOvenState(steamCfg);
    const sSteamRef = refKernels.createInitialOvenState(steamCfg);

    stepOvenKernel(sDry, dryCfg, 180);
    stepOvenKernel(sSteam, steamCfg, 180);
    refKernels.stepOvenKernel(sSteamRef, steamCfg, 180);

    expect(sSteam.oven_spring_ratio).toBeGreaterThan(sDry.oven_spring_ratio);
    expect(Math.abs(sSteam.oven_spring_ratio - sSteamRef.oven_spring_ratio)).toBeLessThan(1e-12);
    expect(Math.abs(sSteam.crust_temp_c - sSteamRef.crust_temp_c)).toBeLessThan(1e-12);
  });

  test("walkin kernel 25 hot deep pans trigger HACCP 41F/5C excursion and matches reference <1e-12", () => {
    const cfg = { panCount: 25, panDepthIn: 6.0, doorOpenFrac: 0.25, kitchenAmbientC: 34 };
    const s = createInitialWalkinState();
    const sRef = refKernels.createInitialWalkinState();
    stepWalkinKernel(s, cfg, 1800);
    refKernels.stepWalkinKernel(sRef, cfg, 1800);

    expect(s.box_excursion_41f).toBe(true);
    expect(s.peak_box_air_c).toBeGreaterThan(5.0);
    expect(Math.abs(s.box_air_c - sRef.box_air_c)).toBeLessThan(1e-12);
    expect(Math.abs(s.frost_mass_kg - sRef.frost_mass_kg)).toBeLessThan(1e-12);
  });

  test("carrier surrogate detects out-of-domain inputs and matches reference <1e-12", () => {
    const inDom = evaluateCarrierSurrogate({ ventFraction: 0.1, foodMassKg: 4.5, ambientC: 22, initialFoodC: 84 });
    const ood = evaluateCarrierSurrogate({ ventFraction: 0.1, foodMassKg: 11.5, ambientC: 42, initialFoodC: 96 });
    const oodRef = refKernels.evaluateCarrierSurrogate({ ventFraction: 0.1, foodMassKg: 11.5, ambientC: 42, initialFoodC: 96 });
    expect(inDom.out_of_domain).toBe(false);
    expect(ood.out_of_domain).toBe(true);
    expect(ood.max_bound_min).toBeGreaterThan(inDom.max_bound_min);
    expect(Math.abs(ood.minutes_until_soggy - oodRef.minutes_until_soggy)).toBeLessThan(1e-12);

    const sSealed = createInitialCarrierState({ ventFraction: 0.0 });
    const sVented = createInitialCarrierState({ ventFraction: 0.9 });
    stepCarrierKernel(sSealed, { ventFraction: 0.0 }, 1200);
    stepCarrierKernel(sVented, { ventFraction: 0.9 }, 1200);
    expect(sVented.crispness_index).toBeGreaterThan(sSealed.crispness_index);
  });
});
