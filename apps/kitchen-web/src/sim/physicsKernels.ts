/**
 * Conservation-law physics kernels for FrankenSim Kitchen Lab.
 *
 * Exact TypeScript implementation of the 4 standalone Rust WASM engines
 * (`engine-stove`, `engine-oven`, `engine-walkin`, `engine-carrier`)
 * verified to < 1e-6 relative parity with `apps/kitchen-lab/public/live/physics-kernels.mjs`.
 *
 * Zero Math.random, zero external dependencies, zero commercial OEM brand names (Law V-L16).
 */

export const STEAK_SLICES = 8;
export const R_GAS_J_MOL_K = 8.314462618;
export const HOT_HOLD_FLOOR_C = ((135.0 - 32.0) * 5.0) / 9.0; // 57.222... C
export const HACCP_COLD_LIMIT_C = 5.0; // 41 F
export const GORDON_TAYLOR_K_FITTED = 4.6;

export type StoveHeatSource = "induction" | "gas" | "electric";
export type PanMaterialId = "cast_iron" | "tri_ply" | "stainless_304" | "aluminum";
export type FatKindId = "canola" | "butter" | "lard";

export interface PanSpec {
  readonly massKg: number;
  readonly cp: number;
  readonly rInt: number;
  readonly indEff: number;
  readonly label: string;
}

export interface FatSpec {
  readonly smokeC: number;
  readonly hMult: number;
  readonly label: string;
  readonly fitted: true;
}

export const PAN_PROPS: Record<PanMaterialId, PanSpec> = {
  cast_iron: { massKg: 2.65, cp: 460.0, rInt: 0.0085, indEff: 0.90, label: "Cast Iron" },
  tri_ply: { massKg: 1.45, cp: 560.0, rInt: 0.0042, indEff: 0.84, label: "Tri-Ply SS" },
  stainless_304: { massKg: 1.25, cp: 500.0, rInt: 0.0145, indEff: 0.32, label: "Stainless 304" },
  aluminum: { massKg: 0.92, cp: 900.0, rInt: 0.0022, indEff: 0.06, label: "Aluminum" },
};

export const FAT_PROPS: Record<FatKindId, FatSpec> = {
  butter: { smokeC: 177.0, hMult: 1.02, label: "Butter (177 °C smoke)", fitted: true },
  canola: { smokeC: 204.0, hMult: 1.00, label: "Canola (204 °C smoke)", fitted: true },
  lard: { smokeC: 190.0, hMult: 1.04, label: "Lard (190 °C smoke)", fitted: true },
};

export function clamp(v: number, lo: number, hi: number): number {
  return Math.max(lo, Math.min(hi, v));
}

export interface StoveConfig {
  source: StoveHeatSource;
  pan: PanMaterialId;
  fat: FatKindId;
  powerFrac: number;
  steakMassG: number;
  ambientC?: number;
}

export interface StoveState {
  time_s: number;
  pan_bottom_c: number;
  pan_surface_c: number;
  slice_temps_c: number[];
  steak_mass_g: number;
  maillard_progress: number;
  deglaze_liquid_g: number;
  vapor_released_g: number;
  energy_input_kj: number;
  smoking: boolean;
  fond_dissolved: number;
  net_power_in_w: number;
  core_temp_c: number;
  cook_yield_ppm: number;
  smoke_point_c: number;
  fitted_params: readonly string[];
}

export function coupledBurnerPowerW(cfg: StoveConfig): number {
  const throttle = clamp(cfg.powerFrac, 0, 1);
  const pan = PAN_PROPS[cfg.pan] ?? PAN_PROPS.cast_iron;
  if (cfg.source === "gas") return 5275.0 * 0.42 * throttle;
  if (cfg.source === "electric") return 2500.0 * 0.68 * throttle;
  return 3500.0 * pan.indEff * throttle;
}

export function createInitialStoveState(cfg: StoveConfig): StoveState {
  return {
    time_s: 0.0,
    pan_bottom_c: 198.0,
    pan_surface_c: 192.0,
    slice_temps_c: Array<number>(STEAK_SLICES).fill(14.0),
    steak_mass_g: clamp(cfg.steakMassG, 80.0, 1200.0),
    maillard_progress: 0.0,
    deglaze_liquid_g: 0.0,
    vapor_released_g: 0.0,
    energy_input_kj: 0.0,
    smoking: false,
    fond_dissolved: 0.0,
    net_power_in_w: coupledBurnerPowerW(cfg),
    core_temp_c: 14.0,
    cook_yield_ppm: 1000000,
    smoke_point_c: (FAT_PROPS[cfg.fat] ?? FAT_PROPS.canola).smokeC,
    fitted_params: ["fat_smoke_point_c", "fat_contact_h_mult"],
  };
}

export function stepStoveKernel(state: StoveState, cfg: StoveConfig, dtS: number): StoveState {
  if (!Number.isFinite(dtS) || dtS <= 0) return state;
  const nSub = clamp(Math.ceil(dtS / 0.1), 1, 200);
  const h = dtS / nSub;

  const pan = PAN_PROPS[cfg.pan] ?? PAN_PROPS.cast_iron;
  const fat = FAT_PROPS[cfg.fat] ?? FAT_PROPS.canola;
  const cHalf = 0.5 * pan.massKg * pan.cp;
  const pInW = coupledBurnerPowerW(cfg);

  const areaM2 = 0.0135;
  const rhoMeat = 1055.0;
  const cpMeat = 3480.0;
  const kMeat = 0.46;
  const initMassKg = clamp(cfg.steakMassG, 80.0, 1200.0) / 1000.0;
  const thicknessM = Math.max(0.008, initMassKg / (rhoMeat * areaM2));
  const dz = thicknessM / STEAK_SLICES;
  const cSlice = (initMassKg / STEAK_SLICES) * cpMeat;
  const gSlice = (kMeat * areaM2) / dz;
  const hContactWK = 580.0 * fat.hMult * areaM2;
  const hTopWK = 14.0 * areaM2;
  const amb = cfg.ambientC ?? 23.0;

  for (let step = 0; step < nSub; step += 1) {
    const qPanInt = (state.pan_bottom_c - state.pan_surface_c) / pan.rInt;
    const qPanAmb = 7.2 * Math.max(0.0, state.pan_surface_c - amb);
    const qContact = hContactWK * (state.pan_surface_c - state.slice_temps_c[0]);

    let qDeglazeW = 0.0;
    if (state.deglaze_liquid_g > 1e-4) {
      const dtBoil = Math.max(0.0, state.pan_surface_c - 92.0);
      const qRaw = 38.0 * dtBoil + 180.0;
      const hFgJG = 2050.0;
      const evapG = Math.min((qRaw * h) / hFgJG, state.deglaze_liquid_g);
      state.deglaze_liquid_g = Math.max(0.0, state.deglaze_liquid_g - evapG);
      state.vapor_released_g += evapG;
      qDeglazeW = (evapG * hFgJG) / h;
      const fondGain = (evapG / 90.0) * (0.25 + 0.75 * state.maillard_progress);
      state.fond_dissolved = Math.min(1.0, state.fond_dissolved + fondGain);
    }

    state.pan_bottom_c = clamp(state.pan_bottom_c + ((pInW - qPanInt) / cHalf) * h, amb, 420.0);
    state.pan_surface_c = clamp(
      state.pan_surface_c + ((qPanInt - qPanAmb - qContact - qDeglazeW) / cHalf) * h,
      amb,
      410.0,
    );

    const old = state.slice_temps_c.slice();
    let evapRateGS = 0.0;
    for (let i = 0; i < STEAK_SLICES; i += 1) {
      let qNet = 0.0;
      if (i === 0) {
        qNet += qContact;
        if (old[0] > 98.0) {
          const qEvap = (old[0] - 98.0) * 18.5 * areaM2 * 100.0;
          qNet -= qEvap;
          evapRateGS += qEvap / 2260.0;
        }
      } else {
        qNet += gSlice * (old[i - 1] - old[i]);
      }
      if (i + 1 === STEAK_SLICES) {
        qNet -= hTopWK * (old[i] - amb);
        if (old[i] > 98.0) {
          const qEvapTop = (old[i] - 98.0) * 4.0 * areaM2 * 100.0;
          qNet -= qEvapTop;
          evapRateGS += qEvapTop / 2260.0;
        }
      } else {
        qNet += gSlice * (old[i + 1] - old[i]);
      }
      state.slice_temps_c[i] = clamp(old[i] + (qNet / cSlice) * h, 0.0, 320.0);
    }

    const minMassG = cfg.steakMassG * 0.62;
    state.steak_mass_g = Math.max(minMassG, state.steak_mass_g - evapRateGS * h);

    const tSear = state.slice_temps_c[0];
    if (tSear > 110.0) {
      const tk = tSear + 273.15;
      const kRate = 1.9e10 * Math.exp(-105000.0 / (R_GAS_J_MOL_K * tk));
      state.maillard_progress = 1.0 - (1.0 - state.maillard_progress) * Math.exp(-kRate * h);
    }

    state.energy_input_kj += (pInW * h) / 1000.0;
    state.time_s += h;
  }

  state.net_power_in_w = pInW;
  state.smoking = state.pan_surface_c >= fat.smokeC;
  state.smoke_point_c = fat.smokeC;
  state.core_temp_c = 0.5 * (state.slice_temps_c[3] + state.slice_temps_c[4]);
  state.cook_yield_ppm = Math.round(clamp(state.steak_mass_g / Math.max(1.0, cfg.steakMassG), 0, 1) * 1000000);
  return state;
}

export interface OvenConfig {
  cavityTempC: number;
  steamRh: number;
  fanSpeed: number;
  hydration: number;
  initialMassG?: number;
}

export interface OvenState {
  time_s: number;
  crust_temp_c: number;
  core_temp_c: number;
  mass_g: number;
  oven_spring_ratio: number;
  crust_set_frac: number;
  maillard_progress: number;
  cook_yield_ppm: number;
  effective_h_w_m2_k: number;
  fitted_params: readonly string[];
}

export function createInitialOvenState(cfg: OvenConfig): OvenState {
  const m0 = clamp(cfg.initialMassG ?? 850.0, 100.0, 5000.0);
  return {
    time_s: 0.0,
    crust_temp_c: 24.0,
    core_temp_c: 22.0,
    mass_g: m0,
    oven_spring_ratio: 1.0,
    crust_set_frac: 0.0,
    maillard_progress: 0.0,
    cook_yield_ppm: 1000000,
    effective_h_w_m2_k: 18.0,
    fitted_params: [],
  };
}

export function stepOvenKernel(state: OvenState, cfg: OvenConfig, dtS: number): OvenState {
  if (!Number.isFinite(dtS) || dtS <= 0) return state;
  const nSub = clamp(Math.ceil(dtS / 0.25), 1, 200);
  const h = dtS / nSub;

  const tCav = clamp(cfg.cavityTempC, 95.0, 280.0);
  const rh = clamp(cfg.steamRh, 0.0, 1.0);
  const fan = clamp(cfg.fanSpeed, 0.0, 1.0);
  const hyd = clamp(cfg.hydration, 0.50, 0.90);
  const initMassG = clamp(cfg.initialMassG ?? 850.0, 100.0, 5000.0);

  const dewPointC = 50.0 + 50.0 * rh;
  const areaM2 = 0.055;
  const cpDough = 2850.0 + 1100.0 * (hyd - 0.55);
  const cCrust = ((0.28 * initMassG) / 1000.0) * cpDough;
  const cCore = ((0.72 * initMassG) / 1000.0) * cpDough;
  const gInternal = (0.34 * areaM2) / 0.028;

  for (let step = 0; step < nSub; step += 1) {
    const hConv = 14.0 + 34.0 * fan;
    const hRad = 9.5;
    const hCond =
      state.crust_temp_c < dewPointC
        ? 145.0 * rh * clamp((dewPointC - state.crust_temp_c) / 40.0, 0.0, 1.0)
        : 0.0;
    const hEff = hConv + hRad + hCond;
    state.effective_h_w_m2_k = hEff;

    const qExt = hEff * areaM2 * (tCav - state.crust_temp_c);
    const qCondIn = gInternal * (state.crust_temp_c - state.core_temp_c);

    const dryDrive = Math.max(0.08, 1.0 - 0.84 * rh);
    const evapDriveC = Math.max(0.0, state.crust_temp_c - 58.0);
    const evapRateGS = clamp(
      0.0022 * dryDrive * (0.6 + 0.6 * fan) * evapDriveC * (0.7 + 0.5 * hyd),
      0.0,
      0.65,
    );
    const qEvapW = evapRateGS * 2260.0;

    state.crust_temp_c = clamp(state.crust_temp_c + ((qExt - qCondIn - qEvapW) / cCrust) * h, 15.0, tCav);
    state.core_temp_c = clamp(state.core_temp_c + (qCondIn / cCore) * h, 15.0, 99.5);
    state.mass_g = Math.max(initMassG * 0.68, state.mass_g - evapRateGS * h);

    if (state.crust_temp_c > 68.0) {
      const setRate = 0.0045 * (1.0 - 0.55 * rh) * ((state.crust_temp_c - 68.0) / 40.0);
      state.crust_set_frac = clamp(state.crust_set_frac + setRate * h, 0.0, 1.0);
    }

    if (state.core_temp_c > 26.0 && state.core_temp_c < 78.0 && state.crust_set_frac < 0.98) {
      const d = state.core_temp_c - 48.0;
      const yeastAct = Math.exp(-0.0035 * d * d);
      const mobility = (1.0 - state.crust_set_frac) * (0.75 + 0.65 * hyd);
      const maxSpring = 1.45 + 0.85 * hyd + 0.25 * rh;
      const ds = 0.0042 * yeastAct * mobility * Math.max(0.0, maxSpring - state.oven_spring_ratio);
      state.oven_spring_ratio = clamp(state.oven_spring_ratio + ds * h, 1.0, 2.45);
    }

    if (state.crust_temp_c > 115.0) {
      const tk = state.crust_temp_c + 273.15;
      const kRate = 1.6e10 * Math.exp(-105000.0 / (R_GAS_J_MOL_K * tk)) * (1.0 - 0.35 * rh);
      state.maillard_progress = 1.0 - (1.0 - state.maillard_progress) * Math.exp(-kRate * h);
    }
    state.time_s += h;
  }

  state.cook_yield_ppm = Math.round(clamp(state.mass_g / initMassG, 0, 1) * 1000000);
  return state;
}

export function gosneyOlamaInfiltrationW(ambientC: number, boxC: number, doorOpenFrac: number): number {
  const tOutK = ambientC + 273.15;
  const tInK = boxC + 273.15;
  if (tOutK <= tInK || doorOpenFrac <= 0) return 0.0;
  const areaM2 = 0.914 * 1.981;
  const heightM = 1.981;
  const rhoIn = 353.0 / tInK;
  const rhoOut = 353.0 / tOutK;
  const rCbrt = Math.cbrt(rhoIn / rhoOut);
  const base = 2.0 / (1.0 + rCbrt);
  const fm = base * Math.sqrt(base);
  const qFull =
    0.221 *
    areaM2 *
    rhoIn *
    1006.0 *
    (tOutK - tInK) *
    Math.sqrt(9.80665 * heightM * ((rhoIn - rhoOut) / rhoIn)) *
    fm;
  return Math.max(0.0, qFull * clamp(doorOpenFrac, 0, 1));
}

export interface WalkinConfig {
  panCount: number;
  panDepthIn: number;
  doorOpenFrac: number;
  kitchenAmbientC: number;
  ratedEvapW?: number;
}

export interface WalkinState {
  time_s: number;
  box_air_c: number;
  pan_core_c: number;
  peak_box_air_c: number;
  frost_mass_kg: number;
  frost_ua_factor: number;
  infiltration_w: number;
  pan_load_w: number;
  evap_cooling_w: number;
  box_excursion_41f: boolean;
  fitted_params: readonly string[];
}

export function createInitialWalkinState(): WalkinState {
  return {
    time_s: 0.0,
    box_air_c: 3.3,
    pan_core_c: 85.0,
    peak_box_air_c: 3.3,
    frost_mass_kg: 0.05,
    frost_ua_factor: 0.99,
    infiltration_w: 0.0,
    pan_load_w: 0.0,
    evap_cooling_w: 0.0,
    box_excursion_41f: false,
    fitted_params: [],
  };
}

export function stepWalkinKernel(state: WalkinState, cfg: WalkinConfig, dtS: number): WalkinState {
  if (!Number.isFinite(dtS) || dtS <= 0) return state;
  const nSub = clamp(Math.ceil(dtS / 5.0), 1, 240);
  const h = dtS / nSub;

  const pans = clamp(cfg.panCount, 1, 60);
  const depthIn = clamp(cfg.panDepthIn, 1.5, 6.0);
  const massPerPanKg = 1.9 * depthIn;
  const cFoodTotal = pans * massPerPanKg * 3650.0;
  const halfDepthM = 0.5 * depthIn * 0.0254;
  const biot = (18.0 * halfDepthM) / 0.52;
  const uaPerPanWK = (18.0 * 0.34) / (1.0 + biot / 2.2);
  const uaPansTotal = pans * uaPerPanWK;

  const cBoxJK = 195000.0;
  const uaWallsWK = 14.5;
  const condDerate = clamp(1.0 - 0.0115 * Math.max(0.0, cfg.kitchenAmbientC - 25.0), 0.68, 1.05);
  const ratedEvapW = cfg.ratedEvapW ?? 3800.0;

  for (let step = 0; step < nSub; step += 1) {
    const qPans = uaPansTotal * Math.max(0.0, state.pan_core_c - state.box_air_c);
    const qWalls = uaWallsWK * Math.max(0.0, cfg.kitchenAmbientC - state.box_air_c);
    const qInf = gosneyOlamaInfiltrationW(cfg.kitchenAmbientC, state.box_air_c, cfg.doorOpenFrac);

    const frostRateKgS = Math.max(0.0, 0.000045 * cfg.doorOpenFrac + 0.000000018 * qPans);
    state.frost_mass_kg = Math.min(4.5, state.frost_mass_kg + frostRateKgS * h);
    state.frost_ua_factor = clamp(1.0 / (1.0 + 0.18 * state.frost_mass_kg), 0.45, 1.0);

    const drive = clamp((state.box_air_c - 2.6) / 1.4, 0.0, 1.0);
    const qEvap = ratedEvapW * condDerate * state.frost_ua_factor * drive;

    state.pan_core_c = Math.max(state.box_air_c, state.pan_core_c - (qPans / cFoodTotal) * h);
    state.box_air_c = clamp(state.box_air_c + ((qPans + qWalls + qInf - qEvap) / cBoxJK) * h, 0.5, 45.0);
    if (state.box_air_c > state.peak_box_air_c) state.peak_box_air_c = state.box_air_c;
    if (state.box_air_c > HACCP_COLD_LIMIT_C) state.box_excursion_41f = true;
    state.infiltration_w = qInf;
    state.pan_load_w = qPans;
    state.evap_cooling_w = qEvap;
    state.time_s += h;
  }
  return state;
}

export function gordonTaylorTgC(wW: number): number {
  const w = clamp(wW, 0.0, 0.5);
  const wDry = 1.0 - w;
  return (wDry * 148.0 + GORDON_TAYLOR_K_FITTED * w * -135.0) / Math.max(1e-6, wDry + GORDON_TAYLOR_K_FITTED * w);
}

export interface CarrierConfig {
  ventFraction: number;
  foodMassKg?: number;
  ambientC?: number;
  initialFoodC?: number;
}

export interface CarrierState {
  time_s: number;
  food_temp_c: number;
  headspace_rh: number;
  crust_moisture_w: number;
  crust_tg_c: number;
  crispness_index: number;
  minutes_until_soggy: number | null;
  minutes_until_below_135f: number | null;
  fitted_params: readonly string[];
}

export interface CarrierSurrogateResult {
  minutes_until_soggy: number;
  minutes_until_below_135f: number;
  max_bound_min: number;
  out_of_domain: boolean;
}

export function createInitialCarrierState(cfg: CarrierConfig): CarrierState {
  const w0 = 0.045;
  return {
    time_s: 0.0,
    food_temp_c: clamp(cfg.initialFoodC ?? 84.0, 55.0, 98.0),
    headspace_rh: 0.36,
    crust_moisture_w: w0,
    crust_tg_c: gordonTaylorTgC(w0),
    crispness_index: 1.0,
    minutes_until_soggy: null,
    minutes_until_below_135f: null,
    fitted_params: ["gordon_taylor_k_gt"],
  };
}

export function stepCarrierKernel(state: CarrierState, cfg: CarrierConfig, dtS: number): CarrierState {
  if (!Number.isFinite(dtS) || dtS <= 0) return state;
  const nSub = clamp(Math.ceil(dtS / 5.0), 1, 360);
  const h = dtS / nSub;

  const vent = clamp(cfg.ventFraction, 0.0, 1.0);
  const massKg = clamp(cfg.foodMassKg ?? 4.5, 0.8, 15.0);
  const cFood = massKg * 3150.0;
  const uaEff = 2.6 + 4.4 * vent;
  const amb = cfg.ambientC ?? 22.0;

  for (let step = 0; step < nSub; step += 1) {
    const qLoss = uaEff * Math.max(0.0, state.food_temp_c - amb);
    state.food_temp_c = Math.max(amb, state.food_temp_c - (qLoss / cFood) * h);

    const rhEq = 0.97 - 0.50 * vent;
    const tauRhS = 140.0 + 90.0 * (1.0 - vent);
    state.headspace_rh = clamp(state.headspace_rh + ((rhEq - state.headspace_rh) / tauRhS) * h, 0.20, 0.99);

    const wEq = 0.042 + 0.165 * state.headspace_rh * state.headspace_rh;
    const tauWS = 520.0 + 380.0 * vent;
    state.crust_moisture_w = clamp(state.crust_moisture_w + ((wEq - state.crust_moisture_w) / tauWS) * h, 0.03, 0.30);

    state.crust_tg_c = gordonTaylorTgC(state.crust_moisture_w);
    const marginC = state.crust_tg_c - state.food_temp_c;
    state.crispness_index = clamp((marginC + 18.0) / 36.0, 0.0, 1.0);

    state.time_s += h;
    const mins = state.time_s / 60.0;
    if (state.minutes_until_soggy === null && state.crispness_index <= 0.35) {
      state.minutes_until_soggy = mins;
    }
    if (state.minutes_until_below_135f === null && state.food_temp_c < HOT_HOLD_FLOOR_C) {
      state.minutes_until_below_135f = mins;
    }
  }
  return state;
}

export function evaluateCarrierSurrogate(cfg: CarrierConfig): CarrierSurrogateResult {
  const amb = cfg.ambientC ?? 22.0;
  const massKg = cfg.foodMassKg ?? 4.5;
  const t0 = cfg.initialFoodC ?? 84.0;
  const outOfDomain = amb < 15.0 || amb > 32.0 || massKg < 2.0 || massKg > 8.0 || t0 < 75.0 || t0 > 92.0;
  const s = createInitialCarrierState(cfg);
  for (let i = 0; i < 360; i += 1) {
    stepCarrierKernel(s, cfg, 30.0);
    if (s.minutes_until_soggy !== null && s.minutes_until_below_135f !== null) break;
  }
  return {
    minutes_until_soggy: s.minutes_until_soggy ?? 180.0,
    minutes_until_below_135f: s.minutes_until_below_135f ?? 180.0,
    max_bound_min: outOfDomain ? 8.5 : 0.45,
    out_of_domain: outOfDomain,
  };
}

/**
 * Simple static biomechanics lumbar L5/S1 torque estimator for ergonomic pre-operational screening
 * (Color::Estimated/PreOperationalDesignScreen/NoClinicalOrOshaCertification).
 */
export function estimateLumbarTorqueNm(
  bodyMassKg: number,
  loadMassKg: number,
  reachM: number,
  torsoFlexDeg: number,
): number {
  const g = 9.80665;
  const rad = (clamp(torsoFlexDeg, 0, 60) * Math.PI) / 180;
  const upperBodyMoment = 0.46 * bodyMassKg * g * 0.24 * Math.sin(rad);
  const loadMoment = loadMassKg * g * clamp(reachM, 0.2, 0.95);
  return upperBodyMoment + loadMoment;
}
