// Typed layout of the 176-element f64 snapshot buffer produced by
// `PocketWatchSim::snapshot()` in `apps/pocket-watch/engine-watch/src/lib.rs`.

export const SNAPSHOT_HEADER_LEN = 32;
export const SNAPSHOT_NODE_STRIDE = 12;
export const NODE_COUNT = 12;

export interface NodeMetadata {
  readonly index: number;
  readonly shortName: string;
  readonly fullName: string;
  readonly materialLabel: string;
  readonly matdbPack: string;
  readonly meltingPointK: number;
}

export const WATCH_NODES_METADATA: readonly NodeMetadata[] = [
  {
    index: 0,
    shortName: "Case 12h (Au)",
    fullName: "18k Gold Case (12h Crown/Bow)",
    materialLabel: "750 yellow gold (Au-Ag-Cu)",
    matdbPack: "gold/silver/copper-pure-nbs-c447",
    meltingPointK: 1158.15,
  },
  {
    index: 1,
    shortName: "Case 9h (Au)",
    fullName: "18k Gold Case (9h Band)",
    materialLabel: "750 yellow gold (Au-Ag-Cu)",
    matdbPack: "gold/silver/copper-pure-nbs-c447",
    meltingPointK: 1158.15,
  },
  {
    index: 2,
    shortName: "Case 3h (Au)",
    fullName: "18k Gold Case (3h Band)",
    materialLabel: "750 yellow gold (Au-Ag-Cu)",
    matdbPack: "gold/silver/copper-pure-nbs-c447",
    meltingPointK: 1158.15,
  },
  {
    index: 3,
    shortName: "Case 6h (Au)",
    fullName: "18k Gold Case (6h Tourbillon Rim)",
    materialLabel: "750 yellow gold (Au-Ag-Cu)",
    matdbPack: "gold/silver/copper-pure-nbs-c447",
    meltingPointK: 1158.15,
  },
  {
    index: 4,
    shortName: "Crystals (Al₂O₃)",
    fullName: "Sapphire Exhibition Crystals",
    materialLabel: "Sapphire",
    matdbPack: "sapphire-al2o3-crystran-reference",
    meltingPointK: 2326.0,
  },
  {
    index: 5,
    shortName: "Mainplate (Brass)",
    fullName: "Brass Mainplate & Geneva Bridges",
    materialLabel: "C36000 Brass",
    matdbPack: "brass-c36000-cda-reference",
    meltingPointK: 1188.0,
  },
  {
    index: 6,
    shortName: "Barrel & Spring",
    fullName: "Mainspring Barrel & Steel Ribbon",
    materialLabel: "Music-wire steel",
    matdbPack: "steel-music-wire-string-heptonstall-ligo-p1400229",
    meltingPointK: 1725.0,
  },
  {
    index: 7,
    shortName: "Going Train",
    fullName: "Center, 3rd & 4th Wheels",
    materialLabel: "C36000 brass",
    matdbPack: "brass-c36000-cda-reference",
    meltingPointK: 1188.0,
  },
  {
    index: 8,
    shortName: "Tourbillon Cage",
    fullName: "1-Minute Tourbillon Carriage & Well",
    materialLabel: "Brass seat (Ti cage inertia)",
    matdbPack: "brass-c36000-cda-reference",
    meltingPointK: 1188.0,
  },
  {
    index: 9,
    shortName: "Escapement & Ruby",
    fullName: "Swiss Lever Fork, Escape Wheel & Rubies",
    materialLabel: "Steel + ruby",
    matdbPack: "steel-1018-onlinemetals-reference",
    meltingPointK: 1725.0,
  },
  {
    index: 10,
    shortName: "Balance Wheel",
    fullName: "Glucydur BeCu Rim & 8 Gold Timing Screws",
    materialLabel: "BeCu C17200",
    matdbPack: "beryllium-copper-c17200-ngk-reference",
    meltingPointK: 1255.0,
  },
  {
    index: 11,
    shortName: "Hairspring",
    fullName: "13-Turn Hairspring",
    materialLabel: "Carbon steel / Nivarox-type",
    matdbPack: "invar-fe36ni-nist-thermal-reference",
    meltingPointK: 1720.0,
  },
];

export interface NodeSnapshot {
  readonly index: number;
  readonly temperatureK: number;
  readonly enthalpyJKg: number;
  readonly meltFraction: number;
  readonly plasticSlump: number;
  readonly moltenGoldG: number;
  readonly frozenGoldG: number;
  readonly ln2FilmG: number;
  readonly frostThicknessUm: number;
  readonly oilFilmUm: number;
  readonly incandescenceRgb: [number, number, number];
}

export interface WatchSnapshot {
  readonly labTimeS: number;
  readonly watchDisplayedTimeS: number;
  readonly mainspringWind: number;
  readonly mainspringTorqueNm: number;
  readonly barrelArborAngleRad: number;
  readonly barrelWheelAngleRad: number;
  readonly centerWheelAngleRad: number;
  readonly thirdWheelAngleRad: number;
  readonly tourbillonCageAngleRad: number;
  readonly escapeWheelVisualRad: number;
  readonly palletForkAngleRad: number;
  readonly balanceAngleRad: number;
  readonly balanceOmegaRadS: number;
  readonly balanceAmplitudeDeg: number;
  readonly totalBeats: number;
  readonly measuredRateSPerDay: number;
  readonly theoreticalRateSPerDay: number;
  readonly measuredBeatErrorMs: number;
  readonly moonphaseAngleRad: number;
  readonly statusCode: number;
  readonly balanceInertiaRatio: number;
  readonly hairspringStiffnessRatio: number;
  readonly pivotFrictionMultiplier: number;
  readonly minPivotClearanceUm: number;
  readonly bFieldEscapementMt: number;
  readonly remanentMagnetizationMt: number;
  readonly goldFoulingJamFactor: number;
  readonly ln2BoilingPowerW: number;
  readonly blowtorchHeatingPowerW: number;
  readonly bluedSteelRgb: [number, number, number];
  readonly nodes: NodeSnapshot[];
}

export function decodeWatchSnapshot(raw: Float64Array | number[]): WatchSnapshot {
  const nodes: NodeSnapshot[] = [];
  for (let i = 0; i < NODE_COUNT; i++) {
    const base = SNAPSHOT_HEADER_LEN + i * SNAPSHOT_NODE_STRIDE;
    nodes.push({
      index: i,
      temperatureK: raw[base] ?? 293.15,
      enthalpyJKg: raw[base + 1] ?? 0,
      meltFraction: raw[base + 2] ?? 0,
      plasticSlump: raw[base + 3] ?? 0,
      moltenGoldG: raw[base + 4] ?? 0,
      frozenGoldG: raw[base + 5] ?? 0,
      ln2FilmG: raw[base + 6] ?? 0,
      frostThicknessUm: raw[base + 7] ?? 0,
      oilFilmUm: raw[base + 8] ?? 0,
      incandescenceRgb: [
        raw[base + 9] ?? 0,
        raw[base + 10] ?? 0,
        raw[base + 11] ?? 0,
      ],
    });
  }

  return {
    labTimeS: raw[0] ?? 0,
    watchDisplayedTimeS: raw[1] ?? 0,
    mainspringWind: raw[2] ?? 0.88,
    mainspringTorqueNm: raw[3] ?? 0.0082,
    barrelArborAngleRad: raw[4] ?? 0,
    barrelWheelAngleRad: raw[5] ?? 0,
    centerWheelAngleRad: raw[6] ?? 0,
    thirdWheelAngleRad: raw[7] ?? 0,
    tourbillonCageAngleRad: raw[8] ?? 0,
    escapeWheelVisualRad: raw[9] ?? 0,
    palletForkAngleRad: raw[10] ?? 0,
    balanceAngleRad: raw[11] ?? 0,
    balanceOmegaRadS: raw[12] ?? 0,
    balanceAmplitudeDeg: raw[13] ?? 285,
    totalBeats: raw[14] ?? 0,
    measuredRateSPerDay: raw[15] ?? 0,
    theoreticalRateSPerDay: raw[16] ?? 0,
    measuredBeatErrorMs: raw[17] ?? 0,
    moonphaseAngleRad: raw[18] ?? 0,
    statusCode: raw[19] ?? 0,
    balanceInertiaRatio: raw[20] ?? 1,
    hairspringStiffnessRatio: raw[21] ?? 1,
    pivotFrictionMultiplier: raw[22] ?? 1,
    minPivotClearanceUm: raw[23] ?? 4.5,
    bFieldEscapementMt: raw[24] ?? 0,
    remanentMagnetizationMt: raw[25] ?? 0,
    goldFoulingJamFactor: raw[26] ?? 0,
    ln2BoilingPowerW: raw[27] ?? 0,
    blowtorchHeatingPowerW: raw[28] ?? 0,
    bluedSteelRgb: [raw[29] ?? 0.25, raw[30] ?? 0.48, raw[31] ?? 0.85],
    nodes,
  };
}

export function formatHorologicalTime(secondsTotal: number): string {
  const s = ((secondsTotal % 43200) + 43200) % 43200;
  const hours = Math.floor(s / 3600) || 12;
  const mins = Math.floor((s % 3600) / 60);
  const secs = s % 60;
  const hh = String(hours).padStart(2, "0");
  const mm = String(mins).padStart(2, "0");
  const ss = secs.toFixed(1).padStart(4, "0");
  return `${hh}:${mm}:${ss}`;
}
