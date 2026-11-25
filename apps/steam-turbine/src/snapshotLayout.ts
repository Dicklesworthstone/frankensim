// Typed decoders for the flat WebAssembly ABI exported by fs-steam-turbine (lib.rs).

export const SNAPSHOT_LEN = 81;
export const BLADE_ROW_STRIDE = 9;
export const STAGE_ROW_STRIDE = 16;
export const NUM_STAGES = 17;
export const CONTINUUM_STRIDE = 8;
export const MOLLIER_PAIR_STRIDE = 2;
export const CAMPBELL_STRIDE = 6;

export const STATUS_LABELS: Record<number, { text: string; tone: "nominal" | "island" | "alarm" | "trip" }> = {
  0: { text: "SYNC NOMINAL", tone: "nominal" },
  1: { text: "ISLAND / RUN-UP", tone: "island" },
  2: { text: "ALARM", tone: "alarm" },
  3: { text: "TRIP: OVERSPEED", tone: "trip" },
  4: { text: "TRIP: THERMAL", tone: "trip" },
  5: { text: "TRIP: VIBRATION", tone: "trip" },
  6: { text: "TRIP: MANUAL ESV", tone: "trip" },
  7: { text: "TRIP: POLE SLIP", tone: "trip" },
};

/** Bit meanings of snapshot slot 80 (`alarmFlags`), in bit order. */
export const ALARM_FLAG_TEXT: readonly string[] = [
  "VIBRATION",
  "BEARING TEMP",
  "WHIRL UNSTABLE",
  "OVERSPEED TRIP",
  "THERMAL TRIP",
  "VIBRATION TRIP",
  "MANUAL TRIP",
  "POLE-SLIP TRIP",
  "DERATING",
];

/** Names of the alarms (bits 0..2) and derating (bit 8) that are set; trips have their own status. */
export function activeAlarmNames(flags: number): string[] {
  const out: string[] = [];
  for (const bit of [0, 1, 2, 8]) {
    if (flags & (1 << bit)) out.push(ALARM_FLAG_TEXT[bit]!);
  }
  return out;
}

/** Protection limits the engine applies (`SteamTurbineSim::limits`). */
export interface EngineLimits {
  vibZoneAbUm: number;
  vibAlarmUm: number;
  vibTripUm: number;
  bearingDerateC: number;
  bearingAlarmC: number;
  bearingTripC: number;
  overspeedTripRpm: number;
}

export function decodeLimits(f: Float64Array): EngineLimits {
  const g = (i: number): number => f[i] ?? Number.NaN;
  return {
    vibZoneAbUm: g(0),
    vibAlarmUm: g(1),
    vibTripUm: g(2),
    bearingDerateC: g(3),
    bearingAlarmC: g(4),
    bearingTripC: g(5),
    overspeedTripRpm: g(6),
  };
}

export interface TurbineSnapshot {
  timeS: number;
  steps: number;
  rpm: number;
  thetaRad: number;
  omegaRadS: number;
  elecMw: number;
  shaftMw: number;
  elecEffPct: number;
  chpEffPct: number;
  statusCode: number;

  liveSteamBar: number;
  liveSteamC: number;
  liveFlowKgS: number;
  chestBar: number;
  chestC: number;
  extHeaderBar: number;
  extHeaderC: number;
  extFlowKgS: number;
  lpFlowKgS: number;
  condKpa: number;
  condSatC: number;
  lpQualityX: number;
  hpEtaPct: number;
  iplpEtaPct: number;
  hpMw: number;
  iplpMw: number;
  axialThrustKn: number;
  windageKw: number;

  condDutyMw: number;
  hCondWm2K: number;
  hCwWm2K: number;
  uCondWm2K: number;
  cwOutC: number;
  condAirPct: number;

  cvLift: number;
  esvOpen: number;
  extLift: number;
  deratingFrac: number;

  loadAngleDeg: number;
  statorCurrentA: number;
  reactiveMvar: number;
  genLossKw: number;
  breakerClosed: boolean;

  frontVibUmPp: number;
  rearVibUmPp: number;
  sealExcursionUm: number;
  rubForceKn: number;
  spectralAbscissaPerS: number;
  dominantWhirlHz: number;
  crit1Rpm: number;
  crit2Rpm: number;
  alfordKxyMnM: number;
  isoLogE: number;
  rotorCenterXUm: number;
  rotorCenterYUm: number;

  frontFilmUm: number;
  rearFilmUm: number;
  frontEccentricity: number;
  oilViscMpaS: number;
  bearingHeatKw: number;

  tHpRotorC: number;
  tHpCasingC: number;
  tIpRotorC: number;
  tIpCasingC: number;
  tLpHoodC: number;
  tFrontBrgC: number;
  tRearBrgC: number;
  tLubeOilC: number;
  diffExpansionMm: number;
  hpThermalStressMpa: number;
  hpCreepRatePerS: number;

  ledgerEnergyResJ: number;
  ledgerMassResKg: number;
  ledgerEntropyGenJK: number;
  ledgerBoundaryMj: number;
  lastNewtonIters: number;

  splDb: number;
  loudnessSone: number;
  sharpnessAcum: number;
  bpfLpHz: number;
  /** One bit per alarm and trip (see `ALARM_FLAG_TEXT`). */
  alarmFlags: number;
}

export function decodeSnapshot(s: Float64Array): TurbineSnapshot {
  const g = (i: number): number => s[i] ?? 0;
  return {
    timeS: g(0),
    steps: g(1),
    rpm: g(2),
    thetaRad: g(3),
    omegaRadS: g(4),
    elecMw: g(5),
    shaftMw: g(6),
    elecEffPct: g(7),
    chpEffPct: g(8),
    statusCode: Math.round(g(9)),

    liveSteamBar: g(10),
    liveSteamC: g(11),
    liveFlowKgS: g(12),
    chestBar: g(13),
    chestC: g(14),
    extHeaderBar: g(15),
    extHeaderC: g(16),
    extFlowKgS: g(17),
    lpFlowKgS: g(18),
    condKpa: g(19),
    condSatC: g(20),
    lpQualityX: g(21),
    hpEtaPct: g(22),
    iplpEtaPct: g(23),
    hpMw: g(24),
    iplpMw: g(25),
    axialThrustKn: g(26),
    windageKw: g(27),

    condDutyMw: g(28),
    hCondWm2K: g(29),
    hCwWm2K: g(30),
    uCondWm2K: g(31),
    cwOutC: g(32),
    condAirPct: g(33),

    cvLift: g(34),
    esvOpen: g(35),
    extLift: g(36),
    deratingFrac: g(37),

    loadAngleDeg: g(38),
    statorCurrentA: g(39),
    reactiveMvar: g(40),
    genLossKw: g(41),
    breakerClosed: g(42) >= 0.5,

    frontVibUmPp: g(43),
    rearVibUmPp: g(44),
    sealExcursionUm: g(45),
    rubForceKn: g(46),
    spectralAbscissaPerS: g(47),
    dominantWhirlHz: g(48),
    crit1Rpm: g(49),
    crit2Rpm: g(50),
    alfordKxyMnM: g(51),
    isoLogE: g(52),
    rotorCenterXUm: g(53),
    rotorCenterYUm: g(54),

    frontFilmUm: g(55),
    rearFilmUm: g(56),
    frontEccentricity: g(57),
    oilViscMpaS: g(58),
    bearingHeatKw: g(59),

    tHpRotorC: g(60),
    tHpCasingC: g(61),
    tIpRotorC: g(62),
    tIpCasingC: g(63),
    tLpHoodC: g(64),
    tFrontBrgC: g(65),
    tRearBrgC: g(66),
    tLubeOilC: g(67),
    diffExpansionMm: g(68),
    hpThermalStressMpa: g(69),
    hpCreepRatePerS: g(70),

    ledgerEnergyResJ: g(71),
    ledgerMassResKg: g(72),
    ledgerEntropyGenJK: g(73),
    ledgerBoundaryMj: g(74),
    lastNewtonIters: g(75),

    splDb: g(76),
    loudnessSone: g(77),
    sharpnessAcum: g(78),
    bpfLpHz: g(79),
    alarmFlags: g(80),
  };
}

export interface BladeRowSpec {
  stageIndex: number;
  zCenterM: number;
  rHubM: number;
  rTipM: number;
  axialChordM: number;
  bladeCount: number;
  airfoilId: number; // 0..6 -> maps to mesh_id 17 + airfoilId
  isRotating: boolean;
  hasShroud: boolean;
}

export function decodeBladeRows(flat: Float64Array): BladeRowSpec[] {
  const rows: BladeRowSpec[] = [];
  for (let i = 0; i + BLADE_ROW_STRIDE <= flat.length; i += BLADE_ROW_STRIDE) {
    rows.push({
      stageIndex: Math.round(flat[i] ?? 0),
      airfoilId: Math.round(flat[i + 1] ?? 0),
      isRotating: (flat[i + 2] ?? 1) >= 0.5,
      zCenterM: flat[i + 3] ?? 0,
      rHubM: flat[i + 4] ?? 0.2,
      rTipM: flat[i + 5] ?? 0.3,
      bladeCount: Math.round(flat[i + 6] ?? 32),
      axialChordM: flat[i + 7] ?? 0.04,
      hasShroud: (flat[i + 8] ?? 0) >= 0.5,
    });
  }
  return rows;
}

export interface StageTelemetry {
  stageIndex: number;
  zM: number;
  rHubM: number;
  rTipM: number;
  pOutBar: number;
  tOutC: number;
  qualityX: number;
  machMax: number;
  etaPct: number;
  powerMw: number;
  cAxialMS: number;
  cSwirlInMS: number;
  cSwirlOutMS: number;
  reaction: number;
  rootCfMpa: number;
  rootVmMpa: number;
  axialThrustKn: number;
}

export function decodeStageProfile(flat: Float64Array): StageTelemetry[] {
  const out: StageTelemetry[] = [];
  for (let i = 0; i + STAGE_ROW_STRIDE <= flat.length; i += STAGE_ROW_STRIDE) {
    out.push({
      stageIndex: out.length,
      zM: flat[i] ?? 0,
      rHubM: flat[i + 1] ?? 0.28,
      rTipM: flat[i + 2] ?? 0.4,
      pOutBar: flat[i + 3] ?? 1,
      tOutC: flat[i + 4] ?? 200,
      qualityX: flat[i + 5] ?? 1,
      machMax: flat[i + 6] ?? 0.3,
      etaPct: flat[i + 7] ?? 85,
      powerMw: flat[i + 8] ?? 2.5,
      cAxialMS: flat[i + 9] ?? 60,
      cSwirlInMS: flat[i + 10] ?? 180,
      cSwirlOutMS: flat[i + 11] ?? 10,
      reaction: flat[i + 12] ?? 0.48,
      rootCfMpa: flat[i + 13] ?? 120,
      rootVmMpa: flat[i + 14] ?? 160,
      axialThrustKn: flat[i + 15] ?? 20,
    });
  }
  return out;
}

export interface GeometryReceipt {
  shaftMassKg: number;
  bladesMassKg: number;
  totalMassKg: number;
  comZM: number;
  polarInertiaKgM2: number;
  transverseInertiaKgM2: number;
  shaftSurfaceAreaM2: number;
  rotorVertexCount: number;
  rotorTriangleCount: number;
  rotorHausdorffM: number;
  rotorEulerChar: number;
  rotorWatertight: boolean;
  totalBaseTriangles: number;
  totalSceneTriangles: number;
}

export function decodeGeometryReceipt(r: Float64Array): GeometryReceipt {
  const g = (i: number): number => r[i] ?? 0;
  return {
    shaftMassKg: g(0),
    bladesMassKg: g(1),
    totalMassKg: g(2),
    comZM: g(3),
    polarInertiaKgM2: g(4),
    transverseInertiaKgM2: g(5),
    shaftSurfaceAreaM2: g(6),
    rotorVertexCount: Math.round(g(7)),
    rotorTriangleCount: Math.round(g(8)),
    rotorHausdorffM: g(9),
    rotorEulerChar: Math.round(g(10)),
    rotorWatertight: g(11) >= 0.5,
    totalBaseTriangles: Math.round(g(12)),
    totalSceneTriangles: Math.round(g(13)),
  };
}

export interface FemShellCertificate {
  f0Hz: number;
  f1Hz: number;
  f2Hz: number;
  stiffnessTraceNM: number;
  diffuserRadResistance: number;
  certifiedFromFem: boolean;
}

export function decodeFemCertificate(f: Float64Array): FemShellCertificate {
  const g = (i: number): number => f[i] ?? 0;
  return {
    f0Hz: g(0),
    f1Hz: g(1),
    f2Hz: g(2),
    stiffnessTraceNM: g(3),
    diffuserRadResistance: g(4),
    certifiedFromFem: g(5) >= 0.5,
  };
}

export interface DesignParams {
  tipClearanceMm: number;
  hpReaction: number;
  lpLastBladeHM: number;
  bearingClearanceUm: number;
  bearingPreload: number;
  casingWallMm: number;
}

export function decodeDesignParams(d: Float64Array): DesignParams {
  const g = (i: number, fb: number): number => d[i] ?? fb;
  return {
    tipClearanceMm: g(0, 0.75),
    hpReaction: g(1, 0.48),
    lpLastBladeHM: g(2, 0.527),
    bearingClearanceUm: g(3, 220),
    bearingPreload: g(4, 0.32),
    casingWallMm: g(5, 95),
  };
}

export interface DesignMetricsDto {
  elecEffPct: number;
  chpEffPct: number;
  elecPowerMw: number;
  maxBladeStressMpa: number;
  casingStressMpa: number;
  spectralAbscissaPerS: number;
  bearingLossKw: number;
  lpLeavingLossKw: number;
  minQualityX: number;
  minGoodmanSf: number;
}

export interface DesignOptimizationReport {
  optimizedParams: DesignParams;
  baseline: DesignMetricsDto;
  optimized: DesignMetricsDto;
  totalEvals: number;
  cmaGenerations: number;
  historyObj: number[];
}

export function decodeOptReport(flat: Float64Array): DesignOptimizationReport {
  const g = (i: number): number => flat[i] ?? 0;
  const parseMetrics = (off: number): DesignMetricsDto => ({
    elecEffPct: g(off),
    chpEffPct: g(off + 1),
    elecPowerMw: g(off + 2),
    maxBladeStressMpa: g(off + 3),
    casingStressMpa: g(off + 4),
    spectralAbscissaPerS: g(off + 5),
    bearingLossKw: g(off + 6),
    lpLeavingLossKw: g(off + 7),
    minQualityX: g(off + 8),
    minGoodmanSf: g(off + 9),
  });
  const historyObj: number[] = [];
  for (let i = 28; i < flat.length; i++) historyObj.push(flat[i] ?? 0);
  return {
    optimizedParams: {
      tipClearanceMm: g(0),
      hpReaction: g(1),
      lpLastBladeHM: g(2),
      bearingClearanceUm: g(3),
      bearingPreload: g(4),
      casingWallMm: g(5),
    },
    baseline: parseMetrics(6),
    optimized: parseMetrics(16),
    totalEvals: Math.round(g(26)),
    cmaGenerations: Math.round(g(27)),
    historyObj,
  };
}
