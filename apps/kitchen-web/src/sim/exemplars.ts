/**
 * Canonical 5 Commercial Kitchen Archetypes & 67-Asset Equipment Catalog metadata.
 *
 * All dimensions, fixture counts, triangle counts, electrical service ratings, and
 * aisle clearances come directly from the verified `.scene3d.v1.json` and
 * `.site_metadata.json` manifests under `assets/kitchens/exemplars/`.
 *
 * Zero commercial OEM brand names (Law V-L16 clean).
 */

export type ExemplarId =
  | "mall_qsr_400"
  | "freestanding_qsr_1450"
  | "fast_casual_1100"
  | "fine_dining_950"
  | "commissary_catering_3200";

export type CameraPresetId = "isometric" | "cookline" | "pass" | "hood" | "walkin";

export interface CameraPose {
  readonly pos: readonly [number, number, number]; // Three.js Y-up meters
  readonly target: readonly [number, number, number];
}

export interface StaffWaypointRoute {
  readonly cohort: "p05" | "p50" | "p95";
  readonly role: string;
  /** Floor path in the manifests' Z-up [x_m, y_m] coordinates (the hero maps y to three.js -z). */
  readonly waypoints: readonly (readonly [number, number])[];
  readonly speedMps: number;
  readonly accentHex: string;
}

export interface ExemplarSpec {
  readonly id: ExemplarId;
  readonly shortName: string;
  readonly fullName: string;
  readonly areaSqFt: number;
  readonly roomM: readonly [number, number, number]; // [width, depth, ceiling] in Blender meters
  readonly fixtures: number;
  readonly triangles: number;
  readonly staffCount: number;
  readonly minAisleClearM: number;
  readonly minAisleNote: string;
  readonly baseGasThermalKw: number;
  readonly baseElectricKw: number;
  readonly baseExhaustCfm: number;
  readonly electricalService: string;
  readonly nominalOrdersPerHr: number;
  readonly maxOrdersPerHr: number;
  readonly tagline: string;
  readonly cameras: Record<CameraPresetId, CameraPose>;
  readonly staffRoutes: readonly StaffWaypointRoute[];
}

export const EXEMPLARS: readonly ExemplarSpec[] = [
  {
    id: "mall_qsr_400",
    shortName: "Mall counter",
    fullName: "Mall Food-Court QSR · 400 sq ft",
    areaSqFt: 400,
    roomM: [4.88, 7.62, 2.743],
    fixtures: 14,
    triangles: 5908,
    staffCount: 3,
    minAisleClearM: 0.965,
    minAisleNote: '38.0" cookline pass (tight door-conflict pair)',
    baseGasThermalKw: 48.0,
    baseElectricKw: 26.5,
    baseExhaustCfm: 1870,
    electricalService: "200A · 208V 3φ",
    nominalOrdersPerHr: 75,
    maxOrdersPerHr: 140,
    tagline: "Single-aisle high-density food-court pod with tight 38\" cookline door overlap.",
    cameras: {
      isometric: { pos: [7.6, 5.8, 2.4], target: [1.35, 0.9, -3.5] },
      cookline: { pos: [2.6, 1.75, -0.4], target: [0.8, 1.05, -2.8] },
      pass: { pos: [3.2, 1.8, -6.9], target: [1.8, 1.0, -4.0] },
      hood: { pos: [3.8, 2.45, -1.2], target: [0.6, 2.1, -2.5] },
      walkin: { pos: [2.1, 1.9, -1.6], target: [4.4, 1.1, -2.2] },
    },
    staffRoutes: [
      { cohort: "p50", role: "Fry / Griddle", waypoints: [[1.28, 1.55], [1.28, 2.85], [1.28, 4.15], [1.28, 2.2]], speedMps: 0.52, accentHex: "#f59e0b" },
      { cohort: "p05", role: "Expo / Makeline", waypoints: [[3.1, 5.6], [1.5, 5.6], [1.45, 4.2], [3.1, 4.4]], speedMps: 0.56, accentHex: "#38bdf8" },
      { cohort: "p95", role: "Cold / Dish", waypoints: [[3.25, 1.8], [3.25, 4.6], [3.25, 3.0]], speedMps: 0.44, accentHex: "#34d399" },
    ],
  },
  {
    id: "freestanding_qsr_1450",
    shortName: "Drive-thru",
    fullName: "Freestanding Drive-Thru QSR · 1,450 sq ft",
    areaSqFt: 1450,
    roomM: [9.754, 13.716, 3.048],
    fixtures: 22,
    triangles: 11344,
    staffCount: 6,
    minAisleClearM: 1.118,
    minAisleNote: '44.0" dual cook aisles · 48.0" primary corridor',
    baseGasThermalKw: 118.0,
    baseElectricKw: 54.0,
    baseExhaustCfm: 5440,
    electricalService: "400A · 208V 3φ",
    nominalOrdersPerHr: 165,
    maxOrdersPerHr: 280,
    tagline: "Dual mirrored cooklines feeding drive-thru expo with dedicated walk-in cooler & freezer.",
    cameras: {
      isometric: { pos: [14.2, 10.2, 4.4], target: [3.2, 1.0, -6.5] },
      cookline: { pos: [3.4, 2.0, -1.0], target: [2.2, 1.1, -5.6] },
      pass: { pos: [5.5, 2.1, -12.6], target: [2.8, 1.1, -8.5] },
      hood: { pos: [6.5, 2.8, -2.5], target: [2.4, 2.2, -5.5] },
      walkin: { pos: [5.8, 2.3, -2.0], target: [8.4, 1.2, -3.5] },
    },
    staffRoutes: [
      { cohort: "p50", role: "West Fry Line", waypoints: [[1.85, 2.8], [1.85, 6.8], [1.85, 4.5]], speedMps: 0.62, accentHex: "#f59e0b" },
      { cohort: "p95", role: "East Grill Line", waypoints: [[4.7, 2.8], [4.7, 7.0], [4.7, 4.9]], speedMps: 0.58, accentHex: "#fb923c" },
      { cohort: "p05", role: "Drive-Thru Expo", waypoints: [[1.6, 9.2], [4.4, 9.2], [4.4, 11.4], [1.6, 11.0]], speedMps: 0.66, accentHex: "#38bdf8" },
      { cohort: "p50", role: "Walk-In Runner", waypoints: [[6.65, 1.5], [6.65, 10.5], [6.65, 5.2]], speedMps: 0.55, accentHex: "#34d399" },
    ],
  },
  {
    id: "fast_casual_1100",
    shortName: "Fast casual",
    fullName: "Exhibition Fast-Casual · 1,100 sq ft",
    areaSqFt: 1100,
    roomM: [7.315, 14.021, 3.048],
    fixtures: 24,
    triangles: 9956,
    staffCount: 5,
    minAisleClearM: 1.067,
    minAisleNote: '42.0" exhibition grill-to-assembly pass',
    baseGasThermalKw: 86.0,
    baseElectricKw: 44.5,
    baseExhaustCfm: 3680,
    electricalService: "300A · 208V 3φ",
    nominalOrdersPerHr: 120,
    maxOrdersPerHr: 210,
    tagline: "Open exhibition charbroiler & flattop theater fronting high-velocity cold-rail assembly.",
    cameras: {
      isometric: { pos: [12.8, 9.2, 3.6], target: [2.2, 1.0, -6.8] },
      cookline: { pos: [3.8, 1.9, -1.6], target: [1.3, 1.1, -5.2] },
      pass: { pos: [6.4, 2.0, -11.5], target: [3.2, 1.1, -8.0] },
      hood: { pos: [5.4, 2.7, -3.0], target: [1.2, 2.1, -5.6] },
      walkin: { pos: [3.8, 2.1, -9.8], target: [6.1, 1.2, -12.0] },
    },
    staffRoutes: [
      { cohort: "p95", role: "Exhibition Grill", waypoints: [[2.1, 2.5], [2.1, 6.2], [2.1, 4.1]], speedMps: 0.55, accentHex: "#f97316" },
      { cohort: "p05", role: "Assembly Line", waypoints: [[4.1, 3.0], [4.1, 8.2], [4.1, 5.4]], speedMps: 0.64, accentHex: "#22d3ee" },
      { cohort: "p50", role: "BOH Prep / Dish", waypoints: [[3.2, 10.2], [5.4, 10.2], [5.4, 12.2], [2.4, 11.8]], speedMps: 0.50, accentHex: "#34d399" },
    ],
  },
  {
    id: "fine_dining_950",
    shortName: "Fine dining",
    fullName: "European-Island Brigade · 950 sq ft",
    areaSqFt: 950,
    roomM: [6.706, 13.106, 3.048],
    fixtures: 26,
    triangles: 10140,
    staffCount: 6,
    minAisleClearM: 1.092,
    minAisleNote: '43.0" perimeter around central 144" Euro island',
    baseGasThermalKw: 74.0,
    baseElectricKw: 62.0,
    baseExhaustCfm: 4200,
    electricalService: "400A · 208V 3φ",
    nominalOrdersPerHr: 55,
    maxOrdersPerHr: 95,
    tagline: "Center-suite cooking island under 144\" island canopy with garde-manger, pastry, and pass.",
    cameras: {
      isometric: { pos: [11.6, 8.4, 3.2], target: [2.0, 1.0, -6.2] },
      cookline: { pos: [5.6, 1.9, -2.4], target: [3.1, 1.0, -5.2] },
      pass: { pos: [3.3, 2.0, -10.8], target: [3.3, 1.05, -7.2] },
      hood: { pos: [5.8, 2.75, -3.0], target: [3.2, 2.1, -5.4] },
      walkin: { pos: [4.2, 2.0, -8.5], target: [1.4, 1.1, -11.2] },
    },
    staffRoutes: [
      { cohort: "p50", role: "Saucier (South Island)", waypoints: [[1.5, 3.8], [4.9, 3.8], [3.2, 3.8]], speedMps: 0.48, accentHex: "#f59e0b" },
      { cohort: "p95", role: "Rôtisseur (West Island)", waypoints: [[1.4, 3.8], [1.4, 7.6], [1.4, 5.4]], speedMps: 0.50, accentHex: "#eab308" },
      { cohort: "p05", role: "Garde-Manger / Pass", waypoints: [[5.1, 3.9], [5.1, 6.8], [3.0, 7.8]], speedMps: 0.54, accentHex: "#c084fc" },
    ],
  },
  {
    id: "commissary_catering_3200",
    shortName: "Commissary",
    fullName: "Production Catering Commissary · 3,200 sq ft",
    areaSqFt: 3200,
    roomM: [15.24, 19.507, 3.658],
    fixtures: 72,
    triangles: 27480,
    staffCount: 11,
    minAisleClearM: 1.524,
    minAisleNote: '60.0" speed-rack spine · 69.5 mm door margin',
    baseGasThermalKw: 195.0,
    baseElectricKw: 233.33,
    baseExhaustCfm: 9800,
    electricalService: "600A · 480V/208V 3φ",
    nominalOrdersPerHr: 420,
    maxOrdersPerHr: 750,
    tagline: "72-fixture high-volume central kitchen: 20-pan roll-in combis, 40-gal kettles, blast chillers & speed-rack spine.",
    cameras: {
      isometric: { pos: [22.8, 15.2, 6.8], target: [5.2, 1.1, -9.5] },
      cookline: { pos: [11.2, 2.5, -4.5], target: [13.5, 1.2, -10.5] },
      pass: { pos: [8.4, 2.4, -17.5], target: [7.6, 1.2, -11.5] },
      hood: { pos: [8.8, 3.2, -6.5], target: [13.2, 2.4, -9.8] },
      walkin: { pos: [6.8, 2.5, -2.5], target: [2.5, 1.3, -4.2] },
    },
    staffRoutes: [
      { cohort: "p95", role: "Hot Production Bank", waypoints: [[11.6, 3.5], [11.6, 12.2], [11.6, 7.6]], speedMps: 0.60, accentHex: "#f59e0b" },
      { cohort: "p50", role: "Speed-Rack Spine", waypoints: [[7.8, 2.5], [7.8, 16.5], [7.8, 9.2]], speedMps: 0.68, accentHex: "#38bdf8" },
      { cohort: "p05", role: "Cold Chain / Blast Chill", waypoints: [[4.5, 2.8], [4.5, 11.0], [4.5, 6.2]], speedMps: 0.56, accentHex: "#34d399" },
      { cohort: "p50", role: "Warewash / Dish Pit", waypoints: [[10.8, 16.8], [13.8, 16.8], [12.0, 15.2]], speedMps: 0.50, accentHex: "#c084fc" },
    ],
  },
];

export interface HeroMetrics {
  readonly thermalKw: number;
  readonly electricalKw: number;
  readonly hvacCoolingKw: number;
  readonly hoodCapturePct: number;
  readonly exhaustCfm: number;
  readonly minAisleIn: number;
  readonly minAisleM: number;
  readonly thermalComfortC: number;
  readonly gasSharePct: number;
  readonly elecSharePct: number;
  readonly hvacSharePct: number;
}

/**
 * Compute live building-level energy, ventilation, and comfort tally for the hero & fleet matrix.
 * Pure deterministic first-principles energy balance:
 *   - Electrification converts gas thermal input (~42% pan efficiency) to induction (~86% efficiency),
 *     dramatically lowering sensible/radiant heat spill into the kitchen and HVAC cooling demand!
 */
export function computeKitchenTelemetry(
  spec: ExemplarSpec,
  ordersPerHr: number,
  electrificationPct: number,
  ambientC: number,
): HeroMetrics {
  const loadRatio = Math.max(0.25, Math.min(1.6, ordersPerHr / spec.nominalOrdersPerHr));
  const elecFrac = Math.max(0, Math.min(1, electrificationPct / 100));

  const gasInputKw = spec.baseGasThermalKw * (1 - elecFrac) * (0.45 + 0.55 * loadRatio);
  // Induction delivers equal pan flux at ~0.52x the primary burner input (86% vs 42% coupling)
  const convertedIndKw = spec.baseGasThermalKw * elecFrac * 0.52 * (0.45 + 0.55 * loadRatio);
  const baseElecKw = spec.baseElectricKw * (0.65 + 0.35 * loadRatio);
  const electricalKw = baseElecKw + convertedIndKw;

  // Convective/radiant heat spilled to room (gas ~58% spill; induction ~14% spill; ambient derates refrig)
  const roomSpillKw =
    gasInputKw * 0.58 +
    convertedIndKw * 0.14 +
    baseElecKw * 0.28 * (1 + Math.max(0, ambientC - 24) * 0.014);

  // Hood exhaust CFM reduces smoothly as gas plume load drops under electrification
  const exhaustCfm = spec.baseExhaustCfm * (0.72 + 0.28 * (1 - elecFrac * 0.42)) * (0.85 + 0.15 * loadRatio);
  const captureBase = 96.4 - loadRatio * 3.1 + elecFrac * 2.4 - Math.max(0, ambientC - 28) * 0.22;
  const hoodCapturePct = Math.max(82.0, Math.min(99.4, captureBase));

  // Uncaptured plume + envelope gain determines HVAC cooling load
  const uncapturedFrac = 1 - hoodCapturePct / 100;
  const hvacCoolingKw = roomSpillKw * (0.38 + uncapturedFrac * 2.2) + (spec.areaSqFt / 400) * Math.max(0, ambientC - 20) * 0.65;

  const totalBar = gasInputKw + electricalKw + hvacCoolingKw + 1e-6;
  const thermalComfortC = ambientC + (roomSpillKw * (0.45 + uncapturedFrac)) / Math.max(12, spec.areaSqFt * 0.028);

  return {
    thermalKw: gasInputKw + electricalKw,
    electricalKw,
    hvacCoolingKw,
    hoodCapturePct,
    exhaustCfm,
    minAisleIn: spec.minAisleClearM / 0.0254,
    minAisleM: spec.minAisleClearM,
    thermalComfortC: Math.min(46.5, thermalComfortC),
    gasSharePct: (gasInputKw / totalBar) * 100,
    elecSharePct: (electricalKw / totalBar) * 100,
    hvacSharePct: (hvacCoolingKw / totalBar) * 100,
  };
}

export type EquipmentTabId = "hot_line" | "ovens_fryers" | "cold_chain" | "ventilation" | "cookware_food";

export interface FlagshipEquipmentItem {
  readonly id: string;
  readonly label: string;
  readonly tab: EquipmentTabId;
  readonly thermalNote: string;
  readonly nominalTempC: number;
}

export const FLAGSHIP_EQUIPMENT: readonly FlagshipEquipmentItem[] = [
  // Hot Line
  { id: "NSF-RANGE-GAS-6B", label: "6-Burner Gas Range", tab: "hot_line", thermalNote: "Open flame grates · 42% pan coupling", nominalTempC: 310 },
  { id: "NSF-INDUCTION-SUITE-4Z", label: "4-Zone Induction Suite", tab: "hot_line", thermalNote: "Eddy-current coil · 86% coupling", nominalTempC: 245 },
  { id: "NSF-GRIDDLE-48", label: '48" Flat-Top Griddle', tab: "hot_line", thermalNote: "25 mm carbon steel platen", nominalTempC: 205 },
  { id: "NSF-CHARBROILER-36", label: '36" Radiant Charbroiler', tab: "hot_line", thermalNote: "High radiant IR + grease plume", nominalTempC: 340 },
  { id: "NSF-TILT-SKILLET-40GAL", label: "40-Gal Tilt Braising Pan", tab: "hot_line", thermalNote: "High-volume sear & braise deck", nominalTempC: 185 },
  // Ovens & Fryers
  { id: "NSF-COMBI-10PAN", label: "10-Pan Combi Oven", tab: "ovens_fryers", thermalNote: "Steam + convection + radiant cavity", nominalTempC: 235 },
  { id: "NSF-FRYER-TWIN-50LB", label: "Twin 50-lb Deep Fryer", tab: "ovens_fryers", thermalNote: "Two-zone cold-well oil vat · 180 °C", nominalTempC: 180 },
  { id: "NSF-COMBI-20PAN-ROLLIN", label: "20-Pan Roll-In Combi", tab: "ovens_fryers", thermalNote: "Commissary mobile speed-rack cavity", nominalTempC: 230 },
  { id: "NSF-KETTLE-40GAL", label: "40-Gal Steam Kettle", tab: "ovens_fryers", thermalNote: "Pressurized steam jacket stockpot", nominalTempC: 105 },
  // Cold Chain
  { id: "NSF-WALKIN-COOLER", label: "Walk-In Cooler Box", tab: "cold_chain", thermalNote: "3.3 °C box · Gosney-Olama door load", nominalTempC: 3.3 },
  { id: "NSF-BLAST-CHILLER-15PAN", label: "15-Pan Blast Chiller", tab: "cold_chain", thermalNote: "Rapid 135→41 °F HACCP drawdown", nominalTempC: -4.0 },
  { id: "NSF-SPEEDRACK-20PAN", label: "20-Pan Speed Rack", tab: "cold_chain", thermalNote: "Mobile aluminum pan transport frame", nominalTempC: 22 },
  { id: "NSF-CARRIER-INSULATED", label: "Insulated Food Carrier", tab: "cold_chain", thermalNote: "Hot-hold crust Tg & vent transport", nominalTempC: 76 },
  // Ventilation
  { id: "NSF-HOOD-TYPE1-CANOPY", label: "Type-I Grease Canopy Hood", tab: "ventilation", thermalNote: "Baffle grease capture + Heskestad plume", nominalTempC: 68 },
  { id: "NSF-HOOD-TYPE2-CONDENSATE", label: "Type-II Condensate Hood", tab: "ventilation", thermalNote: "Steam & warewash vapor exhaust", nominalTempC: 54 },
  { id: "NSF-PREPTABLE-72", label: '72" Stainless Worktable', tab: "ventilation", thermalNote: "304 brushed stainless prep deck", nominalTempC: 22 },
  // Cookware & Food
  { id: "NSF-SKILLET-12IN", label: '12" Cast-Iron Skillet', tab: "cookware_food", thermalNote: "High thermal mass · 2.65 kg seasoned", nominalTempC: 215 },
  { id: "NSF-STOCKPOT-40QT", label: "40-Qt Heavy Stockpot", tab: "cookware_food", thermalNote: "Aqueous boil & Bikerman foam deck", nominalTempC: 100 },
  { id: "NSF-CHEF-KNIFE-10IN", label: '10" Forged Chef Knife', tab: "cookware_food", thermalNote: "High-carbon blade cutting wedge", nominalTempC: 22 },
  { id: "ARCH-FOOD-STEAK-RIBEYE", label: "Ribeye Steak (Morph)", tab: "cookware_food", thermalNote: "8-slice conductive sear + Maillard", nominalTempC: 62 },
  { id: "ARCH-FOOD-SOURDOUGH-BOULE", label: "Sourdough Boule (Morph)", tab: "cookware_food", thermalNote: "Yeast oven-spring + Gordon-Taylor crust", nominalTempC: 96 },
  { id: "NSF-DEGLAZE-BOTTLE", label: "120 mL Deglaze Bottle", tab: "cookware_food", thermalNote: "Ethanol-water NRTL flash vaporizer", nominalTempC: 18 },
];

export interface CatalogEntry {
  readonly archetype_id: string;
  readonly asset_class: string;
  readonly family: string;
  readonly glb_path: string;
  readonly glb_bytes: number;
  readonly measured_extents_m: readonly [number, number, number];
  readonly triangle_count: number;
  readonly joint_count: number;
  readonly socket_ids: readonly string[];
  readonly matdb_card: string;
  readonly morph_targets: readonly string[];
}

export interface CatalogIndex {
  readonly asset_count: number;
  readonly assets: readonly CatalogEntry[];
}

export function categorizeCatalogEntry(entry: CatalogEntry): EquipmentTabId {
  const id = entry.archetype_id;
  if (entry.asset_class === "food" || entry.family === "cookware_and_food" || id.includes("SKILLET-8") || id.includes("SKILLET-10") || id.includes("SKILLET-12") || id.includes("SKILLET-14") || id.includes("STOCKPOT") || id.includes("HOTEL-PAN") || id.includes("KNIFE") || id.includes("BOARD") || id.includes("TONGS") || id.includes("SPATULA") || id.includes("WHISK") || id.includes("LADLE") || id.includes("DEGLAZE") || id.includes("SLICER") || id.includes("GRATER") || id.includes("PROCESSOR")) {
    return "cookware_food";
  }
  if (id.includes("HOOD") || id.includes("DISHWASHER") || id.includes("SINK") || id.includes("BOOSTER") || id.includes("PREPTABLE")) {
    return "ventilation";
  }
  if (id.includes("WALKIN") || id.includes("BLAST") || id.includes("REACHIN") || id.includes("UNDERCOUNTER") || id.includes("MAKELINE") || id.includes("ICEMACHINE") || id.includes("SPEEDRACK") || id.includes("CARRIER") || id.includes("HOLDING")) {
    return "cold_chain";
  }
  if (id.includes("COMBI") || id.includes("FRYER") || id.includes("KETTLE") || id.includes("MICROWAVE") || id.includes("TOASTER") || id.includes("AIRFRYER") || id.includes("SOUSVIDE")) {
    return "ovens_fryers";
  }
  return "hot_line";
}

export function prettyEquipmentName(archetypeId: string): string {
  const flagship = FLAGSHIP_EQUIPMENT.find((f) => f.id === archetypeId);
  if (flagship) return flagship.label;
  return archetypeId
    .replace(/^NSF-/, "")
    .replace(/^ARCH-FOOD-/, "")
    .split("-")
    .map((w) => (w.length <= 2 ? w : w[0] + w.slice(1).toLowerCase()))
    .join(" ");
}
