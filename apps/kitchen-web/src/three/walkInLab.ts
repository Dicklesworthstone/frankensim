/**
 * Order it, keep it, use it (chapter 04, after the cold chain lab): a week of a commissary's meat and fish through its
 * walk-in cooler and walk-in freezer, side by side, their roofs off and their fronts cut low. The week is the model's
 * (sim/walkInWeek.ts): what each truck brings for the lunches until the next, how long each case keeps, what each
 * morning's prep takes and what each night's check throws out. Every case it moves, a cook carries: in off the truck,
 * from the freezer to the walk-in to thaw, out to the kitchen for prep or to the bin, played in time-lapse.
 *
 * The cases stand on the wire shelving along each box's back wall, the hottest-cooked lowest: chicken at the bottom,
 * ground beef above, whole beef and fish at the top. Each carries a day-colour label: the day it must be used by, its
 * band green with three days or more to go, amber with two, orange on its last day; blue frozen, violet thawing.
 */
import * as THREE from "three";
import {
  DAY,
  DAY3,
  FREEZER_PLACES,
  ITEM,
  ITEMS,
  MENU,
  THIS_WEEK,
  TRUCK_DAYS,
  WALKIN_PLACES,
  layoutWeek,
  manyWeeks,
  manyWeeksRun,
  runWeeks,
  thisWeekCounts,
  type CaseAt,
  type End,
  type ItemKey,
  type Layout,
  type Lot,
  type ManyWeeks,
  type Policy,
  type Spot,
  type Staged,
  type Trucks,
  type WeekResult,
  type Where,
} from "../sim/walkInWeek";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { buildFixture, type Model, type Piece } from "./models/all";
import { BACK_SHELF } from "./models/cold";
import { Routine, makeStaff, panSeconds, type RoutineStep, type WorkOpts } from "./models/people";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

type XZ = readonly [number, number];

/* Where everything stands. The two 8 x 10 ft boxes side by side, their fronts at z = d/2 (the front walls' outer faces
 * 87 mm back of that), as one combination box: the freezer's left wall is the cooler's right. The freezer's insulated
 * floor is sunk into the kitchen's slab, flush with it, so the cook walks straight in. */
const COOLER = { id: "NSF-WALKIN-COOLER-8X10-BACK", w: 2.438, d: 3.048, h: 2.438 } as const;
const FREEZER = { id: "NSF-WALKIN-FREEZER-8X10-BACK", w: 2.438, d: 3.048, h: 2.438 } as const;
const PANEL = 0.1016;
const GAP = 0.004;
const X0: Record<Where, number> = { walkin: -GAP / 2 - COOLER.w, freezer: GAP / 2 };
const BOX_Z = COOLER.d / 2;
const FRONT_OUT = BOX_Z - 0.087;
const FRONT_IN = FRONT_OUT - PANEL;
/** The front walls are cut away above this. */
const CUT_Y = 0.3;
const FZ_SINK = PANEL - 0.0015;
/** Each box's floor, world. */
const FLOOR_Y: Record<Where, number> = { walkin: 0, freezer: PANEL - FZ_SINK };
/** Each door's opening: its middle (0.76 m from the box's left side, 0.914 m wide). */
const DOOR0 = 0.76;
const DOOR_W = 0.914;
const DOOR_X: Record<Where, number> = { walkin: X0.walkin + DOOR0 + DOOR_W / 2, freezer: X0.freezer + DOOR0 + DOOR_W / 2 };
/** How far from a door's middle the cook is clear of its leaf's swing (half the opening to the hinge, the leaf the
 *  opening's width, and a body's half width). */
const CLEAR_M = 1.7;
/** Each door's hinge pin, on the opening's left, 5 cm out from the front; and how far from it the leaf (1.01 m with
 *  its lap) reaches, with a foot beside it. The leaf swings out and round to the left, 110 degrees. */
const HINGE_X: Record<Where, number> = { walkin: X0.walkin + DOOR0 - 0.049, freezer: X0.freezer + DOOR0 - 0.049 };
const HINGE_Z_OUT = 0.05;
const REACH_M = 1.18;

/** A case: 17.3 x 14 in, 7.5 in deep for 40 lb, 5 in for 20 lb of salmon; its frame at the top's middle (where the
 *  hands' grip is taken from), the long side across the body as it is carried, its label on the side toward the body. */
const CL = 0.44;
const CD = 0.36;
const caseH = (item: ItemKey) => (item === "salmon" ? 0.13 : 0.19);
const HALF = CL / 2;

/** The back shelving's cases: each unit holds two, 6 cm apart, along its middle; the columns in the model's order,
 *  nearest the door first. Their middles 4 cm behind the shelves' fronts plus half a case. */
const COL_BOX_X: readonly number[] = (() => {
  const xs = BACK_SHELF.x0.flatMap((x) => [x + BACK_SHELF.L / 2 - 0.25, x + BACK_SHELF.L / 2 + 0.25]);
  const door = DOOR0 + DOOR_W / 2;
  return xs.sort((a, b) => Math.abs(a - door) - Math.abs(b - door));
})();
const BACK_FRONT_Z = BOX_Z + (-COOLER.d + PANEL) + BACK_SHELF.gap + BACK_SHELF.D;
const CASE_Z = BACK_FRONT_Z - 0.04 - CD / 2;
/** How far in front of a case's middle the cook stands to take it (the pick and place reach). */
const STANCE = 0.43;
/** The floor places: along the cooler's right side (where its own right wall's shelving would stand) and the freezer's
 *  right of the way in, each case end-on to the cook, who faces +x to take it. */
const FLOOR_BOX: Record<Where, { x: number; z: readonly number[] }> = {
  walkin: { x: 2.15, z: [-0.45, -0.95, -1.45] },
  freezer: { x: 1.68, z: [-0.5, -1.0] },
};
/** Outside: the lane along the fronts to the dock (off to the right), and the one further out to the kitchen (off to
 *  the left), a foot's width clear of an open door's tip (0.99 m out, propped open; the door only swings with the
 *  cook further off than its reach). Inside: just in, past the side shelving, and level with the back. */
const LANE_DOCK = FRONT_OUT + 0.47;
const LANE_KITCHEN = FRONT_OUT + 1.25;
const DOCK: XZ = [4.4, LANE_DOCK];
const KITCHEN: XZ = [-4.4, LANE_KITCHEN];
const AT_DOOR = FRONT_OUT + 0.35;
const IN1 = BOX_Z - 0.45;
const IN2 = BOX_Z - 1.62;
const IN3 = BOX_Z - 1.95;
const FREE_MPS = 1.15;
const CARRY_MPS = 0.95;
/** The stage's quickened time while a job is under way: at least this. A job plays for no longer than a few seconds
 *  and a little over one for each case carried, so a big one plays at much the pace of a small one. */
const LAPSE_MIN = 5;
const JOB_BASE_S = 6;
const JOB_PER_CASE_S = 1.1;
/** The clock's hours a second between jobs. */
const IDLE_HPS = 5;
/** How long a door takes to swing shut on its closer (s, the stage's time), and to be pulled open (s, the job's own
 *  time, as a door is), done this long before the cook is within its leaf's reach. */
const DOOR_S = 0.35;
const OPEN_S = 1.2;
const OPEN_LEAD_S = 0.3;

/** The label's colours: three days or more to go, two, the last day, past it; frozen; thawing. */
const GOOD = "#4ade80";
const SOON = "#fbbf24";
const LAST = "#fb923c";
const PAST = "#f43f5e";
const FROZEN = "#60a5fa";
const THAWING = "#c4b5fd";
/** The items' own colours, in the chart. */
const ITEM_C: Record<ItemKey, string> = { salmon: "#fda4af", ribs: "#f87171", beef: "#fb923c", thighs: "#fde68a" };
const SHORT_NAME: Record<ItemKey, string> = { salmon: "salmon", ribs: "short ribs", beef: "ground beef", thighs: "chicken thighs" };

let detailNow: 1 | 2 | null = null;
const roomDetail = (): 1 | 2 => (detailNow ??= Math.min(window.innerWidth, window.innerHeight) < 600 ? 1 : 2);
export function walkInPieces(): readonly Piece[] {
  const d = roomDetail();
  return [
    [COOLER.id, COOLER, d],
    [FREEZER.id, FREEZER, d],
  ];
}

const lb = (v: number) => `${Math.round(v).toLocaleString("en-US")} lb`;
const smooth = (x: number) => {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
};
const clockText = (hour: number) => {
  const h = Math.min(hour, 7 * 24 - 1 / 60);
  const day = Math.min(6, Math.floor(h / 24));
  const m = Math.floor((h - day * 24) * 60 + 1e-6);
  const hh = Math.floor(m / 60) % 24;
  const mm = m % 60;
  const ap = hh < 12 ? "am" : "pm";
  const h12 = hh % 12 === 0 ? 12 : hh % 12;
  return `${DAY[day]}, ${h12}:${String(mm).padStart(2, "0")} ${ap}`;
};
const cases = (n: number, what: string) => `${n} case${n === 1 ? "" : "s"} of ${what}`;
const lunchesN = (n: number) => `${n.toLocaleString("en-US")} lunch${n === 1 ? "" : "es"}`;
const listAnd = (xs: string[]) => (xs.length <= 1 ? (xs[0] ?? "") : `${xs.slice(0, -1).join(", ")} and ${xs[xs.length - 1]}`);

/* ------------------------------------------------------------------------------------------------------------ */
/* The cases' cardboard and labels, drawn once.                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */
function kraftCanvas(w: number, h: number, seed: number): [HTMLCanvasElement, CanvasRenderingContext2D] {
  const c = document.createElement("canvas");
  c.width = w;
  c.height = h;
  const g = c.getContext("2d")!;
  g.fillStyle = "#b88a5c";
  g.fillRect(0, 0, w, h);
  let s = seed;
  const r = () => ((s = (s * 16807) % 2147483647) / 2147483647);
  // The board's mottle and its fibres, and the flutes faintly through the liner.
  for (let i = 0; i < (w * h) / 60; i++) {
    const v = r();
    g.fillStyle = v < 0.5 ? `rgba(90,58,30,${0.05 + 0.05 * r()})` : `rgba(240,210,170,${0.04 + 0.05 * r()})`;
    g.fillRect(r() * w, r() * h, 1 + r() * 3, 1 + r() * 2);
  }
  g.strokeStyle = "rgba(80,50,25,0.05)";
  g.lineWidth = 1;
  for (let x = 0; x < w; x += 5) {
    g.beginPath();
    g.moveTo(x + 0.5, 0);
    g.lineTo(x + 0.5, h);
    g.stroke();
  }
  // Darker at the creased edges.
  const e = g.createLinearGradient(0, 0, 0, h);
  e.addColorStop(0, "rgba(70,45,22,0.22)");
  e.addColorStop(0.05, "rgba(70,45,22,0)");
  e.addColorStop(0.95, "rgba(70,45,22,0)");
  e.addColorStop(1, "rgba(70,45,22,0.22)");
  g.fillStyle = e;
  g.fillRect(0, 0, w, h);
  return [c, g];
}
const texOf = (c: HTMLCanvasElement) => {
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = 4;
  return t;
};
/** Each item's case, fresh or frozen: its six faces' materials (+x, -x, +y, -y, +z, -z), the print on the long sides. */
const caseMats = new Map<string, THREE.Material[]>();
function caseMaterials(item: ItemKey, frozen: boolean): THREE.Material[] {
  const key = `${item}|${frozen}`;
  let m = caseMats.get(key);
  if (m) return m;
  const it = ITEM[item];
  const H = caseH(item);
  const [sc, sg] = kraftCanvas(512, Math.round((512 * H) / CL), item.length * 31 + (frozen ? 7 : 3));
  const sh = sc.height;
  // Clear of the label, which takes the side's right end (see makeCase).
  const room = 512 * (1 - (labelW(item) + 0.045) / CL) - 26;
  const fit = (text: string, weight: number, px: number) => {
    sg.font = `${weight} ${px}px system-ui, sans-serif`;
    const w = sg.measureText(text).width;
    if (w > room) sg.font = `${weight} ${Math.floor((px * room) / w)}px system-ui, sans-serif`;
  };
  sg.fillStyle = "rgba(48,30,18,0.88)";
  sg.textBaseline = "middle";
  fit(it.name.toUpperCase(), 700, Math.round(sh * 0.2));
  sg.fillText(it.name.toUpperCase(), 26, sh * 0.3);
  fit(`${it.caseLb} LB NET WT · RAW`, 600, Math.round(sh * 0.12));
  sg.fillText(`${it.caseLb} LB NET WT · RAW`, 26, sh * 0.55);
  sg.fillStyle = frozen ? "rgba(30,64,140,0.9)" : "rgba(48,30,18,0.88)";
  const keep = frozen ? "KEEP FROZEN 0 °F OR BELOW" : "KEEP REFRIGERATED 28–40 °F";
  fit(keep, 600, Math.round(sh * 0.12));
  sg.fillText(keep, 26, sh * 0.76);
  const side = new THREE.MeshStandardMaterial({ map: texOf(sc), roughness: 0.92, metalness: 0 });
  const [ec, eg] = kraftCanvas(256, Math.round((256 * H) / CD), 91 + item.length);
  eg.fillStyle = "rgba(48,30,18,0.7)";
  eg.font = `600 ${Math.round(ec.height * 0.14)}px system-ui, sans-serif`;
  eg.textAlign = "center";
  eg.textBaseline = "middle";
  eg.fillText("↑ THIS SIDE UP ↑", 128, ec.height * 0.35);
  const end = new THREE.MeshStandardMaterial({ map: texOf(ec), roughness: 0.92, metalness: 0 });
  const [tc, tg] = kraftCanvas(256, Math.round((256 * CD) / CL), 57 + item.length);
  // The flaps' seam down the middle, taped.
  tg.fillStyle = "rgba(60,38,18,0.55)";
  tg.fillRect(0, tc.height / 2 - 1, 256, 2);
  tg.fillStyle = "rgba(205,160,100,0.85)";
  tg.fillRect(0, tc.height / 2 - 11, 256, 22);
  tg.fillStyle = "rgba(255,240,210,0.18)";
  tg.fillRect(0, tc.height / 2 - 11, 256, 3);
  const top = new THREE.MeshStandardMaterial({ map: texOf(tc), roughness: 0.75, metalness: 0 });
  const [bc] = kraftCanvas(128, 104, 13);
  const bottom = new THREE.MeshStandardMaterial({ map: texOf(bc), roughness: 0.95, metalness: 0 });
  m = [end, end, top, bottom, side, side];
  caseMats.set(key, m);
  return m;
}
/** A label: the use-by day on its coloured band, the item and the day it came in below. */
const labelTex = new Map<string, THREE.CanvasTexture>();
function labelTexture(item: ItemKey, band: string, head: string, foot: string): THREE.CanvasTexture {
  const key = `${item}|${band}|${head}|${foot}`;
  let t = labelTex.get(key);
  if (t) return t;
  const c = document.createElement("canvas");
  c.width = 256;
  c.height = 176;
  const g = c.getContext("2d")!;
  g.fillStyle = "#f6f4ee";
  g.fillRect(0, 0, 256, 176);
  g.fillStyle = band;
  g.fillRect(0, 0, 256, 70);
  g.fillStyle = "#111";
  g.textAlign = "center";
  g.textBaseline = "middle";
  g.font = "800 40px system-ui, sans-serif";
  g.fillText(head, 128, 37);
  g.font = "700 26px system-ui, sans-serif";
  g.fillText(SHORT_NAME[item].toUpperCase(), 128, 104);
  g.font = "600 22px system-ui, sans-serif";
  g.fillStyle = "#444";
  g.fillText(foot, 128, 143);
  g.strokeStyle = "rgba(0,0,0,0.25)";
  g.lineWidth = 3;
  g.strokeRect(1.5, 1.5, 253, 173);
  t = texOf(c);
  labelTex.set(key, t);
  return t;
}

interface CaseMesh {
  group: THREE.Group;
  box: THREE.Mesh;
  label: THREE.Mesh;
  item: ItemKey;
  frozen: boolean;
  key: string;
}
const boxGeos = new Map<ItemKey, THREE.BoxGeometry>();
/** The label: as big as the case's side takes, 256 x 176 like its picture. */
const labelGeos = new Map<ItemKey, THREE.PlaneGeometry>();
const labelW = (item: ItemKey) => Math.min(0.2, (caseH(item) - 0.03) * (256 / 176));
function makeCase(item: ItemKey, frozen: boolean): CaseMesh {
  const H = caseH(item);
  let g = boxGeos.get(item);
  if (!g) boxGeos.set(item, (g = new THREE.BoxGeometry(CL, H, CD)));
  const group = new THREE.Group();
  const box = new THREE.Mesh(g, caseMaterials(item, frozen));
  box.position.y = -H / 2;
  box.castShadow = true;
  box.receiveShadow = true;
  group.add(box);
  let lg = labelGeos.get(item);
  if (!lg) labelGeos.set(item, (lg = new THREE.PlaneGeometry(labelW(item), (labelW(item) * 176) / 256)));
  const label = new THREE.Mesh(lg, new THREE.MeshStandardMaterial({ roughness: 0.55, metalness: 0, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -2 }));
  // On the side toward the carrier's body (out of the shelf once it stands there), at the end the print leaves clear.
  label.position.set(-(CL / 2 - 0.025 - labelW(item) / 2), -H * 0.5, -CD / 2 - 0.0007);
  label.rotation.y = Math.PI;
  label.receiveShadow = true;
  group.add(label);
  return { group, box, label, item, frozen, key: "" };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Places, and the ways between them.                                                                             */
/* ------------------------------------------------------------------------------------------------------------ */
/** A case place's top middle and the way the case faces (its label out of the shelf, or toward the aisle). */
function spotPose(s: Spot, item: ItemKey, out: THREE.Vector3): number {
  const H = caseH(item);
  if (!s.floor) {
    out.set(X0[s.box] + COL_BOX_X[s.col]!, FLOOR_Y[s.box] + BACK_SHELF.tiers[s.tier]! + (s.level + 1) * H, CASE_Z);
    return Math.PI;
  }
  const f = FLOOR_BOX[s.box];
  const row = Math.floor(s.col / f.z.length);
  out.set(X0[s.box] + f.x - row * 0.42, FLOOR_Y[s.box] + (s.level + 1) * H, BOX_Z + f.z[s.col % f.z.length]!);
  return Math.PI / 2;
}
/** Where the cook stands to take a case from a place, or put one there, and which way they face. */
function stanceOf(s: Spot, item: ItemKey): { at: XZ; yaw: number } {
  const p = new THREE.Vector3();
  const yaw = spotPose(s, item, p);
  return s.floor ? { at: [p.x - STANCE, p.z], yaw } : { at: [p.x, p.z + STANCE], yaw };
}
/** From just outside a box's door to a stance in it. */
function inPath(box: Where, s: Spot, stance: XZ): XZ[] {
  const dx = DOOR_X[box];
  const pts: XZ[] = [[dx, AT_DOOR], [dx, IN1]];
  if (!s.floor) pts.push([dx, IN2], [stance[0], IN3]);
  pts.push(stance);
  return pts;
}
type Loc = { kind: "dock" } | { kind: "kitchen" } | { kind: "spot"; spot: Spot; stance: XZ; yaw: number };
const near = (a: XZ, b: XZ) => Math.hypot(a[0] - b[0], a[1] - b[1]) < 0.02;
const ptOf = (l: Loc): XZ => (l.kind === "dock" ? DOCK : l.kind === "kitchen" ? KITCHEN : l.stance);
/** The way from one place to another: through the doors, along the dock's lane if either end is the dock, else the
 *  kitchen's. Within a box from one stance to another: back level with the back's front, across, in. */
function route(a: Loc, b: Loc): XZ[] {
  let pts: XZ[];
  if (a.kind === "spot" && b.kind === "spot" && a.spot.box === b.spot.box) {
    pts = a.spot.floor || b.spot.floor ? [a.stance, [DOOR_X[a.spot.box], IN1], b.stance] : [a.stance, [a.stance[0], IN3], [b.stance[0], IN3], b.stance];
  } else {
    const lane = a.kind === "dock" || b.kind === "dock" ? LANE_DOCK : LANE_KITCHEN;
    const out: XZ[] = a.kind === "spot" ? [...inPath(a.spot.box, a.spot, a.stance)].reverse() : [ptOf(a)];
    const inn: XZ[] = b.kind === "spot" ? inPath(b.spot.box, b.spot, b.stance) : [ptOf(b)];
    const ax = a.kind === "spot" ? DOOR_X[a.spot.box] : null;
    const bx = b.kind === "spot" ? DOOR_X[b.spot.box] : null;
    pts = [...out];
    if (ax !== null) pts.push([ax, lane]);
    if (bx !== null) pts.push([bx, lane]);
    pts.push(...inn);
  }
  const outp: XZ[] = [];
  for (const p of pts) if (!outp.length || !near(outp[outp.length - 1]!, p)) outp.push(p);
  return outp;
}
const headingOf = (from: XZ, to: XZ) => Math.atan2(to[0] - from[0], to[1] - from[1]);

/* ------------------------------------------------------------------------------------------------------------ */
/* The film: the week's jobs, each a routine for the cook, and the clock between them.                            */
/* ------------------------------------------------------------------------------------------------------------ */
interface Carry {
  caseId: number;
  item: ItemKey;
  lot: number;
  from: End;
  to: End;
  /** Routine seconds the case is in the cook's hands (from the pick, or off the truck, to the place, or out of view). */
  t0: number;
  t1: number;
}
interface Job {
  staged: Staged;
  /** The cases on hand before it. */
  before: CaseAt[];
  routine: Routine | null;
  carries: Carry[];
  /** Routine seconds each box's door is open (propped for the job), from it starting to open to it starting to shut. */
  doors: Partial<Record<Where, [number, number]>>;
  /** Routine seconds the job lasts: its routine, and on (the cook waiting at the dock or the kitchen) until the doors
   *  are shut. */
  span: number;
  /** Film seconds it starts and ends, the lapse it plays at. */
  t0: number;
  t1: number;
  lapse: number;
  hour: number;
}
interface Film {
  week: WeekResult;
  layout: Layout;
  jobs: Job[];
  length: number;
  /** Every lot seen in the week, by id. */
  lots: Map<number, Lot>;
}

export function mountWalkInLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-wo-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-wo-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  renderer.localClippingEnabled = true;
  const phone = window.matchMedia("(max-width: 720px)").matches;
  const detail = roomDetail();

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 9, 18);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.55;
  const lights = studioLights(scene, new THREE.Vector3(0, 1.0, 0), 4.4, phone ? 1024 : 2048);
  lights.key.position.set(-2.2, 7.2, 4.6);
  lights.key.intensity = 1.45;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.6, 4.0, -3.4);
  lights.rim.intensity = 0.6;
  lights.hemi.intensity = 0.5;
  {
    const floor = tilePlane(16, 9);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0, 0, 2.6);
    scene.add(floor);
  }

  // The front walls cut above CUT_Y: a fragment goes if it is both above the cut and in the front wall's thickness.
  const frontCut = [new THREE.Plane(new THREE.Vector3(0, -1, 0), CUT_Y), new THREE.Plane(new THREE.Vector3(0, 0, -1), FRONT_IN - 0.001)];
  const clipCopies = new Map<THREE.Material, THREE.Material>();
  function clip(o: THREE.Object3D): void {
    const copy = (m: THREE.Material) => {
      let c = clipCopies.get(m);
      if (!c) {
        c = m.clone();
        c.clippingPlanes = frontCut;
        c.clipIntersection = true;
        c.clipShadows = true;
        clipCopies.set(m, c);
      }
      return c;
    };
    o.traverse((x) => {
      const mesh = x as THREE.Mesh;
      if (!mesh.isMesh) return;
      mesh.material = Array.isArray(mesh.material) ? mesh.material.map(copy) : copy(mesh.material);
    });
  }
  const foam = new THREE.MeshStandardMaterial({ color: 0xd8b75e, roughness: 0.93 });
  const doorCaps: { door: THREE.Object3D; caps: THREE.Mesh[] }[] = [];
  const capAt = new THREE.Vector3();
  /** Keep the door's caps at the cut as it rises on its cam. */
  const placeDoorCaps = () => {
    for (const { door, caps } of doorCaps) {
      door.getWorldPosition(capAt);
      for (const c of caps) c.position.y = CUT_Y - 0.0016 - capAt.y;
    }
  };
  function box(spec: typeof COOLER | typeof FREEZER, where: Where, frost = 0): Model {
    const m = buildFixture(spec.id, spec, detail);
    m.group.position.set(X0[where], where === "freezer" ? -FZ_SINK : 0, BOX_Z);
    m.hooks.set?.("roof", 1);
    if (frost) m.hooks.set?.("frost", frost);
    // The strip curtain hangs from the door's head, which is cut away.
    for (const [name, j] of m.joints) if (name.startsWith("pvc_strip_")) j.visible = false;
    clip(m.group);
    scene.add(m.group);
    // The cut: the panels' foamed-in-place core between their skins, either side of the door (which has its own).
    for (const [a, b] of [[0.001, DOOR0], [DOOR0 + DOOR_W, spec.w - 0.001]] as const) {
      const cap = new THREE.Mesh(new THREE.BoxGeometry(b - a, 0.002, PANEL - 0.0018), foam);
      cap.position.set(X0[where] + (a + b) / 2, CUT_Y - 0.0012, FRONT_OUT - PANEL / 2);
      cap.receiveShadow = true;
      scene.add(cap);
    }
    // The door cut at the same height, swinging with it: a cap over each part of it that stands full height (the
    // slab and the plug), from its outline in the door's frame.
    const door = m.joints.get("cam_lift_door");
    if (door) {
      m.group.updateMatrixWorld(true);
      const inv = door.matrixWorld.clone().invert();
      const bb = new THREE.Box3();
      const caps: THREE.Mesh[] = [];
      door.traverse((x) => {
        const mesh = x as THREE.Mesh;
        if (!mesh.isMesh) return;
        mesh.geometry.computeBoundingBox();
        bb.copy(mesh.geometry.boundingBox!).applyMatrix4(mesh.matrixWorld).applyMatrix4(inv);
        if (bb.max.y - bb.min.y < 1.2 || bb.min.y > CUT_Y) return;
        const c = new THREE.Mesh(new THREE.BoxGeometry(Math.max(0.001, bb.max.x - bb.min.x - 0.002), 0.002, Math.max(0.001, bb.max.z - bb.min.z - 0.002)), foam);
        c.position.set((bb.min.x + bb.max.x) / 2, 0, (bb.min.z + bb.max.z) / 2);
        caps.push(c);
      });
      for (const c of caps) door.add(c);
      doorCaps.push({ door, caps });
    }
    const lamp = new THREE.PointLight(0xeef4ff, 5.5, 5.5, 1.6);
    lamp.position.set(X0[where] + spec.w / 2, spec.h - 0.15, -0.2);
    scene.add(lamp);
    return m;
  }
  const boxes: Record<Where, Model> = { walkin: box(COOLER, "walkin"), freezer: box(FREEZER, "freezer", 0.45) };
  {
    const j = boxes.walkin.joints.get("part_wall_right");
    if (j) j.visible = false;
  }
  boxes.walkin.hooks.set?.("temp", 3);
  boxes.freezer.hooks.set?.("temp", -18);

  /* The floor's warning: a red ring round any case on it. */
  const ringGeo = (() => {
    const s = new THREE.Shape();
    const w = CD / 2 + 0.05;
    const l = CL / 2 + 0.05;
    s.moveTo(-w, -l);
    s.lineTo(w, -l);
    s.lineTo(w, l);
    s.lineTo(-w, l);
    s.lineTo(-w, -l);
    const hole = new THREE.Path();
    const t = 0.025;
    hole.moveTo(-w + t, -l + t);
    hole.lineTo(-w + t, l - t);
    hole.lineTo(w - t, l - t);
    hole.lineTo(w - t, -l + t);
    hole.lineTo(-w + t, -l + t);
    s.holes.push(hole);
    const g = new THREE.ShapeGeometry(s);
    g.rotateX(-Math.PI / 2);
    return g;
  })();
  const ringMat = new THREE.MeshBasicMaterial({ color: PAST, transparent: true, opacity: 0.85, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -4, polygonOffsetUnits: -4, toneMapped: false });
  const rings: THREE.Mesh[] = [];
  const ringOf = (i: number) => {
    while (rings.length <= i) {
      const r = new THREE.Mesh(ringGeo, ringMat);
      r.renderOrder = 2;
      r.visible = false;
      scene.add(r);
      rings.push(r);
    }
    return rings[i]!;
  };

  /* The cook. */
  const cook = makeStaff("p50", phone ? 1 : 2, { seed: 37, apron: 0x3b4a3f, sleeves: "cuffed" });
  scene.add(cook.group);
  cook.group.visible = false;

  /* The cases. */
  const caseMeshes = new Map<number, CaseMesh>();
  const caseOf = (id: number, item: ItemKey, frozen: boolean): CaseMesh => {
    let c = caseMeshes.get(id);
    if (!c || c.item !== item || c.frozen !== frozen) {
      if (c) {
        c.group.removeFromParent();
        (c.label.material as THREE.Material).dispose();
      }
      c = makeCase(item, frozen);
      caseMeshes.set(id, c);
    }
    return c;
  };

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The week, as a film.                                                                                        */
  /* ---------------------------------------------------------------------------------------------------------- */
  function locOf(e: End, item: ItemKey): Loc {
    if (e === "dock") return { kind: "dock" };
    if (e === "kitchen") return { kind: "kitchen" };
    const s = stanceOf(e, item);
    return { kind: "spot", spot: e, stance: s.at, yaw: s.yaw };
  }
  const rimOf = (s: Spot, item: ItemKey) => {
    const p = new THREE.Vector3();
    spotPose(s, item, p);
    return p;
  };
  function makeJob(staged: Staged, before: CaseAt[]): Omit<Job, "t0" | "t1" | "lapse" | "span"> {
    const ev = staged.ev;
    const job = { staged, before, routine: null as Routine | null, carries: [] as Carry[], doors: {} as Job["doors"], hour: ev.hour };
    if (!staged.moves.length) return job;
    const steps: RoutineStep[] = [];
    const startL: Loc = ev.kind === "truck" ? { kind: "dock" } : { kind: "kitchen" };
    let here = startL as Loc;
    let yaw = headingOf(ptOf(startL), [0, ptOf(startL)[1]]);
    const startYaw = yaw;
    const marks: { i0: number; i1: number; m: Staged["moves"][number] }[] = [];
    const boxSteps: Partial<Record<Where, [number, number]>> = {};
    const touch = (b: Where, i0: number, i1: number) => {
      const was = boxSteps[b];
      boxSteps[b] = was ? [Math.min(was[0], i0), Math.max(was[1], i1)] : [i0, i1];
    };
    const carryOpts: WorkOpts = { halfWidth: HALF };
    const walks: { i: number; path: XZ[] }[] = [];
    const walk = (to: Loc, carrying: boolean): XZ[] => {
      const path = route(here, to);
      if (path.length < 2) return [];
      // Face the way off before setting out (a turn about, at the dock or the kitchen).
      const h0 = headingOf(path[0]!, path[1]!);
      if (Math.abs(Math.atan2(Math.sin(h0 - yaw), Math.cos(h0 - yaw))) > 0.35) steps.push({ kind: "turn", yaw: h0, upper: carrying ? "carry" : "free", opts: carrying ? carryOpts : undefined });
      steps.push({ kind: "walk", path: path.slice(1), speed: carrying ? CARRY_MPS : FREE_MPS, upper: carrying ? "carry" : "free", opts: carrying ? carryOpts : undefined });
      walks.push({ i: steps.length - 1, path });
      yaw = headingOf(path[path.length - 2]!, path[path.length - 1]!);
      here = to;
      return path;
    };
    const face = (y: number, carrying: boolean) => {
      if (Math.abs(Math.atan2(Math.sin(y - yaw), Math.cos(y - yaw))) > 0.02) steps.push({ kind: "turn", yaw: y, upper: carrying ? "carry" : "free", opts: carrying ? carryOpts : undefined });
      yaw = y;
    };
    for (const m of staged.moves) {
      const to = locOf(m.to, m.item);
      if (m.from === "dock") {
        if (here.kind !== "dock") walk({ kind: "dock" }, false);
        const i0 = steps.length;
        walk(to, true);
        if (to.kind !== "spot") continue;
        face(to.yaw, true);
        const at = rimOf(to.spot, m.item);
        const o: WorkOpts = { at, surfaceY: at.y - caseH(m.item), halfWidth: HALF };
        steps.push({ kind: "work", work: "place", duration: panSeconds("place", o), opts: o });
        marks.push({ i0, i1: steps.length, m });
        touch(to.spot.box, i0, steps.length);
        continue;
      }
      const from = locOf(m.from, m.item);
      if (from.kind !== "spot") continue;
      const wasIn = here.kind === "spot" ? here.spot.box : null;
      const iw = steps.length;
      if (here.kind !== "spot" || !near(here.stance, from.stance)) walk(from, false);
      face(from.yaw, false);
      if (wasIn !== from.spot.box) touch(from.spot.box, iw, steps.length);
      const at0 = rimOf(from.spot, m.item);
      const o0: WorkOpts = { at: at0, surfaceY: at0.y - caseH(m.item), halfWidth: HALF };
      const i0 = steps.length;
      steps.push({ kind: "work", work: "pick", duration: panSeconds("pick", o0), opts: o0 });
      if (to.kind === "kitchen") {
        walk(to, true);
        marks.push({ i0, i1: steps.length, m });
        touch(from.spot.box, i0, steps.length);
        continue;
      }
      if (to.kind !== "spot") continue;
      if (!near(here.kind === "spot" ? here.stance : ptOf(here), to.stance)) walk(to, true);
      face(to.yaw, true);
      const at1 = rimOf(to.spot, m.item);
      const o1: WorkOpts = { at: at1, surfaceY: at1.y - caseH(m.item), halfWidth: HALF };
      steps.push({ kind: "work", work: "place", duration: panSeconds("place", o1), opts: o1 });
      marks.push({ i0, i1: steps.length, m });
      touch(from.spot.box, i0, steps.length);
      touch(to.spot.box, i0, steps.length);
    }
    // And away: back to the dock after a truck, else to the kitchen.
    walk(startL, false);
    const s0 = ptOf(startL);
    const routine = new Routine(cook, { x: s0[0], z: s0[1], yaw: startYaw }, steps, false);
    const T = routine.stepTimes;
    const at = (i: number) => T[i] ?? routine.duration;
    job.routine = routine;
    job.carries = marks.map(({ i0, i1, m }) => ({ caseId: m.caseId, item: m.item, lot: m.lot, from: m.from, to: m.to, t0: at(i0), t1: at(i1) }));
    // Each door the job uses is pulled open as the cook comes to it, done before they are within its leaf's reach,
    // stays propped while they come and go, and starts to shut on its closer once they are out and clear of its swing
    // after their last time in that box. A walk's times go by its length, at its pace.
    const along = (path: XZ[]) => {
      const pts: { p: XZ; d: number }[] = [{ p: path[0]!, d: 0 }];
      let run = 0;
      for (let i = 1; i < path.length; i++) {
        const [ax, az] = path[i - 1]!;
        const [bx, bz] = path[i]!;
        const len = Math.hypot(bx - ax, bz - az);
        const n = Math.max(1, Math.ceil(len / 0.05));
        for (let k = 1; k <= n; k++) pts.push({ p: [ax + ((bx - ax) * k) / n, az + ((bz - az) * k) / n], d: run + (len * k) / n });
        run += len;
      }
      return pts;
    };
    const timeOn = (wk: { i: number; path: XZ[] }, d: number, len: number) => at(wk.i) + ((at(wk.i + 1) - at(wk.i)) * d) / Math.max(1e-6, len);
    const other = (b: Where): Where => (b === "walkin" ? "freezer" : "walkin");
    const inBox = (b: Where, [x, z]: XZ) => z < FRONT_OUT && Math.abs(x - DOOR_X[b]) < Math.abs(x - DOOR_X[other(b)]);
    const inReach = (b: Where, [x, z]: XZ) => z > FRONT_OUT - 0.05 && Math.hypot(x - HINGE_X[b], z - FRONT_OUT - HINGE_Z_OUT) < REACH_M;
    const clear = (b: Where, [x, z]: XZ) => z > FRONT_OUT + 0.1 && Math.hypot(x - DOOR_X[b], z - FRONT_OUT) >= CLEAR_M;
    for (const b of ["walkin", "freezer"] as const) {
      const s = boxSteps[b];
      if (!s) continue;
      let open = at(s[0]);
      reach: for (const wk of walks) {
        const pts = along(wk.path);
        for (const q of pts)
          if (inReach(b, q.p)) {
            open = timeOn(wk, q.d, pts[pts.length - 1]!.d) - OPEN_LEAD_S;
            break reach;
          }
      }
      let shut = at(s[1]);
      for (let k = walks.length - 1; k >= 0; k--) {
        const wk = walks[k]!;
        const pts = along(wk.path);
        let last = -1;
        for (let n = pts.length - 1; n >= 0 && last < 0; n--) if (inBox(b, pts[n]!.p)) last = n;
        if (last < 0) continue;
        const out = pts.findIndex((q, n) => n > last && clear(b, q.p));
        const len = pts[pts.length - 1]!.d;
        shut = timeOn(wk, out >= 0 ? pts[out]!.d : len, len);
        break;
      }
      job.doors[b] = [Math.max(0, open - OPEN_S), Math.max(open, shut)];
    }
    return job;
  }
  function makeFilm(p: Policy): Film {
    const week = runWeeks(p, thisWeekCounts());
    const layout = layoutWeek(week);
    const lots = new Map<number, Lot>();
    for (const l of week.start) lots.set(l.id, l);
    for (const e of week.events) for (const l of e.stock) if (!lots.has(l.id)) lots.set(l.id, l);
    const jobs: Job[] = [];
    let before = layout.start;
    let t = 0;
    let hour = 0;
    for (const st of layout.events) {
      t += Math.max(0, st.ev.hour - hour) / IDLE_HPS;
      const j = makeJob(st, before);
      const dur = j.routine?.duration ?? 0;
      const most = JOB_BASE_S + JOB_PER_CASE_S * st.moves.length;
      const lapse = dur > 0 ? Math.max(LAPSE_MIN, dur / most) : 1;
      // On until the last door is shut, its swing DOOR_S of film.
      let span = dur;
      for (const w of Object.values(j.doors)) if (w) span = Math.max(span, w[1] + DOOR_S * lapse);
      const len = dur > 0 ? span / lapse : 1.4;
      jobs.push({ ...j, t0: t, t1: t + len, lapse, span });
      t += len;
      hour = st.ev.hour + span / 3600;
      before = st.after;
    }
    t += Math.max(0, 7 * 24 - hour) / IDLE_HPS;
    return { week, layout, jobs, length: t, lots };
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* What the stage shows at a moment of the film.                                                               */
  /* ---------------------------------------------------------------------------------------------------------- */
  interface Moment {
    hour: number;
    job: Job | null;
    /** Routine seconds into the job. */
    r: number;
    cases: { id: number; item: ItemKey; lot: number; at: Spot | null; thrown: boolean }[];
    inHand: number | null;
  }
  function momentAt(f: Film, ft: number): Moment {
    let job: Job | null = null;
    let lastDone: Job | null = null;
    for (const j of f.jobs) {
      if (ft < j.t0) break;
      if (ft < j.t1) {
        job = j;
        break;
      }
      lastDone = j;
    }
    if (!job) {
      // Between jobs: the clock runs on from the last one's end.
      const base = lastDone ? lastDone.hour + lastDone.span / 3600 : 0;
      const from = lastDone ? lastDone.t1 : 0;
      const nextJob = f.jobs.find((j) => j.t0 > ft);
      const hour = Math.min(nextJob ? nextJob.hour : 7 * 24, base + (ft - from) * IDLE_HPS);
      const now = lastDone ? lastDone.staged.after : f.layout.start;
      return { hour, job: null, r: 0, cases: now.map((c) => ({ id: c.caseId, item: c.item, lot: c.lot, at: c.at, thrown: false })), inHand: null };
    }
    const r = (ft - job.t0) * job.lapse;
    const hour = job.hour + r / 3600;
    const map = new Map<number, { id: number; item: ItemKey; lot: number; at: Spot | null; thrown: boolean }>();
    for (const c of job.before) map.set(c.caseId, { id: c.caseId, item: c.item, lot: c.lot, at: c.at, thrown: false });
    let inHand: number | null = null;
    const thrown = job.staged.ev.kind === "check";
    for (const k of job.carries) {
      if (r >= k.t1) {
        if (k.to === "kitchen") map.delete(k.caseId);
        else if (k.to !== "dock") map.set(k.caseId, { id: k.caseId, item: k.item, lot: k.lot, at: k.to, thrown: false });
      } else if (r >= k.t0) {
        inHand = k.caseId;
        map.set(k.caseId, { id: k.caseId, item: k.item, lot: k.lot, at: null, thrown });
      } else if (k.from === "dock") map.delete(k.caseId);
    }
    if (!job.carries.length) for (const c of job.staged.after) map.set(c.caseId, { id: c.caseId, item: c.item, lot: c.lot, at: c.at, thrown: false });
    return { hour, job, r, cases: [...map.values()], inHand };
  }

  /** The label for a lot on an (absolute) day. */
  function labelFor(l: Lot, absDay: number, thrown: boolean): { band: string; head: string; foot: string; state: string } {
    const came = `IN ${DAY3[((l.inDay % 7) + 7) % 7]!.toUpperCase()}`;
    if (l.where === "freezer") return { band: FROZEN, head: "FROZEN", foot: came, state: "frozen" };
    const by = `USE BY ${DAY3[((l.lastDay % 7) + 7) % 7]!.toUpperCase()}`;
    if (thrown) return { band: PAST, head: by, foot: "LAST DAY OVER", state: "past" };
    if (l.readyDay > absDay) return { band: THAWING, head: by, foot: `THAWING · ${came}`, state: "thawing" };
    const n = l.lastDay - absDay + 1;
    return { band: n >= 3 ? GOOD : n === 2 ? SOON : n === 1 ? LAST : PAST, head: by, foot: came, state: n >= 3 ? "good" : n === 2 ? "soon" : n === 1 ? "last" : "past" };
  }

  const tmpP = new THREE.Vector3();
  let shown: Moment | null = null;
  function showMoment(f: Film, m: Moment): void {
    shown = m;
    const absDay = 7 + Math.min(6, Math.floor(m.hour / 24));
    const live = new Set<number>();
    let floorN = 0;
    for (const c of m.cases) {
      const lot = f.lots.get(c.lot);
      if (!lot) continue;
      // A frozen case is still a frozen case's box once it is in the walk-in, thawing.
      const frozenBox = f.week.policy.frozen;
      const cm = caseOf(c.id, c.item, frozenBox);
      live.add(c.id);
      const lab = labelFor(lot, absDay, c.thrown);
      const key = `${lab.band}|${lab.head}|${lab.foot}`;
      if (cm.key !== key) {
        cm.key = key;
        const mat = cm.label.material as THREE.MeshStandardMaterial;
        mat.map = labelTexture(c.item, lab.band, lab.head, lab.foot);
        mat.needsUpdate = true;
      }
      if (c.at) {
        if (cm.group.parent !== scene) scene.add(cm.group);
        const yaw = spotPose(c.at, c.item, tmpP);
        cm.group.position.copy(tmpP);
        cm.group.rotation.set(0, yaw, 0);
        if (c.at.floor && c.at.level === 0) {
          const ring = ringOf(floorN++);
          ring.visible = true;
          ring.position.set(tmpP.x, FLOOR_Y[c.at.box] + 0.002, tmpP.z);
        }
      } else if (m.inHand === c.id) {
        if (cm.group.parent !== cook.item) cook.item.add(cm.group);
        cm.group.position.set(0, 0, 0);
        cm.group.rotation.set(0, 0, 0);
      }
    }
    for (let i = floorN; i < rings.length; i++) rings[i]!.visible = false;
    for (const [id, cm] of caseMeshes) if (!live.has(id) && cm.group.parent) cm.group.removeFromParent();
    // The cook, at work or away; the doors pulled open as the cook comes, propped, and shut on their closers.
    if (m.job?.routine) {
      cook.group.visible = true;
      m.job.routine.apply(Math.min(m.r, m.job.routine.duration));
    } else cook.group.visible = false;
    for (const b of ["walkin", "freezer"] as const) {
      const w = m.job?.doors[b];
      const swing = DOOR_S * (m.job?.lapse ?? 1);
      const k = w ? smooth((m.r - w[0]) / OPEN_S) * (1 - smooth((m.r - w[1]) / swing)) : 0;
      if (doorK[b] !== k) {
        doorK[b] = k;
        boxes[b].hooks.set?.("door", k);
        boxes[b].group.updateMatrixWorld(true);
        placeDoorCaps();
      }
    }
  }
  const doorK: Record<Where, number> = { walkin: -1, freezer: -1 };

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The camera: from in front, raised, both boxes' back shelving in view over their cut fronts.                 */
  /* ---------------------------------------------------------------------------------------------------------- */
  const camera = new THREE.PerspectiveCamera(30, 16 / 9, 0.05, 60);
  const TARGET = new THREE.Vector3(0, 0.75, -0.2);
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    VIEW.set(0, wide ? 0.3 : 0.42, 1).normalize();
    keep.length = 0;
    // Both back shelvings top to bottom, the doors' cut fronts, and the cook's feet on the near lane outside them; on a
    // tall stage, on the kitchen's lane too (on a wide one that would cost the boxes a quarter of their size).
    const l = X0.walkin + BACK_SHELF.x0[0]!;
    const r = X0.freezer + BACK_SHELF.x0[1]! + BACK_SHELF.L;
    for (const x of [l, r]) keep.push(new THREE.Vector3(x, BACK_SHELF.H + 0.28, CASE_Z), new THREE.Vector3(x, BACK_SHELF.tiers[0]!, BACK_FRONT_Z));
    for (const b of ["walkin", "freezer"] as const) {
      keep.push(new THREE.Vector3(X0[b] + DOOR0 - 0.15, CUT_Y, FRONT_OUT), new THREE.Vector3(X0[b] + DOOR0 + DOOR_W + 0.15, CUT_Y, FRONT_OUT));
      keep.push(new THREE.Vector3(DOOR_X[b], 0, LANE_DOCK + 0.2));
      if (!wide) keep.push(new THREE.Vector3(DOOR_X[b], 0, LANE_KITCHEN + 0.15));
    }
  };
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 1.2, maxDistance: 16, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const hudEl = stageEl.querySelector<HTMLElement>(".wo-hud");
  const legendEl = stageEl.querySelector<HTMLElement>(".wo-legend");
  const hudRoom = (): number => {
    if (!hudEl) return 0;
    const parts = ["[data-wo-hud-clock]", "[data-wo-hud-doing]", "[data-wo-hud-note]"].map((k) => hudEl.querySelector<HTMLElement>(k));
    const was = parts.map((e) => e?.textContent ?? "");
    if (parts[0]) parts[0].textContent = "Wednesday, 10:00 pm";
    // Longer than any line the week says in any setting (the longest, 153 characters, a frozen truck's four meats).
    if (parts[1]) parts[1].textContent = "Next, Wednesday at 10:00 pm, the truck: 10 cases of salmon, 11 cases of short ribs, 10 cases of ground beef and 11 cases of chicken thighs, into the freezer";
    if (parts[2]) parts[2].textContent = "1,234 lunches short: bought on the day";
    const most = hudEl.offsetTop + hudEl.offsetHeight;
    parts.forEach((e, i) => {
      if (e) e.textContent = was[i]!;
    });
    return most;
  };
  const ndc = new THREE.Vector3();
  let fitW = 0;
  let fitH = 0;
  const frame = () => {
    if (userMoved) return;
    aim();
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    fitW = w;
    fitH = h;
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    const top = Math.min(0.9, 1 - (2 * (hudRoom() + 8)) / h);
    const bot = -1 + (2 * ((legendEl?.offsetHeight ?? 0) + 16)) / h;
    let dist = 2;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 30; dist += 0.05) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      lo = left = Infinity;
      hi = right = -Infinity;
      for (const k of keep) {
        ndc.copy(k).project(camera);
        lo = Math.min(lo, ndc.y);
        hi = Math.max(hi, ndc.y);
        left = Math.min(left, ndc.x);
        right = Math.max(right, ndc.x);
      }
      if (right - left < 1.92 && hi - lo < top - bot) break;
    }
    const shiftY = (top + bot) / 2 - (hi + lo) / 2;
    const shiftX = (left + right) / 2;
    if (Math.abs(shiftY) > 0.005 || Math.abs(shiftX) > 0.005) {
      camera.setViewOffset(w, h, (shiftX * w) / 2, (shiftY * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.6;
    fog.far = dist * 3.6;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */
  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const truckBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-wo-trucks]"));
  const keepBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-wo-keep]"));
  const dayBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-wo-day]"));
  const playBtn = q<HTMLButtonElement>("[data-wo-play]");
  const cushIn = q<HTMLInputElement>("[data-wo-cushion]");
  const cushVal = q("[data-wo-cushion-val]");
  const statusEl = q("[data-wo-status]");
  const chartEl = q("[data-wo-chart]");
  const hudClock = q("[data-wo-hud-clock]");
  const hudDoing = q("[data-wo-hud-doing]");
  const hudNote = q<HTMLElement>("[data-wo-hud-note]");
  const tiles = {
    thrown: { v: q("[data-wo-thrown]"), s: q("[data-wo-thrown-s]") },
    short: { box: q("[data-wo-short-box]"), v: q("[data-wo-short]"), s: q("[data-wo-short-s]") },
    full: { box: q("[data-wo-full-box]"), v: q("[data-wo-full]"), s: q("[data-wo-full-s]") },
    trucks: { v: q("[data-wo-trucks-t]"), s: q("[data-wo-trucks-s]") },
  };
  const tagCooler = q<HTMLElement>('[data-wo-tag="walkin"]');
  const tagFreezer = q<HTMLElement>('[data-wo-tag="freezer"]');
  const tagFloor = q<HTMLElement>('[data-wo-tag="floor"]');

  let policy: Policy = { trucks: 3, frozen: false, cushion: 0.1 };
  const keyOf = (p: Policy) => `${p.trucks}|${p.frozen}|${Math.round(p.cushion * 100)}`;
  const films = new Map<string, Film>();
  const filmFor = (p: Policy) => {
    let f = films.get(keyOf(p));
    if (!f) films.set(keyOf(p), (f = makeFilm(p)));
    return f;
  };
  let film = filmFor(policy);
  const manys = new Map<string, ManyWeeks>();
  let manyTimer = 0;
  let manyRun: { key: string; run: ReturnType<typeof manyWeeksRun> } | null = null;
  // The setting a reader arrives to, worked out now: the status is whole from the first, and the page below never
  // moves when a sentence lands later.
  manys.set(keyOf(policy), manyWeeks(policy, 300));
  const reduce = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ?? false;
  let playing = !reduce;
  let held = false;
  /** Film seconds. With motion reduced, it waits at Wednesday's prep done. */
  let ft = 0;
  const jobAt = (f: Film, day: number, kind: string) => f.jobs.find((j) => j.staged.ev.day === day && j.staged.ev.kind === kind);
  if (!playing) ft = (jobAt(film, 2, "prep")?.t1 ?? 0) + 0.05;

  /** The film second a clock hour falls at (the start of the job then under way, or between jobs). */
  function filmAt(f: Film, hour: number): number {
    let t = 0;
    let h = 0;
    for (const j of f.jobs) {
      if (hour <= j.hour) return t + (hour - h) / IDLE_HPS;
      const end = j.hour + j.span / 3600;
      if (hour < end) return j.t0 + ((hour - j.hour) * 3600) / j.lapse;
      t = j.t1;
      h = end;
    }
    return Math.min(f.length, t + (hour - h) / IDLE_HPS);
  }
  function setPolicy(p: Policy): void {
    const hour = shown?.hour ?? 0;
    policy = { ...p };
    film = filmFor(policy);
    ft = filmAt(film, hour);
    renderAll();
  }
  /** The 300 weeks for the setting, if worked out yet; else (a setting the reader picked) they are, from a moment
   *  later (a slider's drag settled), a few milliseconds at a time between frames, and then the tiles and status are
   *  said again. */
  function manyNow(): ManyWeeks | null {
    const m = manys.get(keyOf(policy));
    if (m) return m;
    if (!manyTimer) manyTimer = window.setTimeout(manySlice, 80);
    return null;
  }
  function manySlice(): void {
    manyTimer = 0;
    const k = keyOf(policy);
    if (manys.has(k)) return;
    if (manyRun?.key !== k) manyRun = { key: k, run: manyWeeksRun(policy, 300) };
    const t0 = performance.now();
    let done = false;
    while (!done && performance.now() - t0 < 8) done = manyRun.run.step(1);
    if (!done) {
      manyTimer = window.setTimeout(manySlice, 0);
      return;
    }
    manys.set(k, manyRun.run.result());
    manyRun = null;
    renderTiles();
    renderStatus();
  }

  function renderControls(): void {
    for (const b of truckBtns) b.setAttribute("aria-pressed", String(Number(b.dataset.woTrucks) === policy.trucks));
    for (const b of keepBtns) b.setAttribute("aria-pressed", String((b.dataset.woKeep === "frozen") === policy.frozen));
    if (cushIn) {
      cushIn.value = String(Math.round(policy.cushion * 100));
      cushIn.style.setProperty("--p", `${((policy.cushion / 0.25) * 100).toFixed(2)}%`);
    }
    if (cushVal) cushVal.textContent = `${Math.round(policy.cushion * 100)}% over`;
    if (playBtn) playBtn.textContent = playing ? "Pause" : "Play";
  }
  /** The most of a box's shelf places in use at once over the week (its floor's cases not among them). */
  const peakOf = (L: Layout, b: Where) => {
    let most = 0;
    for (const st of [{ after: L.start }, ...L.events]) most = Math.max(most, st.after.filter((c) => c.at.box === b && !c.at.floor).length);
    return most;
  };
  /** The most cases on the floors at once, said: "+ 2 on the walk-in's floor", or "+ 2 on the walk-in's floor, 4 on the
   *  freezer's". */
  const floorSaid = (L: Layout) => {
    const { walkin, freezer } = L.floorMost;
    return walkin && freezer ? `+ ${walkin} on the walk-in's floor, ${freezer} on the freezer's` : `+ ${walkin || freezer} on the ${walkin ? "walk-in's" : "freezer's"} floor`;
  };
  /** The share of the 300 weeks with a lunch short, said: "6% of weeks", or "1 in 300 weeks" under half a percent. */
  const shareSaid = (a: number) => (a > 0 && a < 0.005 ? `${Math.round(a * 300)} in 300 weeks` : `${Math.round(a * 100)}% of weeks`);
  /** The day a box first has cases on its floor, what it holds then, and why they are there. */
  const crowdSaid = (b: Where): string | null => {
    const st = film.layout.events.find((x) => x.after.some((c) => c.at.box === b && c.at.floor));
    if (!st) return null;
    const day = st.ev.day;
    const n = film.layout.floorMost[b];
    const onFloor = `${n === 1 ? "1 case stands" : `up to ${n} cases stand`} on its floor`;
    if (b === "freezer") {
      const k = st.after.filter((c) => c.at.box === b).length;
      const next = TRUCK_DAYS[policy.trucks].find((d) => d > day) ?? TRUCK_DAYS[policy.trucks][0]!;
      return `On ${DAY[day]} the freezer holds ${k} cases, all the frozen meat until ${DAY[next]}'s truck, and ${onFloor}.`;
    }
    // The lunches the walk-in's meat is for: each lot's lunches from then, or its ready day, to its last day.
    const lunches = new Set<number>();
    for (const id of new Set(st.after.filter((c) => c.at.box === b).map((c) => c.lot))) {
      const l = film.lots.get(id);
      if (!l) continue;
      for (let u = Math.max(day + 7, l.readyDay); u <= l.lastDay; u++) if (MENU.some((mm) => mm.day === u % 7 && mm.item === l.item)) lunches.add(u - 7);
    }
    const said = [...lunches].sort((x, y) => x - y).map((u) => `${u >= 7 ? "next " : ""}${DAY[u % 7]}'s`);
    return `On ${DAY[day]} the walk-in holds the meat for ${said.length ? `${listAnd(said)} lunch${said.length === 1 ? "" : "es"}` : "the week"} at once, and ${onFloor}.`;
  };
  function renderTiles(): void {
    const w = film.week;
    const m = manyNow();
    if (tiles.thrown.v) tiles.thrown.v.textContent = lb(w.thrownLb);
    if (tiles.thrown.s) tiles.thrown.s.textContent = `this week\n1 week in 10: over ${m ? lb(m.thrownLb.p90) : "…"}`;
    const short = w.shortLunches;
    tiles.short.box?.setAttribute("data-ok", String(short === 0));
    if (tiles.short.v) tiles.short.v.textContent = short ? `${short.toLocaleString("en-US")} lunches` : "None";
    if (tiles.short.s) tiles.short.s.textContent = `${short ? "bought on the day" : "short this week"}\n${m ? shareSaid(m.anyShort) : "…"} short`;
    const wi = peakOf(film.layout, "walkin");
    const fz = peakOf(film.layout, "freezer");
    const floor = film.layout.floorMost.walkin + film.layout.floorMost.freezer;
    tiles.full.box?.setAttribute("data-ok", String(floor === 0));
    if (tiles.full.v) tiles.full.v.textContent = policy.frozen ? `${fz} of ${FREEZER_PLACES}` : `${wi} of ${WALKIN_PLACES}`;
    const which = policy.frozen ? `the freezer's, at its fullest\nthe walk-in's: ${wi} of ${WALKIN_PLACES}` : "the walk-in's, at its fullest";
    if (tiles.full.s) tiles.full.s.textContent = floor ? `${which}\n${floorSaid(film.layout)}` : which;
    const days = TRUCK_DAYS[policy.trucks].map((d) => DAY3[d]).join(", ");
    if (tiles.trucks.v) tiles.trucks.v.textContent = `${policy.trucks} a week`;
    if (tiles.trucks.s) tiles.trucks.s.textContent = `${days}\nordered by 2 pm the day before`;
  }

  /** The week, said. */
  function renderStatus(): void {
    if (!statusEl) return;
    const w = film.week;
    const m = manyNow();
    const parts: string[] = [];
    const shortDays = w.events.filter((e) => e.kind === "prep" && (e.shortLunches ?? 0) > 0);
    const trucks = TRUCK_DAYS[policy.trucks];
    if (!policy.frozen) {
      if (shortDays.length) {
        const d = shortDays[0]!;
        const it = ITEM[d.item!];
        const before = [...trucks].filter((t) => t < d.day).pop();
        parts.push(
          `${DAY[d.day]}'s ${d.lunches!.toLocaleString("en-US")} lunches of ${it.dish.toLowerCase()} come up short by ${d.shortLunches!.toLocaleString("en-US")}: ${SHORT_NAME[d.item!]} keep ${it.freshDays} days in the walk-in, counting the day they come, so ${before === undefined ? "no truck this week" : `${DAY[before]}'s truck`} can't carry them that far, and the next truck comes too late. They are bought on the day.`,
        );
      } else if (policy.trucks === 5) {
        parts.push(`A truck every weekday brings each day's meat on the morning it is cooked, ordered at 2 pm the day before with the final count in hand: nothing short, and the cushion never comes into it.`);
      } else {
        parts.push(`Every lunch is covered off the trucks: each brings the day's final count and the plan, plus a ${Math.round(policy.cushion * 100)}% cushion, for any lunch before the next truck that its meat will last to.`);
      }
    } else {
      const lead = (n: number) => listAnd(ITEMS.filter((i) => i.thawDays === n).map((i) => SHORT_NAME[i.key]));
      const thaw = `The thaw starts at 8 am, a day ahead of the lunch for ${lead(1)} and two for the 10 lb packs of ${lead(2)}, before the final count is in`;
      parts.push(
        shortDays.length
          ? `${thaw}, so it goes by the plan: ${listAnd(shortDays.map((d) => `${DAY[d.day]} comes up ${lunchesN(d.shortLunches!)} short`))} with a ${Math.round(policy.cushion * 100)}% cushion.`
          : `Frozen keeps for months, so the trucks only decide how full the freezer gets. ${thaw}, so it goes by the plan and its ${Math.round(policy.cushion * 100)}% cushion, a whole case at a time; once thawed it keeps no longer than fresh.`,
      );
    }
    const thrownBy = ITEMS.filter((i) => w.thrownBy[i.key] > 0.5).map((i) => `${lb(w.thrownBy[i.key])} of ${SHORT_NAME[i.key]}`);
    parts.push(thrownBy.length ? `Thrown out this week, each at the end of its last day: ${listAnd(thrownBy)}, ${lb(w.thrownLb)} in all.` : "Nothing was thrown out this week.");
    if (m) parts.push(`Over 300 weeks of counts scattered around the plan, two in three within 10% of it, half throw out more than ${lb(m.thrownLb.mid)} and 1 in 10 more than ${lb(m.thrownLb.p90)}; ${m.anyShort === 0 ? "none comes up short" : m.anyShort < 0.005 ? `${Math.round(m.anyShort * 300)} in 300 come up short somewhere` : `${Math.round(m.anyShort * 100)}% come up short somewhere`}.`);
    const floor = film.layout.floorMost.walkin + film.layout.floorMost.freezer;
    if (floor) {
      for (const b of ["walkin", "freezer"] as const) {
        const said = crowdSaid(b);
        if (said) parts.push(said);
      }
      parts.push(`No shelf place is left that keeps chicken lowest and fish over beef, and a case goes only on its own delivery's, so those stand on the bare floor, where the Food Code allows no case of meat: it wants food at least 6 inches up.`);
    }
    statusEl.textContent = parts.join(" ");
    statusEl.dataset.state = shortDays.length || floor ? "warn" : "ok";
  }

  /* The chart: the week's lots, each a bar from the day it comes in to its last day, item by item; the lunches. */
  function renderChart(): void {
    if (!chartEl) return;
    const W = Math.max(300, Math.round(chartEl.clientWidth));
    const narrow = W < 540;
    const L = narrow ? 78 : 112;
    const R = W - 8;
    const T = 22;
    const laneH = narrow ? 9 : 11;
    const dw = (R - L) / 7;
    const x = (d: number) => L + d * dw;
    const w = film.week;
    /** The (absolute) day a lot's last case leaves, by the week's events; 14 if it is there at the week's end. */
    const goneOn = (l: Lot) => {
      let seen = w.start.some((x) => x.id === l.id);
      for (const e of w.events) {
        if (e.stock.some((x) => x.id === l.id)) seen = true;
        else if (seen) return e.day + 7;
      }
      return 14;
    };
    const rows = ITEMS.map((it) => {
      const lots = [...film.lots.values()].filter((l) => l.item === it.key && (l.lastDay >= 7 || l.where === "freezer") && l.inDay < 14);
      return { it, lots };
    });
    // Lanes in each row: a lot in the first lane it fits.
    let y = T;
    let svg = "";
    const body: string[] = [];
    for (const { it, lots } of rows) {
      const lanes: number[] = [];
      const placed: { l: Lot; lane: number }[] = [];
      for (const l of lots.sort((a, b) => a.inDay - b.inDay)) {
        const a = Math.max(7, l.inDay);
        let lane = lanes.findIndex((end) => end < a);
        if (lane < 0) lane = lanes.push(0) - 1;
        lanes[lane] = l.where === "freezer" ? goneOn(l) : l.lastDay;
        placed.push({ l, lane });
      }
      const n = Math.max(1, lanes.length);
      const h = n * (laneH + 3) + 16;
      body.push(`<text class="ya it" x="${L - 8}" y="${(y + 12).toFixed(1)}" text-anchor="end">${narrow ? SHORT_NAME[it.key].replace("chicken ", "") : it.name}</text>`);
      body.push(`<line class="grid" x1="${L}" x2="${R}" y1="${(y + h).toFixed(1)}" y2="${(y + h).toFixed(1)}" />`);
      for (const { l, lane } of placed) {
        const a = Math.max(0, l.inDay - 7);
        const b = Math.min(7, (l.where === "freezer" ? goneOn(l) : l.lastDay) - 7 + 1);
        const yy = y + 6 + lane * (laneH + 3);
        const cls = l.where === "freezer" ? "lot frozen" : l.state === "thawing" || l.state === "thawed" ? "lot thawed" : "lot";
        body.push(`<rect class="${cls}" x="${(x(a) + 1).toFixed(1)}" y="${yy.toFixed(1)}" width="${Math.max(2, x(b) - x(a) - 2).toFixed(1)}" height="${laneH}" rx="2" style="--c:${ITEM_C[it.key]}" />`);
        const out = w.events.find((e) => e.kind === "check" && e.lots.some((x2) => x2.id === l.id));
        if (out) {
          const lbOut = out.lots.find((x2) => x2.id === l.id)!.lb;
          body.push(`<text class="out" x="${(x(out.day + 1) - 3).toFixed(1)}" y="${(yy + laneH - 1.5).toFixed(1)}" text-anchor="end">✕ ${narrow ? Math.round(lbOut) : lb(lbOut)}</text>`);
        }
      }
      // The lunches: the day's final count on its day, and what came up short.
      for (const mm of MENU.filter((mm) => mm.item === it.key)) {
        const ev = w.events.find((e) => e.kind === "prep" && e.day === mm.day);
        const cx = x(mm.day) + dw / 2;
        const sy = y + h - 4;
        const shortN = ev?.shortLunches ?? 0;
        body.push(`<text class="lunch${shortN ? " short" : ""}" x="${cx.toFixed(1)}" y="${sy.toFixed(1)}" text-anchor="middle">${shortN ? `${narrow ? "" : "short "}${shortN}` : `${THIS_WEEK[mm.day]} lunches`.replace(" lunches", narrow ? "" : " lunches")}</text>`);
      }
      y += h;
    }
    const H = y + 22;
    svg += `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg wo-svg">`;
    for (let d = 0; d < 7; d++) {
      svg += `<g class="day" data-wo-col="${d}"><rect class="hit" x="${x(d).toFixed(1)}" y="${T - 18}" width="${dw.toFixed(1)}" height="${(H - T + 18).toFixed(1)}" />`;
      svg += `<text class="xa${trucksOn(d) ? " truck" : ""}" x="${(x(d) + dw / 2).toFixed(1)}" y="${T - 7}" text-anchor="middle">${narrow ? DAY3[d]!.slice(0, 2) : DAY3[d]}${trucksOn(d) ? " ▲" : ""}</text></g>`;
      svg += `<line class="grid v" x1="${x(d).toFixed(1)}" x2="${x(d).toFixed(1)}" y1="${T - 2}" y2="${y.toFixed(1)}" />`;
    }
    svg += body.join("");
    svg += `<line class="now" data-wo-now x1="${L}" x2="${L}" y1="${T - 2}" y2="${y.toFixed(1)}" />`;
    svg += `<text class="xa sub" x="${(L + (R - L) / 2).toFixed(1)}" y="${H - 6}" text-anchor="middle">${narrow ? "▲ truck · bars: in to last day" : "▲ a truck comes · each bar: a lot, from the day it comes in to its last day · ✕ thrown out"}</text>`;
    svg += `</svg>`;
    chartEl.innerHTML = svg;
    nowLine = chartEl.querySelector("[data-wo-now]");
    chartX = { L, dw };
    for (const g of chartEl.querySelectorAll<SVGGElement>("[data-wo-col]")) {
      g.addEventListener("click", () => seekDay(Number(g.dataset.woCol)));
    }
  }
  const trucksOn = (d: number) => TRUCK_DAYS[policy.trucks].includes(d);
  let nowLine: SVGLineElement | null = null;
  let chartX = { L: 0, dw: 0 };
  function renderAll(): void {
    renderControls();
    renderTiles();
    renderStatus();
    renderChart();
    lastHud = "";
  }
  function seekDay(d: number): void {
    // Just before the day's first job.
    const first = film.jobs.find((j) => j.staged.ev.day === d);
    ft = first ? Math.max(0, first.t0 - 0.4) : filmAt(film, d * 24 + 5);
    if (reduce) playing = false;
    renderControls();
  }

  /* The readouts on the stage. */
  let lastHud = "";
  function doingOf(m: Moment): { doing: string; note: string; kind: string } {
    const j = m.job;
    if (!j) {
      const next = film.jobs.find((x) => x.hour > m.hour);
      if (!next) return { doing: "The week is done", note: "", kind: "" };
      const time = clockText(next.hour).replace(/^[A-Za-z]+, /, "");
      const when = Math.floor(next.hour / 24) === Math.floor(m.hour / 24) ? `at ${time}` : `${DAY[next.staged.ev.day]} at ${time}`;
      const what = whatOf(next);
      return { doing: `Next, ${when}, ${what[0]!.toLowerCase()}${what.slice(1)}`, note: "", kind: "" };
    }
    const note = j.routine ? `Time-lapse ×${Math.round(j.lapse)}` : "";
    const ev = j.staged.ev;
    if (ev.kind === "prep" && (ev.shortLunches ?? 0) > 0) return { doing: whatOf(j), note: `${lunchesN(ev.shortLunches!)} short: bought on the day`, kind: "short" };
    if (ev.kind === "check") return { doing: whatOf(j), note: "Last day over: to the bin", kind: "past" };
    return { doing: whatOf(j), note, kind: "" };
  }
  function whatOf(j: Job): string {
    const ev = j.staged.ev;
    const byItem = (lots: { item: ItemKey; lb: number }[], count: (it: ItemKey) => number) =>
      listAnd(ITEMS.filter((i) => lots.some((l) => l.item === i.key)).map((i) => cases(count(i.key), SHORT_NAME[i.key])));
    const shelfMove = (mv: Staged["moves"][number]) => typeof mv.from === "object" && typeof mv.to === "object" && mv.from.box === mv.to.box;
    const moved = (it: ItemKey) =>
      j.staged.moves.filter((mv) => mv.item === it && !shelfMove(mv) && (mv.from === "dock" || mv.to === "kitchen" || (typeof mv.from === "object" && mv.from.box === "freezer"))).length;
    const shifted = j.staged.moves.filter(shelfMove).length;
    const room = shifted ? `, with ${shifted} case${shifted === 1 ? "" : "s"} moved a shelf to make room` : "";
    if (ev.kind === "truck") return ev.lots.length ? `The truck: ${byItem(ev.lots, moved)}${policy.frozen ? ", into the freezer" : ""}${room}` : `The truck: nothing it needs to bring`;
    if (ev.kind === "thaw") {
      const each = ITEMS.filter((i) => ev.lots.some((l) => l.item === i.key)).map((i) => {
        const l = film.lots.get(ev.lots.find((x) => x.item === i.key)!.id);
        return `${cases(moved(i.key), SHORT_NAME[i.key])}${l ? ` for ${DAY[l.readyDay % 7]}` : ""}`;
      });
      return `Into the walk-in to thaw: ${listAnd(each)}${room}`;
    }
    if (ev.kind === "prep") {
      const n = j.staged.moves.length;
      return `Prep for ${ev.lunches!.toLocaleString("en-US")} lunches of ${ITEM[ev.item!].dish.toLowerCase()}${n ? `: ${n} case${n === 1 ? "" : "s"} out` : ""}`;
    }
    const n = j.staged.moves.length;
    return `The night's check: ${n ? `${n} case${n === 1 ? "" : "s"} on ${n === 1 ? "its" : "their"} last day, ${lb(ev.lots.reduce((s, l) => s + l.lb, 0))}, out` : "nothing on its last day"}`;
  }
  function showHud(m: Moment): void {
    const { doing, note, kind } = doingOf(m);
    const clock = clockText(m.hour);
    const sig = `${clock}|${doing}|${note}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = clock;
    if (hudDoing) hudDoing.textContent = doing;
    if (hudNote) {
      hudNote.textContent = note;
      hudNote.hidden = !note;
      hudNote.dataset.kind = kind;
    }
    for (const b of dayBtns) b.setAttribute("aria-pressed", String(Number(b.dataset.woDay) === Math.min(6, Math.floor(m.hour / 24))));
  }
  const tagAt = new THREE.Vector3();
  const tagSize = new Map<HTMLElement, { w: number; h: number }>();
  const STEM = 22;
  function placeTag(el: HTMLElement | null, at: THREE.Vector3 | null, text: [string, string], warn: boolean, colour: string, W: number, H: number): void {
    if (!el) return;
    let on = !!at;
    if (at) {
      tagAt.copy(at).project(camera);
      on = tagAt.z < 1 && Math.abs(tagAt.x) < 0.96 && Math.abs(tagAt.y) < 0.96;
    }
    if (el.hidden === on) el.hidden = !on;
    if (!on) return;
    const sig = text.join("|") + warn;
    if (el.dataset.last !== sig || !tagSize.has(el)) {
      el.dataset.last = sig;
      const i = el.querySelector("i");
      const b = el.querySelector("b");
      if (i) i.textContent = text[0];
      if (b) b.textContent = text[1];
      el.toggleAttribute("data-boiled", warn);
      el.style.setProperty("--c", colour);
      el.style.setProperty("--stem", `${STEM}px`);
      tagSize.set(el, { w: el.offsetWidth, h: el.offsetHeight });
    }
    const sz = tagSize.get(el)!;
    const px = ((tagAt.x + 1) / 2) * W;
    const py = ((1 - tagAt.y) / 2) * H;
    const x = Math.min(W - sz.w / 2 - 4, Math.max(sz.w / 2 + 4, px));
    const tf = `translate(${x.toFixed(1)}px, ${(py - STEM).toFixed(1)}px)`;
    if (el.dataset.tf !== tf) {
      el.dataset.tf = tf;
      el.style.transform = tf;
    }
    const sx = `${Math.round(px - x)}px`;
    if (el.dataset.sx !== sx) {
      el.dataset.sx = sx;
      el.style.setProperty("--sx", sx);
    }
  }
  const coolerTop = new THREE.Vector3(X0.walkin + 0.4, CUT_Y + 0.01, FRONT_OUT - PANEL / 2);
  const freezerTop = new THREE.Vector3(X0.freezer + COOLER.w - 0.4, CUT_Y + 0.01, FRONT_OUT - PANEL / 2);
  const floorAt = new THREE.Vector3();
  function showTags(m: Moment): void {
    const W = stageEl!.clientWidth;
    const H = stageEl!.clientHeight;
    const count = (b: Where) => m.cases.filter((c) => c.at?.box === b && !c.at.floor).length;
    // On a narrow stage the fronts' outer corners are out of view, and the tiles say as much.
    placeTag(tagCooler, coolerTop, ["Walk-in, 38 °F", `${count("walkin")} of ${WALKIN_PLACES} places`], false, "#7dd3fc", W, H);
    placeTag(tagFreezer, freezerTop, ["Freezer, 0 °F", `${count("freezer")} of ${FREEZER_PLACES} places`], false, FROZEN, W, H);
    const onFloor = m.cases.find((c) => c.at?.floor);
    if (onFloor?.at) {
      spotPose(onFloor.at, onFloor.item, floorAt);
      floorAt.y += 0.05;
    }
    const n = m.cases.filter((c) => c.at?.floor).length;
    // Its first line is hidden on a phone, so the second says where they are.
    placeTag(tagFloor, onFloor?.at ? floorAt : null, ["No shelf place", `${n} case${n === 1 ? "" : "s"} on the floor`], true, PAST, W, H);
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* Running it.                                                                                                 */
  /* ---------------------------------------------------------------------------------------------------------- */
  renderAll();
  frame();
  new ResizeObserver(() => {
    if (stageEl.clientWidth === fitW && stageEl.clientHeight === fitH) return;
    tagSize.clear();
    frame();
  }).observe(stageEl);
  void document.fonts?.ready.then(() => {
    tagSize.clear();
    frame();
  });
  if (chartEl) {
    let lastW = chartEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(chartEl.clientWidth - lastW) < 4) return;
      lastW = chartEl.clientWidth;
      renderChart();
    }).observe(chartEl);
  }
  for (const b of truckBtns) b.addEventListener("click", () => setPolicy({ ...policy, trucks: Number(b.dataset.woTrucks) as Trucks }));
  for (const b of keepBtns) b.addEventListener("click", () => setPolicy({ ...policy, frozen: b.dataset.woKeep === "frozen" }));
  cushIn?.addEventListener("input", () => {
    const c = Math.max(0, Math.min(25, Math.round(Number(cushIn.value) / 5) * 5)) / 100;
    if (c !== policy.cushion) setPolicy({ ...policy, cushion: c });
  });
  for (const b of dayBtns) b.addEventListener("click", () => seekDay(Number(b.dataset.woDay)));
  playBtn?.addEventListener("click", () => {
    playing = !playing;
    if (playing && ft >= film.length - 0.05) ft = 0;
    renderControls();
  });
  void cook.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    if (playing && !held) {
      ft += dt;
      // Round again from Monday.
      if (ft >= film.length) ft = 0;
    }
    const m = momentAt(film, ft);
    showMoment(film, m);
    showHud(m);
    if (nowLine) {
      const nx = (chartX.L + (m.hour / 24) * chartX.dw).toFixed(1);
      if (nowLine.getAttribute("x1") !== nx) {
        nowLine.setAttribute("x1", nx);
        nowLine.setAttribute("x2", nx);
      }
    }
    controls.update();
    fitClip();
    renderer.render(scene, camera);
    showTags(m);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    walkin: {
      info: () => {
        const m = shown;
        const byBox = (b: Where) => (m ? m.cases.filter((c) => c.at?.box === b && !c.at.floor).length : 0);
        return {
          mountMs: Math.round(mountMs),
          running: inView && !document.hidden && gate.ready(),
          policy: { ...policy },
          playing,
          ft: +ft.toFixed(2),
          length: +film.length.toFixed(1),
          hour: m ? +m.hour.toFixed(3) : null,
          job: m?.job ? { kind: m.job.staged.ev.kind, day: m.job.staged.ev.day, moves: m.job.staged.moves.length, lapse: +m.job.lapse.toFixed(1) } : null,
          cook: cook.group.visible ? { x: +cook.placement.x.toFixed(2), z: +cook.placement.z.toFixed(2) } : null,
          inHand: m?.inHand ?? null,
          shelved: { walkin: byBox("walkin"), freezer: byBox("freezer") },
          floor: m ? m.cases.filter((c) => c.at?.floor).length : 0,
          week: { thrownLb: Math.round(film.week.thrownLb), shortLunches: film.week.shortLunches, floorMost: film.layout.floorMost, shelvedOut: film.layout.shelvedOut, jobs: film.jobs.length },
          doors: { ...doorK },
          tags: [tagCooler, tagFreezer, tagFloor].map((el) => (el && !el.hidden ? el.textContent : null)),
          hud: hudEl?.textContent ?? "",
          status: statusEl?.textContent ?? "",
        };
      },
      /** Set the trucks, fresh or frozen, the cushion. */
      set(p: Partial<Policy>) {
        setPolicy({ ...policy, ...p });
      },
      /** Show the film at a clock hour of the week (0 = Monday midnight), held there. */
      seek(hour: number) {
        ft = filmAt(film, hour);
      },
      /** Hold the film where it is (true), or let it run. */
      hold(on: boolean) {
        held = on;
      },
      /** The cases' places in the job under way at a clock hour, without moving the film. */
      layout: () => film.layout.events.map((st) => ({ kind: st.ev.kind, day: st.ev.day, moves: st.moves.length, after: st.after.length })),
      view(eye: [number, number, number] | null, at?: [number, number, number]) {
        userMoved = !!eye;
        if (eye && at) {
          camera.clearViewOffset();
          camera.position.set(...eye);
          controls.target.set(...at);
          controls.update();
        } else frame();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
