/**
 * Cooked the day before, cooled in time (chapter 06, second figure): 30 litres of chili off the stove at 2 pm, cooled
 * four ways in a walk-in cooler and the walk-in freezer beside it. Left in its 40-quart pot on a dunnage rack; in
 * 3-inch hotel pans on a speed rack; in 2-inch pans on another; and in 2-inch pans started in the freezer, wheeled to
 * the walk-in once their middle is down to 41 °F.
 *
 * The front walls are cut away low, and everything holding chili is cut through its middle: each cut face is coloured
 * by the temperature sim/coolDown.ts computes there, through the afternoon and the night. Beside the stage, each way's
 * hours against the Food Code's two limits; under it, the four cuts side by side to scale, and each way's middle
 * against the clock that starts when it passes 135 °F.
 */
import * as THREE from "three";
import { RoundedBoxGeometry } from "three/examples/jsm/geometries/RoundedBoxGeometry.js";
import {
  COLD_C,
  DAY_MIN,
  HOT_C,
  MAX_LITRES,
  MIN_LITRES,
  PAN_2,
  PAN_3,
  SAMPLE_MIN,
  SNAP_MIN,
  STAGE1_MIN,
  TOTAL_MIN,
  WARM_C,
  WAYS,
  cToF,
  centreAt,
  coolAll,
  coolKey,
  pansFor,
  potFill,
  type CoolRun,
  type WayKey,
  type WayRun,
} from "../sim/coolDown";
import type { CoolAsk, CoolDone } from "../sim/coolRuns.worker";
import { sheetPan } from "./coldChainLab";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { buildFixture, type Model, type Piece } from "./models/all";
import { speedRackSlides } from "./models/cold";
import { palette } from "./models/materials";
import { makeStaff } from "./models/people";
import { Shop } from "./models/shop";
import { hotelPan, hotelPanInside } from "./models/smallwares";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

const COOLER = { id: "NSF-WALKIN-COOLER", w: 2.438, d: 3.048, h: 2.438 };
const FREEZER_BOX = { id: "NSF-WALKIN-FREEZER-8X10", w: 2.438, d: 3.048, h: 2.438 };
const RACK = { id: "NSF-SPEEDRACK-20PAN", w: 0.521, d: 0.66, h: 1.753 };
const POT = { id: "NSF-STOCKPOT-40QT", w: 0.38, d: 0.46, h: 0.375 };

/* Where everything stands. The two boxes side by side, their fronts at z = d/2 (the front walls' outer faces 87 mm back
 * of that), as one combination box: the freezer's left wall is the cooler's right. Shelving runs down both sides of
 * each, stocked; the racks and the pot stand in the aisles. */
const GAP = 0.004;
const PANEL = 0.1016;
const COOLER_X0 = -GAP / 2 - COOLER.w;
const FREEZER_X0 = GAP / 2;
const FRONT_OUT = COOLER.d / 2 - 0.087;
const FRONT_IN = FRONT_OUT - PANEL;
/** The freezer stands on its own insulated floor. */
const FZ_FLOOR = PANEL;
/** The front walls are cut away above this. */
const CUT_Y = 0.3;
/** The racks' row: each rack's pans are cut through their middles, on this plane. */
const ROW_Z = -0.55;
/** The racks' x. The freezer pans' rack stands in line with the freezer's door, and once moved, at the cooler's right. */
const RACK_X = { deep: -1.56, shallow: -0.98, moved: -0.4, freezer: FREEZER_X0 + 0.76 + 0.457 } as const;
/** The pot, on an aluminium dunnage rack 8 in off the floor, its handles fore and aft; cut through its axis. */
const POT_AT = { x: -1.45, z: 0.6 };
const DUN = { w: 0.61, d: 0.46, h: 0.2 };
const RACK_SLIDES = speedRackSlides(RACK.h);
/** A sheet pan's rim rides the slide; the hotel pan's floor stands this far below it, on the sheet pan's floor. */
const PAN_DROP = -0.0236;
/** Slides apart, from the top: a 2-inch pan on every other, a 3-inch on every third, air between. */
const EVERY: Record<PanWay, number> = { deep: 3, shallow: 2, freezer: 2 };
type PanWay = Exclude<WayKey, "pot">;
const PAN_OF: Record<PanWay, typeof PAN_2> = { deep: PAN_3, shallow: PAN_2, freezer: PAN_2 };
/** The cook, outside, between the two doors. */
const COOK = { x: 0, z: 2.05 };

/** Each way's colour, on the labels, the tiles and the chart. Kept off the temperature colours. */
export const WAY_COLOR: Record<WayKey, string> = { pot: "#f472b6", deep: "#c084fc", shallow: "#4ade80", freezer: "#e2e8f0" };

/** The boxes, the racks and the pot, for main.ts to build ahead. */
export function coolPieces(): readonly Piece[] {
  const d = roomDetail();
  return [[COOLER.id, COOLER, d], [FREEZER_BOX.id, FREEZER_BOX, d], ...Array.from({ length: 3 }, (): Piece => [RACK.id, RACK, d]), [POT.id, POT, 2]];
}
/** The rooms close up on a laptop, a step coarser on a phone. */
const roomDetail = () => (Math.min(window.innerWidth, window.innerHeight) < 600 ? 1 : 2);

/** The day: the pans go in at 2 pm (minutes after midnight); the stage plays 24 hours from then. */
const START_TOD = 14 * 60;
/** The loop: going in, the day, the move from the freezer (the clock held), and a hold on the end (seconds). */
const IN_S = 1.6;
const DAY_S = 34;
const MOVE_S = 2.4;
const END_S = 4;
/** The chart's hours from 135 °F. */
const CHART_H = 24;

/* ------------------------------------------------------------------------------------------------------------ */
/* Words and colours.                                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

const degF = (c: number) => `${Math.round(cToF(c))} °F`;
/** Metres in inches: "5.3 in". */
const inches = (m: number) => `${(m / 0.0254).toFixed(m < 0.1 ? 1 : 0)} in`;
/** Each way under its cut on a phone. */
const TINY: Record<WayKey, string> = { pot: "Pot", deep: "3 in", shallow: "2 in", freezer: "Freezer" };
/** Each way in a sentence. */
const PHRASE: Record<WayKey, string> = { pot: "the chili left in its pot", deep: "3-inch pans", shallow: "2-inch pans", freezer: "2-inch pans started in the freezer" };
/** "9 h 22" or "45 min". */
const span = (m: number) => {
  const r = Math.round(m);
  return r >= 60 ? `${Math.floor(r / 60)} h${r % 60 ? ` ${String(r % 60).padStart(2, "0")}` : ""}` : `${r} min`;
};
/** The time of day `m` minutes after the pans went in: "2:00 pm". */
const clockAt = (m: number) => {
  const a = (START_TOD + Math.round(m)) % 1440;
  const h = Math.floor(a / 60);
  const mm = String(a % 60).padStart(2, "0");
  return `${h % 12 === 0 ? 12 : h % 12}:${mm} ${h < 12 ? "am" : "pm"}`;
};
const smooth = (x: number) => {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
};

/** The temperature colours (°F): deep blue frozen, sky at 41 °F, yellow at 70 °F, orange at 135 °F, red off the stove. */
const HEAT: readonly (readonly [number, number, number, number])[] = [
  [-10, 0x1e, 0x3a, 0x8a],
  [32, 0x25, 0x63, 0xeb],
  [41, 0x38, 0xbd, 0xf8],
  [70, 0xfa, 0xcc, 0x15],
  [135, 0xf9, 0x73, 0x16],
  [190, 0xdc, 0x26, 0x26],
];
function heatRGB(c: number, out: Uint8Array | Uint8ClampedArray, k: number): void {
  const f = cToF(c);
  let i = 0;
  while (i < HEAT.length - 2 && f > HEAT[i + 1]![0]) i++;
  const a = HEAT[i]!;
  const b = HEAT[i + 1]!;
  const t = Math.min(1, Math.max(0, (f - a[0]) / (b[0] - a[0])));
  out[k] = a[1] + (b[1] - a[1]) * t;
  out[k + 1] = a[2] + (b[2] - a[2]) * t;
  out[k + 2] = a[3] + (b[3] - a[3]) * t;
  out[k + 3] = 255;
}
const heatCss = (f: number) => {
  const px = new Uint8Array(4);
  heatRGB((f - 32) / 1.8, px, 0);
  return `rgb(${px[0]} ${px[1]} ${px[2]})`;
};

/** A run's field at minute `m`, straight between its snapshots, into `out`. */
function fieldAt(run: CoolRun, m: number, out: Float32Array): void {
  const f = Math.max(0, m) / SNAP_MIN;
  const i = Math.min(run.snaps.length - 1, Math.floor(f));
  const j = Math.min(run.snaps.length - 1, i + 1);
  const t = Math.min(1, f - i);
  const a = run.snaps[i]!;
  const b = run.snaps[j]!;
  for (let k = 0; k < out.length; k++) out[k] = a[k]! + (b[k]! - a[k]!) * t;
}

/** The middle's minute at a temperature line, from the run's own marks. */
const lineMin = (r: CoolRun, line: "135" | "70" | "41") => (line === "135" ? r.t135 : line === "70" ? r.t70 : r.t41);

export function mountCoolLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-cd-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-cd-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  renderer.localClippingEnabled = true;
  const phone = window.matchMedia("(max-width: 720px)").matches;
  const detail = roomDetail();

  /* The scene: the walk-in cooler and the walk-in freezer side by side, their roofs off and their fronts cut low. */
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
    const floor = tilePlane(12, 8);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0.2, 0, 2.2);
    scene.add(floor);
  }

  // The front walls cut above CUT_Y: a fragment goes if it is both above the cut and in the front wall's thickness.
  const frontCut = [new THREE.Plane(new THREE.Vector3(0, -1, 0), CUT_Y), new THREE.Plane(new THREE.Vector3(0, 0, -1), FRONT_IN - 0.001)];
  // The racks' pans and the pot: their near halves.
  const rowCut = [new THREE.Plane(new THREE.Vector3(0, 0, -1), ROW_Z)];
  const potCut = [new THREE.Plane(new THREE.Vector3(0, 0, -1), POT_AT.z)];
  /** Every material under `o`, its own clipped copy: one copy per material and set of planes. */
  const clipCopies = new Map<THREE.Plane[], Map<THREE.Material, THREE.Material>>();
  function clip(o: THREE.Object3D, planes: THREE.Plane[], intersect = false): void {
    let made = clipCopies.get(planes);
    if (!made) clipCopies.set(planes, (made = new Map()));
    const copy = (m: THREE.Material) => {
      let c = made!.get(m);
      if (!c) {
        c = m.clone();
        c.clippingPlanes = planes;
        c.clipIntersection = intersect;
        c.clipShadows = true;
        made!.set(m, c);
      }
      return c;
    };
    o.traverse((x) => {
      const mesh = x as THREE.Mesh;
      // Steam stays whole.
      if (!mesh.isMesh && !(x as THREE.Points).isPoints) return;
      mesh.material = Array.isArray(mesh.material) ? mesh.material.map(copy) : copy(mesh.material);
    });
  }

  function box(spec: typeof COOLER, x0: number): Model {
    const m = buildFixture(spec.id, spec, detail);
    m.group.position.set(x0, 0, spec.d / 2);
    m.hooks.set?.("roof", 1);
    // The strip curtain hangs from the door's head, which is cut away.
    for (const [name, j] of m.joints) if (name.startsWith("pvc_strip_")) j.visible = false;
    clip(m.group, frontCut, true);
    scene.add(m.group);
    // The cut: the panels' foamed-in-place core between their skins.
    const cap = new THREE.Mesh(new THREE.BoxGeometry(spec.w - 0.002, 0.002, PANEL - 0.0018), foam);
    cap.position.set(x0 + spec.w / 2, CUT_Y - 0.0012, FRONT_OUT - PANEL / 2);
    cap.receiveShadow = true;
    scene.add(cap);
    // A cold white light over each box, where its ceiling's fixtures were.
    const lamp = new THREE.PointLight(0xeef4ff, 5.5, 5.5, 1.6);
    lamp.position.set(x0 + spec.w / 2, spec.h - 0.15, -0.2);
    scene.add(lamp);
    return m;
  }
  const foam = new THREE.MeshStandardMaterial({ color: 0xd8b75e, roughness: 0.93 });
  const cooler = box(COOLER, COOLER_X0);
  const freezer = box(FREEZER_BOX, FREEZER_X0);
  freezer.hooks.set?.("frost", 0.45);
  // One box of two: the cooler's own right wall, and the shelving on it, give way to the freezer's left wall.
  {
    const j = cooler.joints.get("part_wall_right");
    if (j) j.visible = false;
  }

  /* The racks, each loaded from the top. */
  interface PanSlot {
    group: THREE.Group;
    film: THREE.Mesh;
    cut: THREE.Mesh;
  }
  interface Rack {
    way: PanWay;
    model: Model;
    slots: PanSlot[];
    /** The rack's centre now (it moves, for the freezer's). */
    x: number;
    y: number;
  }
  const racks: Rack[] = (["deep", "shallow", "freezer"] as const).map((way) => {
    const model = buildFixture(RACK.id, RACK, detail);
    model.hooks.set?.("pans", 0);
    model.hooks.set?.("swivel", 0.25);
    scene.add(model.group);
    const r: Rack = { way, model, slots: [], x: RACK_X[way], y: way === "freezer" ? FZ_FLOOR : 0 };
    return r;
  });
  const rackOf = (w: PanWay) => racks.find((r) => r.way === w)!;
  function placeRack(r: Rack, x: number, y: number, z: number): void {
    r.x = x;
    r.y = y;
    r.model.group.position.set(x - RACK.w / 2, y, z + RACK.d / 2);
    r.slots.forEach((s, i) => s.group.position.set(x, y + RACK_SLIDES[RACK_SLIDES.length - 1 - i * EVERY[r.way]]!, z));
  }

  /** The pans' food: its depth's fill of the drawn pan, a pan's load (sheet pan, hotel pan, chili) and its cut. */
  const inside = new Map<number, ReturnType<typeof hotelPanInside>>();
  const insideOf = (panM: number) => {
    let v = inside.get(panM);
    if (!v) inside.set(panM, (v = hotelPanInside("1/1", panM)));
    return v;
  };
  const fullFill = (pan: typeof PAN_2) => pan.foodM / (pan.panM - insideOf(pan.panM).floorY - 0.008);
  const loads = new Map<string, THREE.Group>();
  function load(pan: typeof PAN_2, fill: number): THREE.Group {
    const key = `${pan.panM}|${fill.toFixed(3)}`;
    let g = loads.get(key);
    if (!g) {
      const s = new Shop(detail, `cool_pan_${key}`);
      sheetPan(s);
      s.at([0, PAN_DROP, 0], [0, Math.PI / 2, 0], () => hotelPan(s, "1/1", pan.panM, fill, 4));
      g = s.build(true).group;
      clip(g, rowCut);
      loads.set(key, g);
    }
    return g;
  }
  /** The cut face of a pan's food: the pan's inside across its middle, from its floor to the food's top. */
  const cutGeos = new Map<string, THREE.BufferGeometry>();
  function panCut(pan: typeof PAN_2, fill: number): THREE.BufferGeometry {
    const key = `${pan.panM}|${fill.toFixed(3)}`;
    let g = cutGeos.get(key);
    if (g) return g;
    const ins = insideOf(pan.panM);
    const top = ins.foodY(fill);
    const ys = [ins.floorY, 0.0012, 0.0024, 0.0042, 0.0065, 0.0095, 0.013].filter((y, i) => i === 0 || y < top - 2e-4);
    ys.push(top);
    const pts: THREE.Vector2[] = [];
    for (const y of ys) pts.push(new THREE.Vector2(ins.halfW(y), y));
    for (let i = ys.length - 1; i >= 0; i--) pts.push(new THREE.Vector2(-ins.halfW(ys[i]!), ys[i]!));
    g = new THREE.ShapeGeometry(new THREE.Shape(pts));
    const p = g.attributes.position as THREE.BufferAttribute;
    const uv = g.attributes.uv as THREE.BufferAttribute;
    for (let i = 0; i < p.count; i++) uv.setXY(i, 0.5, (p.getY(i) - ins.floorY) / (top - ins.floorY));
    cutGeos.set(key, g);
    return g;
  }
  // Cling film over each covered pan, beaded with what condenses under it.
  const filmM = new THREE.MeshPhysicalMaterial({ color: 0xf2f6f8, transparent: true, opacity: 0.32, roughness: 0.42, metalness: 0, clearcoat: 1, clearcoatRoughness: 0.18, depthWrite: false, side: THREE.DoubleSide, clippingPlanes: rowCut, clipShadows: true });
  const filmGeos = new Map<number, THREE.BufferGeometry>();
  const filmGeo = (panM: number) => {
    let g = filmGeos.get(panM);
    if (!g) {
      g = new THREE.PlaneGeometry(0.321, 0.526);
      g.rotateX(-Math.PI / 2);
      g.translate(0, PAN_DROP + panM + 0.0014, 0);
      filmGeos.set(panM, g);
    }
    return g;
  };

  /** The heat colours each cut shows: a pan's column of cells, bottom to top; the pot's rows, mirrored across its axis. */
  function heatTex(w: number, h: number): THREE.DataTexture {
    const t = new THREE.DataTexture(new Uint8Array(w * h * 4), w, h, THREE.RGBAFormat);
    t.colorSpace = THREE.SRGBColorSpace;
    t.magFilter = THREE.LinearFilter;
    t.minFilter = THREE.LinearFilter;
    t.wrapS = t.wrapT = THREE.ClampToEdgeWrapping;
    t.needsUpdate = true;
    return t;
  }
  // Unseen until the runs come.
  const cutMats = new Map<WayKey, THREE.MeshBasicMaterial>(WAYS.map((w) => [w.key, new THREE.MeshBasicMaterial({ toneMapped: false, visible: false })]));

  /** Each rack's pans for the batch: each filled to its depth, the last with what is left. */
  function loadRacks(litres: number, covered: boolean): void {
    for (const r of racks) {
      for (const s of r.slots) scene.remove(s.group);
      r.slots = [];
      const pan = PAN_OF[r.way];
      const n = pansFor(r.way, litres);
      for (let i = 0; i < n; i++) {
        const left = litres - i * pan.litres;
        const fill = fullFill(pan) * Math.min(1, left / pan.litres);
        const group = new THREE.Group();
        group.add(load(pan, fill).clone());
        const film = new THREE.Mesh(filmGeo(pan.panM), filmM);
        film.visible = covered;
        film.renderOrder = 3;
        group.add(film);
        const cut = new THREE.Mesh(panCut(pan, fill), cutMats.get(r.way)!);
        cut.position.set(0, PAN_DROP, 0.0004);
        group.add(cut);
        scene.add(group);
        r.slots.push({ group, film, cut });
      }
      placeRack(r, r.x, r.y, ROW_Z);
    }
  }

  /* The pot on its dunnage rack. */
  {
    const alu = palette().alu;
    const dun = new THREE.Group();
    // Five extruded slats across the top, on three inverted channels.
    for (let i = 0; i < 5; i++) {
      const g = new RoundedBoxGeometry(DUN.w, 0.024, 0.052, 2, 0.004);
      const m = new THREE.Mesh(g, alu);
      m.position.set(0, DUN.h - 0.012, -DUN.d / 2 + 0.026 + (i * (DUN.d - 0.052)) / 4);
      dun.add(m);
    }
    for (const x of [-DUN.w / 2 + 0.03, 0, DUN.w / 2 - 0.03]) {
      const top = new RoundedBoxGeometry(0.05, 0.012, DUN.d - 0.01, 2, 0.003);
      const t = new THREE.Mesh(top, alu);
      t.position.set(x, DUN.h - 0.03, 0);
      dun.add(t);
      for (const s of [-1, 1]) {
        const leg = new RoundedBoxGeometry(0.004, DUN.h - 0.036, DUN.d - 0.03, 1, 0.0015);
        const l = new THREE.Mesh(leg, alu);
        l.position.set(x + s * 0.023, (DUN.h - 0.036) / 2, 0);
        dun.add(l);
      }
    }
    dun.traverse((o) => {
      o.castShadow = true;
      o.receiveShadow = true;
    });
    dun.position.set(POT_AT.x, 0, POT_AT.z);
    scene.add(dun);
  }
  const pot = buildFixture(POT.id, POT, 2);
  const potHolder = new THREE.Group();
  pot.group.position.set(-POT.w / 2, 0, POT.d / 2);
  potHolder.add(pot.group);
  potHolder.position.set(POT_AT.x, DUN.h + 0.0005, POT_AT.z);
  scene.add(potHolder);
  pot.hooks.set?.("chili", 1);
  pot.hooks.set?.("boil", 0);
  const potLid = pot.joints.get("pot_lid_lift");
  clip(pot.group, potCut);
  /** The pot as drawn (smallwares.ts stockpot): its wall's outside, the steel's gauge, the bottom's underside over the
   *  bonded base, the bend's outside radius, and the rolled rim's height. The inside: its radius, floor and bend. */
  const POT_RO = POT.w / 2 - 0.0005;
  const POT_T = 0.0014;
  const POT_DISC = 0.0085;
  const POT_RB = 0.012;
  const POT_TOP = POT.h - 0.036 - 0.004;
  const POT_R = POT_RO - POT_T;
  const POT_FLOOR = POT_DISC + POT_T;
  const POT_BEND = POT_RB - POT_T;
  let potCutMesh: THREE.Mesh | null = null;
  function fillPot(litres: number): void {
    const { fill } = potFill(litres);
    pot.hooks.set?.("level", fill / (POT.h - 0.036 - POT_FLOOR - 0.03));
    const sh = new THREE.Shape();
    sh.moveTo(POT_R, POT_FLOOR + fill);
    sh.lineTo(POT_R, POT_FLOOR + POT_BEND);
    sh.absarc(POT_R - POT_BEND, POT_FLOOR + POT_BEND, POT_BEND, 0, -Math.PI / 2, true);
    sh.lineTo(-POT_R + POT_BEND, POT_FLOOR);
    sh.absarc(-POT_R + POT_BEND, POT_FLOOR + POT_BEND, POT_BEND, -Math.PI / 2, -Math.PI, true);
    sh.lineTo(-POT_R, POT_FLOOR + fill);
    sh.closePath();
    const g = new THREE.ShapeGeometry(sh, 6);
    const p = g.attributes.position as THREE.BufferAttribute;
    const uv = g.attributes.uv as THREE.BufferAttribute;
    for (let i = 0; i < p.count; i++) uv.setXY(i, (p.getX(i) + POT_R) / (2 * POT_R), (p.getY(i) - POT_FLOOR) / fill);
    if (potCutMesh) {
      potCutMesh.geometry.dispose();
      potCutMesh.geometry = g;
    } else {
      potCutMesh = new THREE.Mesh(g, cutMats.get("pot")!);
      potCutMesh.position.set(0, 0, 0.0004);
      potHolder.add(potCutMesh);
    }
  }
  {
    // The steel where the cut goes through it, in the cut's plane: the wall and bottom in one U, round the bends, up
    // to the rolled rim; and under it the bonded base, its chamfered edge as turned.
    const u = new THREE.Shape();
    const c = POT_RO - POT_RB;
    u.moveTo(POT_RO, POT_TOP);
    u.lineTo(POT_RO, POT_DISC + POT_RB);
    u.absarc(c, POT_DISC + POT_RB, POT_RB, 0, -Math.PI / 2, true);
    u.lineTo(-c, POT_DISC);
    u.absarc(-c, POT_DISC + POT_RB, POT_RB, -Math.PI / 2, -Math.PI, true);
    u.lineTo(-POT_RO, POT_TOP);
    u.lineTo(-POT_R, POT_TOP);
    u.lineTo(-POT_R, POT_FLOOR + POT_BEND);
    u.absarc(-c, POT_FLOOR + POT_BEND, POT_BEND, -Math.PI, -Math.PI / 2, false);
    u.lineTo(c, POT_FLOOR);
    u.absarc(c, POT_FLOOR + POT_BEND, POT_BEND, -Math.PI / 2, 0, false);
    u.lineTo(POT_R, POT_TOP);
    u.closePath();
    const rd = POT_RO - 0.008;
    const base = new THREE.Shape();
    base.moveTo(-(rd - 0.003), 0);
    base.lineTo(rd - 0.003, 0);
    base.lineTo(rd, 0.002);
    base.lineTo(rd, POT_DISC - 0.0025);
    base.lineTo(rd - 0.0035, POT_DISC);
    base.lineTo(-(rd - 0.0035), POT_DISC);
    base.lineTo(-rd, POT_DISC - 0.0025);
    base.lineTo(-rd, 0.002);
    base.closePath();
    const steel = new THREE.MeshStandardMaterial({ color: 0xd6dadd, metalness: 0.6, roughness: 0.3 });
    const alu = new THREE.MeshStandardMaterial({ color: 0xaeb3b7, metalness: 0.45, roughness: 0.5 });
    for (const [sh, m] of [[u, steel], [base, alu]] as const) {
      const g = new THREE.ShapeGeometry(sh, 8);
      g.translate(0, 0, 0.0004);
      potHolder.add(new THREE.Mesh(g, m));
    }
  }
  function coverPot(covered: boolean): void {
    if (potLid) potLid.visible = covered;
    pot.hooks.open?.(covered ? 0 : 1);
  }

  // The cook, outside, by the freezer, an eye on the boxes.
  const cook = makeStaff("p50", 2, { seed: 23, apron: 0x2b3a55, sleeves: "cuffed" });
  {
    const look = new THREE.Vector3(-1.2, 0, 0.6);
    cook.place(COOK.x, COOK.z, Math.atan2(look.x - COOK.x, look.z - COOK.z));
  }
  scene.add(cook.group);

  /* The camera: from in front, raised, both boxes and the cook in view. */
  const camera = new THREE.PerspectiveCamera(30, 16 / 9, 0.05, 40);
  const TARGET = new THREE.Vector3(0.1, 0.85, 0);
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    VIEW.set(0.06, wide ? 0.5 : 0.72, 1).normalize();
    const L = COOLER_X0;
    const R = FREEZER_X0 + FREEZER_BOX.w;
    keep.length = 0;
    // The aisles' contents and the back walls over them; the boxes' outer sides may run off the edges.
    const l = L + 0.5;
    const r = R - 0.55;
    keep.push(new THREE.Vector3(l, CUT_Y, FRONT_OUT), new THREE.Vector3(r, CUT_Y, FRONT_OUT), new THREE.Vector3(l, 2.05, -COOLER.d / 2), new THREE.Vector3(r, 2.05, -COOLER.d / 2));
    for (const x of Object.values(RACK_X)) keep.push(new THREE.Vector3(x, RACK.h + 0.1, ROW_Z));
    keep.push(new THREE.Vector3(COOK.x, 1.85, COOK.z), new THREE.Vector3(COOK.x, 0, COOK.z + 0.2));
  };
  const controls = createScrollSafeOrbit(camera, canvas, [0, 1, 0], { minDistance: 1.0, maxDistance: 14, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const hudEl = stageEl.querySelector<HTMLElement>(".bo-hud");
  /** The HUD's bottom (px from the stage's top) with the longest words it shows. */
  const hudRoom = (): number => {
    if (!hudEl) return 0;
    const parts = ["[data-cd-hud-clock]", "[data-cd-hud-doing]", "[data-cd-hud-note]"].map((k) => hudEl.querySelector<HTMLElement>(k));
    const was = parts.map((e) => e?.textContent ?? "");
    let most = 0;
    for (const d of longestDoings()) {
      if (parts[0]) parts[0].textContent = "12:00 pm, next day";
      if (parts[1]) parts[1].textContent = d;
      if (parts[2]) parts[2].textContent = "4 of 4 at 41 °F";
      most = Math.max(most, hudEl.offsetTop + hudEl.offsetHeight);
    }
    parts.forEach((e, i) => {
      if (e) e.textContent = was[i]!;
    });
    return most;
  };
  /** Room over the racks for their labels (px). */
  const LABEL_ROOM = 26;
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
    const top = Math.min(0.88, 1 - (2 * (hudRoom() + 8 + LABEL_ROOM)) / h);
    const bot = -0.95;
    let dist = 2;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 20; dist += 0.05) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      lo = left = Infinity;
      hi = right = -Infinity;
      for (const p of keep) {
        ndc.copy(p).project(camera);
        lo = Math.min(lo, ndc.y);
        hi = Math.max(hi, ndc.y);
        left = Math.min(left, ndc.x);
        right = Math.max(right, ndc.x);
      }
      if (right - left < 1.9 && hi - lo < top - bot) break;
    }
    const shiftY = (top + bot) / 2 - (hi + lo) / 2;
    const shiftX = (left + right) / 2;
    if (Math.abs(shiftY) > 0.005 || Math.abs(shiftX) > 0.005) {
      camera.setViewOffset(w, h, (shiftX * w) / 2, (shiftY * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.4;
    fog.far = dist * 3.2;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */
  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const litresIn = q<HTMLInputElement>("[data-cd-litres]");
  const litresVal = q("[data-cd-litres-val]");
  const coverBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-cd-cover]"));
  const statusEl = q("[data-cd-status]");
  const chartEl = q("[data-cd-chart]");
  const cutsEl = q<HTMLCanvasElement>("[data-cd-cuts]");
  const hudClock = q("[data-cd-hud-clock]");
  const hudDoing = q("[data-cd-hud-doing]");
  const hudNote = q("[data-cd-hud-note]");
  const tagEls = new Map(WAYS.map((w) => [w.key, q<HTMLElement>(`[data-cd-tag="${w.key}"]`)] as const));
  const tileEls = new Map(WAYS.map((w) => [w.key, { el: q(`[data-cd-tile-box="${w.key}"]`), v: q(`[data-cd-tile="${w.key}"]`), s: q(`[data-cd-tile-s="${w.key}"]`) }] as const));

  let litres = Math.max(MIN_LITRES, Math.min(MAX_LITRES, Number(litresIn?.value ?? 30)));
  let covered = false;
  /** The batch's runs once worked out (null until they are), and the batch they are for. */
  let runs: WayRun[] | null = null;
  let runsKey = "";
  const runOf = (k: WayKey) => runs?.find((r) => r.way.key === k) ?? null;
  /** The loop's clock (seconds), a minute held by the hook, the minute on show. */
  let loopT = 0;
  let frozen: number | null = null;
  let playMin = 0;

  /** The doings the HUD may show, longest first, for its room. */
  function longestDoings(): string[] {
    const out = ["Into the walk-in and the freezer at 190 °F", "Freezer first: at 41 °F, wheeled to the walk-in"];
    for (const w of WAYS) out.push(`${w.name}: 6 hours, middle at 188 °F`, `${w.name}: middle down to 135 °F, the clock starts`);
    return out;
  }

  /* The runs: off the main thread, kept by batch. */
  const cache = new Map<string, WayRun[]>();
  let worker: Worker | null = null;
  try {
    worker = new Worker(new URL("../sim/coolRuns.worker.ts", import.meta.url), { type: "module" });
    worker.onmessage = (e: MessageEvent<CoolDone>) => {
      cache.set(e.data.key, e.data.runs);
      if (e.data.key === coolKey(litres, covered)) take(e.data.runs, e.data.key);
    };
    worker.onerror = () => {
      worker?.terminate();
      worker = null;
      ask();
    };
  } catch {
    worker = null;
  }
  function ask(): void {
    const key = coolKey(litres, covered);
    const hit = cache.get(key);
    if (hit) return take(hit, key);
    if (statusEl) statusEl.textContent = "Working out the night…";
    if (worker) worker.postMessage({ litres, covered } satisfies CoolAsk);
    else {
      const r = coolAll(litres, covered);
      cache.set(key, r);
      take(r, key);
    }
  }

  /** Heat textures and the fields they show, made to each run's shape. */
  const fields = new Map<WayKey, { tex: THREE.DataTexture; f: Float32Array; w: number; h: number }>();
  function shapeFields(): void {
    if (!runs) return;
    for (const r of runs) {
      const s = r.run.shape;
      const [w, h, n] = s.kind === "slab" ? [1, s.n, s.n] : [2 * s.nr, s.nz, s.nr * s.nz];
      const had = fields.get(r.way.key);
      if (had && had.w === w && had.h === h) continue;
      had?.tex.dispose();
      const tex = heatTex(w, h);
      fields.set(r.way.key, { tex, f: new Float32Array(n), w, h });
      const m = cutMats.get(r.way.key)!;
      m.map = tex;
      m.visible = true;
      m.needsUpdate = true;
    }
  }
  function paintFields(m: number): void {
    if (!runs) return;
    for (const r of runs) {
      const fd = fields.get(r.way.key);
      if (!fd) continue;
      fieldAt(r.run, m, fd.f);
      const px = fd.tex.image.data as Uint8Array;
      const s = r.run.shape;
      if (s.kind === "slab") for (let i = 0; i < s.n; i++) heatRGB(fd.f[i]!, px, i * 4);
      else
        for (let j = 0; j < s.nz; j++)
          for (let c = 0; c < 2 * s.nr; c++) {
            const i = c < s.nr ? s.nr - 1 - c : c - s.nr;
            heatRGB(fd.f[j * s.nr + i]!, px, (j * 2 * s.nr + c) * 4);
          }
      fd.tex.needsUpdate = true;
    }
  }

  function take(r: WayRun[], key: string): void {
    runs = r;
    runsKey = key;
    shapeFields();
    loopT = 0;
    renderTiles();
    renderStatus();
    renderChart();
    lastHud = "";
    for (const t of tagEls.values()) if (t) t.dataset.last = "";
    void gate.refresh(scene);
  }

  /** The batch on the stage: the pans and the pot's fill, as soon as it changes; the runs when they come. */
  function apply(): void {
    renderControls();
    loadRacks(litres, covered);
    fillPot(litres);
    coverPot(covered);
    ask();
  }
  function renderControls(): void {
    if (litresIn) litresIn.value = String(litres);
    if (litresVal) litresVal.textContent = `${litres} litres`;
    for (const b of coverBtns) b.setAttribute("aria-pressed", String((b.dataset.cdCover === "on") === covered));
  }
  function renderTiles(): void {
    if (!runs) return;
    for (const r of runs) {
      const t = tileEls.get(r.way.key);
      if (!t) continue;
      const ok = r.ok1 && r.ok2;
      t.el?.setAttribute("data-ok", String(ok));
      if (t.v) t.v.textContent = r.total === null ? "Over 3 days" : span(r.total);
      const lines: string[] = ["to 41 °F"];
      lines.push(r.stage1 === null ? "70 °F: over 3 days" : `70 °F in ${span(r.stage1)}`);
      lines.push(r.way.key === "pot" ? `${inches(potFill(litres).fill)} deep` : `${r.pans} pan${r.pans === 1 ? "" : "s"}, ${r.way.key === "deep" ? 3 : 2} in deep`);
      lines.push(ok ? "In time" : r.ok1 ? "Too slow to 41 °F" : "Too slow to 70 °F");
      if (t.s) t.s.textContent = lines.join("\n");
    }
  }
  function renderStatus(): void {
    if (!statusEl || !runs) return;
    const okWays = runs.filter((r) => r.ok1 && r.ok2);
    const parts: string[] = [];
    const limits = "to 70 °F within 2 hours and to 41 °F within 6";
    if (okWays.length === 0) parts.push(`None of the four ways cools ${litres} litres${covered ? " under cover" : ""} ${limits}.`);
    else if (okWays.length === 1) {
      const w = okWays[0]!;
      parts.push(`Only ${PHRASE[w.way.key]} ${w.way.key === "pot" ? "cools" : "cool"} ${litres} litres ${limits}: ${span(w.stage1!)} and ${span(w.total!)}.`);
    } else parts.push(`${okWays.map((w) => PHRASE[w.way.key]).join(" and ").replace(/^./, (c) => c.toUpperCase())} cool ${litres} litres ${limits}.`);
    const p = runOf("pot")!;
    const sh = runOf("shallow")!;
    const fill = potFill(litres).fill;
    const far = Math.min(fill / 2, potFill(litres).r);
    const day = centreAt(p.run, DAY_MIN);
    parts.push(
      `In its pot, the middle of the chili is ${inches(far)} from the nearest face, against 1 in in a 2-inch pan: it takes ${span(p.run.t135!)} just to fall to 135 °F, when the clock starts, and is ${cToF(day) > cToF(COLD_C) + 0.5 ? "still " : ""}${degF(day)} a day after it went in.`,
    );
    if (!covered && sh.run.heatOut > 0) {
      const steamShare = (sh.run.steamKg * 2.39e6) / sh.run.heatOut;
      parts.push(`Uncovered, about ${Math.round(steamShare * 20) * 5}% of the heat leaving the 2-inch pans goes off as steam. Covered hot, the steam condenses on the film instead and the pans cool more slowly.`);
    } else if (covered) parts.push(`Covered, the steam condenses under the film and runs back, and every way cools more slowly than uncovered.`);
    statusEl.textContent = parts.join(" ");
    statusEl.dataset.state = okWays.length ? "ok" : "warn";
  }

  /* The chart: each way's middle from 135 °F, the Food Code's two limits shaded where it is too warm too long. */
  interface Live {
    key: WayKey;
    clip: SVGRectElement | null;
    dot: SVGCircleElement | null;
  }
  let lives: Live[] = [];
  let chartX: ((h: number) => number) | null = null;
  let chartY: ((f: number) => number) | null = null;
  let chartL = 0;
  function renderChart(): void {
    if (!chartEl || !runs) return;
    const W = Math.max(280, Math.round(chartEl.clientWidth));
    const narrow = W < 520;
    const H = narrow ? 200 : 230;
    const L = narrow ? 34 : 44;
    const R = W - (narrow ? 12 : 18);
    const T = 12;
    const B = H - 24;
    const f0 = 30;
    const f1 = 140;
    const x = (h: number) => L + (h / CHART_H) * (R - L);
    const y = (f: number) => B - ((f - f0) / (f1 - f0)) * (B - T);
    chartX = x;
    chartY = y;
    chartL = L;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg bo-svg cd-svg"><defs>`;
    for (const r of runs) svg += `<clipPath id="cd-clip-${r.way.key}"><rect x="0" y="0" width="0" height="${H}" data-cd-clip="${r.way.key}" /></clipPath>`;
    svg += `</defs>`;
    svg += `<path class="late" d="M${x(2).toFixed(1)},${y(cToF(WARM_C)).toFixed(1)}V${T}H${R}V${y(cToF(COLD_C)).toFixed(1)}H${x(6).toFixed(1)}V${y(cToF(WARM_C)).toFixed(1)}Z" />`;
    svg += `<text class="late-t" x="${(x(2) + 6).toFixed(1)}" y="${(T + 14).toFixed(1)}">Too warm, too long</text>`;
    for (const f of [41, 70, 135]) {
      svg += `<line class="grid line" x1="${L}" x2="${R}" y1="${y(f).toFixed(1)}" y2="${y(f).toFixed(1)}" />`;
      svg += `<text class="ya" x="${L - 6}" y="${(y(f) + 4).toFixed(1)}" text-anchor="end">${f}°</text>`;
    }
    for (const h of narrow ? [2, 6, 12, 18, 24] : [0, 2, 6, 12, 18, 24]) svg += `<text class="xa${h === 0 ? " st" : h === CHART_H ? " en" : ""}" x="${x(h).toFixed(1)}" y="${H - 6}">${h === 0 ? "0" : `${h} h`}</text>`;
    const pathOf = (r: WayRun) => {
      const t0 = r.run.t135!;
      const pts: string[] = [`${x(0).toFixed(1)},${y(cToF(HOT_C)).toFixed(1)}`];
      for (let i = Math.ceil(t0 / SAMPLE_MIN); i < r.run.centre.length; i++) {
        const h = (i * SAMPLE_MIN - t0) / 60;
        if (h > CHART_H) break;
        pts.push(`${x(h).toFixed(1)},${y(cToF(r.run.centre[i]!)).toFixed(1)}`);
      }
      return `M${pts.join("L")}`;
    };
    const drawn = runs.filter((r) => r.run.t135 !== null);
    for (const r of drawn) svg += `<path class="ln faint" d="${pathOf(r)}" style="--c:${WAY_COLOR[r.way.key]}" />`;
    for (const r of drawn) svg += `<g clip-path="url(#cd-clip-${r.way.key})"><path class="ln" d="${pathOf(r)}" style="--c:${WAY_COLOR[r.way.key]}" /></g>`;
    // Where each crosses 70 °F and 41 °F; a line that runs off the right, where it is then.
    for (const r of drawn)
      for (const line of ["70", "41"] as const) {
        const t = lineMin(r.run, line);
        if (t === null) continue;
        const h = (t - r.run.t135!) / 60;
        if (h <= CHART_H) svg += `<circle class="bp" cx="${x(h).toFixed(1)}" cy="${y(+line).toFixed(1)}" r="3.2" style="--c:${WAY_COLOR[r.way.key]}" />`;
      }
    for (const r of drawn) {
      const tEnd = r.run.t135! + CHART_H * 60;
      if (r.run.t41 !== null && r.run.t41 <= tEnd) continue;
      const c = centreAt(r.run, tEnd);
      svg += `<text class="bpl" x="${(R - 2).toFixed(1)}" y="${(y(cToF(c)) - 7).toFixed(1)}" text-anchor="end" style="--c:${WAY_COLOR[r.way.key]}">${degF(c)} at ${CHART_H} h</text>`;
    }
    for (const r of drawn) svg += `<circle class="now" r="4.2" cx="-20" cy="-20" data-cd-dot="${r.way.key}" style="--c:${WAY_COLOR[r.way.key]}" />`;
    svg += `</svg>`;
    chartEl.innerHTML = svg;
    lives = drawn.map((r) => ({ key: r.way.key, clip: chartEl.querySelector<SVGRectElement>(`[data-cd-clip="${r.way.key}"]`), dot: chartEl.querySelector<SVGCircleElement>(`[data-cd-dot="${r.way.key}"]`) }));
  }
  function showChart(m: number): void {
    if (!runs || !chartX || !chartY) return;
    for (const lv of lives) {
      const r = runOf(lv.key);
      if (!r || r.run.t135 === null) continue;
      const h = (m - r.run.t135) / 60;
      lv.clip?.setAttribute("width", String(h <= 0 ? 0 : Math.max(chartL, chartX(Math.min(CHART_H, h)) + 1)));
      if (lv.dot) {
        const on = h >= 0 && h <= CHART_H;
        lv.dot.setAttribute("cx", on ? chartX(h).toFixed(1) : "-20");
        lv.dot.setAttribute("cy", on ? chartY(cToF(centreAt(r.run, m))).toFixed(1) : "-20");
      }
    }
  }

  /* The four cuts side by side, to scale, with the colours' key. */
  const cutCanvases = new Map<WayKey, HTMLCanvasElement>();
  function drawCuts(m: number): void {
    if (!cutsEl || !runs) return;
    const cssW = cutsEl.clientWidth;
    const cssH = cutsEl.clientHeight;
    if (cssW < 10 || cssH < 10) return;
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    if (cutsEl.width !== Math.round(cssW * dpr) || cutsEl.height !== Math.round(cssH * dpr)) {
      cutsEl.width = Math.round(cssW * dpr);
      cutsEl.height = Math.round(cssH * dpr);
    }
    const g = cutsEl.getContext("2d");
    if (!g) return;
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, cssW, cssH);
    const narrow = cssW < 520;
    const font = narrow ? 10 : 11;
    const keyH = 26;
    const labelH = 2 * font + 10;
    // One scale for all four, fixed over the batches offered: the pot's deepest fill and the four widths fit.
    const gapPx = narrow ? 10 : 22;
    const topW = (pan: typeof PAN_2) => 2 * insideOf(pan.panM).halfW(insideOf(pan.panM).foodY(fullFill(pan)));
    const widthOf = (k: WayKey) => (k === "pot" ? 2 * potFill(1).r : topW(PAN_OF[k]));
    const sumW = WAYS.reduce((a, w) => a + widthOf(w.key), 0);
    const sc = Math.min((cssW - 3 * gapPx - 8) / sumW, (cssH - keyH - labelH - 10) / potFill(MAX_LITRES).fill);
    const base = cssH - keyH - labelH;
    let x0 = (cssW - sumW * sc - 3 * gapPx) / 2;
    for (const r of runs) {
      const fd = fields.get(r.way.key);
      if (!fd) continue;
      // The field, a pixel a cell, drawn smooth into the cut's outline.
      let cv = cutCanvases.get(r.way.key);
      if (!cv) cutCanvases.set(r.way.key, (cv = document.createElement("canvas")));
      if (cv.width !== fd.w || cv.height !== fd.h) {
        cv.width = fd.w;
        cv.height = fd.h;
      }
      const cg = cv.getContext("2d")!;
      const img = cg.createImageData(fd.w, fd.h);
      const src = fd.tex.image.data as Uint8Array;
      // The texture's rows run bottom up; the canvas's top down.
      for (let j = 0; j < fd.h; j++) img.data.set(src.subarray(j * fd.w * 4, (j + 1) * fd.w * 4), (fd.h - 1 - j) * fd.w * 4);
      cg.putImageData(img, 0, 0);
      const s = r.run.shape;
      const wM = widthOf(r.way.key);
      const hM = s.kind === "pot" ? s.fill : s.depth;
      const w = wM * sc;
      const h = hM * sc;
      const outline = () => {
        g.beginPath();
        if (s.kind === "pot") g.rect(x0, base - h, w, h);
        else {
          // The pan's taper: narrower at its floor.
          const ins = insideOf(PAN_OF[r.way.key as PanWay].panM);
          const floorW = 2 * ins.halfW(ins.floorY + 0.012) * sc;
          g.moveTo(x0 + (w - floorW) / 2, base);
          g.lineTo(x0 + (w + floorW) / 2, base);
          g.lineTo(x0 + w, base - h);
          g.lineTo(x0, base - h);
          g.closePath();
        }
      };
      g.save();
      outline();
      g.clip();
      g.imageSmoothingEnabled = true;
      g.drawImage(cv, x0, base - h, w, h);
      g.restore();
      // The steel round it: the pot's or the pan's walls and floor, open at the top.
      g.save();
      outline();
      g.lineWidth = 1;
      g.strokeStyle = "rgb(214 211 209 / 0.55)";
      g.stroke();
      g.restore();
      // Its middle.
      g.beginPath();
      g.arc(x0 + w / 2, base - h / 2, narrow ? 2.6 : 3.2, 0, Math.PI * 2);
      g.fillStyle = "#0b0b0c";
      g.fill();
      g.lineWidth = 1.4;
      g.strokeStyle = "#fff";
      g.stroke();
      // Its name and its middle's temperature.
      g.textAlign = "center";
      g.font = `600 ${font}px ui-monospace, SFMono-Regular, Menlo, monospace`;
      g.fillStyle = WAY_COLOR[r.way.key];
      g.fillText(narrow ? TINY[r.way.key] : r.way.short, x0 + w / 2, base + font + 4);
      g.fillStyle = "#e7e5e4";
      g.fillText(`${narrow ? "" : "middle "}${degF(centreAt(r.run, m))}`, x0 + w / 2, base + 2 * font + 7);
      x0 += w + gapPx;
    }
    // The key.
    const kw = Math.min(cssW - 16, 360);
    const kx = (cssW - kw) / 2;
    const ky = cssH - keyH + 4;
    const fMin = 0;
    const fMax = 200;
    const grad = g.createLinearGradient(kx, 0, kx + kw, 0);
    for (let f = fMin; f <= fMax; f += 10) grad.addColorStop((f - fMin) / (fMax - fMin), heatCss(f));
    g.fillStyle = grad;
    g.fillRect(kx, ky, kw, 7);
    g.font = `500 ${font - 1}px ui-monospace, SFMono-Regular, Menlo, monospace`;
    g.fillStyle = "#a8a29e";
    for (const f of [41, 70, 135, 190]) {
      const tx = kx + ((f - fMin) / (fMax - fMin)) * kw;
      g.fillRect(tx - 0.5, ky - 2, 1, 11);
      g.fillText(`${f} °F`, tx, ky + 19);
    }
  }

  /* The readouts. */
  let lastHud = "";
  /** What has just happened: the latest of each way's moments, within the last hour. */
  function doingAt(m: number, moving: boolean): string {
    if (!runs) return "Working out the night…";
    if (moving) return "Freezer first: at 41 °F, wheeled to the walk-in";
    if (m < 20) return "Into the walk-in and the freezer at 190 °F";
    const best = { t: -Infinity, s: "" };
    const note = (t: number | null, s: () => string) => {
      if (t === null || t > m || m - t > 60 || t <= best.t) return;
      best.t = t;
      best.s = s();
    };
    for (const r of runs) {
      const n = r.way.name;
      const t135 = r.run.t135;
      if (t135 !== null && t135 > 30) note(t135, () => `${n}: middle down to 135 °F, the clock starts`);
      note(r.run.t70, () => `${n}: middle at 70 °F`);
      note(r.run.t41, () => `${n}: middle at 41 °F`);
      if (t135 !== null) {
        note(t135 + STAGE1_MIN, () => `${n}: 2 hours, middle at ${degF(centreAt(r.run, t135 + STAGE1_MIN))}`);
        note(t135 + TOTAL_MIN, () => `${n}: 6 hours, middle at ${degF(centreAt(r.run, t135 + TOTAL_MIN))}`);
      }
    }
    if (best.s) return best.s;
    // Else the slowest's middle.
    const slow = runs.reduce((a, b) => (centreAt(b.run, m) > centreAt(a.run, m) ? b : a));
    return `${slow.way.name}: middle at ${degF(centreAt(slow.run, m))}`;
  }
  function showHud(m: number, moving: boolean): void {
    const cold = runs ? runs.filter((r) => r.run.t41 !== null && r.run.t41 <= m + 1e-9).length : 0;
    const doing = doingAt(m, moving);
    const note = cold ? `${cold} of 4 at 41 °F` : "Cooling";
    const clock = clockAt(m);
    const sig = `${clock}|${doing}|${note}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = START_TOD + m >= 1440 ? `${clock}, next day` : clock;
    if (hudDoing) hudDoing.textContent = doing;
    if (hudNote) {
      hudNote.textContent = note;
      hudNote.dataset.kind = cold ? "cold" : "hot";
    }
  }
  /** The labels over each way: its middle's temperature. One that would cover another is lifted clear on a stem. */
  const tagAt = new THREE.Vector3();
  const tagSize = new Map<HTMLElement, { w: number; h: number }>();
  const tagStem = new Map<HTMLElement, string>();
  const tagItems: { el: HTMLElement; x: number; y: number; w: number; h: number }[] = [];
  const tagPlaced: { l: number; r: number; t: number; b: number }[] = [];
  function anchorOf(k: WayKey): THREE.Vector3 {
    if (k === "pot") return tagAt.set(POT_AT.x, DUN.h + POT.h + 0.06, POT_AT.z);
    const r = rackOf(k);
    return tagAt.set(r.x, r.y + RACK.h + 0.06, ROW_Z);
  }
  function showTags(m: number): void {
    const W = stageEl!.clientWidth;
    const H = stageEl!.clientHeight;
    tagItems.length = 0;
    for (const w of WAYS) {
      const el = tagEls.get(w.key);
      if (!el) continue;
      const r = runOf(w.key);
      anchorOf(w.key).project(camera);
      const on = !!r && tagAt.z < 1 && Math.abs(tagAt.x) < 0.98 && Math.abs(tagAt.y) < 0.98;
      el.hidden = !on;
      if (!on || !r) continue;
      const c = centreAt(r.run, m);
      const text = degF(c);
      if (el.dataset.last !== text || !tagSize.has(el)) {
        el.dataset.last = text;
        const v = el.querySelector("b");
        if (v) v.textContent = text;
        el.toggleAttribute("data-boiled", c <= COLD_C + 1e-6);
        tagSize.set(el, { w: el.offsetWidth, h: el.offsetHeight });
      }
      const sz = tagSize.get(el)!;
      tagItems.push({ el, x: ((tagAt.x + 1) / 2) * W, y: ((1 - tagAt.y) / 2) * H, w: sz.w, h: sz.h });
    }
    tagItems.sort((a, b) => b.y - a.y);
    tagPlaced.length = 0;
    for (const it of tagItems) {
      const l = it.x - it.w / 2 - 3;
      const r = it.x + it.w / 2 + 3;
      let b = it.y - 6;
      for (let moved = true, n = 0; moved && n < 8; n++) {
        moved = false;
        for (const p of tagPlaced)
          if (l < p.r && r > p.l && b > p.t - 3 && b - it.h < p.b + 3) {
            b = p.t - 3;
            moved = true;
          }
      }
      tagPlaced.push({ l, r, t: b - it.h, b });
      const x = Math.min(W - it.w / 2 - 4, Math.max(it.w / 2 + 4, it.x));
      const stem = `${Math.round(it.y - b)}|${Math.round(it.x - x)}`;
      it.el.style.transform = `translate(${x.toFixed(1)}px, ${b.toFixed(1)}px)`;
      if (tagStem.get(it.el) !== stem) {
        tagStem.set(it.el, stem);
        it.el.style.setProperty("--stem", `${Math.round(it.y - b)}px`);
        it.el.style.setProperty("--sx", `${Math.round(it.x - x)}px`);
      }
    }
  }

  /** The pot's steam, as the run has it off its open top. The stockpot draws 1 for a hard boil's 1500 W in a warm
   *  kitchen; in the walk-in's cold, damp air the vapour fogs at once and hangs, so a couple of hundred watts off hot
   *  chili draw as a hard boil's: thick while the top is hot, thinning as it cools. */
  let potSteam = -1;
  function showPot(m: number, t: number): void {
    const p = runOf("pot");
    // The first sample is before any step (no steam yet): from the second on.
    const i = p ? Math.min(p.run.steam.length - 1, Math.max(1, Math.round(m / SAMPLE_MIN))) : 0;
    const watts = p && !covered ? p.run.steam[i]! * Math.PI * potFill(litres).r ** 2 : 0;
    const steam = watts < 20 ? 0 : Math.min(3, watts / 120);
    if (Math.abs(steam - potSteam) > 0.01 || (steam === 0) !== (potSteam === 0)) {
      potSteam = steam;
      pot.hooks.set?.("steam", steam);
    }
    pot.hooks.tick?.(t, 0, 0);
  }

  /** Where the freezer pans' rack is at the loop's minute: in the freezer, or once moved, in the walk-in; and the two
   *  doors, open as it goes through them (`k` 0..1 through the move). */
  function showMove(moved: boolean, k: number | null): void {
    const r = rackOf("freezer");
    const inCooler = k !== null ? k >= 0.5 : moved;
    const x = inCooler ? RACK_X.moved : RACK_X.freezer;
    const y = inCooler ? 0 : FZ_FLOOR;
    if (r.x !== x || r.y !== y) placeRack(r, x, y, ROW_Z);
    const fz = k === null ? 0 : smooth(k / 0.15) * (1 - smooth((k - 0.35) / 0.15));
    const co = k === null ? 0 : smooth((k - 0.4) / 0.15) * (1 - smooth((k - 0.8) / 0.15));
    freezer.hooks.open?.(fz);
    cooler.hooks.open?.(co);
  }

  apply();
  frame();
  new ResizeObserver(() => {
    if (stageEl.clientWidth === fitW && stageEl.clientHeight === fitH) return;
    tagSize.clear();
    frame();
  }).observe(stageEl);
  void document.fonts?.ready.then(() => frame());
  if (chartEl) {
    let lastW = chartEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(chartEl.clientWidth - lastW) < 4) return;
      lastW = chartEl.clientWidth;
      renderChart();
    }).observe(chartEl);
  }
  let pending = 0;
  litresIn?.addEventListener("input", () => {
    litres = Math.max(MIN_LITRES, Math.min(MAX_LITRES, Math.round(Number(litresIn.value))));
    renderControls();
    window.clearTimeout(pending);
    pending = window.setTimeout(apply, 160);
  });
  for (const b of coverBtns)
    b.addEventListener("click", () => {
      const c = b.dataset.cdCover === "on";
      if (c === covered) return;
      covered = c;
      apply();
    });
  void cook.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  /** The loop's minute, and how far through the move it is (null when not moving). */
  function minuteAt(s: number): { m: number; move: number | null } {
    const fz = runOf("freezer");
    const mv = fz?.run.movedMin ?? null;
    if (s < IN_S) return { m: 0, move: null };
    const a = s - IN_S;
    if (mv === null || mv >= DAY_MIN) return { m: Math.min(DAY_MIN, (a / DAY_S) * DAY_MIN), move: null };
    const before = (mv / DAY_MIN) * DAY_S;
    if (a < before) return { m: (a / DAY_S) * DAY_MIN, move: null };
    if (a < before + MOVE_S) return { m: mv, move: (a - before) / MOVE_S };
    return { m: Math.min(DAY_MIN, mv + ((a - before - MOVE_S) / DAY_S) * DAY_MIN), move: null };
  }
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    const loopS = IN_S + DAY_S + MOVE_S + END_S;
    loopT = (loopT + dt) % loopS;
    const t = now / 1000;
    const at = frozen !== null ? { m: Math.min(DAY_MIN, Math.max(0, frozen)), move: null } : minuteAt(loopT);
    playMin = at.m;
    const mv = runOf("freezer")?.run.movedMin ?? null;
    showMove(mv !== null && playMin >= mv, at.move);
    paintFields(playMin);
    showPot(playMin, t);
    for (const b of [cooler, freezer]) b.hooks.tick?.(t, dt, 1);
    cook.stand(t);
    showHud(playMin, at.move !== null);
    showChart(playMin);
    drawCuts(playMin);
    controls.update();
    fitClip();
    renderer.render(scene, camera);
    showTags(playMin);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    cool: {
      info: () => ({
        mountMs: Math.round(mountMs),
        litres,
        covered,
        potSteam: +potSteam.toFixed(3),
        ready: runs !== null && runsKey === coolKey(litres, covered),
        t: +playMin.toFixed(1),
        runs: (runs ?? []).map((r) => ({
          key: r.way.key,
          pans: r.pans,
          t135: r.run.t135 === null ? null : +r.run.t135.toFixed(1),
          stage1: r.stage1 === null ? null : +r.stage1.toFixed(1),
          total: r.total === null ? null : +r.total.toFixed(1),
          ok1: r.ok1,
          ok2: r.ok2,
          movedMin: r.run.movedMin,
          dayF: +cToF(centreAt(r.run, DAY_MIN)).toFixed(1),
        })),
        tags: WAYS.map((w) => tagEls.get(w.key)?.textContent ?? null),
        status: statusEl?.textContent ?? "",
      }),
      hold(min: number | null) {
        frozen = min;
      },
      /** Look from `eye` at `at` (world), or back to the framed view. */
      view(eye: [number, number, number] | null, at?: [number, number, number]) {
        userMoved = !!eye;
        if (eye && at) {
          camera.clearViewOffset();
          camera.position.set(...eye);
          controls.target.set(...at);
          controls.update();
        } else frame();
      },
      set(s: { litres?: number; covered?: boolean }) {
        if (s.litres !== undefined) litres = Math.max(MIN_LITRES, Math.min(MAX_LITRES, Math.round(s.litres)));
        if (s.covered !== undefined) covered = s.covered;
        apply();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
