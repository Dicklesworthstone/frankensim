/**
 * The week (chapter 06): a commissary's packing line, its cooking line behind it, and the week laid out beside them
 * as six speed racks, one per production day.
 *
 * Everything shown follows `runCommissary` (sim/commissaryWeek.ts) for the settings on the page:
 *   - a rack is one day's paid shift, ten tiers of two half pans; one half pan stands for 72 crew-minutes, so a full
 *     rack is exactly the crew's 24 hours. As the day's clock runs the rack fills with the work done: open pans for
 *     the day's own meals, film-covered pans for what is cooked ahead and chilled. Work past the shift (overtime)
 *     stacks on the rack's top under a red glow; empty tiers at the end of the day are crew time nobody needed.
 *   - the cooks at the packing table work while the day's work lasts, then step back and stand: a kinematic
 *     illustration (people.ts Routines, pure functions of the week's clock), timed by the computed hours.
 *   - the kettle, tilting skillet, combi and blast chiller behind them run by the minutes each station cooks that day.
 */
import * as THREE from "three";
import { DAYS, KITCHEN, MENU, SWAPS, runCommissary, type CommissaryRun, type DayResult, type Settings, type Strategy } from "../sim/commissaryWeek";
import { bindCanvasResize, createScrollSafeOrbit, createStudioRenderer, getSoftParticleTexture, observeInView, shaderGate } from "./kit";
import { buildFixture, type Model } from "./models/all";
import { speedRackSlides } from "./models/cold";
import { palette } from "./models/materials";
import { Routine, TOOLS, makeStaff, type Cohort, type RoutineStep, type Staff } from "./models/people";
import { studioEnvironment, studioLights } from "./models/stage";

/* ------------------------------------------------------------------------------------------------------------ */
/* The room.                                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

const SHIFT_H = KITCHEN.shiftMin / 60;
/** Animation seconds per hour of the working day, and between one day's end and the next day's start. */
const SEC_PER_H = 0.8;
const NIGHT_S = 1.2;
/** Half pans per rack, and the crew minutes one stands for. */
const RACK_PANS = 20;
const PAN_MIN = (KITCHEN.cooks * KITCHEN.shiftMin) / RACK_PANS;

const TABLE = { id: "NSF-PREPTABLE-96", w: 2.438, d: 0.762, h: 0.914 };
/** The packing table's cook-side edge (it faces -z, the cooks stand behind it facing +z) and its left end. */
const TABLE_Z = 0.3;
const TABLE_X0 = -2.8;
const COOK_Z = 0.05;
const BACK_Z = -0.12;
const COOKS: { cohort: Cohort; x: number; seed: number; apron: number; share: number; tempo: number }[] = [
  { cohort: "p50", x: -2.32, seed: 41, apron: 0x26344d, share: 1.08, tempo: 0.95 },
  { cohort: "p95", x: -1.58, seed: 57, apron: 0x3b3d43, share: 1.0, tempo: 1.05 },
  { cohort: "p05", x: -0.84, seed: 63, apron: 0x45563f, share: 0.92, tempo: 1.0 },
];

/** The cooking line along the back, fronts facing +z. */
const LINE_Z = -1.45;
const LINE = [
  { id: "NSF-KETTLE-40GAL", w: 0.914, d: 0.965, h: 1.067, x0: -3.3 },
  { id: "NSF-TILT-SKILLET-40GAL", w: 1.219, d: 0.991, h: 1.067, x0: -2.28 },
  { id: "NSF-COMBI-20PAN-ROLLIN", w: 0.879, d: 0.842, h: 1.872, x0: -0.95 },
  { id: "NSF-BLAST-CHILLER-15PAN", w: 0.84, d: 0.88, h: 1.85, x0: 0.02 },
] as const;
const HOLD = { id: "NSF-HOLDING-CABINET", w: 0.66, d: 0.813, h: 1.778, x0: -3.62 };
const CARRIER = { id: "NSF-CARRIER-INSULATED", w: 0.46, d: 0.635, h: 0.625 };
const CARRIER_X0 = -0.26;
const CARRIER_Z = 1.0;

/** The week: six racks, fronts facing +z. */
const RACK = { id: "NSF-SPEEDRACK-20PAN", w: 0.521, d: 0.66, h: 1.753 };
const RACK_X0 = 1.32;
const RACK_STEP = 0.66;
const RACK_Z = 0.62;
const HALF = { id: "NSF-HOTEL-PAN-HALF-4IN", w: 0.325, d: 0.265, h: 0.1 };
const SLIDES = speedRackSlides(RACK.h);
const TRAY = { w: 0.457, d: 0.66, h: 0.025 };

const DAY_NAMES = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"] as const;

/* ------------------------------------------------------------------------------------------------------------ */
/* Helpers.                                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

/** Place a fixture with its front face at z = `z`, its left end (as seen from the front) at x0, facing +z or -z. */
function place(m: Model, x0: number, z: number, facePlusZ = true): THREE.Group {
  const g = m.group;
  if (facePlusZ) g.position.set(x0, 0, z);
  else {
    g.rotation.y = Math.PI;
    g.position.set(x0 + m.ext.w, 0, z);
  }
  return g;
}

function tilePlane(w: number, h: number): THREE.Mesh {
  const g = new THREE.PlaneGeometry(w, h);
  const uv = g.attributes.uv as THREE.BufferAttribute;
  for (let i = 0; i < uv.count; i++) uv.setXY(i, (uv.getX(i) * w) / 0.608, (uv.getY(i) * h) / 0.608);
  const m = new THREE.Mesh(g, palette().tile);
  m.receiveShadow = true;
  return m;
}

function backdrop(): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = 4;
  c.height = 256;
  const x = c.getContext("2d")!;
  const g = x.createLinearGradient(0, 0, 0, 256);
  g.addColorStop(0, "#0a0807");
  g.addColorStop(0.5, "#151510");
  g.addColorStop(0.78, "#1d1a13");
  g.addColorStop(1, "#100e0a");
  x.fillStyle = g;
  x.fillRect(0, 0, 4, 256);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

/** A plating spoon in the grip frame (handle along +z to the bowl), as the plate pose holds it. */
function makeSpoon(): THREE.Object3D {
  const p = palette();
  const K = TOOLS.spoon;
  const len = K.reach + K.handleBack - 0.03;
  const handle = new THREE.Mesh(new THREE.CylinderGeometry(0.0045, 0.006, len, 8).rotateX(Math.PI / 2).translate(0, 0, (K.reach - 0.03 - K.handleBack) / 2), p.ssBright);
  const bowl = new THREE.Mesh(new THREE.SphereGeometry(0.024, 14, 8, 0, Math.PI * 2, Math.PI / 2, Math.PI / 2), p.ssBright);
  bowl.scale.set(1, 1.5, 0.45);
  bowl.rotation.x = Math.PI / 2;
  bowl.position.set(0, 0, K.reach);
  const g = new THREE.Group();
  g.add(handle, bowl);
  g.traverse((o) => {
    if ((o as THREE.Mesh).isMesh) o.castShadow = true;
  });
  return g;
}

/** A kraft meal tray with a portion in it, its base centre at the origin. */
function makeTray(food: THREE.Material): THREE.Group {
  const kraft = new THREE.MeshStandardMaterial({ color: 0xb08a5a, roughness: 0.85 });
  const g = new THREE.Group();
  const w = 0.2;
  const d = 0.15;
  const h = 0.045;
  const t = 0.003;
  const parts: [number, number, number, number, number, number][] = [
    [w, t, d, 0, t / 2, 0],
    [w, h, t, 0, h / 2, d / 2 - t / 2],
    [w, h, t, 0, h / 2, -d / 2 + t / 2],
    [t, h, d, w / 2 - t / 2, h / 2, 0],
    [t, h, d, -w / 2 + t / 2, h / 2, 0],
  ];
  for (const [a, b, c, x, y, z] of parts) {
    const m = new THREE.Mesh(new THREE.BoxGeometry(a, b, c), kraft);
    m.position.set(x, y, z);
    m.castShadow = m.receiveShadow = true;
    g.add(m);
  }
  const mound = new THREE.Mesh(new THREE.SphereGeometry(0.06, 16, 8, 0, Math.PI * 2, 0, Math.PI / 2), food);
  mound.scale.set(1.2, 0.35, 0.85);
  mound.position.y = t;
  g.add(mound);
  return g;
}

/** A sheet-pan tray as the racks carry the half pans: an aluminium floor and a rolled rim, its rim top at y = 0. */
function trayTemplate(): THREE.Group {
  const alu = palette().alu;
  const g = new THREE.Group();
  const t = 0.0016;
  const floor = new THREE.Mesh(new THREE.BoxGeometry(TRAY.w - 0.01, t, TRAY.d - 0.01), alu);
  floor.position.y = -TRAY.h + t / 2;
  g.add(floor);
  const rims: [number, number, number, number][] = [
    [TRAY.w, 0.006, 0, TRAY.d / 2 - 0.003],
    [TRAY.w, 0.006, 0, -TRAY.d / 2 + 0.003],
    [0.006, TRAY.d, TRAY.w / 2 - 0.003, 0],
    [0.006, TRAY.d, -TRAY.w / 2 + 0.003, 0],
  ];
  for (const [a, c, x, z] of rims) {
    const m = new THREE.Mesh(new THREE.BoxGeometry(a, TRAY.h, c), alu);
    m.position.set(x, -TRAY.h / 2, z);
    g.add(m);
  }
  return g;
}

/** Many copies of a template drawn as one instanced mesh per part. */
class Instances {
  readonly meshes: THREE.InstancedMesh[] = [];
  private rel: THREE.Matrix4[] = [];
  private n = 0;
  private readonly m = new THREE.Matrix4();
  constructor(template: THREE.Object3D, readonly max: number, parent: THREE.Object3D, tint?: (m: THREE.Material) => THREE.Material) {
    template.updateMatrixWorld(true);
    const inv = new THREE.Matrix4().copy(template.matrixWorld).invert();
    template.traverseVisible((o) => {
      const mesh = o as THREE.Mesh;
      if (!mesh.isMesh) return;
      const mat = Array.isArray(mesh.material) ? mesh.material.map((x) => (tint ? tint(x) : x)) : tint ? tint(mesh.material) : mesh.material;
      const im = new THREE.InstancedMesh(mesh.geometry, mat, max);
      im.castShadow = true;
      im.receiveShadow = true;
      im.count = 0;
      im.frustumCulled = false;
      this.meshes.push(im);
      this.rel.push(new THREE.Matrix4().multiplyMatrices(inv, mesh.matrixWorld));
      parent.add(im);
    });
  }
  begin(): void {
    this.n = 0;
  }
  push(world: THREE.Matrix4): void {
    if (this.n >= this.max) return;
    this.meshes.forEach((im, i) => im.setMatrixAt(this.n, this.m.multiplyMatrices(world, this.rel[i]!)));
    this.n++;
  }
  end(): void {
    for (const im of this.meshes) {
      im.count = this.n;
      im.instanceMatrix.needsUpdate = true;
    }
  }
}

/** A rack's card: the day's name and its meals, white on a dark card, for a plane 0.46 by 0.17 m. */
function dayLabel(day: string, meals: number): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = 512;
  c.height = 190;
  const x = c.getContext("2d")!;
  x.fillStyle = "#15110e";
  x.beginPath();
  x.roundRect(4, 4, 504, 182, 22);
  x.fill();
  x.strokeStyle = "rgba(255,228,204,0.35)";
  x.lineWidth = 4;
  x.stroke();
  x.fillStyle = "#f4eee8";
  x.textAlign = "center";
  x.font = "700 84px Inter Variable, Inter, system-ui, sans-serif";
  x.fillText(day.toUpperCase(), 256, 98);
  x.font = "500 50px JetBrains Mono Variable, JetBrains Mono, monospace";
  x.fillStyle = "#c9bcae";
  x.fillText(`${meals.toLocaleString("en-US")} meals`, 256, 162);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = 4;
  return t;
}

const clamp01 = (v: number) => Math.max(0, Math.min(1, v));
const fmtH = (min: number) => {
  const h = min / 60;
  return h >= 9.95 ? `${Math.round(h)} h` : `${h.toFixed(1)} h`;
};
const clock = (h: number) => {
  const m = Math.round((6 + h) * 60);
  return `${String(Math.floor(m / 60) % 24).padStart(2, "0")}:${String(m % 60).padStart(2, "0")}`;
};

/* ------------------------------------------------------------------------------------------------------------ */
/* The week's timeline: when each day runs, how long each cook works, what each rack holds.                       */
/* ------------------------------------------------------------------------------------------------------------ */

interface DaySpan {
  /** Animation seconds: the day's start, and its working length (shift and any overtime). */
  t0: number;
  len: number;
  /** Hours each cook works. */
  busyH: number[];
  /** Half pans on its rack at the end of the day: for the day itself, cooked ahead, and past the shift. */
  dayOf: number;
  ahead: number;
  over: number;
  r: DayResult;
}

function timeline(run: CommissaryRun): { days: DaySpan[]; total: number } {
  const days: DaySpan[] = [];
  let t = 0;
  for (const r of run.plan.days) {
    const worked = r.dayOfMin + r.aheadMin;
    const hours = SHIFT_H + r.otMin / KITCHEN.cooks / 60;
    const each = worked / KITCHEN.cooks / 60;
    // Shared a little unevenly (nobody stops at the same minute), never past the day's hours: what one cook cannot
    // take, the others do.
    const busyH = COOKS.map((c) => each * (r.otMin > 0 ? 1 : c.share));
    for (let pass = 0; pass < 2; pass++) {
      const spill = busyH.reduce((a, b) => a + Math.max(0, b - hours), 0);
      const room = busyH.filter((b) => b < hours).length;
      for (let i = 0; i < busyH.length; i++) busyH[i] = busyH[i]! >= hours ? hours : Math.min(hours, busyH[i]! + (room ? spill / room : 0));
    }
    const total = Math.round(worked / PAN_MIN);
    const dayOf = Math.min(total, Math.round(r.dayOfMin / PAN_MIN));
    const over = Math.max(0, total - RACK_PANS);
    days.push({ t0: t, len: hours * SEC_PER_H, busyH, dayOf: dayOf - over, ahead: total - dayOf, over, r });
    t += hours * SEC_PER_H + NIGHT_S;
  }
  return { days, total: t };
}

/** The cooks' week: plating while the day's work lasts, a step back from the table and standing after. */
function cookSteps(spans: DaySpan[], i: number, at: THREE.Vector3, shuffle?: { back: number; fwd: number }): RoutineStep[] {
  const c = COOKS[i]!;
  const steps: RoutineStep[] = [];
  const back: [number, number] = [c.x, BACK_Z];
  const front: [number, number] = [c.x, COOK_Z];
  const sb = shuffle?.back ?? 1;
  const sf = shuffle?.fwd ?? 1;
  spans.forEach((s) => {
    const work = s.busyH[i]! * SEC_PER_H;
    const free = s.len + NIGHT_S - work;
    const busy = work >= 1.2;
    if (busy) steps.push({ kind: "work", work: "plate", duration: work, opts: { at, surfaceY: TABLE.h, tempo: c.tempo } });
    const rest = busy ? free : s.len + NIGHT_S;
    if (rest > sb + sf + 1.2) {
      steps.push({ kind: "shuffle", to: back });
      steps.push({ kind: "stand", duration: rest - sb - sf });
      steps.push({ kind: "shuffle", to: front });
    } else steps.push({ kind: "stand", duration: rest });
  });
  return steps;
}

/* ------------------------------------------------------------------------------------------------------------ */

export function mountWeekLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-week-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-week-canvas]")!;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 9, 17);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.6;
  const lights = studioLights(scene, new THREE.Vector3(0.6, 0.8, -0.3), 5.2, 2048);
  lights.key.position.set(-2.2, 6.4, 4.6);
  lights.key.intensity = 2.0;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(4.5, 3.2, -3.4);
  lights.rim.intensity = 0.8;
  lights.hemi.intensity = 0.5;

  const floor = tilePlane(12, 7);
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(0.8, 0, -0.4);
  scene.add(floor);

  /* The packing table, the holding cabinet and the carriers. */
  const table = buildFixture(TABLE.id, TABLE, 2);
  scene.add(place(table, TABLE_X0, TABLE_Z, false));
  const hold = buildFixture(HOLD.id, HOLD, 2);
  scene.add(place(hold, HOLD.x0, TABLE_Z + TABLE.d, true));
  const carrierTpl = buildFixture(CARRIER.id, CARRIER, 1);
  const carriers: THREE.Object3D[] = [];
  for (let k = 0; k < 6; k++) {
    const g = k === 0 ? carrierTpl.group : carrierTpl.group.clone();
    const col = Math.floor(k / 2);
    g.position.set(CARRIER_X0 + col * (CARRIER.w + 0.03), (k % 2) * CARRIER.h, CARRIER_Z);
    g.visible = false;
    scene.add(g);
    carriers.push(g);
  }

  /* The cooking line. */
  const line = LINE.map((f) => {
    const m = buildFixture(f.id, f, 2);
    scene.add(place(m, f.x0, LINE_Z, true));
    return m;
  });
  const [kettle, skillet, combi, chiller] = line as [Model, Model, Model, Model];

  /* The racks, their labels, and the light on today's. */
  const forecast = KITCHEN.forecast;
  const rackX = (d: number) => RACK_X0 + d * RACK_STEP;
  const labelMats: THREE.MeshBasicMaterial[] = [];
  const labelGeo = new THREE.PlaneGeometry(0.5, 0.186);
  for (let d = 0; d < DAYS.length; d++) {
    const m = buildFixture(RACK.id, RACK, 1);
    scene.add(place(m, rackX(d), RACK_Z, true));
    // Hung on the rack's bottom front rail, under the lowest pan.
    const mat = new THREE.MeshBasicMaterial({ map: dayLabel(DAYS[d]!, forecast[d]!), transparent: true, opacity: 0.62, toneMapped: false });
    labelMats.push(mat);
    const lab = new THREE.Mesh(labelGeo, mat);
    lab.position.set(rackX(d) + RACK.w / 2, 0.165, RACK_Z + 0.014);
    scene.add(lab);
  }
  const glowTex = getSoftParticleTexture();
  const nowGlow = new THREE.Mesh(new THREE.PlaneGeometry(1.0, 1.3), new THREE.MeshBasicMaterial({ map: glowTex, color: 0xf2b134, transparent: true, opacity: 0.55, blending: THREE.AdditiveBlending, depthWrite: false, toneMapped: false }));
  nowGlow.rotation.x = -Math.PI / 2;
  nowGlow.position.y = 0.002;
  scene.add(nowGlow);
  const otGlow = new THREE.Mesh(new THREE.PlaneGeometry(1.1, 1.1), new THREE.MeshBasicMaterial({ map: glowTex, color: 0xfb7185, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false, toneMapped: false }));
  otGlow.position.y = RACK.h + 0.25;
  scene.add(otGlow);

  /* The pans: one template per food, one for chilled (under film), the overflow tinted red; and the trays. */
  const panGroup = new THREE.Group();
  scene.add(panGroup);
  const panTemplate = (food: number, fill: number) => {
    const m = buildFixture(HALF.id, HALF, 1);
    m.hooks.set?.("food", food);
    m.hooks.set?.("fill", fill);
    // To the pan's base centre.
    m.group.position.set(-HALF.w / 2, 0, HALF.d / 2);
    const g = new THREE.Group();
    g.add(m.group);
    return g;
  };
  const MAXP = DAYS.length * RACK_PANS + 24;
  const open = [0, 1, 2, 3].map((f) => new Instances(panTemplate(f, 0.62), MAXP, panGroup));
  // Cling film over a chilled pan, frosted where the cold has misted it: it reads as the chart's ice blue.
  const film = new THREE.MeshPhysicalMaterial({ color: 0x9adcf6, emissive: 0x1d6f93, emissiveIntensity: 0.55, roughness: 0.28, metalness: 0, transparent: true, opacity: 0.82, clearcoat: 1, clearcoatRoughness: 0.15, depthWrite: false });
  const chilledTpl = panTemplate(0, 0.55);
  // Wrapped over the top and down the sides, so a chilled pan shows as one from the front as well as from above.
  const wrap = new THREE.Mesh(new THREE.BoxGeometry(HALF.w + 0.008, HALF.h * 0.94, HALF.d + 0.008), film);
  wrap.position.y = HALF.h * 0.53;
  chilledTpl.add(wrap);
  const chilled = new Instances(chilledTpl, MAXP, panGroup);
  const red = (m: THREE.Material) => {
    const c = (m as THREE.MeshStandardMaterial).clone();
    if ("emissive" in c) {
      c.emissive = new THREE.Color(0xfb3b4a);
      c.emissiveIntensity = 0.32;
    }
    return c;
  };
  const overflow = new Instances(panTemplate(0, 0.62), 40, panGroup, red);
  const trays = new Instances(trayTemplate(), DAYS.length * 10 + 10, panGroup);

  /* The cooks, each with a tray and a pan of food on the table in front, and a spoon. */
  const foodMat = new THREE.MeshStandardMaterial({ color: 0x8a5230, roughness: 0.7 });
  const cooks = COOKS.map((c, i) => {
    const st: Staff = makeStaff(c.cohort, 2, { seed: c.seed, apron: c.apron, sleeves: i === 1 ? "rolled" : "cuffed" });
    st.place(c.x, COOK_Z, 0);
    scene.add(st.group);
    const at = new THREE.Vector3(c.x, TABLE.h + 0.004, TABLE_Z + 0.15);
    const tray = makeTray(foodMat);
    tray.position.copy(at);
    scene.add(tray);
    const src = buildFixture(HALF.id, HALF, 1);
    src.hooks.set?.("food", [0, 1, 2][i]!);
    src.hooks.set?.("fill", 0.5);
    src.group.position.set(c.x + 0.26 - HALF.w / 2, TABLE.h, TABLE_Z + 0.5 + HALF.d / 2);
    scene.add(src.group);
    const spoon = makeSpoon();
    st.hand("R").add(spoon);
    spoon.visible = false;
    return { c, st, at, spoon, routine: null as Routine | null };
  });

  /* Camera: from the front, raised, holding the line, the table and all six racks. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 10, 0.05, 60);
  // A wide stage looks at the line and the racks from the front; a narrow one from beyond the racks' end, so that
  // the row recedes into the picture instead of shrinking to fit across it.
  const TARGET = new THREE.Vector3();
  const VIEW = new THREE.Vector3();
  /** The boxes the view must hold: on a wide stage the whole room; on a narrow one the table and the racks. */
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    keep.length = 0;
    const wide = stageEl.clientWidth >= 1.5 * stageEl.clientHeight;
    if (wide) {
      TARGET.set(0.7, 0.92, 0.1);
      VIEW.set(0.06, 0.3, 1).normalize();
      for (const x of [HOLD.x0, rackX(DAYS.length - 1) + RACK.w]) for (const y of [0.05, RACK.h + 0.2]) for (const z of [LINE_Z, RACK_Z]) keep.push(new THREE.Vector3(x, y, z));
    } else {
      // From in front of the crew's end of the room: the cooks near, the week running away to the right.
      TARGET.set(0.9, 0.95, 0.2);
      VIEW.set(-0.62, 0.36, 1).normalize();
      for (const x of [TABLE_X0 + 0.3, rackX(DAYS.length - 1) + RACK.w]) for (const y of [0.1, RACK.h + 0.15]) for (const z of [RACK_Z - RACK.d, TABLE_Z + TABLE.d]) keep.push(new THREE.Vector3(x, y, z));
    }
  };
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 1.2, maxDistance: 18, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const ndc = new THREE.Vector3();
  const frame = () => {
    if (userMoved) return;
    aim();
    // Fit for the stage's shape now, even while the canvas's own resize waits on the shaders.
    camera.aspect = Math.max(0.2, stageEl.clientWidth / Math.max(1, stageEl.clientHeight));
    camera.updateProjectionMatrix();
    let dist = 3;
    for (; dist < 26; dist += 0.08) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      if (keep.every((p) => (ndc.copy(p).project(camera), Math.abs(ndc.x) < 0.95 && Math.abs(ndc.y) < 0.86))) break;
    }
    // The haze begins past the far racks and closes in well beyond them.
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 0.95;
    fog.far = dist * 2.1;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  frame();
  new ResizeObserver(() => frame()).observe(stageEl);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page: plan, menu, racks, forecast error; the week chart; the tiles; the status line.                    */
  /* ---------------------------------------------------------------------------------------------------------- */

  const planBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-wk-plan]"));
  const menuEl = root.querySelector<HTMLElement>("[data-wk-menu]");
  const racksIn = root.querySelector<HTMLInputElement>("[data-wk-racks]");
  const racksVal = root.querySelector<HTMLElement>("[data-wk-racks-val]");
  const errIn = root.querySelector<HTMLInputElement>("[data-wk-err]");
  const errVal = root.querySelector<HTMLElement>("[data-wk-err-val]");
  const statusEl = root.querySelector<HTMLElement>("[data-wk-status]");
  const chartEl = root.querySelector<HTMLElement>("[data-week-chart]");
  const hudDay = root.querySelector<HTMLElement>("[data-wk-hud-day]");
  const hudClock = root.querySelector<HTMLElement>("[data-wk-hud-clock]");
  const hudNote = root.querySelector<HTMLElement>("[data-wk-hud-note]");
  const tile = (k: string) => ({
    v: root.querySelector<HTMLElement>(`[data-wk-${k}]`),
    bad: root.querySelector<HTMLElement>(`[data-wk-${k}-p90]`),
    delta: root.querySelector<HTMLElement>(`[data-wk-${k}-d]`),
  });
  const tiles = { ot: tile("ot"), idle: tile("idle"), late: tile("late"), waste: tile("waste") };

  const settings: Settings & { swaps: string[] } = { strategy: "day", swaps: [], sigma: Number(errIn?.value ?? 10) / 100, racks: Number(racksIn?.value ?? 3) };

  /* The menu card: every dish and how long it keeps; the three that have a variant that keeps are switches. */
  const dishRows = new Map<string, HTMLElement>();
  if (menuEl) {
    menuEl.innerHTML = "";
    for (const dish of MENU) {
      const li = document.createElement("li");
      const alt = SWAPS[dish.key];
      if (alt) {
        const b = document.createElement("button");
        b.type = "button";
        b.className = "wk-dish";
        b.dataset.wkSwap = dish.key;
        b.setAttribute("aria-pressed", "false");
        li.appendChild(b);
        dishRows.set(dish.key, b);
        b.addEventListener("click", () => {
          const on = !settings.swaps.includes(dish.key);
          settings.swaps = on ? [...settings.swaps, dish.key] : settings.swaps.filter((k) => k !== dish.key);
          update();
        });
      } else {
        const s = document.createElement("span");
        s.className = "wk-dish";
        li.appendChild(s);
        dishRows.set(dish.key, s);
      }
      menuEl.appendChild(li);
    }
  }
  const keepsHtml = (n: number) =>
    `<span class="wk-keeps" aria-label="${n === 0 ? "served the day it is cooked" : `keeps ${n} days chilled`}">${[1, 2, 3, 4].map((k) => `<i${k <= n ? ' class="on"' : ""}></i>`).join("")}<small>${n === 0 ? "same day" : n === 1 ? "1 day" : `${n} days`}</small></span>`;
  function renderMenu(): void {
    for (const dish of MENU) {
      const el = dishRows.get(dish.key);
      if (!el) continue;
      const swapped = settings.swaps.includes(dish.key);
      const d = swapped ? SWAPS[dish.key]! : dish;
      const hint = SWAPS[dish.key] ? `<em>${swapped ? "undo" : "swap"}</em>` : "";
      el.innerHTML = `<b>${d.name}</b>${keepsHtml(d.keeps)}${hint}`;
      if (el instanceof HTMLButtonElement) {
        el.setAttribute("aria-pressed", String(swapped));
        el.title = swapped ? `Back to ${dish.name.toLowerCase()}` : `Swap for ${SWAPS[dish.key]!.name.toLowerCase()}, which keeps ${SWAPS[dish.key]!.keeps} days`;
      }
    }
  }

  /* The chart: a column per day, the crew's 24 hours as a line across. */
  const CAP = KITCHEN.cooks * KITCHEN.shiftMin;
  const TOP = CAP * 1.45;
  type Col = { el: HTMLButtonElement; dayOf: HTMLElement; ahead: HTMLElement; ot: HTMLElement; idle: HTMLElement; tag: HTMLElement };
  const cols: Col[] = [];
  let arcs: SVGSVGElement | null = null;
  if (chartEl) {
    chartEl.innerHTML = "";
    const plot = document.createElement("div");
    plot.className = "wk-plot";
    const cap = document.createElement("div");
    cap.className = "wk-cap";
    cap.style.bottom = `${(CAP / TOP) * 100}%`;
    plot.appendChild(cap);
    arcs = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    arcs.setAttribute("class", "wk-arcs");
    arcs.setAttribute("viewBox", "0 0 600 100");
    arcs.setAttribute("preserveAspectRatio", "none");
    arcs.setAttribute("aria-hidden", "true");
    plot.appendChild(arcs);
    DAYS.forEach((day, d) => {
      const el = document.createElement("button");
      el.type = "button";
      el.className = "wk-col";
      el.innerHTML = `<span class="wk-bar"><i class="wk-idle"></i><i class="wk-dayof"></i><i class="wk-ahead"></i><i class="wk-ot"></i></span><span class="wk-tag"></span><span class="wk-day">${day}<small>${forecast[d]!.toLocaleString("en-US")}</small></span>`;
      el.addEventListener("click", () => jumpTo(d));
      plot.appendChild(el);
      cols.push({ el, dayOf: el.querySelector(".wk-dayof")!, ahead: el.querySelector(".wk-ahead")!, ot: el.querySelector(".wk-ot")!, idle: el.querySelector(".wk-idle")!, tag: el.querySelector(".wk-tag")! });
    });
    chartEl.appendChild(plot);
    const keys = document.createElement("p");
    keys.className = "wk-keys";
    keys.innerHTML = `<span class="dayof">For the day</span><span class="ahead">Cooked ahead</span><span class="ot">Overtime</span><span class="idle">Idle</span><span class="cap">Shift: ${KITCHEN.cooks} cooks × ${SHIFT_H} h</span><span class="note">Crew-hours, the forecast week</span>`;
    chartEl.appendChild(keys);
  }
  const pct = (min: number) => `${((Math.max(0, min) / TOP) * 100).toFixed(2)}%`;
  function renderChart(run: CommissaryRun): void {
    run.plan.days.forEach((r, d) => {
      const c = cols[d];
      if (!c) return;
      const worked = r.dayOfMin + r.aheadMin;
      const base = Math.min(r.dayOfMin, CAP);
      c.dayOf.style.height = pct(base);
      c.ahead.style.bottom = pct(base);
      c.ahead.style.height = pct(r.aheadMin);
      c.ot.style.bottom = pct(CAP);
      c.ot.style.height = pct(r.otMin);
      c.idle.style.bottom = pct(worked);
      c.idle.style.height = pct(r.idleMin);
      c.ot.classList.toggle("late", r.lateMeals > 0.5);
      c.tag.innerHTML = r.otMin > 1 ? `+${fmtH(r.otMin)}` : r.idleMin > 30 ? `${fmtH(r.idleMin)}<span class="w"> idle</span>` : "";
      c.tag.dataset.kind = r.otMin > 1 ? "ot" : "idle";
      c.tag.style.bottom = `calc(${pct(Math.max(worked, CAP))} + 0.3rem)`;
      const parts = [`${DAYS[d]}: ${Math.round(r.meals)} meals`, `${fmtH(r.dayOfMin)} for the day`];
      if (r.aheadMin > 0) parts.push(`${fmtH(r.aheadMin)} cooked ahead`);
      if (r.otMin > 0) parts.push(`${fmtH(r.otMin)} overtime`);
      if (r.idleMin > 0) parts.push(`${fmtH(r.idleMin)} idle`);
      c.el.setAttribute("aria-label", `${parts.join(", ")}. Show this day.`);
    });
    // Arcs from the day cooked for to the day it was cooked on.
    if (arcs) {
      const flows = new Map<string, number>();
      run.plan.days.forEach((r, d) => {
        for (const b of r.batches) flows.set(`${d}>${b.forDay}`, (flows.get(`${d}>${b.forDay}`) ?? 0) + b.minutes);
      });
      const colX = (d: number) => ((d + 0.5) / DAYS.length) * 600;
      let html = "";
      for (const [k, min] of flows) {
        const [a, b] = k.split(">").map(Number) as [number, number];
        const r = run.plan.days[a]!;
        const yA = 100 - ((Math.min(r.dayOfMin, CAP) + r.aheadMin / 2) / TOP) * 100;
        const rb = run.plan.days[b]!;
        const yB = 100 - ((Math.max(rb.dayOfMin + rb.aheadMin, CAP) / TOP) * 100) - 6;
        const xA = colX(a);
        const xB = colX(b);
        const lift = 14 + 6 * (b - a);
        const w = Math.max(1.2, Math.min(2.6, min / 120));
        html += `<path d="M ${xB.toFixed(1)} ${yB.toFixed(1)} C ${xB.toFixed(1)} ${(yB - lift).toFixed(1)}, ${xA.toFixed(1)} ${(Math.min(yA, yB) - lift).toFixed(1)}, ${xA.toFixed(1)} ${yA.toFixed(1)}" stroke-width="${w.toFixed(1)}" />`;
      }
      arcs.innerHTML = html;
    }
  }

  function renderTiles(run: CommissaryRun, base: CommissaryRun | null): void {
    const set = (t: (typeof tiles)[keyof typeof tiles], v: string, bad: string, delta: number | null, unit: string, digits = 0) => {
      if (t.v) t.v.textContent = v;
      if (t.bad) t.bad.textContent = bad;
      if (t.delta) {
        t.delta.hidden = delta === null || Math.abs(delta) < 0.5 * 10 ** -digits;
        if (delta !== null) {
          t.delta.textContent = `${delta > 0 ? "+" : "−"}${Math.abs(delta).toFixed(digits)}${unit}`;
          t.delta.dataset.good = String(delta < 0);
        }
      }
    };
    const h = (x: number) => (x >= 9.95 ? x.toFixed(0) : x.toFixed(1));
    set(tiles.ot, `${h(run.ot.p50)} h`, `${h(run.ot.p90)} h`, base ? run.ot.p50 - base.ot.p50 : null, " h", 1);
    set(tiles.idle, `${h(run.idle.p50)} h`, `${h(run.idle.p90)} h`, base ? run.idle.p50 - base.idle.p50 : null, " h", 1);
    set(tiles.late, `${Math.round(run.late.p50)}`, `${Math.round(run.late.p90)}`, base ? Math.round(run.late.p50 - base.late.p50) : null, "");
    set(tiles.waste, `${Math.round(run.wasted.p50)}`, `${Math.round(run.wasted.p90)}`, base ? Math.round(run.wasted.p50 - base.wasted.p50) : null, "");
  }

  function renderStatus(run: CommissaryRun): void {
    if (!statusEl) return;
    const days = run.plan.days;
    const worst = days.reduce((a, r, i) => (r.otMin > days[a]!.otMin ? i : a), 0);
    const slack = days.reduce((a, r, i) => (r.idleMin > days[a]!.idleMin ? i : a), 0);
    const w = days[worst]!;
    const s = days[slack]!;
    const util = Math.round(run.plan.utilization * 100);
    let text: string;
    let state: "ok" | "warn" | "danger";
    if (settings.strategy === "day") {
      text = `${DAY_NAMES[slack]} leaves ${fmtH(s.idleMin)} of crew time idle. ${w.otMin > 1 ? `${DAY_NAMES[worst]} needs ${fmtH(w.otMin)} of overtime${w.lateMeals > 0.5 ? ` and still sends ${Math.round(w.lateMeals)} meals out late` : ""}.` : "No day runs over."}`;
      state = w.lateMeals > 0.5 ? "danger" : w.otMin > 1 ? "warn" : "ok";
    } else {
      const big = days.reduce((a, r, i) => (r.aheadMin > days[a]!.aheadMin ? i : a), 0);
      const b = days[big]!;
      const pans = b.batches.reduce((x, y) => x + y.pans, 0);
      const names = [...new Set(b.batches.map((x) => x.forDay))].sort((x, y) => x - y).map((x) => DAY_NAMES[x]!);
      const fors = names.length > 1 ? `${names.slice(0, -1).join(", ")} and ${names[names.length - 1]}` : (names[0] ?? "");
      text = b.aheadMin > 0 ? `${DAY_NAMES[big]} cooks ${pans} pans ahead for ${fors}. ` : "Nothing on this menu keeps long enough to cook ahead. ";
      text += w.otMin > 1 ? `${DAY_NAMES[worst]} still needs ${fmtH(w.otMin)} of overtime.` : `No day runs over; the crew works ${util}% of its paid hours.`;
      state = w.lateMeals > 0.5 ? "danger" : w.otMin > 1 ? "warn" : "ok";
    }
    statusEl.textContent = text;
    statusEl.dataset.state = state;
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* Running it.                                                                                                 */
  /* ---------------------------------------------------------------------------------------------------------- */

  let run: CommissaryRun = runCommissary(settings);
  let spans = timeline(run);
  let weekT = 0;
  function buildRoutines(): void {
    for (let i = 0; i < cooks.length; i++) {
      const k = cooks[i]!;
      // Plan once to learn how long the step back and forward take, then again with the standing filled to fit.
      const trial = new Routine(k.st, { x: k.c.x, z: COOK_Z, yaw: 0 }, cookSteps(spans.days, i, k.at), true);
      const steps0 = cookSteps(spans.days, i, k.at);
      let back = 1;
      let fwd = 1;
      steps0.forEach((s, j) => {
        if (s.kind !== "shuffle") return;
        const dur = (trial.stepTimes[j + 1] ?? trial.duration) - trial.stepTimes[j]!;
        if (s.to[1] === BACK_Z) back = dur;
        else fwd = dur;
      });
      k.routine = new Routine(k.st, { x: k.c.x, z: COOK_Z, yaw: 0 }, cookSteps(spans.days, i, k.at, { back, fwd }), true);
    }
  }

  function update(): void {
    settings.strategy = (planBtns.find((b) => b.getAttribute("aria-pressed") === "true")?.dataset.wkPlan as Strategy) ?? settings.strategy;
    run = runCommissary(settings);
    const base = settings.strategy === "ahead" ? runCommissary({ ...settings, strategy: "day" }) : null;
    const frac = spans.total > 0 ? weekT / spans.total : 0;
    spans = timeline(run);
    buildRoutines();
    weekT = frac * spans.total;
    renderMenu();
    renderChart(run);
    renderTiles(run, base);
    renderStatus(run);
    if (racksVal) racksVal.textContent = `${settings.racks} ${settings.racks === 1 ? "rack" : "racks"} · ${(settings.racks ?? 3) * KITCHEN.rackPans} pans`;
    if (errVal) errVal.textContent = settings.sigma > 0 ? `±${Math.round(settings.sigma * 100)}%` : "Exact";
    lastDay = -1;
  }

  for (const b of planBtns)
    b.addEventListener("click", () => {
      for (const o of planBtns) o.setAttribute("aria-pressed", String(o === b));
      update();
    });
  racksIn?.addEventListener("input", () => {
    settings.racks = Number(racksIn.value);
    update();
  });
  errIn?.addEventListener("input", () => {
    settings.sigma = Number(errIn.value) / 100;
    update();
  });

  function jumpTo(d: number): void {
    weekT = spans.days[d]!.t0 + 0.01;
  }

  /* One moment of the week: the racks' pans, the cooks, the line, the carriers, the clock. */
  const M = new THREE.Matrix4();
  const Q = new THREE.Quaternion();
  const S = new THREE.Vector3(1, 1, 1);
  const V3 = new THREE.Vector3();
  let lastPans = "";
  let lastDay = -1;
  let lastClock = "";
  function show(t: number, dt: number): void {
    const days = spans.days;
    let d = days.length - 1;
    for (let i = 0; i < days.length; i++)
      if (t < days[i]!.t0 + days[i]!.len + NIGHT_S) {
        d = i;
        break;
      }
    const span = days[d]!;
    const local = t - span.t0;
    const h = Math.min(local, span.len) / SEC_PER_H;

    // Pans on each rack: the days gone by full, today's as far as its work has gone, the days to come empty.
    const counts = days.map((s, i) => {
      if (i < d) return [s.dayOf, s.ahead, s.over];
      if (i > d) return [0, 0, 0];
      // Work done so far, as a share of the day's whole: the rack is full when the last cook stops.
      const doneH = s.busyH.reduce((a, b) => a + Math.min(h, b), 0);
      const allH = s.busyH.reduce((a, b) => a + b, 0);
      const total = s.dayOf + s.ahead + s.over;
      const done = allH > 0 ? Math.min(total, Math.floor((doneH / allH) * total + 1e-6)) : 0;
      const nd = Math.min(done, s.dayOf);
      const na = Math.min(done - nd, s.ahead);
      return [nd, na, done - nd - na];
    });
    const key = counts.map((c) => c.join(",")).join("|");
    if (key !== lastPans) {
      lastPans = key;
      for (const o of open) o.begin();
      chilled.begin();
      overflow.begin();
      trays.begin();
      counts.forEach(([nd, na, no], i) => {
        const x = rackX(i) + RACK.w / 2;
        const z0 = RACK_Z - RACK.d / 2;
        const n = nd! + na!;
        for (let k = 0; k < n; k++) {
          const tier = Math.floor(k / 2);
          const y = SLIDES[1 + 2 * tier]!;
          if (k % 2 === 0) trays.push(M.compose(V3.set(x, y, z0), Q.identity(), S));
          const z = z0 + (k % 2 === 0 ? 1 : -1) * (HALF.d / 2 + 0.006);
          M.compose(V3.set(x, y - TRAY.h + 0.0018, z), Q.identity(), S);
          if (k < nd!) open[(i * 7 + k * 3) % 4]!.push(M);
          else chilled.push(M);
        }
        if (no! > 0) {
          trays.push(M.compose(V3.set(x, RACK.h + TRAY.h + 0.004, z0), Q.identity(), S));
          for (let k = 0; k < no!; k++) {
            const z = z0 + (k % 2 === 0 ? 1 : -1) * (HALF.d / 2 + 0.006);
            overflow.push(M.compose(V3.set(x, RACK.h + 0.006 + Math.floor(k / 2) * (HALF.h + 0.006), z), Q.identity(), S));
          }
        }
      });
      for (const o of open) o.end();
      chilled.end();
      overflow.end();
      trays.end();
    }
    // Today's rack in the light; the overflow glowing red while there is any.
    nowGlow.position.set(rackX(d) + RACK.w / 2, 0.002, RACK_Z - RACK.d / 2 + 0.1);
    const ot = counts[d]![2]!;
    const otMat = otGlow.material as THREE.MeshBasicMaterial;
    otMat.opacity = ot > 0 ? 0.6 : Math.max(0, otMat.opacity - dt * 1.5);
    if (ot > 0) otGlow.position.x = rackX(d) + RACK.w / 2;
    otGlow.position.z = RACK_Z - RACK.d / 2;
    otGlow.quaternion.copy(camera.quaternion);
    if (d !== lastDay) {
      labelMats.forEach((m, i) => (m.opacity = i === d ? 1 : 0.62));
      lastDay = d;
    }

    // The cooks, and their spoons while they plate.
    cooks.forEach((k, i) => {
      if (!k.routine) return;
      const state = k.routine.apply(t);
      const work = span.busyH[i]! * SEC_PER_H;
      k.spoon.visible = state.activity === "work" && local > 0.45 && local < work - 0.45;
    });

    // The line: each station by the minutes it cooks today, while the crew is at work.
    const busy = local < Math.max(...span.busyH) * SEC_PER_H ? 1 : 0;
    const sm = span.r.stationMin;
    const k1 = busy * clamp01(sm.kettle / 240);
    kettle.hooks.set?.("heat", k1);
    kettle.hooks.set?.("level", 0.25 + 0.6 * k1);
    skillet.hooks.set?.("heat", busy * clamp01(sm.skillet / 240));
    const c1 = busy * clamp01(sm.combi / 300);
    combi.hooks.set?.("heat", c1);
    combi.hooks.set?.("fan", c1 > 0 ? 0.4 + 0.6 * c1 : 0);
    combi.hooks.set?.("steam", 0.6 * c1);
    combi.hooks.tick?.(t, dt, c1);
    chiller.hooks.tick?.(t, dt, span.r.pansChilled > 0 && local < span.len ? 1 : 0.25);
    hold.hooks.set?.("heat", busy ? 0.8 : 0.3);
    hold.hooks.set?.("humidity", busy ? 0.5 : 0.2);
    hold.hooks.tick?.(t, dt, busy);

    // Carriers fill and leave: the day's meals packed so far, a carrier per 200.
    const packed = span.r.meals * clamp01(h / Math.max(0.1, Math.max(...span.busyH)));
    const nCar = local > span.len + NIGHT_S * 0.5 ? 0 : Math.min(carriers.length, Math.round(packed / 200));
    carriers.forEach((g, i) => (g.visible = i < nCar));

    // The clock.
    const c = `${DAYS[d]}|${clock(h)}|${h > SHIFT_H + 1e-6 ? "ot" : local > span.len ? "done" : "on"}`;
    if (c !== lastClock) {
      lastClock = c;
      if (hudDay) hudDay.textContent = DAYS[d]!;
      if (hudClock) hudClock.textContent = clock(h);
      if (hudNote) {
        const idle = span.busyH.filter((b) => h >= b - 1e-6).length;
        const note = local > span.len ? "Shift over" : h > SHIFT_H ? "Overtime" : idle === KITCHEN.cooks ? "Crew idle" : idle > 0 ? `${idle} of ${KITCHEN.cooks} cooks idle` : "Crew working";
        hudNote.textContent = note;
        hudNote.dataset.kind = h > SHIFT_H && local <= span.len ? "ot" : idle > 0 && local <= span.len ? "idle" : "on";
      }
      cols.forEach((col, i) => col.el.classList.toggle("now", i === d));
    }
  }

  update();
  show(0, 0);
  for (const k of cooks)
    void k.st.ready.then(() => {
      // The detailed body has its own hand: the spoon goes to it.
      k.st.hand("R").add(k.spoon);
      void gate.refresh(scene);
    });

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let frozen: number | null = null;
  let last = performance.now();
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    if (frozen === null) weekT = (weekT + dt) % spans.total;
    else weekT = frozen;
    show(weekT, dt);
    controls.update();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    week: {
      info: () => ({
        mountMs: Math.round(mountMs),
        tris: renderer.info.render.triangles,
        calls: renderer.info.render.calls,
        weekS: +spans.total.toFixed(2),
        loopS: cooks.map((k) => +(k.routine?.duration ?? 0).toFixed(2)),
        t: +weekT.toFixed(2),
        plan: run.plan.days.map((r) => ({ dayOf: Math.round(r.dayOfMin), ahead: Math.round(r.aheadMin), ot: Math.round(r.otMin), idle: Math.round(r.idleMin) })),
      }),
      hold(at: number | null) {
        frozen = at;
      },
      day(d: number) {
        jumpTo(d);
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.position.set(...eye);
        controls.target.set(...at);
        controls.update();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
