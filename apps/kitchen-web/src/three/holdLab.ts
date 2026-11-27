/**
 * Kept hot until lunch (chapter 06, third figure): yesterday's chili, reheated and packed hot into four insulated
 * carriers, six pans in two of them and a single pan in the other two; one of each pair primed first with a pan of
 * near-boiling water, the other loaded cold. The day runs from half an hour before the chili is done (the priming) to
 * 1 pm, through the van at 10:30 and lunch at 11:30.
 *
 * Each carrier is cut through the pans' middles, front to back: the chili's cut faces are coloured by the temperature
 * sim/holdHot.ts computes there, and the shell's cut, from its inner skin through the foam to its outer skin, by the
 * temperature across it. Beside the stage, what each carrier's coldest pan reads an inch in from its corner at lunch;
 * under it, that reading against the clock.
 */
import * as THREE from "three";
import { cToF, fToC } from "../sim/coolDown";
import {
  HOLD_WAYS,
  HOT_LINE_C,
  NX,
  NY,
  PAN,
  PRIME_MIN,
  RIM_Y,
  SLICE_MIN,
  WALL_NODES,
  WALL_POS,
  WALL_T,
  holdAll,
  sampleAt,
  sampleOf,
  type HoldKey,
  type HoldRun,
} from "../sim/holdHot";
import type { HoldAsk, HoldDone } from "../sim/holdRuns.worker";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { buildFixture, type Model, type Piece } from "./models/all";
import { palette } from "./models/materials";
import { makeStaff } from "./models/people";
import { Shop, rbox } from "./models/shop";
import { hotelPan, hotelPanInside } from "./models/smallwares";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

const CARRIER = { id: "NSF-CARRIER-INSULATED", w: 0.46, d: 0.635, h: 0.625 } as const;
/** In the carrier's own frame (three/models/support.ts insulatedCarrier): its shell's middle across, the pans' middle
 *  front to back (its front is at z = 0), and the shell's and the cavity's sections there. */
const MID_X = (0.016 + CARRIER.w - 0.014) / 2;
const PAN_Z = -CARRIER.d + 0.034 + 0.004 + 0.265;
const SHELL = { w: CARRIER.w - 0.03, h: CARRIER.h - 0.02, r: 0.034, cy: (0.02 + CARRIER.h) / 2 };
const CAV = { w: 0.36, h: 0.525, r: 0.026, cy: 0.06 + 0.525 / 2 };
/** A carrier stacked on another stands on its lugs this high. */
const STACK_Y = CARRIER.h - 0.01;
/** The hotel pan's floor below its rim, as it hangs on the rails. */
const PAN_SEAT = PAN.panM - 0.001;

/** Each carrier's colour, on the labels, the tiles and the chart. Kept off the temperature colours. */
export const HOLD_COLOR: Record<HoldKey, string> = { fullPrimed: "#4ade80", fullCold: "#e2e8f0", onePrimed: "#c084fc", oneCold: "#f472b6" };

/** The carriers close up on a laptop, a step coarser on a phone: chosen once, so what main.ts builds ahead is used. */
let detailNow: 1 | 2 | null = null;
const holdDetail = (): 1 | 2 => (detailNow ??= Math.min(window.innerWidth, window.innerHeight) < 600 ? 1 : 2);
export function holdPieces(): readonly Piece[] {
  const d = holdDetail();
  return HOLD_WAYS.map((): Piece => [CARRIER.id, CARRIER, d]);
}

/** The day (minutes after midnight): the van, lunch, and the stage's end. */
const VAN = [10 * 60 + 30, 11 * 60 + 15] as const;
const LUNCH = 11 * 60 + 30;
const DAY_END = 13 * 60;
const DONE_MIN = 6 * 60;
const DONE_MAX = 10 * 60 + 30;
/** The chart's clock: from the earliest priming to the stage's end. */
const CHART_T0 = DONE_MIN - PRIME_MIN;
/** The loop: the priming, the swap (the water out, the chili in, slowed to be seen), the hours (seconds an hour), a
 *  hold on the end (seconds). */
const PRIME_S = 2.4;
const SWAP_S = 1.8;
const HOUR_S = 3.6;
const END_S = 4;
/** The pans slide in over the first minutes after packing (and the water out over the last of priming). */
const SLIDE_MIN = 4;

/** A temperature as shown: rounded, but never up to 135 °F from under it, where the verdicts are drawn. */
const degF = (c: number) => `${c < HOT_LINE_C ? Math.min(134, Math.round(cToF(c))) : Math.round(cToF(c))} °F`;
const clockOf = (m: number) => {
  const a = ((Math.round(m) % 1440) + 1440) % 1440;
  const h = Math.floor(a / 60);
  return `${((h + 11) % 12) + 1}:${String(a % 60).padStart(2, "0")} ${h < 12 ? "am" : "pm"}`;
};
const span = (m: number) => {
  const r = Math.round(m);
  const h = Math.floor(r / 60);
  return h ? `${h} h${r % 60 ? ` ${r % 60} min` : ""}` : `${r} min`;
};
const smooth = (x: number) => {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
};

/** The temperature colours (°F), as the cooling figure's: sky at 41 °F, yellow at 70 °F, orange at 135 °F, red at 190. */
const HEAT: readonly (readonly [number, number, number, number])[] = [
  [-10, 0x1e, 0x3a, 0x8a],
  [32, 0x25, 0x63, 0xeb],
  [41, 0x38, 0xbd, 0xf8],
  [70, 0xfa, 0xcc, 0x15],
  [135, 0xf9, 0x73, 0x16],
  [190, 0xdc, 0x26, 0x26],
];
function heatRGB(c: number, out: Uint8Array, k: number): void {
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

/** A rounded rectangle's outline, `n` points a corner, centred on (0, cy), anticlockwise from the right side's middle. */
function roundRectPts(w: number, h: number, r: number, cy: number, n: number): [number, number][] {
  const out: [number, number][] = [];
  const cs: [number, number, number][] = [
    [w / 2 - r, h / 2 - r, 0],
    [-w / 2 + r, h / 2 - r, Math.PI / 2],
    [-w / 2 + r, -h / 2 + r, Math.PI],
    [w / 2 - r, -h / 2 + r, (3 * Math.PI) / 2],
  ];
  for (const [x, y, a0] of cs)
    for (let i = 0; i <= n; i++) {
      const a = a0 + (i / n) * (Math.PI / 2);
      out.push([x + r * Math.cos(a), cy + y + r * Math.sin(a)]);
    }
  return out;
}

export function mountHoldLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-hd-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-hd-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  renderer.localClippingEnabled = true;
  const phone = window.matchMedia("(max-width: 720px)").matches;
  const detail = holdDetail();
  /** On a phone the four stand in two stacks of two, the full carriers on top (so the lone pans' labels sit in their
   *  empty carriers, clear of the pans); else in a row. */
  const stacked = phone;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 6, 14);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.6;
  const lights = studioLights(scene, new THREE.Vector3(0, 0.4, -0.3), 2.4, phone ? 1024 : 2048);
  lights.key.position.set(-1.6, 4.2, 3.2);
  lights.key.intensity = 1.4;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(2.4, 2.8, -2.6);
  lights.rim.intensity = 0.55;
  lights.hemi.intensity = 0.5;
  {
    const floor = tilePlane(8, 6);
    floor.rotation.x = -Math.PI / 2;
    floor.position.set(0, 0, 0.6);
    scene.add(floor);
  }

  // Everything is cut through the pans' middles, its near half gone.
  const CUT_Z = PAN_Z;
  const cut = [new THREE.Plane(new THREE.Vector3(0, 0, -1), CUT_Z)];
  const clipCopies = new Map<THREE.Material, THREE.Material>();
  function clip(o: THREE.Object3D): void {
    const copy = (m: THREE.Material) => {
      let c = clipCopies.get(m);
      if (!c) {
        c = m.clone();
        c.clippingPlanes = cut;
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

  /* The pans: a 2-inch hotel pan of chili under its lid; the priming water's pan, open. Each turned to run front to
   * back, as the carrier's rails take it. */
  const ins = hotelPanInside("1/1", PAN.panM);
  const fullFill = Math.min(1, PAN.foodM / (PAN.panM - ins.floorY - 0.008));
  const foodTop = ins.foodY(fullFill);
  function panTemplate(water: boolean, cut: boolean): THREE.Group {
    const p = palette();
    const s = new Shop(detail, water ? "hold_water_pan" : "hold_chili_pan");
    hotelPan(s, "1/1", PAN.panM, water ? 0 : fullFill, 4);
    if (!water) {
      s.add(p.ss, rbox(s, PAN.lidL, 0.0018, PAN.lidW, 0.0008), [0, PAN.panM + 0.0012, 0]);
      if (detail > 1) s.add(p.ss, rbox(s, PAN.lidL - 0.06, 0.004, PAN.lidW - 0.05, 0.003), [0, PAN.panM + 0.0032, 0]);
    }
    const g = s.build(true).group;
    if (water) {
      const top = new THREE.Mesh(new THREE.PlaneGeometry(2 * ins.halfL(foodTop), 2 * ins.halfW(foodTop)), p.water);
      top.rotation.x = -Math.PI / 2;
      top.position.y = foodTop;
      g.add(top);
    }
    g.rotation.y = Math.PI / 2;
    const holder = new THREE.Group();
    holder.add(g);
    if (cut) clip(holder);
    return holder;
  }
  /** Each pan cut through on its rails, and whole while it slides in or out of the door. */
  const chiliPan = panTemplate(false, true);
  const chiliPanFree = panTemplate(false, false);
  const waterPan = panTemplate(true, true);
  const waterPanFree = panTemplate(true, false);
  /** The cut face of a pan's contents across its middle: its width at each height, from the floor to the top; the
   *  texture runs across it and up it. */
  const cutGeo = (() => {
    const ys = [ins.floorY, 0.0012, 0.0024, 0.0042, 0.0065, 0.0095, 0.013, 0.02, 0.03, 0.04].filter((y, i) => i === 0 || (y > ins.floorY && y < foodTop - 2e-4));
    ys.push(foodTop);
    const pts: THREE.Vector2[] = [];
    for (const y of ys) pts.push(new THREE.Vector2(ins.halfW(y), y));
    for (let i = ys.length - 1; i >= 0; i--) pts.push(new THREE.Vector2(-ins.halfW(ys[i]!), ys[i]!));
    const g = new THREE.ShapeGeometry(new THREE.Shape(pts));
    const pos = g.attributes.position as THREE.BufferAttribute;
    const uv = g.attributes.uv as THREE.BufferAttribute;
    for (let i = 0; i < pos.count; i++) {
      const y = pos.getY(i);
      const hw = ins.halfW(y);
      uv.setXY(i, (pos.getX(i) + hw) / (2 * hw), (y - ins.floorY) / (foodTop - ins.floorY));
    }
    return g;
  })();
  /** The shell's cut: rings from the cavity's outline out to the shell's, the texture running across the wall. */
  const RING_N = 8;
  const ringGeo = (() => {
    const inner = roundRectPts(CAV.w, CAV.h, CAV.r, CAV.cy, 6);
    const outer = roundRectPts(SHELL.w, SHELL.h, SHELL.r, SHELL.cy, 6);
    const P = inner.length;
    const pos: number[] = [];
    const uv: number[] = [];
    const idx: number[] = [];
    for (let l = 0; l <= RING_N; l++) {
      const f = l / RING_N;
      for (let i = 0; i < P; i++) {
        pos.push(inner[i]![0] + (outer[i]![0] - inner[i]![0]) * f, inner[i]![1] + (outer[i]![1] - inner[i]![1]) * f, 0);
        uv.push(f, 0.5);
      }
    }
    for (let l = 0; l < RING_N; l++)
      for (let i = 0; i < P; i++) {
        const a = l * P + i;
        const b = l * P + ((i + 1) % P);
        const c = (l + 1) * P + i;
        const d = (l + 1) * P + ((i + 1) % P);
        idx.push(a, b, c, b, d, c);
      }
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
    g.setIndex(idx);
    g.computeVertexNormals();
    return g;
  })();
  function heatTex(w: number, h: number): THREE.DataTexture {
    const t = new THREE.DataTexture(new Uint8Array(w * h * 4), w, h, THREE.RGBAFormat);
    t.colorSpace = THREE.SRGBColorSpace;
    t.magFilter = THREE.LinearFilter;
    t.minFilter = THREE.LinearFilter;
    t.wrapS = t.wrapT = THREE.ClampToEdgeWrapping;
    t.needsUpdate = true;
    return t;
  }

  /* The four carriers. */
  interface Carrier {
    key: HoldKey;
    model: Model;
    x: number;
    y: number;
    pans: { group: THREE.Group; free: THREE.Group; cutMat: THREE.MeshBasicMaterial; tex: THREE.DataTexture }[];
    water: { group: THREE.Group; free: THREE.Group; cutMat: THREE.MeshBasicMaterial } | null;
    ring: { mat: THREE.MeshBasicMaterial; tex: THREE.DataTexture };
  }
  const WALL_TEX = 48;
  const place = (i: number): [number, number] => (stacked ? [(i % 2 === 0 ? -1 : 1) * 0.255, i < 2 ? STACK_Y : 0] : [(i - 1.5) * 0.62, 0]);
  const carriers: Carrier[] = HOLD_WAYS.map((w, i) => {
    const [x, y] = place(i);
    const model = buildFixture(CARRIER.id, CARRIER, detail);
    model.hooks.set?.("pans", 0);
    model.hooks.open?.(0);
    model.group.position.set(x - MID_X, y, 0);
    clip(model.group);
    scene.add(model.group);
    const ringTex = heatTex(WALL_TEX, 1);
    const ringMat = new THREE.MeshBasicMaterial({ map: ringTex, toneMapped: false, side: THREE.DoubleSide, visible: false });
    const ring = new THREE.Mesh(ringGeo, ringMat);
    ring.position.set(x, y, CUT_Z + 0.0003);
    scene.add(ring);
    const pans = Array.from({ length: w.pans }, (_, s) => {
      const group = chiliPan.clone();
      const tex = heatTex(2 * NX, NY);
      const cutMat = new THREE.MeshBasicMaterial({ map: tex, toneMapped: false, visible: false });
      const face = new THREE.Mesh(cutGeo, cutMat);
      face.position.z = 0.0004;
      group.add(face);
      group.position.set(x, y + RIM_Y[s]! - PAN_SEAT, CUT_Z);
      group.visible = false;
      scene.add(group);
      const free = chiliPanFree.clone();
      free.position.copy(group.position);
      free.visible = false;
      scene.add(free);
      return { group, free, cutMat, tex };
    });
    let water: Carrier["water"] = null;
    if (w.primed) {
      const group = waterPan.clone();
      const cutMat = new THREE.MeshBasicMaterial({ color: 0xffffff, toneMapped: false, transparent: true, opacity: 0.85 });
      const face = new THREE.Mesh(cutGeo, cutMat);
      face.position.z = 0.0004;
      group.add(face);
      group.position.set(x, y + RIM_Y[0]! - PAN_SEAT, CUT_Z);
      group.visible = false;
      scene.add(group);
      const free = waterPanFree.clone();
      free.position.copy(group.position);
      free.visible = false;
      scene.add(free);
      water = { group, free, cutMat };
    }
    return { key: w.key, model, x, y, pans, water, ring: { mat: ringMat, tex: ringTex } };
  });

  // The cook, beside the carriers, an eye on them.
  const cook = makeStaff("p50", 2, { seed: 31, apron: 0x2b3a55, sleeves: "cuffed" });
  const COOK = stacked ? { x: 0.84, z: -0.42 } : { x: 1.42, z: -0.36 };
  cook.place(COOK.x, COOK.z, Math.atan2(-0.6 - COOK.x, 0.3 - COOK.z));
  scene.add(cook.group);

  /* The camera: from in front, a little raised, square on to the cuts. */
  const camera = new THREE.PerspectiveCamera(30, 16 / 9, 0.05, 40);
  const TARGET = new THREE.Vector3(0, stacked ? 0.62 : 0.32, CUT_Z);
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    VIEW.set(0.14, stacked ? 0.5 : 0.62, 1).normalize();
    keep.length = 0;
    for (const c of carriers) {
      keep.push(new THREE.Vector3(c.x - SHELL.w / 2, c.y, CUT_Z), new THREE.Vector3(c.x + SHELL.w / 2, c.y, CUT_Z));
      keep.push(new THREE.Vector3(c.x - SHELL.w / 2, c.y + CARRIER.h, -CARRIER.d), new THREE.Vector3(c.x + SHELL.w / 2, c.y + CARRIER.h, -CARRIER.d));
    }
    keep.push(new THREE.Vector3(COOK.x, 1.8, COOK.z));
  };
  const controls = createScrollSafeOrbit(camera, canvas, [0, 0.4, CUT_Z], { minDistance: 0.6, maxDistance: 9, maxPolarAngle: Math.PI * 0.48 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const hudEl = stageEl.querySelector<HTMLElement>(".bo-hud");
  const hudRoom = (): number => {
    if (!hudEl) return 0;
    const parts = ["[data-hd-hud-clock]", "[data-hd-hud-doing]", "[data-hd-hud-note]"].map((k) => hudEl.querySelector<HTMLElement>(k));
    const was = parts.map((e) => e?.textContent ?? "");
    let most = 0;
    for (const d of longestDoings()) {
      if (parts[0]) parts[0].textContent = "10:45 am";
      if (parts[1]) parts[1].textContent = d;
      if (parts[2]) parts[2].textContent = "4 of 4 at 135 °F or more";
      most = Math.max(most, hudEl.offsetTop + hudEl.offsetHeight);
    }
    parts.forEach((e, i) => {
      if (e) e.textContent = was[i]!;
    });
    return most;
  };
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
    let dist = 1;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 20; dist += 0.04) {
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
  const doneIn = q<HTMLInputElement>("[data-hd-done]");
  const doneVal = q("[data-hd-done-val]");
  const packBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-hd-pack]"));
  const statusEl = q("[data-hd-status]");
  const chartEl = q("[data-hd-chart]");
  const hudClock = q("[data-hd-hud-clock]");
  const hudDoing = q("[data-hd-hud-doing]");
  const hudNote = q("[data-hd-hud-note]");
  const figEl = stageEl.closest("figure");
  const tagEls = new Map(HOLD_WAYS.map((w) => [w.key, q<HTMLElement>(`[data-hd-tag="${w.key}"]`)] as const));
  const tileEls = new Map(HOLD_WAYS.map((w) => [w.key, { el: q(`[data-hd-tile-box="${w.key}"]`), v: q(`[data-hd-tile="${w.key}"]`), s: q(`[data-hd-tile-s="${w.key}"]`) }] as const));

  let done = Math.max(DONE_MIN, Math.min(DONE_MAX, Number(doneIn?.value ?? 420)));
  let startF = 185;
  let runs: HoldRun[] | null = null;
  let runsF = 0;
  const runOf = (k: HoldKey) => runs?.find((r) => r.way.key === k) ?? null;
  let loopT = 0;
  let frozen: number | null = null;
  /** The minute on show, after packing (negative while priming). */
  let playMin = 0;
  /** The lunch verdict: a carrier's coldest pan an inch in from its corner, at 11:30. */
  const atLunch = (r: HoldRun) => sampleOf(r.coldest, LUNCH - done);
  const okAt = (r: HoldRun) => atLunch(r) >= HOT_LINE_C;

  function longestDoings(): string[] {
    return ["Priming two carriers with near-boiling water", "Shut in the carriers until the van", "1 pan, primed: its corner under 135 °F", "Lunch: 1 pan, primed needs reheating"];
  }

  /* The runs: off the main thread, kept by packing temperature. */
  const cache = new Map<number, HoldRun[]>();
  let worker: Worker | null = null;
  try {
    worker = new Worker(new URL("../sim/holdRuns.worker.ts", import.meta.url), { type: "module" });
    worker.onmessage = (e: MessageEvent<HoldDone>) => {
      cache.set(e.data.startF, e.data.runs);
      // A second answer for what is already shown (asked again while the first was on its way) changes nothing.
      if (e.data.startF === startF && runsF !== startF) take(e.data.runs, e.data.startF);
    };
    worker.onerror = () => {
      worker = null;
      ask();
    };
  } catch {
    worker = null;
  }
  function ask(): void {
    const hit = cache.get(startF);
    if (hit) return take(hit, startF);
    figEl?.setAttribute("data-stale", "");
    if (statusEl) statusEl.textContent = "Working out the morning…";
    if (worker) worker.postMessage({ startF } satisfies HoldAsk);
    else {
      const r = holdAll(fToC(startF));
      cache.set(startF, r);
      take(r, startF);
    }
  }
  function take(r: HoldRun[], f: number): void {
    runs = r;
    runsF = f;
    figEl?.removeAttribute("data-stale");
    for (const c of carriers) {
      for (const p of c.pans) p.cutMat.visible = true;
      c.ring.mat.visible = true;
    }
    loopT = 0;
    renderAll();
    for (const t of tagEls.values()) if (t) t.dataset.last = "";
    void gate.refresh(scene);
    // The other packing temperature, ready for a press.
    const other = f === 185 ? 165 : 185;
    if (worker && !cache.has(other)) worker.postMessage({ startF: other } satisfies HoldAsk);
  }
  function renderAll(): void {
    renderControls();
    renderTiles();
    renderStatus();
    renderChart();
    lastHud = "";
  }
  function renderControls(): void {
    if (doneIn) doneIn.value = String(done);
    if (doneVal) doneVal.textContent = clockOf(done);
    for (const b of packBtns) b.setAttribute("aria-pressed", String(Number(b.dataset.hdPack) === startF));
  }
  function renderTiles(): void {
    if (!runs) return;
    for (const r of runs) {
      const t = tileEls.get(r.way.key);
      if (!t) continue;
      const ok = okAt(r);
      t.el?.setAttribute("data-ok", String(ok));
      if (t.v) t.v.textContent = degF(atLunch(r));
      const mid = Math.min(...r.centre.map((c) => sampleOf(c, LUNCH - done)));
      const lines = ["corner, at lunch", `middle ${degF(mid)}`];
      lines.push(ok ? (r.hotMin === null ? "Hot all day" : `Hot for ${span(r.hotMin)}`) : `Too cool at ${clockOf(done + r.hotMin!)}`);
      if (t.s) t.s.textContent = lines.join("\n");
    }
  }
  /** Each carrier, as a sentence names it. */
  const SAID: Record<HoldKey, string> = {
    fullPrimed: "the full carrier that was primed",
    fullCold: "the full carrier loaded cold",
    onePrimed: "the lone pan in a primed carrier",
    oneCold: "the lone pan in a cold carrier",
  };
  const cap = (t: string) => t.charAt(0).toUpperCase() + t.slice(1);
  function renderStatus(): void {
    // While another packing temperature is worked out, its "Working out" stays up.
    if (!statusEl || !runs || runsF !== startF) return;
    const ok = runs.filter(okAt);
    const late = runs.filter((r) => !okAt(r));
    const list = (rs: HoldRun[]) => {
      const has = (k: HoldKey) => rs.some((r) => r.way.key === k);
      const n: string[] = [];
      if (has("fullPrimed") && has("fullCold")) n.push("both full carriers");
      else for (const k of ["fullPrimed", "fullCold"] as const) if (has(k)) n.push(SAID[k]);
      if (has("onePrimed") && has("oneCold")) n.push("both lone pans");
      else for (const k of ["onePrimed", "oneCold"] as const) if (has(k)) n.push(SAID[k]);
      return n.length < 2 ? n.join("") : `${n.slice(0, -1).join(", ")} and ${n[n.length - 1]}`;
    };
    const parts: string[] = [];
    const when = `With the chili done at ${clockOf(done)} and packed at ${runsF} °F`;
    if (!late.length) parts.push(`${when}, all four carriers hold the chili at 135 °F or more until lunch.`);
    else if (!ok.length) parts.push(`${when}, no carrier holds all its chili at 135 °F or more until lunch.`);
    else parts.push(`${when}, ${list(ok)} ${ok.length === 1 ? "holds" : "hold"} the chili at 135 °F or more until lunch.`);
    for (const r of late) parts.push(`${cap(SAID[r.way.key])} falls under 135 °F at ${clockOf(done + r.hotMin!)}, and would need reheating to 165 °F before it is served.`);
    const fp = runOf("fullPrimed")!;
    const fc = runOf("fullCold")!;
    const op = runOf("onePrimed")!;
    const oc = runOf("oneCold")!;
    parts.push(`Six pans keep each other warm: the middle ones lose heat only sideways, through half a foot of chili, so a full carrier keeps every pan at 135 °F or more for ${Math.floor(Math.min(fp.hotMin!, fc.hotMin!) / 60)} hours or longer.`);
    const drop = (r: HoldRun) => Math.round(cToF(r.startC) - cToF(sampleOf(r.coldest, 60)));
    parts.push(`A lone pan has a whole carrier's walls to warm. In a cold carrier it is ${drop(oc)} °F cooler an inch in from its corner after the first hour; in a primed one, ${drop(op)} °F.`);
    if (runsF === 165) parts.push("Packed at 165 °F, the reheating minimum, the chili has 20 °F less to lose before 135 °F than at 185 °F.");
    statusEl.textContent = parts.join(" ");
    statusEl.dataset.state = late.length ? "warn" : "ok";
  }

  /* The chart: each carrier's coolest pan, an inch in from its corner, against the clock. */
  interface Live {
    key: HoldKey;
    clip: SVGRectElement | null;
    dot: SVGCircleElement | null;
  }
  let lives: Live[] = [];
  let chartX: ((m: number) => number) | null = null;
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
    const f0 = 110;
    const f1 = 190;
    const x = (m: number) => L + ((m - CHART_T0) / (DAY_END - CHART_T0)) * (R - L);
    const y = (f: number) => B - ((Math.min(f1, Math.max(f0, f)) - f0) / (f1 - f0)) * (B - T);
    chartX = x;
    chartY = y;
    chartL = L;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg bo-svg cd-svg hd-svg"><defs>`;
    for (const r of runs) svg += `<clipPath id="hd-clip-${r.way.key}"><rect x="0" y="0" width="0" height="${H}" data-hd-clip="${r.way.key}" /></clipPath>`;
    svg += `</defs>`;
    svg += `<rect class="van" x="${x(VAN[0]).toFixed(1)}" y="${T}" width="${(x(VAN[1]) - x(VAN[0])).toFixed(1)}" height="${B - T}" />`;
    svg += `<text class="van-t" x="${((x(VAN[0]) + x(VAN[1])) / 2).toFixed(1)}" y="${(T + 12).toFixed(1)}" text-anchor="middle">Van</text>`;
    svg += `<path class="late" d="M${L},${y(135).toFixed(1)}H${R}V${B}H${L}Z" />`;
    svg += `<text class="late-t" x="${(L + 6).toFixed(1)}" y="${(B - 6).toFixed(1)}">Under 135 °F: reheat before serving</text>`;
    svg += `<line class="lunch" x1="${x(LUNCH).toFixed(1)}" x2="${x(LUNCH).toFixed(1)}" y1="${T}" y2="${B}" />`;
    svg += `<text class="lunch-t" x="${(x(LUNCH) + 4).toFixed(1)}" y="${(T + 12).toFixed(1)}">Lunch</text>`;
    for (const f of [135, 165, 185]) {
      svg += `<line class="grid line" x1="${L}" x2="${R}" y1="${y(f).toFixed(1)}" y2="${y(f).toFixed(1)}" />`;
      svg += `<text class="ya" x="${L - 6}" y="${(y(f) + 4).toFixed(1)}" text-anchor="end">${f}°</text>`;
    }
    for (const h of narrow ? [6, 8, 10, 12] : [6, 7, 8, 9, 10, 11, 12, 13]) svg += `<text class="xa" x="${x(h * 60).toFixed(1)}" y="${H - 6}" text-anchor="middle">${((h + 11) % 12) + 1}${h < 12 ? " am" : " pm"}</text>`;
    const pathOf = (r: HoldRun) => {
      const pts: string[] = [];
      for (let m = 0; done + m <= DAY_END; m += 5) pts.push(`${x(done + m).toFixed(1)},${y(cToF(sampleOf(r.coldest, m))).toFixed(1)}`);
      return `M${pts.join("L")}`;
    };
    for (const r of runs) svg += `<path class="ln faint" d="${pathOf(r)}" style="--c:${HOLD_COLOR[r.way.key]}" />`;
    for (const r of runs) svg += `<g clip-path="url(#hd-clip-${r.way.key})"><path class="ln" d="${pathOf(r)}" style="--c:${HOLD_COLOR[r.way.key]}" /></g>`;
    for (const r of runs) if (r.hotMin !== null && done + r.hotMin <= DAY_END) svg += `<circle class="bp" cx="${x(done + r.hotMin).toFixed(1)}" cy="${y(135).toFixed(1)}" r="3.2" style="--c:${HOLD_COLOR[r.way.key]}" />`;
    for (const r of runs) svg += `<circle class="now" r="4.2" cx="-20" cy="-20" data-hd-dot="${r.way.key}" style="--c:${HOLD_COLOR[r.way.key]}" />`;
    svg += `</svg>`;
    chartEl.innerHTML = svg;
    lives = runs.map((r) => ({ key: r.way.key, clip: chartEl.querySelector<SVGRectElement>(`[data-hd-clip="${r.way.key}"]`), dot: chartEl.querySelector<SVGCircleElement>(`[data-hd-dot="${r.way.key}"]`) }));
    // Drawn to the minute on show now: the loop that keeps it there does not run while the stage is out of view.
    showChart(playMin);
  }
  function showChart(m: number): void {
    if (!runs || !chartX || !chartY) return;
    for (const lv of lives) {
      const r = runOf(lv.key);
      if (!r) continue;
      const on = m >= 0 && done + m <= DAY_END;
      lv.clip?.setAttribute("width", String(m <= 0 ? 0 : Math.max(chartL, chartX(Math.min(DAY_END, done + m)) + 1)));
      if (lv.dot) {
        lv.dot.setAttribute("cx", on ? chartX(done + m).toFixed(1) : "-20");
        lv.dot.setAttribute("cy", on ? chartY(cToF(sampleOf(r.coldest, m))).toFixed(1) : "-20");
      }
    }
  }

  /* The stage: the pans in or out, the cut faces and the shell's cut painted for minute `m` after packing. */
  /** How far out of the door a pan starts (or ends), its middle from the cut (m). */
  const SLIDE_OUT = 0.6;
  /** A pan `off` in front of its place: whole while it moves, cut through once on its rails. */
  function slide(p: { group: THREE.Group; free: THREE.Group }, off: number, shown: boolean): void {
    const moving = off > 1e-3;
    p.group.visible = shown && !moving;
    p.free.visible = shown && moving;
    p.free.position.z = CUT_Z + off;
  }
  const px4 = new Uint8Array(4);
  function paint(m: number): void {
    if (!runs) return;
    const x = Math.max(0, Math.min(sampleAt(SLICE_MIN) - 1e-6, sampleAt(m)));
    const k = Math.floor(x);
    const f = x - k;
    for (const c of carriers) {
      const r = runOf(c.key);
      if (!r) continue;
      // The pans slide in from the front as the chili is packed; the priming water slides out just before.
      const inK = smooth(m / SLIDE_MIN);
      c.pans.forEach((p, s) => slide(p, (1 - smooth(inK * 1.4 - s * 0.08)) * SLIDE_OUT, m >= 0));
      if (c.water) {
        const wv = m < 0;
        slide(c.water, smooth((m + SLIDE_MIN) / SLIDE_MIN) * SLIDE_OUT, wv);
        if (wv) {
          heatRGB(sampleOf(r.water, Math.max(-PRIME_MIN, Math.min(-SAMPLE_GUARD, m))), px4, 0);
          c.water.cutMat.color.setRGB(px4[0]! / 255, px4[1]! / 255, px4[2]! / 255, THREE.SRGBColorSpace);
        }
      }
      // The chili's cut: the run's middle slice, mirrored across the pan.
      const [nx, ny] = r.grid;
      const per = nx * ny;
      if (m >= 0)
        c.pans.forEach((p, s) => {
          const data = p.tex.image.data as Uint8Array;
          const a = (k * r.way.pans + s) * per;
          const b = ((k + 1) * r.way.pans + s) * per;
          for (let j = 0; j < ny; j++)
            for (let col = 0; col < 2 * nx; col++) {
              const i = col < nx ? nx - 1 - col : col - nx;
              const v0 = r.slice[a + j * nx + i]!;
              const v1 = r.slice[b + j * nx + i]!;
              heatRGB(Number.isNaN(v1) ? v0 : v0 + (v1 - v0) * f, data, (j * 2 * nx + col) * 4);
            }
          p.tex.needsUpdate = true;
        });
      // The shell: its nodes' temperatures by depth.
      const xw = Math.max(0, sampleAt(Math.max(-PRIME_MIN, m)));
      const kw = Math.min(Math.floor(xw), r.wall.length / WALL_NODES - 2);
      const fw = Math.min(1, xw - kw);
      const data = c.ring.tex.image.data as Uint8Array;
      for (let t = 0; t < WALL_TEX; t++) {
        const s = ((t + 0.5) / WALL_TEX) * WALL_T;
        let n = 0;
        while (n < WALL_NODES - 2 && s > WALL_POS[n + 1]!) n++;
        const g = Math.min(1, Math.max(0, (s - WALL_POS[n]!) / (WALL_POS[n + 1]! - WALL_POS[n]!)));
        const at = (kk: number) => r.wall[kk * WALL_NODES + n]! + (r.wall[kk * WALL_NODES + n + 1]! - r.wall[kk * WALL_NODES + n]!) * g;
        heatRGB(at(kw) + (at(kw + 1) - at(kw)) * fw, data, t * 4);
      }
      c.ring.tex.needsUpdate = true;
    }
  }
  /** Priming's last sample is 5 minutes before packing: hold the water's colour there through the slide out. */
  const SAMPLE_GUARD = 5;

  /* The readouts. */
  let lastHud = "";
  function doingAt(m: number): string {
    if (!runs) return "Working out the morning…";
    const clock = done + m;
    if (m < 0) return "Priming two carriers with near-boiling water";
    if (m < SLIDE_MIN || (m < 20 && clock < VAN[0])) return `Packed at ${runsF} °F, the doors shut`;
    if (clock >= VAN[0] && clock < VAN[1]) return "In the van, 10:30 to 11:15";
    if (clock >= LUNCH && clock < LUNCH + 30) {
      const cool = runs.filter((r) => !okAt(r));
      return !cool.length ? "Lunch: served from the carriers" : cool.length === 1 ? `Lunch: ${cool[0]!.way.short} needs reheating` : `Lunch: ${cool.length} carriers need reheating`;
    }
    // The latest carrier to fall under 135 °F within the last hour.
    let best = { t: -Infinity, s: "" };
    for (const r of runs) if (r.hotMin !== null && r.hotMin <= m && m - r.hotMin < 60 && r.hotMin > best.t) best = { t: r.hotMin, s: `${r.way.short}: its corner under 135 °F` };
    if (best.s) return best.s;
    return clock < VAN[0] ? "Shut in the carriers until the van" : clock < LUNCH ? "At the lunch, waiting to be served" : "Served";
  }
  function showHud(m: number): void {
    const hot = runs ? runs.filter((r) => sampleOf(r.coldest, Math.max(0, m)) >= HOT_LINE_C).length : 0;
    const doing = doingAt(m);
    const note = !runs || m < 0 ? "Priming" : `${hot} of 4 at 135 °F or more`;
    const clock = clockOf(done + m);
    const sig = `${clock}|${doing}|${note}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = clock;
    if (hudDoing) hudDoing.textContent = doing;
    if (hudNote) {
      hudNote.textContent = note;
      hudNote.dataset.kind = !runs || m < 0 || hot === 4 ? "hot" : "cold";
    }
  }
  /** The labels over each carrier: its coolest pan, an inch in from its corner. */
  const tagAt = new THREE.Vector3();
  const tagSize = new Map<HTMLElement, { w: number; h: number }>();
  const tagStem = new Map<HTMLElement, string>();
  const tagItems: { el: HTMLElement; x: number; y: number; w: number; h: number }[] = [];
  const tagPlaced: { l: number; r: number; t: number; b: number }[] = [];
  function showTags(m: number): void {
    const W = stageEl!.clientWidth;
    const H = stageEl!.clientHeight;
    tagItems.length = 0;
    for (const c of carriers) {
      const el = tagEls.get(c.key);
      if (!el) continue;
      const r = runOf(c.key);
      // Over the carrier's top; or, with another stacked on it, in its own cavity, over its one pan.
      const onTop = !stacked || c.y > 0;
      tagAt.set(c.x, c.y + (onTop ? CARRIER.h + 0.03 : CAV.cy + 0.06), onTop ? -CARRIER.d / 2 : CUT_Z).project(camera);
      const on = !!r && m >= 0 && tagAt.z < 1 && Math.abs(tagAt.x) < 0.98 && Math.abs(tagAt.y) < 0.98;
      el.hidden = !on;
      if (!on || !r) continue;
      const v = sampleOf(r.coldest, m);
      const text = degF(v);
      if (el.dataset.last !== text || !tagSize.has(el)) {
        el.dataset.last = text;
        const b = el.querySelector("b");
        if (b) b.textContent = text;
        el.toggleAttribute("data-boiled", v < HOT_LINE_C);
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

  renderControls();
  ask();
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
  doneIn?.addEventListener("input", () => {
    done = Math.max(DONE_MIN, Math.min(DONE_MAX, Math.round(Number(doneIn.value) / 15) * 15));
    loopT = 0;
    renderAll();
  });
  for (const b of packBtns)
    b.addEventListener("click", () => {
      const f = Number(b.dataset.hdPack);
      if (f === startF) return;
      startF = f;
      renderControls();
      ask();
    });
  void cook.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  /** The loop's minute after packing: the priming, the swap, then the hours to 1 pm, then a hold. */
  const loopS = () => PRIME_S + SWAP_S + ((DAY_END - done - SLIDE_MIN) / 60) * HOUR_S + END_S;
  function minuteAt(s: number): number {
    if (s < PRIME_S) return -PRIME_MIN + (s / PRIME_S) * (PRIME_MIN - SLIDE_MIN);
    const w = s - PRIME_S;
    if (w < SWAP_S) return -SLIDE_MIN + (w / SWAP_S) * 2 * SLIDE_MIN;
    return Math.min(DAY_END - done, SLIDE_MIN + ((w - SWAP_S) / HOUR_S) * 60);
  }
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    loopT = (loopT + dt) % loopS();
    const t = now / 1000;
    playMin = frozen !== null ? Math.min(DAY_END - done, Math.max(-PRIME_MIN, frozen)) : minuteAt(loopT);
    paint(playMin);
    cook.stand(t);
    showHud(playMin);
    showChart(playMin);
    controls.update();
    fitClip();
    renderer.render(scene, camera);
    showTags(playMin);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    hold: {
      info: () => ({
        mountMs: Math.round(mountMs),
        done,
        startF,
        ready: runs !== null && runsF === startF,
        t: +playMin.toFixed(1),
        runs: (runs ?? []).map((r) => ({
          key: r.way.key,
          hotMin: r.hotMin === null ? null : +r.hotMin.toFixed(1),
          lunchF: +cToF(atLunch(r)).toFixed(1),
          ok: okAt(r),
        })),
        tags: HOLD_WAYS.map((w) => {
          const el = tagEls.get(w.key);
          return el && !el.hidden ? el.textContent : null;
        }),
        stale: figEl?.hasAttribute("data-stale") ?? false,
        status: statusEl?.textContent ?? "",
      }),
      /** Hold the stage at minute `min` after packing (null: play). */
      hold(min: number | null) {
        frozen = min;
      },
      view(eye: [number, number, number] | null, at?: [number, number, number]) {
        userMoved = !!eye;
        if (eye && at) {
          camera.clearViewOffset();
          camera.position.set(...eye);
          controls.target.set(...at);
          controls.update();
        } else frame();
      },
      set(s: { done?: number; startF?: number }) {
        if (s.done !== undefined) done = Math.max(DONE_MIN, Math.min(DONE_MAX, Math.round(s.done / 15) * 15));
        if (s.startF !== undefined) startF = s.startF === 165 ? 165 : 185;
        renderControls();
        ask();
        if (runs && runsF === startF) renderAll();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
