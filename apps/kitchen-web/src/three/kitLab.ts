/**
 * Buy one thing (chapter 07, fourth figure): the morning's kitchen with room at the ends of its line for a third combi
 * and a second six-burner range, and a fourth cook to hire. Bought equipment rolls in along the line, a cook hired
 * joins the crew, and the morning is planned again with it (planMorning, run in a worker; sim/buyOne.ts).
 *
 * The floor in front of a station glows red while every place at it is taken: both combis running, both halves of
 * the range worked (or all six burners lit), all three table places. A cook with nothing to do while work waits for
 * a station with no room gets a pause badge. Under the kitchen, how full each station is minute by minute, and what
 * each purchase does to the most meals the van takes at 10:30 with every portion at its best.
 */
import * as THREE from "three";
import { DOOR, KIT, MENUS, RANGE_BURNERS, RANGE_PLACES, SLOTS, TABLE_SPOTS, clockOf, kitOf, madeIt, menuDishes, type Kit, type Morning, type Station } from "../sim/cateringMorning";
import { BUYS, buyKey, limitKey, planKey, type Buy, type KitSettings } from "../sim/buyOne";
import type { KitAsk, KitDone } from "../sim/kitRuns.worker";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { CARRIER, CARRIER_COLS, CARRIER_X0, CUTOFF, END_HOLD_S, LINE_Z, MorningKitchen, OVEN, RANGE, ROLL_MIN, SEC_PER_MIN, SLOT_COLOR, TABLE, TABLE_Z, fmtLate, kitchenPieces, ovenX, ovenX0, rangeSideX, rangeX0, type CookSpec, type KitchenMoment } from "./morningKitchen";
import type { Piece } from "./models/all";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

/** The most the kitchen can hold: one more combi and one more range than it was built with. */
const ROOM: Kit = { ovens: KIT.ovens + 1, ranges: KIT.ranges + 1 };
/** The pieces the figure builds, for main.ts to build ahead. */
export const KIT_PIECES: readonly Piece[] = kitchenPieces(ROOM);

/** The morning's three cooks, and a fourth to hire. */
const COOKS: CookSpec[] = [
  { cohort: "p50", seed: 41, apron: 0x26344d },
  { cohort: "p95", seed: 57, apron: 0x3b3d43 },
  { cohort: "p05", seed: 63, apron: 0x45563f },
  { cohort: "p50", seed: 71, apron: 0x4a3a2e },
];

/** The stations a purchase can add room at, and the table and the cooks' hands, as the page names them. */
type Kind = "oven" | "range" | "table" | "hands";
const KIND_NAME: Record<Kind, string> = { oven: "The combis", range: "The range", table: "The table", hands: "The cooks" };
const KIND_IS: Record<Kind, string> = { oven: "are", range: "is", table: "is", hands: "are" };
/** The name of a kind, and its verb, for a kitchen with `ranges` ranges. */
const kindName = (kind: Kind, ranges: number) => (kind === "range" && ranges > 1 ? "The ranges" : KIND_NAME[kind]);
const kindIs = (kind: Kind, ranges: number) => (kind === "range" && ranges > 1 ? "are" : KIND_IS[kind]);
/** What full means for each: every place taken, or for the cooks, every one at work. */
const KIND_FULL: Record<Kind, string> = { oven: "full", range: "full", table: "full", hands: "all busy" };
const WAIT_FOR: Record<Station, string> = { oven: "a combi", range: "a place at the range", table: "a place at the table" };
const waitFor = (at: Station, ranges: number) => (at === "range" && ranges > 1 ? "a place at a range" : WAIT_FOR[at]);
const BUY_NAME: Record<keyof Buy, string> = { combi: "a third combi", range: "a second range", cook: "a fourth cook" };
const BUY_SHORT: Record<keyof Buy, string> = { combi: "the combi", range: "the range", cook: "the cook" };
/** For the ladder's row of what is bought now, short enough for its label column. */
const BUY_WORD: Record<keyof Buy, string> = { combi: "combi", range: "range", cook: "cook" };
const BUY_KEYS = ["combi", "range", "cook"] as const;
const num = (v: number) => Math.round(v).toLocaleString("en-US");
/** "a", "a and b", "a, b and c". */
const listOf = (xs: string[]) => (xs.length < 2 ? (xs[0] ?? "") : `${xs.slice(0, -1).join(", ")} and ${xs[xs.length - 1]}`);
const bought = (b: Buy) => BUY_KEYS.filter((k) => b[k]);
const sameBuy = (a: Buy, b: Buy) => buyKey(a) === buyKey(b);
/** Plans kept on the page for settings seen before. */
const PLANS_KEPT = 48;
/** The red of a full station, the amber of a cook waiting for room. */
const RED = "#ef4444";
const AMBER = "#fbbf24";

/**
 * How full each station is, minute by minute (sampled at each minute's middle) over a morning: combis running,
 * range places worked and ranges full (both places worked, or every burner lit), table places, the cooks' hands at
 * work, and cooks free while work waits for room.
 */
interface Fill {
  n: number;
  ovens: Uint8Array;
  places: Uint8Array;
  rangesFull: Uint8Array;
  table: Uint8Array;
  hands: Float32Array;
  waiting: Uint8Array;
  /** The station the most pressing waiting work is for, per minute. */
  waitAt: (Station | null)[];
}
function fillOf(m: Morning, endMin: number): Fill {
  const kit = kitOf(m.settings);
  const n = Math.ceil(endMin) + 1;
  const f: Fill = { n, ovens: new Uint8Array(n), places: new Uint8Array(n), rangesFull: new Uint8Array(n), table: new Uint8Array(n), hands: new Float32Array(n), waiting: new Uint8Array(n), waitAt: new Array(n).fill(null) };
  /** Each minute whose middle falls in [t0, t1). */
  const each = (t0: number, t1: number, fn: (k: number) => void) => {
    for (let k = Math.max(0, Math.ceil(t0 - 0.5)); k < n && k + 0.5 < t1; k++) fn(k);
  };
  for (const l of m.loads) each(l.t0, l.t3, (k) => f.ovens[k]!++);
  const sides = new Uint8Array(n);
  const hands = Array.from({ length: m.settings.cooks }, () => new Float32Array(n));
  for (const j of m.jobs) {
    each(j.t0, j.t1, (k) => (hands[j.cook]![k]! += j.hands));
    if (j.at === "range") each(j.t0, j.t1, (k) => (sides[k]! |= 1 << (j.side ?? 0)));
    if (j.at === "table") each(j.t0, j.t1, (k) => f.table[k]!++);
  }
  const lit = Array.from({ length: kit.ranges }, () => new Uint8Array(n));
  for (const b of m.burns) each(b.t0, b.t1, (k) => lit[Math.floor(b.burner / RANGE_BURNERS)]![k]!++);
  for (let k = 0; k < n; k++) {
    let p = 0;
    let full = 0;
    for (let r = 0; r < kit.ranges; r++) {
      const mine = (sides[k]! >> (r * RANGE_PLACES)) & ((1 << RANGE_PLACES) - 1);
      const worked = (mine & 1) + ((mine >> 1) & 1);
      p += worked;
      if (worked === RANGE_PLACES || lit[r]![k]! >= RANGE_BURNERS) full++;
    }
    f.places[k] = p;
    f.rangesFull[k] = full;
    let h = 0;
    for (const x of hands) h += Math.min(1, x[k]!);
    f.hands[k] = h;
  }
  for (const w of m.waits)
    each(w.t0, w.t1, (k) => {
      f.waiting[k] = w.cooks.length;
      f.waitAt[k] = w.at;
    });
  return f;
}
/** Whether a station is full in minute k. */
function fullAt(f: Fill, kit: Kit, cooks: number, kind: Kind, k: number): boolean {
  if (k < 0 || k >= f.n) return false;
  if (kind === "oven") return f.ovens[k]! >= kit.ovens;
  if (kind === "range") return f.rangesFull[k]! >= kit.ranges;
  if (kind === "table") return f.table[k]! >= TABLE_SPOTS;
  return f.hands[k]! >= cooks - 1e-6;
}

export function mountKitLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-kit-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-kit-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 8, 16);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.6;
  const lights = studioLights(scene, new THREE.Vector3(-0.6, 0.9, 0.1), 5.2, phone ? 1024 : 2048);
  lights.key.position.set(-1.4, 6.4, 4.6);
  lights.key.intensity = 1.9;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.4, 3.2, -3.0);
  lights.rim.intensity = 0.75;
  lights.hemi.intensity = 0.5;
  // Wider than the morning's floor: the line runs a combi and a range longer.
  const floor = tilePlane(14, 8);
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(-0.5, 0, 0.2);
  scene.add(floor);

  const kitchen = new MorningKitchen(COOKS, ROOM);
  kitchen.rollOut = true;
  scene.add(kitchen.group);

  /* The floor in front of each station: red while every place at that station is taken. */
  const padTex = (() => {
    const S = 64;
    const cv = document.createElement("canvas");
    cv.width = cv.height = S;
    const x = cv.getContext("2d")!;
    const img = x.createImageData(S, S);
    for (let j = 0; j < S; j++)
      for (let i = 0; i < S; i++) {
        // Soft at the edges: full inside, fading over the outer fifth.
        const e = Math.min(i + 0.5, S - i - 0.5, j + 0.5, S - j - 0.5) / (S * 0.2);
        const a = Math.max(0, Math.min(1, e));
        const o = 4 * (j * S + i);
        img.data[o] = img.data[o + 1] = img.data[o + 2] = 255;
        img.data[o + 3] = Math.round(255 * a * a * (3 - 2 * a));
      }
    x.putImageData(img, 0, 0);
    return new THREE.CanvasTexture(cv);
  })();
  type Pad = { mesh: THREE.Mesh; kind: Kind; piece: number; on: number };
  const pads: Pad[] = [];
  const pad = (kind: Kind, piece: number, x: number, z: number, w: number, d: number) => {
    const mat = new THREE.MeshBasicMaterial({ color: RED, map: padTex, toneMapped: false, transparent: true, opacity: 0, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -2 });
    const mesh = new THREE.Mesh(new THREE.PlaneGeometry(w, d).rotateX(-Math.PI / 2), mat);
    mesh.position.set(x, 0.004, z);
    mesh.renderOrder = 1;
    mesh.visible = false;
    scene.add(mesh);
    pads.push({ mesh, kind, piece, on: 0 });
  };
  for (let o = 0; o < ROOM.ovens; o++) pad("oven", o, ovenX(o), LINE_Z + 0.36, OVEN.w * 0.96, 0.62);
  for (let s = 0; s < ROOM.ranges * RANGE_PLACES; s++) pad("range", Math.floor(s / RANGE_PLACES), rangeSideX(s), LINE_Z + 0.36, (RANGE.w / RANGE_PLACES) * 0.96, 0.62);
  for (let k = 0; k < TABLE_SPOTS; k++) pad("table", 0, TABLE.x0 + (TABLE.w * (k + 0.5)) / TABLE_SPOTS, TABLE_Z - 0.32, (TABLE.w / TABLE_SPOTS) * 0.96, 0.56);

  /* A cook free while work waits for a station with no room: a pause badge over their head. */
  const badgeMat = (() => {
    const cv = document.createElement("canvas");
    cv.width = cv.height = 96;
    const x = cv.getContext("2d")!;
    x.fillStyle = AMBER;
    x.strokeStyle = "rgba(255,255,255,0.92)";
    x.lineWidth = 6;
    x.beginPath();
    x.arc(48, 48, 40, 0, Math.PI * 2);
    x.fill();
    x.stroke();
    x.fillStyle = "#1a1206";
    x.fillRect(34, 30, 9, 36);
    x.fillRect(53, 30, 9, 36);
    const t = new THREE.CanvasTexture(cv);
    t.colorSpace = THREE.SRGBColorSpace;
    return new THREE.SpriteMaterial({ map: t, toneMapped: false, depthTest: false, depthWrite: false });
  })();
  const badges = kitchen.cooks.map(() => {
    const b = new THREE.Sprite(badgeMat);
    b.scale.setScalar(0.26);
    b.renderOrder = 5;
    b.visible = false;
    scene.add(b);
    return b;
  });

  /* Camera: from behind the worktable, raised, the whole line across the aisle, room at its ends included. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 9, 0.05, 60);
  const TARGET = new THREE.Vector3();
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const LEFT = ovenX0(ROOM.ovens - 1);
  const RIGHT = Math.max(rangeX0(ROOM.ranges - 1) + RANGE.w, CARRIER_X0 + CARRIER_COLS * (CARRIER.w + 0.03));
  const aim = () => {
    keep.length = 0;
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    TARGET.set((LEFT + RIGHT) / 2, 0.75, 0.0);
    VIEW.set(0.06, wide ? 1.15 : 1.35, 1).normalize();
    // The line's back top edge (the combis are its tallest), and the floor in front of the table.
    for (const x of [LEFT, RIGHT]) keep.push(new THREE.Vector3(x, OVEN.h + 0.04, LINE_Z - OVEN.d), new THREE.Vector3(x, 0.2, TABLE_Z + TABLE.d));
    keep.push(new THREE.Vector3(RIGHT, 3 * CARRIER.h, TABLE_Z));
  };
  const controls = createScrollSafeOrbit(camera, canvas, [0, 0.75, 0], { minDistance: 1.2, maxDistance: 16, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const ndc = new THREE.Vector3();
  const hudEl = stageEl.querySelector<HTMLElement>(".kt-hud");
  const legendEl = stageEl.querySelector<HTMLElement>(".kt-legend");
  /** The HUD's lowest edge so far at this size: a HUD that wraps for a moment and back must not refit the view. */
  let hudMax = 0;
  let fitW = 0;
  let fitH = 0;
  const frame = () => {
    if (userMoved) return;
    aim();
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    if (w !== fitW || h !== fitH) hudMax = 0;
    fitW = w;
    fitH = h;
    hudMax = Math.max(hudMax, hudEl ? hudEl.offsetTop + hudEl.offsetHeight : 0);
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    // The kitchen fits between the HUD and the legend: fitted to that band's height and the stage's width, then
    // shifted into the band and across to the middle.
    const top = Math.min(0.86, 1 - (2 * (hudMax + 6)) / h);
    const bot = -Math.min(0.92, 1 - (2 * ((legendEl ? h - legendEl.offsetTop : 0) + 4)) / h);
    let dist = 2.5;
    let lo = 0;
    let hi = 0;
    let left = 0;
    let right = 0;
    for (; dist < 22; dist += 0.06) {
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
    fog.near = dist * 1.1;
    fog.far = dist * 2.6;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  frame();
  new ResizeObserver(() => {
    const w = stageEl.clientWidth;
    const h = stageEl.clientHeight;
    if (w === fitW && h === fitH && (hudEl ? hudEl.offsetTop + hudEl.offsetHeight : 0) <= hudMax) return;
    frame();
  }).observe(stageEl);
  if (hudEl)
    new ResizeObserver(() => {
      if (hudEl.offsetTop + hudEl.offsetHeight > hudMax) frame();
    }).observe(hudEl);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */

  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const buyBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-kt-buy]"));
  const menuEl = q("[data-kt-menu]");
  const mealsIn = q<HTMLInputElement>("[data-kt-meals]");
  const mealsVal = q("[data-kt-meals-val]");
  const statusEl = q("[data-kt-status]");
  const chartEl = q("[data-kit-chart]");
  const ladderEl = q("[data-kit-ladder]");
  const ladderNote = q("[data-kit-ladder-note]");
  const hudClock = q("[data-kt-hud-clock]");
  const hudPans = q("[data-kt-hud-pans]");
  const hudNote = q("[data-kt-hud-note]");
  const tile = (k: string) => ({ v: q(`[data-kt-${k}]`), d: q(`[data-kt-${k}-d]`), s: q(`[data-kt-${k}-s]`) });
  const tiles = { van: tile("van"), limit: tile("limit"), full: tile("full"), wait: tile("wait") };
  const tilesEl = q(".kt-lab .wk-tiles");
  const setTile = (t: (typeof tiles)[keyof typeof tiles], v: string, s: string, d: string | null = null, good = true) => {
    if (t.v) t.v.textContent = v;
    if (t.s) t.s.textContent = s;
    if (t.d) {
      t.d.hidden = d === null;
      if (d !== null) {
        t.d.textContent = d;
        t.d.dataset.good = String(good);
      }
    }
  };

  // Every dish swapped: the menu that leans on the combis.
  const settings: KitSettings = { meals: Number(mealsIn?.value ?? 1400), menu: MENUS - 1, combi: false, range: false, cook: false };
  const plans = new Map<string, Extract<KitDone, { kind: "plan" }>>();
  const limits = new Map<string, number>();
  /** The plan on screen, its settings' key, and how full its stations are (it stays while the next is planned). */
  let shown: Extract<KitDone, { kind: "plan" }> | null = null;
  let shownKey = "";
  let fill: Fill | null = null;
  const ready = () => shown !== null && shownKey === planKey(settings);
  const limitOf = (menu: number, b: Buy) => limits.get(limitKey(menu, b));

  /* The menu card: each place's dish, where it cooks and how long it holds; a press swaps it. */
  const slotBtns: HTMLButtonElement[] = [];
  if (menuEl) {
    menuEl.innerHTML = "";
    SLOTS.forEach((_, k) => {
      const li = document.createElement("li");
      const b = document.createElement("button");
      b.type = "button";
      b.className = "wk-dish mo-dish";
      b.style.setProperty("--dish", SLOT_COLOR[k]!);
      b.addEventListener("click", () => {
        settings.menu ^= 1 << k;
        changed(true);
      });
      li.appendChild(b);
      menuEl.appendChild(li);
      slotBtns.push(b);
    });
  }
  const holdBars = (min: number) => {
    const n = min >= 210 ? 4 : min >= 140 ? 3 : min >= 85 ? 2 : 1;
    return `<span class="wk-keeps" aria-hidden="true">${[1, 2, 3, 4].map((k) => `<i${k <= n ? ' class="on"' : ""}></i>`).join("")}<small>${min >= 120 ? `${min / 60} h` : `${min} min`}</small></span>`;
  };
  /** Where a dish cooks: what the plan needs room for. */
  const whereOf = (k: number) => {
    const st = menuDishes(settings.menu)[k]!.stages;
    const range = st.some((x) => x.at === "range");
    const oven = st.some((x) => x.at === "oven");
    if (range && oven) return ["Range, then combi", "both"];
    if (range) return st.some((x) => x.at === "range" && (x.hands ?? 1) >= 1) ? ["Range, never left", "range"] : ["Range, 2 pans a cook", "range"];
    return ["Combi, walk away", "oven"];
  };
  function renderControls(): void {
    const dishes = menuDishes(settings.menu);
    slotBtns.forEach((b, k) => {
      const d = dishes[k]!;
      const swapped = (settings.menu & (1 << k)) !== 0;
      const other = swapped ? SLOTS[k]!.written : SLOTS[k]!.swap;
      const [where, at] = whereOf(k);
      b.innerHTML = `<b>${d.name}</b>${holdBars(d.holdMin)}<em data-at="${at}">${where}</em>`;
      b.setAttribute("aria-pressed", String(swapped));
      b.title = `Swap for ${other.name.toLowerCase()}`;
    });
    for (const b of buyBtns) b.setAttribute("aria-pressed", String(settings[b.dataset.ktBuy as keyof Buy] === true));
    if (mealsVal) mealsVal.textContent = `${num(settings.meals)} meals`;
  }

  /* The chart: how full each station is, minute by minute, for the plan on screen. */
  let X0 = 120;
  let X1 = 600;
  let T1 = 330;
  const xOf = (m: number) => X0 + (Math.max(0, Math.min(m, T1)) / T1) * (X1 - X0);
  let cursor: SVGLineElement | null = null;
  const ROW = 26;
  function renderChart(): void {
    if (!chartEl) return;
    const stale = !ready();
    const W = Math.max(300, Math.round(chartEl.clientWidth - 2 * 12));
    const narrow = W < 480;
    // Room for "Range places · 4" in the 10.5 px mono labels.
    X0 = narrow ? 112 : 132;
    X1 = W - 8;
    const m = shown?.m ?? null;
    // On to half an hour past the van, and further while its "Leaves" label would not fit beside its line.
    T1 = Math.max(300, Math.ceil(((m?.depart ?? DOOR) + 8) / 30) * 30);
    while (m && m.late > 0.5 && xOf(m.depart) + 82 > X1) T1 += 30;
    const kit = m ? kitOf(m.settings) : KIT;
    const cooks = m?.settings.cooks ?? 3;
    const rows: { kind: Kind; label: string; cap: number; used: (k: number) => number }[] = [
      { kind: "oven", label: `Combis · ${kit.ovens}`, cap: kit.ovens, used: (k) => fill!.ovens[k]! },
      { kind: "range", label: `Range places · ${kit.ranges * RANGE_PLACES}`, cap: kit.ranges * RANGE_PLACES, used: (k) => (fullAt(fill!, kit, cooks, "range", k) ? kit.ranges * RANGE_PLACES : fill!.places[k]!) },
      { kind: "table", label: `Table places · ${TABLE_SPOTS}`, cap: TABLE_SPOTS, used: (k) => fill!.table[k]! },
      { kind: "hands", label: `Cooks · ${cooks}`, cap: cooks, used: (k) => fill!.hands[k]! },
    ];
    const y0 = 22;
    const axisY = y0 + rows.length * (ROW + 6);
    const H = axisY + 22;
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg mo-svg kt-svg${stale ? " stale" : ""}" aria-hidden="true">`;
    if (m && m.late > 0.5) svg += `<rect class="late" x="${xOf(DOOR)}" y="${y0 - 4}" width="${xOf(m.depart) - xOf(DOOR)}" height="${axisY - y0 + 4}"/>`;
    const step = W < 460 ? 60 : 30;
    for (let t = 0; t <= T1; t += step) {
      const label = t % (W < 400 ? 120 : 60) === 0 && !(xOf(t) > X1 - 14);
      svg += `<line class="g v" x1="${xOf(t)}" x2="${xOf(t)}" y1="${y0 - 4}" y2="${axisY}"/>${label ? `<text class="xa" x="${xOf(t)}" y="${axisY + 15}">${clockOf(t)}</text>` : ""}`;
    }
    rows.forEach((row, i) => {
      const y = y0 + i * (ROW + 6);
      svg += `<rect class="lane" x="${X0}" y="${y}" width="${X1 - X0}" height="${ROW}"/><text class="ya" x="${X0 - 6}" y="${y + ROW / 2 + 3.5}">${row.label}</text>`;
      if (!m || !fill) return;
      // Each run of minutes alike: how much of the station is in use (to the lane's top: all of it), red when full,
      // and on the cooks' lane, above their work, the cooks free while work waits for room.
      const n = Math.min(fill.n, Math.ceil(T1));
      let k0 = 0;
      for (let k = 1; k <= n; k++) {
        const a = row.used(k0);
        const wa = row.kind === "hands" ? fill.waiting[k0]! : 0;
        if (k < n && Math.abs(row.used(k) - a) < 1e-6 && (row.kind !== "hands" || fill.waiting[k]! === wa)) continue;
        const x = xOf(k0);
        const w = Math.max(0.6, xOf(k) - x);
        const full = fullAt(fill, kit, cooks, row.kind, k0);
        const hh = (Math.min(row.cap, a) / row.cap) * ROW;
        if (hh > 0.2) svg += `<rect class="${full ? "full" : "use"}" x="${x.toFixed(1)}" y="${(y + ROW - hh).toFixed(1)}" width="${w.toFixed(1)}" height="${hh.toFixed(1)}"/>`;
        if (wa > 0) {
          const wh = (Math.min(row.cap, wa) / row.cap) * ROW;
          svg += `<rect class="wait" x="${x.toFixed(1)}" y="${(y + ROW - hh - wh).toFixed(1)}" width="${w.toFixed(1)}" height="${wh.toFixed(1)}"/>`;
        }
        k0 = k;
      }
    });
    // The van.
    const vanR = (m?.late ?? 0) > 0.5 || xOf(DOOR) > X1 - 100;
    svg += `<line class="van" x1="${xOf(DOOR)}" x2="${xOf(DOOR)}" y1="${y0 - 8}" y2="${axisY}"/><text class="vanl${vanR ? " r" : ""}" x="${xOf(DOOR) + (vanR ? -4 : 4)}" y="${y0 - 9}">Van due 10:30</text>`;
    if (m && m.late > 0.5) svg += `<line class="van late" x1="${xOf(m.depart)}" x2="${xOf(m.depart)}" y1="${y0 - 8}" y2="${axisY}"/><text class="vanl late" x="${xOf(m.depart) + 4}" y="${y0 - 9}">${m.finished ? `Leaves ${clockOf(m.depart)}` : "Still waiting"}</text>`;
    if (stale) svg += `<text class="planning" x="${((X0 + X1) / 2).toFixed(1)}" y="${(y0 + axisY) / 2 + 4}">Planning…</text>`;
    svg += `<line class="cur" x1="${X0}" x2="${X0}" y1="${y0 - 6}" y2="${axisY + 2}"/></svg>`;
    chartEl.innerHTML = svg;
    cursor = chartEl.querySelector<SVGLineElement>(".cur");
  }

  /* The ladder: the most meals each purchase gets onto the van with every portion at its best, for this menu. */
  /** The rows: nothing, each purchase alone, and what is bought now when that is more than one thing. */
  const ladderRows = (): Buy[] => {
    const rows = BUYS.slice(0, 4);
    if (bought(settings).length > 1) rows.push({ combi: settings.combi, range: settings.range, cook: settings.cook });
    return rows;
  };
  const rowLabel = (b: Buy) => {
    const ks = bought(b);
    if (!ks.length) return "As built";
    if (ks.length === 1) return `+ ${BUY_NAME[ks[0]!]}`;
    if (ks.length === 3) return "+ all three";
    return `+ ${listOf(ks.map((k) => BUY_WORD[k]))}`;
  };
  function renderLadder(): void {
    if (!ladderEl) return;
    const W = Math.max(280, Math.round(ladderEl.clientWidth));
    const stack = W < 440;
    const rows = ladderRows();
    const rowH = stack ? 32 : 22;
    const barY = stack ? 15 : 4;
    const barH = stack ? 12 : rowH - 8;
    const L = stack ? 2 : Math.min(190, Math.round(W * 0.36));
    const R = W - 62;
    const T = 6;
    const H = T + rows.length * rowH + 20;
    const vals = rows.map((b) => limitOf(settings.menu, b) ?? 0);
    const top = Math.max(settings.meals * 1.1, 400, ...vals);
    const nice = [400, 600, 800, 1000, 1200, 1500, 2000, 2500].find((v) => v >= top) ?? 2500;
    const X = (v: number) => L + (Math.min(v, nice) / nice) * (R - L);
    const base = limitOf(settings.menu, BUYS[0]!);
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg mo-ladder kt-ladder">`;
    for (const v of [0, nice / 2, nice]) svg += `<line class="g v" x1="${X(v)}" x2="${X(v)}" y1="${T}" y2="${T + rows.length * rowH}"/><text class="xa" x="${X(v)}" y="${T + rows.length * rowH + 14}">${v.toLocaleString("en-US")}</text>`;
    // The values go over the order line, on a halo, so a value next to it stays legible.
    let values = "";
    rows.forEach((b, row) => {
      const y = T + row * rowH;
      const v = limitOf(settings.menu, b);
      const cur = sameBuy(b, settings);
      svg += `<g class="row${cur ? " on" : ""}" data-buy="${buyKey(b)}"><rect class="hit" x="0" y="${y}" width="${W}" height="${rowH}"/>`;
      svg += stack ? `<text class="ml st" x="${L}" y="${y + 11}">${rowLabel(b)}</text>` : `<text class="ml" x="${L - 8}" y="${y + rowH / 2 + 4}">${rowLabel(b)}</text>`;
      if (v === undefined) values += `<text class="v pend" x="${X(0) + 5}" y="${y + barY + barH / 2 + 4}">…</text>`;
      else {
        svg += `<rect class="${cur ? "today" : "other"}" x="${X(0)}" y="${y + barY}" width="${Math.max(1, X(v) - X(0))}" height="${barH}"/>`;
        const gain = base !== undefined && row > 0 ? v - base : null;
        values += `<text class="v" x="${X(v) + 5}" y="${y + barY + barH / 2 + 4}">${num(v)}${gain === null ? "" : `<tspan class="${gain > 0 ? "up" : "nil"}"> ${gain > 0 ? `+${num(gain)}` : "+0"}</tspan>`}</text>`;
      }
      svg += `</g>`;
    });
    if (stack) for (let row = 0; row < rows.length; row++) svg += `<line class="order" x1="${X(settings.meals)}" x2="${X(settings.meals)}" y1="${T + row * rowH + barY - 2}" y2="${T + row * rowH + barY + barH + 2}"/>`;
    else svg += `<line class="order" x1="${X(settings.meals)}" x2="${X(settings.meals)}" y1="${T - 4}" y2="${T + rows.length * rowH + 2}"/>`;
    ladderEl.innerHTML = `${svg}<g class="vals">${values}</g></svg>`;
    if (ladderNote) {
      const singles = BUYS.slice(1, 4);
      const all = [BUYS[0]!, ...singles].map((b) => limitOf(settings.menu, b));
      if (all.some((v) => v === undefined)) ladderNote.textContent = "Finding the most meals for each purchase…";
      else {
        const b0 = all[0]!;
        const helps = singles.filter((_, i) => all[i + 1]! > b0);
        const none = singles.filter((_, i) => all[i + 1]! <= b0);
        const bestOf = (xs: Buy[]) => xs.reduce((a, b) => ((limitOf(settings.menu, b) ?? 0) > (limitOf(settings.menu, a) ?? 0) ? b : a));
        let text = `As built, ${num(b0)} meals make the van with every portion at its best.`;
        if (helps.length) {
          const best = bestOf(helps);
          text += ` ${BUY_NAME[bought(best)[0]!].replace(/^./, (c) => c.toUpperCase())} raises that to ${num(limitOf(settings.menu, best)!)}`;
          text += none.length ? `; ${listOf(none.map((b) => BUY_NAME[bought(b)[0]!]))} ${none.length > 1 ? "add" : "adds"} nothing.` : ".";
        } else {
          const pairs = BUYS.slice(4).filter((b) => (limitOf(settings.menu, b) ?? 0) > b0);
          text += ` No one purchase raises it on this menu${pairs.length ? `; ${listOf(bought(bestOf(pairs)).map((k) => BUY_SHORT[k]))} together take it to ${num(limitOf(settings.menu, bestOf(pairs))!)}` : ""}.`;
        }
        ladderNote.textContent = text + " Pick a row to buy it.";
      }
    }
  }
  ladderEl?.addEventListener("click", (e) => {
    const g = (e.target as Element).closest<SVGGElement>("[data-buy]");
    if (!g) return;
    const k = g.dataset.buy ?? "000";
    settings.combi = k[0] === "1";
    settings.range = k[1] === "1";
    settings.cook = k[2] === "1";
    changed(true);
  });

  /* Tiles and status: the plan on screen, against the kitchen as built. */
  /** Over the hour before the van: the share of minutes each station was full, and the cook-minutes spent waiting
   * for room, by the station waited for. */
  /** "waiting for a combi", or "mostly waiting for a combi" when some of the minutes were for another station. */
  const waitingFor = (wait: Record<Station, number>, ranges: number) => {
    const at = (Object.keys(wait) as Station[]).reduce((a, b) => (wait[b] > wait[a] ? b : a), "oven" as Station);
    const total = wait.oven + wait.range + wait.table;
    return `${wait[at] < total - 0.5 ? "mostly " : ""}waiting for ${waitFor(at, ranges)}`;
  };
  const lastHour = () => {
    const m = shown!.m;
    const kit = kitOf(m.settings);
    const k1 = Math.min(fill!.n, Math.ceil(m.depart));
    const k0 = Math.max(0, k1 - 60);
    const share: Record<Kind, number> = { oven: 0, range: 0, table: 0, hands: 0 };
    const wait: Record<Station, number> = { oven: 0, range: 0, table: 0 };
    for (let k = k0; k < k1; k++) {
      for (const kind of Object.keys(share) as Kind[]) if (fullAt(fill!, kit, m.settings.cooks, kind, k)) share[kind] += 1 / Math.max(1, k1 - k0);
      const at = fill!.waitAt[k];
      if (at) wait[at] += fill!.waiting[k]!;
    }
    return { share, wait };
  };
  function renderTiles(): void {
    tilesEl?.toggleAttribute("data-stale", !ready());
    if (!shown) for (const t of [tiles.van, tiles.full, tiles.wait]) setTile(t, "…", "");
    else if (ready()) {
      const { m, built } = shown;
      const any = bought(settings).length > 0;
      // Between the clocks as shown, so the difference reads off them.
      const d = Math.round(m.depart) - Math.round(built.depart);
      // A morning not finished by the cut-off has no departure to show, nor to measure another from.
      if (!m.finished) setTile(tiles.van, "Not yet", `still waiting at ${CUTOFF}, some dishes not packed`);
      else setTile(tiles.van, clockOf(m.depart), m.late > 0.5 ? fmtLate(m.late).toLowerCase() : m.pastBest ? `${num(m.pastBest)} portions past their best` : "on time, every portion at its best", any && d !== 0 && built.finished ? `${d > 0 ? "+" : "−"}${Math.abs(d)} min` : null, d < 0);
      const { share, wait } = lastHour();
      const top = (Object.keys(share) as Kind[]).reduce((a, b) => (share[b] > share[a] + 1e-9 ? b : a), "oven" as Kind);
      if (share[top] < 0.25) setTile(tiles.full, "None", "no station full for long");
      else setTile(tiles.full, kindName(top, kitOf(m.settings).ranges), `${KIND_FULL[top]} ${Math.round(100 * share[top])}% of the hour`);
      const total = wait.oven + wait.range + wait.table;
      if (total < 0.5) setTile(tiles.wait, "None", "no cook stood waiting for room");
      else setTile(tiles.wait, `${num(total)} min`, `all cooks, last hour, ${waitingFor(wait, kitOf(m.settings).ranges)}`);
    }
    const v = limitOf(settings.menu, settings);
    const b0 = limitOf(settings.menu, BUYS[0]!);
    if (v === undefined) setTile(tiles.limit, "…", "");
    else if (!bought(settings).length) setTile(tiles.limit, num(v), "with nothing bought");
    else setTile(tiles.limit, num(v), b0 === undefined ? "" : `as built: ${num(b0)}`, b0 === undefined ? null : v > b0 ? `+${num(v - b0)}` : "+0", v > (b0 ?? v));
  }

  function renderStatus(): void {
    if (!statusEl) return;
    if (failed) {
      statusEl.textContent = failed;
      statusEl.dataset.state = "warn";
      return;
    }
    if (!ready()) {
      statusEl.textContent = "Planning the morning…";
      statusEl.dataset.state = "warn";
      return;
    }
    const { m, built } = shown!;
    const ks = bought(settings);
    const made = madeIt(m);
    const outcome = (x: Morning) => (!x.finished ? `the van is still waiting at ${CUTOFF}, some dishes not yet packed` : x.late > 0.5 ? `the van leaves at ${clockOf(x.depart)}, ${fmtLate(x.late).toLowerCase()}` : x.pastBest ? `the van leaves at 10:30 with ${num(x.pastBest)} portions past their best` : "the van leaves at 10:30 with every portion at its best");
    const { share, wait } = lastHour();
    const top = (Object.keys(share) as Kind[]).reduce((a, b) => (share[b] > share[a] + 1e-9 ? b : a), "oven" as Kind);
    const ranges = kitOf(m.settings).ranges;
    const fullText = share[top] >= 0.25 ? `${kindName(top, ranges).replace(/^The/, "the")} ${kindIs(top, ranges)} ${KIND_FULL[top]} ${Math.round(100 * share[top])}% of the last hour` : "";
    // Minutes summed over the cooks: three cooks free for half the hour stand free 90 minutes between them.
    const total = wait.oven + wait.range + wait.table;
    const waiting = total >= 1 ? `, and the cooks stand free ${num(total)} minutes between them, ${waitingFor(wait, ranges)}` : "";
    let text = `${outcome(m).replace(/^./, (c) => c.toUpperCase())}.`;
    if (!made && fullText) text = `${text.slice(0, -1)}: ${fullText}${waiting}.`;
    if (ks.length && built !== m) {
      const was = outcome(built);
      const same = built.finished && m.finished && madeIt(built) === made && Math.abs(built.depart - m.depart) < 0.5 && built.pastBest === m.pastBest;
      text += !built.finished && !m.finished ? ` As built, the van is still waiting at ${CUTOFF} too.` : same ? ` As built, the same: ${listOf(ks.map((k) => BUY_NAME[k]))} ${ks.length > 1 ? "buy" : "buys"} nothing for this order.` : ` As built, ${was}.`;
    }
    // What would help: the one purchase that raises the most meals most, from the ladder.
    const v0 = limitOf(settings.menu, BUYS[0]!);
    if (!made && v0 !== undefined) {
      const singles = BUYS.slice(1, 4).filter((b) => !sameBuy(b, settings));
      const best = singles.reduce<Buy | null>((a, b) => ((limitOf(settings.menu, b) ?? 0) > (a ? (limitOf(settings.menu, a) ?? 0) : v0) ? b : a), null);
      if (best && !ks.length) text += ` ${BUY_NAME[bought(best)[0]!].replace(/^./, (c) => c.toUpperCase())} gets up to ${num(limitOf(settings.menu, best)!)} meals out on time.`;
      else if (!ks.length) text += " No one purchase raises the most meals on this menu: lower the order, or swap a dish.";
    }
    statusEl.textContent = text;
    statusEl.dataset.state = made ? "ok" : m.late <= 30 ? "warn" : "danger";
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The plans, from the worker.                                                                                 */
  /* ---------------------------------------------------------------------------------------------------------- */

  let worker: Worker | null = null;
  /** Why the plans cannot come, once the worker could not start or has failed. */
  let failed = "";
  function ask(): void {
    if (failed) return;
    if (!worker) {
      try {
        worker = new Worker(new URL("../sim/kitRuns.worker.ts", import.meta.url), { type: "module" });
      } catch {
        failed = "This browser could not start the planner for the kitchen.";
        return;
      }
      worker.onerror = () => {
        failed = "The planner for the kitchen stopped; reload the page to start it again.";
        renderStatus();
      };
      worker.onmessage = (e: MessageEvent<KitDone>) => {
        const d = e.data;
        if (d.kind === "plan") {
          plans.set(d.key, d);
          // The oldest go first, never the plan on screen or the one asked for.
          for (const k of plans.keys()) {
            if (plans.size <= PLANS_KEPT) break;
            if (k !== shownKey && k !== planKey(settings)) plans.delete(k);
          }
          if (d.key === planKey(settings)) apply();
        } else {
          limits.set(d.key, d.meals);
          renderLadder();
          renderTiles();
          renderStatus();
        }
      };
    }
    const key = planKey(settings);
    if (!plans.has(key)) worker.postMessage({ kind: "plan", s: { ...settings } } satisfies KitAsk);
    if (BUYS.some((b) => !limits.has(limitKey(settings.menu, b)))) worker.postMessage({ kind: "limits", menu: settings.menu } satisfies KitAsk);
  }

  let endMin = DOOR + 6;
  let loopS = endMin * SEC_PER_MIN + END_HOLD_S;
  let tMin = 0;
  let loopT = 0;
  /** Put the plan for the settings on screen, if it has come, and everything that shows it. */
  function apply(): void {
    const key = planKey(settings);
    const r = plans.get(key);
    if (r) {
      // Kept as the newest.
      plans.delete(key);
      plans.set(key, r);
      const frac = endMin > 0 ? tMin / endMin : 0;
      shown = r;
      shownKey = key;
      const m = r.m;
      let last = 0;
      for (const j of m.jobs) last = Math.max(last, j.t1);
      // On past the last work, and past the van until its carriers are out of the door.
      endMin = Math.max(DOOR + 6, last + 6, m.depart + ROLL_MIN + 2);
      fill = fillOf(m, endMin);
      kitchen.setKit(kitOf(m.settings));
      loopS = kitchen.plan(m, m.settings.cooks, endMin, endMin * SEC_PER_MIN + END_HOLD_S);
      tMin = frac * endMin;
      loopT = tMin * SEC_PER_MIN;
    }
    renderControls();
    renderChart();
    renderLadder();
    renderTiles();
    renderStatus();
    lastHud = "";
  }
  let pending = 0;
  /** The settings changed: ask for their plan (a slider's steps only once it rests a moment). */
  function changed(now: boolean): void {
    window.clearTimeout(pending);
    const go = () => {
      ask();
      apply();
    };
    if (now) go();
    else {
      renderControls();
      renderLadder();
      pending = window.setTimeout(go, 150);
    }
  }
  for (const b of buyBtns)
    b.addEventListener("click", () => {
      const k = b.dataset.ktBuy as keyof Buy;
      settings[k] = !settings[k];
      changed(true);
    });
  mealsIn?.addEventListener("input", () => {
    settings.meals = Number(mealsIn.value);
    changed(false);
  });

  /* The readouts on the stage, each minute. */
  let lastHud = "";
  function showHud(t: number, a: KitchenMoment): void {
    const m = shown?.m ?? null;
    let note = "";
    let kind = "idle";
    const k = Math.floor(t);
    if (!m || !fill) note = "Planning…";
    else if (t >= m.depart && m.finished) {
      note = `Van gone at ${clockOf(m.depart)}`;
      kind = m.late > 0.5 ? "late" : "held";
    } else if (t >= DOOR) {
      note = "The van is waiting";
      kind = "late";
    } else if (k < fill.n && fill.waiting[k]! > 0) {
      const n = fill.waiting[k]!;
      note = `${n === 1 ? "A cook" : `${n} cooks`} waiting for ${waitFor(fill.waitAt[k]!, kitOf(m.settings).ranges)}`;
      kind = "wait";
    } else {
      const kit = kitOf(m.settings);
      const full = (["oven", "range", "table"] as Kind[]).filter((x) => fullAt(fill!, kit, m.settings.cooks, x, k));
      if (full.length) {
        note = `${listOf(full.map((x) => kindName(x, kit.ranges).replace(/^The /, "")))} full`.replace(/^./, (c) => c.toUpperCase());
        kind = "full";
      }
    }
    const pans = a.hot ? `${a.packed} of ${a.hot} hot pans packed` : "…";
    const sig = `${k}|${note}|${pans}`;
    if (sig === lastHud) return;
    lastHud = sig;
    if (hudClock) hudClock.textContent = clockOf(t);
    if (hudPans) hudPans.textContent = pans;
    if (hudNote) {
      hudNote.textContent = note;
      hudNote.hidden = !note;
      hudNote.dataset.kind = kind;
    }
  }

  /** The pads and badges at minute t, `dt` seconds after the last frame. */
  function showMarks(t: number, dt: number): void {
    const m = shown?.m ?? null;
    const k = Math.floor(t);
    const kit = m ? kitOf(m.settings) : KIT;
    const placed = kitchen.inPlace();
    for (const p of pads) {
      const there = p.kind === "oven" ? (placed.ovens[p.piece] ?? 0) : p.kind === "range" ? (placed.ranges[p.piece] ?? 0) : 1;
      const want = m && fill && there >= 1 && fullAt(fill, kit, m.settings.cooks, p.kind, k) ? 1 : 0;
      p.on += (want - p.on) * Math.min(1, dt * 6);
      if (p.on < 0.01 && want === 0) p.on = 0;
      p.mesh.visible = p.on > 0;
      (p.mesh.material as THREE.MeshBasicMaterial).opacity = 0.62 * p.on;
    }
    // A cook free while work waits for room: their badge, where the routine has them now.
    const waits = m ? m.waits.filter((w) => w.t0 <= t && w.t1 > t) : [];
    kitchen.cooks.forEach((c, i) => {
      const b = badges[i]!;
      const on = c.st.group.visible && !!c.routine && waits.some((w) => w.cooks.includes(kitchen.whoAt(i, loopT)));
      b.visible = on;
      if (on) {
        const at = c.routine!.rootAt(loopT);
        b.position.set(at.x, c.st.dims.standing_height_m + 0.3, at.z);
      }
    });
  }

  renderControls();
  apply();
  ask();
  for (const el of [chartEl, ladderEl]) {
    if (!el) continue;
    let lastW = el.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(el.clientWidth - lastW) < 4) return;
      lastW = el.clientWidth;
      if (el === chartEl) renderChart();
      else renderLadder();
    }).observe(el);
  }
  for (const k of kitchen.cooks) void k.st.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let frozen: number | null = null;
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    if (frozen === null) loopT = (loopT + dt) % loopS;
    else loopT = frozen * SEC_PER_MIN;
    tMin = Math.min(endMin, loopT / SEC_PER_MIN);
    kitchen.pose(loopT);
    const a = kitchen.show(tMin, dt);
    showMarks(tMin, dt);
    if (a.minute || !lastHud) showHud(tMin, a);
    if (cursor) {
      const x = xOf(tMin).toFixed(1);
      cursor.setAttribute("x1", x);
      cursor.setAttribute("x2", x);
    }
    controls.update();
    fitClip();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    kit: {
      info: () => {
        const m = shown?.m ?? null;
        return {
          mountMs: Math.round(mountMs),
          key: planKey(settings),
          ready: ready(),
          depart: m ? clockOf(m.depart) : null,
          late: m ? Math.round(m.late) : null,
          pastBest: m?.pastBest ?? null,
          made: m ? madeIt(m) : null,
          built: shown ? { depart: clockOf(shown.built.depart), late: Math.round(shown.built.late), made: madeIt(shown.built) } : null,
          kit: m ? kitOf(m.settings) : null,
          cooks: kitchen.cooks.filter((k) => k.st.group.visible).length,
          waits: m ? Math.round(m.waits.reduce((s, w) => s + (w.t1 - w.t0) * w.cooks.length, 0)) : null,
          limits: Object.fromEntries(BUYS.map((b) => [buyKey(b), limitOf(settings.menu, b) ?? null])),
          pads: pads.filter((p) => p.mesh.visible).map((p) => `${p.kind}${p.piece}`),
          badges: badges.filter((b) => b.visible).length,
          inPlace: kitchen.inPlace(),
          loopS: +loopS.toFixed(2),
          t: +tMin.toFixed(1),
        };
      },
      hold(min: number | null) {
        frozen = min;
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.position.set(...eye);
        controls.target.set(...at);
        camera.clearViewOffset();
        camera.updateProjectionMatrix();
        controls.update();
      },
      set(s: Partial<KitSettings>) {
        Object.assign(settings, s);
        if (mealsIn) mealsIn.value = String(settings.meals);
        changed(true);
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
