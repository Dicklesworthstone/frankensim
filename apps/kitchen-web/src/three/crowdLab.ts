/**
 * How many cooks fit (chapter 07, second figure). One kitchen, one hood, a crew of one to eight cooks turning out
 * hotel pans: reach-in, prep island, line, pack-out, and back. The runs are computed off the main thread
 * (sim/crewCrowd.ts through sim/crowdRuns.worker.ts): every walk is a leg the crowd planner let go only when it kept
 * clear of everyone. The figure replays the run on screen leg by leg on the cooks' own routines (the same footsteps the
 * planner checked), rings whoever is standing waiting, and marks the floor where the waits piled up; under it, the
 * crew's pans an hour against what the same crew would turn out if nobody were ever in anybody's way.
 */
import * as THREE from "three";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { makeTool, type ToolKind } from "./heroKitchen";
import { buildFixture, type Model, type Piece as FixturePiece } from "./models/all";
import { Routine, makeStaff, type RoutineStep, type Staff, type WorkKind } from "./models/people";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";
import { CREW_COHORTS, DWELL, HOOD, KITCHENS, LINE_Z, MAX_CREW, MORNINGS, RUN_S, START_GAP_S, footprintsOf, legSteps, stanceOf, type CrowdSet, type Kitchen, type Leg, type Piece } from "../sim/crewCrowd";
import type { CrowdAsk, CrowdDone } from "../sim/crowdRuns.worker";

/** Seconds of the run played a second on screen. */
const PLAY_SPEED = 3;
/** The run's clock starts at 9:00. */
const CLOCK0_MIN = 9 * 60;
/** The hood's canopy bottom, metres off the floor. */
const HOOD_Y = 2.03;
const KEYS = ["built", "moved"] as const;
type Key = (typeof KEYS)[number];

const ext = (p: Piece) => ({ w: p.w, d: p.d, h: p.h });
/** The fixtures the figure builds: every real piece of the kitchen and the hood. */
export const CROWD_PIECES: readonly FixturePiece[] = [...KITCHENS.built!.pieces.filter((p) => p.id).map((p): FixturePiece => [p.id, ext(p), 2]), [HOOD.id, HOOD, 1]];

const TOOL_FOR: Partial<Record<WorkKind, ToolKind>> = { chop: "knife", stir: "ladle", saute: "skillet" };
const APRONS = [0x26344d, 0x3a2b22, 0x2c3b2f, 0x4a2430, 0x2a2a2a, 0x3b3424, 0x23343a, 0x402c1c];
/** Chart colours: as built, pulled back, and nobody in the way. */
const LINE_COLOR: Record<Key, string> = { built: "#f59e0b", moved: "#7dd3fc" };

/** Pose a fixture by its front's middle and the way it faces (sim/crewCrowd.ts Piece). */
function poseFixture(g: THREE.Object3D, p: Piece): void {
  const c = Math.cos(p.yaw);
  const s = Math.sin(p.yaw);
  g.rotation.set(0, p.yaw, 0);
  g.position.set(p.fx - (p.w / 2) * c, 0, p.fz + (p.w / 2) * s);
}

/** A cook's replay: their day's legs and each leg's routine, made when first needed. */
interface Play {
  legs: Leg[];
  init: Routine;
  routines: (Routine | null)[];
}

export function mountCrowdLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-crowd-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-crowd-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 12, 24);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.6;
  const lights = studioLights(scene, new THREE.Vector3(0, 0.9, 2.2), 7.5, phone ? 1024 : 2048);
  lights.key.position.set(-2.0, 8.0, 7.0);
  lights.key.intensity = 1.9;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(4.0, 4.0, -3.0);
  lights.rim.intensity = 0.7;
  lights.hemi.intensity = 0.55;

  const K0 = KITCHENS.built!;
  const floor = tilePlane(K0.bounds.maxX - K0.bounds.minX + 0.6, K0.bounds.maxZ - K0.bounds.minZ + 0.8);
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(0, 0, (K0.bounds.minZ + K0.bounds.maxZ) / 2 + 0.2);
  scene.add(floor);

  // Every piece once; the island's tables move with the layout.
  const models: (Model | null)[] = K0.pieces.map((p) => {
    if (!p.id) return null;
    const m = buildFixture(p.id, ext(p), 2);
    poseFixture(m.group, p);
    scene.add(m.group);
    return m;
  });
  models.forEach((m) => m?.hooks.flames?.setAll?.(0));
  /** Which pieces have a cook working at them, and which ranges' burners are on. */
  const lit: boolean[] = models.map(() => false);
  const burning: boolean[] = models.map(() => false);
  function setBurners(): void {
    models.forEach((m, i) => {
      const f = m?.hooks.flames;
      if (!f || lit[i] === burning[i]) return;
      burning[i] = lit[i]!;
      f.setFlow(1, lit[i] ? 0.8 : 0);
      f.setFlow(2, lit[i] ? 0.55 : 0);
    });
  }
  const hood = buildFixture(HOOD.id, HOOD, 1);
  hood.group.position.set(-HOOD.w / 2, HOOD_Y, HOOD.d);
  scene.add(hood.group);

  /*
   * Whoever stands waiting: a ring on the floor round their feet and a pause badge over their head (the badge shows
   * when the island or another cook hides their feet). Red, waiting on someone in the way; amber, for a free place.
   */
  const ringGeo = new THREE.RingGeometry(0.3, 0.38, 40).rotateX(-Math.PI / 2);
  const RED = "#ef4444";
  const AMBER = "#fbbf24";
  const ringMat = (c: string) => new THREE.MeshBasicMaterial({ color: c, toneMapped: false, transparent: true, opacity: 0.95, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -3, polygonOffsetUnits: -3 });
  const badgeMat = (c: string) => {
    const cv = document.createElement("canvas");
    cv.width = cv.height = 96;
    const x = cv.getContext("2d")!;
    x.fillStyle = c;
    x.strokeStyle = "rgba(255,255,255,0.92)";
    x.lineWidth = 6;
    x.beginPath();
    x.arc(48, 48, 40, 0, Math.PI * 2);
    x.fill();
    x.stroke();
    x.fillStyle = "#fff";
    x.fillRect(34, 30, 9, 36);
    x.fillRect(53, 30, 9, 36);
    const t = new THREE.CanvasTexture(cv);
    t.colorSpace = THREE.SRGBColorSpace;
    return new THREE.SpriteMaterial({ map: t, toneMapped: false, depthTest: false, depthWrite: false });
  };
  const MARK = { body: { ring: ringMat(RED), badge: badgeMat(RED) }, place: { ring: ringMat(AMBER), badge: badgeMat(AMBER) } };
  const marks: { ring: THREE.Mesh; badge: THREE.Sprite }[] = [];
  const hideMarks = () =>
    marks.forEach((m) => {
      m.ring.visible = false;
      m.badge.visible = false;
    });
  /* The floor where the waits piled up: a soft glow for each wait, wider for longer. */
  const MAX_HEAT = 400;
  /** Waits closer than this (m) are one place on the floor. */
  const HEAT_CELL = 0.2;
  const glow = (() => {
    const cv = document.createElement("canvas");
    cv.width = cv.height = 64;
    const x = cv.getContext("2d")!;
    const g = x.createRadialGradient(32, 32, 0, 32, 32, 32);
    g.addColorStop(0, "rgba(255,255,255,1)");
    g.addColorStop(0.45, "rgba(255,255,255,0.55)");
    g.addColorStop(1, "rgba(255,255,255,0)");
    x.fillStyle = g;
    x.fillRect(0, 0, 64, 64);
    return new THREE.CanvasTexture(cv);
  })();
  const heatMat = new THREE.MeshBasicMaterial({ color: 0xff5c8a, map: glow, toneMapped: false, transparent: true, opacity: 0.78, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -2 });
  const heat = new THREE.InstancedMesh(new THREE.PlaneGeometry(2, 2).rotateX(-Math.PI / 2), heatMat, MAX_HEAT);
  heat.count = 0;
  heat.frustumCulled = false;
  heat.renderOrder = 1;
  scene.add(heat);

  /* The crew, made as the slider asks for them. */
  const staff: Staff[] = [];
  const tools: Map<ToolKind, THREE.Object3D>[] = [];
  const held: (ToolKind | null)[] = [];
  function cookAt(i: number): Staff {
    while (staff.length <= i) {
      const j = staff.length;
      const st = makeStaff(CREW_COHORTS[j % CREW_COHORTS.length]!, phone ? 1 : 2, { seed: 7 + j * 13, apron: APRONS[j % APRONS.length], sleeves: "cuffed", shadows: !phone });
      scene.add(st.group);
      staff.push(st);
      tools.push(new Map());
      held.push(null);
      const ring = new THREE.Mesh(ringGeo, MARK.body.ring);
      ring.renderOrder = 2;
      ring.visible = false;
      const badge = new THREE.Sprite(MARK.body.badge);
      badge.scale.setScalar(0.3);
      badge.renderOrder = 5;
      badge.visible = false;
      scene.add(ring, badge);
      marks.push({ ring, badge });
      void st.ready.then(() => void gate.refresh(scene));
    }
    return staff[i]!;
  }

  /* Camera: from in front of the island, the whole room across the view. */
  const camera = new THREE.PerspectiveCamera(30, 16 / 9, 0.05, 80);
  const TARGET = new THREE.Vector3(0, 0.7, 2.0);
  const VIEW = new THREE.Vector3(0, phone ? 1.9 : 1.45, 1).normalize();
  const B = K0.bounds;
  const keep = [B.minX, B.maxX].flatMap((x) => [new THREE.Vector3(x, 0, B.maxZ), new THREE.Vector3(x, 0, B.minZ + 0.3), new THREE.Vector3(x * 0.82, HOOD_Y + HOOD.h, 0.2)]);
  // The reach-ins' tops against the left wall.
  keep.push(new THREE.Vector3(B.minX, 2.06, 1.95), new THREE.Vector3(B.minX, 2.06, 1.95 + 2 * 1.372 + 0.04));
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 2, maxDistance: 22, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const hudEls = Array.from(stageEl.querySelectorAll<HTMLElement>(".cr-hud"));
  const legendEl = stageEl.querySelector<HTMLElement>(".cr-legend");
  const ndc = new THREE.Vector3();
  const frame = () => {
    if (userMoved) return;
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    // The room fits between the HUDs and the legend, then the view shifts to sit in that band.
    const hudB = hudEls.reduce((m, e) => Math.max(m, e.offsetTop + e.offsetHeight), 0);
    const top = Math.min(0.9, 1 - (2 * (hudB + 6)) / h);
    const bot = -Math.min(0.92, 1 - (2 * ((legendEl ? h - legendEl.offsetTop : 0) + 4)) / h);
    let dist = 3;
    let lo = 0;
    let hi = 0;
    for (; dist < 30; dist += 0.08) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      lo = Infinity;
      hi = -Infinity;
      let wide = false;
      for (const p of keep) {
        ndc.copy(p).project(camera);
        lo = Math.min(lo, ndc.y);
        hi = Math.max(hi, ndc.y);
        if (Math.abs(ndc.x) >= 0.97) wide = true;
      }
      if (!wide && hi - lo < top - bot) break;
    }
    const shift = (top + bot) / 2 - (hi + lo) / 2;
    if (Math.abs(shift) > 0.005) {
      camera.setViewOffset(w, h, 0, (shift * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.2;
    fog.far = dist * 2.8;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);

  /* The page. */
  const q = <T extends Element = HTMLElement>(sel: string) => root.querySelector<T>(sel);
  const keyBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-cr-key]"));
  const crewIn = q<HTMLInputElement>("[data-cr-crew]");
  const crewVal = q("[data-cr-crew-val]");
  const hudName = q("[data-cr-hud-name]");
  const hudClock = q("[data-cr-hud-clock]");
  const hudWait = q("[data-cr-hud-wait]");
  const chartEl = q("[data-crowd-chart]");
  const statusEl = q("[data-cr-status]");
  const tile = (k: string) => ({ v: q(`[data-cr-${k}]`), d: q(`[data-cr-${k}-d]`), s: q(`[data-cr-${k}-s]`) });
  const tiles = { pans: tile("pans"), lost: tile("lost"), wait: tile("wait"), other: tile("other") };
  const num = (v: number) => Math.round(v).toLocaleString("en-US");

  let key: Key = "built";
  let crew = MAX_CREW;
  const runs = new Map<string, CrowdSet>();
  const runKey = (k: string, n: number, ghosts: boolean) => `${k}/${n}/${ghosts ? "g" : "b"}`;
  const runOf = (k: string, n: number, ghosts: boolean) => runs.get(runKey(k, n, ghosts)) ?? null;

  /* The runs: the one on screen first, then the rest of the chart. */
  let worker: Worker | null = null;
  function ask(): void {
    if (!worker) {
      try {
        worker = new Worker(new URL("../sim/crowdRuns.worker.ts", import.meta.url), { type: "module" });
      } catch {
        if (statusEl) statusEl.textContent = "This browser can't run the crew off the page.";
        return;
      }
      worker.onmessage = (e: MessageEvent<CrowdDone>) => {
        const d = e.data;
        runs.set(runKey(d.kitchen, d.crew, d.ghosts), d.set);
        if (d.kitchen === key && d.crew === crew && !d.ghosts) play();
        renderPage();
      };
    }
    const list: CrowdAsk[] = [];
    const push = (k: string, n: number, ghosts: boolean) => {
      if (!runs.has(runKey(k, n, ghosts)) && !list.some((a) => a.kitchen === k && a.crew === n && a.ghosts === ghosts)) list.push({ kitchen: k, crew: n, ghosts });
    };
    const other: Key = key === "built" ? "moved" : "built";
    push(key, crew, false);
    push(key, crew, true);
    push(other, crew, false);
    for (let n = MAX_CREW; n >= 1; n--) {
      push(key, n, false);
      push(key, n, true);
      push(other, n, false);
      push(other, n, true);
    }
    worker.postMessage(list);
  }

  /* Replay. */
  let plays: Play[] = [];
  let playing: CrowdSet | null = null;
  let t0 = 0;
  const atVec = (a: [number, number, number]) => new THREE.Vector3(a[0], a[1], a[2]);
  function play(): void {
    const set = runOf(key, crew, false);
    if (!set || set === playing) return;
    playing = set;
    const run = set.shown;
    plays = run.cooks.map((day, i) => {
      const st = cookAt(i);
      const legs = day.legs;
      const first = legs[0];
      const startS = i * START_GAP_S;
      const { from, at, surfaceY } = prepStance(i);
      // Before the first leg: standing at their prep place, then prepping the pan they started with, then waiting.
      const steps: RoutineStep[] = [];
      if (startS > 0) steps.push({ kind: "stand", duration: startS });
      steps.push({ kind: "work", work: "chop", duration: DWELL.prep, opts: { at, surfaceY } });
      const gap = (first ? first.go : RUN_S) - startS - DWELL.prep;
      if (gap > 0.01) steps.push({ kind: "stand", duration: gap });
      return { legs, init: new Routine(st, from, steps, false), routines: legs.map(() => null) };
    });
    // The floor under the waits, all three mornings': the waits gathered by where they stood (cooks wait at the same
    // few places), a patch for each place, wider for more minutes.
    const spots = new Map<string, { x: number; z: number; s: number }>();
    for (const [x, z, s] of set.waits) {
      const key = `${Math.round(x / HEAT_CELL)},${Math.round(z / HEAT_CELL)}`;
      const h = spots.get(key);
      if (h) {
        h.x += x * s;
        h.z += z * s;
        h.s += s;
      } else spots.set(key, { x: x * s, z: z * s, s });
    }
    const m4 = new THREE.Matrix4();
    const places = [...spots.values()].sort((a, b) => b.s - a.s).slice(0, MAX_HEAT);
    places.forEach((h, k) => {
      const r = Math.min(0.9, 0.2 + 0.04 * Math.sqrt(h.s / MORNINGS.length));
      m4.makeScale(r, 1, r).setPosition(h.x / h.s, 0.003, h.z / h.s);
      heat.setMatrixAt(k, m4);
    });
    heat.count = places.length;
    heat.instanceMatrix.needsUpdate = true;
    t0 = clock;
    staff.forEach((s, i) => {
      s.group.visible = i < crew;
      marks[i]!.ring.visible = false;
      marks[i]!.badge.visible = false;
    });
  }
  /** Where cook i stands at their own prep place, the planner's stance for their body. */
  function prepStance(i: number): { from: { x: number; z: number; yaw: number }; at: THREE.Vector3; surfaceY: number } {
    const kit = KITCHENS[key]!;
    const sp = kit.spots[kit.prepOrder[i]!]!;
    const g = stanceOf(kit, sp, footprintsOf(cookAt(i)).rWalk);
    return { from: { x: g.x, z: g.z, yaw: g.yaw }, at: atVec(g.at), surfaceY: sp.surfaceY };
  }
  function routineOf(i: number, j: number): Routine {
    const p = plays[i]!;
    let r = p.routines[j];
    if (!r) {
      const leg = p.legs[j]!;
      const next = p.legs[j + 1];
      const steps = legSteps(leg, atVec(leg.at));
      const gap = (next ? next.go : RUN_S) - (leg.arrive + leg.dwell);
      if (gap > 0.01) steps.push({ kind: "stand", duration: gap });
      r = p.routines[j] = new Routine(staff[i]!, leg.from, steps, false);
    }
    return r;
  }

  /** Pose every cook at run time t; returns how many stand waiting on someone and how many for a free place. */
  function show(t: number): { body: number; place: number } {
    let body = 0;
    let place = 0;
    lit.fill(false);
    const kit = KITCHENS[key]!;
    plays.forEach((p, i) => {
      if (i >= crew) return;
      const st = staff[i]!;
      // The leg under way: the last to have set off by t.
      let lo = 0;
      let hi = p.legs.length;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (p.legs[mid]!.go <= t) lo = mid + 1;
        else hi = mid;
      }
      const j = lo - 1;
      const r = j < 0 ? p.init : routineOf(i, j);
      const lt = j < 0 ? t : t - p.legs[j]!.go;
      const s = r.apply(lt);
      // The tool for the work in hand.
      const want = s.activity === "work" && s.work ? (TOOL_FOR[s.work] ?? null) : null;
      if (want !== held[i]) {
        if (held[i]) tools[i]!.get(held[i]!)?.removeFromParent();
        held[i] = want;
        if (want) {
          let tool = tools[i]!.get(want);
          if (!tool) tools[i]!.set(want, (tool = makeTool(want)));
          st.hand("R").add(tool);
          tool.position.set(0, 0, 0);
          tool.quaternion.identity();
        }
      }
      // A range's burners are lit while a cook works at it.
      const cur = j >= 0 ? p.legs[j]! : null;
      if (cur && cur.kind === "hot" && s.activity === "work") lit[kit.spots[cur.spot]!.piece] = true;
      // Waiting: ready for the next leg and not yet gone.
      const next = p.legs[j + 1];
      const waiting = next && t >= next.ready && t < next.go ? (next.blocked ? "body" : "place") : null;
      const mk = marks[i]!;
      mk.ring.visible = mk.badge.visible = waiting !== null;
      if (waiting) {
        const root = r.rootAt(lt);
        mk.ring.position.set(root.x, 0.006, root.z);
        mk.badge.position.set(root.x, st.dims.standing_height_m + 0.32, root.z);
        mk.ring.material = MARK[waiting].ring;
        mk.badge.material = MARK[waiting].badge;
        if (waiting === "body") body++;
        else place++;
      }
    });
    return { body, place };
  }

  /* Chart: pans an hour against cooks, both layouts and nobody in the way. */
  function renderChart(): void {
    if (!chartEl) return;
    const W = Math.max(280, Math.round(chartEl.clientWidth - 2 * 12));
    const H = W < 440 ? 190 : 210;
    const X0 = 44;
    const X1 = W - 14;
    const Y0 = 14;
    const Y1 = H - 30;
    const ghostTop = Math.max(...KEYS.flatMap((k) => Array.from({ length: MAX_CREW }, (_, i) => runOf(k, i + 1, true)?.perHour ?? 0)), 26.4 * MAX_CREW);
    const yMax = Math.ceil((ghostTop * 1.05) / 50) * 50;
    const xOf = (n: number) => X0 + ((n - 1) / (MAX_CREW - 1)) * (X1 - X0);
    const yOf = (v: number) => Y1 - (v / yMax) * (Y1 - Y0);
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg cr-svg" aria-hidden="true">`;
    for (let v = 0; v <= yMax; v += 50) svg += `<line class="g" x1="${X0}" x2="${X1}" y1="${yOf(v).toFixed(1)}" y2="${yOf(v).toFixed(1)}"/><text class="ya" x="${X0 - 6}" y="${(yOf(v) + 3.5).toFixed(1)}">${v}</text>`;
    for (let n = 1; n <= MAX_CREW; n++) svg += `<text class="xa" x="${xOf(n).toFixed(1)}" y="${Y1 + 15}">${n}</text>`;
    svg += `<text class="xl" x="${X1}" y="${H - 2}">cooks</text><text class="yl" x="${X0}" y="${Y0 - 4}">pans an hour</text>`;
    svg += `<line class="sel" x1="${xOf(crew).toFixed(1)}" x2="${xOf(crew).toFixed(1)}" y1="${Y0}" y2="${Y1}"/>`;
    const path = (k: Key, ghosts: boolean) => {
      let d = "";
      for (let n = 1; n <= MAX_CREW; n++) {
        const r = runOf(k, n, ghosts);
        if (!r) continue;
        d += `${d ? "L" : "M"}${xOf(n).toFixed(1)},${yOf(r.perHour).toFixed(1)}`;
      }
      return d;
    };
    const gd = path(key, true);
    if (gd) svg += `<path class="ghost" d="${gd}"/>`;
    for (const k of KEYS) {
      const d = path(k, false);
      if (d) svg += `<path class="run${k === key ? " on" : ""}" d="${d}" stroke="${LINE_COLOR[k]}"/>`;
      for (let n = 1; n <= MAX_CREW; n++) {
        const r = runOf(k, n, false);
        if (r) svg += `<circle cx="${xOf(n).toFixed(1)}" cy="${yOf(r.perHour).toFixed(1)}" r="${k === key && n === crew ? 4.5 : 2.6}" fill="${LINE_COLOR[k]}"/>`;
      }
    }
    svg += "</svg>";
    chartEl.innerHTML = svg;
  }

  function setTile(t: { v: HTMLElement | null; d: HTMLElement | null; s: HTMLElement | null }, value: string, sub: string, delta: string | null = null, good = true): void {
    if (t.v) t.v.textContent = value;
    if (t.s) t.s.textContent = sub;
    if (t.d) {
      t.d.hidden = delta === null;
      t.d.textContent = delta ?? "";
      t.d.dataset.good = String(good);
    }
  }

  /** The share of waiting done in the hot aisle (between the line and the island). */
  function hotAisleShare(run: Pick<CrowdSet, "waits">, k: Kitchen): number {
    const islandFront = Math.min(...k.pieces.filter((p) => p.id && p.yaw === Math.PI).map((p) => p.fz));
    let all = 0;
    let hot = 0;
    for (const [x, z, s] of run.waits) {
      all += s;
      if (z > LINE_Z && z < islandFront && Math.abs(x) < 3.2) hot += s;
    }
    return all > 0 ? hot / all : 0;
  }

  function renderPage(): void {
    keyBtns.forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.crKey === key)));
    if (crewVal) crewVal.textContent = `${crew} cook${crew === 1 ? "" : "s"}`;
    if (hudName) hudName.textContent = `${KITCHENS[key]!.name} · ${crew} cook${crew === 1 ? "" : "s"}`;
    renderChart();
    const b = runOf(key, crew, false);
    const g = runOf(key, crew, true);
    const other: Key = key === "built" ? "moved" : "built";
    const o = runOf(other, crew, false);
    if (!b || !g) {
      setTile(tiles.pans, "…", "running the crew");
      setTile(tiles.lost, "…", "");
      setTile(tiles.wait, "…", "");
      setTile(tiles.other, "…", "");
      if (statusEl) {
        statusEl.textContent = `Running ${crew} cook${crew === 1 ? "" : "s"} through three mornings…`;
        statusEl.dataset.state = "warn";
      }
      return;
    }
    const lost = g.perHour > 0 ? Math.max(0, 1 - b.perHour / g.perHour) : 0;
    const perCook = g.perHour / crew;
    const equiv = perCook > 0 ? b.perHour / perCook : crew;
    const waitMin = b.waitPerHour / 60;
    setTile(tiles.pans, num(b.perHour), `${num(g.perHour)} with nobody in the way`);
    setTile(tiles.lost, `${Math.round(lost * 100)}%`, `${crew} cook${crew === 1 ? "" : "s"} turn out what ${equiv.toFixed(1)} would`, null);
    setTile(tiles.wait, `${waitMin.toFixed(1)} min`, "an hour, on someone in the way");
    if (o) {
      const dv = o.perHour - b.perHour;
      setTile(tiles.other, num(o.perHour), `pans an hour ${other === "moved" ? "with the island pulled back" : "as built"}`, Math.abs(dv) >= 0.5 ? `${dv > 0 ? "+" : "−"}${num(Math.abs(dv))}` : null, dv >= 0);
    } else setTile(tiles.other, "…", other === "moved" ? "with the island pulled back" : "as built");
    if (statusEl) {
      const share = hotAisleShare(b, KITCHENS[key]!);
      const where = share > 0.5 ? `, and ${Math.round(share * 100)}% of that waiting is in the hot aisle` : "";
      const lostTxt = lost < 0.005 ? "They lose nothing to each other" : `They lose ${Math.round(lost * 100)}% waiting on each other${where}`;
      let then = "";
      if (o) {
        const og = runOf(other, crew, true);
        const ol = og && og.perHour > 0 ? Math.max(0, 1 - o.perHour / og.perHour) : null;
        if (ol !== null) then = other === "moved" ? ` Pull the island back 18 in, leaving the hood where it is, and they lose ${Math.round(ol * 100)}%.` : ` As built, with the island 42 in from the line, they lose ${Math.round(ol * 100)}%.`;
      }
      statusEl.textContent = `${crew} cook${crew === 1 ? "" : "s"} turn out ${num(b.perHour)} pans an hour; with nobody ever in the way they would turn out ${num(g.perHour)}. ${lostTxt}.${then}`;
      statusEl.dataset.state = lost > 0.08 ? "warn" : "ok";
    }
  }

  function setKitchen(k: Key): void {
    key = k;
    const kit = KITCHENS[k]!;
    kit.pieces.forEach((p, i) => {
      const m = models[i];
      if (m) poseFixture(m.group, p);
    });
    playing = null;
    plays = [];
    heat.count = 0;
    staff.forEach((s) => (s.group.visible = false));
    hideMarks();
    showStill();
    ask();
    play();
    renderPage();
  }
  /** Before a run is in: the crew standing at their prep places. */
  function showStill(): void {
    for (let i = 0; i < crew; i++) {
      cookAt(i).group.visible = true;
      new Routine(staff[i]!, prepStance(i).from, [{ kind: "stand", duration: 1 }], false).apply(0);
    }
    for (let i = crew; i < staff.length; i++) staff[i]!.group.visible = false;
  }
  for (const btn of keyBtns)
    btn.addEventListener("click", () => {
      const k = btn.dataset.crKey === "moved" ? "moved" : "built";
      if (k !== key) setKitchen(k);
    });
  crewIn?.addEventListener("input", () => {
    const n = Math.max(1, Math.min(MAX_CREW, Math.round(Number(crewIn.value) || MAX_CREW)));
    if (n === crew) return;
    crew = n;
    playing = null;
    plays = [];
    heat.count = 0;
    hideMarks();
    showStill();
    ask();
    play();
    renderPage();
  });

  showStill();
  renderPage();
  frame();
  const ro = new ResizeObserver(() => {
    frame();
    renderChart();
  });
  // The HUD and the legend too: the HUD wraps to a second line once its text is in, and the view must clear it.
  for (const el of [stageEl, ...hudEls, legendEl]) if (el) ro.observe(el);
  const mountMs = performance.now() - mountT0;

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
    if (v && !worker) ask();
  });
  let clock = 0;
  let frozen: number | null = null;
  let last = performance.now();
  let lastHud = "";
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    clock += dt;
    const run = playing?.shown ?? null;
    const t = frozen ?? (run ? ((clock - t0) * PLAY_SPEED) % run.seconds : 0);
    const w = run ? show(t) : { body: 0, place: 0 };
    if (!run) lit.fill(false);
    setBurners();
    const m = CLOCK0_MIN + Math.floor(t / 60);
    const hud = run ? `${Math.floor(m / 60)}:${String(m % 60).padStart(2, "0")} · ${PLAY_SPEED}× speed` : "…";
    const waitTxt = run ? (w.body || w.place ? `${w.body} waiting on someone${w.place ? ` · ${w.place} for a free place` : ""}` : "Nobody waiting") : "…";
    const line = hud + waitTxt;
    if (line !== lastHud) {
      lastHud = line;
      if (hudClock) hudClock.textContent = hud;
      if (hudWait) {
        hudWait.textContent = waitTxt;
        hudWait.dataset.kind = w.body ? "ot" : w.place ? "wait" : "";
      }
    }
    models.forEach((m) => m?.hooks.tick?.(clock, dt, 0.5));
    controls.update();
    fitClip();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    crowd: {
      info: () => ({
        mountMs: Math.round(mountMs),
        tris: renderer.info.render.triangles,
        calls: renderer.info.render.calls,
        key,
        crew,
        runs: [...runs.entries()].map(([k, r]) => ({ k, perHour: +r.perHour.toFixed(1), minGap: r.minGap === null ? null : +r.minGap.toFixed(3), gridlock: r.gridlock, waitPerHour: Math.round(r.waitPerHour) })),
        playing: playing ? `${playing.kitchen}/${playing.crew}` : null,
      }),
      hold(t: number | null) {
        frozen = t;
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.clearViewOffset();
        camera.position.set(...eye);
        controls.target.set(...at);
        controls.update();
      },
    },
  });
}
