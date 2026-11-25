/**
 * The floor: where the commissary's blast chiller should stand.
 *
 * The room is the commissary exemplar, built from its own manifest (every fixture at its catalogue size and its
 * place in the plan). Along the walls a strip shows, for every spot a chiller fits, how far the week's cooked-ahead
 * pans would be carried to reach it (commissaryFloor: distance fields over the free floor). The routes drawn are A*
 * paths with the chiller in place, as wide as their traffic, and a cook walks them with a full hotel pan, at the pace
 * the walking time is quoted at. The cook's motion is a kinematic illustration: posed for reach and stride.
 */
import * as THREE from "three";
import { mergeGeometries } from "three/examples/jsm/utils/BufferGeometryUtils.js";
import {
  CARRY_MPS,
  WALKIN,
  buildFloor,
  fieldsOf,
  placeChiller,
  routesFor,
  trafficOf,
  wallSpots,
  type FloorAnswer,
  type Fields,
  type Floor,
  type FloorManifest,
  type FloorMeta,
  type Route,
  type Spot,
  type Traffic,
} from "../sim/commissaryFloor";
import { menuFor, runCommissary } from "../sim/commissaryWeek";
import { ASSET_BASE, bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { buildFixture } from "./models/all";
import { palette } from "./models/materials";
import { Routine, makeStaff, panSeconds, type RoutineStep, type Staff } from "./models/people";
import { placeFixture } from "./models/shop";
import { studioEnvironment, studioLights } from "./models/stage";

const ID = "commissary_catering_3200";
const SMALL = /KNIFE|BOARD|HOTEL-PAN|GRATER|SCALE|SLICER|PROCESSOR|BOOSTER|HOOD/;
const HOT = /KETTLE|SKILLET|COMBI|BLAST/;
const CHILLER_ID = "NSF-BLAST-CHILLER-15PAN";
const PAN_EXT = { w: 0.53, d: 0.325, h: 0.1 };
const SWAPS = ["chicken", "pasta", "salad"];
/** Where a pan is taken from or set down at each kind of fixture: the height, and how far in from the fixture's front. */
const RAIL = { combi: { y: 0.93, inset: 0.32 }, kettle: { y: 0.9, inset: 0.32 }, skillet: { y: 0.9, inset: 0.31 }, prep: { y: 0.914, inset: 0.3 }, chiller: { y: 0.95, inset: 0.3 } };

function tilePlane(w: number, h: number): THREE.Mesh {
  const g = new THREE.PlaneGeometry(w, h);
  const uv = g.attributes.uv as THREE.BufferAttribute;
  for (let i = 0; i < uv.count; i++) uv.setXY(i, (uv.getX(i) * w) / 0.608, (uv.getY(i) * h) / 0.608);
  const m = new THREE.Mesh(g, palette().tile);
  m.receiveShadow = true;
  return m;
}

/** A dash that runs along a ribbon (u along its length), for traffic that flows one way. */
function dashTexture(): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = 64;
  c.height = 4;
  const x = c.getContext("2d")!;
  const g = x.createLinearGradient(0, 0, 64, 0);
  g.addColorStop(0, "rgba(255,255,255,0.25)");
  g.addColorStop(0.55, "rgba(255,255,255,1)");
  g.addColorStop(0.62, "rgba(255,255,255,0.25)");
  g.addColorStop(1, "rgba(255,255,255,0.25)");
  x.fillStyle = g;
  x.fillRect(0, 0, 64, 4);
  const t = new THREE.CanvasTexture(c);
  t.wrapS = THREE.RepeatWrapping;
  return t;
}

/** A flat ribbon on the floor along a path: u is metres / `period` along it. */
function ribbon(path: readonly (readonly [number, number])[], width: number, y: number, period: number): THREE.BufferGeometry {
  const pos: number[] = [];
  const uv: number[] = [];
  const idx: number[] = [];
  let s = 0;
  for (let i = 0; i < path.length; i++) {
    const p = path[i]!;
    const a = path[Math.max(0, i - 1)]!;
    const b = path[Math.min(path.length - 1, i + 1)]!;
    let tx = b[0] - a[0];
    let tz = b[1] - a[1];
    const tl = Math.hypot(tx, tz) || 1;
    tx /= tl;
    tz /= tl;
    if (i > 0) s += Math.hypot(p[0] - path[i - 1]![0], p[1] - path[i - 1]![1]);
    const nx = -tz * (width / 2);
    const nz = tx * (width / 2);
    pos.push(p[0] + nx, y, p[1] + nz, p[0] - nx, y, p[1] - nz);
    uv.push(s / period, 0, s / period, 1);
    if (i > 0) {
      const k = 2 * i;
      idx.push(k - 2, k - 1, k, k - 1, k + 1, k);
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  g.setIndex(idx);
  return g;
}

/** Every plain mesh under the groups, merged into one mesh per material in world space; the rest left in place. */
function mergeStatic(groups: THREE.Object3D[], into: THREE.Object3D): void {
  // One bin per material and shadow setting, so a mesh that cast no shadow (a lamp, a glass door) still casts none.
  const bins = new Map<string, { mat: THREE.Material; cast: boolean; receive: boolean; geos: THREE.BufferGeometry[] }>();
  const taken: THREE.Mesh[] = [];
  for (const g of groups) {
    g.updateMatrixWorld(true);
    g.traverse((o) => {
      const m = o as THREE.Mesh;
      if (!m.isMesh || (m as THREE.InstancedMesh).isInstancedMesh || (m as THREE.SkinnedMesh).isSkinnedMesh || Array.isArray(m.material)) return;
      if ((m.material as THREE.ShaderMaterial).isShaderMaterial) return;
      for (let q: THREE.Object3D | null = m; q; q = q.parent) if (!q.visible) return;
      const src = m.geometry as THREE.BufferGeometry;
      if (!src.getAttribute("position")) return;
      const geo = src.clone().applyMatrix4(m.matrixWorld);
      for (const name of Object.keys(geo.attributes)) if (!["position", "normal", "uv"].includes(name)) geo.deleteAttribute(name);
      geo.morphAttributes = {};
      if (!geo.getAttribute("normal")) geo.computeVertexNormals();
      if (!geo.getAttribute("uv")) geo.setAttribute("uv", new THREE.BufferAttribute(new Float32Array(geo.getAttribute("position").count * 2), 2));
      if (!geo.index) geo.setIndex([...Array(geo.getAttribute("position").count).keys()]);
      geo.clearGroups();
      const mat = m.material as THREE.Material;
      const key = `${mat.uuid}|${m.castShadow ? 1 : 0}${m.receiveShadow ? 1 : 0}`;
      if (!bins.has(key)) bins.set(key, { mat, cast: m.castShadow, receive: m.receiveShadow, geos: [] });
      bins.get(key)!.geos.push(geo);
      taken.push(m);
    });
  }
  for (const m of taken) m.parent?.remove(m);
  for (const { mat, cast, receive, geos } of bins.values()) {
    const g = geos.length === 1 ? geos[0]! : mergeGeometries(geos, false);
    if (!g) continue;
    const mesh = new THREE.Mesh(g, mat);
    mesh.castShadow = cast;
    mesh.receiveShadow = receive;
    mesh.matrixAutoUpdate = false;
    into.add(mesh);
  }
  for (const g of groups) into.add(g);
}

const fwd = (yaw: number): [number, number] => [Math.sin(yaw), Math.cos(yaw)];

export function mountFloorLab(root: HTMLElement): void {
  const stageEl = root.querySelector<HTMLElement>("[data-floor-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-floor-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  const phone = window.matchMedia("(max-width: 720px)").matches;
  const statusEl = root.querySelector<HTMLElement>("[data-fl-status]");
  const base = `${ASSET_BASE}/exemplars/${ID}/${ID}`;
  if (statusEl) statusEl.textContent = "Building the commissary…";
  void Promise.all([fetch(`${base}.scene3d.v1.json`).then((r) => (r.ok ? r.json() : null)), fetch(`${base}.site_metadata.json`).then((r) => (r.ok ? r.json() : null))])
    .then(([manifest, meta]) => {
      if (!manifest) throw new Error("manifest unavailable");
      start(root, stageEl, canvas, renderer, manifest as FloorManifest & { fixtures: (FloorManifest["fixtures"][number] & { local_extents_m: [number, number, number] })[] }, meta as FloorMeta | null, phone);
    })
    .catch(() => {
      if (statusEl) {
        statusEl.textContent = "The commissary's plan could not be loaded.";
        statusEl.dataset.state = "warn";
      }
    });
}

function start(root: HTMLElement, stageEl: HTMLElement, canvas: HTMLCanvasElement, renderer: THREE.WebGLRenderer, manifest: FloorManifest, meta: FloorMeta | null, phone: boolean): void {
  const mountT0 = performance.now();
  /* The model: the floor, its distance fields and the spots on its walls. */
  const floor: Floor = buildFloor(manifest, meta);
  const fields: Fields = fieldsOf(floor);
  const spots: Spot[] = wallSpots(floor);
  const traffic = new Map<string, Traffic>();
  const trafficFor = (swapped: boolean): Traffic => {
    const k = swapped ? "swapped" : "written";
    if (!traffic.has(k)) {
      const swaps = swapped ? SWAPS : [];
      const run = runCommissary({ strategy: "ahead", swaps, sigma: 0.1 });
      traffic.set(k, trafficOf(run.plan.days, menuFor(swaps)));
    }
    return traffic.get(k)!;
  };

  const W = manifest.room_shell.width_m;
  const D = manifest.room_shell.depth_m;
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x12100c);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.55;
  const centre = new THREE.Vector3(W / 2, 0, -D * 0.4);
  const lights = studioLights(scene, centre, 11.5, phone ? 1024 : 2048);
  lights.key.position.set(W * 0.75, 15, 4);
  lights.key.intensity = 1.9;
  lights.hemi.intensity = 0.6;

  const fl = tilePlane(W, D);
  fl.rotation.x = -Math.PI / 2;
  fl.position.set(W / 2, 0, -D / 2);
  scene.add(fl);
  // The walls cut down to a kerb, so the room reads as a plan from above.
  const kerbMat = new THREE.MeshStandardMaterial({ color: 0x3a342c, roughness: 0.8 });
  const kerb = (x0: number, x1: number, z0: number, z1: number) => {
    const m = new THREE.Mesh(new THREE.BoxGeometry(x1 - x0, 0.16, z1 - z0), kerbMat);
    m.position.set((x0 + x1) / 2, 0.08, (z0 + z1) / 2);
    m.receiveShadow = true;
    scene.add(m);
  };
  kerb(-0.15, W + 0.15, 0, 0.15);
  kerb(-0.15, W + 0.15, -D - 0.15, -D);
  kerb(-0.15, 0, -D, 0);
  kerb(W, W + 0.15, -D, 0);

  /* The fixtures, from the manifest, merged per material (the chiller is built on its own: it moves). */
  const statics: THREE.Object3D[] = [];
  for (const f of manifest.fixtures) {
    if (SMALL.test(f.archetype_id) || f.archetype_id === CHILLER_ID) continue;
    const [w, d, h] = f.local_extents_m;
    const m = buildFixture(f.archetype_id, { w, d, h }, (phone ? 0 : HOT.test(f.archetype_id) ? 1 : 0) as 0 | 1, true);
    placeFixture(m.group, f.pose);
    m.hooks.set?.("lights", 1);
    if (/WALKIN-COOLER/.test(f.archetype_id)) m.hooks.set?.("roof", 1);
    if (/WIRE-SHELVING/.test(f.archetype_id)) m.hooks.set?.("stock", 0.34);
    statics.push(m.group);
  }
  const room = new THREE.Group();
  scene.add(room);
  mergeStatic(statics, room);

  const chillerFx = manifest.fixtures.find((f) => f.archetype_id === CHILLER_ID)!;
  const [cw, cd, ch] = chillerFx.local_extents_m;
  const chiller = buildFixture(CHILLER_ID, { w: cw, d: cd, h: ch }, phone ? 1 : 2, true);
  chiller.hooks.set?.("lights", 1);
  scene.add(chiller.group);
  const ghost = buildFixture(CHILLER_ID, { w: cw, d: cd, h: ch }, 0, false);
  const ghostMat = new THREE.MeshBasicMaterial({ color: 0xa7f3c0, transparent: true, opacity: 0.22, depthWrite: false });
  ghost.group.traverse((o) => {
    const m = o as THREE.Mesh;
    if (m.isMesh) m.material = ghostMat;
  });
  scene.add(ghost.group);

  /** Stand the chiller (or its ghost) at a spot: its front on the spot's front edge, facing the room. */
  const standAt = (g: THREE.Object3D, s: { stance: readonly [number, number]; yaw: number }) => {
    const [fx, fz] = fwd(s.yaw);
    const frontX = s.stance[0] - 0.5 * fx;
    const frontZ = s.stance[1] - 0.5 * fz;
    g.rotation.set(0, s.yaw, 0);
    g.position.set(frontX - (cw / 2) * Math.cos(s.yaw), 0, frontZ + (cw / 2) * Math.sin(s.yaw));
  };
  const drawn = { rect: floor.chiller.rect, stance: floor.chiller.stance, yaw: floor.chiller.yaw };
  standAt(ghost.group, drawn);

  /* The strip along the walls: a week's walking for a chiller at each spot. */
  const strip = new THREE.Mesh(new THREE.BufferGeometry(), new THREE.MeshBasicMaterial({ vertexColors: true, transparent: true, opacity: 1, depthWrite: false, side: THREE.DoubleSide, toneMapped: false }));
  strip.renderOrder = 2;
  scene.add(strip);
  const good = new THREE.Color(0x5fe08a);
  const mid = new THREE.Color(0xffc23d);
  const bad = new THREE.Color(0xff4d6a);
  const paintStrip = (cost: number[]) => {
    const sorted = [...cost].sort((a, b) => a - b);
    const lo = sorted[0]!;
    const hi = sorted[Math.floor(sorted.length * 0.9)]!;
    const pos: number[] = [];
    const col: number[] = [];
    const c = new THREE.Color();
    spots.forEach((s, i) => {
      const t = Math.min(1, Math.max(0, (cost[i]! - lo) / Math.max(1e-6, hi - lo)));
      if (t < 0.5) c.copy(good).lerp(mid, t * 2);
      else c.copy(mid).lerp(bad, (t - 0.5) * 2);
      // A 0.1 m band at the spot's centre, across the chiller's depth, on the floor.
      const [fx, fz] = fwd(s.yaw);
      const ax = Math.cos(s.yaw);
      const az = -Math.sin(s.yaw);
      const cx = s.stance[0] - fx * (0.5 + cd / 2);
      const cz = s.stance[1] - fz * (0.5 + cd / 2);
      const hw = 0.055;
      const hd = cd / 2;
      const q = [
        [cx - ax * hw - fx * hd, cz - az * hw - fz * hd],
        [cx + ax * hw - fx * hd, cz + az * hw - fz * hd],
        [cx + ax * hw + fx * hd, cz + az * hw + fz * hd],
        [cx - ax * hw + fx * hd, cz - az * hw + fz * hd],
      ];
      for (const k of [0, 1, 2, 0, 2, 3]) {
        pos.push(q[k]![0]!, 0.014, q[k]![1]!);
        col.push(c.r, c.g, c.b);
      }
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute("color", new THREE.Float32BufferAttribute(col, 3));
    strip.geometry.dispose();
    strip.geometry = g;
  };

  /* The routes. */
  const dash = dashTexture();
  // Drawn over the brick as they are (not added to it), so the ice blue stays blue.
  const hotMat = new THREE.MeshBasicMaterial({ color: 0xffb35c, map: dash, transparent: true, opacity: 1, depthWrite: false, side: THREE.DoubleSide, toneMapped: false });
  const coldMat = new THREE.MeshBasicMaterial({ color: 0x6fd3ff, map: dash, transparent: true, opacity: 1, depthWrite: false, side: THREE.DoubleSide, toneMapped: false });
  const routeGroup = new THREE.Group();
  routeGroup.renderOrder = 3;
  scene.add(routeGroup);
  const drawRoutes = (routes: Route[]) => {
    for (const c of [...routeGroup.children]) {
      (c as THREE.Mesh).geometry.dispose();
      routeGroup.remove(c);
    }
    const maxTrips = Math.max(...routes.map((r) => r.trips), 1);
    for (const r of routes) {
      if (r.trips <= 0) continue;
      const cold = r.to === WALKIN;
      const w = cold ? 0.16 : 0.07 + 0.2 * Math.sqrt(r.trips / maxTrips);
      // Hot pans flow towards the chiller, loads from it to the walk-in: the dash runs from the path's start.
      const m = new THREE.Mesh(ribbon(r.path, w, cold ? 0.022 : 0.02, 0.9), cold ? coldMat : hotMat);
      m.renderOrder = 3;
      routeGroup.add(m);
    }
  };

  /* The cook and the pan. */
  const st: Staff = makeStaff("p50", phone ? 1 : 2, { seed: 11, apron: "navy", sleeves: "rolled" });
  scene.add(st.group);
  const pan = buildFixture("NSF-HOTEL-PAN-FULL-4IN", PAN_EXT, 1);
  pan.hooks.set?.("food", 2);
  pan.group.position.set(-PAN_EXT.w / 2, -PAN_EXT.h, PAN_EXT.d / 2);
  st.item.add(pan.group);
  let routine: Routine | null = null;
  /** Time windows in the loop when the pan is in the cook's hands, and when the chiller's door is open. */
  let carrying: [number, number][] = [];
  let doorOpen: [number, number][] = [];

  const kindOf = (id: string): Exclude<keyof typeof RAIL, "chiller"> => (/combi/.test(id) ? "combi" : /kettle/.test(id) ? "kettle" : /skillet/.test(id) ? "skillet" : "prep");
  /** The way a cook faces to work a fixture: towards its front. */
  const sourceYaw = (id: string) => floor.stations.get(id)!.yaw + Math.PI;

  /** The cook's loop: a pan from each station in turn, as often as the week sends one, to the chiller and back. */
  const buildRoutine = (routes: Route[], spot: { stance: readonly [number, number]; yaw: number }) => {
    // Six trips a loop, shared between the stations by their pans (largest remainder), each station's fixtures taken in
    // turn, busiest first, and the stations interleaved.
    const n = 6;
    const byKind = new Map<string, Route[]>();
    for (const r of routes) {
      if (r.to !== "chiller" || r.trips <= 0) continue;
      const k = kindOf(r.from);
      if (!byKind.has(k)) byKind.set(k, []);
      byKind.get(k)!.push(r);
    }
    const groups = [...byKind.values()].map((rs) => rs.sort((a, b) => b.trips - a.trips));
    const weight = groups.map((rs) => rs.reduce((s, r) => s + r.trips, 0));
    const total = weight.reduce((s, w) => s + w, 0);
    const raw = weight.map((w) => (w / Math.max(1e-9, total)) * n);
    const counts = raw.map(Math.floor);
    let left = total > 0 ? n - counts.reduce((s, c) => s + c, 0) : 0;
    raw.map((v, i) => [v - Math.floor(v), i] as const).sort((a, b) => b[0] - a[0]).forEach(([, i]) => {
      if (left-- > 0) counts[i]!++;
    });
    const order: Route[] = [];
    for (let k = 0; k < Math.max(0, ...counts); k++) groups.forEach((rs, i) => k < counts[i]! && order.push(rs[k % rs.length]!));
    if (!order.length) {
      routine = null;
      return;
    }
    const steps: RoutineStep[] = [];
    // The chiller's runners: in from its front face, which is half a metre past the stance (wallSpots, and the plan's own
    // stance socket).
    const [sfx, sfz] = fwd(spot.yaw);
    const reach = 0.5 + RAIL.chiller.inset;
    const chillAt = new THREE.Vector3(spot.stance[0] - sfx * reach, RAIL.chiller.y, spot.stance[1] - sfz * reach);
    const first = order[0]!;
    const s0 = floor.stations.get(first.from)!.stance;
    const startYaw = sourceYaw(first.from);
    carrying = [];
    doorOpen = [];
    const push = (s: RoutineStep) => steps.push(s);
    order.forEach((r, i) => {
      const station = floor.stations.get(r.from)!;
      const src = station.stance;
      const yaw = sourceYaw(r.from);
      const [fx, fz] = fwd(yaw);
      const rail = RAIL[kindOf(r.from)];
      // The stance is `out` in front of the fixture's front, and the pan `inset` in from it.
      const at = new THREE.Vector3(src[0] + fx * (station.out + rail.inset), rail.y, src[1] + fz * (station.out + rail.inset));
      if (i > 0) {
        // Back from the chiller to this source, empty-handed: the reverse of its own route.
        const back = [...r.path].reverse();
        push({ kind: "walk", path: back.slice(1), speed: CARRY_MPS, upper: "free" });
        push({ kind: "turn", yaw, upper: "free" });
      }
      const pickS = panSeconds("pick", {});
      const t0 = steps.length;
      push({ kind: "work", work: "pick", duration: pickS, opts: { at, surfaceY: rail.y } });
      push({ kind: "walk", path: r.path.slice(1), speed: CARRY_MPS, upper: "carry" });
      push({ kind: "turn", yaw: spot.yaw + Math.PI, upper: "carry" });
      const placeS = panSeconds("place", {});
      push({ kind: "work", work: "place", duration: placeS, opts: { at: chillAt, surfaceY: RAIL.chiller.y } });
      carrying.push([t0, steps.length - 1]);
    });
    // And back to the first source to begin again.
    push({ kind: "walk", path: [...first.path].reverse().slice(1), speed: CARRY_MPS, upper: "free" });
    push({ kind: "turn", yaw: startYaw, upper: "free" });
    routine = new Routine(st, { x: s0[0], z: s0[1], yaw: startYaw }, steps, true);
    const T = routine.stepTimes;
    // Step indices to times: the pan is held from its pick's start to its place's end; the door is open round a place.
    carrying = carrying.map(([a, b]) => [T[a]!, T[b + 1]!]);
    doorOpen = carrying.map(([, b]) => [b - panSeconds("place", {}) - 1.0, b + 0.6]);
  };

  /* Camera: from in front of the room, high, over the cold store, the hot bank and the chiller's walls. */
  const camera = new THREE.PerspectiveCamera(30, 4 / 3, 0.1, 120);
  const TARGET = new THREE.Vector3(W / 2 + 0.6, 0, -7.2);
  const VIEW = new THREE.Vector3(0.1, 1.05, 1).normalize();
  // The cold store, the hot bank, and the chiller's walls between them; the dish room at the back may crop.
  const keep: THREE.Vector3[] = [];
  for (const x of [0.6, W]) for (const z of [0.2, -13.6]) for (const y of [0, 1.2]) keep.push(new THREE.Vector3(x, y, z));
  const gate = shaderGate(renderer, camera);
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 4, maxDistance: 60, maxPolarAngle: Math.PI * 0.45 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const ndc = new THREE.Vector3();
  const frame = () => {
    if (userMoved) return;
    camera.aspect = Math.max(0.2, stageEl.clientWidth / Math.max(1, stageEl.clientHeight));
    camera.updateProjectionMatrix();
    // A narrower stage wraps its legend onto two lines: aim nearer the front so the room sits above it.
    TARGET.z = camera.aspect < 1.5 ? -5.9 : -7.2;
    let dist = 8;
    for (; dist < 80; dist += 0.25) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      const bottom = camera.aspect < 1.5 ? -0.74 : -0.9;
      if (keep.every((p) => (ndc.copy(p).project(camera), Math.abs(ndc.x) < 0.97 && ndc.y < 0.95 && ndc.y > bottom))) break;
    }
    controls.target.copy(TARGET);
    controls.update();
  };
  bindCanvasResize(stageEl, renderer, camera, gate);
  frame();
  new ResizeObserver(() => frame()).observe(stageEl);

  /* The page. */
  const menuBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-fl-menu]"));
  const spotBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-fl-spot]"));
  const statusEl = root.querySelector<HTMLElement>("[data-fl-status]");
  const hudSpot = root.querySelector<HTMLElement>("[data-fl-hud-spot]");
  const hudPan = root.querySelector<HTMLElement>("[data-fl-hud-pan]");
  const out = (k: string) => root.querySelector<HTMLElement>(`[data-fl-${k}]`);
  let swapped = false;
  let choice: "drawn" | "best" | number = "drawn";
  const answers = new Map<boolean, FloorAnswer>();
  const answerFor = (sw: boolean): FloorAnswer => {
    if (!answers.has(sw)) answers.set(sw, placeChiller(floor, fields, trafficFor(sw), spots));
    return answers.get(sw)!;
  };
  let answer = answerFor(false);
  let current = answer.drawn;

  const km = (m: number) => `${(m / 1000).toFixed(m < 10000 ? 2 : 1)} km`;
  const mins = (m: number) => {
    const s = m / CARRY_MPS;
    const h = Math.floor(s / 3600);
    const mm = Math.round((s - h * 3600) / 60);
    return h ? `${h} h ${String(mm).padStart(2, "0")} min` : `${mm} min`;
  };
  /** The walk to the chiller, metres per pan, averaged over the week's pans. */
  const perPan = (routes: Route[]) => {
    const hot = routes.filter((r) => r.to === "chiller");
    const n = hot.reduce((s, r) => s + r.trips, 0);
    return n ? hot.reduce((s, r) => s + r.trips * r.metres, 0) / n : 0;
  };
  const where = (s: { stance: readonly [number, number] }) => {
    // Name a spot by what it stands nearest, among the fixtures people know.
    const near: [string, string][] = [
      ["combi-rollin-2", "the combis"],
      ["kettle-40gal-1", "the kettle"],
      ["tilt-skillet-40gal-1", "the tilting skillet"],
      [WALKIN, "the walk-in"],
      ["walkin-freezer-10x14", "the freezer"],
      ["preptable-prep-1", "the prep tables"],
      ["dishwasher-conveyor-1", "the dish room"],
    ];
    let best = near[0]![1];
    let bd = Infinity;
    for (const [id, name] of near) {
      const p = floor.stations.get(id)?.stance;
      if (!p) continue;
      const d = Math.hypot(p[0] - s.stance[0], p[1] - s.stance[1]);
      if (d < bd) {
        bd = d;
        best = name;
      }
    }
    return best;
  };

  const update = () => {
    const t = trafficFor(swapped);
    answer = answerFor(swapped);
    paintStrip(answer.cost);
    const spot = choice === "drawn" ? drawn : choice === "best" ? answer.spots[answer.best]! : answer.spots[choice]!;
    current = choice === "drawn" ? answer.drawn : choice === "best" ? answer.atBest : routesFor(floor, t, spot);
    standAt(chiller.group, spot);
    // The ghost stands where the chiller is not: at the shortest walk while it stands as drawn, and as drawn otherwise.
    standAt(ghost.group, choice === "drawn" ? answer.spots[answer.best]! : drawn);
    ghostMat.color.set(choice === "drawn" ? 0x7cf0a4 : 0xd9d4cc);
    drawRoutes(current.routes);
    buildRoutine(current.routes, spot);
    menuBtns.forEach((b) => b.setAttribute("aria-pressed", String((b.dataset.flMenu === "swapped") === swapped)));
    spotBtns.forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.flSpot === choice)));
    const m = current.metres;
    const d0 = answer.drawn.metres;
    const bestSpot = answer.spots[answer.best]!;
    const set = (k: string, v: string) => {
      const el = out(k);
      if (el) el.textContent = v;
    };
    set("km", km(m));
    set("pan", `${perPan(current.routes).toFixed(1)} m`);
    set("time", mins(m));
    set("pans", `${Math.round(t.totalPans)}`);
    const dEl = out("km-d");
    if (dEl) {
      const dm = m - d0;
      dEl.hidden = choice === "drawn" || Math.abs(dm) < 5;
      dEl.textContent = `${dm < 0 ? "−" : "+"}${Math.round(Math.abs(dm / d0) * 100)}%`;
      dEl.dataset.good = String(dm < 0);
    }
    if (hudSpot) hudSpot.textContent = choice === "drawn" ? "Chiller as drawn" : choice === "best" ? "Chiller at the best spot" : "Chiller where you put it";
    if (hudPan) hudPan.textContent = `${perPan(current.routes).toFixed(1)} m per pan`;
    if (statusEl) {
      const p0 = perPan(answer.drawn.routes);
      const pb = perPan(answer.atBest.routes);
      statusEl.textContent =
        `As drawn, the chiller stands by ${where(drawn)} and a pan is carried ${p0.toFixed(1)} m to it. ` +
        `With this menu the shortest walk is beside ${where(bestSpot)}: ${pb.toFixed(1)} m a pan, ${km(d0 - answer.atBest.metres)} a week less.`;
      statusEl.dataset.state = "ok";
    }
  };
  menuBtns.forEach((b) =>
    b.addEventListener("click", () => {
      swapped = b.dataset.flMenu === "swapped";
      update();
    }),
  );
  spotBtns.forEach((b) =>
    b.addEventListener("click", () => {
      choice = b.dataset.flSpot === "best" ? "best" : "drawn";
      update();
    }),
  );
  // A click (not a drag) on the strip stands the chiller at the spot under it.
  const ray = new THREE.Raycaster();
  const ground = new THREE.Plane(new THREE.Vector3(0, 1, 0), 0);
  let down: [number, number] | null = null;
  canvas.addEventListener("pointerdown", (e) => {
    down = [e.clientX, e.clientY];
  });
  canvas.addEventListener("pointerup", (e) => {
    if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 6) return;
    const r = canvas.getBoundingClientRect();
    ray.setFromCamera(new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1), camera);
    const hit = new THREE.Vector3();
    if (!ray.ray.intersectPlane(ground, hit)) return;
    let bi = -1;
    let bd = 0.9;
    spots.forEach((s, i) => {
      const [fx, fz] = fwd(s.yaw);
      const d = Math.hypot(s.stance[0] - fx * (0.5 + cd / 2) - hit.x, s.stance[1] - fz * (0.5 + cd / 2) - hit.z);
      if (d < bd) {
        bd = d;
        bi = i;
      }
    });
    if (bi < 0) return;
    choice = bi === answer.best ? "best" : bi;
    update();
  });

  update();
  const buildMs = Math.round(performance.now() - mountT0);
  void st.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let clock = 0;
  let frozen: number | null = null;
  let last = performance.now();
  const inWin = (w: [number, number][], t: number) => w.some(([a, b]) => t >= a && t < b);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    clock = frozen ?? clock + dt;
    dash.offset.x = -clock * (CARRY_MPS / 0.9);
    if (routine) {
      const t = clock % routine.duration;
      routine.apply(t);
      pan.group.visible = inWin(carrying, t);
      const open = doorOpen.find(([a, b]) => t >= a && t < b);
      chiller.hooks.set?.("door", open ? Math.min(1, (t - open[0]) / 0.6, (open[1] - t) / 0.6) : 0);
    }
    controls.update();
    fitClip();
    renderer.render(scene, camera);
  }
  // Depth precision for the high, distant view: the walk-ins' roof panels lie under a millimetre off their frames.
  const fitClip = clipFitter(camera, scene, 0.1);
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  Object.assign(stageEl, {
    floor: {
      info: () => ({
        buildMs,
        tris: renderer.info.render.triangles,
        calls: renderer.info.render.calls,
        spots: spots.length,
        best: answer.spots[answer.best]!.stance,
        drawnM: Math.round(answer.drawn.metres),
        bestM: Math.round(answer.atBest.metres),
        loopS: routine?.duration ?? 0,
      }),
      hold: (t: number | null) => {
        frozen = t;
      },
      view: (eye: [number, number, number], at: [number, number, number]) => {
        userMoved = true;
        camera.position.set(...eye);
        controls.target.set(...at);
        controls.update();
      },
    },
  });
}
