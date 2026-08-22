// Kitty Hawk scene dressing, three.js layer (bead guzez.13). THIN
// consumer of the tested math in dressing.ts: this file only builds
// meshes and applies poses. Procedural canvas textures keep the app
// self-contained (no asset fetches; the CSP and the size budget both
// stay honest).

import * as THREE from "three";
import {
  campLayout,
  gullFleet,
  gullPose,
  lcg,
  orvillePose,
  railTies,
  type GullPath,
} from "./dressing.ts";

/* ---------- procedural textures (canvas, deterministic) ---------- */

function canvasTexture(
  size: number,
  draw: (ctx: CanvasRenderingContext2D, size: number) => void,
): THREE.CanvasTexture {
  const canvas = document.createElement("canvas");
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext("2d");
  if (ctx === null) {
    // Headless/degraded: a 1x1 texture keeps materials valid.
    canvas.width = 1;
    canvas.height = 1;
  } else {
    draw(ctx, size);
  }
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/** Speckled dune sand with faint wind-ripple banding. */
function sandTexture(): THREE.CanvasTexture {
  const rand = lcg(1903);
  return canvasTexture(512, (ctx, s) => {
    ctx.fillStyle = "#c2b088";
    ctx.fillRect(0, 0, s, s);
    // Wind ripples: soft diagonal bands.
    for (let i = 0; i < 42; i += 1) {
      const y = rand() * s;
      ctx.strokeStyle = rand() < 0.5 ? "rgba(160,140,100,0.16)" : "rgba(230,215,180,0.14)";
      ctx.lineWidth = 2 + rand() * 5;
      ctx.beginPath();
      ctx.moveTo(0, y);
      ctx.bezierCurveTo(s * 0.3, y + 14 * (rand() - 0.5), s * 0.7, y + 14 * (rand() - 0.5), s, y);
      ctx.stroke();
    }
    // Grain speckle.
    for (let i = 0; i < 9000; i += 1) {
      const g = 150 + Math.floor(rand() * 90);
      ctx.fillStyle = `rgba(${g},${g - 16},${g - 46},0.28)`;
      ctx.fillRect(rand() * s, rand() * s, 1.4, 1.4);
    }
  });
}

/** Weathered plank siding for the camp buildings. */
function plankTexture(): THREE.CanvasTexture {
  const rand = lcg(1911);
  return canvasTexture(256, (ctx, s) => {
    ctx.fillStyle = "#8a7354";
    ctx.fillRect(0, 0, s, s);
    const plank = s / 8;
    for (let i = 0; i < 8; i += 1) {
      const shade = 108 + Math.floor(rand() * 40);
      ctx.fillStyle = `rgb(${shade},${shade - 22},${shade - 46})`;
      ctx.fillRect(0, i * plank + 1, s, plank - 2);
      // Wood grain streaks.
      for (let k = 0; k < 6; k += 1) {
        ctx.strokeStyle = `rgba(60,44,26,${0.12 + rand() * 0.15})`;
        ctx.lineWidth = 1;
        const y = i * plank + rand() * plank;
        ctx.beginPath();
        ctx.moveTo(0, y);
        ctx.lineTo(s, y + (rand() - 0.5) * 3);
        ctx.stroke();
      }
    }
  });
}

/** Sky: horizon-to-zenith gradient with a low December sun glow. */
function skyTexture(): THREE.CanvasTexture {
  return canvasTexture(1024, (ctx, s) => {
    // Sphere v = 1 - y/s and the HORIZON is the equator at v = 0.5,
    // i.e. canvas y = 0.5s — the pale sea haze belongs THERE, not at
    // the bottom pole (which points underground).
    const grad = ctx.createLinearGradient(0, s, 0, 0);
    grad.addColorStop(0, "#cfd8de"); // below the horizon (rarely seen)
    grad.addColorStop(0.5, "#d8dfe6"); // sea haze AT the horizon
    grad.addColorStop(0.68, "#b6c9dd");
    grad.addColorStop(0.85, "#7fa3cd");
    grad.addColorStop(1, "#4f7cb4"); // zenith
    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, s, s);
    // Low winter sun (south-east), just above the dome's equator —
    // canvas y maps to sphere v = 1 - y/s, and the horizon sits at
    // v = 0.5, so y must be < 0.5s or the sun paints underground.
    const sun = ctx.createRadialGradient(s * 0.68, s * 0.43, 8, s * 0.68, s * 0.43, s * 0.28);
    sun.addColorStop(0, "rgba(255,244,214,0.95)");
    sun.addColorStop(0.12, "rgba(255,238,200,0.5)");
    sun.addColorStop(1, "rgba(255,238,200,0)");
    ctx.fillStyle = sun;
    ctx.fillRect(0, 0, s, s);
  });
}

/** One soft cloud billboard. */
function cloudTexture(seed: number): THREE.CanvasTexture {
  const rand = lcg(seed);
  return canvasTexture(256, (ctx, s) => {
    ctx.clearRect(0, 0, s, s);
    for (let i = 0; i < 14; i += 1) {
      const x = s * (0.2 + rand() * 0.6);
      const y = s * (0.35 + rand() * 0.3);
      const r = s * (0.08 + rand() * 0.16);
      const puff = ctx.createRadialGradient(x, y, r * 0.1, x, y, r);
      puff.addColorStop(0, "rgba(255,255,255,0.55)");
      puff.addColorStop(1, "rgba(255,255,255,0)");
      ctx.fillStyle = puff;
      ctx.fillRect(0, 0, s, s);
    }
  });
}

/* ---------------------- shared materials ------------------------- */

const WOOD = new THREE.MeshStandardMaterial({ color: 0x6f5231, roughness: 0.9 });
const WOOD_DARK = new THREE.MeshStandardMaterial({ color: 0x4c3a22, roughness: 0.95 });
const IRON = new THREE.MeshStandardMaterial({ color: 0x3a3f45, roughness: 0.6, metalness: 0.55 });
const SUIT = new THREE.MeshStandardMaterial({ color: 0x2d2c33, roughness: 0.9 });
const SKIN = new THREE.MeshStandardMaterial({ color: 0xc99b78, roughness: 0.8 });
const GULL_BODY = new THREE.MeshStandardMaterial({ color: 0xe8e8e4, roughness: 0.85 });
const GULL_WING = new THREE.MeshStandardMaterial({
  color: 0xd6d6d0,
  roughness: 0.85,
  side: THREE.DoubleSide,
});

function box(
  w: number,
  h: number,
  d: number,
  mat: THREE.Material,
  x = 0,
  y = 0,
  z = 0,
): THREE.Mesh {
  const m = new THREE.Mesh(new THREE.BoxGeometry(w, h, d), mat);
  m.position.set(x, y, z);
  return m;
}

/* ------------------------- environment --------------------------- */

/** Sky dome (BackSide sphere) — replaces the flat clear color. */
export function buildSky(): THREE.Mesh {
  const geo = new THREE.SphereGeometry(2600, 32, 18);
  const mat = new THREE.MeshBasicMaterial({
    map: skyTexture(),
    side: THREE.BackSide,
    fog: false,
    depthWrite: false,
  });
  const sky = new THREE.Mesh(geo, mat);
  sky.renderOrder = -10;
  return sky;
}

/** Drifting cloud billboards (updated by animateDressing). */
export function buildClouds(): THREE.Group {
  const group = new THREE.Group();
  const rand = lcg(1908);
  for (let i = 0; i < 9; i += 1) {
    const mat = new THREE.MeshBasicMaterial({
      map: cloudTexture(200 + i),
      transparent: true,
      depthWrite: false,
      opacity: 0.85,
      fog: false,
    });
    const w = 260 + rand() * 320;
    const cloud = new THREE.Mesh(new THREE.PlaneGeometry(w, w * 0.42), mat);
    mat.side = THREE.DoubleSide;
    cloud.rotation.x = -Math.PI / 2 + 0.35; // near-horizontal, tipped at the camera
    cloud.position.set(0, 190 + rand() * 160, -900 + rand() * 1800);
    cloud.userData["baseX"] = -900 + rand() * 1800;
    cloud.userData["driftMps"] = 1.5 + rand() * 2.2;
    group.add(cloud);
  }
  return group;
}

/** Sand skirt far beyond the surveyed tile so the horizon is never
 * void, plus (Kill Devil Hills only) the Atlantic to the EAST —
 * Huffman Prairie is landlocked Ohio pasture and gets NO ocean. */
export function buildOuterGround(tileExtentM: number, withOcean: boolean): THREE.Group {
  const group = new THREE.Group();
  const sand = sandTexture();
  sand.wrapS = THREE.RepeatWrapping;
  sand.wrapT = THREE.RepeatWrapping;
  sand.repeat.set(160, 160);
  const skirt = new THREE.Mesh(
    new THREE.CircleGeometry(2400, 48),
    new THREE.MeshStandardMaterial({ map: sand, color: 0xcabb95, roughness: 1 }),
  );
  skirt.rotation.x = -Math.PI / 2;
  skirt.position.y = -0.35; // tucked under the surveyed tile
  group.add(skirt);
  if (!withOcean) {
    return group;
  }
  const ocean = new THREE.Mesh(
    new THREE.PlaneGeometry(2600, 5200),
    new THREE.MeshStandardMaterial({
      color: 0x2e5468,
      roughness: 0.35,
      metalness: 0.1,
      transparent: true,
      opacity: 0.96,
    }),
  );
  ocean.rotation.x = -Math.PI / 2;
  ocean.position.set(tileExtentM / 2 + 1300, -0.15, 0);
  group.add(ocean);
  return group;
}

/** The surveyed tile gets the same sand texture via a second, textured
 * material — vertex colors still tint water/sand/dune classes. */
export function sandTileMaterial(): THREE.MeshStandardMaterial {
  const sand = sandTexture();
  sand.wrapS = THREE.RepeatWrapping;
  sand.wrapT = THREE.RepeatWrapping;
  sand.repeat.set(130, 130);
  return new THREE.MeshStandardMaterial({ map: sand, vertexColors: true, roughness: 1 });
}

/* --------------------------- the rail ---------------------------- */

/** The 60 ft monorail: 2x4 on edge, half-buried ties, and the small
 * starting trestle. Positions are launch-relative; caller places. */
export function buildRail(railLengthM: number): THREE.Group {
  const group = new THREE.Group();
  const rail = box(railLengthM, 0.09, 0.04, WOOD_DARK, railLengthM / 2, 0.16, 0);
  group.add(rail);
  const cap = box(railLengthM, 0.015, 0.09, IRON, railLengthM / 2, 0.21, 0);
  group.add(cap);
  for (const x of railTies(railLengthM)) {
    group.add(box(0.35, 0.06, 1.15, WOOD, x, 0.05, 0));
  }
  // Starting trestle: two A-frames + a cross bench at the tail.
  for (const dz of [-0.45, 0.45]) {
    group.add(box(0.07, 0.5, 0.07, WOOD, -0.6, 0.25, dz));
  }
  group.add(box(0.09, 0.07, 1.15, WOOD, -0.6, 0.52, 0));
  return group;
}

/* --------------------------- the camp ---------------------------- */

function buildBuilding(widthM: number, depthM: number, heightM: number): THREE.Group {
  const group = new THREE.Group();
  const planks = plankTexture();
  const wall = new THREE.MeshStandardMaterial({ map: planks, roughness: 0.95 });
  const body = box(widthM, heightM, depthM, wall, 0, heightM / 2, 0);
  group.add(body);
  // Gable roof: two pitched slabs.
  const half = widthM / 2;
  const rise = heightM * 0.45;
  const slabLen = Math.hypot(half, rise) + 0.15;
  for (const side of [-1, 1]) {
    const slab = box(slabLen, 0.06, depthM + 0.4, WOOD_DARK);
    slab.position.set((side * half) / 2, heightM + rise / 2, 0);
    slab.rotation.z = -side * Math.atan2(rise, half);
    group.add(slab);
  }
  // Dark doorway on the +x gable end.
  const door = box(0.04, heightM * 0.62, 1.15, WOOD_DARK, widthM / 2 + 0.01, heightM * 0.31, 0);
  door.material = new THREE.MeshStandardMaterial({ color: 0x17120b, roughness: 1 });
  group.add(door);
  return group;
}

export interface Campfire {
  group: THREE.Group;
  light: THREE.PointLight;
  flame: THREE.Mesh;
}

function buildCampfire(): Campfire {
  const group = new THREE.Group();
  const rand = lcg(1917);
  // Stone ring + charred logs.
  for (let i = 0; i < 9; i += 1) {
    const a = (i / 9) * Math.PI * 2;
    const stone = new THREE.Mesh(
      new THREE.DodecahedronGeometry(0.14 + rand() * 0.07),
      new THREE.MeshStandardMaterial({ color: 0x8f8a80, roughness: 1 }),
    );
    stone.position.set(0.62 * Math.cos(a), 0.08, 0.62 * Math.sin(a));
    group.add(stone);
  }
  for (let i = 0; i < 4; i += 1) {
    const log = new THREE.Mesh(new THREE.CylinderGeometry(0.055, 0.07, 0.8, 6), WOOD_DARK);
    log.rotation.set(Math.PI / 2.3, (i / 4) * Math.PI, 0);
    log.position.set(0, 0.12, 0);
    group.add(log);
  }
  const flame = new THREE.Mesh(
    new THREE.ConeGeometry(0.16, 0.5, 7),
    new THREE.MeshBasicMaterial({ color: 0xff9a33, transparent: true, opacity: 0.85 }),
  );
  flame.position.y = 0.38;
  group.add(flame);
  const light = new THREE.PointLight(0xff9440, 14, 22, 2);
  light.position.set(0, 0.9, 0);
  group.add(light);
  return { group, light, flame };
}

function buildChair(): THREE.Group {
  const g = new THREE.Group();
  g.add(box(0.42, 0.04, 0.42, WOOD, 0, 0.45, 0));
  g.add(box(0.42, 0.5, 0.04, WOOD, 0, 0.73, -0.19));
  for (const [dx, dz] of [
    [-0.18, -0.18],
    [0.18, -0.18],
    [-0.18, 0.18],
    [0.18, 0.18],
  ] as const) {
    g.add(box(0.04, 0.45, 0.04, WOOD, dx, 0.225, dz));
  }
  return g;
}

function buildBarrel(): THREE.Mesh {
  return new THREE.Mesh(new THREE.CylinderGeometry(0.3, 0.26, 0.85, 12), WOOD);
}

function buildWorkbench(): THREE.Group {
  const g = new THREE.Group();
  g.add(box(1.8, 0.07, 0.7, WOOD, 0, 0.85, 0));
  for (const [dx, dz] of [
    [-0.8, -0.28],
    [0.8, -0.28],
    [-0.8, 0.28],
    [0.8, 0.28],
  ] as const) {
    g.add(box(0.07, 0.85, 0.07, WOOD_DARK, dx, 0.42, dz));
  }
  // Tools on top: a saw blade, a hammer, a can.
  g.add(box(0.5, 0.015, 0.1, IRON, -0.4, 0.9, 0.05));
  g.add(box(0.06, 0.06, 0.28, WOOD_DARK, 0.25, 0.92, -0.1));
  const can = new THREE.Mesh(new THREE.CylinderGeometry(0.07, 0.07, 0.16, 10), IRON);
  can.position.set(0.6, 0.97, 0.15);
  g.add(can);
  return g;
}

function buildToolchest(): THREE.Group {
  const g = new THREE.Group();
  g.add(box(0.9, 0.4, 0.5, WOOD_DARK, 0, 0.2, 0));
  g.add(box(0.94, 0.05, 0.54, WOOD, 0, 0.44, 0));
  return g;
}

/* -------------------------- the people --------------------------- */

export interface Figure {
  group: THREE.Group;
  leftLeg: THREE.Mesh;
  rightLeg: THREE.Mesh;
  leftArm: THREE.Group;
  rightArm: THREE.Group;
  glasses: THREE.Mesh;
}

/** A standing 1900s figure (dark suit, flat cap) with poseable limbs
 * — Orville on the ground. ~1.78 m tall. */
export function buildFigure(): Figure {
  const group = new THREE.Group();
  const torso = box(0.34, 0.62, 0.22, SUIT, 0, 1.18, 0);
  group.add(torso);
  const head = new THREE.Mesh(new THREE.SphereGeometry(0.115, 12, 10), SKIN);
  head.position.set(0, 1.66, 0);
  group.add(head);
  const cap = new THREE.Mesh(new THREE.CylinderGeometry(0.13, 0.135, 0.055, 12), SUIT);
  cap.position.set(0, 1.755, 0);
  group.add(cap);
  const brim = box(0.12, 0.015, 0.1, SUIT, 0, 1.73, 0.13);
  group.add(brim);
  const mkLeg = (side: number): THREE.Mesh => {
    const leg = box(0.11, 0.85, 0.13, SUIT, side * 0.09, 0.445, 0);
    leg.geometry.translate(0, -0.36, 0);
    leg.position.y = 0.85;
    return leg;
  };
  const leftLeg = mkLeg(-1);
  const rightLeg = mkLeg(1);
  group.add(leftLeg, rightLeg);
  const mkArm = (side: number): THREE.Group => {
    const arm = new THREE.Group();
    arm.position.set(side * 0.22, 1.44, 0);
    const upper = box(0.09, 0.52, 0.1, SUIT, 0, -0.24, 0);
    arm.add(upper);
    const hand = new THREE.Mesh(new THREE.SphereGeometry(0.05, 8, 6), SKIN);
    hand.position.set(0, -0.52, 0);
    arm.add(hand);
    return arm;
  };
  const leftArm = mkArm(-1);
  const rightArm = mkArm(1);
  group.add(leftArm, rightArm);
  // Field glasses: hidden until the pose raises them.
  const glasses = box(0.14, 0.05, 0.07, IRON, 0, 1.62, 0.16);
  glasses.visible = false;
  group.add(glasses);
  return { group, leftLeg, rightLeg, leftArm, rightArm, glasses };
}

/** Wilbur prone on the lower wing (the pilot the player embodies —
 * visible in chase/wingtip/Daniels views). Lies along +x, head fore. */
export function buildProneWilbur(): THREE.Group {
  const group = new THREE.Group();
  const torso = box(0.62, 0.16, 0.3, SUIT, 0, 0.08, 0);
  group.add(torso);
  const head = new THREE.Mesh(new THREE.SphereGeometry(0.1, 12, 10), SKIN);
  head.position.set(0.4, 0.12, 0);
  group.add(head);
  const cap = new THREE.Mesh(new THREE.CylinderGeometry(0.11, 0.115, 0.045, 12), SUIT);
  cap.position.set(0.4, 0.2, 0);
  group.add(cap);
  for (const side of [-1, 1]) {
    const leg = box(0.5, 0.11, 0.12, SUIT, -0.55, 0.06, side * 0.08);
    group.add(leg);
    const arm = box(0.34, 0.09, 0.09, SUIT, 0.28, 0.06, side * 0.2);
    arm.rotation.y = side * 0.25;
    group.add(arm);
  }
  return group;
}

/* ----------------------------- gulls ----------------------------- */

interface GullRig {
  group: THREE.Group;
  leftWing: THREE.Mesh;
  rightWing: THREE.Mesh;
  path: GullPath;
}

function buildGull(path: GullPath): GullRig {
  const group = new THREE.Group();
  const body = new THREE.Mesh(new THREE.CapsuleGeometry(0.09, 0.36, 3, 6), GULL_BODY);
  body.rotation.z = Math.PI / 2;
  group.add(body);
  const beak = new THREE.Mesh(
    new THREE.ConeGeometry(0.03, 0.12, 6),
    new THREE.MeshStandardMaterial({ color: 0xd9a13a, roughness: 0.7 }),
  );
  beak.rotation.z = -Math.PI / 2;
  beak.position.set(0.28, 0.02, 0);
  group.add(beak);
  // Wing geometry is baked flat in the xz-plane with its ROOT at the
  // body, span along ±z — so the mesh's rotation.x is a clean flap
  // hinge (tips beat up and down about the body line).
  const mkWing = (side: number): THREE.Mesh => {
    const geo = new THREE.PlaneGeometry(0.5, 0.95);
    geo.rotateX(-Math.PI / 2);
    geo.translate(0, 0, side * 0.475);
    return new THREE.Mesh(geo, GULL_WING);
  };
  const leftWing = mkWing(-1);
  const rightWing = mkWing(1);
  group.add(leftWing, rightWing);
  return { group, leftWing, rightWing, path };
}

/* ------------------------ assembled diorama ---------------------- */

export interface Dressing {
  group: THREE.Group;
  /** Advance every animated element to scene time t. */
  animate(
    t: number,
    orville: { onRail: boolean; aircraftX: number; releaseX: number | null; releaseT: number | null },
  ): void;
}

/** Build the full diorama around the launch point. `groundY` samples
 * terrain height for prop placement (launch-relative x/z in metres). */
export function buildDressing(
  launch: readonly [number, number, number],
  railLengthM: number,
  tileExtentM: number,
  withOcean: boolean,
  groundY: (xRel: number, zRel: number) => number,
): Dressing {
  const group = new THREE.Group();
  group.add(buildSky());
  const clouds = buildClouds();
  group.add(clouds);
  group.add(buildOuterGround(tileExtentM, withOcean));
  const rail = buildRail(railLengthM);
  rail.position.set(launch[0], launch[1], launch[2]);
  group.add(rail);
  const fire = buildCampfire();
  // Camp props, each seated on the sampled terrain.
  for (const p of campLayout()) {
    let obj: THREE.Object3D;
    switch (p.kind) {
      case "hangar":
        obj = buildBuilding(4.8, 12.5, 2.6);
        break;
      case "shack":
        obj = buildBuilding(4.2, 6.8, 2.4);
        break;
      case "campfire":
        obj = fire.group;
        break;
      case "chair":
        obj = buildChair();
        break;
      case "barrel":
        obj = buildBarrel();
        (obj as THREE.Mesh).position.y = 0.425;
        break;
      case "workbench":
        obj = buildWorkbench();
        break;
      case "toolchest":
        obj = buildToolchest();
        break;
    }
    obj.position.set(launch[0] + p.x, launch[1] + groundY(p.x, p.z), launch[2] + p.z);
    obj.rotation.y = p.rotY;
    group.add(obj);
  }
  const orvilleFig = buildFigure();
  group.add(orvilleFig.group);
  const fleet = gullFleet(14, 1903);
  const gulls = fleet.map((path) => {
    const rig = buildGull(path);
    group.add(rig.group);
    return rig;
  });
  return {
    group,
    animate(t, orville): void {
      // Clouds drift east — position is a pure function of t (frame
      // rate never changes the weather).
      for (const cloud of clouds.children) {
        const baseX = cloud.userData["baseX"] as number;
        const drift = cloud.userData["driftMps"] as number;
        const span = 2800;
        cloud.position.x = ((((baseX + drift * t + 1400) % span) + span) % span) - 1400;
      }
      // Campfire flicker (deterministic in t).
      const flick = 0.82 + 0.18 * Math.sin(11 * t) * Math.sin(17.3 * t + 1.2);
      fire.light.intensity = 14 * flick;
      fire.flame.scale.setScalar(0.9 + 0.2 * flick);
      // Orville.
      const pose = orvillePose(t, orville.onRail, orville.aircraftX, orville.releaseX, orville.releaseT);
      orvilleFig.group.position.set(
        launch[0] + pose.x,
        launch[1] + groundY(pose.x, pose.z),
        launch[2] + pose.z,
      );
      orvilleFig.group.rotation.y = pose.headingRad;
      const swing = pose.gaitRad === 0 ? 0 : 0.7 * Math.sin(pose.gaitRad);
      orvilleFig.leftLeg.rotation.x = swing;
      orvilleFig.rightLeg.rotation.x = -swing;
      if (pose.glassesUp) {
        orvilleFig.leftArm.rotation.x = -2.4;
        orvilleFig.rightArm.rotation.x = -2.4;
        orvilleFig.glasses.visible = true;
        // Face the machine he is watching: Ry(phi) sends +x to
        // (cos phi, -sin phi), so phi = atan2(-dz, dx) for target
        // direction (dx, dz) = (aircraftX - x, 0 - z).
        orvilleFig.group.rotation.y = Math.atan2(pose.z, orville.aircraftX - pose.x);
      } else {
        orvilleFig.leftArm.rotation.x = pose.gaitRad === 0 ? 0 : swing * 0.8;
        orvilleFig.rightArm.rotation.x = pose.gaitRad === 0 ? 0 : -swing * 0.8;
        orvilleFig.glasses.visible = false;
      }
      // Gulls. Ry(phi) sends body-forward +x to (cos phi, -sin phi) in
      // the xz-plane, so facing the heading direction needs phi = -h.
      for (const rig of gulls) {
        const p = gullPose(rig.path, t);
        rig.group.position.set(launch[0] + p.x, launch[1] + p.y, launch[2] + p.z);
        rig.group.rotation.y = -p.headingRad;
        rig.leftWing.rotation.x = p.flapRad;
        rig.rightWing.rotation.x = -p.flapRad;
      }
    },
  };
}
