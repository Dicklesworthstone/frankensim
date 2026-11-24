/**
 * CHAPTER 2 - "THE FLAME & THE PAN": a ribeye seared over gas, induction or an electric coil, and a wine deglaze
 * (#stove).
 *
 * Every visual quantity is read from `stepStoveKernel`'s state each frame (Rule 2: behaviour is computed, never put
 * in). The chapter describes the device, sets the kernel's inputs (source, pan, fat, throttle, the wine that lands in
 * the pan) and draws what comes back:
 *
 *   gas        NSF-RANGE-GAS-6B, framed on its front-left star burner: the valve flow is the throttle (an input), so
 *              the port flames and the knob follow it; their blue light under the pan follows the same flow.
 *   induction  NSF-INDUCTION-HOB-2Z: the knob shows the throttle; the zone's LEDs and lamp show the power the kernel
 *              couples into this pan (net_power_in_w over the zone's 3.5 kW rating), so aluminium barely lights it.
 *   electric   NSF-HOTPLATE-ELEC-2B (models/electricHob.ts): the coil glows at the temperature it must reach to
 *              radiate the kernel's coupled power into the pan bottom at the kernel's pan_bottom_c.
 *   pan        smallwares skillet() of the chosen alloy; the fat, fond, wine and steak sit on its real floor.
 *   fat        a film on the pan floor, its colour and shimmer from pan_surface_c (butter foams, browns, burns);
 *              smoke rises while the kernel reports smoking, denser the further past the smoke point.
 *   steam      emitted at the kernel's vapour rate: the steak's mass loss plus the deglaze's vapor_released_g.
 *   steak      ARCH-FOOD-STEAK-RIBEYE: browning = maillard_progress, shrink = 1 - cook_yield_ppm/1e6,
 *              core = core_temp_c; its raw thickness is the kernel's slab (mass over contact area and density).
 *   deglaze    NSF-DEGLAZE-BOTTLE tilts in and pours; the stream is ballistic parcels leaving the lip, and each
 *              parcel's mass enters deglaze_liquid_g only when it lands. The pool's depth is that mass over the
 *              floor; the fond fades with fond_dissolved; the boil and the steam burst follow the kernel.
 *
 * The kernel runs as a time-lapse (SIM_SPEED simulated seconds per wall second); the readouts and chart are its state.
 */

import * as THREE from "three";
import {
  FAT_PROPS,
  createInitialStoveState,
  stepStoveKernel,
  type FatKindId,
  type PanMaterialId,
  type StoveConfig,
  type StoveHeatSource,
  type StoveState,
} from "../sim/physicsKernels";
import { bindCanvasResize, createScrollSafeOrbit, createStudioRenderer, drawStripChart, observeInView } from "./kit";
import { buildFixture, type Ext, type Model } from "./models/all";
import "./models/electricHob";
import { palette } from "./models/materials";
import { Shop, rbox, cyl, turnJoint } from "./models/shop";
import { skillet, type SkilletInfo, type SkilletKind } from "./models/smallwares";
import { studioEnvironment, studioLights } from "./models/stage";

/* ------------------------------------------------------------------------------------------------------------ */
/* The described assembly                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/** Archetype extents (assets/kitchens/specs/equipment_dimensions.tsv). */
const RANGE: Ext = { w: 0.914, d: 0.864, h: 1.473 };
const HOB: Ext = { w: 0.4, d: 0.71, h: 0.145 };
const PLATE: Ext = { w: 0.305, d: 0.61, h: 0.125 };
const BOTTLE: Ext = { w: 0.076, d: 0.076, h: 0.295 };
/** The skillet: 11 in across the rim, so that a 165 mm ribeye leaves a margin of fat round it. */
const PAN_D = 0.28;
/** The kernel's steak slab: 0.0135 m² of contact at 1055 kg/m³ (stepStoveKernel's areaM2 and rhoMeat). */
const SLAB_AREA_M2 = 0.0135;
const MEAT_RHO = 1055;
/** The ribeye's cut-face box (ARCH-FOOD-STEAK-RIBEYE). */
const STEAK_W = 0.165;
const STEAK_D = 0.115;
const STEAK_YAW = 0.35;
/** The induction zone's rating, and the hotplate element's: what the kernel's coupled power is a fraction of. */
const ZONE_RATED_W = 3500;
/** Wine: 120 g per press of the button, from a 750 mL bottle; ~0.99 g/mL. */
const POUR_G = 120;
const BOTTLE_ML = 750;
const WINE_RHO = 990;
/** Simulated seconds per wall second: a time-lapse, so that a sear of several minutes plays out in about one. */
const SIM_SPEED = 6;

const SKILLET_KIND: Record<PanMaterialId, SkilletKind> = {
  cast_iron: "cast_iron",
  tri_ply: "tri_ply",
  stainless_304: "stainless_304",
  aluminum: "aluminum",
};

const clamp01 = (v: number) => Math.max(0, Math.min(1, v));
const smooth = (e0: number, e1: number, x: number) => {
  const t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3 - 2 * t);
};
const lerp = (a: number, b: number, t: number) => a + (b - a) * t;

function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Textures drawn once                                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

function valueNoise(seed: number, n = 64): (x: number, y: number) => number {
  const r = rng(seed);
  const v = new Float32Array(n * n);
  for (let i = 0; i < v.length; i++) v[i] = r();
  return (x, y) => {
    const xi = Math.floor(x);
    const yi = Math.floor(y);
    const fx = x - xi;
    const fy = y - yi;
    const ux = fx * fx * (3 - 2 * fx);
    const uy = fy * fy * (3 - 2 * fy);
    const x0 = ((xi % n) + n) % n;
    const y0 = ((yi % n) + n) % n;
    const x1 = (x0 + 1) % n;
    const y1 = (y0 + 1) % n;
    const a = v[y0 * n + x0]!;
    const b = v[y0 * n + x1]!;
    const c = v[y1 * n + x0]!;
    const d = v[y1 * n + x1]!;
    return a + (b - a) * ux + (c - a) * uy + (a - b - c + d) * ux * uy;
  };
}

/** A tileable normal map of slow convection cells and fine ripples, for the hot fat's surface. */
function shimmerNormals(): THREE.DataTexture {
  const S = 256;
  const n1 = valueNoise(71, 16);
  const n2 = valueNoise(72, 64);
  const h = new Float32Array(S * S);
  for (let j = 0; j < S; j++)
    for (let i = 0; i < S; i++) {
      const u = (i / S) * 16;
      const v = (j / S) * 16;
      const cell = 1 - Math.abs(2 * n1(u, v) - 1);
      h[j * S + i] = 0.65 * cell + 0.35 * n2((i / S) * 64, (j / S) * 64);
    }
  const data = new Uint8Array(S * S * 4);
  for (let j = 0; j < S; j++)
    for (let i = 0; i < S; i++) {
      const dx = h[j * S + ((i + 1) % S)]! - h[j * S + ((i + S - 1) % S)]!;
      const dy = h[((j + 1) % S) * S + i]! - h[((j + S - 1) % S) * S + i]!;
      const l = Math.hypot(dx * 4, dy * 4, 1);
      const k = (j * S + i) * 4;
      data[k] = Math.round(((-dx * 4) / l) * 127 + 128);
      data[k + 1] = Math.round(((-dy * 4) / l) * 127 + 128);
      data[k + 2] = Math.round((1 / l) * 127 + 128);
      data[k + 3] = 255;
    }
  const t = new THREE.DataTexture(data, S, S, THREE.RGBAFormat);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(3, 3);
  t.magFilter = THREE.LinearFilter;
  t.minFilter = THREE.LinearMipmapLinearFilter;
  t.generateMipmaps = true;
  t.needsUpdate = true;
  return t;
}

/** Fond: browned juices stuck to the floor round the steak's footprint, drawn in the steak's own frame. */
function fondTexture(r: number): THREE.CanvasTexture {
  const S = 512;
  const c = document.createElement("canvas");
  c.width = c.height = S;
  const x = c.getContext("2d")!;
  const rand = rng(8831);
  const nz = valueNoise(8832, 64);
  const px = S / (2 * r);
  const a = STEAK_W / 2;
  const b = STEAK_D / 2;
  for (let i = 0; i < 9000; i++) {
    const th = rand() * Math.PI * 2;
    // Most of it in a band just outside the edge, where the juices run out and brown; a little scattered further.
    const k = rand() < 0.8 ? 0.92 + Math.pow(rand(), 1.6) * 0.38 : 0.9 + rand() * 0.9;
    const ex = Math.cos(th) * a * k;
    const ez = Math.sin(th) * b * k;
    if (Math.hypot(ex, ez) > r * 0.97) continue;
    const u = S / 2 + ex * px;
    const v = S / 2 - ez * px;
    const dens = nz(u / 18, v / 18);
    if (rand() > 0.35 + dens * 0.65) continue;
    const s = (0.4 + rand() * rand() * 3.2) * (S / 512);
    const l = 18 + rand() * 26;
    x.fillStyle = `hsla(${(18 + rand() * 14).toFixed(0)},${(55 + rand() * 25).toFixed(0)}%,${l.toFixed(0)}%,${(0.35 + rand() * 0.55).toFixed(2)})`;
    x.beginPath();
    x.ellipse(u, v, s * (1 + rand()), s, rand() * Math.PI, 0, Math.PI * 2);
    x.fill();
  }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = 8;
  return t;
}

/** A soft, uneven puff for steam and smoke: white, its alpha a noisy radial falloff. */
function puffTexture(): THREE.CanvasTexture {
  const S = 128;
  const c = document.createElement("canvas");
  c.width = c.height = S;
  const x = c.getContext("2d")!;
  const img = x.createImageData(S, S);
  const n = valueNoise(991, 32);
  for (let j = 0; j < S; j++)
    for (let i = 0; i < S; i++) {
      const u = (i + 0.5) / S - 0.5;
      const v = (j + 0.5) / S - 0.5;
      const r = Math.hypot(u, v) * 2;
      const f = 0.5 * n(i / 10, j / 10) + 0.3 * n(i / 5 + 7, j / 5) + 0.2 * n(i / 2.5, j / 2.5 + 3);
      const a = Math.pow(Math.max(0, 1 - r), 1.8) * (0.35 + 0.9 * f);
      const k = (j * S + i) * 4;
      img.data[k] = img.data[k + 1] = img.data[k + 2] = 255;
      img.data[k + 3] = Math.round(Math.min(1, a) * 255);
    }
  x.putImageData(img, 0, 0);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

/** The backdrop: a dark warm studio sweep, brightest just above the horizon. */
function backdrop(): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = 4;
  c.height = 256;
  const x = c.getContext("2d")!;
  const g = x.createLinearGradient(0, 0, 0, 256);
  g.addColorStop(0, "#0b0807");
  g.addColorStop(0.45, "#1a120d");
  g.addColorStop(0.72, "#2a1c13");
  g.addColorStop(1, "#140e0a");
  x.fillStyle = g;
  x.fillRect(0, 0, 4, 256);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Steam and smoke: camera-facing puffs in one instanced draw                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

const enum Puff {
  Steam = 0,
  Smoke = 1,
}

class Wisps {
  readonly mesh: THREE.Mesh;
  private geo: THREE.InstancedBufferGeometry;
  private aPos: THREE.InstancedBufferAttribute;
  private aData: THREE.InstancedBufferAttribute;
  private n: number;
  private live = 0;
  private p: Float32Array;
  private v: Float32Array;
  private age: Float32Array;
  private life: Float32Array;
  private s0: Float32Array;
  private s1: Float32Array;
  private kind: Uint8Array;
  private seed: Float32Array;
  private peak: Float32Array;
  private rand = rng(5521);

  constructor(n: number) {
    this.n = n;
    const plane = new THREE.PlaneGeometry(1, 1);
    this.geo = new THREE.InstancedBufferGeometry();
    this.geo.index = plane.index;
    this.geo.setAttribute("position", plane.attributes.position!);
    this.geo.setAttribute("uv", plane.attributes.uv!);
    this.aPos = new THREE.InstancedBufferAttribute(new Float32Array(n * 3), 3);
    this.aData = new THREE.InstancedBufferAttribute(new Float32Array(n * 4), 4);
    this.aPos.setUsage(THREE.DynamicDrawUsage);
    this.aData.setUsage(THREE.DynamicDrawUsage);
    this.geo.setAttribute("aPos", this.aPos);
    this.geo.setAttribute("aData", this.aData);
    this.geo.instanceCount = 0;
    const mat = new THREE.ShaderMaterial({
      uniforms: { uMap: { value: puffTexture() } },
      vertexShader: /* glsl */ `
        attribute vec3 aPos;
        attribute vec4 aData; // size, alpha, spin, kind
        varying vec2 vUv;
        varying float vAlpha;
        varying float vKind;
        void main() {
          vUv = uv;
          vAlpha = aData.y;
          vKind = aData.w;
          vec4 mv = modelViewMatrix * vec4(aPos, 1.0);
          float c = cos(aData.z);
          float s = sin(aData.z);
          // Drawn out upward: the plume stretches as it rises.
          mv.xy += mat2(c, s, -s, c) * position.xy * aData.x * vec2(1.0, 1.45);
          gl_Position = projectionMatrix * mv;
        }`,
      fragmentShader: /* glsl */ `
        uniform sampler2D uMap;
        varying vec2 vUv;
        varying float vAlpha;
        varying float vKind;
        void main() {
          float a = texture2D(uMap, vUv).a * vAlpha;
          // Steam is condensed droplets, white; oil smoke scatters blue and shows grey-blue.
          vec3 c = mix(vec3(0.93, 0.92, 0.9), vec3(0.5, 0.55, 0.64), vKind);
          gl_FragColor = vec4(c, a);
          #include <colorspace_fragment>
        }`,
      transparent: true,
      depthWrite: false,
    });
    this.mesh = new THREE.Mesh(this.geo, mat);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = 5;
    this.p = new Float32Array(n * 3);
    this.v = new Float32Array(n * 3);
    this.age = new Float32Array(n);
    this.life = new Float32Array(n);
    this.s0 = new Float32Array(n);
    this.s1 = new Float32Array(n);
    this.kind = new Uint8Array(n);
    this.seed = new Float32Array(n);
    this.peak = new Float32Array(n);
  }

  /** Release one puff at (x, y, z) rising at `up` m/s. */
  emit(kind: Puff, x: number, y: number, z: number, up: number, life: number, s0: number, s1: number, peak: number): void {
    if (this.live >= this.n) return;
    const i = this.live++;
    const r = this.rand;
    this.p[i * 3] = x;
    this.p[i * 3 + 1] = y;
    this.p[i * 3 + 2] = z;
    this.v[i * 3] = (r() - 0.5) * 0.03;
    this.v[i * 3 + 1] = up * (0.75 + 0.5 * r());
    this.v[i * 3 + 2] = (r() - 0.5) * 0.03;
    this.age[i] = 0;
    this.life[i] = life * (0.8 + 0.4 * r());
    this.s0[i] = s0;
    this.s1[i] = s1 * (0.8 + 0.4 * r());
    this.kind[i] = kind;
    this.seed[i] = r() * 100;
    this.peak[i] = peak;
  }

  clear(): void {
    this.live = 0;
    this.geo.instanceCount = 0;
  }

  update(dt: number, t: number): void {
    let k = 0;
    for (let i = 0; i < this.live; i++) {
      const age = this.age[i]! + dt;
      if (age >= this.life[i]!) continue;
      // Compact the survivors into the first k slots.
      if (k !== i) {
        for (let c = 0; c < 3; c++) {
          this.p[k * 3 + c] = this.p[i * 3 + c]!;
          this.v[k * 3 + c] = this.v[i * 3 + c]!;
        }
        this.life[k] = this.life[i]!;
        this.s0[k] = this.s0[i]!;
        this.s1[k] = this.s1[i]!;
        this.kind[k] = this.kind[i]!;
        this.seed[k] = this.seed[i]!;
        this.peak[k] = this.peak[i]!;
      }
      this.age[k] = age;
      const sd = this.seed[k]!;
      const y = this.p[k * 3 + 1]!;
      // A buoyant plume: it slows as it mixes, and eddies push it about, more the higher it has risen.
      const sw = this.kind[k] === Puff.Smoke ? 0.11 : 0.07;
      this.v[k * 3] = this.v[k * 3]! * Math.exp(-dt * 1.2) + Math.sin(t * 1.7 + sd + y * 31) * sw * dt;
      this.v[k * 3 + 1] = this.v[k * 3 + 1]! * Math.exp(-dt * 0.45);
      this.v[k * 3 + 2] = this.v[k * 3 + 2]! * Math.exp(-dt * 1.2) + Math.cos(t * 1.3 + sd * 1.3 + y * 27) * sw * dt;
      for (let c = 0; c < 3; c++) this.p[k * 3 + c] = this.p[k * 3 + c]! + this.v[k * 3 + c]! * dt;
      const q = age / this.life[k]!;
      this.aPos.setXYZ(k, this.p[k * 3]!, this.p[k * 3 + 1]!, this.p[k * 3 + 2]!);
      const alpha = this.peak[k]! * smooth(0, 0.18, q) * (1 - smooth(0.45, 1, q));
      this.aData.setXYZW(k, lerp(this.s0[k]!, this.s1[k]!, Math.sqrt(q)), alpha, sd + t * (0.25 + (sd % 1) * 0.3) * (sd % 2 > 1 ? 1 : -1), this.kind[k]!);
      k++;
    }
    this.live = k;
    this.geo.instanceCount = k;
    if (k) {
      this.aPos.needsUpdate = true;
      this.aData.needsUpdate = true;
    }
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Sizzle and boil bubbles                                                                                        */
/* ------------------------------------------------------------------------------------------------------------ */

class Bubbles {
  readonly mesh: THREE.InstancedMesh;
  private n: number;
  private live = 0;
  private x: Float32Array;
  private y: Float32Array;
  private z: Float32Array;
  private age: Float32Array;
  private life: Float32Array;
  private r: Float32Array;
  private m = new THREE.Matrix4();
  private q = new THREE.Quaternion();
  private s = new THREE.Vector3();
  private v = new THREE.Vector3();

  constructor(n: number, mat: THREE.Material) {
    this.n = n;
    this.mesh = new THREE.InstancedMesh(new THREE.SphereGeometry(1, 12, 6, 0, Math.PI * 2, 0, Math.PI / 2), mat, n);
    this.mesh.count = 0;
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = 4;
    this.x = new Float32Array(n);
    this.y = new Float32Array(n);
    this.z = new Float32Array(n);
    this.age = new Float32Array(n);
    this.life = new Float32Array(n);
    this.r = new Float32Array(n);
  }
  emit(x: number, y: number, z: number, r: number, life: number): void {
    if (this.live >= this.n) return;
    const i = this.live++;
    this.x[i] = x;
    this.y[i] = y;
    this.z[i] = z;
    this.r[i] = r;
    this.age[i] = 0;
    this.life[i] = life;
  }
  clear(): void {
    this.live = 0;
    this.mesh.count = 0;
  }
  update(dt: number): void {
    let k = 0;
    for (let i = 0; i < this.live; i++) {
      const age = this.age[i]! + dt;
      if (age >= this.life[i]!) continue;
      this.x[k] = this.x[i]!;
      this.y[k] = this.y[i]!;
      this.z[k] = this.z[i]!;
      this.r[k] = this.r[i]!;
      this.life[k] = this.life[i]!;
      this.age[k] = age;
      // Grows as it fills, then bursts.
      const q = age / this.life[k]!;
      const rr = this.r[k]! * Math.sqrt(smooth(0, 0.8, q)) * (q > 0.94 ? 0.2 : 1);
      this.m.compose(this.v.set(this.x[k]!, this.y[k]!, this.z[k]!), this.q, this.s.set(rr, rr * 0.8, rr));
      this.mesh.setMatrixAt(k, this.m);
      k++;
    }
    this.live = k;
    this.mesh.count = k;
    if (k) this.mesh.instanceMatrix.needsUpdate = true;
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The wine stream: ballistic parcels leaving the bottle's lip                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

const G = 9.81;
const STREAM_RINGS = 40;
const STREAM_SIDES = 12;

interface Parcel {
  x: number;
  y: number;
  z: number;
  vx: number;
  vy: number;
  vz: number;
  /** Volume flow when it left the lip, m³/s; its mass, g. */
  q: number;
  g: number;
}

class Stream {
  readonly mesh: THREE.Mesh;
  parcels: Parcel[] = [];
  private pos: THREE.BufferAttribute;
  private nrm: THREE.BufferAttribute;
  private curve = new THREE.CatmullRomCurve3([new THREE.Vector3(), new THREE.Vector3()], false, "centripetal");
  private tmp = new THREE.Vector3();
  private tan = new THREE.Vector3();
  private nA = new THREE.Vector3();
  private nB = new THREE.Vector3();

  constructor(mat: THREE.Material) {
    const g = new THREE.BufferGeometry();
    this.pos = new THREE.BufferAttribute(new Float32Array(STREAM_RINGS * STREAM_SIDES * 3), 3);
    this.nrm = new THREE.BufferAttribute(new Float32Array(STREAM_RINGS * STREAM_SIDES * 3), 3);
    this.pos.setUsage(THREE.DynamicDrawUsage);
    this.nrm.setUsage(THREE.DynamicDrawUsage);
    g.setAttribute("position", this.pos);
    g.setAttribute("normal", this.nrm);
    const idx: number[] = [];
    for (let i = 0; i < STREAM_RINGS - 1; i++)
      for (let j = 0; j < STREAM_SIDES; j++) {
        const a = i * STREAM_SIDES + j;
        const b = i * STREAM_SIDES + ((j + 1) % STREAM_SIDES);
        const c = a + STREAM_SIDES;
        const d = b + STREAM_SIDES;
        // Counter-clockwise seen from outside: the ring runs anticlockwise about the tangent.
        idx.push(a, b, c, b, d, c);
      }
    g.setIndex(idx);
    this.mesh = new THREE.Mesh(g, mat);
    this.mesh.frustumCulled = false;
    this.mesh.visible = false;
    this.mesh.castShadow = false;
    this.mesh.renderOrder = 3;
  }

  /**
   * Advance every parcel ballistically; those that fall through the surface at `surfaceY` land, and `land` is told
   * where and how much. Then rebuild the tube through the parcels, its radius from continuity (r = √(q / π|v|)).
   */
  update(dt: number, surfaceY: number, land: (x: number, z: number, g: number) => void, t: number): void {
    const kept: Parcel[] = [];
    for (const p of this.parcels) {
      p.vy -= G * dt;
      p.x += p.vx * dt;
      p.y += p.vy * dt;
      p.z += p.vz * dt;
      if (p.y <= surfaceY) land(p.x, p.z, p.g);
      else kept.push(p);
    }
    this.parcels = kept;
    const n = kept.length;
    if (n < 2) {
      this.mesh.visible = false;
      return;
    }
    // Newest (at the lip) first; finish at the oldest's landing point, found from its own trajectory.
    const pts: THREE.Vector3[] = [];
    const rad: number[] = [];
    for (let i = n - 1; i >= 0; i--) {
      const p = kept[i]!;
      pts.push(new THREE.Vector3(p.x, p.y, p.z));
      rad.push(Math.max(0.0006, Math.sqrt(p.q / (Math.PI * Math.max(0.05, Math.hypot(p.vx, p.vy, p.vz))))));
    }
    const o = kept[0]!;
    const disc = o.vy * o.vy + 2 * G * (o.y - surfaceY);
    const tl = (o.vy + Math.sqrt(Math.max(0, disc))) / G;
    pts.push(new THREE.Vector3(o.x + o.vx * tl, surfaceY, o.z + o.vz * tl));
    const vEnd = Math.hypot(o.vx, o.vy - G * tl, o.vz);
    rad.push(Math.max(0.0006, Math.sqrt(o.q / (Math.PI * Math.max(0.05, vEnd)))));
    this.curve.points = pts;
    const P = this.pos.array as Float32Array;
    const N = this.nrm.array as Float32Array;
    const segs = pts.length - 1;
    for (let i = 0; i < STREAM_RINGS; i++) {
      const u = i / (STREAM_RINGS - 1);
      this.curve.getPoint(u, this.tmp);
      this.curve.getTangent(u, this.tan);
      const f = u * segs;
      const k = Math.min(segs - 1, Math.floor(f));
      let r = lerp(rad[k]!, rad[k + 1]!, f - k);
      // Capillary varicose waves grow as the column thins and lengthens.
      r *= 1 + 0.12 * u * Math.sin(u * 46 - t * 38);
      // A frame round the tangent.
      const ref = Math.abs(this.tan.y) < 0.95 ? this.nA.set(0, 1, 0) : this.nA.set(1, 0, 0);
      this.nB.crossVectors(this.tan, ref).normalize();
      this.nA.crossVectors(this.nB, this.tan).normalize();
      for (let j = 0; j < STREAM_SIDES; j++) {
        const a = (j / STREAM_SIDES) * Math.PI * 2;
        const cx = Math.cos(a);
        const sx = Math.sin(a);
        const nx = this.nA.x * cx + this.nB.x * sx;
        const ny = this.nA.y * cx + this.nB.y * sx;
        const nz = this.nA.z * cx + this.nB.z * sx;
        const v = (i * STREAM_SIDES + j) * 3;
        P[v] = this.tmp.x + nx * r;
        P[v + 1] = this.tmp.y + ny * r;
        P[v + 2] = this.tmp.z + nz * r;
        N[v] = nx;
        N[v + 1] = ny;
        N[v + 2] = nz;
      }
    }
    this.pos.needsUpdate = true;
    this.nrm.needsUpdate = true;
    this.mesh.geometry.computeBoundingSphere();
    this.mesh.visible = true;
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The fat's look at temperature                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

interface FatLook {
  /** The fat's own colour (sRGB). */
  color: THREE.Color;
  /** How much of that colour hides the pan metal under it: a clear film 0, foam or burnt solids toward 1. */
  body: number;
  rough: number;
}
type Stop = readonly [number, number, number, number, number, number]; // °C, r, g, b (sRGB), body, roughness

function ramp(stops: readonly Stop[], tc: number, out: FatLook): void {
  let a = stops[0]!;
  let b = stops[stops.length - 1]!;
  for (let i = 1; i < stops.length; i++)
    if (tc <= stops[i]![0]) {
      a = stops[i - 1]!;
      b = stops[i]!;
      break;
    }
  const t = b[0] === a[0] ? 1 : clamp01((tc - a[0]) / (b[0] - a[0]));
  out.color.setRGB(lerp(a[1], b[1], t), lerp(a[2], b[2], t), lerp(a[3], b[3], t), THREE.SRGBColorSpace);
  out.body = lerp(a[4], b[4], t);
  out.rough = lerp(a[5], b[5], t);
}

/**
 * The film's colour at the pan surface temperature. Butter is about 16 % water: it foams until that has boiled off,
 * then its milk solids brown (beurre noisette) and, past the smoke point, burn. Canola and lard are refined fats that
 * barely colour until they smoke.
 */
function fatLook(fat: FatKindId, tc: number, smokeC: number, out: FatLook): void {
  if (fat === "butter")
    ramp([
      [20, 0.97, 0.86, 0.46, 0.95, 0.5],
      [100, 0.96, 0.87, 0.52, 0.92, 0.58],
      [122, 0.9, 0.7, 0.26, 0.4, 0.12],
      [150, 0.72, 0.45, 0.14, 0.5, 0.08],
      [smokeC, 0.46, 0.25, 0.08, 0.62, 0.08],
      [smokeC + 45, 0.16, 0.08, 0.03, 0.85, 0.14],
    ], tc, out);
  else if (fat === "lard")
    ramp([
      [20, 0.96, 0.93, 0.84, 0.04, 0.06],
      [smokeC, 0.92, 0.84, 0.62, 0.06, 0.04],
      [smokeC + 60, 0.62, 0.45, 0.2, 0.2, 0.06],
    ], tc, out);
  else
    ramp([
      [20, 0.92, 0.74, 0.28, 0.06, 0.05],
      [smokeC, 0.86, 0.6, 0.18, 0.08, 0.04],
      [smokeC + 60, 0.55, 0.32, 0.08, 0.25, 0.06],
    ], tc, out);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The chapter                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

interface SourceRig {
  root: THREE.Group;
  model: Model;
  /** Height of the counter or table the bottle stands on, relative to the cooking plane. */
  counterY: number;
}

export function mountStoveLab(root: HTMLElement): void {
  const stageEl = root.querySelector<HTMLElement>("[data-stove-stage]")!;
  const canvas = root.querySelector<HTMLCanvasElement>("[data-stove-canvas]")!;
  const chartCanvas = root.querySelector<HTMLCanvasElement>("[data-stove-chart]");

  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x1a120d, 1.1, 3.2);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.5;
  // The cooking plane is y = 0 and the burner, zone or element under the pan is at the origin.
  const lights = studioLights(scene, new THREE.Vector3(0, 0.04, 0), 0.5, 1024);
  lights.key.intensity = 0.55;
  lights.key.castShadow = false;
  lights.key.position.set(0.34, 0.62, 0.42);
  lights.hemi.intensity = 0.22;
  lights.rim.intensity = 1.25;
  lights.rim.position.set(-0.55, 0.38, -0.6);
  // The key: a warm spot from high front right, pooled on the pan so that the room round it falls away dark.
  const spot = new THREE.SpotLight(0xffd9b0, 6.5, 3, 0.42, 0.75, 1.6);
  spot.position.set(0.55, 0.95, 0.62);
  spot.target.position.set(-0.01, 0, -0.01);
  spot.castShadow = true;
  spot.shadow.mapSize.set(2048, 2048);
  spot.shadow.camera.near = 0.4;
  spot.shadow.camera.far = 2.4;
  spot.shadow.bias = -0.00015;
  spot.shadow.normalBias = 0.0015;
  spot.shadow.radius = 4;
  scene.add(spot, spot.target);

  const camera = new THREE.PerspectiveCamera(32, 16 / 10, 0.01, 12);
  const TARGET = new THREE.Vector3(0, 0.004, 0.01);
  const AZ = 0.38;
  const EL = 0.31;
  const VIEW_DIR = new THREE.Vector3(Math.sin(AZ) * Math.cos(EL), Math.sin(EL), Math.cos(AZ) * Math.cos(EL)).normalize();
  camera.position.copy(TARGET).addScaledVector(VIEW_DIR, 0.6);
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], {
    minDistance: 0.14,
    maxDistance: 1.7,
    maxPolarAngle: Math.PI * 0.485,
  });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  bindCanvasResize(stageEl, renderer, camera);
  /**
   * Frame the pan across the view, whatever the stage's shape; while the bottle pours (`wide` toward 1) ease back and
   * up so that the bottle and its stream are in the picture. Once the viewer orbits, the camera is theirs.
   */
  let wide = 0;
  const target = new THREE.Vector3();
  const frame = () => {
    if (userMoved) return;
    const halfV = THREE.MathUtils.degToRad(camera.fov / 2);
    const halfH = Math.atan(Math.tan(halfV) * camera.aspect);
    const w = smooth(0, 1, wide);
    const dist = Math.max(0.18 / Math.tan(halfH), 0.11 / Math.tan(halfV)) * (1 + 0.3 * w);
    target.copy(TARGET).add(new THREE.Vector3(-0.035 * w, 0.075 * w, -0.02 * w));
    controls.target.copy(target);
    camera.position.copy(target).addScaledVector(VIEW_DIR, dist);
    controls.update();
  };
  new ResizeObserver(frame).observe(stageEl);

  const p = palette();

  /* ---- Heat sources, each built the first time it is chosen ---- */
  const rigs = new Map<StoveHeatSource, SourceRig>();
  const counterShop = (top: number, x0: number, x1: number, z0: number, z1: number, splash: boolean) => {
    const s = new Shop(1, "counter");
    const t = 0.04;
    s.add(p.ss, rbox(s, x1 - x0, t, z1 - z0, 0.002), [(x0 + x1) / 2, top - t / 2, (z0 + z1) / 2]);
    s.add(p.ss, cyl(s, 0.012, x1 - x0, "x", 24), [(x0 + x1) / 2, top - 0.012, z1]);
    if (splash) s.add(p.ss, rbox(s, x1 - x0, 0.26, 0.0016, 0.0006), [(x0 + x1) / 2, top + 0.13, z0 + 0.001]);
    return s.build(true).group;
  };
  function rig(src: StoveHeatSource): SourceRig {
    const hit = rigs.get(src);
    if (hit) return hit;
    const group = new THREE.Group();
    let model: Model;
    let counterY: number;
    if (src === "gas") {
      // The front-left burner of a 36 in range: centre at (0.18, -0.22), its grates' tops at 0.914 m.
      model = buildFixture("NSF-RANGE-GAS-6B", RANGE, 2);
      model.group.position.set(-0.18, -0.914, 0.22);
      model.hooks.flames?.setAll(0);
      model.hooks.flames?.setFlow(model.hooks.burners ?? 6, 0.12);
      counterY = 0;
      // A work table at cooking height against the range's left side.
      group.add(counterShop(0, -0.95, -0.183, -0.62, 0.2, false));
    } else if (src === "induction") {
      // The front zone of the countertop hob, centre (0.2, -0.2) on its glass at 0.145 m.
      model = buildFixture("NSF-INDUCTION-HOB-2Z", HOB, 2);
      model.group.position.set(-HOB.w / 2, -HOB.h, 0.2);
      counterY = -HOB.h;
      group.add(counterShop(counterY, -0.95, 0.95, -0.62, 0.3, false));
    } else {
      model = buildFixture("NSF-HOTPLATE-ELEC-2B", PLATE, 2);
      model.group.position.set(-PLATE.w / 2, -PLATE.h, 0.19);
      counterY = -PLATE.h;
      group.add(counterShop(counterY, -0.95, 0.95, -0.62, 0.3, false));
    }
    group.add(model.group);
    // Only what is near the pan casts into the key's shadow map; the rest of the appliance just receives.
    group.updateMatrixWorld(true);
    const near = new THREE.Box3(new THREE.Vector3(-0.32, -0.3, -0.32), new THREE.Vector3(0.32, 0.4, 0.32));
    const bb = new THREE.Box3();
    group.traverse((o) => {
      if (!(o as THREE.Mesh).isMesh || !o.castShadow) return;
      bb.setFromObject(o);
      if (!bb.intersectsBox(near)) o.castShadow = false;
    });
    group.visible = false;
    scene.add(group);
    const r = { root: group, model, counterY };
    rigs.set(src, r);
    return r;
  }

  // The flame's light on the grate and the pan's underside, and the coil's.
  const flameLight = new THREE.PointLight(0x4d7bff, 0, 0.45, 2);
  flameLight.position.set(0, -0.022, 0);
  scene.add(flameLight);
  const coilLight = new THREE.PointLight(0xff5a1c, 0, 0.5, 2);
  coilLight.position.set(0, -0.012, 0);
  scene.add(coilLight);

  /* ---- The pan, the fat, the fond, the wine, the steak ---- */
  const pans = new Map<PanMaterialId, { group: THREE.Group; info: SkilletInfo }>();
  const panRoot = new THREE.Group();
  scene.add(panRoot);
  function pan(id: PanMaterialId): { group: THREE.Group; info: SkilletInfo } {
    const hit = pans.get(id);
    if (hit) return hit;
    const s = new Shop(2, `skillet-${id}`);
    const info = skillet(s, SKILLET_KIND[id], PAN_D);
    const group = s.build(true).group;
    // Handle to the front left, away from the range's other burners.
    group.rotation.y = -0.62;
    group.visible = false;
    panRoot.add(group);
    const out = { group, info };
    pans.set(id, out);
    return out;
  }

  const FILM = 0.0009; // the fat film's depth on the floor
  const fatLookNow: FatLook = { color: new THREE.Color(), body: 0, rough: 0.1 };
  const shimmer = shimmerNormals();
  // A film this thin is seen as the pan's own floor gone wet: the metal under it, tinted by the fat, its gloss the
  // fat's surface. Foam and browned solids hide the metal by their `body`.
  const oilMat = new THREE.MeshPhysicalMaterial({
    color: 0x222018,
    metalness: 0.5,
    roughness: 0.06,
    clearcoat: 1,
    clearcoatRoughness: 0.035,
    clearcoatNormalMap: shimmer,
    clearcoatNormalScale: new THREE.Vector2(0.1, 0.1),
    envMapIntensity: 1.3,
  });
  const oil = new THREE.Mesh(new THREE.CircleGeometry(1, 120).rotateX(-Math.PI / 2), oilMat);
  oil.receiveShadow = true;
  scene.add(oil);
  /** The floor each pan's film lies on (sRGB colour, metalness): seasoned iron, polished clad steel, bare aluminium. */
  const FLOOR: Record<PanMaterialId, readonly [number, number, number, number]> = {
    cast_iron: [0.13, 0.125, 0.12, 0.45],
    tri_ply: [0.8, 0.81, 0.82, 1],
    stainless_304: [0.78, 0.79, 0.8, 1],
    aluminum: [0.8, 0.81, 0.82, 0.92],
  };
  const floorCol = new THREE.Color();
  const tintCol = new THREE.Color();

  const fondMat = new THREE.MeshStandardMaterial({
    map: fondTexture(0.13),
    transparent: true,
    opacity: 0,
    roughness: 0.75,
    metalness: 0,
    depthWrite: false,
    polygonOffset: true,
    polygonOffsetFactor: -1,
    polygonOffsetUnits: -1,
  });
  const fond = new THREE.Mesh(new THREE.CircleGeometry(0.13, 96).rotateX(-Math.PI / 2), fondMat);
  fond.receiveShadow = true;
  fond.renderOrder = 1;

  const wineMat = new THREE.MeshPhysicalMaterial({
    color: 0x5c0c18,
    metalness: 0,
    roughness: 0.04,
    transparent: true,
    opacity: 0,
    clearcoat: 1,
    clearcoatRoughness: 0.03,
    normalMap: shimmer,
    normalScale: new THREE.Vector2(0.25, 0.25),
    depthWrite: false,
    envMapIntensity: 1.5,
  });
  const pool = new THREE.Mesh(new THREE.CircleGeometry(1, 120).rotateX(-Math.PI / 2), wineMat);
  pool.renderOrder = 3;
  pool.visible = false;
  scene.add(pool);

  const steakHolder = new THREE.Group();
  steakHolder.rotation.y = STEAK_YAW;
  scene.add(steakHolder);
  steakHolder.add(fond);

  const cfg: StoveConfig = {
    source: "induction",
    pan: "cast_iron",
    fat: "canola",
    powerFrac: 0.72,
    steakMassG: 340,
    ambientC: 23,
  };
  // The raw steak's thickness is the kernel's slab: mass over contact area and density.
  const slabT = cfg.steakMassG / 1000 / (MEAT_RHO * SLAB_AREA_M2);
  const steak = buildFixture("ARCH-FOOD-STEAK-RIBEYE", { w: STEAK_W, d: STEAK_D, h: slabT * (1 + 0.35 * 0.3) }, 2);
  steak.group.position.set(-STEAK_W / 2, 0, STEAK_D / 2);
  steakHolder.add(steak.group);

  /* ---- The bottle and its pour ---- */
  const bottle = buildFixture("NSF-DEGLAZE-BOTTLE", BOTTLE, 2);
  const BOTTLE_MID = 0.12;
  const LIP = new THREE.Vector3(0, BOTTLE.h - BOTTLE_MID, 0);
  bottle.group.position.set(-BOTTLE.w / 2, -BOTTLE_MID, BOTTLE.d / 2);
  const bottlePivot = new THREE.Group();
  bottlePivot.add(bottle.group);
  scene.add(bottlePivot);
  let bottleMl = 675;
  bottle.hooks.set?.("level", bottleMl / BOTTLE_ML);

  // Red wine in a thin column: dark ruby, glowing where the light comes through it, a wet gloss.
  const streamMat = new THREE.MeshPhysicalMaterial({
    color: 0x4a0712,
    emissive: 0x5a0a18,
    emissiveIntensity: 0.55,
    metalness: 0,
    roughness: 0.05,
    clearcoat: 0.7,
    clearcoatRoughness: 0.04,
    transparent: true,
    opacity: 0.94,
    envMapIntensity: 0.7,
  });
  const stream = new Stream(streamMat);
  scene.add(stream.mesh);
  const ripples: { mesh: THREE.Mesh; age: number }[] = [];
  const rippleGeo = new THREE.RingGeometry(0.82, 1, 48).rotateX(-Math.PI / 2);
  for (let i = 0; i < 6; i++) {
    const m = new THREE.Mesh(rippleGeo, new THREE.MeshBasicMaterial({ color: 0xe8a6a0, transparent: true, opacity: 0, depthWrite: false, blending: THREE.AdditiveBlending }));
    m.visible = false;
    m.renderOrder = 4;
    scene.add(m);
    ripples.push({ mesh: m, age: 9 });
  }
  let rippleNext = 0;
  let rippleClock = 0;

  const wisps = new Wisps(360);
  scene.add(wisps.mesh);
  const bubbleMat = new THREE.MeshPhysicalMaterial({ color: 0xfff4dc, metalness: 0, roughness: 0.03, transparent: true, opacity: 0.45, clearcoat: 1, depthWrite: false, envMapIntensity: 1.8 });
  const bubbles = new Bubbles(160, bubbleMat);
  scene.add(bubbles.mesh);

  /* ---- Kernel state and history ---- */
  let state: StoveState = createInitialStoveState(cfg);
  const histPan: number[] = [];
  const histCore: number[] = [];
  function prewarmHistory(): void {
    histPan.length = 0;
    histCore.length = 0;
    histPan.push(state.pan_surface_c);
    histCore.push(state.core_temp_c);
    for (let s = 0; s < 16; s += 1) {
      stepStoveKernel(state, cfg, 1.5);
      histPan.push(state.pan_surface_c);
      histCore.push(state.core_temp_c);
    }
  }
  prewarmHistory();

  /* ---- The scene follows the configuration ---- */
  let floorY = 0;
  let floorR = 0.11;
  function syncSource(): void {
    for (const [k, r] of rigs) r.root.visible = k === cfg.source;
    rig(cfg.source).root.visible = true;
  }
  function syncPan(): void {
    for (const [k, v] of pans) v.group.visible = k === cfg.pan;
    const cur = pan(cfg.pan);
    cur.group.visible = true;
    floorY = cur.info.floorY;
    // The fat spreads into the floor-to-wall bend; the wall hides the film's edge.
    floorR = cur.info.floorR;
    const rFilm = floorR + 0.006;
    oil.scale.set(rFilm, 1, rFilm);
    oil.position.y = floorY + FILM;
    pool.scale.set(rFilm + 0.004, 1, rFilm + 0.004);
    steakHolder.position.y = floorY + FILM * 0.5;
    fond.position.y = -FILM * 0.5 + 0.0002;
    fond.scale.setScalar(Math.min(1, (floorR + 0.004) / 0.13));
  }

  /* ---- DOM ---- */
  const srcBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-stove-source]"));
  const panBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-stove-pan]"));
  const fatBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-stove-fat]"));
  const powerInput = root.querySelector<HTMLInputElement>("[data-stove-power]");
  const powerVal = root.querySelector<HTMLElement>("[data-stove-power-val]");
  const deglazeBtn = root.querySelector<HTMLButtonElement>("[data-stove-deglaze]");
  const resetBtn = root.querySelector<HTMLButtonElement>("[data-stove-reset]");

  const outPan = root.querySelector<HTMLElement>("[data-st-pan]");
  const outCore = root.querySelector<HTMLElement>("[data-st-core]");
  const outMaillard = root.querySelector<HTMLElement>("[data-st-maillard]");
  const outFond = root.querySelector<HTMLElement>("[data-st-fond]");
  const outVapor = root.querySelector<HTMLElement>("[data-st-vapor]");
  const outPower = root.querySelector<HTMLElement>("[data-st-power]");
  const smokeBadge = root.querySelector<HTMLElement>("[data-st-smoke]");

  function syncStoveReadouts(): void {
    if (outPan) outPan.textContent = `${state.pan_surface_c.toFixed(1)} °C`;
    if (outCore) outCore.textContent = `${state.core_temp_c.toFixed(1)} °C`;
    if (outMaillard) outMaillard.textContent = `${(state.maillard_progress * 100).toFixed(1)}%`;
    if (outFond) outFond.textContent = `${(state.fond_dissolved * 100).toFixed(1)}%`;
    if (outVapor) outVapor.textContent = `${state.vapor_released_g.toFixed(1)} g`;
    if (outPower) outPower.textContent = `${Math.round(state.net_power_in_w)} W`;
    if (smokeBadge) {
      const fat = FAT_PROPS[cfg.fat];
      if (state.smoking) {
        const fatName: Record<string, string> = { butter: "butter", canola: "canola oil", lard: "lard" };
        smokeBadge.textContent = `Smoking: the ${fatName[cfg.fat] ?? fat.label.split(" ")[0]!.toLowerCase()} is past its ${fat.smokeC.toFixed(0)} °C smoke point`;
        smokeBadge.dataset.state = "warn";
      } else {
        smokeBadge.textContent = `Clean sear · smoke point ${fat.smokeC.toFixed(0)} °C`;
        smokeBadge.dataset.state = "ok";
      }
    }
    if (chartCanvas) {
      drawStripChart(
        chartCanvas,
        [
          { label: "Pan °C", color: "#f97316", values: histPan },
          { label: "Core °C", color: "#38bdf8", values: histCore },
        ],
        {
          minY: 10,
          maxY: 310,
          unit: "°C",
          threshold: { value: state.smoke_point_c, color: "#fb7185", label: `Smoke ${state.smoke_point_c}°C` },
        },
      );
    }
  }
  function pushHistory(): void {
    histPan.push(state.pan_surface_c);
    histCore.push(state.core_temp_c);
    if (histPan.length > 90) {
      histPan.shift();
      histCore.shift();
    }
  }

  for (const btn of srcBtns) {
    btn.addEventListener("click", () => {
      cfg.source = btn.dataset.stoveSource as StoveHeatSource;
      for (const b of srcBtns) b.setAttribute("aria-pressed", String(b.dataset.stoveSource === cfg.source));
      syncSource();
      syncStoveReadouts();
    });
  }
  for (const btn of panBtns) {
    btn.addEventListener("click", () => {
      cfg.pan = btn.dataset.stovePan as PanMaterialId;
      for (const b of panBtns) b.setAttribute("aria-pressed", String(b.dataset.stovePan === cfg.pan));
      syncPan();
      syncStoveReadouts();
    });
  }
  for (const btn of fatBtns) {
    btn.addEventListener("click", () => {
      cfg.fat = btn.dataset.stoveFat as FatKindId;
      for (const b of fatBtns) b.setAttribute("aria-pressed", String(b.dataset.stoveFat === cfg.fat));
      syncStoveReadouts();
    });
  }
  powerInput?.addEventListener("input", () => {
    cfg.powerFrac = Number(powerInput.value) / 100;
    if (powerVal) powerVal.textContent = `${Math.round(cfg.powerFrac * 100)}%`;
  });

  /* ---- The pour: a choreographed hand (an input), the wine itself ballistic ---- */
  type PourPhase = "rest" | "lift" | "pour" | "back";
  let phase: PourPhase = "rest";
  let phaseT = 0;
  let pourTarget = 0; // grams still to pour
  const TILT = 2.0; // rad from upright at the pour
  const TILT_ON = 1.62; // the neck dips below the wine's surface: flow starts
  const Q_MAX_G = 52; // g/s, a steady glugging pour from a 750 mL neck
  const JET_AREA = 0.4 * Math.PI * 0.0094 * 0.0094; // the part of the 18.8 mm bore the wine leaves through
  const POUR_DIR = new THREE.Vector3(1, 0, 0.18).normalize();
  const IMPACT = new THREE.Vector3(-0.058, 0, -0.06);
  const restPose = { pos: new THREE.Vector3(), quat: new THREE.Quaternion() };
  const pourPose = { pos: new THREE.Vector3(), quat: new THREE.Quaternion() };
  const tiltAxis = new THREE.Vector3(POUR_DIR.z, 0, -POUR_DIR.x);
  const lipNow = new THREE.Vector3();
  const lipPrev = new THREE.Vector3();
  const axisNow = new THREE.Vector3();
  let haveLip = false;
  let tilt = 0;

  /** Where the pour must hold the bottle for its stream to land at IMPACT: solved from the nominal jet's flight. */
  function solvePourPose(): void {
    const v0 = Q_MAX_G / 1e6 / (WINE_RHO / 1000) / JET_AREA;
    const ax = Math.sin(TILT);
    const ay = Math.cos(TILT);
    const surface = floorY + FILM;
    const lipH = 0.115 - surface; // the lip 11.5 cm over the pan's floor, clear of the rim
    // Falling lipH from a downward start of v0·|ay|: lipH = v0|ay| t + g t²/2.
    const vy0 = -v0 * ay;
    const tf = (-vy0 + Math.sqrt(vy0 * vy0 + 2 * G * lipH)) / G;
    const run = v0 * ax * tf;
    const lip = IMPACT.clone().addScaledVector(POUR_DIR, -run);
    lip.y = surface + lipH;
    pourPose.quat.setFromAxisAngle(tiltAxis, TILT);
    pourPose.pos.copy(lip).sub(LIP.clone().applyQuaternion(pourPose.quat));
  }
  function setRest(): void {
    const counter = rig(cfg.source).counterY;
    // Standing on the counter behind the pan, to the left (for gas, on the side table beside the range).
    restPose.pos.set(-0.235, counter + BOTTLE_MID, -0.23);
    restPose.quat.identity();
  }
  const qTmp = new THREE.Quaternion();
  function poseBottle(dt: number): void {
    setRest();
    solvePourPose();
    let u = 0;
    if (phase === "lift") u = smooth(0, 1, phaseT / 1.1);
    else if (phase === "pour") u = 1;
    else if (phase === "back") u = 1 - smooth(0, 1, phaseT / 1.0);
    // Up and over: rise first, then reach across; tilt lags the reach.
    const lift = Math.sin(u * Math.PI) * 0.07;
    bottlePivot.position.lerpVectors(restPose.pos, pourPose.pos, smooth(0.1, 1, u));
    bottlePivot.position.y = lerp(restPose.pos.y, pourPose.pos.y, smooth(0, 0.7, u)) + lift;
    const tu = smooth(0.25, 1, u);
    tilt = TILT * tu;
    qTmp.slerpQuaternions(restPose.quat, pourPose.quat, tu);
    bottlePivot.quaternion.copy(qTmp);
    bottlePivot.updateMatrixWorld(true);
    lipPrev.copy(lipNow);
    lipNow.copy(LIP).applyMatrix4(bottlePivot.matrixWorld);
    if (!haveLip) {
      lipPrev.copy(lipNow);
      haveLip = true;
    }
    axisNow.set(0, 1, 0).applyQuaternion(bottlePivot.quaternion);
    void dt;
  }

  /** Emit this frame's wine at the lip: flow rises as the neck dips under the surface, glugs, stops at the target. */
  function pourFlow(dt: number, t: number): void {
    const open = smooth(TILT_ON, TILT, tilt);
    if (open <= 0 || pourTarget <= 0 || bottleMl <= 0) return;
    const glug = 0.82 + 0.18 * Math.sin(t * 2 * Math.PI * 5.5);
    const g = Math.min(pourTarget, bottleMl * (WINE_RHO / 1000), Q_MAX_G * open * glug * dt);
    if (g <= 0) return;
    pourTarget -= g;
    bottleMl -= g / (WINE_RHO / 1000);
    const q = g / 1e6 / (WINE_RHO / 1000) / dt;
    const v0 = q / JET_AREA;
    const lv = lipNow.clone().sub(lipPrev).divideScalar(Math.max(1e-3, dt));
    stream.parcels.push({
      x: lipNow.x,
      y: lipNow.y,
      z: lipNow.z,
      vx: axisNow.x * v0 + lv.x * 0.6,
      vy: axisNow.y * v0 + lv.y * 0.6,
      vz: axisNow.z * v0 + lv.z * 0.6,
      q,
      g,
    });
    bottle.hooks.set?.("level", bottleMl / BOTTLE_ML);
  }

  deglazeBtn?.addEventListener("click", () => {
    // A fresh bottle when this one cannot give a full measure.
    if (phase === "rest" && bottleMl < POUR_G / (WINE_RHO / 1000) + 1) {
      bottleMl = BOTTLE_ML;
      bottle.hooks.set?.("level", 1);
    }
    pourTarget += POUR_G;
    if (phase === "rest" || phase === "back") {
      phase = "lift";
      phaseT = 0;
    }
  });
  function resetAll(): void {
    state = createInitialStoveState(cfg);
    phase = "rest";
    phaseT = 0;
    pourTarget = 0;
    stream.parcels = [];
    stream.mesh.visible = false;
    bottleMl = 675;
    bottle.hooks.set?.("level", bottleMl / BOTTLE_ML);
    wisps.clear();
    bubbles.clear();
    histPan.length = 0;
    histCore.length = 0;
    histPan.push(state.pan_surface_c);
    histCore.push(state.core_temp_c);
    lastMass = state.steak_mass_g;
    lastVapor = state.vapor_released_g;
    syncStoveReadouts();
  }
  resetBtn?.addEventListener("click", resetAll);

  syncSource();
  syncPan();
  syncStoveReadouts();
  frame();

  let inView = true;
  observeInView(stageEl, (visible) => {
    inView = visible;
  });

  /* ---- Per-frame coupling ---- */
  let lastTs = performance.now();
  let tAnim = 0;
  let chartAccum = 0;
  let lastMass = state.steak_mass_g;
  let lastVapor = state.vapor_released_g;
  let searRate = 0; // g/s of water leaving the steak (kernel mass loss), smoothed for display
  let boilRate = 0; // g/s of deglaze liquid boiling off (kernel vapor_released_g)
  let shownBrowning = -1;
  let emitSteam = 0;
  let emitBoil = 0;
  let emitSmoke = 0;
  let emitSizzle = 0;
  let emitBubble = 0;
  const rand = rng(1207);
  const v3 = new THREE.Vector3();

  /** A point on the floor in the ring between the steak's edge (scaled by k) and the pan wall. */
  function floorPoint(out: THREE.Vector3, shrink: number, edgeOnly: boolean): THREE.Vector3 {
    const s = 1 - 0.45 * shrink;
    const a = (STEAK_W / 2) * s;
    const b = (STEAK_D / 2) * s;
    for (let tries = 0; tries < 8; tries++) {
      const th = rand() * Math.PI * 2;
      const k = edgeOnly ? 1.02 + rand() * 0.08 : 1.05 + rand() * 0.9;
      out.set(Math.cos(th) * a * k, 0, Math.sin(th) * b * k);
      // The steak's frame to the pan's: the holder's turn about y.
      const c = Math.cos(STEAK_YAW);
      const sn = Math.sin(STEAK_YAW);
      out.set(out.x * c + out.z * sn, 0, -out.x * sn + out.z * c);
      if (Math.hypot(out.x, out.z) < floorR - 0.004) return out;
    }
    return out.set(floorR * 0.8, 0, 0);
  }

  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - lastTs) / 1000));
    lastTs = now;
    if (!inView || document.hidden) return;
    tAnim += dt;

    // The pour lands its wine in the pan before the kernel steps, so the kernel sees it this frame.
    phaseT += dt;
    if (phase === "lift" && phaseT >= 1.1) {
      phase = "pour";
      phaseT = 0;
    } else if (phase === "pour" && (pourTarget <= 0 || bottleMl <= 0)) {
      phase = "back";
      phaseT = 0;
    } else if (phase === "back" && phaseT >= 1.0) {
      phase = "rest";
      phaseT = 0;
    }
    poseBottle(dt);
    pourFlow(dt, tAnim);
    const wideTo = phase === "rest" ? 0 : 1;
    if (wide !== wideTo) {
      wide = wideTo > wide ? Math.min(1, wide + dt / 0.9) : Math.max(0, wide - dt / 1.4);
      frame();
    }
    const surface = floorY + FILM + Math.max(0, state.deglaze_liquid_g) / 1000 / WINE_RHO / Math.max(1e-4, Math.PI * floorR * floorR - SLAB_AREA_M2);
    stream.update(dt, surface, (x, z, g) => {
      state.deglaze_liquid_g += g;
      rippleClock += 1;
      if (rippleClock % 3 === 0) {
        const r = ripples[rippleNext++ % ripples.length]!;
        r.age = 0;
        r.mesh.position.set(x, surface + 0.0004, z);
        r.mesh.visible = true;
      }
    }, tAnim);

    const simDt = dt * SIM_SPEED;
    stepStoveKernel(state, cfg, simDt);
    chartAccum += dt;
    if (chartAccum >= 0.25) {
      chartAccum = 0;
      pushHistory();
      syncStoveReadouts();
    }

    // Vapour rates, from the kernel's own bookkeeping: the steak's lost mass and the deglaze's boiled-off mass.
    const dm = Math.max(0, lastMass - state.steak_mass_g);
    const dv = Math.max(0, state.vapor_released_g - lastVapor);
    lastMass = state.steak_mass_g;
    lastVapor = state.vapor_released_g;
    const k = 1 - Math.exp(-dt / 0.25);
    searRate += (dm / Math.max(1e-6, simDt) - searRate) * k;
    boilRate += (dv / Math.max(1e-6, simDt) - boilRate) * k;

    /* Heat source */
    const r = rig(cfg.source);
    const hooks = r.model.hooks;
    if (cfg.source === "gas") {
      hooks.flames?.setFlow(0, cfg.powerFrac);
      hooks.tick?.(tAnim, dt, cfg.powerFrac);
      flameLight.intensity = 0.05 * cfg.powerFrac * (0.92 + 0.08 * Math.sin(tAnim * 23));
    } else flameLight.intensity = 0;
    if (cfg.source === "induction") {
      hooks.set?.("zone0", state.net_power_in_w / ZONE_RATED_W);
      // The knob is the cook's setting; setZone turned it to the coupled power, so turn it back to the setting.
      turnJoint(r.model.joints.get("knob_0"), "z", -(0.25 + cfg.powerFrac * 1.1 * Math.PI));
    }
    if (cfg.source === "electric") {
      // The sheath temperature that radiates the coupled power into the pan bottom: εσA(Te⁴ - Tpb⁴) = P, with the
      // element's upper surface (≈1.75 m of 8.4 mm tube, half facing up) and an oxidised-sheath / pan emissivity 0.8.
      const sigmaEA = 5.670e-8 * 0.8 * (Math.PI * 0.0084 * 1.75 * 0.5);
      const tpb = state.pan_bottom_c + 273.15;
      const te = Math.pow(tpb ** 4 + Math.max(0, state.net_power_in_w) / sigmaEA, 0.25) - 273.15;
      // Visible from the Draper point (~525 °C): dull red, cherry near 750 °C, orange near 900 °C, yellow by 1300 °C.
      const heat = te < 525 ? 0.2 * clamp01((te - 400) / 125) : 0.2 + 0.8 * clamp01((te - 525) / 800);
      hooks.set?.("heat0", heat);
      hooks.set?.("knob0", cfg.powerFrac);
      coilLight.intensity = 0.06 * Math.pow(clamp01((heat - 0.2) / 0.8), 1.5);
    } else coilLight.intensity = 0;

    /* Fat */
    const tc = state.pan_surface_c;
    const smokeC = FAT_PROPS[cfg.fat].smokeC;
    fatLook(cfg.fat, tc, smokeC, fatLookNow);
    const fl = FLOOR[cfg.pan];
    floorCol.setRGB(fl[0], fl[1], fl[2], THREE.SRGBColorSpace);
    // Light crosses the film twice: the metal's reflection takes the fat's colour, softened.
    tintCol.setRGB(1, 1, 1).lerp(fatLookNow.color, 0.7);
    floorCol.multiply(tintCol);
    oilMat.color.copy(floorCol).lerp(fatLookNow.color, fatLookNow.body);
    oilMat.metalness = fl[3] * (1 - fatLookNow.body);
    oilMat.roughness = lerp(0.07, fatLookNow.rough, fatLookNow.body);
    oilMat.clearcoat = 1 - 0.6 * smooth(0.6, 1, fatLookNow.body);
    // Convection cells in the film strengthen and quicken as it heats and thins.
    const sh = smooth(60, smokeC + 20, tc);
    oilMat.clearcoatNormalScale.setScalar(0.05 + 0.55 * sh);
    shimmer.offset.x = (shimmer.offset.x + dt * (0.004 + 0.03 * sh)) % 1;
    shimmer.offset.y = (shimmer.offset.y + dt * (0.003 + 0.022 * sh)) % 1;

    /* Fond and wine */
    fondMat.opacity = clamp01(state.maillard_progress * 1.6) * (1 - state.fond_dissolved) * 0.95;
    const liquid = Math.max(0, state.deglaze_liquid_g);
    pool.visible = liquid > 0.05;
    if (pool.visible) {
      pool.position.y = surface;
      const depth = surface - floorY - FILM;
      // Red wine absorbs strongly: a film shows the pan through it, a few millimetres hide the floor.
      wineMat.opacity = 0.45 + 0.5 * smooth(0, 0.0012, depth);
      // Dark garnet, turning to a brown pan sauce as the fond dissolves into it and as it reduces.
      const brown = clamp01(state.fond_dissolved * 0.8 + (1 - smooth(0, 60, liquid)) * 0.35);
      wineMat.color.setRGB(lerp(0.2, 0.2, brown), lerp(0.025, 0.085, brown), lerp(0.045, 0.03, brown), THREE.SRGBColorSpace);
      wineMat.normalScale.setScalar(0.1 + 1.2 * clamp01(boilRate / 1.5));
    }
    for (const rp of ripples) {
      if (rp.age > 0.6) {
        rp.mesh.visible = false;
        continue;
      }
      rp.age += dt;
      const q = rp.age / 0.6;
      rp.mesh.scale.setScalar(0.004 + 0.03 * Math.sqrt(q));
      (rp.mesh.material as THREE.MeshBasicMaterial).opacity = 0.32 * (1 - q);
    }

    /* Steak */
    const shrink = 1 - state.cook_yield_ppm / 1e6;
    if (Math.abs(state.maillard_progress - shownBrowning) >= 0.01) {
      shownBrowning = state.maillard_progress;
      steak.hooks.set?.("browning", shownBrowning);
    }
    steak.hooks.set?.("core", state.core_temp_c);
    steak.hooks.set?.("shrink", shrink);

    /* Steam, smoke, sizzle, boil */
    const shrinkS = 1 - 0.45 * Math.min(0.3, shrink);
    // Sizzle steam from the steak's edges, where the water driven out of the seared face boils off.
    emitSteam += 9 * Math.pow(searRate, 0.8) * simDt;
    while (emitSteam >= 1) {
      emitSteam -= 1;
      floorPoint(v3, Math.min(0.3, shrink), true);
      wisps.emit(Puff.Steam, v3.x, floorY + 0.004 + rand() * slabT * shrinkS, v3.z, 0.15, 1.8, 0.01, 0.05, 0.14);
    }
    // The deglaze: a burst off the whole pool, as fast as the kernel boils it.
    emitBoil += 16 * Math.pow(boilRate, 0.8) * simDt;
    while (emitBoil >= 1) {
      emitBoil -= 1;
      floorPoint(v3, Math.min(0.3, shrink), false);
      wisps.emit(Puff.Steam, v3.x, surface + 0.002, v3.z, 0.26, 1.9, 0.014, 0.085, 0.11);
    }
    if (state.smoking) {
      const over = clamp01((tc - state.smoke_point_c) / 35);
      emitSmoke += (8 + 40 * over) * dt;
      while (emitSmoke >= 1) {
        emitSmoke -= 1;
        floorPoint(v3, Math.min(0.3, shrink), false);
        wisps.emit(Puff.Smoke, v3.x, floorY + 0.003, v3.z, 0.09, 3.0, 0.012, 0.12, 0.22 + 0.25 * over);
      }
    }
    // Frying bubbles at the steak's edge, at the rate its water boils into the fat; boiling bubbles in the wine.
    emitSizzle += 60 * searRate * simDt;
    while (emitSizzle >= 1) {
      emitSizzle -= 1;
      floorPoint(v3, Math.min(0.3, shrink), true);
      bubbles.emit(v3.x, floorY + FILM, v3.z, 0.0007 + rand() * 0.0012, 0.12 + rand() * 0.3);
    }
    if (pool.visible) {
      emitBubble += 12 * boilRate * simDt;
      while (emitBubble >= 1) {
        emitBubble -= 1;
        floorPoint(v3, Math.min(0.3, shrink), false);
        bubbles.emit(v3.x, surface, v3.z, 0.001 + rand() * 0.002, 0.15 + rand() * 0.35);
      }
    }
    wisps.update(dt, tAnim);
    bubbles.update(dt);

    controls.update();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
}
