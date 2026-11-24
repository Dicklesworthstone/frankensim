/**
 * Premixed gas flames on a burner's ports, all of one fixture's in one instanced draw.
 *
 * Each port carries the two cones a premixed hydrocarbon flame shows: a short bright inner cone (the reaction zone,
 * where the primary air burns) inside a faint violet-blue outer mantle (the diffusion zone, where secondary air
 * finishes it). Both lengthen with the gas flow the valve passes; at low flow the inner cone sits close to the port
 * and the mantle almost vanishes. Flicker is turbulence at the flame tips, not the base, which stays anchored.
 *
 * Flow is an input (a valve position), never an output: whatever burns the fuel is the kernel's business.
 */
import * as THREE from "three";
import type { V3 } from "./shop";

/** One port: its position, and the direction its jet leaves in (normalised by the set). */
export interface Port {
  at: V3;
  dir: V3;
  /** Which burner (valve) feeds it. */
  burner: number;
}

function flameGeometry(): THREE.BufferGeometry {
  // Inner cone and outer mantle, unit height along +y, base at 0; `aLayer` tells the shader which is which.
  const parts: THREE.BufferGeometry[] = [];
  const layers: number[] = [];
  for (const [layer, w, h] of [
    [0, 0.42, 0.55],
    [1, 1, 1],
  ] as const) {
    const pts: THREE.Vector2[] = [];
    for (let i = 0; i <= 7; i++) {
      const u = i / 7;
      // A teardrop: full width a little above the base, narrowing to a point.
      const r = w * Math.pow(Math.sin(Math.PI * Math.min(1, u * 0.92 + 0.08)), 0.8) * (1 - u * 0.55);
      pts.push(new THREE.Vector2(Math.max(0, r), u * h));
    }
    pts[pts.length - 1]!.x = 0;
    const g = new THREE.LatheGeometry(pts, 8);
    parts.push(g);
    layers.push(layer, g.attributes.position.count);
  }
  const merged = new THREE.BufferGeometry();
  const total = parts.reduce((n, g) => n + g.attributes.position.count, 0);
  const pos = new Float32Array(total * 3);
  const lay = new Float32Array(total);
  const idx: number[] = [];
  let off = 0;
  for (const [k, g] of parts.entries()) {
    pos.set(g.attributes.position.array as Float32Array, off * 3);
    lay.fill(k, off, off + g.attributes.position.count);
    const gi = g.index!.array;
    for (let i = 0; i < gi.length; i++) idx.push(gi[i]! + off);
    off += g.attributes.position.count;
  }
  merged.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  merged.setAttribute("aLayer", new THREE.BufferAttribute(lay, 1));
  merged.setIndex(idx);
  return merged;
}

const vert = /* glsl */ `
  attribute float aLayer;
  attribute float aFlow;
  attribute float aSeed;
  uniform float uTime;
  uniform float uSize;
  varying float vLayer;
  varying float vH;
  varying float vFlow;
  void main() {
    vLayer = aLayer;
    vFlow = aFlow;
    vec3 p = position;
    vH = p.y / mix(0.55, 1.0, aLayer);
    float on = step(0.02, aFlow);
    // The jet lengthens with flow; the mantle more than the inner cone.
    float len = uSize * mix(0.35, 1.0, aFlow) * mix(1.0, mix(0.4, 1.25, aFlow), aLayer);
    float wid = uSize * 0.32 * mix(0.75, 1.0, aFlow);
    float t = uTime * (11.0 + aSeed * 7.0) + aSeed * 40.0;
    float flick = 1.0 + 0.09 * sin(t) * vH + 0.07 * sin(t * 2.3 + 1.7) * vH * aLayer;
    p.xz *= wid;
    p.y *= len * flick * on;
    // Tips waver; bases stay on the port.
    p.x += sin(t * 0.7 + aSeed * 9.0) * 0.12 * wid * vH * vH;
    p.z += cos(t * 0.6 + aSeed * 5.0) * 0.12 * wid * vH * vH;
    gl_Position = projectionMatrix * modelViewMatrix * instanceMatrix * vec4(p, 1.0);
  }
`;
const frag = /* glsl */ `
  varying float vLayer;
  varying float vH;
  varying float vFlow;
  uniform float uGain;
  void main() {
    // Inner cone: bright blue-cyan, white-blue at its base. Mantle: faint violet-blue, gone at the tip.
    vec3 inner = mix(vec3(0.55, 0.8, 1.0), vec3(0.12, 0.38, 1.0), smoothstep(0.0, 0.9, vH));
    vec3 outer = mix(vec3(0.18, 0.32, 1.0), vec3(0.42, 0.18, 0.85), smoothstep(0.2, 1.0, vH));
    // A hint of sodium-yellow at the mantle tips when the burner is wide open.
    outer = mix(outer, vec3(1.0, 0.55, 0.18), smoothstep(0.75, 1.0, vH) * smoothstep(0.7, 1.0, vFlow) * 0.35);
    vec3 c = mix(inner, outer, vLayer);
    float a = mix(0.85, 0.32, vLayer) * (1.0 - smoothstep(0.55, 1.0, vH)) * smoothstep(0.0, 0.08, vH + 0.02);
    gl_FragColor = vec4(c * a * uGain, 1.0);
  }
`;

let geometry: THREE.BufferGeometry | null = null;

export class FlameSet {
  readonly mesh: THREE.InstancedMesh;
  private flow: THREE.InstancedBufferAttribute;
  private byBurner: number[][] = [];
  private mat: THREE.ShaderMaterial;

  /** `size` is the jet length at full flow, in metres (a 30 000 BTU/h star burner's is about 25 mm). */
  constructor(ports: readonly Port[], size = 0.024) {
    geometry ??= flameGeometry();
    this.mat = new THREE.ShaderMaterial({
      vertexShader: vert,
      fragmentShader: frag,
      uniforms: { uTime: { value: 0 }, uSize: { value: size }, uGain: { value: 1.25 } },
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    });
    const g = geometry.clone();
    this.flow = new THREE.InstancedBufferAttribute(new Float32Array(ports.length), 1);
    const seed = new THREE.InstancedBufferAttribute(new Float32Array(ports.length), 1);
    g.setAttribute("aFlow", this.flow);
    g.setAttribute("aSeed", seed);
    this.mesh = new THREE.InstancedMesh(g, this.mat, ports.length);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = 3;
    const m = new THREE.Matrix4();
    const q = new THREE.Quaternion();
    const up = new THREE.Vector3(0, 1, 0);
    const d = new THREE.Vector3();
    for (const [i, p] of ports.entries()) {
      d.set(p.dir[0], p.dir[1], p.dir[2]).normalize();
      q.setFromUnitVectors(up, d);
      m.compose(new THREE.Vector3(p.at[0], p.at[1], p.at[2]), q, new THREE.Vector3(1, 1, 1));
      this.mesh.setMatrixAt(i, m);
      seed.setX(i, ((i * 0.6180339887) % 1) + p.burner * 0.13);
      (this.byBurner[p.burner] ??= []).push(i);
    }
    this.mesh.instanceMatrix.needsUpdate = true;
  }
  /** A valve position, 0 shut to 1 wide open. */
  setFlow(burner: number, flow: number): void {
    for (const i of this.byBurner[burner] ?? []) this.flow.setX(i, flow);
    this.flow.needsUpdate = true;
  }
  setAll(flow: number): void {
    for (let b = 0; b < this.byBurner.length; b++) this.setFlow(b, flow);
  }
  tick(t: number): void {
    this.mat.uniforms.uTime!.value = t;
  }
  get burners(): number {
    return this.byBurner.length;
  }
}

/**
 * The ports of a round star or ring burner: `n` drilled round a head of radius `r`, at height `y` above (x, z),
 * angled up from the horizontal by `lift`.
 */
export function ringPorts(burner: number, at: V3, r: number, n: number, lift = 0.55, phase = 0): Port[] {
  const out: Port[] = [];
  for (let i = 0; i < n; i++) {
    const a = phase + (i / n) * Math.PI * 2;
    const [c, s] = [Math.cos(a), Math.sin(a)];
    out.push({ at: [at[0] + c * r, at[1], at[2] + s * r], dir: [c * Math.cos(lift), Math.sin(lift), s * Math.cos(lift)], burner });
  }
  return out;
}

/** Ports along a straight pipe burner (griddles, charbroilers, salamanders): along x from x0 to x1, jets up. */
export function linePorts(burner: number, y: number, z: number, x0: number, x1: number, n: number, tilt = 0.2): Port[] {
  const out: Port[] = [];
  for (let i = 0; i < n; i++) {
    const x = x0 + ((i + 0.5) / n) * (x1 - x0);
    const side = i % 2 ? 1 : -1;
    out.push({ at: [x, y, z + side * 0.004], dir: [0, 1, side * tilt], burner });
  }
  return out;
}
