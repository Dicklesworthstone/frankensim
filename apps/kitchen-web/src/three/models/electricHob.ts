/**
 * Model family: the open-coil electric hotplate. A countertop two-element unit, the elements one behind the other on
 * a stainless case: each a sheathed tubular element wound as a flat double spiral, the S-bend at its
 * centre and both cold ends leaving side by side for the terminal block, resting on a three-armed spider in a chrome
 * drip bowl, a polished trim ring round the opening. The pan stands on the coil itself: the coil's top is the cooking
 * plane, at the archetype's height.
 *
 * Archetype: NSF-HOTPLATE-ELEC-2B (about 305 x 610 x 125 mm: two 8 in, 203 mm, elements).
 * Front element centre at (w/2, h, -0.19), back element at (w/2, h, -0.45) in the fixture frame (see shop.ts).
 *
 * Joints: "knob_0" (front element), "knob_1" (back element), each turning about its own z.
 * Live names, `Hooks.set(name, v)`:
 *   heat0, heat1  0..1  element glow, as materials.ts setGlow reads it (0.2 is the first dull red, 1 yellow-white).
 *   knob0, knob1  0..1  the control's setting, off to high: turns the knob and lights the element's pilot lamp.
 * `hooks.glows` holds the two sheath materials, front first.
 */
import * as THREE from "three";
import { register, type Builder, type Hooks } from "./index";
import { glowMaterial, palette, setGlow } from "./materials";
import { bulletFoot, circleHole, cyl, knob, lathe, plate, rbox, screw, sweep, turnJoint, type Shop, type V3 } from "./shop";

const TAU = Math.PI * 2;

/** The element's centre line: a flat double spiral of pitch `pitch` out to `rOut`, its leads dipping to the block. */
function coilPath(rOut: number, pitch: number, y: number): THREE.Vector3[] {
  const k = pitch / Math.PI; // each arm's radius grows by 2 * pitch per turn; the other arm sits between
  const r0 = 0.012; // the S-bend at the centre: two semicircles of r0 / 2
  const d = 0.075; // half the angle between the two leads
  const theta = (rOut - r0) / k - Math.PI + d;
  const pts: [number, number, number][] = [];
  const at = (r: number, a: number, dy = 0): [number, number, number] => [r * Math.cos(a), y + dy, r * Math.sin(a)];
  // Lead B, from the terminal block in to the outer turn.
  const aB = theta - d;
  pts.push(at(0.132, aB, -0.024), at(0.118, aB, -0.021), at(0.109, aB, -0.008), at(0.104, aB, -0.001));
  // Arm B, outside in: the point reflection of arm A.
  const nB = Math.ceil((theta + Math.PI - d) / 0.07);
  for (let i = nB; i >= 0; i--) {
    const u = ((theta + Math.PI - d) * i) / nB;
    const [x, yy, z] = at(r0 + k * u, u);
    pts.push([-x, yy, -z]);
  }
  // The S: over the top round (-r0/2, 0), under the bottom round (r0/2, 0).
  for (let i = 1; i <= 8; i++) {
    const a = Math.PI - (Math.PI * i) / 8;
    pts.push([-r0 / 2 + (r0 / 2) * Math.cos(a), y, (r0 / 2) * Math.sin(a)]);
  }
  for (let i = 1; i < 8; i++) {
    const a = Math.PI + (Math.PI * i) / 8;
    pts.push([r0 / 2 + (r0 / 2) * Math.cos(a), y, (r0 / 2) * Math.sin(a)]);
  }
  // Arm A, inside out.
  const nA = Math.ceil((theta + d) / 0.07);
  for (let i = 0; i <= nA; i++) {
    const u = ((theta + d) * i) / nA;
    pts.push(at(r0 + k * u, u));
  }
  const aA = theta + d;
  pts.push(at(0.104, aA, -0.001), at(0.109, aA, -0.008), at(0.118, aA, -0.021), at(0.132, aA, -0.024));
  return pts.map(([x, yy, z]) => new THREE.Vector3(x, yy, z));
}

/** The angle midway between the two leads (in the x-z plane, from +x toward +z), so that the block can face back. */
function leadAngle(rOut: number, pitch: number): number {
  const k = pitch / Math.PI;
  return (rOut - 0.012) / k - Math.PI + 0.075;
}

interface Element {
  glow: THREE.MeshStandardMaterial;
}

/**
 * The sheath's section: a round tube rolled flat on top after winding, so that the pan sits on a band, not a line.
 * Counter-clockwise from the bottom, for `sweep`; y = 0 is the tube's axis, its flat top at 0.685 r.
 */
function sheathSection(r: number, n: number): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < n; i++) {
    const a = (i / n) * TAU - Math.PI / 2;
    let y = r * Math.sin(a);
    if (y > 0.55 * r) y = 0.55 * r + (y - 0.55 * r) * 0.3;
    out.push([r * 1.06 * Math.cos(a), y]);
  }
  return out;
}
const SHEATH_TOP = 0.685;

const glowTex = new Map<string, THREE.CanvasTexture>();
/**
 * How brightly each stretch of the element glows, along its length (u): the cold ends dark where the resistance wire
 * stops short of the terminals, a shade cooler where the coil lies on the spider's arms, and a few per cent of
 * winding variation. `pts` are the sheath's centre points, evenly spaced, in the element's plan; `arms` the spider's
 * arm angles.
 */
function glowProfile(key: string, pts: readonly THREE.Vector3[], arms: readonly number[], rOut: number): THREE.CanvasTexture {
  const hit = glowTex.get(key);
  if (hit) return hit;
  const W = 1024;
  const c = document.createElement("canvas");
  c.width = W;
  c.height = 2;
  const x = c.getContext("2d")!;
  const n = pts.length;
  // Distance along the sheath from either end to where it first reaches the outer turn's radius.
  const run: number[] = [0];
  for (let i = 1; i < n; i++) run.push(run[i - 1]! + pts[i]!.distanceTo(pts[i - 1]!));
  const total = run[n - 1]!;
  const rr = (i: number) => Math.hypot(pts[i]!.x, pts[i]!.z);
  let enterA = 0;
  while (enterA < n - 1 && (rr(enterA) > rOut + 0.004 || pts[enterA]!.y < -0.002)) enterA++;
  let enterB = n - 1;
  while (enterB > 0 && (rr(enterB) > rOut + 0.004 || pts[enterB]!.y < -0.002)) enterB--;
  const coldA = run[enterA]! + 0.03;
  const coldB = run[enterB]! - 0.03;
  let seed = 4111;
  const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  const wob = Array.from({ length: 24 }, () => rnd());
  for (let px = 0; px < W; px++) {
    const u = (px + 0.5) / W;
    const at = u * total;
    // The sample nearest this distance.
    let i = Math.min(n - 1, Math.max(0, Math.round(u * (n - 1))));
    while (i > 0 && run[i]! > at) i--;
    const p = pts[i]!;
    let v = 1;
    // Cold ends: dark at the terminals, warming over 30 mm into the coil.
    v *= Math.min(1, Math.max(0, (at - (coldA - 0.03)) / 0.03), Math.max(0, (coldB + 0.03 - at) / 0.03));
    // The spider's arms draw heat where the coil rests on them.
    const rho = Math.hypot(p.x, p.z);
    const al = Math.atan2(p.z, p.x);
    for (const a of arms) {
      const da = Math.atan2(Math.sin(al - a), Math.cos(al - a));
      if (Math.cos(da) <= 0) continue;
      const lat = Math.abs(rho * Math.sin(da));
      v *= 1 - 0.42 * Math.exp(-((lat / 0.006) ** 2));
    }
    // Winding variation: a slow ripple of a few per cent.
    const k = u * (wob.length - 1);
    const k0 = Math.floor(k);
    const fr = k - k0;
    v *= 0.95 + 0.08 * (wob[k0]! * (1 - fr) + wob[Math.min(wob.length - 1, k0 + 1)]! * fr);
    const lin = Math.max(0, Math.min(1, v));
    const g = Math.round(255 * Math.pow(lin, 1 / 2.2));
    x.fillStyle = `rgb(${g},${g},${g})`;
    x.fillRect(px, 0, 1, 2);
  }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.wrapS = THREE.ClampToEdgeWrapping;
  glowTex.set(key, t);
  return t;
}

/** One element in its bowl and ring, centred at (x, z); its coil's top at `top`, the cooktop sheet's top at `yT`. */
function element(s: Shop, x: number, z: number, top: number, yT: number): Element {
  const p = palette();
  const tubeR = 0.0042;
  const pitch = 0.0145;
  const rOut = 0.097;
  const yc = top - SHEATH_TOP * tubeR;
  // The sheath: oxidised nickel-chrome alloy, dark grey with a dull metallic sheen until it glows.
  const glow = glowMaterial(0x45413c);
  glow.metalness = 0.6;
  glow.roughness = 0.46;
  // Turn the spiral so that its leads run to the back (-z), where the terminal block sits under the sheet. A turn
  // about y by φ carries the plan angle a (from +x toward +z) to a - φ.
  const back = -Math.PI / 2;
  const arms = [0, 1, 2].map((i) => back + Math.PI / 3 + (i * TAU) / 3);
  const turn = new THREE.Matrix4().makeRotationY(leadAngle(rOut, pitch) + Math.PI / 2);
  const curve = new THREE.CatmullRomCurve3(coilPath(rOut, pitch, 0).map((v) => v.applyMatrix4(turn)), false, "centripetal");
  const N = s.seg(640, 120);
  const pts = curve.getSpacedPoints(N);
  const sec = sheathSection(tubeR, s.detail === 2 ? 14 : s.detail === 1 ? 7 : 6);
  const g = sweep(pts.map((v) => [v.x, v.y, v.z] as V3), sec);
  // The glow map runs along the sheath, u from end to end (the end caps, after the rings, are cold).
  const uv = g.attributes.uv as THREE.BufferAttribute;
  const ring = (N + 1) * sec.length;
  for (let i = 0; i < uv.count; i++) uv.setXY(i, i < ring ? Math.floor(i / sec.length) / N : 0, 0.5);
  glow.emissiveMap = glowProfile(`${N}`, pts, arms, rOut);
  // A glowing sheath is brightest where it faces the eye and dims toward its silhouette, where the oxide's
  // emissivity falls off: the glowing core of the tube against darker edges.
  glow.onBeforeCompile = (sh) => {
    sh.fragmentShader = sh.fragmentShader.replace(
      "#include <emissivemap_fragment>",
      "#include <emissivemap_fragment>\n\ttotalEmissiveRadiance *= mix(0.5, 1.0, pow(saturate(dot(normal, normalize(vViewPosition))), 0.6));",
    );
  };
  glow.customProgramCacheKey = () => "kit-sheath-core";
  s.add(glow, g, [x, yc, z]);
  // The cold ends' crimped collars where the sheath meets the terminal block, and the block itself (steatite).
  const [bx, bz] = [Math.cos(back), Math.sin(back)];
  s.add(p.ceramic, rbox(s, 0.034, 0.016, 0.018, 0.002), [x + bx * 0.142, yc - 0.024, z + bz * 0.142]);
  // The spider: three stamped arms on edge under the coil, joined at a small hub.
  for (const a of arms) {
    const [cx, cz] = [Math.cos(a), Math.sin(a)];
    const len = 0.092;
    s.add(p.ssHeat, rbox(s, len, 0.009, 0.0018, 0.0006), [x + cx * (0.012 + len / 2), yc - tubeR - 0.0045, z + cz * (0.012 + len / 2)], [0, -a, 0]);
    // The arm's outer foot, turned down onto the bowl.
    s.add(p.ssHeat, rbox(s, 0.006, 0.018, 0.0018, 0.0006), [x + cx * 0.103, yc - tubeR - 0.012, z + cz * 0.103], [0, -a, 0]);
  }
  s.add(p.ssHeat, cyl(s, 0.014, 0.0018, "y", 24), [x, yc - tubeR - 0.0095, z]);
  // The drip bowl: chrome-plated steel, a rolled lip round its centre hole, its flange under the trim ring.
  const bowl = new THREE.MeshStandardMaterial({ color: 0xc9c4bb, metalness: 1, roughness: 0.2, map: bowlTint() });
  bowl.userData.uv = "own";
  const t = 0.0008;
  const prof: [number, number][] = [
    [0.03, -0.031],
    [0.034, -0.0315],
    [0.06, -0.0305],
    [0.088, -0.0262],
    [0.104, -0.016],
    [0.111, -0.0045],
    [0.116, -0.0011],
    [0.123, -0.0011],
    [0.123, -0.0011 + t],
    [0.116, -0.0011 + t],
    [0.111 - t, -0.0045 + t * 0.4],
    [0.104 - t, -0.016 + t * 0.5],
    [0.088 - t * 0.4, -0.0262 + t],
    [0.06, -0.0305 + t],
    [0.034, -0.0315 + t],
    [0.0305, -0.031 + t * 1.6],
    [0.03, -0.031],
  ];
  s.add(bowl, lathe(s, prof, 72), [x, yT, z]);
  // The trim ring: polished, domed, sitting on the sheet round the opening.
  s.add(p.chrome, lathe(s, [
    [0.1112, 0.0002],
    [0.1345, 0.0002],
    [0.1348, 0.0011],
    [0.129, 0.0029],
    [0.119, 0.0043],
    [0.1135, 0.0041],
    [0.1112, 0.0026],
    [0.1112, 0.0002],
  ], 96), [x, yT, z]);
  return { glow };
}

let bowlTex: THREE.CanvasTexture | null = null;
/** Heat tint on the drip bowl: the plating straw and blue in rings where the element radiates hardest. */
function bowlTint(): THREE.CanvasTexture {
  if (bowlTex) return bowlTex;
  const c = document.createElement("canvas");
  c.width = 4;
  c.height = 128;
  const x = c.getContext("2d")!;
  const g = x.createLinearGradient(0, 0, 0, 128);
  // v runs along the lathe profile, from the centre hole out to the flange and back.
  g.addColorStop(0, "#9c8a78");
  g.addColorStop(0.12, "#b7a48a");
  g.addColorStop(0.25, "#c9b99d");
  g.addColorStop(0.34, "#a7b2c4");
  g.addColorStop(0.44, "#e8e6e2");
  g.addColorStop(0.56, "#f2f0ec");
  g.addColorStop(0.7, "#e8e6e2");
  g.addColorStop(0.8, "#b9b3ae");
  g.addColorStop(1, "#8e8378");
  x.fillStyle = g;
  x.fillRect(0, 0, 4, 128);
  bowlTex = new THREE.CanvasTexture(c);
  bowlTex.colorSpace = THREE.SRGBColorSpace;
  return bowlTex;
}

function lamp(): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({ color: 0x1a0d0b, metalness: 0, roughness: 0.35, emissive: 0xffa12e, emissiveIntensity: 0 });
  m.userData.uv = "own";
  return m;
}

const hotplate: Builder = (s, e) => {
  const p = palette();
  const { w, d, h } = e;
  const footH = 0.024;
  const yT = h - 0.0094; // the sheet's top; the trim rings stand 4.5 mm on it, the coils 5 mm above them
  const zs = [-0.19, -Math.min(0.45, d - 0.16)];
  // Feet, the case's sides, back and bottom, and the front fascia.
  for (const fx of [0.035, w - 0.035]) for (const fz of [-0.04, -d + 0.04]) bulletFoot(s, p.rubber, [fx, 0, fz], 0.013, footH);
  const bodyH = yT - 0.0012 - footH;
  const yMid = footH + bodyH / 2;
  // Sides the full depth; the back and the fascia fit between them, so that no two panels share a face.
  for (const sx of [0.0008, w - 0.0008]) s.add(p.ss, rbox(s, 0.0016, bodyH, d, 0.0006), [sx, yMid, -d / 2]);
  s.add(p.ss, rbox(s, w - 0.0032, bodyH, 0.0016, 0.0006), [w / 2, yMid, -d + 0.0008]);
  s.add(p.ssDull, rbox(s, w - 0.003, 0.0016, d - 0.003, 0.0006), [w / 2, footH + 0.0008, -d / 2]);
  const zf = -0.0008;
  s.add(p.ss, rbox(s, w - 0.0032, bodyH, 0.0016, 0.0006), [w / 2, yMid, zf]);
  // The cooktop sheet with its two openings, its front edge turned down over the fascia.
  const sheet = new THREE.Shape();
  sheet.moveTo(0, 0);
  sheet.lineTo(w, 0);
  sheet.lineTo(w, d);
  sheet.lineTo(0, d);
  sheet.lineTo(0, 0);
  for (const z of zs) sheet.holes.push(circleHole(w / 2, -z, 0.1125));
  const top = plate(s, sheet, 0.0012, 0, 64);
  top.rotateX(-Math.PI / 2);
  s.add(p.ss, top, [0, yT - 0.0012, 0]);
  // The fold: in front of the sheet's edge and the fascia, its top flush with the sheet's but beside it, not over it.
  s.add(p.ss, rbox(s, w, 0.012, 0.0014, 0.0006), [w / 2, yT - 0.006, 0.0007]);
  // Under the openings a dark well, so that the case's inside never shows through the bowls' centre holes.
  for (const z of zs) s.add(p.enamel, cyl(s, 0.04, 0.0015, "y", 24), [w / 2, yT - 0.05, z]);
  const els = zs.map((z) => element(s, w / 2, z, h, yT));
  // Controls: an infinite switch per element, its knob on a chrome bezel, a pilot lamp above it.
  const ky = footH + bodyH * 0.45;
  const knobs: string[] = [];
  const lamps: THREE.MeshStandardMaterial[] = [];
  for (const [i, kx] of [[0, w * 0.3], [1, w * 0.7]] as const) {
    const name = `knob_${i}`;
    const j = s.joint(name, [kx, ky, zf + 0.0008]);
    knob(j, p.bakelite, p.ssBright, [0, 0, 0], 0.017);
    knobs.push(name);
    const l = lamp();
    lamps.push(l);
    s.add(l, cyl(s, 0.0034, 0.002, "z", 12), [kx, ky + 0.03, zf + 0.0018]);
    if (s.fine(2)) {
      s.add(p.ssBright, cyl(s, 0.0046, 0.0012, "z", 16), [kx, ky + 0.03, zf + 0.0012]);
      // The dial's printed ticks, off at the top and round clockwise to high.
      for (let k = 0; k < 7; k++) {
        const a = Math.PI / 2 - 0.25 - (k / 6) * 1.1 * Math.PI;
        s.add(p.enamel, rbox(s, 0.0012, k === 0 || k === 6 ? 0.005 : 0.003, 0.0004, 0.0001), [kx + Math.cos(a) * 0.026, ky + Math.sin(a) * 0.026, zf + 0.0011], [0, 0, a - Math.PI / 2]);
      }
    }
  }
  if (s.fine(2)) for (const sx of [0.014, w - 0.014]) screw(s, p.ssBright, [sx, ky, zf + 0.0008], 0.0024);
  return (b) => {
    const hooks: Hooks = {
      glows: els.map((x) => x.glow),
      set: (name, v) => {
        const m = /^(heat|knob)(\d)$/.exec(name);
        if (!m) return;
        const i = Number(m[2]);
        if (m[1] === "heat") {
          const g = els[i]?.glow;
          if (g) setGlow(g, v);
        } else {
          const k = Math.min(1, Math.max(0, v));
          turnJoint(b.joints.get(knobs[i] ?? ""), "z", k <= 0.001 ? 0 : -(0.25 + k * 1.1 * Math.PI));
          const l = lamps[i];
          if (l) l.emissiveIntensity = k > 0.001 ? 2.4 : 0;
        }
      },
    };
    return hooks;
  };
};

register({
  "NSF-HOTPLATE-ELEC-2B": hotplate,
});
