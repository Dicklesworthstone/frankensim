/**
 * The hot line: ranges, burners and grates. Dimensions are the archetypes' (assets/kitchens/specs), the layouts those
 * of the equipment class. A 36 in restaurant range is six open star burners under three cast-iron grates (one per
 * front-and-back pair), a stainless manifold panel of knobs over a standard oven, a flue riser and high shelf behind,
 * all on 6 in adjustable legs.
 *
 * NSF-RANGE-GAS-6B, as built here, front to back:
 *   legs       6 in stainless tube legs on mounting plates, threaded studs into adjustable bullet feet;
 *   kick       a recessed kick plate under the door, slotted for the oven burner's combustion air;
 *   face       the face frame round the oven opening, slotted where the door's spring hinge arms run in;
 *   door       a 1 in stainless door pan on a bottom hinge, a tube handle on cast brackets, and inside it a porcelain
 *              liner pressed into a raised panel that fits the opening, vent slots along its top;
 *   oven       a porcelain-enamel cavity (blue-grey, white-flecked) with rounded corners, wire side ladders with four
 *              runner positions, two chrome racks, the thermostat's capillary bulb, the flue slot at the back;
 *   manifold   the valve panel folding over a rounded nose into a landing ledge; seven knobs (six burner valves and the
 *              oven thermostat) on printed dial plates, polished bezels and valve stems;
 *   cooktop    side and back rails the grates stand on, a full-width drip tray, cast-iron venturi tubes rising to
 *              lift-off star burner heads, a standing pilot between each front and back burner with flash tubes to
 *              both, and three cast grates of drafted, crowned bars with tapered fingers;
 *   back       a flue riser with a louvred flue outlet, and a high shelf with a hemmed, turned-down front edge.
 *
 * Joints: "oven_door" (about x at its bottom hinge), "knob_0".."knob_5" (burner valves, front row then back row, left
 * to right), "knob_6" (oven thermostat). Hooks: `flames` (burners 0..5, pilots as burner 6), `open` swings the door.
 * Sockets the manifest names (handle at y 0.704, z -0.01; the pan seat over the front-centre burner) stay where they were.
 */
import * as THREE from "three";
import { mergeVertices } from "three/examples/jsm/utils/BufferGeometryUtils.js";
import { ringPorts, FlameSet, type Port } from "./flames";
import { register, type Builder, type Hooks } from "./index";
import { MAX_ANISOTROPY, palette } from "./materials";
import { castSection, cyl, lathe, plate, rbox, rectPath, roundRect, screw, sweep, turnJoint, wire, type Shop, type V3 } from "./shop";

const TAU = Math.PI * 2;
const smooth = (a: number, b: number, x: number) => {
  const t = Math.max(0, Math.min(1, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Parts                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/** The outline of a six-lobed star burner head of mean radius `r`: points and outward normals, counter-clockwise. */
function starOutline(r: number, lobes: number, depth: number, n: number): { p: [number, number][]; nrm: [number, number][] } {
  const p: [number, number][] = [];
  const nrm: [number, number][] = [];
  for (let i = 0; i < n; i++) {
    const a = (i / n) * Math.PI * 2;
    const rr = r + depth * Math.cos(lobes * a);
    const dr = -depth * lobes * Math.sin(lobes * a);
    const [x, y] = [rr * Math.cos(a), rr * Math.sin(a)];
    // Tangent d/da (r cos a, r sin a); the outward normal is it turned clockwise.
    const [tx, ty] = [dr * Math.cos(a) - rr * Math.sin(a), dr * Math.sin(a) + rr * Math.cos(a)];
    const l = Math.hypot(tx, ty);
    p.push([x, y]);
    nrm.push([ty / l, -tx / l]);
  }
  return { p, nrm };
}

/**
 * A slab cut from `shape` (its own xy), `len` long along +z from 0, its edges eased by `bevel` without growing past the
 * outline (unlike `plate`, whose bevel bulges outward).
 */
function slab(s: Shop, shape: THREE.Shape, len: number, bevel = 0, curve = 12): THREE.BufferGeometry {
  const b = s.fine(1) ? bevel : 0;
  const g = new THREE.ExtrudeGeometry(shape, {
    depth: Math.max(1e-5, len - 2 * b),
    bevelEnabled: b > 0,
    bevelThickness: b,
    bevelSize: b,
    bevelOffset: -b,
    bevelSegments: s.detail === 2 ? 2 : 1,
    curveSegments: s.seg(curve, 3),
  });
  if (b > 0) g.translate(0, 0, b);
  return g;
}

/**
 * A lift-off cast-iron star burner centred at (x, z) with its port line at height y, on its venturi throat. Returns
 * the ports, numbered for `burner`. `r` is the head's mean radius (a 30 000 BTU/h head is about 55 mm). The head is a
 * star-shaped wall 14 mm tall with the ports drilled round it at mid-height, a lid with a raised centre boss over it,
 * and a bell underneath narrowing into the throat the venturi delivers its premixed gas into.
 */
export function starBurner(s: Shop, x: number, y: number, z: number, burner: number, r = 0.052): Port[] {
  const p = palette();
  const lobes = 6;
  const depth = r * 0.2;
  const n = s.detail === 2 ? 120 : s.detail === 1 ? 48 : 24;
  const { p: outline, nrm } = starOutline(r, lobes, depth, n);
  const v2 = (pts: [number, number][]) => pts.map(([a, b]) => new THREE.Vector2(a, b));
  const ease = s.fine(2) ? 0.0012 : 0;
  // The wall (shape y maps to -z after the turn).
  const shape = new THREE.Shape(v2(outline));
  shape.holes.push(new THREE.Path(v2(starOutline(r - 0.0065, lobes, depth * 0.9, n).p.reverse())));
  const wall = slab(s, shape, 0.014, ease, 12);
  wall.rotateX(-Math.PI / 2);
  s.add(p.castIron, wall, [x, y - 0.008, z]);
  // The lid: a star a little inside the wall's outline, its edge eased, a raised boss at its centre.
  const lid = slab(s, new THREE.Shape(v2(starOutline(r - 0.0012, lobes, depth * 0.97, n).p)), 0.0034, s.fine(2) ? 0.0011 : 0);
  lid.rotateX(-Math.PI / 2);
  s.add(p.castIron, lid, [x, y + 0.0055, z]);
  s.add(p.castIron, lathe(s, [
    [0.019, 0.0085],
    [0.017, 0.0105],
    [0.012, 0.0118],
    [0, 0.012],
  ], 28), [x, y, z]);
  // The bell and throat beneath: climbing from the throat up to the wall's foot, so that it faces outward.
  s.add(p.castIron, lathe(s, [
    [0.0165, -0.072],
    [0.0172, -0.05],
    [0.019, -0.034],
    [r * 0.56, -0.02],
    [r * 0.74, -0.012],
    [r * 0.8, -0.0075],
  ], 32), [x, y, z]);
  // The ports: a hole every few millimetres round the outline. Close up each is a dark bore just proud of the wall.
  const nPorts = s.fine(1) ? 36 : 24;
  const ports: Port[] = [];
  for (let i = 0; i < nPorts; i++) {
    const k = Math.floor(((i + 0.5) / nPorts) * outline.length);
    const [px, pz] = outline[k]!;
    const [nx, nz] = nrm[k]!;
    ports.push({ at: [x + px, y - 0.001, z - pz], dir: [nx * 0.85, 0.52, -nz * 0.85], burner });
    if (s.fine(2)) s.add(p.enamel, cyl(s, 0.0008, 0.0005, "x", 8), [x + px + nx * 0.00015, y - 0.001, z - pz - nz * 0.00015], [0, Math.atan2(nz, nx), 0]);
  }
  return ports;
}

/** A rounded rectangle's centre line in plan at height y: `w` by `d` about (cx, cz), corners of radius r. */
function planLoop(cx: number, cz: number, w: number, d: number, r: number, y: number, arc: number): V3[] {
  const out: V3[] = [];
  const c: [number, number][] = [
    [cx + w / 2 - r, cz + d / 2 - r],
    [cx - w / 2 + r, cz + d / 2 - r],
    [cx - w / 2 + r, cz - d / 2 + r],
    [cx + w / 2 - r, cz - d / 2 + r],
  ];
  for (let k = 0; k < 4; k++)
    for (let i = 0; i <= arc; i++) {
      const a = ((k + i / arc) * Math.PI) / 2;
      out.push([c[k]![0] + r * Math.cos(a), y, c[k]![1] + r * Math.sin(a)]);
    }
  return out;
}

/**
 * A cast-iron grate `w` by `d` centred at (x, z), its crowns at `top`, standing `h` tall on cast feet. The frame is one
 * casting of drafted, crowned bars with rounded corners; a cross bar parts each burner's cell from the next; round
 * each burner (centres at depths `bzs`) eight fingers reach in from the frame and cross bars, deep at the root and
 * tapering in plan and in depth to a rounded tip. The tips' crowns, level with the frame's, carry the pan over the
 * flame. Crowned bars meet along lines, never over coplanar areas.
 */
function castGrate(s: Shop, x: number, z: number, w: number, d: number, top: number, bzs: readonly number[], h = 0.034): void {
  const p = palette();
  const bw = 0.0128;
  const fh = 0.024;
  const arc = s.fine(2) ? 6 : s.fine(1) ? 3 : 2;
  const sec = castSection(bw, fh, 0.0032, 1.22, arc);
  // The frame, swept round its centre line.
  s.add(p.castIron, sweep(planLoop(x, z, w - bw, d - bw, 0.014, top, s.fine(2) ? 5 : s.fine(1) ? 2 : 1), sec, { closed: true }));
  const xl = x - w / 2 + bw / 2;
  const xr = x + w / 2 - bw / 2;
  const zFront = z + d / 2 - bw / 2;
  const zBack = z - d / 2 + bw / 2;
  // Cross bars between neighbouring burners.
  const order = [...bzs].sort((a, b) => b - a);
  const cuts: number[] = [];
  for (let i = 1; i < order.length; i++) {
    const zc = (order[i - 1]! + order[i]!) / 2;
    cuts.push(zc);
    s.add(p.castIron, sweep([[xl, top, zc], [xr, top, zc]], sec));
  }
  // Fingers: tapered in plan (1.25 to 0.75 of a bar), in depth (1 to 0.45), rounding off at the tip.
  const st: [number, number, number][] = s.fine(2)
    ? [[0, 1.22, 1], [0.3, 1.08, 0.86], [0.6, 0.92, 0.66], [0.82, 0.8, 0.52], [0.93, 0.72, 0.47], [0.985, 0.48, 0.42], [1, 0.16, 0.36]]
    : s.fine(1)
      ? [[0, 1.2, 1], [0.6, 0.9, 0.62], [0.93, 0.72, 0.46], [1, 0.3, 0.38]]
      : [[0, 1.15, 1], [1, 0.6, 0.45]];
  const rest = 0.056;
  for (const [bi, bz] of order.entries()) {
    const zHi = bi === 0 ? zFront : cuts[bi - 1]!;
    const zLo = bi === order.length - 1 ? zBack : cuts[bi]!;
    for (let k = 0; k < 8; k++) {
      const a = (k / 8) * TAU + Math.PI / 4;
      const [c, sn] = [Math.cos(a), Math.sin(a)];
      // Where the ray from the burner meets the cell's bars (their centre lines).
      const tx = c > 1e-6 ? (xr - x) / c : c < -1e-6 ? (xl - x) / c : Infinity;
      const tz = sn > 1e-6 ? (zHi - bz) / sn : sn < -1e-6 ? (zLo - bz) / sn : Infinity;
      const reach = Math.min(tx, tz);
      if (reach <= rest + 0.01) continue;
      const root: V3 = [x + c * reach, top, bz + sn * reach];
      const path: V3[] = st.map(([t]) => [root[0] + (x + c * rest - root[0]) * t, top, root[2] + (bz + sn * rest - root[2]) * t]);
      s.add(p.castIron, sweep(path, sec, { scale: (i) => [st[i]![1], st[i]![2]] }));
    }
  }
  // Cast feet under the frame's corners and the cross bars' ends, down to the rails.
  const footH = h - fh + 0.003;
  const feet: [number, number][] = [];
  for (const fz of [zFront - 0.004, ...cuts, zBack + 0.004]) for (const fx of [xl + 0.004, xr - 0.004]) feet.push([fx, fz]);
  for (const [fx, fz] of feet) s.add(p.castIron, cyl(s, 0.0072, footH, "y", 14, 0.0098), [fx, top - h + footH / 2, fz]);
}

/**
 * A cast-iron pot grate, `w` by `d`, its top at `top` and its feet `h` below, centred at (x, z): a frame with eight
 * fingers reaching in toward the one burner under it, their tips the rests a pan stands on over the flame.
 */
export function potGrate(s: Shop, x: number, z: number, w: number, d: number, top: number, h = 0.034): void {
  castGrate(s, x, z, w, d, top, [z], h);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Controls                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A restaurant range knob facing +z from a base at z = 0, radius `r`: a black phenolic moulding with a flared skirt,
 * a grip ringed with 18 rounded flutes, a softly domed face and a raised pointer, with a polished pointer insert close up.
 */
function rangeKnob(j: Shop, r: number): void {
  const p = palette();
  const prof: [number, number][] = [
    [0, 0.004],
    [r * 0.93, 0.004],
    [r, 0.0056],
    [r * 0.985, 0.0082],
    [r * 0.86, 0.0102],
    [r * 0.8, 0.0118],
    [r * 0.79, 0.019],
    [r * 0.78, 0.0255],
    [r * 0.74, 0.0284],
    [r * 0.52, 0.0296],
    [0, 0.03],
  ];
  let g = lathe(j, prof, j.fine(2) ? 144 : 48);
  g.rotateX(Math.PI / 2);
  if (j.fine(1)) {
    g.deleteAttribute("uv");
    g.deleteAttribute("normal");
    g = mergeVertices(g, 1e-7);
    const P = g.attributes.position as THREE.BufferAttribute;
    for (let i = 0; i < P.count; i++) {
      const [vx, vy, vz] = [P.getX(i), P.getY(i), P.getZ(i)];
      const rr = Math.hypot(vx, vy);
      const band = smooth(0.0118, 0.0136, vz) * (1 - smooth(0.0246, 0.0268, vz));
      if (band <= 0 || rr < 1e-6) continue;
      const f = Math.pow(0.5 + 0.5 * Math.cos(18 * Math.atan2(vy, vx)), 2.2);
      const k = 1 - (0.0017 * band * f) / rr;
      P.setXY(i, vx * k, vy * k);
    }
    g.computeVertexNormals();
  }
  j.add(p.bakelite, g);
  // The pointer: a raised rib from the centre to the rim, pointing up at rest.
  j.add(p.bakelite, rbox(j, 0.0044, r * 0.68, 0.0032, 0.0012), [0, r * 0.36, 0.0302]);
  if (j.fine(2)) j.add(p.ssBright, rbox(j, 0.0014, r * 0.5, 0.001, 0.0004), [0, r * 0.42, 0.0316]);
}

type Dial = "burner" | "oven";
const dialCache = new Map<Dial, THREE.MeshStandardMaterial>();

/**
 * The printed stainless dial plate behind a knob, as a material for a disc whose UVs span it. Burner valves: OFF at
 * the top, HI a quarter turn anticlockwise, a flame band tapering round to LO at the bottom (see `valveAngle`). The
 * oven thermostat: OFF at the top, then 150 to 500 °F round anticlockwise.
 */
function dialMaterial(kind: Dial): THREE.MeshStandardMaterial {
  const hit = dialCache.get(kind);
  if (hit) return hit;
  const S = 512;
  const c = document.createElement("canvas");
  c.width = c.height = S;
  const x = c.getContext("2d")!;
  const R = S / 2;
  // Brushed ground, a fine rolled edge.
  const g = x.createLinearGradient(0, 0, 0, S);
  g.addColorStop(0, "#d4d7db");
  g.addColorStop(0.5, "#c3c7cc");
  g.addColorStop(1, "#d0d3d7");
  x.fillStyle = g;
  x.fillRect(0, 0, S, S);
  let seed = 911;
  const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  for (let i = 0; i < 900; i++) {
    const v = 175 + Math.floor(rnd() * 70);
    x.fillStyle = `rgba(${v},${v},${v + 4},0.18)`;
    x.fillRect(rnd() * S, rnd() * S, 30 + rnd() * 200, 1);
  }
  x.strokeStyle = "rgba(40,42,46,0.55)";
  x.lineWidth = 4;
  x.beginPath();
  x.arc(R, R, R - 5, 0, TAU);
  x.stroke();
  // Angle φ is the pointer's turn anticlockwise from straight up, as the front sees it.
  const at = (phi: number, rad: number): [number, number] => [R - rad * Math.sin(phi), R - rad * Math.cos(phi)];
  x.fillStyle = "#16181b";
  x.textAlign = "center";
  x.textBaseline = "middle";
  const label = (t: string, phi: number, rad: number, size: number) => {
    x.font = `700 ${size}px "Helvetica Neue", Helvetica, Arial, sans-serif`;
    const [px, py] = at(phi, rad);
    x.fillText(t, px, py);
  };
  if (kind === "burner") {
    label("OFF", 0, R * 0.83, 46);
    // The flame band: wide at HI (φ = 90°), tapering to LO (φ = 180°).
    x.beginPath();
    const steps = 48;
    for (let i = 0; i <= steps; i++) {
      const phi = Math.PI / 2 + (i / steps) * (Math.PI / 2);
      const [px, py] = at(phi, R * 0.8 + (1 - i / steps) * 22);
      i ? x.lineTo(px, py) : x.moveTo(px, py);
    }
    for (let i = steps; i >= 0; i--) {
      const phi = Math.PI / 2 + (i / steps) * (Math.PI / 2);
      const [px, py] = at(phi, R * 0.8 - (1 - i / steps) * 18);
      x.lineTo(px, py);
    }
    x.closePath();
    x.fill();
    label("HI", Math.PI / 2 - 0.32, R * 0.8, 40);
    label("LO", Math.PI + 0.3, R * 0.8, 40);
    // A tick at OFF.
    const [t0x, t0y] = at(0, R * 0.69);
    const [t1x, t1y] = at(0, R * 0.74);
    x.lineWidth = 6;
    x.strokeStyle = "#16181b";
    x.beginPath();
    x.moveTo(t0x, t0y);
    x.lineTo(t1x, t1y);
    x.stroke();
  } else {
    label("OFF", 0, R * 0.84, 40);
    const temps = [150, 200, 250, 300, 350, 400, 450, 500];
    temps.forEach((t, i) => {
      const phi = 0.62 + (i / (temps.length - 1)) * (TAU - 1.24);
      label(String(t), phi, R * 0.83, 34);
      const [a0x, a0y] = at(phi, R * 0.68);
      const [a1x, a1y] = at(phi, R * 0.72);
      x.lineWidth = 5;
      x.strokeStyle = "#16181b";
      x.beginPath();
      x.moveTo(a0x, a0y);
      x.lineTo(a1x, a1y);
      x.stroke();
    });
    label("°F", Math.PI, R * 0.6, 30);
  }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = MAX_ANISOTROPY;
  const m = new THREE.MeshStandardMaterial({ color: 0xffffff, map: t, metalness: 0.55, roughness: 0.34 });
  m.userData.uv = "own";
  dialCache.set(kind, m);
  return m;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Ranges                                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

interface RangeLayout {
  /** Burner centres across, and their depths behind the front edge. */
  xs: number[];
  zs: number[];
}

/** Knob angle for a valve opening: off points up; turned 90° counter-clockwise is high, on round to low at 180°. */
export function valveAngle(flow: number): number {
  return flow <= 0.01 ? 0 : Math.PI / 2 + (1 - flow) * (Math.PI / 2);
}

/** A restaurant range: open burners over a standard oven. The 36 in archetype has six burners. */
function openRange(L: RangeLayout): Builder {
  return (s, e) => {
    const p = palette();
    const { w, d, h } = e;
    const t = 0.0012; // 18 gauge sheet
    const legH = 0.152;
    const zf = -0.056; // the face frame round the oven opening
    const zp = zf + 0.022; // the manifold panel's and the door's front faces
    const deck = 0.88; // the rails the grates stand on
    const top = 0.914; // the grates' crowns: the cooking height
    const back = -d;
    const riserD = 0.095;
    const riserFront = back + riserD;
    const doorY0 = legH + 0.064;
    const doorY1 = 0.735;
    const panelY0 = 0.742;
    const kickTop = doorY0 - 0.004;
    const railW = 0.064;
    const ledgeBack = -0.095;
    // The oven: opening and cavity.
    const cw = 0.67;
    const ch = 0.37;
    const cx = w / 2;
    const cy0 = doorY0 + 0.095; // the hearth, over the oven burner
    const czB = -0.73;
    const cz0 = zf - t;
    const ovenR = 0.024;

    /* --- Legs: 6 in stainless tube on a bolted plate, a threaded stud into an adjustable bullet foot. ---------- */
    for (const lx of [0.055, w - 0.055])
      for (const lz of [zf - 0.06, back + 0.06]) {
        s.add(p.ssDull, rbox(s, 0.064, 0.005, 0.064, 0.0015), [lx, legH - 0.0025, lz]);
        if (s.fine(2)) for (const [bx, bz] of [[-1, -1], [1, -1], [-1, 1], [1, 1]] as const) s.add(p.ssBright, cyl(s, 0.0042, 0.0032, "y", 6), [lx + bx * 0.024, legH - 0.0066, lz + bz * 0.024]);
        s.add(p.ss, lathe(s, [
          [0, 0.046],
          [0.0185, 0.046],
          [0.019, 0.05],
          [0.019, legH - 0.018],
          [0.021, legH - 0.012],
          [0.024, legH - 0.0055],
          [0.024, legH - 0.005],
        ], 28), [lx, 0, lz]);
        // The stud, the adjusting hex and the foot.
        s.add(p.chrome, cyl(s, 0.0075, 0.022, "y", 12), [lx, 0.036, lz]);
        if (s.fine(1)) s.add(p.ssBright, cyl(s, 0.012, 0.007, "y", 6), [lx, 0.0275, lz]);
        s.add(p.ssBright, lathe(s, [
          [0, 0],
          [0.0235, 0],
          [0.025, 0.003],
          [0.0235, 0.008],
          [0.017, 0.017],
          [0.0105, 0.024],
          [0, 0.0245],
        ], 28), [lx, 0, lz]);
      }

    /* --- The body shell: side panels, back, bottom. Panels abut; none shares a face with another. ------------- */
    const sideTop = deck - 0.03;
    for (const sx of [t / 2, w - t / 2]) s.add(p.ss, rbox(s, t, sideTop - legH, zf - t - back, 0.0005), [sx, (legH + sideTop) / 2, (zf - t + back) / 2]);
    s.add(p.ssDull, rbox(s, w - 2 * t, h - t - legH, t, 0.0005), [w / 2, (legH + h - t) / 2, back + t / 2]);
    s.add(p.ssDull, rbox(s, w - 2 * t, t, zf - 0.006 - t - back - t, 0.0005), [w / 2, legH + t / 2, (zf - 0.006 - t + back + t) / 2]);
    // The gas inlet out of the back, low on the right: a 3/4 in black iron nipple through a flanged boss, capped for
    // shipping, kept within 45 mm of the back so the range still stands against a wall.
    if (s.fine(1)) {
      const gx = w - 0.11;
      const gy = 0.3;
      s.add(p.castIron, lathe(s, [
        [0.0134, 0],
        [0.024, 0],
        [0.024, 0.003],
        [0.0175, 0.008],
        [0.0134, 0.008],
      ], 20).rotateX(-Math.PI / 2), [gx, gy, back]);
      s.add(p.castIron, cyl(s, 0.0134, 0.032, "z", 20), [gx, gy, back - 0.016]);
      s.add(p.castIron, lathe(s, [
        [0.0148, 0],
        [0.0168, 0.0015],
        [0.0168, 0.011],
        [0.0158, 0.0125],
        [0, 0.013],
      ], 6).rotateX(-Math.PI / 2), [gx, gy, back - 0.03]);
    }

    /* --- The kick plate: recessed under the door, slotted for the oven burner's air, dark behind. -------------- */
    {
      const kw = w - 2 * t;
      const kh = kickTop - legH;
      const kick = new THREE.Shape();
      rectPath(kick, kw, kh, 0, 0, 0);
      if (s.fine(1)) {
        const n = 13;
        for (let i = 0; i < n; i++) kick.holes.push(rectPath(new THREE.Path(), 0.038, 0.009, 0.0045, -kw / 2 + 0.09 + ((kw - 0.18) * i) / (n - 1), 0.004, true) as THREE.Path);
        s.add(p.enamel, rbox(s, kw - 0.12, kh - 0.02, 0.004), [w / 2, legH + kh / 2 + 0.004, zf - 0.006 - t - 0.025]);
      }
      s.add(p.ss, plate(s, kick, t, 0, 4), [w / 2, legH + kh / 2, zf - 0.006 - t]);
    }

    /* --- The face frame round the oven opening, slotted for the door's hinge arms. ----------------------------- */
    const hingeX = cw / 2 + 0.018;
    {
      const fh = panelY0 - kickTop;
      const frame = new THREE.Shape();
      rectPath(frame, w, fh, 0, 0, 0);
      const oy = cy0 + ch / 2 - (kickTop + fh / 2);
      frame.holes.push(rectPath(new THREE.Path(), cw, ch, ovenR, 0, oy, true) as THREE.Path);
      if (s.fine(1)) for (const sx of [-1, 1]) frame.holes.push(rectPath(new THREE.Path(), 0.006, 0.026, 0.002, sx * hingeX, doorY0 + 0.016 - (kickTop + fh / 2), true) as THREE.Path);
      s.add(p.ss, plate(s, frame, t), [w / 2, kickTop + fh / 2, zf - t]);
      if (s.fine(1)) for (const sx of [-1, 1]) s.add(p.enamel, rbox(s, 0.02, 0.04, 0.004), [cx + sx * hingeX, doorY0 + 0.016, zf - 0.03]);
    }

    /* --- The oven door on its bottom hinge. --------------------------------------------------------------------- */
    const dw = w - 0.05;
    const dh = doorY1 - doorY0;
    const door = s.joint("oven_door", [w / 2, doorY0, zf]);
    // The outer pan, 20 mm deep, standing 1.5 mm off the face frame.
    door.add(p.ss, rbox(door, dw, dh, 0.0205, 0.0035), [0, dh / 2, 0.0015 + 0.01025]);
    // The handle: a 7/8 in tube on two cast brackets, its centre line at the manifest's handle socket.
    const hl = 0.32;
    const hy = 0.704 - doorY0;
    door.add(p.ssBright, wire(door, [[-hl, hy, 0.022], [-hl, hy, 0.046], [hl, hy, 0.046], [hl, hy, 0.022]], 0.0105, 0.017, 16));
    for (const sx of [-1, 1]) {
      door.add(p.ssBright, lathe(door, [
        [0, 0],
        [0.0175, 0],
        [0.0175, 0.0012],
        [0.0145, 0.0042],
        [0.0118, 0.0062],
        [0, 0.0064],
      ], 28).rotateX(Math.PI / 2), [sx * hl, hy, 0.022]);
    }
    // The liner: a flat porcelain flange over the pan's back, and the raised pressing that fits into the opening.
    if (door.fine(1)) {
      door.add(p.ovenEnamel, plate(door, roundRect(dw - 0.014, dh - 0.014, 0.008), t), [0, dh / 2, 0.0002]);
      const bev = 0.007;
      const pw = cw - 0.01 - 2 * bev;
      const ph = ch - 0.01 - 2 * bev;
      const py = cy0 + ch / 2 - doorY0;
      door.add(p.ovenEnamel, plate(door, roundRect(pw, ph, 0.022), 0.02, bev), [0, py, -0.0198]);
      // Vent slots across the top of the pressing, and two rows of louvres low down.
      if (door.fine(2)) {
        for (let i = 0; i < 9; i++) door.add(p.enamel, plate(door, roundRect(0.048, 0.0065, 0.003), 0.0008), [-0.24 + i * 0.06, py + ph / 2 - 0.035, -0.0201 - 0.0008]);
        for (let i = 0; i < 4; i++) door.add(p.ssBright, cyl(door, 0.004, 0.0012, "z", 10), [(i < 2 ? -1 : 1) * (pw / 2 - 0.012), py + (i % 2 ? 1 : -1) * (ph / 2 - 0.012), -0.0202]);
      }
      // The spring hinge arms, running from the door's corners in through the face frame's slots.
      for (const sx of [-1, 1]) door.add(p.ss, rbox(door, 0.003, 0.11, 0.014, 0.0008), [sx * hingeX, -0.005, -0.014]);
    }

    /* --- The oven cavity: porcelain enamel with rounded corners, ladders, racks, thermostat bulb, flue. ------- */
    if (s.fine(1)) {
      const tl = 0.0012;
      const R = ovenR;
      const cd = cz0 - czB;
      const ring = roundRect(cw + 2 * tl, ch + 2 * tl, R + tl);
      ring.holes.push(rectPath(new THREE.Path(), cw, ch, R, 0, 0, true) as THREE.Path);
      s.add(p.ovenEnamel, new THREE.ExtrudeGeometry(ring, { depth: cd, bevelEnabled: false, curveSegments: s.seg(8, 3) }), [cx, cy0 + ch / 2, czB]);
      s.add(p.ovenEnamel, plate(s, roundRect(cw + 2 * tl, ch + 2 * tl, R + tl), tl), [cx, cy0 + ch / 2, czB - tl]);
      // The flue slot at the back top, and the heat slots along the floor's sides.
      s.add(p.enamel, plate(s, roundRect(0.24, 0.018, 0.009), 0.0008), [cx, cy0 + ch - 0.04, czB]);
      if (s.fine(2))
        for (const sx of [-1, 1])
          for (let i = 0; i < 6; i++) s.add(p.enamel, plate(s, roundRect(0.012, 0.07, 0.006), 0.0008).rotateX(-Math.PI / 2), [cx + sx * (cw / 2 - 0.035), cy0, czB + 0.07 + i * 0.1]);
      // Wire side ladders: two posts each, four runners, hooked into the walls.
      const levels = [0.19, 0.4, 0.61, 0.82].map((f) => cy0 + ch * f);
      const zA = cz0 - 0.045;
      const zB = czB + 0.045;
      const lx = cw / 2 - 0.012;
      for (const sx of [-1, 1]) {
        const xw = cx + sx * (cw / 2);
        const xr = cx + sx * lx;
        for (const zz of [zA, zB]) {
          s.add(p.chrome, wire(s, [[xw, cy0 + 0.03, zz], [xr, cy0 + 0.04, zz], [xr, cy0 + ch - 0.04, zz], [xw, cy0 + ch - 0.03, zz]], 0.0026, 0.008, 8));
        }
        for (const y of levels) {
          s.add(p.chrome, wire(s, [[xr, y, zA + 0.004], [xr - sx * 0.004, y, zA - 0.01], [xr - sx * 0.004, y, zB + 0.01], [xr, y, zB - 0.004]], 0.0024, 0.01, 8));
        }
      }
      // Two racks on the middle runners: a framed grid of wires front to back, two cross wires under them, a stop
      // turned up at the back.
      const rw = 2 * lx - 0.008;
      const rd = zA - zB + 0.004;
      for (const y0 of [levels[1]!, levels[2]!]) {
        const y = y0 + 0.0024 + 0.0032;
        const zc = (zA + zB) / 2;
        const loop = planLoop(cx, zc, rw, rd, 0.016, y, s.fine(2) ? 4 : 2);
        s.add(p.chrome, sweep(loop, circleSection(0.0032, s.fine(2) ? 8 : 5), { closed: true, caps: false }));
        const n = s.fine(2) ? 26 : 13;
        for (let i = 1; i < n; i++) s.add(p.chrome, cyl(s, 0.0021, rd - 0.006, "z", s.fine(2) ? 8 : 5), [cx - rw / 2 + (rw * i) / n, y + 0.0004, zc]);
        for (const f of [0.3, 0.7]) s.add(p.chrome, cyl(s, 0.0026, rw - 0.026, "x", 8), [cx, y - 0.0043, zc - rd / 2 + rd * f]);
        s.add(p.chrome, wire(s, [[cx - rw / 2 + 0.03, y, zc - rd / 2], [cx - rw / 2 + 0.03, y + 0.03, zc - rd / 2], [cx + rw / 2 - 0.03, y + 0.03, zc - rd / 2], [cx + rw / 2 - 0.03, y, zc - rd / 2]], 0.0026, 0.008, 8));
      }
      // The thermostat's capillary bulb, clipped high on the left wall.
      if (s.fine(2)) {
        s.add(p.copper, cyl(s, 0.0042, 0.16, "z", 14), [cx - cw / 2 + 0.008, cy0 + ch - 0.05, czB + 0.2]);
        s.add(p.copper, wire(s, [[cx - cw / 2 + 0.008, cy0 + ch - 0.05, czB + 0.12], [cx - cw / 2 + 0.006, cy0 + ch - 0.05, czB + 0.06], [cx - cw / 2 + 0.002, cy0 + ch - 0.03, czB + 0.02], [cx - cw / 2, cy0 + ch - 0.03, czB]], 0.0012, 0.01, 6));
        for (const zz of [czB + 0.15, czB + 0.25]) s.add(p.ssDull, rbox(s, 0.003, 0.014, 0.008, 0.0008), [cx - cw / 2 + 0.0015, cy0 + ch - 0.05, zz]);
      }
    }

    /* --- The manifold panel: one folded sheet from the panel face over a rounded nose into the landing ledge. --- */
    {
      const qp = -zp;
      const qb = -ledgeBack;
      const rn = 0.011;
      const prof = new THREE.Shape();
      prof.moveTo(qp, panelY0);
      prof.lineTo(qp, deck - rn);
      prof.absarc(qp + rn, deck - rn, rn, Math.PI, Math.PI / 2, true);
      prof.lineTo(qb, deck);
      prof.lineTo(qb, sideTop);
      prof.lineTo(-zf + 0.003, sideTop);
      prof.lineTo(-zf + 0.003, panelY0);
      prof.lineTo(qp, panelY0);
      const g = slab(s, prof, w - 2 * t, 0.0008, 10);
      g.rotateY(Math.PI / 2);
      s.add(p.ss, g, [t, 0, 0]);
      // End caps closing the fold, standing a hair proud of the side panels.
      for (const ex of [-0.0005, w - t]) {
        const cap = slab(s, prof, t + 0.0005, 0.0004, 10);
        cap.rotateY(Math.PI / 2);
        s.add(p.ss, cap, [ex, 0, 0]);
      }
      if (s.fine(2)) for (const sx of [0.022, w - 0.022]) for (const sy of [panelY0 + 0.016, deck - 0.03]) screw(s, p.ssBright, [sx, sy, zp], 0.0026);
    }
    // Seven controls: back and front valve of each column side by side, the oven thermostat at the right.
    const knobY = (panelY0 + deck - 0.012) / 2;
    const valves = L.xs.length * L.zs.length;
    const slots = [3, 0, 4, 1, 5, 2, 6];
    const knobJoints: string[] = [];
    for (const [pos, id] of slots.entries()) {
      const kx = 0.07 + ((w - 0.14) * pos) / (slots.length - 1);
      const oven = id === valves;
      const kr = oven ? 0.0265 : 0.0245;
      const j = s.joint(`knob_${id}`, [kx, knobY, zp]);
      rangeKnob(j, kr);
      knobJoints[id] = `knob_${id}`;
      // Fixed on the panel: the dial plate, a polished bezel round the stem, the D-shaft valve stem.
      if (s.fine(2)) {
        const dial = new THREE.CircleGeometry(oven ? 0.044 : 0.04, s.seg(64, 16));
        s.add(dialMaterial(oven ? "oven" : "burner"), dial, [kx, knobY, zp + 0.0004]);
      } else if (s.fine(1)) s.add(p.ssBright, cyl(s, oven ? 0.044 : 0.04, 0.0008, "z", 32), [kx, knobY, zp + 0.0004]);
      if (s.fine(1)) {
        s.add(p.ssBright, lathe(s, [
          [0.0042, 0],
          [0.0095, 0],
          [0.0095, 0.0008],
          [0.0078, 0.0024],
          [0.0042, 0.0026],
        ], 28).rotateX(Math.PI / 2), [kx, knobY, zp + 0.0006]);
        s.add(p.chrome, cyl(s, 0.0031, 0.009, "z", 12), [kx, knobY, zp + 0.0045]);
      }
    }

    /* --- The cooktop: rails, drip tray, venturis, burners, pilots and flash tubes, grates. --------------------- */
    // Side rails the outer grates stand on, from the riser to the ledge; the back rail between them.
    for (const sx of [-0.0006, w - railW]) s.add(p.ss, rbox(s, railW + 0.0006, deck - sideTop, ledgeBack - riserFront - t, 0.0015), [sx + (railW + 0.0006) / 2, (deck + sideTop) / 2, (ledgeBack + riserFront + t) / 2]);
    s.add(p.ss, rbox(s, w - 2 * railW, 0.02, -0.7 - riserFront - t, 0.0015), [w / 2, deck - 0.01, (-0.7 + riserFront + t) / 2]);
    // Grate supports under the seams between grates.
    const pitch = (v: number[]) => (v.length > 1 ? Math.min(...v.slice(1).map((x, i) => x - v[i]!)) : 0.3);
    const gw = pitch(L.xs) - 0.004;
    for (let i = 1; i < L.xs.length; i++) s.add(p.ssHeat, rbox(s, 0.03, 0.018, ledgeBack + 0.7, 0.002), [(L.xs[i - 1]! + L.xs[i]!) / 2, deck - 0.009, (ledgeBack - 0.7) / 2]);
    // The drip tray: a full-width pan with rolled rims.
    {
      const tz0 = ledgeBack - 0.006;
      const tz1 = riserFront + t + 0.006;
      const tw = w - 2 * t - 0.012;
      const ty = 0.752;
      s.add(p.ssHeat, rbox(s, tw, 0.0016, tz0 - tz1, 0.0006), [w / 2, ty, (tz0 + tz1) / 2]);
      if (s.fine(1)) {
        for (const zz of [tz0, tz1]) s.add(p.ssHeat, rbox(s, tw, 0.024, 0.0016, 0.0006), [w / 2, ty + 0.012, zz]);
        for (const xx of [w / 2 - tw / 2, w / 2 + tw / 2]) s.add(p.ssHeat, rbox(s, 0.0016, 0.0232, tz0 - tz1 - 0.0004, 0.0006), [xx, ty + 0.0116, (tz0 + tz1) / 2]);
        s.add(p.ssHeat, cyl(s, 0.003, tw - 0.0012, "x", 10), [w / 2, ty + 0.024, tz0]);
      }
    }
    // Burners on their venturis, the pilots and flash tubes, and the grates.
    const ports: Port[] = [];
    const burnerY = deck + 0.008;
    const pz = -(L.zs[0]! + L.zs[L.zs.length - 1]!) / 2;
    const gd = -0.07 - -0.725;
    const gz = (-0.07 + -0.725) / 2;
    let burner = 0;
    for (const bz of L.zs)
      for (const bx of L.xs) {
        ports.push(...starBurner(s, bx, burnerY, -bz, burner++));
        if (s.fine(1)) {
          // The venturi: a cast tube from its air shutter behind the manifold, back and up into the throat.
          const off = bz < 0.4 ? -0.036 : 0.036;
          const vy = burnerY - 0.082;
          const path: V3[] = [[bx + off, vy, -0.075], [bx + off, vy, -bz + (bz < 0.4 ? 0.07 : 0.08)], [bx, vy, -bz]];
          s.add(p.castIron, wire(s, path, 0.0145, 0.05, s.fine(2) ? 14 : 8));
          s.add(p.castIron, cyl(s, 0.0168, 0.02, "y", 16), [bx, vy + 0.008, -bz]);
          if (s.fine(2)) {
            s.add(p.castIron, cyl(s, 0.019, 0.028, "z", 18, 0.0165), [bx + off, vy, -0.087]);
            s.add(p.brass, cyl(s, 0.0062, 0.016, "z", 6), [bx + off, vy, -0.068]);
          }
        }
      }
    for (const bx of L.xs) {
      castGrate(s, bx, gz, gw, gd, top, L.zs.map((z) => -z));
      // The standing pilot between the front and back burner: its head under a hood, fed by an aluminium tube.
      if (s.fine(1)) {
        s.add(p.ssDull, cyl(s, 0.0042, 0.016, "y", 12), [bx, deck - 0.0055, pz]);
        s.add(p.ssHeat, lathe(s, [
          [0.0045, 0],
          [0.009, 0.002],
          [0.0095, 0.009],
          [0.0085, 0.0095],
          [0.0085, 0.0025],
          [0.0045, 0.001],
        ], 16), [bx, deck - 0.004, pz]);
        s.add(p.ssDull, wire(s, [[bx, deck - 0.013, pz], [bx, deck - 0.05, pz], [bx - 0.065, deck - 0.075, pz + 0.04], [bx - 0.065, deck - 0.085, -0.075]], 0.0016, 0.012, 6));
        // Flash tubes: U-channels from the pilot rising to each burner's ports.
        for (const bz of L.zs) {
          const z1 = -bz + (bz < 0.4 ? -0.062 : 0.062);
          const z0 = pz + (bz < 0.4 ? 0.012 : -0.012);
          const y0 = deck - 0.004;
          const y1 = burnerY - 0.006;
          const len = Math.hypot(z1 - z0, y1 - y0);
          const pitchA = Math.atan2(y1 - y0, Math.abs(z1 - z0)) * (z1 > z0 ? -1 : 1);
          s.at([bx, (y0 + y1) / 2, (z0 + z1) / 2], [pitchA, 0, 0], () => {
            s.add(p.ssHeat, rbox(s, 0.014, 0.0012, len - 0.0006, s.fine(2) ? 0.0004 : 0), [0, -0.005, 0]);
            for (const sx of [-1, 1]) s.add(p.ssHeat, rbox(s, 0.0012, 0.01, len, s.fine(2) ? 0.0004 : 0), [sx * 0.0066, 0, 0]);
          });
        }
      }
      ports.push(...ringPorts(valves, [bx, deck + 0.006, pz], 0.0035, 3, 1.25));
    }

    /* --- Flue riser with its louvred outlet, and the high shelf. ----------------------------------------------- */
    const shelfD = 0.27;
    const shelfLip = 0.038;
    for (const sx of [t / 2, w - t / 2]) s.add(p.ss, rbox(s, t, h - shelfLip - sideTop, riserD + t, 0.0005), [sx, (sideTop + h - shelfLip) / 2, back + (riserD + t) / 2]);
    {
      const rw = w - 2 * t;
      const y0 = 0.752;
      const rh = h - t - y0;
      const riser = new THREE.Shape();
      rectPath(riser, rw, rh, 0, 0, 0);
      const fy = deck + 0.022 - (y0 + rh / 2);
      if (s.fine(1))
        for (let row = 0; row < 3; row++)
          for (let i = 0; i < 6; i++) riser.holes.push(rectPath(new THREE.Path(), 0.118, 0.011, 0.0055, -rw / 2 + 0.085 + i * ((rw - 0.17) / 5), fy + row * 0.022, true) as THREE.Path);
      s.add(p.ss, plate(s, riser, t, 0, 4), [w / 2, y0 + rh / 2, riserFront]);
      s.add(p.enamel, rbox(s, rw - 0.04, 0.09, 0.004), [w / 2, deck + 0.044, riserFront - 0.03]);
    }
    // The shelf: a sheet with its front edge turned down and hemmed back, its ends closed, on two gussets.
    {
      const prof = new THREE.Shape();
      prof.moveTo(0, h);
      prof.lineTo(shelfD, h);
      prof.lineTo(shelfD, h - shelfLip);
      prof.lineTo(shelfD - 0.012, h - shelfLip);
      prof.lineTo(shelfD - 0.012, h - shelfLip + t);
      prof.lineTo(shelfD - t, h - shelfLip + t);
      prof.lineTo(shelfD - t, h - t);
      prof.lineTo(0, h - t);
      prof.lineTo(0, h);
      // Shape x runs forward from the back; turned so, the sheet runs from x = w - t back to t.
      const g = slab(s, prof, w - 2 * t, 0, 4);
      g.rotateY(-Math.PI / 2);
      s.add(p.ss, g, [w - t, 0, back]);
      const end = new THREE.Shape();
      rectPath(end, shelfD, shelfLip, 0.002, shelfD / 2, h - shelfLip / 2);
      for (const ex of [0, w - t]) {
        const eg = slab(s, end, t, 0.0004, 4);
        eg.rotateY(-Math.PI / 2);
        s.add(p.ss, eg, [ex + t, 0, back]);
      }
      for (const sx of [0.03, w - 0.03]) {
        const tri = new THREE.Shape([new THREE.Vector2(0, 0), new THREE.Vector2(shelfD - riserD - 0.012, 0), new THREE.Vector2(0, -0.12)]);
        const gg = plate(s, tri, t);
        gg.rotateY(-Math.PI / 2);
        s.add(p.ss, gg, [sx + t / 2, h - t, riserFront + t]);
      }
    }

    return (b) => {
      const flames = new FlameSet(ports, 0.024);
      b.group.add(flames.mesh);
      flames.setAll(0);
      flames.setFlow(valves, 0.12);
      const hooks: Hooks = {
        flames,
        burners: valves,
        open: (f) => turnJoint(b.joints.get("oven_door"), "x", f * (Math.PI / 2) * 0.98),
        tick: (tt) => flames.tick(tt),
      };
      // Knobs follow their valves.
      const set = flames.setFlow.bind(flames);
      flames.setFlow = (i: number, f: number) => {
        set(i, f);
        if (i < valves) turnJoint(b.joints.get(knobJoints[i]!), "z", valveAngle(f));
      };
      return hooks;
    };
  };
}

/** A round wire's section, radius r, for `sweep`. */
function circleSection(r: number, n: number): [number, number][] {
  return Array.from({ length: n }, (_, i) => [r * Math.cos((i / n) * TAU), r * Math.sin((i / n) * TAU)] as [number, number]);
}

register({
  "NSF-RANGE-GAS-6B": openRange({ xs: [0.18, 0.457, 0.734], zs: [0.22, 0.58] }),
});

export type { V3 };
