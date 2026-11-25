/**
 * Ovens and countertop appliances. Dimensions are the archetypes' (assets/kitchens/specs); the construction is that of
 * the equipment class.
 *
 * A combi-steamer is a convection oven that can also flood its cavity with steam: a radial fan wheel at the back draws
 * cavity air in through the eye of a perforated air baffle and throws it out past a ring element round its rim, while a
 * nozzle sprays water onto the hub to flash it into steam. Pans ride on GN runners of hanging rack frames (or of a
 * roll-in trolley on the 20-pan floor model). The door is triple glazed, hinged left, sealed by a gasket on the cavity
 * frame, and drains its condensate into a drip channel; the controls are a touch screen and a central dial on the
 * right-hand column, which also houses the electronics.
 *
 * Live names (`set`):
 *   combi ovens      fan 0..1 wheel speed · heat 0..1 element glow · steam 0..1 cavity haze · cutaway 0|1 hide the
 *                    right-hand column and cavity wall · dial 0..1 the central dial · trolley 0..1 (20-pan) roll the
 *                    rack out of the cavity
 *   holding cabinet  heat 0..1 temperature knob turns, heat lamp lights · humidity 0..1 humidity knob turns, its
 *                    lamp lights, vapour rises off the water pan
 *   microwave        power 0..1 magnetron on: the cavity lamp lights
 *   toaster          heat 0..1 quartz glow, pilot lamp and heat knob · speed 0..1 belt speed and its knob
 *   air fryer        heat 0..1 element glow · fan 0..1 impeller speed
 *   sous vide        heat 0..1 heater lamp · flow 0..1 the pump's churn of the water surface
 *   steam table      steam 0..1 wisps off the pans · heat 0..1 well knobs turn and pilot lamps light
 *   booster heater   temp 0..1 thermometer needle (full scale) · pressure 0..1 pressure needle · heat 0..1 element lamp
 *
 * Joints: combis cavity_door, fan (fan_2), dial, cutaway, rack_trolley (20-pan); holding cabinet cabinet_door,
 * temp_knob, humidity_knob; microwave door; toaster heat_knob, speed_knob, sprocket_front, sprocket_rear; air fryer
 * basket_slide, fan; sous vide set_wheel; steam table knob_0..3; booster heater temp_needle, pressure_needle.
 * `open(frac)` swings or slides the main door or drawer; `tick` runs fans, belts, steam and water.
 */
import * as THREE from "three";
import { register, type Builder, type Ext, type Hooks } from "./index";
import { glowMaterial, palette, perforation, setGlow } from "./materials";
import {
  ball,
  box,
  caster,
  circleHole,
  cyl,
  lathe,
  plate,
  rbox,
  rectPath,
  roundLeg,
  roundRect,
  screw,
  boxUVs,
  slideJoint,
  torus,
  turnJoint,
  type Shop,
  type V3,
} from "./shop";

type Mat = THREE.Material;

/* ------------------------------------------------------------------------------------------------------------ */
/* The family's own finishes: shared by every model that uses them, made the first time they are asked for.        */
/* ------------------------------------------------------------------------------------------------------------ */

function canvas2d(w: number, h: number): [HTMLCanvasElement, CanvasRenderingContext2D] {
  const c = document.createElement("canvas");
  c.width = w;
  c.height = h;
  return [c, c.getContext("2d")!];
}
function ctex(c: HTMLCanvasElement, srgb = true): THREE.CanvasTexture {
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.anisotropy = 8;
  t.needsUpdate = true;
  return t;
}
function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}
function rrect(x: CanvasRenderingContext2D, x0: number, y0: number, w: number, h: number, r: number): void {
  x.beginPath();
  x.moveTo(x0 + r, y0);
  x.arcTo(x0 + w, y0, x0 + w, y0 + h, r);
  x.arcTo(x0 + w, y0 + h, x0, y0 + h, r);
  x.arcTo(x0, y0 + h, x0, y0, r);
  x.arcTo(x0, y0, x0 + w, y0, r);
  x.closePath();
}

/** A combi's touch screen mid-programme: the cavity set point on a ring gauge, the three modes, the time line. */
function combiUi(): HTMLCanvasElement {
  // Drawn at twice the layout's size, so the type stays sharp from a hand's breadth away.
  const [c, x] = canvas2d(512, 640);
  x.scale(2, 2);
  const g = x.createLinearGradient(0, 0, 0, 320);
  g.addColorStop(0, "#0d1822");
  g.addColorStop(1, "#05090d");
  x.fillStyle = g;
  x.fillRect(0, 0, 256, 320);
  x.fillStyle = "#1b2a36";
  x.fillRect(12, 14, 232, 2);
  for (let i = 0; i < 4; i++) {
    x.fillStyle = i === 1 ? "#f2a23a" : "#3a4c5a";
    rrect(x, 14 + i * 60, 24, 48, 8, 4);
    x.fill();
  }
  // The ring gauge: set point amber, the core probe's progress blue inside it.
  x.lineCap = "round";
  x.lineWidth = 12;
  x.strokeStyle = "#23323e";
  x.beginPath();
  x.arc(128, 132, 72, Math.PI * 0.75, Math.PI * 2.25);
  x.stroke();
  x.strokeStyle = "#f2a23a";
  x.beginPath();
  x.arc(128, 132, 72, Math.PI * 0.75, Math.PI * 1.82);
  x.stroke();
  x.lineWidth = 6;
  x.strokeStyle = "#3aa6f2";
  x.beginPath();
  x.arc(128, 132, 56, Math.PI * 0.75, Math.PI * 1.35);
  x.stroke();
  x.fillStyle = "#f4f7fa";
  x.font = "bold 44px sans-serif";
  x.textAlign = "center";
  x.fillText("185°", 128, 142);
  x.font = "16px sans-serif";
  x.fillStyle = "#8fa3b3";
  x.fillText("72°", 128, 172);
  // Mode keys: steam, convection, combination.
  const modes: [string, number][] = [
    ["#3aa6f2", 52],
    ["#f2633a", 128],
    ["#b05cf2", 204],
  ];
  for (const [col, cx] of modes) {
    x.fillStyle = "#142230";
    x.beginPath();
    x.arc(cx, 238, 26, 0, Math.PI * 2);
    x.fill();
    x.strokeStyle = col;
    x.lineWidth = cx === 204 ? 4 : 2;
    x.stroke();
    x.fillStyle = col;
    if (cx === 52) for (const [dx, dy, r] of [[-7, 3, 7], [3, -3, 9], [11, 4, 6]] as const) {
      x.beginPath();
      x.arc(cx + dx, 238 + dy, r, 0, Math.PI * 2);
      x.fill();
    }
    else
      for (let k = 0; k < 4; k++) {
        x.save();
        x.translate(cx, 238);
        x.rotate((k * Math.PI) / 2 + (cx === 204 ? 0.4 : 0));
        x.beginPath();
        x.ellipse(8, 0, 9, 4, 0.3, 0, Math.PI * 2);
        x.fill();
        x.restore();
      }
  }
  // The time line.
  x.fillStyle = "#23323e";
  rrect(x, 20, 288, 216, 8, 4);
  x.fill();
  x.fillStyle = "#4ade80";
  rrect(x, 20, 288, 132, 8, 4);
  x.fill();
  x.fillStyle = "#c9d6e0";
  x.font = "14px sans-serif";
  x.textAlign = "left";
  x.fillText("12:40", 20, 282);
  x.textAlign = "right";
  x.fillText("20:00", 236, 282);
  return c;
}

/** A membrane keypad: the number pad and the start and stop keys, printed on a dark polyester overlay. */
function keypadArt(): HTMLCanvasElement {
  const [c, x] = canvas2d(192, 384);
  x.fillStyle = "#1c1f23";
  x.fillRect(0, 0, 192, 384);
  const keys = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "×2", "0", "+30"];
  x.textAlign = "center";
  x.textBaseline = "middle";
  for (let i = 0; i < 12; i++) {
    const cx = 18 + (i % 3) * 58;
    const cy = 16 + Math.floor(i / 3) * 62;
    x.fillStyle = "#2f343a";
    rrect(x, cx, cy, 48, 50, 7);
    x.fill();
    x.strokeStyle = "#6b737c";
    x.lineWidth = 1.5;
    x.stroke();
    x.fillStyle = "#eef1f4";
    x.font = i === 9 || i === 11 ? "bold 17px sans-serif" : "bold 24px sans-serif";
    x.fillText(keys[i]!, cx + 24, cy + 26);
  }
  // Stop and start.
  for (const [col, cx, sym] of [
    ["#b8322c", 18, "stop"],
    ["#2f9a52", 104, "start"],
  ] as const) {
    x.fillStyle = col;
    rrect(x, cx, 268, 70, 96, 9);
    x.fill();
    x.fillStyle = "#ffffff";
    if (sym === "stop") x.fillRect(cx + 23, 304, 24, 24);
    else {
      x.beginPath();
      x.moveTo(cx + 25, 300);
      x.lineTo(cx + 50, 316);
      x.lineTo(cx + 25, 332);
      x.closePath();
      x.fill();
    }
  }
  return c;
}

/** Seven-segment style readout on a dark window. */
function readout(text: string, col: string, w = 256, h = 96, sub?: string): HTMLCanvasElement {
  const [c, x] = canvas2d(w, h);
  x.fillStyle = "#06090b";
  x.fillRect(0, 0, w, h);
  x.fillStyle = col;
  x.shadowColor = col;
  x.shadowBlur = 8;
  x.font = `bold ${Math.round(h * 0.62)}px monospace`;
  x.textAlign = "center";
  x.textBaseline = "middle";
  x.fillText(text, w / 2, h * (sub ? 0.42 : 0.52));
  if (sub) {
    x.shadowBlur = 3;
    x.font = `${Math.round(h * 0.2)}px sans-serif`;
    x.fillText(sub, w / 2, h * 0.84);
  }
  return c;
}

/** A gauge dial face: white enamel, scale from 0 to `max` in `major` steps, a red band above `red`. */
function dialFace(max: number, major: number, red: number, unit: string): HTMLCanvasElement {
  const S = 256;
  const [c, x] = canvas2d(S, S);
  x.fillStyle = "#f4f2ec";
  x.fillRect(0, 0, S, S);
  const cx = S / 2;
  const a0 = Math.PI * 0.75;
  const span = Math.PI * 1.5;
  const ang = (v: number) => a0 + (v / max) * span;
  x.strokeStyle = "#c8302a";
  x.lineWidth = 10;
  x.beginPath();
  x.arc(cx, cx, 100, ang(red), ang(max));
  x.stroke();
  x.strokeStyle = "#16181b";
  const minor = major / 5;
  for (let v = 0; v <= max + 1e-6; v += minor) {
    const a = ang(v);
    const big = Math.abs(v / major - Math.round(v / major)) < 1e-6;
    x.lineWidth = big ? 3 : 1.4;
    x.beginPath();
    x.moveTo(cx + Math.cos(a) * (big ? 86 : 94), cx + Math.sin(a) * (big ? 86 : 94));
    x.lineTo(cx + Math.cos(a) * 106, cx + Math.sin(a) * 106);
    x.stroke();
    if (big) {
      x.fillStyle = "#16181b";
      x.font = "bold 19px sans-serif";
      x.textAlign = "center";
      x.textBaseline = "middle";
      x.fillText(String(Math.round(v)), cx + Math.cos(a) * 68, cx + Math.sin(a) * 68);
    }
  }
  x.font = "bold 18px sans-serif";
  x.fillText(unit, cx, cx + 52);
  x.strokeStyle = "#16181b";
  x.lineWidth = 4;
  x.beginPath();
  x.arc(cx, cx, 124, 0, Math.PI * 2);
  x.stroke();
  return c;
}

/** Smooth value noise on a seeded 64 × 64 lattice (wrapping), about 0..1. */
function noiseField(seed: number): (u: number, v: number) => number {
  const r = rng(seed);
  const lat = new Float32Array(64 * 64);
  for (let i = 0; i < lat.length; i++) lat[i] = r();
  return (u, v) => {
    const i = Math.floor(u);
    const j = Math.floor(v);
    const fu = u - i;
    const fv = v - j;
    const su = fu * fu * (3 - 2 * fu);
    const sv = fv * fv * (3 - 2 * fv);
    const at = (a: number, b: number) => lat[(a & 63) + (b & 63) * 64]!;
    const [a, b, c, d] = [at(i, j), at(i + 1, j), at(i, j + 1), at(i + 1, j + 1)];
    return a + (b - a) * su + (c - a) * sv + (a - b - c + d) * su * sv;
  };
}
/** Fractal sum of `oct` octaves of a noise field, normalised to about 0..1. */
function fbm(n: (u: number, v: number) => number, u: number, v: number, oct = 5): number {
  let a = 0;
  let amp = 0.5;
  let f = 1;
  let norm = 0;
  for (let k = 0; k < oct; k++) {
    a += amp * n(u * f + k * 17.3, v * f + k * 9.1);
    norm += amp;
    f *= 2.03;
    amp *= 0.5;
  }
  return a / norm;
}
function sstep(a: number, b: number, t: number): number {
  const k = Math.max(0, Math.min(1, (t - a) / (b - a)));
  return k * k * (3 - 2 * k);
}

/**
 * A wisp of steam, twice as tall as it is wide: a plume that leaves its source narrow and spreads as it rises, torn
 * into strands by domain-warped noise, fading out at its foot and its crown. White; the alpha carries it all.
 */
function wisp(seed: number): HTMLCanvasElement {
  const W = 96;
  const H = 192;
  const [c, x] = canvas2d(W, H);
  const img = x.createImageData(W, H);
  const n = noiseField(seed);
  for (let py = 0; py < H; py++) {
    const v = py / (H - 1);
    const up = 1 - v;
    const cx = 0.5 + 0.34 * (fbm(n, 1.3, v * 1.8, 3) - 0.5) * up;
    const half = 0.09 + 0.24 * up;
    for (let px = 0; px < W; px++) {
      const u = px / (W - 1);
      const dx = (u - cx) / half;
      const env = Math.exp(-dx * dx * 1.8) * sstep(0, 0.3, v) * sstep(0, 0.22, up);
      const warp = fbm(n, u * 2.6 + 4.1, v * 1.4 + 2.7, 3);
      const f = fbm(n, u * 3.8 + warp * 2.4, v * 2.1 - warp * 1.3, 5);
      const strands = sstep(0.42, 0.8, f);
      const fine = 0.55 + 0.45 * sstep(0.35, 0.85, fbm(n, u * 9 + 31, v * 5 + 7, 3));
      const a = Math.min(1, env * strands * fine * 1.25);
      const o = (py * W + px) * 4;
      img.data[o] = 255;
      img.data[o + 1] = 255;
      img.data[o + 2] = 255;
      img.data[o + 3] = Math.round(255 * a);
    }
  }
  x.putImageData(img, 0, 0);
  return c;
}
let wispTexs: THREE.CanvasTexture[] | null = null;
/** Four different wisps, made once. */
function wisps(): THREE.CanvasTexture[] {
  if (!wispTexs) wispTexs = [301, 302, 303, 304].map((sd) => ctex(wisp(sd)));
  return wispTexs;
}

/**
 * The black ceramic frit fired onto the back of a door's outer pane, for a pane W × H centred on (cx, cy): solid out
 * to the pane's edges, with a clear window `win` (absolute, in the pane's plane) whose corners are rounded by `r`.
 * Built as geometry, so the window's edge stays exact however close one looks.
 */
function fritGeom(s: Shop, W: number, H: number, cx: number, cy: number, win: { x0: number; x1: number; y0: number; y1: number }, r: number): THREE.BufferGeometry {
  const sh = roundRect(W, H, 0.003, cx, cy);
  rectPath(sh, win.x1 - win.x0, win.y1 - win.y0, r, (win.x0 + win.x1) / 2, (win.y0 + win.y1) / 2, true);
  return new THREE.ShapeGeometry(sh, s.seg(10, 3));
}

/**
 * A sheet of glass w × h × t centred on `at`, facing ±z: its two faces in `face` and its four edges in `edge` (glass
 * seen edge on is green and nearly opaque, which is how a pane's thickness reads).
 */
function pane(s: Shop, face: Mat, edge: Mat, w: number, h: number, t: number, at: V3): void {
  s.add(face, new THREE.PlaneGeometry(w, h), [at[0], at[1], at[2] + t / 2]);
  s.add(face, new THREE.PlaneGeometry(w, h).rotateY(Math.PI), [at[0], at[1], at[2] - t / 2]);
  if (!s.fine(1)) return;
  for (const sg of [-1, 1]) {
    s.add(edge, new THREE.PlaneGeometry(t, h).rotateY((sg * Math.PI) / 2), [at[0] + (sg * w) / 2, at[1], at[2]]);
    s.add(edge, new THREE.PlaneGeometry(w, t).rotateX((-sg * Math.PI) / 2), [at[0], at[1] + (sg * h) / 2, at[2]]);
  }
}

/**
 * Draw a door's glass from the inside out (inner panes, their edges, then the tinted outer pane), so that seen from the
 * front each layer is blended over the one behind it. Every transparent part a shop makes shares one render order.
 */
function layerGlass(root: THREE.Object3D): void {
  const fm = mats();
  root.traverse((o) => {
    const m = o as THREE.Mesh;
    if (!m.isMesh) return;
    if (m.material === fm.paneClear) m.renderOrder = 3;
    else if (m.material === fm.paneEdge) m.renderOrder = 4;
    else if (m.material === fm.paneSmoke) m.renderOrder = 5;
  });
}

/** A section point: x, y and the unit normal (nx, ny) of the face that is to be seen. */
type SecPt = readonly [number, number, number, number];
/** An arc of a concave cove about (cx, cy), radius r, from angle a0 to a1, its normals toward the centre. */
function covePts(cx: number, cy: number, r: number, a0: number, a1: number, n: number): SecPt[] {
  const out: SecPt[] = [];
  for (let i = 0; i <= n; i++) {
    const a = a0 + ((a1 - a0) * i) / n;
    out.push([cx + r * Math.cos(a), cy + r * Math.sin(a), -Math.cos(a), -Math.sin(a)]);
  }
  return out;
}
/**
 * A sheet `t` thick swept along z from z0 to z1 through the open section `pts`: its face along the given normals,
 * its back `t` behind, and its two long edges closed (so a cutaway shows the sheet's thickness). Smooth across the
 * section, as a deep-drawn liner's coves are.
 */
function sweepZ(pts: readonly SecPt[], z0: number, z1: number, t: number): THREE.BufferGeometry {
  const pos: number[] = [];
  const nrm: number[] = [];
  const idx: number[] = [];
  const quad = (a: number, b: number, c: number, d: number, nx: number, ny: number, nz: number) => {
    // a, b, c, d round the quad; keep the winding whose face normal agrees with (nx, ny, nz).
    const P = (i: number) => new THREE.Vector3(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);
    const fn = new THREE.Vector3().subVectors(P(b), P(a)).cross(new THREE.Vector3().subVectors(P(c), P(a)));
    if (fn.x * nx + fn.y * ny + fn.z * nz >= 0) idx.push(a, b, c, a, c, d);
    else idx.push(a, c, b, a, d, c);
  };
  for (const [off, flip] of [
    [0, 1],
    [t, -1],
  ] as const) {
    const o = pos.length / 3;
    for (const [x, y, nx, ny] of pts)
      for (const z of [z0, z1]) {
        pos.push(x - nx * off, y - ny * off, z);
        nrm.push(nx * flip, ny * flip, 0);
      }
    for (let i = 0; i < pts.length - 1; i++) {
      const a = o + i * 2;
      const [nx, ny] = [pts[i]![2] + pts[i + 1]![2], pts[i]![3] + pts[i + 1]![3]];
      quad(a, a + 2, a + 3, a + 1, nx * flip, ny * flip, 0);
    }
  }
  for (const [i, j] of [
    [0, 1],
    [pts.length - 1, pts.length - 2],
  ] as const) {
    const [x, y, nx, ny] = pts[i]!;
    const tx = x - pts[j]![0];
    const ty = y - pts[j]![1];
    const tl = Math.hypot(tx, ty) || 1;
    const o = pos.length / 3;
    for (const [px, py] of [
      [x, y],
      [x - nx * t, y - ny * t],
    ] as const)
      for (const z of [z0, z1]) {
        pos.push(px, py, z);
        nrm.push(tx / tl, ty / tl, 0);
      }
    quad(o, o + 2, o + 3, o + 1, tx, ty, 0);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("normal", new THREE.Float32BufferAttribute(nrm, 3));
  g.setIndex(idx);
  return g;
}

/**
 * UVs projected along z onto `bounds` (the geometry's own xy bounds if none), 0..1: a slice or a steak painted as one
 * picture, its sides taking the colour at the picture's edge.
 */
function planarUV(g: THREE.BufferGeometry, bounds?: THREE.Box2): THREE.BufferGeometry {
  g.computeBoundingBox();
  const b = bounds ?? new THREE.Box2(new THREE.Vector2(g.boundingBox!.min.x, g.boundingBox!.min.y), new THREE.Vector2(g.boundingBox!.max.x, g.boundingBox!.max.y));
  const p = g.attributes.position as THREE.BufferAttribute;
  const uv = new Float32Array(p.count * 2);
  const [sx, sy] = [b.max.x - b.min.x || 1, b.max.y - b.min.y || 1];
  for (let i = 0; i < p.count; i++) {
    uv[i * 2] = (p.getX(i) - b.min.x) / sx;
    uv[i * 2 + 1] = (p.getY(i) - b.min.y) / sy;
  }
  g.setAttribute("uv", new THREE.BufferAttribute(uv, 2));
  return g;
}

/** A path's outline, sampled, in the canvas frame that maps the bounds [bx0, bx1] × [by0, by1] onto S × S. */
function tracePath(x: CanvasRenderingContext2D, pts: THREE.Vector2[], b: THREE.Box2, S: number): void {
  const sx = S / (b.max.x - b.min.x);
  const sy = S / (b.max.y - b.min.y);
  x.beginPath();
  for (const [i, q] of pts.entries()) {
    const [cx, cy] = [(q.x - b.min.x) * sx, S - (q.y - b.min.y) * sy];
    if (i) x.lineTo(cx, cy);
    else x.moveTo(cx, cy);
  }
  x.closePath();
}

interface FamilyMats {
  /** Perforated stainless sheet: 5 mm holes at 8 mm stagger (air baffles, ducts, baskets). */
  perf: THREE.MeshStandardMaterial;
  /** The microwave door's screen: black steel with 1.5 mm holes, the choke that keeps the waves in. */
  screenMesh: THREE.MeshStandardMaterial;
  combiUi: THREE.MeshBasicMaterial;
  keypad: THREE.MeshStandardMaterial;
  mwDisplay: THREE.MeshBasicMaterial;
  fryDisplay: THREE.MeshBasicMaterial;
  holdDisplay: THREE.MeshBasicMaterial;
  sousDisplay: THREE.MeshBasicMaterial;
  dialTemp: THREE.MeshStandardMaterial;
  dialPress: THREE.MeshStandardMaterial;
  /** Clear polycarbonate (tubs, lids). */
  poly: THREE.MeshPhysicalMaterial;
  /** Vacuum-bag film. */
  bag: THREE.MeshPhysicalMaterial;
  /** A sous-vide bath's water: clean tap water, nearly clear, faintly blue-green in depth. */
  bath: THREE.MeshPhysicalMaterial;
  /** Red gauge needles, valve handles. */
  red: THREE.MeshStandardMaterial;
  /** Grey electrical enclosure paint. */
  greyBox: THREE.MeshStandardMaterial;
  /** Mica sheet: the microwave's waveguide and stirrer covers. */
  mica: THREE.MeshStandardMaterial;
  /** Mineral-wool insulation, seen only in a cutaway. */
  wool: THREE.MeshStandardMaterial;
  /** A halogen lamp seen through frosted glass. */
  lamp: THREE.MeshBasicMaterial;
  /** A basket's heap of fries seen from above (the fill under the loose ones). */
  fries: THREE.MeshStandardMaterial;
  /** Single fries, pale to well browned; each is painted darker at its tips. */
  fry: THREE.MeshStandardMaterial[];
  dough: THREE.MeshStandardMaterial;
  /** A ribeye's cut face, marbled, its fat cap and the seam round the spinalis; sides take the edge's colour. */
  beef: THREE.MeshStandardMaterial;
  /** A salmon portion: coral flesh with its white myosepta in chevrons. */
  salmon: THREE.MeshStandardMaterial;
  gravy: THREE.MeshStandardMaterial;
  /** Door glass: the outer pane's grey-bronze tint, clear inner panes, and edges seen edge on. */
  paneSmoke: THREE.MeshPhysicalMaterial;
  paneClear: THREE.MeshPhysicalMaterial;
  paneEdge: THREE.MeshStandardMaterial;
  /** The black glass fascia of a combi's control column. */
  blackGlass: THREE.MeshPhysicalMaterial;
  /** Black ceramic frit, seen through the glass it is fired onto. */
  frit: THREE.MeshStandardMaterial;
  /** A microwave cavity's punched lamp grille and air inlet: steel with 2.5 mm holes at 4 mm. */
  grille: THREE.MeshStandardMaterial;
}

/** The outline of a slice of sandwich bread (the pan loaf's crown flaring over its sides), centred, in metres. */
function breadShape(): THREE.Shape {
  const s = new THREE.Shape();
  s.moveTo(-0.036, -0.052);
  s.lineTo(0.036, -0.052);
  s.quadraticCurveTo(0.047, -0.052, 0.047, -0.041);
  s.lineTo(0.049, 0.014);
  s.bezierCurveTo(0.05, 0.021, 0.056, 0.023, 0.0555, 0.032);
  s.bezierCurveTo(0.054, 0.051, 0.031, 0.06, 0, 0.06);
  s.bezierCurveTo(-0.031, 0.06, -0.054, 0.051, -0.0555, 0.032);
  s.bezierCurveTo(-0.056, 0.023, -0.05, 0.021, -0.049, 0.014);
  s.lineTo(-0.047, -0.041);
  s.quadraticCurveTo(-0.047, -0.052, -0.036, -0.052);
  return s;
}

/**
 * A slice's face as one picture over the outline's bounds: the open crumb (its pores drawn into a bump map as well),
 * the crust round the edge, and when toasted the browning, deepest toward the edges and blotched where the quartz
 * found it. The area outside the outline is crust, so that a slice's sides, projected from it, are crust too.
 */
function breadArt(toasted: boolean): { map: THREE.CanvasTexture; bump: THREE.CanvasTexture } {
  const S = 512;
  const shape = breadShape();
  const pts = shape.getPoints(30);
  const b = new THREE.Box2().setFromPoints(pts);
  const [c, x] = canvas2d(S, S);
  const [cb, xb] = canvas2d(S, S);
  const crust = toasted ? "#6a3816" : "#b87a3e";
  x.fillStyle = crust;
  x.fillRect(0, 0, S, S);
  xb.fillStyle = "#6e6e6e";
  xb.fillRect(0, 0, S, S);
  for (const ctx of [x, xb]) {
    ctx.save();
    tracePath(ctx, pts, b, S);
    ctx.clip();
  }
  x.fillStyle = toasted ? "#d9a35c" : "#efe3c9";
  x.fillRect(0, 0, S, S);
  xb.fillStyle = "#9a9a9a";
  xb.fillRect(0, 0, S, S);
  const r = rng(toasted ? 61 : 62);
  // The crumb's cells: small, elongated in the direction the loaf rose, a few big holes.
  for (let i = 0; i < 3200; i++) {
    const px = r() * S;
    const py = r() * S;
    const big = r() < 0.04;
    const rx = big ? 3 + r() * 5 : 0.8 + r() * 2.2;
    const ry = rx * (1.1 + r() * 0.9);
    const rot = (r() - 0.5) * 0.8;
    const shade = r();
    x.fillStyle = toasted ? `rgba(${110 + shade * 40 | 0},${58 + shade * 30 | 0},${22 + shade * 14 | 0},${0.35 + r() * 0.35})` : `rgba(${176 + shade * 30 | 0},${150 + shade * 26 | 0},${112 + shade * 20 | 0},${0.3 + r() * 0.35})`;
    x.beginPath();
    x.ellipse(px, py, rx, ry, rot, 0, Math.PI * 2);
    x.fill();
    xb.fillStyle = `rgba(30,30,30,${0.5 + r() * 0.4})`;
    xb.beginPath();
    xb.ellipse(px, py, rx, ry, rot, 0, Math.PI * 2);
    xb.fill();
    // The lit lip of the cell.
    x.fillStyle = toasted ? "rgba(238,190,120,0.35)" : "rgba(255,250,236,0.45)";
    x.beginPath();
    x.ellipse(px - rx * 0.35, py - ry * 0.5, rx * 0.7, ry * 0.35, rot, 0, Math.PI * 2);
    x.fill();
  }
  if (toasted) {
    // Browning: blotches where the elements and the belt's wires were hottest, and the edges darkest.
    for (let i = 0; i < 90; i++) {
      const px = r() * S;
      const py = r() * S;
      const rad = 18 + r() * 60;
      const g = x.createRadialGradient(px, py, 0, px, py, rad);
      g.addColorStop(0, `rgba(122,62,20,${0.18 + r() * 0.22})`);
      g.addColorStop(1, "rgba(122,62,20,0)");
      x.fillStyle = g;
      x.fillRect(px - rad, py - rad, rad * 2, rad * 2);
    }
    x.globalCompositeOperation = "multiply";
    const g = x.createRadialGradient(S / 2, S * 0.52, S * 0.12, S / 2, S * 0.52, S * 0.62);
    g.addColorStop(0, "rgb(255,246,226)");
    g.addColorStop(1, "rgb(150,96,56)");
    x.fillStyle = g;
    x.fillRect(0, 0, S, S);
    x.globalCompositeOperation = "source-over";
  }
  // The crust band, soft on its inner side.
  x.strokeStyle = toasted ? "rgba(92,48,18,0.35)" : "rgba(196,144,84,0.4)";
  x.lineWidth = 34;
  tracePath(x, pts, b, S);
  x.stroke();
  x.strokeStyle = crust;
  x.lineWidth = 18;
  tracePath(x, pts, b, S);
  x.stroke();
  xb.strokeStyle = "#5a5a5a";
  xb.lineWidth = 18;
  tracePath(xb, pts, b, S);
  xb.stroke();
  x.restore();
  xb.restore();
  return { map: ctex(c), bump: ctex(cb, false) };
}

/** One fry along u (its length) and round v: golden, its tips and edges browned, flecked where it blistered. */
function fryArt(): HTMLCanvasElement {
  const [c, x] = canvas2d(128, 32);
  x.fillStyle = "#e8b45e";
  x.fillRect(0, 0, 128, 32);
  const g = x.createLinearGradient(0, 0, 128, 0);
  g.addColorStop(0, "rgba(160,96,40,0.6)");
  g.addColorStop(0.1, "rgba(196,126,52,0.3)");
  g.addColorStop(0.3, "rgba(196,126,52,0)");
  g.addColorStop(0.72, "rgba(196,126,52,0)");
  g.addColorStop(0.9, "rgba(196,126,52,0.3)");
  g.addColorStop(1, "rgba(160,96,40,0.6)");
  x.fillStyle = g;
  x.fillRect(0, 0, 128, 32);
  // The four long edges (the section's rounded corners, a band of v each) brown first.
  for (const v of [6, 14, 22, 30]) {
    const g2 = x.createLinearGradient(0, v - 3, 0, v + 3);
    g2.addColorStop(0, "rgba(184,110,40,0)");
    g2.addColorStop(0.5, "rgba(184,110,40,0.32)");
    g2.addColorStop(1, "rgba(184,110,40,0)");
    x.fillStyle = g2;
    x.fillRect(0, v - 3, 128, 6);
  }
  const r = rng(93);
  for (let i = 0; i < 160; i++) {
    x.fillStyle = r() < 0.7 ? `rgba(170,96,34,${0.15 + r() * 0.3})` : `rgba(255,236,180,${0.2 + r() * 0.3})`;
    x.beginPath();
    x.ellipse(r() * 128, r() * 32, 1 + r() * 3, 0.6 + r() * 1.2, 0, 0, Math.PI * 2);
    x.fill();
  }
  return c;
}

/**
 * One fry along x, `L` long and `a` across: a square section with rounded edges, narrowing toward its ends, which are
 * cut square with soft corners; bowed by `bow` across its length and rolled by `roll`. u runs along it, v round it.
 */
function fryGeom(L: number, a: number, bow: number, roll: number): THREE.BufferGeometry {
  const nL = 5;
  const sec: [number, number, number, number][] = [];
  const c = a * 0.13;
  for (let k = 0; k < 4; k++) {
    const f = (k * Math.PI) / 2 + roll;
    const [ny, nz] = [Math.cos(f), Math.sin(f)];
    const [ty, tz] = [-nz, ny];
    for (const sg of [-1, 1]) {
      const m = Math.hypot(1, 0.22);
      sec.push([(ny * a) / 2 + sg * ty * (a / 2 - c), (nz * a) / 2 + sg * tz * (a / 2 - c), (ny + sg * 0.22 * ty) / m, (nz + sg * 0.22 * tz) / m]);
    }
  }
  sec.push(sec[0]!);
  const pos: number[] = [];
  const nrm: number[] = [];
  const uv: number[] = [];
  const idx: number[] = [];
  const ring = sec.length;
  for (let i = 0; i <= nL; i++) {
    const u = i / nL;
    const e = 2 * u - 1;
    const sc = 1 - 0.16 * e ** 4;
    const x = e * (L / 2);
    const y = bow * (1 - e * e);
    for (const [j, [py, pz, qy, qz]] of sec.entries()) {
      pos.push(x, y + py * sc, pz * sc);
      nrm.push(0, qy, qz);
      uv.push(0.04 + u * 0.92, j / (ring - 1));
    }
  }
  for (let i = 0; i < nL; i++)
    for (let j = 0; j < ring - 1; j++) {
      const a0 = i * ring + j;
      const b0 = a0 + ring;
      // The section runs anticlockwise from +y toward +z, so (along the section) × (along the fry) faces out.
      idx.push(a0, a0 + 1, b0, a0 + 1, b0 + 1, b0);
    }
  // The cut ends: a fan from a centre pushed a little out.
  for (const sg of [-1, 1]) {
    const sc = 1 - 0.16;
    const base = pos.length / 3;
    pos.push(sg * (L / 2) + sg * a * 0.12, 0, 0);
    nrm.push(sg, 0, 0);
    uv.push(sg < 0 ? 0.16 : 0.84, 0.45);
    for (const [py, pz] of sec) {
      pos.push(sg * (L / 2), py * sc, pz * sc);
      nrm.push(sg * 0.8, py / a, pz / a);
      uv.push(sg < 0 ? 0.1 : 0.9, 0.45);
    }
    for (let j = 0; j < ring - 1; j++) sg > 0 ? idx.push(base, base + 1 + j, base + 2 + j) : idx.push(base, base + 2 + j, base + 1 + j);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("normal", new THREE.Float32BufferAttribute(nrm, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  g.setIndex(idx);
  return g;
}

/** A basket of fries from above: a jumble of sticks, the ones beneath in shadow. */
function friesHeapArt(): HTMLCanvasElement {
  const S = 256;
  const [c, x] = canvas2d(S, S);
  x.fillStyle = "#8a5a26";
  x.fillRect(0, 0, S, S);
  const r = rng(94);
  for (let layer = 0; layer < 3; layer++)
    for (let i = 0; i < 150; i++) {
      const len = 26 + r() * 26;
      const wid = 5 + r() * 2;
      x.save();
      x.translate(r() * S, r() * S);
      x.rotate(r() * Math.PI);
      x.shadowColor = "rgba(40,20,5,0.55)";
      x.shadowBlur = 3;
      const k = 0.55 + layer * 0.2 + r() * 0.12;
      const g = x.createLinearGradient(-len / 2, 0, len / 2, 0);
      const base = [246 * k, 205 * k, 128 * k].map((v) => Math.min(255, v | 0)).join(",");
      const tip = [190 * k, 112 * k, 44 * k].map((v) => Math.min(255, v | 0)).join(",");
      g.addColorStop(0, `rgb(${tip})`);
      g.addColorStop(0.2, `rgb(${base})`);
      g.addColorStop(0.8, `rgb(${base})`);
      g.addColorStop(1, `rgb(${tip})`);
      x.fillStyle = g;
      rrect(x, -len / 2, -wid / 2, len, wid, 2);
      x.fill();
      x.restore();
    }
  return c;
}

/**
 * A flat food cut from `shape` (a slice, a steak, a fillet), `t` thick and centred on y = 0, lying in xz with the
 * shape's +y toward -z; its edges rounded by `bev`, and its whole surface painted from one picture over the shape's
 * bounds.
 */
function flatFood(s: Shop, shape: THREE.Shape, t: number, bev: number): THREE.BufferGeometry {
  const b = new THREE.Box2().setFromPoints(shape.getPoints(30));
  const fine = s.fine(1);
  const bv = fine ? Math.min(bev, t * 0.45) : 0;
  const g = new THREE.ExtrudeGeometry(shape, {
    depth: Math.max(1e-4, t - 2 * bv),
    bevelEnabled: bv > 0,
    bevelThickness: bv,
    bevelSize: bv * 0.85,
    bevelSegments: s.fine(2) ? 3 : 1,
    curveSegments: s.fine(2) ? 8 : fine ? 3 : 2,
  });
  planarUV(g, b);
  g.translate(0, 0, -(t - 2 * bv) / 2);
  g.rotateX(-Math.PI / 2);
  return g;
}

/** `pts` (a closed outline) pushed out by `o` along its vertex normals. */
function offsetOutline(pts: THREE.Vector2[], o: number): THREE.Vector2[] {
  const q = pts.length > 2 && pts[0]!.distanceTo(pts[pts.length - 1]!) < 1e-7 ? pts.slice(0, -1) : pts;
  const sg = THREE.ShapeUtils.isClockWise(q) ? -1 : 1;
  const n = q.length;
  return q.map((p, i) => {
    const a = q[(i + n - 1) % n]!;
    const b = q[(i + 1) % n]!;
    const e0 = new THREE.Vector2(p.x - a.x, p.y - a.y).normalize();
    const e1 = new THREE.Vector2(b.x - p.x, b.y - p.y).normalize();
    const nn = new THREE.Vector2(e0.y + e1.y, -(e0.x + e1.x)).normalize().multiplyScalar(sg);
    return new THREE.Vector2(p.x + nn.x * o, p.y + nn.y * o);
  });
}

/** Signed distance from (x, y) to a closed polygon: negative inside. */
function polyDistance(q: THREE.Vector2[], x: number, y: number): number {
  let d2 = Infinity;
  let inside = false;
  for (let i = 0, j = q.length - 1; i < q.length; j = i++) {
    const [a, b] = [q[j]!, q[i]!];
    const [ex, ey] = [b.x - a.x, b.y - a.y];
    const k = Math.max(0, Math.min(1, ((x - a.x) * ex + (y - a.y) * ey) / (ex * ex + ey * ey || 1)));
    const [dx, dy] = [x - a.x - ex * k, y - a.y - ey * k];
    d2 = Math.min(d2, dx * dx + dy * dy);
    if (a.y > y !== b.y > y && x < a.x + ((y - a.y) * ex) / ey) inside = !inside;
  }
  return inside ? -Math.sqrt(d2) : Math.sqrt(d2);
}

interface BagSpec {
  /** The bag's rectangle in the food's frame (the food lies in xz, centred on y = 0). */
  x0: number;
  x1: number;
  z0: number;
  z1: number;
  /** The heat seal: a band across the bag, at x in [a, b] or at z in [a, b]. */
  seal: { axis: "x" | "z"; a: number; b: number };
  /** A sideways set (in y) of the bag by z, for a bag hung over a wall. */
  bend?: (z: number) => number;
  /** Lying on a floor: the film underneath is pressed flat and needs no shape. */
  flat?: boolean;
}

/**
 * A vacuum bag round a food that flatFood cut from `shape` (`t` thick, edges rounded by `bev`). Close to, the bag's two
 * films are height fields over its rectangle: drawn tight over the food and down round its edge to meet a few
 * millimetres outside it, then running out together to the bag's cut edges, creased where the vacuum drew them in
 * (the creases radiating from the food and dying away from it) and cockled everywhere else; flat and doubled across
 * the heat seal. Further off, a shell round the food and a sheet.
 */
function vacBag(s: Shop, mat: THREE.Material, shape: THREE.Shape, t: number, bev: number, spec: BagSpec, seed: number): void {
  const { x0, x1, z0, z1, seal, bend } = spec;
  const bendAt = (z: number) => (bend ? bend(z) : 0);
  const inSeal = (x: number, z: number) => {
    const u = seal.axis === "x" ? x : z;
    return sstep(seal.a - 0.004, seal.a, u) * (1 - sstep(seal.b, seal.b + 0.004, u));
  };
  const sealBox = () => {
    const [w, l] = seal.axis === "x" ? [seal.b - seal.a, z1 - z0] : [x1 - x0, seal.b - seal.a];
    const [ax, az] = seal.axis === "x" ? [(seal.a + seal.b) / 2, (z0 + z1) / 2] : [(x0 + x1) / 2, (seal.a + seal.b) / 2];
    s.add(mat, new THREE.PlaneGeometry(w, l).rotateX(-Math.PI / 2), [ax, bendAt(az) + 0.0006, az]);
  };
  if (!s.fine(2)) {
    const pts = shape.getPoints(3);
    const bvF = Math.min(bev, t * 0.45) * 0.85;
    s.add(mat, flatFood(s, new THREE.Shape(offsetOutline(pts, bvF + 0.0012)), t + 0.0018, Math.min(0.006, t * 0.4)));
    const nz = bend ? 8 : 1;
    const g = new THREE.PlaneGeometry(x1 - x0, z1 - z0, 1, nz).rotateX(-Math.PI / 2);
    const pa = g.attributes.position as THREE.BufferAttribute;
    for (let i = 0; i < pa.count; i++) pa.setY(i, bendAt(pa.getZ(i) + (z0 + z1) / 2));
    s.add(mat, g, [(x0 + x1) / 2, 0, (z0 + z1) / 2]);
    sealBox();
    return;
  }
  const q = shape.getPoints(30);
  const qb = new THREE.Box2().setFromPoints(q);
  const [ocx, ocy] = [(qb.min.x + qb.max.x) / 2, (qb.min.y + qb.max.y) / 2];
  const nf = noiseField(seed);
  const H = t / 2 + 0.0009;
  const bvF = Math.min(bev, t * 0.45) * 0.85;
  const r = bvF + Math.max(0.0045, 0.5 * H);
  const sep = 0.0002;
  const cell = 0.004;
  const nx = Math.max(2, Math.round((x1 - x0) / cell));
  const nz = Math.max(2, Math.round((z1 - z0) / cell));
  const top: number[] = [];
  const bot: number[] = [];
  const uv: number[] = [];
  for (let j = 0; j <= nz; j++)
    for (let i = 0; i <= nx; i++) {
      const x = x0 + ((x1 - x0) * i) / nx;
      const z = z0 + ((z1 - z0) * j) / nz;
      // The food's frame puts the shape's +y toward -z.
      const d = polyDistance(q, x, -z);
      let h = 0;
      if (d <= 0) h = H + 0.0003 * fbm(nf, x * 90, z * 90, 2);
      else if (d < r) h = H * Math.pow(1 - Math.pow(d / r, 2.4), 1 / 2.4);
      // The films together: the creases radiate from the food, the cockle is everywhere off it.
      let wv = 0;
      if (d > 0) {
        const a = Math.atan2(-z - ocy, x - ocx);
        const ridge = 1 - Math.abs(2 * fbm(nf, Math.cos(a) * 7 + d * 9, Math.sin(a) * 7 + d * 6, 3) - 1);
        const near = sstep(r * 0.5, r + 0.004, d) * Math.exp(-Math.max(0, d - r) / 0.03);
        wv = 0.0017 * ridge ** 4 * near + 0.0006 * (fbm(nf, x * 26 + 3, z * 26 + 7, 3) - 0.5) * sstep(r, r + 0.012, d);
        wv *= 1 - inSeal(x, z);
      }
      const b = bendAt(z);
      top.push(x, b + wv + h + sep, z);
      bot.push(x, b + wv - h - sep, z);
      uv.push(x / 0.04, z / 0.04);
    }
  const idx: number[] = [];
  const idxR: number[] = [];
  for (let j = 0; j < nz; j++)
    for (let i = 0; i < nx; i++) {
      const a = j * (nx + 1) + i;
      const [b, c, e] = [a + 1, a + nx + 1, a + nx + 2];
      idx.push(a, c, b, b, c, e);
      idxR.push(a, b, c, b, e, c);
    }
  const film = (pos: number[], ix: number[]) => {
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
    g.setIndex(ix);
    g.computeVertexNormals();
    return g;
  };
  s.add(mat, film(top, idx));
  if (spec.flat) s.add(mat, new THREE.PlaneGeometry(x1 - x0, z1 - z0).rotateX(Math.PI / 2), [(x0 + x1) / 2, -H - sep, (z0 + z1) / 2]);
  else s.add(mat, film(bot, idxR));
  sealBox();
}

/** The outline of a ribeye steak lying flat: the eye, its cap round the top, the fat tail at the left; centred, metres. */
function steakShape(): THREE.Shape {
  const s = new THREE.Shape();
  s.moveTo(-0.072, -0.03);
  s.bezierCurveTo(-0.05, -0.058, 0.02, -0.06, 0.058, -0.04);
  s.bezierCurveTo(0.084, -0.024, 0.088, 0.012, 0.072, 0.034);
  s.bezierCurveTo(0.054, 0.058, 0.012, 0.062, -0.022, 0.054);
  s.bezierCurveTo(-0.05, 0.048, -0.074, 0.034, -0.084, 0.012);
  s.bezierCurveTo(-0.09, -0.004, -0.084, -0.02, -0.072, -0.03);
  return s;
}

/**
 * A ribeye's face over its outline's bounds, as it looks part way through its bath: the eye a rosy brown-pink, finely
 * marbled; the spinalis cap round its upper edge, parted from it by a seam of fat; the fat cap along the top and the
 * tail's fat at the left, both continued past the outline so that the sides there read as fat.
 */
function steakArt(): HTMLCanvasElement {
  const S = 512;
  const pts = steakShape().getPoints(30);
  const b = new THREE.Box2().setFromPoints(pts);
  const [c, x] = canvas2d(S, S);
  const r = rng(71);
  x.fillStyle = "#a2544c";
  x.fillRect(0, 0, S, S);
  // The fat cap along the top and round to the tail, drawn wide so that it runs past the outline.
  x.strokeStyle = "#eadcc6";
  x.lineCap = "round";
  x.lineWidth = 40;
  x.beginPath();
  x.moveTo(S * 0.88, S * 0.18);
  x.bezierCurveTo(S * 0.62, -S * 0.02, S * 0.22, S * 0.0, S * 0.03, S * 0.36);
  x.stroke();
  x.save();
  tracePath(x, pts, b, S);
  x.clip();
  // Flesh: mottled so it is not a flat colour.
  for (let i = 0; i < 260; i++) {
    const px = r() * S;
    const py = r() * S;
    const rad = 10 + r() * 40;
    const g = x.createRadialGradient(px, py, 0, px, py, rad);
    const k = r();
    g.addColorStop(0, k < 0.5 ? "rgba(178,92,84,0.35)" : "rgba(128,52,50,0.3)");
    g.addColorStop(1, "rgba(150,70,64,0)");
    x.fillStyle = g;
    x.fillRect(px - rad, py - rad, rad * 2, rad * 2);
  }
  // The seam round the spinalis, and the cap's fat.
  x.strokeStyle = "rgba(238,224,204,0.7)";
  x.shadowColor = "rgba(238,224,204,0.8)";
  x.shadowBlur = 6;
  x.lineWidth = 9;
  x.beginPath();
  x.moveTo(S * 0.06, S * 0.52);
  x.bezierCurveTo(S * 0.2, S * 0.3, S * 0.6, S * 0.22, S * 0.95, S * 0.42);
  x.stroke();
  x.strokeStyle = "#eadcc6";
  x.lineWidth = 30;
  x.beginPath();
  x.moveTo(S * 0.98, S * 0.12);
  x.bezierCurveTo(S * 0.62, -S * 0.05, S * 0.2, -S * 0.03, -S * 0.02, S * 0.36);
  x.stroke();
  // Marbling: intramuscular fat as a fine, soft web of flecks and short streaks over the whole eye, a little denser
  // toward its lower half, sparse in the cap.
  x.shadowColor = "rgba(238,224,204,0.7)";
  x.shadowBlur = 4;
  for (let i = 0; i < 1700; i++) {
    const px = S * (0.08 + r() * 0.88);
    const py = S * (0.3 + Math.sqrt(r()) * 0.68);
    x.fillStyle = `rgba(240,227,208,${0.14 + r() * 0.32})`;
    x.beginPath();
    x.ellipse(px, py, 0.8 + r() * 2.4, 0.5 + r() * 1.1, r() * Math.PI, 0, Math.PI * 2);
    x.fill();
  }
  for (let i = 0; i < 170; i++) {
    let px = S * (0.1 + r() * 0.84);
    let py = S * (0.32 + Math.sqrt(r()) * 0.64);
    let a = r() * Math.PI * 2;
    const n = 2 + Math.floor(r() * 5);
    x.strokeStyle = `rgba(240,226,206,${0.18 + r() * 0.3})`;
    x.lineWidth = 0.9 + r() * 1.8;
    x.beginPath();
    x.moveTo(px, py);
    for (let k = 0; k < n; k++) {
      a += (r() - 0.5) * 1.1;
      px += Math.cos(a) * (3 + r() * 7);
      py += Math.sin(a) * (3 + r() * 7);
      x.lineTo(px, py);
    }
    x.stroke();
  }
  x.shadowBlur = 0;
  // Fat tail at the left.
  const g = x.createRadialGradient(S * 0.02, S * 0.55, 0, S * 0.02, S * 0.55, S * 0.16);
  g.addColorStop(0, "rgba(236,224,204,1)");
  g.addColorStop(0.75, "rgba(236,224,204,0.85)");
  g.addColorStop(1, "rgba(236,224,204,0)");
  x.fillStyle = g;
  x.fillRect(0, 0, S, S);
  x.restore();
  return c;
}

/** The outline of a salmon portion cut across a side: rounded at the thick back, tapering to the belly. Centred, metres. */
function salmonShape(): THREE.Shape {
  const s = new THREE.Shape();
  s.moveTo(-0.08, -0.03);
  s.bezierCurveTo(-0.03, -0.04, 0.05, -0.036, 0.078, -0.022);
  s.bezierCurveTo(0.09, -0.012, 0.09, 0.016, 0.076, 0.026);
  s.bezierCurveTo(0.04, 0.042, -0.04, 0.04, -0.078, 0.03);
  s.bezierCurveTo(-0.088, 0.018, -0.088, -0.022, -0.08, -0.03);
  return s;
}

/** A salmon portion's face: coral, paler toward the belly, its myosepta in white chevrons, the brown lateral line. */
function salmonArt(): HTMLCanvasElement {
  const W = 512;
  const H = 256;
  const [c, x] = canvas2d(W, H);
  const g = x.createLinearGradient(0, 0, 0, H);
  g.addColorStop(0, "#e2734c");
  g.addColorStop(0.55, "#ec8a60");
  g.addColorStop(1, "#f2a27e");
  x.fillStyle = g;
  x.fillRect(0, 0, W, H);
  x.strokeStyle = "rgba(168,84,58,0.6)";
  x.lineWidth = 7;
  x.beginPath();
  x.moveTo(0, H * 0.42);
  x.bezierCurveTo(W * 0.3, H * 0.46, W * 0.7, H * 0.4, W, H * 0.44);
  x.stroke();
  const r = rng(72);
  for (let i = 0; i < 15; i++) {
    const x0 = -W * 0.04 + i * W * 0.075 + (r() - 0.5) * 6;
    x.strokeStyle = `rgba(250,232,218,${0.7 + r() * 0.25})`;
    x.lineWidth = 2.5 + r() * 2.5;
    x.beginPath();
    x.moveTo(x0 + 22, -4);
    x.quadraticCurveTo(x0 - 6, H * 0.2, x0 + 14, H * 0.42);
    x.quadraticCurveTo(x0 + 34, H * 0.66, x0 + 6, H + 4);
    x.stroke();
  }
  return c;
}

/** Vacuum-bag film's embossed channels: a diamond lattice, as a bump map. */
function bagEmboss(): THREE.CanvasTexture {
  const [c, x] = canvas2d(64, 64);
  x.fillStyle = "#808080";
  x.fillRect(0, 0, 64, 64);
  x.strokeStyle = "#b8b8b8";
  x.lineWidth = 3;
  for (const k of [-64, 0, 64]) {
    x.beginPath();
    x.moveTo(k, 0);
    x.lineTo(k + 64, 64);
    x.moveTo(k + 64, 0);
    x.lineTo(k, 64);
    x.stroke();
  }
  const t = ctex(c, false);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(14, 14);
  return t;
}
/* ------------------------------------------------------------------------------------------------------------ */
/* Hot-table food and its board, each painted on one tile that repeats seamlessly, laid on in metres              */
/* ------------------------------------------------------------------------------------------------------------ */

/** Run `f` at (x, y) and at each wrapped copy that reaches into an S-pixel tile within `reach`, so the tile repeats. */
function wrapped(S: number, x: number, y: number, reach: number, f: (x: number, y: number) => void): void {
  for (const ox of [-S, 0, S])
    for (const oy of [-S, 0, S]) {
      const [px, py] = [x + ox, y + oy];
      if (px > -reach && px < S + reach && py > -reach && py < S + reach) f(px, py);
    }
}
const grey = (v: number) => {
  const k = Math.max(0, Math.min(255, Math.round(v)));
  return `rgb(${k},${k},${k})`;
};
const tint = (rgb: readonly number[], k: number) => `rgb(${rgb.map((v) => Math.max(0, Math.min(255, Math.round(v * k)))).join(",")})`;
const mixRgb = (a: readonly number[], b: readonly number[], t: number) => a.map((v, i) => v + (b[i]! - v) * t);

/** Steamed long-grain rice, a tile 80 mm across (0.156 mm a pixel): three layers of grains, those beneath in shadow. */
function riceArt(): [HTMLCanvasElement, HTMLCanvasElement] {
  const S = 512;
  const [c, x] = canvas2d(S, S);
  const [b, y] = canvas2d(S, S);
  x.fillStyle = "#958b77";
  x.fillRect(0, 0, S, S);
  y.fillStyle = "#000";
  y.fillRect(0, 0, S, S);
  y.globalCompositeOperation = "lighten";
  const r = rng(73);
  for (let layer = 0; layer < 3; layer++)
    for (let i = 0; i < 560; i++) {
      const [gx, gy, a] = [r() * S, r() * S, r() * Math.PI];
      // A cooked grain is 7 to 9 mm long and about 2.3 mm across.
      const len = 45 + r() * 13;
      const wid = 14 + r() * 3;
      const k = 0.84 + layer * 0.07 + r() * 0.04;
      const hgt = 80 + layer * 60 + r() * 30;
      wrapped(S, gx, gy, len, (px, py) => {
        for (const ctx of [x, y]) {
          ctx.save();
          ctx.translate(px, py);
          ctx.rotate(a);
          ctx.scale(1, wid / len);
          const g = ctx.createRadialGradient(-len * 0.06, -len * 0.1, 0, 0, 0, len / 2);
          if (ctx === y) {
            g.addColorStop(0, grey(hgt));
            g.addColorStop(1, grey(hgt * 0.45));
          } else {
            g.addColorStop(0, tint([255, 254, 249], k));
            g.addColorStop(0.6, tint([243, 239, 228], k));
            g.addColorStop(1, tint([200, 192, 172], k));
          }
          ctx.fillStyle = g;
          ctx.beginPath();
          ctx.arc(0, 0, len / 2, 0, Math.PI * 2);
          ctx.fill();
          ctx.restore();
        }
      });
    }
  return [c, b];
}

/** An irregular crumb of `n` corners round (px, py), its radii `rs` times `rad`, turned by `a0`. */
function crumb(ctx: CanvasRenderingContext2D, px: number, py: number, rad: number, rs: readonly number[], a0: number): void {
  ctx.beginPath();
  for (const [k, q] of rs.entries()) {
    const a = a0 + (k / rs.length) * Math.PI * 2;
    if (k) ctx.lineTo(px + Math.cos(a) * rad * q, py + Math.sin(a) * rad * q);
    else ctx.moveTo(px + Math.cos(a) * rad * q, py + Math.sin(a) * rad * q);
  }
  ctx.closePath();
  ctx.fill();
}

/** Chili con carne, a tile 120 mm across (0.234 mm a pixel): mince in a red-brown sauce, onion and tomato through it. */
function chiliArt(): [HTMLCanvasElement, HTMLCanvasElement] {
  const S = 512;
  const [c, x] = canvas2d(S, S);
  const [b, y] = canvas2d(S, S);
  x.fillStyle = "#47190b";
  x.fillRect(0, 0, S, S);
  y.fillStyle = grey(40);
  y.fillRect(0, 0, S, S);
  const r = rng(37);
  const pool = (n: number, alpha: number) => {
    for (let i = 0; i < n; i++) {
      const [px0, py0, rad] = [r() * S, r() * S, 30 + r() * 70];
      const col = r() < 0.5 ? "116,40,14" : "50,17,7";
      wrapped(S, px0, py0, rad, (px, py) => {
        const g = x.createRadialGradient(px, py, 0, px, py, rad);
        g.addColorStop(0, `rgba(${col},${alpha})`);
        g.addColorStop(1, `rgba(${col},0)`);
        x.fillStyle = g;
        x.fillRect(px - rad, py - rad, rad * 2, rad * 2);
      });
    }
  };
  pool(40, 0.45);
  y.globalCompositeOperation = "lighten";
  // Mince: browned crumbs 2 to 5 mm, lit from above, sauce dark in the cracks between them.
  const cols = [
    [86, 36, 16],
    [104, 46, 22],
    [70, 27, 12],
    [120, 56, 28],
  ];
  for (let i = 0; i < 1700; i++) {
    const [px0, py0, a0] = [r() * S, r() * S, r() * 6.28];
    const rad = 5 + r() * 8;
    const rs = Array.from({ length: 7 }, () => 0.6 + r() * 0.5);
    const col = cols[Math.floor(r() * cols.length)]!;
    const hgt = 110 + r() * 110;
    wrapped(S, px0, py0, rad * 1.2, (px, py) => {
      const g = x.createRadialGradient(px - rad * 0.3, py - rad * 0.3, 0, px, py, rad * 1.1);
      g.addColorStop(0, tint(col, 1.2));
      g.addColorStop(0.7, tint(col, 0.9));
      g.addColorStop(1, tint(col, 0.5));
      x.fillStyle = g;
      crumb(x, px, py, rad, rs, a0);
      const gh = y.createRadialGradient(px, py, 0, px, py, rad * 1.1);
      gh.addColorStop(0, grey(hgt));
      gh.addColorStop(1, grey(hgt * 0.5));
      y.fillStyle = gh;
      crumb(y, px, py, rad, rs, a0);
    });
  }
  // The sauce runs back over some of it.
  pool(45, 0.4);
  // Onion gone soft and translucent, tomato, flecks of chili.
  for (const [n, col, r0, r1] of [
    [70, "rgba(200,150,100,0.3)", 4, 8],
    [140, "rgba(150,30,12,0.55)", 3, 7],
    [260, "rgba(40,10,4,0.7)", 1, 2.2],
  ] as const)
    for (let i = 0; i < n; i++) {
      const [px0, py0, a0] = [r() * S, r() * S, r() * 6.28];
      const rad = r0 + r() * (r1 - r0);
      const rs = Array.from({ length: 5 }, () => 0.7 + r() * 0.35);
      wrapped(S, px0, py0, rad * 1.2, (px, py) => {
        x.fillStyle = col;
        crumb(x, px, py, rad, rs, a0);
      });
    }
  return [c, b];
}

/**
 * Baked macaroni and cheese, a tile 100 mm across (0.195 mm a pixel): elbows in a golden cheese sauce, three deep,
 * browned in patches where the top caught under the heat, each elbow's open end a dark hollow.
 */
function macArt(): [HTMLCanvasElement, HTMLCanvasElement] {
  const S = 512;
  const [c, x] = canvas2d(S, S);
  const [b, y] = canvas2d(S, S);
  x.fillStyle = "#a8641a";
  x.fillRect(0, 0, S, S);
  y.fillStyle = grey(30);
  y.fillRect(0, 0, S, S);
  const r = rng(53);
  const patches = Array.from({ length: 24 }, () => [r() * S, r() * S, 30 + r() * 70] as const);
  const brown = (px: number, py: number) => {
    let v = 0;
    for (const [bx, by, br] of patches) {
      let dx = Math.abs(px - bx) % S;
      let dy = Math.abs(py - by) % S;
      dx = Math.min(dx, S - dx);
      dy = Math.min(dy, S - dy);
      v = Math.max(v, 1 - Math.hypot(dx, dy) / br);
    }
    return sstep(0, 0.8, v);
  };
  const sauce = (n: number, alpha: number) => {
    for (let i = 0; i < n; i++) {
      const [px0, py0, rad] = [r() * S, r() * S, 14 + r() * 30];
      const col = tint(mixRgb([222, 146, 40], [160, 86, 28], brown(px0, py0)), 1);
      wrapped(S, px0, py0, rad, (px, py) => {
        const g = x.createRadialGradient(px, py, 0, px, py, rad);
        g.addColorStop(0, col.replace("rgb", "rgba").replace(")", `,${alpha})`));
        g.addColorStop(1, col.replace("rgb", "rgba").replace(")", ",0)"));
        x.fillStyle = g;
        x.fillRect(px - rad, py - rad, rad * 2, rad * 2);
      });
    }
  };
  x.lineCap = "round";
  y.lineCap = "round";
  for (let layer = 0; layer < 3; layer++) {
    for (let i = 0; i < 170; i++) {
      const [ex, ey, a0] = [r() * S, r() * S, r() * 6.28];
      // An elbow's centreline bends at 5 to 6.5 mm; the tube is about 4.3 mm across.
      const R = 26 + r() * 8;
      const lw = 20 + r() * 4;
      const sweep = 1.7 + r() * 0.9;
      const k = 0.78 + layer * 0.11 + r() * 0.05;
      const col = mixRgb([238, 180, 72], [176, 100, 36], brown(ex, ey) * (0.5 + layer * 0.25));
      const hgt = 90 + layer * 55;
      const hole = r() < 0.2;
      wrapped(S, ex, ey, R + lw, (px, py) => {
        const arc = (ctx: CanvasRenderingContext2D, rr: number, w: number, style: string) => {
          ctx.strokeStyle = style;
          ctx.lineWidth = w;
          ctx.beginPath();
          ctx.arc(px, py, rr, a0, a0 + sweep);
          ctx.stroke();
        };
        arc(x, R, lw + 5, "rgba(104,56,18,0.55)");
        arc(x, R, lw, tint(col, k));
        arc(x, R - lw * 0.18, lw * 0.26, "rgba(255,226,150,0.3)");
        // An open end turned up, seen into.
        if (hole) {
          x.fillStyle = "rgba(122,66,20,0.45)";
          x.beginPath();
          x.ellipse(px + Math.cos(a0) * R, py + Math.sin(a0) * R, lw * 0.22, lw * 0.13, a0, 0, Math.PI * 2);
          x.fill();
        }
        arc(y, R, lw, grey(hgt));
        arc(y, R, lw * 0.45, grey(hgt + 60));
      });
    }
    sauce(layer < 2 ? 70 : 30, layer < 2 ? 0.5 : 0.3);
  }
  // Where it caught: the crust darker still at the patches' hearts.
  for (const [bx, by, br] of patches)
    wrapped(S, bx, by, br, (px, py) => {
      const g = x.createRadialGradient(px, py, 0, px, py, br * 0.6);
      g.addColorStop(0, "rgba(110,52,14,0.45)");
      g.addColorStop(1, "rgba(110,52,14,0)");
      x.fillStyle = g;
      x.fillRect(px - br, py - br, br * 2, br * 2);
    });
  return [c, b];
}

/** Whole cooked green beans heaped, a tile 200 mm across (0.39 mm a pixel): three layers, those beneath in shadow. */
function greenBeanArt(): [HTMLCanvasElement, HTMLCanvasElement] {
  const S = 512;
  const [c, x] = canvas2d(S, S);
  const [b, y] = canvas2d(S, S);
  x.fillStyle = "#18230a";
  x.fillRect(0, 0, S, S);
  y.fillStyle = "#000";
  y.fillRect(0, 0, S, S);
  x.lineCap = "round";
  y.lineCap = "round";
  const r = rng(61);
  for (let layer = 0; layer < 3; layer++)
    for (let i = 0; i < 64; i++) {
      // A bean is 60 to 100 mm long and 7 to 8 mm across.
      const [px0, py0, a] = [r() * S, r() * S, r() * Math.PI];
      const L = 150 + r() * 110;
      const wd = 17 + r() * 4;
      const bow = (r() - 0.5) * 34;
      const k = 0.5 + layer * 0.24 + r() * 0.08;
      const [ca, sa] = [Math.cos(a), Math.sin(a)];
      wrapped(S, px0, py0, L / 2 + wd, (px, py) => {
        const path = (ctx: CanvasRenderingContext2D, off: number) => {
          ctx.beginPath();
          ctx.moveTo(px - (ca * L) / 2 - sa * off, py - (sa * L) / 2 + ca * off);
          ctx.quadraticCurveTo(px - sa * (bow + off), py + ca * (bow + off), px + (ca * L) / 2 - sa * off, py + (sa * L) / 2 + ca * off);
        };
        for (const [ctx, w, style] of [
          [x, wd + 4, "rgba(14,26,6,0.75)"],
          [x, wd, tint([86, 132, 46], k)],
          [y, wd, grey(70 + layer * 60)],
          [y, wd * 0.5, grey(110 + layer * 60)],
        ] as const) {
          ctx.strokeStyle = style;
          ctx.lineWidth = w;
          path(ctx, 0);
          ctx.stroke();
        }
        x.strokeStyle = "rgba(176,210,120,0.32)";
        x.lineWidth = 4;
        path(x, -wd * 0.22);
        x.stroke();
        // The stem end, trimmed.
        x.fillStyle = "rgba(120,112,52,0.8)";
        x.beginPath();
        x.arc(px - (ca * L) / 2, py - (sa * L) / 2, wd * 0.28, 0, Math.PI * 2);
        x.fill();
      });
    }
  return [c, b];
}

/** Corn kernels, a tile 80 mm across (0.156 mm a pixel): plump and glossy, pale where they were cut from the cob. */
function cornArt(): [HTMLCanvasElement, HTMLCanvasElement] {
  const S = 512;
  const [c, x] = canvas2d(S, S);
  const [b, y] = canvas2d(S, S);
  x.fillStyle = "#6a3f04";
  x.fillRect(0, 0, S, S);
  y.fillStyle = "#000";
  y.fillRect(0, 0, S, S);
  y.globalCompositeOperation = "lighten";
  const r = rng(29);
  for (let layer = 0; layer < 3; layer++)
    for (let i = 0; i < 150; i++) {
      // A kernel is about 7 mm across and 8 long.
      const [px0, py0, a] = [r() * S, r() * S, r() * 6.28];
      const [kw, kh] = [40 + r() * 8, 46 + r() * 8];
      const k = 0.8 + layer * 0.1 + r() * 0.05;
      wrapped(S, px0, py0, kh, (px, py) => {
        for (const ctx of [x, y]) {
          ctx.save();
          ctx.translate(px, py);
          ctx.rotate(a);
          const g = ctx.createRadialGradient(-kw * 0.12, -kh * 0.12, 0, 0, 0, kh * 0.62);
          if (ctx === y) {
            g.addColorStop(0, grey(110 + layer * 55));
            g.addColorStop(1, grey(30 + layer * 40));
          } else {
            g.addColorStop(0, tint([255, 212, 62], k));
            g.addColorStop(0.55, tint([240, 178, 22], k));
            g.addColorStop(1, tint([192, 124, 8], k));
          }
          ctx.fillStyle = g;
          rrect(ctx, -kw / 2, -kh / 2, kw, kh, kw * 0.38);
          ctx.fill();
          if (ctx === x) {
            x.strokeStyle = "rgba(120,72,4,0.5)";
            x.lineWidth = 2;
            x.stroke();
            // The cut tip, pale; a glint.
            x.fillStyle = "rgba(250,230,160,0.35)";
            x.fillRect(-kw * 0.3, kh / 2 - 7, kw * 0.6, 5);
            x.fillStyle = "rgba(255,250,220,0.4)";
            x.beginPath();
            x.ellipse(-kw * 0.14, -kh * 0.18, 4, 2.5, 0.5, 0, Math.PI * 2);
            x.fill();
          }
          ctx.restore();
        }
      });
    }
  return [c, b];
}

/**
 * A cutting board's HDPE, a tile 300 mm across: a warm off-white, a faint stain or two, crossed by knife scores that
 * catch a little grey, roughen it where they run and are cut into it.
 */
function boardArt(): [HTMLCanvasElement, HTMLCanvasElement, HTMLCanvasElement] {
  const S = 512;
  const [c, x] = canvas2d(S, S);
  const [rc, rx] = canvas2d(S, S);
  const [bc, bx] = canvas2d(S, S);
  const r = rng(97);
  x.fillStyle = "#ddd7c8";
  x.fillRect(0, 0, S, S);
  for (let i = 0; i < 6; i++) {
    const [px0, py0, rad] = [r() * S, r() * S, 60 + r() * 120];
    wrapped(S, px0, py0, rad, (px, py) => {
      const g = x.createRadialGradient(px, py, 0, px, py, rad);
      g.addColorStop(0, "rgba(196,186,166,0.12)");
      g.addColorStop(1, "rgba(196,186,166,0)");
      x.fillStyle = g;
      x.fillRect(px - rad, py - rad, rad * 2, rad * 2);
    });
  }
  rx.fillStyle = grey(160);
  rx.fillRect(0, 0, S, S);
  bx.fillStyle = grey(128);
  bx.fillRect(0, 0, S, S);
  for (let i = 0; i < 900; i++) {
    // Most cuts run along the board, a few across it.
    const a = (r() < 0.6 ? 0 : Math.PI / 2) + (r() - 0.5) * 0.7;
    const len = 10 + r() * 70;
    const [px0, py0] = [r() * S, r() * S];
    const alpha = 0.08 + r() * 0.16;
    wrapped(S, px0, py0, len, (px, py) => {
      for (const [ctx, col, lw] of [
        [x, `rgba(112,104,90,${alpha})`, 0.8],
        [rx, "rgba(240,240,240,0.5)", 1.2],
        [bx, "rgba(40,40,40,0.45)", 0.9],
      ] as const) {
        ctx.strokeStyle = col;
        ctx.lineWidth = lw;
        ctx.beginPath();
        ctx.moveTo(px, py);
        ctx.lineTo(px + Math.cos(a) * len, py + Math.sin(a) * len);
        ctx.stroke();
      }
    });
  }
  return [c, rc, bc];
}

function std(p: THREE.MeshStandardMaterialParameters, uv: "box" | "own" = "box", uvScale?: number): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial(p);
  m.userData.uv = uv;
  if (uvScale) m.userData.uvScale = uvScale;
  return m;
}

/** A 256 px tile of `base` scattered with `n` small dots of the `dots` colours. */
function speck(base: string, dots: string[], n: number, seed: number): THREE.CanvasTexture {
  const [c, x] = canvas2d(256, 256);
  const r = rng(seed);
  x.fillStyle = base;
  x.fillRect(0, 0, 256, 256);
  for (let i = 0; i < n; i++) {
    x.fillStyle = dots[Math.floor(r() * dots.length)]!;
    x.beginPath();
    x.ellipse(r() * 256, r() * 256, 1 + r() * 3.5, 1 + r() * 2.5, r() * 3, 0, Math.PI * 2);
    x.fill();
  }
  const t = ctex(c);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  return t;
}

interface FoodMats {
  /** Food: what the pans hold. Baked mac and cheese, its elbows in the sauce and browned where the top caught. */
  mac: THREE.MeshStandardMaterial;
  /** Loose elbows standing out of the mac: cheese-coated, and browned. */
  elbow: THREE.MeshStandardMaterial;
  elbowBrown: THREE.MeshStandardMaterial;
  /** Chili con carne: a glossy red-brown sauce full of mince, onion and tomato. */
  chili: THREE.MeshStandardMaterial;
  /** Kidney beans in the chili. */
  bean: THREE.MeshStandardMaterial;
  /** Whole cooked green beans: the heap seen from above, and loose pods in two greens. */
  greensHeap: THREE.MeshStandardMaterial;
  greens: THREE.MeshStandardMaterial;
  greensDark: THREE.MeshStandardMaterial;
  /** Steamed long-grain rice, grain by grain. */
  rice: THREE.MeshStandardMaterial;
  mash: THREE.MeshStandardMaterial;
  /** Corn kernels. */
  corn: THREE.MeshStandardMaterial;
  /** A cutting board's HDPE: off-white, satin, crossed by fine knife scores. UVs in metres (a tile is 0.3 m). */
  board: THREE.MeshStandardMaterial;
}

let food: FoodMats | null = null;

/**
 * The hot table's food and its board, painted the first time a hot table is built rather than with the rest of the
 * family's materials: a kitchen without one never paints them. Food tops carry their own UVs in metres; each map
 * repeats every `tile` m.
 */
function foodMats(): FoodMats {
  if (food) return food;
  const tiled = (cv: HTMLCanvasElement, srgb: boolean, tile: number) => {
    const t = ctex(cv, srgb);
    t.wrapS = t.wrapT = THREE.RepeatWrapping;
    t.repeat.set(1 / tile, 1 / tile);
    return t;
  };
  const art = ([c, b]: [HTMLCanvasElement, HTMLCanvasElement], tile: number) => ({ map: tiled(c, true, tile), bumpMap: tiled(b, false, tile) });
  const mashMap = speck("#e4d5ad", ["#efe3c2", "#d8c79c", "#ece0bc"], 1800, 21);
  mashMap.repeat.set(1 / 0.06, 1 / 0.06);
  const mashBump = tiled(mashMap.image as HTMLCanvasElement, false, 0.06);
  const chili = std({ ...art(chiliArt(), 0.12), bumpScale: 4, roughness: 0.32, metalness: 0 }, "own");
  const [bm, br, bb] = boardArt();
  return (food = {
    mac: std({ ...art(macArt(), 0.1), bumpScale: 1.6, roughness: 0.5, metalness: 0 }, "own"),
    // Open tubes: both faces drawn, so an end that stands clear of the sauce shows the pasta's inside.
    elbow: std({ color: 0xe9b04c, roughness: 0.42, metalness: 0, side: THREE.DoubleSide }),
    elbowBrown: std({ color: 0xc27f36, roughness: 0.5, metalness: 0, side: THREE.DoubleSide }),
    chili,
    bean: std({ color: 0x4c1610, roughness: 0.3, metalness: 0 }),
    greensHeap: std({ ...art(greenBeanArt(), 0.2), bumpScale: 2, roughness: 0.48, metalness: 0 }, "own"),
    greens: std({ color: 0x5a8a2e, roughness: 0.42, metalness: 0 }),
    greensDark: std({ color: 0x41661f, roughness: 0.46, metalness: 0 }),
    rice: std({ ...art(riceArt(), 0.08), bumpScale: 2.4, roughness: 0.55, metalness: 0 }, "own"),
    mash: std({ map: mashMap, bumpMap: mashBump, bumpScale: 0.8, roughness: 0.62, metalness: 0 }, "own"),
    corn: std({ ...art(cornArt(), 0.08), bumpScale: 2.2, roughness: 0.34, metalness: 0 }, "own"),
    board: std({ map: tiled(bm, true, 0.3), roughnessMap: tiled(br, false, 0.3), bumpMap: tiled(bb, false, 0.3), bumpScale: 0.5, roughness: 0.85, metalness: 0 }, "own"),
  });
}

/** A slice of sandwich bread, its face and crust as one picture: raw, and toasted golden-brown. */
interface BreadMats {
  breadRaw: THREE.MeshStandardMaterial;
  breadToast: THREE.MeshStandardMaterial;
}

let breads: BreadMats | null = null;

/** The toaster's bread, painted the first time a toaster is built. */
function breadMats(): BreadMats {
  if (breads) return breads;
  const bread = (toasted: boolean) => {
    const a = breadArt(toasted);
    const m = new THREE.MeshStandardMaterial({ map: a.map, bumpMap: a.bump, bumpScale: 2.2, roughness: 0.92, metalness: 0 });
    m.userData.uv = "own";
    return m;
  };
  return (breads = { breadRaw: bread(false), breadToast: bread(true) });
}

let fam: FamilyMats | null = null;

function mats(): FamilyMats {
  if (fam) return fam;
  const basic = (c: HTMLCanvasElement) => {
    const m = new THREE.MeshBasicMaterial({ map: ctex(c) });
    m.userData.uv = "own";
    return m;
  };
  const perfTex = perforation(16, 0.62);
  perfTex.repeat.set(1, 1);
  const meshTex = perforation(32, 0.6);
  // Clear polycarbonate. Everything made of it is a closed solid, so its front faces alone draw both of its surfaces.
  const poly = new THREE.MeshPhysicalMaterial({ color: 0xe6eef2, metalness: 0, roughness: 0.08, transparent: true, opacity: 0.13, clearcoat: 1, envMapIntensity: 1.4, depthWrite: false });
  poly.userData.uv = "own";
  // Vacuum-bag film: thin, glossy and slightly milky; its embossed channels catch the light.
  const bag = new THREE.MeshPhysicalMaterial({
    color: 0xf4f7f8,
    metalness: 0,
    roughness: 0.16,
    transparent: true,
    opacity: 0.11,
    clearcoat: 1,
    clearcoatRoughness: 0.06,
    envMapIntensity: 1.8,
    bumpMap: bagEmboss(),
    bumpScale: 0.6,
    depthWrite: false,
    side: THREE.DoubleSide,
  });
  bag.userData.uv = "own";
  const ui = basic(combiUi());
  const own = <T extends THREE.Material>(m: T): T => {
    m.userData.uv = "own";
    return m;
  };
  const fryTex = ctex(fryArt());
  const grilleTex = perforation(15, 0.62);
  grilleTex.repeat.set(1, 2 / 3);
  fam = {
    perf: std({ color: 0xb9bdc2, metalness: 0.85, roughness: 0.42, alphaMap: perfTex, alphaTest: 0.5, side: THREE.DoubleSide }, "box", 125 / 16),
    screenMesh: std({ color: 0x101113, metalness: 0.5, roughness: 0.5, alphaMap: meshTex, alphaTest: 0.5, side: THREE.DoubleSide }, "box", 400 / 32),
    combiUi: ui,
    keypad: std({ map: ctex(keypadArt()), roughness: 0.45, metalness: 0 }, "own"),
    mwDisplay: basic(readout("1:30", "#62e0ff", 256, 96)),
    fryDisplay: basic(readout("2:45", "#ff9f3a", 256, 96, "FRY · 200°")),
    holdDisplay: basic(readout("71° 65%", "#ff6a3a", 320, 96)),
    sousDisplay: basic(readout("57.5°", "#ffffff", 256, 96)),
    dialTemp: std({ map: ctex(dialFace(120, 20, 82, "°C")), roughness: 0.5 }, "own"),
    dialPress: std({ map: ctex(dialFace(10, 2, 8.6, "bar")), roughness: 0.5 }, "own"),
    poly,
    bag,
    bath: new THREE.MeshPhysicalMaterial({ color: 0xa8ccd4, metalness: 0, roughness: 0.03, transparent: true, opacity: 0.15, clearcoat: 1, clearcoatRoughness: 0.02, envMapIntensity: 1.3, depthWrite: false }),
    red: std({ color: 0xc8302a, metalness: 0.1, roughness: 0.45 }),
    greyBox: std({ color: 0x8b9096, metalness: 0.2, roughness: 0.55 }),
    mica: std({ color: 0xcdb98e, metalness: 0, roughness: 0.5, side: THREE.DoubleSide }, "box"),
    wool: std({ color: 0xffffff, metalness: 0, roughness: 1, map: speck("#cdc29e", ["#c2b68e", "#d8cfae", "#bfb38a"], 2400, 15) }, "box", 5),
    lamp: (() => {
      const m = new THREE.MeshBasicMaterial({ color: 0xffe0a8 });
      m.userData.uv = "own";
      return m;
    })(),
    fries: own(new THREE.MeshStandardMaterial({ map: ctex(friesHeapArt()), roughness: 0.62 })),
    fry: [0xfff8ec, 0xf8dcae, 0xe8b07a].map((col) => own(new THREE.MeshStandardMaterial({ map: fryTex, color: col, roughness: 0.62 }))),
    dough: std({ color: 0xead8b4, roughness: 0.88 }),
    beef: own(new THREE.MeshStandardMaterial({ map: ctex(steakArt()), roughness: 0.42 })),
    salmon: own(new THREE.MeshStandardMaterial({ map: ctex(salmonArt()), roughness: 0.4 })),
    gravy: std({ color: 0x6a4224, roughness: 0.22 }),
    paneSmoke: own(new THREE.MeshPhysicalMaterial({ color: 0x2a2b2a, metalness: 0, roughness: 0.02, transparent: true, opacity: 0.4, envMapIntensity: 1.7, clearcoat: 1, clearcoatRoughness: 0.02, depthWrite: false })),
    paneClear: own(new THREE.MeshPhysicalMaterial({ color: 0xd8e6e2, metalness: 0, roughness: 0.02, transparent: true, opacity: 0.1, envMapIntensity: 1.5, depthWrite: false })),
    paneEdge: own(new THREE.MeshStandardMaterial({ color: 0x4f7c6c, metalness: 0, roughness: 0.08, transparent: true, opacity: 0.72, depthWrite: false })),
    blackGlass: own(new THREE.MeshPhysicalMaterial({ color: 0x060708, metalness: 0, roughness: 0.08, clearcoat: 1, clearcoatRoughness: 0.03, envMapIntensity: 1.2 })),
    frit: new THREE.MeshStandardMaterial({ color: 0x0b0c0d, roughness: 0.35, metalness: 0 }),
    grille: own(new THREE.MeshStandardMaterial({ color: 0xa9adb2, metalness: 0.85, roughness: 0.45, alphaMap: grilleTex, alphaTest: 0.5, side: THREE.DoubleSide })),
  };
  return fam;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Generic parts                                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

/** A box by its extents [x0, x1] × [y0, y1] × [z0, z1], its edges rounded by `r`. */
export function slab(s: Shop, m: Mat, x0: number, x1: number, y0: number, y1: number, z0: number, z1: number, r = 0): void {
  const [w, h, d] = [Math.abs(x1 - x0), Math.abs(y1 - y0), Math.abs(z1 - z0)];
  s.add(m, r > 0 ? rbox(s, w, h, d, r) : box(w, h, d), [(x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2]);
}
/** A small part by its extents, rounded only close up, where the rounding can be seen. */
function bit(s: Shop, m: Mat, x0: number, x1: number, y0: number, y1: number, z0: number, z1: number, r = 0.001): void {
  slab(s, m, x0, x1, y0, y1, z0, z1, s.fine(2) ? r : 0);
}

/** k + 1 points per corner round a w × d rectangle centred on the origin in (x, z), its corners rounded by r. */
function rrRing(w: number, d: number, r: number, k: number): [number, number][] {
  const rr = Math.max(0, Math.min(r, w / 2, d / 2));
  const [cx, cz] = [w / 2 - rr, d / 2 - rr];
  const out: [number, number][] = [];
  for (const [x, z, a0] of [
    [cx, cz, 0],
    [-cx, cz, Math.PI / 2],
    [-cx, -cz, Math.PI],
    [cx, -cz, Math.PI * 1.5],
  ] as const)
    for (let i = 0; i <= k; i++) {
      const a = a0 + (i / k) * (Math.PI / 2);
      out.push([x + rr * Math.cos(a), z + rr * Math.sin(a)]);
    }
  return out;
}

/**
 * A drawn vessel (a hotel pan, a well, a tub): a w × d rectangle with corners rounded by r, swept through `profile`
 * [inset, y] from the outer bottom up and over the rim and down to the inner floor; each ring is the outline shrunk by
 * its inset (negative grows it, for a flange). Repeating a point makes a crisp fold. Both ends are capped flat.
 */
export function vessel(s: Shop, w: number, d: number, r: number, profile: readonly (readonly [number, number])[], caps: [boolean, boolean] = [true, true]): THREE.BufferGeometry {
  const k = s.fine(2) ? 5 : s.fine(1) ? 2 : 1;
  const rings = profile.map(([i, y]) => rrRing(w - 2 * i, d - 2 * i, r - i, k).map(([x, z]) => [x, y, z] as const));
  const N = rings[0]!.length;
  const pos: number[] = [];
  const idx: number[] = [];
  for (const ring of rings) for (const q of ring) pos.push(q[0], q[1], q[2]);
  for (let j = 0; j < rings.length - 1; j++)
    for (let i = 0; i < N; i++) {
      const a = j * N + i;
      const b = j * N + ((i + 1) % N);
      const c = (j + 1) * N + ((i + 1) % N);
      const dd = (j + 1) * N + i;
      idx.push(a, c, b, a, dd, c);
    }
  const cap = (ring: readonly (readonly [number, number, number])[], up: boolean) => {
    const o = pos.length / 3;
    pos.push(0, ring[0]![1], 0);
    for (const q of ring) pos.push(q[0], q[1], q[2]);
    for (let i = 0; i < N; i++) {
      const [p0, p1] = [o + 1 + i, o + 1 + ((i + 1) % N)];
      if (up) idx.push(o, p1, p0);
      else idx.push(o, p0, p1);
    }
  };
  if (caps[0]) cap(rings[0]!, false);
  if (caps[1]) cap(rings[rings.length - 1]!, true);
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/**
 * A food surface: a heightfield w × d (corners rounded by r) at y = 0, lumped by `amp` at a wavelength of about
 * `lambda` metres; the edge is pulled down a little where the food meets the pan wall.
 */
function foodTop(s: Shop, w: number, d: number, r: number, amp: number, lambda: number, seed: number): THREE.BufferGeometry {
  const n = s.fine(2) ? 1 : s.fine(1) ? 0.5 : 0.2;
  const nx = Math.max(3, Math.round((w / lambda) * 3 * n));
  const nz = Math.max(3, Math.round((d / lambda) * 3 * n));
  const g = new THREE.PlaneGeometry(w, d, nx, nz);
  g.rotateX(-Math.PI / 2);
  const p = g.attributes.position as THREE.BufferAttribute;
  const rr = Math.min(r, w / 2, d / 2);
  const [cx, cz] = [w / 2 - rr, d / 2 - rr];
  const k = (Math.PI * 2) / lambda;
  for (let i = 0; i < p.count; i++) {
    let [x, z] = [p.getX(i), p.getZ(i)];
    const [ax, az] = [Math.abs(x) - cx, Math.abs(z) - cz];
    if (ax > 0 && az > 0) {
      const l = Math.hypot(ax, az);
      if (l > rr) {
        x = Math.sign(x) * (cx + (ax / l) * rr);
        z = Math.sign(z) * (cz + (az / l) * rr);
      }
    }
    const edge = Math.min(w / 2 - Math.abs(x), d / 2 - Math.abs(z));
    const lump =
      Math.sin(x * k + seed) * Math.cos(z * k * 0.8 + seed * 1.7) +
      0.55 * Math.sin((x * 0.6 - z) * k * 1.9 + seed * 3.1) +
      0.3 * Math.cos((x + z * 0.7) * k * 3.3 + seed * 0.7);
    const y = amp * (0.5 + 0.32 * lump) - amp * 0.6 * Math.exp(-edge / (lambda * 0.25));
    p.setXYZ(i, x, y, z);
  }
  g.computeVertexNormals();
  return g;
}

const UP = new THREE.Vector3(0, 1, 0);

/**
 * A round tube along a polyline, bent as a fitter bends pipe: one straight run between bends, each bend a torus arc of
 * radius `bend`. Cheap however long its runs (a handle, a cable, a copper line).
 */
export function tube(s: Shop, m: Mat, pts: readonly V3[], r: number, bend: number, sides = 12, caps = true): void {
  const P = pts.map((q) => new THREE.Vector3(q[0], q[1], q[2]));
  const n = s.seg(sides, 5);
  let from = P[0]!.clone();
  for (let i = 1; i < P.length; i++) {
    const q = P[i]!;
    let to = q.clone();
    let next: THREE.Vector3 | null = null;
    if (i < P.length - 1) {
      const a = q.clone().sub(P[i - 1]!).normalize();
      const b = P[i + 1]!.clone().sub(q).normalize();
      const c = THREE.MathUtils.clamp(a.dot(b), -1, 1);
      const th = Math.acos(c);
      if (th > 1e-3) {
        const t = Math.min(bend * Math.tan(th / 2), q.distanceTo(P[i - 1]!) * 0.49, P[i + 1]!.distanceTo(q) * 0.49);
        const R = t / Math.tan(th / 2);
        const A = q.clone().addScaledVector(a, -t);
        const nrm = b.clone().addScaledVector(a, -c).normalize();
        const C = A.clone().addScaledVector(nrm, R);
        const ex = nrm.clone().negate();
        const ez = new THREE.Vector3().crossVectors(ex, a);
        const g = new THREE.TorusGeometry(R, r, n, Math.max(2, Math.ceil((th / (Math.PI / 2)) * s.seg(8, 3))), th);
        g.applyMatrix4(new THREE.Matrix4().makeBasis(ex, a, ez).setPosition(C));
        to = A;
        next = q.clone().addScaledVector(b, t);
        s.add(m, g);
      }
    }
    const len = from.distanceTo(to);
    if (len > 1e-5) {
      const g = new THREE.CylinderGeometry(r, r, len, n, 1, true);
      g.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(UP, to.clone().sub(from).normalize()));
      g.translate((from.x + to.x) / 2, (from.y + to.y) / 2, (from.z + to.z) / 2);
      s.add(m, g);
    }
    from = next ?? to;
  }
  if (caps)
    for (const [e, o] of [
      [P[0]!, P[1]!],
      [P[P.length - 1]!, P[P.length - 2]!],
    ] as const) {
      const g = new THREE.CircleGeometry(r, n);
      g.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 0, 1), e.clone().sub(o).normalize()));
      g.translate(e.x, e.y, e.z);
      s.add(m, g);
    }
}

/** A swivel caster under a plate at `at` (its underside), the wheel on y = 0: the kit's close up, a cheap one further off. */
function wheelSet(s: Shop, wheel: Mat, metal: Mat, at: V3, r: number): void {
  if (s.fine(2)) return caster(s, wheel, metal, at, r);
  const [x, y, z] = at;
  s.add(metal, box(0.05, 0.005, 0.05), [x, y - 0.0025, z]);
  s.add(metal, box(0.04, y - r, r * 1.1), [x, (y + r) / 2, z - r * 0.3]);
  s.add(wheel, cyl(s, r, 0.028, "x", 12), [x, r, z - r * 0.4]);
}

/** A thin open rod (no end caps: its ends always meet something), `sides` round, along `axis`. */
function rod(r: number, len: number, axis: "x" | "y" | "z", sides = 6): THREE.BufferGeometry {
  const g = new THREE.CylinderGeometry(r, r, len, sides, 1, true);
  if (axis === "x") g.rotateZ(-Math.PI / 2);
  else if (axis === "z") g.rotateX(Math.PI / 2);
  return g;
}

/**
 * A round-section loop round a rectangle in the xy plane facing +z (a door gasket, a bumper): four straight runs and
 * four quarter bends of radius `bend`, the tube `r` thick.
 */
function rectLoop(s: Shop, m: Mat, x0: number, x1: number, y0: number, y1: number, z: number, r: number, bend: number, sides = 10): void {
  const [cx0, cx1, cy0, cy1] = [x0 + bend, x1 - bend, y0 + bend, y1 - bend];
  s.add(m, rod(r, cx1 - cx0, "x", sides), [(cx0 + cx1) / 2, y0, z]);
  s.add(m, rod(r, cx1 - cx0, "x", sides), [(cx0 + cx1) / 2, y1, z]);
  s.add(m, rod(r, cy1 - cy0, "y", sides), [x0, (cy0 + cy1) / 2, z]);
  s.add(m, rod(r, cy1 - cy0, "y", sides), [x1, (cy0 + cy1) / 2, z]);
  const k = s.seg(10, 3);
  for (const [cx, cy, a] of [
    [cx1, cy1, 0],
    [cx0, cy1, Math.PI / 2],
    [cx0, cy0, Math.PI],
    [cx1, cy0, Math.PI * 1.5],
  ] as const)
    s.add(m, new THREE.TorusGeometry(bend, r, sides, k, Math.PI / 2), [cx, cy, z], [0, 0, a]);
}

/** A ring of `n` hex-ish bolts or rivets round a circle, facing +z. */
function boltCircle(s: Shop, m: Mat, at: V3, R: number, n: number, r = 0.003, phase = 0, rot?: V3): void {
  if (!s.fine(2)) return;
  for (let i = 0; i < n; i++) {
    const a = phase + (i / n) * Math.PI * 2;
    screw(s, m, [at[0] + Math.cos(a) * R, at[1] + Math.sin(a) * R, at[2]], r, rot);
  }
}

/**
 * A bank of pressed louvres in a sheet facing +z, centred at `at` in the frame turned by `rot`: n slots `w` wide,
 * `pitch` apart, each a dark opening under a lip pressed out of the sheet.
 */
function louvres(s: Shop, at: V3, rot: V3, n: number, w: number, pitch: number): void {
  const p = palette();
  s.at(at, rot, () => {
    for (let i = 0; i < n; i++) {
      const y = (i - (n - 1) / 2) * pitch;
      s.add(p.enamel, box(w, pitch * 0.42, 0.0012), [0, y - pitch * 0.06, 0.0009]);
      s.add(p.ss, s.fine(2) ? rbox(s, w + 0.004, pitch * 0.5, 0.0012, 0.0005) : box(w + 0.004, pitch * 0.5, 0.0012), [0, y + pitch * 0.02, 0.0034], [-0.55, 0, 0]);
      if (s.fine(2)) for (const sx of [-1, 1]) s.add(p.ss, box(0.0012, pitch * 0.42, 0.004), [sx * (w / 2 + 0.0014), y, 0.0024]);
    }
  });
}

/** A data plate facing +z: a brushed plate on four rivets, its engraved rating lines dark. No maker's name. */
function dataPlate(s: Shop, at: V3, rot: V3 = [0, 0, 0], w = 0.09, h = 0.055): void {
  if (!s.fine(1)) return;
  const p = palette();
  s.at(at, rot, () => {
    s.add(p.ssBright, rbox(s, w, h, 0.0008, 0.0003), [0, 0, 0.0004]);
    if (!s.fine(2)) return;
    for (const [sx, sy] of [
      [-1, -1],
      [1, -1],
      [-1, 1],
      [1, 1],
    ] as const)
      s.add(p.ss, lathe(s, [[0, 0], [0.0022, 0], [0.002, 0.0008], [0, 0.0011]], 10), [sx * (w / 2 - 0.005), sy * (h / 2 - 0.005), 0.0008], [Math.PI / 2, 0, 0]);
    const lines = [0.8, 0.55, 0.62, 0.48, 0.7, 0.4];
    for (const [i, f] of lines.entries()) s.add(p.enamel, box(w * 0.75 * f, h * 0.045, 0.0004), [-w * 0.36 + (w * 0.75 * f) / 2, h * 0.3 - i * h * 0.12, 0.001]);
    s.add(p.enamel, box(w * 0.2, h * 0.24, 0.0004), [w * 0.3, -h * 0.2, 0.001]);
  });
}

/** An adjustable foot: a stainless sleeve under the base at height `y1`, a threaded stem and a domed pad on the floor. */
function adjFoot(s: Shop, x: number, z: number, y1: number, r = 0.019): void {
  const p = palette();
  const sl = y1 * 0.5;
  s.add(p.ss, cyl(s, r, sl, "y", 20), [x, y1 - sl / 2, z]);
  s.add(p.ssBright, cyl(s, r * 0.42, y1 - sl - 0.008, "y", 12), [x, (y1 - sl + 0.008) / 2, z]);
  s.add(p.ssBright, lathe(s, [[0, 0], [r * 1.25, 0], [r * 1.3, 0.003], [r * 1.1, 0.008], [r * 0.5, 0.011], [0, 0.011]], 20), [x, 0, z]);
  if (s.fine(2)) s.add(p.ssBright, cyl(s, r * 0.75, 0.005, "y", 6), [x, y1 - sl - 0.004, z]);
}

/**
 * The steam a pan or a cavity gives off: torn, noise-textured wisps (view-aligned sprites, so they read from every
 * side) drifting in a box, their strength live. Each is a plume taller than it is wide; off a surface they rise, swell
 * and thin out as they go, swaying; in a cavity the fan churns them round, turning slowly as they are carried.
 */
class Haze {
  readonly group = new THREE.Group();
  private items: { s: THREE.Sprite; m: THREE.SpriteMaterial; x: number; y: number; z: number; ph: number; sp: number; w: number; h: number; rot: number; spin: number }[] = [];
  private level = 0;
  constructor(
    private bx: { x0: number; x1: number; y0: number; y1: number; z0: number; z1: number },
    n: number,
    size: number,
    seed: number,
    private rise = 0,
    tint = 0xffffff,
    private peak = 0.5,
  ) {
    const r = rng(seed);
    const tex = wisps();
    for (let i = 0; i < n; i++) {
      const rot = rise > 0 ? (r() - 0.5) * 0.5 : r() * Math.PI * 2;
      const m = new THREE.SpriteMaterial({ map: tex[i % tex.length], color: tint, transparent: true, depthWrite: false, opacity: 0, rotation: rot });
      m.userData.uv = "own";
      const sp = new THREE.Sprite(m);
      sp.renderOrder = 1;
      const w = size * (0.45 + r() * 0.4);
      const h = w * (1.7 + r() * 0.7);
      sp.scale.set(w, h, 1);
      this.group.add(sp);
      this.items.push({ s: sp, m, x: r(), y: r(), z: r(), ph: r() * 6.28, sp: 0.6 + r() * 0.8, w, h, rot, spin: (r() - 0.5) * 0.5 });
    }
    this.group.visible = false;
  }
  get tris(): number {
    return this.items.length * 2;
  }
  set(v: number): void {
    this.level = Math.max(0, Math.min(1, v));
    for (const it of this.items) it.m.opacity = this.peak * this.level;
    this.group.visible = this.level > 0.01;
  }
  tick(t: number): void {
    if (!this.group.visible) return;
    const b = this.bx;
    for (const it of this.items) {
      if (this.rise > 0) {
        // Rising off a surface: born small and dense, swelling and thinning as it climbs, leaning as it sways.
        const u = (it.y + t * this.rise * it.sp) % 1;
        const y = b.y0 + u * (b.y1 - b.y0);
        const k = Math.sin(Math.PI * Math.min(1, u * 1.25)) * (1 - u * 0.35);
        const g = 0.55 + 0.9 * u;
        it.s.scale.set(it.w * g, it.h * g, 1);
        it.m.opacity = this.peak * this.level * Math.max(0, k);
        it.m.rotation = it.rot + 0.25 * Math.sin(t * 0.7 + it.ph) * u;
        it.s.position.set(b.x0 + (b.x1 - b.x0) * (it.x + 0.06 * Math.sin(t * 0.9 + it.ph) * u), y + (it.h * g) * 0.3, b.z0 + (b.z1 - b.z0) * it.z);
      } else {
        // Churned round by the fan: a slow loop in the cavity, each wisp turning as it is carried.
        const a = t * 0.35 * it.sp + it.ph;
        const fx = 0.5 + (it.x - 0.5) * 0.7 + 0.14 * Math.cos(a);
        const fy = 0.5 + (it.y - 0.5) * 0.7 + 0.14 * Math.sin(a);
        it.m.rotation = it.rot + t * it.spin;
        const br = 0.75 + 0.25 * Math.sin(t * 0.6 * it.sp + it.ph * 2);
        it.m.opacity = this.peak * this.level * br;
        it.s.position.set(b.x0 + (b.x1 - b.x0) * fx, b.y0 + (b.y1 - b.y0) * fy, b.z0 + (b.z1 - b.z0) * (0.2 + 0.6 * it.z));
      }
    }
  }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Combi-steamers                                                                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

interface Box3 {
  x0: number;
  x1: number;
  y0: number;
  y1: number;
  z0: number;
  z1: number;
}

interface CombiLayout {
  pans: 6 | 10 | 20;
  w: number;
  d: number;
  h: number;
  /** Body underside and top. */
  base: number;
  top: number;
  /** Door front and inner faces, the cavity frame's face behind the door, the back panel. */
  zDoorFront: number;
  zDoorBack: number;
  zFrame: number;
  zBack: number;
  door: { x0: number; x1: number; y0: number; y1: number };
  /** The control column, right of the door, and its face. */
  panel: { x0: number; x1: number; z: number };
  /** The working cavity: z0 the baffle's face, z1 the opening; `wall` the back wall behind the fan. */
  cav: Box3 & { wall: number };
  cx: number;
  rails: number[];
  fans: number[];
  fanR: number;
  fanZ: number;
  /** Top of the stand the oven sits on (0: none). */
  stand: number;
}

function combiLayout(pans: 6 | 10 | 20, e: Ext): CombiLayout {
  const { w, d, h } = e;
  const big = pans === 20;
  const stand = pans === 10 ? Math.max(0.3, h - 1.062) : 0;
  const base = pans === 10 ? stand + 0.05 : big ? 0.008 : 0.06;
  const top = h - 0.046;
  const zDoorFront = -0.038;
  const zDoorBack = zDoorFront - (big ? 0.055 : 0.05);
  const zFrame = zDoorBack - 0.006;
  const zBack = -d + 0.024;
  const doorX1 = w - (big ? 0.215 : 0.209);
  const door = { x0: 0.006, x1: doorX1, y0: big ? base + 0.004 : base + 0.032, y1: top - 0.014 };
  const cx = (door.x0 + door.x1) / 2;
  const cw = big ? 0.49 : 0.46;
  const wall = zBack + (big ? 0.032 : 0.03);
  const z0 = wall + (big ? 0.07 : 0.058);
  const y0 = big ? 0.035 : base + 0.085;
  const y1 = top - 0.1;
  const cav = { x0: cx - cw / 2, x1: cx + cw / 2, y0, y1, z0, z1: zFrame - 0.0015, wall };
  let rails: number[];
  if (big) rails = Array.from({ length: 20 }, (_, i) => y0 + 0.19 + i * 0.068);
  else {
    const n = pans;
    const topRail = y1 - 0.078;
    rails = Array.from({ length: n }, (_, i) => topRail - (n - 1 - i) * 0.068);
  }
  const H = y1 - y0;
  const fans = pans === 6 ? [(y0 + y1) / 2] : pans === 10 ? [y0 + 0.25 * H, y0 + 0.75 * H] : [y0 + 0.28 * H, y0 + 0.72 * H];
  const fanR = pans === 6 ? 0.15 : pans === 10 ? 0.145 : 0.18;
  return {
    pans,
    w,
    d,
    h,
    base,
    top,
    zDoorFront,
    zDoorBack,
    zFrame,
    zBack,
    door,
    panel: { x0: doorX1 + 0.006, x1: w - 0.005, z: zDoorFront - 0.004 },
    cav,
    cx,
    rails,
    fans,
    fanR,
    fanZ: wall + 0.024,
    stand,
  };
}

/** Half the spacing of the GN runner rods a pan's rim rests on, and of the frame rods that carry them. */
const RUNNER = 0.155;
const FRAME = 0.205;
/** How far in from the door's left edge its hinge pins stand (in the plane of its outer face). */
const HINGE_IN = 0.004;

export interface CombiCavity {
  /** The working cavity's clear box (fixture frame): wall to wall in x, floor to ceiling in y, baffle face (z0) to the opening (z1). */
  box: Box3;
  /** Rail heights, bottom level first: the y where a pan's rim rests on its runners. */
  rails: number[];
  /** The x of the shelf centre line, and of the left and right runner rods under a pan's rim. */
  cx: number;
  runners: [number, number];
  /** The z of a pan's centre when pushed home on its runners. */
  zc: number;
  /**
   * The door (fixture frame): its hinge pins' line at x, z (in the plane of its outer face, hinged left); the leaf's
   * length from the pins to its free edge; how far its inner face lies behind the pins' line; the pull's grip along
   * the leaf from the pins and out from its face, and where along the leaf the pull begins to stand out.
   */
  door: { x: number; z: number; len: number; back: number; hAlong: number; hOut: number; pullFrom: number };
}

/** Where pans go in a combi oven: its cavity box and rail heights in the fixture frame, or null for another archetype. */
export function combiCavity(id: string, e: Ext): CombiCavity | null {
  const pans = id === "NSF-COMBI-6PAN" ? 6 : id === "NSF-COMBI-10PAN" ? 10 : id === "NSF-COMBI-20PAN-ROLLIN" ? 20 : 0;
  if (!pans) return null;
  const L = combiLayout(pans, e);
  const { x0, x1, y0, y1, z0, z1 } = L.cav;
  // The pull (see the door in combi()): a grip on two standoffs 30 mm in from the free edge, its axis 26 mm out.
  const len = L.door.x1 - L.door.x0 - HINGE_IN;
  const hAlong = len - 0.03;
  const door = { x: L.door.x0 + HINGE_IN, z: L.zDoorFront, len, back: L.zDoorFront - L.zDoorBack, hAlong, hOut: 0.026, pullFrom: hAlong - 0.022 };
  return { box: { x0, x1, y0, y1, z0, z1 }, rails: [...L.rails], cx: L.cx, runners: [L.cx - RUNNER, L.cx + RUNNER], zc: z0 + (pans === 20 ? 0.042 : 0.024) + 0.265, door };
}

/** A backward-curved fan blade in the wheel's plane: from radius r0 to r1, swept `sweep` radians, `t` thick, `depth` deep along +z. */
function fanBlade(s: Shop, r0: number, r1: number, sweep: number, t: number, depth: number): THREE.BufferGeometry {
  const n = s.fine(2) ? 8 : 3;
  const a: THREE.Vector2[] = [];
  const b: THREE.Vector2[] = [];
  for (let i = 0; i <= n; i++) {
    const u = i / n;
    const r = r0 + (r1 - r0) * u;
    const th = sweep * Math.sqrt(u);
    a.push(new THREE.Vector2(r * Math.cos(th), r * Math.sin(th)));
    const th2 = th + t / r;
    b.push(new THREE.Vector2(r * Math.cos(th2), r * Math.sin(th2)));
  }
  const shape = new THREE.Shape([...a, ...b.reverse()]);
  return new THREE.ExtrudeGeometry(shape, { depth, bevelEnabled: false, curveSegments: 1 });
}

/**
 * A combi-steamer of `pans` GN 1/1 levels: countertop (6), on its stand (10), or the floor model with a roll-in rack
 * trolley (20).
 */
function combi(pans: 6 | 10 | 20): Builder {
  return (s, e) => {
    const p = palette();
    const fm = mats();
    const L = combiLayout(pans, e);
    const { w, base, top, zFrame, zBack, cav, cx } = L;
    const big = pans === 20;
    const glow = glowMaterial(0x2b2724);
    const fine1 = s.fine(1);
    const fine2 = s.fine(2);
    /** The radius of the cavity liner's coves, and of the front frame's opening that matches them. */
    const rc = big ? 0.03 : 0.025;

    /* Feet, or the stand. */
    const footX = [0.05, w - 0.05];
    const footZ = [zFrame - 0.05, zBack + 0.05];
    if (pans === 10) {
      const st = L.stand;
      const sx = [0.012, w - 0.012];
      const sz = [zFrame - 0.012, zBack + 0.012];
      const leg = 0.04;
      for (const x of sx)
        for (const z of sz) {
          const lx = x + (x < w / 2 ? leg / 2 : -leg / 2);
          const lz = z + (z < -0.3 ? leg / 2 : -leg / 2);
          bit(s, p.ss, lx - leg / 2, lx + leg / 2, 0.075, st - 0.001, lz - leg / 2, lz + leg / 2, 0.006);
          s.add(p.ssBright, lathe(s, [[0, 0], [0.024, 0], [0.025, 0.006], [0.02, 0.03], [0.015, 0.075], [0, 0.075]], 20), [lx, 0, lz]);
          if (fine2) s.add(p.enamel, box(leg - 0.008, 0.002, leg - 0.008), [lx, st + 0.0005, lz]);
        }
      // Top frame: square tube round the perimeter, the oven's feet standing on it.
      bit(s, p.ss, sx[0]!, sx[1]!, st - 0.03, st, sz[0]! - 0.03, sz[0]!, 0.004);
      bit(s, p.ss, sx[0]!, sx[1]!, st - 0.03, st, sz[1]!, sz[1]! + 0.03, 0.004);
      for (const x of sx) bit(s, p.ss, x < w / 2 ? x : x - 0.03, x < w / 2 ? x + 0.03 : x, st - 0.03, st, sz[1]! + 0.03, sz[0]! - 0.03, 0.004);
      // The undershelf: a sheet with its edges turned down, and a lower stretcher.
      bit(s, p.ss, sx[0]! + 0.04, sx[1]! - 0.04, 0.16, 0.1615, sz[1]! + 0.04, sz[0]! - 0.04, 0.0006);
      for (const z of [sz[0]! - 0.04, sz[1]! + 0.04]) bit(s, p.ss, sx[0]! + 0.04, sx[1]! - 0.04, 0.13, 0.1615, z - 0.0008, z + 0.0008, 0.0005);
      for (const x of [sx[0]! + 0.04, sx[1]! - 0.04]) bit(s, p.ss, x - 0.0008, x + 0.0008, 0.13, 0.1615, sz[1]! + 0.04, sz[0]! - 0.04, 0.0005);
      // Pan slides: two side frames of angle runners, a GN 1/1 lengthwise between them.
      const mid = w / 2;
      const levels: number[] = [];
      for (let y = st - 0.075; y > 0.24; y -= 0.065) levels.push(y);
      for (const sg of [-1, 1]) {
        const xf = mid + sg * (RUNNER + 0.022);
        for (const z of [sz[0]! - 0.05, sz[1]! + 0.05]) bit(s, p.ss, xf - 0.0015, xf + 0.0015, 0.1615, st - 0.03, z - 0.015, z + 0.015, 0.0006);
        for (const y of levels) {
          slab(s, p.ss, Math.min(xf, xf - sg * 0.024), Math.max(xf, xf - sg * 0.024), y - 0.0015, y, sz[1]! + 0.035, sz[0]! - 0.035, fine2 ? 0.0005 : 0);
          if (fine1) slab(s, p.ss, xf - 0.0015 + sg * 0.0015, xf + 0.0015 + sg * 0.0015, y - 0.016, y, sz[1]! + 0.035, sz[0]! - 0.035, fine2 ? 0.0005 : 0);
        }
      }
      // Two empty GN 1/1 pans stored on the slides.
      if (fine1)
        for (const [k, y] of [levels[1], levels[4]].entries()) {
          if (y === undefined) continue;
          const pd = k ? 0.065 : 0.1;
          s.add(p.ss, gnPan(s, 0.325, 0.53, pd), [mid, y + 0.0005, (sz[0]! + sz[1]!) / 2]);
        }
    }
    if (pans === 10) {
      const fx = [0.027, w - 0.027];
      const fz = [zFrame - 0.027, zBack + 0.027];
      s.at([0, L.stand, 0], [0, 0, 0], () => {
        for (const x of fx) for (const z of fz) adjFoot(s, x, z, base - L.stand, 0.017);
      });
    } else if (big)
      for (const x of footX) for (const z of footZ) s.add(p.ssBright, cyl(s, 0.03, base, "y", 20), [x, base / 2, z]);
    else for (const x of footX) for (const z of footZ) adjFoot(s, x, z, base, 0.019);

    /* The housing: a skin of 1.2 mm sheet over the insulated cavity, its front frame round the opening. */
    const sk = 0.0012;
    const xR = L.panel.x1;
    slab(s, p.ss, 0, sk, base, top - sk, zBack, zFrame - 0.0015, 0.0005);
    slab(s, p.ss, 0, xR, top - sk, top, zBack, L.zDoorFront, 0.0005);
    // The top's front return, stopped short of the upper hinge bracket at the left.
    slab(s, p.ss, 0.033, xR, top - 0.012, top - sk, L.zDoorFront - sk, L.zDoorFront, 0.0005);
    slab(s, p.ssDull, 0.001, xR - 0.001, base, base + sk, zBack, zFrame - 0.0015);
    slab(s, p.ssDull, 0, xR, base, top - sk, zBack - sk, zBack);
    // Folded side returns at the back corners, and the skin's screws down the rear edge.
    if (fine2)
      for (const x of [0, w])
        for (let y = base + 0.06; y < top - 0.03; y += (top - base - 0.1) / 4) screw(s, p.ssBright, [x < w / 2 ? 0 : xR, y, zBack + 0.012], 0.0026, [0, x < w / 2 ? -Math.PI / 2 : Math.PI / 2, 0]);
    {
      const fw = L.panel.x0;
      const fh = top - sk - base;
      const sh = roundRect(fw, fh, 0.002, fw / 2, base + fh / 2);
      rectPath(sh, cav.x1 - cav.x0, cav.y1 - cav.y0, rc, cx, (cav.y0 + cav.y1) / 2, true);
      s.add(p.ss, plate(s, sh, 0.0015), [0, 0, zFrame - 0.0015]);
    }
    // The gasket: a hollow EPDM profile on the cavity frame, pressed by the door's inner pane.
    {
      const g = 0.014;
      const [gx0, gx1, gy0, gy1] = [cav.x0 - g, cav.x1 + g, cav.y0 - g, cav.y1 + g];
      const z = zFrame + 0.0035;
      if (fine1) rectLoop(s, p.gasket, gx0, gx1, gy0, gy1, z, 0.0042, 0.02, fine2 ? 10 : 6);
    }
    // Door drip channel: a formed trough under the door that leads its condensate into the cavity drain. The roll-in's
    // door reaches the floor, so its channel is a gutter across the cavity floor just inside the opening.
    if (big) {
      if (fine1) {
        const z0g = cav.z1 - 0.035;
        const z1g = cav.z1 - 0.008;
        slab(s, p.enamel, cav.x0, cav.x1, cav.y0 - 0.012, cav.y0 + 0.0003, z0g, z1g);
        const n = fine2 ? 40 : 0;
        for (let i = 0; i < n; i++) {
          const x = cav.x0 + 0.01 + ((cav.x1 - cav.x0 - 0.02) * i) / (n - 1);
          s.add(p.ss, box(0.004, 0.004, z1g - z0g), [x, cav.y0 - 0.0018, (z0g + z1g) / 2]);
        }
        if (!fine2) s.add(fm.perf, new THREE.PlaneGeometry(cav.x1 - cav.x0, z1g - z0g).rotateX(-Math.PI / 2), [cx, cav.y0 - 0.0015, (z0g + z1g) / 2]);
        slab(s, p.ss, cav.x0, cav.x1, cav.y0 - 0.006, cav.y0 - 0.0045, z0g, z1g);
      }
    } else {
      // A deep-drawn trough under the door, its rim rolled, falling to an outlet at its right end that drains into the
      // cavity's drain; it stands a little proud of the shut door so that the door's condensate falls into it when it
      // opens.
      const y1 = L.door.y0 - 0.004;
      const dep = Math.min(0.02, y1 - base - 0.006);
      const [x0, x1] = [L.door.x0 + 0.03, L.door.x1];
      const zf = L.zDoorFront + 0.006;
      const tw = x1 - x0;
      const td = zf - zFrame;
      const trough = fine2
        ? vessel(s, tw, td, 0.008, [
            [0.006, -dep],
            [0.0015, -dep + 0.004],
            [0, -0.004],
            [-0.0012, -0.0012],
            [0, 0],
            [0.0016, -0.0008],
            [0.0026, -0.004],
            [0.0074, -dep + 0.0012],
          ])
        : vessel(s, tw, td, 0.006, [
            [0.006, -dep],
            [0, -0.002],
            [0, 0],
            [0.0022, -0.002],
            [0.0074, -dep + 0.0012],
          ]);
      s.add(p.ss, trough, [(x0 + x1) / 2, y1, (zFrame + zf) / 2]);
      if (fine2) {
        // The outlet: a dark hole in a pressed dimple, its pipe stub down through the base.
        const zc = (zFrame + zf) / 2;
        s.add(p.enamel, cyl(s, 0.0055, 0.0008, "y", 16), [x1 - 0.022, y1 - dep + 0.0016, zc]);
        s.add(p.ss, torus(s, 0.0062, 0.0011, "y", Math.PI * 2, 18, 6), [x1 - 0.022, y1 - dep + 0.0016, zc]);
        s.add(p.ss, cyl(s, 0.006, y1 - dep - base - 0.001, "y", 14, 0.006, true), [x1 - 0.022, (y1 - dep + base + 0.001) / 2, zc]);
      }
    }
    // Hinges. The door turns on two pins standing in the plane of its outer face, just in from its left edge, so that it
    // swings forward and its leaf never passes behind the front into the space beside the oven (a griddle or a
    // toaster may stand flush there). Each pin is carried by a cast bracket fixed to the cavity frame, one over the
    // door and one under it, ending in a knuckle round the pin; the door's own knuckles, on straps along its left
    // edge, turn on the same pins between them, a bronze washer in each joint.
    const kr = 0.0075;
    const pivot: V3 = [L.door.x0 + HINGE_IN, 0, L.zDoorFront];
    {
      const [px, pz] = [pivot[0], pivot[2]];
      const ns = s.seg(20, 10);
      const yLo = Math.max(0.0015, L.door.y0 - 0.0105);
      for (const [b0, b1] of [[L.door.y1 + 0.0015, top - sk - 0.0003], [yLo, L.door.y0 - 0.0015]] as const) {
        bit(s, p.ss, 0.0006, 0.031, b0 + 0.0004, b1 - 0.0004, zFrame, pz - 0.0025, 0.0012);
        s.add(p.ssBright, cyl(s, kr, b1 - b0, "y", ns), [px, (b0 + b1) / 2, pz]);
        if (fine2) {
          // The bracket's two screws into the frame, on its front.
          screw(s, p.ssBright, [0.024, (b0 + b1) / 2, pz - 0.0025], 0.0019);
          // The washer between this knuckle and the door's.
          const yw = b0 > L.door.y1 ? L.door.y1 + 0.00075 : L.door.y0 - 0.00075;
          s.add(p.brass, cyl(s, kr * 0.92, 0.0011, "y", ns), [px, yw, pz]);
        }
      }
      // The lower pin's head under its knuckle, where there is room above the floor.
      if (fine1 && yLo > 0.004) s.add(p.ssBright, lathe(s, [[0, 0], [kr * 0.5, 0], [kr * 0.62, 0.0008], [kr * 0.68, 0.0018], [kr * 0.68, 0.003], [0, 0.003]], s.seg(16, 8)), [px, yLo - 0.0025, pz]);
    }

    /* The cavity: deep-drawn stainless, its runners, the fan and element behind the air baffle. */
    const ct = 0.0015;
    // The liner is drawn from one sheet, so its walls meet the ceiling and the floor in generous coves. The left half
    // and the ceiling stay; the right wall and its two coves belong to the cutaway.
    const nArc = s.seg(8, 3);
    const linerL: SecPt[] = [
      ...covePts(cav.x0 + rc, cav.y0 + rc, rc, -Math.PI / 2, -Math.PI, nArc),
      ...covePts(cav.x0 + rc, cav.y1 - rc, rc, Math.PI, Math.PI / 2, nArc),
      [cav.x1 - rc, cav.y1, 0, -1],
    ];
    const linerR: SecPt[] = [...covePts(cav.x1 - rc, cav.y1 - rc, rc, Math.PI / 2, 0, nArc), ...covePts(cav.x1 - rc, cav.y0 + rc, rc, 0, -Math.PI / 2, nArc)];
    if (fine1) {
      s.add(p.ss, sweepZ(linerL, cav.wall, cav.z1, ct));
      {
        const bw = roundRect(cav.x1 - cav.x0, cav.y1 - cav.y0, rc, cx, (cav.y0 + cav.y1) / 2);
        s.add(p.ssHeat, plate(s, bw, ct, 0, 8), [0, 0, cav.wall - ct]);
      }
      // The floor, its drain at the back under the baffle's foot.
      const dz = cav.z0 + 0.04;
      const fz1 = big ? cav.z1 - 0.035 : cav.z1;
      const fl = roundRect(cav.x1 - cav.x0 - 2 * rc, fz1 - cav.wall, 0, cx, (fz1 + cav.wall) / 2);
      fl.holes.push(circleHole(cx, dz, 0.024));
      const flg = plate(s, fl, ct);
      flg.rotateX(Math.PI / 2);
      s.add(p.ss, flg, [0, cav.y0, 0]);
      if (big) slab(s, p.ss, cav.x0 + rc, cav.x1 - rc, cav.y0 - ct, cav.y0, cav.z1 - 0.008, cav.z1);
      s.add(p.enamel, lathe(s, [[0.024, 0], [0.023, -0.012], [0.016, -0.03], [0, -0.032]], 20), [cx, cav.y0 - ct, dz]);
      // The drain's strainer: a slotted stainless cap.
      if (fine2) {
        const sh = new THREE.Shape().absarc(0, 0, 0.028, 0, Math.PI * 2, false);
        for (let i = 0; i < 6; i++) {
          const a = (i / 6) * Math.PI * 2;
          const hole = new THREE.Path();
          const [c, sn] = [Math.cos(a), Math.sin(a)];
          const pts = [
            [0.006, -0.0018],
            [0.021, -0.0018],
            [0.021, 0.0018],
            [0.006, 0.0018],
          ].map(([u, v]) => new THREE.Vector2(u! * c - v! * sn, u! * sn + v! * c));
          hole.setFromPoints(pts.reverse());
          sh.holes.push(hole);
        }
        const g = plate(s, sh, 0.002, 0.0005);
        g.rotateX(-Math.PI / 2);
        s.add(p.ssBright, g, [cx, cav.y0 + 0.0006, dz]);
      } else s.add(p.ssBright, cyl(s, 0.028, 0.002, "y", 16), [cx, cav.y0 + 0.0016, dz]);
      // The trolley's floor guides (roll-in): two angles that steer its wheels home.
      if (big)
        for (const sg of [-1, 1]) {
          const x = cx + sg * (FRAME + 0.028);
          slab(s, p.ss, x - 0.0015, x + 0.0015, cav.y0, cav.y0 + 0.03, cav.z0 + 0.02, cav.z1 - 0.01, 0.0005);
          slab(s, p.ss, Math.min(x, x + sg * 0.025), Math.max(x, x + sg * 0.025), cav.y0, cav.y0 + 0.0025, cav.z0 + 0.02, cav.z1 - 0.01, 0.0005);
        }
    }

    if (!fine1) slab(s, p.ssHeat, cav.x0, cav.x1, cav.y0, cav.y1, cav.wall, cav.z1 - 0.002);

    // Mineral wool packed between the cavity and the skin: only the cutaway shows it.
    if (fine1) {
      const g = 0.004;
      slab(s, fm.wool, sk + g, cav.x0 - ct - g, base + sk + g, top - sk - g, cav.wall - ct, zFrame - 0.0015 - g);
      slab(s, fm.wool, cav.x0 - ct - g, cav.x1 + 0.03, cav.y1 + ct + g, top - sk - g, cav.wall - ct, zFrame - 0.0015 - g);
      if (!big) slab(s, fm.wool, cav.x0 - ct - g, cav.x1 + 0.03, base + sk + g, cav.y0 - 0.034, cav.wall - ct, zFrame - 0.0015 - g);
    }

    // The lamp: a halogen behind toughened glass in a sealed bezel on the left wall, one per 300 mm of height.
    if (fine1) {
      const H = cav.y1 - cav.y0;
      const nl = H > 1.2 ? 3 : H > 0.65 ? 2 : 1;
      for (let i = 0; i < nl; i++) {
        const y = nl > 1 ? cav.y1 - 0.075 - (i * (cav.y1 - cav.y0 - 0.15)) / (nl - 1) : cav.y1 - 0.075;
        const z = cav.z1 - 0.13;
        const x = cav.x0;
        slab(s, p.ss, x, x + 0.005, y - 0.032, y + 0.032, z - 0.05, z + 0.05, 0.002);
        slab(s, fm.lamp, x + 0.0053, x + 0.0058, y - 0.02, y + 0.02, z - 0.035, z + 0.035);
        slab(s, p.glass, x + 0.0062, x + 0.0074, y - 0.026, y + 0.026, z - 0.043, z + 0.043, 0.0006);
        if (fine2) {
          for (const sy of [-1, 1]) slab(s, p.ssBright, x + 0.005, x + 0.0085, y + sy * 0.029 - 0.003, y + sy * 0.029 + 0.003, z - 0.05, z + 0.05, 0.001);
          for (const sz of [-1, 1]) slab(s, p.ssBright, x + 0.005, x + 0.0085, y - 0.026, y + 0.026, z + sz * 0.047 - 0.003, z + sz * 0.047 + 0.003, 0.001);
          for (const sg of [-1, 1]) screw(s, p.ssBright, [x + 0.0085, y + sg * 0.029, z + sg * 0.047], 0.0022, [0, Math.PI / 2, 0]);
        }
      }
    }

    // Hanging rack frames (6 and 10): rod frames hooked on wall pins, a runner rod for each level.
    const rodR = 0.003;
    const runnerSet = (sh: Shop, xc: number, zf: number, zb: number, ys: number[], y0: number, y1: number, sides: number[]) => {
      for (const sg of sides) {
        const xf = xc + sg * FRAME;
        const xr = xc + sg * RUNNER;
        for (const z of [zf, zb]) {
          if (fine2) sh.add(p.ss, rod(0.004, y1 - y0, "y", 10), [xf, (y0 + y1) / 2, z]);
          else sh.add(p.ss, box(0.007, y1 - y0, 0.007), [xf, (y0 + y1) / 2, z]);
        }
        for (const y of ys) {
          const yc = y - rodR;
          if (fine2) {
            sh.add(p.ss, rod(rodR, zf - zb, "z", 7), [xr, yc, (zf + zb) / 2]);
            for (const z of [zf, zb]) sh.add(p.ss, rod(rodR, Math.abs(xf - xr) + rodR, "x", 6), [(xf + xr) / 2, yc, z]);
            // A stop bent up at the back end, so that a pan cannot be pushed into the baffle.
            sh.add(p.ss, rod(rodR, 0.018, "y", 6), [xr, yc + 0.009, zb]);
          } else {
            sh.add(p.ss, box(rodR * 2, rodR * 2, zf - zb), [xr, yc, (zf + zb) / 2]);
            if (fine1) sh.add(p.ss, box(Math.abs(xf - xr), rodR * 2, rodR * 2), [(xf + xr) / 2, yc, zf]);
          }
        }
        // Top and bottom ties.
        for (const y of [y0 + 0.004, y1 - 0.004]) sh.add(p.ss, fine2 ? rod(0.004, zf - zb, "z", 10) : box(0.007, 0.007, zf - zb), [xf, y, (zf + zb) / 2]);
      }
    };
    if (fine1 && !big) {
      const y0 = cav.y0 + 0.02;
      const y1 = cav.y1 - 0.03;
      runnerSet(s, cx, cav.z1 - 0.03, cav.z0 + 0.022, L.rails, y0, y1, [-1, 1]);
      // Hooks over pins on the walls.
      if (fine2)
        for (const sg of [-1, 1])
          for (const z of [cav.z1 - 0.03, cav.z0 + 0.022]) {
            const xf = cx + sg * FRAME;
            const xw = sg < 0 ? cav.x0 : cav.x1;
            tube(s, p.ss, [[xf, y1, z], [xf, y1 + 0.014, z], [xw - sg * 0.008, y1 + 0.014, z], [xw - sg * 0.008, y1 - 0.004, z]], 0.0035, 0.005, 6);
            s.add(p.ssBright, cyl(s, 0.005, 0.012, "x", 10), [xw - sg * 0.006, y1 + 0.006, z]);
          }
    }

    // The air baffle: a stainless plate across the back wall, the fan's intake a perforated eye in it; air returns to
    // the cavity through the slots down its sides. It swings open on its left edge for cleaning.
    if (fine1) {
      const bx0 = cav.x0 + 0.03;
      const bx1 = cav.x1 - 0.03;
      const by0 = cav.y0 + 0.028;
      const by1 = cav.y1 - 0.028;
      const sh = roundRect(bx1 - bx0, by1 - by0, 0.006, (bx0 + bx1) / 2, (by0 + by1) / 2);
      const eye = L.fanR * 0.62;
      for (const fy of L.fans) sh.holes.push(circleHole(cx, fy, eye));
      if (fine2) {
        const ns = Math.floor((by1 - by0 - 0.06) / 0.055);
        for (const x of [bx0 + 0.016, bx1 - 0.016])
          for (let i = 0; i < ns; i++) {
            const y = by0 + 0.04 + i * 0.055;
            const hole = new THREE.Path();
            rectPath(hole, 0.008, 0.036, 0, x, y + 0.018, true);
            sh.holes.push(hole);
          }
      }
      s.add(p.ss, plate(s, sh, 0.0015, 0, 24), [0, 0, cav.z0 - 0.0015]);
      for (const fy of L.fans) {
        s.add(fm.perf, new THREE.CircleGeometry(eye + 0.003, s.seg(40, 12)), [cx, fy, cav.z0 - 0.0008]);
        // A pressed rim round the eye.
        if (fine2) s.add(p.ss, torus(s, eye + 0.002, 0.0018, "z", Math.PI * 2, 40, 4), [cx, fy, cav.z0 - 0.0005]);
      }
      // Hinge knuckles on the left edge and quarter-turn fasteners on the right.
      if (fine2) {
        for (const y of [by0 + 0.05, by1 - 0.05]) s.add(p.ssBright, cyl(s, 0.004, 0.03, "y", 10), [bx0 - 0.004, y, cav.z0 - 0.004]);
        for (const y of [by0 + 0.04, by1 - 0.04]) {
          s.add(p.ssBright, cyl(s, 0.007, 0.003, "z", 14), [bx1 - 0.012, y, cav.z0 + 0.0015]);
          s.add(p.ss, box(0.009, 0.0016, 0.0016), [bx1 - 0.012, y, cav.z0 + 0.0032]);
        }
      }
      // Brackets holding the baffle off the wall.
      if (fine2) for (const y of [by0 + 0.02, by1 - 0.02]) for (const x of [bx0 + 0.01, bx1 - 0.01]) slab(s, p.ss, x - 0.006, x + 0.006, y - 0.006, y + 0.006, cav.wall, cav.z0 - 0.0015);
    }

    // Behind it: the radial fan wheel(s), the ring element round each, the shaft seal, and the steam nozzle that
    // sprays the hub.
    const fanJoints: string[] = [];
    if (fine1)
      for (const [fi, fy] of L.fans.entries()) {
        const R = L.fanR;
        const name = fi ? `fan_${fi + 1}` : "fan";
        fanJoints.push(name);
        const f = s.joint(name, [cx, fy, L.fanZ]);
        const dep = 0.036;
        // Back plate, blades, and the front shroud ring whose open eye faces the baffle.
        f.add(p.ssDull, cyl(f, R, 0.002, "z", 48), [0, 0, -dep / 2 - 0.001]);
        const nb = fine2 ? 14 : 10;
        for (let i = 0; i < nb; i++) f.add(p.ssDull, fanBlade(f, R * 0.42, R * 0.98, 0.55, 0.0016, dep), [0, 0, -dep / 2], [0, 0, (i / nb) * Math.PI * 2]);
        f.add(p.ssDull, lathe(f, [[R * 0.6, 0], [R, 0], [R, 0.002], [R * 0.6, 0.004]], 48), [0, 0, dep / 2], [Math.PI / 2, 0, 0]);
        // Hub, the nut that holds it, and the slinger.
        f.add(p.ssBright, lathe(f, [[0, -0.004], [0.034, -0.004], [0.034, -0.002], [0.026, 0.004], [0.012, 0.014], [0, 0.016]], 24), [0, 0, -dep / 2], [Math.PI / 2, 0, 0]);
        if (fine2) f.add(p.ssBright, cyl(f, 0.009, 0.008, "z", 6), [0, 0, -dep / 2 + 0.02]);
        // The element: a sheathed resistance tube bent in two turns round the wheel, its ends taken back through the
        // wall at the top.
        const gap = 0.36;
        for (const [k, Re] of [R + 0.022, R + 0.036].entries()) {
          s.add(glow, torus(s, Re, 0.0042, "z", Math.PI * 2 - gap, 56, 6), [cx, fy, L.fanZ + (k ? 0.006 : -0.004)], [0, 0, Math.PI / 2 + gap / 2]);
          for (const sg of [-1, 1]) {
            const a = Math.PI / 2 + sg * (gap / 2);
            const [ex, ey] = [cx + Math.cos(a) * Re, fy + Math.sin(a) * Re];
            tube(s, glow, [[ex, ey, L.fanZ + (k ? 0.006 : -0.004)], [ex, ey, L.fanZ - 0.012], [ex, ey + 0.004, cav.wall]], 0.0042, 0.006, 8, false);
            if (fine2) s.add(p.ceramic, cyl(s, 0.007, 0.012, "z", 12), [ex, ey + 0.004, cav.wall + 0.006]);
          }
        }
        // Brackets clipping the element's turns to the wall.
        if (fine2)
          for (let i = 0; i < 3; i++) {
            const a = -Math.PI / 2 + ((i - 1) * Math.PI) / 2.2;
            const rr = R + 0.029;
            s.add(p.ss, box(0.008, 0.024, 0.002), [cx + Math.cos(a) * rr, fy + Math.sin(a) * rr, L.fanZ - 0.002], [0, 0, a + Math.PI / 2]);
            s.add(p.ss, box(0.008, 0.002, L.fanZ - cav.wall), [cx + Math.cos(a) * rr, fy + Math.sin(a) * rr, (L.fanZ + cav.wall) / 2]);
          }
        // Shaft and its seal flange on the wall.
        s.add(p.ssBright, cyl(s, 0.009, L.fanZ - dep / 2 - cav.wall, "z", 14), [cx, fy, (cav.wall + L.fanZ - dep / 2) / 2]);
        s.add(p.ssBright, cyl(s, 0.03, 0.004, "z", 24), [cx, fy, cav.wall + 0.002]);
        boltCircle(s, p.ssBright, [cx, fy, cav.wall + 0.004], 0.024, 4, 0.0024, Math.PI / 4);
        // The steam nozzle: a 6 mm pipe in through the wall above the wheel, down in front of it to the hub.
        {
          const zn = L.fanZ + dep / 2 + 0.007;
          const a = 0.6;
          const [ex, ey] = [cx + (R + 0.052) * Math.cos(a), fy + (R + 0.052) * Math.sin(a)];
          tube(s, p.ssBright, [[ex, ey, cav.wall], [ex, ey, zn], [cx + 0.02, fy + 0.014, zn]], 0.003, 0.012, 8);
          s.add(p.brass, cyl(s, 0.0045, 0.012, "x", 10, 0.003), [cx + 0.016, fy + 0.012, zn], [0, 0, 0.6]);
          if (fine2) s.add(p.brass, cyl(s, 0.007, 0.006, "z", 6), [ex, ey, cav.wall + 0.003]);
        }
      }
    // The motor behind the back wall: its cover bulges out of the back panel.
    for (const fy of L.fans) {
      s.add(p.ssDull, cyl(s, 0.07, 0.02, "z", 32), [cx, fy, zBack - 0.01]);
      if (fine1) s.add(p.ssDull, torus(s, 0.068, 0.003, "z", Math.PI * 2, 32, 6), [cx, fy, zBack - 0.019]);
      boltCircle(s, p.ssBright, [cx, fy, zBack - 0.02], 0.058, 4, 0.0026, Math.PI / 4, [0, Math.PI, 0]);
    }

    // The roll-in rack trolley (20): a stainless frame on high-temperature casters with 20 levels of runners, a drip
    // tray at its foot and a handle across its front.
    const trolleyRest: V3 = [0, 0, 0];
    if (big && fine1) {
      const t = s.joint("rack_trolley", trolleyRest);
      const zf = cav.z1 - 0.02;
      const zb = cav.z0 + 0.02;
      const xl = cx - FRAME - 0.012;
      const xr = cx + FRAME + 0.012;
      const yb = cav.y0 + 0.1;
      const yt = L.rails[L.rails.length - 1]! + 0.05;
      const post = 0.022;
      for (const x of [xl, xr])
        for (const z of [zf - post / 2, zb + post / 2]) {
          slab(t, p.ss, x - post / 2, x + post / 2, yb, yt, z - post / 2, z + post / 2, fine2 ? 0.003 : 0);
          t.at([0, cav.y0, 0], [0, 0, 0], () => wheelSet(t, p.poly, p.ss, [x, yb - cav.y0, z + 0.012], 0.038));
        }
      for (const y of [yb, yt - 0.012])
        for (const z of [zf - post / 2, zb + post / 2]) slab(t, p.ss, xl, xr, y, y + 0.012, z - 0.01, z + 0.01, fine2 ? 0.002 : 0);
      for (const y of [yb, yt - 0.012]) for (const x of [xl, xr]) slab(t, p.ss, x - 0.01, x + 0.01, y, y + 0.012, zb + post, zf - post, fine2 ? 0.002 : 0);
      if (fine1) slab(t, p.ss, xl + 0.012, xr - 0.012, yb + 0.012, yb + 0.0135, zb + post, zf - post, 0.0005);
      runnerSet(t, cx, zf - post / 2, zb + post / 2, L.rails, yb, yt, [-1, 1]);
      tube(t, p.ssBright, [[xl + 0.03, yt - 0.06, zf], [xl + 0.03, yt - 0.06, zf + 0.012], [xr - 0.03, yt - 0.06, zf + 0.012], [xr - 0.03, yt - 0.06, zf]], 0.009, 0.008, 12);
    }

    /*
     * The door, hinged left: a stainless chassis glazed three times. The outer pane, 6 mm toughened glass with a grey
     * tint, covers the whole front; a black frit fired onto its back hides the chassis and leaves a clear window, its
     * edge broken into the printer's halftone. The middle pane sits in the chassis behind a dark spacer; the inner
     * pane, which the gasket seals on, rides on a stainless carrier that swings out for cleaning between the panes.
     * Through the window one sees down into the door: the reveal, the spacer, the carrier's frame set well back, and
     * the cavity beyond it.
     */
    const door = s.joint("cavity_door", pivot);
    // Built from the door's back face at its left edge: x along the door, y absolute, z out through its thickness.
    door.push([-HINGE_IN, 0, L.zDoorBack - L.zDoorFront]);
    {
      const W = L.door.x1 - L.door.x0;
      const T = L.zDoorFront - L.zDoorBack;
      const [y0, y1] = [L.door.y0, L.door.y1];
      const H = y1 - y0;
      const ym = (y0 + y1) / 2;
      const dx = L.door.x0;
      // The clear window (door frame: x from the hinge pivot, y absolute): the cavity opening and a margin.
      const win = { x0: cav.x0 - dx - 0.024, x1: cav.x1 - dx + 0.024, y0: Math.max(cav.y0 - 0.03, y0 + 0.03), y1: cav.y1 + 0.024 };
      const [wcx, wcy, ww, wh] = [(win.x0 + win.x1) / 2, (win.y0 + win.y1) / 2, win.x1 - win.x0, win.y1 - win.y0];
      pane(door, fm.paneSmoke, fm.paneEdge, W, H, 0.006, [W / 2, ym, T - 0.003]);
      door.add(fm.frit, fritGeom(door, W, H, W / 2, ym, win, 0.022), [0, 0, T - 0.0066]);
      // The chassis: formed stainless round the window, its reveal seen through it at an angle.
      const fz0 = 0.0068;
      const fz1 = T - 0.0082;
      {
        const sh = roundRect(W - 0.002, H - 0.002, 0.004, W / 2, ym);
        rectPath(sh, ww + 0.012, wh + 0.012, 0.025, wcx, wcy, true);
        door.add(p.ss, plate(door, sh, fz1 - fz0, fine1 ? 0.0012 : 0), [0, 0, fz0]);
      }
      // The middle pane in its dark spacer.
      if (fine1) {
        const zm = fz0 + (fz1 - fz0) * 0.58;
        pane(door, fm.paneClear, fm.paneEdge, ww + 0.03, wh + 0.03, 0.004, [wcx, wcy, zm]);
        const sp = roundRect(ww + 0.016, wh + 0.016, 0.027, wcx, wcy);
        rectPath(sp, ww + 0.003, wh + 0.003, 0.0215, wcx, wcy, true);
        door.add(p.enamel, plate(door, sp, 0.012, fine2 ? 0.0008 : 0), [0, 0, zm - 0.006]);
      }
      // The inner pane, and the carrier that holds it: a formed frame whose opening is a little larger than the cavity's.
      const ix0 = cav.x0 - 0.034 - dx;
      const ix1 = cav.x1 + 0.034 - dx;
      const iy0 = Math.max(cav.y0 - 0.034, y0 + 0.012);
      const iy1 = cav.y1 + 0.034;
      const [icx, icy] = [(ix0 + ix1) / 2, (iy0 + iy1) / 2];
      pane(door, fm.paneClear, fm.paneEdge, ix1 - ix0, iy1 - iy0, 0.005, [icx, icy, 0.0025]);
      if (fine1) {
        const cr = roundRect(ix1 - ix0 - 0.002, iy1 - iy0 - 0.002, 0.006, icx, icy);
        rectPath(cr, cav.x1 - cav.x0 + 0.016, Math.min(cav.y1 - cav.y0 + 0.016, iy1 - iy0 - 0.03), 0.02, cx - dx, (cav.y0 + cav.y1) / 2, true);
        door.add(p.ss, plate(door, cr, 0.006, fine2 ? 0.001 : 0), [0, 0, 0.0069]);
      }
      // The door's back sheet round the inner pane (heat-tinted from the cavity side).
      {
        const bs = roundRect(W - 0.002, H - 0.002, 0.003, W / 2, ym);
        rectPath(bs, ix1 - ix0 - 0.03, iy1 - iy0 - 0.03, 0.006, icx, icy, true);
        door.add(p.ssDull, plate(door, bs, 0.0012), [0, 0, 0.0055]);
      }
      if (fine2) {
        // The carrier's hinge knuckles at the left, its quarter-turn catches at the right, and the condensate lip
        // along the inner pane's foot.
        for (const yy of [iy0 + 0.1, iy1 - 0.1]) {
          door.add(p.ssBright, cyl(door, 0.004, 0.034, "y", 12), [ix0 + 0.006, yy, 0.0005 - 0.004]);
          slab(door, p.ssBright, ix0 + 0.004, ix0 + 0.03, yy - 0.012, yy + 0.012, -0.0012, 0.0055, 0.001);
        }
        for (const yy of [iy0 + 0.08, iy1 - 0.08]) {
          slab(door, p.ssBright, ix1 - 0.012, ix1 + 0.008, yy - 0.012, yy + 0.012, -0.0012, 0.0055, 0.001);
          door.add(p.ssBright, cyl(door, 0.0065, 0.004, "z", 16), [ix1 - 0.004, yy, -0.0032]);
          door.add(p.enamel, box(0.009, 0.0015, 0.0012), [ix1 - 0.004, yy, -0.0058], [0, 0, 0.5]);
        }
        slab(door, p.ss, ix0 + 0.02, ix1 - 0.02, iy0 - 0.006, iy0 + 0.002, -0.004, 0.0055, 0.001);
        // The latch: a sprung tongue with a roller that the frame's striker catches.
        slab(door, p.ssBright, W - 0.03, W - 0.012, ym - 0.03, ym + 0.03, -0.005, 0.0068, 0.002);
        door.add(p.rubber, cyl(door, 0.006, 0.014, "x", 14), [W - 0.021, ym, -0.0075]);
      }
      // The door's hinge knuckles round the pins, at its outer face, each on a strap screwed along its left edge.
      for (const [k0, k1] of [[y1 - 0.034, y1 + 0.0001], [y0 - 0.0001, y0 + 0.034]] as const) {
        door.add(p.ssBright, cyl(door, 0.0075, k1 - k0, "y", s.seg(20, 10)), [HINGE_IN, (k0 + k1) / 2, T]);
        if (fine1) bit(door, p.ssBright, -0.003, 0.0015, k0 + 0.002, k1 - 0.002, 0.008, T - 0.002, 0.0008);
        if (fine2) for (const z of [0.016, T * 0.55]) screw(door, p.ssBright, [-0.003, (k0 + k1) / 2, z], 0.0018, [0, -Math.PI / 2, 0]);
      }
      // The handle: a heavy vertical grip on cast standoffs at the door's right edge.
      const hl = big ? 0.5 : pans === 10 ? 0.36 : 0.26;
      const hx = W - 0.03;
      const hz = T;
      const hy = big ? Math.min(ym, 1.0) : ym;
      for (const sg of [-1, 1]) {
        door.add(p.ssBright, lathe(door, [[0, 0], [0.019, 0], [0.019, 0.004], [0.015, 0.012], [0.013, 0.026], [0, 0.026]], 24), [hx, hy + sg * (hl / 2), hz], [Math.PI / 2, 0, 0]);
        if (fine2) screw(door, p.ssBright, [hx, hy + sg * (hl / 2 - 0.012), hz + 0.0005], 0.002);
      }
      tube(door, p.bakelite, [[hx, hy + hl / 2, hz + 0.018], [hx, hy + hl / 2 - 0.035, hz + 0.026], [hx, hy - hl / 2 + 0.035, hz + 0.026], [hx, hy - hl / 2, hz + 0.018]], 0.0115, 0.03, 16);
    }
    door.pop();

    /* The right-hand column: side panel, cavity wall and controls, all hidden by the cutaway. */
    const cut = s.joint("cutaway", [0, 0, 0]);
    {
      slab(cut, p.ss, xR - sk, xR, base, top - sk, zBack, zFrame, 0.0005);
      if (fine1) cut.add(p.ss, sweepZ(linerR, cav.wall, cav.z1, ct));
      // Ventilation louvres for the electronics at the side's lower rear, the data plate at its lower front.
      if (fine1) louvres(cut, [xR + 0.0001, base + 0.12, zBack + 0.12], [0, Math.PI / 2, 0], 6, 0.11, 0.016);
      dataPlate(cut, [xR, base + 0.09, zFrame - 0.12], [0, Math.PI / 2, 0]);
      // The control column: a stainless case, a black glass fascia with the touch screen and the central dial.
      const P = L.panel;
      const [py0, py1] = [L.door.y0, L.door.y1];
      slab(cut, p.ss, P.x0, P.x1, py0, py1, zFrame, P.z - 0.002, 0.004);
      const fx0 = P.x0 + 0.008;
      const fx1 = P.x1 - 0.008;
      const fw = fx1 - fx0;
      const sw = fw - 0.024;
      const shh = sw * 1.25;
      const sTop = big ? Math.min(py1 - 0.06, 1.72) : pans === 10 ? Math.min(py1 - 0.05, 1.68) : py1 - 0.03;
      const fy1 = Math.min(py1 - 0.008, sTop + 0.028);
      const dialY = sTop - shh - 0.075;
      const fy0 = pans === 6 ? py0 + 0.008 : dialY - 0.14;
      slab(cut, fm.blackGlass, fx0, fx1, fy0, fy1, P.z - 0.003, P.z, 0.003);
      const pcx = (fx0 + fx1) / 2;
      slab(cut, p.screen, pcx - sw / 2, pcx + sw / 2, sTop - shh, sTop, P.z - 0.001, P.z + 0.0008, 0.0015);
      if (fine1) cut.add(fm.combiUi, new THREE.PlaneGeometry(sw - 0.006, shh - 0.007), [pcx, sTop - shh / 2, P.z + 0.0014]);
      // The dial: a knurled stainless ring round a black grip, pushed to confirm; a lit ring round it.
      const dial = cut.joint("dial", [pcx, dialY, P.z]);
      dial.add(p.ssBright, lathe(dial, [[0, 0], [0.034, 0], [0.035, 0.004], [0.034, 0.02], [0.031, 0.023], [0, 0.023]], 48), [0, 0, 0], [Math.PI / 2, 0, 0]);
      dial.add(p.bakelite, lathe(dial, [[0, 0.0235], [0.027, 0.0235], [0.028, 0.025], [0.02, 0.027], [0, 0.0255]], 40), [0, 0, 0], [Math.PI / 2, 0, 0]);
      if (fine2) {
        for (let i = 0; i < 36; i++) {
          const a = (i / 36) * Math.PI * 2;
          dial.add(p.ssBright, box(0.002, 0.0026, 0.014), [Math.cos(a) * 0.0348, Math.sin(a) * 0.0348, 0.012], [0, 0, a]);
        }
        dial.add(p.ssBright, cyl(dial, 0.0022, 0.001, "z", 10), [0, 0.019, 0.0262]);
      }
      if (fine1) cut.add(p.ledAmber, torus(cut, 0.04, 0.0012, "z", Math.PI * 2, 48, 4), [pcx, dialY, P.z + 0.0004]);
      if (fine2) {
        // A USB port under a rubber flap, and the status light.
        slab(cut, p.rubber, pcx - 0.012, pcx + 0.012, dialY - 0.07, dialY - 0.058, P.z, P.z + 0.003, 0.0015);
        slab(cut, p.ledGreen, pcx - 0.018, pcx + 0.018, dialY - 0.085, dialY - 0.0835, P.z, P.z + 0.0007);
      }
      // Below the fascia (tall units): the column's face, a louvred service panel at its foot.
      if (fine1 && pans !== 6) louvres(cut, [pcx, py0 + 0.12, P.z - 0.002], [0, 0, 0], 8, fw - 0.06, 0.016);
      if (fine2 && pans !== 6) for (const sy of [py0 + 0.03, fy0 - 0.03]) for (const sx of [P.x0 + 0.02, P.x1 - 0.02]) screw(cut, p.ssBright, [sx, sy, P.z - 0.002], 0.0024);
      // The hand shower: a spray gun parked nozzle up in a moulded cradle, its hose drawn off a retracting reel
      // inside the column through a grommet beside it, hanging in a loop.
      if (fine2) {
        const hy = pans === 6 ? py0 + 0.13 : py0 + 0.25;
        const hz = P.z - 0.002;
        cut.add(p.plastic, rbox(cut, 0.046, 0.074, 0.006, 0.0025), [pcx, hy - 0.004, hz + 0.003]);
        cut.add(p.plastic, rbox(cut, 0.036, 0.012, 0.024, 0.004), [pcx, hy - 0.036, hz + 0.014]);
        // The gun: its grip, the head bent forward, a stainless nozzle with its spray face, and the trigger.
        const gz = hz + 0.024;
        cut.add(p.bakelite, rbox(cut, 0.026, 0.07, 0.03, 0.011), [pcx, hy - 0.006, gz]);
        cut.at([pcx, hy + 0.026, gz + 0.002], [0.38, 0, 0], () => {
          cut.add(p.bakelite, lathe(cut, [[0, 0], [0.0135, 0], [0.0135, 0.008], [0.0115, 0.018], [0.0092, 0.026], [0, 0.026]], 24));
          cut.add(p.ssBright, lathe(cut, [[0, 0.025], [0.0092, 0.025], [0.0096, 0.03], [0.0088, 0.038], [0.0082, 0.039], [0, 0.039]], 24));
          cut.add(p.enamel, cyl(cut, 0.0062, 0.0006, "y", 16), [0, 0.0392, 0]);
        });
        tube(cut, p.charcoal, [[pcx, hy + 0.018, gz + 0.016], [pcx, hy + 0.012, gz + 0.025], [pcx, hy - 0.02, gz + 0.023], [pcx, hy - 0.03, gz + 0.016]], 0.0035, 0.008, 8);
        cut.add(p.ssBright, cyl(cut, 0.0072, 0.012, "y", 16), [pcx, hy - 0.046, gz]);
        // The hose, and the grommet it runs through.
        const gx = pcx + 0.048;
        const gy = hy - 0.055;
        cut.add(p.ssBright, torus(cut, 0.0095, 0.0024, "z", Math.PI * 2, 24, 8), [gx, gy, hz + 0.0016]);
        cut.add(p.rubber, cyl(cut, 0.0072, 0.0025, "z", 16), [gx, gy, hz + 0.00125]);
        tube(cut, p.rubber, [[pcx, hy - 0.052, gz], [pcx + 0.003, hy - 0.088, gz + 0.006], [pcx + 0.026, hy - 0.104, gz + 0.004], [gx, hy - 0.09, hz + 0.016], [gx, gy, hz + 0.008], [gx, gy, hz - 0.004]], 0.0058, 0.014, 12);
      }
      // The core-temperature probe: a needle on a silicone grip, parked in its holder on the right wall at the front,
      // its cord hanging in a loop down the wall and up to its socket near the top.
      if (fine2) {
        const px = cav.x1 - 0.011;
        const pz = cav.z1 - 0.013;
        const yh = cav.y0 + (cav.y1 - cav.y0) * 0.42;
        const sy = cav.y1 - 0.07;
        const sz = cav.z1 - 0.05;
        slab(cut, p.ss, cav.x1 - 0.003, cav.x1, yh - 0.034, yh + 0.002, pz - 0.009, pz + 0.009, 0.001);
        slab(cut, p.ss, px + 0.003, cav.x1 - 0.003, yh - 0.028, yh - 0.004, pz - 0.0012, pz + 0.0012, 0.0004);
        cut.add(p.ss, cyl(cut, 0.0042, 0.03, "y", 16, 0.0042, true), [px, yh - 0.015, pz]);
        cut.add(p.ssBright, torus(cut, 0.0042, 0.0008, "y", Math.PI * 2, 16, 4), [px, yh, pz]);
        cut.add(p.ssBright, cyl(cut, 0.0016, 0.124, "y", 10), [px, yh + 0.008 - 0.062, pz]);
        cut.add(p.ssBright, cyl(cut, 0.0002, 0.008, "y", 10, 0.0016), [px, yh - 0.12, pz]);
        cut.add(p.gasket, lathe(cut, [[0, 0], [0.0042, 0], [0.0062, 0.006], [0.0066, 0.034], [0.0058, 0.046], [0.004, 0.054], [0.0026, 0.068], [0, 0.07]], 20), [px, yh + 0.006, pz]);
        for (let i = 0; i < 5; i++) cut.add(p.gasket, torus(cut, 0.0066, 0.0007, "y", Math.PI * 2, 18, 4), [px, yh + 0.016 + i * 0.005, pz]);
        tube(
          cut,
          p.gasket,
          [
            [px, yh + 0.074, pz],
            [px, yh + 0.098, pz - 0.006],
            [px - 0.001, yh + 0.082, pz - 0.024],
            [cav.x1 - 0.013, cav.y0 + 0.09, pz - 0.032],
            [cav.x1 - 0.014, cav.y0 + 0.062, sz + 0.004],
            [cav.x1 - 0.015, cav.y0 + 0.1, sz],
            [cav.x1 - 0.016, sy - 0.03, sz],
            [cav.x1 - 0.018, sy, sz],
          ],
          0.0021,
          0.018,
          8,
          false,
        );
        cut.add(p.ssBright, cyl(cut, 0.012, 0.004, "x", 24), [cav.x1 - 0.002, sy, sz]);
        cut.add(p.plastic, cyl(cut, 0.0068, 0.014, "x", 16, 0.0085), [cav.x1 - 0.011, sy, sz]);
        for (let i = 0; i < 6; i++) {
          const a = (i / 6) * Math.PI * 2;
          cut.add(p.plastic, box(0.01, 0.0016, 0.0016), [cav.x1 - 0.011, sy + Math.cos(a) * 0.0082, sz + Math.sin(a) * 0.0082], [a, 0, 0]);
        }
      }
      // The door's striker on the frame.
      if (fine2) slab(cut, p.ssBright, L.door.x1 - 0.034, L.door.x1 - 0.008, (py0 + py1) / 2 - 0.04, (py0 + py1) / 2 + 0.04, zFrame, zFrame + 0.004, 0.0015);
    }

    /* On top: the vapour outlet, the air intake beside it, and the louvred vent over the electronics. */
    {
      const vx = cav.x0 + 0.07;
      const vz = cav.wall + 0.05;
      const vh = L.h - 0.004 - top;
      s.add(p.ss, lathe(s, [[0.026, 0], [0.03, 0], [0.03, vh - 0.006], [0.033, vh - 0.003], [0.03, vh], [0.026, vh - 0.002], [0.026, 0.01]], 32), [vx, top, vz]);
      s.add(p.enamel, cyl(s, 0.026, 0.002, "y", 24), [vx, top + 0.012, vz]);
      s.add(p.ss, lathe(s, [[0.014, 0], [0.017, 0], [0.017, vh * 0.6], [0.014, vh * 0.6], [0.014, 0.006]], 24), [vx + 0.085, top, vz]);
      s.add(p.enamel, cyl(s, 0.014, 0.002, "y", 16), [vx + 0.085, top + 0.008, vz]);
      if (fine1) louvres(s, [w - 0.12, top + 0.0001, zBack + 0.11], [-Math.PI / 2, 0, 0], 7, 0.13, 0.017);
      if (fine2)
        for (const x of [0.03, w - 0.03])
          for (const z of [zBack + 0.02, L.zDoorFront - 0.02]) screw(s, p.ssBright, [x, top, z], 0.0026, [-Math.PI / 2, 0, 0]);
    }

    /* At the back: water inlets, the drain, the supply cable's gland and the bonding stud. */
    if (fine1) {
      const zb = zBack - sk;
      for (const [i, x] of [w - 0.1, w - 0.16].entries()) {
        s.add(p.brass, cyl(s, 0.011, 0.018, "z", 6), [x, base + 0.08, zb - 0.009]);
        s.add(p.brass, cyl(s, 0.008, 0.008, "z", 16), [x, base + 0.08, zb - 0.022]);
        if (fine2) s.add(i ? p.ledGreen : fm.red, cyl(s, 0.0115, 0.002, "z", 16), [x, base + 0.08, zb - 0.001]);
      }
      // The drain: down through the base on the stand and countertop models; the floor model's leaves at the back.
      if (base > 0.03) {
        s.add(p.poly, cyl(s, 0.025, base - 0.014, "y", 20), [cx, (base + 0.014) / 2, zBack + 0.06]);
        s.add(p.ssBright, torus(s, 0.026, 0.003, "y", Math.PI * 2, 24, 6), [cx, base - 0.02, zBack + 0.06]);
      } else {
        tube(s, p.poly, [[cx, 0.06, zBack + 0.004], [cx, 0.06, -L.d + 0.012], [cx, 0.001, -L.d + 0.012]], 0.011, 0.014, 16);
        s.add(p.ssBright, torus(s, 0.012, 0.002, "z", Math.PI * 2, 24, 6), [cx, 0.06, zb - 0.002]);
      }
      s.add(p.bakelite, cyl(s, 0.014, 0.016, "z", 6), [w - 0.08, top - 0.1, zb - 0.008]);
      tube(s, p.rubber, [[w - 0.08, top - 0.1, zb - 0.016], [w - 0.08, top - 0.12, zb - 0.014], [w - 0.06, base + 0.02, zb - 0.014], [w - 0.06, 0.006, zb - 0.014]], 0.009, 0.03, 10);
      if (fine2) s.add(p.brass, cyl(s, 0.004, 0.02, "z", 8), [0.06, base + 0.06, zb - 0.01]);
      dataPlate(s, [w - 0.2, top - 0.08, zb], [0, Math.PI, 0]);
    }

    /* Steam in the cavity. */
    const haze = new Haze({ x0: cav.x0 + 0.1, x1: cav.x1 - 0.1, y0: cav.y0 + 0.1, y1: cav.y1 - 0.1, z0: cav.z0 + 0.06, z1: cav.z1 - 0.08 }, pans === 6 ? 22 : pans === 10 ? 34 : 50, 0.36, pans, 0, 0xf4f6f8, 0.42);
    if (fine1) {
      s.extra(haze.group);
      s.tris += haze.tris;
    }

    return (b) => {
      layerGlass(b.group);
      const fans = fanJoints.map((n) => b.joints.get(n));
      const doorJ = b.joints.get("cavity_door");
      const cutJ = b.joints.get("cutaway");
      const trolley = b.joints.get("rack_trolley");
      let speed = 0;
      let angle = 0;
      setGlow(glow, 0);
      // The trolley's wheels: front and rear contact points at rest; past the opening the floor drops to the kitchen's.
      const wf = cav.z1 - 0.02 - 0.011 + 0.012 - 0.038 * 0.4;
      const wr = cav.z0 + 0.02 + 0.011 + 0.012 - 0.038 * 0.4;
      const floorAt = (z: number) => cav.y0 * (1 - THREE.MathUtils.smoothstep(z, cav.z1 - 0.005, cav.z1 + 0.04));
      const X = new THREE.Vector3(1, 0, 0);
      const pivotRest = new THREE.Vector3(cx, cav.y0, wr);
      const rollTrolley = (j: THREE.Object3D, dz: number) => {
        const [hf, hr] = [floorAt(wf + dz), floorAt(wr + dz)];
        const th = Math.asin(THREE.MathUtils.clamp((hr - hf) / (wf - wr), -0.2, 0.2));
        j.quaternion.setFromAxisAngle(X, th);
        const pv = pivotRest.clone().add(new THREE.Vector3(0, hr - cav.y0, dz));
        j.position.copy(pv.sub(pivotRest.clone().applyQuaternion(j.quaternion)));
      };
      const hooks: Hooks = {
        glows: [glow],
        open: (f) => turnJoint(doorJ, "y", -Math.max(0, Math.min(1, f)) * 2.356),
        set: (name, v) => {
          if (name === "fan") speed = Math.max(0, Math.min(1, v));
          else if (name === "heat") setGlow(glow, v);
          else if (name === "steam") haze.set(v);
          else if (name === "cutaway" && cutJ) cutJ.visible = v < 0.5;
          else if (name === "dial") turnJoint(b.joints.get("dial"), "z", -v * Math.PI * 1.5);
          else if (name === "trolley" && trolley) rollTrolley(trolley, Math.max(0, Math.min(1, v)) * 0.62);
        },
        tick: (t, dt) => {
          angle += Math.max(0, Math.min(0.1, dt)) * speed * 14;
          for (const [i, f] of fans.entries()) turnJoint(f, "z", i % 2 ? -angle : angle);
          haze.tick(t);
        },
      };
      return hooks;
    };
  };
}

/** A GN 1/1 or US hotel pan, rim on y = 0, `depth` deep: rolled rim, drawn walls with a draft, rounded corners. */
function gnPan(s: Shop, w: number, d: number, depth: number): THREE.BufferGeometry {
  const t = 0.0008;
  const lip = 0.008;
  const dr = depth * 0.08;
  // The outer wall starts at the curl's foot, `lip` in: a sheet t thick, its inner face never outside its outer one.
  if (!s.fine(2))
    return vessel(s, w, d, 0.028, [
      [lip + dr + 0.008, -depth],
      [lip, -0.006],
      [0, -0.002],
      [0.0015, 0],
      [lip + t, -0.002],
      [lip + dr + t + 0.008, -depth + t],
    ]);
  return vessel(s, w, d, 0.028, [
    [lip + dr + 0.012, -depth],
    [lip + dr + 0.003, -depth + 0.004],
    [lip, -0.006],
    [lip - 0.002, -0.0035],
    [0.0015, -0.003],
    [0, -0.0015],
    [0.0015, 0],
    [lip - 0.002, 0],
    [lip + t, -0.002],
    [lip + dr + t + 0.003, -depth + 0.004 + t],
    [lip + dr + t + 0.012, -depth + t],
  ]);
}

/** An 18 × 26 in (or half-size) sheet pan, rim on y = 0: a wire-beaded rim over a 25 mm wall. */
function sheetPan(s: Shop, w: number, d: number): THREE.BufferGeometry {
  if (!s.fine(2))
    return vessel(s, w, d, 0.012, [
      [0.006, -0.025],
      [0.0012, -0.005],
      [-0.0026, -0.0018],
      [0.0012, 0.0004],
      [0.0036, -0.0236],
      [0.0066, -0.0246],
    ]);
  return vessel(s, w, d, 0.012, [
    [0.006, -0.025],
    [0.003, -0.0235],
    [0.0012, -0.005],
    [-0.0012, -0.0038],
    [-0.0026, -0.0018],
    [-0.0016, 0.0004],
    [0.0012, 0.0004],
    [0.0026, -0.002],
    [0.0036, -0.0236],
    [0.0066, -0.0246],
  ]);
}

/** Dough rolls proofing on a sheet pan's floor (y = 0 at the floor): domes close up, a lumpy sheet of them further off. */
function rolls(s: Shop, m: Mat, w: number, d: number, nx: number, nz: number, r: number, seed: number): void {
  const rr = rng(seed);
  const px = w / nx;
  const pz = d / nz;
  if (s.fine(2)) {
    for (let i = 0; i < nx; i++)
      for (let j = 0; j < nz; j++) {
        const k = 0.9 + rr() * 0.2;
        const g = lathe(s, [[r * k * 0.9, 0], [r * k * 1.02, r * 0.3], [r * k * 0.82, r * 0.8], [r * k * 0.4, r * 1.08], [0, r * 1.12]], 10);
        s.add(m, g, [-w / 2 + px * (i + 0.5) + (rr() - 0.5) * 0.006, 0, -d / 2 + pz * (j + 0.5) + (rr() - 0.5) * 0.006], [0, rr() * 3, 0], [1, k, 1 + (rr() - 0.5) * 0.08]);
      }
    return;
  }
  const g = new THREE.PlaneGeometry(w, d, nx * 2, nz * 2);
  g.rotateX(-Math.PI / 2);
  const p = g.attributes.position as THREE.BufferAttribute;
  for (let i = 0; i < p.count; i++) {
    const [x, z] = [p.getX(i) + w / 2, p.getZ(i) + d / 2];
    const bump = Math.max(0, Math.sin((x / px) * Math.PI)) * Math.max(0, Math.sin((z / pz) * Math.PI));
    p.setY(i, r * 1.1 * Math.sqrt(bump));
  }
  g.computeVertexNormals();
  s.add(m, g);
}

/**
 * An appliance knob facing +z from z = 0: the kit's skirted profile on a bezel, with its pointer ridge and (close up)
 * the grip flutes, cheap enough to use by the row; a plain tapered puck from far off.
 */
export function dial(s: Shop, body: Mat, bezel: Mat, at: V3, r = 0.022): void {
  s.at(at, [0, 0, 0], () => {
    if (s.detail === 0) {
      s.add(body, new THREE.CylinderGeometry(r * 0.7, r, 0.027, 8).rotateX(Math.PI / 2), [0, 0, 0.0155]);
      return;
    }
    s.add(bezel, cyl(s, r * 1.18, 0.002, "z", 32), [0, 0, 0.001]);
    const g = lathe(s, [[0, 0.002], [r, 0.002], [r, 0.006], [r * 0.82, 0.009], [r * 0.72, 0.012], [r * 0.7, 0.024], [r * 0.6, 0.028], [0, 0.029]], s.fine(2) ? 36 : 12);
    g.rotateX(Math.PI / 2);
    s.add(body, g);
    s.add(body, s.fine(2) ? rbox(s, r * 0.28, r * 1.3, 0.012, 0.002) : box(r * 0.28, r * 1.3, 0.012), [0, r * 0.2, 0.022]);
    if (s.fine(2))
      for (let i = 0; i < 12; i++) {
        const a = (i / 12) * Math.PI * 2;
        s.add(body, box(r * 0.16, r * 0.16, 0.011), [Math.cos(a) * r * 0.7, Math.sin(a) * r * 0.7, 0.018], [0, 0, a]);
      }
  });
}

/** A round stainless leg with its bullet foot; a plain post from far off. */
function leg(s: Shop, m: Mat, x: number, z: number, y0: number, y1: number, r = 0.0205): void {
  if (s.detail === 0) s.add(m, cyl(s, r, y1 - y0, "y", 8), [x, (y0 + y1) / 2, z]);
  else roundLeg(s, m, x, z, y0, y1, r);
}

/** A lamp whose light is live: one material per model, dark when off. */
function liveLamp(on: number): THREE.MeshBasicMaterial {
  const m = new THREE.MeshBasicMaterial({ color: on });
  m.userData.uv = "own";
  m.userData.on = on;
  return m;
}
function setLamp(m: THREE.MeshBasicMaterial, v: number): void {
  const c = new THREE.Color(m.userData.on as number);
  const off = new THREE.Color(0x2a2724);
  m.color.copy(off.lerp(c, Math.max(0, Math.min(1, v))));
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Heated holding and proofing cabinet                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A full-height insulated holding and proofing cabinet: a double-walled body on 5 in casters inside a perimeter
 * bumper; a full-height door on three barrel hinges with a double-glazed window and a magnetic gasket; universal
 * slides for sheet or hotel pans; a top-mounted module with the heater, fan and controls, its air returned down a
 * perforated duct on the back wall; a water pan on the floor for humidity.
 */
const holdingCabinet: Builder = (s, e) => {
  const p = palette();
  const fm = mats();
  const { w, d, h } = e;
  const fine1 = s.fine(1);
  const fine2 = s.fine(2);
  const [bx0, bx1] = [0.02, w - 0.02];
  const zf = -0.068;
  const zb = -d + 0.048;
  const yb = 0.17;
  const yMod = h - 0.178;
  const doorY: [number, number] = [yb + 0.045, yMod - 0.006];
  const zDoorBack = zf + 0.006;
  const zDoorFront = -0.016;
  const cav = { x0: 0.07, x1: w - 0.07, y0: yb + 0.065, y1: yMod - 0.03, z0: zb + 0.035, z1: zf - 0.002 };
  const cx = w / 2;

  // Casters and the base pan; the bumper round it.
  for (const x of [bx0 + 0.07, bx1 - 0.07])
    for (const [i, z] of [zf - 0.07, zb + 0.07].entries()) {
      wheelSet(s, p.poly, p.ss, [x, yb, z], 0.0635);
      if (fine2 && i === 0) s.add(fm.red, box(0.03, 0.005, 0.022), [x, 0.14, z + 0.045], [-0.3, 0, 0]);
    }
  slab(s, p.ss, bx0, bx1, yb, yb + 0.05, zb, zf, fine2 ? 0.003 : 0);
  s.at([0, yb + 0.025, 0], [Math.PI / 2, 0, 0], () => rectLoop(s, p.rubber, 0.02, w - 0.02, zb, zf, 0, 0.02, 0.03, fine2 ? 12 : 6));
  // The body: double skin, its front frame round the opening.
  slab(s, p.ss, bx0, bx0 + 0.0012, yb + 0.05, yMod, zb, zf, fine2 ? 0.0005 : 0);
  slab(s, p.ss, bx1 - 0.0012, bx1, yb + 0.05, yMod, zb, zf, fine2 ? 0.0005 : 0);
  slab(s, p.ssDull, bx0, bx1, yb + 0.05, yMod, zb - 0.0012, zb);
  {
    const sh = roundRect(bx1 - bx0, yMod - yb - 0.05, 0.002, cx, (yMod + yb + 0.05) / 2);
    rectPath(sh, cav.x1 - cav.x0, cav.y1 - cav.y0, 0.01, cx, (cav.y0 + cav.y1) / 2, true);
    s.add(p.ss, plate(s, sh, 0.002), [0, 0, zf - 0.002]);
  }
  // The cavity liner (aluminium in most, stainless here) and its universal slides.
  if (fine1) {
    slab(s, p.ssDull, cav.x0 - 0.0012, cav.x0, cav.y0, cav.y1, cav.z0, cav.z1);
    slab(s, p.ssDull, cav.x1, cav.x1 + 0.0012, cav.y0, cav.y1, cav.z0, cav.z1);
    slab(s, p.ssDull, cav.x0, cav.x1, cav.y1, cav.y1 + 0.0012, cav.z0, cav.z1);
    slab(s, p.ssDull, cav.x0, cav.x1, cav.y0 - 0.0012, cav.y0, cav.z0, cav.z1);
    slab(s, p.ssHeat, cav.x0, cav.x1, cav.y0, cav.y1, cav.z0 - 0.0012, cav.z0);
    // The return duct down the back wall: perforated sheet in a formed frame.
    const [dy0, dy1] = [cav.y0 + 0.08, cav.y1 - 0.04];
    s.add(fm.perf, new THREE.PlaneGeometry(cav.x1 - cav.x0 - 0.1, dy1 - dy0), [cx, (dy0 + dy1) / 2, cav.z0 + 0.03]);
    for (const x of [cav.x0 + 0.045, cav.x1 - 0.045]) slab(s, p.ssDull, x - 0.006, x + 0.006, dy0, dy1, cav.z0, cav.z0 + 0.031);
    for (const y of [dy0, dy1]) slab(s, p.ssDull, cav.x0 + 0.04, cav.x1 - 0.04, y - 0.006, y + 0.006, cav.z0, cav.z0 + 0.031);
  }
  const rails: number[] = [];
  for (let y = cav.y0 + 0.1; y < cav.y1 - 0.06; y += 0.0762) rails.push(y);
  if (fine1)
    for (const sg of [-1, 1]) {
      const xu = cx + sg * 0.248;
      const xi = cx + sg * 0.22;
      for (const z of [cav.z1 - 0.04, cav.z0 + 0.07]) slab(s, p.ss, xu - 0.008, xu + 0.008, cav.y0 + 0.03, cav.y1 - 0.02, z - 0.0125, z + 0.0125, fine2 ? 0.001 : 0);
      for (const y of rails) {
        slab(s, p.ss, Math.min(xu, xi), Math.max(xu, xi), y - 0.0016, y, cav.z0 + 0.05, cav.z1 - 0.025);
        if (fine2) slab(s, p.ss, xu - 0.0008 + sg * 0.008, xu + 0.0008 + sg * 0.008, y - 0.018, y, cav.z0 + 0.05, cav.z1 - 0.025);
      }
    }
  // Its load: rolls proofing on two sheet pans, an empty pan above (none in the "-EMPTY" variant, which a view fills).
  if (fine1 && s.name !== "NSF-HOLDING-CABINET-EMPTY") {
    const pz = (cav.z0 + 0.05 + cav.z1 - 0.025) / 2;
    const loads = fine2 ? [5, 9, 12] : [5, 9];
    for (const [k, li] of loads.entries()) {
      const y = rails[li]!;
      if (y === undefined) continue;
      s.add(p.alu, sheetPan(s, 0.46, 0.66), [cx, y + 0.0038, pz]);
      if (k < 2) s.at([cx, y + 0.0038 - 0.0246 + 0.0006, pz], [0, 0, 0], () => rolls(s, fm.dough, 0.42, 0.62, 4, 7, 0.03, 40 + k));
    }
  }
  // The water pan on the floor, its immersion element along its bottom.
  if (fine1) {
    const wz = (cav.z0 + cav.z1) / 2;
    s.add(p.ss, gnPan(s, 0.3, 0.46, 0.065), [cx, cav.y0 + 0.066, wz]);
    s.add(p.water, box(0.27, 0.0016, 0.43), [cx, cav.y0 + 0.035, wz]);
    if (fine2) tube(s, p.ssHeat, [[cx - 0.08, cav.y0 + 0.07, wz - 0.24], [cx - 0.08, cav.y0 + 0.012, wz - 0.18], [cx - 0.08, cav.y0 + 0.012, wz + 0.15], [cx + 0.08, cav.y0 + 0.012, wz + 0.15], [cx + 0.08, cav.y0 + 0.012, wz - 0.18], [cx + 0.08, cav.y0 + 0.07, wz - 0.24]], 0.004, 0.03, 8);
  }
  const vapour = new Haze({ x0: cx - 0.12, x1: cx + 0.12, y0: cav.y0 + 0.06, y1: cav.y0 + 0.45, z0: (cav.z0 + cav.z1) / 2 - 0.18, z1: (cav.z0 + cav.z1) / 2 + 0.18 }, 8, 0.16, 81, 0.25, 0xffffff, 0.3);
  if (fine1) {
    s.extra(vapour.group);
    s.tris += vapour.tris;
  }

  // The top module: heater, fan and humidity control behind a louvred case, its panel above the door.
  slab(s, p.ss, bx0, bx1, yMod, h - 0.002, zb, zf, 0.004);
  slab(s, p.ss, bx0 + 0.004, bx1 - 0.004, yMod + 0.004, h - 0.006, zf, zDoorFront - 0.003, 0.004);
  const pz = zDoorFront - 0.003;
  const pY = (yMod + h) / 2;
  slab(s, p.bakelite, bx0 + 0.03, bx1 - 0.03, yMod + 0.025, h - 0.028, pz - 0.002, pz + 0.001, 0.003);
  s.add(fm.holdDisplay, new THREE.PlaneGeometry(0.15, 0.045), [cx, pY + 0.012, pz + 0.0018]);
  slab(s, p.screen, cx - 0.08, cx + 0.08, pY - 0.016, pY + 0.04, pz, pz + 0.0012, 0.002);
  const heatLamp = liveLamp(0xffa040);
  const humLamp = liveLamp(0x4ade80);
  const knobs: string[] = [];
  for (const [i, x] of [bx0 + 0.11, bx1 - 0.11].entries()) {
    const name = i ? "humidity_knob" : "temp_knob";
    const j = s.joint(name, [x, pY, pz + 0.001]);
    dial(j, p.bakelite, p.ssBright, [0, 0, 0], 0.022);
    knobs.push(name);
    if (fine1) s.add(i ? humLamp : heatLamp, cyl(s, 0.0035, 0.003, "z", 12), [x + (i ? -0.05 : 0.05), pY + 0.03, pz + 0.002]);
    if (fine2) s.add(p.enamel, cyl(s, 0.032, 0.0006, "z", 32), [x, pY, pz + 0.0017]);
  }
  // Power switch.
  if (fine1) slab(s, p.plastic, cx + 0.12, cx + 0.145, pY - 0.035, pY - 0.005, pz, pz + 0.006, fine2 ? 0.002 : 0);
  if (fine1) for (const sg of [-1, 1]) louvres(s, [sg < 0 ? bx0 : bx1, (yMod + h) / 2, (zb + zf) / 2], [0, sg * Math.PI / 2, 0], 5, 0.25, 0.022);
  // Push handles on the back corners.
  for (const x of [bx0 + 0.07, bx1 - 0.07]) {
    tube(s, p.ss, [[x, 0.85, zb], [x, 0.85, -d + 0.016], [x, 1.45, -d + 0.016], [x, 1.45, zb]], 0.0135, 0.03, 14);
    if (fine1) for (const y of [0.85, 1.45]) s.add(p.ss, cyl(s, 0.022, 0.004, "z", 20), [x, y, zb - 0.002]);
  }
  dataPlate(s, [bx1, yb + 0.15, zf - 0.12], [0, Math.PI / 2, 0]);
  if (fine1) tube(s, p.rubber, [[bx1 - 0.1, h - 0.1, zb], [bx1 - 0.1, h - 0.1, zb - 0.02], [bx1 - 0.1, yb + 0.06, zb - 0.02]], 0.005, 0.03, 8);

  // The door, hinged on its left edge: insulated panel, double-glazed window, magnetic gasket on its back.
  const door = s.joint("cabinet_door", [bx0, 0, zDoorBack]);
  {
    const W = bx1 - bx0;
    const T = zDoorFront - zDoorBack;
    const [y0, y1] = doorY;
    const H = y1 - y0;
    const [wy0, wy1] = [y0 + 0.2, y1 - 0.13];
    const [wx0, wx1] = [0.075, W - 0.075];
    const skin = roundRect(W, H, 0.006, W / 2, (y0 + y1) / 2);
    rectPath(skin, wx1 - wx0, wy1 - wy0, 0.012, (wx0 + wx1) / 2, (wy0 + wy1) / 2, true);
    door.add(p.ss, plate(door, skin, T - 0.004, fine1 ? 0.002 : 0), [0, 0, 0.004]);
    door.add(p.ssDull, plate(door, skin, 0.004), [0, 0, 0]);
    // The window: a sealed double-glazed unit, two 4 mm panes held apart by a dark spacer that runs round the
    // reveal, so that the gap between them reads; their edges sit in the door's rebate.
    const [wcx, wcy, ww, wh] = [(wx0 + wx1) / 2, (wy0 + wy1) / 2, wx1 - wx0, wy1 - wy0];
    pane(door, fm.paneClear, fm.paneEdge, ww + 0.012, wh + 0.012, 0.004, [wcx, wcy, T - 0.006]);
    if (fine1) {
      pane(door, fm.paneClear, fm.paneEdge, ww + 0.012, wh + 0.012, 0.004, [wcx, wcy, 0.008]);
      const sp = roundRect(ww + 0.004, wh + 0.004, 0.013, wcx, wcy);
      rectPath(sp, ww - 0.018, wh - 0.018, 0.006, wcx, wcy, true);
      door.add(p.charcoal, plate(door, sp, T - 0.0186, fine2 ? 0.0006 : 0), [0, 0, 0.0103]);
    }
    // The window's trim ring.
    if (fine1) {
      const tr = roundRect(wx1 - wx0 + 0.024, wy1 - wy0 + 0.024, 0.016, (wx0 + wx1) / 2, (wy0 + wy1) / 2);
      rectPath(tr, wx1 - wx0 - 0.004, wy1 - wy0 - 0.004, 0.01, (wx0 + wx1) / 2, (wy0 + wy1) / 2, true);
      door.add(p.ssBright, plate(door, tr, 0.0016, fine2 ? 0.0005 : 0), [0, 0, T]);
    }
    // Magnetic gasket on the back, against the frame.
    if (fine1) rectLoop(door, p.gasket, 0.02, W - 0.02, y0 + 0.02, y1 - 0.02, -0.0025, 0.0045, 0.02, fine2 ? 8 : 5);
    // The latch handle: an edge pull on two standoffs; its latch tongue behind.
    const hx = W - 0.04;
    const hy = 0.98;
    for (const sg of [-1, 1]) door.add(p.ssBright, cyl(door, 0.012, 0.022, "z", 16), [hx, hy + sg * 0.16, T + 0.011]);
    tube(door, p.ssBright, [[hx, hy + 0.16, T + 0.016], [hx, hy - 0.16, T + 0.016]], 0.012, 0.01, 16);
    if (fine2) slab(door, p.ssBright, W - 0.03, W - 0.012, hy - 0.025, hy + 0.025, -0.006, 0.002, 0.002);
    // Barrel hinges: knuckles on the pivot, a leaf on the door's face.
    for (const y of [y0 + 0.12, (y0 + y1) / 2, y1 - 0.12]) {
      door.add(p.ssBright, cyl(door, 0.0075, 0.09, "y", 16), [0, y, 0]);
      if (fine1) slab(door, p.ssBright, 0.004, 0.034, y - 0.04, y + 0.04, T, T + 0.002, fine2 ? 0.001 : 0);
      if (fine2) for (const sy of [-1, 1]) screw(door, p.ssBright, [0.022, y + sy * 0.025, T + 0.002], 0.0024);
    }
  }

  return (b) => {
    layerGlass(b.group);
    setLamp(heatLamp, 0);
    setLamp(humLamp, 0);
    const doorJ = b.joints.get("cabinet_door");
    return {
      open: (f) => turnJoint(doorJ, "y", -Math.max(0, Math.min(1, f)) * 3.05),
      set: (name, v) => {
        if (name === "heat") {
          setLamp(heatLamp, v > 0.05 ? 1 : 0);
          turnJoint(b.joints.get(knobs[0]!), "z", -v * Math.PI * 1.5);
        } else if (name === "humidity") {
          setLamp(humLamp, v > 0.05 ? 1 : 0);
          turnJoint(b.joints.get(knobs[1]!), "z", -v * Math.PI * 1.5);
          vapour.set(v);
        }
      },
      tick: (t) => vapour.tick(t),
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Heavy-duty commercial microwave                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A heavy-duty commercial microwave: a stainless wrapper over a welded cavity with a ceramic floor (the waves come up
 * through it from a bottom feed, so there is no turntable), mica waveguide covers, a lamp behind a perforated
 * patch. The door's window is a perforated steel screen between two panes: its holes, far smaller than the 122 mm
 * wavelength, keep the microwaves in. A membrane keypad and display on the right, a grab handle on the door.
 */
const microwave: Builder = (s, e) => {
  const p = palette();
  const fm = mats();
  const { w, d, h } = e;
  const fine1 = s.fine(1);
  const fine2 = s.fine(2);
  const yb = 0.014;
  const zf = -0.04;
  const zb = -d + 0.006;
  const doorX: [number, number] = [0.006, 0.312];
  const doorY: [number, number] = [yb + 0.008, h - 0.01];
  const cav = { x0: 0.034, x1: 0.284, y0: 0.05, y1: h - 0.055, z0: -0.43, z1: zf - 0.002 };
  const lamp = liveLamp(0xffd9a0);

  for (const x of [0.03, w - 0.03]) for (const z of [zf - 0.03, zb + 0.03]) s.add(p.rubber, cyl(s, 0.014, yb, "y", 16, 0.011), [x, yb / 2, z]);
  // Wrapper: one sheet bent over the top and down both sides, its corners on a 6 mm radius.
  // (The flats stop short of the bends; each bend is a quarter cylinder, its axis along z.)
  const rb = 0.006;
  slab(s, p.ss, rb, w - rb, h - 0.0012, h, zb, zf);
  slab(s, p.ss, 0, 0.0012, yb, h - rb, zb, zf);
  slab(s, p.ss, w - 0.0012, w, yb, h - rb, zb, zf);
  for (const [x, a0] of [
    [rb, Math.PI],
    [w - rb, Math.PI / 2],
  ] as const) {
    s.add(p.ss, new THREE.CylinderGeometry(rb, rb, zf - zb, s.seg(12, 3), 1, true, a0, Math.PI / 2).rotateX(Math.PI / 2), [x, h - rb, (zf + zb) / 2]);
  }
  // The bottom pan sits between the sides and behind the front frame, not flush with their faces.
  slab(s, p.ssDull, 0.0012, w - 0.0012, yb, yb + 0.0012, zb, zf - 0.002);
  slab(s, p.ssDull, 0.0012, w - 0.0012, yb, h - 0.0012, zb, zb + 0.0012);
  // Vent louvres: sides near the back, and the back's exhaust grille.
  if (fine1) {
    for (const sg of [-1, 1]) louvres(s, [sg < 0 ? 0 : w, h * 0.55, zb + 0.1], [0, sg * Math.PI / 2, 0], 8, 0.12, 0.016);
    louvres(s, [w * 0.6, h * 0.5, zb], [0, Math.PI, 0], 9, 0.2, 0.016);
  }
  // The front frame round the cavity opening and the control strip's bezel.
  {
    const sh = roundRect(w - 0.001, h - yb - 0.001, 0.0058, w / 2, (h + yb) / 2);
    rectPath(sh, cav.x1 - cav.x0, cav.y1 - cav.y0, 0.008, (cav.x0 + cav.x1) / 2, (cav.y0 + cav.y1) / 2, true);
    s.add(p.ss, plate(s, sh, 0.002), [0, 0, zf - 0.002]);
    // Interlock slots the door's latch hooks enter.
    if (fine1) for (const y of [cav.y0 + 0.04, cav.y1 - 0.04]) slab(s, p.enamel, cav.x1 + 0.012, cav.x1 + 0.018, y - 0.012, y + 0.012, zf, zf + 0.0005);
  }
  // Controls: display, membrane keypad, the strip proud of the frame.
  const kx0 = doorX[1] + 0.006;
  const kx1 = w - 0.006;
  const kz = -0.014;
  slab(s, p.ss, kx0, kx1, doorY[0], doorY[1], zf, kz - 0.001, 0.003);
  slab(s, p.charcoal, kx0 + 0.006, kx1 - 0.006, doorY[0] + 0.012, doorY[1] - 0.012, kz - 0.0015, kz, 0.002);
  const kcx = (kx0 + kx1) / 2;
  slab(s, p.screen, kcx - 0.038, kcx + 0.038, h - 0.07, h - 0.035, kz - 0.0005, kz + 0.0008, 0.001);
  s.add(fm.mwDisplay, new THREE.PlaneGeometry(0.07, 0.028), [kcx, h - 0.0525, kz + 0.0012]);
  s.add(fm.keypad, new THREE.PlaneGeometry(0.076, 0.152), [kcx, h - 0.17, kz + 0.0006]);
  if (fine2) slab(s, p.ssBright, kcx - 0.04, kcx + 0.04, doorY[0] + 0.03, doorY[0] + 0.034, kz - 0.001, kz + 0.001, 0.001);
  dataPlate(s, [w * 0.5, h * 0.35, zb], [0, Math.PI, 0], 0.08, 0.045);
  if (fine1) tube(s, p.rubber, [[w - 0.06, 0.07, zb], [w - 0.06, 0.07, -d + 0.0005], [w - 0.06, yb + 0.004, -d + 0.0005]], 0.004, 0.02, 8);

  // The cavity: welded stainless, its side walls each a sheet with a punched patch (the lamp's grille at the upper
  // front left, the cooling air's inlet at the upper rear right), stamped with two shelf ribs; the back pressed with a
  // stiffening panel.
  if (fine1) {
    const t = 0.0012;
    const lampAt = { y: cav.y1 - 0.045, z: cav.z1 - 0.11, w: 0.06, h: 0.04 };
    const inletAt = { y: cav.y1 - 0.06, z: cav.z0 + 0.1, w: 0.08, h: 0.0533 };
    const sideWall = (x: number, hole: { y: number; z: number; w: number; h: number }) => {
      // In the plate's own frame x runs along -z and y up; it is turned to stand in the yz plane, t thick along +x.
      const sh = roundRect(cav.z1 - cav.z0, cav.y1 - cav.y0, 0, -(cav.z0 + cav.z1) / 2, (cav.y0 + cav.y1) / 2);
      sh.holes.push(rectPath(new THREE.Path(), hole.w, hole.h, fine2 ? 0.004 : 0, -hole.z, hole.y, true) as THREE.Path);
      const g = plate(s, sh, t);
      g.rotateY(Math.PI / 2);
      s.add(p.ssDull, g, [x, 0, 0]);
    };
    sideWall(cav.x0 - t, lampAt);
    sideWall(cav.x1, inletAt);
    slab(s, p.ssDull, cav.x0, cav.x1, cav.y1, cav.y1 + t, cav.z0, cav.z1);
    slab(s, p.ssDull, cav.x0, cav.x1, cav.y0, cav.y1, cav.z0 - t, cav.z0);
    // The punched patches: the lamp glowing behind its grille in a dark reflector box; darkness behind the inlet.
    s.add(fm.grille, new THREE.PlaneGeometry(lampAt.w, lampAt.h).rotateY(Math.PI / 2), [cav.x0 - t / 2, lampAt.y, lampAt.z]);
    s.add(lamp, new THREE.PlaneGeometry(lampAt.w, lampAt.h).rotateY(Math.PI / 2), [cav.x0 - 0.007, lampAt.y, lampAt.z]);
    {
      const [ly0, ly1, lz0, lz1] = [lampAt.y - lampAt.h / 2, lampAt.y + lampAt.h / 2, lampAt.z - lampAt.w / 2, lampAt.z + lampAt.w / 2];
      slab(s, p.ssBright, cav.x0 - 0.0075, cav.x0 - t, ly0 - 0.001, ly1 + 0.001, lz0 - 0.001, lz0);
      slab(s, p.ssBright, cav.x0 - 0.0075, cav.x0 - t, ly0 - 0.001, ly1 + 0.001, lz1, lz1 + 0.001);
      slab(s, p.ssBright, cav.x0 - 0.0075, cav.x0 - t, ly0 - 0.001, ly0, lz0, lz1);
      slab(s, p.ssBright, cav.x0 - 0.0075, cav.x0 - t, ly1, ly1 + 0.001, lz0, lz1);
    }
    s.add(fm.grille, new THREE.PlaneGeometry(inletAt.w, inletAt.h).rotateY(-Math.PI / 2), [cav.x1 + t / 2, inletAt.y, inletAt.z]);
    s.add(p.enamel, new THREE.PlaneGeometry(inletAt.w, inletAt.h).rotateY(-Math.PI / 2), [cav.x1 + 0.008, inletAt.y, inletAt.z]);
    // Shelf ribs stamped into both side walls, standing 3 mm proud.
    for (const x of [cav.x0, cav.x1])
      for (const y of [cav.y0 + 0.075, cav.y0 + 0.145]) s.add(p.ssDull, fine2 ? rbox(s, 0.006, 0.011, cav.z1 - cav.z0 - 0.07, 0.0029) : box(0.004, 0.009, cav.z1 - cav.z0 - 0.07), [x, y, (cav.z0 + cav.z1) / 2 - 0.01]);
    // The back's pressed panel: a shallow raised field with soft edges.
    if (fine2) s.add(p.ssDull, rbox(s, cav.x1 - cav.x0 - 0.06, cav.y1 - cav.y0 - 0.07, 0.006, 0.0028), [(cav.x0 + cav.x1) / 2, (cav.y0 + cav.y1) / 2, cav.z0 - 0.0015]);
    // The ceramic floor over the bottom feed, sealed into its frame.
    slab(s, p.ceramic, cav.x0, cav.x1, cav.y0 - 0.006, cav.y0, cav.z0, cav.z1, fine2 ? 0.001 : 0);
    if (fine2) for (const x of [cav.x0 + 0.0015, cav.x1 - 0.0015]) s.add(p.gasket, cyl(s, 0.0018, cav.z1 - cav.z0, "z", 8), [x, cav.y0, (cav.z0 + cav.z1) / 2]);
    // Mica waveguide cover on the right wall, the stirrer cover above, each held by plastic rivets.
    const wg: V3 = [cav.x1 - 0.0008, cav.y1 - 0.042, (cav.z0 + cav.z1) / 2 + 0.02];
    s.add(fm.mica, new THREE.PlaneGeometry(0.12, 0.072).rotateY(-Math.PI / 2), wg);
    s.add(fm.mica, new THREE.PlaneGeometry(0.16, 0.2).rotateX(Math.PI / 2), [(cav.x0 + cav.x1) / 2, cav.y1 - 0.0008, (cav.z0 + cav.z1) / 2]);
    if (fine2) {
      for (const sz of [-1, 1]) s.add(p.hdpe, lathe(s, [[0, 0], [0.004, 0], [0.0036, 0.0012], [0, 0.0016]], 12), [wg[0], wg[1], wg[2] + sz * 0.052], [0, 0, Math.PI / 2]);
      for (const sx of [-1, 1]) for (const sz of [-1, 1]) s.add(p.hdpe, lathe(s, [[0, 0], [0.004, 0], [0.0036, 0.0012], [0, 0.0016]], 12), [(cav.x0 + cav.x1) / 2 + sx * 0.07, cav.y1 - 0.0008, (cav.z0 + cav.z1) / 2 + sz * 0.09], [Math.PI, 0, 0]);
    }
  } else slab(s, p.charcoal, cav.x0, cav.x1, cav.y0, cav.y1, cav.z0, cav.z1 - 0.003);

  // The door.
  const door = s.joint("door", [doorX[0], 0, zf]);
  {
    const W = doorX[1] - doorX[0];
    const T = 0.026;
    const [y0, y1] = doorY;
    const [wx0, wx1, wy0, wy1] = [0.032, W - 0.048, y0 + 0.042, y1 - 0.042];
    const fr = roundRect(W, y1 - y0, 0.004, W / 2, (y0 + y1) / 2);
    rectPath(fr, wx1 - wx0, wy1 - wy0, 0.006, (wx0 + wx1) / 2, (wy0 + wy1) / 2, true);
    door.add(p.ss, plate(door, fr, T - 0.006, fine1 ? 0.0015 : 0), [0, 0, 0.006]);
    door.add(p.glass, box(wx1 - wx0 + 0.008, wy1 - wy0 + 0.008, 0.003), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, T - 0.004]);
    door.add(fm.screenMesh, new THREE.PlaneGeometry(wx1 - wx0 + 0.004, wy1 - wy0 + 0.004), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, T - 0.0085]);
    if (fine1) door.add(p.glass, box(wx1 - wx0 + 0.008, wy1 - wy0 + 0.008, 0.002), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, 0.009]);
    // The choke: a dark seal band round the back, and the interlock hooks.
    door.add(p.charcoal, plate(door, fr, 0.006), [0, 0, 0]);
    if (fine2) for (const y of [cav.y0 + 0.04, cav.y1 - 0.04]) slab(door, p.ssBright, W - 0.022 + 0.006, W - 0.016 + 0.006, y - 0.01, y + 0.01, -0.014, 0, 0.001);
    // The grab handle.
    const hx = W - 0.022;
    for (const y of [y0 + 0.06, y1 - 0.06]) door.add(p.ssBright, cyl(door, 0.009, 0.016, "z", 14), [hx, y, T + 0.008]);
    tube(door, p.ssBright, [[hx, y1 - 0.06, T + 0.01], [hx, y1 - 0.08, T + 0.02], [hx, y0 + 0.08, T + 0.02], [hx, y0 + 0.06, T + 0.01]], 0.0085, 0.02, 14);
  }

  return (b) => {
    setLamp(lamp, 0);
    const doorJ = b.joints.get("door");
    return {
      open: (f) => turnJoint(doorJ, "y", -Math.max(0, Math.min(1, f)) * 1.745),
      set: (name, v) => {
        if (name === "power") setLamp(lamp, v > 0.02 ? 1 : 0);
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Conveyor toaster                                                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * The conveyor belt's chain: cross rods on side links, carried over two sprockets. Each rod is an instance whose place
 * round the loop the tick advances; the sprockets turn with it.
 */
class Chain {
  readonly mesh: THREE.InstancedMesh;
  private L: number;
  private off = 0;
  private m = new THREE.Matrix4();
  constructor(
    private cx: number,
    private zF: number,
    private zR: number,
    private yC: number,
    private R: number,
    private pitch: number,
    width: number,
    mat: Mat,
  ) {
    const run = zF - zR;
    this.L = 2 * run + 2 * Math.PI * R;
    const n = Math.floor(this.L / pitch);
    this.pitch = this.L / n;
    const parts: THREE.BufferGeometry[] = [new THREE.CylinderGeometry(0.0011, 0.0011, width, 5, 1, true).rotateZ(Math.PI / 2)];
    for (const sx of [-1, 1]) parts.push(new THREE.BoxGeometry(0.0016, 0.0042, this.pitch * 1.15).translate(sx * (width / 2 + 0.0008), 0, -this.pitch * 0.5));
    const g = mergeGeoms(parts);
    this.mesh = new THREE.InstancedMesh(g, mat, n);
    this.mesh.castShadow = true;
    this.mesh.receiveShadow = true;
    this.place();
  }
  get tris(): number {
    return ((this.mesh.geometry.index?.count ?? 0) / 3) * this.mesh.count;
  }
  private at(u: number, m: THREE.Matrix4): void {
    const { zF, zR, yC, R } = this;
    const run = zF - zR;
    let y: number;
    let z: number;
    let a: number;
    if (u < run) [z, y, a] = [zF - u, yC + R, 0];
    else if (u < run + Math.PI * R) {
      const f = (u - run) / R;
      [z, y, a] = [zR - R * Math.sin(f), yC + R * Math.cos(f), f];
    } else if (u < 2 * run + Math.PI * R) [z, y, a] = [zR + (u - run - Math.PI * R), yC - R, Math.PI];
    else {
      const f = (u - 2 * run - Math.PI * R) / R;
      [z, y, a] = [zF + R * Math.sin(f), yC - R * Math.cos(f), Math.PI + f];
    }
    m.makeRotationX(a).setPosition(this.cx, y, z);
  }
  private place(): void {
    const { m } = this;
    for (let i = 0; i < this.mesh.count; i++) {
      this.at((((i * this.pitch + this.off) % this.L) + this.L) % this.L, m);
      this.mesh.setMatrixAt(i, m);
    }
    this.mesh.instanceMatrix.needsUpdate = true;
    // The links only move round the same loop, so its bounds, taken once, hold however far the belt has run.
    if (!this.mesh.boundingSphere) this.mesh.computeBoundingSphere();
  }
  advance(dist: number): void {
    this.off = (this.off + dist) % this.L;
    this.place();
  }
  get travel(): number {
    return this.off;
  }
}
function mergeGeoms(parts: THREE.BufferGeometry[]): THREE.BufferGeometry {
  const pos: number[] = [];
  const nrm: number[] = [];
  const idx: number[] = [];
  for (const g0 of parts) {
    const g = g0;
    const o = pos.length / 3;
    const pa = g.attributes.position as THREE.BufferAttribute;
    const na = g.attributes.normal as THREE.BufferAttribute;
    for (let i = 0; i < pa.count; i++) {
      pos.push(pa.getX(i), pa.getY(i), pa.getZ(i));
      nrm.push(na.getX(i), na.getY(i), na.getZ(i));
    }
    const ia = g.index!;
    for (let i = 0; i < ia.count; i++) idx.push(ia.getX(i) + o);
  }
  const out = new THREE.BufferGeometry();
  out.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  out.setAttribute("normal", new THREE.Float32BufferAttribute(nrm, 3));
  out.setIndex(idx);
  return out;
}

/**
 * A countertop conveyor toaster: bread laid on the wire belt at the front rides through between quartz elements
 * (two banks above, one between the belt's runs), drops off the back onto a chute that returns it to the front tray
 * under the feed. A crumb tray slides out underneath; the controls are a heat knob, the belt-speed knob and a rocker.
 */
const conveyorToaster: Builder = (s, e) => {
  const p = palette();
  const fm = breadMats();
  const { w, d, h } = e;
  const fine1 = s.fine(1);
  const fine2 = s.fine(2);
  const glow = glowMaterial(0x30201a);
  const yb = 0.026;
  const zf = -0.07;
  const zb = -d + 0.004;
  const [x0, x1] = [0.004, w - 0.004];
  const belt = { x0: 0.054, x1: w - 0.054, yC: 0.181, R: 0.014, zF: -0.016, zR: -0.425 };
  const yTop = belt.yC + belt.R;
  const cx = w / 2;

  // Legs.
  for (const x of [x0 + 0.03, x1 - 0.03]) for (const z of [zf - 0.03, zb + 0.03]) s.add(p.ssBright, lathe(s, [[0, 0], [0.012, 0], [0.013, 0.004], [0.009, 0.01], [0.008, yb], [0, yb]], 14), [x, 0, z]);
  // Housing: a wrapper over top and sides, back and bottom pans.
  slab(s, p.ss, x0, x1, h - 0.0012, h, zb, zf, 0.0006);
  for (const x of [x0, x1 - 0.0012]) slab(s, p.ss, x, x + 0.0012, yb, h, zb, zf, 0.0006);
  // The pans sit between the side sheets and under the top, not flush with their outer faces.
  slab(s, p.ssDull, x0 + 0.0012, x1 - 0.0012, yb, yb + 0.0012, zb, zf - 0.0015);
  slab(s, p.ssDull, x0 + 0.0012, x1 - 0.0012, yb, h - 0.0012, zb, zb + 0.0012);
  if (fine1) {
    for (const sg of [-1, 1]) louvres(s, [sg < 0 ? x0 : x1, h - 0.07, zb + 0.12], [0, sg * Math.PI / 2, 0], 4, 0.14, 0.018);
    louvres(s, [cx, h - 0.0001 + 0.0001, zb + 0.1], [-Math.PI / 2, 0, 0], 5, 0.22, 0.02);
  }
  // Front panel: the feed opening over the belt, the discharge below it, the controls above.
  const feed = { y0: yTop - 0.006, y1: yTop + 0.07 };
  const dis = { y0: 0.048, y1: 0.125 };
  {
    const sh = roundRect(x1 - x0, h - yb, 0.006, cx, (h + yb) / 2);
    rectPath(sh, belt.x1 - belt.x0 + 0.02, feed.y1 - feed.y0, 0.006, cx, (feed.y0 + feed.y1) / 2, true);
    rectPath(sh, belt.x1 - belt.x0 + 0.02, dis.y1 - dis.y0, 0.006, cx, (dis.y0 + dis.y1) / 2, true);
    rectPath(sh, belt.x1 - belt.x0 + 0.004, 0.016, 0.002, cx, yb + 0.011, true);
    s.add(p.ss, plate(s, sh, 0.0015, fine1 ? 0.0006 : 0), [0, 0, zf - 0.0015]);
    // A heat deflector lip over the feed.
    s.add(p.ss, box(belt.x1 - belt.x0 + 0.04, 0.0012, 0.022), [cx, feed.y1 + 0.004, zf + 0.008], [-0.25, 0, 0]);
  }
  // Controls.
  const kY = (feed.y1 + h) / 2 + 0.005;
  const heatKnob = s.joint("heat_knob", [cx - 0.09, kY, zf]);
  dial(heatKnob, p.bakelite, p.ssBright, [0, 0, 0], 0.018);
  const speedKnob = s.joint("speed_knob", [cx + 0.03, kY, zf]);
  dial(speedKnob, p.bakelite, p.ssBright, [0, 0, 0], 0.018);
  if (fine2) for (const x of [cx - 0.09, cx + 0.03]) s.add(p.enamel, cyl(s, 0.027, 0.0006, "z", 32), [x, kY, zf + 0.0007]);
  slab(s, p.plastic, cx + 0.1, cx + 0.125, kY - 0.016, kY + 0.016, zf, zf + 0.008, 0.002);
  const pilot = liveLamp(0xffa040);
  s.add(pilot, cyl(s, 0.004, 0.004, "z", 12), [cx + 0.1125, kY + 0.028, zf + 0.002]);
  // The tunnel: heat-tinted walls round the belt, reflectors over the top bank.
  if (fine1) {
    slab(s, p.ssHeat, belt.x0 - 0.008, belt.x1 + 0.008, yTop + 0.085, yTop + 0.087, belt.zR - 0.02, zf - 0.002);
    for (const x of [belt.x0 - 0.009, belt.x1 + 0.009]) slab(s, p.ssHeat, x - 0.001, x + 0.001, 0.13, yTop + 0.087, belt.zR - 0.02, zf - 0.002);
    s.add(p.ssBright, box(belt.x1 - belt.x0, 0.001, belt.zF - belt.zR - 0.08), [cx, yTop + 0.075, (belt.zR + zf) / 2 - 0.01]);
  } else slab(s, p.charcoal, belt.x0, belt.x1, feed.y0, feed.y1, belt.zR, zf - 0.003);
  // Quartz elements: a coiled filament in a quartz tube, ceramic end caps in the side walls.
  const quartz = (y: number, z: number) => {
    const L = belt.x1 - belt.x0 + 0.012;
    s.add(glow, cyl(s, 0.0022, L, "x", 8), [cx, y, z]);
    if (fine1) s.add(p.glass, cyl(s, 0.0055, L, "x", 14), [cx, y, z]);
    if (fine1) for (const sg of [-1, 1]) s.add(p.ceramic, cyl(s, 0.0065, 0.008, "x", 12), [cx + sg * (L / 2 + 0.002), y, z]);
  };
  for (const z of [-0.12, -0.2, -0.28, -0.36]) quartz(yTop + 0.05, z);
  for (const z of [-0.16, -0.3]) quartz(belt.yC, z);
  // The belt: sprockets on shafts at each end, the chain over them, side guides over the loading end.
  const sprockets: string[] = [];
  for (const [i, z] of [belt.zF, belt.zR].entries()) {
    const name = i ? "sprocket_rear" : "sprocket_front";
    sprockets.push(name);
    const j = s.joint(name, [cx, belt.yC, z]);
    j.add(p.ssBright, cyl(j, 0.004, belt.x1 - belt.x0 + 0.03, "x", 10), [0, 0, 0]);
    for (const sx of [belt.x0 + 0.004, belt.x1 - 0.004]) {
      j.add(p.ssDull, cyl(j, belt.R - 0.0015, 0.003, "x", 16), [sx - cx, 0, 0]);
      if (fine2)
        for (let k = 0; k < 10; k++) {
          const a = (k / 10) * Math.PI * 2;
          j.add(p.ssDull, box(0.003, 0.003, 0.004), [sx - cx, Math.cos(a) * (belt.R - 0.0005), Math.sin(a) * (belt.R - 0.0005)], [a, 0, 0]);
        }
    }
  }
  // The loading end: the belt runs out over the front on two side plates bolted to the front panel, each carrying
  // the front shaft's bearing and standing a little above the belt as a guide; a crumb sheet under the belt between
  // them. Nothing here moves but the belt.
  {
    const t = 0.0015;
    const [py0, py1] = [belt.yC - 0.02, yTop + 0.012];
    const [pz0, pz1] = [zf - 0.002, -0.0008];
    for (const [sg, xo] of [
      [-1, belt.x0 - 0.007],
      [1, belt.x1 + 0.007 + t],
    ] as const) {
      if (!fine1) {
        slab(s, p.ss, xo - t, xo, py0, py1, pz0, pz1);
        continue;
      }
      const g = plate(s, roundRect(pz1 - pz0, py1 - py0, 0.011, (pz0 + pz1) / 2, (py0 + py1) / 2), t, fine2 ? 0.0004 : 0);
      g.rotateY(-Math.PI / 2);
      s.add(p.ss, g, [xo, 0, 0]);
      const xb = xo - (sg < 0 ? t : 0) + sg * 0.003;
      s.add(p.ssBright, cyl(s, 0.0075, 0.006, "x", 18), [xb, belt.yC, belt.zF]);
      if (fine2) {
        s.add(p.ssBright, cyl(s, 0.0042, 0.008, "x", 6), [xb + sg * 0.004, belt.yC, belt.zF]);
        for (const y of [py0 + 0.01, py1 - 0.01]) screw(s, p.ssBright, [xo - (sg < 0 ? t : 0), y, zf + 0.012], 0.0022, [0, (sg * Math.PI) / 2, 0]);
      }
    }
    slab(s, p.ss, belt.x0 - 0.007, belt.x1 + 0.007, py0 + 0.001, py0 + 0.002, zf, -0.004, 0.0003);
  }
  let chain: Chain | null = null;
  if (fine1) {
    chain = new Chain(cx, belt.zF, belt.zR, belt.yC, belt.R + 0.0012, fine2 ? 0.0127 : 0.0254, belt.x1 - belt.x0, p.chrome);
    s.extra(chain.mesh);
    s.tris += chain.tris;
  } else s.add(p.chrome, box(belt.x1 - belt.x0, 0.002, belt.zF - belt.zR), [cx, yTop, (belt.zF + belt.zR) / 2]);
  // Bread: slices of a sandwich loaf, their crowns toward the back; two waiting on the belt, half fed into the
  // opening, one well inside, and one toasted coming down onto the front tray. Toast curls a little as it dries.
  const slice = (at: V3, rot: V3, done: boolean) => {
    const g = flatFood(s, breadShape(), 0.0125, fine2 ? 0.0022 : 0.0016);
    if (done) {
      const pa = g.attributes.position as THREE.BufferAttribute;
      for (let i = 0; i < pa.count; i++) pa.setY(i, pa.getY(i) + 0.0016 * ((pa.getX(i) / 0.055) ** 2 + (pa.getZ(i) / 0.06) ** 2));
    }
    s.add(done ? fm.breadToast : fm.breadRaw, g, at, rot);
  };
  const ty = yTop + 0.0012 + 0.00625;
  if (fine1) {
    slice([cx - 0.0655, ty, -0.072], [0, 0.03, 0], false);
    slice([cx + 0.0655, ty, -0.088], [0, -0.025, 0], false);
    slice([cx + 0.05, ty, -0.27], [0, 0.04, 0], true);
  }
  // The return chute from the back, and the front toast tray: a drawn stainless tray with a rolled rim, standing out
  // of the discharge and falling slightly toward the front.
  if (fine1) {
    const g = box(belt.x1 - belt.x0 + 0.01, 0.0012, Math.hypot(belt.zR - zf, 0.1));
    s.add(p.ss, g, [cx, (0.152 + dis.y0 + 0.008) / 2, (belt.zR + zf) / 2], [Math.atan2(0.152 - dis.y0 - 0.008, zf - belt.zR), 0, 0]);
  }
  {
    const tw = belt.x1 - belt.x0 + 0.012;
    const td = 0.13;
    const dep = 0.018;
    const tz = zf - 0.08 + td / 2;
    const prof: [number, number][] = fine2
      ? [
          [0.009, -dep],
          [0.003, -dep + 0.004],
          [0.0008, -0.004],
          [-0.0016, -0.0014],
          [0, 0],
          [0.0016, -0.0009],
          [0.0022, -0.004],
          [0.0042, -dep + 0.0042],
          [0.0102, -dep + 0.0012],
        ]
      : [
          [0.008, -dep],
          [0, -0.002],
          [0, 0],
          [0.0022, -0.002],
          [0.0092, -dep + 0.0012],
        ];
    s.at([cx, dis.y0 + dep + 0.002, tz], [0.07, 0, 0], () => {
      if (fine1) {
        s.add(p.ss, vessel(s, tw, td, 0.01, prof));
        slice([-0.03, -dep + 0.0012 + 0.0065, -0.002], [0, 0.04, 0], true);
      } else s.add(p.ss, box(tw, dep, td), [0, -dep / 2, 0]);
    });
  }
  // The crumb tray: a shallow drawn pan the full length of the base, its front plate closing a slot in the front
  // panel, a pull tab on it. open() draws it out.
  const trayRest: V3 = [0, 0, 0];
  const tray = s.joint("crumb_tray", trayRest);
  {
    const tw = belt.x1 - belt.x0;
    const z0 = zb + 0.02;
    const z1 = zf - 0.0005;
    if (fine1)
      tray.add(p.ssDull, vessel(tray, tw, z1 - z0, 0.006, [[0.004, -0.0135], [0, -0.0005], [0, 0], [0.001, 0], [0.005, -0.0125]]), [cx, yb + 0.0165, (z0 + z1) / 2]);
    tray.add(p.ss, rbox(tray, tw + 0.012, 0.021, 0.0015, 0.0006), [cx, yb + 0.011, zf + 0.00075]);
    tray.add(p.ssBright, rbox(tray, 0.06, 0.009, 0.007, 0.002), [cx, yb + 0.011, zf + 0.0015 + 0.0035]);
    if (fine2) tray.add(p.enamel, rbox(tray, 0.044, 0.003, 0.0012, 0.0005), [cx, yb + 0.011, zf + 0.0015 + 0.0072]);
  }
  dataPlate(s, [x1, 0.1, zf - 0.08], [0, Math.PI / 2, 0], 0.07, 0.04);
  if (fine1) tube(s, p.rubber, [[cx + 0.1, 0.08, zb], [cx + 0.1, 0.08, -d + 0.0005], [cx + 0.1, yb + 0.004, -d + 0.0005]], 0.004, 0.02, 8);

  return (b) => {
    let speed = 0.5;
    setGlow(glow, 0);
    const sj = sprockets.map((n) => b.joints.get(n));
    turnJoint(b.joints.get("speed_knob"), "z", -speed * Math.PI * 1.5);
    setLamp(pilot, 0);
    const trayJ = b.joints.get("crumb_tray");
    return {
      glows: [glow],
      // The only thing on a conveyor toaster that opens is its crumb tray.
      open: (f) => slideJoint(trayJ, "z", Math.max(0, Math.min(1, f)) * 0.14, trayRest),
      set: (name, v) => {
        if (name === "heat") {
          setGlow(glow, v);
          setLamp(pilot, v > 0.05 ? 1 : 0);
          turnJoint(b.joints.get("heat_knob"), "z", -v * Math.PI * 1.5);
        } else if (name === "speed") {
          speed = Math.max(0, Math.min(1, v));
          turnJoint(b.joints.get("speed_knob"), "z", -speed * Math.PI * 1.5);
        }
      },
      tick: (_t, dt) => {
        if (!chain) return;
        chain.advance(Math.max(0, Math.min(0.1, dt)) * (0.004 + 0.03 * speed));
        for (const j of sj) turnJoint(j, "x", chain.travel / (belt.R + 0.0012));
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* High-velocity air fryer                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A commercial countertop air fryer: an impeller over a spiral element blows air down hard through a perforated basket
 * in a drawer that runs on telescopic slides; the exhaust leaves through a catalytic filter under the louvred top, so
 * it needs no hood. Touch keys and a display above the drawer.
 */
const airFryer: Builder = (s, e) => {
  const p = palette();
  const fm = mats();
  const { w, d, h } = e;
  const fine1 = s.fine(1);
  const fine2 = s.fine(2);
  const glow = glowMaterial(0x2a2420);
  const yb = 0.016;
  const zf = -0.06;
  const zb = -d + 0.004;
  const cx = w / 2;
  const cav = { x0: 0.045, x1: w - 0.045, y0: 0.05, y1: 0.33, z0: -0.52, z1: zf };
  const op = { x0: 0.045, x1: w - 0.045, y0: 0.05, y1: 0.29 };

  for (const x of [0.04, w - 0.04]) for (const z of [zf - 0.04, zb + 0.05]) s.add(p.rubber, cyl(s, 0.013, yb, "y", 16, 0.011), [x, yb / 2, z]);
  // The case: a wrapper with rounded front corners, back and bottom pans, the front frame round the drawer.
  slab(s, p.ss, 0, w, h - 0.0012, h, zb, zf, 0.0006);
  for (const x of [0, w - 0.0012]) slab(s, p.ss, x, x + 0.0012, yb, h, zb, zf, 0.0006);
  slab(s, p.ssDull, 0.0012, w - 0.0012, yb, yb + 0.0012, zb, zf - 0.002);
  slab(s, p.ssDull, 0.0012, w - 0.0012, yb, h - 0.0012, zb, zb + 0.0012);
  {
    const sh = roundRect(w - 0.002, h - yb - 0.002, 0.008, w / 2, (h + yb) / 2);
    rectPath(sh, op.x1 - op.x0, op.y1 - op.y0, 0.01, cx, (op.y0 + op.y1) / 2, true);
    s.add(p.ss, plate(s, sh, 0.002, fine1 ? 0.0008 : 0), [0, 0, zf - 0.002]);
  }
  if (fine1) {
    louvres(s, [cx, h, zb + 0.12], [-Math.PI / 2, 0, 0], 6, 0.3, 0.02);
    for (const sg of [-1, 1]) louvres(s, [sg < 0 ? 0 : w, 0.12, zb + 0.12], [0, sg * Math.PI / 2, 0], 4, 0.14, 0.018);
  }
  // Controls: a black glass strip with the display and the touch keys, lit.
  const pz = -0.036;
  const py = (op.y1 + h) / 2 + 0.004;
  slab(s, p.ss, 0.03, w - 0.03, op.y1 + 0.012, h - 0.012, zf, pz - 0.0015, 0.004);
  slab(s, p.bakelite, 0.036, w - 0.036, op.y1 + 0.018, h - 0.018, pz - 0.002, pz, 0.003);
  slab(s, p.screen, cx - 0.06, cx + 0.06, py - 0.022, py + 0.022, pz - 0.001, pz + 0.0008, 0.002);
  s.add(fm.fryDisplay, new THREE.PlaneGeometry(0.112, 0.04), [cx, py, pz + 0.0014]);
  for (const [i, x] of [0.075, 0.12, w - 0.12, w - 0.075].entries()) {
    s.add(i === 3 ? p.ledGreen : p.ledWhite, torus(s, 0.011, 0.0011, "z", Math.PI * 2, 32, 4), [x, py, pz + 0.0004]);
    if (fine2) s.add(i === 3 ? p.ledGreen : p.ledWhite, box(0.008, 0.0016, 0.001), [x, py, pz + 0.0004], [0, 0, i * 0.6]);
  }
  dataPlate(s, [w * 0.35, h * 0.5, zb], [0, Math.PI, 0], 0.08, 0.045);
  if (fine1) tube(s, p.rubber, [[w * 0.7, 0.07, zb], [w * 0.7, 0.07, -d + 0.0005], [w * 0.7, yb + 0.004, -d + 0.0005]], 0.004, 0.02, 8);

  // The cavity: heat-tinted, the spiral element and the impeller over the basket, the slides' fixed members.
  if (fine1) {
    slab(s, p.ssHeat, cav.x0 - 0.0012, cav.x0, cav.y0, cav.y1, cav.z0, cav.z1 - 0.002);
    slab(s, p.ssHeat, cav.x1, cav.x1 + 0.0012, cav.y0, cav.y1, cav.z0, cav.z1 - 0.002);
    slab(s, p.ssHeat, cav.x0, cav.x1, cav.y1, cav.y1 + 0.0012, cav.z0, cav.z1 - 0.002);
    slab(s, p.ssHeat, cav.x0, cav.x1, cav.y0 - 0.0012, cav.y0, cav.z0, cav.z1 - 0.002);
    slab(s, p.ssHeat, cav.x0, cav.x1, cav.y0, cav.y1, cav.z0 - 0.0012, cav.z0);
    for (const x of [cav.x0 + 0.003, cav.x1 - 0.003]) slab(s, p.ssBright, x - 0.003, x + 0.003, 0.17, 0.2, cav.z0 + 0.02, cav.z1 - 0.01, fine2 ? 0.001 : 0);
  } else slab(s, p.charcoal, op.x0, op.x1, op.y0, op.y1, cav.z0, zf - 0.003);
  const ez = (cav.z0 + cav.z1) / 2 - 0.02;
  if (fine1) {
    for (const R of [0.05, 0.085, 0.12]) s.add(glow, torus(s, R, 0.004, "y", Math.PI * 1.85, 48, 6), [cx, 0.3, ez], [0, R * 9, 0]);
    s.add(glow, box(0.004, 0.006, 0.08), [cx + 0.07, 0.3, ez - 0.02], [0, 0.4, 0]);
    // The guard under the impeller.
    for (let i = -3; i <= 3; i++) s.add(p.ssDull, rod(0.0016, 0.3, "x", 5), [cx, 0.312, ez + i * 0.04]);
  }
  const fan = s.joint("fan", [cx, 0.322, ez]);
  if (fine1) {
    fan.add(p.ssDull, cyl(fan, 0.13, 0.002, "y", 40), [0, 0.006, 0]);
    for (let i = 0; i < 12; i++) fan.add(p.ssDull, fanBlade(fan, 0.05, 0.128, 0.5, 0.0014, 0.012).rotateX(Math.PI / 2), [0, 0.006, 0], [0, (i / 12) * Math.PI * 2, 0]);
    fan.add(p.ssBright, cyl(fan, 0.018, 0.012, "y", 16), [0, 0, 0]);
  }

  // The drawer: front, handle, pan, basket of fries, the slides' running members.
  const rest: V3 = [0, 0, 0];
  const dr = s.joint("basket_slide", rest);
  {
    const fz0 = zf + 0.001;
    const fz1 = -0.034;
    slab(dr, p.ss, op.x0 - 0.005, op.x1 + 0.005, op.y0 - 0.005, op.y1 + 0.005, fz0, fz1, 0.006);
    if (fine1) slab(dr, p.gasket, op.x0 + 0.004, op.x1 - 0.004, op.y0 + 0.004, op.y1 - 0.004, zf - 0.004, fz0, 0.002);
    const hy = 0.15;
    for (const x of [cx - 0.11, cx + 0.11]) dr.add(p.ssBright, cyl(dr, 0.011, 0.024, "z", 16), [x, hy, fz1 + 0.012]);
    tube(dr, p.bakelite, [[cx - 0.11, hy, fz1 + 0.018], [cx - 0.09, hy, fz1 + 0.026], [cx + 0.09, hy, fz1 + 0.026], [cx + 0.11, hy, fz1 + 0.018]], 0.0095, 0.02, 14);
    // The drawer pan, its walls rising to the opening's top.
    const pw = op.x1 - op.x0 - 0.018;
    const pd = 0.42;
    const pzc = zf - 0.004 - pd / 2;
    if (fine1) {
      dr.add(p.ssDull, vessel(dr, pw, pd, 0.02, [[0.006, 0.058], [0.002, 0.062], [0, 0.28], [0.0012, 0.2812], [0.0024, 0.28], [0.0024, 0.064], [0.007, 0.0598]]), [cx, 0, pzc]);
      for (const sg of [-1, 1]) slab(dr, p.ssBright, cx + sg * (pw / 2) - 0.002, cx + sg * (pw / 2) + 0.002, 0.175, 0.195, pzc - pd / 2 + 0.02, pzc + pd / 2 - 0.01, fine2 ? 0.0008 : 0);
      // The basket: perforated walls and floor in a rolled wire frame, a fold-down handle bracket at its front.
      const bw = pw - 0.03;
      const bd = pd - 0.05;
      const by0 = 0.085;
      const by1 = 0.26;
      dr.add(fm.perf, new THREE.PlaneGeometry(bw, bd).rotateX(-Math.PI / 2), [cx, by0, pzc]);
      for (const sg of [-1, 1]) {
        dr.add(fm.perf, new THREE.PlaneGeometry(bd, by1 - by0).rotateY(Math.PI / 2), [cx + sg * (bw / 2), (by0 + by1) / 2, pzc]);
        dr.add(fm.perf, new THREE.PlaneGeometry(bw, by1 - by0), [cx, (by0 + by1) / 2, pzc + sg * (bd / 2)]);
      }
      // Rolled wire rims round its floor and its top (rectLoop builds in xy; the frame turns it flat).
      for (const y of [by0, by1])
        dr.at([cx, y, pzc], [Math.PI / 2, 0, 0], () => rectLoop(dr, p.ss, -bw / 2, bw / 2, -bd / 2, bd / 2, 0, 0.0025, 0.012, 6));
      // Fries: a heap close up, a lumpy fill further off.
      const fr = rng(91);
      const hy0 = by0 + 0.072;
      dr.add(fm.fries, foodTop(dr, bw - 0.01, bd - 0.01, 0.01, 0.016, 0.07, 4), [cx, hy0, pzc]);
      if (fine2) {
        // The loose fries over the heap, dropped one by one so that each comes to rest on the heap or on the fries
        // already there: a height map in 1 cm cells, starting from the heap's own surface.
        const k = (Math.PI * 2) / 0.07;
        const [hw, hd] = [(bw - 0.01) / 2, (bd - 0.01) / 2];
        const heapAt = (x: number, z: number) => {
          const edge = Math.min(hw - Math.abs(x), hd - Math.abs(z));
          const lump = Math.sin(x * k + 4) * Math.cos(z * k * 0.8 + 6.8) + 0.55 * Math.sin((x * 0.6 - z) * k * 1.9 + 12.4) + 0.3 * Math.cos((x + z * 0.7) * k * 3.3 + 2.8);
          return hy0 + 0.016 * (0.5 + 0.32 * lump) - 0.016 * 0.6 * Math.exp(-Math.max(0, edge) / (0.07 * 0.25));
        };
        const hc = 0.01;
        const [gx, gz] = [Math.ceil((2 * hw) / hc), Math.ceil((2 * hd) / hc)];
        const hm = new Float32Array((gx + 1) * (gz + 1));
        for (let j = 0; j <= gz; j++) for (let i = 0; i <= gx; i++) hm[j * (gx + 1) + i] = heapAt(-hw + i * hc, -hd + j * hc);
        const cellOf = (x: number, z: number) => {
          const i = Math.max(0, Math.min(gx, Math.round((x + hw) / hc)));
          const j = Math.max(0, Math.min(gz, Math.round((z + hd) / hc)));
          return j * (gx + 1) + i;
        };
        for (let n = 0; n < 104; n++) {
          const L = 0.042 + fr() * 0.05;
          const a = 0.0072 + fr() * 0.0026;
          const yaw = fr() * Math.PI;
          const [ca, sa] = [Math.cos(yaw), Math.sin(yaw)];
          const ex = Math.abs(ca) * (L / 2) + a;
          const ez = Math.abs(sa) * (L / 2) + a;
          const px = (fr() * 2 - 1) * Math.max(0, bw / 2 - 0.006 - ex);
          const pz = (fr() * 2 - 1) * Math.max(0, bd / 2 - 0.006 - ez);
          // Rotated by yaw about y, the fry's +x runs along (cos yaw, -sin yaw) in xz.
          const along = (t: number): [number, number] => [px + ca * t, pz - sa * t];
          const ts = [-0.5, -0.25, 0, 0.25, 0.5].map((f) => f * L);
          const hs = ts.map((t) => hm[cellOf(...along(t))]!);
          let slope = Math.max(-0.35, Math.min(0.35, (hs[4]! - hs[0]!) / L));
          if (fr() < 0.3) slope *= 0.3;
          let yc = -Infinity;
          for (const [q, t] of ts.entries()) yc = Math.max(yc, hs[q]! - slope * t + a * 0.42);
          const tint = fr();
          dr.add(fm.fry[tint < 0.42 ? 0 : tint < 0.84 ? 1 : 2]!, fryGeom(L, a, (fr() - 0.5) * 0.007, (fr() - 0.5) * 0.5), [cx + px, yc, pzc + pz], [0, yaw, Math.atan(slope)]);
          for (const t of ts) {
            const c0 = cellOf(...along(t));
            hm[c0] = Math.max(hm[c0]!, yc + slope * t + a * 0.6);
          }
        }
      }
    }
  }
  dr.add(p.ss, box(0.001, 0.001, 0.001), [cx, 0.1, zf - 0.01]);

  return (b) => {
    let speed = 0;
    let angle = 0;
    setGlow(glow, 0);
    const drawer = b.joints.get("basket_slide");
    const fanJ = b.joints.get("fan");
    return {
      glows: [glow],
      open: (f) => slideJoint(drawer, "z", Math.max(0, Math.min(1, f)) * 0.32, rest),
      set: (name, v) => {
        if (name === "heat") setGlow(glow, v);
        else if (name === "fan") speed = Math.max(0, Math.min(1, v));
      },
      tick: (_t, dt) => {
        angle += Math.max(0, Math.min(0.1, dt)) * speed * 18;
        turnJoint(fanJ, "y", angle);
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Sous-vide bath                                                                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A sous-vide bath: a clear polycarbonate tub of water, an immersion circulator clamped to its back wall (a PID head
 * with its display over a stainless skirt that hides the heater and the pump, which draws water in through the skirt's
 * slots and jets it out at the bottom), a lid notched round it, and the food sealed in vacuum bags.
 */
const sousVide: Builder = (s, e) => {
  const p = palette();
  const fm = mats();
  const { w, d, h } = e;
  const fine1 = s.fine(1);
  const fine2 = s.fine(2);
  const tw = w - 0.024;
  const td = d - 0.022;
  const tzc = -td / 2;
  const cx = w / 2;
  const rim = 0.206;
  const yw = 0.165;
  // The tub: drawn polycarbonate with a draft, a flanged rim with hand-holds at the sides.
  const prof: [number, number][] = [
    [0.032, 0],
    [0.026, 0.004],
    [0.012, 0.188],
    [0.012, 0.192],
    [0, 0.197],
    [0, 0.203],
    [0.003, rim],
    [0.012, rim],
    [0.015, 0.2],
    [0.029, 0.007],
    [0.035, 0.003],
  ];
  s.add(fm.poly, vessel(s, tw, td, 0.03, prof), [cx, 0, tzc]);
  for (const sg of [-1, 1]) {
    const hs = roundRect(0.012, 0.12, 0.004, 0, 0);
    hs.holes.push(rectPath(new THREE.Path(), 0.004, 0.08, 0.002, -sg * 0.001, 0, true) as THREE.Path);
    const g = plate(s, hs, 0.005);
    g.rotateX(-Math.PI / 2);
    s.add(fm.poly, g, [cx + sg * (tw / 2 + 0.006), rim - 0.005, tzc]);
  }
  const insetAt = (y: number) => 0.029 - ((y - 0.007) * (0.029 - 0.015)) / (0.2 - 0.007);
  // The water: its body, and a surface the circulator's jet keeps moving.
  {
    const i0 = insetAt(0.004) + 0.001;
    const i1 = insetAt(yw) + 0.001;
    s.add(fm.bath, vessel(s, tw, td, 0.03, [[i0 + 0.006, 0.0035], [i0, 0.006], [i1, yw]], [true, false]), [cx, 0, tzc]);
    // The meniscus climbing the wall: a narrow band that reads as the water line.
    if (fine1) s.add(fm.bath, vessel(s, tw, td, 0.03, [[i1 + 0.0012, yw - 0.0035], [i1 + 0.0006, yw + 0.001], [i1 - 0.0004, yw + 0.0022]], [false, false]), [cx, 0, tzc]);
  }
  const sw = tw - 2 * (insetAt(yw) + 0.001);
  const sd = td - 2 * (insetAt(yw) + 0.001);
  const nx = fine2 ? 28 : 12;
  const nz = fine2 ? 44 : 18;
  const surfG = new THREE.PlaneGeometry(sw, sd, nx, nz).rotateX(-Math.PI / 2);
  {
    const pa = surfG.attributes.position as THREE.BufferAttribute;
    const rr = 0.03 - insetAt(yw);
    const [ccx, ccz] = [sw / 2 - rr, sd / 2 - rr];
    for (let i = 0; i < pa.count; i++) {
      let [x, z] = [pa.getX(i), pa.getZ(i)];
      const [ax, az] = [Math.abs(x) - ccx, Math.abs(z) - ccz];
      if (ax > 0 && az > 0) {
        const l = Math.hypot(ax, az);
        if (l > rr) [x, z] = [Math.sign(x) * (ccx + (ax / l) * rr), Math.sign(z) * (ccz + (az / l) * rr)];
      }
      pa.setXYZ(i, x, 0, z);
    }
  }
  const surf = new THREE.Mesh(surfG, fm.bath);
  surf.renderOrder = 2;
  const base = Float32Array.from(surfG.attributes.position!.array as Float32Array);
  s.extra(surf, [cx, yw, tzc]);
  s.tris += (surfG.index!.count / 3) | 0;

  // The circulator, against the back wall.
  const ccz = -td + 0.062;
  const heatLamp = liveLamp(0xff7a2a);
  {
    // Skirt: slotted stainless, the pump's outlet near its foot.
    s.add(p.ss, cyl(s, 0.029, 0.172, "y", 28, 0.029, true), [cx, 0.03 + 0.086, ccz]);
    s.add(p.ssDull, cyl(s, 0.0275, 0.004, "y", 24), [cx, 0.032, ccz]);
    if (fine2)
      for (let i = 0; i < 12; i++) {
        const a = (i / 12) * Math.PI * 2;
        s.add(p.enamel, box(0.005, 0.06, 0.0012), [cx + Math.cos(a) * 0.0293, 0.09, ccz + Math.sin(a) * 0.0293], [0, -a + Math.PI / 2, 0]);
      }
    s.add(p.ssBright, cyl(s, 0.006, 0.02, "z", 12), [cx, 0.045, ccz + 0.035]);
    // Head: a plastic body with a stainless band, the display pod on its front, the set wheel on top.
    const ht = h - 0.014;
    s.add(p.plastic, lathe(s, [[0, 0.2], [0.033, 0.2], [0.035, 0.21], [0.035, ht - 0.014], [0.031, ht - 0.002], [0, ht]], 32), [cx, 0, ccz]);
    s.add(p.ssBright, cyl(s, 0.0355, 0.018, "y", 32), [cx, 0.232, ccz]);
    s.add(p.charcoal, rbox(s, 0.05, 0.07, 0.018, 0.006), [cx, 0.3, ccz + 0.03]);
    slab(s, p.screen, cx - 0.02, cx + 0.02, 0.297, 0.325, ccz + 0.039, ccz + 0.0405, 0.001);
    s.add(fm.sousDisplay, new THREE.PlaneGeometry(0.036, 0.016), [cx, 0.312, ccz + 0.0412]);
    s.add(heatLamp, torus(s, 0.006, 0.0009, "z", Math.PI * 2, 20, 4), [cx, 0.282, ccz + 0.0393]);
    const wheel = s.joint("set_wheel", [cx, ht, ccz]);
    wheel.add(p.ssBright, cyl(wheel, 0.026, 0.008, "y", 32), [0, 0.004, 0]);
    if (fine2) for (let i = 0; i < 30; i++) wheel.add(p.ssBright, box(0.002, 0.007, 0.002), [Math.cos((i / 30) * 6.283) * 0.026, 0.004, Math.sin((i / 30) * 6.283) * 0.026], [0, -(i / 30) * 6.283, 0]);
    // The clamp: a C over the tub's back wall, its thumbscrew pressing a pad on the outside.
    const zw = -td + 0.012;
    slab(s, p.plastic, cx - 0.022, cx + 0.022, 0.2, 0.236, ccz - 0.036, ccz - 0.028, 0.003);
    slab(s, p.plastic, cx - 0.022, cx + 0.022, rim + 0.001, rim + 0.012, -d + 0.008, ccz - 0.028, 0.003);
    slab(s, p.plastic, cx - 0.022, cx + 0.022, 0.14, rim + 0.012, -d + 0.008, -d + 0.016, 0.003);
    s.add(p.ssBright, cyl(s, 0.004, zw - (-d + 0.016) + 0.002, "z", 10), [cx, 0.165, (zw - 0.006 + (-d + 0.016)) / 2]);
    s.add(p.plastic, cyl(s, 0.012, 0.007, "z", 16), [cx, 0.165, -d + 0.0045]);
    s.add(p.rubber, cyl(s, 0.009, 0.003, "z", 14), [cx, 0.165, zw - 0.004]);
    if (fine1) tube(s, p.rubber, [[cx, ht - 0.01, ccz - 0.02], [cx, ht - 0.004, ccz - 0.045], [cx + 0.06, 0.34, -d + 0.006], [cx + 0.09, 0.006, -d + 0.006]], 0.0035, 0.03, 8);
  }
  // The lid, notched round the circulator, a moulded grip at its front.
  if (fine1) {
    const ls = roundRect(tw, td, 0.03, cx, tzc);
    ls.holes.push(circleHole(cx, ccz, 0.04));
    const g = plate(s, ls, 0.004);
    g.rotateX(Math.PI / 2);
    s.add(fm.poly, g, [0, rim + 0.0045, 0]);
    s.add(fm.poly, rbox(s, 0.08, 0.012, 0.02, 0.004), [cx, rim + 0.011, tzc + td / 2 - 0.04]);
  }
  // The bags: a ribeye hanging in its bag from a clip on the left wall, a salmon portion lying on the floor. Each is
  // built lying flat in its own frame (the food in xz); the hanging bag is stood on edge, its open end set over to
  // pass up the inside of the wall under the clip.
  if (fine1) {
    const steakT = 0.026;
    const steakBev = fine2 ? 0.0045 : 0.003;
    // Stood up: the frame's x runs toward the front of the tub (-z), its y toward the right (+x), its z down.
    s.at([0.058, 0.1, -0.2], [Math.PI / 2, 0, -Math.PI / 2], () => {
      s.add(fm.beef, flatFood(s, steakShape(), steakT, steakBev));
      vacBag(
        s,
        fm.bag,
        steakShape(),
        steakT,
        steakBev,
        {
          x0: -0.106,
          x1: 0.104,
          z0: -0.113,
          z1: 0.082,
          seal: { axis: "z", a: -0.104, b: -0.092 },
          bend: (z) => -0.022 * sstep(-0.066, -0.082, z),
        },
        71,
      );
    });
    // Its clip over the tub's wall, pinching the bag's open end against the wall.
    slab(s, p.enamel, 0.004, 0.022, rim - 0.035, rim + 0.007, -0.212, -0.188, 0.002);
    slab(s, p.enamel, 0.004, 0.05, rim + 0.001, rim + 0.007, -0.212, -0.188, 0.002);
    slab(s, p.enamel, 0.044, 0.05, rim - 0.03, rim + 0.007, -0.212, -0.188, 0.002);
    const salT = 0.024;
    const salBev = fine2 ? 0.006 : 0.004;
    s.at([0.18, 0.0065 + salT / 2 + 0.0011, -0.25], [0, 0.1, 0], () => {
      s.add(fm.salmon, flatFood(s, salmonShape(), salT, salBev));
      vacBag(s, fm.bag, salmonShape(), salT, salBev, { x0: -0.1, x1: 0.1, z0: -0.058, z1: 0.06, seal: { axis: "x", a: 0.08, b: 0.092 }, flat: true }, 83);
    });
  }

  return () => {
    let flow = 0.6;
    setLamp(heatLamp, 0);
    const pa = surfG.attributes.position as THREE.BufferAttribute;
    return {
      set: (name, v) => {
        if (name === "heat") setLamp(heatLamp, v > 0.05 ? 1 : 0);
        else if (name === "flow") flow = Math.max(0, Math.min(1, v));
      },
      tick: (t) => {
        const A = 0.0022 * flow;
        for (let i = 0; i < pa.count; i++) {
          const [x, z] = [base[i * 3]!, base[i * 3 + 2]!];
          const r = Math.hypot(x - 0, z - (ccz - tzc + 0.04));
          const y = A * (Math.sin(r * 90 - t * 7) * Math.exp(-r * 6) + 0.35 * Math.sin(x * 40 + t * 2.1) * Math.cos(z * 31 - t * 1.7));
          pa.setY(i, y);
        }
        pa.needsUpdate = true;
        surfG.computeVertexNormals();
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Four-well hot-food table                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/** How far in from its outline a gnPan's inner wall stands at height y below the rim (negative), at this detail. */
function panWallIn(s: Shop, depth: number, y: number): number {
  const lip = 0.008;
  const t = 0.0008;
  const dr = depth * 0.08;
  const [i1, y1] = s.fine(2) ? [lip + dr + t + 0.003, -depth + 0.004 + t] : [lip + dr + t + 0.008, -depth + t];
  const f = Math.max(0, Math.min(1, (y + 0.002) / (y1 + 0.002)));
  return lip + t + f * (i1 - lip - t);
}

/** Signed distance from (x, z) to a w × d rectangle centred on the origin, corners rounded by r: negative inside. */
function rrSdf(x: number, z: number, w: number, d: number, r: number): number {
  const qx = Math.abs(x) - (w / 2 - r);
  const qz = Math.abs(z) - (d / 2 - r);
  return Math.hypot(Math.max(qx, 0), Math.max(qz, 0)) + Math.min(Math.max(qx, qz), 0) - r;
}

/** How a food lies in its pan, in metres, heights from the pan's rim. */
interface Lie {
  /** Where the food meets the wall, below the rim. */
  edge: number;
  /** How far it climbs the steel above that (a sauce's meniscus, rice packed to the wall), and over what distance in. */
  climb: number;
  climbIn: number;
  /** The crown of its heap above the edge (its peak off centre, where the spoon left it), and how far in from the wall the food takes to rise to it. */
  heap: number;
  rise: number;
  /** Lumps: [height, wavelength] pairs. */
  lumps: readonly (readonly [number, number])[];
  /** The depth of the spoon's swirl cut into a soft heap. */
  swirl?: number;
  /** Mesh spacing relative to the usual, for a food with finer relief. */
  step?: number;
  seed: number;
}

/**
 * A food's surface height at (x, z) over its outline w × d (corners r): at the wall it stands at `edge` plus its climb
 * up the steel, and it rises away from the wall into its heap, lumps and swirls. Lumps finer than the mesh can carry
 * at this detail are left to the food's bump map.
 */
function lieHeight(s: Shop, w: number, d: number, r: number, L: Lie): (x: number, z: number) => number {
  const n = noiseField(L.seed);
  const q = rng(L.seed * 7 + 3);
  const [px, pz] = [(q() - 0.5) * w * 0.3, (q() - 0.5) * d * 0.3];
  const [ax, az] = [w / 2 + Math.abs(px), d / 2 + Math.abs(pz)];
  const [ox, oz] = [q() * 30, q() * 30];
  const minLam = s.fine(2) ? 0.015 : s.fine(1) ? 0.04 : 1;
  const lumps = L.lumps.filter(([, lam]) => lam >= minLam);
  return (x, z) => {
    const dIn = Math.max(0, -rrSdf(x, z, w, d, r));
    const crown = (1 - Math.min(1, ((x - px) / ax) ** 2 + ((z - pz) / az) ** 2)) ** 1.3;
    let body = L.heap * crown;
    for (const [k, [a, lam]] of lumps.entries()) body += a * n((x + ox + k * 7.7) / lam, (z + oz + k * 3.1) / lam);
    if (L.swirl && s.fine(2)) {
      // A two-armed spiral about the peak, ridges 40 mm apart, as the back of a spoon leaves a soft heap.
      const ang = Math.atan2(z - pz, x - px);
      body += L.swirl * crown * (0.5 + 0.5 * Math.sin(Math.hypot(x - px, z - pz) * 157 - ang * 2));
    }
    return L.edge + L.climb * Math.exp(-dIn / L.climbIn) + sstep(0, L.rise, dIn) * body;
  };
}

/**
 * The top of a pan's food over its outline w × d (corners r), heights from `h`: meshed in rings that follow the wall
 * in to a line down the middle, where the last ring's two sides meet and are welded, closest together at the wall,
 * where the food meets the steel. UVs are (x, z) in metres.
 */
function foodMesh(s: Shop, w: number, d: number, r: number, h: (x: number, z: number) => number, step = 1): THREE.BufferGeometry {
  const fine2 = s.fine(2);
  const fine1 = s.fine(1);
  const sp = (fine2 ? 0.013 : fine1 ? 0.04 : 1) * step;
  const kc = fine2 ? 5 : fine1 ? 2 : 1;
  const [hx, hz] = [w / 2, d / 2];
  const dMax = Math.min(hx, hz);
  const near = fine2 ? [0, 0.0015, 0.004] : [0];
  const last = near[near.length - 1]!;
  const k = Math.max(1, Math.ceil((dMax - last) / sp));
  const rings = [...near, ...Array.from({ length: k }, (_, i) => last + ((i + 1) / k) * (dMax - last))];
  // Each ring has the same points: kc + 1 round each corner, and the long and short sides split as the outline's are.
  const nx = Math.max(0, Math.round((2 * (hx - r)) / sp) - 1);
  const nz = Math.max(0, Math.round((2 * (hz - r)) / sp) - 1);
  const pos: number[] = [];
  const uv: number[] = [];
  const idx: number[] = [];
  let N = 0;
  for (const dl of rings) {
    const [ex, ez] = [hx - dl, hz - dl];
    // The rings round their corners more as they close in, so none folds to a sharp corner until the last, a line.
    const rj = Math.max(0, Math.min(ex, ez, Math.max(r - dl, Math.min(ex, ez) * 0.35 * (dl / dMax))));
    const [cx, cz] = [ex - rj, ez - rj];
    const corners = [
      [cx, cz, 0, nx],
      [-cx, cz, Math.PI / 2, nz],
      [-cx, -cz, Math.PI, nx],
      [cx, -cz, Math.PI * 1.5, nz],
    ] as const;
    const ring: [number, number][] = [];
    for (const [c, [ox, oz, a0, ns]] of corners.entries()) {
      for (let i = 0; i <= kc; i++) {
        const a = a0 + (i / kc) * (Math.PI / 2);
        ring.push([ox + rj * Math.cos(a), oz + rj * Math.sin(a)]);
      }
      // The straight run to the next corner.
      const [qx, qz] = ring[ring.length - 1]!;
      const [nx2, nz2, na] = corners[(c + 1) % 4]!;
      const [tx, tz] = [nx2 + rj * Math.cos(na), nz2 + rj * Math.sin(na)];
      for (let i = 1; i <= ns; i++) ring.push([qx + ((tx - qx) * i) / (ns + 1), qz + ((tz - qz) * i) / (ns + 1)]);
    }
    N = ring.length;
    for (const [x, z] of ring) {
      pos.push(x, h(x, z), z);
      uv.push(x, z);
    }
  }
  const M = rings.length;
  // The last ring lies on the middle line, each of its points there twice or more: weld them, so the surface runs
  // smoothly across the line.
  const weld = new Map<string, number>();
  const at = (k: number) => {
    if (k < (M - 1) * N) return k;
    const key = Math.round(pos[k * 3]! * 1e6) + "," + Math.round(pos[k * 3 + 2]! * 1e6);
    const w0 = weld.get(key);
    if (w0 !== undefined) return w0;
    weld.set(key, k);
    return k;
  };
  for (let j = 0; j < M - 1; j++)
    for (let i = 0; i < N; i++) {
      const a = j * N + i;
      const b = j * N + ((i + 1) % N);
      const c = at(a + N);
      const e = at(b + N);
      idx.push(a, c, b);
      if (c !== e) idx.push(b, c, e);
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/** A whole green bean `L` long, `r` thick, bowed by `bow`, along x: a slightly flattened pod tapering to its tips. */
function podGeom(L: number, r: number, bow: number, n = 6, m = 8): THREE.BufferGeometry {
  const pos: number[] = [];
  const idx: number[] = [];
  for (let i = 0; i <= m; i++) {
    const e = (2 * i) / m - 1;
    const rr = r * Math.sqrt(Math.max(0, 1 - Math.abs(e) ** 5)) * (1 + 0.05 * Math.sin(e * Math.PI * 4));
    for (let j = 0; j < n; j++) {
      const a = (j / n) * Math.PI * 2;
      pos.push((e * L) / 2, bow * (1 - e * e) + rr * Math.cos(a), rr * Math.sin(a) * 0.86);
    }
  }
  // The section runs anticlockwise from +y toward +z, so (along the section) × (along the pod) faces out.
  for (let i = 0; i < m; i++)
    for (let j = 0; j < n; j++) {
      const a = i * n + j;
      const b = i * n + ((j + 1) % n);
      idx.push(a, b, a + n, b, b + n, a + n);
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/** A kidney bean, 17 × 10 × 7 mm, lying flat along x, its hilum side (+z) drawn in. */
function kidneyBean(): THREE.BufferGeometry {
  const g = new THREE.SphereGeometry(1, 9, 6);
  const p = g.attributes.position as THREE.BufferAttribute;
  for (let i = 0; i < p.count; i++) {
    const [x, y, z] = [p.getX(i), p.getY(i), p.getZ(i)];
    p.setXYZ(i, x * 0.0085, y * 0.0036, z * 0.005 * (z > 0 ? 1 - 0.4 * (1 - x * x) : 1) + 0.0022 * x * x);
  }
  g.computeVertexNormals();
  return g;
}

/**
 * A 60 in, four-well electric hot-food table: a 14-gauge top with four wells (each a drawn liner over its own element,
 * filled with water), hotel pans of food in them, a full-length poly cutting board on brackets along the server's
 * side over the control panel (an infinite switch and a pilot for each well), a drain manifold under the wells with
 * its valve, on four legs with an undershelf.
 */
const steamTable: Builder = (s, e) => {
  const p = palette();
  const fm = mats();
  const food = foodMats();
  const { w, d, h } = e;
  const fine1 = s.fine(1);
  const fine2 = s.fine(2);
  const top = 0.895;
  const zf = -0.192;
  const zb = -d;
  const yb = 0.62;
  const zc = (zf + zb) / 2;
  const wells = [0, 1, 2, 3].map((i) => 0.19 + i * ((w - 0.38) / 3));
  const lamp = liveLamp(0xffa040);

  // Legs and undershelf.
  const legX = [0.045, w - 0.045];
  const legZ = [zf - 0.04, zb + 0.04];
  const r2 = (r: number) => (fine2 ? r : 0);
  for (const x of legX) for (const z of legZ) leg(s, p.ss, x, z, 0, yb);
  slab(s, p.ssDull, legX[0]!, legX[1]!, 0.2, 0.2015, legZ[1]!, legZ[0]!, r2(0.0006));
  if (fine1) for (const z of legZ) slab(s, p.ssDull, legX[0]!, legX[1]!, 0.17, 0.2015, z - 0.0008, z + 0.0008, r2(0.0005));
  for (const x of legX) slab(s, p.ss, x - 0.012, x + 0.012, 0.18, 0.2, legZ[1]!, legZ[0]!, fine2 ? 0.002 : 0);
  // The body: sides, back, bottom, front panel; the top with its wells, a marine edge at the back.
  slab(s, p.ss, 0, w, yb, top, zf - 0.0012, zf, r2(0.0006));
  slab(s, p.ss, 0, w, yb, top, zb, zb + 0.0012, r2(0.0006));
  for (const x of [0, w - 0.0012]) slab(s, p.ss, x, x + 0.0012, yb, top, zb, zf, r2(0.0006));
  slab(s, p.ssDull, 0.0012, w - 0.0012, yb, yb + 0.0012, zb + 0.0012, zf - 0.0012);
  {
    const sh = roundRect(w, zf - zb, 0.004, w / 2, zc);
    for (const x of wells) rectPath(sh, 0.31, 0.51, 0.01, x, zc, true);
    const g = plate(s, sh, 0.0019, r2(0.0007));
    g.rotateX(Math.PI / 2);
    s.add(p.ss, g, [0, top, 0]);
    slab(s, p.ss, 0, w, top, h - 0.004, zb, zb + 0.0015, r2(0.0006));
    s.add(p.ss, cyl(s, 0.004, w, "x", 12), [w / 2, h - 0.004, zb + 0.003]);
    // The bull-nosed front edge.
    if (fine1) s.add(p.ss, cyl(s, 0.012, w, "x", 20, 0.012, true), [w / 2, top - 0.012, zf + 0.0005]);
  }
  // Wells: drawn liners with water standing in them.
  if (fine1)
    for (const x of wells) {
      s.add(p.ssDull, vessel(s, 0.31, 0.51, 0.01, [[0.03, -0.205], [0.012, -0.2], [0, -0.004], [0, -0.0018]], [true, false]), [x, top, zc]);
      s.add(p.water, box(0.27, 0.002, 0.47), [x, top - 0.17, zc]);
    }
  // The pans and their food: mac and cheese; green beans and rice in halves; chili; mash, corn and gravy in thirds.
  // Each food meets its pan's drafted wall at its own level, 0.2 mm into the 0.8 mm steel (crossing it, never lying
  // along it), and rises from there into its heap, below the rim.
  const py = top + 0.003;
  const depth = 0.1;
  const fill = (x: number, z: number, pw: number, pd: number, m: Mat, L: Lie) => {
    if (fine1) s.add(p.ss, gnPan(s, pw, pd, depth), [x, py, z]);
    const iw = fine1 ? panWallIn(s, depth, L.edge + L.climb) - 0.0002 : 0.004;
    const [fw, fd, fr] = [pw - 2 * iw, pd - 2 * iw, Math.max(0.002, 0.028 - iw)];
    const h = lieHeight(s, fw, fd, fr, L);
    s.add(m, foodMesh(s, fw, fd, fr, h, L.step), [x, py, z]);
    return { x, z, h, fw, fd };
  };
  // Baked mac: a low dome, lumpy where the elbows pile. Green beans: a loose heap. Rice: heaped up to the middle and
  // packed against the steel. Chili: a thick sauce, nearly level, its meniscus climbing the wall.
  const mac = fill(wells[0]!, zc, 0.325, 0.53, food.mac, { edge: -0.033, climb: 0.0015, climbIn: 0.003, heap: 0.01, rise: 0.03, lumps: [[0.005, 0.045], [0.002, 0.022]], seed: 1 });
  const greens = fill(wells[1]!, zc - 0.133, 0.325, 0.265, food.greensHeap, { edge: -0.045, climb: 0, climbIn: 0.003, heap: 0.016, rise: 0.04, lumps: [[0.007, 0.05], [0.003, 0.025]], seed: 2 });
  fill(wells[1]!, zc + 0.133, 0.325, 0.265, food.rice, { edge: -0.036, climb: 0.002, climbIn: 0.004, heap: 0.022, rise: 0.045, lumps: [[0.004, 0.04], [0.0018, 0.02]], seed: 3 });
  const chili = fill(wells[2]!, zc, 0.325, 0.53, food.chili, { edge: -0.031, climb: 0.003, climbIn: 0.0025, heap: 0.006, rise: 0.02, lumps: [[0.0035, 0.05], [0.0015, 0.028]], seed: 4 });
  // Mash swirled by the spoon; corn; gravy, level but for its meniscus.
  const thirds: [Mat, Lie][] = [
    [food.mash, { edge: -0.032, climb: 0.001, climbIn: 0.003, heap: 0.016, rise: 0.035, lumps: [[0.004, 0.04]], swirl: 0.008, step: 0.7, seed: 5 }],
    [food.corn, { edge: -0.033, climb: 0.0012, climbIn: 0.003, heap: 0.01, rise: 0.03, lumps: [[0.003, 0.03], [0.0015, 0.016]], seed: 6 }],
    [fm.gravy, { edge: -0.036, climb: 0.002, climbIn: 0.002, heap: 0, rise: 0.01, lumps: [[0.0004, 0.06]], seed: 7 }],
  ];
  for (const [k, [m, L]] of thirds.entries()) fill(wells[3]!, zc + (k - 1) * 0.177, 0.325, 0.176, m, L);
  // Adaptor bars across the wells between the fractional pans.
  for (const z of [zc]) slab(s, p.ss, wells[1]! - 0.16, wells[1]! + 0.16, top, top + 0.003, z - 0.0015, z + 0.0015, r2(0.0005));
  for (const z of [zc - 0.0885, zc + 0.0885]) slab(s, p.ss, wells[3]! - 0.16, wells[3]! + 0.16, top, top + 0.003, z - 0.0015, z + 0.0015, r2(0.0005));
  if (fine2) {
    const r = rng(5);
    // Elbows standing out of the mac: each end's tube sunk at least half into the sauce, the arch well clear of it.
    const eu = new THREE.Euler();
    const v = new THREE.Vector3();
    const [R, tr] = [0.0058, 0.0023];
    for (let tries = 0, placed = 0; tries < 400 && placed < 70; tries++) {
      const arc = 2.1 + r() * 0.6;
      const t0 = (Math.PI - arc) / 2;
      const [ex, ez] = [(r() - 0.5) * (mac.fw - 0.03), (r() - 0.5) * (mac.fd - 0.03)];
      const rot: V3 = [(r() - 0.5) * 1, r() * Math.PI * 2, (r() - 0.5) * 1];
      eu.set(rot[0], rot[1], rot[2]);
      const at = (a: number) => v.set(R * Math.cos(a), R * Math.sin(a), 0).applyEuler(eu);
      let cy = Infinity;
      for (const a of [t0, Math.PI - t0]) {
        at(a);
        cy = Math.min(cy, mac.h(ex + v.x, ez + v.z) - v.y - tr * 0.5);
      }
      at(Math.PI / 2);
      if (cy + v.y + tr < mac.h(ex + v.x, ez + v.z) + 0.003) continue;
      const g = new THREE.TorusGeometry(R, tr, 5, 8, arc);
      g.rotateZ(t0);
      s.add(r() < 0.4 ? food.elbowBrown : food.elbow, g, [mac.x + ex, py + cy, mac.z + ez], rot);
      placed++;
    }
    // Kidney beans in the chili, mostly under the sauce, their backs showing.
    for (let i = 0; i < 80; i++) {
      const [bx, bz] = [(r() - 0.5) * (chili.fw - 0.03), (r() - 0.5) * (chili.fd - 0.03)];
      s.add(food.bean, kidneyBean(), [chili.x + bx, py + chili.h(bx, bz) - 0.0012 - r() * 0.0012, chili.z + bz], [(r() - 0.5) * 0.3, r() * Math.PI * 2, (r() - 0.5) * 0.3]);
    }
    // Green beans dropped one by one onto the heap, each coming to rest on it or on the beans already there: a height
    // map in 1 cm cells, starting from the heap's own surface.
    {
      const [hw, hd] = [greens.fw / 2 - 0.006, greens.fd / 2 - 0.006];
      const hc = 0.01;
      const [gx, gz] = [Math.ceil((2 * hw) / hc), Math.ceil((2 * hd) / hc)];
      const hm = new Float32Array((gx + 1) * (gz + 1));
      for (let j = 0; j <= gz; j++) for (let i = 0; i <= gx; i++) hm[j * (gx + 1) + i] = greens.h(-hw + i * hc, -hd + j * hc);
      const cell = (x: number, z: number) => Math.max(0, Math.min(gz, Math.round((z + hd) / hc))) * (gx + 1) + Math.max(0, Math.min(gx, Math.round((x + hw) / hc)));
      for (let n = 0; n < 64; n++) {
        const L = 0.06 + r() * 0.04;
        const rad = 0.0036 + r() * 0.0006;
        const yaw = r() * Math.PI;
        const [ca, sa] = [Math.cos(yaw), Math.sin(yaw)];
        const bx = (r() * 2 - 1) * Math.max(0, hw - Math.abs(ca) * (L / 2) - rad);
        const bz = (r() * 2 - 1) * Math.max(0, hd - Math.abs(sa) * (L / 2) - rad);
        // Turned by yaw about y, the pod's +x runs along (cos yaw, -sin yaw) in xz.
        const along = (t: number) => cell(bx + ca * t, bz - sa * t);
        const ts = [-0.5, -0.25, 0, 0.25, 0.5].map((f) => f * L);
        const hs = ts.map((t) => hm[along(t)]!);
        const slope = Math.max(-0.3, Math.min(0.3, (hs[4]! - hs[0]!) / L));
        let yc = -Infinity;
        for (const [k, t] of ts.entries()) yc = Math.max(yc, hs[k]! - slope * t + rad * 0.55);
        // Nothing heaped up over the rim.
        if (yc > -0.012) continue;
        const tone = r() < 0.45 ? food.greensDark : food.greens;
        s.add(tone, podGeom(L, rad, (r() - 0.35) * 0.008), [greens.x + bx, py + yc, greens.z + bz], [0, yaw, Math.atan(slope)]);
        for (const t of ts) hm[along(t)] = Math.max(hm[along(t)]!, yc + slope * t + rad * 0.8);
      }
    }
    // A serving spoon in the mac, a ladle in the chili.
    const bowl: [number, number][] = [];
    for (let i = 0; i <= 8; i++) bowl.push([Math.sin((i / 8) * Math.PI / 2) * 0.042, -Math.cos((i / 8) * Math.PI / 2) * 0.034]);
    bowl.push([0.0435, 0.0012]);
    for (let i = 8; i >= 0; i--) bowl.push([Math.sin((i / 8) * Math.PI / 2) * 0.0405, -Math.cos((i / 8) * Math.PI / 2) * 0.0325]);
    // The spoon: a shallow oval bowl half buried, its handle up over the pan's back rim.
    const sx = wells[0]! + 0.02;
    const sz = zc + 0.1;
    const sy = py + mac.h(sx - mac.x, sz - mac.z);
    s.add(p.ssBright, lathe(s, bowl, 24), [sx, sy - 0.002, sz], [-0.12, 0, 0], [0.75, 0.45, 1.05]);
    tube(s, p.ssBright, [[sx, sy, sz + 0.042], [sx + 0.006, py - 0.005, sz + 0.1], [sx + 0.012, py + 0.016, zc + 0.262], [sx + 0.013, py + 0.008, zc + 0.3]], 0.0042, 0.02, 8);
    // The ladle: its bowl full, sitting a little proud of the chili.
    const lx = wells[2]! - 0.06;
    const lz = zc - 0.1;
    const ly = py + chili.h(lx - chili.x, lz - chili.z) + 0.006;
    s.add(p.ssBright, lathe(s, bowl, 24), [lx, ly, lz], [0.18, 0, -0.12]);
    tube(s, p.ssBright, [[lx + 0.042, ly + 0.006, lz], [wells[2]! + 0.12, py + 0.026, lz - 0.03], [wells[2]! + 0.163, py + 0.03, lz - 0.035], [wells[2]! + 0.168, py + 0.004, lz - 0.035]], 0.0045, 0.008, 8);
  }
  // The server's side: the control panel, the cutting board on its brackets.
  const kY = 0.73;
  const knobs: string[] = [];
  for (const [i, x] of wells.entries()) {
    const name = `knob_${i}`;
    const j = s.joint(name, [x, kY, zf]);
    dial(j, p.bakelite, p.ssBright, [0, 0, 0], 0.02);
    knobs.push(name);
    if (fine1) s.add(lamp, cyl(s, 0.005, 0.005, "z", 14), [x + 0.06, kY + 0.02, zf + 0.0025]);
    if (fine2) s.add(p.ssBright, torus(s, 0.006, 0.0012, "z", Math.PI * 2, 16, 4), [x + 0.06, kY + 0.02, zf + 0.001]);
    if (fine2) s.add(p.enamel, cyl(s, 0.029, 0.0006, "z", 32), [x, kY, zf + 0.0007]);
  }
  dataPlate(s, [w - 0.1, 0.68, zf], [0, 0, 0], 0.08, 0.045);
  {
    // The board: 12.7 mm HDPE, its edges chamfered 2 mm.
    const [bw, bd, bt] = [w - 0.01, -0.004 - (zf + 0.002), 0.0127];
    let g: THREE.BufferGeometry;
    if (fine1) {
      const ch = 0.002;
      g = new THREE.ExtrudeGeometry(roundRect(bw - 2 * ch, bd - 2 * ch, 0.006), { depth: bt - 2 * ch, bevelEnabled: true, bevelThickness: ch, bevelSize: ch, bevelSegments: 1, curveSegments: fine2 ? 4 : 2 });
      g.translate(0, 0, ch);
      g.rotateX(-Math.PI / 2);
    } else {
      g = box(bw, bt, bd).translate(0, bt / 2, 0);
      boxUVs(g, 1);
    }
    s.add(food.board, g, [w / 2, 0.8, (zf + 0.002 - 0.004) / 2]);
  }
  if (fine1) for (const x of [0.15, w / 2, w - 0.15]) {
    const tri = new THREE.Shape([new THREE.Vector2(0, 0), new THREE.Vector2(0.16, 0), new THREE.Vector2(0, -0.08)]);
    const g = plate(s, tri, 0.003);
    g.rotateY(-Math.PI / 2);
    s.add(p.ss, g, [x - 0.0015, 0.8, zf]);
    s.add(p.ss, box(0.03, 0.002, 0.17), [x, 0.799, zf + 0.085]);
  }
  // Under the wells: a drain from each to a manifold, its ball valve at the front right; the supply's junction box.
  if (fine1) {
    const my = yb - 0.03;
    const mz = zc + 0.05;
    if (fine2) for (const x of wells) tube(s, p.copper, [[x, yb, zc + 0.15], [x, my, zc + 0.15], [x, my, mz]], 0.008, 0.02, 10, false);
    else for (const x of wells) s.add(p.copper, rod(0.008, yb - my, "y", 5), [x, (yb + my) / 2, mz]);
    tube(s, p.copper, [[wells[0]! - 0.02, my, mz], [w - 0.12, my, mz], [w - 0.12, my, zf - 0.08]], 0.011, 0.03, 12);
    s.add(p.brass, cyl(s, 0.018, 0.05, "z", 14), [w - 0.12, my, zf - 0.08]);
    s.add(fm.red, box(0.012, 0.006, 0.09), [w - 0.12, my + 0.024, zf - 0.07], [0, 0.4, 0]);
    s.add(fm.greyBox, rbox(s, 0.12, 0.08, 0.05, 0.004), [0.3, yb - 0.04, zb + 0.06]);
    tube(s, fm.greyBox, [[0.3, yb - 0.08, zb + 0.06], [0.3, 0.25, zb + 0.06], [0.3, 0.205, zb + 0.08]], 0.009, 0.04, 8);
  }
  const steam = new Haze({ x0: 0.08, x1: w - 0.08, y0: top + 0.0, y1: top + 0.3, z0: zb + 0.12, z1: zf - 0.1 }, 18, 0.16, 61, 0.12, 0xffffff, 0.35);
  if (fine1) {
    s.extra(steam.group);
    s.tris += steam.tris;
  }

  return (b) => {
    setLamp(lamp, 0);
    return {
      set: (name, v) => {
        if (name === "steam") steam.set(v);
        else if (name === "heat") {
          setLamp(lamp, v > 0.05 ? 1 : 0);
          for (const k of knobs) turnJoint(b.joints.get(k), "z", -Math.max(0, Math.min(1, v)) * Math.PI * 1.6);
        }
      },
      tick: (t) => steam.tick(t),
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Sanitising-rinse booster heater                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A pressure gauge or thermometer facing +z from its back at z = 0: a chrome case, the dial face, a needle (a joint
 * named `needle`, zero at the scale's start), a glass lens, a stem out of its back along -z.
 */
function gauge(s: Shop, face: Mat, at: V3, r: number, needle: string): void {
  const p = palette();
  s.at(at, [0, 0, 0], () => {
    if (s.detail === 0) s.add(p.chrome, cyl(s, r * 1.08, 0.023, "z", 12), [0, 0, 0.0115]);
    else s.add(p.chrome, lathe(s, [[0, 0], [r * 1.02, 0], [r * 1.06, 0.004], [r * 1.06, 0.018], [r * 1.12, 0.021], [r * 1.08, 0.025], [r * 0.98, 0.023], [0, 0.023]], 40), [0, 0, 0], [Math.PI / 2, 0, 0]);
    s.add(face, new THREE.CircleGeometry(r, s.seg(40, 12)), [0, 0, 0.0235]);
    const n = s.joint(needle, [0, 0, 0.0245]);
    n.add(p.enamel, taperBar(r * 0.9, 0.0026, 0.0008, 0.0008), [-r * 0.12, 0, 0]);
    n.add(p.enamel, cyl(n, 0.003, 0.002, "z", 12), [0, 0, 0.0008]);
    s.add(p.glass, cyl(s, r * 1.0, 0.0015, "z", 32), [0, 0, 0.0262]);
  });
}
/** A flat needle `len` long along +x from x = 0, tapering from w0 to w1, t thick. */
function taperBar(len: number, w0: number, w1: number, t: number): THREE.BufferGeometry {
  const sh = new THREE.Shape([new THREE.Vector2(0, -w0 / 2), new THREE.Vector2(len, -w1 / 2), new THREE.Vector2(len, w1 / 2), new THREE.Vector2(0, w0 / 2)]);
  return new THREE.ExtrudeGeometry(sh, { depth: t, bevelEnabled: false });
}
/** Needle turn for a fraction of full scale: the 270° scale starts at 7:30 and runs clockwise. */
function needleAngle(f: number): number {
  return -(Math.PI * 0.75 + Math.max(0, Math.min(1, f)) * Math.PI * 1.5);
}

/**
 * An electric booster heater: an insulated stainless tank in a cabinet on legs, raising the dish machine's rinse water
 * to 82 °C with immersion elements. Cold water comes up through a ball valve past a pressure gauge into its foot; the
 * hot outlet rises past a thermometer to a tee carrying the pressure-relief valve (its discharge piped down toward the
 * floor) and runs on to the machine. Supply cable enters a junction box through flexible conduit.
 */
const boosterHeater: Builder = (s, e) => {
  const p = palette();
  const fm = mats();
  const { w, d, h } = e;
  const fine1 = s.fine(1);
  const fine2 = s.fine(2);
  const yb = 0.1;
  const yt = 0.39;
  const zf = -0.09;
  const zb = -d + 0.01;
  const [x0, x1] = [0.025, w - 0.025];
  const pz = -0.045;
  const pr = 0.011;
  const lamp = liveLamp(0xffa040);

  for (const x of [x0 + 0.035, x1 - 0.035]) for (const z of [zf - 0.035, zb + 0.035]) leg(s, p.ss, x, z, 0, yb, 0.016);
  // The cabinet: wrapper and front panel on screws.
  slab(s, p.ss, x0, x1, yb, yt, zb, zf - 0.0015, 0.006);
  slab(s, p.ss, x0 + 0.01, x1 - 0.01, yb + 0.01, yt - 0.01, zf - 0.002, zf, 0.002);
  if (fine2) for (const x of [x0 + 0.02, x1 - 0.02]) for (const y of [yb + 0.02, yt - 0.02]) screw(s, p.ssBright, [x, y, zf], 0.0025);
  if (fine1) louvres(s, [w / 2, yt, zb + 0.08], [-Math.PI / 2, 0, 0], 4, 0.18, 0.018);
  // Controls: power rocker, element pilot, the low-water reset; data plate.
  slab(s, p.charcoal, w / 2 - 0.07, w / 2 + 0.07, yt - 0.09, yt - 0.03, zf, zf + 0.0015, 0.003);
  slab(s, p.plastic, w / 2 - 0.05, w / 2 - 0.025, yt - 0.075, yt - 0.045, zf + 0.0015, zf + 0.009, 0.002);
  s.add(lamp, cyl(s, 0.006, 0.006, "z", 14), [w / 2 + 0.005, yt - 0.06, zf + 0.003]);
  s.add(fm.red, cyl(s, 0.007, 0.008, "z", 14), [w / 2 + 0.04, yt - 0.06, zf + 0.004]);
  dataPlate(s, [w / 2 - 0.06, yb + 0.09, zf + 0.0015], [0, 0, 0], 0.09, 0.055);

  // Inlet: up from the floor through a ball valve and past the pressure gauge into the tank's foot.
  const ix = x0 + 0.06;
  const iy = yb + 0.06;
  tube(s, p.copper, [[ix, 0.002, pz], [ix, iy, pz], [ix, iy, zf - 0.02]], pr, 0.02, 14);
  if (fine1) s.add(p.brass, cyl(s, 0.018, 0.012, "z", 6), [ix, iy, zf + 0.006]);
  s.add(p.brass, cyl(s, 0.017, 0.05, "y", 6), [ix, 0.045, pz]);
  if (fine1) s.add(p.brass, ball(s, 0.02, 16), [ix, 0.045, pz]);
  s.add(fm.red, box(0.07, 0.006, 0.014), [ix + 0.035, 0.045, pz + 0.022], [0, 0, 0]);
  if (fine1) s.add(p.brass, cyl(s, 0.006, 0.02, "z", 10), [ix, 0.045, pz + 0.015]);
  // Pressure gauge on a tee in the riser.
  s.add(p.brass, cyl(s, 0.015, 0.03, "y", 6), [ix, 0.1, pz]);
  if (fine1) s.add(p.brass, cyl(s, 0.006, 0.02, "z", 10), [ix, 0.1, pz + 0.018]);
  gauge(s, fm.dialPress, [ix, 0.1, pz + 0.026], 0.026, "pressure_needle");
  // Outlet: out of the tank's top front, past the thermometer, up to a tee; the run turns back over the top to the
  // machine, the relief valve on the tee's other side, its discharge piped down to an air gap over the floor.
  const ox = x1 - 0.07;
  const oy = yt - 0.05;
  const ry = h - 0.05;
  const lx = x0 + 0.03;
  tube(s, p.copper, [[ox, oy, zf - 0.01], [ox, oy, pz], [ox, ry - 0.016, pz]], pr, 0.022, 14, false);
  tube(s, p.copper, [[ox - 0.02, ry, pz], [lx, ry, pz], [lx, ry, -d + 0.05]], pr, 0.03, 14, false);
  if (fine1) s.add(p.brass, cyl(s, 0.018, 0.012, "z", 6), [ox, oy, zf + 0.006]);
  s.add(p.brass, cyl(s, 0.015, 0.032, "y", 6), [ox, oy, pz]);
  if (fine1) s.add(p.brass, cyl(s, 0.006, 0.02, "z", 10), [ox, oy, pz + 0.018]);
  gauge(s, fm.dialTemp, [ox, oy, pz + 0.026], 0.026, "temp_needle");
  // The tee, and the union where the run leaves for the machine.
  s.add(p.brass, cyl(s, 0.0145, 0.044, "x", 16), [ox, ry, pz]);
  s.add(p.brass, cyl(s, 0.0145, 0.02, "y", 16), [ox, ry - 0.012, pz]);
  if (fine1) for (const sx of [-1, 1]) s.add(p.brass, cyl(s, 0.016, 0.006, "x", 6), [ox + sx * 0.022, ry, pz]);
  s.add(p.brass, cyl(s, 0.017, 0.022, "z", 6), [lx, ry, -d + 0.06]);
  // The relief valve: inlet from the side, spring bonnet and test lever above, outlet down.
  const vx = w - 0.042;
  tube(s, p.brass, [[ox + 0.022, ry, pz], [vx - 0.012, ry, pz]], 0.009, 0.01, 12, false);
  s.at([vx, ry, pz], [0, 0, 0], () => {
    s.add(p.brass, cyl(s, 0.0145, 0.014, "x", 6), [-0.012, 0, 0]);
    s.add(p.brass, lathe(s, [[0.012, -0.026], [0.016, -0.02], [0.017, 0.006], [0.013, 0.012], [0.011, 0.034], [0.008, 0.038], [0, 0.039]], 20), [0, 0, 0]);
    s.add(p.brass, cyl(s, 0.0125, 0.012, "y", 6), [0, -0.026, 0]);
    // The test lever on its pin, its stem through the cap.
    if (fine1) s.add(p.brass, cyl(s, 0.0022, 0.012, "z", 8), [0.006, 0.042, 0]);
    s.add(p.brass, rbox(s, 0.042, 0.004, 0.009, 0.0015), [-0.012, 0.043, 0], [0, 0, -0.08]);
    if (fine2) s.add(p.brass, cyl(s, 0.0018, 0.01, "y", 8), [0, 0.04, 0]);
  });
  tube(s, p.copper, [[vx, ry - 0.03, pz], [vx, 0.035, pz]], 0.009, 0.01, 12);
  if (fine1) s.add(p.brass, cyl(s, 0.0105, 0.012, "y", 16), [vx, 0.035, pz]);
  // The supply: a junction box at the front, flexible conduit down to the floor.
  {
    const jx = w / 2 + 0.08;
    slab(s, fm.greyBox, jx - 0.045, jx + 0.045, yb + 0.03, yb + 0.12, zf, zf + 0.036, 0.003);
    slab(s, fm.greyBox, jx - 0.043, jx + 0.043, yb + 0.032, yb + 0.118, zf + 0.036, zf + 0.039, 0.002);
    if (fine2) for (const sx of [-1, 1]) for (const sy of [-1, 1]) screw(s, p.ssBright, [jx + sx * 0.036, yb + 0.075 + sy * 0.036, zf + 0.039], 0.0022);
    s.add(p.ssBright, cyl(s, 0.012, 0.014, "y", 6), [jx, yb + 0.023, zf + 0.018]);
    if (fine1) {
      const n = fine2 ? 18 : 0;
      tube(s, fm.greyBox, [[jx, yb + 0.016, zf + 0.018], [jx, 0.002, zf + 0.018]], 0.01, 0.02, 10);
      for (let i = 0; i < n; i++) s.add(fm.greyBox, torus(s, 0.0102, 0.0012, "y", Math.PI * 2, 14, 4), [jx, 0.006 + (i * (yb + 0.01)) / n, zf + 0.018]);
    }
  }

  return (b) => {
    setLamp(lamp, 0);
    turnJoint(b.joints.get("temp_needle"), "z", needleAngle(0.68));
    turnJoint(b.joints.get("pressure_needle"), "z", needleAngle(0.35));
    return {
      set: (name, v) => {
        if (name === "temp") turnJoint(b.joints.get("temp_needle"), "z", needleAngle(v));
        else if (name === "pressure") turnJoint(b.joints.get("pressure_needle"), "z", needleAngle(v));
        else if (name === "heat") setLamp(lamp, v > 0.05 ? 1 : 0);
      },
    };
  };
};

register({
  "NSF-COMBI-6PAN": combi(6),
  "NSF-COMBI-10PAN": combi(10),
  "NSF-COMBI-20PAN-ROLLIN": combi(20),
  "NSF-HOLDING-CABINET": holdingCabinet,
  "NSF-HOLDING-CABINET-EMPTY": holdingCabinet,
  "NSF-MICROWAVE-COMM": microwave,
  "NSF-TOASTER-CONVEYOR": conveyorToaster,
  "NSF-AIRFRYER-COMM": airFryer,
  "NSF-SOUSVIDE-BATH": sousVide,
  "NSF-STEAMTABLE-HOT-60": steamTable,
  "NSF-BOOSTER-HEATER-36KW": boosterHeater,
});
