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
  const [c, x] = canvas2d(256, 320);
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

/** A soft puff for steam and vapour. */
function puff(): HTMLCanvasElement {
  const S = 128;
  const [c, x] = canvas2d(S, S);
  const r = rng(77);
  for (let i = 0; i < 9; i++) {
    const px = S / 2 + (r() - 0.5) * 34;
    const py = S / 2 + (r() - 0.5) * 34;
    const rad = 28 + r() * 22;
    const g = x.createRadialGradient(px, py, 0, px, py, rad);
    g.addColorStop(0, "rgba(255,255,255,0.32)");
    g.addColorStop(1, "rgba(255,255,255,0)");
    x.fillStyle = g;
    x.fillRect(0, 0, S, S);
  }
  return c;
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
  /** Food: what the pans hold. */
  mac: THREE.MeshStandardMaterial;
  chili: THREE.MeshStandardMaterial;
  bean: THREE.MeshStandardMaterial;
  greens: THREE.MeshStandardMaterial;
  rice: THREE.MeshStandardMaterial;
  mash: THREE.MeshStandardMaterial;
  corn: THREE.MeshStandardMaterial;
  fries: THREE.MeshStandardMaterial;
  dough: THREE.MeshStandardMaterial;
  beef: THREE.MeshStandardMaterial;
  fat: THREE.MeshStandardMaterial;
  salmon: THREE.MeshStandardMaterial;
  gravy: THREE.MeshStandardMaterial;
  toast: THREE.MeshStandardMaterial;
  crust: THREE.MeshStandardMaterial;
  crumb: THREE.MeshStandardMaterial;
  crustRaw: THREE.MeshStandardMaterial;
  steamTex: THREE.CanvasTexture;
}
let fam: FamilyMats | null = null;

function mats(): FamilyMats {
  if (fam) return fam;
  const std = (p: THREE.MeshStandardMaterialParameters, uv: "box" | "own" = "box", uvScale?: number) => {
    const m = new THREE.MeshStandardMaterial(p);
    m.userData.uv = uv;
    if (uvScale) m.userData.uvScale = uvScale;
    return m;
  };
  const basic = (c: HTMLCanvasElement) => {
    const m = new THREE.MeshBasicMaterial({ map: ctex(c) });
    m.userData.uv = "own";
    return m;
  };
  const perfTex = perforation(16, 0.62);
  perfTex.repeat.set(1, 1);
  const meshTex = perforation(32, 0.6);
  const speck = (base: string, dots: string[], n: number, seed: number) => {
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
  };
  const poly = new THREE.MeshPhysicalMaterial({ color: 0xe6eef2, metalness: 0, roughness: 0.12, transparent: true, opacity: 0.2, clearcoat: 1, depthWrite: false, side: THREE.DoubleSide });
  poly.userData.uv = "own";
  const bag = new THREE.MeshPhysicalMaterial({ color: 0xf2f4f5, metalness: 0, roughness: 0.18, transparent: true, opacity: 0.16, clearcoat: 1, clearcoatRoughness: 0.2, depthWrite: false });
  bag.userData.uv = "own";
  const ui = basic(combiUi());
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
    red: std({ color: 0xc8302a, metalness: 0.1, roughness: 0.45 }),
    greyBox: std({ color: 0x8b9096, metalness: 0.2, roughness: 0.55 }),
    mica: std({ color: 0xcdb98e, metalness: 0, roughness: 0.5, side: THREE.DoubleSide }, "box"),
    wool: std({ color: 0xffffff, metalness: 0, roughness: 1, map: speck("#cdc29e", ["#c2b68e", "#d8cfae", "#bfb38a"], 2400, 15) }, "box", 5),
    lamp: (() => {
      const m = new THREE.MeshBasicMaterial({ color: 0xffe0a8 });
      m.userData.uv = "own";
      return m;
    })(),
    mac: std({ color: 0xe8b443, roughness: 0.5, map: speck("#e8b443", ["#f6cf6a", "#c98a25", "#f2c25a"], 900, 11) }, "box", 6),
    chili: std({ color: 0xffffff, roughness: 0.42, map: speck("#6a2414", ["#4a160c", "#8a3a1c", "#3a120a", "#a54a22"], 1400, 12) }, "box", 5),
    bean: std({ color: 0x4a140e, roughness: 0.35 }),
    greens: std({ color: 0x4f8030, roughness: 0.45 }),
    rice: std({ color: 0xffffff, roughness: 0.8, map: speck("#eeeadf", ["#ffffff", "#d9d3c4", "#f8f6ef"], 2600, 13) }, "box", 6),
    mash: std({ color: 0xe9dcb8, roughness: 0.72 }),
    corn: std({ color: 0xffffff, roughness: 0.4, map: speck("#f0c02c", ["#ffd84a", "#d9a218", "#fbe07a"], 2000, 14) }, "box", 8),
    fries: std({ color: 0xe3ac4a, roughness: 0.6 }),
    dough: std({ color: 0xead8b4, roughness: 0.88 }),
    beef: std({ color: 0x8a2a2c, roughness: 0.48 }),
    fat: std({ color: 0xe8d8c8, roughness: 0.6 }),
    salmon: std({ color: 0xe07a52, roughness: 0.45 }),
    gravy: std({ color: 0x6a4224, roughness: 0.22 }),
    toast: std({ color: 0xffffff, roughness: 0.85, map: speck("#b77a3c", ["#8e5426", "#c98d4a", "#6e3e1c"], 1800, 16) }, "box", 6),
    crust: std({ color: 0x6e3e1e, roughness: 0.7 }),
    crumb: std({ color: 0xffffff, roughness: 0.9, map: speck("#efe2c4", ["#e0cfa8", "#f7eedb", "#d8c49a"], 2200, 17) }, "box", 6),
    crustRaw: std({ color: 0xc8a572, roughness: 0.75 }),
    steamTex: ctex(puff()),
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

/** The steam a pan or a cavity gives off: soft sprites drifting in a box, their strength live. */
class Haze {
  readonly group = new THREE.Group();
  private items: { s: THREE.Sprite; m: THREE.SpriteMaterial; x: number; y: number; z: number; ph: number; sp: number; sc: number }[] = [];
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
    for (let i = 0; i < n; i++) {
      const m = new THREE.SpriteMaterial({ map: mats().steamTex, color: tint, transparent: true, depthWrite: false, opacity: 0, rotation: r() * 6.28 });
      m.userData.uv = "own";
      const sp = new THREE.Sprite(m);
      sp.renderOrder = 1;
      const sc = size * (0.7 + r() * 0.6);
      sp.scale.set(sc, sc, sc);
      this.group.add(sp);
      this.items.push({ s: sp, m, x: r(), y: r(), z: r(), ph: r() * 6.28, sp: 0.6 + r() * 0.8, sc });
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
      let u: number;
      let y: number;
      if (this.rise > 0) {
        // Wisps rising off a surface, growing and fading as they go.
        u = (it.y + t * this.rise * it.sp) % 1;
        y = b.y0 + u * (b.y1 - b.y0);
        const k = Math.sin(Math.PI * u);
        it.s.scale.setScalar(it.sc * (0.5 + u));
        it.m.opacity = this.peak * this.level * k;
        it.s.position.set(b.x0 + (b.x1 - b.x0) * (it.x + 0.06 * Math.sin(t * 0.9 + it.ph) * u), y, b.z0 + (b.z1 - b.z0) * it.z);
      } else {
        // Churned round by the fan: a slow loop in the cavity.
        const a = t * 0.35 * it.sp + it.ph;
        const fx = 0.5 + (it.x - 0.5) * 0.7 + 0.14 * Math.cos(a);
        const fy = 0.5 + (it.y - 0.5) * 0.7 + 0.14 * Math.sin(a);
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
}

/** Where pans go in a combi oven: its cavity box and rail heights in the fixture frame, or null for another archetype. */
export function combiCavity(id: string, e: Ext): CombiCavity | null {
  const pans = id === "NSF-COMBI-6PAN" ? 6 : id === "NSF-COMBI-10PAN" ? 10 : id === "NSF-COMBI-20PAN-ROLLIN" ? 20 : 0;
  if (!pans) return null;
  const L = combiLayout(pans, e);
  const { x0, x1, y0, y1, z0, z1 } = L.cav;
  return { box: { x0, x1, y0, y1, z0, z1 }, rails: [...L.rails], cx: L.cx, runners: [L.cx - RUNNER, L.cx + RUNNER], zc: z0 + (pans === 20 ? 0.042 : 0.024) + 0.265 };
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
    slab(s, p.ss, 0, xR, top - 0.012, top - sk, L.zDoorFront - sk, L.zDoorFront, 0.0005);
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
      rectPath(sh, cav.x1 - cav.x0, cav.y1 - cav.y0, 0.006, cx, (cav.y0 + cav.y1) / 2, true);
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
      const y0 = base + 0.005;
      const y1 = L.door.y0 - 0.004;
      const [x0, x1] = [L.door.x0 + 0.03, L.door.x1];
      const zf = L.zDoorFront + 0.003;
      slab(s, p.ss, x0, x1, y0, y0 + 0.0012, zFrame, zf, 0.0005);
      slab(s, p.ss, x0, x1, y0, y1, zf - 0.0012, zf, 0.0005);
      for (const x of [x0, x1 - 0.0012]) slab(s, p.ss, x, x + 0.0012, y0, y1, zFrame, zf, 0.0004);
      if (fine2) s.add(p.enamel, cyl(s, 0.005, 0.0008, "y", 12), [x1 - 0.02, y0 + 0.0016, (zFrame + zf) / 2]);
    }
    // Hinge brackets: a pin top and bottom on the frame's left edge, the door's knuckles turning on them.
    const pivot: V3 = [L.door.x0, 0, L.zDoorBack];
    for (const y of [L.door.y1 + 0.0015, Math.max(0.0005, L.door.y0 - 0.0105)]) {
      slab(s, p.ssBright, 0, 0.034, y, y + 0.009, zFrame, L.zDoorBack + 0.009, 0.0015);
      if (fine2) screw(s, p.ssBright, [0.022, y + 0.0045, L.zDoorBack + 0.009], 0.0022);
    }

    /* The cavity: deep-drawn stainless, its runners, the fan and element behind the air baffle. */
    const ct = 0.0015;
    if (fine1) {
      slab(s, p.ss, cav.x0 - ct, cav.x0, cav.y0, cav.y1, cav.wall, cav.z1);
      slab(s, p.ss, cav.x0 - ct, cav.x1 + ct, cav.y1, cav.y1 + ct, cav.wall, cav.z1);
      slab(s, p.ssHeat, cav.x0 - ct, cav.x1 + ct, cav.y0, cav.y1, cav.wall - ct, cav.wall);
      // The floor, its drain at the back under the baffle's foot.
      const dz = cav.z0 + 0.04;
      const fz1 = big ? cav.z1 - 0.035 : cav.z1;
      const fl = roundRect(cav.x1 - cav.x0 + 2 * ct, fz1 - cav.wall, 0, cx, (fz1 + cav.wall) / 2);
      fl.holes.push(circleHole(cx, dz, 0.024));
      const flg = plate(s, fl, ct);
      flg.rotateX(Math.PI / 2);
      s.add(p.ss, flg, [0, cav.y0, 0]);
      if (big) slab(s, p.ss, cav.x0 - ct, cav.x1 + ct, cav.y0 - ct, cav.y0, cav.z1 - 0.008, cav.z1);
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

    /* The door: a stainless frame carrying three panes, its outer one smoked with a black border; hinged left. */
    const door = s.joint("cavity_door", pivot);
    {
      const W = L.door.x1 - L.door.x0;
      const T = L.zDoorFront - L.zDoorBack;
      const [y0, y1] = [L.door.y0, L.door.y1];
      const H = y1 - y0;
      const ym = (y0 + y1) / 2;
      door.add(p.glassSmoke, rbox(door, W, H, 0.006, 0.004), [W / 2, ym, T - 0.003]);
      // Black enamel border fired onto the outer pane's back, hiding the frame.
      if (fine1) {
        const fr = roundRect(W - 0.002, H - 0.002, 0.004, W / 2, ym);
        rectPath(fr, W - 0.085, H - 0.1, 0.01, W / 2 - 0.008, ym, true);
        door.add(p.enamel, plate(door, fr, 0.001), [0, 0, T - 0.0076]);
      }
      // Frame: stiles and rails of formed stainless, the door's visible edges.
      const fz0 = 0.0068;
      const fz1 = T - 0.0082;
      slab(door, p.ss, 0, 0.032, y0, y1, fz0, fz1, 0.002);
      slab(door, p.ss, W - 0.048, W, y0, y1, fz0, fz1, 0.002);
      slab(door, p.ss, 0.03, W - 0.046, y1 - 0.038, y1, fz0, fz1, 0.002);
      slab(door, p.ss, 0.03, W - 0.046, y0, y0 + 0.048, fz0, fz1, 0.002);
      // The middle pane, and the inner pane on its own carrier: it swings out for cleaning between the panes.
      if (fine1) door.add(p.glass, box(W - 0.05, H - 0.06, 0.004), [W / 2 - 0.006, ym, T * 0.5]);
      const ix0 = cav.x0 - 0.034 - L.door.x0;
      const ix1 = cav.x1 + 0.034 - L.door.x0;
      const iy0 = Math.max(cav.y0 - 0.034, y0 + 0.012);
      const iy1 = cav.y1 + 0.034;
      door.add(p.glass, rbox(door, ix1 - ix0, iy1 - iy0, 0.005, 0.003), [(ix0 + ix1) / 2, (iy0 + iy1) / 2, 0.0025]);
      // The door's back sheet round the inner pane (heat-tinted from the cavity side).
      {
        const bs = roundRect(W - 0.002, H - 0.002, 0.003, W / 2, ym);
        rectPath(bs, ix1 - ix0 - 0.03, iy1 - iy0 - 0.03, 0.006, (ix0 + ix1) / 2, (iy0 + iy1) / 2, true);
        door.add(p.ssDull, plate(door, bs, 0.0012), [0, 0, 0.0055]);
      }
      if (fine2) {
        // Clips holding the inner pane, and the condensate lip along its foot.
        for (const yy of [iy0 + 0.08, iy1 - 0.08]) slab(door, p.ssBright, ix1 - 0.012, ix1 + 0.008, yy - 0.012, yy + 0.012, -0.0012, 0.0055, 0.001);
        slab(door, p.ss, ix0 + 0.02, ix1 - 0.02, iy0 - 0.006, iy0 + 0.002, -0.004, 0.0055, 0.001);
        // The latch tongue that the frame's striker catches.
        slab(door, p.ssBright, W - 0.03, W - 0.012, ym - 0.03, ym + 0.03, -0.005, 0.0068, 0.002);
      }
      // Hinge knuckles round the pins.
      for (const yy of [y1 - 0.02, y0 + 0.02]) door.add(p.ssBright, cyl(door, 0.0062, 0.04, "y", 16), [0.0001, yy, 0.0062 - 0.006]);
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

    /* The right-hand column: side panel, cavity wall and controls, all hidden by the cutaway. */
    const cut = s.joint("cutaway", [0, 0, 0]);
    {
      slab(cut, p.ss, xR - sk, xR, base, top - sk, zBack, zFrame, 0.0005);
      if (fine1) slab(cut, p.ss, cav.x1, cav.x1 + ct, cav.y0, cav.y1, cav.wall, cav.z1);
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
      slab(cut, p.bakelite, fx0, fx1, fy0, fy1, P.z - 0.003, P.z, 0.003);
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
      // The hand shower: a spray gun in a holder at the column's foot, its hose reeled inside.
      if (fine2) {
        const hy = pans === 6 ? py0 + 0.06 : py0 + 0.25;
        const hz = P.z - 0.002;
        slab(cut, p.plastic, pcx - 0.024, pcx + 0.024, hy - 0.045, hy + 0.035, hz, hz + 0.005, 0.004);
        tube(cut, p.bakelite, [[pcx, hy - 0.036, hz + 0.016], [pcx, hy + 0.012, hz + 0.016], [pcx, hy + 0.024, hz + 0.012]], 0.0095, 0.012, 14);
        cut.add(p.ssBright, cyl(cut, 0.0085, 0.012, "y", 16, 0.0075), [pcx, hy + 0.029, hz + 0.011]);
        cut.add(p.charcoal, box(0.006, 0.022, 0.008), [pcx, hy - 0.008, hz + 0.027], [0.15, 0, 0]);
        tube(cut, p.rubber, [[pcx, hy - 0.036, hz + 0.016], [pcx, hy - 0.05, hz + 0.008], [pcx, hy - 0.05, hz - 0.002]], 0.005, 0.008, 8);
      }
      // The column's left cheek, where it meets the door's right edge, and the door's striker.
      slab(cut, p.ss, P.x0, P.x0 + sk, py0, py1, zFrame, P.z - 0.004);
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
      s.add(p.poly, cyl(s, 0.025, base - 0.014, "y", 20), [cx, (base + 0.014) / 2, zBack + 0.06]);
      s.add(p.ssBright, torus(s, 0.026, 0.003, "y", Math.PI * 2, 24, 6), [cx, base - 0.02, zBack + 0.06]);
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
  if (!s.fine(2))
    return vessel(s, w, d, 0.028, [
      [lip + dr + 0.008, -depth],
      [lip + 0.0015, -0.006],
      [0, -0.002],
      [0.0015, 0],
      [lip + t, -0.002],
      [lip + dr + t + 0.008, -depth + t],
    ]);
  return vessel(s, w, d, 0.028, [
    [lip + dr + 0.012, -depth],
    [lip + dr + 0.003, -depth + 0.004],
    [lip + 0.0015, -0.006],
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
  // Its load: rolls proofing on two sheet pans, an empty pan above.
  if (fine1) {
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
    if (fine2) s.add(p.enamel, cyl(s, 0.032, 0.0005, "z", 32), [x, pY, pz + 0.0013]);
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
    door.add(p.glass, box(wx1 - wx0 + 0.02, wy1 - wy0 + 0.02, 0.004), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, T - 0.008]);
    if (fine1) door.add(p.glass, box(wx1 - wx0 + 0.02, wy1 - wy0 + 0.02, 0.004), [(wx0 + wx1) / 2, (wy0 + wy1) / 2, 0.01]);
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
  slab(s, p.ssDull, 0.001, w - 0.001, yb, yb + 0.0012, zb, zf);
  slab(s, p.ssDull, 0.001, w - 0.001, yb, h - 0.001, zb, zb + 0.0012);
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

  // The cavity.
  if (fine1) {
    slab(s, p.ssDull, cav.x0 - 0.0012, cav.x0, cav.y0, cav.y1, cav.z0, cav.z1);
    slab(s, p.ssDull, cav.x1, cav.x1 + 0.0012, cav.y0, cav.y1, cav.z0, cav.z1);
    slab(s, p.ssDull, cav.x0, cav.x1, cav.y1, cav.y1 + 0.0012, cav.z0, cav.z1);
    slab(s, p.ssDull, cav.x0, cav.x1, cav.y0, cav.y1, cav.z0 - 0.0012, cav.z0);
    // The ceramic floor over the bottom feed, sealed into its frame.
    slab(s, p.ceramic, cav.x0, cav.x1, cav.y0 - 0.006, cav.y0, cav.z0, cav.z1, fine2 ? 0.001 : 0);
    // Mica waveguide cover on the right wall, the stirrer cover above.
    s.add(fm.mica, new THREE.PlaneGeometry(0.12, 0.09).rotateY(-Math.PI / 2), [cav.x1 - 0.0008, (cav.y0 + cav.y1) / 2 + 0.03, (cav.z0 + cav.z1) / 2 - 0.04]);
    s.add(fm.mica, new THREE.PlaneGeometry(0.16, 0.2).rotateX(Math.PI / 2), [(cav.x0 + cav.x1) / 2, cav.y1 - 0.0008, (cav.z0 + cav.z1) / 2]);
    // The lamp behind its perforated patch, top left.
    s.add(lamp, new THREE.PlaneGeometry(0.05, 0.035).rotateY(Math.PI / 2), [cav.x0 + 0.0004, cav.y1 - 0.04, cav.z1 - 0.12]);
    s.add(fm.perf, new THREE.PlaneGeometry(0.05, 0.035).rotateY(Math.PI / 2), [cav.x0 + 0.0012, cav.y1 - 0.04, cav.z1 - 0.12]);
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
    const m = new THREE.Matrix4();
    for (let i = 0; i < this.mesh.count; i++) {
      this.at((((i * this.pitch + this.off) % this.L) + this.L) % this.L, m);
      this.mesh.setMatrixAt(i, m);
    }
    this.mesh.instanceMatrix.needsUpdate = true;
    this.mesh.computeBoundingSphere();
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
  const fm = mats();
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
  slab(s, p.ssDull, x0, x1, yb, yb + 0.0012, zb, zf);
  slab(s, p.ssDull, x0, x1, yb, h, zb, zb + 0.0012);
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
  if (fine2) for (const x of [cx - 0.09, cx + 0.03]) s.add(p.enamel, cyl(s, 0.027, 0.0005, "z", 32), [x, kY, zf + 0.0004]);
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
  for (const x of [belt.x0 - 0.006, belt.x1 + 0.006]) slab(s, p.ss, x - 0.0008, x + 0.0008, belt.yC - 0.022, yTop + 0.018, belt.zF - 0.03, -0.004, 0.0004);
  slab(s, p.ss, belt.x0 - 0.007, belt.x1 + 0.007, belt.yC - 0.024, belt.yC - 0.0225, zf, -0.004, 0.0004);
  let chain: Chain | null = null;
  if (fine1) {
    chain = new Chain(cx, belt.zF, belt.zR, belt.yC, belt.R + 0.0012, fine2 ? 0.0127 : 0.0254, belt.x1 - belt.x0, p.chrome);
    s.extra(chain.mesh);
    s.tris += chain.tris;
  } else s.add(p.chrome, box(belt.x1 - belt.x0, 0.002, belt.zF - belt.zR), [cx, yTop, (belt.zF + belt.zR) / 2]);
  // Bread riding the belt, toast arriving on the tray: a crust round a crumb face that browns.
  const slice = (at: V3, rot: V3, done: boolean) =>
    s.at(at, rot, () => {
      s.add(done ? fm.crust : fm.crustRaw, rbox(s, 0.1, 0.012, 0.095, 0.005), [0, 0, 0]);
      if (fine2) for (const sy of [-1, 1]) s.add(done ? fm.toast : fm.crumb, rbox(s, 0.09, 0.001, 0.085, 0.0004), [0, sy * 0.0058, 0.002]);
    });
  if (fine1) {
    slice([cx - 0.06, yTop + 0.0075, -0.04], [0, 0.05, 0], false);
    slice([cx + 0.055, yTop + 0.0075, -0.2], [0, -0.04, 0], true);
  }
  // The return chute from the back, and the front discharge tray with its raised sides.
  if (fine1) {
    const g = box(belt.x1 - belt.x0 + 0.01, 0.0012, Math.hypot(belt.zR - zf, 0.1));
    s.add(p.ss, g, [cx, (0.152 + dis.y0 + 0.008) / 2, (belt.zR + zf) / 2], [Math.atan2(0.152 - dis.y0 - 0.008, zf - belt.zR), 0, 0]);
  }
  slab(s, p.ss, belt.x0 - 0.005, belt.x1 + 0.005, dis.y0, dis.y0 + 0.0015, zf - 0.03, -0.004, 0.0005);
  slab(s, p.ss, belt.x0 - 0.005, belt.x1 + 0.005, dis.y0, dis.y0 + 0.02, -0.0055, -0.004, 0.0005);
  for (const x of [belt.x0 - 0.005, belt.x1 + 0.0035]) slab(s, p.ss, x, x + 0.0015, dis.y0, dis.y0 + 0.02, zf, -0.004, 0.0005);
  if (fine1) slice([cx - 0.02, dis.y0 + 0.0085, -0.05], [0.04, 0.3, 0], true);
  // The crumb tray, its pull tab under the discharge.
  slab(s, p.ss, belt.x0, belt.x1, yb + 0.002, yb + 0.016, zb + 0.02, zf + 0.0015, 0.0006);
  slab(s, p.ssBright, cx - 0.03, cx + 0.03, yb + 0.004, yb + 0.014, zf + 0.0015, zf + 0.009, 0.002);
  dataPlate(s, [x1, 0.1, zf - 0.08], [0, Math.PI / 2, 0], 0.07, 0.04);
  if (fine1) tube(s, p.rubber, [[cx + 0.1, 0.08, zb], [cx + 0.1, 0.08, -d + 0.0005], [cx + 0.1, yb + 0.004, -d + 0.0005]], 0.004, 0.02, 8);

  return (b) => {
    let speed = 0.5;
    setGlow(glow, 0);
    const sj = sprockets.map((n) => b.joints.get(n));
    turnJoint(b.joints.get("speed_knob"), "z", -speed * Math.PI * 1.5);
    setLamp(pilot, 0);
    return {
      glows: [glow],
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
  slab(s, p.ssDull, 0.001, w - 0.001, yb, yb + 0.0012, zb, zf);
  slab(s, p.ssDull, 0.001, w - 0.001, yb, h - 0.001, zb, zb + 0.0012);
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
      dr.add(fm.fries, foodTop(dr, bw - 0.01, bd - 0.01, 0.01, 0.016, 0.07, 4), [cx, by0 + 0.072, pzc]);
      if (fine2)
        for (let i = 0; i < 70; i++) {
          const L = 0.05 + fr() * 0.04;
          dr.add(fm.fries, box(0.0085, 0.0085, L), [cx + (fr() - 0.5) * (bw - 0.06), by0 + 0.07 + fr() * 0.035, pzc + (fr() - 0.5) * (bd - 0.08)], [(fr() - 0.5) * 0.9, fr() * 3.14, (fr() - 0.5) * 0.6]);
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
    s.add(p.water, vessel(s, tw, td, 0.03, [[i0 + 0.006, 0.0035], [i0, 0.006], [i1, yw]], [true, false]), [cx, 0, tzc]);
    // The meniscus climbing the wall: a narrow band that reads as the water line.
    if (fine1) s.add(p.water, vessel(s, tw, td, 0.03, [[i1 + 0.0012, yw - 0.0035], [i1 + 0.0006, yw + 0.001], [i1 - 0.0004, yw + 0.0022]], [false, false]), [cx, 0, tzc]);
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
  const surf = new THREE.Mesh(surfG, p.water);
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
    s.add(fm.sousDisplay, new THREE.PlaneGeometry(0.036, 0.016), [cx, 0.312, ccz + 0.0408]);
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
  // The bags: a steak standing in its bag clipped to the side, a salmon fillet on the floor.
  if (fine1) {
    const bag = (food: () => void, bw: number, bh: number, bt: number) => {
      food();
      s.add(fm.bag, rbox(s, bw, bh, bt, Math.min(0.01, bt / 2 - 1e-4)), [0, 0, 0]);
      s.add(fm.bag, box(bw + 0.004, 0.035, 0.0012), [0, bh / 2 + 0.016, 0]);
    };
    s.at([0.064, 0.115, -0.2], [0, Math.PI / 2, 0.06], () =>
      bag(
        () => {
          const st = new THREE.Shape();
          st.moveTo(-0.07, -0.045);
          st.bezierCurveTo(-0.02, -0.06, 0.05, -0.055, 0.075, -0.02);
          st.bezierCurveTo(0.085, 0.02, 0.06, 0.05, 0.02, 0.055);
          st.bezierCurveTo(-0.03, 0.06, -0.08, 0.03, -0.07, -0.045);
          const g = plate(s, st, 0.026, s.fine(2) ? 0.004 : 0);
          g.translate(0, 0, -0.013);
          s.add(fm.beef, g, [0, 0, 0]);
          s.add(fm.fat, rbox(s, 0.11, 0.01, 0.024, 0.004), [0.005, 0.052, 0], [0, 0, -0.1]);
        },
        0.18,
        0.14,
        0.034,
      ),
    );
    // Its clip over the tub's wall, pinching the bag's seal.
    slab(s, p.enamel, 0.004, 0.022, rim - 0.035, rim + 0.007, -0.212, -0.188, 0.002);
    slab(s, p.enamel, 0.004, 0.05, rim + 0.001, rim + 0.007, -0.212, -0.188, 0.002);
    slab(s, p.enamel, 0.044, 0.05, rim - 0.03, rim + 0.007, -0.212, -0.188, 0.002);
    s.at([cx + 0.03, 0.03, -0.24], [0, 0.5, 0], () =>
      bag(
        () => {
          s.add(fm.salmon, rbox(s, 0.16, 0.024, 0.075, 0.01), [0, 0, 0]);
          if (s.fine(2)) for (let i = 0; i < 6; i++) s.add(fm.fat, box(0.003, 0.0006, 0.06), [-0.06 + i * 0.024, 0.0121, 0], [0, 0.5, 0]);
        },
        0.2,
        0.036,
        0.11,
      ),
    );
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

/** A pan's food at `level` below its rim (negative y), lumped: the food meets the drafted wall where it is. */
function panFood(s: Shop, w: number, d: number, depth: number, level: number, amp: number, lambda: number, seed: number): THREE.BufferGeometry {
  const lip = 0.008;
  const t = 0.0008;
  const dr = depth * 0.08;
  const u = Math.min(1, Math.max(0, (-level - 0.002) / (depth - 0.006)));
  const inset = s.fine(1) ? lip + t + u * (dr + 0.003) + 0.0015 : 0.004;
  const g = foodTop(s, w - 2 * inset, d - 2 * inset, Math.max(0.004, 0.028 - inset), amp, lambda, seed);
  g.translate(0, level, 0);
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
  slab(s, p.ssDull, 0.001, w - 0.001, yb, yb + 0.0012, zb, zf);
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
  const py = top + 0.003;
  const fill = (x: number, z: number, pw: number, pd: number, depth: number, m: Mat, level: number, amp: number, lambda: number, seed: number) => {
    if (fine1) s.add(p.ss, gnPan(s, pw, pd, depth), [x, py, z]);
    s.add(m, panFood(s, pw, pd, depth, level, amp, lambda, seed), [x, py, z]);
  };
  fill(wells[0]!, zc, 0.325, 0.53, 0.1, fm.mac, -0.028, 0.012, 0.05, 1);
  fill(wells[1]!, zc - 0.133, 0.325, 0.265, 0.1, fm.greens, -0.03, 0.01, 0.04, 2);
  fill(wells[1]!, zc + 0.133, 0.325, 0.265, 0.1, fm.rice, -0.026, 0.008, 0.06, 3);
  fill(wells[2]!, zc, 0.325, 0.53, 0.1, fm.chili, -0.026, 0.004, 0.08, 4);
  for (const [k, m] of [fm.mash, fm.corn, fm.gravy].entries()) fill(wells[3]!, zc + (k - 1) * 0.177, 0.325, 0.176, 0.1, m, k === 2 ? -0.03 : -0.025, k === 0 ? 0.018 : k === 1 ? 0.006 : 0.001, k === 0 ? 0.05 : 0.025, 5 + k);
  // Adaptor bars across the wells between the fractional pans.
  for (const z of [zc]) slab(s, p.ss, wells[1]! - 0.16, wells[1]! + 0.16, top, top + 0.003, z - 0.0015, z + 0.0015, r2(0.0005));
  for (const z of [zc - 0.0885, zc + 0.0885]) slab(s, p.ss, wells[3]! - 0.16, wells[3]! + 0.16, top, top + 0.003, z - 0.0015, z + 0.0015, r2(0.0005));
  if (fine2) {
    // Green beans, kidney beans in the chili.
    const r = rng(5);
    for (let i = 0; i < 46; i++) s.add(fm.greens, cyl(s, 0.0035, 0.05 + r() * 0.03, "x", 6), [wells[1]! + (r() - 0.5) * 0.25, py - 0.022 + r() * 0.006, zc - 0.133 + (r() - 0.5) * 0.18], [0, r() * 3.14, (r() - 0.5) * 0.3]);
    for (let i = 0; i < 40; i++) s.add(fm.bean, new THREE.SphereGeometry(0.006, 6, 4), [wells[2]! + (r() - 0.5) * 0.26, py - 0.021 + r() * 0.003, zc + (r() - 0.5) * 0.44], [0, r() * 3, 0], [1.5, 0.7, 0.9]);
    // A serving spoon in the mac, a ladle in the chili.
    const bowl: [number, number][] = [];
    for (let i = 0; i <= 8; i++) bowl.push([Math.sin((i / 8) * Math.PI / 2) * 0.042, -Math.cos((i / 8) * Math.PI / 2) * 0.034]);
    bowl.push([0.0435, 0.0012]);
    for (let i = 8; i >= 0; i--) bowl.push([Math.sin((i / 8) * Math.PI / 2) * 0.0405, -Math.cos((i / 8) * Math.PI / 2) * 0.0325]);
    // The spoon: a shallow oval bowl half buried, its handle up over the pan's back rim.
    const sx = wells[0]! + 0.02;
    const sz = zc + 0.1;
    s.add(p.ssBright, lathe(s, bowl, 24), [sx, py - 0.024, sz], [-0.12, 0, 0], [0.75, 0.45, 1.05]);
    tube(s, p.ssBright, [[sx, py - 0.022, sz + 0.042], [sx + 0.006, py - 0.005, sz + 0.1], [sx + 0.012, py + 0.016, zc + 0.262], [sx + 0.013, py + 0.008, zc + 0.3]], 0.0042, 0.02, 8);
    const lx = wells[2]! - 0.06;
    const lz = zc - 0.1;
    s.add(p.ssBright, lathe(s, bowl, 24), [lx, py - 0.018, lz], [0.18, 0, -0.12]);
    tube(s, p.ssBright, [[lx + 0.042, py - 0.012, lz], [wells[2]! + 0.12, py + 0.026, lz - 0.03], [wells[2]! + 0.163, py + 0.03, lz - 0.035], [wells[2]! + 0.168, py + 0.004, lz - 0.035]], 0.0045, 0.008, 8);
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
    if (fine2) s.add(p.enamel, cyl(s, 0.029, 0.0005, "z", 32), [x, kY, zf + 0.0004]);
  }
  dataPlate(s, [w - 0.1, 0.68, zf], [0, 0, 0], 0.08, 0.045);
  slab(s, p.hdpe, 0.005, w - 0.005, 0.8, 0.8127, zf + 0.002, -0.004, fine1 ? 0.003 : 0);
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
  "NSF-MICROWAVE-COMM": microwave,
  "NSF-TOASTER-CONVEYOR": conveyorToaster,
  "NSF-AIRFRYER-COMM": airFryer,
  "NSF-SOUSVIDE-BATH": sousVide,
  "NSF-STEAMTABLE-HOT-60": steamTable,
  "NSF-BOOSTER-HEATER-36KW": boosterHeater,
});
