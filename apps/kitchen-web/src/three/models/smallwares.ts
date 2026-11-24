/**
 * Model family: smallwares and food. Cookware a line cook reaches for (skillets, stockpots, hotel pans), the utensils
 * beside it (a forged chef's knife, tongs, a fish turner, a balloon whisk, a ladle, the deglazing bottle), and the food
 * the page cooks, whose colour and shape are live.
 *
 * Every vessel is a real hollow body: a profile with a wall thickness turned or lofted round its axis, a rim, a flat
 * base. Long utensils lie along the depth axis, handle toward the front (+z), business end toward the back.
 *
 * Live quantities, `hooks.set(name, value)` (default in brackets):
 *   NSF-STOCKPOT-*            level 0..1 liquid height [0.7]; boil 0..1 [0] (tick rolls the surface, raises bubbles,
 *                             steam above about 0.6). open(frac) lifts and tips the lid (joint "pot_lid_lift").
 *   NSF-HOTEL-PAN-*           fill 0..1 food surface height [0.6]; food 0 braised meat in sauce, 1 rice, 2 vegetables,
 *                             3 soup [0]. Each food rides its own joint ("food_0".."food_3").
 *   NSF-TONGS-12IN            open(frac) opens the upper arm up to 30° about the hinge pin (joint "tongs_pinch").
 *   NSF-DEGLAZE-BOTTLE        level 0..1 wine in the bottle, 1 = a full 750 mL [0.7]. The lip is at (w/2, h, -d/2).
 *   ARCH-FOOD-STEAK-RIBEYE    browning 0..1 Maillard crust on both faces (stove kernel maillard_progress) [0];
 *                             shrink 0..0.3 (1 - cook_yield_ppm/1e6): contracts across, thickens [0]; core °C
 *                             (optional): the grey band on the cut edges and the colour inside it. Joint "steak".
 *   ARCH-FOOD-SOURDOUGH-BOULE spring 1..1.85 (oven kernel oven_spring_ratio): volume, mostly height, score opening into
 *                             an ear [1.85]; crust 0..1 (crust_set_frac, or maillard_progress): pale dough to deep
 *                             brown [0.8]. Joint "boule".
 *   ARCH-FOOD-EGG-FRIED       set 0..1: albumen clear to opaque white, lacy brown edges near 1 [0.9].
 *   ARCH-FOOD-FRY-LOAD        browning 0..1: blanched pale to deep golden, tips and corners first [0.75].
 *   ARCH-FOOD-BRAISED-PAN-PORTION  no live quantities: a fixed portion of the hotel pans' braise.
 *
 * Detail: 2 is for close-ups (rivets, coil springs, steam, bubbles, herbs); 1 drops the smallest parts and halves the
 * facets; 0 keeps silhouettes and the foods' painted colour only.
 *
 * Exports for other families and chapters: skillet() builds a pan centred on its axis, its handle along +z;
 * hotelPan() a gastronorm pan with optional food.
 */
import * as THREE from "three";
import { mergeVertices } from "three/examples/jsm/utils/BufferGeometryUtils.js";
import { register, type Builder, type Hooks } from "./index";
import { palette } from "./materials";
import { bar, cyl, lathe, plate, rbox, torus, ball, turnJoint, Shop, type V3, type Built } from "./shop";

/* ------------------------------------------------------------------------------------------------------------ */
/* Small maths, seeded noise, canvases                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

type P2 = [number, number];
const PI = Math.PI;
const TAU = Math.PI * 2;
const clamp = (v: number, a: number, b: number) => Math.max(a, Math.min(b, v));
const clamp01 = (v: number) => Math.max(0, Math.min(1, v));
const lerp = (a: number, b: number, t: number) => a + (b - a) * t;
const smooth = (e0: number, e1: number, x: number) => {
  const t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3 - 2 * t);
};

function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** Smooth value noise on a 256-periodic lattice, 0..1. */
function noise2(seed: number): (x: number, y: number) => number {
  const r = rng(seed);
  const N = 256;
  const v = new Float32Array(N * N);
  for (let i = 0; i < v.length; i++) v[i] = r();
  return (x, y) => {
    const xi = Math.floor(x);
    const yi = Math.floor(y);
    const fx = x - xi;
    const fy = y - yi;
    const ux = fx * fx * (3 - 2 * fx);
    const uy = fy * fy * (3 - 2 * fy);
    const x0 = xi & 255;
    const y0 = yi & 255;
    const x1 = (x0 + 1) & 255;
    const y1 = (y0 + 1) & 255;
    const a = v[y0 * N + x0];
    const b = v[y0 * N + x1];
    const c = v[y1 * N + x0];
    const d = v[y1 * N + x1];
    return a + (b - a) * ux + (c - a) * uy + (a - b - c + d) * ux * uy;
  };
}
/** Fractal sum of `oct` octaves of `n`, 0..1. */
function fbm(n: (x: number, y: number) => number, x: number, y: number, oct = 4): number {
  let s = 0;
  let a = 0.5;
  let f = 1;
  let norm = 0;
  for (let i = 0; i < oct; i++) {
    s += a * n(x * f + i * 17.3, y * f - i * 11.1);
    norm += a;
    a *= 0.5;
    f *= 2.03;
  }
  return s / norm;
}

function canvas(w: number, h: number): [HTMLCanvasElement, CanvasRenderingContext2D] {
  const c = document.createElement("canvas");
  c.width = w;
  c.height = h;
  return [c, c.getContext("2d", { willReadFrequently: true })!];
}
function canvasTex(c: HTMLCanvasElement, srgb: boolean, wrap = false): THREE.CanvasTexture {
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.wrapS = t.wrapT = wrap ? THREE.RepeatWrapping : THREE.ClampToEdgeWrapping;
  t.anisotropy = 8;
  t.needsUpdate = true;
  return t;
}

/** A canvas filled pixel by pixel: `f(x, y)` returns [r, g, b] 0..255. */
function paint(w: number, h: number, f: (x: number, y: number) => readonly [number, number, number]): HTMLCanvasElement {
  const [c, x] = canvas(w, h);
  const img = x.createImageData(w, h);
  const d = img.data;
  for (let j = 0; j < h; j++)
    for (let i = 0; i < w; i++) {
      const [r, g, b] = f(i, j);
      const k = (j * w + i) * 4;
      d[k] = r;
      d[k + 1] = g;
      d[k + 2] = b;
      d[k + 3] = 255;
    }
  x.putImageData(img, 0, 0);
  return c;
}

function std(p: THREE.MeshStandardMaterialParameters, uv: "box" | "own" = "box", uvScale?: number): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial(p);
  m.userData.uv = uv;
  if (uvScale) m.userData.uvScale = uvScale;
  return m;
}
function phys(p: THREE.MeshPhysicalMaterialParameters, uv: "box" | "own" = "box", uvScale?: number): THREE.MeshPhysicalMaterial {
  const m = new THREE.MeshPhysicalMaterial(p);
  m.userData.uv = uv;
  if (uvScale) m.userData.uvScale = uvScale;
  return m;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Finishes of the smallwares, made once and shared                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

interface Finishes {
  /** Seasoned cast iron: black polymerised oil over sand-cast grain, a satin sheen. */
  iron: THREE.MeshStandardMaterial;
  /** Seasoned carbon steel: blue-black with the bronze and brown patches seasoning leaves. */
  carbon: THREE.MeshStandardMaterial;
  /** PTFE non-stick: charcoal with a fine reinforcing fleck. */
  ptfe: THREE.MeshStandardMaterial;
  /** Hard-anodised aluminium: dark grey satin. */
  anod: THREE.MeshStandardMaterial;
  /** Red silicone handle sleeve. */
  silicone: THREE.MeshStandardMaterial;
  /** Black POM (acetal) knife scales and handles: dense, slightly glossy. */
  pom: THREE.MeshStandardMaterial;
  /** Knife steel: satin flats ground lengthwise. */
  blade: THREE.MeshPhysicalMaterial;
  /** Polished edge bevel and bolster. */
  edge: THREE.MeshStandardMaterial;
  /** Dark green bottle glass. */
  bottle: THREE.MeshPhysicalMaterial;
  /** Red wine as seen in the bottle: nearly black. */
  wine: THREE.MeshPhysicalMaterial;
  /** Foil capsule, deep burgundy. */
  foil: THREE.MeshStandardMaterial;
  /** Plain cream paper label. */
  label: THREE.MeshStandardMaterial;
  /** Boiling-water bubbles. */
  bubble: THREE.MeshPhysicalMaterial;
  /** Stock in a pot: water tinted by what is cooking in it. */
  stock: THREE.MeshPhysicalMaterial;
  /** Spun satin stainless, the inside of a pan or pot: fine concentric lines from the spinning lathe. */
  spun: THREE.MeshStandardMaterial;
  /** The same, spun brighter, for polished outsides and lids. */
  spunBright: THREE.MeshStandardMaterial;
  /** Spun aluminium, natural finish. */
  spunAlu: THREE.MeshStandardMaterial;
  /** Deep-drawn satin stainless (hotel pans): faint lines round the walls, duller than spun. */
  satin: THREE.MeshStandardMaterial;
}
let fin: Finishes | null = null;

function finishes(): Finishes {
  if (fin) return fin;
  // Cast iron: dark speckled colour with a sandy bump, repeating every 80 mm.
  const ni = noise2(4401);
  const ironC = paint(256, 256, (x, y) => {
    const n = fbm(ni, x / 9, y / 9, 3);
    const g = fbm(ni, x / 1.3, y / 1.3, 2);
    const v = 28 + n * 12 + (g - 0.5) * 14;
    return [v + 2, v + 1, v];
  });
  const ironB = paint(256, 256, (x, y) => {
    const g = fbm(ni, x / 1.1 + 40, y / 1.1, 2);
    const v = 128 + (g - 0.5) * 200;
    return [v, v, v];
  });
  const iron = std({ color: 0xffffff, metalness: 0.42, roughness: 0.52, map: canvasTex(ironC, true, true), bumpMap: canvasTex(ironB, false, true), bumpScale: 0.9 }, "box", 1 / 0.08);
  // Carbon steel: blue-black, bronze-brown where the seasoning is thinner, 160 mm repeat.
  const nc = noise2(9127);
  const carbC = paint(256, 256, (x, y) => {
    const n = fbm(nc, x / 40, y / 40, 4);
    const b = smooth(0.52, 0.78, n) * 0.8;
    const f = (fbm(nc, x / 2, y / 2, 2) - 0.5) * 8;
    return [lerp(30, 66, b) + f, lerp(32, 48, b) + f, lerp(40, 34, b) + f];
  });
  const carbon = std({ color: 0xffffff, metalness: 0.68, roughness: 0.42, map: canvasTex(carbC, true, true) }, "box", 1 / 0.16);
  // PTFE: charcoal with the light fleck of the reinforced coating.
  const r = rng(77);
  const [pc, px] = canvas(256, 256);
  px.fillStyle = "#2a2b2e";
  px.fillRect(0, 0, 256, 256);
  for (let i = 0; i < 2600; i++) {
    const v = 60 + Math.floor(r() * 90);
    px.fillStyle = `rgba(${v},${v},${v + 4},${(0.2 + r() * 0.5).toFixed(2)})`;
    px.fillRect(r() * 256, r() * 256, 0.6 + r(), 0.6 + r());
  }
  const ptfe = std({ color: 0xffffff, metalness: 0.12, roughness: 0.36, map: canvasTex(pc, true, true) }, "box", 1 / 0.05);
  // Blade steel: fine lengthwise grind lines, its anisotropy along them (v in the box projection of the flats).
  const rb = rng(3133);
  const [bc, bx] = canvas(256, 256);
  bx.fillStyle = "#dfe2e6";
  bx.fillRect(0, 0, 256, 256);
  for (let i = 0; i < 1800; i++) {
    const v = 190 + Math.floor(rb() * 60);
    bx.fillStyle = `rgba(${v},${v + 1},${v + 3},${(0.1 + rb() * 0.2).toFixed(2)})`;
    bx.fillRect(rb() * 256, rb() * 256, 1, 20 + rb() * 120);
  }
  const blade = phys({ color: 0xd4d8dc, metalness: 1, roughness: 0.24, map: canvasTex(bc, true, true), anisotropy: 0.5, anisotropyRotation: PI / 2 }, "box", 1 / 0.2);
  const bottle = phys({ color: 0x2e5a24, metalness: 0, roughness: 0.06, transparent: true, opacity: 0.5, clearcoat: 1, clearcoatRoughness: 0.04, envMapIntensity: 1.5, side: THREE.DoubleSide, depthWrite: false }, "own");
  const wine = phys({ color: 0x1c0507, metalness: 0, roughness: 0.08, clearcoat: 1, clearcoatRoughness: 0.05 }, "own");
  // Label paper: faint fibre.
  const nl = noise2(555);
  const labC = paint(128, 128, (x, y) => {
    const v = (fbm(nl, x / 3, y / 3, 3) - 0.5) * 14;
    return [236 + v, 226 + v, 200 + v];
  });
  // Spun finish: lines running all the way round (u), varying across them (v), 50 mm of profile per repeat.
  const rs = rng(4242);
  const [sc, sx] = canvas(16, 512);
  const [sr, srx] = canvas(16, 512);
  for (let y = 0; y < 512; y++) {
    const v = 214 + Math.floor(rs() * 32);
    sx.fillStyle = `rgb(${v},${v + 1},${v + 3})`;
    sx.fillRect(0, y, 16, 1);
    const q = 118 + Math.floor(rs() * 48);
    srx.fillStyle = `rgb(${q},${q},${q})`;
    srx.fillRect(0, y, 16, 1);
  }
  const spunMap = canvasTex(sc, true, true);
  const spunRough = canvasTex(sr, false, true);
  const spun = std({ color: 0xc9cdd2, metalness: 1, roughness: 0.3, map: spunMap, roughnessMap: spunRough }, "own", 1 / 0.05);
  const spunBright = std({ color: 0xdde0e4, metalness: 1, roughness: 0.16, map: spunMap }, "own", 1 / 0.05);
  const spunAlu = std({ color: 0xc8ccd0, metalness: 0.92, roughness: 0.4, map: spunMap, roughnessMap: spunRough }, "own", 1 / 0.05);
  const satin = std({ color: 0xc2c6cb, metalness: 0.95, roughness: 0.36, map: spunMap, roughnessMap: spunRough }, "own", 1 / 0.12);
  fin = {
    satin,
    spun,
    spunBright,
    spunAlu,
    iron,
    carbon,
    ptfe,
    anod: std({ color: 0x45474b, metalness: 0.55, roughness: 0.48 }),
    silicone: std({ color: 0x8e1914, metalness: 0, roughness: 0.6 }),
    pom: std({ color: 0x151617, metalness: 0, roughness: 0.32 }),
    blade,
    edge: std({ color: 0xe4e7ea, metalness: 1, roughness: 0.1 }),
    bottle,
    wine,
    foil: std({ color: 0x5a1220, metalness: 0.85, roughness: 0.32 }),
    label: std({ color: 0xffffff, metalness: 0, roughness: 0.86, map: canvasTex(labC, true, true) }, "box", 1 / 0.04),
    bubble: phys({ color: 0xf4f8fa, metalness: 0, roughness: 0.04, transparent: true, opacity: 0.42, clearcoat: 1, depthWrite: false }, "own"),
    stock: phys({ color: 0x8a5a1e, metalness: 0, roughness: 0.04, transparent: true, opacity: 0.86, clearcoat: 1, clearcoatRoughness: 0.03, envMapIntensity: 1.8, depthWrite: false }, "own"),
  };
  return fin;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Geometry: rings, revolves, lofts                                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

interface RingOpts {
  /** Join each ring's last point to its first (default true). */
  closed?: boolean;
  /** Closed, but with the first point repeated at the end so that UVs can wrap; its normals are averaged back. */
  seam?: boolean;
  /** Close the first and last rings with flat fans on their own vertices. */
  caps?: boolean;
  /** UVs per ring point. */
  uv?: (j: number, i: number) => P2;
}

/**
 * Quads between consecutive rings of points (each a flat [x, y, z, ...] list of the same length). A profile turned
 * counter-clockwise in (radius, height), each ring running from +x toward +z, faces outward; `outward` fixes any
 * other closed solid.
 */
function ringMesh(rings: number[][], o: RingOpts = {}): THREE.BufferGeometry {
  if (o.seam) {
    const g = ringMesh(
      rings.map((r) => [...r, r[0], r[1], r[2]]),
      { ...o, seam: false, closed: false },
    );
    const nn = rings[0].length / 3 + 1;
    const nr = g.attributes.normal as THREE.BufferAttribute;
    const v = new THREE.Vector3();
    for (let j = 0; j < rings.length; j++) {
      const a = j * nn;
      const b = j * nn + nn - 1;
      v.set(nr.getX(a) + nr.getX(b), nr.getY(a) + nr.getY(b), nr.getZ(a) + nr.getZ(b)).normalize();
      nr.setXYZ(a, v.x, v.y, v.z);
      nr.setXYZ(b, v.x, v.y, v.z);
    }
    return g;
  }
  const closed = o.closed ?? true;
  const m = rings.length;
  const n = rings[0].length / 3;
  const capV = o.caps ? 2 * (n + 1) : 0;
  const pos = new Float32Array((m * n + capV) * 3);
  const uvs = new Float32Array((m * n + capV) * 2);
  for (let j = 0; j < m; j++) {
    pos.set(rings[j], j * n * 3);
    for (let i = 0; i < n; i++) {
      const [u, v] = o.uv ? o.uv(j, i) : [i / n, j / Math.max(1, m - 1)];
      uvs[(j * n + i) * 2] = u;
      uvs[(j * n + i) * 2 + 1] = v;
    }
  }
  const idx: number[] = [];
  const lim = closed ? n : n - 1;
  for (let j = 0; j < m - 1; j++)
    for (let i = 0; i < lim; i++) {
      const i1 = (i + 1) % n;
      const a = j * n + i;
      const b = j * n + i1;
      const c = (j + 1) * n + i1;
      const d = (j + 1) * n + i;
      idx.push(a, d, b, b, d, c);
    }
  if (o.caps) {
    // Each cap: a centre and its ring copied, so that its flat normal does not bend the sides'.
    let base = m * n;
    for (const [j, flip] of [
      [0, false],
      [m - 1, true],
    ] as const) {
      const ring = rings[j];
      let cx = 0;
      let cy = 0;
      let cz = 0;
      for (let i = 0; i < n; i++) {
        cx += ring[i * 3] / n;
        cy += ring[i * 3 + 1] / n;
        cz += ring[i * 3 + 2] / n;
      }
      pos.set([cx, cy, cz], base * 3);
      pos.set(ring, (base + 1) * 3);
      for (let i = 0; i < lim; i++) {
        const i1 = (i + 1) % n;
        if (!flip) idx.push(base, base + 1 + i, base + 1 + i1);
        else idx.push(base, base + 1 + i1, base + 1 + i);
      }
      base += n + 1;
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  g.setAttribute("uv", new THREE.BufferAttribute(uvs, 2));
  g.setIndex(idx);
  g.computeVertexNormals();
  return g;
}

/** Make a closed solid face outward: flip every triangle if its signed volume is negative. */
function outward(g: THREE.BufferGeometry): THREE.BufferGeometry {
  const p = g.attributes.position;
  const ix = g.index!;
  let vol = 0;
  const a = new THREE.Vector3();
  const b = new THREE.Vector3();
  const c = new THREE.Vector3();
  for (let t = 0; t < ix.count; t += 3) {
    a.fromBufferAttribute(p, ix.getX(t));
    b.fromBufferAttribute(p, ix.getX(t + 1));
    c.fromBufferAttribute(p, ix.getX(t + 2));
    vol += a.dot(b.cross(c));
  }
  if (vol < 0) {
    for (let t = 0; t < ix.count; t += 3) {
      const k = ix.getX(t + 1);
      ix.setX(t + 1, ix.getX(t + 2));
      ix.setX(t + 2, k);
    }
    ix.needsUpdate = true;
    g.computeVertexNormals();
  }
  return g;
}

/**
 * A solid of revolution about +y through the profile [radius, height], turned counter-clockwise for outward faces.
 * `warp(theta, r, y)` may move any point (a pour spout, a dent). Repeat a profile point for a crisp edge.
 */
function revolve(prof: readonly P2[], n: number, warp?: (th: number, r: number, y: number) => P2, uvOf?: THREE.Material, polar = false): THREE.BufferGeometry {
  const rings = prof.map(([r, y]) => {
    const ring: number[] = [];
    for (let i = 0; i < n; i++) {
      const th = (i / n) * TAU;
      const [rr, yy] = warp ? warp(th, r, y) : [r, y];
      ring.push(rr * Math.cos(th), yy, rr * Math.sin(th));
    }
    return ring;
  });
  if (!uvOf) return ringMesh(rings);
  // Turned-finish UVs in metres: u round the axis (a whole number of texture repeats, so that it wraps), v along
  // the profile. Brushed grain then runs round the vessel as a spun finish does.
  const scale = (uvOf.userData.uvScale as number | undefined) ?? 1;
  const rep = ((uvOf as THREE.MeshStandardMaterial).map?.repeat.x ?? 1) || 1;
  const rRef = Math.max(...prof.map(([r]) => r));
  const U = Math.max(1, Math.round(TAU * rRef * scale * rep)) / rep;
  const arcLen = [0];
  for (let k = 1; k < prof.length; k++) arcLen.push(arcLen[k - 1] + Math.hypot(prof[k][0] - prof[k - 1][0], prof[k][1] - prof[k - 1][1]));
  if (polar) {
    // Isotropic finishes: unwrapped flat by arc length from the nearer end of the profile (the base centre outside,
    // the floor centre inside), which keeps distances true on a shallow vessel and has no seam.
    const total = arcLen[arcLen.length - 1];
    const startsOnAxis = prof[0][0] < 1e-4;
    const endsOnAxis = prof[prof.length - 1][0] < 1e-4;
    return ringMesh(rings, {
      uv: (j, i) => {
        const fromEnd = endsOnAxis && (!startsOnAxis || arcLen[j] > total / 2);
        const sd = (fromEnd ? total - arcLen[j] : arcLen[j]) * scale;
        const th = (i / n) * TAU;
        return [sd * Math.cos(th) + (fromEnd ? 0.37 : 0), sd * Math.sin(th)];
      },
    });
  }
  return ringMesh(rings, { seam: true, uv: (j, i) => [(i / n) * U, arcLen[j] * scale] });
}

const ownCache = new Map<THREE.Material, THREE.Material>();
/** The same finish, keeping the UVs its geometry brings (turned parts), sharing every texture. */
function own<M extends THREE.Material>(m: M): M {
  if (m.userData.uv === "own") return m;
  let o = ownCache.get(m);
  if (!o) {
    o = m.clone();
    o.userData = { ...m.userData, uv: "own" };
    ownCache.set(m, o);
  }
  return o as M;
}

/** Add a turned part with its own spun UVs. */
function turned(s: Shop, mat: THREE.Material, prof0: readonly P2[], n0: number, warp?: (th: number, r: number, y: number) => P2, at: V3 = [0, 0, 0]): void {
  // Farther off, fewer facets round and fewer points down the section.
  const prof = lodProf(s, prof0);
  const n = s.detail === 2 ? n0 : Math.max(s.detail === 1 ? 12 : 10, Math.round(n0 * (s.detail === 1 ? 0.5 : 0.45)));
  // Spun finishes carry turned UVs; isotropic ones (iron, patina, coatings) unwrap flat from the axis.
  if (mat.userData.uv === "own") s.add(mat, revolve(prof, n, warp, mat), at);
  else {
    const m = own(mat);
    s.add(m, revolve(prof, n, warp, m, true), at);
  }
}

/** An arc of `n` segments about (cx, cy), radius r, from angle a0 to a1, start point included when `first`. */
function arc(out: P2[], cx: number, cy: number, r: number, a0: number, a1: number, n: number, first = false): void {
  for (let k = first ? 0 : 1; k <= n; k++) {
    const a = a0 + ((a1 - a0) * k) / n;
    out.push([cx + r * Math.cos(a), cy + r * Math.sin(a)]);
  }
}

/**
 * A rounded rectangle's outline in the xz plane at height y, half-sizes A (x) and B (z), corner radius rc: k + 1
 * points per corner and mx, mz segments along the x and z sides, from +x toward +z, the same count whatever rc is.
 */
function rrRing(A: number, B: number, rc: number, y: number, k: number, mx: number, mz: number): number[] {
  const out: number[] = [];
  const r = Math.max(0, Math.min(rc, A, B));
  const ax = A - r;
  const bz = B - r;
  const corners: [number, number, number][] = [
    [ax, bz, 0],
    [-ax, bz, PI / 2],
    [-ax, -bz, PI],
    [ax, -bz, (3 * PI) / 2],
  ];
  const sideStart: P2[] = [
    [A, -bz],
    [ax, B],
    [-A, bz],
    [-ax, -B],
  ];
  const sideEnd: P2[] = [
    [A, bz],
    [-ax, B],
    [-A, -bz],
    [ax, -B],
  ];
  for (let c = 0; c < 4; c++) {
    const segs = c % 2 === 0 ? mz : mx;
    for (let s = 1; s < segs; s++) {
      const t = s / segs;
      out.push(lerp(sideStart[c][0], sideEnd[c][0], t), y, lerp(sideStart[c][1], sideEnd[c][1], t));
    }
    const [ox, oz, a0] = corners[c];
    for (let j = 0; j <= k; j++) {
      const a = a0 + (j / k) * (PI / 2);
      out.push(ox + r * Math.cos(a), y, oz + r * Math.sin(a));
    }
  }
  return out;
}

/** A rounded rectangle traced in 2D (x, y), for loft sections: from +x toward +y. */
function rr2(A: number, B: number, rc: number, k: number, mx = 1, my = 1): P2[] {
  const ring = rrRing(A, B, rc, 0, k, mx, my);
  const out: P2[] = [];
  for (let i = 0; i < ring.length; i += 3) out.push([ring[i], ring[i + 2]]);
  return out;
}

/**
 * A body lofted through the offsets of one rounded rectangle: each profile point [d, y] is the outline grown by d
 * (half-sizes A + d, B + d, corner radius rc + d). Hotel pans, food beds: their walls stay parallel at every corner.
 */
function rrRevolve(prof: readonly P2[], A: number, B: number, rc: number, k: number, mx: number, mz: number, spunOf?: THREE.Material): THREE.BufferGeometry {
  const rings = prof.map(([d, y]) => rrRing(A + d, B + d, rc + d, y, k, mx, mz));
  if (!spunOf) return ringMesh(rings);
  // A spun or satin finish's lines run round the pan (v constant along a ring), so u may be anything continuous.
  const scale = (spunOf.userData.uvScale as number | undefined) ?? 1;
  const arcLen = [0];
  for (let j = 1; j < prof.length; j++) arcLen.push(arcLen[j - 1] + Math.hypot(prof[j][0] - prof[j - 1][0], prof[j][1] - prof[j - 1][1]));
  const n = rings[0].length / 3;
  return ringMesh(rings, { seam: true, uv: (j, i) => [(i / n) * 8, arcLen[j] * scale] });
}

/**
 * Sweep a closed 2D section along a path: section x along `side` (kept square to the path), section y along
 * tangent × side. `sec(j)` gives the section at path point j (the same point count each time). Capped and outward.
 */
function sweep(path: readonly V3[], sec: (j: number) => readonly P2[], side: V3 = [1, 0, 0]): THREE.BufferGeometry {
  const rings: number[][] = [];
  const T = new THREE.Vector3();
  const S = new THREE.Vector3();
  const U = new THREE.Vector3();
  const sv = new THREE.Vector3(...side);
  for (let j = 0; j < path.length; j++) {
    const p0 = path[Math.max(0, j - 1)];
    const p1 = path[Math.min(path.length - 1, j + 1)];
    T.set(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]).normalize();
    S.copy(sv).addScaledVector(T, -sv.dot(T)).normalize();
    U.crossVectors(T, S);
    const ring: number[] = [];
    for (const [x, y] of sec(j)) {
      const p = path[j];
      ring.push(p[0] + S.x * x + U.x * y, p[1] + S.y * x + U.y * y, p[2] + S.z * x + U.z * y);
    }
    rings.push(ring);
  }
  return outward(ringMesh(rings, { caps: true }));
}

/** Sample a polyline (with quadratic-rounded corners of radius `bend`) evenly into about `n` points. */
function bentPath(pts: readonly V3[], bend: number, n: number): V3[] {
  const v = pts.map((p) => new THREE.Vector3(p[0], p[1], p[2]));
  const path = new THREE.CurvePath<THREE.Vector3>();
  let from = v[0].clone();
  for (let i = 1; i < v.length; i++) {
    const p = v[i];
    if (i === v.length - 1) {
      path.add(new THREE.LineCurve3(from, p.clone()));
      break;
    }
    const next = v[i + 1];
    const inDir = p.clone().sub(v[i - 1]).normalize();
    const outDir = next.clone().sub(p).normalize();
    const cut = Math.min(bend, p.distanceTo(v[i - 1]) * 0.49, next.distanceTo(p) * 0.49);
    const a = p.clone().addScaledVector(inDir, -cut);
    const b = p.clone().addScaledVector(outDir, cut);
    path.add(new THREE.LineCurve3(from, a));
    path.add(new THREE.QuadraticBezierCurve3(a, p.clone(), b));
    from = b;
  }
  return path.getSpacedPoints(n).map((p) => [p.x, p.y, p.z] as V3);
}

/** Point on a polyline profile at height y (the first span that brackets it), with the profile's unit tangent. */
function profileAt(prof: readonly P2[], y: number): { r: number; tr: number; ty: number } {
  for (let k = 0; k < prof.length - 1; k++) {
    const [r0, y0] = prof[k];
    const [r1, y1] = prof[k + 1];
    if ((y - y0) * (y - y1) <= 0 && y0 !== y1) {
      const t = (y - y0) / (y1 - y0);
      const l = Math.hypot(r1 - r0, y1 - y0);
      return { r: lerp(r0, r1, t), tr: (r1 - r0) / l, ty: (y1 - y0) / l };
    }
  }
  const last = prof[prof.length - 1];
  return { r: last[0], tr: 0, ty: 1 };
}

/** A section thinned for distance: every other point at detail 1, every third at detail 0 (ends kept). */
function lodProf(s: Shop, prof: readonly P2[]): P2[] {
  if (s.detail === 2 || prof.length <= (s.detail === 1 ? 8 : 6)) return prof.slice();
  return decimate(prof, s.detail === 1 ? 2 : 3);
}

/** A round rivet head of radius r and height h, domed, its base at `at` and its axis along `dir`. */
function rivet(s: Shop, mat: THREE.Material, at: V3, dir: V3, r = 0.003, h = 0.0012): void {
  if (s.detail === 0) return;
  const g = lathe(s, [
    [0, -0.0004],
    [r, -0.0004],
    [r, h * 0.25],
    [r * 0.75, h * 0.8],
    [0, h],
  ], s.detail === 2 ? 10 : 6);
  g.applyQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), new THREE.Vector3(...dir).normalize()));
  s.add(mat, g, at);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Skillets                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

export type SkilletKind = "cast_iron" | "tri_ply" | "stainless_304" | "aluminum" | "carbon_steel" | "nonstick";

/** What a skillet's builder reports: the cooking floor's height and radius, the inner radius at the rim, the rim's height. */
export interface SkilletInfo {
  floorY: number;
  floorR: number;
  innerR: number;
  rimY: number;
}

interface BowlSpec {
  /** Outer radius at the rim, outer height of the body. */
  R: number;
  H: number;
  /** Wall and base thickness. */
  t: number;
  tb: number;
  /** Wall angle from horizontal leaving the base bend, and at the rim (smaller flares). */
  a0: number;
  a1: number;
  /** Base-to-wall bend radius. */
  rb: number;
  rim: "round" | "flat";
  /** Raise the whole vessel (over an encapsulated disc). */
  lift?: number;
}

interface BowlProfile {
  outer: P2[];
  cap: P2[];
  inner: P2[];
  floorY: number;
  floorR: number;
  innerR: number;
  rimY: number;
}

/** The section of a pressed or cast vessel: outer surface from the base centre to the rim, the rim, the inside back. */
function bowlProfile(s: Shop, b: BowlSpec): BowlProfile {
  const { R, H, t, tb, a0, a1, rb } = b;
  const y0 = b.lift ?? 0;
  const Ht = b.rim === "round" ? H - t / 2 : H;
  const Sy = rb * (1 - Math.cos(a0));
  const run = (Ht - Sy) * (1 / Math.tan(a0) + 1 / Math.tan(a1)) * 0.5;
  const Sx = R - run;
  const bx = Sx - rb * Math.sin(a0);
  const outer: P2[] = [];
  outer.push([0, 0], [bx * 0.5, 0], [bx, 0]);
  arc(outer, bx, rb, rb, -PI / 2, -PI / 2 + a0, s.seg(7, 3));
  // The wall: a quadratic from the bend (leaving at a0) to the rim (arriving at a1).
  const S: P2 = outer[outer.length - 1];
  const T: P2 = [R, Ht];
  const d0: P2 = [Math.cos(a0), Math.sin(a0)];
  const d1: P2 = [Math.cos(a1), Math.sin(a1)];
  const det = d0[0] * d1[1] - d1[0] * d0[1];
  let C: P2 = [(S[0] + T[0]) / 2, (S[1] + T[1]) / 2];
  if (Math.abs(det) > 1e-4) {
    const dx = T[0] - S[0];
    const dy = T[1] - S[1];
    const u = (dx * d1[1] - d1[0] * dy) / det;
    if (u > 0) C = [S[0] + d0[0] * u, S[1] + d0[1] * u];
  }
  const nw = s.seg(9, 3);
  for (let k = 1; k <= nw; k++) {
    const q = k / nw;
    outer.push([(1 - q) * (1 - q) * S[0] + 2 * q * (1 - q) * C[0] + q * q * T[0], (1 - q) * (1 - q) * S[1] + 2 * q * (1 - q) * C[1] + q * q * T[1]]);
  }
  // The inside: the outside offset inward, base thickness tb blending to wall thickness t round the bend.
  const inner: P2[] = outer.map((p, k) => {
    const pa = outer[Math.max(0, k - 1)];
    const pb = outer[Math.min(outer.length - 1, k + 1)];
    const l = Math.hypot(pb[0] - pa[0], pb[1] - pa[1]);
    const tx = (pb[0] - pa[0]) / l;
    const ty = (pb[1] - pa[1]) / l;
    const ang = Math.atan2(ty, tx);
    const th = lerp(tb, t, smooth(0, a0 * 0.9, ang));
    return [p[0] - ty * th, p[1] + tx * th];
  });
  inner[0][0] = 0;
  const Ot = outer[outer.length - 1];
  const It = inner[inner.length - 1];
  const cap: P2[] = [];
  if (b.rim === "round") {
    const cx = (Ot[0] + It[0]) / 2;
    const cy = (Ot[1] + It[1]) / 2;
    const rr = Math.hypot(Ot[0] - It[0], Ot[1] - It[1]) / 2;
    const aO = Math.atan2(Ot[1] - cy, Ot[0] - cx);
    let aI = Math.atan2(It[1] - cy, It[0] - cx);
    if (aI < aO) aI += TAU;
    arc(cap, cx, cy, rr, aO, aI, s.seg(6, 3), true);
  } else cap.push([Ot[0], Ot[1]], [It[0], It[1]]);
  const lift = (p: P2[]) => p.map(([r, y]) => [r, y + y0] as P2);
  return {
    outer: lift(outer),
    cap: lift(cap),
    inner: lift(inner.slice().reverse()),
    floorY: y0 + tb,
    floorR: bx,
    innerR: It[0],
    rimY: y0 + H,
  };
}

/** Pour spouts on the rim at angles `at`: the upper wall pulled out and dipped into a lip. */
function spoutWarp(at: number[], A: number, H: number, y0 = 0): (th: number, r: number, y: number) => P2 {
  return (th, r, y) => {
    let lobe = 0;
    for (const a of at) lobe += Math.pow(Math.max(0, Math.cos(th - a)), 70);
    const k = Math.pow(smooth(y0 + H * 0.35, y0 + H, y), 2);
    return [r + A * lobe * k, y - A * 0.35 * lobe * k];
  };
}

/** A hanging eye at a handle's end: a flattened ring about the vertical through `at`, tipped by `tilt` about x. */
function handleEye(s: Shop, mat: THREE.Material, at: V3, tilt: number, ro: number, ri: number, th: number): void {
  const e = Math.min(th * 0.45, (ro - ri) * 0.45);
  const prof: P2[] = [];
  const k = s.seg(3, 1);
  arc(prof, ri + e, -th / 2 + e, e, PI, 1.5 * PI, k, true);
  arc(prof, ro - e, -th / 2 + e, e, -PI / 2, 0, k, true);
  arc(prof, ro - e, th / 2 - e, e, 0, PI / 2, k, true);
  arc(prof, ri + e, th / 2 - e, e, PI / 2, PI, k, true);
  prof.push(prof[0]);
  s.add(mat, revolve(prof, s.fine(2) ? 30 : 10), at, [tilt, 0, 0]);
}

/**
 * Build a skillet centred on its axis, its base at the origin, its handle along +z. `diameterM` is the rim's
 * diameter. Options: `handle`, the handle's reach beyond the rim; `height`, the overall height the handle may rise
 * to; `helper`, a helper handle opposite (default on from 12 in, and always on cast iron).
 */
export function skillet(s: Shop, kind: SkilletKind, diameterM: number, opts: { handle?: number; height?: number; helper?: boolean } = {}): SkilletInfo {
  const p = palette();
  const f = finishes();
  const D = diameterM;
  const R = D / 2;
  const iron = kind === "cast_iron";
  const total = opts.height ?? 0.032 + 0.085 * D;
  const H = Math.min(0.026 + 0.075 * D, total - 0.006);
  const L = opts.handle ?? (iron ? 0.07 + 0.3 * D : 0.1 + 0.27 * D);
  const helper = opts.helper ?? (iron || D >= 0.3);
  const n = s.seg(D > 0.32 ? 80 : 96, 12);

  // Each pan's section: cast iron thick and steep, carbon steel pressed thin with flared walls, the clad and
  // aluminium pans between. The disc-base stainless pan stands on its encapsulated disc.
  const disc = kind === "stainless_304" ? 0.005 : 0;
  const spec: BowlSpec =
    kind === "cast_iron"
      ? { R, H, t: 0.0045, tb: 0.0052, a0: 1.32, a1: 1.18, rb: 0.016, rim: "round" }
      : kind === "carbon_steel"
        ? { R, H, t: 0.0024, tb: 0.0026, a0: 1.02, a1: 0.72, rb: 0.028, rim: "round" }
        : kind === "tri_ply"
          ? { R, H, t: 0.0026, tb: 0.0026, a0: 1.25, a1: 0.86, rb: 0.022, rim: "flat" }
          : kind === "stainless_304"
            ? { R, H: H - disc, t: 0.0012, tb: 0.0012, a0: 1.25, a1: 0.9, rb: 0.02, rim: "round", lift: disc }
            : { R: R - 0.002, H, t: 0.0038, tb: 0.0045, a0: 1.2, a1: 0.86, rb: 0.024, rim: "round" };
  const b = bowlProfile(s, spec);
  const warp = iron ? spoutWarp([0, PI], 0.0065 * (D / 0.2), H) : undefined;
  const outMat = iron ? f.iron : kind === "carbon_steel" ? f.carbon : kind === "nonstick" ? f.anod : kind === "aluminum" ? f.spunAlu : kind === "tri_ply" ? f.spunBright : f.spun;
  const inMat = iron ? f.iron : kind === "carbon_steel" ? f.carbon : kind === "nonstick" ? f.ptfe : kind === "aluminum" ? f.spunAlu : f.spun;
  if (outMat === inMat) turned(s, outMat, [...b.outer, ...b.cap.slice(1, -1), ...b.inner], n, warp);
  else {
    turned(s, outMat, b.outer, n, warp);
    turned(s, inMat, b.inner, n, warp);
    if (kind === "tri_ply") {
      // The cut rim shows the cladding: stainless, the aluminium core a shade lower where it polishes back, stainless.
      const [o, i] = [b.cap[0], b.cap[1]];
      const at = (q: number): P2 => [lerp(o[0], i[0], q), lerp(o[1], i[1], q)];
      turned(s, outMat, [at(0), at(0.2)], n);
      turned(s, f.spunAlu, [[at(0.2)[0], at(0.2)[1]], [at(0.2)[0], at(0.2)[1] - 0.0002], [at(0.8)[0], at(0.8)[1] - 0.0002], [at(0.8)[0], at(0.8)[1]]], n);
      turned(s, inMat, [at(0.8), at(1)], n);
    } else if (kind !== "aluminum" && kind !== "nonstick") turned(s, outMat, b.cap, n, warp);
    else turned(s, inMat, [b.cap[0], b.cap[b.cap.length - 1]], n);
  }
  if (kind === "aluminum" || kind === "nonstick") {
    // Rolled rim: the edge curled outward into a bead, stiffening the thin aluminium.
    const top = b.cap[Math.floor(b.cap.length / 2)];
    s.add(outMat, torus(s, b.outer[b.outer.length - 1][0] + 0.0006, 0.0021, "y", TAU, s.seg(96, 10), s.seg(6, 3)), [0, top[1] - 0.0024, 0]);
  }
  if (disc) {
    // The encapsulated base: an aluminium disc sealed under a stainless cap, a shallow step under the pan.
    const rd = b.floorR + 0.006;
    turned(s, f.spunBright, [
      [0, 0],
      [rd - 0.002, 0],
      [rd, 0.0015],
      [rd, disc - 0.0012],
      [rd - 0.003, disc + 0.0004],
      [0, disc + 0.0004],
    ], n);
  }

  // Where the handle meets the wall.
  const yRoot = b.rimY - (iron ? 0.006 : 0.011);
  const wallAt = profileAt(b.outer, yRoot);
  const zRoot = wallAt.r;
  const yEnd = total - (iron ? 0.0068 : 0.0075);
  const zEnd = R + L;
  const rise = Math.atan2(yEnd - yRoot, zEnd - zRoot);

  if (iron) {
    // One casting: the handle beam flares into the wall like a fillet, tapers, and ends in a hanging eye.
    const ro = 0.068 * D;
    const ri = 0.028 * D;
    const eyeTh = 0.05 * D;
    const zEye = zEnd - ro;
    const nu = s.detail === 0 ? 5 : s.seg(26, 8);
    const path: V3[] = [];
    const secs: P2[][] = [];
    for (let j = 0; j <= nu; j++) {
      const u = j / nu;
      const z = lerp(zRoot - 0.003, zEye - ro * 0.35, u);
      const y = lerp(yRoot, yEnd, (z - zRoot) / (zEnd - zRoot));
      path.push([0, y, z]);
      const fl = Math.exp(-u * 9);
      const A = lerp(0.07 * D, 0.05 * D, u) + 0.05 * D * fl;
      const Bh = lerp(0.034 * D, eyeTh / 2, smooth(0.5, 1, u)) + 0.02 * D * fl;
      secs.push(s.detail === 0 ? rr2(A, Bh, Math.min(A, Bh) * 0.7, 1, 1, 1) : rr2(A, Bh, Math.min(A, Bh) * 0.7, s.seg(5, 2)));
    }
    s.add(f.iron, sweep(path, (j) => secs[j]));
    if (s.fine(1)) handleEye(s, f.iron, [0, lerp(yRoot, yEnd, (zEye - zRoot) / (zEnd - zRoot)), zEye], -rise, ro, ri, eyeTh);
  } else {
    // A riveted handle. Its bracket is a plate shaped to the wall, the rivets through it: peened tails flush on the
    // plate, round heads inside the pan. Carbon steel's is a flat forged strap rising to an eye; the others a long
    // cast-stainless (or plated steel) stick with a channel along its top, which aluminium pans sleeve in silicone.
    const strap = kind === "carbon_steel";
    const sleeve = kind === "aluminum" || kind === "nonstick";
    const hMat = strap ? f.carbon : sleeve ? p.chrome : kind === "tri_ply" ? p.ss : p.ssBright;
    const ro = strap ? 0.016 : 0.015;
    const zEye = zEnd - ro;
    const yB = yRoot - 0.004;
    const wB = profileAt(b.outer, yB);
    const nB = s.seg(12, 4);
    const half = strap ? 0.021 : 0.02;
    const bh = strap ? 0.0145 : 0.0125;
    const bpath: V3[] = [];
    for (let k = 0; k <= nB; k++) {
      const a = PI / 2 + ((k / nB - 0.5) * 2 * half) / wB.r;
      const rr = wB.r + 0.0012 * wB.ty;
      bpath.push([Math.cos(a) * rr, yB - 0.0012 * wB.tr, Math.sin(a) * rr]);
    }
    if (s.fine(1)) s.add(hMat, sweep(bpath, (j) => {
      const u = j / nB;
      const e = Math.sqrt(Math.max(0.08, 1 - Math.pow(2 * u - 1, 4)));
      return rr2(bh * e + 0.001, 0.0012, 0.001, 3, 2, 1);
    }, [0, wB.ty, wB.tr]));
    const pins: P2[] = kind === "tri_ply" || kind === "stainless_304" ? [[-0.0085, yB], [0.0085, yB]] : [[-0.0105, yB - 0.004], [0.0105, yB - 0.004], [0, yB + 0.0075]];
    for (const [x, y] of s.fine(1) ? pins : []) {
      const wp = profileAt(b.outer, y);
      s.add(hMat, cyl(s, strap ? 0.0036 : 0.0032, 0.0008, "z", 14), [x, y - wp.tr * 0.0026, wp.r + wp.ty * 0.0026], [Math.atan2(wp.tr, wp.ty), 0, 0]);
      const pin = profileAt(b.inner.slice().reverse(), y);
      rivet(s, strap ? f.carbon : kind === "nonstick" ? p.alu : p.ssBright, [x, y, pin.r - 0.0002], [0, 0.3, -1], strap ? 0.0038 : 0.0036, strap ? 0.0017 : 0.0011);
    }
    const pts = bentPath([[0, yB + (strap ? 0.008 : 0.002), wB.r + 0.0012], [0, yB + 0.01, zRoot + 0.026], [0, yEnd, zEye - 0.006]], 0.026, s.detail === 0 ? 5 : s.detail === 1 ? 9 : s.seg(26, 8));
    const groove = kind === "tri_ply" ? 0.0022 : 0;
    const k0 = Math.floor(pts.length * 0.3);
    // Under a sleeve only the exposed root of the stick is drawn, and a short tang into the eye.
    const stick = sleeve ? pts.slice(0, k0 + 2) : pts;
    s.add(hMat, sweep(stick, (j) => {
      const u = j / (pts.length - 1);
      const fl = Math.exp(-u * 16);
      if (strap) return rr2(lerp(0.0125, 0.0095, smooth(0, 0.4, u)) + 0.004 * fl, 0.0021 + 0.002 * fl, 0.0016, s.seg(3, 1), 2, 1);
      const A = lerp(0.0105, 0.0085, u) + 0.0035 * fl;
      const Bh = 0.0058 + 0.0025 * fl;
      const sec = s.fine(2) ? rr2(A, Bh, Math.min(A, Bh) * 0.8, 5, 4, 2) : rr2(A, Bh, Math.min(A, Bh) * 0.8, 1, s.detail, 1);
      return groove ? sec.map(([x, y]) => [x, y > 0 ? y - groove * Math.max(0, 1 - (x / (A * 0.78)) ** 2) * smooth(0.02, 0.2, u) * (y / Bh) : y] as P2) : sec;
    }));
    if (s.fine(1)) handleEye(s, hMat, [0, yEnd, zEye], -rise, ro, strap ? 0.0066 : 0.0068, strap ? 0.0042 : 0.0085);
    if (sleeve) {
      // The silicone sleeve, slid over the handle from the end: a soft oval, swelling a little under the palm.
      const sp = pts.slice(k0, pts.length - 1);
      s.add(f.silicone, sweep(sp, (j) => {
        const u = j / (sp.length - 1);
        const sw = Math.sin(u * PI) * (1 - 0.3 * u);
        return s.fine(2) ? rr2(0.0122 + 0.0014 * sw, 0.0078 + 0.0009 * sw, 0.0066, 6, 2, 1) : rr2(0.0122 + 0.0014 * sw, 0.0078 + 0.0009 * sw, 0.0066, s.detail === 1 ? 2 : 1, 1, 1);
      }));
    }
  }

  if (helper && s.fine(1)) {
    const yH = b.rimY - (iron ? 0.0055 : 0.01);
    const zH = profileAt(b.outer, yH).r;
    if (iron) {
      // A cast helper lug opposite the handle, its underside scooped for the fingers.
      const tab = new THREE.Shape();
      const hw = 0.11 * D;
      const reach = 0.12 * D;
      tab.moveTo(-hw, -0.004);
      tab.bezierCurveTo(-hw, reach * 0.8, -hw * 0.5, reach, 0, reach);
      tab.bezierCurveTo(hw * 0.5, reach, hw, reach * 0.8, hw, -0.004);
      tab.lineTo(-hw, -0.004);
      const g = plate(s, tab, 0.0062, 0.0022, 16);
      g.rotateX(-PI / 2);
      s.add(f.iron, g, [0, yH - 0.0028, -zH + 0.002], [-0.12, 0, 0]);
    } else {
      // A riveted helper loop: a bent bar on a small plate.
      const hm = kind === "carbon_steel" ? f.carbon : kind === "nonstick" || kind === "aluminum" ? p.chrome : p.ssBright;
      const w2 = 0.032;
      const loop: V3[] = [[-w2, yH - 0.002, -zH + 0.003], [-w2 * 0.85, yH + 0.002, -zH - 0.03], [w2 * 0.85, yH + 0.002, -zH - 0.03], [w2, yH - 0.002, -zH + 0.003]];
      if (s.fine(2)) s.add(hm, bar(s, loop, 0.0042, 0.012, 12));
      else s.add(hm, sweep(bentPath(loop, 0.012, 10), () => circleSec(0.0042, 6), [0, 1, 0]));
      for (const x of [-w2, w2]) {
        const pin = profileAt(b.inner.slice().reverse(), yH - 0.004);
        rivet(s, hm === f.carbon ? f.carbon : p.ssBright, [x * 0.9, yH - 0.004, -pin.r + 0.0002], [0, 0.3, 1], 0.003, 0.001);
      }
    }
  }
  return { floorY: b.floorY, floorR: b.floorR, innerR: b.innerR, rimY: b.rimY };
}

/** A skillet archetype: the bowl toward the back of its box, the handle running forward to the front. */
function skilletFixture(kind: SkilletKind): Builder {
  return (s, e) => {
    // The rim's diameter leaves room in the box for the spouts (cast iron) or the rolled bead, and the helper.
    const iron = kind === "cast_iron";
    const D = iron ? e.w / (1 + 0.0145 / 0.2) : e.w - (kind === "nonstick" || kind === "aluminum" ? 0.0016 : 0.001);
    const helper = kind !== "tri_ply";
    const hr = helper ? (iron ? 0.12 * D + 0.001 : 0.0265) : 0;
    s.at([e.w / 2, 0, -e.d + hr + D / 2], [0, 0, 0], () => {
      skillet(s, kind, D, { handle: e.d - hr - D - 0.001, height: e.h, helper });
    });
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Stockpots                                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

/** A loop handle riveted to a vertical wall of radius R at angle th (from +x toward +z), its plate's centre at y. */
function potHandle(s: Shop, R: number, th: number, y: number, span: number, reach: number, r: number, inner: number): void {
  const p = palette();
  const [c, sn] = [Math.cos(th), Math.sin(th)];
  const out = (rad: number, x: number, yy: number): V3 => [c * rad - sn * x, yy, sn * rad + c * x];
  // The plate: a strip of the wall's curvature, 2.5 mm thick, rounded at its ends.
  const half = span / 2 + 0.012;
  const nP = s.seg(16, 6);
  const path: V3[] = [];
  for (let k = 0; k <= nP; k++) {
    const a = th + ((k / nP - 0.5) * 2 * half) / R;
    path.push([Math.cos(a) * (R + 0.0011), y, Math.sin(a) * (R + 0.0011)]);
  }
  if (s.fine(1))
    s.add(p.ss, sweep(path, (j) => {
      const u = j / nP;
      const e = Math.sqrt(Math.max(0.05, 1 - Math.pow(2 * u - 1, 6)));
      return s.fine(2) ? rr2(0.017 * e, 0.0014, 0.0012, 3, 2, 1) : rr2(0.017 * e, 0.0014, 0.0012, 1, 1, 1);
    }, [0, 1, 0]));
  // The loop: a 304 bar bent square, from the plate out and across.
  const yb = y + 0.004;
  const loop: V3[] = [out(R + 0.001, -span / 2, yb - 0.006), out(R + reach, -span / 2 * 0.92, yb + 0.004), out(R + reach, span / 2 * 0.92, yb + 0.004), out(R + 0.001, span / 2, yb - 0.006)];
  if (s.fine(2)) s.add(p.ssBright, bar(s, loop, r, reach * 0.55, 14));
  else s.add(p.ssBright, sweep(bentPath(loop, reach * 0.55, s.detail === 1 ? 14 : 6), () => circleSec(r, s.detail === 1 ? 6 : 4), [0, 1, 0]));
  // Two rivets each end: domed heads on the plate, flat heads inside the pot (one each end from farther off).
  for (const x of [-span / 2 - 0.004, span / 2 + 0.004])
    for (const dy of s.fine(2) ? [-0.007, 0.007] : [0]) {
      const a = th + x / R;
      const [ca, sa] = [Math.cos(a), Math.sin(a)];
      rivet(s, p.ssBright, [ca * (R + 0.0024), y + dy, sa * (R + 0.0024)], [ca, 0, sa], 0.0027, 0.0013);
      if (s.fine(2)) rivet(s, p.ss, [ca * (inner - 0.0001), y + dy, sa * (inner - 0.0001)], [-ca, 0, -sa], 0.0029, 0.0007);
    }
}

/** A heavy stainless stockpot with an encapsulated base, riveted loop handles front and back, and a looped lid. */
const stockpot: Builder = (s, e) => {
  const p = palette();
  const f = finishes();
  const R = e.w / 2 - 0.0005;
  const disc = 0.0085;
  const lidH = 0.012;
  const knobH = 0.024;
  const H = e.h - lidH - knobH;
  const t = 0.0014;
  const n = s.detail === 1 ? 48 : s.seg(112, 12);
  const cx = e.w / 2;
  const cz = -e.d / 2;
  // The body: drawn 304, straight walls, a tight bottom bend, its rim rolled outward into a bead.
  const rb = 0.012;
  const outer: P2[] = [];
  outer.push([0, disc], [(R - rb) * 0.5, disc], [R - rb, disc]);
  arc(outer, R - rb, disc + rb, rb, -PI / 2, 0, s.seg(7, 3));
  for (let k = 1; k <= 2; k++) outer.push([R, lerp(disc + rb, H - 0.004, k / 2)]);
  const inner: P2[] = [];
  for (let k = 2; k >= 1; k--) inner.push([R - t, lerp(disc + rb, H - 0.004, k / 2)]);
  arc(inner, R - rb, disc + rb, rb - t, 0, -PI / 2, s.seg(7, 3), true);
  inner.push([(R - rb) * 0.5, disc + t], [0, disc + t]);
  const cap: P2[] = [];
  arc(cap, R - t / 2, H - 0.004, t / 2, 0, PI, 4);
  s.at([cx, 0, cz], [0, 0, 0], () => {
    turned(s, f.spun, [...outer, ...cap, ...inner], n);
    if (s.fine(1)) s.add(p.ss, torus(s, R + 0.0012, 0.0024, "y", TAU, s.seg(112, 10), s.fine(2) ? 8 : 3), [0, H - 0.0035, 0]);
    // The impact-bonded base: an aluminium slab sealed in a stainless cap, a step under the pot with a turned chamfer.
    const rd = R - 0.008;
    if (s.fine(1)) turned(s, f.spunBright, [
      [0, 0],
      [rd - 0.003, 0],
      [rd, 0.002],
      [rd, disc - 0.0025],
      [rd - 0.0035, disc + 0.0006],
      [0, disc + 0.0006],
    ], n);
    // Riveted loop handles at the front and back (the grips the scene's sockets name).
    const yh = H - 0.03;
    const span = Math.min(0.11, R * 0.75);
    const reach = Math.max(0.026, (e.d - e.w) / 2 - 0.006);
    potHandle(s, R, PI / 2, yh, span, reach, 0.0042 + R * 0.006, R - t);
    potHandle(s, R, -PI / 2, yh, span, reach, 0.0042 + R * 0.006, R - t);
  });

  // The lid: a shallow dome with a rolled edge that sits on the rim, a riveted loop handle on top.
  const lid = s.joint("pot_lid_lift", [cx, H, cz]);
  const Rl = R + 0.003;
  const lprof: P2[] = [[0, lidH - 0.0015]];
  for (let k = 1; k <= 6; k++) {
    const q = k / 6;
    lprof.push([Rl * 0.93 * q, (lidH - 0.0015) * (1 - q * q) + 0.0015]);
  }
  arc(lprof, Rl * 0.93 + 0.0016, 0.0016, 0.0024, PI * 0.75, -PI * 0.1, 6);
  lprof.push([Rl * 0.93, 0.0006]);
  for (let k = 5; k >= 0; k--) {
    const q = k / 6;
    lprof.push([Rl * 0.93 * q, (lidH - 0.0027) * (1 - q * q) + 0.0006]);
  }
  turned(lid, f.spunBright, lprof, n);
  // A short skirt inside the rim keeps the lid centred (hidden in the pot unless the lid is up close).
  if (s.fine(2)) lid.add(p.ss, lathe(lid, [
    [R - t - 0.0035, 0.002],
    [R - t - 0.0022, 0.0015],
    [R - t - 0.0022, -0.006],
    [R - t - 0.0032, -0.006],
    [R - t - 0.0045, 0.0015],
  ], n));
  const lw = Math.min(0.045, R * 0.32);
  const lloop: V3[] = [[-lw, lidH - 0.002, 0], [-lw * 0.8, lidH + knobH - 0.006, 0], [lw * 0.8, lidH + knobH - 0.006, 0], [lw, lidH - 0.002, 0]];
  if (s.fine(2)) lid.add(p.ssBright, bar(lid, lloop, 0.0045, 0.012, 14));
  else lid.add(p.ssBright, sweep(bentPath(lloop, 0.012, s.detail === 1 ? 12 : 6), () => circleSec(0.0045, s.detail === 1 ? 6 : 4), [0, 0, 1]));
  for (const x of s.fine(1) ? [-lw, lw] : []) {
    lid.add(p.ss, rbox(lid, 0.016, 0.003, 0.014, 0.0012), [x, lidH - 0.0012, 0]);
    rivet(lid, p.ssBright, [x + Math.sign(x) * 0.004, lidH + 0.0003, 0], [0, 1, 0], 0.0024, 0.001);
  }

  // The liquid: a disc inside the wall, its surface rolled by the boil in `tick`; bubbles and steam on top.
  const floor = disc + t;
  const Rin = R - t - 0.0005;
  const K = s.seg(14, 4);
  const N = s.fine(2) ? 64 : s.seg(48, 16);
  const rings: number[][] = [];
  for (let j = 0; j <= K; j++) {
    const ring: number[] = [];
    const rr = (Rin * j) / K;
    for (let i = 0; i < N; i++) ring.push(rr * Math.cos((i / N) * TAU), 0, rr * Math.sin((i / N) * TAU));
    rings.push(ring);
  }
  const surf = ringMesh(rings, { uv: (j, i) => [0.5 + (0.5 * j * Math.cos((i / N) * TAU)) / K, 0.5 + (0.5 * j * Math.sin((i / N) * TAU)) / K] });
  // The inner-to-outer ring order makes it face down; turn it to face up.
  const ix = surf.index!;
  for (let k = 0; k < ix.count; k += 3) {
    const a = ix.getX(k + 1);
    ix.setX(k + 1, ix.getX(k + 2));
    ix.setX(k + 2, a);
  }
  surf.computeVertexNormals();
  // Per pot: its own ripple bump (a shared image), scrolled by the boil.
  const stockM = f.stock.clone();
  stockM.bumpMap = rippleTexture().clone();
  stockM.bumpMap.needsUpdate = true;
  stockM.bumpScale = 0;
  const liquid = new THREE.Mesh(surf, stockM);
  liquid.renderOrder = 2;
  const liqHolder = new THREE.Group();
  liqHolder.add(liquid);
  s.extra(liqHolder, [cx, 0, cz]);
  s.tris += ix.count / 3;
  const base = Float32Array.from(surf.attributes.position.array as Float32Array);

  const maxB = s.fine(2) ? 40 : s.fine(1) ? 16 : 0;
  const bubG = new THREE.SphereGeometry(1, s.fine(2) ? 8 : 6, s.fine(2) ? 5 : 4);
  const bubbles = new THREE.InstancedMesh(bubG, f.bubble, Math.max(1, maxB));
  bubbles.count = 0;
  bubbles.renderOrder = 3;
  bubbles.frustumCulled = false;
  liqHolder.add(bubbles);
  s.tris += (bubG.index!.count / 3) * maxB;
  const br = rng(2024);
  const bub = Array.from({ length: maxB }, () => ({ r: Math.sqrt(br()) * Rin * 0.92, a: br() * TAU, ph: br(), sp: 0.6 + br() * 0.8, sz: 0.4 + br() * 0.6 }));

  const steamMat = new THREE.SpriteMaterial({ map: steamTexture(), color: 0xffffff, transparent: true, opacity: 0, depthWrite: false });
  const puffs = Array.from({ length: s.fine(1) ? 10 : 0 }, (_, i) => {
    const sp = new THREE.Sprite(steamMat);
    sp.renderOrder = 4;
    liqHolder.add(sp);
    return { sp, ph: i / 10, a: br() * TAU, r: br() * Rin * 0.6 };
  });

  let level = 0.7;
  let boil = 0;
  let lidOpen = 0;
  const m4 = new THREE.Matrix4();
  const q4 = new THREE.Quaternion();
  const v4 = new THREE.Vector3();
  const sc4 = new THREE.Vector3();
  const surfaceY = () => floor + level * (H - floor - 0.03);
  const place = (t: number) => {
    const y = surfaceY();
    liquid.visible = level > 0.005;
    liquid.position.y = y;
    // Rolling surface: travelling ripples, and a slow convection dome welling up off-centre at a hard boil.
    const pos = surf.attributes.position as THREE.BufferAttribute;
    const amp = boil * boil * 0.0035;
    const cxw = Math.cos(t * 0.7) * Rin * 0.25;
    const czw = Math.sin(t * 0.9) * Rin * 0.25;
    for (let i = 0; i < pos.count; i++) {
      const x = base[i * 3];
      const z = base[i * 3 + 2];
      const rr = Math.hypot(x, z) / Rin;
      const edge = 1 - rr * rr * rr * rr;
      let dy = 0;
      if (amp > 0) {
        dy += amp * 0.35 * Math.sin(x * 160 + t * 7.1) * Math.sin(z * 140 - t * 6.3);
        dy += amp * 0.25 * Math.sin((x + z) * 230 - t * 11);
        dy += amp * 1.2 * Math.exp(-((x - cxw) ** 2 + (z - czw) ** 2) / (Rin * Rin * 0.12)) * (0.6 + 0.4 * Math.sin(t * 3.3));
      }
      // The meniscus climbs the wall a little.
      pos.setY(i, dy * edge + 0.0012 * Math.pow(rr, 12));
    }
    pos.needsUpdate = true;
    surf.computeVertexNormals();
    // Bubbles: rising from the floor, then domes at the surface that grow and burst.
    const active = Math.round(maxB * clamp01(boil * 1.1));
    bubbles.count = active;
    for (let i = 0; i < active; i++) {
      const b = bub[i];
      const ph = (b.ph + t * b.sp * (0.35 + boil)) % 1;
      const rmax = (0.0012 + 0.0035 * boil) * b.sz;
      let yy: number;
      let rad: number;
      let flat = 1;
      if (ph < 0.72) {
        yy = floor - y + (y - floor) * (ph / 0.72);
        rad = rmax * (0.35 + 0.65 * (ph / 0.72));
      } else {
        const q = (ph - 0.72) / 0.28;
        yy = 0.0004;
        rad = rmax * (1 + 0.6 * q) * (q > 0.85 ? (1 - q) / 0.15 : 1);
        flat = 0.55;
      }
      const a = b.a + ph * 0.6;
      v4.set(Math.cos(a) * b.r, y + yy, Math.sin(a) * b.r);
      sc4.set(rad, rad * flat, rad);
      m4.compose(v4, q4, sc4);
      bubbles.setMatrixAt(i, m4);
    }
    bubbles.instanceMatrix.needsUpdate = true;
    // Fine ripples and bursting cells, scrolling faster the harder it boils.
    stockM.bumpScale = 1.2 * boil;
    stockM.bumpMap!.offset.set(Math.sin(t * 0.8) * 0.05 + t * 0.07 * boil, Math.cos(t * 0.6) * 0.05 - t * 0.05 * boil);
    // Steam: faint wisps above a hard boil, fainter with the lid on.
    const st = smooth(0.55, 1, boil) * (level > 0.01 ? 1 : 0) * (0.25 + 0.75 * lidOpen);
    steamMat.opacity = 0.16 * st;
    for (const pf of puffs) {
      const ph = (pf.ph + t * 0.12) % 1;
      pf.sp.visible = st > 0.01;
      pf.sp.position.set(Math.cos(pf.a + ph) * pf.r, y + 0.02 + ph * R * 1.6, Math.sin(pf.a + ph) * pf.r);
      const sz = R * (0.5 + ph * 1.4);
      pf.sp.scale.set(sz, sz, 1);
      (pf.sp.material as THREE.SpriteMaterial).rotation = ph * 2;
    }
  };
  place(0);
  let clock = 0;
  return (built: Built): Hooks => {
    const j = built.joints.get("pot_lid_lift");
    const rest = j ? j.position.clone() : new THREE.Vector3();
    return {
      set: (name, v) => {
        if (name === "level") level = clamp01(v);
        else if (name === "boil") boil = clamp01(v);
        place(clock);
      },
      open: (fr) => {
        lidOpen = clamp01(fr);
        if (!j) return;
        turnJoint(j, "x", -0.32 * lidOpen);
        j.position.set(rest.x, rest.y + 0.07 * lidOpen, rest.z + 0.03 * lidOpen);
      },
      tick: (t) => {
        clock = t;
        if (boil > 0 || level > 0) place(t);
      },
    };
  };
};

let rippleTex: THREE.CanvasTexture | null = null;
/** Ripples and burst bubble cells on boiling liquid: a tileable bump, 120 mm per repeat. */
function rippleTexture(): THREE.CanvasTexture {
  if (rippleTex) return rippleTex;
  const n = noise2(1717);
  const c = paint(256, 256, (x, y) => {
    const a = fbm(n, x / 18, y / 18, 3);
    const cell = Math.abs(Math.sin(a * 22));
    const v = 128 + (a - 0.5) * 120 + (cell - 0.5) * 50;
    return [v, v, v];
  });
  rippleTex = canvasTex(c, false, true);
  rippleTex.repeat.set(3, 3);
  return rippleTex;
}

let steamTex: THREE.CanvasTexture | null = null;
/** A soft round puff for steam sprites. */
function steamTexture(): THREE.CanvasTexture {
  if (steamTex) return steamTex;
  const [c, x] = canvas(128, 128);
  const g = x.createRadialGradient(64, 64, 4, 64, 64, 62);
  g.addColorStop(0, "rgba(255,255,255,0.9)");
  g.addColorStop(0.45, "rgba(255,255,255,0.35)");
  g.addColorStop(1, "rgba(255,255,255,0)");
  x.fillStyle = g;
  x.fillRect(0, 0, 128, 128);
  const r = rng(91);
  for (let i = 0; i < 40; i++) {
    x.fillStyle = `rgba(255,255,255,${(r() * 0.12).toFixed(3)})`;
    x.beginPath();
    x.arc(30 + r() * 68, 30 + r() * 68, 6 + r() * 16, 0, TAU);
    x.fill();
  }
  steamTex = canvasTex(c, true);
  return steamTex;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Hotel (gastronorm) pans                                                                                        */
/* ------------------------------------------------------------------------------------------------------------ */

/** EN 631-1 gastronorm sizes, outer flange width × depth across (metres). */
const GN: Record<string, P2> = {
  "1/1": [0.53, 0.325],
  "2/3": [0.354, 0.325],
  "1/2": [0.325, 0.265],
  "1/3": [0.325, 0.176],
  "1/4": [0.265, 0.162],
  "1/6": [0.176, 0.162],
  "1/9": [0.176, 0.108],
};
export type GNSize = keyof typeof GN;
/** What a hotel pan holds: braised meat in sauce, rice, vegetables, soup. */
export type PanFood = 0 | 1 | 2 | 3;

interface PanSection {
  prof: P2[];
  A0: number;
  B0: number;
  rc0: number;
  /** Offset of the inner wall at height y. */
  innerD(y: number): number;
  tS: number;
}

/**
 * The section of a pan drawn from 22-gauge (0.8 mm) 304: the bottom bend, the tapered wall (9° for stacking), the flat
 * flange and the bead rolled under its edge. In offsets of the bottom's footprint.
 */
function panSection(s: Shop, W: number, Dp: number, h: number): PanSection {
  const tS = 0.0008;
  const fl = 0.0125;
  const rr = 0.0022;
  const rb = Math.min(0.014, h * 0.2);
  const rf = 0.0022;
  const tau = 0.155;
  const [ct, st] = [Math.cos(tau), Math.sin(tau)];
  const nrm: P2 = [ct, -st];
  // The wall top, under the flange, and its bottom, from the top's half-sizes.
  const At = W / 2 - fl;
  const Bt = Dp / 2 - fl;
  const C2y = h - tS - rf;
  const drop = C2y + rf * st - rb * (1 - st);
  // How far the wall's top stands out from the base's flat: so that the flange's edge lands on the nominal size.
  const inset = drop * Math.tan(tau) + rf * ct - rf - rb * (1 - ct);
  const A0 = At - inset - rf;
  const B0 = Bt - inset - rf;
  const rc0 = Math.max(0.004, 0.026 - inset - rf);
  const prof: P2[] = [];
  for (const q of [1, 0.5]) prof.push([-B0 * q, 0]);
  prof.push([-rb, 0]);
  arc(prof, -rb, rb, rb, -PI / 2, -tau, s.seg(5, 2));
  const Pb = prof[prof.length - 1];
  const len = (C2y + rf * st - Pb[1]) / ct;
  const Pw: P2 = [Pb[0] + st * len, Pb[1] + ct * len];
  prof.push(Pw);
  const C2: P2 = [Pw[0] + nrm[0] * rf, Pw[1] + nrm[1] * rf];
  arc(prof, C2[0], C2[1], rf, PI - tau, PI / 2, s.seg(4, 2));
  const df = C2[0] + fl - rr;
  const Cb: P2 = [df, h - rr];
  const a1 = PI - Math.asin((rr - tS) / rr);
  prof.push([Cb[0] + rr * Math.cos(a1), h - tS]);
  arc(prof, Cb[0], Cb[1], rr, a1, PI * 2.5, s.seg(9, 4));
  prof.push([C2[0], h]);
  arc(prof, C2[0], C2[1], rf + tS, PI / 2, PI - tau, s.seg(4, 2));
  const Iw: P2 = [C2[0] - nrm[0] * (rf + tS), C2[1] - nrm[1] * (rf + tS)];
  const Ib: P2 = [-rb + (rb - tS) * ct, rb - (rb - tS) * st];
  void Iw;
  arc(prof, -rb, rb, rb - tS, -tau, -PI / 2, s.seg(5, 2), true);
  for (const q of [0.5, 1]) prof.push([-B0 * q, tS]);
  return {
    prof,
    A0,
    B0,
    rc0,
    tS,
    innerD: (y: number) => Ib[0] + Math.max(0, y - Ib[1]) * Math.tan(tau),
  };
}

let foodCache: Record<string, THREE.Material> | null = null;
/** The hotel-pan foods' finishes: braising sauce and meat, rice, a vegetable medley, soup. Static, shared. */
function panFoods(): Record<string, THREE.Material> {
  if (foodCache) return foodCache;
  const n = noise2(808);
  // Sauce: a reduced braising liquid, glossy and nearly opaque, fat beads catching the light.
  const sauceC = paint(256, 256, (x, y) => {
    const v = fbm(n, x / 50, y / 50, 3);
    // Beads of rendered fat, amber, on the reduced liquor.
    const fat = smooth(0.8, 0.86, fbm(n, x / 4 + 50, y / 4, 2)) * 0.45;
    return [lerp(48, 66, v) + fat * 50, lerp(21, 28, v) + fat * 30, lerp(9, 12, v) + fat * 8];
  });
  // Braised meat: dark mahogany under its glaze, lighter where fibres split open, near black on the edges.
  // Fibres bundled and wandering (the coordinates warped), split open here and there to a redder brown.
  const warp = (x: number, y: number): P2 => [x + 14 * (fbm(n, x / 60 + 3, y / 60, 2) - 0.5), y + 10 * (fbm(n, x / 50, y / 50 + 8, 2) - 0.5)];
  const meatC = paint(256, 256, (x, y) => {
    const [wx, wy] = warp(x, y);
    const f = fbm(n, wx / 1.6, wy / 16, 3);
    const v = fbm(n, x / 26 + 9, y / 26, 3);
    const split = smooth(0.6, 0.78, f);
    const glaze = smooth(0.42, 0.25, v);
    return [lerp(62, 84, v) * (1 - 0.4 * glaze) + split * 46, lerp(30, 40, v) * (1 - 0.4 * glaze) + split * 18, lerp(15, 19, v) * (1 - 0.35 * glaze) + split * 8];
  });
  const meatB = paint(256, 256, (x, y) => {
    const [wx, wy] = warp(x, y);
    const f = 0.7 * fbm(n, wx / 1.6, wy / 16, 3) + 0.3 * fbm(n, x / 5 + 40, y / 5, 2);
    return [f * 255, f * 255, f * 255];
  });
  // Steamed long-grain rice: grains drawn one by one, slightly translucent at their ends.
  const r = rng(616);
  const [rc, rx] = canvas(512, 512);
  const [rbc, rbx] = canvas(512, 512);
  rx.fillStyle = "#cfc7b4";
  rx.fillRect(0, 0, 512, 512);
  rbx.fillStyle = "#202020";
  rbx.fillRect(0, 0, 512, 512);
  for (let i = 0; i < 5200; i++) {
    const x = r() * 512;
    const y = r() * 512;
    const a = r() * PI;
    const v = 228 + Math.floor(r() * 26);
    for (const [ctx, col] of [
      [rx, `rgb(${v},${v - 2},${v - 10})`],
      [rbx, `rgb(${170 + Math.floor(r() * 80)},0,0)`],
    ] as const) {
      ctx.save();
      ctx.translate(x, y);
      ctx.rotate(a);
      ctx.fillStyle = col;
      ctx.beginPath();
      ctx.ellipse(0, 0, 7.5, 2.1, 0, 0, TAU);
      ctx.fill();
      ctx.restore();
    }
  }
  // Vegetables seen in the gaps: a green and orange medley under its butter gloss.
  const vegC = paint(256, 256, (x, y) => {
    const v = fbm(n, x / 6, y / 6, 3);
    const o = fbm(n, x / 9 + 30, y / 9, 3) > 0.6;
    return o ? [212, 112, 32] : [lerp(40, 96, v), lerp(80, 140, v), lerp(20, 40, v)];
  });
  // Soup: a vegetable broth, golden, fat eyes on its surface.
  const soupC = paint(256, 256, (x, y) => {
    const v = fbm(n, x / 24, y / 24, 3);
    const eye = fbm(n, x / 5 + 70, y / 5, 2) > 0.7 ? 1 : 0;
    return [lerp(168, 196, v) + eye * 30, lerp(104, 128, v) + eye * 26, lerp(34, 46, v) + eye * 10];
  });
  foodCache = {
    sauce: phys({ color: 0xffffff, map: canvasTex(sauceC, true, true), roughness: 0.12, clearcoat: 1, clearcoatRoughness: 0.04, envMapIntensity: 1.3 }, "box", 1 / 0.16),
    meat: phys({ color: 0xffffff, map: canvasTex(meatC, true, true), bumpMap: canvasTex(meatB, false, true), bumpScale: 3.2, roughness: 0.55, clearcoat: 0.3, clearcoatRoughness: 0.3 }, "box", 1 / 0.05),
    carrot: phys({ color: 0xc85e1a, roughness: 0.45, clearcoat: 0.35, clearcoatRoughness: 0.25, sheen: 0.3, sheenColor: new THREE.Color(1, 0.6, 0.3) }, "box"),
    onion: phys({ color: 0xc08e58, roughness: 0.2, clearcoat: 1, clearcoatRoughness: 0.08 }, "box"),
    rice: std({ color: 0xffffff, map: canvasTex(rc, true, true), bumpMap: canvasTex(rbc, false, true), bumpScale: 3, roughness: 0.55 }, "box", 1 / 0.07),
    veg: phys({ color: 0xffffff, map: canvasTex(vegC, true, true), roughness: 0.35, clearcoat: 0.5 }, "box", 1 / 0.1),
    broccoli: std({ color: 0x3f6b26, roughness: 0.75 }, "box"),
    bean: phys({ color: 0x5a8a2e, roughness: 0.32, clearcoat: 0.5 }, "box"),
    soup: phys({ color: 0xffffff, map: canvasTex(soupC, true, true), roughness: 0.06, clearcoat: 1, transparent: true, opacity: 0.94 }, "box", 1 / 0.15),
    herb: std({ color: 0x3e6a22, roughness: 0.6, side: THREE.DoubleSide }, "box"),
  };
  return foodCache;
}

/**
 * An irregular braised chunk, w × h × d: a superellipsoid (a rounded block, as braised short rib sets) pushed about
 * by noise and slumped a little at the base, so that no two look cast from one mould.
 */
function chunk(s: Shop, w: number, h: number, d: number, seed: number): THREE.BufferGeometry {
  // An icosphere has no poles to pinch; merged so that its normals are smooth.
  const ico = new THREE.IcosahedronGeometry(1, s.detail === 2 ? 3 : s.detail === 1 ? 1 : 0);
  ico.deleteAttribute("normal");
  ico.deleteAttribute("uv");
  const g = mergeVertices(ico);
  const n = noise2(seed);
  const p = g.attributes.position as THREE.BufferAttribute;
  const e = 0.68;
  const sp = (v: number) => Math.sign(v) * Math.pow(Math.abs(v), e);
  for (let i = 0; i < p.count; i++) {
    const [x, y, z] = [p.getX(i), p.getY(i), p.getZ(i)];
    const k = 1 + 0.5 * (fbm(n, x * 1.4 + 5, z * 1.4 + y * 1.2, 3) - 0.5) + 0.1 * (n(x * 6 + 30, y * 6 + z * 5) - 0.5);
    const yy = sp(y) * (y < 0 ? 0.75 : 1) - (y > 0.6 ? (y - 0.6) * 0.3 : 0);
    p.setXYZ(i, (sp(x) * k * w) / 2, (yy * k * h) / 2, (sp(z) * k * d) / 2);
  }
  g.computeVertexNormals();
  return g;
}

/**
 * The food in a pan, its surface at y = 0 in the current frame (so a joint can carry it up and down): a surface
 * the shape of the pan's inside at that height, half-sunk pieces on it. `A`, `B`, `rc` are the surface's half-sizes
 * and corner radius; pieces keep inside `inA`, `inB`.
 */
function foodBed(s: Shop, kind: PanFood, A: number, B: number, rc: number, inA: number, inB: number, seed: number): void {
  const m = panFoods();
  const k = s.seg(6, 2);
  const mx = s.seg(24, 4);
  const mz = s.seg(14, 3);
  const surfMat = kind === 0 ? m.sauce : kind === 1 ? m.rice : kind === 2 ? m.veg : m.soup;
  const prof: P2[] = [];
  const steps = kind === 1 || kind === 2 ? s.seg(8, 2) : s.seg(4, 2);
  for (let j = 0; j <= steps; j++) prof.push([-B * (1 - j / steps) + 1e-5, 0]);
  const g = rrRevolve(prof, A, B, rc, k, mx, mz);
  // The inner-to-outer order faces down; turn it up. Rice mounds a little and is never quite level.
  const ix = g.index!;
  for (let t = 0; t < ix.count; t += 3) {
    const a = ix.getX(t + 1);
    ix.setX(t + 1, ix.getX(t + 2));
    ix.setX(t + 2, a);
  }
  if (kind === 1 || kind === 2) {
    const n = noise2(seed);
    const pos = g.attributes.position as THREE.BufferAttribute;
    for (let i = 0; i < pos.count; i++) {
      const x = pos.getX(i);
      const z = pos.getZ(i);
      const e = Math.min(1, (A - Math.abs(x)) / 0.03, (B - Math.abs(z)) / 0.03);
      pos.setY(i, Math.max(0, e) * (0.006 * (1 - (x / A) ** 2) + 0.004 * fbm(n, x * 40, z * 40, 3)));
    }
  }
  g.computeVertexNormals();
  s.add(surfMat, g);
  const r = rng(seed);
  const at = (): [number, number] => [(r() * 2 - 1) * (inA - 0.03), (r() * 2 - 1) * (inB - 0.03)];
  if (kind === 0) {
    // Braised short rib, half-sunk in its sauce, with carrot and pearl onion between.
    const cols = Math.max(2, Math.round((inA * 2) / 0.095));
    const rows = Math.max(1, Math.round((inB * 2) / 0.095));
    for (let i = 0; i < cols && s.fine(1); i++)
      for (let j = 0; j < rows; j++) {
        if (!s.fine(2) && (i + j * cols) % 4 === 3) continue;
        const x = -inA + 0.04 + ((inA * 2 - 0.08) * (i + 0.5)) / cols + (r() - 0.5) * 0.028;
        const z = -inB + 0.04 + ((inB * 2 - 0.08) * (j + 0.5)) / rows + (r() - 0.5) * 0.028;
        const w = 0.045 + r() * 0.035;
        const h = 0.026 + r() * 0.016;
        s.add(m.meat, chunk(s, w, h, w * (0.55 + r() * 0.35), 900 + i * 7 + j), [x, -0.002, z], [(r() - 0.5) * 0.35, r() * PI, (r() - 0.5) * 0.35]);
        // A few thyme leaves caught on top.
        if (s.fine(2))
          for (let l = 0; l < 1; l++) s.add(m.herb, cyl(s, 0.0022, 0.0004, "y", 6), [x + (r() - 0.5) * 0.02, h * 0.42, z + (r() - 0.5) * 0.015], [(r() - 0.5) * 0.5, r() * PI, (r() - 0.5) * 0.5], [1, 1, 0.45]);
      }
    if (s.fine(2))
      for (let i = 0; i < cols * rows * 1.5; i++) {
        const [x, z] = at();
        if (r() < 0.6) s.add(m.carrot, cyl(s, 0.009 + r() * 0.004, 0.008, "y", 12), [x, 0.001, z], [(r() - 0.5) * 0.6, 0, (r() - 0.5) * 0.6]);
        else s.add(m.onion, ball(s, 0.011, 10), [x, 0.0055, z], undefined, [1, 0.85, 1]);
      }
  } else if (kind === 2) {
    // Steamed vegetables: broccoli florets, carrot coins, green beans.
    const count = s.fine(2) ? 46 : s.fine(1) ? 18 : 0;
    for (let i = 0; i < count; i++) {
      const [x, z] = at();
      const q = r();
      if (q < 0.35) {
        s.add(m.broccoli, ball(s, 0.012 + r() * 0.005, 9), [x, 0.008, z], undefined, [1, 0.75, 1]);
        s.add(m.bean, cyl(s, 0.0035, 0.016, "y", 6), [x, 0.0, z]);
      } else if (q < 0.7) s.add(m.carrot, cyl(s, 0.011, 0.005, "y", 12), [x, 0.007, z], [(r() - 0.5) * 0.9, 0, (r() - 0.5) * 0.9]);
      else s.add(m.bean, cyl(s, 0.0035, 0.05 + r() * 0.03, "x", 6), [x, 0.008, z], [0, r() * PI, (r() - 0.5) * 0.3]);
    }
  } else if (kind === 3 && s.fine(1)) {
    // Soup: carrot coins and herb leaves floating.
    for (let i = 0; i < (s.fine(2) ? 40 : 14); i++) {
      const [x, z] = at();
      if (r() < 0.5) s.add(m.carrot, cyl(s, 0.008, 0.003, "y", 10), [x, 0.0003, z], [(r() - 0.5) * 0.3, 0, (r() - 0.5) * 0.3]);
      else s.add(m.herb, cyl(s, 0.004, 0.0004, "y", 6), [x, 0.0008, z], undefined, [1, 1, 0.5]);
    }
  }
}

/**
 * A gastronorm pan centred at the current frame's origin, its base at y = 0: `size` an EN 631-1 fraction, `depthM`
 * its depth, `fill` 0..1 how full and `food` of what. Static; the pan archetypes carry the live version.
 */
export function hotelPan(s: Shop, size: GNSize, depthM: number, fill = 0, food: PanFood = 0): void {
  const [W, Dp] = GN[size] ?? GN["1/1"];
  const sec = panBody(s, W, Dp, depthM);
  if (fill > 0.005) {
    const y = sec.tS + fill * (depthM - sec.tS - 0.008);
    const d = sec.innerD(y) + 0.001;
    s.at([0, y, 0], [0, 0, 0], () => foodBed(s, food, sec.A0 + d, sec.B0 + d, sec.rc0 + d, sec.A0 + sec.innerD(0), sec.B0 + sec.innerD(0), 31 + Math.round(W * 100)));
  }
}

/** The pan itself (body, flange, bead, anti-jam lugs) at the current frame's origin. */
function panBody(s: Shop, W: number, Dp: number, h: number): PanSection {
  const p = palette();
  const f = finishes();
  const sec = panSection(s, W, Dp, h);
  const k = s.seg(6, 2);
  const mx = s.seg(8, 2);
  const mz = s.seg(6, 2);
  s.add(f.satin, rrRevolve(lodProf(s, sec.prof), sec.A0, sec.B0, sec.rc0, k, mx, mz, f.satin));
  if (s.fine(2)) {
    // Anti-jam lugs: dimples pressed out of each corner below the flange, so that stacked pans part.
    const y = h - 0.02;
    const d = sec.innerD(y) + sec.tS;
    const rc = sec.rc0 + d;
    for (const [sx, sz] of [
      [1, 1],
      [-1, 1],
      [-1, -1],
      [1, -1],
    ] as const) {
      const cxl = sx * (sec.A0 + d - rc + rc * Math.SQRT1_2);
      const czl = sz * (sec.B0 + d - rc + rc * Math.SQRT1_2);
      s.add(p.ssDull, ball(s, 0.005, 10), [cxl, y, czl], [0, Math.atan2(-sz, sx), 0], [0.3, 0.6, 1.2]);
    }
  }
  return sec;
}

/** A hotel pan archetype: live fill and food. The shown food is built at once, any other when first asked for. */
function panFixture(): Builder {
  return (s, e) => {
    const cx = e.w / 2;
    const cz = -e.d / 2;
    let sec!: PanSection;
    s.at([cx, 0, cz], [0, 0, 0], () => {
      sec = panBody(s, e.w, e.d, e.h);
    });
    const y0 = sec.tS;
    const yMax = e.h - 0.008;
    const inA = sec.A0 + sec.innerD(0);
    const inB = sec.B0 + sec.innerD(0);
    const dm = sec.innerD(yMax) + 0.001;
    // Each food is a bed sized for a full pan, on its own joint; set("fill") carries it and fits it to the wall.
    const bed = (sh: Shop, kd: PanFood) => foodBed(sh, kd, sec.A0 + dm, sec.B0 + dm, sec.rc0 + dm, inA, inB, 71 + kd * 13);
    bed(s.joint("food_0", [cx, 0, cz]), 0);
    const detail = s.detail;
    return (b) => {
      let fill = 0.6;
      let food: PanFood = 0;
      const beds = new Map<PanFood, THREE.Object3D>();
      const j0 = b.joints.get("food_0");
      if (j0) beds.set(0, j0);
      const get = (kd: PanFood): THREE.Object3D => {
        let o = beds.get(kd);
        if (!o) {
          const sh = new Shop(detail, `food_${kd}`);
          bed(sh, kd);
          o = sh.build().group;
          o.position.set(cx, 0, cz);
          b.group.add(o);
          b.joints.set(`food_${kd}`, o as THREE.Group);
          beds.set(kd, o);
        }
        return o;
      };
      const show = () => {
        const y = y0 + fill * (yMax - y0);
        const d = sec.innerD(y) + 0.001;
        if (fill > 0.005) get(food);
        for (const [kd, o] of beds) {
          o.visible = kd === food && fill > 0.005;
          o.position.y = y;
          // Fit the bed to the tapered wall at this height (a scale is near enough to the true offset).
          o.scale.set((sec.A0 + d) / (sec.A0 + dm), 1, (sec.B0 + d) / (sec.B0 + dm));
        }
      };
      show();
      return {
        set: (name, v) => {
          if (name === "fill") fill = clamp01(v);
          else if (name === "food") food = clamp(Math.round(v), 0, 3) as PanFood;
          show();
        },
      };
    };
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Food shared parts                                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

type RGB = readonly [number, number, number];
/** Piecewise-linear colour ramp: stops at `at`, colours `c`. */
function ramp(at: readonly number[], c: readonly RGB[], k: number): [number, number, number] {
  if (k <= at[0]) return [c[0][0], c[0][1], c[0][2]];
  for (let i = 1; i < at.length; i++)
    if (k <= at[i]) {
      const t = (k - at[i - 1]) / (at[i] - at[i - 1]);
      return [lerp(c[i - 1][0], c[i][0], t), lerp(c[i - 1][1], c[i][1], t), lerp(c[i - 1][2], c[i][2], t)];
    }
  const l = c[c.length - 1];
  return [l[0], l[1], l[2]];
}
const adiff = (a: number, b: number) => {
  let d = (a - b) % TAU;
  if (d > PI) d -= TAU;
  if (d < -PI) d += TAU;
  return d;
};

/** A live canvas texture: an ImageData written pixel by pixel and pushed to the GPU when repainted. */
class LiveCanvas {
  readonly c: HTMLCanvasElement;
  readonly x: CanvasRenderingContext2D;
  readonly img: ImageData;
  readonly tex: THREE.CanvasTexture;
  constructor(readonly w: number, readonly h: number) {
    [this.c, this.x] = canvas(w, h);
    this.img = this.x.createImageData(w, h);
    this.img.data.fill(255);
    this.tex = canvasTex(this.c, true);
  }
  push(): void {
    this.x.putImageData(this.img, 0, 0);
    this.tex.needsUpdate = true;
  }
}

/** Build index-sharing parts of one closed mesh as separate geometries, each with its own UVs, the normals shared. */
function subGeometry(pos: ArrayLike<number>, nrm: ArrayLike<number>, tris: number[], uv: (v: number) => P2): THREE.BufferGeometry {
  const map = new Map<number, number>();
  const P: number[] = [];
  const N: number[] = [];
  const U: number[] = [];
  const I: number[] = [];
  for (const v of tris) {
    let k = map.get(v);
    if (k === undefined) {
      k = map.size;
      map.set(v, k);
      P.push(pos[v * 3], pos[v * 3 + 1], pos[v * 3 + 2]);
      N.push(nrm[v * 3], nrm[v * 3 + 1], nrm[v * 3 + 2]);
      U.push(...uv(v));
    }
    I.push(k);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(P, 3));
  g.setAttribute("normal", new THREE.Float32BufferAttribute(N, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(U, 2));
  g.setIndex(I);
  return g;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* A boneless ribeye                                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * The rib eye in section, as the butcher cuts it across the grain: the eye (longissimus dorsi) filling the middle,
 * the cap (spinalis dorsi) a crescent round its back and one end, a kernel of fat in the seam between them, the
 * small complexus at the cap's end, intramuscular fat flecking both, and the fat edge outside the cap. Laid out in
 * polar form about the steak's middle: the outline R(θ), and inward from it the fat edge, the cap, the seam.
 */
class Ribeye {
  readonly a: number;
  readonly b: number;
  private sx = 1;
  private sz = 1;
  private ox = 0;
  private oz = 0;
  private nT = 1024;
  private Rt = new Float32Array(1024);
  private fatT = new Float32Array(1024);
  private capT = new Float32Array(1024);
  readonly thC = 1.35 * PI;
  readonly thK = 1.35 * PI - 0.95;
  readonly thX = 1.35 * PI - 1.3;
  constructor(a: number, b: number, seed: number) {
    this.a = a;
    this.b = b;
    const n = noise2(seed);
    const R0 = (th: number) => {
      const [c, sn] = [Math.cos(th), Math.sin(th)];
      const re = 1 / Math.sqrt((c / a) ** 2 + (sn / b) ** 2);
      return re * (1 + 0.035 * Math.cos(2 * th + 0.4) + 0.03 * Math.cos(3 * th - 0.8) + 0.018 * Math.cos(5 * th + 1.3) + 0.01 * Math.cos(7 * th - 0.2)) * (1 - 0.06 * Math.max(0, sn) ** 2);
    };
    let [x0, x1, z0, z1] = [1, -1, 1, -1];
    for (let i = 0; i < this.nT; i++) {
      const th = (i / this.nT) * TAU;
      const r = R0(th);
      this.Rt[i] = r;
      x0 = Math.min(x0, r * Math.cos(th));
      x1 = Math.max(x1, r * Math.cos(th));
      z0 = Math.min(z0, r * Math.sin(th));
      z1 = Math.max(z1, r * Math.sin(th));
    }
    // Fit the outline to the box exactly.
    this.sx = (2 * a) / (x1 - x0);
    this.sz = (2 * b) / (z1 - z0);
    this.ox = (x0 + x1) / 2;
    this.oz = (z0 + z1) / 2;
    for (let i = 0; i < this.nT; i++) {
      const th = (i / this.nT) * TAU;
      // The fat edge: thick outside the cap, a thin silverskin round the eye's free side.
      const wc = Math.max(0, Math.cos(clamp(adiff(th, this.thC) / 1.9, -1, 1) * (PI / 2)));
      this.fatT[i] = 0.0011 + 0.0058 * Math.pow(wc, 0.8) * (0.55 + 0.9 * n(th * 3, 7.7));
      // The cap: a crescent, widest past its middle, tapering to nothing at both ends.
      const u = adiff(th, this.thC) / 1.55;
      const cap = Math.abs(u) < 1 ? Math.pow(Math.cos(u * (PI / 2)), 0.75) * (1 + 0.25 * u) : 0;
      this.capT[i] = 0.025 * cap * (0.85 + 0.3 * n(th * 2.3 + 9, 1.1));
    }
  }
  private look(t: Float32Array, th: number): number {
    const f = ((((th / TAU) % 1) + 1) % 1) * this.nT;
    const i = Math.floor(f);
    return lerp(t[i % this.nT], t[(i + 1) % this.nT], f - i);
  }
  /** The outline point at angle th, in the steak's frame (x, z), centred on the box. */
  outline(th: number): P2 {
    const r = this.look(this.Rt, th);
    return [(r * Math.cos(th) - this.ox) * this.sx, (r * Math.sin(th) - this.oz) * this.sz];
  }
  /** Polar position (angle, depth inside the outline) of a point (x, z) in the steak's frame. */
  polar(x: number, z: number): [number, number] {
    const X = x / this.sx + this.ox;
    const Z = z / this.sz + this.oz;
    const th = Math.atan2(Z, X);
    return [th, this.look(this.Rt, th) - Math.hypot(X, Z)];
  }
  fatW(th: number): number {
    return this.look(this.fatT, th);
  }
  capW(th: number): number {
    return this.look(this.capT, th);
  }
}

/** Raw-to-seared colour stops of lean and of fat, against local browning k (0 raw, 1 deep crust, beyond: charring). */
const SEAR_K = [0, 0.2, 0.5, 0.8, 1.05, 1.35] as const;
const LEAN_SEAR: RGB[] = [
  [0, 0, 0],
  [132, 96, 80],
  [122, 72, 40],
  [98, 52, 24],
  [74, 37, 16],
  [44, 23, 10],
];
const FAT_SEAR: RGB[] = [
  [0, 0, 0],
  [232, 212, 172],
  [214, 166, 94],
  [180, 118, 54],
  [140, 82, 34],
  [80, 44, 16],
];
/** The colour of beef at temperature tc on a surface not in contact: raw red, rare, pink, grey. */
const BEEF_T = [38, 50, 57, 63, 70, 85] as const;
const BEEF_C: RGB[] = [
  [142, 26, 34],
  [158, 38, 44],
  [180, 76, 74],
  [172, 110, 100],
  [142, 106, 92],
  [124, 92, 78],
];

const steakBuilder: Builder = (s, e) => {
  const a = e.w / 2;
  const b = e.d / 2;
  // Raw thickness: thickening as it cooks (set "shrink") brings it to the box's height.
  const T = e.h / (1 + 0.35 * 0.3);
  const rib = new Ribeye(a, b, 5150);
  const nRel = noise2(611);
  const nFine = noise2(612);
  const nFib = noise2(613);
  const nEdge = noise2(614);
  const relief = (x: number, z: number) => 0.65 * fbm(nRel, x * 150 + 3, z * 150 + 7, 3) + 0.35 * fbm(nRel, x * 55 + 31, z * 55, 2);
  const spotOf = (x: number, z: number) => smooth(0.38, 0.66, relief(x, z));

  // --- The mesh: two cut faces and the rounded side, one closed body so that their normals meet smoothly. ---
  const N = [28, 64, 160][s.detail];
  const K = [4, 9, 22][s.detail];
  const M = [4, 6, 14][s.detail];
  const er = 0.0036;
  const O: P2[] = [];
  for (let i = 0; i < N; i++) O.push(rib.outline((i / N) * TAU));
  const nrm: P2[] = O.map((_, i) => {
    const p0 = O[(i + N - 1) % N];
    const p1 = O[(i + 1) % N];
    const [tx, tz] = [p1[0] - p0[0], p1[1] - p0[1]];
    const l = Math.hypot(tx, tz);
    return [tz / l, -tx / l];
  });
  // Point the normals outward (the outline runs from +x toward +z, so (tz, -tx) points out).
  const Q: P2[] = O.map((p, i) => [p[0] - nrm[i][0] * er, p[1] - nrm[i][1] * er]);
  let cxs = 0;
  let czs = 0;
  for (const [x, z] of Q) {
    cxs += x / N;
    czs += z / N;
  }
  const pos: number[] = [];
  const kind: number[] = []; // 0 face, 1 side
  const perim: number[] = [];
  const rowV: number[] = [];
  const push = (x: number, y: number, z: number, kd: number, u = 0, v = 0) => {
    pos.push(x, y, z);
    kind.push(kd);
    perim.push(u);
    rowV.push(v);
    return pos.length / 3 - 1;
  };
  // Cumulative perimeter, for the side's u.
  const plen = [0];
  for (let i = 1; i <= N; i++) plen.push(plen[i - 1] + Math.hypot(O[i % N][0] - O[i - 1][0], O[i % N][1] - O[i - 1][1]));
  const L = plen[N];
  const faceY = (x: number, z: number, k: number, top: boolean) => {
    const fade = 1 - Math.pow(k / K, 3);
    if (!top) return 0.00025 * relief(x, z) * fade;
    return T + (0.0011 * (1 - (k / K) ** 2) + 0.0009 * (relief(x, z) - 0.5)) * fade;
  };
  // Face normals from the height field itself, so that the fan of the polar grid never shows in the shading.
  const faceN = new Map<number, [number, number, number]>();
  const ringOf = (k: number, top: boolean) => {
    const out: number[] = [];
    for (let i = 0; i < N; i++) {
      const x = cxs + (Q[i][0] - cxs) * (k / K);
      const z = czs + (Q[i][1] - czs) * (k / K);
      const v = push(x, faceY(x, z, k, top), z, 0);
      const ep = 0.0004;
      const rr = Math.hypot(x - cxs, z - czs) || 1e-6;
      const kf = (dx: number, dz: number) => (k * Math.hypot(x + dx - cxs, z + dz - czs)) / rr;
      const gx = (faceY(x + ep, z, kf(ep, 0), top) - faceY(x - ep, z, kf(-ep, 0), top)) / (2 * ep);
      const gz = (faceY(x, z + ep, kf(0, ep), top) - faceY(x, z - ep, kf(0, -ep), top)) / (2 * ep);
      const l = Math.hypot(gx, 1, gz);
      faceN.set(v, top ? [-gx / l, 1 / l, -gz / l] : [gx / l, -1 / l, gz / l]);
      out.push(v);
    }
    return out;
  };
  const topC = push(cxs, faceY(cxs, czs, 0, true), czs, 0);
  const topR: number[][] = [];
  for (let k = 1; k < K; k++) topR.push(ringOf(k, true));
  // Side rows from the top junction (s = 1) to the bottom junction (s = -1), denser toward the rounded edges.
  const q = 3;
  const side: number[][] = [];
  for (let j = 0; j <= M; j++) {
    const sj = Math.cos((PI * j) / M);
    const off = er * (Math.pow(Math.max(0, 1 - Math.pow(Math.abs(sj), q)), 1 / q) - 1);
    const y = T / 2 + (sj * T) / 2;
    const row: number[] = [];
    for (let i = 0; i < N; i++) row.push(push(O[i][0] + nrm[i][0] * off, y, O[i][1] + nrm[i][1] * off, 1, plen[i] / L, (sj + 1) / 2));
    side.push(row);
  }
  const botR: number[][] = [];
  for (let k = K - 1; k >= 1; k--) botR.push(ringOf(k, false));
  const botC = push(cxs, 0, czs, 0);
  const faceTris: number[] = [];
  const sideTris: number[] = [];
  const quad = (A: number[], B: number[], out: number[]) => {
    for (let i = 0; i < N; i++) {
      const i1 = (i + 1) % N;
      out.push(A[i], B[i], A[i1], A[i1], B[i], B[i1]);
    }
  };
  const fan = (c: number, A: number[], out: number[], flip: boolean) => {
    for (let i = 0; i < N; i++) {
      const i1 = (i + 1) % N;
      if (flip) out.push(c, A[i1], A[i]);
      else out.push(c, A[i], A[i1]);
    }
  };
  const topRings = [...topR, side[0]];
  fan(topC, topRings[0], faceTris, false);
  for (let k = 0; k < topRings.length - 1; k++) quad(topRings[k], topRings[k + 1], faceTris);
  for (let j = 0; j < M; j++) quad(side[j], side[j + 1], sideTris);
  const botRings = [side[M], ...botR];
  for (let k = 0; k < botRings.length - 1; k++) quad(botRings[k], botRings[k + 1], faceTris);
  fan(botC, botRings[botRings.length - 1], faceTris, true);
  const all = new THREE.BufferGeometry();
  all.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  all.setIndex([...faceTris, ...sideTris]);
  outward(all);
  all.computeVertexNormals();
  const ix = Array.from(all.index!.array as ArrayLike<number>);
  const nr = all.attributes.normal.array as Float32Array;
  // Inner rings take the field's normal; the ring next to the side blends it with the mesh's.
  const lastInner = new Set([...topR[topR.length - 1], ...botR[0]]);
  for (const [v, n3] of faceN) {
    const w = lastInner.has(v) ? 0.5 : 1;
    const x = lerp(nr[v * 3], n3[0], w);
    const y = lerp(nr[v * 3 + 1], n3[1], w);
    const z = lerp(nr[v * 3 + 2], n3[2], w);
    const l = Math.hypot(x, y, z);
    nr[v * 3] = x / l;
    nr[v * 3 + 1] = y / l;
    nr[v * 3 + 2] = z / l;
  }
  const faceGeo = subGeometry(pos, nr, ix.slice(0, faceTris.length), (v) => [(pos[v * 3] + a) / (2 * a), (pos[v * 3 + 2] + b) / (2 * b)]);
  // The side's u wraps: give it a seam column by hand.
  const sideIx = ix.slice(faceTris.length);
  const seamFix = new Map<number, number>();
  const extraPos = pos.slice();
  const fixTris: number[] = [];
  for (let t = 0; t < sideIx.length; t += 3) {
    const tri = sideIx.slice(t, t + 3);
    const us = tri.map((v) => perim[v]);
    if (Math.max(...us) - Math.min(...us) > 0.5) {
      // A triangle straddling u = 0: its low-u corners move to copies at u = 1.
      for (let c = 0; c < 3; c++)
        if (us[c] < 0.5) {
          const v = tri[c];
          let w = seamFix.get(v);
          if (w === undefined) {
            w = extraPos.length / 3;
            extraPos.push(pos[v * 3], pos[v * 3 + 1], pos[v * 3 + 2]);
            perim.push(perim[v] + 1);
            rowV.push(rowV[v]);
            seamFix.set(v, w);
          }
          tri[c] = w;
        }
    }
    fixTris.push(...tri);
  }
  const nrX = Array.from(nr);
  for (const [v, w] of seamFix) nrX.splice(w * 3, 3, nr[v * 3], nr[v * 3 + 1], nr[v * 3 + 2]);
  const sideGeo = subGeometry(extraPos, nrX, fixTris, (v) => [perim[v], rowV[v]]);
  void kind;

  // --- The cut face's tissue map: fat fraction, which muscle, fibre, relief, crust grain, depth from the edge. ---
  const FW = s.detail === 2 ? 1024 : s.detail === 1 ? 512 : 192;
  const FH = Math.round((FW * b) / a);
  const pxm = FW / (2 * a);
  const marb = (() => {
    const [, mx] = canvas(FW, FH);
    mx.fillStyle = "#000";
    mx.fillRect(0, 0, FW, FH);
    mx.lineCap = "round";
    mx.lineJoin = "round";
    const r = rng(9091);
    const toPx = (x: number, z: number): P2 => [(x + a) * pxm, (1 - (z + b) / (2 * b)) * FH];
    const region = (x: number, z: number) => {
      const [th, d] = rib.polar(x, z);
      const fw = rib.fatW(th);
      const cw = rib.capW(th);
      if (d < fw + 0.0005) return -1;
      return d < fw + cw ? 1 : 0;
    };
    // Intramuscular fat: flecks and short streaks, heavier in the cap and toward the eye's rim, the eye's heart leaner.
    const count = s.detail === 2 ? 2300 : s.detail === 1 ? 1100 : 320;
    for (let i = 0; i < count; i++) {
      let x = (r() * 2 - 1) * a;
      let z = (r() * 2 - 1) * b;
      const reg = region(x, z);
      if (reg < 0) continue;
      const [, d] = rib.polar(x, z);
      if (reg === 0 && r() > 0.25 + 0.75 * smooth(0.036, 0.008, d)) continue;
      mx.fillStyle = mx.strokeStyle = `rgba(255,255,255,${(0.6 + r() * 0.4).toFixed(2)})`;
      if (r() < 0.55) {
        // A fleck: an irregular blob, two or three overlapping ellipses.
        const rr = (0.00025 + r() * r() * 0.0011) * pxm;
        for (let k = 0; k < 2 + Math.floor(r() * 2); k++) {
          mx.beginPath();
          const [px, py] = toPx(x + (r() - 0.5) * rr / pxm, z + (r() - 0.5) * rr / pxm);
          mx.ellipse(px, py, rr * (0.6 + r()), rr * (0.3 + r() * 0.5), r() * PI, 0, TAU);
          mx.fill();
        }
        continue;
      }
      let ang = r() * TAU;
      const steps = 2 + Math.floor(r() * 5);
      let wid = (reg === 1 ? 0.0002 + r() * 0.00045 : 0.00012 + r() * 0.00035) * pxm;
      mx.lineWidth = Math.max(0.7, wid);
      let [px, py] = toPx(x, z);
      for (let k = 0; k < steps; k++) {
        ang += (r() - 0.5) * 1.1;
        const st = 0.0005 + r() * 0.0012;
        x += Math.cos(ang) * st;
        z += Math.sin(ang) * st;
        if (region(x, z) < 0) break;
        const [qx, qy] = toPx(x, z);
        mx.beginPath();
        mx.moveTo(px, py);
        mx.lineTo(qx, qy);
        mx.stroke();
        [px, py] = [qx, qy];
        wid *= 0.6 + r() * 0.6;
        mx.lineWidth = Math.max(0.6, wid);
      }
    }
    for (let i = 0; i < 4; i++) {
      const th = rib.thK + (r() - 0.5) * 0.5;
      const [ox, oz] = rib.outline(th);
      const d0 = rib.fatW(th) + rib.capW(th) + 0.001;
      let x = ox * (1 - d0 / Math.hypot(ox, oz));
      let z = oz * (1 - d0 / Math.hypot(ox, oz));
      let ang = Math.atan2(-z, -x) + (r() - 0.5) * 1.2;
      mx.strokeStyle = "rgba(255,255,255,0.95)";
      mx.lineWidth = (0.0005 + r() * 0.0006) * pxm;
      mx.beginPath();
      mx.moveTo(...toPx(x, z));
      const len = 8 + Math.floor(r() * 14);
      for (let k = 0; k < len; k++) {
        ang += (r() - 0.5) * 0.8;
        x += Math.cos(ang) * 0.0011;
        z += Math.sin(ang) * 0.0011;
        mx.lineTo(...toPx(x, z));
        if (k > 3 && r() < 0.12) mx.lineWidth *= 0.7;
      }
      mx.stroke();
    }
    return mx.getImageData(0, 0, FW, FH).data;
  })();
  const fF = new Uint8Array(FW * FH); // fat 0..255
  const fM = new Uint8Array(FW * FH); // 0 eye, 1 cap, 2 complexus
  const fFib = new Uint8Array(FW * FH);
  const fSpot = new Uint8Array(FW * FH);
  const fFine = new Uint8Array(FW * FH);
  const fEdge = new Uint8Array(FW * FH);
  const fEdgeFat = new Uint8Array(FW * FH);
  for (let j = 0; j < FH; j++)
    for (let i = 0; i < FW; i++) {
      const x = -a + (i + 0.5) / pxm;
      const z = -b + (1 - (j + 0.5) / FH) * 2 * b;
      const [th, d0] = rib.polar(x, z);
      // Organic boundaries: each line between tissues wanders by a millimetre or so.
      const d = d0 + 0.0011 * (fbm(nEdge, x * 300, z * 300, 2) - 0.5);
      const fw = rib.fatW(th);
      const cw = rib.capW(th);
      const k = j * FW + i;
      let fat = 0;
      let mus = 0;
      let edgeFat = 0;
      if (d < fw) {
        // The ventral side's thin silverskin is a translucent membrane, the cap's fat edge solid.
        fat = fw < 0.0016 ? 0.55 : 1;
        edgeFat = 1;
      } else if (cw > 0.0008 && d < fw + cw) mus = 1;
      else if (cw > 0.0008 && d < fw + cw + 0.0009 + 0.0007 * nEdge(th * 9, 3)) fat = 0.95;
      // The kernel: a lens of fat in the seam; the complexus: a small muscle at the cap's end, rimmed with fat.
      const [kx, kz] = rib.outline(rib.thK);
      const kd = rib.fatW(rib.thK) + rib.capW(rib.thK) + 0.0008;
      const kl = Math.hypot(kx, kz);
      const kcx = kx * (1 - kd / kl);
      const kcz = kz * (1 - kd / kl);
      const ta = Math.atan2(kz, kx) + PI / 2;
      const [ux, uz] = [x - kcx, z - kcz];
      const along = ux * Math.cos(ta) + uz * Math.sin(ta);
      const across = -ux * Math.sin(ta) + uz * Math.cos(ta);
      const wob = 0.75 + 0.5 * fbm(nEdge, x * 160 + 40, z * 160, 2);
      const lens = (along / 0.0115) ** 2 + (across / 0.0052) ** 2;
      if (lens < wob) fat = 1;
      const [xx, xz] = rib.outline(rib.thX);
      const xd = rib.fatW(rib.thX) + 0.0065;
      const xl = Math.hypot(xx, xz);
      const cxx = xx * (1 - xd / xl);
      const cxz = xz * (1 - xd / xl);
      const cpx = (((x - cxx) / 0.0078) ** 2 + ((z - cxz) / 0.0056) ** 2) / (0.7 + 0.6 * fbm(nEdge, x * 140, z * 140 + 70, 2));
      if (cpx < 1) {
        mus = 2;
        fat = Math.max(0, fat - 0.6);
      } else if (cpx < 1.2 && fat < 1) fat = Math.max(fat, 0.85 * smooth(0.35, 0.6, fbm(nEdge, x * 400, z * 400, 2)));
      const m = marb[k * 4] / 255;
      if (!edgeFat) fat = Math.max(fat, m * 0.95);
      fF[k] = Math.round(fat * 255);
      fM[k] = mus;
      fEdgeFat[k] = edgeFat ? 255 : 0;
      fFib[k] = Math.round(fbm(nFib, x * 1700, z * 1700, 2) * 255);
      fSpot[k] = Math.round(spotOf(x, z) * 255);
      fFine[k] = Math.round(fbm(nFine, x * 850, z * 850, 3) * 255);
      fEdge[k] = Math.round(clamp01(d0 / 0.006) * 255);
    }
  const face = new LiveCanvas(FW, FH);
  // Static bump: the high spots, the crust's grain, fat standing a little proud of the lean.
  const bumpC = paint(FW, FH, (i, j) => {
    const k = j * FW + i;
    const v = 70 + fSpot[k] * 0.35 + (fFine[k] - 128) * 0.55 + fF[k] * 0.12;
    return [v, v, v];
  });
  const faceMat = phys({ map: face.tex, bumpMap: canvasTex(bumpC, false), bumpScale: 0.6, roughness: 0.4, clearcoat: 0.4, clearcoatRoughness: 0.32, metalness: 0 }, "own");
  const SW = s.detail === 2 ? 1024 : s.detail === 1 ? 512 : 192;
  const SH = s.detail === 2 ? 128 : s.detail === 1 ? 64 : 24;
  const sideC = new LiveCanvas(SW, SH);
  const sideMat = phys({ map: sideC.tex, roughness: 0.42, clearcoat: 0.4, clearcoatRoughness: 0.35, metalness: 0 }, "own");
  // The side's tissue round the perimeter: fat outside the cap, lean (with its silverskin) elsewhere.
  const sFat = new Float32Array(SW);
  const sNoise = new Float32Array(SW * SH);
  for (let i = 0; i < SW; i++) {
    const u = (i + 0.5) / SW;
    let lo = 0;
    while (lo < N - 1 && plen[lo + 1] / L < u) lo++;
    const th = ((lo + (u * L - plen[lo]) / (plen[lo + 1] - plen[lo])) / N) * TAU;
    sFat[i] = smooth(0.0016, 0.0028, rib.fatW(th));
    // Fibres run through the steak's thickness, so its side is striated vertically; the membrane over it is smooth.
    for (let j = 0; j < SH; j++) sNoise[j * SW + i] = 0.35 * fbm(nFib, u * 1400, (j / SH) * 3, 2) + 0.65 * fbm(nEdge, u * 120, (j / SH) * 4, 3);
  }

  const leanRaw: RGB[] = [
    [152, 24, 32],
    [134, 20, 30],
    [146, 30, 38],
  ];
  let browning = 0;
  let core = Number.NaN;
  let shrink = 0;
  const paintFace = () => {
    const d = face.img.data;
    const bb = browning;
    for (let k = 0; k < FW * FH; k++) {
      const fat = fF[k] / 255;
      const fib = fFib[k] / 255;
      const spot = fSpot[k] / 255;
      const fine = fFine[k] / 255;
      const edge = 0.86 + 0.14 * (fEdge[k] / 255);
      const kk = clamp((bb * (0.8 + 0.42 * spot) + (fine - 0.5) * 0.3 * bb) * edge, 0, 1.25);
      const lr = leanRaw[fM[k]];
      const sh = 0.8 + 0.28 * fib;
      LEAN_SEAR[0] = [lr[0] * sh, lr[1] * sh, lr[2] * sh];
      const ef = fEdgeFat[k] ? 1 : 0;
      FAT_SEAR[0] = [236 - ef * 4, 226 - ef * 8, 210 - ef * 18];
      const lc = ramp(SEAR_K, LEAN_SEAR, kk);
      const fc = ramp(SEAR_K, FAT_SEAR, kk * 0.92);
      // Thin flecks of fat show the lean's pink through them.
      // Seared flecks take colour with the lean round them.
      const fa = (fat < 1 ? Math.pow(fat, 0.8) : 1) * (ef ? 1 : 1 - 0.4 * Math.min(1, kk));
      const grain = 1 + (fine - 0.5) * 0.3 * Math.min(1, kk);
      d[k * 4] = clamp(lerp(lc[0], fc[0], fa) * grain, 0, 255);
      d[k * 4 + 1] = clamp(lerp(lc[1], fc[1], fa) * grain, 0, 255);
      d[k * 4 + 2] = clamp(lerp(lc[2], fc[2], fa) * grain, 0, 255);
    }
    face.push();
  };
  const paintSide = () => {
    const d = sideC.img.data;
    const bb = browning;
    const tc = Number.isFinite(core) ? core : 14 + 46 * Math.pow(bb, 1.2);
    const crustD = 0.0007 + 0.0018 * bb;
    for (let j = 0; j < SH; j++) {
      const v = 1 - (j + 0.5) / SH;
      const dv = Math.min(v, 1 - v) * T;
      // Temperature in from the seared faces: the grey band is where it passed about 66 °C.
      const tau = bb > 0.01 ? tc + (Math.max(tc + 1, 102) - tc) * Math.exp(-dv / 0.0035) * smooth(0, 0.15, bb) : tc;
      const cooked = ramp(BEEF_T, BEEF_C, tau);
      for (let i = 0; i < SW; i++) {
        const k = j * SW + i;
        const nz = sNoise[k];
        const fat = sFat[i];
        let lean = cooked;
        // Silverskin: a pale translucent film over the lean.
        const sk = smooth(0.5, 0.7, nz) * 0.18;
        lean = [lerp(lean[0], 214, sk), lerp(lean[1], 190, sk), lerp(lean[2], 186, sk)];
        const fatTone = smooth(55, 80, tau);
        const fatC: [number, number, number] = [lerp(234, 226, fatTone), lerp(224, 204, fatTone), lerp(206, 160, fatTone)];
        let c: [number, number, number] = [lerp(lean[0], fatC[0], fat), lerp(lean[1], fatC[1], fat), lerp(lean[2], fatC[2], fat)];
        if (dv < crustD && bb > 0.01) {
          const kk = bb * (0.6 + 0.5 * nz) * smooth(crustD, crustD * 0.3, dv) * 1.05;
          LEAN_SEAR[0] = c;
          FAT_SEAR[0] = c;
          const lc = ramp(SEAR_K, LEAN_SEAR, kk);
          const fc = ramp(SEAR_K, FAT_SEAR, kk * 0.9);
          c = [lerp(lc[0], fc[0], fat), lerp(lc[1], fc[1], fat), lerp(lc[2], fc[2], fat)];
        }
        const g = 0.94 + 0.1 * nz;
        d[k * 4] = clamp(c[0] * g, 0, 255);
        d[k * 4 + 1] = clamp(c[1] * g, 0, 255);
        d[k * 4 + 2] = clamp(c[2] * g, 0, 255);
      }
    }
    sideC.push();
  };
  const dress = () => {
    faceMat.bumpScale = 0.5 + 2.4 * browning;
    faceMat.roughness = lerp(0.36, 0.62, browning);
    sideMat.roughness = lerp(0.38, 0.52, browning);
    faceMat.clearcoat = lerp(0.42, 0.28, browning);
  };
  paintFace();
  paintSide();
  dress();
  const j = s.joint("steak", [a, 0, -b]);
  j.add(faceMat, faceGeo);
  j.add(sideMat, sideGeo);
  return (built) => {
    const sj = built.joints.get("steak");
    let lastB = -1;
    let lastC = Number.NaN;
    const apply = () => {
      if (sj) sj.scale.set(1 - 0.45 * shrink, 1 + 0.35 * shrink, 1 - 0.45 * shrink);
    };
    apply();
    return {
      set: (name, v) => {
        if (name === "browning") {
          browning = clamp01(v);
          if (Math.abs(browning - lastB) > 0.004) {
            lastB = browning;
            paintFace();
            paintSide();
            dress();
          }
        } else if (name === "core") {
          core = v;
          if (!(Math.abs(core - lastC) < 0.4)) {
            lastC = core;
            paintSide();
          }
        } else if (name === "shrink") {
          shrink = clamp(v, 0, 0.3);
          apply();
        }
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* A sourdough boule                                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

interface BouleProf {
  L: number;
  r: Float32Array;
  y: Float32Array;
  nr: Float32Array;
  ny: Float32Array;
}
const PROF_N = 1024;
/**
 * The loaf's meridian, crown to the middle of its base, resampled at even arc length: a dome down to its widest
 * point at y0, then a tucked-under edge rolling into a flat base.
 */
function bouleProf(R: number, H: number, y0: number): BouleProf {
  const S = 4096;
  const raw: P2[] = [];
  for (let i = 0; i <= S; i++) {
    const al = (i / S) * PI;
    const c = Math.cos(al);
    const sn = Math.max(0, Math.sin(al));
    raw.push(c >= 0 ? [R * Math.pow(sn, 0.9), y0 + (H - y0) * c] : [R * Math.pow(sn, 0.5), y0 * (1 - Math.pow(-c, 0.25))]);
  }
  const cum = [0];
  for (let i = 1; i <= S; i++) cum.push(cum[i - 1] + Math.hypot(raw[i][0] - raw[i - 1][0], raw[i][1] - raw[i - 1][1]));
  const L = cum[S];
  const r = new Float32Array(PROF_N + 1);
  const y = new Float32Array(PROF_N + 1);
  const nr = new Float32Array(PROF_N + 1);
  const ny = new Float32Array(PROF_N + 1);
  let k = 0;
  for (let i = 0; i <= PROF_N; i++) {
    const t = (i / PROF_N) * L;
    while (k < S - 1 && cum[k + 1] < t) k++;
    const f = clamp01((t - cum[k]) / Math.max(1e-12, cum[k + 1] - cum[k]));
    r[i] = lerp(raw[k][0], raw[k + 1][0], f);
    y[i] = lerp(raw[k][1], raw[k + 1][1], f);
  }
  for (let i = 0; i <= PROF_N; i++) {
    const a = Math.max(0, i - 2);
    const b = Math.min(PROF_N, i + 2);
    const tr = r[b] - r[a];
    const ty = y[b] - y[a];
    const l = Math.hypot(tr, ty) || 1;
    nr[i] = -ty / l;
    ny[i] = tr / l;
  }
  return { L, r, y, nr, ny };
}

/** Crust colour against local browning: raw dough, pale gold, golden, chestnut, deep caramel, nearly burnt. */
const CRUST_K = [0, 0.2, 0.45, 0.7, 0.95, 1.2] as const;
const CRUST_C: RGB[] = [
  [212, 192, 158],
  [218, 180, 120],
  [196, 136, 64],
  [142, 80, 36],
  [90, 50, 24],
  [50, 28, 13],
];

const bouleBuilder: Builder = (s, e) => {
  // Material coordinates: s, arc length from the crown as a fraction of the meridian, and the azimuth θ. They are
  // fixed in the dough, so the crust's marks stretch with it as it springs. The texture is their polar map (a disc
  // of radius 0.42 about the middle) with two strips at the sides for the ear's crust and its underside.
  const D = 0.42;
  const nLump = noise2(7101);
  const nFine = noise2(7102);
  const nBl = noise2(7103);
  const nFl = noise2(7104);
  const nSt = noise2(7105);
  // The score: an arc in the material plane (p, q) = s (cos θ, sin θ), the ear on the side toward its centre,
  // the opening (and the camera's usual view into it) on the far side, toward the front.
  const BETA = 0.4;
  const cdx = -Math.sin(BETA);
  const cdz = -Math.cos(BETA);
  const C0 = 0.4;
  const RS = 0.46;
  const PHI = 0.5;
  const PHIE = 0.44;
  const score = (p: number, q: number): P2 => {
    const wx = p - C0 * cdx;
    const wz = q - C0 * cdz;
    return [Math.hypot(wx, wz) - RS, Math.atan2(-cdx * wz + cdz * wx, -cdx * wx - cdz * wz)];
  };
  const gEnd = (phi: number) => {
    const t = Math.abs(phi) / PHI;
    return t >= 1 ? 0 : Math.sqrt(1 - t * t);
  };
  const ringPhase = (s0: number, th: number) => s0 / 0.042 + th / TAU;
  const ringFade = (s0: number) => smooth(0.02, 0.07, s0) * (1 - smooth(0.56, 0.64, s0));

  let spring = 1.85;
  let crust = 0.8;
  const fOf = () => clamp01((spring - 1) / 0.85);
  let prof = bouleProf(0.0893, 0.0795, 0.024);
  const reprofile = () => {
    const f = fOf();
    const H = lerp(0.05, 0.0738, f);
    prof = bouleProf(lerp(0.0855, 0.0893, f), H, lerp(0.2, 0.3, f) * H);
  };

  /** Displacement along the surface normal at a material point, in metres. */
  const disp = (p: number, q: number, s0: number, th: number, yb: number): number => {
    const f = fOf();
    const L = prof.L;
    let h = 0.0011 * (fbm(nLump, p * 7 + 2, q * 7 + 5, 3) - 0.5) * smooth(0, 0.008, yb);
    const rf = ringFade(s0);
    if (rf > 0) h -= 0.00042 * (1 - 0.45 * f) * rf * (0.5 + 0.5 * Math.cos(TAU * ringPhase(s0, th)));
    const [ds, phi] = score(p, q);
    const g = gEnd(phi);
    if (g > 0 && ds > -0.1 && ds < 0.2) {
      const d = ds * L;
      const W = (0.0016 + 0.022 * f) * g;
      const Dp = (0.0012 + 0.011 * f) * g;
      const Hr = 0.003 * f * g;
      if (d < 0) h += Hr * smooth(-0.014, -0.0015, d);
      else if (d < W) {
        const u = d / W;
        h += Hr * (1 - smooth(0, 0.3, u)) - Dp * Math.sin(PI * Math.pow(u, 0.8)) * (1 - 0.3 * u);
      } else h += 0.0007 * f * g * smooth(W + 0.012, W, d);
    }
    return h;
  };
  /** The surface point at material (θ, s0), displaced. */
  const P3 = (th: number, s0: number, out: number[], o: number) => {
    const fi = clamp(s0, 0, 1) * PROF_N;
    const i = Math.min(PROF_N - 1, Math.floor(fi));
    const t = fi - i;
    const r = lerp(prof.r[i], prof.r[i + 1], t);
    const y = lerp(prof.y[i], prof.y[i + 1], t);
    const nr = lerp(prof.nr[i], prof.nr[i + 1], t);
    const ny = lerp(prof.ny[i], prof.ny[i + 1], t);
    const [c, sn] = [Math.cos(th), Math.sin(th)];
    const h = disp(s0 * c, s0 * sn, s0, th, y);
    out[o] = (r + nr * h) * c;
    out[o + 1] = Math.max(0, y + ny * h);
    out[o + 2] = (r + nr * h) * sn;
  };
  const Pq = (p: number, q: number, out: number[], o = 0) => P3(Math.atan2(q, p), Math.hypot(p, q), out, o);

  // --- The mesh: a polar grid in material coordinates (denser over the crown where the score is), and the ear. ---
  const N = [24, 56, 120][s.detail];
  const M = [14, 30, 84][s.detail];
  const nE = [0, 16, 48][s.detail];
  const KT = s.detail === 2 ? 7 : 4;
  const rowS: number[] = [];
  {
    const w = (x: number) => lerp(1.7, 0.55, smooth(0.3, 0.7, x));
    const cum = [0];
    const Z = 2000;
    for (let i = 1; i <= Z; i++) cum.push(cum[i - 1] + w((i - 0.5) / Z));
    let k = 0;
    for (let j = 0; j <= M; j++) {
      const t = (j / M) * cum[Z];
      while (k < Z - 1 && cum[k + 1] < t) k++;
      rowS.push((k + (t - cum[k]) / (cum[k + 1] - cum[k])) / Z);
    }
  }
  const mTh: number[] = [];
  const mS: number[] = [];
  const uv: number[] = [];
  const vtx = (th: number, s0: number) => {
    mTh.push(th);
    mS.push(s0);
    uv.push(0.5 + D * s0 * Math.cos(th), 0.5 + D * s0 * Math.sin(th));
    return mS.length - 1;
  };
  const top = vtx(0, 0);
  const rows: number[][] = [];
  for (let j = 1; j < M; j++) {
    const row: number[] = [];
    for (let i = 0; i < N; i++) row.push(vtx((i / N) * TAU, rowS[j]));
    rows.push(row);
  }
  const bot = vtx(0, 1);
  const domeV = mS.length;
  const idx: number[] = [];
  for (let i = 0; i < N; i++) idx.push(top, rows[0][(i + 1) % N], rows[0][i]);
  for (let j = 0; j < rows.length - 1; j++)
    for (let i = 0; i < N; i++) {
      const i1 = (i + 1) % N;
      const [a, b, c, d] = [rows[j][i], rows[j][i1], rows[j + 1][i1], rows[j + 1][i]];
      idx.push(a, b, d, b, c, d);
    }
  const last = rows[rows.length - 1];
  for (let i = 0; i < N; i++) idx.push(last[i], last[(i + 1) % N], bot);
  const domeI = idx.length;
  // The ear: a thin flap, sections along the score, each a loop up its crust face to the lip and back beneath.
  const LOOP = 2 * KT;
  const flap0 = mS.length;
  for (let e = 0; e <= nE; e++)
    for (let k = 0; k < LOOP; k++) {
      const topSheet = k < KT;
      const uu = topSheet ? k / (KT - 1) : (LOOP - 1 - k) / (KT - 1);
      mTh.push(0);
      mS.push(0);
      uv.push(topSheet ? 0.93 + 0.06 * uu : 0.01 + 0.06 * uu, 0.02 + (0.96 * e) / nE);
    }
  for (let e = 0; e < nE; e++)
    for (let k = 0; k < LOOP - 1; k++) {
      const a = flap0 + e * LOOP + k;
      const b = a + 1;
      const c = a + LOOP + 1;
      const d = a + LOOP;
      idx.push(a, d, b, b, d, c);
    }
  const nV = mS.length;
  const pos = new Float32Array(nV * 3);
  const geo = new THREE.BufferGeometry();
  const posA = new THREE.BufferAttribute(pos, 3);
  posA.setUsage(THREE.DynamicDrawUsage);
  geo.setAttribute("position", posA);
  geo.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  geo.setIndex(idx);
  let flipDome = false;
  let flipFlap = false;
  const tmp: number[] = [0, 0, 0, 0, 0, 0];
  const shape = () => {
    reprofile();
    const out = pos as unknown as number[];
    for (let v = 0; v < domeV; v++) P3(mTh[v], mS[v], out, v * 3);
    const f = fOf();
    const L = prof.L;
    const V = new THREE.Vector3();
    const A = new THREE.Vector3();
    const T = new THREE.Vector3();
    const Nn = new THREE.Vector3();
    const ep = 0.004;
    for (let e = 0; e <= nE && nE > 0; e++) {
      const phi = lerp(-PHIE, PHIE, e / nE);
      const te = Math.sqrt(Math.max(0, 1 - (phi / PHIE) ** 2));
      const g = gEnd(phi);
      const W = (0.0016 + 0.022 * f) * g;
      const Hr = 0.003 * f * g;
      // The material direction away from the score's centre at this point, and along it.
      const ca = Math.cos(phi);
      const sa = Math.sin(phi);
      const dx = -cdx * ca - -cdz * sa;
      const dz = -cdx * sa + -cdz * ca;
      const rr = RS - 0.0035 / L;
      const px = C0 * cdx + rr * dx;
      const pz = C0 * cdz + rr * dz;
      Pq(px, pz, tmp, 0);
      V.set(tmp[0], tmp[1], tmp[2]);
      Pq(px + dx * ep, pz + dz * ep, tmp, 0);
      Pq(px - dx * ep, pz - dz * ep, tmp, 3);
      A.set(tmp[0] - tmp[3], tmp[1] - tmp[4], tmp[2] - tmp[5]).normalize();
      Pq(px - dz * ep, pz + dx * ep, tmp, 0);
      Pq(px + dz * ep, pz - dx * ep, tmp, 3);
      T.set(tmp[0] - tmp[3], tmp[1] - tmp[4], tmp[2] - tmp[5]).normalize();
      Nn.crossVectors(T, A).normalize();
      if (Nn.y < 0) Nn.negate();
      A.addScaledVector(Nn, -A.dot(Nn)).normalize();
      // The section: a curved lip from its root (buried in the crust) up and over the opening.
      const earH = ((0.0004 + 0.0078 * f) * Math.pow(g, 0.7) - 0.003 * (1 - smooth(0, 0.12, f))) * te;
      const bury = -0.0008 - 0.0025 * (1 - smooth(0, 0.12, f));
      const b0x = -0.002;
      const c1x = 0.004 * te;
      const b2x = (0.0035 + 0.55 * W) * te;
      const b0y = bury;
      const c1y = 0.6 * earH + Hr * 0.1;
      const b2y = earH;
      for (let k = 0; k < LOOP; k++) {
        const topSheet = k < KT;
        const u = topSheet ? k / (KT - 1) : (LOOP - 1 - k) / (KT - 1);
        const x = (1 - u) ** 2 * b0x + 2 * u * (1 - u) * c1x + u * u * b2x;
        const y = (1 - u) ** 2 * b0y + 2 * u * (1 - u) * c1y + u * u * b2y;
        let tx = 2 * (1 - u) * (c1x - b0x) + 2 * u * (b2x - c1x);
        let ty = 2 * (1 - u) * (c1y - b0y) + 2 * u * (b2y - c1y);
        const tl = Math.hypot(tx, ty) || 1;
        tx /= tl;
        ty /= tl;
        const th = lerp(0.0045, 0.0012, u) * (0.4 + 0.6 * te) * smooth(0, 0.1, f + 0.02);
        const sg = topSheet ? 0.5 : -0.5;
        const sx = x - ty * th * sg;
        const sy = y + tx * th * sg;
        const o = (flap0 + e * LOOP + k) * 3;
        pos[o] = V.x + A.x * sx + Nn.x * sy;
        pos[o + 1] = V.y + A.y * sx + Nn.y * sy;
        pos[o + 2] = V.z + A.z * sx + Nn.z * sy;
      }
    }
    posA.needsUpdate = true;
    geo.computeVertexNormals();
    geo.computeBoundingSphere();
  };
  shape();
  // Wind the dome outward (by its volume) and the ear's crust face up.
  {
    let vol = 0;
    const a = new THREE.Vector3();
    const b = new THREE.Vector3();
    const c = new THREE.Vector3();
    const get = (v: THREE.Vector3, i: number) => v.set(pos[i * 3], pos[i * 3 + 1], pos[i * 3 + 2]);
    for (let t = 0; t < domeI; t += 3) vol += get(a, idx[t]).dot(get(b, idx[t + 1]).cross(get(c, idx[t + 2])));
    flipDome = vol < 0;
    if (nE > 0) {
      const t = domeI + (Math.floor(nE / 2) * (LOOP - 1) + Math.floor(KT / 2)) * 6;
      get(a, idx[t]);
      get(b, idx[t + 1]).sub(a);
      get(c, idx[t + 2]).sub(a);
      flipFlap = b.cross(c).y < 0;
    }
    const flipRange = (t0: number, t1: number) => {
      for (let t = t0; t < t1; t += 3) [idx[t + 1], idx[t + 2]] = [idx[t + 2], idx[t + 1]];
    };
    if (flipDome) flipRange(0, domeI);
    if (flipFlap) flipRange(domeI, idx.length);
    geo.setIndex(idx);
    geo.computeVertexNormals();
  }

  // --- The texture: static layers in material space, then the colour of crust, flour and opening, repainted. ---
  const TW = [256, 512, 1024][s.detail];
  const NP = TW * TW;
  const lS = new Float32Array(NP);
  const lD = new Float32Array(NP);
  const lG = new Float32Array(NP);
  const lRing = new Uint8Array(NP);
  const lFine = new Uint8Array(NP);
  const lBl = new Uint8Array(NP);
  const lFl = new Uint8Array(NP);
  const lSt = new Uint8Array(NP);
  const zone = new Uint8Array(NP); // 0 loaf, 1 ear crust, 2 ear underside
  for (let j = 0; j < TW; j++)
    for (let i = 0; i < TW; i++) {
      const k = j * TW + i;
      const u = (i + 0.5) / TW;
      const v = 1 - (j + 0.5) / TW;
      let p = (u - 0.5) / D;
      let q = (v - 0.5) / D;
      if (u > 0.925 || u < 0.075) {
        // The ear's strips: across from its root (0) to its lip (1), along the score.
        const uu = u > 0.925 ? (u - 0.93) / 0.06 : (u - 0.01) / 0.06;
        const vv = clamp01((v - 0.02) / 0.96);
        zone[k] = u > 0.925 ? 1 : 2;
        lS[k] = clamp01(uu);
        lG[k] = lerp(-PHIE, PHIE, vv);
        p = clamp01(uu) * 3;
        q = vv * 9;
        lFine[k] = Math.round(fbm(nFine, p * 9, q * 9, 3) * 255);
        lBl[k] = Math.round(smooth(0.62, 0.8, fbm(nBl, p * 40, q * 40, 2)) * 255);
        lSt[k] = Math.round(fbm(nSt, p * 2, q * 60, 2) * 255);
        continue;
      }
      const s0 = Math.min(1, Math.hypot(p, q));
      const th = Math.atan2(q, p);
      lS[k] = s0;
      const [ds, phi] = score(p, q);
      lD[k] = ds;
      lG[k] = phi;
      const rp = ringPhase(s0, th);
      const fr = rp - Math.floor(rp);
      // The cane's line wanders a little and its width varies.
      const fw = 0.1 + 0.08 * fbm(nFl, p * 40 + 9, q * 40, 2);
      lRing[k] = Math.round(Math.exp(-(((fr - 0.5 + 0.06 * (fbm(nLump, p * 30, q * 30, 2) - 0.5)) / fw) ** 2)) * ringFade(s0) * 255);
      lFine[k] = Math.round(fbm(nFine, p * 60, q * 60, 3) * 255);
      lBl[k] = Math.round(smooth(0.64, 0.8, fbm(nBl, p * 210, q * 210, 2)) * 255);
      lFl[k] = Math.round(fbm(nFl, p * 22, q * 22, 3) * 255);
      // Stretch marks in the opening run along the cut.
      lSt[k] = Math.round(fbm(nSt, ds * 140, phi * 3, 2) * 255);
    }
  const tex = new LiveCanvas(TW, TW);
  const bumpC = paint(TW, TW, (i, j) => {
    const k = j * TW + i;
    const b = 150 - lRing[k] * 0.22 + (lFine[k] - 128) * 0.35 - lBl[k] * 0.16 + (lFl[k] - 128) * 0.08;
    return [b, b, b];
  });
  const mat = phys({ map: tex.tex, bumpMap: canvasTex(bumpC, false), bumpScale: 1.3, roughness: 0.66, clearcoat: 0, clearcoatRoughness: 0.5, metalness: 0 }, "own");
  const repaint = () => {
    const f = fOf();
    const c = crust;
    const L = prof.L;
    const d = tex.img.data;
    const flourT = 0.42 + 0.2 * f;
    const flourC: RGB = [lerp(242, 226, c), lerp(238, 206, c), lerp(228, 166, c)];
    for (let k = 0; k < NP; k++) {
      const fine = lFine[k] / 255;
      const bl = lBl[k] / 255;
      let col: [number, number, number];
      const z = zone[k];
      if (z) {
        const uu = lS[k];
        const g = gEnd(lG[k]);
        const kk = z === 1 ? c * (1.0 + 0.24 * uu * g + 0.16 * (fine - 0.5) + 0.08 * bl) : c * (0.5 + 0.36 * uu * uu + 0.12 * (lSt[k] / 255 - 0.5));
        col = ramp(CRUST_K, CRUST_C, kk);
      } else {
        const s0 = lS[k];
        const ring = lRing[k] / 255;
        const rf = ringFade(s0);
        const fl = lFl[k] / 255;
        let kk = c * (0.9 + 0.2 * (fine - 0.5) + 0.07 * bl + 0.14 * (fl - 0.5) + 0.08 * (1 - ring) * rf + 0.12 * smooth(0.6, 0.7, s0));
        let open = 0;
        const g = gEnd(lG[k]);
        if (g > 0) {
          const dm = lD[k] * L;
          const W = (0.0016 + 0.022 * f) * g;
          if (dm < 0) kk += c * 0.34 * g * smooth(-0.007, 0, dm);
          else if (dm < W) {
            const u = dm / W;
            open = smooth(0, 0.12, u) * smooth(1, 0.8, u);
            const ko = c * (0.56 + 0.3 * smooth(0.5, 1, u) + 0.12 * (lSt[k] / 255 - 0.5));
            kk = lerp(kk, ko, open);
          } else kk += c * 0.2 * g * smooth(W + 0.004, W, dm);
        }
        col = ramp(CRUST_K, CRUST_C, kk);
        // Rice flour from the banneton: in the cane's grooves and in patches, thinning as the dough stretches.
        const fa = smooth(flourT, flourT + 0.25, ring * (0.35 + 0.9 * fl)) * 0.85 + smooth(0.68, 0.85, fl) * 0.16 * smooth(0.04, 0.12, s0) * (1 - smooth(0.5, 0.62, s0));
        const a = clamp01(fa) * (1 - open) * (1 - 0.3 * c);
        col = [lerp(col[0], flourC[0], a), lerp(col[1], flourC[1], a), lerp(col[2], flourC[2], a)];
      }
      d[k * 4] = clamp(col[0], 0, 255);
      d[k * 4 + 1] = clamp(col[1], 0, 255);
      d[k * 4 + 2] = clamp(col[2], 0, 255);
    }
    tex.push();
    mat.roughness = lerp(0.72, 0.56, c);
    mat.clearcoat = 0.16 * c;
  };
  repaint();
  const j = s.joint("boule", [e.w / 2, 0, -e.d / 2]);
  j.add(mat, geo);
  return () => {
    let lastS = spring;
    let lastP = spring;
    let lastC = crust;
    return {
      set: (name, v) => {
        if (name === "spring") {
          spring = clamp(v, 1, 1.85);
          if (Math.abs(spring - lastS) > 0.002) {
            lastS = spring;
            shape();
          }
          if (Math.abs(spring - lastP) > 0.015) {
            lastP = spring;
            repaint();
          }
        } else if (name === "crust") {
          crust = clamp01(v);
          if (Math.abs(crust - lastC) > 0.004) {
            lastC = crust;
            repaint();
          }
        }
      },
    };
  };
};

/* ------------------------------------------------------------------------------------------------------------ */
/* A fried egg                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

const eggBuilder: Builder = (s, e) => {
  const nO = noise2(4401);
  const nT = noise2(4402);
  const nL = noise2(4403);
  const nB = noise2(4404);
  // The white's outline about the yolk, irregular as an egg spreads in a pan; fitted to the box.
  const Rof = (th: number) =>
    0.045 * (1 + 0.08 * Math.cos(2 * th + 0.7) + 0.05 * Math.cos(3 * th - 1.1) + 0.035 * Math.cos(5 * th + 0.4) + 0.012 * Math.cos(9 * th + 2) + 0.07 * (fbm(nO, Math.cos(th) * 2 + 3, Math.sin(th) * 2 + 5, 3) - 0.5));
  const rtOf = (th: number) => 0.025 * (1 + 0.16 * (fbm(nT, Math.cos(th) * 1.6 + 7, Math.sin(th) * 1.6, 2) - 0.5));
  let [x0, x1, z0, z1] = [1, -1, 1, -1];
  for (let i = 0; i < 720; i++) {
    const th = (i / 720) * TAU;
    const r = Rof(th);
    x0 = Math.min(x0, r * Math.cos(th));
    x1 = Math.max(x1, r * Math.cos(th));
    z0 = Math.min(z0, r * Math.sin(th));
    z1 = Math.max(z1, r * Math.sin(th));
  }
  const sx = e.w / (x1 - x0);
  const sz = e.d / (z1 - z0);
  const toX = (x: number) => (x - x0) * sx;
  const toZ = (z: number) => -e.d + (z - z0) * sz;
  const yx = toX(0);
  const yz = toZ(0);
  /** The white's thickness at (r, θ) about the yolk: a thin spread, the thick white standing round the yolk. */
  const thick = (r: number, th: number) => smooth(rtOf(th) + 0.007, rtOf(th) - 0.003, r);
  const hOf = (r: number, th: number) => {
    const R = Rof(th);
    let h = 0.0021 + 0.0042 * thick(r, th) + 0.0009 * smooth(0.022, 0.012, r);
    const er = 0.0028;
    if (r > R - er) h *= Math.sqrt(Math.max(0, 1 - ((r - (R - er)) / er) ** 2));
    // A crisp edge ruffles as it fries.
    const edge = smooth(0.009, 0, R - r);
    h += edge * 0.0007 * (fbm(nL, Math.cos(th) * 20 + r * 300, Math.sin(th) * 20, 2) - 0.4) * smooth(0, 0.0015, R - r);
    return 0.0003 + Math.max(0, h);
  };

  // --- The white: a polar grid from under the yolk to its edge. ---
  const N = [16, 40, 96][s.detail];
  const K = [8, 18, 40][s.detail];
  const rIn = 0.011;
  const pos: number[] = [];
  const uvs: number[] = [];
  for (let k = 0; k <= K; k++) {
    const t = Math.pow(k / K, 0.85);
    for (let i = 0; i <= N; i++) {
      const th = (i / N) * TAU;
      const r = lerp(rIn, Rof(th), t);
      const X = toX(r * Math.cos(th));
      const Z = toZ(r * Math.sin(th));
      pos.push(X, k === K ? 0.0003 : hOf(r, th), Z);
      uvs.push(X / e.w, (Z + e.d) / e.d);
    }
  }
  const idx: number[] = [];
  for (let k = 0; k < K; k++)
    for (let i = 0; i < N; i++) {
      const a = k * (N + 1) + i;
      const b = a + 1;
      const c = a + N + 2;
      const d = a + N + 1;
      idx.push(a, b, d, b, c, d);
    }
  const white = new THREE.BufferGeometry();
  white.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  white.setAttribute("uv", new THREE.Float32BufferAttribute(uvs, 2));
  white.setIndex(idx);
  white.computeVertexNormals();
  // Face it up, and weld the seam's normals.
  {
    const nn = white.attributes.normal as THREE.BufferAttribute;
    if (nn.getY(Math.floor(K / 2) * (N + 1)) < 0) {
      for (let t = 0; t < idx.length; t += 3) [idx[t + 1], idx[t + 2]] = [idx[t + 2], idx[t + 1]];
      white.setIndex(idx);
      white.computeVertexNormals();
    }
    for (let k = 0; k <= K; k++) {
      const a = k * (N + 1);
      const b = a + N;
      const v = new THREE.Vector3(nn.getX(a) + nn.getX(b), nn.getY(a) + nn.getY(b), nn.getZ(a) + nn.getZ(b)).normalize();
      nn.setXYZ(a, v.x, v.y, v.z);
      nn.setXYZ(b, v.x, v.y, v.z);
    }
  }

  // --- Its texture: translucent raw albumen setting to white, the edge browning into lace with holes. ---
  const TW = [128, 384, 768][s.detail];
  const TH = Math.round((TW * e.d) / e.w);
  const NP = TW * TH;
  const lR = new Float32Array(NP);
  const lTk = new Uint8Array(NP);
  const lLace = new Uint8Array(NP);
  const lHole = new Uint8Array(NP);
  const lBub = new Uint8Array(NP);
  for (let j = 0; j < TH; j++)
    for (let i = 0; i < TW; i++) {
      const k = j * TW + i;
      const X = ((i + 0.5) / TW) * e.w;
      const Z = -e.d + (1 - (j + 0.5) / TH) * e.d;
      const x = X / sx + x0;
      const z = (Z + e.d) / sz + z0;
      const r = Math.hypot(x, z);
      const th = Math.atan2(z, x);
      lR[k] = Rof(th) - r;
      lTk[k] = Math.round(thick(r, th) * 255);
      lLace[k] = Math.round(fbm(nL, x * 260, z * 260, 3) * 255);
      lHole[k] = Math.round(fbm(nB, x * 520, z * 520, 2) * 255);
      lBub[k] = Math.round(smooth(0.66, 0.74, fbm(nB, x * 900 + 40, z * 900, 2)) * 255);
    }
  const tex = new LiveCanvas(TW, TH);
  const bumpC = paint(TW, TH, (i, j) => {
    const k = j * TW + i;
    const b = 140 + lBub[k] * 0.25 + (lLace[k] - 128) * 0.3 * smooth(0.01, 0, lR[k]);
    return [b, b, b];
  });
  const whiteM = phys({ map: tex.tex, bumpMap: canvasTex(bumpC, false), bumpScale: 0.6, roughness: 0.22, clearcoat: 0.7, clearcoatRoughness: 0.14, metalness: 0, transparent: true, alphaTest: 0.08 }, "own");
  const LACE_K = [0, 0.35, 0.7, 1] as const;
  const LACE_C: RGB[] = [
    [250, 249, 243],
    [222, 178, 104],
    [172, 108, 46],
    [104, 60, 26],
  ];
  let set = 0.9;
  const repaint = () => {
    const d = tex.img.data;
    const lace = smooth(0.62, 1, set);
    for (let k = 0; k < NP; k++) {
      const dr = lR[k];
      if (dr < -0.0004) {
        d[k * 4 + 3] = 0;
        continue;
      }
      const tk = lTk[k] / 255;
      const sl = clamp01(set * 1.35 - 0.35 * tk - 0.04 * (lLace[k] / 255 - 0.5));
      const op = smooth(0, 0.85, sl);
      let col: [number, number, number] = [lerp(206, 250, op), lerp(210, 249, op), lerp(190, 243, op)];
      let a = lerp(0.26, 1, Math.pow(op, 0.8));
      // Small bubbles set as pale specks.
      const bub = lBub[k] / 255;
      col = [lerp(col[0], 255, bub * 0.3 * op), lerp(col[1], 255, bub * 0.3 * op), lerp(col[2], 255, bub * 0.3 * op)];
      const eb = smooth(0.011, 0.0, dr);
      if (lace > 0 && eb > 0) {
        const lk = clamp01(lace * eb * (0.55 + 0.9 * (lLace[k] / 255 - 0.25)));
        col = ramp(LACE_K, LACE_C, lk);
        a = Math.max(a, smooth(0, 0.3, lk));
        // Lace: the thinnest white fries through into holes.
        if (lHole[k] / 255 > 1 - 0.38 * lace * smooth(0.0055, 0.0005, dr)) a = 0;
      }
      if (dr < 0) a *= smooth(-0.0004, 0, dr);
      d[k * 4] = clamp(col[0], 0, 255);
      d[k * 4 + 1] = clamp(col[1], 0, 255);
      d[k * 4 + 2] = clamp(col[2], 0, 255);
      d[k * 4 + 3] = clamp(a * 255, 0, 255);
    }
    tex.push();
    whiteM.roughness = lerp(0.12, 0.34, set);
    whiteM.clearcoat = lerp(0.9, 0.5, set);
  };
  repaint();
  s.add(whiteM, white);

  // --- The yolk: a dome standing in the thick white, glossy under its membrane. ---
  const yolkC = paint(128, 128, (x, y) => {
    const v = fbm(nB, x / 18, y / 18, 3);
    return [lerp(232, 248, v), lerp(122, 146, v), lerp(6, 20, v)];
  });
  const yolkM = phys({ color: 0xffffff, map: canvasTex(yolkC, true, true), roughness: 0.16, clearcoat: 1, clearcoatRoughness: 0.05, sheen: 0.4, sheenColor: new THREE.Color(1, 0.75, 0.35), metalness: 0 }, "own", 1 / 0.03);
  const ry = 0.0166;
  const yc = 0.0026;
  const yprof: P2[] = [];
  const nY = [4, 9, 18][s.detail];
  const aMax = Math.acos((0.0004 - yc) / (0.86 * ry));
  for (let k = 0; k <= nY; k++) {
    const a = aMax * (1 - k / nY);
    yprof.push([ry * Math.sin(a), yc + 0.86 * ry * Math.cos(a)]);
  }
  s.add(yolkM, revolve(yprof, [10, 24, 48][s.detail], undefined, yolkM, true), [yx, 0, yz]);
  return () => ({
    set: (name, v) => {
      if (name === "set") {
        const nv = clamp01(v);
        if (Math.abs(nv - set) > 0.004) {
          set = nv;
          repaint();
        }
      }
    },
  });
};

/* ------------------------------------------------------------------------------------------------------------ */
/* A basket load of fries                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

const FRY_K = [0, 0.25, 0.5, 0.75, 1, 1.25] as const;
const FRY_C: RGB[] = [
  [232, 222, 186],
  [238, 214, 150],
  [222, 162, 70],
  [194, 122, 40],
  [150, 84, 28],
  [98, 56, 20],
];

const friesBuilder: Builder = (s, e) => {
  const r = rng(3131);
  const count = [24, 40, 110][s.detail];
  const P = [4, 8, 12][s.detail];
  const LS = [1, 2, 3][s.detail];
  const half = 0.00476;
  // The pile: fries dropped one at a time onto a height field, each coming to rest on what is below it.
  const cell = 0.003;
  const nx = Math.ceil(e.w / cell);
  const nz = Math.ceil(e.d / cell);
  const hf = new Float32Array(nx * nz);
  const hAt = (x: number, z: number) => {
    const i = clamp(Math.floor(x / cell), 0, nx - 1);
    const j = clamp(Math.floor((z + e.d) / cell), 0, nz - 1);
    return hf[j * nx + i];
  };
  const gauss = () => {
    let u = 0;
    for (let k = 0; k < 4; k++) u += r();
    return (u - 2) / 0.58;
  };
  const pos: number[] = [];
  const uv: number[] = [];
  const idx: number[] = [];
  const up = new THREE.Vector3();
  const side = new THREE.Vector3();
  const ax = new THREE.Vector3();
  const Y = new THREE.Vector3(0, 1, 0);
  // Section: a rounded square, P points, a corner every quarter (u = 0.125, 0.375, ...).
  const sec: [number, number, number][] = [];
  for (let k = 0; k < P; k++) {
    const c = Math.floor((k * 4) / P);
    const f = (k * 4) / P - c;
    const ang = PI / 4 + (c * PI) / 2 + (P === 4 ? 0 : (f - 0.5 + 0.5 / (P / 4)) * (PI / 2) * 0.9);
    sec.push([Math.cos(ang), Math.sin(ang), (k + 0.5) / P]);
  }
  for (let n = 0; n < count; n++) {
    const L = 0.05 + r() * 0.045;
    const a0 = half * (0.9 + 0.18 * r());
    // The rounded square's corner reaches a little past its half-size.
    const a = a0 * 1.2;
    let best: { cx: number; cz: number; psi: number; ya: number; yb: number } | null = null;
    let bestH = Infinity;
    for (let t = 0; t < 6; t++) {
      const psi = r() * PI;
      const [dx, dz] = [Math.cos(psi) * (L / 2), Math.sin(psi) * (L / 2)];
      // Most fall toward the middle of the basket and heap there; the rest anywhere in it.
      const wide = r() < 0.35;
      const cx = clamp(e.w / 2 + (wide ? (r() - 0.5) * e.w : gauss() * e.w * 0.17), Math.abs(dx) + a + 0.0004, e.w - Math.abs(dx) - a - 0.0004);
      const cz = clamp(-e.d / 2 + (wide ? (r() - 0.5) * e.d : gauss() * e.d * 0.17), -e.d + Math.abs(dz) + a + 0.0004, -Math.abs(dz) - a - 0.0004);
      // Rest on the highest support in each half, then clear anything between.
      const S = 9;
      let ya = 0;
      let yb = 0;
      const hs: number[] = [];
      for (let k = 0; k < S; k++) {
        const q = k / (S - 1) - 0.5;
        let h = 0;
        for (const o of [-a, 0, a]) h = Math.max(h, hAt(cx + 2 * q * dx - (o * dz) / (L / 2), cz + 2 * q * dz + (o * dx) / (L / 2)));
        hs.push(h);
        if (q <= 0) ya = Math.max(ya, h);
        if (q >= 0) yb = Math.max(yb, h);
      }
      for (let k = 0; k < S; k++) {
        const q = k / (S - 1);
        const need = hs[k] - lerp(ya, yb, q);
        if (need > 0) {
          ya += need;
          yb += need;
        }
      }
      const slope = Math.atan2(yb - ya, L);
      if (Math.abs(slope) > 0.5) {
        const m = (ya + yb) / 2;
        const lim = Math.tan(0.5) * (L / 2);
        ya = clamp(ya, m - lim, m + lim);
        yb = clamp(yb, m - lim, m + lim);
      }
      const top = Math.max(ya, yb) + 2 * a;
      if (top > e.h - 0.0005) continue;
      const score = (ya + yb) / 2 + r() * 0.004;
      if (score < bestH) {
        bestH = score;
        best = { cx, cz, psi, ya: ya + a, yb: yb + a };
      }
    }
    if (!best) continue;
    const { cx, cz, psi, ya, yb } = best;
    const A = new THREE.Vector3(cx - (Math.cos(psi) * L) / 2, ya, cz - (Math.sin(psi) * L) / 2);
    const B = new THREE.Vector3(cx + (Math.cos(psi) * L) / 2, yb, cz + (Math.sin(psi) * L) / 2);
    ax.subVectors(B, A).normalize();
    side.crossVectors(Y, ax).normalize();
    up.crossVectors(ax, side).normalize();
    const roll = (r() - 0.5) * 0.3 + (r() < 0.5 ? 0 : PI / 4) * (P === 4 ? 0 : 0);
    const [cr, sr] = [Math.cos(roll), Math.sin(roll)];
    const band = Math.floor(r() * 4);
    // A fry bridging a gap droops; one lying on the bottom cannot.
    const sag = Math.min(0.0012 + r() * 0.0025, Math.max(0, Math.min(ya, yb) - a - 0.0002));
    const wob = (r() - 0.5) * 0.002;
    const taper = 0.94 + r() * 0.08;
    const base = pos.length / 3;
    for (let k = 0; k <= LS; k++) {
      const q = k / LS;
      const ctr = new THREE.Vector3().lerpVectors(A, B, q).addScaledVector(up, -sag * 4 * q * (1 - q)).addScaledVector(side, wob * Math.sin(PI * q));
      const aa = a0 * lerp(1, taper, q);
      for (let j = 0; j <= P; j++) {
        const [c0, s0, u] = sec[j % P];
        // A square of half-size aa with corners rounded: push the circle's points out toward the square.
        const sq = Math.max(Math.abs(c0), Math.abs(s0));
        const kx = (c0 / sq) * 0.78 + c0 * 0.3;
        const ky = (s0 / sq) * 0.78 + s0 * 0.3;
        const lx = (kx * cr - ky * sr) * aa;
        const ly = (kx * sr + ky * cr) * aa;
        pos.push(ctr.x + side.x * lx + up.x * ly, ctr.y + side.y * lx + up.y * ly, ctr.z + side.z * lx + up.z * ly);
        uv.push(j === P ? 1 + u - 1 / P / 2 + 1 / P / 2 : u, (band + 0.04 + 0.92 * q) / 4);
      }
    }
    for (let k = 0; k < LS; k++)
      for (let j = 0; j < P; j++) {
        const a0 = base + k * (P + 1) + j;
        const b0 = a0 + 1;
        const c0 = a0 + P + 2;
        const d0 = a0 + P + 1;
        idx.push(a0, b0, d0, b0, c0, d0);
      }
    // The cut ends, browned most.
    for (const end of [0, LS]) {
      const ring = base + end * (P + 1);
      const c = pos.length / 3;
      let mx = 0;
      let my = 0;
      let mz = 0;
      for (let j = 0; j < P; j++) {
        mx += pos[(ring + j) * 3] / P;
        my += pos[(ring + j) * 3 + 1] / P;
        mz += pos[(ring + j) * 3 + 2] / P;
      }
      const capBase = pos.length / 3 + 1;
      pos.push(mx, my, mz);
      uv.push(0.5, (band + (end ? 0.985 : 0.015)) / 4);
      for (let j = 0; j < P; j++) {
        pos.push(pos[(ring + j) * 3], pos[(ring + j) * 3 + 1], pos[(ring + j) * 3 + 2]);
        uv.push(0.5 + 0.08 * Math.cos((j / P) * TAU), (band + (end ? 0.985 : 0.015)) / 4);
      }
      for (let j = 0; j < P; j++) {
        const j1 = (j + 1) % P;
        if (end) idx.push(c, capBase + j, capBase + j1);
        else idx.push(c, capBase + j1, capBase + j);
      }
    }
    // The fry now lies on the field.
    const tl = 2 * a + 0.0004;
    for (let k = 0; k <= 16; k++) {
      const q = k / 16;
      const px = lerp(A.x, B.x, q);
      const pz = lerp(A.z, B.z, q);
      const py = lerp(ya, yb, q) - a + tl;
      for (const o of [-a, 0, a]) {
        const i = clamp(Math.floor((px - side.x * o) / cell), 0, nx - 1);
        const j = clamp(Math.floor((pz - side.z * o + e.d) / cell), 0, nz - 1);
        hf[j * nx + i] = Math.max(hf[j * nx + i], py);
      }
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("uv", new THREE.Float32BufferAttribute(uv, 2));
  g.setIndex(idx);
  g.computeVertexNormals();
  {
    // Wound outward fry by fry is easiest checked once: the first side face should face away from its axis.
    const p0 = new THREE.Vector3().fromArray(pos, idx[0] * 3);
    const p1 = new THREE.Vector3().fromArray(pos, idx[1] * 3);
    const p2 = new THREE.Vector3().fromArray(pos, idx[2] * 3);
    const nrm = new THREE.Vector3().crossVectors(p1.clone().sub(p0), p2.clone().sub(p0));
    const ctr = new THREE.Vector3();
    for (let j = 0; j < P; j++) ctr.add(new THREE.Vector3().fromArray(pos, j * 3));
    ctr.divideScalar(P);
    if (nrm.dot(p0.clone().sub(ctr)) < 0) {
      for (let t = 0; t < idx.length; t += 3) [idx[t + 1], idx[t + 2]] = [idx[t + 2], idx[t + 1]];
      g.setIndex(idx);
      g.computeVertexNormals();
    }
  }
  // The shared skin: four bands of fries browning at slightly different rates, corners and cut ends first.
  const FU = [32, 64, 128][s.detail];
  const FV = [256, 512, 1024][s.detail];
  const nF = noise2(3132);
  const lN = new Uint8Array(FU * FV);
  const lBl = new Uint8Array(FU * FV);
  for (let j = 0; j < FV; j++)
    for (let i = 0; i < FU; i++) {
      lN[j * FU + i] = Math.round(fbm(nF, (i / FU) * 6, (j / FV) * 90, 3) * 255);
      lBl[j * FU + i] = Math.round(smooth(0.64, 0.78, fbm(nF, (i / FU) * 30 + 50, (j / FV) * 500, 2)) * 255);
    }
  const tex = new LiveCanvas(FU, FV);
  tex.tex.wrapS = THREE.RepeatWrapping;
  const bump = paint(FU, FV, (i, j) => {
    const b = 128 + lBl[j * FU + i] * 0.35 + (lN[j * FU + i] - 128) * 0.25;
    return [b, b, b];
  });
  const mat = phys({ map: tex.tex, bumpMap: canvasTex(bump, false), bumpScale: 0.8, roughness: 0.5, clearcoat: 0.3, clearcoatRoughness: 0.3, metalness: 0 }, "own");
  let browning = 0.75;
  const repaint = () => {
    const d = tex.img.data;
    for (let j = 0; j < FV; j++) {
      const vb = 1 - (j + 0.5) / FV;
      const band = Math.min(3, Math.floor(vb * 4));
      const q = clamp01((vb * 4 - band - 0.04) / 0.92);
      const mult = [0.86, 0.96, 1.05, 1.16][band];
      const endness = Math.exp(-((Math.min(q, 1 - q) / 0.07) ** 2));
      for (let i = 0; i < FU; i++) {
        const u = (i + 0.5) / FU;
        const cu = Math.abs(((u + 0.125) % 0.25) - 0.125);
        const corner = Math.exp(-((cu / 0.03) ** 2));
        const k = j * FU + i;
        const kk = browning * mult * (0.82 + 0.28 * corner + 0.42 * endness + 0.2 * (lN[k] / 255 - 0.5) + 0.1 * (lBl[k] / 255));
        const c = ramp(FRY_K, FRY_C, kk);
        d[k * 4] = c[0];
        d[k * 4 + 1] = c[1];
        d[k * 4 + 2] = c[2];
      }
    }
    tex.push();
    mat.roughness = lerp(0.6, 0.42, browning);
    mat.clearcoat = lerp(0.05, 0.22, browning);
  };
  repaint();
  s.add(mat, g);
  return () => ({
    set: (name, v) => {
      if (name === "browning") {
        const nv = clamp01(v);
        if (Math.abs(nv - browning) > 0.004) {
          browning = nv;
          repaint();
        }
      }
    },
  });
};

/* ------------------------------------------------------------------------------------------------------------ */
/* A braise, as it comes out of a full-size pan                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

const braiseBuilder: Builder = (s, e) => {
  const m = panFoods();
  // The block takes the pan's taper (0.155 rad) and corners; its sauce stands 22 mm below the box's top, the meat
  // proud of it.
  const level = e.h - 0.022;
  const grow = Math.tan(0.155) * level;
  const A = e.w / 2 - grow;
  const B = e.d / 2 - grow;
  const rc = 0.022;
  s.at([e.w / 2, 0, -e.d / 2], [0, 0, 0], () => {
    const prof: P2[] = [];
    const steps = s.seg(4, 2);
    for (let j = 0; j <= steps; j++) prof.push([-B * (1 - j / steps) + 1e-5, 0]);
    const ws = s.seg(6, 2);
    for (let j = 1; j <= ws; j++) prof.push([(grow * j) / ws, (level * j) / ws]);
    const g = rrRevolve(prof, A, B, rc, s.seg(6, 2), s.seg(24, 4), s.seg(14, 3));
    g.computeVertexNormals();
    s.add(m.sauce, outwardOpen(g, [0, level / 2, 0]));
    s.at([0, level, 0], [0, 0, 0], () => foodBed(s, 0, A + grow, B + grow, rc + grow, A + grow - 0.012, B + grow - 0.012, 4242));
  });
  return () => ({});
};

/** Turn an open shell so that its faces point away from `c` (by their majority). */
function outwardOpen(g: THREE.BufferGeometry, c: V3): THREE.BufferGeometry {
  const p = g.attributes.position;
  const ix = g.index!;
  let vote = 0;
  const a = new THREE.Vector3();
  const b = new THREE.Vector3();
  const d = new THREE.Vector3();
  const C = new THREE.Vector3(c[0], c[1], c[2]);
  for (let t = 0; t < ix.count; t += 3) {
    a.fromBufferAttribute(p, ix.getX(t));
    b.fromBufferAttribute(p, ix.getX(t + 1)).sub(a);
    d.fromBufferAttribute(p, ix.getX(t + 2)).sub(a);
    vote += Math.sign(b.cross(d).dot(a.sub(C)));
  }
  if (vote < 0) {
    for (let t = 0; t < ix.count; t += 3) {
      const k = ix.getX(t + 1);
      ix.setX(t + 1, ix.getX(t + 2));
      ix.setX(t + 2, k);
    }
    ix.needsUpdate = true;
  }
  g.computeVertexNormals();
  return g;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Utensils                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

const circleSec = (r: number, n: number): P2[] => Array.from({ length: n }, (_, k) => [r * Math.cos((k / n) * TAU), r * Math.sin((k / n) * TAU)] as P2);
/** A superellipse section, half-sizes a × b, exponent p (2 an ellipse, larger toward a rectangle). */
const superSec = (a: number, b: number, p: number, n: number, cx = 0, cy = 0): P2[] =>
  Array.from({ length: n }, (_, k) => {
    const t = (k / n) * TAU;
    const c = Math.cos(t);
    const sn = Math.sin(t);
    return [cx + a * Math.sign(c) * Math.pow(Math.abs(c), 2 / p), cy + b * Math.sign(sn) * Math.pow(Math.abs(sn), 2 / p)] as P2;
  });
/** Keep every `step`-th point of a profile, and its ends. */
function decimate<T>(pts: readonly T[], step: number): T[] {
  if (step <= 1) return pts.slice();
  const out = pts.filter((_, i) => i % step === 0);
  if ((pts.length - 1) % step !== 0) out.push(pts[pts.length - 1]);
  return out;
}
/** A heightfield patch over a grid of rows (i) and columns (k), faced toward `upY` (+1 up, -1 down). */
function gridMesh(nu: number, nv: number, at: (i: number, k: number) => V3, upY: number): THREE.BufferGeometry {
  const pos: number[] = [];
  for (let i = 0; i <= nu; i++) for (let k = 0; k <= nv; k++) pos.push(...at(i, k));
  const idx: number[] = [];
  for (let i = 0; i < nu; i++)
    for (let k = 0; k < nv; k++) {
      const a = i * (nv + 1) + k;
      const b = a + 1;
      const c = a + nv + 2;
      const d = a + nv + 1;
      idx.push(a, d, b, b, d, c);
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setIndex(idx);
  g.computeVertexNormals();
  const n = g.attributes.normal as THREE.BufferAttribute;
  let sy = 0;
  for (let i = 0; i < n.count; i++) sy += n.getY(i);
  if (sy * upY < 0) {
    for (let t = 0; t < idx.length; t += 3) [idx[t + 1], idx[t + 2]] = [idx[t + 2], idx[t + 1]];
    g.setIndex(idx);
    g.computeVertexNormals();
  }
  return g;
}

/**
 * A forged chef's knife lying on its side: x across the blade (edge at x = 0, spine at x = w), z along it (tip at
 * the back, butt at the front), y its thickness. The blade: a full flat grind from a distal-tapered spine to a
 * polished edge bevel, the edge straight from the heel then sweeping up through the belly to the point. A full
 * bolster, the tang running the handle's length between two acetal scales on three rivets. It rests on its handle
 * and its point, so it lies a degree or two off level.
 */
function chefKnife(inches: 8 | 10): Builder {
  return (s, e) => {
    const p = palette();
    const f = finishes();
    const big = inches === 10;
    const Lh = big ? 0.131 : 0.125;
    const Lbol = 0.014;
    const zB = -Lh;
    const zH = zB - Lbol;
    const Lb = e.d + zH - 0.0005;
    const yc = 0.0095;
    const xSp = e.w - 0.0012;
    const xEd = 0.0012;
    const xT = xEd + 0.55 * (xSp - xEd);
    const spine = (t: number) => (t < 0.55 ? xSp : xSp - (xSp - xT) * ((t - 0.55) / 0.45) ** 1.8);
    const edgeX = (t: number) => (t < 0.5 ? xEd : xEd + (xT - xEd) * ((t - 0.5) / 0.5) ** 2.4);
    const zm = -Lh / 2;
    const Z = (z: number) => z - zm;
    const tsU = (u: number) => (0.0071 + 0.0012 * Math.sin(PI * u)) * (u > 0.93 ? Math.sqrt(Math.max(0, 1 - ((u - 0.93) / 0.07) ** 2)) : 1);
    // At rest it bears on its point and on the low side of its handle or bolster: find the tilt that has both down.
    const support: P2[] = [];
    for (let k = 0; k <= 40; k++) support.push([yc - 0.0012 - tsU(k / 40), Z(lerp(-Lh, -0.0004, k / 40))]);
    for (let k = 0; k <= 10; k++) support.push([yc - lerp(0.0083, 0.0019, Math.pow(smooth(0.1, 1, k / 10), 0.7)), Z(lerp(-Lh + 0.002, -Lh - Lbol, k / 10))]);
    const tipP: P2 = [yc - 0.00045, Z(-e.d + 0.0005)];
    const low = (al: number) => Math.min(...support.map(([y, z]) => y * Math.cos(al) - z * Math.sin(al)));
    const tipY = (al: number) => tipP[0] * Math.cos(al) - tipP[1] * Math.sin(al);
    let tilt = 0;
    for (let al = 0; al > -0.1; al -= 0.0002) {
      tilt = al;
      if (tipY(al) <= low(al)) break;
    }
    s.push([0, 0.0002 - Math.min(low(tilt), tipY(tilt)), zm], [tilt, 0, 0]);
    // --- Blade ---
    const nL = [8, 24, 60][s.detail];
    const nX = [1, 3, 8][s.detail];
    const main: number[][] = [];
    const bevel: number[][] = [];
    for (let j = 0; j <= nL; j++) {
      const t = Math.pow(j / nL, 0.9);
      const z = Z(zH - t * Lb);
      const xs = spine(t);
      const xe = edgeX(t);
      const hh = Math.max(0, xs - xe);
      const ts = lerp(big ? 0.0027 : 0.0025, 0.0009, Math.pow(t, 1.2));
      const wb = Math.min(0.0012, 0.35 * hh);
      const tb = 0.00035 * Math.min(1, hh / 0.004) + 0.0001;
      const face = (k: number, sg: number): number[] => {
        const q = k / nX;
        const x = lerp(xs, xe + wb, q);
        // A faint hollow in the grind.
        const th = lerp(ts, tb, q) - 0.00015 * Math.sin(PI * q) * Math.min(1, hh / 0.01);
        return [x, yc + (sg * th) / 2, z];
      };
      const m: number[] = [];
      for (let k = nX; k >= 0; k--) m.push(...face(k, -1));
      m.push(xs + 0.00022, yc - ts / 4, z, xs + 0.00022, yc + ts / 4, z);
      for (let k = 0; k <= nX; k++) m.push(...face(k, 1));
      main.push(m);
      bevel.push([...face(nX, 1), xe, yc, z, ...face(nX, -1)]);
    }
    s.add(f.blade, ringMesh(main, { closed: false }));
    s.add(f.edge, ringMesh(bevel, { closed: false }));
    // --- Handle outline: its top in line with the spine, its butt hooking down a little. ---
    const xHi = (u: number) => xSp + 0.0002 - 0.0012 * Math.sin(PI * u);
    const hh = (u: number) => 0.0215 + 0.0035 * Math.sin(PI * u * 0.9) + 0.002 * smooth(0.85, 1, u);
    const xLo = (u: number) => xHi(u) - hh(u) - 0.004 * smooth(0.75, 1, u);
    const zU = (u: number) => lerp(zB, -0.0004, u);
    // --- Bolster: a loft from the handle's section down over the heel to the blade. ---
    const nb = [3, 6, 10][s.detail];
    const nr = [8, 16, 28][s.detail];
    const bol: number[][] = [];
    for (let j = 0; j <= nb; j++) {
      const q = j / nb;
      const lo = lerp(xLo(0) - 0.0003, 0.0002, smooth(0.25, 1, q));
      const hi = lerp(xHi(0) + 0.0003, xSp + 0.0006, q);
      const hy = lerp(0.0083, 0.0019, Math.pow(smooth(0.1, 1, q), 0.7));
      const ring: number[] = [];
      for (const [x, y] of superSec((hi - lo) / 2, hy, 3, nr, (hi + lo) / 2, yc)) ring.push(x, y, Z(lerp(zB + 0.002, zH, q)));
      bol.push(ring);
    }
    s.add(p.ssBright, outward(ringMesh(bol, { caps: true })));
    // --- Scales, top and bottom, domed across and rounded off at the butt; the tang between them. ---
    const nu = [3, 10, 28][s.detail];
    const nv = [2, 6, 12][s.detail];
    const ts = tsU;
    const dome = (v: number) => Math.pow(Math.max(0, 1 - (2 * v - 1) ** 2), 0.38);
    for (const sg of [1, -1]) {
      s.add(
        f.pom,
        gridMesh(
          nu,
          nv,
          (i, k) => {
            const u = i / nu;
            const v = k / nv;
            return [lerp(xLo(u), xHi(u), v), yc + sg * (0.0012 + ts(u) * dome(v)), Z(zU(u))];
          },
          sg,
        ),
      );
    }
    const top: number[] = [];
    const bot: number[] = [];
    const per: [number, number][] = [];
    for (let i = 0; i <= nu; i++) per.push([xHi(i / nu), i / nu]);
    for (let i = nu; i >= 0; i--) per.push([xLo(i / nu), i / nu]);
    for (const [x, u] of per) {
      top.push(x, yc + 0.0012, Z(zU(u)));
      bot.push(x, yc - 0.0012, Z(zU(u)));
    }
    s.add(f.blade, outwardOpen(ringMesh([top, bot], { closed: true }), [(xSp + 0.02) / 2, yc, Z(-Lh / 2)]));
    if (s.fine(1))
      for (const u of [0.2, 0.5, 0.8])
        for (const sg of s.fine(2) ? [1, -1] : [1]) {
          const y = yc + sg * (0.0012 + ts(u)) - sg * 0.0001;
          s.add(p.ssBright, cyl(s, 0.0024, 0.0006, "y", 16), [(xLo(u) + xHi(u)) / 2, y, Z(zU(u))]);
        }
    s.pop();
    return () => ({});
  };
}

/**
 * Locking tongs lying with one arm on the bench and the other above it: two pressed stainless channels on a hinge
 * pin, a torsion spring coiled about the pin, a pull ring through the hinge to lock them, scalloped cups at the
 * working end. The upper arm turns on joint "tongs_pinch"; open(frac) springs it 0.52 rad apart.
 */
const tongsBuilder: Builder = (s, e) => {
  const p = palette();
  const cx = e.w / 2;
  const zPin = -0.022;
  const yPin = 0.019;
  const zHeadC = -e.d + 0.03;
  const A = 0.023;
  const B = 0.028;
  const dc = 0.005;
  const nPath = [5, 12, 36][s.detail];
  const t = 0.0008;
  // The channel's flanges run out where the arm is pressed into its cup.
  const chan = (sg: number) => (j: number) => {
    const fd = lerp(0.0045, 0.0013, smooth(0.8, 0.95, j / nPath));
    const pts: P2[] = [
      [-0.0095, fd],
      [-0.0095, 0.0006],
      [-0.0089, 0],
      [0.0089, 0],
      [0.0095, 0.0006],
      [0.0095, fd],
      [0.0095 - t, fd],
      [0.0095 - t, 0.0012],
      [0.0089 - t, t],
      [-0.0089 + t, t],
      [-0.0095 + t, 0.0012],
      [-0.0095 + t, fd],
    ];
    const use = s.detail === 0 ? [pts[0], pts[2], pts[3], pts[5], pts[6], pts[8], pts[9], pts[11]] : pts;
    return use.map(([x, y]) => [x, y * sg] as P2);
  };
  /** A scalloped cup, its rim at yRim, dished toward `sg` (+1 opening up). */
  const head = (sh: Shop, y0: number, sg: number, zc: number, phase: number) => {
    const K = [1, 3, 7][sh.detail];
    const M = [10, 18, 48][sh.detail];
    const teeth = 9;
    const rings: number[][] = [];
    const surf = (k: number, off: number) => {
      const ring: number[] = [];
      for (let j = 0; j < M; j++) {
        const ph = (j / M) * TAU;
        const tooth = Math.cos(teeth * ph + phase);
        const q = k / K;
        const rr = q * (1 + 0.045 * tooth * q * q);
        const a = A * rr * Math.cos(ph);
        const b = B * rr * Math.sin(ph);
        const y = y0 + sg * (dc * (rr * rr) + 0.0012 * Math.max(0, tooth) * q ** 6) + off;
        ring.push(cx + a, y, zc + b);
      }
      return ring;
    };
    for (let k = 0; k <= K; k++) rings.push(surf(k, 0));
    for (let k = K; k >= 0; k--) rings.push(surf(k, -sg * t));
    sh.add(p.ssBright, outward(ringMesh(rings, { closed: true })));
  };
  // Lower arm: from under the pin down to the bench, along it to its cup.
  const lowerPath = bentPath(
    [
      [cx, yPin - 0.0058, zPin + 0.002],
      [cx, 0.0, zPin - 0.04],
      [cx, 0.0, zHeadC + B * 0.35],
    ],
    0.03,
    nPath,
  );
  s.add(p.ss, sweep(lowerPath, chan(1), [-1, 0, 0]));
  head(s, 0.0009, 1, zHeadC, 0);
  // The hinge: rolled knuckle, pin heads, the torsion coil, the pull ring.
  s.add(p.ss, cyl(s, 0.0058, 0.0205, "x", 20), [cx, yPin, zPin]);
  if (s.fine(1)) for (const sx of [-1, 1]) s.add(p.ssBright, ball(s, 0.0034, 12, PI / 2), [cx + sx * 0.0103, yPin, zPin], [0, 0, (-sx * PI) / 2], [1, 0.45, 1]);
  if (s.fine(2)) {
    const coil: V3[] = [];
    const turns = 4;
    const nC = turns * [6, 10, 16][s.detail];
    for (let k = 0; k <= nC; k++) {
      const a = (k / nC) * turns * TAU;
      coil.push([cx - 0.0175 + (0.0062 * k) / nC, yPin + 0.0068 * Math.cos(a), zPin + 0.0068 * Math.sin(a)]);
    }
    s.add(p.chrome, sweep(coil, () => circleSec(0.0007, s.fine(2) ? 6 : 4), [1, 0, 0]));
  }
  if (s.fine(1)) s.add(p.chrome, torus(s, 0.0098, 0.0012, "x", PI * 2, 40, 8), [cx, yPin, zPin + 0.0088]);
  // Upper arm and its cup, on the joint.
  const j = s.joint("tongs_pinch", [cx, yPin, zPin]);
  const upperPath = bentPath(
    [
      [0, 0.0058, 0.002],
      [0, 0.0315 - yPin, -0.04],
      [0, 0.0118 - yPin, zHeadC - zPin + B * 0.35],
    ].map(([x, y, z]) => [x + cx, y, z] as V3),
    0.03,
    nPath,
  ).map(([x, y, z]) => [x - cx, y, z] as V3);
  // The joint's frame has the pivot at its origin: lay the upper parts out about it.
  j.push([-cx, -yPin, -zPin]);
  j.add(
    p.ss,
    sweep(
      upperPath.map(([x, y, z]) => [x + cx, y + yPin, z + zPin] as V3),
      chan(-1),
      [-1, 0, 0],
    ),
  );
  head(j, 0.0118, -1, zHeadC, PI);
  j.pop();
  return (b) => {
    const jj = b.joints.get("tongs_pinch");
    return {
      open: (fr) => turnJoint(jj, "x", 0.52 * clamp01(fr)),
    };
  };
};

/**
 * A slotted fish turner: a thin flexible blade with an angled leading edge and five slots, an offset neck rising
 * to a riveted acetal handle that slopes down to rest its butt on the bench.
 */
const spatulaBuilder: Builder = (s, e) => {
  const p = palette();
  const f = finishes();
  const d = e.d;
  const tb = 0.0009;
  const sh = new THREE.Shape();
  // Outline in (x, z): the leading edge at the back, angled; one side straight, the other curved in to the neck.
  sh.moveTo(0.0035, -d + 0.016);
  sh.lineTo(0.0805, -d + 0.0012);
  sh.quadraticCurveTo(0.0815, -d + 0.07, 0.066, -d + 0.13);
  sh.quadraticCurveTo(0.058, -d + 0.158, 0.05, -d + 0.168);
  sh.lineTo(0.033, -d + 0.168);
  sh.quadraticCurveTo(0.024, -d + 0.15, 0.012, -d + 0.1);
  sh.quadraticCurveTo(0.002, -d + 0.06, 0.0035, -d + 0.016);
  if (s.fine(1)) {
    for (let k = 0; k < 5; k++) {
      const x0 = 0.017 + k * 0.0108;
      const za = -d + 0.032 - k * 0.002;
      const zb = -d + 0.112 - k * 0.007;
      const xb = x0 + 0.003 + k * 0.001;
      const hole = new THREE.Path();
      const w = 0.0021;
      hole.moveTo(x0 - w, za);
      hole.lineTo(xb - w, zb);
      hole.absarc(xb, zb, w, PI, 0, true);
      hole.lineTo(x0 + w, za);
      hole.absarc(x0, za, w, 0, PI, true);
      sh.holes.push(hole);
    }
  }
  const blade = plate(s, sh, tb);
  blade.rotateX(PI / 2);
  blade.translate(0, tb, 0);
  // The shape's y became z and its thickness went down; the turn above maps (x, y) to (x, z) with y flipped.
  const pa = blade.attributes.position as THREE.BufferAttribute;
  for (let i = 0; i < pa.count; i++) {
    // A faint lengthwise flex, rising toward the neck, and a ground leading edge.
    const z = pa.getZ(i);
    const q = clamp01((z + d) / 0.168);
    pa.setY(i, pa.getY(i) + 0.0012 * q * q);
  }
  blade.computeVertexNormals();
  s.add(p.ss, blade);
  // The neck: a flat bar from the blade's root, cranked up to the handle.
  const zA = -d + 0.205;
  const yA = 0.019;
  const neck = bentPath(
    [
      [0.0415, 0.0014, -d + 0.15],
      [0.0415, 0.0016, -d + 0.172],
      [0.0415, yA - 0.002, zA - 0.004],
      [0.0415, yA - 0.0025, zA + 0.03],
    ],
    0.012,
    [6, 12, 24][s.detail],
  );
  s.add(p.ss, sweep(neck, (jj) => superSec(lerp(0.0075, 0.006, jj / neck.length), 0.0011, 4, 12), [1, 0, 0]));
  // The handle: a rounded bar sloping to the bench, three rivets through it.
  const yB = 0.0085;
  const zEnd = -0.0012;
  const nh = [4, 10, 22][s.detail];
  const hPath: V3[] = [];
  for (let k = 0; k <= nh; k++) hPath.push([0.0415, lerp(yA, yB, k / nh), lerp(zA, zEnd, k / nh)]);
  const hw = (q: number) => lerp(0.0118, 0.0135, q) * (q > 0.94 ? Math.sqrt(Math.max(0.05, 1 - ((q - 0.94) / 0.06) ** 2)) : 1);
  const ht = (q: number) => lerp(0.0076, 0.0085, q) * (q > 0.94 ? Math.sqrt(Math.max(0.05, 1 - ((q - 0.94) / 0.06) ** 2)) : 1) * (q < 0.04 ? 0.92 + 2 * q : 1);
  s.add(f.pom, sweep(hPath, (jj) => superSec(hw(jj / nh), ht(jj / nh), 3, [8, 16, 24][s.detail]), [1, 0, 0]));
  if (s.fine(1)) {
    const slope = Math.atan2(yA - yB, zEnd - zA);
    for (const q of [0.2, 0.5, 0.8]) {
      const y = lerp(yA, yB, q) + ht(q) - 0.0001;
      s.add(p.ssBright, cyl(s, 0.0026, 0.0006, "y", 16), [0.0415, y, lerp(zA, zEnd, q)], [slope, 0, 0]);
    }
  }
  return () => ({});
};

/**
 * A twelve-inch balloon whisk at rest on the bench: six loops of spring wire crossing at the tip (twelve strands),
 * gathered into a ferrule on a sealed tubular handle with a hanging ring. Its axis tilts down to the handle's butt.
 */
const whiskBuilder: Builder = (s, e) => {
  const p = palette();
  const f = finishes();
  const R = e.h / 2 - 0.0012;
  const La = 0.18;
  const aW = 0.07;
  const aEnd = 0.29;
  const rF = 0.0068;
  const cx = e.w / 2;
  const zTip = -e.d + 0.0012;
  const slope = (0.0116 - (R + 0.0009)) / (aEnd - aW);
  const ang = Math.atan(slope);
  const dirY = Math.sin(ang);
  const dirZ = Math.cos(ang);
  const axis = (a: number): V3 => [cx, R + 0.0011 + (a - aW) * dirY, zTip + a * dirZ];
  // Unit vectors across the axis: x, and the one perpendicular to both (mostly up).
  const ux: V3 = [1, 0, 0];
  const uy: V3 = [0, dirZ, -dirY];
  const rOf = (a: number) => (a <= aW ? R * Math.sqrt(Math.max(0, 1 - ((aW - a) / aW) ** 2)) : lerp(R, rF, Math.pow((a - aW) / (La - aW), 1.5)));
  const loops = s.detail === 0 ? 4 : 6;
  const nA = [5, 12, 26][s.detail];
  const sec = circleSec(0.0008, [3, 5, 6][s.detail]);
  for (let k = 0; k < loops; k++) {
    const ph = (k / loops) * PI;
    const sc = 1 - 0.018 * (k % 3);
    const dir: V3 = [Math.cos(ph) * ux[0] + Math.sin(ph) * uy[0], Math.cos(ph) * ux[1] + Math.sin(ph) * uy[1], Math.cos(ph) * ux[2] + Math.sin(ph) * uy[2]];
    // Down one side from inside the ferrule to the tip, round it, and back up the other: the tip by angle, so
    // that the wire turns over in an even curve, the shank by length.
    const n1 = Math.max(3, Math.round(nA * 0.5));
    const half: [number, number][] = [];
    for (let i = 0; i <= n1; i++) {
      const ps = (i / n1) * (PI / 2);
      half.push([aW * (1 - Math.cos(ps)), R * Math.sin(ps)]);
    }
    for (let i = 1; i <= nA - n1; i++) {
      const a = aW + ((La + 0.008 - aW) * i) / (nA - n1);
      half.push([a, rOf(Math.min(a, La)) * (a > La ? 0.9 : 1)]);
    }
    const pts: V3[] = [];
    const put = (a0: number, r0: number, sg: number) => {
      const a = a0 + 0.0009 * k;
      const c = axis(a);
      const r = r0 * sc * sg;
      pts.push([c[0] + dir[0] * r, c[1] + dir[1] * r, c[2] + dir[2] * r]);
    };
    for (let i = half.length - 1; i >= 1; i--) put(half[i][0], half[i][1], 1);
    put(0, 0, 1);
    for (let i = 1; i < half.length; i++) put(half[i][0], half[i][1], -1);
    // The loop's plane normal keeps the section square to the wire.
    const nrm = new THREE.Vector3(...dir).cross(new THREE.Vector3(0, dirY, dirZ)).normalize();
    s.add(p.chrome, sweep(pts, () => sec, [nrm.x, nrm.y, nrm.z]));
  }
  // Ferrule and handle, turned about the axis (built about y, then laid along it).
  const n = [8, 16, 36][s.detail];
  const hp: P2[] = [
    [0, -0.002],
    [0.0062, -0.002],
    [0.0072, 0.0],
    [0.0098, 0.012],
    [0.0104, 0.0135],
    [0.0108, 0.016],
    [0.0116, 0.04],
    [0.0118, 0.075],
    [0.0114, 0.1],
    [0.0108, aEnd - La - 0.006],
    [0.0094, aEnd - La - 0.0015],
    [0.006, aEnd - La],
    [0, aEnd - La + 0.0006],
  ];
  const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), new THREE.Vector3(0, dirY, dirZ));
  const eul = new THREE.Euler().setFromQuaternion(q);
  const prof = decimate(hp, s.detail === 0 ? 2 : 1);
  s.add(f.spunBright, revolve(prof, n, undefined, f.spunBright), axis(La), [eul.x, eul.y, eul.z]);
  // The hanging ring through a lug at the butt.
  const ringAt = axis(aEnd + 0.0055);
  if (s.fine(1)) s.add(p.chrome, torus(s, 0.0058, 0.0011, "x", PI * 2, 28, 6), ringAt, [eul.x, 0, 0]);
  return () => ({});
};

/**
 * An eight-ounce ladle, drawn in one piece: a deep spun bowl with a rolled rim, upright on the bench, and its
 * handle rising from the rim in a swan neck, running out nearly level and ending in a hook.
 */
const ladleBuilder: Builder = (s, e) => {
  const p = palette();
  const f = finishes();
  const cx = e.w / 2;
  const cz = -e.d + e.w / 2;
  const Rr = e.w / 2 - 0.0012;
  const prof: P2[] = [
    [0, 0],
    [0.016, 0],
    [0.03, 0.0045],
    [0.042, 0.0142],
    [0.0482, 0.0272],
    [Rr - 0.0004, 0.04],
    [Rr - 0.0006, 0.0462],
    [Rr - 0.0002, 0.0478],
    [Rr + 0.0008, 0.0486],
    [Rr + 0.0012, 0.0476],
    [Rr + 0.0004, 0.0466],
    [Rr - 0.0014, 0.0466],
    [Rr - 0.0014, 0.04],
    [0.0474, 0.0272],
    [0.0412, 0.0148],
    [0.0296, 0.0054],
    [0.0158, 0.0008],
    [0, 0.0008],
  ];
  turned(s, f.spunBright, prof, s.seg(72, 10), undefined, [cx, 0, cz]);
  const top = e.h - 0.0013;
  const path = bentPath(
    [
      [cx, 0.0445, cz + Rr - 0.003],
      [cx, 0.06, cz + Rr + 0.012],
      [cx, top - 0.002, cz + Rr + 0.05],
      [cx, top, -0.03],
      [cx, top, -0.012],
      [cx, top - 0.012, -0.0018],
      [cx, top - 0.022, -0.009],
    ],
    0.01,
    [10, 24, 56][s.detail],
  );
  const nP = path.length;
  s.add(
    p.ss,
    sweep(
      path,
      (j) => {
        const q = j / (nP - 1);
        const w = q < 0.08 ? lerp(0.016, 0.0115, q / 0.08) : q < 0.75 ? lerp(0.0115, 0.0165, (q - 0.08) / 0.67) : lerp(0.0165, 0.011, (q - 0.75) / 0.25);
        return superSec(w / 2, 0.0011, 5, [8, 12, 16][s.detail]);
      },
      [1, 0, 0],
    ),
  );
  return () => ({});
};

/* ------------------------------------------------------------------------------------------------------------ */
/* A Bordeaux bottle of red, for deglazing                                                                        */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * 750 mL in dark green glass: a deep punt, straight shoulders, a long neck with its collar under a foil capsule, a
 * plain cream label. set("level", 0..1) is the fraction of 750 mL left, its surface at the height that holds it.
 */
const bottleBuilder: Builder = (s, e) => {
  const f = finishes();
  const cx = e.w / 2;
  const cz = -e.d / 2;
  const Rb = e.w / 2 - 0.0007;
  const Ri = Rb - 0.0025;
  const outer: P2[] = [
    [0, 0.022],
    [0.012, 0.019],
    [0.02, 0.0108],
    [0.0262, 0.0032],
    [0.0292, 0.0004],
    [0.031, 0],
    [0.0338, 0.0006],
    [Rb - 0.0011, 0.0028],
    [Rb, 0.0068],
    [Rb, 0.19],
    [Rb - 0.0008, 0.197],
    [0.0332, 0.2075],
    [0.0262, 0.2205],
    [0.0188, 0.2312],
    [0.0148, 0.2398],
    [0.0142, 0.25],
    [0.014, 0.2778],
    [0.0153, 0.2792],
    [0.0155, 0.2872],
    [0.0146, 0.2892],
    [0.0145, 0.2928],
  ];
  const inner: P2[] = [
    [0.0094, 0.2928],
    [0.0094, 0.276],
    [0.0095, 0.25],
    [0.0108, 0.2398],
    [0.0158, 0.2305],
    [0.0236, 0.219],
    [0.0302, 0.2065],
    [Ri - 0.0005, 0.197],
    [Ri, 0.19],
    [Ri, 0.0095],
    [Ri - 0.0022, 0.0052],
    [0.0262, 0.0058],
    [0.0206, 0.0132],
    [0.0124, 0.0222],
    [0, 0.0252],
  ];
  const st = [3, 2, 1][s.detail];
  const n = [7, 16, 48][s.detail];
  const glassProf = [...decimate(outer, st), ...decimate(inner, st)];
  s.add(f.bottle, revolve(glassProf, n, undefined, f.bottle, true), [cx, 0, cz]);
  // Foil capsule, a hair proud of the glass, down over the collar.
  const cap: P2[] = [
    [0.01445, 0.2455],
    [0.01465, 0.2455],
    [0.01465, 0.2776],
    [0.01585, 0.2792],
    [0.01605, 0.2872],
    [0.01515, 0.2896],
    [0.01505, 0.2938],
    [0.0128, 0.2948],
    [0, 0.2949],
    [0, 0.2944],
    [0.0127, 0.2943],
    [0.0149, 0.2934],
    [0.0149, 0.2896],
    [0.0158, 0.2872],
    [0.0156, 0.2794],
    [0.01445, 0.2778],
    [0.01445, 0.2455],
  ];
  const capL = s.detail === 0 ? [cap[1], cap[2], cap[4], cap[7], cap[8], cap[9], cap[10], cap[15]] : cap;
  s.add(f.foil, outward(revolve(capL.slice().reverse(), n)), [cx, 0, cz]);
  // Label: a band round the front of the body.
  const lab: number[][] = [];
  const nl = [6, 12, 28][s.detail];
  for (const y of [0.055, 0.148]) {
    const ring: number[] = [];
    for (let k = 0; k <= nl; k++) {
      const th = PI / 2 - 1.25 + (2.5 * k) / nl;
      ring.push(cx + (Rb + 0.0004) * Math.cos(th), y, cz + (Rb + 0.0004) * Math.sin(th));
    }
    lab.push(ring);
  }
  s.add(f.label, outwardOpen(ringMesh(lab, { closed: false }), [cx, 0.1, cz]));
  // The wine: the inside of the glass filled to its level; the level found from the volume it holds.
  const winner = inner
    .slice()
    .reverse()
    .map(([r, y]) => [Math.max(0, r - 0.0002), y + 0.0002] as P2);
  const vol: number[] = [0];
  for (let k = 1; k < winner.length; k++) {
    const [r0, y0] = winner[k - 1];
    const [r1, y1] = winner[k];
    vol.push(vol[k - 1] + (PI * (y1 - y0) * (r0 * r0 + r0 * r1 + r1 * r1)) / 3);
  }
  const yFull = (() => {
    // 750 mL: find the height that holds it (the punt displaces its own volume, already outside the profile).
    for (let k = 1; k < winner.length; k++)
      if (vol[k] >= 0.00075) {
        const t = (0.00075 - vol[k - 1]) / (vol[k] - vol[k - 1]);
        return lerp(winner[k - 1][1], winner[k][1], t);
      }
    return 0.25;
  })();
  const yAt = (frac: number) => {
    const v = frac * Math.min(0.00075, vol[vol.length - 1]);
    for (let k = 1; k < winner.length; k++)
      if (vol[k] >= v) {
        const t = (v - vol[k - 1]) / Math.max(1e-12, vol[k] - vol[k - 1]);
        return lerp(winner[k - 1][1], winner[k][1], clamp01(t));
      }
    return yFull;
  };
  const wineMesh = new THREE.Mesh(new THREE.BufferGeometry(), f.wine);
  wineMesh.castShadow = false;
  const wineProf = (y: number): P2[] => {
    const out: P2[] = [];
    for (let k = 0; k < winner.length; k++) {
      const [r, yy] = winner[k];
      if (yy < y) out.push([r, yy]);
      else {
        const [r0, y0] = winner[k - 1];
        out.push([lerp(r0, r, (y - y0) / (yy - y0)), y]);
        break;
      }
    }
    out.push([0, y]);
    // Inside to out along the punt, up the wall, in across the surface: counter-clockwise.
    return out.reverse();
  };
  let level = 0.7;
  const rebuild = () => {
    wineMesh.geometry.dispose();
    const y = Math.max(winner[0][1] + 0.0008, yAt(level));
    wineMesh.visible = level > 0.002;
    wineMesh.geometry = level > 0.002 ? outward(revolve(decimate(wineProf(y).reverse(), st), n)) : new THREE.BufferGeometry();
  };
  rebuild();
  s.extra(wineMesh, [cx, 0, cz]);
  s.tris += decimate(winner, st).length * n * 2;
  return () => ({
    set: (name, v) => {
      if (name === "level") {
        const nv = clamp01(v);
        if (Math.abs(nv - level) > 0.002) {
          level = nv;
          rebuild();
        }
      }
    },
  });
};

/* ------------------------------------------------------------------------------------------------------------ */
/* Registry                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

register({
  "NSF-SKILLET-8IN": skilletFixture("cast_iron"),
  "NSF-SKILLET-10IN": skilletFixture("tri_ply"),
  "NSF-SKILLET-12IN": skilletFixture("carbon_steel"),
  "NSF-SKILLET-14IN": skilletFixture("nonstick"),
  "NSF-STOCKPOT-8QT": stockpot,
  "NSF-STOCKPOT-20QT": stockpot,
  "NSF-STOCKPOT-40QT": stockpot,
  "NSF-HOTEL-PAN-FULL-2IN": panFixture(),
  "NSF-HOTEL-PAN-FULL-4IN": panFixture(),
  "NSF-HOTEL-PAN-FULL-6IN": panFixture(),
  "NSF-HOTEL-PAN-HALF-4IN": panFixture(),
  "ARCH-FOOD-STEAK-RIBEYE": steakBuilder,
  "ARCH-FOOD-SOURDOUGH-BOULE": bouleBuilder,
  "ARCH-FOOD-EGG-FRIED": eggBuilder,
  "ARCH-FOOD-FRY-LOAD": friesBuilder,
  "ARCH-FOOD-BRAISED-PAN-PORTION": braiseBuilder,
  "NSF-CHEF-KNIFE-8IN": chefKnife(8),
  "NSF-CHEF-KNIFE-10IN": chefKnife(10),
  "NSF-TONGS-12IN": tongsBuilder,
  "NSF-SPATULA-SLOTTED": spatulaBuilder,
  "NSF-WHISK-12IN": whiskBuilder,
  "NSF-LADLE-8OZ": ladleBuilder,
  "NSF-DEGLAZE-BOTTLE": bottleBuilder,
});

