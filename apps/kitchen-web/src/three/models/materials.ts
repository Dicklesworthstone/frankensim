/**
 * The finishes a commercial kitchen is made of, each one material shared by every model that uses it so that a
 * whole room merges into one mesh per finish. Textures are drawn once on canvases from a seeded generator, so every
 * load looks the same.
 *
 * Stainless is the AISI 304 No. 4 finish NSF equipment ships in: a fine directional grain whose highlights stretch
 * across it (MeshPhysicalMaterial anisotropy along the texture's u, which `boxUVs` lays level on every face).
 */
import * as THREE from "three";

function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}
function canvas(w: number, h: number): [HTMLCanvasElement, CanvasRenderingContext2D] {
  const c = document.createElement("canvas");
  c.width = w;
  c.height = h;
  return [c, c.getContext("2d")!];
}
/**
 * Anisotropic filtering asked of every texture here. three clamps it to the renderer's own maximum when it uploads
 * (`capabilities.getMaxAnisotropy()`, 16 on every current desktop and phone GPU), so this is "the most there is".
 */
export const MAX_ANISOTROPY = 16;

function tex(c: HTMLCanvasElement, srgb: boolean, repeat = 1): THREE.CanvasTexture {
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(repeat, repeat);
  t.anisotropy = MAX_ANISOTROPY;
  t.needsUpdate = true;
  return t;
}

/** A standard normal deviate from a uniform source (Box–Muller). */
function gauss(r: () => number): number {
  return Math.sqrt(-2 * Math.log(Math.max(1e-9, r()))) * Math.cos(2 * Math.PI * r());
}

/** Scale a field to zero mean and unit spread. */
function standardise(f: Float32Array): Float32Array {
  let m = 0;
  for (let i = 0; i < f.length; i++) m += f[i]!;
  m /= f.length;
  let v = 0;
  for (let i = 0; i < f.length; i++) v += (f[i]! - m) ** 2;
  const s = Math.sqrt(v / f.length) || 1;
  for (let i = 0; i < f.length; i++) f[i] = (f[i]! - m) / s;
  return f;
}

/**
 * The No. 4 grain as a field W texels along the grain (u) by H across it (v), zero mean and unit spread. Every row of
 * texels is a line of polishing scratches: a few overlapping strokes, each from a few centimetres to most of a tile
 * long with tapered ends and a random depth, the odd one cut deep, over a fine pitting noise. Rows are independent,
 * so the grain averages away within a few rows instead of leaving stripes behind in the coarser mip levels.
 */
function grainField(r: () => number, W: number, H: number): Float32Array {
  const f = new Float32Array(W * H);
  for (let y = 0; y < H; y++) {
    const row = y * W;
    const base = 0.3 * gauss(r);
    for (let x = 0; x < W; x++) f[row + x] = base + 0.25 * gauss(r);
    const strokes = 3 + Math.floor(r() * 4);
    for (let k = 0; k < strokes; k++) {
      const len = Math.max(8, Math.round(W * (0.05 + 0.9 * r() * r())));
      const x0 = Math.floor(r() * W);
      const a = gauss(r) * (r() < 0.02 ? 1.7 : 0.75);
      for (let i = 0; i < len; i++) f[row + ((x0 + i) % W)] += a * Math.sqrt(Math.sin((Math.PI * (i + 0.5)) / len));
    }
  }
  return standardise(f);
}

/**
 * Soft, tiling cloudiness W by H: value noise on wrapped lattices of `cells` (one octave each, halving in weight), the
 * slow variation a rolled and polished sheet shows across a panel. Zero mean, unit spread.
 */
function mottleField(r: () => number, W: number, H: number, cells: readonly (readonly [number, number])[], stretch = 1): Float32Array {
  const f = new Float32Array(W * H);
  let weight = 1;
  for (const [cx, cy] of cells) {
    const lat = Array.from({ length: cx * cy }, () => gauss(r));
    const at = (i: number, j: number) => lat[(((j % cy) + cy) % cy) * cx + (((i % cx) + cx) % cx)]!;
    const ease = (t: number) => t * t * (3 - 2 * t);
    for (let y = 0; y < H; y++) {
      const v = (y / H) * cy;
      const j = Math.floor(v);
      const ty = ease(v - j);
      for (let x = 0; x < W; x++) {
        const u = (x / W) * cx;
        const i = Math.floor(u);
        const tx = ease(u - i);
        const a = at(i, j) + (at(i + 1, j) - at(i, j)) * tx;
        const b = at(i, j + 1) + (at(i + 1, j + 1) - at(i, j + 1)) * tx;
        f[y * W + x] += weight * (a + (b - a) * ty);
      }
    }
    weight *= 0.5 * stretch;
  }
  return standardise(f);
}

/**
 * Box-filter a W×H field one mip level down, to max(1, ⌊W/2⌋) × max(1, ⌊H/2⌋) as GL sizes its levels. Each texel
 * below averages the area of the level above that it covers, weighting the texels it only partly covers, so that an
 * odd size's spare row or column is shared out rather than read past or dropped.
 */
function halve(f: Float32Array, W: number, H: number): [Float32Array, number, number] {
  const w = Math.max(1, W >> 1);
  const h = Math.max(1, H >> 1);
  const sx = W / w;
  const sy = H / h;
  const o = new Float32Array(w * h);
  for (let y = 0; y < h; y++) {
    const y0 = y * sy;
    const y1 = y0 + sy;
    for (let x = 0; x < w; x++) {
      const x0 = x * sx;
      const x1 = x0 + sx;
      let s = 0;
      for (let j = Math.floor(y0); j < Math.min(H, Math.ceil(y1 - 1e-9)); j++) {
        const wy = Math.min(j + 1, y1) - Math.max(j, y0);
        for (let i = Math.floor(x0); i < Math.min(W, Math.ceil(x1 - 1e-9)); i++) s += wy * (Math.min(i + 1, x1) - Math.max(i, x0)) * f[j * W + i]!;
      }
      o[y * w + x] = s / (sx * sy);
    }
  }
  return [o, w, h];
}

/**
 * A grey tiling texture with its own mip chain. Level 0 is `mean + grain·g + mottle·m` (0..1); each level down
 * box-filters both fields, and the grain's share fades by `fade` per level on top of that. Seen from across the room a
 * brushed sheet is an even satin: the scratches are far below a pixel and only their average roughness is left. The
 * GPU's own mipmaps would keep each texel row's contrast through several levels, and those rows alias into regular
 * bands across every big stainless face.
 */
function mipTexture(W: number, H: number, grain: Float32Array, mottle: Float32Array | null, o: { mean: number; grain: number; mottle: number; fade: number; srgb: boolean; repeat: [number, number] }): THREE.DataTexture {
  const levels: { data: Uint8Array; width: number; height: number }[] = [];
  let g = grain;
  let m = mottle;
  let [w, h] = [W, H];
  for (let k = 0; ; k++) {
    const data = new Uint8Array(w * h * 4);
    const ga = o.grain * Math.pow(o.fade, k);
    for (let i = 0; i < w * h; i++) {
      const v = o.mean + ga * g[i]! + (m ? o.mottle * m[i]! : 0);
      const b = Math.max(0, Math.min(255, Math.round(v * 255)));
      data[i * 4] = data[i * 4 + 1] = data[i * 4 + 2] = b;
      data[i * 4 + 3] = 255;
    }
    levels.push({ data, width: w, height: h });
    if (w === 1 && h === 1) break;
    const [g2, w2, h2] = halve(g, w, h);
    if (m) m = halve(m, w, h)[0];
    g = g2;
    [w, h] = [w2, h2];
  }
  const t = new THREE.DataTexture(levels[0]!.data, W, H, THREE.RGBAFormat, THREE.UnsignedByteType);
  t.mipmaps = levels as unknown as typeof t.mipmaps;
  t.generateMipmaps = false;
  t.minFilter = THREE.LinearMipmapLinearFilter;
  t.magFilter = THREE.LinearFilter;
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(o.repeat[0], o.repeat[1]);
  t.colorSpace = o.srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.anisotropy = MAX_ANISOTROPY;
  t.needsUpdate = true;
  return t;
}

/** The mean of the brushed roughness map; a material's `roughness` is its effective roughness over this. */
const GRAIN_MEAN = 0.6;

/**
 * AISI 304 No. 4 brushed finish, carried mostly in roughness.
 *   rough  The grain: 256 texels along it by 1024 across over a 0.25 m tile (about 1 mm by 0.25 mm a texel), mean
 *          `GRAIN_MEAN`, ±0.08 a standard deviation, fading to 0.35 a mip level beyond the box filter's own averaging.
 *   map    Near white (the colour is the material's): a slow cloudiness across a 1 m tile, ±1.6 %, with the grain's
 *          long streaks laid faintly over it at ±0.5 % and faded the same way.
 */
function brushed(seed: number): { map: THREE.DataTexture; rough: THREE.DataTexture } {
  const r = rng(seed);
  const rough = mipTexture(256, 1024, grainField(r, 256, 1024), null, { mean: GRAIN_MEAN, grain: 0.08, mottle: 0, fade: 0.35, srgb: false, repeat: [4, 4] });
  const map = mipTexture(256, 256, grainField(r, 256, 256), mottleField(r, 256, 256, [[3, 5], [7, 11], [16, 23]], 0.9), {
    mean: 0.915,
    grain: 0.005,
    mottle: 0.016,
    fade: 0.35,
    srgb: true,
    repeat: [1, 1],
  });
  return { map, rough };
}

/** The same textures (one upload) shifted, so that two finishes sharing them never line their grain up. */
function shifted<T extends THREE.Texture>(t: T, du: number, dv: number): T {
  const c = t.clone();
  c.offset.set(du, dv);
  c.needsUpdate = true;
  return c;
}

/** Grain and pitting for cast iron and seasoned steel: dark speckle, rough. 0.25 m. */
function speckle(seed: number, base: string, dots: number, light: number): THREE.CanvasTexture {
  const r = rng(seed);
  const S = 256;
  const [c, x] = canvas(S, S);
  x.fillStyle = base;
  x.fillRect(0, 0, S, S);
  for (let i = 0; i < dots; i++) {
    const v = Math.floor(r() * light);
    x.fillStyle = `rgba(${v},${v},${v},${(0.1 + r() * 0.3).toFixed(2)})`;
    const s = 0.5 + r() * 1.8;
    x.fillRect(r() * S, r() * S, s, s);
  }
  return tex(c, true, 4);
}

/**
 * Oven porcelain enamel: a dark blue-grey glass fired onto steel, flecked with white and pale blue where the frit's
 * opacifier gathered, a fleck every half millimetre or so. 0.12 m a repeat; flecks are isotropic, so its own mips
 * average them cleanly into the enamel's mean blue-grey.
 */
function enamelFleck(): THREE.CanvasTexture {
  const r = rng(60013);
  const S = 256;
  const [c, x] = canvas(S, S);
  x.fillStyle = "#2f3845";
  x.fillRect(0, 0, S, S);
  // A slow cloud in the ground coat, then the flecks: many tiny pale ones, fewer larger white ones, some dark pits.
  // Every mark is drawn at its wrapped copies too, so the tile repeats without seams.
  const dot = (cx: number, cy: number, rad: number) => {
    for (const ox of [-S, 0, S])
      for (const oy of [-S, 0, S]) {
        if (cx + ox < -rad || cx + ox > S + rad || cy + oy < -rad || cy + oy > S + rad) continue;
        x.beginPath();
        x.arc(cx + ox, cy + oy, rad, 0, Math.PI * 2);
        x.fill();
      }
  };
  for (let i = 0; i < 60; i++) {
    x.fillStyle = r() < 0.5 ? "rgba(62,74,90,0.05)" : "rgba(30,34,42,0.06)";
    dot(r() * S, r() * S, 12 + r() * 26);
  }
  for (let i = 0; i < 2600; i++) {
    const big = r() < 0.05;
    const l = 140 + Math.floor(r() * 90);
    x.fillStyle = r() < 0.1 ? `rgba(14,16,20,0.55)` : `rgba(${l - 10},${l},${Math.min(255, l + 18)},${(big ? 0.6 : 0.18 + r() * 0.3).toFixed(2)})`;
    const sz = big ? 0.8 + r() * 0.7 : 0.4 + r() * 0.45;
    dot(r() * S, r() * S, sz);
  }
  return tex(c, true, 1);
}

/** Red-brown quarry tile, 6 in (152 mm) square, on 3 mm grout lines: colour and bump. 4 tiles per repeat. */
function quarry(): { map: THREE.CanvasTexture; bump: THREE.CanvasTexture; rough: THREE.CanvasTexture } {
  const r = rng(22913);
  const S = 512;
  const n = 4;
  const T = S / n;
  const g = Math.max(2, Math.round((0.003 / 0.152) * T));
  const [c, x] = canvas(S, S);
  const [cb, xb] = canvas(S, S);
  const [cr, xr] = canvas(S, S);
  x.fillStyle = "#5b544c";
  x.fillRect(0, 0, S, S);
  xb.fillStyle = "#202020";
  xb.fillRect(0, 0, S, S);
  xr.fillStyle = "#d8d8d8";
  xr.fillRect(0, 0, S, S);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++) {
      const hue = 9 + r() * 8;
      const lum = 20 + r() * 6;
      x.fillStyle = `hsl(${hue.toFixed(1)},${(46 + r() * 14).toFixed(0)}%,${lum.toFixed(1)}%)`;
      x.fillRect(i * T + g / 2, j * T + g / 2, T - g, T - g);
      for (let k = 0; k < 260; k++) {
        const v = r() < 0.5 ? 0 : 255;
        x.fillStyle = `rgba(${v},${v * 0.8},${v * 0.6},${(r() * 0.06).toFixed(3)})`;
        const s = 1 + r() * 3;
        x.fillRect(i * T + g / 2 + r() * (T - g), j * T + g / 2 + r() * (T - g), s, s);
      }
      xb.fillStyle = "#e0e0e0";
      xb.fillRect(i * T + g / 2 + 1, j * T + g / 2 + 1, T - g - 2, T - g - 2);
      xr.fillStyle = `rgb(${(110 + r() * 40) | 0},${(110 + r() * 40) | 0},${(110 + r() * 40) | 0})`;
      xr.fillRect(i * T + g / 2, j * T + g / 2, T - g, T - g);
    }
  return { map: tex(c, true), bump: tex(cb, false), rough: tex(cr, false) };
}

/** Pebbled fibreglass-reinforced plastic wall panel, the kind every kitchen is lined with. */
function pebble(): THREE.CanvasTexture {
  const r = rng(5501);
  const S = 256;
  const [c, x] = canvas(S, S);
  x.fillStyle = "#808080";
  x.fillRect(0, 0, S, S);
  for (let i = 0; i < 2600; i++) {
    const v = 90 + Math.floor(r() * 90);
    x.fillStyle = `rgba(${v},${v},${v},0.55)`;
    x.beginPath();
    x.arc(r() * S, r() * S, 0.8 + r() * 2.2, 0, Math.PI * 2);
    x.fill();
  }
  return tex(c, false, 8);
}

/** A round-hole perforation alpha map: `holes` across the texture, holes `frac` of the pitch wide. */
export function perforation(holes: number, frac: number, stagger = true): THREE.CanvasTexture {
  const S = 512;
  const [c, x] = canvas(S, S);
  x.fillStyle = "#fff";
  x.fillRect(0, 0, S, S);
  x.fillStyle = "#000";
  const p = S / holes;
  for (let j = 0; j <= holes; j++)
    for (let i = -1; i <= holes; i++) {
      const cx = (i + (stagger && j % 2 ? 0.5 : 0) + 0.5) * p;
      const cy = (j + 0.5) * p;
      x.beginPath();
      x.arc(cx, cy, (p * frac) / 2, 0, Math.PI * 2);
      x.fill();
    }
  return tex(c, false);
}

type Std = THREE.MeshStandardMaterial;
type Phys = THREE.MeshPhysicalMaterial;

export interface Palette {
  /** AISI 304 No. 4 brushed: cabinets, tops, splashes, hoods. */
  ss: Phys;
  /** The same sheet rougher and greyer: table undershelves, galvanised interior panels. */
  ssDull: Std;
  /** Bright polished trim: rails, bezels, handles. */
  ssBright: Std;
  /** Darkened heat-tinted steel inside ovens and under burners. */
  ssHeat: Std;
  /** Satin anodised aluminium: speed racks, sheet pans. */
  alu: Std;
  /** Chrome wire. */
  chrome: Std;
  /** Black cast iron: grates, burner heads, skillets. */
  castIron: Std;
  /** Carbon steel griddle plate, seasoned dark and smooth. */
  griddle: Std;
  /** Black powder coat and enamel. */
  enamel: Std;
  /** Charcoal powder coat. */
  charcoal: Std;
  /** Black rubber and EPDM gasket. */
  rubber: Std;
  /** Grey door gasket and the like. */
  gasket: Std;
  /** Black phenolic knobs and handles. */
  bakelite: Std;
  /** Black ABS and nylon plastics. */
  plastic: Std;
  /** White HDPE: cutting boards, carrier shells. */
  hdpe: Std;
  /** Grey polycarbonate and polyethylene: caster wheels, bins. */
  poly: Std;
  /** Clear tempered glass. */
  glass: Phys;
  /** Smoked oven-door glass. */
  glassSmoke: Phys;
  /** Brass, for gas fittings and valve bodies. */
  brass: Std;
  /** Copper refrigerant lines. */
  copper: Std;
  /** Cream insulated panel skin (walk-ins). */
  panel: Std;
  /** Pebbled white FRP wall lining. */
  frp: Std;
  /** Quarry tile floor. */
  tile: Std;
  /** Painted drywall and ceiling. */
  paint: Std;
  /** Indicator lamps and lit displays. */
  ledAmber: THREE.MeshBasicMaterial;
  ledGreen: THREE.MeshBasicMaterial;
  ledRed: THREE.MeshBasicMaterial;
  ledWhite: THREE.MeshBasicMaterial;
  /** A dark LCD with a faint glow. */
  screen: Std;
  /** Frying oil, amber and glossy. */
  oil: Phys;
  /** Water. */
  water: Phys;
  /** Fine heat-resistant ceramic (salamander and broiler radiants, element insulators). */
  ceramic: Std;
  /** Oven porcelain enamel: dark blue-grey, white-flecked, glossy. Oven cavities, door liners, broiler pans. */
  ovenEnamel: Std;
}

let cached: Palette | null = null;

/** The shared palette, made the first time it is asked for. */
export function palette(): Palette {
  if (cached) return cached;
  const b = brushed(104729);
  const bd = { map: shifted(b.map, 0.37, 0.61), rough: shifted(b.rough, 0.43, 0.29) };
  const ci = speckle(311, "#555", 2200, 120);
  const q = quarry();
  const peb = pebble();
  const std = (p: THREE.MeshStandardMaterialParameters, uv: "box" | "own" = "box"): Std => {
    const m = new THREE.MeshStandardMaterial(p);
    m.userData.uv = uv;
    return m;
  };
  const basic = (color: number) => {
    const m = new THREE.MeshBasicMaterial({ color });
    m.userData.uv = "own";
    return m;
  };
  // The grain runs along u (level on every vertical face, see shop.ts `boxUVs`). A real No. 4 sheet spreads its
  // reflection across the grain, averaging the room above and below the horizon into a light satin grey. The
  // renderer's image-based anisotropy does not integrate a stretched lobe: it bends the normal toward the stretch, so a
  // vertical stretch on a wall panel samples one band above or below the horizon and turns the steel gunmetal. Laid
  // along the grain, the bend keeps the reflection on the horizon, which matches the photographed tone. The base
  // colour, with the map's mean, gives stainless's reflectance (about 0.6).
  const ss = new THREE.MeshPhysicalMaterial({
    color: 0xd2d5d9,
    metalness: 1,
    roughness: 0.24 / GRAIN_MEAN,
    map: b.map,
    roughnessMap: b.rough,
    anisotropy: 0.45,
    anisotropyRotation: 0,
  });
  ss.userData.uv = "box";
  const glass = new THREE.MeshPhysicalMaterial({ color: 0xdfe8ee, metalness: 0, roughness: 0.04, transparent: true, opacity: 0.18, envMapIntensity: 1.4, depthWrite: false });
  glass.userData.uv = "own";
  const glassSmoke = new THREE.MeshPhysicalMaterial({ color: 0x1b1d20, metalness: 0, roughness: 0.05, transparent: true, opacity: 0.62, envMapIntensity: 1.3, clearcoat: 1, depthWrite: false });
  glassSmoke.userData.uv = "own";
  const oil = new THREE.MeshPhysicalMaterial({ color: 0xc58a1c, metalness: 0, roughness: 0.06, transmission: 0, transparent: true, opacity: 0.86, clearcoat: 1, clearcoatRoughness: 0.03 });
  oil.userData.uv = "own";
  const water = new THREE.MeshPhysicalMaterial({ color: 0x8fb4c4, metalness: 0, roughness: 0.03, transparent: true, opacity: 0.45, clearcoat: 1, depthWrite: false });
  water.userData.uv = "own";
  cached = {
    ss,
    ssDull: std({ color: 0xa9adb2, metalness: 0.85, roughness: 0.34 / GRAIN_MEAN, map: bd.map, roughnessMap: bd.rough }),
    ssBright: std({ color: 0xdfe2e6, metalness: 1, roughness: 0.14 }),
    ssHeat: std({ color: 0x3b3630, metalness: 0.75, roughness: 0.68, map: bd.map }),
    alu: std({ color: 0xc8ccd0, metalness: 0.9, roughness: 0.42, map: bd.map }),
    chrome: std({ color: 0xe8eaec, metalness: 1, roughness: 0.1 }),
    castIron: std({ color: 0x26262a, metalness: 0.55, roughness: 0.78, map: ci }),
    griddle: std({ color: 0x2c2b2a, metalness: 0.7, roughness: 0.42, map: ci }),
    enamel: std({ color: 0x141518, metalness: 0.2, roughness: 0.42 }),
    charcoal: std({ color: 0x34363a, metalness: 0.3, roughness: 0.55 }),
    rubber: std({ color: 0x18191b, metalness: 0, roughness: 0.92 }),
    gasket: std({ color: 0x45474b, metalness: 0, roughness: 0.85 }),
    bakelite: std({ color: 0x111213, metalness: 0.05, roughness: 0.3 }),
    plastic: std({ color: 0x1d1e21, metalness: 0, roughness: 0.55 }),
    hdpe: std({ color: 0xe9e8e2, metalness: 0, roughness: 0.62 }),
    poly: std({ color: 0x7d8186, metalness: 0, roughness: 0.6 }),
    glass,
    glassSmoke,
    brass: std({ color: 0xc9a258, metalness: 1, roughness: 0.3 }),
    copper: std({ color: 0xc27a4f, metalness: 1, roughness: 0.32 }),
    panel: std({ color: 0xe8e4da, metalness: 0.15, roughness: 0.55, bumpMap: peb, bumpScale: 0.3 }),
    frp: std({ color: 0xeeece6, metalness: 0, roughness: 0.68, bumpMap: peb, bumpScale: 0.6 }),
    tile: std({ color: 0xffffff, metalness: 0, roughness: 0.78, map: q.map, bumpMap: q.bump, bumpScale: 1.2, roughnessMap: q.rough }),
    paint: std({ color: 0xd8d6d0, metalness: 0, roughness: 0.9 }),
    ledAmber: basic(0xffb347),
    ledGreen: basic(0x4ade80),
    ledRed: basic(0xff4d4d),
    ledWhite: basic(0xf7fbff),
    screen: std({ color: 0x0b1118, metalness: 0.2, roughness: 0.15, emissive: 0x0d3a4a, emissiveIntensity: 0.55 }),
    oil,
    water,
    ceramic: std({ color: 0xd9d2c4, metalness: 0, roughness: 0.85 }),
    ovenEnamel: std({ color: 0xffffff, metalness: 0, roughness: 0.17, map: enamelFleck() }),
  };
  // The floor tile's texture repeats every 4 tiles: 0.608 m; the enamel's fleck every 0.12 m.
  cached.tile.userData.uvScale = 1 / 0.608;
  cached.ovenEnamel.userData.uvScale = 1 / 0.12;
  return cached;
}

/** Glowing parts whose heat is live: an element, a radiant, an ember bed. `heat` 0 is cold, 1 is bright orange (see `setGlow`). */
export function glowMaterial(base = 0x2a2522): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({ color: base, metalness: 0.3, roughness: 0.6, emissive: 0xff5a14, emissiveIntensity: 0 });
  m.userData.uv = "own";
  return m;
}
/**
 * Set a glow material's heat, 0 to 1: nothing below 0.2, then the colour of a body heating from about 900 K (the first
 * dull red) through cherry red (0.6) and orange (0.85) to orange-yellow at 1. The colour follows the Planck locus
 * (a close fit of it in sRGB); the brightness is held low enough that the renderer's filmic tone curve keeps the hue
 * saturated rather than bleaching it to cream, as hot metal looks to the eye.
 */
export function setGlow(m: THREE.MeshStandardMaterial, heat: number): void {
  const h = Math.max(0, Math.min(1, heat));
  // The surface's own colour and reflections, kept from the first call, fade as the glow takes over, so that a hot
  // element in a bright room reads as glowing metal, not as glow washed with the room's light.
  const base = (m.userData.glowBase ??= { color: m.color.getHex(), env: m.envMapIntensity }) as { color: number; env: number };
  const f = h <= 0.2 ? 0 : (h - 0.2) / 0.8;
  const keep = 1 - 0.7 * Math.min(1, f * 1.4);
  m.color.setHex(base.color).multiplyScalar(keep);
  m.envMapIntensity = base.env * keep;
  if (h <= 0.2) {
    m.emissiveIntensity = 0;
    return;
  }
  const T = 900 + 1100 * f;
  const c01 = (v: number) => Math.max(0, Math.min(1, v));
  const g = c01((99.4708 * Math.log(T / 100) - 161.1196) / 255);
  const b = T <= 1900 ? 0 : c01((138.5177 * Math.log(T / 100 - 10) - 305.0448) / 255);
  m.emissive.setRGB(1, g, b, THREE.SRGBColorSpace);
  m.emissiveIntensity = 0.06 + 1.5 * Math.pow(f, 1.7);
}
