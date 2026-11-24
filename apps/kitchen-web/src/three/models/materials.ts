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
function tex(c: HTMLCanvasElement, srgb: boolean, repeat = 1): THREE.CanvasTexture {
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(repeat, repeat);
  t.anisotropy = 8;
  t.needsUpdate = true;
  return t;
}

/** No. 4 brushed grain: a colour map near white (the colour comes from the material) and a roughness map. 0.5 m. */
function brushed(seed: number, coarse: number): { map: THREE.CanvasTexture; rough: THREE.CanvasTexture } {
  const r = rng(seed);
  const S = 512;
  const [c, x] = canvas(S, S);
  const [cr, xr] = canvas(S, S);
  x.fillStyle = "#e9ebee";
  x.fillRect(0, 0, S, S);
  xr.fillStyle = "#7a7a7a";
  xr.fillRect(0, 0, S, S);
  for (let i = 0; i < 5200 * coarse; i++) {
    const y = r() * S;
    const x0 = r() * S;
    const len = 30 + r() * 260;
    const v = 200 + Math.floor(r() * 55);
    x.fillStyle = `rgba(${v},${v + 1},${v + 3},${(0.05 + r() * 0.12).toFixed(3)})`;
    x.fillRect(x0, y, len, r() < 0.8 ? 1 : 2);
    x.fillRect(x0 - S, y, len, 1);
    const rv = 90 + Math.floor(r() * 90);
    xr.fillStyle = `rgba(${rv},${rv},${rv},0.16)`;
    xr.fillRect(x0, y, len, 1);
    xr.fillRect(x0 - S, y, len, 1);
  }
  return { map: tex(c, true, 2), rough: tex(cr, false, 2) };
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
}

let cached: Palette | null = null;

/** The shared palette, made the first time it is asked for. */
export function palette(): Palette {
  if (cached) return cached;
  const b = brushed(104729, 1);
  const bd = brushed(7919, 0.7);
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
  const ss = new THREE.MeshPhysicalMaterial({
    color: 0xc9cdd2,
    metalness: 1,
    roughness: 0.34,
    map: b.map,
    roughnessMap: b.rough,
    anisotropy: 0.55,
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
    ssDull: std({ color: 0xa9adb2, metalness: 0.85, roughness: 0.5, map: bd.map, roughnessMap: bd.rough }),
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
  };
  // The floor tile's texture repeats every 4 tiles: 0.608 m.
  cached.tile.userData.uvScale = 1 / 0.608;
  return cached;
}

/** Glowing parts whose heat is live: an element, a radiant, an ember bed. `heat` 0 is cold, 1 is bright orange. */
export function glowMaterial(base = 0x2a2522): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({ color: base, metalness: 0.3, roughness: 0.6, emissive: 0xff5a14, emissiveIntensity: 0 });
  m.userData.uv = "own";
  return m;
}
/** Set a glow material's heat: blackbody-ish, dark red at 0.3, orange at 0.7, yellow-white near 1. */
export function setGlow(m: THREE.MeshStandardMaterial, heat: number): void {
  const h = Math.max(0, Math.min(1, heat));
  const r = 1;
  const g = 0.18 + 0.55 * h * h;
  const bl = 0.02 + 0.25 * Math.max(0, h - 0.75) * 4;
  m.emissive.setRGB(r, g, bl);
  m.emissiveIntensity = h <= 0.2 ? 0 : 2.6 * Math.pow((h - 0.2) / 0.8, 1.6);
}
