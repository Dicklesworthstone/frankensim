/**
 * Deterministic procedural Canvas PBR texture generator for FrankenSim Kitchen Lab.
 *
 * Uses a seeded LCG (zero Math.random) to generate tactile commercial-kitchen surfaces:
 * brushed 304 stainless steel, commercial quarry tile floor, seasoned cast iron,
 * marbled ribeye & Maillard sear crust, floured/blistered sourdough crust, diamond-plate
 * walk-in aluminum, perforated rubber anti-fatigue mats, woven chef twill, and soft
 * multi-lobe culinary steam/smoke wisps.
 */

import * as THREE from "three";

function makeRng(seed0: number): () => number {
  let s = seed0 >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

function makeCanvas(w: number, h: number): [HTMLCanvasElement, CanvasRenderingContext2D] {
  const c = document.createElement("canvas");
  c.width = w;
  c.height = h;
  const ctx = c.getContext("2d")!;
  return [c, ctx];
}

function finalizeTex(
  canvas: HTMLCanvasElement,
  opts: { srgb?: boolean; repeat?: [number, number]; anisotropy?: number } = {},
): THREE.CanvasTexture {
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = opts.srgb === false ? THREE.NoColorSpace : THREE.SRGBColorSpace;
  tex.wrapS = THREE.RepeatWrapping;
  tex.wrapT = THREE.RepeatWrapping;
  if (opts.repeat) tex.repeat.set(opts.repeat[0], opts.repeat[1]);
  tex.anisotropy = opts.anisotropy ?? 4;
  tex.needsUpdate = true;
  return tex;
}

export interface PbrMapSet {
  readonly map: THREE.CanvasTexture;
  readonly bumpMap: THREE.CanvasTexture;
  readonly roughnessMap?: THREE.CanvasTexture;
}

let brushedSteelCache: PbrMapSet | null = null;
export function getBrushedSteelMaps(): PbrMapSet {
  if (brushedSteelCache) return brushedSteelCache;
  const rand = makeRng(104729);
  const w = 256;
  const h = 256;
  const [cCol, ctxCol] = makeCanvas(w, h);
  const [cBump, ctxBump] = makeCanvas(w, h);
  const [cRough, ctxRough] = makeCanvas(w, h);

  ctxCol.fillStyle = "#8c939b";
  ctxCol.fillRect(0, 0, w, h);
  ctxBump.fillStyle = "#808080";
  ctxBump.fillRect(0, 0, w, h);
  ctxRough.fillStyle = "#626262";
  ctxRough.fillRect(0, 0, w, h);

  // Fine directional anisotropic brushed grain streaks
  for (let i = 0; i < 1400; i += 1) {
    const y = rand() * h;
    const x = rand() * w;
    const len = 28 + rand() * 110;
    const v = 115 + Math.floor(rand() * 48);
    const a = 0.08 + rand() * 0.14;
    ctxCol.fillStyle = `rgba(${v},${v + 3},${v + 7},${a.toFixed(3)})`;
    ctxCol.fillRect(x, y, len, 1);
    if (x + len > w) ctxCol.fillRect(x - w, y, len, 1);

    const bv = 95 + Math.floor(rand() * 75);
    ctxBump.fillStyle = `rgba(${bv},${bv},${bv},0.28)`;
    ctxBump.fillRect(x, y, len, 1);

    const rv = 70 + Math.floor(rand() * 68);
    ctxRough.fillStyle = `rgba(${rv},${rv},${rv},0.22)`;
    ctxRough.fillRect(x, y, len, 1.5);
  }
  // Subtle panel sheen bands
  for (let i = 0; i < 24; i += 1) {
    const y = rand() * h;
    const grad = ctxCol.createLinearGradient(0, y, 0, y + 14);
    grad.addColorStop(0, "rgba(255,255,255,0)");
    grad.addColorStop(0.5, "rgba(255,248,240,0.05)");
    grad.addColorStop(1, "rgba(255,255,255,0)");
    ctxCol.fillStyle = grad;
    ctxCol.fillRect(0, y, w, 14);
  }

  brushedSteelCache = {
    map: finalizeTex(cCol, { srgb: true, repeat: [2, 2] }),
    bumpMap: finalizeTex(cBump, { srgb: false, repeat: [2, 2] }),
    roughnessMap: finalizeTex(cRough, { srgb: false, repeat: [2, 2] }),
  };
  return brushedSteelCache;
}

let quarryTileCache: PbrMapSet | null = null;
export function getQuarryTileMaps(): PbrMapSet {
  if (quarryTileCache) return quarryTileCache;
  const rand = makeRng(22913);
  const size = 256;
  const tiles = 4;
  const tilePx = size / tiles;
  const grout = 4;

  const [cCol, ctxCol] = makeCanvas(size, size);
  const [cBump, ctxBump] = makeCanvas(size, size);
  const [cRough, ctxRough] = makeCanvas(size, size);

  // Dark commercial charcoal epoxy grout bed
  ctxCol.fillStyle = "#191513";
  ctxCol.fillRect(0, 0, size, size);
  ctxBump.fillStyle = "#222222";
  ctxBump.fillRect(0, 0, size, size);
  ctxRough.fillStyle = "#d8d8d8";
  ctxRough.fillRect(0, 0, size, size);

  for (let ty = 0; ty < tiles; ty += 1) {
    for (let tx = 0; tx < tiles; tx += 1) {
      const x0 = tx * tilePx + grout * 0.5;
      const y0 = ty * tilePx + grout * 0.5;
      const tw = tilePx - grout;
      const th = tilePx - grout;

      // Authentic warm commercial red-charcoal quarry tile with subtle kiln variation
      const r = 64 + Math.floor(rand() * 14);
      const g = 34 + Math.floor(rand() * 9);
      const b = 28 + Math.floor(rand() * 8);
      ctxCol.fillStyle = `rgb(${r},${g},${b})`;
      ctxCol.fillRect(x0, y0, tw, th);

      // Subtle beveled tile highlight along top-left edge
      ctxCol.strokeStyle = "rgba(255,220,195,0.08)";
      ctxCol.lineWidth = 1;
      ctxCol.strokeRect(x0 + 0.5, y0 + 0.5, tw - 1, th - 1);

      ctxBump.fillStyle = "#c8c8c8";
      ctxBump.fillRect(x0, y0, tw, th);

      const rg = 138 + Math.floor(rand() * 28);
      ctxRough.fillStyle = `rgb(${rg},${rg},${rg})`;
      ctxRough.fillRect(x0, y0, tw, th);

      // Fine abrasive slip-resistant mineral grit speckles
      for (let k = 0; k < 42; k += 1) {
        const sx = x0 + 2 + rand() * (tw - 4);
        const sy = y0 + 2 + rand() * (th - 4);
        const dark = rand() > 0.5;
        ctxCol.fillStyle = dark ? "rgba(18,12,10,0.18)" : "rgba(210,165,135,0.09)";
        ctxCol.fillRect(sx, sy, 2, 2);
        ctxBump.fillStyle = dark ? "rgba(160,160,160,0.35)" : "rgba(235,235,235,0.35)";
        ctxBump.fillRect(sx, sy, 2, 2);
      }
    }
  }

  quarryTileCache = {
    map: finalizeTex(cCol, { srgb: true, repeat: [12, 12] }),
    bumpMap: finalizeTex(cBump, { srgb: false, repeat: [12, 12] }),
    roughnessMap: finalizeTex(cRough, { srgb: false, repeat: [12, 12] }),
  };
  return quarryTileCache;
}

let castIronCache: PbrMapSet | null = null;
export function getCastIronMaps(): PbrMapSet {
  if (castIronCache) return castIronCache;
  const rand = makeRng(88417);
  const size = 256;
  const [cCol, ctxCol] = makeCanvas(size, size);
  const [cBump, ctxBump] = makeCanvas(size, size);

  ctxCol.fillStyle = "#1c1918";
  ctxCol.fillRect(0, 0, size, size);
  ctxBump.fillStyle = "#808080";
  ctxBump.fillRect(0, 0, size, size);

  for (let i = 0; i < 6500; i += 1) {
    const x = rand() * size;
    const y = rand() * size;
    const rad = 0.7 + rand() * 1.9;
    const shade = 20 + Math.floor(rand() * 26);
    ctxCol.fillStyle = `rgba(${shade + 3},${shade},${shade},0.34)`;
    ctxCol.beginPath();
    ctxCol.arc(x, y, rad, 0, Math.PI * 2);
    ctxCol.fill();

    const bv = 70 + Math.floor(rand() * 120);
    ctxBump.fillStyle = `rgba(${bv},${bv},${bv},0.45)`;
    ctxBump.beginPath();
    ctxBump.arc(x, y, rad, 0, Math.PI * 2);
    ctxBump.fill();
  }

  castIronCache = {
    map: finalizeTex(cCol, { srgb: true, repeat: [3, 3] }),
    bumpMap: finalizeTex(cBump, { srgb: false, repeat: [3, 3] }),
  };
  return castIronCache;
}

export interface DualFoodMaps {
  readonly rawMap: THREE.CanvasTexture;
  readonly cookedMap: THREE.CanvasTexture;
  readonly bumpMap: THREE.CanvasTexture;
}

let steakMapsCache: DualFoodMaps | null = null;
export function getSteakMaps(): DualFoodMaps {
  if (steakMapsCache) return steakMapsCache;
  const rand = makeRng(51907);
  const size = 256;
  const [cRaw, ctxRaw] = makeCanvas(size, size);
  const [cCook, ctxCook] = makeCanvas(size, size);
  const [cBump, ctxBump] = makeCanvas(size, size);

  // Raw crimson-ruby ribeye muscle base
  const gRaw = ctxRaw.createRadialGradient(size * 0.48, size * 0.5, 16, size * 0.5, size * 0.5, size * 0.55);
  gRaw.addColorStop(0, "#8d2022");
  gRaw.addColorStop(0.7, "#76181b");
  gRaw.addColorStop(1, "#a62b2d");
  ctxRaw.fillStyle = gRaw;
  ctxRaw.fillRect(0, 0, size, size);

  // Rich mahogany caramelized Maillard sear crust base
  const gCook = ctxCook.createRadialGradient(size * 0.48, size * 0.5, 12, size * 0.5, size * 0.5, size * 0.56);
  gCook.addColorStop(0, "#3b1709");
  gCook.addColorStop(0.55, "#4d210d");
  gCook.addColorStop(0.88, "#6a3316");
  gCook.addColorStop(1, "#3e1c0c");
  ctxCook.fillStyle = gCook;
  ctxCook.fillRect(0, 0, size, size);

  ctxBump.fillStyle = "#7e7e7e";
  ctxBump.fillRect(0, 0, size, size);

  // Muscle fiber striations
  for (let i = 0; i < 420; i += 1) {
    const y = rand() * size;
    ctxRaw.strokeStyle = `rgba(110,18,22,${(0.16 + rand() * 0.22).toFixed(2)})`;
    ctxRaw.lineWidth = 1 + rand() * 1.6;
    ctxRaw.beginPath();
    ctxRaw.moveTo(0, y);
    ctxRaw.bezierCurveTo(size * 0.33, y + (rand() - 0.5) * 24, size * 0.66, y + (rand() - 0.5) * 24, size, y);
    ctxRaw.stroke();

    ctxBump.strokeStyle = `rgba(${100 + Math.floor(rand() * 70)},${100 + Math.floor(rand() * 70)},${100 + Math.floor(rand() * 70)},0.25)`;
    ctxBump.lineWidth = 1.2;
    ctxBump.stroke();
  }

  // Organic white intramuscular marbling veins (golden rendered fat when seared)
  for (let i = 0; i < 65; i += 1) {
    let x = rand() * size;
    let y = rand() * size;
    ctxRaw.strokeStyle = `rgba(246,232,218,${(0.38 + rand() * 0.42).toFixed(2)})`;
    ctxRaw.lineWidth = 0.8 + rand() * 2.2;
    ctxCook.strokeStyle = `rgba(155,104,52,${(0.26 + rand() * 0.25).toFixed(2)})`;
    ctxCook.lineWidth = 1.0 + rand() * 2.0;

    ctxRaw.beginPath();
    ctxCook.beginPath();
    ctxRaw.moveTo(x, y);
    ctxCook.moveTo(x, y);
    for (let s = 0; s < 5; s += 1) {
      x += (rand() - 0.45) * 32;
      y += (rand() - 0.5) * 24;
      ctxRaw.lineTo(x, y);
      ctxCook.lineTo(x, y);
    }
    ctxRaw.stroke();
    ctxCook.stroke();
  }

  // Sear micro-blisters & coarse black peppercorn / Maillard char crust speckles on cooked map
  for (let i = 0; i < 1200; i += 1) {
    const x = rand() * size;
    const y = rand() * size;
    const rad = 0.8 + rand() * 2.6;
    ctxCook.fillStyle = rand() > 0.35 ? "rgba(32,12,5,0.34)" : "rgba(135,68,26,0.28)";
    ctxCook.beginPath();
    ctxCook.arc(x, y, rad, 0, Math.PI * 2);
    ctxCook.fill();

    const bv = rand() > 0.5 ? 210 : 52;
    ctxBump.fillStyle = `rgba(${bv},${bv},${bv},0.28)`;
    ctxBump.beginPath();
    ctxBump.arc(x, y, rad, 0, Math.PI * 2);
    ctxBump.fill();
  }

  steakMapsCache = {
    rawMap: finalizeTex(cRaw, { srgb: true }),
    cookedMap: finalizeTex(cCook, { srgb: true }),
    bumpMap: finalizeTex(cBump, { srgb: false }),
  };
  return steakMapsCache;
}

let sourdoughMapsCache: DualFoodMaps | null = null;
export function getSourdoughMaps(): DualFoodMaps {
  if (sourdoughMapsCache) return sourdoughMapsCache;
  const rand = makeRng(73301);
  const size = 256;
  const [cRaw, ctxRaw] = makeCanvas(size, size);
  const [cBake, ctxBake] = makeCanvas(size, size);
  const [cBump, ctxBump] = makeCanvas(size, size);

  ctxRaw.fillStyle = "#e5dac3";
  ctxRaw.fillRect(0, 0, size, size);
  ctxBake.fillStyle = "#9b531b";
  ctxBake.fillRect(0, 0, size, size);
  ctxBump.fillStyle = "#888888";
  ctxBump.fillRect(0, 0, size, size);

  // Warm mottled blistered crust patches + fermentation micro-blisters
  for (let i = 0; i < 950; i += 1) {
    const x = rand() * size;
    const y = rand() * size;
    const r = 1.5 + rand() * 6.5;
    const dark = rand() > 0.45;
    ctxBake.fillStyle = dark ? "rgba(88,36,10,0.22)" : "rgba(206,134,54,0.22)";
    ctxBake.beginPath();
    ctxBake.arc(x, y, r, 0, Math.PI * 2);
    ctxBake.fill();

    const bv = dark ? 195 : 75;
    ctxBump.fillStyle = `rgba(${bv},${bv},${bv},0.28)`;
    ctxBump.beginPath();
    ctxBump.arc(x, y, r * 0.75, 0, Math.PI * 2);
    ctxBump.fill();
  }

  // Artisan rice-flour dusting halo + subtle banneton spiral rings
  for (let ring = 0; ring < 7; ring += 1) {
    const y = 18 + ring * 34;
    ctxRaw.strokeStyle = "rgba(255,252,245,0.45)";
    ctxRaw.lineWidth = 4;
    ctxRaw.beginPath();
    ctxRaw.moveTo(0, y);
    ctxRaw.lineTo(size, y);
    ctxRaw.stroke();

    ctxBake.strokeStyle = "rgba(245,232,210,0.34)";
    ctxBake.lineWidth = 3.2;
    ctxBake.beginPath();
    ctxBake.moveTo(0, y);
    ctxBake.lineTo(size, y);
    ctxBake.stroke();
  }

  // Fine flour speckles
  for (let i = 0; i < 1600; i += 1) {
    const x = rand() * size;
    const y = rand() * size;
    ctxBake.fillStyle = `rgba(250,242,226,${(0.06 + rand() * 0.16).toFixed(2)})`;
    ctxBake.fillRect(x, y, 1.5, 1.5);
  }

  sourdoughMapsCache = {
    rawMap: finalizeTex(cRaw, { srgb: true }),
    cookedMap: finalizeTex(cBake, { srgb: true }),
    bumpMap: finalizeTex(cBump, { srgb: false }),
  };
  return sourdoughMapsCache;
}

let diamondPlateCache: PbrMapSet | null = null;
export function getDiamondPlateMaps(): PbrMapSet {
  if (diamondPlateCache) return diamondPlateCache;
  const size = 256;
  const [cCol, ctxCol] = makeCanvas(size, size);
  const [cBump, ctxBump] = makeCanvas(size, size);

  ctxCol.fillStyle = "#5d6875";
  ctxCol.fillRect(0, 0, size, size);
  ctxBump.fillStyle = "#686868";
  ctxBump.fillRect(0, 0, size, size);

  const step = 32;
  for (let row = 0; row < size / step; row += 1) {
    for (let col = 0; col < size / step; col += 1) {
      const cx = (col + 0.5) * step;
      const cy = (row + 0.5) * step;
      const ang = (row + col) % 2 === 0 ? Math.PI / 4 : -Math.PI / 4;
      for (const [ctx, fill] of [
        [ctxCol, "#8896a6"],
        [ctxBump, "#ececec"],
      ] as const) {
        ctx.save();
        ctx.translate(cx, cy);
        ctx.rotate(ang);
        ctx.fillStyle = fill;
        ctx.beginPath();
        ctx.ellipse(0, 0, 10, 3.2, 0, 0, Math.PI * 2);
        ctx.fill();
        ctx.restore();
      }
    }
  }

  diamondPlateCache = {
    map: finalizeTex(cCol, { srgb: true, repeat: [5, 5] }),
    bumpMap: finalizeTex(cBump, { srgb: false, repeat: [5, 5] }),
  };
  return diamondPlateCache;
}

let antiFatigueMatCache: PbrMapSet | null = null;
export function getAntiFatigueMatMaps(): PbrMapSet {
  if (antiFatigueMatCache) return antiFatigueMatCache;
  const size = 128;
  const [cCol, ctxCol] = makeCanvas(size, size);
  const [cBump, ctxBump] = makeCanvas(size, size);

  ctxCol.fillStyle = "#171413";
  ctxCol.fillRect(0, 0, size, size);
  ctxBump.fillStyle = "#b0b0b0";
  ctxBump.fillRect(0, 0, size, size);

  const step = 16;
  for (let y = step / 2; y < size; y += step) {
    for (let x = step / 2; x < size; x += step) {
      ctxCol.fillStyle = "#090706";
      ctxCol.beginPath();
      ctxCol.arc(x, y, 4.4, 0, Math.PI * 2);
      ctxCol.fill();

      ctxBump.fillStyle = "#222222";
      ctxBump.beginPath();
      ctxBump.arc(x, y, 4.8, 0, Math.PI * 2);
      ctxBump.fill();
    }
  }

  antiFatigueMatCache = {
    map: finalizeTex(cCol, { srgb: true, repeat: [6, 2] }),
    bumpMap: finalizeTex(cBump, { srgb: false, repeat: [6, 2] }),
  };
  return antiFatigueMatCache;
}

let chefTwillBumpCache: THREE.CanvasTexture | null = null;
export function getChefTwillBump(): THREE.CanvasTexture {
  if (chefTwillBumpCache) return chefTwillBumpCache;
  const size = 128;
  const [c, ctx] = makeCanvas(size, size);
  ctx.fillStyle = "#808080";
  ctx.fillRect(0, 0, size, size);
  ctx.strokeStyle = "#9e9e9e";
  ctx.lineWidth = 1;
  for (let i = -size; i < size * 2; i += 4) {
    ctx.beginPath();
    ctx.moveTo(i, 0);
    ctx.lineTo(i + size, size);
    ctx.stroke();
  }
  ctx.strokeStyle = "#686868";
  for (let i = 0; i < size; i += 4) {
    ctx.beginPath();
    ctx.moveTo(0, i);
    ctx.lineTo(size, i);
    ctx.stroke();
  }
  chefTwillBumpCache = finalizeTex(c, { srgb: false, repeat: [8, 8] });
  return chefTwillBumpCache;
}

let wispyCloudCache: THREE.CanvasTexture | null = null;
export function getWispyCloudTexture(): THREE.CanvasTexture {
  if (wispyCloudCache) return wispyCloudCache;
  const rand = makeRng(99371);
  const size = 128;
  const [c, ctx] = makeCanvas(size, size);
  ctx.clearRect(0, 0, size, size);

  // Layered off-center overlapping soft clouds for organic curling steam wisps (zero polka-dot look!)
  for (let i = 0; i < 14; i += 1) {
    const ang = rand() * Math.PI * 2;
    const dist = rand() * size * 0.20;
    const cx = size * 0.5 + Math.cos(ang) * dist;
    const cy = size * 0.5 + Math.sin(ang) * dist;
    const rad = size * (0.20 + rand() * 0.22);
    const g = ctx.createRadialGradient(cx, cy, 1, cx, cy, rad);
    g.addColorStop(0, "rgba(255,255,255,0.22)");
    g.addColorStop(0.45, "rgba(255,255,255,0.11)");
    g.addColorStop(0.82, "rgba(255,255,255,0.03)");
    g.addColorStop(1, "rgba(255,255,255,0)");
    ctx.fillStyle = g;
    ctx.beginPath();
    ctx.arc(cx, cy, rad, 0, Math.PI * 2);
    ctx.fill();
  }
  wispyCloudCache = finalizeTex(c, { srgb: true });
  wispyCloudCache.wrapS = THREE.ClampToEdgeWrapping;
  wispyCloudCache.wrapT = THREE.ClampToEdgeWrapping;
  return wispyCloudCache;
}
