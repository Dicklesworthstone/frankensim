/**
 * The soft steam and smoke sprite used by every chapter's plumes: overlapping radial clouds drawn on a canvas from a
 * seeded generator, so it is the same on every load.
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
