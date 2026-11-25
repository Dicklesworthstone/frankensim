/**
 * The sourdough boule's static texel layers (smallwares' boule): its material coordinates, where the score runs, the
 * banneton's rings and the noise its crust is painted from. They depend on nothing but the texture's size, so a
 * worker (bouleTexels.worker.ts) makes them while the page loads, and the builder repaints the crust from them as the
 * loaf bakes. Built here or in the worker, the arrays are the same.
 */
import { TAU, clamp01, fbm, lerp, noise2, smooth } from "./texels";

/** The crust texture's side, by detail. */
export const BOULE_TW = [256, 512, 1024] as const;
/**
 * Material coordinates: s, arc length from the crown as a fraction of the meridian, and the azimuth θ. They are fixed
 * in the dough, so the crust's marks stretch with it as it springs. The texture is their polar map (a disc of radius
 * BOULE_D about the middle) with two strips at the sides for the ear's crust and its underside.
 */
export const BOULE_D = 0.42;
export const nLump = noise2(7101);
const nFine = noise2(7102);
const nBl = noise2(7103);
const nFl = noise2(7104);
const nSt = noise2(7105);
// The score: an arc in the material plane (p, q) = s (cos θ, sin θ), the ear on the side toward its centre, the
// opening (and the camera's usual view into it) on the far side, toward the front.
const BETA = 0.4;
export const cdx = -Math.sin(BETA);
export const cdz = -Math.cos(BETA);
export const C0 = 0.4;
export const RS = 0.46;
const PHI = 0.5;
export const PHIE = 0.44;
/** The score's (distance across it, angle along it) at material point (p, q), into out[0], out[1]. */
export const score = (p: number, q: number, out: Float64Array): void => {
  const wx = p - C0 * cdx;
  const wz = q - C0 * cdz;
  out[0] = Math.hypot(wx, wz) - RS;
  out[1] = Math.atan2(-cdx * wz + cdz * wx, -cdx * wx - cdz * wz);
};
export const gEnd = (phi: number) => {
  const t = Math.abs(phi) / PHI;
  return t >= 1 ? 0 : Math.sqrt(1 - t * t);
};
export const ringPhase = (s0: number, th: number) => s0 / 0.042 + th / TAU + 1.7 * (fbm(nLump, Math.cos(th) * s0 * 7 + 11, Math.sin(th) * s0 * 7, 3) - 0.5);
export const ringFade = (s0: number) => smooth(0.02, 0.07, s0) * (1 - smooth(0.56, 0.64, s0));

export interface BouleLayers {
  TW: number;
  lS: Float32Array;
  lD: Float32Array;
  lG: Float32Array;
  lRing: Uint8Array;
  lFine: Uint8Array;
  lBl: Uint8Array;
  lFl: Uint8Array;
  lSt: Uint8Array;
  lTear: Uint8Array;
  /** 0 loaf, 1 ear crust, 2 ear underside. */
  zone: Uint8Array;
}

/** The static layers, texel by texel, in material space. */
export function bouleLayers(TW: number): BouleLayers {
  const D = BOULE_D;
  const sc = new Float64Array(2);
  const NP = TW * TW;
  const lS = new Float32Array(NP);
  const lD = new Float32Array(NP);
  const lG = new Float32Array(NP);
  const lRing = new Uint8Array(NP);
  const lFine = new Uint8Array(NP);
  const lBl = new Uint8Array(NP);
  const lFl = new Uint8Array(NP);
  const lSt = new Uint8Array(NP);
  const lTear = new Uint8Array(NP);
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
      score(p, q, sc);
      const ds = sc[0];
      const phi = sc[1];
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
      // Where the flour in the rings tore apart as the dough stretched.
      lTear[k] = Math.round(fbm(nSt, p * 95 + 40, q * 95, 2) * 255);
    }
  return { TW, lS, lD, lG, lRing, lFine, lBl, lFl, lSt, lTear, zone };
}

/** The layers' buffers, to hand from a worker without a copy. */
export const bouleBuffers = (l: BouleLayers): ArrayBuffer[] =>
  [l.lS, l.lD, l.lG, l.lRing, l.lFine, l.lBl, l.lFl, l.lSt, l.lTear, l.zone].map((a) => a.buffer as ArrayBuffer);

/** Layers made ahead (by the worker), by texture side. They are only read, so one set serves every boule that size. */
const ready = new Map<number, BouleLayers>();
export function provideBouleLayers(l: BouleLayers): void {
  ready.set(l.TW, l);
}
export function bouleLayersFor(TW: number): BouleLayers {
  let l = ready.get(TW);
  if (!l) {
    l = bouleLayers(TW);
    ready.set(TW, l);
  }
  return l;
}
