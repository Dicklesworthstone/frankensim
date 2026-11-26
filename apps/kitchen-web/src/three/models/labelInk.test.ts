import { describe, expect, test } from "bun:test";
import { LABEL_TW, labelInk, labelPaper, type LabelPaper } from "./foodTexels";
import { clamp, lerp, rng, smooth } from "./texels";

// The wine label's roughness and normal maps from labelInk (the worker's bare-paper normal map, redone round the ink)
// must be the bytes labelMat made before, when it ran its own ink loop and then greyTex and normalTex over the whole
// label. The ink is drawn with fonts on a canvas, which bun has none of, so these inks are painted from seeded noise
// over the real paper: sparse and dense, every r + g + b from 0 to 765, ink on the edges and in the corners.

const W = LABEL_TW;

/** normalTexels as it was before patchNormalTexels shared its loop (texels.ts). */
function normalTexelsBefore(h: Float32Array, w: number, hh: number, k: number, wrap: boolean): Uint8ClampedArray {
  const d = new Uint8ClampedArray(w * hh * 4);
  const at = (i: number, j: number) => {
    if (wrap) return h[((j + hh) % hh) * w + ((i + w) % w)]!;
    return h[clamp(j, 0, hh - 1) * w + clamp(i, 0, w - 1)]!;
  };
  for (let j = 0; j < hh; j++)
    for (let i = 0; i < w; i++) {
      const dx = (at(i + 1, j) - at(i - 1, j)) * k;
      const dy = (at(i, j - 1) - at(i, j + 1)) * k;
      const l = Math.hypot(dx, dy, 1);
      const o = (j * w + i) * 4;
      d[o] = (0.5 - (0.5 * dx) / l) * 255;
      d[o + 1] = (0.5 - (0.5 * dy) / l) * 255;
      d[o + 2] = (0.5 + 0.5 / l) * 255;
      d[o + 3] = 255;
    }
  return d;
}

/** greyTexels (texels.ts), as greyTex used it. */
function greyTexelsBefore(v: ArrayLike<number>, w: number, hh: number): Uint8ClampedArray {
  const d = new Uint8ClampedArray(w * hh * 4);
  for (let k = 0; k < w * hh; k++) {
    const g = clamp(v[k]! * 255, 0, 255);
    d[k * 4] = d[k * 4 + 1] = d[k * 4 + 2] = g;
    d[k * 4 + 3] = 255;
  }
  return d;
}

/** labelMat's ink before labelInk: its loop over the read-back canvas, then greyTex(rough) and normalTex(hP, 0.5). */
function labelInkBefore(ink8: Uint8ClampedArray, height: Float32Array): { rough: Uint8ClampedArray; normal: Uint8ClampedArray; height: Float32Array } {
  const hP = height.slice();
  const rough = new Float32Array(W * W);
  for (let k = 0; k < W * W; k++) {
    const lum = (ink8[k * 4]! + ink8[k * 4 + 1]! + ink8[k * 4 + 2]!) / 765;
    const inked = smooth(0.72, 0.45, lum);
    rough[k] = lerp(0.88, 0.55, inked);
    hP[k]! -= 1.2 * inked;
  }
  return { rough: greyTexelsBefore(rough, W, W), normal: normalTexelsBefore(hP, W, W, 0.5, false), height: hP };
}

const bare = labelPaper();
const fresh = (): LabelPaper => ({ paper: bare.paper, height: bare.height.slice(), normal: bare.normal.slice() });

/** The paper with ink over it: `spots` of a random colour, `strokes` rectangles blended in at a random strength (as a
 * canvas's antialiased edges and translucent hatching are), and, if `edges`, the four edges and corners inked. Row 512
 * holds every r + g + b from 0 to 765. */
function inked(seed: number, spots: number, strokes: number, edges: boolean): Uint8ClampedArray {
  const r = rng(seed);
  const d = bare.paper.slice();
  const put = (k: number, R: number, G: number, B: number, a: number) => {
    d[k * 4] = lerp(d[k * 4]!, R, a);
    d[k * 4 + 1] = lerp(d[k * 4 + 1]!, G, a);
    d[k * 4 + 2] = lerp(d[k * 4 + 2]!, B, a);
  };
  for (let k = 0; k < W * W; k++) if (r() < spots) put(k, r() * 256, r() * 256, r() * 256, 1);
  const inks: [number, number, number][] = [
    [96, 22, 34],
    [34, 26, 24],
  ];
  for (let n = 0; n < strokes; n++) {
    const [R, G, B] = inks[n % 2]!;
    const x0 = Math.floor(r() * W);
    const y0 = Math.floor(r() * W);
    const x1 = Math.min(W, x0 + 1 + Math.floor(r() * (r() < 0.5 ? 4 : 60)));
    const y1 = Math.min(W, y0 + 1 + Math.floor(r() * (r() < 0.5 ? 4 : 60)));
    const a = r() < 0.4 ? 1 : r();
    for (let j = y0; j < y1; j++) for (let i = x0; i < x1; i++) put(j * W + i, R, G, B, a);
  }
  if (edges)
    for (let q = 0; q < W; q++) {
      const a = 0.2 + 0.8 * r();
      for (const k of [q, (W - 1) * W + q, q * W, q * W + W - 1]) put(k, 34, 26, 24, a);
    }
  for (let s = 0; s <= 765; s++) {
    const R = Math.min(255, s);
    const G = Math.min(255, s - R);
    const k = (512 * W + s) * 4;
    d[k] = R;
    d[k + 1] = G;
    d[k + 2] = s - R - G;
  }
  return d;
}

const same = (a: ArrayLike<number>, b: ArrayLike<number>) => {
  expect(a.length).toBe(b.length);
  let diff = 0;
  for (let k = 0; k < a.length; k++) if (!Object.is(a[k], b[k])) diff++;
  expect(diff).toBe(0);
};

describe("labelInk", () => {
  test("the worker's bare paper normal map is normalTex's of the bare relief", () => {
    same(bare.normal, normalTexelsBefore(bare.height, W, W, 0.5, false));
  });

  const cases: [string, () => Uint8ClampedArray][] = [
    ["bare paper", () => bare.paper.slice()],
    ["sparse ink, like the label's", () => inked(91, 0.002, 1400, false)],
    ["dense ink on the edges and corners", () => inked(92, 0.3, 5000, true)],
    ["ink everywhere", () => inked(93, 1, 0, true)],
  ];
  for (const [name, make] of cases)
    test(`${name}: roughness and normal maps byte for byte`, () => {
      const ink = make();
      const was = labelInkBefore(ink, bare.height);
      const t = fresh();
      const now = labelInk(ink, t);
      same(now.rough, was.rough);
      same(now.normal, was.normal);
      same(t.height, was.height);
    });
});
