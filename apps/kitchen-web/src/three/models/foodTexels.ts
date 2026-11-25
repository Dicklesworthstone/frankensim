/**
 * The pixel layers of the procedural textures that cost most to paint: the pan foods' (smallwares' panFoods) and the
 * wine label's paper (smallwares' labelMat). They are arithmetic on seeded noise and nothing else, so a worker
 * (foodTexels.worker.ts) makes them while the page loads, and the builders only copy them into canvases. Made here or
 * in the worker, the arrays are the same.
 */
import { apnoise, clamp01, fbm, greyTexels, lerp, noise2, normalTexels, paintTexels, smooth, tfbm, TAU } from "./texels";

/** The braised meat's tile, in texels a side. */
export const MEAT_TW = 512;
/** The label's side, in texels. */
export const LABEL_TW = 1024;

/** RGBA pixels of the pan foods' maps: sauce, veg and soup 256 a side, the meat's colour, normal and roughness MEAT_TW. */
export interface PanFoodTexels {
  sauce: Uint8ClampedArray;
  meat: Uint8ClampedArray;
  meatNormal: Uint8ClampedArray;
  meatRough: Uint8ClampedArray;
  veg: Uint8ClampedArray;
  soup: Uint8ClampedArray;
}

export function panFoodTexels(): PanFoodTexels {
  const n = noise2(808);
  // Sauce: a reduced braising liquid, glossy and nearly opaque, fat beads catching the light.
  const sauce = paintTexels(256, 256, (x, y, o) => {
    const v = fbm(n, x / 50, y / 50, 3);
    // Beads of rendered fat, amber, on the reduced liquor.
    const fat = smooth(0.8, 0.86, fbm(n, x / 4 + 50, y / 4, 2)) * 0.45;
    o[0] = lerp(48, 66, v) + fat * 50;
    o[1] = lerp(21, 28, v) + fat * 30;
    o[2] = lerp(9, 12, v) + fat * 8;
  });
  // Braised meat, a 5 cm tile with u along the grain: bundles of fibre a millimetre or two across, wandering, split
  // apart into dark crevices; a mahogany glaze over most of it, torn patches showing the redder, matt strands.
  const MT = MEAT_TW;
  const nWv = apnoise(8101, 4, 4);
  const nBd = apnoise(8102, 5, 36);
  const nSt = apnoise(8103, 8, 140);
  const nTr = apnoise(8104, 6, 6);
  const mC = new Uint8ClampedArray(MT * MT * 4);
  const mH = new Float32Array(MT * MT);
  const mR = new Float32Array(MT * MT);
  for (let j = 0; j < MT; j++)
    for (let i = 0; i < MT; i++) {
      const u = i / MT;
      const v = j / MT;
      const vw = v + 0.035 * (tfbm(nWv, u * 4, v * 4, 2) - 0.5);
      const bd = tfbm(nBd, u * 5, vw * 36, 2);
      const st = tfbm(nSt, u * 8, vw * 140, 2);
      const crev = smooth(0.3, 0.2, bd);
      const torn = 0.75 * smooth(0.64, 0.76, tfbm(nTr, u * 6, v * 6, 2));
      const k = j * MT + i;
      mH[k] = 0.55 * bd + 0.3 * st * (0.4 + 0.6 * torn) - 0.45 * crev;
      const glaze = 1 - torn;
      const lv = 0.92 + 0.16 * st;
      const dk = 1 - 0.5 * crev;
      mC[k * 4] = lerp(lerp(90, 108, st), lerp(56, 74, bd), glaze) * lv * dk;
      mC[k * 4 + 1] = lerp(lerp(46, 56, st), lerp(26, 34, bd), glaze) * lv * dk;
      mC[k * 4 + 2] = lerp(lerp(28, 34, st), lerp(13, 17, bd), glaze) * lv * dk;
      mC[k * 4 + 3] = 255;
      mR[k] = clamp01(lerp(0.42, 0.7, torn) + 0.2 * crev + 0.08 * (st - 0.5));
    }
  // Vegetables seen in the gaps: a green and orange medley under its butter gloss.
  const veg = paintTexels(256, 256, (x, y, o) => {
    const v = fbm(n, x / 6, y / 6, 3);
    const or = fbm(n, x / 9 + 30, y / 9, 3) > 0.6;
    o[0] = or ? 212 : lerp(40, 96, v);
    o[1] = or ? 112 : lerp(80, 140, v);
    o[2] = or ? 32 : lerp(20, 40, v);
  });
  // Soup: a vegetable broth, golden and a little cloudy, small eyes of fat on it, each paler with a dark rim.
  const soup = paintTexels(256, 256, (x, y, o) => {
    const v = fbm(n, x / 24, y / 24, 3);
    const f = fbm(n, x / 2.6 + 70, y / 2.6, 2) * smooth(0.45, 0.65, fbm(n, x / 40 + 9, y / 40, 2));
    const eye = smooth(0.66, 0.69, f);
    const rim = smooth(0.62, 0.66, f) - eye;
    o[0] = lerp(160, 184, v) + eye * 34 - rim * 50;
    o[1] = lerp(96, 116, v) + eye * 30 - rim * 36;
    o[2] = lerp(30, 40, v) + eye * 8 - rim * 12;
  });
  return { sauce, meat: mC, meatNormal: normalTexels(mH, MT, MT, 3, true), meatRough: greyTexels(mR, MT, MT), veg, soup };
}

/** The label's paper before any ink: its RGBA pixels and its relief, LABEL_TW a side. */
export interface LabelPaper {
  paper: Uint8ClampedArray;
  height: Float32Array;
}

export function labelPaper(): LabelPaper {
  const W = LABEL_TW;
  const nP = noise2(5551);
  const paper = new Uint8ClampedArray(W * W * 4);
  const hP = new Float32Array(W * W);
  for (let j = 0; j < W; j++)
    for (let i = 0; i < W; i++) {
      const k = j * W + i;
      // Laid lines about a millimetre apart, chain lines every 25 mm, the fibre's cloud over both.
      const laid = 0.5 + 0.5 * Math.cos((j / W) * 93 * TAU);
      const cd = (((i / W) * 93) % 25) - 12.5;
      const chain = Math.exp(-(cd * cd) / 0.5);
      const fib = fbm(nP, i / 6, j / 6, 3);
      hP[k] = 0.25 * laid + 0.3 * chain + 0.6 * fib;
      const v = (fib - 0.5) * 16 - 4 * laid - 6 * chain;
      paper[k * 4] = 238 + v;
      paper[k * 4 + 1] = 229 + v;
      paper[k * 4 + 2] = 204 + v * 1.1;
      paper[k * 4 + 3] = 255;
    }
  return { paper, height: hP };
}

/** The layers' buffers, to hand from a worker without a copy. */
export const panFoodBuffers = (t: PanFoodTexels): ArrayBuffer[] => [t.sauce, t.meat, t.meatNormal, t.meatRough, t.veg, t.soup].map((a) => a.buffer as ArrayBuffer);
export const labelPaperBuffers = (t: LabelPaper): ArrayBuffer[] => [t.paper.buffer as ArrayBuffer, t.height.buffer as ArrayBuffer];

// Layers made ahead by the worker, each for the one builder that paints from it (its material is made once a page).
let pan: PanFoodTexels | null = null;
let panTaken = false;
let label: LabelPaper | null = null;
let labelTaken = false;

export function providePanFoodTexels(t: PanFoodTexels): void {
  if (!panTaken) pan = t;
}
export function provideLabelPaper(t: LabelPaper): void {
  if (!labelTaken) label = t;
}
/** The pan foods' layers: the worker's if they have arrived, else made now. */
export function panFoodTexelsFor(): PanFoodTexels {
  const t = pan ?? panFoodTexels();
  pan = null;
  panTaken = true;
  return t;
}
/** The label's paper: the worker's if it has arrived, else made now. The caller owns it (labelMat inks its relief). */
export function labelPaperFor(): LabelPaper {
  const t = label ?? labelPaper();
  label = null;
  labelTaken = true;
  return t;
}
