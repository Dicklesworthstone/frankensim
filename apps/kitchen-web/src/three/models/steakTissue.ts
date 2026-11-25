/**
 * The ribeye's tissue, texel by texel: what the cut face and the side show raw (fat, lean, marbling, the cap and its
 * seam), with the relief and the crust's high spots, for smallwares' steak to paint and sear. Free of three.js and the
 * DOM, so a worker computes it while the page loads (steakTissue.worker.ts, steakPrefetch.ts); the same function on
 * the same inputs gives the same arrays wherever it runs.
 */
import { PI, TAU, Cells, clamp, fbm, hash2, lerp, noise2, rng, smooth, type P2 } from "./texels";
import { adiff } from "./texels";

/**
 * The rib eye in section, as the butcher cuts it across the grain: the eye (longissimus dorsi) filling the middle,
 * the cap (spinalis dorsi) a crescent round its back and one end, a kernel of fat in the seam between them, the
 * small complexus at the cap's end, intramuscular fat flecking both, and the fat edge outside the cap. Laid out in
 * polar form about the steak's middle: the outline R(θ), and inward from it the fat edge, the cap, the seam.
 */
export class Ribeye {
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
    // ((v % 1) + 1) % 1, without the float remainders where v is known small: for |v| < 1, v % 1 is v, and w = v + 1
    // (below 2 when v ≤ 0.5) leaves w - 1 (exact) or w itself as the remainder.
    const v = th / TAU;
    let fr: number;
    if (v > -1 && v <= 0.5) {
      const w = v + 1;
      fr = w >= 1 ? w - 1 : w;
    } else fr = ((v % 1) + 1) % 1;
    const f = fr * this.nT;
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
  /** `polar`, written into out[0] (angle) and out[1] (depth). */
  polarInto(x: number, z: number, out: Float64Array): void {
    const X = x / this.sx + this.ox;
    const Z = z / this.sz + this.oz;
    const th = Math.atan2(Z, X);
    out[0] = th;
    out[1] = this.look(this.Rt, th) - Math.hypot(X, Z);
  }
  fatW(th: number): number {
    return this.look(this.fatT, th);
  }
  capW(th: number): number {
    return this.look(this.capT, th);
  }
}

/** The side's profile: a superellipse of this exponent from the top junction to the bottom one. */
export const STEAK_SIDE_Q = 3;

/** The outline at the builder's resolution: N points round, their angles, outward normals, top edge radii, lengths. */
export function steakOutline(rib: Ribeye, detail: 0 | 1 | 2) {
  const N = [28, 64, 160][detail];
  const O: P2[] = [];
  const thO: number[] = [];
  for (let i = 0; i < N; i++) {
    thO.push((i / N) * TAU);
    O.push(rib.outline((i / N) * TAU));
  }
  const nrm: P2[] = O.map((_, i) => {
    const p0 = O[(i + N - 1) % N];
    const p1 = O[(i + 1) % N];
    const [tx, tz] = [p1[0] - p0[0], p1[1] - p0[1]];
    const l = Math.hypot(tx, tz);
    return [tz / l, -tx / l];
  });
  // The edge's rounding: soft over the fat cap, tighter on the lean; the bottom flatter, where it lies on the pan.
  const erT = thO.map((th) => 0.0034 + 0.0026 * smooth(0.0016, 0.0045, rib.fatW(th)));
  const plen = [0];
  for (let i = 1; i <= N; i++) plen.push(plen[i - 1] + Math.hypot(O[i % N][0] - O[i - 1][0], O[i % N][1] - O[i - 1][1]));
  const L = plen[N];
  return { N, thO, O, nrm, erT, plen, L };
}

export interface SteakTissue {
  /** The cut face's map: FW x FH texels, pxm texels a metre. */
  FW: number;
  FH: number;
  pxm: number;
  fF: Uint8Array;
  fEF: Uint8Array;
  fR: Uint8Array;
  fG: Uint8Array;
  fB: Uint8Array;
  fSpot: Uint8Array;
  fFine: Uint8Array;
  hRaw: Float32Array;
  hCrust: Float32Array;
  /** The side's map: SW texels round, SH down. */
  SW: number;
  SH: number;
  sR: Uint8Array;
  sG: Uint8Array;
  sB: Uint8Array;
  sF: Uint8Array;
  sEF: Uint8Array;
  sSpot: Uint8Array;
  sN: Uint8Array;
  sMem: Uint8Array;
  sH: Float32Array;
}

/** The tissue of a ribeye in a w x d x h box (the box allows for the thickening as it cooks) at a detail level. */
export function steakTissue(w: number, d: number, h: number, detail: 0 | 1 | 2): SteakTissue {
  const a = w / 2;
  const b = d / 2;
  const T = h / (1 + 0.35 * 0.3);
  const rib = new Ribeye(a, b, 5150);
  const nRel = noise2(611);
  const nFine = noise2(612);
  const nEdge = noise2(614);
  const nWarp = noise2(615);
  const nOn = noise2(616);
  const relief = (x: number, z: number) => 0.6 * fbm(nRel, x * 150 + 3, z * 150 + 7, 3) + 0.4 * fbm(nRel, x * 55 + 31, z * 55, 2);
  const spotOf = (x: number, z: number) => smooth(0.36, 0.66, relief(x, z));
  const { N, nrm, erT, plen, L } = steakOutline(rib, detail);
  const q = STEAK_SIDE_Q;

  // --- The tissue map of the cut, shared by both faces and (at the outline) by the side. ---
  const FW = detail === 2 ? 1024 : detail === 1 ? 512 : 192;
  const FH = Math.round((FW * b) / a);
  const pxm = FW / (2 * a);
  const NP = FW * FH;
  const fF = new Uint8Array(NP); // fat 0..255
  const fEF = new Uint8Array(NP); // the outer fat edge
  const fR = new Uint8Array(NP); // raw colour
  const fG = new Uint8Array(NP);
  const fB = new Uint8Array(NP);
  const fSpot = new Uint8Array(NP);
  const fFine = new Uint8Array(NP);
  const hRaw = new Float32Array(NP);
  const hCrust = new Float32Array(NP);
  const fine2 = detail >= 1;
  // Bundles of fibres seen end on (fascicles, ~1.6 mm), the secondary bundles whose walls carry the marbling
  // (~5 mm, smaller in the loose cap), fat cells for flecks, and crust blisters.
  const pad = 0.01;
  const cFas = fine2 ? new Cells(-a - pad, -b - pad, a + pad, b + pad, 0.0016, 31) : null;
  const cEye = new Cells(-a - pad, -b - pad, a + pad, b + pad, 0.0056, 37);
  const cCap = new Cells(-a - pad, -b - pad, a + pad, b + pad, 0.0042, 41);
  const cFl = new Cells(-a - pad, -b - pad, a + pad, b + pad, 0.0021, 43);
  const cCr = fine2 ? new Cells(-a - pad, -b - pad, a + pad, b + pad, 0.0011, 47) : null;
  const wo = new Float64Array(4);
  const [kx, kz] = rib.outline(rib.thK);
  const kd = rib.fatW(rib.thK) + rib.capW(rib.thK) + 0.0008;
  const kl = Math.hypot(kx, kz);
  const kcx = kx * (1 - kd / kl);
  const kcz = kz * (1 - kd / kl);
  const ta = Math.atan2(kz, kx) + PI / 2;
  const [xx, xz] = rib.outline(rib.thX);
  const xd = rib.fatW(rib.thX) + 0.0065;
  const xl = Math.hypot(xx, xz);
  const cxx = xx * (1 - xd / xl);
  const cxz = xz * (1 - xd / xl);
  // A few long intermuscular streaks running in from the cap's seam across the eye.
  const streaks: P2[][] = [];
  {
    const r = rng(9091);
    for (let i = 0; i < 5; i++) {
      const th = rib.thK + (r() - 0.5) * 1.4;
      const [ox, oz] = rib.outline(th);
      const d0 = rib.fatW(th) + rib.capW(th) + 0.001;
      let x = ox * (1 - d0 / Math.hypot(ox, oz));
      let z = oz * (1 - d0 / Math.hypot(ox, oz));
      let ang = Math.atan2(-z, -x) + (r() - 0.5) * 1.2;
      const pts: P2[] = [[x, z]];
      const len = 6 + Math.floor(r() * 16);
      for (let k = 0; k < len; k++) {
        ang += (r() - 0.5) * 0.7;
        x += Math.cos(ang) * 0.0012;
        z += Math.sin(ang) * 0.0012;
        pts.push([x, z]);
      }
      streaks.push(pts);
    }
  }
  let sx0 = Infinity;
  let sx1 = -Infinity;
  let sz0 = Infinity;
  let sz1 = -Infinity;
  for (const pts of streaks)
    for (const [x, z] of pts) {
      sx0 = Math.min(sx0, x - 0.002);
      sx1 = Math.max(sx1, x + 0.002);
      sz0 = Math.min(sz0, z - 0.002);
      sz1 = Math.max(sz1, z + 0.002);
    }
  // The streaks' segments, flattened, and each streak's box padded as the whole set's is.
  const segN = streaks.reduce((n, pts) => n + pts.length - 1, 0);
  const sgX = new Float64Array(segN);
  const sgZ = new Float64Array(segN);
  const sgDX = new Float64Array(segN);
  const sgDZ = new Float64Array(segN);
  const sgL2 = new Float64Array(segN);
  const sgK = new Float64Array(segN);
  const sgN = new Float64Array(segN);
  const stEnd: number[] = [];
  const stBox: number[] = [];
  {
    let o = 0;
    for (const pts of streaks) {
      let bx0 = Infinity;
      let bx1 = -Infinity;
      let bz0 = Infinity;
      let bz1 = -Infinity;
      for (let k = 0; k < pts.length - 1; k++, o++) {
        const [x0, z0] = pts[k];
        const [x1, z1] = pts[k + 1];
        sgX[o] = x0;
        sgZ[o] = z0;
        sgDX[o] = x1 - x0;
        sgDZ[o] = z1 - z0;
        sgL2[o] = sgDX[o] * sgDX[o] + sgDZ[o] * sgDZ[o];
        sgK[o] = k;
        sgN[o] = pts.length - 1;
      }
      for (const [x, z] of pts) {
        bx0 = Math.min(bx0, x - 0.002);
        bx1 = Math.max(bx1, x + 0.002);
        bz0 = Math.min(bz0, z - 0.002);
        bz1 = Math.max(bz1, z + 0.002);
      }
      stEnd.push(o);
      stBox.push(bx0, bx1, bz0, bz1);
    }
  }
  /**
   * The nearest streak to (x, z): its distance and how far along it (0..1), into sdOut. Only a distance under the
   * streak's width (at most 1.05 mm) is ever used; a streak whose box, padded by 2 mm, misses the point cannot be that
   * near, nor nearer than one that is, so it is passed over (the nearest found is the same whenever it matters).
   */
  const sdOut = new Float64Array(2);
  const streakD = (x: number, z: number): void => {
    let best = 1;
    let along = 0;
    if (!(x < sx0 || x > sx1 || z < sz0 || z > sz1))
      for (let s = 0, o = 0; s < stEnd.length; s++) {
        const end = stEnd[s];
        if (x < stBox[s * 4] || x > stBox[s * 4 + 1] || z < stBox[s * 4 + 2] || z > stBox[s * 4 + 3]) {
          o = end;
          continue;
        }
        for (; o < end; o++) {
          const x0 = sgX[o];
          const z0 = sgZ[o];
          const dx = sgDX[o];
          const dz = sgDZ[o];
          const t = clamp(((x - x0) * dx + (z - z0) * dz) / sgL2[o], 0, 1);
          const dd = Math.hypot(x - x0 - dx * t, z - z0 - dz * t);
          if (dd < best) {
            best = dd;
            along = (sgK[o] + t) / sgN[o];
          }
        }
      }
    sdOut[0] = best;
    sdOut[1] = along;
  };
  const cosTa = Math.cos(ta);
  const sinTa = Math.sin(ta);
  const po = new Float64Array(2);
  const tissue = (x: number, z: number, k: number) => {
    rib.polarInto(x, z, po);
    const th = po[0];
    const d0 = po[1];
    // Organic boundaries: each line between tissues wanders by a millimetre or so.
    const d = d0 + 0.0011 * (fbm(nEdge, x * 300, z * 300, 2) - 0.5);
    const fw = rib.fatW(th);
    const cw = rib.capW(th);
    let fat = 0;
    let mus = 0;
    let ef = 0;
    if (d < fw) {
      // The cap's fat edge is solid; the eye's free side has only a thin translucent silverskin.
      fat = fw < 0.0016 ? 0.5 : 1;
      ef = 1;
    } else if (cw > 0.0008 && d < fw + cw) mus = 1;
    else if (cw > 0.0008 && d < fw + cw + 0.0011 + 0.0012 * nEdge(th * 9, 3)) fat = 0.97;
    // The kernel: a lens of fat in the seam; the complexus: a small muscle at the cap's end, rimmed with fat.
    const ux = x - kcx;
    const uz = z - kcz;
    const al = ux * cosTa + uz * sinTa;
    const ac = -ux * sinTa + uz * cosTa;
    // (fbm is 0..1, so the lens's wander reaches 1.25 at most and the complexus's divisor 1.3: past those, no test.)
    const kq = (al / 0.0125) ** 2 + (ac / 0.0058) ** 2;
    if (kq < 1.3 && kq < 0.75 + 0.5 * fbm(nEdge, x * 160 + 40, z * 160, 2)) fat = 1;
    const cq = ((x - cxx) / 0.0078) ** 2 + ((z - cxz) / 0.0056) ** 2;
    if (cq < 1.6) {
      const cpx = cq / (0.7 + 0.6 * fbm(nEdge, x * 140, z * 140 + 70, 2));
      if (cpx < 1) {
        mus = 2;
        fat = Math.max(0, fat - 0.6);
      } else if (cpx < 1.22 && fat < 1) fat = Math.max(fat, 0.9 * smooth(0.3, 0.55, fbm(nEdge, x * 400, z * 400, 2)));
    }
    // Marbling: fat laid down in the walls between bundles of fibres, a lacy network, denser in the cap and toward
    // the eye's rim, the eye's heart leaner; flecks of fat cells between; a few long streaks.
    const wx = x + 0.0055 * (fbm(nWarp, x * 70, z * 70, 3) - 0.5) + 0.0016 * (fbm(nWarp, x * 260 + 9, z * 260, 2) - 0.5);
    const wz = z + 0.0055 * (fbm(nWarp, x * 70 + 17, z * 70, 3) - 0.5) + 0.0016 * (fbm(nWarp, x * 260, z * 260 + 9, 2) - 0.5);
    const dens = mus === 1 ? 0.62 : mus === 2 ? 0.3 : 0.2 + 0.34 * smooth(0.034, 0.008, d - fw - cw);
    let marb = 0;
    if (!ef) {
      // Whole walls between two bundles carry fat or none, so the threads run on and stop rather than crazing.
      (mus === 1 ? cCap : cEye).at(wx, wz, wo);
      const edge = (wo[1] - wo[0]) / 2;
      const pair = hash2(Math.round(Math.min(wo[2], wo[3]) * 1e7), Math.round(Math.max(wo[2], wo[3]) * 1e7), 61);
      const along = fbm(nOn, wx * 330 + 50, wz * 330, 2);
      const on = (pair < dens ? 1 : 0) * smooth(0.3, 0.45, along + 0.25 * (dens - pair));
      const wid = (0.00013 + 0.0005 * Math.pow(along, 2) * (0.7 + dens)) * (1 + 1.2 * smooth(0.0006, 0, wo[0] - 0.5 * wo[1]));
      marb = on * smooth(wid * 1.8, wid * 0.25, edge);
      // Fine threads off the main walls, fewer and thinner, and flecks of fat cells between.
      cFl.at(wx * 1.03 + 0.003, wz, wo);
      const e2 = (wo[1] - wo[0]) / 2;
      const p2 = hash2(Math.round(Math.min(wo[2], wo[3]) * 1e7), Math.round(Math.max(wo[2], wo[3]) * 1e7), 67);
      if (p2 < dens * 0.35) marb = Math.max(marb, smooth(0.00014, 0.00004, e2) * smooth(0.35, 0.6, along) * 0.8);
      const rf = 0.00014 + 0.00055 * hash2(Math.round(wo[2] * 1e6), 3, 5) * (0.4 + dens);
      if (wo[2] < dens * 0.3) marb = Math.max(marb, smooth(rf, rf * 0.4, wo[0]) * 0.9);
      if (mus === 0) {
        streakD(x, z);
        const sw = 0.0009 * (1 - sdOut[1]) + 0.00015;
        marb = Math.max(marb, smooth(sw, sw * 0.5, sdOut[0]));
      }
      fat = Math.max(fat, marb);
    }
    // The raw colour: cherry red, the cap a shade deeper and the complexus darker, each bundle a little different.
    let fas = 0;
    let fid = 0.5;
    if (cFas) {
      cFas.at(wx, wz, wo);
      fas = smooth(0.00011, 0.0, (wo[1] - wo[0]) / 2);
      fid = wo[2];
    }
    const lv = fbm(nFine, x * 60, z * 60, 3);
    const base0 = mus === 1 ? 124 : mus === 2 ? 106 : 134;
    const base1 = mus === 1 ? 20 : mus === 2 ? 17 : 22;
    const base2 = mus === 1 ? 30 : mus === 2 ? 26 : 30;
    const sh = (0.9 + 0.14 * lv + 0.04 * (fid - 0.5)) * (1 - 0.035 * fas);
    // Fat: creamy, layered in the cap's edge, pink where a thin fleck lets the lean show through.
    const fatTone = ef ? 1 - 0.04 * Math.sin(d * 2400 + 6 * fbm(nEdge, x * 90, z * 90, 2)) : 1;
    const pink = fw < 0.0016 && ef;
    const fa = fat;
    fR[k] = clamp(lerp(base0 * sh, pink ? 220 : 240 * fatTone, fa), 0, 255);
    fG[k] = clamp(lerp(base1 * sh, pink ? 196 : 231 * fatTone, fa), 0, 255);
    fB[k] = clamp(lerp(base2 * sh, pink ? 196 : 212 * fatTone, fa), 0, 255);
    fF[k] = Math.round(fat * 255);
    fEF[k] = ef ? 255 : 0;
    const fine = fbm(nFine, x * 850, z * 850, 3);
    const spot = spotOf(x, z);
    fSpot[k] = Math.round(spot * 255);
    fFine[k] = Math.round(fine * 255);
    // Relief: fat stands a hair proud of the lean, the bundles' walls are faint grooves. The crust: blistered and
    // craggy, its high spots where the pan touched hardest.
    hRaw[k] = 0.00006 * fat - 0.000018 * fas + 0.00005 * (lv - 0.5);
    if (cCr) cCr.at(x, z, wo);
    const blister = cCr ? smooth(0.0007, 0.0, wo[0]) * (wo[2] < 0.55 ? 1 : 0.3) : 0;
    hCrust[k] = 0.00014 * (fine - 0.5) + 0.00009 * blister + 0.0001 * (spot - 0.5);
  };
  for (let j = 0; j < FH; j++)
    for (let i = 0; i < FW; i++) tissue(-a + (i + 0.5) / pxm, -b + (1 - (j + 0.5) / FH) * 2 * b, j * FW + i);

  // --- The side: the tissue at the outline, carried down through the thickness (the fibres run through it). ---
  const SW = detail === 2 ? 1024 : detail === 1 ? 512 : 192;
  const SH = detail === 2 ? 128 : detail === 1 ? 64 : 24;
  const SP = SW * SH;
  const sR = new Uint8Array(SP);
  const sG = new Uint8Array(SP);
  const sB = new Uint8Array(SP);
  const sF = new Uint8Array(SP);
  const sEF = new Uint8Array(SP);
  const sSpot = new Uint8Array(SP);
  const sN = new Uint8Array(SP);
  const sMem = new Uint8Array(SP);
  const sH = new Float32Array(SP);
  const sideOfV = (v: number) => {
    // The plan offset of the side surface at height fraction v (0 at the bottom junction, 1 at the top).
    const sj = 2 * v - 1;
    return Math.pow(Math.max(0, 1 - Math.pow(Math.abs(sj), q)), 1 / q) - 1;
  };
  for (let i = 0; i < SW; i++) {
    const u = (i + 0.5) / SW;
    let lo = 0;
    while (lo < N - 1 && plen[lo + 1] / L < u) lo++;
    const f = (u * L - plen[lo]) / (plen[lo + 1] - plen[lo]);
    const th = ((lo + f) / N) * TAU;
    const [ox, oz] = rib.outline(th);
    const hi = (lo + 1) % N;
    const nx0 = lerp(nrm[lo][0], nrm[hi][0], f);
    const nz0 = lerp(nrm[lo][1], nrm[hi][1], f);
    const er = lerp(erT[lo], erT[hi], f);
    for (let j = 0; j < SH; j++) {
      const v = 1 - (j + 0.5) / SH;
      const k = j * SW + i;
      // The plan point under this texel, a hair inside the surface; the fibres wave a little through the height.
      const wav = 0.0016 * (fbm(nWarp, u * L * 140, v * 5, 2) - 0.5);
      const fwt = rib.fatW(th);
      const off = er * (v >= 0.5 ? 1 : 0.6) * sideOfV(v) - (fwt < 0.0016 ? fwt + 0.0007 : 0.0005);
      const px = ox + nx0 * off - nz0 * wav;
      const pz = oz + nz0 * off + nx0 * wav;
      const fx = clamp((px + a) * pxm - 0.5, 0, FW - 1.001);
      const fy = clamp((1 - (pz + b) / (2 * b)) * FH - 0.5, 0, FH - 1.001);
      const fi = Math.floor(fx);
      const fj = Math.floor(fy);
      const tx = fx - fi;
      const ty = fy - fj;
      const fk = fj * FW + fi;
      const bl = (A: Uint8Array) => lerp(lerp(A[fk], A[fk + 1], tx), lerp(A[fk + FW], A[fk + FW + 1], tx), ty);
      const stri = fbm(nFine, u * L * 1400, v * T * 140, 2);
      const g = 0.95 + 0.1 * stri;
      const ef = bl(fEF) / 255;
      // Marbling meets the side as soft streaks, not the cut face's crisp threads; the fat edge stays solid.
      const fat0 = bl(fF) / 255;
      const dash = smooth(0.42, 0.62, fbm(nOn, u * L * 300 + 7, v * T * 90, 2));
      const isFat = lerp(fat0 * 0.45 * dash, fat0, ef);
      const lr = (bl(fR) - fat0 * 240) / Math.max(0.05, 1 - fat0);
      const lg = (bl(fG) - fat0 * 231) / Math.max(0.05, 1 - fat0);
      const lb = (bl(fB) - fat0 * 212) / Math.max(0.05, 1 - fat0);
      const lean = fat0 > 0.95 ? [134, 22, 30] : [clamp(lr, 60, 200), clamp(lg, 10, 120), clamp(lb, 15, 120)];
      // The cap's fat in section: creamy, faintly yellow, with the grain of its connective sheets running across it.
      const fn = stri * 0.6 + 0.4 * fbm(nFine, u * L * 90 + 3, v * T * 30, 2);
      sR[k] = clamp(lerp(lean[0] * g, 238 - 10 * fn, isFat), 0, 255);
      sG[k] = clamp(lerp(lean[1] * g, 228 - 14 * fn, isFat), 0, 255);
      sB[k] = clamp(lerp(lean[2] * g, 208 - 26 * fn, isFat), 0, 255);
      sF[k] = Math.round(isFat * 255);
      sEF[k] = ef > 0.5 ? 255 : 0;
      sSpot[k] = fSpot[fk];
      sN[k] = Math.round(fbm(nEdge, u * L * 260, v * T * 260, 3) * 255);
      // Silverskin: a pearly film in patches over the lean of the eye's free side.
      sMem[k] = Math.round(smooth(0.48, 0.7, fbm(nOn, u * L * 90, v * T * 60, 3)) * (1 - isFat) * 255);
      sH[k] = 0.00005 * (stri - 0.5) + 0.00004 * isFat;
    }
  }
  return { FW, FH, pxm, fF, fEF, fR, fG, fB, fSpot, fFine, hRaw, hCrust, SW, SH, sR, sG, sB, sF, sEF, sSpot, sN, sMem, sH };
}

/** The buffers of a tissue, for a worker to hand them over without a copy. */
export const steakBuffers = (t: SteakTissue): ArrayBuffer[] =>
  [t.fF, t.fEF, t.fR, t.fG, t.fB, t.fSpot, t.fFine, t.hRaw, t.hCrust, t.sR, t.sG, t.sB, t.sF, t.sEF, t.sSpot, t.sN, t.sMem, t.sH].map((x) => x.buffer as ArrayBuffer);

const ready = new Map<string, SteakTissue>();
const keyOf = (w: number, d: number, h: number, detail: number) => `${w}|${d}|${h}|${detail}`;
/** A tissue computed elsewhere (in a worker) for these inputs. Read only, as every tissue is once made. */
export function provideSteakTissue(w: number, d: number, h: number, detail: 0 | 1 | 2, t: SteakTissue): void {
  ready.set(keyOf(w, d, h, detail), t);
}
/** The tissue for these inputs: the one a worker made, if it has arrived, or made here and now. */
export function steakTissueFor(w: number, d: number, h: number, detail: 0 | 1 | 2): SteakTissue {
  return ready.get(keyOf(w, d, h, detail)) ?? steakTissue(w, d, h, detail);
}
