import { describe, expect, test } from "bun:test";
import * as THREE from "three";
import { tuckUnder } from "./people";

/**
 * A tube round the y axis: `n` corners round (the first column repeated at the end, as a texture seam splits it),
 * rows from y0 to y1 every `dy`, normals out (or in), every corner on one bone.
 */
function tube(r: number, y0: number, y1: number, dy: number, bone: number, n = 24, inward = false): THREE.BufferGeometry {
  const rows = Math.round((y1 - y0) / dy) + 1;
  const pos: number[] = [];
  const nrm: number[] = [];
  const si: number[] = [];
  const sw: number[] = [];
  for (let j = 0; j < rows; j++)
    for (let i = 0; i <= n; i++) {
      const a = (i / n) * Math.PI * 2;
      pos.push(Math.sin(a) * r, y0 + j * dy, Math.cos(a) * r);
      const s = inward ? -1 : 1;
      nrm.push(Math.sin(a) * s, 0, Math.cos(a) * s);
      si.push(bone, 0, 0, 0);
      sw.push(1, 0, 0, 0);
    }
  const idx: number[] = [];
  for (let j = 0; j + 1 < rows; j++)
    for (let i = 0; i < n; i++) {
      const a = j * (n + 1) + i;
      const b = a + 1;
      const c = a + n + 1;
      const d = c + 1;
      idx.push(a, c, b, b, c, d);
    }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute("normal", new THREE.Float32BufferAttribute(nrm, 3));
  g.setAttribute("skinIndex", new THREE.BufferAttribute(Uint8Array.from(si), 4));
  g.setAttribute("skinWeight", new THREE.Float32BufferAttribute(sw, 4));
  g.setIndex(idx);
  return g;
}

const SKIN = 2;
const CLOTH = 5;
/** A sleeve from y 0 to 0.3, and skin inside it from -0.1 to 0.4. */
const sleeve = () => tube(0.1, 0, 0.3, 0.01, CLOTH);
const radius = (g: THREE.BufferGeometry, v: number) => Math.hypot(g.getAttribute("position").getX(v), g.getAttribute("position").getZ(v));

describe("tuckUnder", () => {
  test("drops only skin deep under the cloth, keeps skin by its openings and outside it", () => {
    const skin = tube(0.095, -0.1, 0.4, 0.01, SKIN);
    const before = skin.getIndex()!.count / 3;
    const r = tuckUnder(skin, sleeve());
    expect(r.dropped).toBeGreaterThan(0);
    expect(skin.getIndex()!.count / 3).toBe(before - r.dropped);
    const pos = skin.getAttribute("position");
    const ix = skin.getIndex()!;
    const kept = new Set<number>();
    for (let k = 0; k < ix.count; k++) kept.add(ix.getX(k));
    // Every row the openings could show is still there; the rows deep inside are gone.
    for (let v = 0; v < pos.count; v++) {
      const y = pos.getY(v);
      if (y < 0.025 || y > 0.275) expect(kept.has(v)).toBe(true);
      if (y > 0.045 && y < 0.255) expect(kept.has(v)).toBe(false);
    }
  });

  test("the texture seam is not an opening: the seam column goes with the rest", () => {
    const skin = tube(0.095, -0.1, 0.4, 0.01, SKIN);
    tuckUnder(skin, sleeve());
    const pos = skin.getAttribute("position");
    const ix = skin.getIndex()!;
    for (let k = 0; k < ix.count; k++) {
      const v = ix.getX(k);
      const y = pos.getY(v);
      // The seam is at angle 0 (x = 0, z > 0); nothing deep inside survives there either.
      if (Math.abs(pos.getX(v)) < 1e-6 && y > 0.045 && y < 0.255) throw new Error(`seam corner kept at y ${y}`);
    }
  });

  test("skin under the cloth takes the cloth's bones, fading back to its own at the openings", () => {
    const skin = tube(0.095, -0.1, 0.4, 0.01, SKIN);
    tuckUnder(skin, sleeve());
    const pos = skin.getAttribute("position");
    const si = skin.getAttribute("skinIndex");
    const sw = skin.getAttribute("skinWeight");
    const share = (v: number, bone: number) => {
      let w = 0;
      for (let k = 0; k < 4; k++) if (si.getComponent(v, k) === bone) w += sw.getComponent(v, k);
      return w;
    };
    for (let v = 0; v < pos.count; v++) {
      const y = pos.getY(v);
      let sum = 0;
      for (let k = 0; k < 4; k++) sum += sw.getComponent(v, k);
      expect(sum).toBeCloseTo(1, 6);
      // Outside the sleeve, and at its very edges, the skin keeps its own bone.
      if (y < -0.005 || y > 0.305) expect(share(v, SKIN)).toBeCloseTo(1, 6);
      // From 1.2 cm in, the cloth's.
      if (y > 0.015 && y < 0.285) expect(share(v, CLOTH)).toBeCloseTo(1, 6);
    }
  });

  test("skin lying on the cloth is set 4 mm under it, away from the openings; skin already deep stays put", () => {
    const close = tube(0.099, -0.1, 0.4, 0.01, SKIN);
    tuckUnder(close, sleeve(), 0.03);
    const deep = tube(0.09, -0.1, 0.4, 0.01, SKIN);
    tuckUnder(deep, sleeve(), 0.03);
    const pos = close.getAttribute("position");
    for (let v = 0; v < pos.count; v++) {
      const y = pos.getY(v);
      if (y > 0.015 && y < 0.285) expect(radius(close, v)).toBeCloseTo(0.096, 3);
      if (y < -0.005 || y > 0.305) expect(radius(close, v)).toBeCloseTo(0.099, 6);
      expect(radius(deep, v)).toBeCloseTo(0.09, 6);
    }
  });

  test("skin facing in (an inner face's partner) is left as it is", () => {
    const skin = tube(0.095, -0.1, 0.4, 0.01, SKIN, 24, true);
    const n = skin.getIndex()!.count;
    const r = tuckUnder(skin, sleeve());
    expect(r).toEqual({ moved: 0, dropped: 0 });
    expect(skin.getIndex()!.count).toBe(n);
  });
});
