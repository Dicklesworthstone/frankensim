import { describe, expect, test } from "bun:test";
import * as THREE from "three";
import { metresPerUv } from "./people";

// metresPerUv over the arrays must give the accessor loop's number to the bit.

/** The loop metresPerUv used before it read the arrays: three.js's accessors and Vector3. */
function reference(g: THREE.BufferGeometry): number {
  const pos = g.getAttribute("position");
  const uv = g.getAttribute("uv");
  if (!uv) return 1;
  const ix = g.index;
  const n = ix ? ix.count : pos.count;
  const at = (k: number) => (ix ? ix.getX(k) : k);
  let a3 = 0;
  let a2 = 0;
  const p = [new THREE.Vector3(), new THREE.Vector3(), new THREE.Vector3()];
  const u = [new THREE.Vector2(), new THREE.Vector2(), new THREE.Vector2()];
  const dp = new THREE.Vector3();
  const dq = new THREE.Vector3();
  for (let t = 0; t + 2 < n; t += 3) {
    for (let k = 0; k < 3; k++) {
      p[k]!.fromBufferAttribute(pos as THREE.BufferAttribute, at(t + k));
      u[k]!.fromBufferAttribute(uv as THREE.BufferAttribute, at(t + k));
    }
    a3 += dp.subVectors(p[1]!, p[0]!).cross(dq.subVectors(p[2]!, p[0]!)).length() / 2;
    a2 += Math.abs((u[1]!.x - u[0]!.x) * (u[2]!.y - u[0]!.y) - (u[2]!.x - u[0]!.x) * (u[1]!.y - u[0]!.y)) / 2;
  }
  return a2 > 0 ? Math.sqrt(a3 / a2) : 1;
}

/** `g` with its positions stored as the cooks' files store them: normalised signed 16-bit integers, interleaved. */
function quantizedPositions(g: THREE.BufferGeometry): THREE.BufferGeometry {
  const pos = g.getAttribute("position");
  const q = new Int16Array(pos.count * 4);
  for (let i = 0; i < pos.count; i++)
    for (let c = 0; c < 3; c++) q[i * 4 + c] = Math.round(THREE.MathUtils.clamp(pos.getComponent(i, c), -1, 1) * 32767);
  const out = g.clone();
  out.setAttribute("position", new THREE.InterleavedBufferAttribute(new THREE.InterleavedBuffer(q, 4), 3, 0, true));
  return out;
}

/** `g` with its UVs stored as a glTF file may store them: normalised 16-bit integers, interleaved. */
function quantizedUvs(g: THREE.BufferGeometry): THREE.BufferGeometry {
  const uv = g.getAttribute("uv");
  const q = new Uint16Array(uv.count * 4);
  for (let i = 0; i < uv.count; i++) {
    q[i * 4 + 1] = Math.round(THREE.MathUtils.clamp(uv.getX(i), 0, 1) * 65535);
    q[i * 4 + 2] = Math.round(THREE.MathUtils.clamp(uv.getY(i), 0, 1) * 65535);
  }
  const out = g.clone();
  out.setAttribute("uv", new THREE.InterleavedBufferAttribute(new THREE.InterleavedBuffer(q, 4), 2, 1, true));
  return out;
}

describe("metresPerUv", () => {
  test("matches the accessor loop", () => {
    const shapes: [string, THREE.BufferGeometry][] = [
      ["box", new THREE.BoxGeometry(0.6, 0.9, 0.4, 2, 3, 2)],
      ["sphere", new THREE.SphereGeometry(0.21, 13, 9)],
      ["cylinder", new THREE.CylinderGeometry(0.05, 0.08, 0.7, 19, 2)],
      ["torus", new THREE.TorusGeometry(0.3, 0.04, 7, 23)],
      ["box unindexed", new THREE.BoxGeometry(0.2, 0.3, 0.4).toNonIndexed()],
    ];
    for (const [name, g] of shapes) {
      g.applyMatrix4(new THREE.Matrix4().makeScale(1.3, 0.7, 2.1));
      for (const [kind, x] of [["float", g], ["quantized UVs", quantizedUvs(g)], ["quantized positions and UVs", quantizedUvs(quantizedPositions(g))]] as const) {
        const want = reference(x);
        const got = metresPerUv(x);
        expect(Object.is(got, want), `${name} ${kind}: ${got} vs ${want}`).toBe(true);
        expect(want).not.toBe(1);
      }
    }
    const noUv = new THREE.BoxGeometry(1, 1, 1);
    noUv.deleteAttribute("uv");
    expect(metresPerUv(noUv)).toBe(1);
    const flat = new THREE.PlaneGeometry(1, 1);
    flat.getAttribute("uv").array.fill(0.5);
    expect(Object.is(metresPerUv(flat), reference(flat))).toBe(true);
  });
});
