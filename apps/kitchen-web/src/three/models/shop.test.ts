import { describe, expect, test } from "bun:test";
import * as THREE from "three";
import { RoundedBoxGeometry } from "three/examples/jsm/geometries/RoundedBoxGeometry.js";
import { applyMatrixInto, applyNormalInto, attributeValues, boundSphere, boxUVs, transformGeometry } from "./shop";

// The shop's fast paths must give three's own results to the bit: every float of positions, normals, UVs and
// bounding spheres, signed zeros included.

function geometries(): [string, THREE.BufferGeometry][] {
  const shape = new THREE.Shape();
  shape.moveTo(0, 0);
  shape.lineTo(0.3, 0);
  shape.quadraticCurveTo(0.4, 0.05, 0.35, 0.2);
  shape.lineTo(0, 0.25);
  const out: [string, THREE.BufferGeometry][] = [
    ["box", new THREE.BoxGeometry(0.6, 0.9, 0.4, 2, 3, 2)],
    ["cylinder", new THREE.CylinderGeometry(0.05, 0.08, 0.7, 19, 2, true)],
    ["sphere", new THREE.SphereGeometry(0.21, 13, 9)],
    ["torus", new THREE.TorusGeometry(0.3, 0.04, 7, 23)],
    ["extrude", new THREE.ExtrudeGeometry(shape, { depth: 0.12, bevelEnabled: true, bevelSize: 0.01, bevelThickness: 0.008, bevelSegments: 2 })],
    ["lathe", new THREE.LatheGeometry([new THREE.Vector2(0, 0), new THREE.Vector2(0.1, 0.02), new THREE.Vector2(0.12, 0.3)], 17)],
    ["plane", new THREE.PlaneGeometry(1.3, 0.7, 3, 2)],
    ["roundedBox", new RoundedBoxGeometry(0.5, 0.3, 0.2, 3, 0.03)],
  ];
  out.push(["box unindexed", new THREE.BoxGeometry(0.2, 0.3, 0.4).toNonIndexed()]);
  return out;
}

function matrices(): [string, THREE.Matrix4][] {
  const q = new THREE.Quaternion().setFromEuler(new THREE.Euler(0.3, -1.1, 2.4));
  const persp = new THREE.Matrix4().makePerspective(-0.4, 0.4, 0.3, -0.3, 0.5, 20);
  return [
    ["identity", new THREE.Matrix4()],
    ["translate", new THREE.Matrix4().makeTranslation(1.25, -0.5, 3.75)],
    ["rotate", new THREE.Matrix4().makeRotationFromQuaternion(q)],
    ["scale", new THREE.Matrix4().makeScale(1.7, 0.3, 2.9)],
    ["mirror", new THREE.Matrix4().makeScale(-1, 1, 1)],
    // Singular: the normal matrix is all zeros, so every normal is zero length.
    ["flatten", new THREE.Matrix4().makeScale(1, 0, 1)],
    ["compose", new THREE.Matrix4().compose(new THREE.Vector3(0.71, 0.92, -2.3), q, new THREE.Vector3(0.4, 1.6, 0.9))],
    ["perspective", persp.multiply(new THREE.Matrix4().makeTranslation(0, 0, -3))],
  ];
}

function same(a: ArrayLike<number>, b: ArrayLike<number>): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (!Object.is(a[i], b[i])) return false;
  return true;
}

function sameBox(a: THREE.Box3 | null, b: THREE.Box3 | null): boolean {
  if (!a || !b) return a === b;
  return same(a.min.toArray(), b.min.toArray()) && same(a.max.toArray(), b.max.toArray());
}

function sameSphere(a: THREE.Sphere | null, b: THREE.Sphere | null): boolean {
  if (!a || !b) return a === b;
  return Object.is(a.radius, b.radius) && Object.is(a.center.x, b.center.x) && Object.is(a.center.y, b.center.y) && Object.is(a.center.z, b.center.z);
}

describe("shop fast paths", () => {
  test("transformGeometry matches applyMatrix4", () => {
    for (const [gn, g] of geometries()) {
      // Plain float data, so the fast path (not its fallback to three) is what is compared.
      for (const name of ["position", "normal", "uv"]) {
        const a = g.attributes[name] as THREE.BufferAttribute;
        expect(a.array instanceof Float32Array && !a.normalized && !(a as unknown as THREE.InterleavedBufferAttribute).isInterleavedBufferAttribute, `${gn}: ${name} plain`).toBe(true);
      }
      expect(g.attributes.tangent, `${gn}: no tangents`).toBeUndefined();
      for (const [mn, m] of matrices()) {
        for (const bounded of [false, true]) {
          const ref = g.clone();
          const fast = g.clone();
          if (bounded) {
            for (const x of [ref, fast]) {
              x.computeBoundingBox();
              x.computeBoundingSphere();
            }
          }
          ref.applyMatrix4(m);
          transformGeometry(fast, m);
          const where = `${gn} × ${mn}${bounded ? " (bounded)" : ""}`;
          expect(same(fast.attributes.position!.array, ref.attributes.position!.array), `${where}: position`).toBe(true);
          expect(same(fast.attributes.normal!.array, ref.attributes.normal!.array), `${where}: normal`).toBe(true);
          expect(fast.attributes.position!.version, `${where}: position version`).toBe(ref.attributes.position!.version);
          expect(fast.attributes.normal!.version, `${where}: normal version`).toBe(ref.attributes.normal!.version);
          expect(sameBox(fast.boundingBox, ref.boundingBox), `${where}: box`).toBe(true);
          expect(sameSphere(fast.boundingSphere, ref.boundingSphere), `${where}: sphere`).toBe(true);
        }
      }
    }
  });

  test("transformGeometry matches applyMatrix4 without normals", () => {
    for (const [gn, g] of geometries()) {
      const bare = new THREE.BufferGeometry();
      bare.setAttribute("position", g.attributes.position!.clone());
      for (const [mn, m] of matrices()) {
        const ref = bare.clone();
        const fast = bare.clone();
        ref.computeBoundingSphere();
        fast.computeBoundingSphere();
        ref.applyMatrix4(m);
        transformGeometry(fast, m);
        expect(same(fast.attributes.position!.array, ref.attributes.position!.array), `${gn} × ${mn}: position`).toBe(true);
        expect(fast.attributes.normal, `${gn} × ${mn}: still no normals`).toBeUndefined();
        expect(sameSphere(fast.boundingSphere, ref.boundingSphere), `${gn} × ${mn}: sphere`).toBe(true);
      }
    }
  });

  test("boundSphere matches computeBoundingSphere", () => {
    for (const [gn, g] of geometries()) {
      for (const [mn, m] of matrices()) {
        const moved = g.clone().applyMatrix4(m);
        const ref = moved.clone();
        ref.computeBoundingSphere();
        const fast = moved.clone();
        boundSphere(fast);
        expect(sameSphere(fast.boundingSphere, ref.boundingSphere), `${gn} × ${mn}`).toBe(true);
      }
    }
    const empty = new THREE.BufferGeometry();
    empty.setAttribute("position", new THREE.BufferAttribute(new Float32Array(0), 3));
    const ref = empty.clone();
    ref.computeBoundingSphere();
    boundSphere(empty);
    expect(sameSphere(empty.boundingSphere, ref.boundingSphere), "empty").toBe(true);
  });

  test("boxUVs matches its accessor loop", () => {
    // The accessor loop boxUVs used before its fast path.
    const reference = (g: THREE.BufferGeometry, scale: number) => {
      const p = g.attributes.position!;
      const n = g.attributes.normal!;
      const uv = g.attributes.uv as THREE.BufferAttribute;
      for (let i = 0; i < p.count; i++) {
        const [ax, ay, az] = [Math.abs(n.getX(i)), Math.abs(n.getY(i)), Math.abs(n.getZ(i))];
        if (ay >= ax && ay >= az) uv.setXY(i, p.getX(i) * scale, p.getZ(i) * scale);
        else if (ax >= az) uv.setXY(i, p.getZ(i) * scale, p.getY(i) * scale);
        else uv.setXY(i, p.getX(i) * scale, p.getY(i) * scale);
      }
      uv.needsUpdate = true;
    };
    for (const [gn, g] of geometries()) {
      for (const [mn, m] of matrices()) {
        for (const scale of [1, 4, 0.3]) {
          const moved = g.clone().applyMatrix4(m);
          const ref = moved.clone();
          reference(ref, scale);
          const fast = moved.clone();
          boxUVs(fast, scale);
          expect(same(fast.attributes.uv!.array, ref.attributes.uv!.array), `${gn} × ${mn} × ${scale}`).toBe(true);
          expect(fast.attributes.uv!.version).toBe(ref.attributes.uv!.version);
        }
      }
    }
  });

  // The kinds of attribute a glTF file can hand three.js: every array type, normalised or not, plain or interleaved.
  const TYPES = [Float32Array, Uint32Array, Uint16Array, Uint8Array, Int32Array, Int16Array, Int8Array] as const;
  function fill(T: (typeof TYPES)[number], n: number): ArrayLike<number> & { [i: number]: number } {
    const a = new T(n);
    const lo = T === Float32Array ? -3 : T === Uint32Array || T === Uint16Array || T === Uint8Array ? 0 : -(2 ** (8 * T.BYTES_PER_ELEMENT - 1));
    const hi = T === Float32Array ? 3 : T === Uint32Array || T === Uint16Array || T === Uint8Array ? 2 ** (8 * T.BYTES_PER_ELEMENT) - 1 : 2 ** (8 * T.BYTES_PER_ELEMENT - 1) - 1;
    let r = 12345;
    for (let i = 0; i < n; i++) {
      r = (Math.imul(r, 1103515245) + 12345) >>> 0;
      // The ends of the range, zero, and values between.
      a[i] = i % 7 === 0 ? lo : i % 7 === 1 ? hi : i % 7 === 2 ? 0 : T === Float32Array ? lo + (r / 2 ** 32) * (hi - lo) : Math.floor(lo + (r / 2 ** 32) * (hi - lo + 1));
    }
    return a;
  }
  function attributes(T: (typeof TYPES)[number], itemSize: number, normalized: boolean): [string, THREE.BufferAttribute | THREE.InterleavedBufferAttribute][] {
    const count = 53;
    const plainA = new THREE.BufferAttribute(fill(T, count * itemSize) as THREE.TypedArray, itemSize, normalized);
    // Interleaved: two attributes in one buffer, this one second, with a spare component after it.
    const stride = itemSize + 3;
    const buf = new THREE.InterleavedBuffer(fill(T, count * stride) as THREE.TypedArray, stride);
    const inter = new THREE.InterleavedBufferAttribute(buf, itemSize, 2, normalized);
    return [
      ["plain", plainA],
      ["interleaved", inter],
    ];
  }

  test("attributeValues reads as three's accessors do", () => {
    const get = (a: THREE.BufferAttribute | THREE.InterleavedBufferAttribute, i: number, c: number) => [a.getX(i), a.getY(i), a.getZ(i), a.getW(i)][c]!;
    for (const T of TYPES)
      for (const normalized of [false, true])
        for (const itemSize of [1, 2, 3, 4])
          for (const [kind, a] of attributes(T, itemSize, normalized))
            for (let size = 1; size <= itemSize; size++) {
              const v = attributeValues(a, size);
              const ref: number[] = [];
              for (let i = 0; i < a.count; i++) for (let c = 0; c < size; c++) ref.push(get(a, i, c));
              expect(same(v, ref), `${T.name} ${normalized ? "normalised" : "raw"} ${kind} itemSize ${itemSize} size ${size}`).toBe(true);
            }
    const half = new THREE.Float16BufferAttribute(new Float32Array([0.5, -1.25, 3, 1e-3, 65504, -0]), 3);
    expect(same(attributeValues(half, 3), [half.getX(0), half.getY(0), half.getZ(0), half.getX(1), half.getY(1), half.getZ(1)]), "half floats").toBe(true);
    // Types three.js cannot denormalise are refused, as three.js refuses them.
    const odd = new THREE.BufferAttribute(new Float64Array(6), 3, true);
    expect(() => odd.getX(0)).toThrow("Invalid component type.");
    expect(() => attributeValues(odd, 3)).toThrow("Invalid component type.");
  });

  test("applyMatrixInto and applyNormalInto match Vector3 over quantized data", () => {
    const v = new THREE.Vector3();
    const k = new THREE.Matrix3();
    for (const T of [Float32Array, Int16Array, Int8Array, Uint16Array] as const)
      for (const [kind, a] of attributes(T, 3, T !== Float32Array))
        for (const [mn, m] of matrices()) {
          const ref = new Float32Array(a.count * 3);
          for (let i = 0; i < a.count; i++) v.fromBufferAttribute(a, i).applyMatrix4(m).toArray(ref, i * 3);
          const out = new Float32Array(a.count * 3);
          applyMatrixInto(attributeValues(a, 3), out, m);
          expect(same(out, ref), `${T.name} ${kind} × ${mn}: points`).toBe(true);
          k.getNormalMatrix(m);
          for (let i = 0; i < a.count; i++) v.fromBufferAttribute(a, i).applyMatrix3(k).normalize().toArray(ref, i * 3);
          applyNormalInto(attributeValues(a, 3), out, k);
          expect(same(out, ref), `${T.name} ${kind} × ${mn}: normals`).toBe(true);
        }
  });
});

