// Hero airframe loader (A6): the Smithsonian NASM 3-D scan
// (flyer-smithsonian.glb — draco glTF 2.0, CC0 1.0, provenance in
// data/wright-flyer/assets/flyer-model-provenance-v1.json) as an
// OPTIONAL presentation skin behind ?hero=1.
//
// ROLE-BOUNDARY LAW (flyer-model-provenance-v1): this mesh is
// PRESENTATION ONLY. Physics never reads it; no dimension from it
// feeds any model; the parametric airframe (FLYER_DIM, dossier-sourced)
// remains the articulation rig and always loads first. The scan is
// scaled to the dossier span so the two agree by construction.
//
// All assets are LOCAL (public/): no runtime fetches outside origin —
// CSP/no-fetch law respected.

import * as THREE from "three";

/** Dossier span [m] (FLYER_DIM.span mirror; flyer-reference.json). */
const DOSSIER_SPAN_M = 12.29;

export interface HeroModel {
  readonly group: THREE.Group;
  /** Dispose GPU resources for every loaded geometry/material. */
  dispose(): void;
}

let cached: Promise<HeroModel | null> | null = null;

/** Lazy-load + normalize the scan. Resolves NULL on any failure (the
 * app keeps the parametric airframe — a missing skin is silent, a
 * broken one would be a bug). Idempotent per page load. */
export function loadHeroAirframe(): Promise<HeroModel | null> {
  if (cached !== null) {
    return cached;
  }
  cached = (async (): Promise<HeroModel | null> => {
    try {
      const [{ GLTFLoader }, { DRACOLoader }] = await Promise.all([
        import("three/examples/jsm/loaders/GLTFLoader.js"),
        import("three/examples/jsm/loaders/DRACOLoader.js"),
      ]);
      const draco = new DRACOLoader();
      draco.setDecoderPath("/draco/");
      const loader = new GLTFLoader();
      loader.setDRACOLoader(draco);
      const gltf = await loader.loadAsync("/models/flyer-smithsonian.glb");
      const root = gltf.scene;
      // Normalize: rotate the span onto the parametric rig's wing axis
      // FIRST (nose is local +z there), THEN measure and center — the
      // centering offset must be computed in the FINAL orientation, or
      // a 90° remap swings the model off its mount by the difference
      // between the pre- and post-rotation bbox centers.
      root.updateMatrixWorld(true);
      const rawBox = new THREE.Box3().setFromObject(root);
      const rawSize = rawBox.getSize(new THREE.Vector3());
      if (![rawSize.x, rawSize.y, rawSize.z].every((v) => Number.isFinite(v) && v > 0)) {
        return null;
      }
      // The scan's span axis is X in Smithsonian convention; verify by
      // picking the LARGEST extent (robust to re-exports).
      const spanAxis =
        rawSize.x >= rawSize.y && rawSize.x >= rawSize.z
          ? "x"
          : rawSize.y >= rawSize.z
            ? "y"
            : "z";
      const scale = DOSSIER_SPAN_M / rawSize[spanAxis];
      root.scale.setScalar(scale);
      if (spanAxis === "x") {
        root.rotation.y = Math.PI / 2; // span x -> z
      } else if (spanAxis === "y") {
        root.rotation.z = Math.PI / 2; // span y -> z
      }
      root.updateMatrixWorld(true);
      const scaled = new THREE.Box3().setFromObject(root);
      const center = scaled.getCenter(new THREE.Vector3());
      root.position.sub(center);
      const group = new THREE.Group();
      group.add(root);
      group.visible = false; // caller reveals after alignment check
      let disposed = false;
      return {
        group,
        dispose(): void {
          if (disposed) {
            return;
          }
          disposed = true;
          // Drop the module cache TOO: a scene rebuild (?hero=1 again
          // after R/N) must re-load fresh GPU resources, never remount
          // these disposed geometries.
          cached = null;
          root.traverse((obj) => {
            const mesh = obj as THREE.Mesh;
            if (mesh.isMesh) {
              mesh.geometry?.dispose();
              const m = mesh.material;
              if (Array.isArray(m)) {
                m.forEach((mm) => mm.dispose());
              } else {
                m?.dispose();
              }
            }
          });
        },
      };
    } catch (err) {
      console.warn(
        JSON.stringify({
          suite: "wright-flyer-app",
          stage: "hero-model-unavailable",
          message: err instanceof Error ? err.message : String(err),
        }),
      );
      return null;
    }
  })();
  return cached;
}
