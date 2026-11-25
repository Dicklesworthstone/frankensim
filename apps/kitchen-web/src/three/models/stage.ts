/**
 * Lighting for every scene on the page: a studio environment for the stainless to reflect (soft overhead banks, a
 * warm key, a cool rim), one shadow-casting key light, and filmic tone mapping. Metals read as metal only when they
 * have something to reflect, so the environment carries most of the look.
 */
import * as THREE from "three";
import { RoomEnvironment } from "three/examples/jsm/environments/RoomEnvironment.js";
import { glOf, viewGl } from "../kit";

export function makeRenderer(canvas: HTMLCanvasElement, opts: { alpha?: boolean; shadows?: boolean } = {}): THREE.WebGLRenderer | null {
  return viewGl(canvas, { alpha: opts.alpha ?? true, exposure: 1.0, shadows: opts.shadows ?? true });
}

const envCache = new WeakMap<THREE.WebGLRenderer, THREE.Texture>();

/** A studio room to reflect: RoomEnvironment's soft boxes, warmed on one side and cooled on the other. */
export function studioEnvironment(view: THREE.WebGLRenderer): THREE.Texture {
  // Built once by the renderer that holds the GPU state, for every view it draws.
  const r = glOf(view);
  const hit = envCache.get(r);
  if (hit) return hit;
  const room = new RoomEnvironment();
  // Tint the room's lamps: warm key side, cool fill side, so steel shows both.
  let k = 0;
  room.traverse((o) => {
    const m = (o as THREE.Mesh).material as THREE.MeshBasicMaterial | THREE.MeshLambertMaterial | undefined;
    if (!m || !("color" in m)) return;
    if ((m as THREE.MeshBasicMaterial).isMeshBasicMaterial) {
      const c = m.color;
      if (k % 3 === 0) c.multiply(new THREE.Color(1.0, 0.9, 0.78));
      else if (k % 3 === 1) c.multiply(new THREE.Color(0.82, 0.9, 1.0));
      k++;
    }
  });
  const pm = new THREE.PMREMGenerator(r);
  const t = pm.fromScene(room, 0.03).texture;
  pm.dispose();
  envCache.set(r, t);
  return t;
}

export interface Lights {
  key: THREE.DirectionalLight;
  hemi: THREE.HemisphereLight;
  rim: THREE.DirectionalLight;
}

/**
 * The scene's lights: a warm key from high front right casting soft shadows over `span` metres round `centre`, a sky
 * and floor bounce, and a cool rim from behind left.
 */
export function studioLights(scene: THREE.Scene, centre: THREE.Vector3, span: number, shadowSize = 2048): Lights {
  const hemi = new THREE.HemisphereLight(0xfff3e6, 0x3a3028, 0.55);
  scene.add(hemi);
  const key = new THREE.DirectionalLight(0xfff0de, 1.9);
  key.position.copy(centre).add(new THREE.Vector3(span * 0.55, span * 1.1, span * 0.75));
  key.target.position.copy(centre);
  key.castShadow = true;
  key.shadow.mapSize.set(shadowSize, shadowSize);
  const sc = key.shadow.camera;
  sc.left = -span;
  sc.right = span;
  sc.top = span;
  sc.bottom = -span;
  sc.near = 0.05;
  sc.far = span * 4;
  key.shadow.bias = -0.0002;
  key.shadow.normalBias = 0.012 * Math.max(0.2, span / 4);
  key.shadow.radius = 3;
  scene.add(key, key.target);
  const rim = new THREE.DirectionalLight(0xa9c8ff, 0.7);
  rim.position.copy(centre).add(new THREE.Vector3(-span * 0.9, span * 0.6, -span * 0.8));
  scene.add(rim);
  return { key, hemi, rim };
}

/** Fit a renderer and camera to their element's size, now and whenever it changes. */
export function fitToElement(el: HTMLElement, r: THREE.WebGLRenderer, cam: THREE.PerspectiveCamera): () => void {
  const fit = () => {
    const w = Math.max(1, el.clientWidth);
    const h = Math.max(1, el.clientHeight);
    r.setSize(w, h, false);
    cam.aspect = w / h;
    cam.updateProjectionMatrix();
  };
  fit();
  const ro = new ResizeObserver(fit);
  ro.observe(el);
  return () => ro.disconnect();
}
