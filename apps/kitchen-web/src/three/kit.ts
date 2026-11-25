/**
 * Shared helpers for the page's 3D chapters: the renderer, scroll-safe orbit controls, the thermal colour ramp, the
 * soft particle sprite, canvas resizing, in-view pausing and the strip charts. Every model on the page is built in
 * code (src/three/models/); the kitchens' manifests come from ASSET_BASE.
 */

import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { getWispyCloudTexture } from "./textures";

export const ASSET_BASE = `${import.meta.env.BASE_URL}kitchen-assets`;

/** Create an ACESFilmic PBR WebGLRenderer with crisp PCFSoftShadowMap shadows. */
export function createStudioRenderer(canvas: HTMLCanvasElement, alpha = true): THREE.WebGLRenderer | null {
  try {
    const renderer = new THREE.WebGLRenderer({
      canvas,
      antialias: true,
      alpha,
      powerPreference: "high-performance",
    });
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.06;
    renderer.shadowMap.enabled = true;
    renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    return renderer;
  } catch {
    // Each chapter's model runs inside its 3D loop, so without WebGL its status line would wait forever: say why.
    canvas
      .closest("section")
      ?.querySelectorAll<HTMLElement>(".status")
      .forEach((s) => {
        s.dataset.state = "warn";
        s.textContent = "This model runs with its 3D view, which needs WebGL.";
      });
    const note = document.createElement("p");
    note.className = "webgl-fallback";
    note.textContent = "This 3D view needs WebGL, which this browser could not start.";
    canvas.replaceWith(note);
    return null;
  }
}

/**
 * Draw a scene only once its shaders are linked. The first draw of a new material otherwise waits on the GPU while
 * its shader compiles, freezing the page (two seconds for the hero's first kitchen). `refresh()` compiles whatever
 * the scene now holds in the background (KHR_parallel_shader_compile, where the browser has it) and `ready()` stays
 * false until that is done, so the frame loop skips its draw and the canvas keeps its last frame meanwhile. A compile
 * that has not finished in 4 s stops holding the draw, so a view can never stay blank.
 */
export interface ShaderGate {
  ready(): boolean;
  /** Compile what `scene` holds now. Resolves when linked or timed out. */
  refresh(scene: THREE.Scene): Promise<void>;
  /**
   * Run `fn` now or, while a compile is in flight, once it is done (a later call with the same key replaces it).
   * For calls that wait on the GPU, such as resizing the canvas, which would otherwise stall behind the compile.
   */
  defer(key: string, fn: () => void): void;
}
export function shaderGate(renderer: THREE.WebGLRenderer, camera: THREE.Camera): ShaderGate {
  let pending = 0;
  const queued = new Map<string, () => void>();
  return {
    ready: () => pending === 0,
    refresh(scene) {
      pending++;
      const compiled = renderer.compileAsync(scene, camera).catch(() => undefined);
      const timeout = new Promise((r) => setTimeout(r, 4000));
      return Promise.race([compiled, timeout]).then(() => {
        if (--pending > 0) return;
        const fns = [...queued.values()];
        queued.clear();
        for (const fn of fns) fn();
      });
    },
    defer(key, fn) {
      if (pending === 0) fn();
      else queued.set(key, fn);
    },
  };
}

/**
 * Scroll-safe OrbitControls: wheel zoom is OFF by default so normal vertical page
 * scrolling is NEVER hijacked; clicking/focusing the canvas unlocks zoom until blur/leave.
 */
export function createScrollSafeOrbit(
  camera: THREE.PerspectiveCamera,
  canvas: HTMLCanvasElement,
  target: readonly [number, number, number],
  opts?: { minDistance?: number; maxDistance?: number; maxPolarAngle?: number },
): OrbitControls {
  const controls = new OrbitControls(camera, canvas);
  controls.target.set(target[0], target[1], target[2]);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.enableZoom = false;
  controls.enablePan = false;
  controls.minDistance = opts?.minDistance ?? 0.5;
  controls.maxDistance = opts?.maxDistance ?? 35;
  controls.maxPolarAngle = opts?.maxPolarAngle ?? Math.PI * 0.49;

  const activateZoom = () => {
    controls.enableZoom = true;
  };
  const deactivateZoom = () => {
    controls.enableZoom = false;
  };
  canvas.addEventListener("pointerdown", activateZoom);
  canvas.addEventListener("focus", activateZoom);
  canvas.addEventListener("mouseleave", deactivateZoom);
  canvas.addEventListener("blur", deactivateZoom);

  controls.update();
  return controls;
}

/** Scientific ironbow / blackbody thermal false-color ramp in linear THREE.Color. */
export function thermalColor(tempC: number, minC = 0, maxC = 320, target = new THREE.Color()): THREE.Color {
  const t = Math.max(0, Math.min(1, (tempC - minC) / Math.max(1e-3, maxC - minC)));
  if (t < 0.25) {
    const u = t / 0.25;
    return target.setRGB(0.06 + u * 0.06, 0.14 + u * 0.24, 0.38 + u * 0.34);
  }
  if (t < 0.52) {
    const u = (t - 0.25) / 0.27;
    return target.setRGB(0.12 + u * 0.68, 0.38 - u * 0.22, 0.72 - u * 0.46);
  }
  if (t < 0.8) {
    const u = (t - 0.52) / 0.28;
    return target.setRGB(0.8 + u * 0.18, 0.16 + u * 0.48, 0.26 - u * 0.22);
  }
  const u = (t - 0.8) / 0.2;
  return target.setRGB(0.98 + u * 0.02, 0.64 + u * 0.32, 0.04 + u * 0.72);
}

/** Soft organic wispy cloud sprite texture for steam, smoke, and exhaust plumes. */
export function getSoftParticleTexture(): THREE.CanvasTexture {
  return getWispyCloudTexture();
}

export interface Scene3dFixture {
  readonly archetype_id: string;
  readonly fixture_id: string;
  readonly gltf_node?: string;
  readonly pose: { x_m: number; y_m: number; z_m: number; yaw_rad: number };
  readonly local_extents_m: readonly [number, number, number];
  readonly world_aabb_m?: { min_m: readonly [number, number, number]; max_m: readonly [number, number, number] };
  readonly swept_volumes?: readonly {
    volume_id: string;
    kind: string;
    aabb_world_m: { min_m: readonly [number, number, number]; max_m: readonly [number, number, number] };
    clearance_required_m: number;
  }[];
  readonly triangle_count?: number;
}

/** Sync a WebGLRenderer and PerspectiveCamera cleanly to container ResizeObserver measurements. */
export function bindCanvasResize(
  container: HTMLElement,
  renderer: THREE.WebGLRenderer,
  camera: THREE.PerspectiveCamera,
  gate?: ShaderGate,
): () => void {
  let lastW = 0;
  let lastH = 0;
  const apply = () => {
    const w = Math.max(1, Math.floor(container.clientWidth));
    const h = Math.max(1, Math.floor(container.clientHeight));
    if (Math.abs(w - lastW) <= 2 && Math.abs(h - lastH) <= 2) return;
    lastW = w;
    lastH = h;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  };
  // Resizing waits on the GPU: behind a shader compile it is held until the compile is done.
  const resize = () => (gate ? gate.defer("resize", apply) : apply());
  resize();
  const ro = new ResizeObserver(resize);
  ro.observe(container);
  return () => ro.disconnect();
}

/** Observe element viewport intersection so offscreen WebGL loops pause immediately. */
export function observeInView(el: Element, cb: (visible: boolean) => void, rootMargin = "40px 0px"): () => void {
  const io = new IntersectionObserver(
    (entries) => {
      const entry = entries[0];
      if (entry) cb(entry.isIntersecting);
    },
    { rootMargin },
  );
  io.observe(el);
  return () => io.disconnect();
}

export interface StripSeries {
  readonly label: string;
  readonly color: string;
  readonly values: readonly number[];
  readonly dashed?: boolean;
}

/** Crisp Retina 2D strip-chart renderer with subtle gradient area fill for live physics readouts. */
export function drawStripChart(
  canvas: HTMLCanvasElement,
  series: readonly StripSeries[],
  opts: { minY: number; maxY: number; unit: string; threshold?: { value: number; color: string; label: string } },
): void {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(220, Math.floor(canvas.clientWidth || 320));
  const h = Math.max(96, Math.floor(canvas.clientHeight || 116));
  if (canvas.width !== Math.floor(w * dpr) || canvas.height !== Math.floor(h * dpr)) {
    canvas.width = Math.floor(w * dpr);
    canvas.height = Math.floor(h * dpr);
  }
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  ctx.save();
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, w, h);

  const padL = 38;
  const padR = 12;
  const padT = 10;
  const padB = 20;
  const pw = w - padL - padR;
  const ph = h - padT - padB;

  // Gridlines
  ctx.strokeStyle = "rgba(255, 224, 196, 0.09)";
  ctx.lineWidth = 1;
  ctx.font = '500 10px "JetBrains Mono Variable", "JetBrains Mono", monospace';
  ctx.fillStyle = "rgba(214, 198, 184, 0.62)";
  for (let i = 0; i <= 3; i += 1) {
    const f = i / 3;
    const y = padT + ph * (1 - f);
    const val = opts.minY + (opts.maxY - opts.minY) * f;
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + pw, y);
    ctx.stroke();
    ctx.fillText(`${Math.round(val)}`, 4, y + 3);
  }

  // Optional threshold line (e.g. smoke point or HACCP 41 F)
  if (opts.threshold) {
    const f = Math.max(0, Math.min(1, (opts.threshold.value - opts.minY) / Math.max(1e-3, opts.maxY - opts.minY)));
    const y = padT + ph * (1 - f);
    ctx.save();
    ctx.setLineDash([4, 4]);
    ctx.strokeStyle = opts.threshold.color;
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + pw, y);
    ctx.stroke();
    ctx.restore();
    ctx.fillStyle = opts.threshold.color;
    ctx.fillText(opts.threshold.label, padL + 4, Math.max(padT + 9, y - 4));
  }

  for (let si = 0; si < series.length; si += 1) {
    const s = series[si];
    if (s.values.length < 2) continue;
    ctx.save();
    if (s.dashed) ctx.setLineDash([5, 3]);
    ctx.strokeStyle = s.color;
    ctx.lineWidth = 2.1;
    ctx.beginPath();
    const n = s.values.length;
    let lastX = padL;
    let lastY = padT + ph;
    for (let i = 0; i < n; i += 1) {
      const x = padL + (i / (n - 1)) * pw;
      const f = Math.max(0, Math.min(1, (s.values[i] - opts.minY) / Math.max(1e-3, opts.maxY - opts.minY)));
      const y = padT + ph * (1 - f);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
      lastX = x;
      lastY = y;
    }
    ctx.stroke();

    // Subtle luminous area fill under primary series
    if (si === 0 && !s.dashed) {
      ctx.lineTo(padL + pw, padT + ph);
      ctx.lineTo(padL, padT + ph);
      ctx.closePath();
      const grad = ctx.createLinearGradient(0, padT, 0, padT + ph);
      grad.addColorStop(0, `${s.color}2e`);
      grad.addColorStop(1, `${s.color}00`);
      ctx.fillStyle = grad;
      ctx.fill();
    }

    // Crisp live endpoint dot
    ctx.setLineDash([]);
    ctx.fillStyle = s.color;
    ctx.beginPath();
    ctx.arc(lastX, lastY, 3, 0, Math.PI * 2);
    ctx.fill();
    ctx.restore();
  }

  // Inline legend at bottom
  let lx = padL;
  for (const s of series) {
    ctx.fillStyle = s.color;
    ctx.fillRect(lx, h - 10, 8, 3);
    ctx.fillStyle = "rgba(247, 241, 235, 0.84)";
    ctx.fillText(s.label, lx + 12, h - 5);
    lx += ctx.measureText(s.label).width + 24;
  }
  ctx.fillStyle = "rgba(214, 198, 184, 0.55)";
  ctx.fillText(opts.unit, w - padR - ctx.measureText(opts.unit).width, h - 5);
  ctx.restore();
}
