/**
 * Shared helpers for the page's 3D chapters: the renderer, scroll-safe orbit controls, the thermal colour ramp, the
 * soft particle sprite, canvas resizing, in-view pausing and the strip charts. Every model on the page is built in
 * code (src/three/models/); the kitchens' manifests come from ASSET_BASE.
 */

import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { getWispyCloudTexture } from "./textures";

export const ASSET_BASE = `${import.meta.env.BASE_URL}kitchen-assets`;

/* ------------------------------------------------------------------------------------------------------------ */
/* One WebGL context for the page                                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Every 3D view on the page draws with one WebGLRenderer. With a context each, the nine views linked their own copies
 * of the shaders they share (355 programs for 220 distinct), uploaded the same textures again (78 of 134 Mpx), built
 * the studio environment nine times and waited on the GPU to set up each context. A view keeps its own canvas, a 2D
 * one: its frame is drawn in the lower-left corner of the shared drawing buffer (viewport and scissor at its size)
 * and copied across unchanged. Its settings (tone mapping, exposure, colour space, shadows, clipping, clear colour,
 * size, pixel ratio) stay its own and are set on the shared renderer before each of its draws.
 *
 * `?gl=each` gives every view its own context again, for comparison.
 */
let shared: THREE.WebGLRenderer | null | undefined;
let bufW = 0;
let bufH = 0;
const ownGl = new WeakMap<object, THREE.WebGLRenderer>();
const EACH = typeof location !== "undefined" && new URLSearchParams(location.search).get("gl") === "each";
/**
 * three reads each shader program's logs and link status back the first time it draws with it, to report a failed
 * compile: four calls per program, each waiting on the GPU process. The built page skips them (as three advises for
 * production); development, and `?shaders=check` (scripts/verify-chrome.mjs), keep them.
 */
const CHECK_SHADERS = import.meta.env.DEV || (typeof location !== "undefined" && new URLSearchParams(location.search).get("shaders") === "check");

function sharedRenderer(): THREE.WebGLRenderer | null {
  if (shared !== undefined) return shared;
  try {
    shared = new THREE.WebGLRenderer({ canvas: document.createElement("canvas"), antialias: true, alpha: true, powerPreference: "high-performance" });
    shared.debug.checkShaderErrors = CHECK_SHADERS;
    shared.setPixelRatio(1);
  } catch {
    shared = null;
  }
  return shared;
}

/** The renderer holding a view's GPU state: for passes of its own (an environment map) that must not be copied out. */
export function glOf(r: THREE.WebGLRenderer): THREE.WebGLRenderer {
  return ownGl.get(r) ?? r;
}

/**
 * A view on `canvas` drawn by the page's shared renderer, standing in for a WebGLRenderer of its own: `render`,
 * sizing, clear colour and the settings above are the view's; anything else is the shared renderer's.
 */
function viewRenderer(canvas: HTMLCanvasElement, alpha: boolean, exposure: number, shadows: boolean): THREE.WebGLRenderer | null {
  const gl = sharedRenderer();
  const ctx = gl && canvas.getContext("2d", { alpha });
  if (!gl || !ctx) return null;
  let cssW = canvas.clientWidth || canvas.width;
  let cssH = canvas.clientHeight || canvas.height;
  let ratio = 1;
  const clear = new THREE.Color(0x000000);
  let clearA = alpha ? 0 : 1;
  const counts = { calls: 0, triangles: 0, points: 0, lines: 0, frame: 0 };
  const view: Record<string | symbol, unknown> = {
    domElement: canvas,
    toneMapping: THREE.ACESFilmicToneMapping,
    toneMappingExposure: exposure,
    outputColorSpace: THREE.SRGBColorSpace,
    localClippingEnabled: false,
    shadowMap: { enabled: shadows, type: THREE.PCFSoftShadowMap },
    info: {
      render: counts,
      get memory() {
        return gl.info.memory;
      },
      get programs() {
        return gl.info.programs;
      },
    },
    setPixelRatio(v: number) {
      ratio = v;
      size(cssW, cssH, false);
    },
    getPixelRatio: () => ratio,
    setSize: (w: number, h: number, style = true) => size(w, h, style),
    getSize: (t: THREE.Vector2) => t.set(cssW, cssH),
    getDrawingBufferSize: (t: THREE.Vector2) => t.set(canvas.width, canvas.height),
    setClearColor(c: THREE.ColorRepresentation, a = 1) {
      clear.set(c);
      clearA = a;
    },
    getClearColor: (t: THREE.Color) => t.copy(clear),
    getClearAlpha: () => clearA,
    setClearAlpha(a: number) {
      clearA = a;
    },
    compile(scene: THREE.Object3D, camera: THREE.Camera, target?: THREE.Scene | null) {
      apply();
      return gl.compile(scene, camera, target);
    },
    compileAsync(scene: THREE.Object3D, camera: THREE.Camera, target?: THREE.Scene | null) {
      apply();
      return gl.compileAsync(scene, camera, target);
    },
    render(scene: THREE.Object3D, camera: THREE.Camera) {
      const w = canvas.width;
      const h = canvas.height;
      if (!w || !h) return;
      // The shared buffer only grows, to the largest view drawn: resizing it on every switch would reallocate it.
      if (w > bufW || h > bufH) {
        bufW = Math.max(bufW, w);
        bufH = Math.max(bufH, h);
        gl.setSize(bufW, bufH, false);
      }
      apply();
      gl.setViewport(0, 0, w, h);
      gl.setScissor(0, 0, w, h);
      gl.setScissorTest(true);
      gl.render(scene, camera);
      const r = gl.info.render;
      counts.calls = r.calls;
      counts.triangles = r.triangles;
      counts.points = r.points;
      counts.lines = r.lines;
      counts.frame = r.frame;
      // "copy": the view's pixels become the frame's exactly (no blending over the last one), alpha and all.
      ctx.globalCompositeOperation = "copy";
      ctx.drawImage(gl.domElement, 0, bufH - h, w, h, 0, 0, w, h);
      // For checks by script: the view has drawn a frame.
      if (canvas.dataset.drawn !== "1") canvas.dataset.drawn = "1";
    },
    dispose() {},
  };
  const shadow = view.shadowMap as { enabled: boolean; type: THREE.ShadowMapType };
  const apply = () => {
    gl.toneMapping = view.toneMapping as THREE.ToneMapping;
    gl.toneMappingExposure = view.toneMappingExposure as number;
    gl.outputColorSpace = view.outputColorSpace as THREE.ColorSpace;
    gl.localClippingEnabled = view.localClippingEnabled as boolean;
    gl.shadowMap.enabled = shadow.enabled;
    gl.shadowMap.type = shadow.type;
    gl.setClearColor(clear, clearA);
  };
  function size(w: number, h: number, style: boolean): void {
    cssW = w;
    cssH = h;
    canvas.width = Math.floor(w * ratio);
    canvas.height = Math.floor(h * ratio);
    if (style) {
      canvas.style.width = `${w}px`;
      canvas.style.height = `${h}px`;
    }
  }
  const proxy = new Proxy(view, {
    get(t, k) {
      if (k in t) return t[k];
      const v = (gl as unknown as Record<string | symbol, unknown>)[k];
      return typeof v === "function" ? v.bind(gl) : v;
    },
    set(t, k, v) {
      if (k in t) t[k] = v;
      else (gl as unknown as Record<string | symbol, unknown>)[k] = v;
      return true;
    },
  }) as unknown as THREE.WebGLRenderer;
  ownGl.set(proxy, gl);
  canvas.dataset.view = "webgl";
  gl.domElement.addEventListener("webglcontextlost", () => (canvas.dataset.view = "lost"));
  return proxy;
}

/**
 * A view's renderer: the page's shared one (above), or with `?gl=each` a WebGLRenderer of its own. ACES filmic tone
 * mapping at `exposure`, sRGB output, PCF soft shadows, pixel ratio up to 2.
 */
export function viewGl(canvas: HTMLCanvasElement, opts: { alpha?: boolean; exposure?: number; shadows?: boolean } = {}): THREE.WebGLRenderer | null {
  const alpha = opts.alpha ?? true;
  const exposure = opts.exposure ?? 1.06;
  const shadows = opts.shadows ?? true;
  let r: THREE.WebGLRenderer | null = null;
  if (!EACH) r = viewRenderer(canvas, alpha, exposure, shadows);
  else
    try {
      r = new THREE.WebGLRenderer({ canvas, antialias: true, alpha, powerPreference: "high-performance" });
      r.debug.checkShaderErrors = CHECK_SHADERS;
      r.outputColorSpace = THREE.SRGBColorSpace;
      r.toneMapping = THREE.ACESFilmicToneMapping;
      r.toneMappingExposure = exposure;
      r.shadowMap.enabled = shadows;
      r.shadowMap.type = THREE.PCFSoftShadowMap;
      canvas.dataset.view = "webgl";
    } catch {
      r = null;
    }
  r?.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  return r;
}

/** Create an ACESFilmic PBR renderer for a view, with crisp PCFSoftShadowMap shadows (see viewGl). */
export function createStudioRenderer(canvas: HTMLCanvasElement, alpha = true): THREE.WebGLRenderer | null {
  const renderer = viewGl(canvas, { alpha });
  if (renderer) return renderer;
  {
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
/**
 * three's compileAsync, except that a material disposed while it compiles stops being waited for. three's own polls
 * each material's program and throws once the material has been disposed (a piece swapped out before its shaders were
 * ready, as picking pieces quickly in the equipment chapter does): its promise then never settled, so the gate counted
 * it pending until its 4 s timeout, and the error went uncaught. After `giveUpMs` it stops asking: a lost context
 * never reports a program ready, and three's would poll it every 10 ms for good.
 */
function compileSettled(renderer: THREE.WebGLRenderer, scene: THREE.Object3D, camera: THREE.Camera, giveUpMs: number): Promise<void> {
  return new Promise((resolve) => {
    const materials = renderer.compile(scene, camera);
    const t0 = performance.now();
    const check = () => {
      for (const m of materials) {
        const program = (renderer.properties.get(m) as { currentProgram?: { isReady(): boolean } }).currentProgram;
        if (!program || program.isReady()) materials.delete(m);
      }
      if (materials.size === 0 || performance.now() - t0 > giveUpMs) resolve();
      else setTimeout(check, 10);
    };
    // Without the extension three counts every program ready at once, and its first draw waits for the link: give the
    // driver a moment's start first, as three's compileAsync does.
    if (renderer.extensions.has("KHR_parallel_shader_compile")) check();
    else setTimeout(check, 10);
  });
}

/** How long a compile may hold a view's draw (ShaderGate). */
const GATE_MS = 4000;

export function shaderGate(renderer: THREE.WebGLRenderer, camera: THREE.Camera): ShaderGate {
  let pending = 0;
  const queued = new Map<string, () => void>();
  return {
    ready: () => pending === 0,
    refresh(scene) {
      pending++;
      const compiled = compileSettled(renderer, scene, camera, GATE_MS).catch(() => undefined);
      const timeout = new Promise((r) => setTimeout(r, GATE_MS));
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
 * Clip planes that follow the camera, for a scene whose drawn parts all lie inside `box`: the near plane at 0.6 of the
 * distance to the box's nearest point (a plane, it cuts nothing nearer than that within 53° of the view axis, wider
 * than any of the page's frusta reach), never under `minNear`, and the far plane just past the box's farthest corner.
 * A fixed near plane of a few centimetres leaves the depth buffer able to part only about z² / (near · 2²⁴): half a
 * millimetre at 20 m, so a pan on a counter or a roof panel on its frame flickers from a high, distant view. Fitted to
 * the box, the same view parts surfaces a few hundredths of a millimetre apart.
 */
/**
 * Clip planes fitted each frame to the box round everything a scene draws (hidden parts included, a margin and a
 * standing person's headroom added for what moves about in it), the box measured again now and then for parts that
 * arrive late.
 */
export function clipFitter(camera: THREE.PerspectiveCamera, scene: THREE.Object3D, minNear: number, headroom = 2.4): () => void {
  const box = new THREE.Box3();
  let age = Infinity;
  return () => {
    if (++age > 90) {
      box.setFromObject(scene).expandByScalar(0.25);
      box.max.y = Math.max(box.max.y, headroom);
      age = 0;
    }
    fitClipPlanes(camera, box, minNear);
  };
}

export function fitClipPlanes(camera: THREE.PerspectiveCamera, box: THREE.Box3, minNear = 0.05): void {
  const p = camera.position;
  const dNear = box.distanceToPoint(p);
  let dFar = 0;
  for (let i = 0; i < 8; i++) {
    const x = i & 1 ? box.max.x : box.min.x;
    const y = i & 2 ? box.max.y : box.min.y;
    const z = i & 4 ? box.max.z : box.min.z;
    dFar = Math.max(dFar, Math.hypot(x - p.x, y - p.y, z - p.z));
  }
  const near = Math.max(minNear, dNear * 0.6);
  const far = dFar + 1;
  if (Math.abs(camera.near - near) > near * 0.05 || Math.abs(camera.far - far) > far * 0.05) {
    camera.near = near;
    camera.far = far;
    camera.updateProjectionMatrix();
  }
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
