// The live soundboard in 3D: the engine's own surface (its quadrature points and the 32 board
// coordinates' mode shapes, interpolated onto a grid) displaced by the modal sum w(x, y) of each
// snapshot, coloured by its sign, with the microphone where the engine listens. Heights are
// exaggerated to a fixed fraction of the board; the scale is stated on the panel.

import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { reducedMotion } from "../util";

const UP: [number, number, number] = [70, 222, 198];
const DOWN: [number, number, number] = [240, 92, 44];
const SPRUCE: [number, number, number] = [176, 136, 84];

export class SoundboardView {
  private renderer: THREE.WebGLRenderer | null = null;
  private scene = new THREE.Scene();
  private camera = new THREE.PerspectiveCamera(32, 1.4, 0.05, 40);
  private controls: OrbitControls | null = null;
  private geom: THREE.BufferGeometry | null = null;
  private mesh: THREE.Mesh | null = null;
  private vshapes = new Float32Array(0); // [vertex][mode]
  private modes = 0;
  private nv = 0;
  private inside = new Uint8Array(0);
  private w = new Float32Array(0);
  private grain = new Float32Array(0);
  private peak = 0;
  private spanX = 1;
  private centre = new THREE.Vector3();
  private mic: THREE.Group | null = null;
  private rings: THREE.Mesh[] = [];
  private lastInteract = -1e9;
  private dirty = true;
  /** Peak |w| at the last update [m], for the scale label. */
  peakM = 0;
  /** Visual height of the largest displacement, as a fraction of the board's width. */
  readonly lift = 0.07;
  ok = false;

  constructor(private readonly host: HTMLElement) {
    try {
      this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "high-performance" });
    } catch {
      host.insertAdjacentHTML("beforeend", '<div class="stage-note">WebGL is unavailable here, so the soundboard cannot be drawn.</div>');
      return;
    }
    const r = this.renderer;
    r.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    r.outputColorSpace = THREE.SRGBColorSpace;
    r.toneMapping = THREE.ACESFilmicToneMapping;
    r.toneMappingExposure = 1.1;
    r.domElement.className = "pp-gl";
    host.prepend(r.domElement);
    this.camera.up.set(0, 0, 1);
    this.scene.add(new THREE.HemisphereLight(0xfff1dc, 0x16110c, 0.9));
    const key = new THREE.DirectionalLight(0xffe2b8, 2.2);
    key.position.set(1.5, -1.6, 3.4);
    const rim = new THREE.DirectionalLight(0x8fd6ca, 0.9);
    rim.position.set(-2.4, 2.8, 1.4);
    this.scene.add(key, rim);
    this.controls = new OrbitControls(this.camera, r.domElement);
    this.controls.enablePan = false;
    this.controls.minDistance = 1.4;
    this.controls.maxDistance = 6;
    this.controls.addEventListener("start", () => { this.lastInteract = performance.now(); });
    // a turn by hand draws at once, playing or not
    this.controls.addEventListener("change", () => { this.lastInteract = performance.now(); this.renderer?.render(this.scene, this.camera); });
    // a finger scrolls the page unless it starts a deliberate rotation on the board
    r.domElement.style.touchAction = "pan-y";
    new ResizeObserver(() => this.resize()).observe(host);
    this.resize();
    this.ok = true;
  }

  /** The engine's board: quadrature points [m] and mode shapes [point][mode]. */
  setBoard(xy: Float32Array, shapes: Float32Array): void {
    if (!this.renderer) return;
    const points = xy.length / 2;
    const modes = Math.round(shapes.length / points);
    this.modes = modes;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (let i = 0; i < points; i++) { const x = xy[2 * i]!, y = xy[2 * i + 1]!; minX = Math.min(minX, x); maxX = Math.max(maxX, x); minY = Math.min(minY, y); maxY = Math.max(maxY, y); }
    const spacing = Math.sqrt(((maxX - minX) * (maxY - minY)) / points);
    this.spanX = maxX - minX;
    this.centre.set((minX + maxX) / 2, (minY + maxY) / 2, 0);
    // grid over the board, each vertex an inverse-distance blend of its nearest quadrature points
    const step = Math.max(spacing * 0.6, (maxX - minX) / 110);
    const nx = Math.ceil((maxX - minX) / step) + 1, ny = Math.ceil((maxY - minY) / step) + 1;
    const cell = spacing * 2;
    const gx = Math.ceil((maxX - minX) / cell) + 1, gy = Math.ceil((maxY - minY) / cell) + 1;
    const buckets: number[][] = Array.from({ length: gx * gy }, () => []);
    for (let i = 0; i < points; i++) buckets[Math.floor((xy[2 * i + 1]! - minY) / cell) * gx + Math.floor((xy[2 * i]! - minX) / cell)]!.push(i);
    const nv = nx * ny;
    this.nv = nv;
    this.vshapes = new Float32Array(nv * modes);
    this.inside = new Uint8Array(nv);
    const pos = new Float32Array(nv * 3);
    const limit = (spacing * 1.25) ** 2;
    // inside = not reachable from the grid's border through points far from every quadrature
    // point, so gaps between sparse points stay part of the board (pianoStage's raster law)
    const nearest: [number, number][][] = new Array(nv);
    const near = new Uint8Array(nv);
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const v = j * nx + i, x = minX + i * step, y = minY + j * step;
      const cx = Math.floor((x - minX) / cell), cy = Math.floor((y - minY) / cell);
      const list: [number, number][] = [];
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        if (cx + dx < 0 || cx + dx >= gx || cy + dy < 0 || cy + dy >= gy) continue;
        for (const p of buckets[(cy + dy) * gx + cx + dx]!) list.push([(xy[2 * p]! - x) ** 2 + (xy[2 * p + 1]! - y) ** 2, p]);
      }
      list.sort((a, b) => a[0] - b[0]);
      nearest[v] = list.slice(0, 4);
      if (list.length && list[0]![0] <= limit) near[v] = 1;
    }
    const outside = new Uint8Array(nv);
    const stack: number[] = [];
    const seed = (v: number): void => { if (!near[v] && !outside[v]) { outside[v] = 1; stack.push(v); } };
    for (let i = 0; i < nx; i++) { seed(i); seed((ny - 1) * nx + i); }
    for (let j = 0; j < ny; j++) { seed(j * nx); seed(j * nx + nx - 1); }
    while (stack.length) {
      const v = stack.pop()!;
      const i = v % nx, j = (v - i) / nx;
      if (i > 0) seed(v - 1);
      if (i < nx - 1) seed(v + 1);
      if (j > 0) seed(v - nx);
      if (j < ny - 1) seed(v + nx);
    }
    const grain = new Float32Array(nv);
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const v = j * nx + i, x = minX + i * step, y = minY + j * step;
      pos[3 * v] = x; pos[3 * v + 1] = y;
      // spruce grain runs along the board's long axis: a few soft stripes, fixed
      grain[v] = 0.93 + 0.035 * Math.sin(x * 47 + Math.sin(y * 2.3) * 1.6) + 0.025 * Math.sin(x * 131 + y * 0.9) + 0.012 * Math.sin(x * 317);
      const use = nearest[v]!;
      if (outside[v] || !use.length) continue;
      this.inside[v] = 1;
      let wsum = 0;
      for (const [d] of use) wsum += 1 / (d + 1e-6);
      for (const [d, p] of use) {
        const wt = 1 / (d + 1e-6) / wsum;
        for (let k = 0; k < modes; k++) this.vshapes[v * modes + k]! += wt * shapes[p * modes + k]!;
      }
    }
    // round off the grid's staircase outline: each vertex on the edge moves to the mean of the
    // board vertices around it (twice); interior vertices and every field value stay put
    for (let pass = 0; pass < 2; pass++) {
      const moved = new Float32Array(pos);
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
        const v = j * nx + i;
        if (!this.inside[v]) continue;
        const edge = i === 0 || j === 0 || i === nx - 1 || j === ny - 1 || !this.inside[v - 1] || !this.inside[v + 1] || !this.inside[v - nx] || !this.inside[v + nx];
        if (!edge) continue;
        let sx = 0, sy = 0, n = 0;
        for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) {
          const u = (j + dj) * nx + (i + di);
          if (i + di < 0 || i + di >= nx || j + dj < 0 || j + dj >= ny || !this.inside[u]) continue;
          sx += pos[3 * u]!; sy += pos[3 * u + 1]!; n++;
        }
        moved[3 * v] = sx / n; moved[3 * v + 1] = sy / n;
      }
      pos.set(moved);
    }
    const idx: number[] = [];
    for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
      const a = j * nx + i, b = a + 1, c = a + nx, d = c + 1;
      if (this.inside[a] && this.inside[b] && this.inside[d]) idx.push(a, b, d);
      if (this.inside[a] && this.inside[d] && this.inside[c]) idx.push(a, d, c);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    g.setAttribute("color", new THREE.BufferAttribute(new Float32Array(nv * 3), 3));
    g.setIndex(idx);
    this.geom = g;
    this.grain = grain;
    this.w = new Float32Array(nv);
    this.mesh = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.55, metalness: 0, side: THREE.DoubleSide }));
    this.scene.add(this.mesh);
    // a soft shadow under the board
    const c = document.createElement("canvas");
    c.width = c.height = 128;
    const cg = c.getContext("2d")!;
    const rg = cg.createRadialGradient(64, 64, 4, 64, 64, 64);
    rg.addColorStop(0, "rgba(0,0,0,0.55)"); rg.addColorStop(1, "rgba(0,0,0,0)");
    cg.fillStyle = rg; cg.fillRect(0, 0, 128, 128);
    const shadow = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), new THREE.MeshBasicMaterial({ map: new THREE.CanvasTexture(c), transparent: true, depthWrite: false }));
    shadow.scale.set((maxX - minX) * 1.5, (maxY - minY) * 1.3, 1);
    shadow.position.set(this.centre.x, this.centre.y, -0.08);
    this.scene.add(shadow);
    this.buildMic();
    this.paint(new Float32Array(modes));
    this.place();
    this.controls?.target.copy(this.centre);
    this.controls?.update();
    this.renderer.render(this.scene, this.camera);
  }

  private buildMic(): void {
    // the engine listens at (0.675, 1.0) m, 1 m above the board (its info: microphoneM)
    const g = new THREE.Group();
    const body = new THREE.Mesh(new THREE.CylinderGeometry(0.018, 0.014, 0.12, 20), new THREE.MeshStandardMaterial({ color: 0x2b2e34, metalness: 0.7, roughness: 0.35 }));
    body.rotation.x = Math.PI / 2;
    const grille = new THREE.Mesh(new THREE.SphereGeometry(0.024, 20, 14), new THREE.MeshStandardMaterial({ color: 0xc9cdd3, metalness: 0.9, roughness: 0.3 }));
    grille.position.z = -0.07;
    g.add(body, grille);
    const lead = new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, 0, 0.06), new THREE.Vector3(0, 0, 0.9)]), new THREE.LineBasicMaterial({ color: 0x5a5d63 }));
    g.add(lead);
    for (let i = 0; i < 3; i++) {
      const ring = new THREE.Mesh(new THREE.RingGeometry(0.03, 0.034, 48), new THREE.MeshBasicMaterial({ color: 0x9ff0e2, transparent: true, opacity: 0, side: THREE.DoubleSide, depthWrite: false }));
      ring.position.z = -0.075;
      g.add(ring);
      this.rings.push(ring);
    }
    g.position.set(0.675, 1.0, 1.0);
    this.scene.add(g);
    this.mic = g;
  }

  /** Board coordinates from a snapshot: w = Σ shapes·q, painted and lifted. */
  paint(coords: Float32Array): void {
    if (!this.geom) return;
    const m = this.modes, w = this.w;
    let peak = 0;
    for (let v = 0; v < this.nv; v++) {
      if (!this.inside[v]) continue;
      let s = 0;
      const o = v * m;
      for (let k = 0; k < m; k++) s += this.vshapes[o + k]! * (coords[k] ?? 0);
      w[v] = s;
      peak = Math.max(peak, Math.abs(s));
    }
    // the colour and height scale follow the motion down slowly, so a decaying board stays visible
    if (peak > this.peak || reducedMotion()) this.peak = peak;
    else this.peak = Math.max(peak, this.peak * 0.97);
    this.peakM = this.peak;
    const scale = this.peak > 1e-15 ? this.peak : 1;
    const pos = this.geom.getAttribute("position") as THREE.BufferAttribute;
    const col = this.geom.getAttribute("color") as THREE.BufferAttribute;
    const P = pos.array as Float32Array, C = col.array as Float32Array;
    const lift = this.lift * this.spanX;
    for (let v = 0; v < this.nv; v++) {
      const t = Math.max(-1, Math.min(1, w[v]! / scale));
      P[3 * v + 2] = t * lift;
      const a = Math.min(1, Math.abs(t) * 1.6) ** 0.7;
      const hue = t >= 0 ? UP : DOWN;
      const gr = this.grain[v] ?? 1;
      C[3 * v] = srgb((SPRUCE[0] + (hue[0] - SPRUCE[0]) * a) * gr);
      C[3 * v + 1] = srgb((SPRUCE[1] + (hue[1] - SPRUCE[1]) * a) * gr);
      C[3 * v + 2] = srgb((SPRUCE[2] + (hue[2] - SPRUCE[2]) * a) * gr);
    }
    pos.needsUpdate = true;
    col.needsUpdate = true;
    this.geom.computeVertexNormals();
    this.dirty = true;
  }

  /** Pulse the microphone's rings with the energy flowing into the board (0..1). */
  pulse(level: number, wall: number): void {
    this.rings.forEach((ring, i) => {
      const phase = ((wall / 900 + i / 3) % 1);
      const s = 1 + phase * 5;
      ring.scale.set(s, s, 1);
      (ring.material as THREE.MeshBasicMaterial).opacity = Math.max(0, Math.min(1, level)) * (1 - phase) * 0.7;
    });
    if (level > 0.001) this.dirty = true;
  }

  private resize(): void {
    if (!this.renderer) return;
    const w = this.host.clientWidth || 600;
    const h = this.host.clientHeight || 400;
    this.renderer.setSize(w, h, false);
    this.renderer.domElement.style.width = "100%";
    this.renderer.domElement.style.height = "100%";
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    if (this.geom) this.place();
    this.renderer.render(this.scene, this.camera);
  }
  /** The opening view, stepped back in a narrow frame so the whole board fits. */
  private place(): void {
    const d = 3.15 * Math.max(1, 1.35 / this.camera.aspect);
    this.camera.position.set(this.centre.x + d * 0.48, this.centre.y - d * 0.74, d * 0.74);
    this.controls?.target.copy(this.centre);
    this.controls?.update();
  }

  /** Draw if anything changed; turn slowly when nobody is steering. */
  frame(wall: number, dt: number, spin: boolean): void {
    if (!this.renderer || !this.controls) return;
    if (spin && !reducedMotion() && wall - this.lastInteract > 5000 && this.geom) {
      const t = this.controls.target;
      const p = this.camera.position.clone().sub(t);
      p.applyAxisAngle(new THREE.Vector3(0, 0, 1), dt * 0.05);
      this.camera.position.copy(t).add(p);
      this.dirty = true;
    }
    this.controls.update();
    if (!this.dirty) return;
    this.dirty = false;
    this.renderer.render(this.scene, this.camera);
  }
}

const srgb = (c: number): number => (c / 255) ** 2.2;
