import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { dpr, fmt, reducedMotion, asset } from "./util";

interface Blob { offset: number; length: number; dtype: string }
interface Body {
  kind: "head" | "shell";
  head?: number;
  modes: number;
  offset_coordinate: number;
  frequencies_hz: number[];
  nodes: Blob;
  shapes: Blob;
  tris: Blob;
  q: Blob;
}
interface VizHeader { bodies: Body[]; frames: number; dt_s: number; coordinates: number }

const DRUM_DEPTH = 0.1651;
const DRUM_RADIUS = 0.1703;

function view(buf: ArrayBuffer, b: Blob): Float32Array | Uint32Array {
  return b.dtype === "uint32" ? new Uint32Array(buf, b.offset, b.length) : new Float32Array(buf, b.offset, b.length);
}

const DOWN = new THREE.Color("#e86a3c");
const UP = new THREE.Color("#8fb9d6");
const REST = new THREE.Color("#3a3631");

interface Mesh3 { body: Body; geo: THREE.BufferGeometry; base: Float32Array; shapes: Float32Array; q: Float32Array; mesh: THREE.Mesh }

/** Replays REAL solver state: mode shapes of the actual mesh times integrated modal coordinates. */
export class ModalViewer {
  private renderer: THREE.WebGLRenderer | null = null;
  private scene = new THREE.Scene();
  private camera = new THREE.PerspectiveCamera(36, 1, 0.005, 10);
  private controls: OrbitControls | null = null;
  private meshes: Mesh3[] = [];
  private header: VizHeader | null = null;
  private frame = 0;
  private playing = false;
  private raf = 0;
  private scale = 1;
  private visible = false;
  private loadToken = 0;
  private framesPerTick = 2;

  constructor(private readonly stage: HTMLElement, private readonly ui: {
    play?: HTMLButtonElement | null;
    scrub?: HTMLInputElement | null;
    clock?: HTMLElement | null;
    scaleLabel?: HTMLElement | null;
    label?: HTMLElement | null;
  }) {
    this.scene.background = null;
    const key = new THREE.DirectionalLight(0xfff1dc, 1.5);
    key.position.set(0.4, -0.6, 1.2);
    const rim = new THREE.DirectionalLight(0x8fb9d6, 0.9);
    rim.position.set(-0.8, 0.7, 0.3);
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.35), key, rim);
    this.camera.up.set(0, 0, 1);
    ui.play?.addEventListener("click", () => this.toggle());
    ui.scrub?.addEventListener("input", () => {
      if (!this.header) return;
      this.pause();
      this.frame = Math.round((Number(ui.scrub!.value) / 1000) * (this.header.frames - 1));
      this.update();
    });
    new IntersectionObserver((e) => {
      this.visible = e.some((x) => x.isIntersecting);
      if (!this.visible) this.pause();
    }).observe(stage);
    document.addEventListener("visibilitychange", () => { if (document.hidden) this.pause(); });
    new ResizeObserver(() => this.resize()).observe(stage);
  }

  private ensureRenderer(): THREE.WebGLRenderer | null {
    if (this.renderer) return this.renderer;
    try {
      this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "high-performance" });
    } catch {
      this.stage.insertAdjacentHTML("beforeend", '<div class="stage-note">WebGL is unavailable in this browser, so the solver-state view cannot be drawn.</div>');
      return null;
    }
    this.renderer.setPixelRatio(dpr());
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.stage.prepend(this.renderer.domElement);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = false;
    this.controls.enablePan = false;
    this.controls.addEventListener("change", () => this.render());
    this.resize();
    return this.renderer;
  }

  async load(name: string, label?: string): Promise<void> {
    const token = ++this.loadToken;
    if (!this.ensureRenderer()) return;
    this.pause();
    const [header, buf] = await Promise.all([
      fetch(asset(`viz/${name}.json`)).then((r) => r.json() as Promise<VizHeader>),
      fetch(asset(`viz/${name}.bin`)).then((r) => r.arrayBuffer()),
    ]);
    if (token !== this.loadToken) return;
    for (const m of this.meshes) {
      this.scene.remove(m.mesh);
      m.geo.dispose();
    }
    this.scene.children.filter((c) => c.userData.decor).forEach((c) => this.scene.remove(c));
    this.meshes = [];
    this.header = header;
    let maxDisp = 1e-12;
    let span = 0.1;
    const isShell = header.bodies.some((b) => b.kind === "shell");
    for (const body of header.bodies) {
      const nodes = view(buf, body.nodes) as Float32Array;
      const shapes = view(buf, body.shapes) as Float32Array;
      const tris = view(buf, body.tris) as Uint32Array;
      const q = view(buf, body.q) as Float32Array;
      const nn = body.kind === "head" ? nodes.length / 2 : nodes.length / 3;
      const base = new Float32Array(nn * 3);
      for (let i = 0; i < nn; i++) {
        if (body.kind === "head") {
          base[3 * i] = nodes[2 * i]!;
          base[3 * i + 1] = nodes[2 * i + 1]!;
          base[3 * i + 2] = body.head === 0 ? DRUM_DEPTH / 2 : -DRUM_DEPTH / 2;
          span = DRUM_RADIUS;
        } else {
          base[3 * i] = nodes[3 * i]!;
          base[3 * i + 1] = nodes[3 * i + 1]!;
          base[3 * i + 2] = nodes[3 * i + 2]!;
          span = Math.max(span, Math.hypot(nodes[3 * i]!, nodes[3 * i + 1]!));
        }
      }
      // Peak displacement over the replay, sampled sparsely, to set one honest exaggeration.
      const stride = Math.max(1, Math.floor(header.frames / 400));
      const per = body.kind === "head" ? 1 : 3;
      for (let f = 0; f < header.frames; f += stride) {
        for (let i = 0; i < nn; i += 3) {
          let s = 0;
          for (let k = 0; k < body.modes; k++) s += shapes[(i * body.modes + k) * per + (per - 1)]! * q[f * body.modes + k]!;
          if (Math.abs(s) > maxDisp) maxDisp = Math.abs(s);
        }
      }
      const geo = new THREE.BufferGeometry();
      geo.setAttribute("position", new THREE.BufferAttribute(base.slice(), 3));
      geo.setAttribute("color", new THREE.BufferAttribute(new Float32Array(nn * 3), 3));
      geo.setIndex(new THREE.BufferAttribute(tris, 1));
      geo.computeVertexNormals();
      const mat = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: isShell ? 0.38 : 0.62, metalness: isShell ? 0.55 : 0.05, side: THREE.DoubleSide, transparent: body.kind === "head" && body.head === 1, opacity: body.head === 1 ? 0.55 : 1 });
      const mesh = new THREE.Mesh(geo, mat);
      this.scene.add(mesh);
      this.meshes.push({ body, geo, base, shapes, q, mesh });
    }
    if (!isShell) {
      const shell = new THREE.Mesh(
        new THREE.CylinderGeometry(DRUM_RADIUS + 0.008, DRUM_RADIUS + 0.008, DRUM_DEPTH, 72, 1, true),
        new THREE.MeshStandardMaterial({ color: 0x3a2a1a, roughness: 0.7, metalness: 0.1, side: THREE.DoubleSide, transparent: true, opacity: 0.28 }),
      );
      shell.rotation.x = Math.PI / 2;
      shell.userData.decor = true;
      this.scene.add(shell);
    }
    this.scale = (0.07 * span) / maxDisp;
    const dist = span * (isShell ? 3.4 : 3.7);
    this.camera.position.set(dist * 0.15, -dist * 0.82, dist * 0.62);
    this.controls?.target.set(0, 0, isShell ? 0.01 : 0);
    this.controls?.update();
    this.frame = 0;
    const slow = 1 / (60 * this.framesPerTick * header.dt_s);
    if (this.ui.scaleLabel) this.ui.scaleLabel.textContent = `displacement ×${fmt(this.scale, 2)} · ${fmt(slow, 2)}× slower than real time`;
    if (this.ui.label && label) this.ui.label.textContent = label;
    this.update();
    if (!reducedMotion() && this.visible) this.play();
  }

  private update(): void {
    if (!this.header) return;
    const f = Math.min(this.frame, this.header.frames - 1);
    for (const m of this.meshes) {
      const { body, base, shapes, q } = m;
      const pos = m.geo.getAttribute("position") as THREE.BufferAttribute;
      const col = m.geo.getAttribute("color") as THREE.BufferAttribute;
      const p = pos.array as Float32Array;
      const c = col.array as Float32Array;
      const nn = base.length / 3;
      const modes = body.modes;
      const qf = f * modes;
      const ref = 0.07 / this.scale;
      for (let i = 0; i < nn; i++) {
        let dx = 0, dy = 0, dz = 0;
        if (body.kind === "head") {
          for (let k = 0; k < modes; k++) dz -= shapes[i * modes + k]! * q[qf + k]!; // positive downward
        } else {
          const o = i * modes * 3;
          for (let k = 0; k < modes; k++) {
            const a = q[qf + k]!;
            dx += shapes[o + 3 * k]! * a;
            dy += shapes[o + 3 * k + 1]! * a;
            dz += shapes[o + 3 * k + 2]! * a;
          }
        }
        p[3 * i] = base[3 * i]! + dx * this.scale;
        p[3 * i + 1] = base[3 * i + 1]! + dy * this.scale;
        p[3 * i + 2] = base[3 * i + 2]! + dz * this.scale;
        const t = Math.max(-1, Math.min(1, dz / (ref * 0.6 || 1)));
        const target = t < 0 ? DOWN : UP;
        const a = Math.abs(t);
        c[3 * i] = REST.r + (target.r - REST.r) * a;
        c[3 * i + 1] = REST.g + (target.g - REST.g) * a;
        c[3 * i + 2] = REST.b + (target.b - REST.b) * a;
      }
      pos.needsUpdate = true;
      col.needsUpdate = true;
      m.geo.computeVertexNormals();
    }
    if (this.ui.clock) this.ui.clock.textContent = `t = ${(f * this.header.dt_s * 1000).toFixed(3)} ms`;
    if (this.ui.scrub) this.ui.scrub.value = String(Math.round((f / Math.max(1, this.header.frames - 1)) * 1000));
    this.ui.scrub?.style.setProperty("--fill", `${(f / Math.max(1, this.header.frames - 1)) * 100}%`);
    this.render();
  }

  private render(): void {
    this.renderer?.render(this.scene, this.camera);
  }

  private resize(): void {
    if (!this.renderer) return;
    const w = this.stage.clientWidth;
    const h = this.stage.clientHeight || w * 0.75;
    this.renderer.setSize(w, h, false);
    this.renderer.domElement.style.width = "100%";
    this.renderer.domElement.style.height = "100%";
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.render();
  }

  toggle(): void {
    if (this.playing) this.pause();
    else this.play();
  }

  play(): void {
    if (!this.header || this.playing) return;
    if (this.frame >= this.header.frames - 1) this.frame = 0;
    this.playing = true;
    if (this.ui.play) this.ui.play.textContent = "Pause";
    const step = (): void => {
      if (!this.playing || !this.header) return;
      this.frame += reducedMotion() ? 1 : this.framesPerTick;
      if (this.frame >= this.header.frames - 1) {
        this.frame = this.header.frames - 1;
        this.update();
        this.pause();
        return;
      }
      this.update();
      this.raf = requestAnimationFrame(step);
    };
    this.raf = requestAnimationFrame(step);
  }

  pause(): void {
    this.playing = false;
    cancelAnimationFrame(this.raf);
    if (this.ui.play) this.ui.play.textContent = this.header && this.frame > 0 && this.frame < this.header.frames - 1 ? "Resume" : "Play slow motion";
  }
}
