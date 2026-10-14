import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { dpr, reducedMotion } from "./util";

export interface ScaleRow { midi: number; unison: number; L: number; mu: number; T: number; EI: number; B: number; f1: number; strike: number; hammer_kg: number; felt_m: number; duplex_m: number; detune_cents: number }
export interface PianoModel {
  scale_et: ScaleRow[];
  scale_stretch: ScaleRow[];
  scale_raw: ScaleRow[];
  board: { nodes: number[]; tris: number[]; thickness_m: number[]; bridge_stations: Record<string, [number, number]> };
  board_obj: { vertices: number[]; groups: Record<string, number[]> };
  board_modes_400: number[];
}

const WOOD: Record<string, number> = { panel: 0xc39a62, rib: 0x8a6a42, bridge: 0xa86f3c, bar: 0x7a5a36 };

/** The Model D soundboard exactly as the solver meshed it, plus schematic string courses. */
export class BoardView {
  private renderer: THREE.WebGLRenderer | null = null;
  private scene = new THREE.Scene();
  private camera = new THREE.PerspectiveCamera(34, 4 / 3, 0.05, 30);
  private controls: OrbitControls | null = null;
  private strikeDots: THREE.Points | null = null;
  private keyIndex = new Map<number, number[]>();
  private visible = false;
  private raf = 0;
  private spinning = false;
  private lastInteraction = 0;
  readonly front: 1 | -1;
  readonly center: THREE.Vector3;

  constructor(private readonly stage: HTMLElement, private readonly model: PianoModel) {
    const st = model.board.bridge_stations;
    const rows = model.scale_et;
    // Which side is the keyboard? The agraffe line (bridge station + speaking length) is
    // nearly straight across the keys; choose the direction that makes it straightest.
    const spread = (sign: number): number => {
      const ys = rows.filter((r) => st[r.midi]).map((r) => st[r.midi]![1] + sign * r.L);
      const m = ys.reduce((a, b) => a + b, 0) / ys.length;
      return ys.reduce((a, b) => a + (b - m) ** 2, 0);
    };
    this.front = spread(1) <= spread(-1) ? 1 : -1;
    const v = model.board_obj.vertices;
    let cx = 0, cy = 0;
    for (let i = 0; i < v.length; i += 3) { cx += v[i]!; cy += v[i + 1]!; }
    this.center = new THREE.Vector3(cx / (v.length / 3), cy / (v.length / 3), 0);
  }

  mount(): boolean {
    try {
      this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "high-performance" });
    } catch {
      this.stage.insertAdjacentHTML("beforeend", '<div class="stage-note">WebGL is unavailable, so the soundboard cannot be drawn.</div>');
      return false;
    }
    const r = this.renderer;
    r.setPixelRatio(dpr());
    r.outputColorSpace = THREE.SRGBColorSpace;
    r.toneMapping = THREE.ACESFilmicToneMapping;
    r.toneMappingExposure = 1.05;
    this.stage.prepend(r.domElement);
    this.camera.up.set(0, 0, 1);
    this.build();
    this.controls = new OrbitControls(this.camera, r.domElement);
    this.controls.target.copy(this.center);
    this.controls.enablePan = false;
    this.controls.minDistance = 1.2;
    this.controls.maxDistance = 6;
    this.controls.addEventListener("start", () => { this.lastInteraction = performance.now(); this.kick(); });
    this.controls.addEventListener("change", () => this.render());
    const d = 2.55;
    this.camera.position.set(this.center.x + d * 0.55, this.center.y + this.front * d * 0.9, d * 1.1);
    this.controls.update();
    new ResizeObserver(() => this.resize()).observe(this.stage);
    new IntersectionObserver((e) => {
      this.visible = e.some((x) => x.isIntersecting);
      this.kick();
    }).observe(this.stage);
    document.addEventListener("visibilitychange", () => this.kick());
    this.resize();
    return true;
  }

  private build(): void {
    const s = this.scene;
    s.add(new THREE.HemisphereLight(0xfff3df, 0x1a1410, 0.75));
    const key = new THREE.DirectionalLight(0xffe7c2, 2.4);
    key.position.set(1.8, -1.2, 3.2);
    const rim = new THREE.DirectionalLight(0x8fb9d6, 1.0);
    rim.position.set(-2.2, 2.4, 1.2);
    s.add(key, rim);

    const { vertices, groups } = this.model.board_obj;
    const pos = new Float32Array(vertices);
    for (const [name, idx] of Object.entries(groups)) {
      const kind = name.includes("panel") ? "panel" : name.includes("rib") ? "rib" : name.includes("bridge") ? "bridge" : "bar";
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
      g.setIndex(idx);
      g.computeVertexNormals();
      const m = new THREE.MeshStandardMaterial({ color: WOOD[kind], roughness: kind === "panel" ? 0.58 : 0.66, metalness: 0.02, side: THREE.DoubleSide, flatShading: kind !== "panel" });
      const mesh = new THREE.Mesh(g, m);
      mesh.userData.kind = kind;
      s.add(mesh);
    }

    // The 88 bridge stations where each key's strings drive the board (from the FSB dump).
    // The model's strings are 1-D, so no string layout is drawn: only these real points.
    const st = this.model.board.bridge_stations;
    const dots: number[] = [];
    for (const row of this.model.scale_et) {
      const p = st[row.midi];
      if (!p) continue;
      this.keyIndex.set(row.midi, [dots.length / 3]);
      dots.push(p[0], p[1], 0.052);
    }
    const dg = new THREE.BufferGeometry();
    dg.setAttribute("position", new THREE.Float32BufferAttribute(dots, 3));
    dg.setAttribute("color", new THREE.Float32BufferAttribute(new Array(dots.length).fill(0.62), 3));
    this.strikeDots = new THREE.Points(dg, new THREE.PointsMaterial({ size: 0.026, vertexColors: true, transparent: true, opacity: 0.95, sizeAttenuation: true, depthTest: false }));
    s.add(this.strikeDots);
  }

  /** intensity[midi] in 0..1: lights the struck courses and their hammer points. */
  setActivity(intensity: Map<number, number>): void {
    if (!this.strikeDots) return;
    const attr = this.strikeDots.geometry.getAttribute("color") as THREE.BufferAttribute;
    const dc = attr.array as Float32Array;
    for (const [midi, idx] of this.keyIndex) {
      const a = Math.min(1, intensity.get(midi) ?? 0);
      const o = idx[0]! * 3;
      dc[o] = 0.62 + a * 1.6;
      dc[o + 1] = 0.6 + a * 0.95;
      dc[o + 2] = 0.56 - a * 0.2;
    }
    attr.needsUpdate = true;
    this.render();
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

  private kick(): void {
    const shouldSpin = this.visible && !document.hidden && !reducedMotion();
    if (shouldSpin && !this.spinning) {
      this.spinning = true;
      let last = performance.now();
      const began = last;
      const tick = (now: number): void => {
        if (!this.spinning) return;
        // An idle page should cost nothing: stop turning 25 s after the last start or touch.
        if (now - Math.max(began, this.lastInteraction) > 25000) {
          this.spinning = false;
          this.render();
          return;
        }
        const dt = Math.min(0.05, (now - last) / 1000);
        last = now;
        if (this.controls && now - this.lastInteraction > 4000) {
          const t = this.controls.target;
          const p = this.camera.position.clone().sub(t);
          p.applyAxisAngle(new THREE.Vector3(0, 0, 1), dt * 0.06);
          this.camera.position.copy(t).add(p);
          this.camera.lookAt(t);
          this.controls.update();
        }
        this.render();
        this.raf = requestAnimationFrame(tick);
      };
      this.raf = requestAnimationFrame(tick);
    } else if (!shouldSpin) {
      this.spinning = false;
      cancelAnimationFrame(this.raf);
      this.render();
    }
  }

  render(): void {
    this.renderer?.render(this.scene, this.camera);
  }
}

/** An 88-key strip (A0..C8) with absolutely positioned keys. */
export function buildKeyboard(host: HTMLElement): Map<number, HTMLElement> {
  const keys = new Map<number, HTMLElement>();
  const isBlack = (m: number): boolean => [1, 3, 6, 8, 10].includes(m % 12);
  const whites: number[] = [];
  for (let m = 21; m <= 108; m++) if (!isBlack(m)) whites.push(m);
  const ww = 100 / whites.length;
  for (let m = 21; m <= 108; m++) {
    const k = document.createElement("div");
    if (isBlack(m)) {
      const left = whites.indexOf(m - 1);
      k.className = "k";
      k.style.left = `${(left + 1) * ww - ww * 0.32}%`;
      k.style.width = `${ww * 0.64}%`;
    } else {
      k.className = "w";
      k.style.left = `${whites.indexOf(m) * ww}%`;
      k.style.width = `${ww}%`;
    }
    host.append(k);
    keys.set(m, k);
  }
  return keys;
}
