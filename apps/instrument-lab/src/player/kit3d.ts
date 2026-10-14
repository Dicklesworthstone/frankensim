// The drum kit in 3D. Each piece the file plays is the fs-vkplate solver in its own Worker; its
// surface here is the solver's own grid, lifted and painted with the computed transverse velocity
// of every node (warm up, cool down), so the wave you see crossing a head is the one you hear.
// Heights are exaggerated to a fixed fraction of each piece's radius, and each piece's colour scale
// follows its own motion down as it decays (stated on the panel).

import * as THREE from "three";
import { RoomEnvironment } from "three/examples/jsm/environments/RoomEnvironment.js";
import type { KitField, KitReady } from "../live/mix";
import { PIECE_NAMES, type PieceId } from "../midi/arrange";
import { reducedMotion } from "../util";

interface Spec { pos: [number, number, number]; tilt: number; r: number; kind: "drum" | "kick" | "cymbal" | "hat"; depth: number; finish?: number }
// metres, y up, +z toward the listener; radii are the kit's real sizes until the solver reports its own
const SPECS: Record<PieceId, Spec> = {
  kick: { pos: [0, 0.3, 0.0], tilt: 0, r: 0.229, kind: "kick", depth: 0.42, finish: 0x5c1d1a },
  snare: { pos: [-0.46, 0.64, 0.3], tilt: 0.16, r: 0.178, kind: "drum", depth: 0.14, finish: 0xc9cdd3 },
  "snare-off": { pos: [-0.46, 0.64, 0.3], tilt: 0.16, r: 0.178, kind: "drum", depth: 0.14, finish: 0xc9cdd3 },
  "tom-12": { pos: [-0.1, 0.84, -0.06], tilt: 0.36, r: 0.152, kind: "drum", depth: 0.2, finish: 0x5c1d1a },
  "tom-14": { pos: [0.52, 0.55, 0.26], tilt: 0.1, r: 0.178, kind: "drum", depth: 0.36, finish: 0x5c1d1a },
  "hihat-closed": { pos: [-0.86, 0.93, 0.08], tilt: 0.08, r: 0.178, kind: "hat", depth: 0 },
  "hihat-open": { pos: [-0.86, 0.93, 0.08], tilt: 0.08, r: 0.178, kind: "hat", depth: 0 },
  crash: { pos: [-0.56, 1.34, -0.36], tilt: 0.42, r: 0.229, kind: "cymbal", depth: 0 },
  ride: { pos: [0.66, 1.18, -0.24], tilt: 0.3, r: 0.254, kind: "cymbal", depth: 0 },
};
const ORDER: PieceId[] = ["kick", "tom-12", "tom-14", "snare", "hihat-closed", "crash", "ride"];
/** Two solver pieces at one place: the one struck last is shown. */
const TWIN: Partial<Record<PieceId, PieceId>> = { "hihat-open": "hihat-closed", "hihat-closed": "hihat-open", "snare-off": "snare", snare: "snare-off" };

interface Surface {
  geom: THREE.BufferGeometry;
  base: Float32Array; // resting height of every vertex (cymbal bow and bell)
  cellOf: Int32Array; // grid vertex of every solver node
  has: Uint8Array;
  m: number;
  R: number;
  h: number;
  /** true: the field is light added on top (kick, cymbals); false: paint multiplied into a white head */
  glow: boolean;
}
interface PieceState {
  id: PieceId;
  group: THREE.Group;
  surface: Surface | null;
  queue: KitField[];
  shown: KitField | null;
  vref: number;
  kpeak: number;
  kinetic: number;
  force: number;
  label: HTMLElement;
  used: boolean;
  flashes: { at: number; u: number; v: number }[];
  ring: THREE.Mesh;
}

export class KitView {
  private renderer: THREE.WebGLRenderer | null = null;
  private scene = new THREE.Scene();
  private camera = new THREE.PerspectiveCamera(30, 1.5, 0.05, 40);
  private pieces = new Map<PieceId, PieceState>();
  private labels: HTMLElement;
  private dirty = true;
  private target = new THREE.Vector3(-0.1, 0.72, 0); // the kit's middle: the hi-hat reaches further out than the ride
  private orbit = 0;
  private scratch = new Float32Array(0);
  ok = false;

  constructor(private readonly host: HTMLElement) {
    this.labels = document.createElement("div");
    this.labels.className = "pp-kit-labels";
    host.append(this.labels);
    try {
      this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "high-performance" });
    } catch {
      host.insertAdjacentHTML("beforeend", '<div class="stage-note">WebGL is unavailable here, so the kit cannot be drawn.</div>');
      return;
    }
    const r = this.renderer;
    r.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    r.outputColorSpace = THREE.SRGBColorSpace;
    r.toneMapping = THREE.ACESFilmicToneMapping;
    r.toneMappingExposure = 1.0;
    r.domElement.className = "pp-gl";
    r.domElement.style.touchAction = "pan-y";
    host.prepend(r.domElement);
    const pmrem = new THREE.PMREMGenerator(r);
    this.scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
    this.scene.add(new THREE.HemisphereLight(0xfff0dc, 0x120c08, 0.35));
    const key = new THREE.DirectionalLight(0xffe6c4, 1.6);
    key.position.set(1.2, 3.2, 2.4);
    this.scene.add(key);
    this.buildRoom();
    for (const id of [...ORDER, "snare-off", "hihat-open"] as PieceId[]) this.buildPiece(id);
    new ResizeObserver(() => this.resize()).observe(host);
    this.resize();
    this.ok = true;
  }

  private buildRoom(): void {
    const c = document.createElement("canvas");
    c.width = c.height = 256;
    const g = c.getContext("2d")!;
    const rg = g.createRadialGradient(128, 128, 10, 128, 128, 128);
    rg.addColorStop(0, "rgba(122, 40, 32, 0.55)"); rg.addColorStop(0.7, "rgba(80, 26, 22, 0.35)"); rg.addColorStop(1, "rgba(0, 0, 0, 0)");
    g.fillStyle = rg;
    g.fillRect(0, 0, 256, 256);
    const rug = new THREE.Mesh(new THREE.CircleGeometry(1.6, 64), new THREE.MeshBasicMaterial({ map: new THREE.CanvasTexture(c), transparent: true, depthWrite: false }));
    rug.rotation.x = -Math.PI / 2;
    rug.position.set(0, 0.001, 0.05);
    this.scene.add(rug);
  }

  private buildPiece(id: PieceId): void {
    const s = SPECS[id];
    const group = new THREE.Group();
    group.position.set(...s.pos);
    // the playing surface lies in the group's local xy plane, its normal along local +z
    if (s.kind === "kick") group.rotation.set(0, 0, 0);
    else group.rotation.set(-Math.PI / 2 + s.tilt, 0, 0);
    const chrome = new THREE.MeshStandardMaterial({ color: 0xdfe3e8, metalness: 1, roughness: 0.22 });
    if (s.kind === "drum" || s.kind === "kick") {
      const shell = new THREE.Mesh(new THREE.CylinderGeometry(s.r * 1.005, s.r * 1.005, s.depth, 48, 1, true),
        new THREE.MeshStandardMaterial({ color: s.finish ?? 0x5c1d1a, metalness: s.finish === 0xc9cdd3 ? 1 : 0.35, roughness: s.finish === 0xc9cdd3 ? 0.2 : 0.38, side: THREE.DoubleSide }));
      shell.rotation.x = Math.PI / 2;
      shell.position.z = -s.depth / 2;
      group.add(shell);
      for (const z of [0, -s.depth]) {
        const hoop = new THREE.Mesh(new THREE.TorusGeometry(s.r * 1.012, s.r * 0.022, 10, 64), chrome);
        hoop.position.z = z + (z === 0 ? 0.004 : 0);
        group.add(hoop);
      }
      // lugs
      for (let k = 0; k < 8; k++) {
        const a = (k / 8) * Math.PI * 2;
        const lug = new THREE.Mesh(new THREE.BoxGeometry(s.r * 0.06, s.r * 0.06, s.depth * 0.35), chrome);
        lug.position.set(Math.cos(a) * s.r * 1.05, Math.sin(a) * s.r * 1.05, -s.depth * 0.5);
        group.add(lug);
      }
      if (s.kind === "kick") {
        // the rear head, and spurs to the floor
        const back = new THREE.Mesh(new THREE.CircleGeometry(s.r, 48), new THREE.MeshStandardMaterial({ color: 0x15151a, roughness: 0.8 }));
        back.position.z = -s.depth;
        back.rotation.y = Math.PI;
        group.add(back);
      }
    } else {
      // stand
      const stand = new THREE.Mesh(new THREE.CylinderGeometry(0.009, 0.011, 1, 10), new THREE.MeshStandardMaterial({ color: 0x8d939b, metalness: 0.9, roughness: 0.3 }));
      const tube = new THREE.Group();
      stand.scale.y = s.pos[1];
      stand.position.y = -s.pos[1] / 2;
      tube.add(stand);
      tube.position.set(...s.pos);
      this.scene.add(tube);
      if (s.kind === "hat") {
        const bottom = new THREE.Mesh(new THREE.CircleGeometry(s.r, 48), new THREE.MeshStandardMaterial({ color: 0x8a6a2e, metalness: 0.95, roughness: 0.32, side: THREE.DoubleSide }));
        bottom.position.z = -0.018;
        group.add(bottom);
      }
    }
    // the surface, until the solver supplies its own grid: a plain head or a bronze disc
    const placeholder = new THREE.Mesh(new THREE.CircleGeometry(s.r, 48), this.surfaceMaterial(s.kind));
    placeholder.name = "placeholder";
    group.add(placeholder);
    const ring = new THREE.Mesh(new THREE.RingGeometry(0.01, 0.014, 32), new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0, depthWrite: false, side: THREE.DoubleSide }));
    ring.position.z = 0.01;
    group.add(ring);
    this.scene.add(group);
    const label = document.createElement("span");
    label.className = "pp-kit-label";
    label.textContent = PIECE_NAMES[id];
    this.labels.append(label);
    const hidden = id === "snare-off" || id === "hihat-open";
    group.visible = !hidden;
    label.hidden = hidden;
    if (id === "hihat-open") group.position.y += 0.012; // the top cymbal lifted off the bottom one
    this.pieces.set(id, { id, group, surface: null, queue: [], shown: null, vref: 1e-6, kpeak: 0, kinetic: 0, force: 0, label, used: false, flashes: [], ring });
  }

  private surfaceMaterial(kind: Spec["kind"]): THREE.Material {
    if (kind === "cymbal" || kind === "hat") return new THREE.MeshStandardMaterial({ color: 0xc9973f, metalness: 0.95, roughness: 0.3, vertexColors: false, side: THREE.DoubleSide });
    if (kind === "kick") return new THREE.MeshStandardMaterial({ color: 0x1e1e24, roughness: 0.7, side: THREE.DoubleSide });
    return new THREE.MeshStandardMaterial({ color: 0xece6d8, roughness: 0.75, side: THREE.DoubleSide });
  }

  /** Which pieces this file plays; the rest of the kit stays, dimmed. */
  setUsed(used: PieceId[]): void {
    for (const p of this.pieces.values()) {
      p.used = used.includes(p.id);
      this.dimPiece(p, !p.used && !(TWIN[p.id] && used.includes(TWIN[p.id]!)));
    }
    // a twin (open hat, snare with its wires off) shares its partner's place: show the main one
    // unless only the twin is played
    for (const [alt, main] of [["hihat-open", "hihat-closed"], ["snare-off", "snare"]] as [PieceId, PieceId][]) {
      this.showTwin(used.includes(main) || !used.includes(alt) ? main : alt);
    }
    this.requestRender();
  }
  private renderQueued = false;
  /** Draw once on the next frame (the player's clock only runs while it plays). */
  private requestRender(): void {
    this.dirty = true;
    if (this.renderQueued || !this.renderer) return;
    this.renderQueued = true;
    requestAnimationFrame(() => {
      this.renderQueued = false;
      if (!this.dirty || !this.renderer) return;
      this.dirty = false;
      this.renderer.render(this.scene, this.camera);
      this.lastLabels = 0;
      this.placeLabels();
    });
  }
  private showTwin(id: PieceId): void {
    const other = TWIN[id];
    if (!other) return;
    const a = this.pieces.get(id)!, b = this.pieces.get(other)!;
    if (a.group.visible && !b.group.visible) return;
    a.group.visible = true; a.label.hidden = false;
    b.group.visible = false; b.label.hidden = true;
    this.dirty = true;
  }
  private dimPiece(p: PieceState, dim: boolean): void {
    p.group.traverse((o) => {
      const m = (o as THREE.Mesh).material as THREE.MeshStandardMaterial | undefined;
      if (!m || !("opacity" in m)) return;
      if (o === p.ring) return;
      m.transparent = dim || m.transparent;
      m.opacity = dim ? 0.28 : 1;
      m.needsUpdate = true;
    });
    p.label.classList.toggle("dim", dim);
  }

  /** The solver's grid for a piece: its surface replaces the placeholder. */
  setReady(id: PieceId, r: KitReady): void {
    const p = this.pieces.get(id);
    if (!p || !this.renderer) return;
    const s = SPECS[p.id];
    const m = 2 * Math.ceil(r.radius / r.h) + 3;
    const nodes = r.xy.length / 2;
    const cellOf = new Int32Array(nodes);
    const has = new Uint8Array(m * m);
    for (let k = 0; k < nodes; k++) {
      const i = Math.round(r.xy[2 * k]! / r.h) + (m - 1) / 2, j = Math.round(r.xy[2 * k + 1]! / r.h) + (m - 1) / 2;
      cellOf[k] = j * m + i;
      has[j * m + i] = 1;
    }
    const pos = new Float32Array(m * m * 3);
    const base = new Float32Array(m * m);
    for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) {
      const v = j * m + i, x = (i - (m - 1) / 2) * r.h, y = (j - (m - 1) / 2) * r.h;
      pos[3 * v] = x; pos[3 * v + 1] = y;
      const rr = Math.hypot(x, y) / r.radius;
      // a cymbal's bow falls toward the edge and its bell rises in the middle
      if (s.kind === "cymbal" || s.kind === "hat") base[v] = r.radius * (0.05 * (1 - Math.min(1, rr) ** 1.4) + (rr < 0.2 ? 0.05 * (1 - rr / 0.2) ** 0.7 : 0));
      pos[3 * v + 2] = base[v]!;
    }
    // the grid's outline is a staircase: move every node on it out to the rim, so the surface
    // meets its hoop (or ends at the cymbal's edge) in a circle; its field values are unchanged
    for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) {
      const v = j * m + i;
      if (!has[v]) continue;
      const edge = !has[v - 1] || !has[v + 1] || !has[v - m] || !has[v + m] || i === 0 || j === 0 || i === m - 1 || j === m - 1;
      if (!edge) continue;
      const x = pos[3 * v]!, y = pos[3 * v + 1]!, rr = Math.hypot(x, y);
      // a drum head ends just inside its hoop (no two surfaces fighting for the same pixels)
      const rim = r.radius * (s.kind === "drum" || s.kind === "kick" ? 0.975 : 1);
      if (rr > 0.5 * r.radius) { pos[3 * v] = (x / rr) * rim; pos[3 * v + 1] = (y / rr) * rim; }
    }
    // and nothing inside reaches past the rim (a fold there shades as a dark wedge)
    const rim = r.radius * (s.kind === "drum" || s.kind === "kick" ? 0.975 : 1);
    for (let v = 0; v < m * m; v++) {
      const x = pos[3 * v]!, y = pos[3 * v + 1]!, rr = Math.hypot(x, y);
      if (rr > rim) { pos[3 * v] = (x / rr) * rim; pos[3 * v + 1] = (y / rr) * rim; }
    }
    const idx: number[] = [];
    for (let j = 0; j < m - 1; j++) for (let i = 0; i < m - 1; i++) {
      const a = j * m + i, b = a + 1, c = a + m, d = c + 1;
      if (has[a] && has[b] && has[d]) idx.push(a, b, d);
      if (has[a] && has[d] && has[c]) idx.push(a, d, c);
    }
    const glow = s.kind !== "drum"; // a dark kick head and bronze glow; white heads take paint
    const geom = new THREE.BufferGeometry();
    geom.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    geom.setAttribute("color", new THREE.BufferAttribute(new Float32Array(m * m * 3).fill(glow ? 0 : 1), 3));
    geom.setIndex(idx);
    geom.computeVertexNormals();
    const mat = this.surfaceMaterial(s.kind) as THREE.MeshStandardMaterial;
    mat.vertexColors = !glow;
    const mesh = new THREE.Mesh(geom, mat);
    mesh.scale.setScalar(s.r / r.radius); // drawn at the kit's size; the field is the solver's
    const old = p.group.getObjectByName("placeholder");
    if (old) p.group.remove(old);
    mesh.name = "surface";
    p.group.add(mesh);
    if (glow) {
      const light = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({ vertexColors: true, transparent: true, blending: THREE.AdditiveBlending, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -2, side: THREE.DoubleSide }));
      light.scale.copy(mesh.scale);
      light.name = "glow";
      light.renderOrder = 2;
      p.group.add(light);
    }
    p.surface = { geom, base, cellOf, has, m, R: r.radius, h: r.h, glow };
    this.dimPiece(p, !p.used);
    this.requestRender();
  }

  /** A field computed for context frame `f.frame`: shown when the audio clock reaches it. */
  field(id: PieceId, f: KitField): void {
    const p = this.pieces.get(id);
    if (!p) return;
    p.queue.push(f);
    if (p.queue.length > 120) p.queue.splice(0, p.queue.length - 120);
  }

  /** A stroke, as it is heard: a ring where the stick lands (and a twin takes its partner's place). */
  stroke(id: PieceId, u: number, v: number, at: number): void {
    const p = this.pieces.get(id);
    if (!p) return;
    if (TWIN[id]) this.showTwin(id);
    p.flashes.push({ at, u, v });
    if (p.flashes.length > 8) p.flashes.shift();
  }

  /** Clear everything queued (after a stop or a jump). */
  clear(): void {
    for (const p of this.pieces.values()) { p.queue.length = 0; p.flashes.length = 0; p.kinetic = 0; p.shown = null; this.applyField(p, null, 0); }
    this.requestRender();
  }

  private applyField(p: PieceState, f: KitField | null, fade: number): void {
    const s = p.surface;
    if (!s) return;
    const spec = SPECS[p.id];
    const pos = s.geom.getAttribute("position") as THREE.BufferAttribute;
    const col = s.geom.getAttribute("color") as THREE.BufferAttribute;
    const P = pos.array as Float32Array, C = col.array as Float32Array;
    const rest = s.glow ? 0 : 1;
    for (let v = 0; v < s.m * s.m; v++) { P[3 * v + 2] = s.base[v]!; C[3 * v] = C[3 * v + 1] = C[3 * v + 2] = rest; }
    if (f && fade > 0.002) {
      let vmax = 0;
      for (let k = 0; k < f.v.length; k++) vmax = Math.max(vmax, Math.abs(f.v[k]!));
      p.vref = Math.max(vmax, p.vref * 0.9, 1e-7);
      const lift = (spec.kind === "drum" || spec.kind === "kick" ? 0.07 : 0.1) * s.R;
      // the field on the grid, then one 3×3 mean over the nodes that exist (the live kit's 2D
      // painter blurs by about half a cell for the same reason: a coarse grid drawn cell by cell
      // looks like noise rather than a wave)
      const m = s.m, raw = (this.scratch.length >= m * m ? this.scratch : (this.scratch = new Float32Array(m * m)));
      raw.fill(0, 0, m * m);
      for (let k = 0; k < f.v.length && k < s.cellOf.length; k++) raw[s.cellOf[k]!] = f.v[k]! / p.vref;
      for (let k = 0; k < f.v.length && k < s.cellOf.length; k++) {
        const v = s.cellOf[k]!;
        let sum = 0, n = 0;
        for (let dj = -m; dj <= m; dj += m) for (let di = -1; di <= 1; di++) { const w = v + dj + di; if (w >= 0 && w < m * m && s.has[w]) { sum += raw[w]! * (di === 0 && dj === 0 ? 2 : 1); n += di === 0 && dj === 0 ? 2 : 1; } }
        const u = sum / n;
        P[3 * v + 2] = s.base[v]! + u * lift * fade;
        const a = Math.min(1, Math.abs(u) * 1.8) ** 0.85 * fade;
        if (s.glow) {
          // light added on bronze or on the dark kick head: gold up, cyan down
          const [r, g, b] = u >= 0 ? [1.0, 0.55, 0.12] : [0.12, 0.8, 0.75];
          const lit = 0.5 * a; // added light: kept below white so the shape stays readable
          C[3 * v] = r * lit; C[3 * v + 1] = g * lit; C[3 * v + 2] = b * lit;
        } else {
          // paint on a coated white head: ember up, deep verdigris down
          const [r, g, b] = u >= 0 ? [0.95, 0.26, 0.08] : [0.04, 0.48, 0.44];
          C[3 * v] = 1 + (r - 1) * a; C[3 * v + 1] = 1 + (g - 1) * a; C[3 * v + 2] = 1 + (b - 1) * a;
        }
      }
    }
    pos.needsUpdate = true;
    col.needsUpdate = true;
    s.geom.computeVertexNormals();
  }

  private resize(): void {
    if (!this.renderer) return;
    const w = this.host.clientWidth || 600, h = this.host.clientHeight || 400;
    this.renderer.setSize(w, h, false);
    this.renderer.domElement.style.width = "100%";
    this.renderer.domElement.style.height = "100%";
    this.camera.aspect = w / h;
    // narrow views step back so the whole kit stays in frame
    this.back = Math.max(1, 1.45 / (w / h)); // keep the kit's full width in a narrow frame
    this.place();
    this.camera.updateProjectionMatrix();
    this.requestRender();
  }
  private back = 1;
  private place(): void {
    const d = 3.05 * this.back;
    this.camera.position.set(this.target.x + Math.sin(this.orbit) * d, 2.05 + 0.3 * (this.back - 1), Math.cos(this.orbit) * d);
    this.camera.lookAt(this.target);
  }

  /** `nowFrame`: the context frame being heard; `wall`: performance.now(). */
  frame(nowFrame: number, wall: number, dt: number, sway: boolean): void {
    if (!this.renderer) return;
    for (const p of this.pieces.values()) {
      let next: KitField | null = null;
      while (p.queue.length && p.queue[0]!.frame <= nowFrame) next = p.queue.shift()!;
      if (next) {
        p.shown = next;
        p.kinetic = next.kinetic;
        p.kpeak = Math.max(p.kpeak * 0.995, next.kinetic);
        p.force = Math.max(p.force * 0.9, next.force);
      } else if (p.shown && nowFrame - p.shown.frame > 48000 * 0.12) {
        p.kinetic *= 0.8; // the piece fell asleep: no more fields
      }
      const fade = p.kpeak > 0 && p.shown ? Math.max(0, Math.min(1, (Math.log10(Math.max(p.kinetic, 1e-30) / p.kpeak) + 4.5) / 4.5)) : 0;
      if (next || (p.shown && fade < 0.01 && p.kinetic > 0)) { this.applyField(p, p.shown, fade); this.dirty = true; }
      // stick marks
      const ringMat = p.ring.material as THREE.MeshBasicMaterial;
      const f = p.flashes.filter((x) => x.at <= wall);
      const last = f[f.length - 1];
      if (last && wall - last.at < 260) {
        const k = (wall - last.at) / 260;
        const R = SPECS[p.id].r;
        p.ring.position.set(last.u * R, last.v * R, 0.012 + (p.surface ? 0 : 0));
        p.ring.scale.setScalar(1 + k * 4);
        ringMat.opacity = 0.85 * (1 - k);
        this.dirty = true;
      } else if (ringMat.opacity > 0) { ringMat.opacity = 0; this.dirty = true; }
      p.flashes = p.flashes.filter((x) => wall - x.at < 400);
    }
    if (sway && !reducedMotion()) {
      this.orbit = Math.sin(wall / 9000) * 0.18;
      this.place();
      this.dirty = true;
    }
    void dt;
    if (!this.dirty) return;
    this.dirty = false;
    this.renderer.render(this.scene, this.camera);
    this.placeLabels();
  }

  private lastLabels = 0;
  private placeLabels(): void {
    const now = performance.now();
    if (now - this.lastLabels < 90) return;
    this.lastLabels = now;
    const w = this.host.clientWidth, h = this.host.clientHeight;
    const v = new THREE.Vector3();
    for (const p of this.pieces.values()) {
      if (!p.group.visible) continue;
      const s = SPECS[p.id];
      // under the front of a drum's shell, below a cymbal's edge
      if (s.kind === "drum") v.set(0, -s.r * 1.02, -s.depth * 0.9);
      else v.set(0, -s.r * 1.15, 0);
      p.group.localToWorld(v);
      if (s.kind === "kick") v.set(s.pos[0], 0.04, s.pos[2] + 0.1);
      v.project(this.camera);
      // labels stay inside the frame (half a label's width from either side)
      const half = (p.label.offsetWidth || 80) / 2 + 6;
      const lx = Math.max(half, Math.min(w - half, ((v.x + 1) / 2) * w));
      p.label.style.transform = `translate(${lx}px, ${Math.min(h - 26, ((1 - v.y) / 2) * h)}px) translate(-50%, 0)`;
      const e = p.kinetic;
      p.label.classList.toggle("on", e > 1e-9 && p.used);
      p.label.dataset.energy = e > 1e-9 ? fmtJ(e) : "";
    }
  }

  /** The piece struck most recently that is still moving, for the readout. */
  energyOf(id: PieceId): { kinetic: number; force: number } {
    const p = this.pieces.get(id);
    return { kinetic: p?.kinetic ?? 0, force: p?.force ?? 0 };
  }
}

export function fmtJ(j: number): string {
  const f = (v: number, u: string): string => `${v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v.toFixed(2)} ${u}`;
  return j >= 1 ? f(j, "J") : j >= 1e-3 ? f(j * 1e3, "mJ") : j >= 1e-6 ? f(j * 1e6, "µJ") : j >= 1e-9 ? f(j * 1e9, "nJ") : "at rest";
}
