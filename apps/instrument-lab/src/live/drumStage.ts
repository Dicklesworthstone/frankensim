// three.js stage for the live drum: the real finite-element head mesh deformed
// by the engine's accepted modal state, a shell and hoops drawn from the real
// radius and depth, strike picking in head-plane metres and a small orbit.
//
// One rAF loop at most. It runs only while the stage is on screen and the tab
// is visible, and it renders only when something changed.

import * as THREE from 'three';

export type StageView = 'batter' | 'resonant';
export type StageMode = 'strike' | 'orbit';

export interface StageCallbacks {
  /** A strike at head-centred metres (engine frame: x right, y away from the viewer at rest). */
  strike(xM: number, yM: number, source: 'pointer' | 'keyboard'): void;
  /** The displayed exaggeration factor changed (0 = the heads are at rest). */
  scale(factor: number): void;
  /** Transient hint for the HUD (empty string clears it). */
  hint(text: string): void;
}

const TAP_SLOP_PX = 6;
const TAP_MS = 250;
/** Bearing edge to hoop: the reference's outer radius is the clear radius + 7.5 mm. */
const SHELL_OFFSET_M = 0.0075;
const HOOP_TUBE_M = 0.0045;
/** Shown amplitude target as a fraction of the radius, and the envelope release. */
const SHOWN_FRACTION = 0.16;
const RELEASE_S = 2.5;
const FLOOR_M = 2e-6;

const C_UP = new THREE.Color('#fff7e4');
const C_REST = new THREE.Color('#a39c8b');
const C_DOWN = new THREE.Color('#e2652c');

function niceBelow(value: number): number {
  if (!(value > 0) || !Number.isFinite(value)) return 1;
  const exponent = Math.floor(Math.log10(value));
  const base = 10 ** exponent;
  for (const m of [5, 2, 1]) if (m * base <= value) return m * base;
  return base;
}

interface PointerTrack {
  x: number;
  y: number;
  startX: number;
  startY: number;
  t0: number;
  type: string;
}

export class DrumStage {
  readonly canvas: HTMLCanvasElement;
  private readonly host: HTMLElement;
  private readonly cb: StageCallbacks;
  private readonly renderer: THREE.WebGLRenderer | null;
  private readonly scene = new THREE.Scene();
  private readonly camera = new THREE.PerspectiveCamera(30, 4 / 3, 0.01, 20);
  private readonly drum = new THREE.Group();
  private readonly reducedMotion: boolean;

  private radius = 0.1703;
  private depth = 0.1651;
  private shell: THREE.Mesh | null = null;
  private collars: THREE.Mesh[] = [];
  private hoops: THREE.Mesh[] = [];
  private heads: [THREE.Mesh, THREE.Mesh] | null = null;
  private wires: THREE.LineSegments[] = [];
  private placeholder: THREE.Mesh[] = [];
  private nodes: Float32Array | null = null;
  private displacement: [Float32Array, Float32Array] | null = null;
  private colors: [Float32Array, Float32Array] | null = null;

  private envelope = FLOOR_M;
  private factor = 1;
  private factorChangedAt = 0;
  private lastPeak = 0;
  private shownFactor = -1;

  private marker: THREE.Mesh;
  private markerBorn = -1e9;
  private aim: THREE.Mesh;

  private yaw = THREE.MathUtils.degToRad(28);
  private pitch = THREE.MathUtils.degToRad(31);
  private view: StageView = 'batter';
  private viewBlend = 0;
  private viewTarget = 0;
  private mode: StageMode = 'strike';
  private started = false;

  private readonly pointers = new Map<number, PointerTrack>();
  private gesture: 'none' | 'strike' | 'pending' | 'orbit' | 'orbit2' = 'none';
  private pendingHit: [number, number] | null = null;

  private visible = false;
  private pageVisible = !document.hidden;
  private raf = 0;
  private lastTime = 0;
  private dirty = true;
  private disposed = false;
  private readonly resizeObserver: ResizeObserver;
  private readonly intersection: IntersectionObserver;
  private readonly onVisibility = (): void => {
    this.pageVisible = !document.hidden;
    this.schedule();
  };

  constructor(host: HTMLElement, callbacks: StageCallbacks) {
    this.host = host;
    this.cb = callbacks;
    this.reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    let renderer: THREE.WebGLRenderer | null = null;
    try {
      renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: 'high-performance' });
      renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
      renderer.setClearColor(0x000000, 0);
    } catch {
      renderer = null;
    }
    this.renderer = renderer;
    // Without WebGL the drum still plays: strikes are raycast in metres on this
    // canvas and the audio engine never needed the GPU. Only the drawing is absent.
    this.canvas = renderer ? renderer.domElement : document.createElement('canvas');
    if (!renderer) {
      const note = document.createElement('p');
      note.className = 'small muted';
      note.style.cssText = 'position:absolute;left:12px;right:12px;bottom:40px;text-align:center;pointer-events:none';
      note.textContent = 'WebGL is unavailable, so the head is not drawn. Tap inside this area to strike; the sound is still computed.';
      host.append(note);
    }
    this.canvas.setAttribute('aria-hidden', 'true');
    host.prepend(this.canvas);

    this.scene.add(new THREE.HemisphereLight(0xf4ead6, 0x0c0e12, 1.15));
    const key = new THREE.DirectionalLight(0xfff4e0, 2.2);
    key.position.set(-0.6, 1.4, 0.8);
    this.scene.add(key);
    const rim = new THREE.DirectionalLight(0x63d2c0, 0.35);
    rim.position.set(0.9, -0.3, -0.9);
    this.scene.add(rim);
    const under = new THREE.DirectionalLight(0xf1eadb, 0.6);
    under.position.set(0.2, -1.4, 0.5);
    this.scene.add(under);
    this.scene.add(this.drum);

    const ring = new THREE.RingGeometry(0.72, 1, 48);
    ring.rotateX(-Math.PI / 2);
    this.marker = new THREE.Mesh(ring, new THREE.MeshBasicMaterial({
      color: 0x63d2c0, transparent: true, opacity: 0, depthTest: false, side: THREE.DoubleSide,
    }));
    this.marker.renderOrder = 10;
    this.scene.add(this.marker);
    const dot = new THREE.RingGeometry(0.0, 1, 24);
    dot.rotateX(-Math.PI / 2);
    this.aim = new THREE.Mesh(dot, new THREE.MeshBasicMaterial({
      color: 0xf1eadb, transparent: true, opacity: 0.0, depthTest: false, side: THREE.DoubleSide,
    }));
    this.aim.renderOrder = 11;
    this.scene.add(this.aim);

    this.setGeometry(this.radius, this.depth);

    this.canvas.addEventListener('pointerdown', this.onPointerDown);
    this.canvas.addEventListener('pointermove', this.onPointerMove);
    this.canvas.addEventListener('pointerup', this.onPointerUp);
    this.canvas.addEventListener('pointercancel', this.onPointerCancel);
    this.canvas.addEventListener('contextmenu', (event) => event.preventDefault());
    document.addEventListener('visibilitychange', this.onVisibility);

    this.resizeObserver = new ResizeObserver(() => this.resize());
    this.resizeObserver.observe(host);
    this.intersection = new IntersectionObserver((entries) => {
      for (const entry of entries) this.visible = entry.isIntersecting;
      this.schedule();
    }, { threshold: 0.01 });
    this.intersection.observe(host);
    this.resize();
  }

  // ------------------------------------------------------------ geometry
  setGeometry(radius: number, depth: number): void {
    if (radius === this.radius && depth === this.depth && this.shell) return;
    this.radius = radius;
    this.depth = depth;
    const outer = radius + SHELL_OFFSET_M;
    for (const object of [this.shell, ...this.collars, ...this.hoops]) {
      if (!object) continue;
      this.drum.remove(object);
      object.geometry.dispose();
    }
    const shellMaterial = new THREE.MeshStandardMaterial({ color: 0x1c2029, metalness: 0.35, roughness: 0.42, side: THREE.DoubleSide });
    this.shell = new THREE.Mesh(new THREE.CylinderGeometry(outer, outer, depth, 96, 1, true), shellMaterial);
    this.drum.add(this.shell);
    const collarMaterial = new THREE.MeshStandardMaterial({ color: 0x6f6959, roughness: 0.8, side: THREE.DoubleSide });
    const hoopMaterial = new THREE.MeshStandardMaterial({ color: 0xd9a94f, metalness: 0.75, roughness: 0.32 });
    this.collars = [];
    this.hoops = [];
    for (const sign of [1, -1]) {
      const collar = new THREE.Mesh(new THREE.RingGeometry(radius, outer, 96, 1), collarMaterial);
      collar.rotation.x = -Math.PI / 2;
      collar.position.y = sign * depth / 2;
      this.collars.push(collar);
      this.drum.add(collar);
      const hoop = new THREE.Mesh(new THREE.TorusGeometry(outer, HOOP_TUBE_M, 10, 128), hoopMaterial);
      hoop.rotation.x = Math.PI / 2;
      hoop.position.y = sign * (depth / 2 + HOOP_TUBE_M * 0.6);
      this.hoops.push(hoop);
      this.drum.add(hoop);
    }
    if (!this.heads) this.buildPlaceholder();
    else for (const [i, head] of this.heads.entries()) head.position.y = i === 0 ? depth / 2 : -depth / 2;
    for (const [i, wire] of this.wires.entries()) wire.position.y = i === 0 ? depth / 2 : -depth / 2;
    this.dirty = true;
    this.schedule();
  }

  private buildPlaceholder(): void {
    for (const disc of this.placeholder) {
      this.drum.remove(disc);
      disc.geometry.dispose();
    }
    const material = new THREE.MeshStandardMaterial({ color: C_REST, roughness: 0.62, side: THREE.DoubleSide });
    this.placeholder = [1, -1].map((sign) => {
      const disc = new THREE.Mesh(new THREE.CircleGeometry(this.radius, 96), material);
      disc.rotation.x = -Math.PI / 2;
      disc.position.y = sign * this.depth / 2;
      this.drum.add(disc);
      return disc;
    });
  }

  /** Install the engine's actual head mesh (shared by both heads). */
  setMesh(nodes: Float32Array, triangles: Uint32Array): void {
    for (const disc of this.placeholder) {
      this.drum.remove(disc);
      disc.geometry.dispose();
    }
    this.placeholder = [];
    if (this.heads) for (const head of this.heads) { this.drum.remove(head); head.geometry.dispose(); }
    for (const wire of this.wires) { this.drum.remove(wire); wire.geometry.dispose(); }
    this.nodes = nodes;
    const count = nodes.length / 2;
    const edgeIndex = uniqueEdges(triangles);
    const material = new THREE.MeshStandardMaterial({
      vertexColors: true, roughness: 0.58, metalness: 0.0, side: THREE.DoubleSide,
      polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1,
    });
    // Thinner ink on denser meshes, so the film still reads as a film.
    const density = Math.min(1, Math.max(0, (triangles.length / 3 - 300) / 1500));
    const wireMaterial = new THREE.LineBasicMaterial({ color: 0x1b1712, transparent: true, opacity: 0.24 - 0.14 * density });
    const heads: THREE.Mesh[] = [];
    const colors: Float32Array[] = [];
    this.wires = [];
    for (const sign of [1, -1]) {
      const geometry = new THREE.BufferGeometry();
      const position = new THREE.BufferAttribute(new Float32Array(count * 3), 3);
      position.setUsage(THREE.DynamicDrawUsage);
      const color = new THREE.BufferAttribute(new Float32Array(count * 3), 3);
      color.setUsage(THREE.DynamicDrawUsage);
      for (let i = 0; i < count; i++) {
        position.setXYZ(i, nodes[2 * i] ?? 0, 0, -(nodes[2 * i + 1] ?? 0));
        color.setXYZ(i, C_REST.r, C_REST.g, C_REST.b);
      }
      geometry.setAttribute('position', position);
      geometry.setAttribute('color', color);
      geometry.setIndex(new THREE.BufferAttribute(triangles, 1));
      geometry.computeVertexNormals();
      geometry.computeBoundingSphere();
      const head = new THREE.Mesh(geometry, material);
      head.position.y = sign * this.depth / 2;
      heads.push(head);
      colors.push(color.array as Float32Array);
      this.drum.add(head);
      const wireGeometry = new THREE.BufferGeometry();
      wireGeometry.setAttribute('position', position);
      wireGeometry.setIndex(new THREE.BufferAttribute(edgeIndex, 1));
      const wire = new THREE.LineSegments(wireGeometry, wireMaterial);
      wire.position.y = sign * this.depth / 2;
      wire.renderOrder = 2;
      this.wires.push(wire);
      this.drum.add(wire);
    }
    this.heads = [heads[0]!, heads[1]!];
    this.colors = [colors[0]!, colors[1]!];
    this.displacement = null;
    this.dirty = true;
    this.schedule();
  }

  /** Latest accepted displacements [m], positive downward, one per mesh node. */
  setDisplacements(batter: Float32Array, resonant: Float32Array): void {
    if (!this.nodes || batter.length !== this.nodes.length / 2 || resonant.length !== batter.length) return;
    let peak = 0;
    for (let i = 0; i < batter.length; i++) {
      const a = Math.abs(batter[i] ?? 0);
      const b = Math.abs(resonant[i] ?? 0);
      if (a > peak) peak = a;
      if (b > peak) peak = b;
    }
    if (peak === 0 && this.lastPeak === 0 && this.displacement) return;
    this.lastPeak = peak;
    this.displacement = [batter, resonant];
    this.dirty = true;
    this.schedule();
  }

  /** Clear motion (a rebuild starts at rest). */
  rest(): void {
    this.displacement = null;
    this.envelope = FLOOR_M;
    this.shownFactor = -1;
    this.lastPeak = 0;
    this.dirty = true;
    this.schedule();
  }

  // ------------------------------------------------------------- markers
  showMarker(xM: number, yM: number): void {
    // Seen from below, mark the point the visitor touched on the resonant
    // head; the stick itself lands on the batter head above it.
    const below = this.viewTarget === 1;
    this.marker.position.set(xM, below ? -this.depth / 2 - 0.002 : this.depth / 2 + 0.002, -yM);
    this.markerBorn = performance.now();
    this.dirty = true;
    this.schedule();
  }

  setAim(xM: number, yM: number, visible: boolean): void {
    const below = this.viewTarget === 1;
    this.aim.position.set(xM, below ? -this.depth / 2 - 0.0025 : this.depth / 2 + 0.0025, -yM);
    (this.aim.material as THREE.MeshBasicMaterial).opacity = visible ? 0.85 : 0;
    this.dirty = true;
    this.schedule();
  }

  setView(view: StageView): void {
    this.view = view;
    this.viewTarget = view === 'resonant' ? 1 : 0;
    const plane = this.viewTarget === 1 ? -this.depth / 2 - 0.0025 : this.depth / 2 + 0.0025;
    this.aim.position.y = plane;
    this.marker.position.y = plane;
    if (this.reducedMotion) this.viewBlend = this.viewTarget;
    this.dirty = true;
    this.schedule();
  }

  setMode(mode: StageMode): void {
    this.mode = mode;
    this.host.dataset.mode = mode;
  }

  setStarted(started: boolean): void {
    this.started = started;
    this.dirty = true;
    this.schedule();
  }

  get exaggeration(): number {
    return this.factor;
  }

  // -------------------------------------------------------------- input
  /** Head-plane metres under a client point, or null when not on the drumhead. */
  pick(clientX: number, clientY: number): [number, number] | 'rim' | null {
    const rect = this.canvas.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) return null;
    const ndc = new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
    const ray = new THREE.Raycaster();
    ray.setFromCamera(ndc, this.camera);
    // Whichever head faces the camera; a tap on the resonant head (seen from
    // below) strikes the batter head at the same x, y: the stick is on top.
    const below = this.camera.position.y < 0;
    const plane = new THREE.Plane(new THREE.Vector3(0, below ? -1 : 1, 0), below ? -this.depth / 2 : -this.depth / 2);
    const hit = new THREE.Vector3();
    if (!ray.ray.intersectPlane(plane, hit)) return null;
    const x = hit.x;
    const y = -hit.z;
    const r = Math.hypot(x, y);
    if (r > this.radius + SHELL_OFFSET_M + HOOP_TUBE_M) return null;
    if (r > this.radius * 0.985) return 'rim';
    return [x, y];
  }

  private readonly onPointerDown = (event: PointerEvent): void => {
    if (!this.started) return;
    this.canvas.setPointerCapture(event.pointerId);
    this.pointers.set(event.pointerId, {
      x: event.clientX, y: event.clientY, startX: event.clientX, startY: event.clientY, t0: performance.now(), type: event.pointerType,
    });
    if (this.pointers.size > 1) {
      this.gesture = 'orbit2';
      this.pendingHit = null;
      return;
    }
    const hit = this.pick(event.clientX, event.clientY);
    if (hit === 'rim') {
      this.cb.hint('Rim shots are not modelled: strike inside the head');
      this.gesture = 'pending';
      this.pendingHit = null;
      return;
    }
    if (hit && this.mode === 'strike') {
      // Strike on press: zero added latency. A drag that begins on the head
      // never orbits, so a stroke is never mistaken for a camera move.
      this.cb.strike(hit[0], hit[1], 'pointer');
      this.gesture = 'strike';
      return;
    }
    this.gesture = 'pending';
    this.pendingHit = hit;
  };

  private readonly onPointerMove = (event: PointerEvent): void => {
    const track = this.pointers.get(event.pointerId);
    if (!track) {
      if (event.pointerType === 'mouse' && this.started) {
        const hit = this.pick(event.clientX, event.clientY);
        this.canvas.style.cursor = Array.isArray(hit) ? 'crosshair' : 'grab';
      }
      return;
    }
    const dx = event.clientX - track.x;
    const dy = event.clientY - track.y;
    track.x = event.clientX;
    track.y = event.clientY;
    if (this.gesture === 'pending' && Math.hypot(event.clientX - track.startX, event.clientY - track.startY) > TAP_SLOP_PX) {
      this.gesture = 'orbit';
      this.pendingHit = null;
      this.canvas.style.cursor = 'grabbing';
    }
    if (this.gesture === 'orbit') this.orbit(dx, dy);
    else if (this.gesture === 'orbit2') this.orbit(dx / this.pointers.size, dy / this.pointers.size);
  };

  private readonly onPointerUp = (event: PointerEvent): void => {
    const track = this.pointers.get(event.pointerId);
    this.pointers.delete(event.pointerId);
    if (!track) return;
    if (this.gesture === 'pending' && this.pendingHit && performance.now() - track.t0 <= TAP_MS) {
      this.cb.strike(this.pendingHit[0], this.pendingHit[1], 'pointer');
    }
    if (this.pointers.size === 0) {
      this.gesture = 'none';
      this.pendingHit = null;
      if (event.pointerType === 'mouse') this.canvas.style.cursor = 'grab';
    }
  };

  private readonly onPointerCancel = (event: PointerEvent): void => {
    this.pointers.delete(event.pointerId);
    if (this.pointers.size === 0) {
      this.gesture = 'none';
      this.pendingHit = null;
    }
  };

  private orbit(dx: number, dy: number): void {
    const rect = this.canvas.getBoundingClientRect();
    const scale = Math.PI / Math.max(240, rect.width);
    this.yaw -= dx * scale;
    this.pitch = THREE.MathUtils.clamp(this.pitch + dy * scale, THREE.MathUtils.degToRad(8), THREE.MathUtils.degToRad(84));
    this.dirty = true;
    this.schedule();
  }

  // ------------------------------------------------------------ the loop
  private resize(): void {
    const width = Math.max(1, this.host.clientWidth);
    const height = Math.max(1, this.host.clientHeight);
    if (this.renderer) {
      this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
      this.renderer.setSize(width, height, false);
    } else {
      this.canvas.style.width = `${width}px`;
      this.canvas.style.height = `${height}px`;
    }
    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.dirty = true;
    this.schedule();
  }

  private schedule(): void {
    if (this.disposed || this.raf !== 0) return;
    if (!this.visible || !this.pageVisible) return;
    this.raf = requestAnimationFrame(this.frame);
  }

  private readonly frame = (now: number): void => {
    this.raf = 0;
    if (this.disposed || !this.visible || !this.pageVisible) return;
    const dt = this.lastTime > 0 ? Math.min(0.1, (now - this.lastTime) / 1000) : 1 / 60;
    this.lastTime = now;
    let animating = false;

    if (!this.started && !this.reducedMotion) {
      this.yaw += dt * 0.12;
      this.dirty = true;
      animating = true;
    }
    if (this.viewBlend !== this.viewTarget) {
      const step = dt / 0.6;
      this.viewBlend = this.viewTarget > this.viewBlend
        ? Math.min(this.viewTarget, this.viewBlend + step)
        : Math.max(this.viewTarget, this.viewBlend - step);
      this.dirty = true;
      animating = true;
    }
    const age = (now - this.markerBorn) / 1000;
    const markerMaterial = this.marker.material as THREE.MeshBasicMaterial;
    if (age < 0.7) {
      const t = age / 0.7;
      markerMaterial.opacity = 0.95 * (1 - t);
      const size = 0.012 + (this.reducedMotion ? 0 : 0.02 * t);
      this.marker.scale.setScalar(size);
      this.dirty = true;
      animating = true;
    } else if (markerMaterial.opacity !== 0) {
      markerMaterial.opacity = 0;
      this.dirty = true;
    }
    this.aim.scale.setScalar(0.0045);

    if (this.dirty) {
      this.updateHeads(dt);
      this.placeCamera();
      this.renderer?.render(this.scene, this.camera);
      this.dirty = false;
    }
    if (animating || this.lastPeak > 0) this.schedule();
  };

  private placeCamera(): void {
    const height = this.host.clientHeight || 1;
    const width = this.host.clientWidth || 1;
    // Fixed framing (sized for a 0.25 m head, 0.3 m deep), so a radius or
    // depth change is seen as a change of the drum, not of the camera.
    // Narrow (phone) stages also carry the controls and HUD: frame wider.
    const extent = 0.25 * Math.max(1, 1.25 / (width / height)) ** 0.6;
    const vertical = THREE.MathUtils.degToRad(this.camera.fov) / 2;
    const horizontal = Math.atan(Math.tan(vertical) * (width / height));
    const distance = extent / Math.sin(Math.min(vertical, horizontal));
    const ease = this.viewBlend * this.viewBlend * (3 - 2 * this.viewBlend);
    const pitch = THREE.MathUtils.lerp(this.pitch, -this.pitch, ease);
    this.camera.position.set(
      distance * Math.cos(pitch) * Math.sin(this.yaw),
      distance * Math.sin(pitch),
      distance * Math.cos(pitch) * Math.cos(this.yaw),
    );
    this.camera.up.set(0, 1, 0);
    this.camera.lookAt(0, 0, 0);
  }

  private updateHeads(dt: number): void {
    if (!this.heads || !this.colors || !this.nodes) return;
    const displacement = this.displacement;
    const peak = this.lastPeak;
    this.envelope = Math.max(FLOOR_M, peak, this.envelope * Math.exp(-dt / RELEASE_S));
    const wanted = niceBelow((SHOWN_FRACTION * this.radius) / this.envelope);
    const now = performance.now();
    if (wanted < this.factor || (wanted > this.factor && now - this.factorChangedAt > 400)) {
      if (wanted !== this.factor) {
        this.factor = wanted;
        this.factorChangedAt = now;
      }
    }
    const shown = this.envelope > FLOOR_M * 1.01 ? this.factor : 0;
    if (shown !== this.shownFactor) {
      this.shownFactor = shown;
      this.cb.scale(shown);
    }
    const k = this.factor;
    const count = this.nodes.length / 2;
    for (let h = 0; h < 2; h++) {
      const head = this.heads[h]!;
      const position = head.geometry.getAttribute('position') as THREE.BufferAttribute;
      const colors = this.colors[h]!;
      const w = displacement ? displacement[h]! : null;
      for (let i = 0; i < count; i++) {
        const value = w ? (w[i] ?? 0) : 0;
        // Positive engine displacement is DOWNWARD for both heads.
        position.setY(i, -value * k);
        const s = Math.max(-1, Math.min(1, value / this.envelope));
        const target = s < 0 ? C_UP : C_DOWN;
        const a = Math.abs(s);
        colors[3 * i] = C_REST.r + (target.r - C_REST.r) * a;
        colors[3 * i + 1] = C_REST.g + (target.g - C_REST.g) * a;
        colors[3 * i + 2] = C_REST.b + (target.b - C_REST.b) * a;
      }
      position.needsUpdate = true;
      (head.geometry.getAttribute('color') as THREE.BufferAttribute).needsUpdate = true;
      head.geometry.computeVertexNormals();
    }
  }

  dispose(): void {
    this.disposed = true;
    if (this.raf) cancelAnimationFrame(this.raf);
    this.raf = 0;
    this.resizeObserver.disconnect();
    this.intersection.disconnect();
    document.removeEventListener('visibilitychange', this.onVisibility);
    this.scene.traverse((object) => {
      const mesh = object as THREE.Mesh;
      if (mesh.geometry) mesh.geometry.dispose();
    });
    this.renderer?.dispose();
  }
}

function uniqueEdges(triangles: Uint32Array): Uint32Array {
  const seen = new Set<number>();
  const edges: number[] = [];
  for (let t = 0; t + 2 < triangles.length; t += 3) {
    const a = triangles[t] ?? 0;
    const b = triangles[t + 1] ?? 0;
    const c = triangles[t + 2] ?? 0;
    for (const [p, q] of [[a, b], [b, c], [c, a]] as const) {
      const lo = Math.min(p, q);
      const hi = Math.max(p, q);
      const key = lo * 1_048_576 + hi;
      if (!seen.has(key)) {
        seen.add(key);
        edges.push(lo, hi);
      }
    }
  }
  return new Uint32Array(edges);
}
