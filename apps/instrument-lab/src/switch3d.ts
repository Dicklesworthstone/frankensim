import * as THREE from "three";
import { RoundedBoxGeometry } from "three/examples/jsm/geometries/RoundedBoxGeometry.js";
import { RoomEnvironment } from "three/examples/jsm/environments/RoomEnvironment.js";
import { reducedMotion } from "./util";

export type SwitchState = "studio" | "live";

const UP = -0.52; // lever tilt (rad) toward the top legend
const DOWN = 0.52;

/** A brushed-aluminium panel with a chrome bat-handle toggle and two lamps. */
export class Switch3D {
  onChange: ((s: SwitchState) => void) | null = null;
  private renderer: THREE.WebGLRenderer;
  private scene = new THREE.Scene();
  private camera = new THREE.PerspectiveCamera(26, 1, 0.1, 100);
  private lever = new THREE.Group();
  private lamps: { mat: THREE.MeshStandardMaterial; glow: THREE.Sprite; light: THREE.PointLight; color: THREE.Color; level: number }[] = [];
  private angle = UP;
  private vel = 0;
  private target = UP;
  private state: SwitchState = "studio";
  private raf = 0;
  private last = 0;
  private drag: { y0: number; a0: number; moved: boolean } | null = null;
  private ro: ResizeObserver;
  private disposed = false;

  constructor(private host: HTMLElement) {
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "low-power" });
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 0.82;
    this.renderer.shadowMap.enabled = true;
    this.renderer.shadowMap.type = THREE.PCFShadowMap;
    const cv = this.renderer.domElement;
    cv.style.cssText = "display:block;width:100%;height:100%;touch-action:none;cursor:grab";
    host.append(cv);
    const pm = new THREE.PMREMGenerator(this.renderer);
    this.scene.environment = pm.fromScene(new RoomEnvironment(), 0.035).texture;
    this.scene.environmentIntensity = 0.62;
    pm.dispose();
    this.build();
    this.camera.position.set(2.4, -6.6, 19.5);
    this.camera.lookAt(0.15, -0.35, 0.9);
    this.ro = new ResizeObserver(() => this.resize());
    this.ro.observe(host);
    this.resize();
    cv.addEventListener("pointerdown", this.down);
    window.addEventListener("pointermove", this.move);
    window.addEventListener("pointerup", this.up);
    window.addEventListener("pointercancel", this.up);
    this.kick();
  }

  get value(): SwitchState {
    return this.state;
  }

  set(s: SwitchState, animate = true): void {
    this.state = s;
    this.target = s === "studio" ? UP : DOWN;
    if (!animate || reducedMotion()) {
      this.angle = this.target;
      this.vel = 0;
    }
    this.kick();
  }

  toggle(): void {
    this.set(this.state === "studio" ? "live" : "studio");
    this.onChange?.(this.state);
  }

  dispose(): void {
    this.disposed = true;
    cancelAnimationFrame(this.raf);
    this.ro.disconnect();
    window.removeEventListener("pointermove", this.move);
    window.removeEventListener("pointerup", this.up);
    window.removeEventListener("pointercancel", this.up);
    this.renderer.dispose();
  }

  // ------------------------------------------------------------------ scene
  private build(): void {
    const s = this.scene;
    const key = new THREE.DirectionalLight(0xfff3e0, 2.2);
    key.position.set(-5, 7, 9);
    key.castShadow = true;
    key.shadow.mapSize.set(1024, 1024);
    key.shadow.camera.left = -5; key.shadow.camera.right = 5; key.shadow.camera.top = 6; key.shadow.camera.bottom = -6;
    key.shadow.radius = 5;
    key.shadow.bias = -0.0004;
    s.add(key);
    const rim = new THREE.DirectionalLight(0x9fd8ff, 0.6);
    rim.position.set(6, -3, 4);
    s.add(rim);

    // the panel: brushed aluminium with engraved legends
    const W = 5.6, H = 8.2, T = 0.5;
    const brushed = brushedTexture();
    const plateMat = new THREE.MeshPhysicalMaterial({ color: 0x8e949c, metalness: 1, roughness: 0.5, roughnessMap: brushed, clearcoat: 0.25, clearcoatRoughness: 0.6 });
    const plate = new THREE.Mesh(new RoundedBoxGeometry(W, H, T, 6, 0.38), plateMat);
    plate.receiveShadow = true;
    s.add(plate);
    const legend = new THREE.Mesh(new THREE.PlaneGeometry(W - 0.5, H - 0.5), new THREE.MeshPhysicalMaterial({ map: legendTexture(), transparent: true, metalness: 0.6, roughness: 0.55, depthWrite: false }));
    legend.position.z = T / 2 + 0.002;
    legend.receiveShadow = true;
    s.add(legend);
    // a darker recessed bed under the panel
    const bed = new THREE.Mesh(new RoundedBoxGeometry(W + 0.5, H + 0.5, 0.3, 4, 0.5), new THREE.MeshStandardMaterial({ color: 0x15171c, metalness: 0.4, roughness: 0.7 }));
    bed.position.z = -0.35;
    bed.receiveShadow = true;
    s.add(bed);

    const chrome = new THREE.MeshPhysicalMaterial({ color: 0xf2f4f7, metalness: 1, roughness: 0.08, clearcoat: 1, clearcoatRoughness: 0.05 });
    const steel = new THREE.MeshPhysicalMaterial({ color: 0x9aa0a8, metalness: 1, roughness: 0.3 });
    // slotted screws in the corners
    for (const [x, y] of [[-1, -1], [1, -1], [-1, 1], [1, 1]] as const) {
      const g = new THREE.Group();
      const head = new THREE.Mesh(new THREE.SphereGeometry(0.24, 24, 12, 0, Math.PI * 2, 0, Math.PI / 2), steel);
      head.rotation.x = Math.PI / 2;
      head.scale.set(1, 0.45, 1);
      const slot = new THREE.Mesh(new THREE.BoxGeometry(0.4, 0.05, 0.12), new THREE.MeshStandardMaterial({ color: 0x23262c, roughness: 0.8 }));
      slot.position.z = 0.09;
      slot.rotation.z = (x * 0.7 + y * 0.35);
      g.add(head, slot);
      g.position.set(x * (W / 2 - 0.45), y * (H / 2 - 0.45), T / 2);
      head.castShadow = true;
      s.add(g);
    }
    // collar: washer, hex nut, threaded boss
    const washer = new THREE.Mesh(new THREE.CylinderGeometry(0.95, 0.95, 0.08, 48), steel);
    washer.rotation.x = Math.PI / 2;
    washer.position.z = T / 2 + 0.04;
    const nut = new THREE.Mesh(new THREE.CylinderGeometry(0.72, 0.72, 0.34, 6), chrome);
    nut.rotation.x = Math.PI / 2;
    nut.rotation.y = Math.PI / 6;
    nut.position.z = T / 2 + 0.25;
    const boss = new THREE.Mesh(new THREE.CylinderGeometry(0.46, 0.5, 0.42, 40), steel);
    boss.rotation.x = Math.PI / 2;
    boss.position.z = T / 2 + 0.55;
    for (const m of [washer, nut, boss]) { m.castShadow = true; m.receiveShadow = true; s.add(m); }

    // the bat handle, pivoting just inside the boss
    this.lever.position.z = T / 2 + 0.62;
    const dome = new THREE.Mesh(new THREE.SphereGeometry(0.4, 32, 16, 0, Math.PI * 2, 0, Math.PI / 2), chrome);
    dome.rotation.x = Math.PI / 2;
    const shaft = new THREE.Mesh(new THREE.CylinderGeometry(0.19, 0.3, 3.6, 40), chrome);
    shaft.rotation.x = Math.PI / 2;
    shaft.position.z = 1.85;
    const tip = new THREE.Mesh(new THREE.SphereGeometry(0.3, 32, 24), chrome);
    tip.scale.set(1, 1, 1.25);
    tip.position.z = 3.7;
    for (const m of [dome, shaft, tip]) { m.castShadow = true; this.lever.add(m); }
    s.add(this.lever);

    // lamps: brass for the studio render (top), verdigris for live (bottom)
    const glowTex = glowTexture();
    for (const [y, hex] of [[H / 2 - 1.05, 0xffb547], [-H / 2 + 1.05, 0x5cf0d6]] as const) {
      const color = new THREE.Color(hex);
      const bezel = new THREE.Mesh(new THREE.TorusGeometry(0.34, 0.07, 16, 40), chrome);
      bezel.position.set(-1.75, y, T / 2 + 0.05);
      const mat = new THREE.MeshStandardMaterial({ color: color.clone().multiplyScalar(0.25), emissive: color, emissiveIntensity: 0, roughness: 0.25, metalness: 0 });
      const lens = new THREE.Mesh(new THREE.SphereGeometry(0.3, 32, 16, 0, Math.PI * 2, 0, Math.PI / 2), mat);
      lens.rotation.x = Math.PI / 2;
      lens.scale.set(1, 0.6, 1);
      lens.position.set(-1.75, y, T / 2 + 0.02);
      const glow = new THREE.Sprite(new THREE.SpriteMaterial({ map: glowTex, color, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false }));
      glow.scale.set(2.6, 2.6, 1);
      glow.position.set(-1.75, y, T / 2 + 0.3);
      const light = new THREE.PointLight(color, 0, 3.5, 2);
      light.position.set(-1.75, y, T / 2 + 0.5);
      s.add(bezel, lens, glow, light);
      this.lamps.push({ mat, glow, light, color, level: 0 });
    }
  }

  private resize(): void {
    const r = this.host.getBoundingClientRect();
    const w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    // keep the whole panel in view whatever the aspect
    const k = w / h < 0.72 ? 0.72 / Math.max(0.45, w / h) : 1;
    this.camera.position.set(2.4 * k, -6.6 * k, 19.5 * k);
    this.camera.lookAt(0.15, -0.35, 0.9);
    this.camera.updateProjectionMatrix();
    this.kick();
  }

  // --------------------------------------------------------------- pointer
  private down = (e: PointerEvent): void => {
    this.drag = { y0: e.clientY, a0: this.angle, moved: false };
    this.renderer.domElement.style.cursor = "grabbing";
    (e.target as Element).setPointerCapture?.(e.pointerId);
    e.preventDefault();
  };
  private move = (e: PointerEvent): void => {
    if (!this.drag) return;
    const dy = e.clientY - this.drag.y0;
    if (Math.abs(dy) > 5) this.drag.moved = true;
    if (!this.drag.moved) return;
    const span = this.host.getBoundingClientRect().height * 0.28;
    this.angle = Math.max(UP - 0.06, Math.min(DOWN + 0.06, this.drag.a0 + (dy / span) * (DOWN - UP)));
    this.target = this.angle;
    this.vel = 0;
    this.kick();
  };
  private up = (): void => {
    if (!this.drag) return;
    const d = this.drag;
    this.drag = null;
    this.renderer.domElement.style.cursor = "grab";
    const next: SwitchState = d.moved ? (this.angle < 0 ? "studio" : "live") : (this.state === "studio" ? "live" : "studio");
    const changed = next !== this.state;
    this.set(next);
    if (changed) {
      if ("vibrate" in navigator) try { navigator.vibrate(12); } catch { /* not allowed */ }
      this.onChange?.(next);
    }
  };

  // --------------------------------------------------------------- animation
  private kick(): void {
    if (this.raf || this.disposed) return;
    this.last = performance.now();
    this.raf = requestAnimationFrame(this.frame);
  }
  private frame = (now: number): void => {
    this.raf = 0;
    const dt = Math.min(0.033, (now - this.last) / 1000);
    this.last = now;
    if (!this.drag) {
      // a stiff detent spring with a little overshoot: the bat snaps over centre
      const k = 520, c = 24;
      const steps = 4;
      for (let i = 0; i < steps; i++) {
        const acc = -k * (this.angle - this.target) - c * this.vel;
        this.vel += (acc * dt) / steps;
        this.angle += (this.vel * dt) / steps;
      }
    }
    this.lever.rotation.x = this.angle;
    let busy = this.drag !== null || Math.abs(this.angle - this.target) > 1e-4 || Math.abs(this.vel) > 1e-3;
    const on = [this.state === "studio", this.state === "live"];
    this.lamps.forEach((l, i) => {
      const goal = on[i] ? 1 : 0;
      l.level += (goal - l.level) * Math.min(1, dt * (goal ? 14 : 7));
      if (Math.abs(goal - l.level) > 1e-3) busy = true;
      l.mat.emissiveIntensity = 0.04 + 1.25 * l.level;
      (l.glow.material as THREE.SpriteMaterial).opacity = 0.7 * l.level;
      l.light.intensity = 1.6 * l.level;
    });
    this.renderer.render(this.scene, this.camera);
    if (busy) this.raf = requestAnimationFrame(this.frame);
  };
}

function brushedTexture(): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = 256; c.height = 512;
  const g = c.getContext("2d")!;
  g.fillStyle = "#8a8a8a";
  g.fillRect(0, 0, c.width, c.height);
  let seed = 7;
  const rnd = (): number => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  for (let i = 0; i < 2600; i++) {
    const y = rnd() * c.height, v = 110 + rnd() * 70 | 0;
    g.strokeStyle = `rgba(${v},${v},${v},${0.25 + rnd() * 0.35})`;
    g.lineWidth = 0.6 + rnd();
    g.beginPath();
    const x0 = rnd() * c.width;
    g.moveTo(x0 - 40 - rnd() * 120, y);
    g.lineTo(x0 + 40 + rnd() * 120, y);
    g.stroke();
  }
  const t = new THREE.CanvasTexture(c);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(1, 1);
  return t;
}

function legendTexture(): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = 1020; c.height = 1540;
  const g = c.getContext("2d")!;
  const engrave = (text: string, x: number, y: number, font: string, spacing: number): void => {
    g.font = font;
    g.textAlign = "left";
    g.textBaseline = "middle";
    const chars = [...text];
    const width = chars.reduce((w, ch) => w + g.measureText(ch).width, 0) + spacing * (chars.length - 1);
    let cx = x - width / 2;
    for (const ch of chars) {
      g.fillStyle = "rgba(255,255,255,0.32)";
      g.fillText(ch, cx, y + 2.2);
      g.fillStyle = "rgba(18,20,24,0.88)";
      g.fillText(ch, cx, y);
      cx += g.measureText(ch).width + spacing;
    }
  };
  const mid = c.width / 2 + 70;
  engrave("STUDIO RENDER", mid, 170, "700 60px 'Inter Variable', Inter, system-ui, sans-serif", 7);
  engrave("OFFLINE · FULL RESOLUTION", mid, 245, "500 32px 'JetBrains Mono', ui-monospace, monospace", 4);
  engrave("LIVE IN THIS TAB", mid, c.height - 245, "700 60px 'Inter Variable', Inter, system-ui, sans-serif", 7);
  engrave("WEBASSEMBLY · REAL TIME", mid, c.height - 170, "500 32px 'JetBrains Mono', ui-monospace, monospace", 4);
  // travel marks either side of the collar
  g.strokeStyle = "rgba(18,20,24,0.8)";
  g.lineWidth = 6;
  for (const sgn of [-1, 1]) {
    const y = c.height / 2 + sgn * 300;
    g.beginPath();
    g.moveTo(c.width / 2 - 40, y + sgn * -24);
    g.lineTo(c.width / 2, y);
    g.lineTo(c.width / 2 + 40, y + sgn * -24);
    g.stroke();
  }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = 4;
  return t;
}

function glowTexture(): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = c.height = 128;
  const g = c.getContext("2d")!;
  const r = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  r.addColorStop(0, "rgba(255,255,255,1)");
  r.addColorStop(0.18, "rgba(255,255,255,0.55)");
  r.addColorStop(0.5, "rgba(255,255,255,0.12)");
  r.addColorStop(1, "rgba(255,255,255,0)");
  g.fillStyle = r;
  g.fillRect(0, 0, 128, 128);
  return new THREE.CanvasTexture(c);
}
