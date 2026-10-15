type Child = Node | string | null | undefined | false;

export function el<K extends keyof HTMLElementTagNameMap>(tag: K, attrs: Record<string, string> = {}, ...children: Child[]): HTMLElementTagNameMap[K] {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") node.className = v;
    else node.setAttribute(k, v);
  }
  for (const c of children) if (c !== null && c !== undefined && c !== false) node.append(c);
  return node;
}

export function fmt(v: number, digits = 3): string {
  if (!Number.isFinite(v)) return "—";
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e7)) return fmtSci(v);
  if (a >= 1000) return Math.round(v).toLocaleString("en-US");
  return Number(v.toPrecision(digits)).toString();
}

export function fmtSci(v: number, digits = 2): string {
  if (!Number.isFinite(v)) return "—";
  if (v === 0) return "0";
  const [m, e] = v.toExponential(digits - 1).split("e");
  const exp = Number(e);
  const sup = String(exp).replace(/-/g, "⁻").replace(/\d/g, (d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(d)]!);
  return `${m}×10${sup}`;
}

const DPR_CAP = 2;
export function dpr(): number {
  return Math.min(window.devicePixelRatio || 1, DPR_CAP);
}

/** Size a canvas backing store to its CSS width at a fixed aspect ratio. */
export function hidpi(canvas: HTMLCanvasElement, aspect: number, minCssHeight = 0): void {
  const w = canvas.clientWidth || canvas.parentElement?.clientWidth || 600;
  const cssH = Math.max(minCssHeight, Math.round(w * aspect));
  canvas.style.height = `${cssH}px`;
  const d = dpr();
  const bw = Math.round(w * d);
  const bh = Math.round(cssH * d);
  if (canvas.width !== bw) canvas.width = bw;
  if (canvas.height !== bh) canvas.height = bh;
}

export function onVisible(target: Element, cb: () => void, rootMargin = "300px"): void {
  const io = new IntersectionObserver((entries) => {
    if (entries.some((e) => e.isIntersecting)) {
      io.disconnect();
      cb();
    }
  }, { rootMargin });
  io.observe(target);
}

export const reducedMotion = (): boolean => window.matchMedia("(prefers-reduced-motion: reduce)").matches;

export const NOTE_NAMES = ["C", "C♯", "D", "D♯", "E", "F", "F♯", "G", "G♯", "A", "A♯", "B"];
export function noteName(midi: number): string {
  return `${NOTE_NAMES[midi % 12]}${Math.floor(midi / 12) - 1}`;
}

export async function json<T>(url: string): Promise<T> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`${url}: HTTP ${res.status}`);
  return (await res.json()) as T;
}

/** Scientific notation as DOM with a real superscript (monospace fonts lack superscript glyphs). */
export function sciNode(v: number, digits = 2): HTMLElement {
  const span = document.createElement("span");
  if (!Number.isFinite(v) || v === 0) {
    span.textContent = Number.isFinite(v) ? "0" : "—";
    return span;
  }
  const [m, e] = v.toExponential(digits - 1).split("e");
  const sup = document.createElement("sup");
  sup.className = "e";
  sup.textContent = String(Number(e)).replace("-", "−");
  span.append(`${m}×10`, sup);
  return span;
}

/** The deployment's base path ("/" at the domain root, "/music/" under frankensim.org/music). */
export const BASE = import.meta.env.BASE_URL;

/** Content hash of each published data, audio, image, replay and MIDI file (a build sets it; vite.config.ts). */
const HASHES: Record<string, string> = import.meta.env.VITE_ASSET_HASHES ?? {};

/**
 * A URL for a file published with the site (data, audio, images, live engines). In a build it carries the
 * file's content hash (`?v=`), so a file that changed in a deploy is never served from an older cache.
 */
export function asset(path: string): string {
  const p = path.replace(/^\/+/, "");
  const v = HASHES[p];
  return BASE + p + (v ? `?v=${v}` : "");
}
