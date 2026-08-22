// Analog HUD dial panel (UI overhaul, root guzez). THIN DOM/SVG writer
// over the tested gauge math in gauges.ts — all angles, pinning, and
// classification come from there; this file only draws. Brass-and-cream
// period aesthetic, walnut panel, red danger arcs.

import {
  IDLE_INPUTS,
  type DialView,
  type GaugeSpec,
  type HudDialInputs,
  type LeverView,
  dialSetFrom,
  leverSetFrom,
  tickMarks,
} from "./gauges.ts";

const SVG_NS = "http://www.w3.org/2000/svg";
const FACE = "#f3e9d2";
const BRASS = "#b08d3f";
const BRASS_DARK = "#7a5f24";
const INK = "#241d12";
const RED = "#a4321f";

function polar(cx: number, cy: number, r: number, deg: number): [number, number] {
  const rad = ((deg - 90) * Math.PI) / 180;
  return [cx + r * Math.cos(rad), cy + r * Math.sin(rad)];
}

function arcPath(cx: number, cy: number, r: number, fromDeg: number, toDeg: number): string {
  const [x0, y0] = polar(cx, cy, r, fromDeg);
  const [x1, y1] = polar(cx, cy, r, toDeg);
  const large = toDeg - fromDeg > 180 ? 1 : 0;
  return `M ${x0.toFixed(2)} ${y0.toFixed(2)} A ${r} ${r} 0 ${large} 1 ${x1.toFixed(2)} ${y1.toFixed(2)}`;
}

function el<K extends keyof SVGElementTagNameMap>(
  tag: K,
  attrs: Record<string, string>,
): SVGElementTagNameMap[K] {
  const node = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) {
    node.setAttribute(k, v);
  }
  return node;
}

interface DialDom {
  needle: SVGLineElement;
  reading: SVGTextElement;
  face: SVGCircleElement;
  /** Pivot coordinate (the rotate() origin) — owned by buildDial so
   * update() never re-derives it from the dial's id. */
  center: number;
}

/** Needle angle for a VALUE position on the face (arc endpoints). */
function valueDeg(value: number, spec: GaugeSpec): number {
  return spec.startDeg + ((value - spec.min) / (spec.max - spec.min)) * spec.sweepDeg;
}

function buildDial(view: DialView, size: number): { svg: SVGSVGElement; dom: DialDom } {
  const c = size / 2;
  const svg = el("svg", {
    viewBox: `0 0 ${size} ${size}`,
    width: String(size),
    height: String(size),
  });
  // Brass bezel + face.
  svg.appendChild(el("circle", { cx: `${c}`, cy: `${c}`, r: `${c - 1}`, fill: BRASS_DARK }));
  svg.appendChild(el("circle", { cx: `${c}`, cy: `${c}`, r: `${c - 3}`, fill: BRASS }));
  const face = el("circle", { cx: `${c}`, cy: `${c}`, r: `${c - 7}`, fill: FACE });
  svg.appendChild(face);
  // Danger arc(s). A symmetric dial (min = −max) whose danger is felt
  // in MAGNITUDE (the pitch dial) gets the mirrored arc too, so the
  // face marks nose-down danger exactly as loudly as nose-up.
  if (view.spec.redline !== null) {
    const [from, to] = view.spec.redline;
    const drawArc = (loVal: number, hiVal: number): void => {
      svg.appendChild(
        el("path", {
          d: arcPath(c, c, c - 12, valueDeg(loVal, view.spec), valueDeg(hiVal, view.spec)),
          stroke: RED,
          "stroke-width": "4",
          fill: "none",
        }),
      );
    };
    drawArc(from, to);
    if (view.spec.min === -view.spec.max && from > 0) {
      drawArc(-to, -from);
    }
  }
  // Ticks + numerals. A full-circle face (360° sweep) would draw its
  // last tick ON its first — skip the duplicate.
  const fullCircle = view.spec.sweepDeg === 360;
  const allTicks = tickMarks(view.spec, view.majorTicks);
  const drawnTicks = fullCircle ? allTicks.slice(0, -1) : allTicks;
  for (const [i, t] of drawnTicks.entries()) {
    const [x0, y0] = polar(c, c, c - 9, t.deg);
    const [x1, y1] = polar(c, c, c - 16, t.deg);
    svg.appendChild(
      el("line", {
        x1: `${x0.toFixed(2)}`,
        y1: `${y0.toFixed(2)}`,
        x2: `${x1.toFixed(2)}`,
        y2: `${y1.toFixed(2)}`,
        stroke: INK,
        "stroke-width": "1.6",
      }),
    );
    // Label every other tick on dense faces.
    if (view.majorTicks <= 9 || i % 2 === 0) {
      const [lx, ly] = polar(c, c, c - 24, t.deg);
      const label = el("text", {
        x: `${lx.toFixed(2)}`,
        y: `${(ly + 3).toFixed(2)}`,
        "text-anchor": "middle",
        "font-size": "9",
        fill: INK,
        "font-family": "Georgia, serif",
      });
      label.textContent = String(Math.round(t.value / view.labelDiv));
      svg.appendChild(label);
    }
  }
  // Instrument name + provenance stamp.
  const name = el("text", {
    x: `${c}`,
    y: `${c - 14}`,
    "text-anchor": "middle",
    "font-size": "8.5",
    fill: INK,
    "font-family": "Georgia, serif",
    "letter-spacing": "1",
  });
  name.textContent = view.label;
  svg.appendChild(name);
  const stamp = el("text", {
    x: `${c}`,
    y: `${size - 14}`,
    "text-anchor": "middle",
    "font-size": "6.5",
    fill: view.provenance === "MODERN" ? BRASS_DARK : "#6b6152",
    "font-family": "Georgia, serif",
    "letter-spacing": "1.5",
  });
  stamp.textContent = view.provenance === "MODERN" ? "MODERN" : "· 1903 ·";
  svg.appendChild(stamp);
  // Needle + hub + reading.
  const needle = el("line", {
    x1: `${c}`,
    y1: `${c + 8}`,
    x2: `${c}`,
    y2: `${14}`,
    stroke: INK,
    "stroke-width": "2.2",
    "stroke-linecap": "round",
  });
  svg.appendChild(needle);
  svg.appendChild(el("circle", { cx: `${c}`, cy: `${c}`, r: "4.5", fill: BRASS_DARK }));
  const reading = el("text", {
    x: `${c}`,
    y: `${c + 22}`,
    "text-anchor": "middle",
    "font-size": "12",
    fill: INK,
    "font-family": "Georgia, serif",
    "font-weight": "bold",
  });
  svg.appendChild(reading);
  return { svg, dom: { needle, reading, face, center: c } };
}

interface LeverDom {
  pointer: HTMLDivElement;
  value: HTMLDivElement;
  wrap: HTMLDivElement;
}

function buildLever(view: LeverView): { root: HTMLDivElement; dom: LeverDom } {
  const root = document.createElement("div");
  root.className = "wf-lever";
  const label = document.createElement("div");
  label.className = "wf-lever-label";
  label.textContent = view.label;
  const track = document.createElement("div");
  track.className = "wf-lever-track";
  const center = document.createElement("div");
  center.className = "wf-lever-center";
  const pointer = document.createElement("div");
  pointer.className = "wf-lever-pointer";
  pointer.style.left = `${(50 + view.fraction * 46).toFixed(1)}%`;
  track.appendChild(center);
  track.appendChild(pointer);
  const value = document.createElement("div");
  value.className = "wf-lever-value";
  value.textContent = view.reading;
  root.appendChild(label);
  root.appendChild(track);
  root.appendChild(value);
  return { root, dom: { pointer, value, wrap: root } };
}

export interface HudDials {
  update(input: HudDialInputs): void;
  dispose(): void;
}

/** Build the cockpit panel into `container` (bottom-center). */
export function createHudDials(container: HTMLElement): HudDials {
  const panel = document.createElement("div");
  panel.id = "wf-dial-panel";
  const dialDom = new Map<string, DialDom>();
  for (const view of dialSetFrom(IDLE_INPUTS)) {
    const cell = document.createElement("div");
    cell.className = "wf-dial";
    const { svg, dom } = buildDial(view, view.id === "anemometer" ? 128 : 108);
    cell.appendChild(svg);
    panel.appendChild(cell);
    dialDom.set(view.id, dom);
  }
  const leverWrap = document.createElement("div");
  leverWrap.id = "wf-lever-stack";
  const leverDom = new Map<string, LeverDom>();
  for (const view of leverSetFrom(IDLE_INPUTS)) {
    const { root, dom } = buildLever(view);
    leverWrap.appendChild(root);
    leverDom.set(view.id, dom);
  }
  panel.appendChild(leverWrap);
  container.appendChild(panel);
  return {
    update(input: HudDialInputs): void {
      for (const view of dialSetFrom(input)) {
        const dom = dialDom.get(view.id);
        if (dom === undefined) {
          continue;
        }
        const c = dom.center;
        dom.needle.setAttribute("transform", `rotate(${view.needleDeg.toFixed(2)} ${c} ${c})`);
        dom.reading.textContent =
          view.unit === ""
            ? view.reading
            : view.unit === "°"
              ? `${view.reading}°`
              : `${view.reading} ${view.unit}`;
        dom.face.setAttribute("fill", view.danger ? "#f0d3c5" : FACE);
      }
      for (const view of leverSetFrom(input)) {
        const dom = leverDom.get(view.id);
        if (dom === undefined) {
          continue;
        }
        dom.pointer.style.left = `${(50 + view.fraction * 46).toFixed(1)}%`;
        dom.pointer.classList.toggle("at-stop", view.atStop);
        dom.value.textContent = view.atStop ? `${view.reading} — AT STOP` : view.reading;
      }
    },
    dispose(): void {
      panel.remove();
    },
  };
}

/** Big top-center phase banner. */
export interface PhaseBannerEl {
  set(text: string | null, tone: "info" | "good" | "end"): void;
  dispose(): void;
}

export function createPhaseBanner(container: HTMLElement): PhaseBannerEl {
  const elDiv = document.createElement("div");
  elDiv.id = "wf-phase-banner";
  elDiv.style.display = "none";
  container.appendChild(elDiv);
  let last = "";
  return {
    set(text: string | null, tone: "info" | "good" | "end"): void {
      if (text === null) {
        elDiv.style.display = "none";
        last = "";
        return;
      }
      elDiv.style.display = "block";
      elDiv.dataset["tone"] = tone;
      if (text !== last) {
        elDiv.textContent = text;
        last = text;
        elDiv.classList.remove("wf-banner-pop");
        void elDiv.offsetWidth; // restart the pop animation
        elDiv.classList.add("wf-banner-pop");
      }
    },
    dispose(): void {
      elDiv.remove();
    },
  };
}
