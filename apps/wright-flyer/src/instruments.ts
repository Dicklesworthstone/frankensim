// FIELD INSTRUMENTS control panel — wires dormant visualization/lesson modules
// into the Wright Flyer app. New file only; touches nothing else.
//
// Adapter law: every dormant module is opt-in via its checkbox row; enabling
// lazy-imports the module and renders its most useful live readout as a tile
// in the right-side stack, updated from the onFrame sample callback. Each
// adapter is fault-isolated: a broken module shows "unavailable", never throws.
//
// Evidence law (bead wf-root-guzez.15): a tile whose inputs are made up in
// the page says so on its first line (ILLUSTRATION) and carries no evidence.
// The evidence tiles (badges, scorecard, applicability) need real receipts;
// no receipt source is connected to this page, so they show NO_RECEIPT_SOURCE
// and build no badge, digest, verdict or receipt link.

import type { ModeFamily, PublishedLabeledPole } from "./eigenView.ts";
import type { FieldArrays } from "./fieldViz.ts";
import type { StripLoadsState } from "./forceOverlay.ts";
import type { YawDecomposition } from "./lateralView.ts";
import type { DesignPoint, SweepRunner } from "./sweeps.ts";

export interface InstrumentSample {
  tick: number;
  xM: number;
  hM: number;
  uMps: number;
  wMps: number;
  qRadS: number;
  thetaRad: number;
  dcRad: number;
  warpRad: number;
  omegaRadS: number;
  gustWMps: number;
  phase: string;
}

export interface InstrumentsPanel {
  dispose(): void;
}

// ---------------------------------------------------------------- plumbing

type Sample = InstrumentSample | null;

interface LoadedAdapter {
  update(sample: Sample): void;
}

interface TileRecord {
  content: HTMLElement;
  adapter: LoadedAdapter | null;
  pending: boolean;
  enabled: boolean;
}

function num(n: number, digits = 2): string {
  return Number.isFinite(n) ? n.toFixed(digits) : "—";
}

function setLines(content: HTMLElement, lines: readonly string[]): void {
  content.textContent = "";
  for (const line of lines) {
    const d = document.createElement("div");
    d.textContent = line;
    content.appendChild(d);
  }
}

interface AdapterDef {
  readonly key: string;
  readonly label: string;
  load(content: HTMLElement): Promise<LoadedAdapter>;
}

/** First line of every tile whose values are made up in the page. */
export const ILLUSTRATION =
  "ILLUSTRATION: made up in the page from the live state, not engine output; carries no evidence";

/** The whole content of an evidence tile while no receipt source is connected. */
export const NO_RECEIPT_SOURCE =
  "unavailable: no receipt source is connected to this page, so no badge, digest or verdict is shown";

/** The augmented engine's poles at the canonical trim, or why they are missing. */
type EngineModes =
  | {
      readonly ok: true;
      readonly trimVMps: number;
      readonly trimAlphaRad: number;
      readonly poles: readonly PublishedLabeledPole[];
    }
  | { readonly ok: false; readonly reason: string };

let engineModes: Promise<EngineModes> | null = null;

/**
 * fs-flyer's augmented linearization (rigid 4 + canard 2 + prop speed, no
 * pilot) at the canonical trim, from the wasm export flyer_modes_json (bead
 * frankensim-byipt). Computed once per page, on the main thread, when the
 * eigenmode tile is first enabled; the poles belong to the trim, not to the
 * live state.
 */
function loadEngineModes(): Promise<EngineModes> {
  engineModes ??= (async (): Promise<EngineModes> => {
    try {
      const pkg = await import("./wasm-pkg/fs_flyer_wasm.js");
      await pkg.default();
      const out = JSON.parse(pkg.flyer_modes_json()) as {
        refusal?: { code: string };
        trim?: { v_mps: number; alpha_rad: number };
        poles?: { re: number; im: number; family: ModeFamily; shift: number }[];
      };
      if (out.refusal) return { ok: false, reason: `engine refused: ${out.refusal.code}` };
      if (!out.trim || !out.poles) return { ok: false, reason: "engine modes envelope malformed" };
      return {
        ok: true,
        trimVMps: out.trim.v_mps,
        trimAlphaRad: out.trim.alpha_rad,
        poles: out.poles.map((p) => ({ re: p.re, im: p.im, family: p.family, attributionShift: p.shift })),
      };
    } catch (error) {
      return { ok: false, reason: `engine modes unavailable: ${String(error)}` };
    }
  })();
  return engineModes;
}

/**
 * One family's poles, each conjugate pair once as `re ± im i`. A part within
 * 1e-9 of the largest pole magnitude (`scale`) is rounding and reads as 0, so
 * a pure integrator shows as neutral; an unstable root names its doubling time.
 */
function poleLine(family: string, poles: readonly PublishedLabeledPole[], scale: number): string {
  const tiny = 1e-9 * scale;
  const shown = poles
    .filter((p) => p.im >= -tiny)
    .map((p) => {
      const re = Math.abs(p.re) <= tiny ? 0 : p.re;
      const root = Math.abs(p.im) <= tiny ? num(re) : `${num(re)} ± ${num(Math.abs(p.im))}i`;
      if (re > 0) return `${root} (unstable, doubles in ${num(Math.LN2 / re)} s)`;
      return re === 0 ? `${root} (neutral)` : root;
    });
  return `${family}: ${shown.join(", ")} /s`;
}

function noReceiptSource(content: HTMLElement): Promise<LoadedAdapter> {
  setLines(content, [NO_RECEIPT_SOURCE]);
  return Promise.resolve({ update: () => setLines(content, [NO_RECEIPT_SOURCE]) });
}

// --------------------------------------------------------------- adapters

// NOTE: module values are loaded via await import() below BY DESIGN: each
// adapter must stay dormant until its checkbox enables it (lazy loading of
// heavy visualization modules), which a static import cannot express. Types
// are still declared via top-level `import type` above.
const ADAPTERS: readonly AdapterDef[] = [
  {
    key: "force-overlay",
    label: "Force overlay",
    async load(content) {
      const fo = await import("./forceOverlay.ts");
      const massKg = 338;
      const areaM2 = 15.1;
      return {
        update(sample) {
          if (!sample) return;
          const q = 0.5 * 1.225 * sample.uMps * sample.uMps;
          const lift = q * areaM2 * 0.9;
          const drag = q * areaM2 * 0.08;
          const thrust: readonly [number, number, number] = [drag, 0, 0];
          const weight: readonly [number, number, number] = [0, 0, -massKg * 9.81];
          const half = lift / 2;
          const state: StripLoadsState = {
            nStrips: 2,
            positions: new Float64Array([0, -2.4, 0, 0, 2.4, 0]),
            forces: new Float64Array([drag, 0, half, drag, 0, half]),
            thrustN: thrust,
            thrustAt: [0, 0, 1.5],
            weightN: weight,
            cgAt: [0, 0, 0.9],
            netN: [2 * drag + thrust[0], 0, lift + weight[2]],
          };
          const built = fo.buildForceOverlay(state);
          if (!built.ok) {
            setLines(content, [`refused: ${built.refusal.code}`]);
            return;
          }
          const div = fo.firstOverlayDivergence(built.value, state);
          setLines(content, [
            "APPROX aero estimate (Cl 0.9 / Cd 0.08 hardcoded) — NOT sim forces",
            `net [${num(state.netN[0])} ${num(state.netN[1])} ${num(state.netN[2])}] N`,
            `thrust [${num(thrust[0])} 0 0] N`,
            `weight [0 0 ${num(weight[2])}] N`,
            `strips 2 · lift ${num(lift)} N · drag ${num(drag)} N`,
            `bit-for-bit audit: ${div < 0 ? "faithful" : `divergence @${div}`}`,
          ]);
        },
      };
    },
  },
  {
    key: "tracers",
    label: "Flow tracers",
    async load(content) {
      const tr = await import("./tracers.ts");
      const svc = new tr.TracerService();
      const releasePoint: readonly [number, number, number] = [-4, 0, 2];
      return {
        update(sample) {
          if (!sample) return;
          const sampler = (p: readonly [number, number, number]) =>
            p[2] < 0 || p[2] > 80
              ? null
              : ([sample.uMps, 0, sample.wMps + sample.gustWMps] as const);
          const tick = svc.checkpoint().tick;
          if (tick % 10 === 0) svc.release(tick, releasePoint);
          const advanced = svc.advance(
            {
              tickA: tick,
              tickB: tick + 1,
              idA: `s${tick}`,
              idB: `s${tick + 1}`,
              samplerA: sampler,
              samplerB: sampler,
            },
            0.02,
          );
          const cp = svc.checkpoint();
          const alive = cp.tracers.filter((t) => t.alive).length;
          const streak = svc.streakline(releasePoint);
          setLines(content, [
            ILLUSTRATION,
            "tracers ride a uniform wind equal to the aircraft's own air velocity",
            `tick ${cp.tick} · advanced ${
              advanced.ok ? String(advanced.value) : advanced.refusal.code
            }`,
            `tracers ${cp.tracers.length} (${alive} alive)`,
            `retained ${svc.retainedPoints()} pts (cap ${tr.MAX_POINTS_PER_TRACER}/tracer)`,
            `streakline ${streak.length / 3} pts`,
          ]);
        },
      };
    },
  },
  {
    key: "porpoises",
    label: "Why it porpoises",
    async load(content) {
      const pv = await import("./porpoisesView.ts");
      const cmd: number[] = [];
      const act: number[] = [];
      const windowCap = 250; // > 2*MAX_LAG_TICKS + 4
      return {
        update(sample) {
          if (!sample) return;
          cmd.push(sample.dcRad);
          act.push(sample.thetaRad);
          if (cmd.length > windowCap) {
            cmd.shift();
            act.shift();
          }
          const lines: string[] = [
            "ILLUSTRATIVE pole (heuristic from live rates — not engine eigenmodes)",
          ];
          const pole = pv.poleIndicator({
            reSigmaPerS: 0.4 + 0.5 * Math.abs(sample.qRadS),
            imOmegaRadPerS: 1.2 + 0.05 * Math.abs(sample.omegaRadS),
          });
          if (pole.ok) {
            lines.push(
              pole.value.timeToDoubleS === null
                ? "t2× — (non-growing)"
                : `t2× ${num(pole.value.timeToDoubleS)} s`,
            );
            lines.push(
              pole.value.periodS === null
                ? "period — (aperiodic)"
                : `period ${num(pole.value.periodS)} s`,
            );
          } else {
            lines.push(`pole refused: ${pole.refusal.code}`);
          }
          if (cmd.length >= 2 * pv.MAX_LAG_TICKS + 4) {
            const lag = pv.estimateDelayTicks(Float64Array.from(cmd), Float64Array.from(act));
            lines.push(
              lag.ok
                ? `cmd→actual lag ${lag.value > 0 ? "+" : ""}${lag.value} ticks`
                : `lag refused: ${lag.refusal.code}`,
            );
          } else {
            lines.push(`lag: warming up (${cmd.length}/${2 * pv.MAX_LAG_TICKS + 4})`);
          }
          const dc = Math.abs(sample.dcRad);
          const gu = Math.abs(sample.gustWMps);
          const tot = dc + gu + 0.05;
          const att = pv.attributionView([
            { component: "elevator", share: dc / tot },
            { component: "gust", share: gu / tot },
            { component: "trim", share: 0.05 / tot },
          ]);
          lines.push(
            att.ok
              ? `dominant ${att.value.dominant} · residual ${num(att.value.residual, 3)}`
              : `attribution refused: ${att.refusal.code}`,
          );
          setLines(content, lines.slice(0, 6));
        },
      };
    },
  },
  {
    key: "eigenmodes",
    label: "Eigenmodes",
    async load(content) {
      const ev = await import("./eigenView.ts");
      const modes = await loadEngineModes();
      const engineLines: string[] = [];
      if (modes.ok) {
        const alphaDeg = (modes.trimAlphaRad * 180) / Math.PI;
        engineLines.push(
          `ENGINE: fs-flyer augmented linearization at trim (${num(modes.trimVMps, 1)} m/s, α ${num(alphaDeg, 1)}°), no pilot; poles of the trim, not of the live state`,
        );
        const scale = Math.max(...modes.poles.map((p) => Math.hypot(p.re, p.im)));
        const families = ev.groupModeFamilies(modes.poles);
        if (families.ok) {
          for (const g of families.value) engineLines.push(poleLine(g.family, g.poles, scale));
        } else {
          engineLines.push(`families refused: ${families.refusal.code}`);
        }
      } else {
        engineLines.push(modes.reason);
      }
      return {
        update(sample) {
          if (!sample) return;
          const labels = ["u", "w", "q", "theta", "dc", "warp"];
          const mags = new Float64Array([
            Math.max(Math.abs(sample.uMps) / 10, 1e-3),
            Math.max(Math.abs(sample.wMps), 1e-3),
            Math.max(Math.abs(sample.qRadS), 1e-3),
            Math.max(Math.abs(sample.thetaRad), 1e-3),
            Math.max(Math.abs(sample.dcRad), 1e-3),
            Math.max(Math.abs(sample.warpRad), 1e-3),
          ]);
          const proj = ev.teachingProjection(labels, mags);
          const lines: string[] = [...engineLines, ILLUSTRATION];
          if (proj.ok) {
            lines.push(`live-state rigid share ${(proj.value.rigidShare * 100).toFixed(1)}%`);
            lines.push(`beyond-4-state ${(proj.value.beyondFourStateContent * 100).toFixed(1)}%`);
            lines.push(proj.value.caption);
          } else {
            lines.push(`projection refused: ${proj.refusal.code}`);
          }
          setLines(content, lines);
        },
      };
    },
  },
  {
    key: "field-viz",
    label: "Field glyphs",
    async load(content) {
      const fv = await import("./fieldViz.ts");
      return {
        update(sample) {
          if (!sample) return;
          const n = 24; // 4 × 3 × 2 probe grid around the plane
          const points = new Float64Array(3 * n);
          const u = new Float64Array(3 * n);
          const divAnalytic = new Float64Array(n);
          const divFd = new Float64Array(n);
          const gradNorm = new Float64Array(n);
          const validity = new Uint8Array(n).fill(1);
          const singularityCore = new Uint8Array(n);
          let k = 0;
          for (let ix = 0; ix < 4; ix += 1) {
            for (let iy = 0; iy < 3; iy += 1) {
              for (let iz = 0; iz < 2; iz += 1) {
                points[3 * k] = sample.xM - 3 + 2 * ix;
                points[3 * k + 1] = iy;
                points[3 * k + 2] = Math.max(sample.hM, 2) - 1 + iz;
                const shear = 1 + 0.02 * ix;
                u[3 * k] = sample.uMps * shear;
                u[3 * k + 1] = 0;
                u[3 * k + 2] = sample.wMps + sample.gustWMps;
                divAnalytic[k] = 0.02 * sample.uMps;
                divFd[k] = 0.02 * sample.uMps;
                gradNorm[k] = k === 0 ? 0 : 0.02 * sample.uMps;
                singularityCore[k] = k === 5 ? 1 : 0;
                k += 1;
              }
            }
          }
          const field: FieldArrays = {
            n,
            points,
            u,
            divAnalytic,
            divFd,
            gradNorm,
            validity,
            singularityCore,
            omittedComponents: ["viscous"],
            forceCoupledSupported: ["pressure", "viscous"],
          };
          const glyphs = fv.buildGlyphInstances(field);
          const overlay = fv.divergenceOverlay(field, "analytic");
          let masked = 0;
          for (let i = 0; i < overlay.masked.length; i += 1) {
            if ((overlay.masked[i] ?? 0) !== 0) masked += 1;
          }
          const legend = fv.legendConfig("velocity", field);
          setLines(content, [
            ILLUSTRATION,
            "field = the aircraft's air velocity with a 2% shear, not the engine's field service",
            glyphs.ok
              ? `glyphs ${glyphs.value.count} / ${n} pts`
              : `glyphs refused: ${glyphs.refusal.code}`,
            `masked (|div| overlay): ${masked}`,
            `legend: ${legend.label} [${legend.units}]`,
            `omitted: ${legend.omittedComponents.join(", ") || "none"}`,
            fv.WAKE_FADE_PRESENTATION_ONLY,
          ]);
        },
      };
    },
  },
  {
    key: "lessons",
    label: "Lessons",
    async load(content) {
      const { curatedLessons, validateLesson } = await import("./lessons.ts");
      const catalog = curatedLessons();
      const guide = document.createElement("p");
      guide.textContent =
        "Reading guide · static teaching material. The named views and references describe each lesson. Reading a step does not run an experiment or produce evidence.";

      const lessonChoice = document.createElement("select");
      const stepChoice = document.createElement("select");
      function labeledChoice(title: string, choice: HTMLSelectElement): HTMLLabelElement {
        const label = document.createElement("label");
        label.style.display = "block";
        choice.style.maxWidth = "100%";
        // Native select navigation must not also steer the aircraft.
        choice.addEventListener("keydown", (event) => event.stopPropagation());
        label.append(title, document.createElement("br"), choice);
        return label;
      }
      const notes = document.createElement("div");
      notes.setAttribute("aria-live", "polite");
      notes.setAttribute("aria-atomic", "true");
      for (const lesson of catalog) {
        lessonChoice.add(new Option(lesson.title, lesson.id));
      }
      content.replaceChildren(guide, labeledChoice("Lesson", lessonChoice),
        labeledChoice("Step", stepChoice), notes);

      function showStep(): void {
        const lesson = catalog[lessonChoice.selectedIndex];
        const step = lesson?.steps[stepChoice.selectedIndex];
        if (lesson === undefined || step === undefined) return;
        const views = document.createElement("p");
        views.textContent = `Named views: ${step.overlays.join(", ")}`;
        const statements = document.createElement("ul");
        for (const claim of step.voicedClaims) {
          const item = document.createElement("li");
          item.textContent = claim;
          statements.appendChild(item);
        }
        notes.replaceChildren(views);
        if (statements.childElementCount > 0) notes.appendChild(statements);
      }
      function selectLesson(): void {
        stepChoice.replaceChildren();
        const lesson = catalog[lessonChoice.selectedIndex];
        if (lesson === undefined) {
          stepChoice.disabled = true;
          setLines(notes, ["No lessons available."]);
          return;
        }
        const check = validateLesson(lesson);
        stepChoice.disabled = !check.ok;
        if (!check.ok) {
          setLines(notes, [check.refusal.message]);
          return;
        }
        lesson.steps.forEach((step, index) => {
          stepChoice.add(new Option(`${index + 1} of ${lesson.steps.length} · ${step.title}`));
        });
        showStep();
      }
      lessonChoice.addEventListener("change", selectLesson);
      stepChoice.addEventListener("change", showStep);
      selectLesson();
      return { update() {} };
    },
  },
  {
    key: "badges",
    label: "Evidence badges",
    load: noReceiptSource,
  },
  {
    key: "scorecard",
    label: "Scorecard",
    load: noReceiptSource,
  },
  {
    key: "sweeps",
    label: "Sweeps",
    async load(content) {
      const sw = await import("./sweeps.ts");
      const grid = sw.makeSweepGrid(
        { name: "u", values: [8, 10, 12] },
        { name: "h", values: [3, 10] },
      );
      if (!grid.ok) {
        return {
          update() {
            setLines(content, [`grid refused: ${grid.refusal.code}`]);
          },
        };
      }
      const engine = sw.makeSweepEngine(grid.value, [11, 23], "wf-inst-1", new Map());
      if (!engine.ok) {
        return {
          update() {
            setLines(content, [`engine refused: ${engine.refusal.code}`]);
          },
        };
      }
      const eng = engine.value;
      const points: DesignPoint[] = grid.value;
      const runner: SweepRunner = (point, seed) => {
        const uu = point.config["u"] ?? 0;
        const hh = point.config["h"] ?? 1;
        return (uu / Math.max(hh, 0.5)) * (1 + 0.001 * seed);
      };
      return {
        update() {
          for (let i = 0; i < 4; i += 1) eng.step(() => true, runner);
          const prog = eng.progress();
          const last = eng.records[eng.records.length - 1];
          const csvRows = sw.exportCsv(eng, points).split("\n").length - 1;
          setLines(content, [
            ILLUSTRATION,
            "each sweep point is u / h, a placeholder, not a flight",
            `progress ${prog.completed}/${prog.total}`,
            last === undefined
              ? "no units yet"
              : `last: point ${last.pointIndex} member ${last.member} → ${num(last.value)}${
                  last.fromCache ? " (cached)" : ""
                }`,
            `csv rows ${csvRows} · QoS gate open`,
          ]);
        },
      };
    },
  },
  {
    key: "lateral",
    label: "Lateral view",
    async load(content) {
      const lv = await import("./lateralView.ts");
      const windowRows: YawDecomposition[] = [];
      return {
        update(sample) {
          if (!sample) return;
          const dyn = 0.5 * 1.225 * sample.uMps * sample.uMps;
          const induced = -0.4 * sample.warpRad * dyn * 0.5; // adverse-yaw sign
          const rudder = 0.3 * sample.dcRad * dyn * 0.5;
          const profile = 0.01 * sample.omegaRadS;
          windowRows.push({
            tick: sample.tick,
            warpCommandRad: sample.warpRad,
            loadedTwistRad: 0.8 * sample.warpRad,
            inducedDragYawNm: induced,
            rudderYawNm: rudder,
            profileYawNm: profile,
            netYawNm: induced + rudder + profile,
          });
          if (windowRows.length > 240) windowRows.shift();
          const verdict = lv.adverseYawVerdict(windowRows);
          const spiral = lv.spiralIndicator(0.05 * sample.omegaRadS);
          const lines = [
            ILLUSTRATION,
            "yaw terms use fixed page coefficients; the spiral sign follows propeller speed",
            verdict.ok
              ? `adverse yaw: ${verdict.value.adverse ? "YES" : "no"} (mean sign product ${num(
                  verdict.value.meanSignProduct,
                  3,
                )})`
              : `verdict: ${verdict.refusal.code} (commanded ticks pending)`,
            verdict.ok
              ? `commanded ticks ${verdict.value.commandedTicks}`
              : `window ${windowRows.length}`,
          ];
          if (spiral.ok) {
            lines.push(
              spiral.value.divergent
                ? `spiral DIVERGENT σ=${num(spiral.value.reSigmaPerS, 3)} t2× ${num(
                    spiral.value.timeToDoubleS ?? Number.NaN,
                  )} s`
                : `spiral convergent σ=${num(spiral.value.reSigmaPerS, 3)}`,
            );
          } else {
            lines.push(`spiral refused: ${spiral.refusal.code}`);
          }
          setLines(content, lines);
        },
      };
    },
  },
  {
    key: "applicability",
    label: "Applicability",
    load: noReceiptSource,
  },
];

// ------------------------------------------------------------ the panel

export function createInstrumentsPanel(
  container: HTMLElement,
  hooks: {
    onFrame(cb: (sample: InstrumentSample | null) => void): void;
  },
): InstrumentsPanel {
  const panel = document.createElement("div");
  panel.className = "wf-inst-panel";

  const head = document.createElement("div");
  head.className = "wf-inst-head";
  head.textContent = "FIELD INSTRUMENTS";
  panel.appendChild(head);

  const tiles = new Map<string, TileRecord>();
  const checkboxes = new Map<string, HTMLInputElement>();

  const view = document.createElement("div");
  view.className = "wf-inst-view wf-inst-view-hidden";

  function ensureTile(def: AdapterDef): TileRecord {
    let rec = tiles.get(def.key);
    if (rec === undefined) {
      const tile = document.createElement("div");
      tile.className = "wf-inst-tile";
      const label = document.createElement("div");
      label.className = "wf-inst-tile-label";
      label.textContent = def.label;
      tile.appendChild(label);
      const content = document.createElement("div");
      tile.appendChild(content);
      view.appendChild(tile);
      rec = { content, adapter: null, pending: false, enabled: false };
      tiles.set(def.key, rec);
    }
    return rec;
  }

  function syncView(): void {
    let any = false;
    for (const rec of tiles.values()) {
      if (rec.enabled) {
        any = true;
        break;
      }
    }
    view.classList.toggle("wf-inst-view-hidden", !any);
  }

  function syncUrl(): void {
    try {
      const enabled = [...tiles.entries()]
        .filter(([, rec]) => rec.enabled)
        .map(([key]) => key)
        .sort();
      const url = new URL(window.location.href);
      if (enabled.length > 0) url.searchParams.set("inst", enabled.join(","));
      else url.searchParams.delete("inst");
      window.history.replaceState(null, "", url);
    } catch {
      /* URL sync is best-effort */
    }
  }

  function setEnabled(key: string, on: boolean): void {
    const def = ADAPTERS.find((a) => a.key === key);
    if (def === undefined) return;
    const cb = checkboxes.get(key);
    if (cb !== undefined) cb.checked = on;
    if (!on) {
      const rec = tiles.get(key);
      if (rec !== undefined) {
        rec.enabled = false;
        rec.adapter = null;
        rec.content.parentElement?.remove();
        tiles.delete(key);
      }
      return;
    }
    const rec = ensureTile(def);
    rec.enabled = true;
    if (rec.adapter !== null || rec.pending) return;
    rec.pending = true;
    setLines(rec.content, ["loading…"]);
    def
      .load(rec.content)
      .then((adapter) => {
        rec.pending = false;
        rec.adapter = adapter;
      })
      .catch(() => {
        rec.pending = false;
        setLines(rec.content, ["unavailable"]);
      });
  }

  for (const def of ADAPTERS) {
    const row = document.createElement("label");
    row.className = "wf-inst-row";
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.className = "wf-inst-check";
    cb.addEventListener("change", () => {
      setEnabled(def.key, cb.checked);
      syncView();
      syncUrl();
    });
    const text = document.createElement("span");
    text.textContent = def.label;
    row.appendChild(cb);
    row.appendChild(text);
    panel.appendChild(row);
    checkboxes.set(def.key, cb);
  }

  panel.appendChild(view);
  container.appendChild(panel);

  // ?inst=csv-keys read on init
  try {
    const raw = new URLSearchParams(window.location.search).get("inst");
    if (raw !== null && raw.trim() !== "") {
      for (const key of raw.split(",")) {
        if (ADAPTERS.some((a) => a.key === key.trim())) setEnabled(key.trim(), true);
      }
    }
  } catch {
    /* location parsing is best-effort */
  }
  syncView();

  hooks.onFrame((sample) => {
    for (const rec of tiles.values()) {
      if (!rec.enabled || rec.adapter === null) continue;
      try {
        rec.adapter.update(sample);
      } catch {
        setLines(rec.content, ["unavailable"]);
      }
    }
  });

  let hidden = false;
  const onKey = (e: KeyboardEvent): void => {
    if (e.key === "i" || e.key === "I") {
      hidden = !hidden;
      panel.style.display = hidden ? "none" : "";
    }
  };
  window.addEventListener("keydown", onKey);

  return {
    dispose() {
      window.removeEventListener("keydown", onKey);
      panel.remove();
    },
  };
}
