/**
 * The morning (chapter 07): a catering commissary's hot side from across its worktable, every meal leaving on one van
 * at 10:30.
 *
 * Everything shown follows `planMorning` (sim/cateringMorning.ts) for the settings on the page, acted out on the
 * kitchen (MorningKitchen, morningKitchen.ts). The chart is the plan itself, minute for minute: each cook's hands,
 * each burner, each oven and each table place.
 */
import * as THREE from "three";
import { BURNERS, DOOR, MENUS, OVENS, SERVE, SLOTS, TABLE_SPOTS, clockOf, madeIt, menuDishes, planMorning, type Crew, type Morning, type MorningSettings } from "../sim/cateringMorning";
import type { LimitAsk, LimitScore } from "../sim/morningLimits.worker";
import { bindCanvasResize, clipFitter, createScrollSafeOrbit, createStudioRenderer, observeInView, shaderGate } from "./kit";
import { CARRIER, CARRIER_COLS, CARRIER_X0, CUTOFF, END_HOLD_S, LINE_Z, MorningKitchen, OVEN_X0, RANGE, SEC_PER_MIN, SLOT_COLOR, TABLE, TABLE_Z, fmtLate, shortName, type CookSpec } from "./morningKitchen";
import { studioEnvironment, studioLights } from "./models/stage";
import { backdrop, tilePlane } from "./weekLab";

export { MORNING_PIECES, SLOT_COLOR } from "./morningKitchen";

/** The cooks, the first two to four of them on shift. */
const COOKS: CookSpec[] = [
  { cohort: "p50", seed: 41, apron: 0x26344d },
  { cohort: "p95", seed: 57, apron: 0x3b3d43 },
  { cohort: "p05", seed: 63, apron: 0x45563f },
  { cohort: "p50", seed: 71, apron: 0x4a3a2e },
];

/* ------------------------------------------------------------------------------------------------------------ */

export function mountMorningLab(root: HTMLElement): void {
  const mountT0 = performance.now();
  const stageEl = root.querySelector<HTMLElement>("[data-morning-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-morning-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = createStudioRenderer(canvas, false);
  if (!maybeRenderer) return;
  const renderer = maybeRenderer;
  renderer.toneMappingExposure = 1.0;
  const phone = window.matchMedia("(max-width: 720px)").matches;

  const scene = new THREE.Scene();
  scene.background = backdrop();
  scene.fog = new THREE.Fog(0x151510, 8, 16);
  scene.environment = studioEnvironment(renderer);
  scene.environmentIntensity = 0.6;
  const lights = studioLights(scene, new THREE.Vector3(-0.4, 0.9, 0.1), 4.6, phone ? 1024 : 2048);
  lights.key.position.set(-1.2, 6.4, 4.6);
  lights.key.intensity = 1.9;
  lights.key.shadow.normalBias = 0.02;
  lights.rim.position.set(3.4, 3.2, -3.0);
  lights.rim.intensity = 0.75;
  lights.hemi.intensity = 0.5;

  const floor = tilePlane(9, 6);
  floor.rotation.x = -Math.PI / 2;
  floor.position.set(0, 0, 0.2);
  scene.add(floor);

  const kitchen = new MorningKitchen(COOKS);
  scene.add(kitchen.group);
  const cooks = kitchen.cooks;

  /* Camera: from behind the worktable, raised, the line across the aisle. */
  const camera = new THREE.PerspectiveCamera(32, 16 / 10, 0.05, 60);
  const TARGET = new THREE.Vector3();
  const VIEW = new THREE.Vector3();
  const keep: THREE.Vector3[] = [];
  const aim = () => {
    keep.length = 0;
    const wide = stageEl.clientWidth >= 1.3 * stageEl.clientHeight;
    if (wide) {
      TARGET.set(-0.1, 0.75, 0.0);
      VIEW.set(0.1, 1.15, 1).normalize();
      for (const x of [OVEN_X0[0], CARRIER_X0 + CARRIER_COLS * (CARRIER.w + 0.03)]) keep.push(new THREE.Vector3(x, 1.85, LINE_Z - 0.4), new THREE.Vector3(x, 0.2, TABLE_Z + TABLE.d));
    } else {
      TARGET.set(-0.6, 0.8, 0.0);
      VIEW.set(0.35, 1.2, 1).normalize();
      for (const x of [OVEN_X0[0] + 0.1, RANGE.x0 + RANGE.w]) for (const y of [0.2, 1.8]) for (const z of [LINE_Z - 0.3, TABLE_Z + TABLE.d]) keep.push(new THREE.Vector3(x, y, z));
    }
  };
  const controls = createScrollSafeOrbit(camera, canvas, [TARGET.x, TARGET.y, TARGET.z], { minDistance: 1.2, maxDistance: 14, maxPolarAngle: Math.PI * 0.47 });
  let userMoved = false;
  controls.addEventListener("start", () => {
    userMoved = true;
  });
  const ndc = new THREE.Vector3();
  const hudEl = stageEl.querySelector<HTMLElement>(".mo-hud");
  const legendEl = stageEl.querySelector<HTMLElement>(".mo-legend");
  const frame = () => {
    if (userMoved) return;
    aim();
    const w = Math.max(1, stageEl.clientWidth);
    const h = Math.max(1, stageEl.clientHeight);
    camera.aspect = Math.max(0.2, w / h);
    camera.clearViewOffset();
    camera.updateProjectionMatrix();
    // The kitchen fits between the HUD and the legend (on a phone they take a big share of the stage): fitted to the
    // band's height, then the view shifted to sit in it.
    const top = Math.min(0.86, 1 - (2 * ((hudEl ? hudEl.offsetTop + hudEl.offsetHeight : 0) + 6)) / h);
    const bot = -Math.min(0.9, 1 - (2 * ((legendEl ? h - legendEl.offsetTop : 0) + 4)) / h);
    let dist = 2.5;
    let lo = 0;
    let hi = 0;
    for (; dist < 20; dist += 0.06) {
      camera.position.copy(TARGET).addScaledVector(VIEW, dist);
      camera.lookAt(TARGET);
      camera.updateMatrixWorld(true);
      lo = Infinity;
      hi = -Infinity;
      let wide = false;
      for (const p of keep) {
        ndc.copy(p).project(camera);
        lo = Math.min(lo, ndc.y);
        hi = Math.max(hi, ndc.y);
        if (Math.abs(ndc.x) >= 0.95) wide = true;
      }
      if (!wide && hi - lo < top - bot) break;
    }
    const shift = (top + bot) / 2 - (hi + lo) / 2;
    if (Math.abs(shift) > 0.005) {
      camera.setViewOffset(w, h, 0, (shift * h) / 2, w, h);
      camera.updateProjectionMatrix();
    }
    const fog = scene.fog as THREE.Fog;
    fog.near = dist * 1.1;
    fog.far = dist * 2.6;
    controls.target.copy(TARGET);
    controls.update();
  };
  const gate = shaderGate(renderer, camera);
  bindCanvasResize(stageEl, renderer, camera, gate);
  frame();
  new ResizeObserver(() => frame()).observe(stageEl);

  /* ---------------------------------------------------------------------------------------------------------- */
  /* The page.                                                                                                   */
  /* ---------------------------------------------------------------------------------------------------------- */

  const q = <T extends Element = HTMLElement>(s: string) => root.querySelector<T>(s);
  const crewBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-mo-crew]"));
  const aheadBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-mo-ahead]"));
  const presetBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-mo-preset]"));
  const menuEl = q("[data-mo-menu]");
  const mealsIn = q<HTMLInputElement>("[data-mo-meals]");
  const mealsVal = q("[data-mo-meals-val]");
  const statusEl = q("[data-mo-status]");
  const chartEl = q("[data-morning-chart]");
  const ladderEl = q("[data-morning-ladder]");
  const ladderNote = q("[data-morning-ladder-note]");
  const hudClock = q("[data-mo-hud-clock]");
  const hudPans = q("[data-mo-hud-pans]");
  const hudNote = q("[data-mo-hud-note]");
  const tile = (k: string) => ({ v: q(`[data-mo-${k}]`), sub: q(`[data-mo-${k}-sub]`), delta: q(`[data-mo-${k}-d]`) });
  const tiles = { van: tile("van"), cost: tile("cost"), hands: tile("hands"), limit: tile("limit") };

  const settings: MorningSettings = { meals: Number(mealsIn?.value ?? 700), cooks: 3, menu: 0, ahead: false };

  /* The menu card: each place's dish, where it cooks and how long it holds; a press swaps it. */
  const slotBtns: HTMLButtonElement[] = [];
  if (menuEl) {
    menuEl.innerHTML = "";
    SLOTS.forEach((_, k) => {
      const li = document.createElement("li");
      const b = document.createElement("button");
      b.type = "button";
      b.className = "wk-dish mo-dish";
      b.style.setProperty("--dish", SLOT_COLOR[k]!);
      b.addEventListener("click", () => {
        settings.menu ^= 1 << k;
        update();
      });
      li.appendChild(b);
      menuEl.appendChild(li);
      slotBtns.push(b);
    });
  }
  const holdBars = (min: number) => {
    const n = min >= 210 ? 4 : min >= 140 ? 3 : min >= 85 ? 2 : 1;
    return `<span class="wk-keeps" aria-hidden="true">${[1, 2, 3, 4].map((k) => `<i${k <= n ? ' class="on"' : ""}></i>`).join("")}<small>${min >= 120 ? `${min / 60} h` : `${min} min`}</small></span>`;
  };
  const whereOf = (k: number) => {
    const d = menuDishes(settings.menu)[k]!;
    const st = settings.ahead && d.reheat ? d.reheat : d.stages;
    const range = st.find((x) => x.at === "range");
    const oven = st.find((x) => x.at === "oven");
    if (settings.ahead && d.reheat) return "Made yesterday";
    if (range && oven) return "Range, then oven";
    if (st.some((x) => x.at === "range" && (x.hands ?? 1) >= 1)) return "Range, never left";
    if (range) return "Range, 2 pans a cook";
    return "Oven, walk away";
  };
  function renderMenu(): void {
    const dishes = menuDishes(settings.menu);
    slotBtns.forEach((b, k) => {
      const d = dishes[k]!;
      const swapped = (settings.menu & (1 << k)) !== 0;
      const other = swapped ? SLOTS[k]!.written : SLOTS[k]!.swap;
      b.innerHTML = `<b>${d.name}</b>${holdBars(d.holdMin)}<em>${whereOf(k)}</em>`;
      b.setAttribute("aria-pressed", String(swapped));
      b.title = `Swap for ${other.name.toLowerCase()}`;
    });
    for (const b of crewBtns) b.setAttribute("aria-pressed", String(Number(b.dataset.moCrew) === settings.cooks));
    for (const b of aheadBtns) b.setAttribute("aria-pressed", String((b.dataset.moAhead === "yes") === settings.ahead));
    for (const b of presetBtns) {
      const want = b.dataset.moPreset === "written" ? 0 : MENUS - 1;
      b.setAttribute("aria-pressed", String(settings.menu === want));
    }
    if (mealsVal) mealsVal.textContent = `${settings.meals.toLocaleString("en-US")} meals`;
  }

  /* The chart: the plan, lane by lane. Each cook's lane is their hands (a pan minded is half a lane); then the six
   * burners, the two ovens (loading and unloading darker) and the table's three places. */
  const X0 = 64;
  let X1 = 600;
  let W = 610;
  let domain = 300;
  const xOf = (m: number) => X0 + (Math.min(m, domain) / domain) * (X1 - X0);
  let cursor: SVGLineElement | null = null;
  const COOK_H = 15;
  const BURN_H = 5;
  const OVEN_H = 13;
  const SPOT_H = 6;
  function renderChart(m: Morning, endMin: number): void {
    if (!chartEl) return;
    W = Math.max(300, Math.round(chartEl.clientWidth - 2 * 12));
    X1 = W - 8;
    domain = Math.max(300, Math.ceil(endMin / 30) * 30);
    const n = settings.cooks;
    const y0 = 20;
    const cookY = (c: number) => y0 + c * (COOK_H + 3);
    const rangeY = cookY(n) + 10;
    const ovenY = rangeY + BURNERS * BURN_H + 10;
    const tableY = ovenY + OVENS * (OVEN_H + 3) + 7;
    const axisY = tableY + TABLE_SPOTS * SPOT_H + 6;
    const H = axisY + 22;
    const r = (x0: number, x1: number, y: number, h: number, fill: string, op = 1, cls = "") =>
      x1 - x0 > 0.05 ? `<rect${cls ? ` class="${cls}"` : ""} x="${xOf(x0).toFixed(1)}" y="${y.toFixed(1)}" width="${Math.max(0.6, xOf(x1) - xOf(x0)).toFixed(1)}" height="${h.toFixed(1)}" fill="${fill}"${op < 1 ? ` fill-opacity="${op}"` : ""}/>` : "";
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg mo-svg" aria-hidden="true">`;
    // Late: from the van's time to when it actually leaves.
    if (m.late > 0.5) svg += `<rect class="late" x="${xOf(DOOR)}" y="${y0 - 4}" width="${xOf(m.depart) - xOf(DOOR)}" height="${axisY - y0 + 4}"/>`;
    const step = W < 460 ? 60 : 30;
    for (let t = 0; t <= domain; t += step) {
      const label = t % (W < 400 ? 120 : 60) === 0 && !(xOf(t) > X1 - 14);
      svg += `<line class="g v" x1="${xOf(t)}" x2="${xOf(t)}" y1="${y0 - 4}" y2="${axisY}"/>${label ? `<text class="xa" x="${xOf(t)}" y="${axisY + 15}">${clockOf(t)}</text>` : ""}`;
    }
    // Lanes.
    for (let c = 0; c < n; c++) svg += `<rect class="lane" x="${X0}" y="${cookY(c)}" width="${X1 - X0}" height="${COOK_H}"/><text class="ya" x="${X0 - 6}" y="${cookY(c) + COOK_H / 2 + 3.5}">Cook ${c + 1}</text>`;
    svg += `<rect class="lane" x="${X0}" y="${rangeY}" width="${X1 - X0}" height="${BURNERS * BURN_H}"/><text class="ya" x="${X0 - 6}" y="${rangeY + (BURNERS * BURN_H) / 2 + 3.5}">Burners</text>`;
    for (let o = 0; o < OVENS; o++) svg += `<rect class="lane" x="${X0}" y="${ovenY + o * (OVEN_H + 3)}" width="${X1 - X0}" height="${OVEN_H}"/><text class="ya" x="${X0 - 6}" y="${ovenY + o * (OVEN_H + 3) + OVEN_H / 2 + 3.5}">Oven ${o + 1}</text>`;
    svg += `<rect class="lane" x="${X0}" y="${tableY}" width="${X1 - X0}" height="${TABLE_SPOTS * SPOT_H}"/><text class="ya" x="${X0 - 6}" y="${tableY + (TABLE_SPOTS * SPOT_H) / 2 + 3.5}">Table</text>`;
    // Cooks' hands: whole-hand work fills the lane; pans minded stack in halves.
    for (let c = 0; c < n; c++) {
      const jobs = m.jobs.filter((j) => j.cook === c).sort((a, b) => a.t0 - b.t0);
      // When each half of the lane is free again: the lower half for a pan minded, the upper for a second beside it.
      const halfFree = [-Infinity, -Infinity];
      for (const j of jobs) {
        if (j.hands >= 1 - 1e-9) svg += r(j.t0, j.t1, cookY(c), COOK_H, SLOT_COLOR[j.slot]!, j.kind === "prep" ? 0.6 : 0.95);
        else {
          const half = halfFree[0]! <= j.t0 + 1e-6 ? 0 : 1;
          halfFree[half] = j.t1;
          svg += r(j.t0, j.t1, cookY(c) + (half ? 0 : COOK_H / 2), COOK_H / 2 - 0.5, SLOT_COLOR[j.slot]!, 0.85);
        }
      }
    }
    for (const b of m.burns) svg += r(b.t0, b.t1, rangeY + b.burner * BURN_H + 0.5, BURN_H - 1, SLOT_COLOR[b.slot]!, 0.9);
    for (const l of m.loads) {
      const y = ovenY + l.oven * (OVEN_H + 3);
      // Loading, the program, done and shut until a cook is free (fainter), unloading.
      svg += r(l.t0, l.t1, y, OVEN_H, SLOT_COLOR[l.slot]!, 1) + r(l.t1, l.t2, y + 2, OVEN_H - 4, SLOT_COLOR[l.slot]!, 0.55) + r(l.t2, l.tu, y + 4, OVEN_H - 8, SLOT_COLOR[l.slot]!, 0.3) + r(l.tu, l.t3, y, OVEN_H, SLOT_COLOR[l.slot]!, 1);
      const w = xOf(l.t2) - xOf(l.t1);
      if (w > 46) svg += `<text class="lbl" x="${((xOf(l.t1) + xOf(l.t2)) / 2).toFixed(1)}" y="${y + OVEN_H / 2 + 3.5}">${l.verb} · ${l.pans} pans</text>`;
    }
    for (const j of m.jobs) if (j.at === "table") svg += r(j.t0, j.t1, tableY + j.spot * SPOT_H + 0.5, SPOT_H - 1, SLOT_COLOR[j.slot]!, j.kind === "prep" ? 0.6 : 0.95);
    // The van.
    // The van's label: to the left of its line when a late departure's label takes the right, or near the edge.
    const vanR = m.late > 0.5 || xOf(DOOR) > X1 - 100;
    svg += `<line class="van" x1="${xOf(DOOR)}" x2="${xOf(DOOR)}" y1="${y0 - 8}" y2="${axisY}"/><text class="vanl${vanR ? " r" : ""}" x="${xOf(DOOR) + (vanR ? -4 : 4)}" y="${y0 - 9}">Van due 10:30</text>`;
    if (m.late > 0.5) svg += `<line class="van late" x1="${xOf(m.depart)}" x2="${xOf(m.depart)}" y1="${y0 - 8}" y2="${axisY}"/><text class="vanl late${xOf(m.depart) > X1 - 90 ? " r" : ""}" x="${xOf(m.depart) + (xOf(m.depart) > X1 - 90 ? -4 : 4)}" y="${y0 - 9}">${m.finished ? `Leaves ${clockOf(m.depart)}` : "Still waiting"}</text>`;
    svg += `<line class="cur" x1="${X0}" x2="${X0}" y1="${y0 - 6}" y2="${axisY + 2}"/></svg>`;
    const dishes = menuDishes(settings.menu);
    const keys = dishes.map((d, k) => `<span style="--k:${SLOT_COLOR[k]}">${d.name}</span>`).join("");
    chartEl.innerHTML = svg + `<p class="wk-keys mo-keys">${keys}<span class="half">Half a lane: a pan minded</span></p>`;
    cursor = chartEl.querySelector<SVGLineElement>(".cur");
  }

  function renderTiles(m: Morning, base: Morning | null): void {
    const set = (t: (typeof tiles)[keyof typeof tiles], v: string, sub: string, delta: number | null, unit: string, upGood = false) => {
      if (t.v) t.v.textContent = v;
      if (t.sub) t.sub.textContent = sub;
      if (t.delta) {
        t.delta.hidden = delta === null || Math.abs(delta) < 0.5;
        if (delta !== null) {
          t.delta.textContent = `${delta > 0 ? "+" : "−"}${Math.abs(Math.round(delta)).toLocaleString("en-US")}${unit}`;
          t.delta.dataset.good = String(upGood ? delta > 0 : delta < 0);
        }
      }
    };
    if (!m.finished) set(tiles.van, "Not yet", `still waiting at ${CUTOFF}`, null, " min");
    else set(tiles.van, clockOf(m.depart), fmtLate(m.late), base?.finished ? m.depart - base.depart : null, " min");
    const cost = m.late > 0.5 ? m.onTimePastBest : 0;
    set(tiles.cost, cost === null ? "Can't" : cost.toLocaleString("en-US"), cost === null ? "not even starting at 6:00" : cost ? "portions past their best" : "nothing past its best", null, "");
    const span = m.depart;
    const pct = Math.round((100 * m.handsMin) / (settings.cooks * span));
    set(tiles.hands, `${pct}%`, settings.ahead && m.aheadMin > 0 ? `${(m.aheadMin / 60).toFixed(1)} h moved to the day before` : `${(m.handsMin / 60).toFixed(1)} of ${((settings.cooks * span) / 60).toFixed(1)} crew-hours`, base ? Math.round((100 * m.handsMin) / (settings.cooks * span) - (100 * base.handsMin) / (settings.cooks * base.depart)) : null, " pts");
  }

  let swapTimer = 0;
  function renderStatus(m: Morning): void {
    window.clearTimeout(swapTimer);
    if (!statusEl) return;
    const dishes = menuDishes(settings.menu);
    const station = m.bottleneck === "range" ? "the range" : m.bottleneck === "oven" ? "the ovens" : m.bottleneck === "hands" ? "the cooks' hands" : "the table";
    let text: string;
    if (m.late <= 0.5) {
      const last = Math.max(...m.batches.map((b) => b.packed));
      text = `Everything is packed by ${clockOf(last)} and the van leaves at 10:30, every dish at its best at lunch.`;
    } else {
      // The place packed last, the sides among them.
      const lateDish = dishes[m.lastPacked.indexOf(Math.max(...m.lastPacked))]!;
      const busy = Math.round(100 * m.use[m.bottleneck]);
      text = `${m.finished ? `The van leaves at ${clockOf(m.depart)}, ${fmtLate(m.late).toLowerCase()}, waiting on` : `The van is still waiting at ${CUTOFF} for`} the ${lateDish.name.toLowerCase()}: ${station} ${m.bottleneck === "hands" ? "are" : "is"} busy ${busy}% of the 90 minutes before ${m.finished ? "it leaves" : CUTOFF}${m.bottleneck === "range" && settings.cooks > 2 ? ", and only two cooks fit at it" : ""}.`;
      if (m.onTimePastBest !== null && m.onTimePastBest > 0) text += ` Started early enough for 10:30, ${m.onTimePastBest.toLocaleString("en-US")} portions would be past their best by lunch.`;
      // The one swap that gets the van away soonest, when one gets it away sooner: up to three more mornings to plan,
      // so once the settings rest rather than on every step of the slider.
      const said = text;
      const s = { ...settings };
      swapTimer = window.setTimeout(() => {
        let swap = -1;
        let soonest = m.depart - 0.5;
        SLOTS.forEach((_, k) => {
          if (s.menu & (1 << k)) return;
          const depart = runFor({ ...s, menu: s.menu | (1 << k) }).depart;
          if (depart < soonest) {
            soonest = depart;
            swap = k;
          }
        });
        if (swap >= 0 && statusEl) statusEl.textContent = `${said} Try ${SLOTS[swap]!.swap.name.toLowerCase()}.`;
      }, 250);
    }
    statusEl.textContent = text;
    statusEl.dataset.state = m.late <= 0.5 ? "ok" : m.late <= 30 ? "warn" : "danger";
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* Where the limit is: each menu's most meals by 10:30, run in a worker for the crew on the page.              */
  /* ---------------------------------------------------------------------------------------------------------- */

  const limits = new Map<string, number>();
  let limitKey = "";
  let worker: Worker | null = null;
  let ladderQueued = false;
  const lkey = (menu: number, ahead: boolean) => `${menu}|${ahead ? 1 : 0}`;
  function askLimits(): void {
    const key = String(settings.cooks);
    if (!ladderEl || key === limitKey) return;
    limitKey = key;
    limits.clear();
    worker?.terminate();
    worker = null;
    let w: Worker;
    try {
      w = new Worker(new URL("../sim/morningLimits.worker.ts", import.meta.url), { type: "module" });
    } catch {
      ladderFailed();
      return;
    }
    worker = w;
    w.onerror = () => {
      w.terminate();
      if (worker === w) {
        worker = null;
        ladderFailed();
      }
    };
    w.onmessage = (e: MessageEvent<LimitScore>) => {
      if (e.data.key !== limitKey) return;
      limits.set(lkey(e.data.menu, e.data.ahead), e.data.meals);
      if (limits.size === 2 * MENUS) {
        w.terminate();
        if (worker === w) worker = null;
      }
      if (!ladderQueued) {
        ladderQueued = true;
        requestAnimationFrame(() => {
          ladderQueued = false;
          renderLadder();
          renderLimitTile();
        });
      }
    };
    const ask: LimitAsk = { key, cooks: settings.cooks };
    w.postMessage(ask);
  }
  const ladderFailed = () => {
    const wrap = ladderEl?.closest<HTMLElement>(".rs-map-wrap");
    if (wrap) wrap.hidden = true;
  };
  /**
   * A menu's limit: the worker's search, raised to the order on screen when that order makes the van. The search
   * halves the gap between an order that makes it and one that doesn't, and the plan is not always easier with fewer
   * meals, so it can stop short of an order that works.
   */
  const limitOf = (menu: number, ahead: boolean): number | undefined => {
    const v = limits.get(lkey(menu, ahead));
    if (v === undefined) return v;
    return menu === settings.menu && ahead === settings.ahead && madeIt(morning) && settings.meals > v ? settings.meals : v;
  };
  function renderLimitTile(): void {
    const v = limitOf(settings.menu, settings.ahead);
    const t = tiles.limit;
    if (t.v) t.v.textContent = v === undefined ? "…" : v.toLocaleString("en-US");
    if (t.sub) t.sub.textContent = v === undefined ? "working it out" : v >= settings.meals ? `${(v - settings.meals).toLocaleString("en-US")} to spare` : `${(settings.meals - v).toLocaleString("en-US")} too many`;
    if (t.delta) t.delta.hidden = true;
  }
  /** Menus by number, the written menu first and the one with every swap last; their labels. */
  const menuLabel = (menu: number) =>
    menuDishes(menu)
      .slice(0, 3)
      .map((d) => shortName(d.name))
      .join(" · ");
  function renderLadder(): void {
    if (!ladderEl) return;
    const W = Math.max(280, Math.round(ladderEl.clientWidth));
    // On a phone each menu's name sits over its bar; wider, to its left.
    const stack = W < 440;
    const rowH = stack ? 32 : 22;
    const barY = stack ? 15 : 4;
    const barH = stack ? 12 : rowH - 8;
    const L = stack ? 2 : Math.min(210, Math.round(W * 0.42));
    const R = W - 34;
    const T = 6;
    const H = T + MENUS * rowH + 20;
    const vals = [...limits.values()];
    const top = Math.max(settings.meals * 1.15, 800, ...vals);
    const nice = [800, 1000, 1200, 1500, 2000, 2500].find((v) => v >= top) ?? 2500;
    const X = (v: number) => L + (Math.min(v, nice) / nice) * (R - L);
    // Menus ordered by how many dishes are swapped, then by number.
    const order = Array.from({ length: MENUS }, (_, k) => k).sort((a, b) => popcount(a) - popcount(b) || a - b);
    let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" class="rs-svg mo-ladder">`;
    for (const v of [0, nice / 2, nice]) svg += `<line class="g v" x1="${X(v)}" x2="${X(v)}" y1="${T}" y2="${T + MENUS * rowH}"/><text class="xa" x="${X(v)}" y="${T + MENUS * rowH + 14}">${v.toLocaleString("en-US")}</text>`;
    order.forEach((menu, row) => {
      const y = T + row * rowH;
      const today = limitOf(menu, false);
      const ahead = limitOf(menu, true);
      const cur = menu === settings.menu;
      svg += `<g class="row${cur ? " on" : ""}" data-menu="${menu}"><rect class="hit" x="0" y="${y}" width="${W}" height="${rowH}"/>`;
      svg += stack ? `<text class="ml st" x="${L}" y="${y + 11}">${menuLabel(menu)}</text>` : `<text class="ml" x="${L - 8}" y="${y + rowH / 2 + 4}">${menuLabel(menu)}</text>`;
      if (ahead !== undefined && today !== undefined && ahead > today) svg += `<rect class="ahead" x="${X(0)}" y="${y + barY}" width="${X(ahead) - X(0)}" height="${barH}"/>`;
      if (today !== undefined) svg += `<rect class="today" x="${X(0)}" y="${y + barY}" width="${Math.max(1, X(today) - X(0))}" height="${barH}"/>`;
      const shown = settings.ahead ? ahead : today;
      if (shown !== undefined) svg += `<text class="v" x="${X(Math.max(shown, settings.ahead && ahead !== undefined ? ahead : 0)) + 5}" y="${y + barY + barH / 2 + 4}">${shown.toLocaleString("en-US")}</text>`;
      svg += `</g>`;
    });
    // This order's line: down the whole chart, or on a phone through each row's bar only, clear of the names.
    if (stack) for (let row = 0; row < MENUS; row++) svg += `<line class="order" x1="${X(settings.meals)}" x2="${X(settings.meals)}" y1="${T + row * rowH + barY - 2}" y2="${T + row * rowH + barY + barH + 2}"/>`;
    else svg += `<line class="order" x1="${X(settings.meals)}" x2="${X(settings.meals)}" y1="${T - 4}" y2="${T + MENUS * rowH + 2}"/>`;
    ladderEl.innerHTML = svg + "</svg>";
    if (ladderNote) {
      if (limits.size < 2 * MENUS) ladderNote.textContent = `Finding the limit for menu ${Math.min(2 * MENUS, limits.size + 1)} of ${2 * MENUS}…`;
      else {
        const best = order.reduce((a, b) => ((limitOf(b, false) ?? 0) > (limitOf(a, false) ?? 0) ? b : a), 0);
        const bt = limitOf(best, false) ?? 0;
        const ba = limitOf(best, true) ?? 0;
        const w0 = limitOf(0, false) ?? 0;
        ladderNote.textContent = `With ${settings.cooks} cooks, the menu as written gets ${w0.toLocaleString("en-US")} meals out by 10:30 at their best; ${menuLabel(best).toLowerCase()} gets ${bt.toLocaleString("en-US")}, and ${ba.toLocaleString("en-US")} with what keeps cooked the day before. Pick a row to plan that menu.`;
      }
    }
  }
  const popcount = (v: number) => (v & 1) + ((v >> 1) & 1) + ((v >> 2) & 1);
  ladderEl?.addEventListener("click", (e) => {
    const g = (e.target as Element).closest<SVGGElement>("[data-menu]");
    if (!g) return;
    settings.menu = Number(g.dataset.menu);
    update();
  });
  if (ladderEl) {
    let lastW = ladderEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(ladderEl.clientWidth - lastW) < 4) return;
      lastW = ladderEl.clientWidth;
      renderLadder();
    }).observe(ladderEl);
  }

  /* ---------------------------------------------------------------------------------------------------------- */
  /* Running it.                                                                                                 */
  /* ---------------------------------------------------------------------------------------------------------- */

  const cache = new Map<string, Morning>();
  const runFor = (s: MorningSettings) => {
    const k = JSON.stringify([s.meals, s.cooks, s.menu, s.ahead]);
    let r = cache.get(k);
    if (!r) {
      r = planMorning(s);
      cache.set(k, r);
    }
    return r;
  };
  let morning: Morning = runFor(settings);
  let endMin = DOOR;
  let loopS = 1;
  let tMin = 0;
  let loopT = 0;

  function update(): void {
    const frac = loopS > 0 ? tMin / Math.max(1, endMin) : 0;
    morning = runFor(settings);
    const base = settings.menu !== 0 || settings.ahead ? runFor({ ...settings, menu: 0, ahead: false }) : null;
    let last = 0;
    for (const j of morning.jobs) last = Math.max(last, j.t1);
    endMin = Math.max(DOOR, morning.depart, last) + 6;
    loopS = kitchen.plan(morning, settings.cooks, endMin, endMin * SEC_PER_MIN + END_HOLD_S);
    tMin = frac * endMin;
    loopT = tMin * SEC_PER_MIN;
    renderMenu();
    renderChart(morning, endMin);
    renderTiles(morning, base);
    renderStatus(morning);
    askLimits();
    renderLadder();
    renderLimitTile();
  }
  let pending = 0;
  const later = () => {
    if (pending) return;
    pending = window.setTimeout(() => {
      pending = 0;
      update();
    }, 90);
  };

  for (const b of crewBtns)
    b.addEventListener("click", () => {
      settings.cooks = Number(b.dataset.moCrew) as Crew;
      update();
    });
  for (const b of aheadBtns)
    b.addEventListener("click", () => {
      settings.ahead = b.dataset.moAhead === "yes";
      update();
    });
  for (const b of presetBtns)
    b.addEventListener("click", () => {
      settings.menu = b.dataset.moPreset === "written" ? 0 : MENUS - 1;
      update();
    });
  mealsIn?.addEventListener("input", () => {
    settings.meals = Number(mealsIn.value);
    if (mealsVal) mealsVal.textContent = `${settings.meals.toLocaleString("en-US")} meals`;
    later();
  });

  /* One moment of the morning: the kitchen, and the readouts when the minute changes. */
  function show(t: number, dt: number): void {
    const now = kitchen.show(t, dt);
    if (now.minute) {
      const gone = morning.finished && t >= morning.depart;
      if (hudClock) {
        hudClock.textContent = clockOf(t);
        hudClock.dataset.kind = t > DOOR && !gone ? "late" : "";
      }
      if (hudPans) hudPans.textContent = `${now.packed} of ${now.hot} hot pans packed`;
      if (hudNote) {
        hudNote.textContent = gone ? (morning.late > 0.5 ? `Van gone at ${clockOf(morning.depart)}` : "Van gone at 10:30") : t >= DOOR ? "The van is waiting" : `${now.burning} of ${BURNERS} burners lit`;
        hudNote.dataset.kind = gone ? (morning.late > 0.5 ? "late" : "held") : t >= DOOR ? "late" : "idle";
      }
    }
    if (cursor) {
      const x = xOf(t).toFixed(1);
      cursor.setAttribute("x1", x);
      cursor.setAttribute("x2", x);
    }
  }

  update();
  show(0, 0);
  if (chartEl) {
    let lastW = chartEl.clientWidth;
    new ResizeObserver(() => {
      if (Math.abs(chartEl.clientWidth - lastW) < 4) return;
      lastW = chartEl.clientWidth;
      renderChart(morning, endMin);
    }).observe(chartEl);
  }
  for (const k of cooks) void k.st.ready.then(() => void gate.refresh(scene));

  let inView = false;
  observeInView(stageEl, (v) => {
    inView = v;
  });
  let frozen: number | null = null;
  let last = performance.now();
  const fitClip = clipFitter(camera, scene, 0.05);
  function animate(now: number): void {
    requestAnimationFrame(animate);
    const dt = Math.max(0, Math.min(0.05, (now - last) / 1000));
    last = now;
    if (!inView || document.hidden || !gate.ready()) return;
    if (frozen === null) loopT = (loopT + dt) % loopS;
    else loopT = frozen * SEC_PER_MIN;
    tMin = Math.min(endMin, loopT / SEC_PER_MIN);
    kitchen.pose(loopT);
    show(tMin, dt);
    controls.update();
    fitClip();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);
  void gate.refresh(scene);

  // A hook for checking the page by script.
  Object.assign(stageEl, {
    morning: {
      info: () => ({
        mountMs: Math.round(mountMs),
        tris: renderer.info.render.triangles,
        calls: renderer.info.render.calls,
        loopS: +loopS.toFixed(2),
        routineS: cooks.slice(0, settings.cooks).map((k) => +(k.routine?.duration ?? 0).toFixed(2)),
        t: +tMin.toFixed(1),
        depart: clockOf(morning.depart),
        late: Math.round(morning.late),
        pastBest: morning.pastBest,
        onTimePastBest: morning.onTimePastBest,
        limits: Object.fromEntries(limits),
        acts: cooks.slice(0, settings.cooks).map((k) => {
          const st = k.routine?.apply(loopT);
          return st ? `${st.activity}${st.work ? `:${st.work}` : ""}` : "none";
        }),
        serve: clockOf(SERVE),
      }),
      hold(min: number | null) {
        frozen = min;
      },
      view(eye: [number, number, number], at: [number, number, number]) {
        userMoved = true;
        camera.position.set(...eye);
        controls.target.set(...at);
        controls.update();
      },
      set(s: Partial<MorningSettings>) {
        Object.assign(settings, s);
        update();
      },
    },
  });
  const mountMs = performance.now() - mountT0;
}
