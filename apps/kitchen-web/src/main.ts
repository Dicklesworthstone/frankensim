import "@fontsource-variable/inter";
import "@fontsource-variable/jetbrains-mono";
import "./styles/main.css";

import { EXEMPLARS, computeKitchenTelemetry, type ExemplarId } from "./sim/exemplars";
import { BRIGADE_PIECES, mountBrigadeStage } from "./three/brigadeStage";
import { coldPieces, mountColdChainLab } from "./three/coldChainLab";
import { EQUIPMENT_PIECES, mountEquipmentTurntable } from "./three/equipmentTurntable";
import { mountHeroKitchen, type HeroController } from "./three/heroKitchen";
import { OVEN_BOULE_DETAIL, OVEN_PIECES, mountOvenLab } from "./three/ovenLab";
import { STOVE_PIECES, STOVE_STEAK, mountStoveLab } from "./three/stoveLab";
import { BOIL_PIECES, mountBoilLab } from "./three/boilLab";
import { coolPieces, mountCoolLab } from "./three/coolLab";
import { holdPieces, mountHoldLab } from "./three/holdLab";
import { hoodPieces, mountHoodLab } from "./three/hoodLab";
import { mountWalkInLab, walkInPieces } from "./three/walkInLab";
import { prefetchSteakTissue } from "./three/models/steakPrefetch";
import { prefetchBouleLayers } from "./three/models/boulePrefetch";
import { prefetchFoodTexels } from "./three/models/foodPrefetch";
import { buildAhead, type Piece } from "./three/models/all";
import { WEEK_PIECES, mountWeekLab } from "./three/weekLab";
import { floorPieces, mountFloorLab } from "./three/floorLab";
import { RUSH_PIECES, mountRushLab } from "./three/rushLab";
import { MORNING_PIECES, mountMorningLab } from "./three/morningLab";
import { PIVOT_PIECES, mountPivotLab } from "./three/pivotLab";
import { CROWD_PIECES, mountCrowdLab } from "./three/crowdLab";
import { PAIR_PIECES, mountPairLab } from "./three/pairLab";
import { KIT_PIECES, mountKitLab } from "./three/kitLab";
import { GLASSES_PIECES, mountGlassesLab } from "./three/glassesLab";
import { TRIP_PIECES, mountTripLab } from "./three/tripLab";

declare global {
  interface Window {
    __KITCHEN_LAB__?: {
      hero: HeroController | null;
      mountedChapters: string[];
      forceMountAll: () => void;
    };
  }
}

type FleetScenario = "gas" | "electric" | "summer";

function initScrollChrome(): void {
  const progressBar = document.querySelector<HTMLElement>("[data-progress-bar]");
  const navLinks = Array.from(document.querySelectorAll<HTMLAnchorElement>(".chapters a[data-spy]"));
  const sections = navLinks
    .map((a) => ({ link: a, el: document.getElementById(a.dataset.spy ?? "") }))
    .filter((x): x is { link: HTMLAnchorElement; el: HTMLElement } => x.el !== null);

  // Once per frame at most, and every layout read before any write, so scrolling never forces a synchronous
  // layout of the whole page.
  let queued = false;
  let current: (typeof sections)[number] | undefined;
  const update = () => {
    queued = false;
    const y = window.scrollY;
    const maxScroll = Math.max(1, document.documentElement.scrollHeight - window.innerHeight);
    // The chapter being read: the last one whose top has passed a third of the way down the screen.
    // Above the first chapter (the hero and the "why" band) no link is lit.
    const probe = window.innerHeight * 0.34;
    let next: (typeof sections)[number] | undefined;
    for (const s of sections) {
      if (s.el.getBoundingClientRect().top <= probe) next = s;
    }

    document.documentElement.classList.toggle("scrolled", y > 24);
    if (progressBar) progressBar.style.transform = `scaleX(${Math.max(0, Math.min(1, y / maxScroll)).toFixed(4)})`;
    if (next !== current) {
      current?.link.removeAttribute("aria-current");
      next?.link.setAttribute("aria-current", "location");
      current = next;
    }
  };
  const onScroll = () => {
    if (queued) return;
    queued = true;
    requestAnimationFrame(update);
  };

  window.addEventListener("scroll", onScroll, { passive: true });
  window.addEventListener("resize", onScroll, { passive: true });
  update();

  // Close the phone menu when one of its links is tapped.
  const menu = document.querySelector<HTMLDetailsElement>(".menu");
  menu?.querySelectorAll("a").forEach((a) => {
    a.addEventListener("click", () => {
      menu.open = false;
    });
  });
}

/** Headings, panels and cards fade up once as they come into view. Without the observer they are simply shown. */
function initReveal(): void {
  const items = Array.from(document.querySelectorAll<HTMLElement>(".reveal"));
  if (!("IntersectionObserver" in window)) {
    for (const el of items) el.classList.add("seen");
    return;
  }
  const io = new IntersectionObserver(
    (entries) => {
      for (const e of entries) {
        if (!e.isIntersecting) continue;
        e.target.classList.add("seen");
        io.unobserve(e.target);
      }
    },
    { rootMargin: "0px 0px -8% 0px", threshold: 0.04 },
  );
  for (const el of items) io.observe(el);
}

function initFleetMatrix(hero: HeroController | null): void {
  const tbody = document.querySelector<HTMLElement>("[data-fleet-body]");
  const scenarioBtns = Array.from(document.querySelectorAll<HTMLButtonElement>("[data-fleet-scenario]"));
  if (!tbody) return;

  let scenario: FleetScenario = "gas";

  const renderRows = () => {
    const elecPct = scenario === "electric" ? 100 : 15;
    const ambC = scenario === "summer" ? 38 : 26;
    const loadMult = scenario === "summer" ? 1.25 : 1.0;

    tbody.innerHTML = "";
    for (const spec of EXEMPLARS) {
      const m = computeKitchenTelemetry(spec, spec.nominalOrdersPerHr * loadMult, elecPct, ambC);
      const amps = spec.electricalService.split(" ")[0] ?? spec.electricalService;
      const tr = document.createElement("tr");
      // data-label names each cell for the phone layout, where every kitchen becomes a card.
      tr.innerHTML = `
        <td class="name"><b>${spec.shortName}</b><small>${spec.areaSqFt.toLocaleString()} ft² · ${amps} service</small></td>
        <td class="num" data-label="Equipment">${spec.fixtures} pieces</td>
        <td class="num" data-label="Cooking load">${m.thermalKw.toFixed(1)} kW</td>
        <td class="num" data-label="Exhaust">${Math.round(m.exhaustCfm).toLocaleString()} cfm</td>
        <td class="num" data-label="Cooling">${m.hvacCoolingKw.toFixed(1)} kW</td>
        <td class="num" data-label="At the line">${m.thermalComfortC.toFixed(1)} °C</td>
        <td class="num" data-label="Tightest aisle">${m.minAisleIn.toFixed(0)} in · ${m.minAisleM.toFixed(2)} m</td>
        <td class="act"><button type="button" class="btn sm" data-load-hero="${spec.id}">View in 3D<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M12 19V5M6 11l6-6 6 6"/></svg></button></td>
      `;
      tbody.appendChild(tr);
    }

    tbody.querySelectorAll<HTMLButtonElement>("[data-load-hero]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const id = btn.dataset.loadHero as ExemplarId;
        if (hero) void hero.setExemplar(id);
        document.getElementById("hero")?.scrollIntoView({ behavior: "smooth" });
      });
    });
  };

  for (const btn of scenarioBtns) {
    btn.addEventListener("click", () => {
      scenario = btn.dataset.fleetScenario as FleetScenario;
      for (const b of scenarioBtns) b.setAttribute("aria-pressed", String(b.dataset.fleetScenario === scenario));
      renderRows();
    });
  }
  renderRows();
}

function initRangeTrackFills(): void {
  const sliders = Array.from(document.querySelectorAll<HTMLInputElement>('input[type="range"]'));
  const syncSlider = (el: HTMLInputElement) => {
    const min = Number(el.min || 0);
    const max = Number(el.max || 100);
    const val = Number(el.value || 0);
    const pct = max > min ? ((val - min) / (max - min)) * 100 : 0;
    el.style.setProperty("--p", `${Math.max(0, Math.min(100, pct)).toFixed(2)}%`);
  };
  for (const el of sliders) {
    syncSlider(el);
    el.addEventListener("input", () => syncSlider(el));
    el.addEventListener("change", () => syncSlider(el));
  }
  // Presets and reset buttons set slider values without an input event; catch up after any click.
  document.addEventListener("click", () => requestAnimationFrame(() => sliders.forEach(syncSlider)));
}

/** A link to a figure's method (a chapter head's "How it's computed", or the page's address) opens it where it is
 * folded under the figure. */
function initHowLinks(): void {
  const open = (id: string) => {
    const el = id ? document.getElementById(id) : null;
    const d = el?.matches(".fig-foot.how") ? el.querySelector<HTMLDetailsElement>(":scope > details") : null;
    if (d) d.open = true;
  };
  document.querySelectorAll<HTMLAnchorElement>('a[href^="#"]').forEach((a) => a.addEventListener("click", () => open(a.hash.slice(1))));
  const fromAddress = () => {
    // A hand-edited address can carry a broken escape: it opens nothing rather than stopping the page's start.
    try {
      if (location.hash) open(decodeURIComponent(location.hash.slice(1)));
    } catch {
      /* not an id on this page */
    }
  };
  window.addEventListener("hashchange", fromAddress);
  fromAddress();
}

/**
 * A jump to a chapter (the nav, a link, an address ending in #) passes figures that mount on the way and grow above
 * its target as they do (their charts, maps and menus are built when they mount), and a smooth scroll ends where its
 * target was when it set off: on a phone the reader landed a screen or more short of the chapter. For a few seconds
 * after a jump, until the reader scrolls themselves, each change in a chapter's or a figure's size puts the target
 * back where the jump puts it.
 */
function initJumpPin(): void {
  if (!("ResizeObserver" in window)) return;
  let target: HTMLElement | null = null;
  let since = 0;
  let until = 0;
  const pin = (id: string) => {
    let el: HTMLElement | null = null;
    try {
      el = id ? document.getElementById(decodeURIComponent(id)) : null;
    } catch {
      /* a broken escape in the address: no target */
    }
    target = el;
    since = performance.now();
    until = since + 2500;
  };
  // The reader's own scroll ends it. (A tap or a key on a link comes before its click, which pins again.)
  for (const ev of ["wheel", "touchstart", "pointerdown", "keydown"]) window.addEventListener(ev, () => (target = null), { passive: true, capture: true });
  document.addEventListener("click", (e) => {
    const a = (e.target as Element | null)?.closest?.<HTMLAnchorElement>('a[href^="#"]');
    if (a && a.hash.length > 1) pin(a.hash.slice(1));
  });
  if (location.hash.length > 1) pin(location.hash.slice(1));
  const ro = new ResizeObserver(() => {
    if (!target) return;
    const now = performance.now();
    if (now > until) {
      target = null;
      return;
    }
    // A figure that mounts late still grows under the pin: each change holds it a little longer, to 12 s at most.
    until = Math.min(since + 12000, now + 2500);
    target.scrollIntoView({ block: "start" });
  });
  for (const el of document.querySelectorAll("section[id], figure.fig")) ro.observe(el);
}

/** Mount one 3D chapter. A chapter that fails (no WebGL, a lost context) leaves its section's text and controls in
 * place and never stops the rest of the page from starting. */
function tryMount<T>(id: string, mount: () => T): T | null {
  try {
    return mount();
  } catch (err) {
    console.error(`[kitchen-lab] the ${id} 3D view could not start`, err);
    const section = document.getElementById(id);
    section?.classList.add("no-3d");
    // Status lines are written by the chapter's loop, which never started.
    section?.querySelectorAll<HTMLElement>(".status").forEach((s) => {
      s.dataset.state = "warn";
      s.textContent = "This model runs with its 3D view, which could not start.";
    });
    const hud = section?.querySelector<HTMLElement>("[data-hero-status]");
    if (hud) hud.textContent = "A still: the live kitchen could not start";
    return null;
  }
}

function bootstrap(): void {
  // The inline script in index.html hides .reveal blocks until they scroll in; it shows them again if this module
  // never runs.
  document.documentElement.classList.add("booted");
  initScrollChrome();
  initReveal();
  // A browser that keeps form values over a reload (Firefox) or a history step puts a control where the chapter
  // behind it does not start: a slider elsewhere, a search the equipment's chips do not show. Every chapter starts at
  // its HTML default, so its controls do too.
  for (const el of document.querySelectorAll<HTMLInputElement>("input")) el.value = el.defaultValue;
  initRangeTrackFills();
  initHowLinks();
  initJumpPin();
  // The stove chapter's ribeye texture is the heaviest single build on the page: start it on another thread now,
  // while the hero kitchen builds, and mount the stove once it has arrived. The oven's loaf crust is the next
  // heaviest; the oven waits for it only to build the loaf, after its own mount.
  const steakReady = prefetchSteakTissue(STOVE_STEAK);
  void prefetchBouleLayers(OVEN_BOULE_DETAIL);
  // The pan foods' textures are the heaviest paint in the hero kitchen's build (its first hotel pan), the wine label's
  // paper in the stove's: both start on another thread now too.
  prefetchFoodTexels();

  const heroRoot = document.getElementById("hero");
  const hero = heroRoot ? tryMount("hero", () => mountHeroKitchen(heroRoot)) : null;
  initFleetMatrix(hero);

  const mountedChapters: string[] = hero ? ["hero"] : [];
  /** A chapter to mount: `done` once asked for, `ran` once its mount has run; `pieces` what it builds (see below). */
  type Pending = {
    id: string;
    el: HTMLElement;
    mount: (el: HTMLElement) => void;
    done: boolean;
    ran?: boolean;
    after?: Promise<void>;
    pieces?: readonly Piece[] | (() => Promise<readonly Piece[]>);
  };
  const pendingMounts: Pending[] = [
    { id: "equipment", el: document.getElementById("equipment")!, mount: mountEquipmentTurntable, done: false, pieces: EQUIPMENT_PIECES },
    { id: "stove", el: document.getElementById("stove")!, mount: mountStoveLab, done: false, after: steakReady, pieces: () => steakReady.then(() => STOVE_PIECES) },
    { id: "boil", el: document.getElementById("boil")!, mount: mountBoilLab, done: false, pieces: BOIL_PIECES },
    { id: "oven", el: document.getElementById("oven")!, mount: mountOvenLab, done: false, pieces: OVEN_PIECES },
    { id: "cold", el: document.getElementById("cold")!, mount: mountColdChainLab, done: false, pieces: coldPieces() },
    { id: "walkin-week", el: document.getElementById("walkin-week")!, mount: mountWalkInLab, done: false, pieces: walkInPieces() },
    { id: "brigade", el: document.getElementById("brigade")!, mount: mountBrigadeStage, done: false, pieces: BRIGADE_PIECES },
    { id: "pivot", el: document.getElementById("pivot")!, mount: mountPivotLab, done: false, pieces: PIVOT_PIECES },
    { id: "week", el: document.getElementById("week")!, mount: mountWeekLab, done: false, pieces: WEEK_PIECES },
    { id: "cool", el: document.getElementById("cool")!, mount: mountCoolLab, done: false, pieces: coolPieces() },
    { id: "hold", el: document.getElementById("hold")!, mount: mountHoldLab, done: false, pieces: holdPieces() },
    { id: "floor", el: document.getElementById("floor")!, mount: mountFloorLab, done: false, pieces: floorPieces },
    { id: "hood", el: document.getElementById("hood")!, mount: mountHoodLab, done: false, after: steakReady, pieces: () => steakReady.then(() => hoodPieces()) },
    { id: "rush", el: document.getElementById("rush")!, mount: mountRushLab, done: false, pieces: RUSH_PIECES },
    { id: "glasses", el: document.getElementById("glasses")!, mount: mountGlassesLab, done: false, pieces: GLASSES_PIECES },
    { id: "morning", el: document.getElementById("morning")!, mount: mountMorningLab, done: false, pieces: MORNING_PIECES },
    { id: "crowd", el: document.getElementById("crowd")!, mount: mountCrowdLab, done: false, pieces: CROWD_PIECES },
    { id: "pair", el: document.getElementById("pair")!, mount: mountPairLab, done: false, pieces: PAIR_PIECES },
    { id: "kit", el: document.getElementById("kit")!, mount: mountKitLab, done: false, pieces: KIT_PIECES },
    { id: "trip", el: document.getElementById("trip")!, mount: mountTripLab, done: false, pieces: TRIP_PIECES },
  ].filter((entry) => entry.el !== null);

  // Each chapter mounts in a task of its own, in the order asked. Two that come into range together, or the stove
  // (whose steak arrives on its own time) behind another, would otherwise build in one long task with no frame
  // drawn between them.
  const queue: Array<() => void> = [];
  let draining = false;
  const drain = () => {
    queue.shift()?.();
    if (queue.length) setTimeout(drain, 0);
    else draining = false;
  };
  const enqueue = (job: () => void) => {
    queue.push(job);
    if (!draining) {
      draining = true;
      setTimeout(drain, 0);
    }
  };
  const mountOne = (entry: Pending) => {
    if (entry.done) return;
    entry.done = true;
    const go = () => {
      entry.ran = true;
      if (tryMount(entry.id, () => (entry.mount(entry.el), true))) mountedChapters.push(entry.id);
    };
    if (entry.after) void entry.after.then(() => enqueue(go));
    else enqueue(go);
  };

  const lazyObserver = new IntersectionObserver(
    (entries) => {
      for (const ioEntry of entries) {
        if (!ioEntry.isIntersecting) continue;
        const target = pendingMounts.find((p) => p.el === ioEntry.target);
        if (target) {
          mountOne(target);
          lazyObserver.unobserve(ioEntry.target);
        }
      }
    },
    { rootMargin: "420px 0px" },
  );

  for (const entry of pendingMounts) {
    lazyObserver.observe(entry.el);
  }

  // Building its pieces of equipment is most of the task that mounts a chapter. Once the hero's first kitchen is up,
  // each chapter within three screens has its pieces built ahead, one to an idle moment, for its mount to take; a
  // chapter that mounts first builds whatever is not ready itself. (The floor's pieces come from the plan it fetches.)
  const aheadQueue: Array<{ entry: Pending; piece: Piece }> = [];
  let aheadBusy = false;
  const idle = (cb: (d: IdleDeadline) => void) => {
    if ("requestIdleCallback" in window) requestIdleCallback(cb, { timeout: 250 });
    else setTimeout(() => cb({ didTimeout: true, timeRemaining: () => 0 }), 50);
  };
  const buildSome = (d: IdleDeadline) => {
    let built = 0;
    // While the browser has 10 ms to spare before its next frame. A chapter on screen that animates leaves little of
    // that on a slow phone, so after a quarter of a second's wait one piece is built anyway: a frame or two late now
    // and then while the reader is mid-chapter, against a second's freeze on reaching the next.
    while (aheadQueue.length && (d.didTimeout ? built === 0 : d.timeRemaining() >= 10)) {
      const { entry, piece } = aheadQueue.shift()!;
      if (entry.ran) continue;
      try {
        buildAhead(piece);
      } catch {
        // Left to the chapter's mount, which builds it itself and reports any failure in its own place.
      }
      built++;
    }
    if (aheadQueue.length) idle(buildSome);
    else aheadBusy = false;
  };
  const queueAhead = (entry: Pending, pieces: readonly Piece[]) => {
    if (entry.ran) return;
    for (const piece of pieces) aheadQueue.push({ entry, piece });
    if (aheadQueue.length && !aheadBusy) {
      aheadBusy = true;
      idle(buildSome);
    }
  };
  const aheadObserver = new IntersectionObserver(
    (entries) => {
      for (const ioEntry of entries) {
        if (!ioEntry.isIntersecting) continue;
        aheadObserver.unobserve(ioEntry.target);
        const entry = pendingMounts.find((p) => p.el === ioEntry.target);
        const pieces = entry?.pieces;
        if (!entry || !pieces) continue;
        if (typeof pieces === "function") void pieces().then((list) => queueAhead(entry, list), () => undefined);
        else queueAhead(entry, pieces);
      }
    },
    { rootMargin: "300% 0px" },
  );
  void (hero?.ready ?? Promise.resolve()).then(() => {
    for (const entry of pendingMounts) if (entry.pieces) aheadObserver.observe(entry.el);
  });

  window.__KITCHEN_LAB__ = {
    hero,
    mountedChapters,
    forceMountAll: () => {
      for (const entry of pendingMounts) {
        mountOne(entry);
        lazyObserver.unobserve(entry.el);
      }
    },
  };
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", bootstrap);
} else {
  bootstrap();
}
