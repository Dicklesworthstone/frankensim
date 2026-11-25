import "@fontsource-variable/inter";
import "@fontsource-variable/jetbrains-mono";
import "./styles/main.css";

import { EXEMPLARS, computeKitchenTelemetry, type ExemplarId } from "./sim/exemplars";
import { mountBrigadeStage } from "./three/brigadeStage";
import { mountColdChainLab } from "./three/coldChainLab";
import { mountEquipmentTurntable } from "./three/equipmentTurntable";
import { mountHeroKitchen, type HeroController } from "./three/heroKitchen";
import { OVEN_BOULE_DETAIL, mountOvenLab } from "./three/ovenLab";
import { STOVE_STEAK, mountStoveLab } from "./three/stoveLab";
import { prefetchSteakTissue } from "./three/models/steakPrefetch";
import { prefetchBouleLayers } from "./three/models/boulePrefetch";
import { mountWeekLab } from "./three/weekLab";
import { mountFloorLab } from "./three/floorLab";
import { mountRushLab } from "./three/rushLab";

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
  initRangeTrackFills();
  // The stove chapter's ribeye texture is the heaviest single build on the page: start it on another thread now,
  // while the hero kitchen builds, and mount the stove once it has arrived. The oven's loaf crust is the next
  // heaviest; the oven waits for it only to build the loaf, after its own mount.
  const steakReady = prefetchSteakTissue(STOVE_STEAK);
  void prefetchBouleLayers(OVEN_BOULE_DETAIL);

  const heroRoot = document.getElementById("hero");
  const hero = heroRoot ? tryMount("hero", () => mountHeroKitchen(heroRoot)) : null;
  initFleetMatrix(hero);

  const mountedChapters: string[] = hero ? ["hero"] : [];
  const pendingMounts: Array<{ id: string; el: HTMLElement; mount: (el: HTMLElement) => void; done: boolean; after?: Promise<void> }> = [
    { id: "equipment", el: document.getElementById("equipment")!, mount: mountEquipmentTurntable, done: false },
    { id: "stove", el: document.getElementById("stove")!, mount: mountStoveLab, done: false, after: steakReady },
    { id: "oven", el: document.getElementById("oven")!, mount: mountOvenLab, done: false },
    { id: "cold", el: document.getElementById("cold")!, mount: mountColdChainLab, done: false },
    { id: "brigade", el: document.getElementById("brigade")!, mount: mountBrigadeStage, done: false },
    { id: "week", el: document.getElementById("week")!, mount: mountWeekLab, done: false },
    { id: "floor", el: document.getElementById("floor")!, mount: mountFloorLab, done: false },
    { id: "rush", el: document.getElementById("rush")!, mount: mountRushLab, done: false },
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
  const mountOne = (entry: (typeof pendingMounts)[number]) => {
    if (entry.done) return;
    entry.done = true;
    const go = () => {
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
