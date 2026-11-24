import "@fontsource-variable/inter";
import "@fontsource-variable/jetbrains-mono";
import "./styles/main.css";

import { EXEMPLARS, computeKitchenTelemetry, type ExemplarId } from "./sim/exemplars";
import { mountBrigadeStage } from "./three/brigadeStage";
import { mountColdChainLab } from "./three/coldChainLab";
import { mountEquipmentTurntable } from "./three/equipmentTurntable";
import { mountHeroKitchen, type HeroController } from "./three/heroKitchen";
import { mountOvenLab } from "./three/ovenLab";
import { mountStoveLab } from "./three/stoveLab";

declare global {
  interface Window {
    __KITCHEN_LAB__?: {
      hero: HeroController;
      mountedChapters: string[];
      forceMountAll: () => void;
    };
  }
}

type FleetScenario = "gas" | "electric" | "summer";

function initScrollChrome(): void {
  const progressBar = document.querySelector<HTMLElement>("[data-progress-bar]");
  const navLinks = Array.from(document.querySelectorAll<HTMLAnchorElement>(".chapter-nav a[data-spy]"));
  const sections = navLinks
    .map((a) => ({ link: a, el: document.getElementById(a.dataset.spy ?? "") }))
    .filter((x): x is { link: HTMLAnchorElement; el: HTMLElement } => x.el !== null);

  const onScroll = () => {
    const y = window.scrollY;
    document.documentElement.classList.toggle("scrolled", y > 24);

    if (progressBar) {
      const maxScroll = Math.max(1, document.documentElement.scrollHeight - window.innerHeight);
      const frac = Math.max(0, Math.min(1, y / maxScroll));
      progressBar.style.transform = `scaleX(${frac.toFixed(4)})`;
    }

    // Active scroll-spy chapter
    let current = sections[0];
    const probeY = y + window.innerHeight * 0.34;
    for (const s of sections) {
      if (s.el.offsetTop <= probeY) current = s;
    }
    for (const s of sections) {
      if (s === current) s.link.setAttribute("aria-current", "location");
      else s.link.removeAttribute("aria-current");
    }
  };

  window.addEventListener("scroll", onScroll, { passive: true });
  window.addEventListener("resize", onScroll, { passive: true });
  onScroll();

  // Close mobile menu drawer cleanly when any link inside it is tapped
  const mobileMenu = document.querySelector<HTMLDetailsElement>(".mobile-menu");
  mobileMenu?.querySelectorAll("a").forEach((a) => {
    a.addEventListener("click", () => {
      mobileMenu.open = false;
    });
  });
}

function initFleetMatrix(hero: HeroController): void {
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
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td><strong>${spec.shortName}</strong></td>
        <td class="num">${spec.areaSqFt.toLocaleString()} ft²</td>
        <td class="num">${spec.fixtures} · ${spec.electricalService}</td>
        <td class="num">${m.thermalKw.toFixed(1)} kW</td>
        <td class="num">${Math.round(m.exhaustCfm).toLocaleString()} CFM</td>
        <td class="num">${m.hvacCoolingKw.toFixed(1)} kW</td>
        <td class="num">${m.thermalComfortC.toFixed(1)} °C</td>
        <td class="num">${m.minAisleIn.toFixed(0)}" (${m.minAisleM.toFixed(2)} m)</td>
        <td><button type="button" class="chip-btn chip-subtle" data-load-hero="${spec.id}">Load in 3D ↑</button></td>
      `;
      tbody.appendChild(tr);
    }

    tbody.querySelectorAll<HTMLButtonElement>("[data-load-hero]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const id = btn.dataset.loadHero as ExemplarId;
        void hero.setExemplar(id);
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
  const syncSlider = (el: HTMLInputElement) => {
    const min = Number(el.min || 0);
    const max = Number(el.max || 100);
    const val = Number(el.value || 0);
    const pct = max > min ? ((val - min) / (max - min)) * 100 : 0;
    el.style.setProperty("--pct", `${Math.max(0, Math.min(100, pct)).toFixed(2)}%`);
  };
  document.querySelectorAll<HTMLInputElement>('input[type="range"]').forEach((el) => {
    syncSlider(el);
    el.addEventListener("input", () => syncSlider(el));
    el.addEventListener("change", () => syncSlider(el));
  });
}

function bootstrap(): void {
  initScrollChrome();
  initRangeTrackFills();

  const heroRoot = document.getElementById("hero");
  if (!heroRoot) return;
  const hero = mountHeroKitchen(heroRoot);
  initFleetMatrix(hero);

  const mountedChapters: string[] = ["hero"];
  const pendingMounts: Array<{ id: string; el: HTMLElement; mount: (el: HTMLElement) => void; done: boolean }> = [
    { id: "equipment", el: document.getElementById("equipment")!, mount: mountEquipmentTurntable, done: false },
    { id: "stove", el: document.getElementById("stove")!, mount: mountStoveLab, done: false },
    { id: "oven", el: document.getElementById("oven")!, mount: mountOvenLab, done: false },
    { id: "cold", el: document.getElementById("cold")!, mount: mountColdChainLab, done: false },
    { id: "brigade", el: document.getElementById("brigade")!, mount: mountBrigadeStage, done: false },
  ].filter((entry) => entry.el !== null);

  const mountOne = (entry: (typeof pendingMounts)[number]) => {
    if (entry.done) return;
    entry.done = true;
    mountedChapters.push(entry.id);
    entry.mount(entry.el);
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
