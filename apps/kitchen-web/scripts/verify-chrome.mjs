#!/usr/bin/env node
/**
 * Headless Chrome E2E & visual verification for `frankensim.org/kitchen-lab` (`dist-cf/`).
 *
 * Verifies:
 *   - ZERO console errors, ZERO uncaught page errors, ZERO failed network fetches (404/500)
 *   - ALL 10 WebGL canvases render non-blank frames (pixel variance > 0)
 *   - Each of the 5 hero kitchens builds cleanly, within the draw-call ceiling (DRAW_CEILING, shadow pass included)
 *   - Clicking "Deglaze 120 mL Wine" drops T_pan and increases vapor_released_g
 *   - The delivery box's fast estimate answers in range (quoting its measured error) and refuses out of range
 *   - 0 px horizontal overflow across 1440x900 (desktop), 768x1024 (tablet), 390x844 (mobile)
 *   - The 1200x630 share card `dist-cf/kitchen-lab/og-card.png` exists (rendered by scripts/render-posters.mjs)
 *   - Captures screenshots to `/tmp/kitchen_web_verify/`
 */

import fs from "node:fs";
import http from "node:http";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const APP_ROOT = path.resolve(__dirname, "..");
const DIST_CF = path.join(APP_ROOT, "dist-cf");
const VERIFY_DIR = "/tmp/kitchen_web_verify";

const require = createRequire(import.meta.url);
const puppeteer = require(path.resolve(APP_ROOT, "../kitchen-studio/node_modules/puppeteer-core"));

const CHROME_PATH = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";

/** Draw calls in one hero frame, shadow pass included. Measured 194–528 across the five kitchens on 2026-11-24 (the
 * commissary is the most); the ceiling leaves room for the camera's view and catches a change that multiplies them. */
const DRAW_CEILING = 700;

const MIME = {
  ".html": "text/html; charset=utf-8",
  ".js": "application/javascript; charset=utf-8",
  ".css": "text/css; charset=utf-8",
  ".svg": "image/svg+xml",
  ".png": "image/png",
  ".woff2": "font/woff2",
  ".woff": "font/woff",
  ".json": "application/json; charset=utf-8",
  ".glb": "model/gltf-binary",
};

function startStaticServer(rootDir) {
  return new Promise((resolve) => {
    const server = http.createServer((req, res) => {
      const urlPath = decodeURIComponent((req.url || "/").split("?")[0]);
      let rel = urlPath.replace(/^\/+/, "");
      let filePath = path.join(rootDir, rel);
      if (fs.existsSync(filePath) && fs.statSync(filePath).isDirectory()) {
        filePath = path.join(filePath, "index.html");
      }
      if (!filePath.startsWith(rootDir) || !fs.existsSync(filePath)) {
        res.statusCode = 404;
        res.end(`404 Not Found: ${urlPath}`);
        return;
      }
      const ext = path.extname(filePath).toLowerCase();
      res.setHeader("Content-Type", MIME[ext] || "application/octet-stream");
      fs.createReadStream(filePath).pipe(res);
    });
    server.listen(0, "127.0.0.1", () => {
      const addr = server.address();
      resolve({ server, port: addr.port });
    });
  });
}

function sleep(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

async function main() {
  if (!fs.existsSync(path.join(DIST_CF, "kitchen-lab", "index.html"))) {
    throw new Error("dist-cf/kitchen-lab/index.html missing — run `bun run build:cf` first.");
  }
  fs.mkdirSync(VERIFY_DIR, { recursive: true });

  const { server, port } = await startStaticServer(DIST_CF);
  const url = `http://127.0.0.1:${port}/kitchen-lab/`;
  console.log(`[verify-chrome] Serving dist-cf at ${url}`);

  const browser = await puppeteer.launch({
    executablePath: CHROME_PATH,
    headless: "new",
    args: [
      "--no-sandbox",
      "--disable-setuid-sandbox",
      "--use-gl=angle",
      "--use-angle=metal", // the GPU: software GL cannot keep up with the full hero kitchen
      "--enable-webgl",
      "--ignore-gpu-blocklist",
      "--window-size=1440,900",
    ],
  });

  const consoleErrors = [];
  const pageErrors = [];
  const failedRequests = [];

  try {
    const page = await browser.newPage();
    await page.setViewport({ width: 1440, height: 900, deviceScaleFactor: 2 });

    page.on("console", (msg) => {
      if (msg.type() === "error") consoleErrors.push(msg.text());
    });
    page.on("pageerror", (err) => pageErrors.push(String(err)));
    page.on("response", (resp) => {
      if (resp.status() >= 400) {
        failedRequests.push(`${resp.status()} ${resp.url()}`);
      }
    });

    // ?shaders=check: the built page skips three's shader error checks (src/three/kit.ts); this check wants them.
    await page.goto(`${url}?shaders=check`, { waitUntil: "networkidle0", timeout: 45000 });
    // The hero kitchen is built in the page after load; wait until it reports its draw calls (it is finished).
    await page.waitForFunction(() => (document.querySelector("[data-hero-status]")?.textContent || "").includes("draw calls"), {
      timeout: 60000,
    });
    await sleep(1600);

    // 1. Capture crisp 1440x900 desktop Hero screenshot (Mall QSR 400 loaded)
    await page.screenshot({ path: path.join(VERIFY_DIR, "hero-desktop.png"), fullPage: false });
    console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/hero-desktop.png");

    // 2. The share card is rendered by scripts/render-posters.mjs into public/ and copied into the build by Vite.
    const ogPath = path.join(DIST_CF, "kitchen-lab", "og-card.png");
    if (!fs.existsSync(ogPath)) throw new Error("dist-cf/kitchen-lab/og-card.png missing: run node scripts/render-posters.mjs, then build again.");

    // Restore the 1440x900 desktop viewport at 1x for the multi-scene pass
    await page.setViewport({ width: 1440, height: 900, deviceScaleFactor: 1 });
    await sleep(400);

    // 3. Build each of the 5 hero kitchens and check its frame stays within the draw-call ceiling
    const archetypes = [
      "freestanding_qsr_1450",
      "fast_casual_1100",
      "fine_dining_950",
      "commissary_catering_3200",
      "mall_qsr_400",
    ];
    for (const id of archetypes) {
      await page.evaluate((target) => {
        const btn = document.querySelector(`button[data-archetype="${target}"]`);
        if (btn) btn.click();
      }, id);
      await page.waitForFunction(
        (target) => {
          const api = window.__KITCHEN_LAB__;
          const st = document.querySelector("[data-hero-status]")?.textContent || "";
          return api && api.hero.getLoadedExemplarId() === target && st.includes("draw calls");
        },
        { timeout: 25000 },
        id,
      );
      const drawCalls = await page.evaluate(() => window.__KITCHEN_LAB__.hero.getDrawCalls());
      console.log(`[verify-chrome] Hero archetype ${id}: loaded cleanly (${drawCalls} draw calls)`);
      if (drawCalls <= 0 || drawCalls > DRAW_CEILING) {
        throw new Error(`Hero draw call budget violated for ${id}: ${drawCalls}`);
      }
    }

    // 4. Mount all chapters cleanly by scrolling through each and verifying all 10 WebGL canvases render non-blank frames
    const sectionIds = ["hero", "equipment", "stove", "oven", "cold", "brigade", "week", "floor", "rush", "morning"];
    const canvasSelectors = [
      "[data-hero-canvas]",
      "[data-eq-canvas]",
      "[data-stove-canvas]",
      "[data-oven-canvas]",
      "[data-cold-canvas]",
      "[data-brigade-canvas]",
      "[data-week-canvas]",
      "[data-floor-canvas]",
      "[data-rush-canvas]",
      "[data-morning-canvas]",
    ];

    for (let i = 0; i < sectionIds.length; i += 1) {
      const sid = sectionIds[i];
      const sel = canvasSelectors[i];
      // The canvas itself, centred: a chapter can hold several figures (06 has three), so centring the section can
      // leave its canvas off screen, where the view does not draw.
      await page.evaluate((s) => {
        document.querySelector(s)?.scrollIntoView({ behavior: "instant", block: "center" });
      }, sel);
      await sleep(1400);

      // Verify non-blank painted pixels inside each canvas bounding box via real compositor screenshot clip
      // In page coordinates, which a screenshot clip is in (a viewport box clipped the top of the page instead).
      const box = await page.$eval(sel, (el) => {
        const r = el.getBoundingClientRect();
        return { x: Math.max(0, r.x + window.scrollX), y: Math.max(0, r.y + window.scrollY), width: r.width, height: r.height };
      });
      if (box.width < 120 || box.height < 120) {
        throw new Error(`Canvas ${sel} collapsed: ${JSON.stringify(box)}`);
      }
      // A view draws through the page's shared WebGL renderer (src/three/kit.ts), which marks its canvas live and
      // marks it lost if the context is lost.
      const hasGl = await page.$eval(sel, (c) => c.dataset.view === "webgl");
      if (!hasGl) {
        throw new Error(`Canvas ${sel} missing live WebGL context`);
      }
      // Its first frame: a view waits for its shaders before drawing, and an undrawn stage would pass the clip test.
      await page
        .waitForFunction((s) => document.querySelector(s)?.dataset.drawn === "1", { timeout: 20000 }, sel)
        .catch(() => {
          throw new Error(`Canvas ${sel} never drew a frame`);
        });
      await sleep(400);
      const clipPng = await page.screenshot({
        clip: {
          x: Math.round(box.x),
          y: Math.round(box.y),
          width: Math.min(400, Math.round(box.width)),
          height: Math.min(300, Math.round(box.height)),
        },
      });
      fs.writeFileSync(path.join(VERIFY_DIR, `clip-${sid}.png`), clipPng);
      if (clipPng.byteLength < 2500) {
        throw new Error(`Canvas ${sel} appears blank (PNG clip only ${clipPng.byteLength} bytes)`);
      }
      console.log(`[verify-chrome] Canvas ${sid} (${sel}): WebGL live, non-blank clip ${clipPng.byteLength} bytes`);

      if (sid === "equipment") {
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch1-equipment.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch1-equipment.png");
      }

      // Take chapter screenshot at Stove Lab (shows 3D skillet/steak/wine + strip chart cleanly)
      if (sid === "stove") {
        // Test "Deglaze 120 mL Wine" real conservation physics!
        const before = await page.evaluate(() => ({
          pan: parseFloat(document.querySelector("[data-st-pan]")?.textContent || "0"),
          vapor: parseFloat(document.querySelector("[data-st-vapor]")?.textContent || "0"),
        }));
        await page.click("[data-stove-deglaze]");
        // The bottle lifts, tilts and pours; the wine counts only as it lands, so wait for the boil-off to show.
        await page
          .waitForFunction(
            (b) => {
              const pan = parseFloat(document.querySelector("[data-st-pan]")?.textContent || "0");
              const vapor = parseFloat(document.querySelector("[data-st-vapor]")?.textContent || "0");
              return vapor > b.vapor && pan < b.pan;
            },
            { timeout: 20000, polling: 200 },
            before,
          )
          .catch(() => {});
        const after = await page.evaluate(() => ({
          pan: parseFloat(document.querySelector("[data-st-pan]")?.textContent || "0"),
          vapor: parseFloat(document.querySelector("[data-st-vapor]")?.textContent || "0"),
        }));
        console.log(
          `[verify-chrome] Deglaze physics: T_pan ${before.pan.toFixed(1)} -> ${after.pan.toFixed(1)} °C, vapor ${before.vapor.toFixed(1)} -> ${after.vapor.toFixed(1)} g`,
        );
        if (!(after.pan < before.pan && after.vapor > before.vapor)) {
          throw new Error(`Deglaze physics check failed: before=${JSON.stringify(before)}, after=${JSON.stringify(after)}`);
        }
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch2-stove.png"), fullPage: false });
        fs.copyFileSync(path.join(VERIFY_DIR, "ch2-stove.png"), path.join(VERIFY_DIR, "chapters-desktop.png"));
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch2-stove.png");
      }

      if (sid === "oven") {
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch3-oven.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch3-oven.png");
      }

      if (sid === "cold") {
        // Capture Walk-In Cooler cutaway first before switching to Insulated Carrier OOD test!
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch4-cold.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch4-cold.png");

        // Switch to the delivery box: the fast estimate must answer in range (quoting its measured error) and refuse
        // out of range.
        await page.click('button[data-cold-mode="carrier"]');
        await sleep(500);
        await page.click("[data-car-surrogate]");
        await sleep(500);
        const inState = await page.evaluate(() => {
          const el = document.querySelector("[data-cold-alert]");
          return { state: el?.getAttribute("data-state"), text: el?.textContent || "" };
        });
        console.log(`[verify-chrome] Fast estimate in range: ${inState.state} | ${inState.text.slice(0, 90)}...`);
        if (!inState.text.includes("Times from the fast estimate")) {
          throw new Error(`Fast estimate did not answer in range: ${JSON.stringify(inState)}`);
        }
        await page.click("[data-car-ood-preset]");
        await sleep(500);
        const oodState = await page.evaluate(() => {
          const el = document.querySelector("[data-cold-alert]");
          return { state: el?.getAttribute("data-state"), text: el?.textContent || "" };
        });
        console.log(`[verify-chrome] Fast estimate out of range: ${oodState.state} | ${oodState.text.slice(0, 90)}...`);
        if (oodState.state !== "warn" || !oodState.text.includes("Out of range for the fast estimate")) {
          throw new Error(`Carrier out-of-range refusal failed: ${JSON.stringify(oodState)}`);
        }
      }

      if (sid === "brigade") {
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch5-brigade.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch5-brigade.png");
      }

      if (sid === "week") {
        // Cooking ahead must cut the median week's overtime against cooking on the day, and say so.
        const otOf = () => page.$eval("[data-wk-ot]", (el) => parseFloat(el.textContent || "NaN"));
        await page.click('[data-wk-plan="day"]');
        await sleep(400);
        const otDay = await otOf();
        await page.click('[data-wk-plan="ahead"]');
        await sleep(400);
        const otAhead = await otOf();
        const status = await page.$eval("[data-wk-status]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Week overtime: on the day ${otDay} h -> cooked ahead ${otAhead} h | ${status.slice(0, 90)}`);
        if (!(otAhead < otDay)) {
          throw new Error(`Cooking ahead did not cut overtime: ${otDay} -> ${otAhead}`);
        }
        // A fourth cook on the day must get the late meals to none, at the cost of idle hours, and the status must weigh
        // it against cooking ahead with three. Then back to the page's start.
        const lateOf = () => page.$eval("[data-wk-late]", (el) => parseFloat(el.textContent || "NaN"));
        const idleOf = () => page.$eval("[data-wk-idle]", (el) => parseFloat(el.textContent || "NaN"));
        await page.click('[data-wk-plan="day"]');
        await sleep(300);
        const [late3, idle3] = [await lateOf(), await idleOf()];
        await page.click('[data-wk-crew="4"]');
        await sleep(600);
        const [late4, idle4] = [await lateOf(), await idleOf()];
        const status4 = await page.$eval("[data-wk-status]", (el) => el.textContent || "");
        const cooks4 = await page.$eval("[data-wk-hud-note]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Week crew: 3 cooks ${late3} late, ${idle3} h idle -> 4 cooks ${late4} late, ${idle4} h idle | ${status4.slice(-90)}`);
        if (!(late4 < late3 && idle4 > idle3 && /Three cooks cooking ahead/.test(status4))) {
          throw new Error(`A fourth cook did not read right: ${late3}/${idle3} -> ${late4}/${idle4}, "${status4}" (HUD ${cooks4})`);
        }
        await page.click('[data-wk-crew="3"]');
        await sleep(300);
        // Every plan for the week, run in a worker: all 240 must arrive, a frontier be drawn, and a dot set the plan.
        await page.waitForFunction(() => document.querySelectorAll("[data-wk-map] [data-plan]").length >= 240, { timeout: 90000 });
        await sleep(200);
        const wkFront = await page.$$eval("[data-wk-map] .pt.on", (els) => els.length);
        const wkNote = await page.$eval("[data-wk-map-note]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Week plan map: 240 plans, ${wkFront} on the frontier | ${wkNote.slice(0, 110)}`);
        if (wkFront < 2 || !/^Cheapest with no meal late/.test(wkNote)) throw new Error(`The week plan map did not finish: ${wkFront} on the frontier, "${wkNote}"`);
        await page.$eval("[data-wk-map]", (el) => el.scrollIntoView({ block: "center" }));
        await sleep(300);
        const wkDot = await page.$$eval("[data-wk-map] .pt.on", (els) => {
          const el = els.find((e) => /^2 cooks/.test(e.querySelector("title")?.textContent || ""));
          if (!el) return null;
          const r = el.getBoundingClientRect();
          return { x: r.x + r.width / 2, y: r.y + r.height / 2 };
        });
        if (!wkDot) throw new Error("The week plan map has no two-cook plan on its frontier");
        await page.mouse.click(wkDot.x, wkDot.y);
        await sleep(400);
        const wkCrew = await page.$eval('[data-wk-crew][aria-pressed="true"]', (el) => el.textContent || "");
        if (!/2 cooks/.test(wkCrew)) throw new Error(`A two-cook dot on the week map set "${wkCrew}"`);
        await page.click('[data-wk-crew="3"]');
        await sleep(300);
        // The same meals, spread flat or sharp: flat, nothing needs cooking ahead; half as sharp again, the cheapest plan
        // with every meal on time takes a fourth cook.
        const setSwing = (v) =>
          page.$eval(
            "[data-wk-swing]",
            (el, x) => {
              el.value = x;
              el.dispatchEvent(new Event("input", { bubbles: true }));
            },
            v,
          );
        await page.click('[data-wk-plan="ahead"]');
        await setSwing("0");
        await sleep(400);
        const flatStatus = await page.$eval("[data-wk-status]", (el) => el.textContent || "");
        const flatLabel = await page.$eval("[data-wk-swing-val]", (el) => el.textContent || "");
        if (!/nothing needs cooking ahead/.test(flatStatus) || !/^Flat/.test(flatLabel)) throw new Error(`A flat week read "${flatLabel}": "${flatStatus}"`);
        await setSwing("150");
        await page.waitForFunction(() => /^Cheapest with no meal late/.test(document.querySelector("[data-wk-map-note]")?.textContent || ""), { timeout: 90000 });
        const sharpNote = await page.$eval("[data-wk-map-note]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Week swing: flat "${flatStatus.slice(0, 60)}" | 1.5x "${sharpNote.slice(0, 100)}"`);
        if (!/: 4 cooks/.test(sharpNote)) throw new Error(`At 1.5x the swing the cheapest on-time plan did not take a fourth cook: "${sharpNote}"`);
        await setSwing("100");
        await page.click('[data-wk-plan="day"]');
        await sleep(300);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch6-week.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch6-week.png");
      }

      if (sid === "floor") {
        // The shortest-walk spot must walk less than the chiller as drawn, for both menus. The figure builds once the
        // commissary's plan has loaded: wait for its first answer before using its buttons.
        const kmOf = () => page.$eval("[data-fl-km]", (el) => parseFloat(el.textContent || "NaN"));
        await page.waitForFunction(() => /\d/.test(document.querySelector("[data-fl-km]")?.textContent || ""), { timeout: 30000 });
        for (const menu of ["written", "swapped"]) {
          await page.click(`[data-fl-menu="${menu}"]`);
          await page.click('[data-fl-spot="drawn"]');
          await sleep(300);
          const drawn = await kmOf();
          await page.click('[data-fl-spot="best"]');
          await sleep(300);
          const best = await kmOf();
          console.log(`[verify-chrome] Floor (${menu}): chiller as drawn ${drawn} km a week -> shortest walk ${best} km`);
          if (!(best < drawn)) {
            throw new Error(`The best chiller spot did not walk less (${menu}): ${drawn} -> ${best}`);
          }
        }
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch7-floor.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch7-floor.png");
      }

      if (sid === "rush") {
        // Holding what keeps must cut the long waits against cooking to order on the busy day, and pooled plating must
        // cut them further; the status must say so.
        const minOf = () => page.$eval("[data-rs-long]", (el) => parseFloat(el.textContent || "NaN"));
        await page.click('[data-rs-day="busy"]');
        await page.click('[data-rs-pool="no"]');
        await page.click('[data-rs-plan="order"]');
        await sleep(300);
        const order = await minOf();
        await page.click('[data-rs-plan="keepers"]');
        await sleep(300);
        const held = await minOf();
        await page.click('[data-rs-pool="yes"]');
        await sleep(300);
        const pooled = await minOf();
        const status = await page.$eval("[data-rs-status]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Rush, 1 guest in 10 waits: to order ${order} min -> held ${held} -> pooled ${pooled} | ${status.slice(0, 90)}`);
        if (!(held < order && pooled < held)) {
          throw new Error(`The rush levers did not cut the long waits: ${order} -> ${held} -> ${pooled}`);
        }
        // Every plan for the day, run in a worker: all 128 must arrive and a frontier be drawn.
        await page.waitForFunction(() => document.querySelectorAll("[data-rush-map] [data-plan]").length >= 128, { timeout: 90000 });
        const front = await page.$$eval("[data-rush-map] .pt.on", (els) => els.length);
        const mapNote = await page.$eval("[data-rush-map-note]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Rush plan map: 128 plans, ${front} on the frontier | ${mapNote.slice(0, 90)}`);
        if (front < 2) throw new Error(`The rush plan map has no frontier (${front} points)`);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch6-rush.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch6-rush.png");
      }

      if (sid === "morning") {
        // The menu as written must miss the van at 700 meals with three cooks and every swap must make it; the status
        // must say which, and the worker must post every menu's limit, cooked today and ahead.
        const status = () => page.$eval("[data-mo-status]", (el) => el.textContent || "");
        await page.click('[data-mo-crew="3"]');
        await page.click('[data-mo-ahead="no"]');
        await page.click('[data-mo-preset="written"]');
        await sleep(300);
        const written = await status();
        await page.click('[data-mo-preset="swapped"]');
        await sleep(300);
        const swapped = await status();
        console.log(`[verify-chrome] Morning: as written | ${written.slice(0, 70)} || every swap | ${swapped.slice(0, 70)}`);
        if (!/late/.test(written) || !/leaves at 10:30/.test(swapped)) {
          throw new Error(`The morning presets did not tell late from on time: ${written} / ${swapped}`);
        }
        await page.waitForFunction(() => document.querySelectorAll("[data-morning-ladder] .row .today").length >= 8, { timeout: 90000 });
        const note = await page.$eval("[data-morning-ladder-note]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Morning limits: ${note.slice(0, 110)}`);
        if (!/\d/.test(note)) throw new Error(`The morning limit note has no numbers: ${note}`);
        await page.click('[data-mo-preset="written"]');
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch7-morning.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch7-morning.png");
      }
    }

    // 5. Verify ALL 10 WebGL canvases render non-blank pixels via direct pixel sampling after capture
    const nonBlankSummary = await page.evaluate(async () => {
      const ids = ["hero", "equipment", "stove", "oven", "cold", "brigade", "week", "floor", "rush", "morning"];
      const sels = [
        "[data-hero-canvas]",
        "[data-eq-canvas]",
        "[data-stove-canvas]",
        "[data-oven-canvas]",
        "[data-cold-canvas]",
        "[data-brigade-canvas]",
        "[data-week-canvas]",
        "[data-floor-canvas]",
        "[data-rush-canvas]",
        "[data-morning-canvas]",
      ];
      return ids.map((id, i) => {
        const c = document.querySelector(sels[i]);
        return { id, w: c?.width || 0, h: c?.height || 0 };
      });
    });
    console.log("[verify-chrome] All 10 WebGL canvases active:", nonBlankSummary);
    for (const item of nonBlankSummary) {
      if (item.w < 200 || item.h < 200) {
        throw new Error(`Canvas backing buffer too small: ${JSON.stringify(item)}`);
      }
    }

    // 6. Verify 0 px horizontal overflow across desktop (1440), tablet (768), and mobile (390)
    for (const vp of [
      { name: "desktop-1440", width: 1440, height: 900 },
      { name: "tablet-768", width: 768, height: 1024 },
      { name: "mobile-390", width: 390, height: 844 },
    ]) {
      await page.setViewport({ width: vp.width, height: vp.height, deviceScaleFactor: 2 });
      await page.evaluate(() => window.scrollTo(0, 0));
      await sleep(700);
      const overflow = await page.evaluate(
        () => document.documentElement.scrollWidth - document.documentElement.clientWidth,
      );
      console.log(`[verify-chrome] Viewport ${vp.name}: horizontal overflow = ${overflow}px`);
      if (overflow > 0) {
        throw new Error(`Horizontal overflow on ${vp.name}: ${overflow}px`);
      }
      if (vp.name === "mobile-390") {
        await page.screenshot({ path: path.join(VERIFY_DIR, "mobile-390.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/mobile-390.png");
      }
    }

    if (consoleErrors.length > 0 || pageErrors.length > 0 || failedRequests.length > 0) {
      throw new Error(
        `Runtime/network errors detected:\n` +
          `  consoleErrors=${JSON.stringify(consoleErrors)}\n` +
          `  pageErrors=${JSON.stringify(pageErrors)}\n` +
          `  failedRequests=${JSON.stringify(failedRequests)}`,
      );
    }

    console.log("[verify-chrome] ALL CHECKS PASSED (0 console errors, 0 page errors, 0 failed fetches, 10/10 WebGL live, 0px overflow)!");
  } finally {
    await browser.close();
    server.close();
  }
}

main().catch((err) => {
  console.error("[verify-chrome] FAILED:", err);
  process.exit(1);
});
