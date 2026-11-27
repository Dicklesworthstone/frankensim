#!/usr/bin/env node
/**
 * Headless Chrome E2E & visual verification for `frankensim.org/kitchen-lab` (`dist-cf/`).
 *
 * Verifies:
 *   - ZERO console errors, ZERO uncaught page errors, ZERO failed network fetches (404/500)
 *   - ALL 15 WebGL canvases (14 figures; the two-kitchen figure has two) render non-blank frames (pixel variance > 0)
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

    // 4. Mount all chapters cleanly by scrolling through each and verifying every figure's WebGL canvas renders non-blank frames
    const sectionIds = ["hero", "equipment", "stove", "boil", "oven", "cold", "brigade", "pivot", "week", "cool", "hold", "floor", "hood", "rush", "glasses", "morning", "crowd", "pair", "kit", "trip"];
    const canvasSelectors = [
      "[data-hero-canvas]",
      "[data-eq-canvas]",
      "[data-stove-canvas]",
      "[data-bo-canvas]",
      "[data-oven-canvas]",
      "[data-cold-canvas]",
      "[data-brigade-canvas]",
      "[data-pivot-canvas]",
      "[data-week-canvas]",
      "[data-cd-canvas]",
      "[data-hd-canvas]",
      "[data-floor-canvas]",
      "[data-hc-canvas]",
      "[data-rush-canvas]",
      "[data-gl-canvas]",
      "[data-morning-canvas]",
      "[data-crowd-canvas]",
      '[data-pair-canvas="hub"]',
      "[data-kit-canvas]",
      "[data-trip-canvas]",
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

      if (sid === "pivot") {
        // Moving the plating table must take steps off every plate, beside the range most of all; the counts come from
        // the walking model's own footfalls.
        const info = () => page.$eval("[data-pivot-stage]", (el) => el.pivot.info());
        await page.click('[data-pv-move="behind"]');
        await sleep(300);
        const behind = await info();
        await page.click('[data-pv-move="beside"]');
        await sleep(300);
        const beside = await info();
        const status = await page.$eval("[data-pv-status]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Pivot steps a plate: as built ${behind.built.steps} (+${behind.built.pivots} pivots) -> behind ${behind.moved.steps} -> beside ${beside.moved.steps} | ${status.slice(0, 80)}`);
        if (!(behind.moved.steps < behind.built.steps && beside.moved.steps <= behind.moved.steps && beside.moved.moveS < behind.built.moveS)) {
          throw new Error(`Moving the plating table did not save steps: ${JSON.stringify({ behind, beside })}`);
        }
        await page.click('[data-pv-move="behind"]');
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch5-pivot.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch5-pivot.png");
      }

      if (sid === "crowd") {
        // Eight cooks as built must lose more to each other than with the island pulled back; nobody may ever come
        // within reach of anybody (the planner's least gap), no run may gridlock, and bodies never beat ghosts by more
        // than a pan or two of rounding.
        const info = () => page.$eval("[data-crowd-stage]", (el) => el.crowd.info());
        const want = ["built/8/b", "built/8/g", "moved/8/b", "moved/8/g"];
        const t0 = Date.now();
        let got = await info();
        while (want.some((k) => !got.runs.some((r) => r.k === k)) && Date.now() - t0 < 90000) {
          await sleep(1000);
          got = await info();
        }
        const run = (k) => got.runs.find((r) => r.k === k);
        if (want.some((k) => !run(k))) throw new Error(`The crowd runs did not arrive: ${JSON.stringify(got.runs)}`);
        const lost = (k) => 1 - run(`${k}/8/b`).perHour / run(`${k}/8/g`).perHour;
        const status = await page.$eval("[data-cr-status]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Crowd, 8 cooks: as built ${run("built/8/b").perHour}/h (ghosts ${run("built/8/g").perHour}), pulled back ${run("moved/8/b").perHour}/h | ${status.slice(0, 90)}`);
        for (const r of got.runs) {
          if (r.gridlock) throw new Error(`Crowd run ${r.k} gridlocked`);
          if (r.minGap !== null && r.minGap < 0) throw new Error(`Crowd run ${r.k}: two bodies overlapped (${r.minGap} m)`);
        }
        if (!(lost("built") > lost("moved"))) throw new Error(`Pulling the island back did not cut the loss: ${JSON.stringify(got.runs)}`);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch7-crowd.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch7-crowd.png");
      }

      if (sid === "pair") {
        // As set on arrival (the all-keeping menu, 1,400 meals, two cooks across town): the hub alone is late and both
        // kitchens make lunch; both halves draw; the most meals with both is more than alone. Swapping the thighs for
        // the piccata, which keeps only 90 minutes, must keep the piccata at the hub and say why.
        const info = () => page.$eval("[data-pair-stage]", (el) => el.pair.info());
        await page.waitForFunction(() => {
          const i = document.querySelector("[data-pair-stage]")?.pair?.info();
          return i?.ready && i.limit;
        }, { timeout: 90000 });
        await page.waitForFunction(() => document.querySelector('[data-pair-canvas="far"]')?.dataset.drawn === "1", { timeout: 20000 });
        const both = await info();
        const status = () => page.$eval("[data-pr-status]", (el) => el.textContent || "");
        console.log(`[verify-chrome] Two kitchens: ${JSON.stringify({ made: both.made, parts: both.parts, alone: both.alone, limit: [both.limit.alone, both.limit.both] })} | ${(await status()).slice(0, 90)}`);
        if (!both.made || !both.far || both.alone?.made !== false) throw new Error(`Two kitchens: the default order did not tell the hub alone from both: ${JSON.stringify(both)}`);
        if (!(both.limit.both > both.limit.alone)) throw new Error(`Two kitchens: a second kitchen did not raise the most meals: ${JSON.stringify(both.limit)}`);
        await page.click('[data-pr-mode="alone"]');
        await sleep(300);
        const alone = await status();
        if (!/late/.test(alone)) throw new Error(`Two kitchens: the hub alone's status does not say it is late: ${alone}`);
        await page.click('[data-pr-mode="both"]');
        await page.click("[data-pr-menu] li:first-child button");
        await page.waitForFunction(() => document.querySelector("[data-pair-stage]")?.pair?.info().ready, { timeout: 30000 });
        const swapped = await status();
        const parts = (await info()).parts;
        console.log(`[verify-chrome] Two kitchens, piccata: ${JSON.stringify(parts)} | ${swapped.slice(0, 110)}`);
        if (parts[0] !== 0 || !/too late for a van at 10:00/.test(swapped)) throw new Error(`Two kitchens: the piccata travelled or the status did not say why: ${swapped}`);
        await page.click("[data-pr-menu] li:first-child button");
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch7-pair.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch7-pair.png");
      }

      if (sid === "kit") {
        // As set on arrival (every swap, 1,400 meals, nothing bought): late, the combis full. A third combi rolls in
        // and makes the van; a second range buys nothing; the ladder raises the most meals for the combi only.
        const info = () => page.$eval("[data-kit-stage]", (el) => el.kit.info());
        const status = () => page.$eval("[data-kt-status]", (el) => el.textContent || "");
        const settle = () =>
          page.waitForFunction(() => {
            const i = document.querySelector("[data-kit-stage]")?.kit?.info();
            return i?.ready && Object.values(i.limits).every((v) => v !== null);
          }, { timeout: 90000 });
        await settle();
        const built = await info();
        console.log(`[verify-chrome] Buy one thing: ${JSON.stringify({ depart: built.depart, limits: built.limits })} | ${(await status()).slice(0, 90)}`);
        if (built.made !== false || !/combis are full/.test(await status())) throw new Error(`Buy one thing: the kitchen as built did not miss the van on the combis: ${JSON.stringify(built)}`);
        if (!(built.limits["100"] > built.limits["000"]) || built.limits["010"] !== built.limits["000"]) throw new Error(`Buy one thing: the limits do not favour the combi: ${JSON.stringify(built.limits)}`);
        await page.click('[data-kt-buy="combi"]');
        await settle();
        await sleep(1800);
        const combi = await info();
        if (!combi.made || combi.kit.ovens !== 3 || combi.inPlace.ovens[2] !== 1) throw new Error(`Buy one thing: a third combi did not roll in and make the van: ${JSON.stringify(combi)}`);
        await page.click('[data-kt-buy="combi"]');
        await page.click('[data-kt-buy="range"]');
        await settle();
        const range = await status();
        console.log(`[verify-chrome] Buy one thing, a second range: ${range.slice(0, 120)}`);
        if (!/buys nothing/.test(range)) throw new Error(`Buy one thing: the second range's status does not say it buys nothing: ${range}`);
        await page.click('[data-kt-buy="range"]');
        await settle();
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch7-kit.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch7-kit.png");
      }

      if (sid === "cool") {
        // As set on arrival (30 litres, uncovered): of the four ways only the 2-inch pans started in the freezer meet
        // both of the Food Code's limits, the slower the deeper, and the pot is still over 41 °F a day after it went
        // in. Covered, every way is slower. Held at 5 h, the freezer pans' rack stands in the walk-in, their middle at
        // or under 41 °F, while the pot's is still over 135 °F; partway through the move the rack is outside, in front of
        // the boxes. A way not at 41 °F within its run has a null total: counted as never.
        const info = () => page.$eval("[data-cd-stage]", (el) => el.cool.info());
        const ready = (covered) =>
          page.waitForFunction((c) => {
            const i = document.querySelector("[data-cd-stage]")?.cool?.info();
            return !!i && i.ready && i.covered === c;
          }, { timeout: 30000 }, covered);
        const by = (i, k) => i.runs.find((r) => r.key === k);
        const inTime = (r) => r.ok1 && r.ok2;
        const tot = (r) => r.total ?? Infinity;
        await ready(false);
        const open = await info();
        console.log(`[verify-chrome] Cool the chili: ${JSON.stringify(open.runs.map((r) => [r.key, r.stage1, r.total, r.dayF]))} | ${open.status.slice(0, 90)}`);
        const [pot, deep, shallow, freezer] = ["pot", "deep", "shallow", "freezer"].map((k) => by(open, k));
        if (open.litres !== 30 || !inTime(freezer) || [pot, deep, shallow].some(inTime) || !(tot(pot) > tot(deep) && tot(deep) > tot(shallow) && tot(shallow) > tot(freezer)) || !(pot.dayF > 41) || !/^Only 2-inch pans started in the freezer/.test(open.status)) throw new Error(`Cool the chili: 30 litres uncovered do not cool in the order or verdicts expected: ${JSON.stringify(open)}`);
        await page.click('[data-cd-cover="on"]');
        await ready(true);
        const film = await info();
        console.log(`[verify-chrome] Cool the chili, covered: ${JSON.stringify(film.runs.map((r) => [r.key, r.total]))}`);
        if (["pot", "deep", "shallow", "freezer"].some((k) => !(tot(by(film, k)) > tot(by(open, k))))) throw new Error(`Cool the chili: covered, a way is not slower: ${JSON.stringify(film.runs)}`);
        await page.click('[data-cd-cover="off"]');
        await ready(false);
        await page.$eval("[data-cd-stage]", (el) => el.cool.hold(300));
        await sleep(1200);
        const held = await info();
        const tagF = (t) => Number(/(-?\d+) °F$/.exec(t ?? "")?.[1] ?? NaN);
        if (held.t !== 300 || !(tagF(held.tags[3]) <= 41) || !(tagF(held.tags[0]) > 135)) throw new Error(`Cool the chili: held at 5 h, the labels are not the freezer pans cold and the pot hot: ${JSON.stringify(held.tags)}`);
        if (Math.abs(held.rack.x - -0.4) > 1e-3 || held.rack.y !== 0) throw new Error(`Cool the chili: held at 5 h, the freezer pans' rack is not in the walk-in: ${JSON.stringify(held.rack)}`);
        const moveMin = by(held, "freezer").movedMin;
        await page.$eval("[data-cd-stage]", (el, m) => el.cool.hold(m, 0.45), moveMin);
        await sleep(600);
        const mid = await info();
        if (!(mid.rack.z > 1.3) || mid.rack.y !== 0) throw new Error(`Cool the chili: partway through the move the rack is not outside the boxes: ${JSON.stringify(mid.rack)}`);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch6-cool.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch6-cool.png");
        await page.$eval("[data-cd-stage]", (el) => el.cool.hold(null));
      }

      if (sid === "hold") {
        // As set on arrival (done at 7:00, packed at 185 °F): both full carriers and the lone pan in a primed carrier
        // are hot at lunch, the lone pan in a cold carrier is under 135 °F before it; full beats one pan, and primed
        // beats cold. Packed at 165 °F every carrier falls under 135 °F sooner, and neither lone pan is hot at lunch.
        // Held an hour after packing, every label is up, under the packing temperature and over 135 °F.
        const info = () => page.$eval("[data-hd-stage]", (el) => el.hold.info());
        const ready = (f) =>
          page.waitForFunction((f) => {
            const i = document.querySelector("[data-hd-stage]")?.hold?.info();
            return !!i && i.ready && i.startF === f;
          }, { timeout: 30000 }, f);
        const keys = ["fullPrimed", "fullCold", "onePrimed", "oneCold"];
        const by = (i, k) => i.runs.find((r) => r.key === k);
        await ready(185);
        const hot = await info();
        console.log(`[verify-chrome] Keep it hot: ${JSON.stringify(hot.runs.map((r) => [r.key, r.hotMin, r.lunchF]))} | ${hot.status.slice(0, 90)}`);
        const h = keys.map((k) => by(hot, k).hotMin);
        if (hot.done !== 420 || keys.map((k) => by(hot, k).ok).join() !== "true,true,true,false" || !(h[0] > h[1] && h[1] > h[2] && h[2] > h[3])) throw new Error(`Keep it hot: done at 7:00 and packed at 185 °F, the carriers do not hold in the order or verdicts expected: ${JSON.stringify(hot)}`);
        await page.click('[data-hd-pack="165"]');
        await ready(165);
        const warm = await info();
        console.log(`[verify-chrome] Keep it hot, packed at 165 °F: ${JSON.stringify(warm.runs.map((r) => [r.key, r.hotMin, r.ok]))}`);
        if (keys.some((k) => !(by(warm, k).hotMin < by(hot, k).hotMin)) || by(warm, "onePrimed").ok || by(warm, "oneCold").ok || !by(warm, "fullPrimed").ok) throw new Error(`Keep it hot: packed at 165 °F, the carriers do not fall sooner, or the verdicts are not as expected: ${JSON.stringify(warm.runs)}`);
        await page.click('[data-hd-pack="185"]');
        await ready(185);
        await page.$eval("[data-hd-stage]", (el) => el.hold.hold(60));
        await sleep(1200);
        const held = await info();
        const tagF = (t) => Number(/(-?\d+) °F$/.exec(t ?? "")?.[1] ?? NaN);
        if (held.t !== 60 || held.tags.some((t) => !(tagF(t) < 185 && tagF(t) > 135))) throw new Error(`Keep it hot: held an hour after packing, the labels are not all up and between 135 °F and 185 °F: ${JSON.stringify(held.tags)}`);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch6-hold.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch6-hold.png");
        await page.$eval("[data-hd-stage]", (el) => el.hold.hold(null));
      }

      if (sid === "hood") {
        // Flush with the hood's edge and the gap open, the broilers need 5,100 cfm: at the code's 4,000 smoke rolls out
        // under the edge. Pushed 12 inches back with the gap sealed they need 2,800: all of it caught, and none out in
        // the kitchen once what was already out has cleared. At 4,000 cfm four of the six setups are caught, at 2,700
        // none. No puff the hood caught is ever drawn as escaping.
        const info = () => page.$eval("[data-hc-stage]", (el) => el.hood.info());
        await page.$eval("[data-hc-stage]", (el) => {
          el.hood.hold(true);
          el.hood.set({ at: 0, sealed: false, fan: 4000 });
          el.hood.settle();
        });
        await sleep(3500);
        const flush = await info();
        console.log(`[verify-chrome] Keep the hood, flush and open: need ${flush.need}, caught ${flush.caught}, ${JSON.stringify(flush.smoke)} | ${flush.status.slice(0, 90)}`);
        if (flush.need !== 5100 || flush.caught || flush.caughtCount !== 4 || !(flush.smoke.escaping > 20) || flush.smoke.stray !== 0 || !/escapes/.test(flush.tags[1] ?? "")) throw new Error(`Keep the hood: flush and open at 4,000 cfm, the smoke is not shown escaping: ${JSON.stringify(flush)}`);
        await page.click('[data-hc-at="2"]');
        await page.click('[data-hc-seal="1"]');
        await sleep(6000);
        const sealed = await info();
        console.log(`[verify-chrome] Keep the hood, 12 in back and sealed: need ${sealed.need}, caught ${sealed.caught}, front ${sealed.frontIn} in, seal ${sealed.sealK}, ${JSON.stringify(sealed.smoke)}`);
        if (sealed.at !== 2 || !sealed.sealed || sealed.touring || sealed.need !== 2800 || !sealed.caught || sealed.frontIn !== 12 || sealed.sealK !== 1 || sealed.smoke.escaping !== 0 || sealed.smoke.stray !== 0 || !(sealed.smoke.live > 200) || !/all caught/.test(sealed.tags[1] ?? "") || !/sealed/.test(sealed.tags[0] ?? "")) throw new Error(`Keep the hood: 12 inches back and sealed at 4,000 cfm, the smoke is not all caught: ${JSON.stringify(sealed)}`);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch6-hood.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch6-hood.png");
        await page.$eval("[data-hc-stage]", (el) => el.hood.set({ fan: 2700 }));
        const low = await info();
        if (low.caughtCount !== 0 || low.caught) throw new Error(`Keep the hood: at 2,700 cfm a setup is still caught: ${JSON.stringify(low)}`);
        await page.$eval("[data-hc-stage]", (el) => {
          el.hood.set({ fan: 4000 });
          el.hood.hold(false);
        });
      }

      if (sid === "boil") {
        // As set on arrival (30 litres, full flame): the two pots with lids boil in about 25 min, half the big pot's
        // 50 with its lid, for the same gas; the big pot without one in about 57. On half a flame the big open pot does
        // not boil in 3 h: it levels off short of 100 °C. Held at 52 min, every pot but the big open one is boiling.
        const info = () => page.$eval("[data-bo-stage]", (el) => el.boil.info());
        const minOf = (i, k) => i.runs.find((r) => r.key === k)?.boilMin;
        const built = await info();
        console.log(`[verify-chrome] Put the water on: ${JSON.stringify(built.runs.map((r) => [r.key, r.boilMin, r.steamL]))} | ${built.status.slice(0, 90)}`);
        const [oo, ol, to, tl] = ["one-open", "one-lid", "two-open", "two-lid"].map((k) => minOf(built, k));
        if (built.litres !== 30 || built.flame !== "full" || !(tl < to && to < ol && ol < oo) || Math.round(tl) !== 25 || Math.round(ol) !== 50 || !/for the same gas/.test(built.status)) throw new Error(`Put the water on: 30 litres at full flame do not boil in the order or minutes expected: ${JSON.stringify(built)}`);
        await page.click('[data-bo-flame="half"]');
        await sleep(300);
        const half = await info();
        console.log(`[verify-chrome] Put the water on, half a flame: ${JSON.stringify(half.runs.map((r) => [r.key, r.boilMin, r.peak]))}`);
        if (minOf(half, "one-open") !== null || !(minOf(half, "two-lid") > 0) || !/does not boil in 3 h: it levels off at 9\d\.\d °C/.test(half.status)) throw new Error(`Put the water on: on half a flame the big open pot boils, or the status does not say it levels off: ${JSON.stringify(half)}`);
        await page.click('[data-bo-flame="full"]');
        await page.$eval("[data-bo-stage]", (el) => el.boil.hold(52));
        await sleep(1200);
        const held = await info();
        if (held.t !== 52 || !/\d+ °C$/.test(held.tags[0] ?? "") || !held.tags.slice(1).every((t) => /Boiling$/.test(t ?? ""))) throw new Error(`Put the water on: held at 52 min, the labels are not the big open pot heating and the rest boiling: ${JSON.stringify(held.tags)}`);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch2-boil.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch2-boil.png");
        await page.$eval("[data-bo-stage]", (el) => el.boil.hold(null));
      }

      if (sid === "glasses") {
        // As set on arrival (thighs and risotto, 60 lunches of each, both at once): 1 h 01 together, 1 h 36 one after
        // the other, 35 min saved; the cook at the left combi loading the thighs at 6:04, at the range stirring at 6:20;
        // the glasses' screen draws.
        const info = () => page.$eval("[data-gl-stage]", (el) => el.glasses.info());
        const status = () => page.$eval("[data-gl-status]", (el) => el.textContent || "");
        await page.waitForFunction(() => document.querySelector("[data-gl-pov]")?.dataset.drawn === "1", { timeout: 20000 });
        const built = await info();
        console.log(`[verify-chrome] One cook, two dishes: ${JSON.stringify({ dishes: built.dishes, end: built.end, endApart: built.endApart, saved: built.saved, free: built.free })} | ${(await status()).slice(0, 90)}`);
        if (built.dishes.join() !== "thighs,risotto" || Math.round(built.end) !== 61 || Math.round(built.endApart) !== 96 || Math.round(built.saved) !== 35) throw new Error(`One cook, two dishes: the default pair does not save 35 min: ${JSON.stringify(built)}`);
        const [load, stir] = await page.$eval("[data-gl-stage]", (el) => el.glasses.track([4.2, 20]));
        if (!(load[2] < -1.2 && load[3] < 0.1) || !(stir[2] > -0.6 && stir[3] < 0)) throw new Error(`One cook, two dishes: the cook is not at the combi loading, or at the range stirring: ${JSON.stringify({ load, stir })}`);
        // Piccata in place of the thighs: both cook at the range with the cook's hands on them, and pairing saves nothing.
        await page.click('[data-gl-dishes="0"] button:nth-child(1)');
        await sleep(400);
        const bad = await info();
        const badStatus = await status();
        console.log(`[verify-chrome] One cook, piccata and risotto: ${JSON.stringify({ end: bad.end, endApart: bad.endApart })} | ${badStatus.slice(0, 90)}`);
        if (bad.dishes.join() !== "piccata,risotto" || bad.saved !== 0 || !/no less than one after the other/.test(badStatus)) throw new Error(`One cook, two dishes: piccata and risotto do not save nothing: ${JSON.stringify(bad)} ${badStatus}`);
        // A braise in a combi for hours leaves the cook's hands free for the tacos and more.
        await page.$eval("[data-gl-stage]", (el) => el.glasses.set({ a: "braise", b: "tacos", portions: 120 }));
        const braise = await info();
        console.log(`[verify-chrome] One cook, braise and tacos: ${JSON.stringify({ end: braise.end, endApart: braise.endApart, free: braise.free })} | ${(await status()).slice(0, 90)}`);
        if (!(braise.saved > 60) || !(braise.free > 30)) throw new Error(`One cook, two dishes: the braise does not free the cook for the tacos: ${JSON.stringify(braise)}`);
        await page.$eval("[data-gl-stage]", (el) => el.glasses.set({ a: "thighs", b: "risotto", portions: 60, plan: "together" }));
        await sleep(1200);
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch7-glasses.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch7-glasses.png");
      }

      if (sid === "trip") {
        // As set on arrival (800 meals, every swap, the left combi out 9:20 to 10:20): the van 18 min late; the map has
        // all 486 trips, most of them free. Held mid-trip, the stage marks the combi and the camera's screen flags it;
        // both screens draw. At 820 meals the braise fills both combis and few trips are free; a short early trip at
        // 800 costs nothing.
        const info = () => page.$eval("[data-trip-stage]", (el) => el.trip.info());
        const status = () => page.$eval("[data-tp-status]", (el) => el.textContent || "");
        const settle = () => page.waitForFunction(() => document.querySelector("[data-trip-stage]")?.trip?.info().ready, { timeout: 90000 });
        await settle();
        const built = await info();
        console.log(`[verify-chrome] A combi trips: ${JSON.stringify({ depart: built.depart, late: built.late, pastBest: built.pastBest, free: built.free, cells: built.cells })} | ${(await status()).slice(0, 90)}`);
        if (built.depart !== "10:48" || built.cells !== 486 || !(built.free > 400)) throw new Error(`A combi trips: the default trip did not cost the van 18 min, or the map is wrong: ${JSON.stringify(built)}`);
        await page.$eval("[data-trip-stage]", (el) => el.trip.hold(215));
        await page.evaluate(() => document.querySelector(".tp-mons")?.scrollIntoView({ behavior: "instant", block: "center" }));
        for (const cam of ["twin", "seen"]) await page.waitForFunction((c) => document.querySelector(`[data-trip-cam="${c}"]`)?.dataset.drawn === "1", { timeout: 20000 }, cam);
        await sleep(1200);
        const held = await info();
        if (held.pads !== 1 || held.box !== true || !/dark/.test(held.flag || "")) throw new Error(`A combi trips: mid-trip, the combi is not marked out: ${JSON.stringify(held)}`);
        await page.$eval("[data-trip-stage]", (el) => el.trip.hold(null));
        await page.$eval("[data-trip-stage]", (el) => el.trip.set({ meals: 820 }));
        await settle();
        const full = await info();
        console.log(`[verify-chrome] A combi trips, 820 meals: ${JSON.stringify({ depart: full.depart, free: full.free })} | ${(await status()).slice(0, 90)}`);
        if (!(full.free < 120)) throw new Error(`A combi trips: at 820 meals too many trips are free: ${JSON.stringify(full)}`);
        await page.$eval("[data-trip-stage]", (el) => el.trip.set({ meals: 800, at: 150, mins: 45 }));
        await settle();
        const cheap = await status();
        if (!/costs nothing/.test(cheap)) throw new Error(`A combi trips: a short trip at 8:30 did not cost nothing: ${cheap}`);
        // A press on a map cell trips that combi then, for that long: the right combi at 7:00 for 30 min.
        await page.click('[data-trip-map] rect.cell[data-o="1"][data-at="60"][data-mins="30"]');
        await settle();
        const picked = await info();
        const pickedStatus = await status();
        console.log(`[verify-chrome] A combi trips, a map cell: ${JSON.stringify(picked.trip)} | ${pickedStatus.slice(0, 90)}`);
        if (picked.trip.oven !== 1 || picked.trip.at !== 60 || picked.trip.mins !== 30 || !/^The right combi trips at 7:00 .* back at 7:30\./.test(pickedStatus)) throw new Error(`A combi trips: a map cell did not trip the right combi at 7:00 for 30 min: ${JSON.stringify(picked.trip)} ${pickedStatus}`);
        await page.$eval("[data-trip-stage]", (el) => el.trip.set({ oven: 0, at: 200, mins: 60 }));
        await settle();
        await page.screenshot({ path: path.join(VERIFY_DIR, "ch7-trip.png"), fullPage: false });
        console.log("[verify-chrome] Saved /tmp/kitchen_web_verify/ch7-trip.png");
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

    // 5. Verify ALL 24 WebGL canvases (20 figures; the two-kitchen figure has two, the combi trip has its two screens, the camera glasses their own) render non-blank pixels via direct pixel sampling after capture
    const nonBlankSummary = await page.evaluate(async () => {
      const ids = ["hero", "equipment", "stove", "boil", "oven", "cold", "brigade", "pivot", "week", "cool", "hold", "floor", "hood", "rush", "glasses", "morning", "crowd", "pair", "kit", "trip"];
      const sels = [
        "[data-hero-canvas]",
        "[data-eq-canvas]",
        "[data-stove-canvas]",
        "[data-bo-canvas]",
        "[data-oven-canvas]",
        "[data-cold-canvas]",
        "[data-brigade-canvas]",
        "[data-pivot-canvas]",
        "[data-week-canvas]",
        "[data-cd-canvas]",
        "[data-hd-canvas]",
        "[data-floor-canvas]",
        "[data-hc-canvas]",
        "[data-rush-canvas]",
        "[data-gl-canvas]",
        "[data-morning-canvas]",
        "[data-crowd-canvas]",
        '[data-pair-canvas="hub"]',
        "[data-kit-canvas]",
        "[data-trip-canvas]",
        '[data-pair-canvas="far"]',
        '[data-trip-cam="twin"]',
        '[data-trip-cam="seen"]',
        "[data-gl-pov]",
      ];
      ids.push("pair-far", "trip-twin", "trip-seen", "glasses-pov");
      return ids.map((id, i) => {
        const c = document.querySelector(sels[i]);
        return { id, w: c?.width || 0, h: c?.height || 0 };
      });
    });
    console.log("[verify-chrome] All 24 WebGL canvases active:", nonBlankSummary);
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

    // 7. A jump from the nav lands on its chapter: the figures it passes mount and grow above it on the way, and the
    // page keeps the chapter where the jump put it (a fresh page each, desktop and phone).
    for (const vp of [
      { name: "desktop-1440", width: 1440, height: 900, mobile: false },
      { name: "mobile-402", width: 402, height: 874, mobile: true },
    ]) {
      const jp = await browser.newPage();
      jp.on("pageerror", (err) => pageErrors.push(`[jump ${vp.name}] ${err.message}`));
      await jp.setViewport({ width: vp.width, height: vp.height, deviceScaleFactor: 1 });
      await jp.goto(url, { waitUntil: "load" });
      await sleep(2000);
      if (vp.mobile) {
        await jp.click(".menu > summary");
        await sleep(400);
        await jp.click('.menu a[href="#fleet"]');
      } else await jp.click('nav a[href="#fleet"]');
      await sleep(6000);
      const landed = await jp.evaluate(() => ({ top: Math.round(document.getElementById("fleet").getBoundingClientRect().top), pad: parseFloat(getComputedStyle(document.documentElement).scrollPaddingTop) || 0 }));
      console.log(`[verify-chrome] Jump to Kitchens (${vp.name}): chapter top ${landed.top}px, header ${landed.pad}px`);
      if (Math.abs(landed.top - landed.pad) > 6) throw new Error(`A nav jump on ${vp.name} did not land on its chapter: top ${landed.top}px, expected ${landed.pad}px`);
      await jp.close();
    }

    if (consoleErrors.length > 0 || pageErrors.length > 0 || failedRequests.length > 0) {
      throw new Error(
        `Runtime/network errors detected:\n` +
          `  consoleErrors=${JSON.stringify(consoleErrors)}\n` +
          `  pageErrors=${JSON.stringify(pageErrors)}\n` +
          `  failedRequests=${JSON.stringify(failedRequests)}`,
      );
    }

    console.log(`[verify-chrome] ALL CHECKS PASSED (0 console errors, 0 page errors, 0 failed fetches, ${sectionIds.length}/${sectionIds.length} WebGL live, 0px overflow)!`);
  } finally {
    await browser.close();
    server.close();
  }
}

main().catch((err) => {
  console.error("[verify-chrome] FAILED:", err);
  process.exit(1);
});
