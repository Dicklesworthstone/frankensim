#!/usr/bin/env node
/**
 * Render the hero's loading posters from the live 3D kitchen, so the first thing a visitor sees matches the scene
 * that replaces it. Run after `bun run build:cf`, whenever the models or the hero camera change:
 *
 *   node scripts/render-posters.mjs
 *
 * Serves dist-cf, opens the page in Chrome with the GPU, waits until the mall kitchen is built and the camera has
 * settled, hides everything on the page but the 3D canvas, and screenshots the hero stage as WebP:
 *   ../../assets/kitchens/web_posters/hero-desktop.webp   (1920 px wide viewport)
 *   ../../assets/kitchens/web_posters/hero-mobile.webp    (390 px wide viewport at 2.77x, about 1080 px wide)
 * and the 1200x630 share card (the top of the page as a visitor first sees it) as public/og-card.png.
 * Then rerun `bun run build:cf` so the new images are staged.
 */

import fs from "node:fs";
import http from "node:http";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const APP_ROOT = path.resolve(__dirname, "..");
const DIST_CF = path.join(APP_ROOT, "dist-cf");
const OUT_DIR = path.resolve(APP_ROOT, "../../assets/kitchens/web_posters");
const OG_CARD = path.join(APP_ROOT, "public", "og-card.png");

const require = createRequire(import.meta.url);
const puppeteer = require(path.resolve(APP_ROOT, "../kitchen-studio/node_modules/puppeteer-core"));
const CHROME_PATH = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";

const MIME = { ".html": "text/html; charset=utf-8", ".js": "application/javascript", ".css": "text/css", ".svg": "image/svg+xml", ".png": "image/png", ".webp": "image/webp", ".woff2": "font/woff2", ".json": "application/json", ".glb": "model/gltf-binary" };

function serve(rootDir) {
  return new Promise((resolve) => {
    const server = http.createServer((req, res) => {
      let p = path.join(rootDir, decodeURIComponent((req.url || "/").split("?")[0]).replace(/^\/+/, ""));
      if (fs.existsSync(p) && fs.statSync(p).isDirectory()) p = path.join(p, "index.html");
      if (!p.startsWith(rootDir) || !fs.existsSync(p)) {
        res.statusCode = 404;
        res.end();
        return;
      }
      res.setHeader("Content-Type", MIME[path.extname(p).toLowerCase()] || "application/octet-stream");
      fs.createReadStream(p).pipe(res);
    });
    server.listen(0, "127.0.0.1", () => resolve({ server, port: server.address().port }));
  });
}

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function shoot(browser, url, viewport, file, { stageOnly = true } = {}) {
  const page = await browser.newPage();
  const errors = [];
  page.on("pageerror", (e) => errors.push(String(e)));
  await page.setViewport(viewport);
  await page.goto(url, { waitUntil: "networkidle0", timeout: 60000 });
  await page.waitForFunction(() => (document.querySelector("[data-hero-status]")?.textContent || "").includes("draw calls"), { timeout: 90000 });
  // The camera settles in about half a second. Shoot soon after, so the poster matches the first live frames it
  // cross-fades into (the cooks have barely moved from where they start).
  await sleep(2000);
  let box = { x: 0, y: 0, width: viewport.width, height: viewport.height };
  if (stageOnly) {
    // Only the 3D canvas: on wide screens the headline and panel sit over the stage and must not be baked in.
    await page.addStyleTag({ content: "body * { visibility: hidden !important; } .hero-stage canvas { visibility: visible !important; opacity: 1 !important; transition: none !important; }" });
    await sleep(400);
    box = await page.$eval(".hero-stage", (el) => {
      const r = el.getBoundingClientRect();
      return { x: r.x, y: r.y + window.scrollY, width: r.width, height: r.height };
    });
  }
  const type = file.endsWith(".png") ? "png" : "webp";
  await page.screenshot({ path: file, type, ...(type === "webp" ? { quality: 82 } : {}), clip: box, captureBeyondViewport: false });
  await page.close();
  if (errors.length) throw new Error(`page errors while rendering ${file}: ${errors.join("; ")}`);
  console.log(`[render-posters] ${path.relative(APP_ROOT, file)}  ${Math.round(box.width * viewport.deviceScaleFactor)}x${Math.round(box.height * viewport.deviceScaleFactor)}  ${(fs.statSync(file).size / 1024).toFixed(0)} KiB`);
}

async function main() {
  if (!fs.existsSync(path.join(DIST_CF, "kitchen-lab", "index.html"))) throw new Error("dist-cf/kitchen-lab/index.html missing: run `bun run build:cf` first.");
  fs.mkdirSync(OUT_DIR, { recursive: true });
  const { server, port } = await serve(DIST_CF);
  const url = `http://127.0.0.1:${port}/kitchen-lab/`;
  const browser = await puppeteer.launch({
    executablePath: CHROME_PATH,
    headless: "new",
    args: ["--use-gl=angle", "--use-angle=metal", "--enable-gpu", "--ignore-gpu-blocklist"],
  });
  try {
    await shoot(browser, url, { width: 1920, height: 1100, deviceScaleFactor: 1 }, path.join(OUT_DIR, "hero-desktop.webp"));
    await shoot(browser, url, { width: 390, height: 844, deviceScaleFactor: 2.77, isMobile: true, hasTouch: true }, path.join(OUT_DIR, "hero-mobile.webp"));
    await shoot(browser, url, { width: 1200, height: 630, deviceScaleFactor: 1 }, OG_CARD, { stageOnly: false });
  } finally {
    await browser.close();
    server.close();
  }
}

main().catch((err) => {
  console.error("[render-posters] FAILED:", err);
  process.exit(1);
});
