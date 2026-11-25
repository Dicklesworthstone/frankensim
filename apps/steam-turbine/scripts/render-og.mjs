#!/usr/bin/env node
/**
 * Render the share cards from the live turbine scene. Run after `npm run build:cf`, whenever the scene,
 * the engine or the card text changes:
 *
 *   node scripts/render-og.mjs
 *
 * Serves dist-cf, opens /turbine/ in Chrome with the GPU, lets the engine run, switches to the half
 * cutaway with the temperature field from the HP-and-valves camera, and screenshots the 3D viewport
 * at 2400 x 1260. The shot is the background of a typeset card. The four numbers on the card and the
 * colour-scale ticks are read from the page at the moment of capture, so they are engine output, not
 * typed in. Writes public/og.png (1200 x 630) and public/twitter.png (1200 x 600). Then rerun
 * `npm run build:cf` so the cards are staged.
 */

import fs from "node:fs";
import http from "node:http";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const APP_ROOT = path.resolve(__dirname, "..");
const DIST_CF = path.join(APP_ROOT, "dist-cf");
const ASSETS = path.join(DIST_CF, "turbine", "assets");
const OUT = [
  { file: path.join(APP_ROOT, "public", "og.png"), width: 1200, height: 630 },
  { file: path.join(APP_ROOT, "public", "twitter.png"), width: 1200, height: 600 },
];

const require = createRequire(import.meta.url);
const puppeteer = require(path.resolve(APP_ROOT, "../kitchen-studio/node_modules/puppeteer-core"));
const CHROME_PATH = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const MIME = { ".html": "text/html; charset=utf-8", ".js": "application/javascript", ".css": "text/css", ".wasm": "application/wasm", ".woff2": "font/woff2", ".woff": "font/woff", ".svg": "image/svg+xml", ".png": "image/png" };
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

const asset = (re) => {
  const name = fs.readdirSync(ASSETS).find((f) => re.test(f));
  if (!name) throw new Error(`no asset matching ${re} in ${ASSETS}`);
  return `/turbine/assets/${name}`;
};

function serve(routes) {
  return new Promise((resolve) => {
    const server = http.createServer((req, res) => {
      const url = (req.url || "/").split("?")[0];
      const route = routes.get(url);
      if (route) {
        res.setHeader("Content-Type", route.type);
        res.end(route.body());
        return;
      }
      let p = path.join(DIST_CF, decodeURIComponent(url));
      if (fs.existsSync(p) && fs.statSync(p).isDirectory()) p = path.join(p, "index.html");
      if (!p.startsWith(DIST_CF) || !fs.existsSync(p)) {
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

const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);

function cardHtml(m, height) {
  const chips = [
    ["Shaft", m.rpm],
    ["Electrical output", m.mw],
    ["Condenser", m.cond],
    ["HP rotor metal", m.temp],
  ]
    .map(([label, value]) => `<div class="chip"><span class="chip-label">${esc(label)}</span><span class="chip-value">${esc(value)}</span></div>`)
    .join("");
  return `<!doctype html><html><head><meta charset="utf-8"><style>
@font-face { font-family: "Inter"; src: url(${m.fontInter}) format("woff2"); font-weight: 100 900; }
@font-face { font-family: "JBM"; src: url(${m.fontMono}) format("woff2"); font-weight: 400; }
* { box-sizing: border-box; margin: 0; padding: 0; }
html, body { width: 1200px; height: ${height}px; overflow: hidden; background: #03060a; }
body { position: relative; font-family: "Inter", system-ui, sans-serif; color: #f4f8fc; }
.render { position: absolute; left: 0; top: -${height === 600 ? 64 : 58}px; width: 1200px; height: 630px; }
.shade { position: absolute; inset: 0; background:
  linear-gradient(to top, rgba(3,6,10,0.97) 0%, rgba(3,6,10,0.90) 25%, rgba(3,6,10,0.45) 47%, rgba(3,6,10,0) 63%),
  linear-gradient(to bottom, rgba(3,6,10,0.78) 0%, rgba(3,6,10,0) 20%),
  linear-gradient(to left, rgba(3,6,10,0.6) 0%, rgba(3,6,10,0) 24%),
  radial-gradient(ellipse 80% 90% at 62% 36%, rgba(0,0,0,0) 55%, rgba(0,0,0,0.5) 100%); }
.glow { position: absolute; left: 520px; top: 120px; width: 520px; height: 260px; background: radial-gradient(ellipse at center, rgba(255,190,90,0.16), rgba(255,190,90,0) 70%); mix-blend-mode: screen; }
.brand { position: absolute; left: 52px; top: 40px; display: flex; align-items: center; gap: 14px; }
.brand svg { width: 40px; height: 40px; filter: drop-shadow(0 0 12px rgba(125,211,252,0.35)); }
.brand-name { font-size: 25px; font-weight: 750; letter-spacing: -0.01em; }
.brand-sep { width: 1px; height: 22px; background: rgba(255,255,255,0.25); }
.brand-tag { font-size: 15px; font-weight: 600; letter-spacing: 0.16em; text-transform: uppercase; color: #9fb4c8; }
.scale { position: absolute; right: 48px; top: 36px; width: 290px; padding: 12px 16px 11px; border-radius: 12px;
  background: rgba(6,11,17,0.62); border: 1px solid rgba(160,200,235,0.18); backdrop-filter: blur(6px); }
.scale-head { display: flex; justify-content: space-between; font-size: 13px; font-weight: 600; color: #c9d7e4; letter-spacing: 0.02em; }
.scale-bar { height: 10px; border-radius: 5px; margin: 9px 0 7px; background: ${m.cbGradient}; box-shadow: 0 0 14px rgba(255,170,80,0.25); }
.scale-ticks { display: flex; justify-content: space-between; font-family: "JBM", monospace; font-size: 12.5px; color: #9fb0c0; }
.copy { position: absolute; left: 52px; bottom: ${height === 600 ? 38 : 46}px; width: 1096px; }
.kicker { font-size: 15px; font-weight: 650; letter-spacing: 0.2em; text-transform: uppercase; color: #7dd3fc; }
.title { margin-top: 8px; font-size: 86px; line-height: 0.98; font-weight: 820; letter-spacing: -0.035em;
  background: linear-gradient(100deg, #ffffff 0%, #ffffff 46%, #ffd79a 74%, #ffb45c 100%); -webkit-background-clip: text; background-clip: text; color: transparent; }
.sub { margin-top: 14px; width: 840px; font-size: 21px; line-height: 1.42; font-weight: 450; color: #bccad8; }
.row { margin-top: 22px; display: flex; align-items: flex-end; justify-content: space-between; }
.chips { display: flex; gap: 10px; }
.chip { display: flex; flex-direction: column; gap: 3px; padding: 9px 14px 10px; border-radius: 10px; background: rgba(125,211,252,0.07); border: 1px solid rgba(125,211,252,0.22); }
.chip-label { font-size: 11.5px; font-weight: 650; letter-spacing: 0.12em; text-transform: uppercase; color: #8aa3b8; }
.chip-value { font-family: "JBM", monospace; font-size: 21px; color: #f1f6fb; }
.url { font-family: "JBM", monospace; font-size: 19px; color: #7dd3fc; padding-bottom: 6px; }
</style></head><body>
<img class="render" src="/__og/render.png" alt="">
<div class="glow"></div>
<div class="shade"></div>
<div class="brand">
  <svg viewBox="0 0 32 32"><rect width="32" height="32" rx="7" fill="#162636"/><circle cx="16" cy="16" r="10" fill="none" stroke="#7dd3fc" stroke-width="2"/><circle cx="16" cy="16" r="3.5" fill="#eab663"/><path d="M16 6v6M16 20v6M6 16h6M20 16h6M9 9l4.2 4.2M18.8 18.8L23 23M23 9l-4.2 4.2M13.2 18.8L9 23" stroke="#7dd3fc" stroke-width="1.8" stroke-linecap="round"/></svg>
  <span class="brand-name">FrankenSim</span><span class="brand-sep"></span><span class="brand-tag">Live simulation</span>
</div>
<div class="scale">
  <div class="scale-head"><span>${esc(m.cbTitle)}</span><span>${esc(m.cbUnit)}</span></div>
  <div class="scale-bar"></div>
  <div class="scale-ticks"><span>${esc(m.cbMin)}</span><span>${esc(m.cbMid)}</span><span>${esc(m.cbMax)}</span></div>
</div>
<div class="copy">
  <div class="kicker">Rust physics, compiled to WebAssembly</div>
  <div class="title">Steam Turbine</div>
  <div class="sub">A 50 MW turbine-generator stepped live in the browser: IAPWS-IF97 steam, a coupled plant graph, rotor critical speeds and a 2-D thermal field.</div>
  <div class="row"><div class="chips">${chips}</div><div class="url">frankensim.org/turbine &#8599;</div></div>
</div>
</body></html>`;
}

async function captureScene(browser, port) {
  const page = await browser.newPage();
  const errors = [];
  page.on("pageerror", (e) => errors.push(String(e)));
  // 2400 x 1260 viewport under the 54 px top bar: twice the card, the card's aspect.
  await page.setViewport({ width: 2400, height: 1314, deviceScaleFactor: 1 });
  await page.goto(`http://127.0.0.1:${port}/turbine/`, { waitUntil: "networkidle0", timeout: 90000 });
  await page.waitForFunction(() => !document.getElementById("boot-overlay"), { timeout: 120000 });
  await page.click("#btn-zen");
  await sleep(6000); // let the plant settle from its initial state
  await page.evaluate(() => {
    document.querySelector('[data-shader="temp"]').click();
    document.querySelector('[data-cam="hp_valves"]').click();
    const q = document.getElementById("btn-quick-cutaway");
    for (let i = 0; i < 3 && !q.textContent.includes("Half"); i++) q.click();
  });
  await sleep(3500); // camera flight
  await page.keyboard.press("-"); // one notch back (1.25x distance) so the whole HP and LP blade path fits
  await sleep(1500);
  const m = await page.evaluate(() => {
    const t = (id) => (document.getElementById(id)?.textContent || "").trim();
    const kpi = (id) => {
      const el = document.getElementById(id);
      const unit = el?.querySelector("small")?.textContent || "";
      const num = (el?.textContent || "").replace(unit, "").trim();
      return `${num} ${unit}`.trim();
    };
    return {
      rpm: kpi("kpi-rpm"), mw: kpi("kpi-mw"), cond: kpi("kpi-cond"), temp: kpi("kpi-temp"),
      cbTitle: t("cb-title"), cbUnit: t("cb-unit"), cbMin: t("cb-min"), cbMid: t("cb-mid"), cbMax: t("cb-max"),
      cbGradient: getComputedStyle(document.getElementById("cb-bar")).backgroundImage,
    };
  });
  await page.addStyleTag({ content: "body * { visibility: hidden !important; } #viewport, #viewport * { visibility: visible !important; }" });
  await sleep(300);
  const box = await page.$eval("#viewport", (el) => {
    const r = el.getBoundingClientRect();
    return { x: r.x, y: r.y, width: r.width, height: r.height };
  });
  const png = await page.screenshot({ type: "png", clip: box });
  await page.close();
  if (errors.length) throw new Error(`page errors during capture: ${errors.join("; ")}`);
  if (!m.cbGradient || m.cbGradient === "none") throw new Error("colour scale gradient not found on the page");
  console.log(`[render-og] scene ${box.width}x${box.height}; ${m.rpm} | ${m.mw} | ${m.cond} | ${m.temp} | scale ${m.cbMin}..${m.cbMax}`);
  return { png: Buffer.from(png), m };
}

async function main() {
  if (!fs.existsSync(path.join(DIST_CF, "turbine", "index.html"))) throw new Error("dist-cf/turbine/index.html missing: run `npm run build:cf` first.");
  const routes = new Map();
  const { server, port } = await serve(routes);
  const browser = await puppeteer.launch({
    executablePath: CHROME_PATH,
    headless: "new",
    args: ["--use-gl=angle", "--use-angle=metal", "--enable-gpu", "--ignore-gpu-blocklist"],
  });
  try {
    const { png, m } = await captureScene(browser, port);
    m.fontInter = asset(/^inter-latin-wght-normal-.*\.woff2$/);
    m.fontMono = asset(/^jetbrains-mono-latin-400-normal-.*\.woff2$/);
    routes.set("/__og/render.png", { type: "image/png", body: () => png });
    for (const { file, width, height } of OUT) {
      routes.set("/__og/card.html", { type: "text/html; charset=utf-8", body: () => cardHtml(m, height) });
      const page = await browser.newPage();
      await page.setViewport({ width, height, deviceScaleFactor: 1 });
      await page.goto(`http://127.0.0.1:${port}/__og/card.html`, { waitUntil: "networkidle0", timeout: 30000 });
      await page.evaluate(() => document.fonts.ready);
      await sleep(200);
      await page.screenshot({ path: file, type: "png", clip: { x: 0, y: 0, width, height } });
      await page.close();
      console.log(`[render-og] ${path.relative(APP_ROOT, file)}  ${width}x${height}  ${(fs.statSync(file).size / 1024).toFixed(0)} KiB`);
    }
  } finally {
    await browser.close();
    server.close();
  }
}

main().catch((err) => {
  console.error("[render-og] FAILED:", err);
  process.exit(1);
});
