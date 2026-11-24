#!/usr/bin/env node
/**
 * Stage Cloudflare Workers Assets bundle for frankensim.org/kitchen-lab:
 *   1. Run `tsc --noEmit`
 *   2. Run `vite build --base /kitchen-lab/ --outDir dist-cf/kitchen-lab --emptyOutDir`
 *   3. Copy the five kitchens' build manifests (`.scene3d.v1.json`, `.site_metadata.json`) and the hero's loading
 *      posters from `../../assets/kitchens/` into `dist-cf/kitchen-lab/kitchen-assets/`. The page builds every room,
 *      fixture and cook in code from those manifests, so no Blender `.glb` is staged.
 *   4. Copy `cloudflare/_headers` into `dist-cf/_headers`
 *   5. Run Law V-L16 zero-OEM-brand audit across staged HTML/JS/CSS
 */

import fs from "node:fs";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const APP_ROOT = path.resolve(__dirname, "..");
const REPO_ASSETS = path.resolve(APP_ROOT, "../../assets/kitchens");
const DIST_CF = path.join(APP_ROOT, "dist-cf");
const DIST_SITE = path.join(DIST_CF, "kitchen-lab");
const DIST_ASSETS = path.join(DIST_SITE, "kitchen-assets");

const EXEMPLAR_DIRS = [
  "mall_qsr_400",
  "freestanding_qsr_1450",
  "fast_casual_1100",
  "fine_dining_950",
  "commissary_catering_3200",
];

function copyFileEnsured(src, dst) {
  fs.mkdirSync(path.dirname(dst), { recursive: true });
  fs.copyFileSync(src, dst);
  return fs.statSync(dst).size;
}

function run(cmd, args) {
  execFileSync(cmd, args, { cwd: APP_ROOT, stdio: "inherit" });
}

console.log("[stage-cf] 1/5 Type-checking TypeScript (`tsc --noEmit`)...");
run("bunx", ["tsc", "--noEmit"]);

console.log("[stage-cf] 2/5 Building Vite production bundle (`--base /kitchen-lab/`)...");
run("bunx", [
  "vite",
  "build",
  "--base",
  "/kitchen-lab/",
  "--outDir",
  "dist-cf/kitchen-lab",
  "--emptyOutDir",
]);

console.log("[stage-cf] 3/5 Staging the kitchens' build manifests and the hero posters into dist-cf/kitchen-lab/kitchen-assets/...");
let copiedFiles = 0;
let copiedBytes = 0;

// The 5 kitchens' manifests: the hero builds each room from these (heroKitchen.ts fetches nothing else).
for (const dir of EXEMPLAR_DIRS) {
  const srcDir = path.join(REPO_ASSETS, "exemplars", dir);
  for (const suffix of [".scene3d.v1.json", ".site_metadata.json"]) {
    const fname = `${dir}${suffix}`;
    copiedBytes += copyFileEnsured(
      path.join(srcDir, fname),
      path.join(DIST_ASSETS, "exemplars", dir, fname),
    );
    copiedFiles += 1;
  }
}

// The hero's loading posters, rendered from the live 3D kitchen by scripts/render-posters.mjs (which needs a staged
// build to render from, so a first build may have none yet). index.html references both; `bun run verify` fails on
// the 404 if one is missing.
const posterSrcDir = path.join(REPO_ASSETS, "web_posters");
for (const entry of ["hero-desktop.webp", "hero-mobile.webp"]) {
  const src = path.join(posterSrcDir, entry);
  if (!fs.existsSync(src)) {
    console.warn(`[stage-cf] WARNING: missing hero poster ${src}; run node scripts/render-posters.mjs, then build again.`);
    continue;
  }
  copiedBytes += copyFileEnsured(src, path.join(DIST_ASSETS, "web_posters", entry));
  copiedFiles += 1;
}

console.log(
  `[stage-cf] Staged ${copiedFiles} asset files (${(copiedBytes / (1024 * 1024)).toFixed(2)} MiB total).`,
);

console.log("[stage-cf] 4/5 Copying cloudflare/_headers -> dist-cf/_headers...");
copyFileEnsured(
  path.join(APP_ROOT, "cloudflare", "_headers"),
  path.join(DIST_CF, "_headers"),
);

console.log("[stage-cf] 5/5 Running Law V-L16 zero-OEM-brand hygiene scan on dist-cf/kitchen-lab/...");
const OEM_REGEX =
  /\b(cambro|vulcan|hobart|pitco|henny\s*penny|blodgett|hoshizaki|manitowoc|vitamix|robot\s*coupe|rational\s*ag|garland|true\s*mfg|frialator|thermomix|turbo\s*air|traulsen|alto-shaam|cleveland\s*range|southbend)\b/gi;

const textFiles = [];
function collectText(dir) {
  for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, ent.name);
    if (ent.isDirectory()) collectText(p);
    else if (/\.(html|js|css|svg)$/i.test(ent.name)) textFiles.push(p);
  }
}
collectText(DIST_SITE);

const oemLeaks = [];
for (const f of textFiles) {
  const txt = fs.readFileSync(f, "utf8");
  const matches = txt.match(OEM_REGEX);
  if (matches) oemLeaks.push({ file: path.relative(APP_ROOT, f), matches });
}
if (oemLeaks.length > 0) {
  console.error("[stage-cf] FAIL Law V-L16 OEM brand leak detected:", oemLeaks);
  process.exit(1);
}

console.log(`[stage-cf] PASS! Staged dist-cf cleanly (${textFiles.length} bundle files scanned, 0 OEM leaks).`);
