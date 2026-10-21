import { createHash } from "node:crypto";
import { copyFileSync, existsSync, readdirSync, readFileSync, statSync } from "node:fs";
import { join, relative, resolve, sep } from "node:path";
import { fileURLToPath } from "node:url";
import { defineConfig, type Plugin } from "vite";

// The instruments-only review build: one page (instruments.html) and only its own data, figures and audio.
// tools/build-review.sh stages that data (tools/build_instruments.py output) in REVIEW_PUBLIC, runs this config with
// REVIEW_OUT as the output directory and VITE_REVIEW_SHA for the footer. vite.config.ts (the full site) is untouched.
const here = fileURLToPath(new URL(".", import.meta.url));
const staged = process.env.REVIEW_PUBLIC;
const outDir = process.env.REVIEW_OUT;
if (!staged || !outDir || !existsSync(join(staged, "data", "instruments", "index.json"))) {
  throw new Error("vite.instruments.config.ts: run it through tools/build-review.sh (REVIEW_PUBLIC must hold data/instruments/index.json, REVIEW_OUT names the output)");
}

/** util.asset() asks for `path?v=<content hash>` so a republished review never serves older data from a cache. */
function assetHashes(): Plugin {
  return {
    name: "struck:review-asset-hashes",
    config() {
      const hashes: Record<string, string> = {};
      const walk = (d: string): void => {
        for (const name of readdirSync(d).sort()) {
          const p = join(d, name);
          if (statSync(p).isDirectory()) walk(p);
          else hashes[relative(staged!, p).split(sep).join("/")] = createHash("sha256").update(readFileSync(p)).digest("hex").slice(0, 10);
        }
      };
      walk(join(staged!, "data"));
      return { define: { "import.meta.env.VITE_ASSET_HASHES": JSON.stringify(hashes) } };
    },
  };
}

/** The host may want index.html at the root: serve the same page under both names. */
function indexAlias(): Plugin {
  let dir = "";
  return {
    name: "struck:review-index-alias",
    configResolved(config) { dir = resolve(config.root, config.build.outDir); },
    closeBundle() { copyFileSync(join(dir, "instruments.html"), join(dir, "index.html")); },
  };
}

export default defineConfig({
  base: "./", // relative URLs: the page works wherever the host puts it
  publicDir: staged,
  plugins: [assetHashes(), indexAlias()],
  build: {
    outDir: resolve(outDir),
    emptyOutDir: true,
    target: "es2022",
    sourcemap: false,
    chunkSizeWarningLimit: 900,
    // The host caps the number of files, not their size: inline the small font subsets (JetBrains Mono, Inter's minor
    // subsets, and the duplicate .woff beside each .woff2) into the stylesheet instead of shipping ~20 more files.
    assetsInlineLimit: 30_000,
    rollupOptions: { input: resolve(here, "instruments.html") },
  },
});
