import { createHash } from "node:crypto";
import { existsSync, mkdirSync, readdirSync, readFileSync, renameSync, statSync } from "node:fs";
import { join, relative, resolve, sep } from "node:path";
import { fileURLToPath } from "node:url";
import { defineConfig, type Plugin } from "vite";

// Public URL for the share cards (%VITE_SITE_URL% in index.html); the Cloudflare build sets
// VITE_SITE_URL=https://frankensim.org/music/ and --base /music/.
process.env.VITE_SITE_URL ??= "https://frankensim-struck.vercel.app/";

const here = fileURLToPath(new URL(".", import.meta.url));

/** SHA-256 over every file under `dir` (relative path and bytes, in path order), first 12 hex digits. */
function treeHash(dir: string): string {
  const h = createHash("sha256");
  const walk = (d: string): void => {
    for (const name of readdirSync(d).sort()) {
      const p = join(d, name);
      if (statSync(p).isDirectory()) walk(p);
      else h.update(relative(dir, p).split(sep).join("/")).update("\0").update(readFileSync(p)).update("\0");
    }
  };
  walk(dir);
  return h.digest("hex").slice(0, 12);
}

/**
 * The page's data, audio, images, replays and MIDI rolls (public/{data,audio,img,viz,midi}) keep their
 * paths, served with a one-hour cache. A build gives the page each file's content hash
 * (import.meta.env.VITE_ASSET_HASHES, path -> 10 hex digits), and util.asset() asks for `path?v=<hash>`:
 * a file that changed gets a new URL, so right after a deploy no visitor pairs new code and engines with
 * an hour-old strandbar.json or recording; an unchanged file keeps its URL and its cache. The dev server
 * serves plain paths.
 */
function assetHashes(): Plugin {
  return {
    name: "struck:asset-hashes",
    apply: "build",
    config() {
      const hashes: Record<string, string> = {};
      const walk = (d: string): void => {
        for (const name of readdirSync(d).sort()) {
          const p = join(d, name);
          if (statSync(p).isDirectory()) walk(p);
          else hashes[relative(join(here, "public"), p).split(sep).join("/")] = createHash("sha256").update(readFileSync(p)).digest("hex").slice(0, 10);
        }
      };
      for (const dir of ["data", "audio", "img", "viz", "midi"]) if (existsSync(join(here, "public", dir))) walk(join(here, "public", dir));
      return { define: { "import.meta.env.VITE_ASSET_HASHES": JSON.stringify(hashes) } };
    },
  };
}

/**
 * The live engines (public/live: module Workers, their wasm-bindgen glue and wasm, the level table,
 * the piano boards, the mixer worklet) are static files the bundler never sees, and each Worker finds
 * its glue, wasm and data relative to its own URL. A build serves the whole tree under one versioned
 * directory, live/<hash of its contents>/, and the page's LIVE_BASE_URL names it
 * (import.meta.env.VITE_LIVE_VERSION). So a page always loads a matching set of Workers, glue and
 * wasm, and every one of those URLs can be cached immutably (cloudflare/_headers, vercel.json).
 * Without the version a stale cached wasm could meet a newer engine.js for as long as it was cached.
 * The dev server keeps serving public/live at live/.
 */
function versionedLive(): Plugin {
  const source = join(here, "public", "live");
  let version = "";
  let outDir = "";
  let log: (msg: string) => void = () => {};
  return {
    name: "struck:versioned-live",
    apply: "build",
    config() {
      version = treeHash(source);
      return { define: { "import.meta.env.VITE_LIVE_VERSION": JSON.stringify(version) } };
    },
    configResolved(config) {
      outDir = resolve(config.root, config.build.outDir);
      log = (msg) => config.logger.info(msg);
    },
    closeBundle() {
      const live = join(outDir, "live");
      if (!existsSync(live)) throw new Error(`versioned-live: ${live} was not copied from public/`);
      const parked = join(outDir, `.live-${version}`);
      renameSync(live, parked);
      mkdirSync(live);
      renameSync(parked, join(live, version));
      // what ships must be exactly what the version was computed from (a file changed during the build
      // would otherwise ship under a name that promises other bytes)
      const shipped = treeHash(join(live, version));
      if (shipped !== version || treeHash(source) !== version) {
        throw new Error(`versioned-live: live/${version}/ does not hold exactly public/live as hashed (shipped ${shipped}): public/live changed during the build, or the output directory was not emptied; build again`);
      }
      log(`live engines served from live/${version}/`);
    },
  };
}

// The live engines run in module Workers that feed an AudioWorklet sink; they
// are plain static files under public/live and never pass through the bundler.
export default defineConfig({
  plugins: [versionedLive(), assetHashes()],
  build: { target: "es2022", sourcemap: false, chunkSizeWarningLimit: 900 },
  // The page imports the drum and cymbal descriptions verbatim from the repo (?raw). The watcher skips the
  // engines' cargo targets (GBs) and the build outputs: watching them kept an idle dev server near a full core.
  server: {
    port: 5199,
    strictPort: false,
    fs: { allow: ["../.."] },
    watch: { ignored: ["**/target/**", "**/dist/**", "**/dist-cf/**", "**/.wrangler/**", "**/.vercel/**"] },
  },
});
