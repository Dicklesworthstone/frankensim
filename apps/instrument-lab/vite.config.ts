import { defineConfig } from "vite";

// Public URL for the share cards (%VITE_SITE_URL% in index.html); the Cloudflare build sets
// VITE_SITE_URL=https://frankensim.org/music/ and --base /music/.
process.env.VITE_SITE_URL ??= "https://frankensim-struck.vercel.app/";

// The live engines run in module Workers that feed an AudioWorklet sink; they
// are plain static files under public/live and never pass through the bundler.
export default defineConfig({
  build: { target: "es2022", sourcemap: false, chunkSizeWarningLimit: 900 },
  // The page imports the drum and cymbal descriptions verbatim from the repo (?raw).
  server: { port: 5199, strictPort: false, fs: { allow: ["../.."] } },
});
