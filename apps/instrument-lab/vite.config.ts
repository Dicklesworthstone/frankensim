import { defineConfig } from "vite";

// The live engines run in module Workers that feed an AudioWorklet sink; they
// are plain static files under public/live and never pass through the bundler.
export default defineConfig({
  build: { target: "es2022", sourcemap: false, chunkSizeWarningLimit: 900 },
  // The page imports the drum and cymbal descriptions verbatim from the repo (?raw).
  server: { port: 5199, strictPort: false, fs: { allow: ["../.."] } },
});
