import { defineConfig } from "vite";

// COOP/COEP: cross-origin isolation is required for SharedArrayBuffer (the
// seqlock state ring + leased field snapshots, plan §4.3/§7.3). The app must
// ALSO work without isolation (degraded transferable-pool mode) — the runtime
// banner in src/capability.ts reports which mode is actually live, so these
// headers are an enhancement, never a hidden requirement.
const crossOriginIsolationHeaders = {
  "Cross-Origin-Opener-Policy": "same-origin",
  "Cross-Origin-Embedder-Policy": "require-corp",
};

export default defineConfig({
  server: { headers: crossOriginIsolationHeaders },
  preview: { headers: crossOriginIsolationHeaders },
  build: {
    target: "es2022",
    sourcemap: true,
  },
});
