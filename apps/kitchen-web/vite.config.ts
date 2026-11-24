import fs from "node:fs";
import path from "node:path";
import { defineConfig, type Plugin } from "vite";

const REPO_KITCHEN_ASSETS = path.resolve(__dirname, "../../assets/kitchens");

const MIME_BY_EXT: Record<string, string> = {
  ".glb": "model/gltf-binary",
  ".json": "application/json; charset=utf-8",
  ".png": "image/png",
  ".svg": "image/svg+xml",
  ".webp": "image/webp",
};

function serveRepoKitchenAssets(): Plugin {
  return {
    name: "serve-repo-kitchen-assets",
    configureServer(server) {
      server.middlewares.use((req, res, next) => {
        const rawUrl = req.url?.split("?")[0] ?? "";
        let rel = "";
        if (rawUrl.startsWith("/kitchen-lab/kitchen-assets/")) {
          rel = rawUrl.slice("/kitchen-lab/kitchen-assets/".length);
        } else if (rawUrl.startsWith("/kitchen-assets/")) {
          rel = rawUrl.slice("/kitchen-assets/".length);
        } else {
          return next();
        }
        const safeRel = path.normalize(rel).replace(/^(\.\.(\/|\\|$))+/, "");
        const abs = path.join(REPO_KITCHEN_ASSETS, safeRel);
        if (!abs.startsWith(REPO_KITCHEN_ASSETS) || !fs.existsSync(abs) || fs.statSync(abs).isDirectory()) {
          res.statusCode = 404;
          res.end("Not found");
          return;
        }
        const ext = path.extname(abs).toLowerCase();
        res.setHeader("Content-Type", MIME_BY_EXT[ext] ?? "application/octet-stream");
        res.setHeader("Cache-Control", "no-cache");
        fs.createReadStream(abs).pipe(res);
      });
    },
  };
}

export default defineConfig({
  plugins: [serveRepoKitchenAssets()],
  build: {
    target: "es2022",
    sourcemap: false,
    chunkSizeWarningLimit: 900,
    rollupOptions: {
      output: {
        manualChunks: {
          three: ["three"],
        },
      },
    },
  },
});
