// Serves the BUILT app for e2e boots (bead frankensim-xsz8b).
// Leaf A uses `vite preview` (COOP/COEP on => SAB transport path).
// Leaf b will add a headerless static server for the degraded row.
import { spawn } from "node:child_process";
import { mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import net from "node:net";

export class ServeRefusal extends Error {
  constructor(code, message) {
    super(message);
    this.code = code;
  }
}

/** Grab an ephemeral free port (closed again immediately; local race acceptable). */
export function freePort() {
  return new Promise((resolve, reject) => {
    const srv = net.createServer();
    srv.unref();
    srv.on("error", reject);
    srv.listen(0, "127.0.0.1", () => {
      const { port } = srv.address();
      srv.close(() => resolve(port));
    });
  });
}

async function waitForHttp(url, deadlineMs, label) {
  const start = Date.now();
  while (Date.now() - start < deadlineMs) {
    try {
      const res = await fetch(url);
      if (res.ok) return;
    } catch {
      /* not up yet */
    }
    await new Promise((r) => setTimeout(r, 150));
  }
  throw new ServeRefusal("SERVE_TIMEOUT", `${label} never became reachable at ${url}`);
}

/**
 * Start `vite preview` for apps/wright-flyer with isolation headers.
 * Resolves { baseUrl, stop() }; stop() kills the child and never throws.
 */
export async function startVitePreview({ timeoutMs = 30000 } = {}) {
  const appRoot = path.resolve(import.meta.dirname, "..");
  const port = await freePort();
  const child = spawn(
    process.execPath,
    [path.join(appRoot, "node_modules", "vite", "bin", "vite.js"), "preview", "--host", "127.0.0.1", "--port", String(port), "--strictPort"],
    { cwd: appRoot, stdio: ["ignore", "pipe", "pipe"] },
  );
  let stderrTail = "";
  child.stderr.on("data", (d) => {
    stderrTail = (stderrTail + d.toString()).slice(-2000);
  });
  const baseUrl = `http://127.0.0.1:${port}/`;
  try {
    await waitForHttp(baseUrl, timeoutMs, "vite preview");
  } catch (error) {
    child.kill("SIGTERM");
    throw error;
  }
  return {
    baseUrl,
    port,
    stop() {
      child.kill("SIGTERM");
    },
    stderrTail: () => stderrTail,
  };
}

const MIME = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript",
  ".mjs": "text/javascript",
  ".css": "text/css",
  ".json": "application/json",
  ".map": "application/json",
  ".wasm": "application/wasm",
  ".webp": "image/webp",
  ".png": "image/png",
  ".svg": "image/svg+xml",
  ".ico": "image/x-icon",
  ".txt": "text/plain; charset=utf-8",
  ".woff2": "font/woff2",
  ".glb": "model/gltf-binary",
  ".gltf": "model/gltf+json",
};

/**
 * Serve the BUILT dist/ with NO isolation headers (leaf frankensim-y473e):
 * the app must degrade to the transferable-pool/postMessage transport and
 * still fly a full deterministic rail run. Deliberately plain node:http.
 */
export async function startStaticDist({ timeoutMs = 15000 } = {}) {
  const http = await import("node:http");
  const { readFile, stat } = await import("node:fs/promises");
  const appRoot = path.resolve(import.meta.dirname, "..");
  const distRoot = path.join(appRoot, "dist");
  const port = await freePort();
  const server = http.createServer(async (req, res) => {
    try {
      const urlPath = decodeURIComponent(new URL(req.url, "http://x").pathname);
      let rel = urlPath === "/" ? "index.html" : urlPath.replace(/^\/+/, "");
      let filePath = path.normalize(path.join(distRoot, rel));
      if (!filePath.startsWith(distRoot + path.sep) && filePath !== distRoot) {
        res.writeHead(403).end("forbidden");
        return;
      }
      try {
        const s = await stat(filePath);
        if (s.isDirectory()) filePath = path.join(filePath, "index.html");
      } catch {
        // SPA fallback: unknown extension-less routes get the shell.
        if (!path.extname(rel)) filePath = path.join(distRoot, "index.html");
      }
      const body = await readFile(filePath);
      const type = MIME[path.extname(filePath).toLowerCase()] ?? "application/octet-stream";
      res.writeHead(200, { "content-type": type });
      res.end(body);
    } catch {
      res.writeHead(404).end("not found");
    }
  });
  await new Promise((resolve, reject) => {
    server.once("error", reject);
    server.listen(port, "127.0.0.1", resolve);
  });
  const baseUrl = `http://127.0.0.1:${port}/`;
  await waitForHttp(baseUrl, timeoutMs, "static dist server");
  return {
    baseUrl,
    port,
    stop() {
      server.close();
    },
    stderrTail: () => "",
  };
}

/** Temporary scratch dir for run artifacts (never inside the repo tree). */
export function makeArtifactDir(prefix = "wf-e2e-") {
  return mkdtempSync(path.join(tmpdir(), prefix));
}
