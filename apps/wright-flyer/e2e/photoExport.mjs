// Focused browser proof for E9.2a. Vite owns an ephemeral child port,
// Chrome exercises the real source module and download path, and the
// retained artifact is checked for a PNG signature.
// Repro: node e2e/photoExport.mjs

import { spawn } from "node:child_process";
import { readFile, stat } from "node:fs/promises";
import path from "node:path";
import puppeteer from "puppeteer-core";
import { makeArtifactDir } from "./serve.mjs";
import { resolveChromeBin } from "./boot.mjs";

const appRoot = path.resolve(import.meta.dirname, "..");
const artifactDir = makeArtifactDir("wf-photo-e2e-");
const filename = "wright-flyer-attract-mode.png";
const artifact = path.join(artifactDir, filename);

function startViteDev(timeoutMs = 30_000) {
  const child = spawn(
    process.execPath,
    [
      path.join(appRoot, "node_modules", "vite", "bin", "vite.js"),
      "--host",
      "127.0.0.1",
      "--port",
      "0",
      "--strictPort",
    ],
    { cwd: appRoot, stdio: ["ignore", "pipe", "pipe"] },
  );
  let output = "";
  const ready = new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      reject(new Error(`vite dev did not bind an ephemeral port: ${output.slice(-2000)}`));
    }, timeoutMs);
    const consume = (chunk) => {
      output = (output + chunk.toString()).slice(-4000);
      const plain = output.replace(/\x1b\[[0-9;]*m/g, "");
      const match = plain.match(/http:\/\/127\.0\.0\.1:(\d+)\//);
      if (match !== null) {
        clearTimeout(timer);
        resolve({
          baseUrl: `http://127.0.0.1:${match[1]}/`,
          stop() {
            child.kill("SIGTERM");
          },
        });
      }
    };
    child.stdout.on("data", consume);
    child.stderr.on("data", consume);
    child.once("error", (error) => {
      clearTimeout(timer);
      reject(error);
    });
    child.once("exit", (code) => {
      clearTimeout(timer);
      reject(new Error(`vite dev exited ${code}: ${output.slice(-2000)}`));
    });
  });
  return ready.catch((error) => {
    child.kill("SIGTERM");
    throw error;
  });
}

async function waitForArtifact(timeoutMs = 15_000) {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    try {
      const info = await stat(artifact);
      if (info.size > 8) return info;
    } catch {
      // Chrome has not committed the download yet.
    }
    await new Promise((resolve) => setTimeout(resolve, 100));
  }
  throw new Error(`photo download did not appear at ${artifact}`);
}

const chrome = resolveChromeBin();
if (chrome === null) {
  console.error(
    JSON.stringify({ suite: "wf-photo-e2e", verdict: "REFUSED", code: "CHROME_NOT_FOUND" }),
  );
  process.exit(5);
}

const server = await startViteDev();
const browser = await puppeteer.launch({
  executablePath: chrome,
  headless: true,
  args: ["--disable-dev-shm-usage", "--window-size=960,640"],
  defaultViewport: { width: 960, height: 640 },
});

try {
  const page = await browser.newPage();
  const rows = [];
  page.on("console", (message) => {
    try {
      rows.push(JSON.parse(message.text()));
    } catch {
      // Only structured app rows participate in the verdict.
    }
  });
  const cdp = await page.createCDPSession();
  await cdp.send("Page.setDownloadBehavior", { behavior: "allow", downloadPath: artifactDir });
  await page.goto(`${server.baseUrl}?demo=1`, { waitUntil: "load", timeout: 30_000 });
  await page.waitForFunction(
    () => {
      const canvas = document.querySelector("#app canvas");
      return canvas instanceof HTMLCanvasElement && canvas.width > 0 && canvas.height > 0;
    },
    { timeout: 30_000 },
  );
  await page.evaluate(() => {
    window.dispatchEvent(new KeyboardEvent("keydown", { code: "KeyI" }));
  });
  const info = await waitForArtifact();
  const bytes = await readFile(artifact);
  assertPng(bytes);
  const success = rows.find((row) => row.stage === "photo-export");
  const refusal = rows.find((row) => row.stage === "photo-export-refused");
  if (success === undefined || refusal !== undefined) {
    throw new Error(`missing success or observed refusal: ${JSON.stringify({ success, refusal })}`);
  }
  if (success.filename !== filename || success.bytes !== info.size) {
    throw new Error(
      `download/log mismatch: ${JSON.stringify({ success, fileBytes: info.size, filename })}`,
    );
  }
  console.log(
    JSON.stringify({
      suite: "wf-photo-e2e",
      verdict: "PASS",
      artifact,
      bytes: info.size,
      width: success.width,
      height: success.height,
    }),
  );
} finally {
  await browser.close();
  server.stop();
}

function assertPng(bytes) {
  const signature = [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a];
  if (signature.some((value, index) => bytes[index] !== value)) {
    throw new Error(`download is not a PNG: ${bytes.subarray(0, 8).toString("hex")}`);
  }
}
