// One headless boot of the REAL built app (bead frankensim-xsz8b).
// Fresh browser per boot: worker/wasm state must be cold for the
// determinism claim to mean anything.
import { accessSync } from "node:fs";
import { JsonlCapture } from "./lib.mjs";
import puppeteer from "puppeteer-core";

const CHROME_CANDIDATES = [
  process.env.WF_CHROME_BIN,
  "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
  "/usr/bin/google-chrome",
  "/usr/bin/google-chrome-stable",
  "/usr/bin/chromium",
  "/usr/bin/chromium-browser",
].filter(Boolean);

export function resolveChromeBin() {
  for (const candidate of CHROME_CANDIDATES) {
    try {
      // accessSync is sync-on-purpose: boot path wants fail-fast.
      accessSync(candidate);
      return candidate;
    } catch {
      /* try next */
    }
  }
  return null;
}

export class BootRefusal extends Error {
  constructor(code, message) {
    super(message);
    this.code = code;
  }
}

/**
 * Boot `baseUrl + query` in a fresh headless Chrome, capture console JSONL,
 * and return when the run terminates (results card visible) or a typed
 * refusal fires. Never leaves a browser alive on failure.
 */
export async function bootOnce({ baseUrl, query, timeoutMs = 60000, chromeBin }) {
  const bin = chromeBin ?? resolveChromeBin();
  if (!bin) {
    throw new BootRefusal(
      "CHROME_NOT_FOUND",
      "No Chrome/Chromium binary found. Set WF_CHROME_BIN to the executable path.",
    );
  }
  const capture = new JsonlCapture();
  const browser = await puppeteer.launch({
    executablePath: bin,
    headless: true,
    args: ["--disable-dev-shm-usage", "--window-size=1280,800"],
    defaultViewport: { width: 1280, height: 800 },
  });
  try {
    const page = await browser.newPage();
    page.on("console", (msg) => {
      capture.push(msg.text());
    });
    page.on("pageerror", (err) => {
      capture.push(JSON.stringify({ stage: "page-error", message: String(err?.message ?? err) }));
    });
    const url = baseUrl.replace(/\/?$/, "/") + "?" + query;
    await page.goto(url, { waitUntil: "load", timeout: timeoutMs });
    // Terminal OR refusal, whichever comes first. The card is the app's own
    // end-of-run signal; the JSONL line is cross-checked by extractReceipts.
    await Promise.race([
      page.waitForFunction(
        () => {
          const el = document.getElementById("wf-results-card");
          return el !== null && el.style.display !== "none" && (el.textContent ?? "").length > 0;
        },
        { timeout: timeoutMs, polling: 250 },
      ),
      new Promise((_, reject) =>
        setTimeout(() => reject(new BootRefusal("RUN_TIMEOUT", `run did not terminate within ${timeoutMs}ms`)), timeoutMs),
      ),
    ]).catch(async (error) => {
      const refusals = capture.lines().filter((r) => r.stage === "sim-refusal");
      if (refusals.length > 0) {
        throw new BootRefusal("SIM_REFUSAL", JSON.stringify(refusals));
      }
      throw error;
    });
    return { lines: capture.lines(), captureRefusals: capture.refusals() };
  } finally {
    await browser.close().catch(() => {});
  }
}
