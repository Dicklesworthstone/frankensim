/**
 * Start the ribeye's tissue in a worker as the page loads, so the stove chapter's mount (smallwares' steak) finds it
 * made instead of spending the main thread on it. Resolves when it has arrived, or at once (or after `waitMs`) if
 * no worker can make it: the builder then makes the tissue itself, the same arrays either way.
 */
import { provideSteakTissue, type SteakTissue } from "./steakTissue";
import type { SteakAsk } from "./steakTissue.worker";

export function prefetchSteakTissue(ask: SteakAsk, waitMs = 10000): Promise<void> {
  return new Promise((resolve) => {
    let worker: Worker;
    try {
      worker = new Worker(new URL("./steakTissue.worker.ts", import.meta.url), { type: "module" });
    } catch {
      resolve();
      return;
    }
    let timer = 0;
    const done = () => {
      clearTimeout(timer);
      worker.terminate();
      resolve();
    };
    timer = window.setTimeout(done, waitMs);
    worker.onmessage = (e: MessageEvent<SteakTissue>) => {
      provideSteakTissue(ask.w, ask.d, ask.h, ask.detail, e.data);
      done();
    };
    worker.onerror = done;
    worker.postMessage(ask);
  });
}
