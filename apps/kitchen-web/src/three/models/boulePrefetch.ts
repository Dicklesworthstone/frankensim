/**
 * Start the boule's texel layers in a worker as the page loads, so the oven chapter's loaf (smallwares' boule) finds
 * them made instead of spending the main thread on them. Resolves when they have arrived, or at once (or after
 * `waitMs`) if no worker can make them: the builder then makes them itself, the same arrays either way. Asking again
 * for the same detail returns the same promise.
 */
import { BOULE_TW, provideBouleLayers, type BouleLayers } from "./bouleTexels";
import type { BouleAsk } from "./bouleTexels.worker";

const asked = new Map<number, Promise<void>>();

export function prefetchBouleLayers(detail: 0 | 1 | 2, waitMs = 10000): Promise<void> {
  const ask: BouleAsk = { TW: BOULE_TW[detail] };
  const known = asked.get(ask.TW);
  if (known) return known;
  const p = new Promise<void>((resolve) => {
    let worker: Worker;
    try {
      worker = new Worker(new URL("./bouleTexels.worker.ts", import.meta.url), { type: "module" });
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
    worker.onmessage = (e: MessageEvent<BouleLayers>) => {
      provideBouleLayers(e.data);
      done();
    };
    worker.onerror = done;
    worker.postMessage(ask);
  });
  asked.set(ask.TW, p);
  return p;
}
