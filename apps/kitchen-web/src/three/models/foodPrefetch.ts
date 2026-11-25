/**
 * Start the pan foods', the knife's and the wine label's texel layers (foodTexels.ts) in a worker as the page loads, so
 * the hero kitchen's first hotel pan and knife and the stove's bottle find them made instead of spending the main
 * thread on them. If no
 * worker can make them, or one is late, the builders make them themselves: the same arrays either way.
 */
import { provideKnifeTexels, provideLabelPaper, providePanFoodTexels } from "./foodTexels";
import type { FoodTexelsReply } from "./foodTexels.worker";

let started = false;

export function prefetchFoodTexels(waitMs = 15000): void {
  if (started) return;
  started = true;
  let worker: Worker;
  try {
    worker = new Worker(new URL("./foodTexels.worker.ts", import.meta.url), { type: "module" });
  } catch {
    return;
  }
  let left = 3;
  const done = () => {
    clearTimeout(timer);
    worker.terminate();
  };
  const timer = window.setTimeout(done, waitMs);
  worker.onmessage = (e: MessageEvent<FoodTexelsReply>) => {
    const r = e.data;
    if (r.kind === "pan") providePanFoodTexels(r.t);
    else if (r.kind === "knife") provideKnifeTexels(r.t);
    else provideLabelPaper(r.t);
    if (--left === 0) done();
  };
  worker.onerror = done;
  worker.postMessage(null);
}
