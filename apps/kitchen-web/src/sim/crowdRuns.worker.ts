/**
 * The crowd figure's runs, off the page's main thread: one crew in one kitchen, with bodies or as ghosts, over the
 * pooled mornings (runMornings).
 * The page posts the runs it wants, the one on screen first; each is posted back when done, and a new list replaces
 * whatever was still waiting (a run under way finishes first).
 */
import { KITCHENS, runMornings, type CrowdSet } from "./crewCrowd";

export interface CrowdAsk {
  kitchen: string;
  crew: number;
  ghosts: boolean;
}
export interface CrowdDone extends CrowdAsk {
  set: CrowdSet;
}

const scope = self as unknown as { onmessage: ((e: MessageEvent<CrowdAsk[]>) => void) | null; postMessage(m: CrowdDone): void };
const keyOf = (a: CrowdAsk) => `${a.kitchen}/${a.crew}/${a.ghosts ? "g" : "b"}`;
const done = new Set<string>();
let queue: CrowdAsk[] = [];
let busy = false;

function pump(): void {
  if (busy) return;
  busy = true;
  // A turn of the event loop first, so a newer list can take this one's place: the run is taken from the queue only
  // then, skipping any a list posted meanwhile still named but that is done already (else it would run twice).
  setTimeout(() => {
    let a = queue.shift();
    while (a && (done.has(keyOf(a)) || !KITCHENS[a.kitchen])) a = queue.shift();
    if (!a) {
      busy = false;
      return;
    }
    const set = runMornings(KITCHENS[a.kitchen]!, a.crew, a.ghosts);
    done.add(keyOf(a));
    scope.postMessage({ ...a, set });
    busy = false;
    pump();
  }, 0);
}

scope.onmessage = (e) => {
  queue = e.data.filter((a) => !done.has(keyOf(a)));
  pump();
};
