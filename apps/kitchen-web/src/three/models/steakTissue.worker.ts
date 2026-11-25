/** The ribeye's tissue (steakTissue.ts), computed off the page's main thread and handed back without a copy. */
import { steakBuffers, steakTissue, type SteakTissue } from "./steakTissue";

export interface SteakAsk {
  w: number;
  d: number;
  h: number;
  detail: 0 | 1 | 2;
}

const scope = self as unknown as {
  onmessage: ((e: MessageEvent<SteakAsk>) => void) | null;
  postMessage(m: SteakTissue, transfer: ArrayBuffer[]): void;
};
scope.onmessage = (e) => {
  const { w, d, h, detail } = e.data;
  const t = steakTissue(w, d, h, detail);
  scope.postMessage(t, steakBuffers(t));
};
