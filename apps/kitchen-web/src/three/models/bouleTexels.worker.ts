/** The boule's static texel layers (bouleTexels.ts), computed off the page's main thread and handed back without a copy. */
import { bouleBuffers, bouleLayers, type BouleLayers } from "./bouleTexels";

export interface BouleAsk {
  TW: number;
}

const scope = self as unknown as {
  onmessage: ((e: MessageEvent<BouleAsk>) => void) | null;
  postMessage(m: BouleLayers, transfer: ArrayBuffer[]): void;
};
scope.onmessage = (e) => {
  const l = bouleLayers(e.data.TW);
  scope.postMessage(l, bouleBuffers(l));
};
