/** The pan foods' layers, then the label's paper (foodTexels.ts), computed off the page's main thread and handed back without a copy. */
import { labelPaper, labelPaperBuffers, panFoodBuffers, panFoodTexels, type LabelPaper, type PanFoodTexels } from "./foodTexels";

export type FoodTexelsReply = { kind: "pan"; t: PanFoodTexels } | { kind: "label"; t: LabelPaper };

const scope = self as unknown as {
  onmessage: ((e: MessageEvent<null>) => void) | null;
  postMessage(m: FoodTexelsReply, transfer: ArrayBuffer[]): void;
};
scope.onmessage = () => {
  // The hero kitchen's pans need theirs first; the label is the stove chapter's.
  const pan = panFoodTexels();
  scope.postMessage({ kind: "pan", t: pan }, panFoodBuffers(pan));
  const label = labelPaper();
  scope.postMessage({ kind: "label", t: label }, labelPaperBuffers(label));
};
