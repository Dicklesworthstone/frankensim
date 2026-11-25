/** The pan foods' layers, the knife's, then the label's paper (foodTexels.ts), computed off the page's main thread and handed back without a copy. */
import { knifeBuffers, knifeTexels, labelPaper, labelPaperBuffers, panFoodBuffers, panFoodTexels, type KnifeTexels, type LabelPaper, type PanFoodTexels } from "./foodTexels";

export type FoodTexelsReply = { kind: "pan"; t: PanFoodTexels } | { kind: "knife"; t: KnifeTexels } | { kind: "label"; t: LabelPaper };

const scope = self as unknown as {
  onmessage: ((e: MessageEvent<null>) => void) | null;
  postMessage(m: FoodTexelsReply, transfer: ArrayBuffer[]): void;
};
scope.onmessage = () => {
  // The hero kitchen's pans need theirs first, then its cooks' knives; the label is the stove chapter's.
  const pan = panFoodTexels();
  scope.postMessage({ kind: "pan", t: pan }, panFoodBuffers(pan));
  const knife = knifeTexels();
  scope.postMessage({ kind: "knife", t: knife }, knifeBuffers(knife));
  const label = labelPaper();
  scope.postMessage({ kind: "label", t: label }, labelPaperBuffers(label));
};
