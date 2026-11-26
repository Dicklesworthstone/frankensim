/**
 * Where each menu's limit is, off the page's main thread: for a crew, every menu (cooked this morning, and with what
 * keeps cooked the day before), the most meals that make the 10:30 van with every portion at its best at lunch
 * (morningLimit). Posted back one at a time, as they come.
 */
import { MENUS, morningLimit, type Crew } from "./cateringMorning";

export interface LimitAsk {
  key: string;
  cooks: Crew;
}
export interface LimitScore {
  key: string;
  menu: number;
  ahead: boolean;
  meals: number;
}

const scope = self as unknown as { onmessage: ((e: MessageEvent<LimitAsk>) => void) | null; postMessage(m: LimitScore): void };
scope.onmessage = (e) => {
  const { key, cooks } = e.data;
  for (const ahead of [false, true])
    for (let menu = 0; menu < MENUS; menu++) scope.postMessage({ key, menu, ahead, meals: morningLimit({ cooks, menu, ahead }) });
};
