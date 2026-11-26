/**
 * One purchase for the morning's kitchen (chapter 07, fourth figure): a third combi, a second six-burner range, or a
 * fourth cook, and any two or all three. Each is the morning's own planner (planMorning) run with that equipment and
 * crew, so what a purchase buys is what the plan does with it: the most meals that make the 10:30 van with every
 * portion at its best (morningLimit), and which station is full while the cooks stand free with no room at it.
 */
import { KIT, type Crew, type MorningSettings } from "./cateringMorning";

/** What can be bought: a third combi, a second range, a fourth cook. */
export interface Buy {
  combi: boolean;
  range: boolean;
  cook: boolean;
}
export interface KitSettings extends Buy {
  meals: number;
  menu: number;
}
/** Every choice of purchases, nothing first, then one at a time, then two, then all three. */
export const BUYS: readonly Buy[] = [
  { combi: false, range: false, cook: false },
  { combi: true, range: false, cook: false },
  { combi: false, range: true, cook: false },
  { combi: false, range: false, cook: true },
  { combi: true, range: true, cook: false },
  { combi: true, range: false, cook: true },
  { combi: false, range: true, cook: true },
  { combi: true, range: true, cook: true },
];
export const buyKey = (b: Buy): string => `${+b.combi}${+b.range}${+b.cook}`;
export const planKey = (s: KitSettings): string => `${s.meals}|${s.menu}|${buyKey(s)}`;
export const limitKey = (menu: number, b: Buy): string => `${menu}|${buyKey(b)}`;
/** The planner's settings for an order and what was bought: the hub's three cooks and its kit, plus the purchases. */
export const morningOf = (meals: number, menu: number, b: Buy): MorningSettings => ({
  meals,
  menu,
  cooks: (b.cook ? 4 : 3) as Crew,
  ahead: false,
  kit: { ovens: KIT.ovens + (b.combi ? 1 : 0), ranges: KIT.ranges + (b.range ? 1 : 0) },
});
