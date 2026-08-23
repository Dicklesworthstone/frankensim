// Landing-menu model (UI overhaul, root guzez). PURE: a menu selection
// maps to exactly the URL query the app already honors — the menu adds
// no third configuration path, it just types the params for the player.
// Repro: node --test test/menu.test.ts

import { flightByIndex } from "./missions/flights.ts";

export interface MenuSelection {
  readonly mode: "human" | "historical" | "fixed";
  readonly site: "kdh" | "huffman";
  readonly assist: boolean;
  /** Mission preset: the 1-based Dec-17 flight number, or undefined
   * for the default ensemble. Missions are Kill Devil Hills facts —
   * a Huffman selection silently drops it (like assist), because the
   * URL must never claim a scenario the engine would not run. */
  readonly flight?: number;
}

export const DEFAULT_SELECTION: MenuSelection = { mode: "human", site: "kdh", assist: false };

/** True when the assist toggle can actually engage for a selection:
 * the menu offers it as a human-pilot aid, and the Huffman scenario
 * factory pins assist OFF (protocol.ts huffmanScenario), so a Huffman
 * URL must never claim it. */
export function assistAvailable(sel: Pick<MenuSelection, "mode" | "site">): boolean {
  return sel.mode === "human" && sel.site === "kdh";
}

/** The query string for a selection. Fixed mode omits `mode` (the
 * app default); Kill Devil Hills omits `site`; `assist` appears only
 * where it can engage (see assistAvailable); `flight` appears only
 * for a valid mission id at KDH — the URL never claims a scenario the
 * engine would silently drop. */
export function menuQuery(sel: MenuSelection): string {
  const parts = ["sim=1"];
  if (sel.mode !== "fixed") {
    parts.push(`mode=${sel.mode}`);
  }
  if (sel.site === "huffman") {
    parts.push("site=huffman");
  }
  if (sel.assist && assistAvailable(sel)) {
    parts.push("assist=1");
  }
  if (sel.site === "kdh" && flightByIndex(sel.flight ?? 0) !== null) {
    parts.push(`flight=${sel.flight}`);
  }
  return `?${parts.join("&")}`;
}

/** Mission selector chips (labels from missions/flights.ts data). */
export const FLIGHT_CHIPS: readonly { readonly id: number; readonly label: string }[] = [
  { id: 1, label: "F1 · Orville" },
  { id: 2, label: "F2 · Wilbur" },
  { id: 3, label: "F3 · Orville" },
  { id: 4, label: "F4 · Wilbur" },
];

/** Mode cards shown on the landing menu. */
export const MODE_CARDS: readonly {
  mode: MenuSelection["mode"];
  title: string;
  blurb: string;
}[] = [
  {
    mode: "human",
    title: "TAKE THE CONTROLS",
    blurb: "Fly it yourself. The 1903 Flyer is unstable in pitch — expect to work.",
  },
  {
    mode: "historical",
    title: "1903 PILOT MODEL",
    blurb: "Ride along with the calibrated historical-pilot family (member 3).",
  },
  {
    mode: "fixed",
    title: "RECORDED CONTROLS",
    blurb: "Watch the fixed-controls reconstruction of the December 17 run.",
  },
];

/** Key bindings card (mirrors input.ts, the ONLY binding authority). */
export const KEY_LINES: readonly string[] = [
  "S / ↓  pull canard (nose UP)   W / ↑  push (nose DOWN)",
  "A / ←  warp left     D / →  warp right",
  "Space  recenter    V  camera toggle    1-6  camera presets",
  "M  sound on/off      H  controls card   T  telemetry",
  "R  replay with ghost N  fresh relaunch  P  photo mode",
  "Drag the view = hip cradle   Gamepad left stick = warp/pull",
  "J  guided journey (watch → assist → authentic)",
];
