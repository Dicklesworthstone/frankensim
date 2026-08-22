// Landing-menu model (UI overhaul, root guzez). PURE: a menu selection
// maps to exactly the URL query the app already honors — the menu adds
// no third configuration path, it just types the params for the player.
// Repro: node --test test/menu.test.ts

export interface MenuSelection {
  readonly mode: "human" | "historical" | "fixed";
  readonly site: "kdh" | "huffman";
  readonly assist: boolean;
}

export const DEFAULT_SELECTION: MenuSelection = { mode: "human", site: "kdh", assist: false };

/** The query string for a selection. Fixed mode omits `mode` (the
 * app default); Kill Devil Hills omits `site`; assist appears only in
 * human mode (it is a human-pilot aid — other modes ignore it, so the
 * URL never claims an aid that cannot engage). */
export function menuQuery(sel: MenuSelection): string {
  const parts = ["sim=1"];
  if (sel.mode !== "fixed") {
    parts.push(`mode=${sel.mode}`);
  }
  if (sel.site === "huffman") {
    parts.push("site=huffman");
  }
  if (sel.assist && sel.mode === "human") {
    parts.push("assist=1");
  }
  return `?${parts.join("&")}`;
}

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
  "S / ↓  pull canard (nose up)",
  "W / ↑  push canard (nose down)",
  "A / ←  warp left     D / →  warp right",
  "Space  recenter      R  replay with ghost",
  "T  telemetry panel",
];
