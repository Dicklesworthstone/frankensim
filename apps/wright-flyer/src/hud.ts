// E2.4 HUD skeleton (bead wf-root-guzez.3.4): the period triad the plan
// asks for (the Wrights' own instruments: anemometer, stopwatch, engine
// rev counter) + modern overlay slots. PURE formatting; the DOM writer is
// a thin consumer.

import type { FlyerPose } from "./airframe/pose.ts";

export interface HudState {
  /** Airspeed at the aircraft [m/s] (anemometer). */
  airspeedMps: number;
  /** Elapsed flight time [s] (stopwatch). */
  elapsedS: number;
  /** Engine speed [rpm] (rev counter). */
  engineRpm: number;
  /** Camera preset name. */
  camera: string;
  /** Pose flags from the rig. */
  pose: FlyerPose;
}

/** Format the period triad + modern slots as display lines (tested pure). */
export function hudLines(s: HudState): string[] {
  const mph = s.airspeedMps * 2.23694;
  const lines = [
    `anemometer ${mph.toFixed(1)} mph (${s.airspeedMps.toFixed(1)} m/s)`,
    `stopwatch ${s.elapsedS.toFixed(1)} s`,
    `engine ${Math.round(s.engineRpm)} rpm`,
    `camera ${s.camera}`,
  ];
  if (s.pose.clamped) {
    lines.push("CONTROL AT STOP");
  }
  if (s.pose.schematicPreview) {
    lines.push("SCHEMATIC PREVIEW (geometry beyond ±25% of 1903)");
  }
  return lines;
}
