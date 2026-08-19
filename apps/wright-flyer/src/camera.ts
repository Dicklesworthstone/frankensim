// E2.4 camera presets (bead wf-root-guzez.3.4): the plan's five —
// chase / wingtip / daniels-tripod / onboard-prone / free (arrival).
// PURE: preset x (time, launch, aircraft position) -> pos/look.

import { arrivalCamera } from "./terrainMesh.ts";

export type CameraPreset = "free" | "chase" | "wingtip" | "daniels" | "onboard";

export interface CameraShot {
  pos: [number, number, number];
  look: [number, number, number];
}

/** Daniels' tripod stood a few metres downwind-left of the rail start —
 * the fixed historical viewpoint of the famous photograph. */
export function cameraFor(
  preset: CameraPreset,
  t: number,
  launch: [number, number, number],
  aircraft: [number, number, number],
): CameraShot {
  const look: [number, number, number] = [aircraft[0], aircraft[1] + 1.0, aircraft[2]];
  switch (preset) {
    case "free":
      return arrivalCamera(t, launch);
    case "chase":
      return { pos: [aircraft[0] - 18, aircraft[1] + 4.5, aircraft[2] + 7], look };
    case "wingtip":
      return { pos: [aircraft[0] + 0.5, aircraft[1] + 1.4, aircraft[2] - 7.5], look };
    case "daniels":
      return { pos: [launch[0] - 6, launch[1] + 1.6, launch[2] + 9], look };
    case "onboard":
      // Prone pilot's eye: just above the lower wing, left of center.
      return {
        pos: [aircraft[0] - 0.6, aircraft[1] + 1.35, aircraft[2] - 0.3],
        look: [aircraft[0] + 30, aircraft[1] + 1.0, aircraft[2]],
      };
  }
}

/** Number-key bindings (Digit1..Digit5). */
export const PRESET_KEYS: Record<string, CameraPreset> = {
  Digit1: "free",
  Digit2: "chase",
  Digit3: "wingtip",
  Digit4: "daniels",
  Digit5: "onboard",
};
