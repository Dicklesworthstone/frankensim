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
      // Directly BEHIND the machine, close enough to read Wilbur on
      // the cradle (the old 18-m back / 7-m side framing read as a
      // detached observer, not a piloting view).
      return { pos: [aircraft[0] - 11, aircraft[1] + 3.2, aircraft[2] + 0.001], look };
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

/** Presentation-plane camera state: smoothed position, look point,
 * and vertical FOV. The scene holds ONE instance and advances it with
 * easeCameraToward every frame. */
export interface CameraState {
  pos: [number, number, number];
  look: [number, number, number];
  fovDeg: number;
}

/** Exponential critically-damped approach toward a target shot:
 * blend factor `1 - exp(-k·dt)` is frame-rate independent, so the
 * glide feels identical at 30 and 144 fps. Pure in its inputs. */
export function easeCameraToward(
  prev: CameraState,
  target: CameraShot,
  dtS: number,
  k = 5.5,
): CameraState {
  if (!Number.isFinite(dtS) || dtS < 0 || dtS > 1) {
    throw new RangeError(`camera dt out of domain: ${dtS}`);
  }
  const t = 1 - Math.exp(-k * dtS);
  const mix = (
    a: readonly [number, number, number],
    b: readonly [number, number, number],
  ): [number, number, number] => [
    a[0] + (b[0] - a[0]) * t,
    a[1] + (b[1] - a[1]) * t,
    a[2] + (b[2] - a[2]) * t,
  ];
  return { pos: mix(prev.pos, target.pos), look: mix(prev.look, target.look), fovDeg: prev.fovDeg };
}

/** Base FOV widened by airspeed — a subtle rush cue, clamped to
 * +8.4° at racing speed so presets keep their framing. */
export function speedFov(baseDeg: number, airspeedMps: number): number {
  const extra = Math.max(0, Math.min(24, airspeedMps - 9)) * 0.35;
  return baseDeg + extra;
}

/** Number-key bindings (Digit1..Digit5). */
export const PRESET_KEYS: Record<string, CameraPreset> = {
  Digit1: "free",
  Digit2: "chase",
  Digit3: "wingtip",
  Digit4: "daniels",
  Digit5: "onboard",
};
