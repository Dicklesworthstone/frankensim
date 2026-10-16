// E2.4 camera presets (bead wf-root-guzez.3.4): the plan's five —
// chase / wingtip / daniels-tripod / onboard-prone / free (arrival).
// PURE: preset x (time, launch, aircraft position) -> pos/look.

import { arrivalCamera } from "./terrainMesh.ts";

export type CameraPreset = "free" | "chase" | "wingtip" | "daniels" | "onboard" | "binoculars";

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
  /** Optional live attitude for the onboard seat: pitch couples the
   * eye and the look direction to the REAL airframe state, gust01 adds
   * deterministic buffet, and the head angles let the pilot look
   * around without moving the machine. Presentation-only. */
  attitude?: {
    pitchRad?: number;
    gust01?: number;
    headYawRad?: number;
    headPitchRad?: number;
  },
): CameraShot {
  const look: [number, number, number] = [aircraft[0], aircraft[1] + 1.0, aircraft[2]];
  switch (preset) {
    case "free":
      return arrivalCamera(t, launch);
    case "chase":
      // Directly BEHIND the machine, close enough to read Wilbur on
      // the cradle and see the wing warping and rudder action.
      return { pos: [aircraft[0] - 11.5, aircraft[1] + 3.2, aircraft[2] + 0.001], look };
    case "wingtip":
      // Starboard wingtip camera looking inward past the spinning props across the wingspan
      return { pos: [aircraft[0] + 0.5, aircraft[1] + 1.4, aircraft[2] - 7.5], look };
    case "daniels":
      return { pos: [launch[0] - 6, launch[1] + 1.6, launch[2] + 9], look };
    case "onboard": {
      // Prone pilot's eye: just above the lower wing, left of center,
      // looking through the canard. The eye OFFSET rides the airframe:
      // pitch theta about world z carries forward (cos, sin) and up
      // (-sin, cos), so climbing tilts the eye back and the horizon
      // dips — the cockpit lives in the machine, not above the map.
      const th = attitude?.pitchRad ?? 0;
      const fwd: [number, number] = [Math.cos(th), Math.sin(th)];
      const up: [number, number] = [-Math.sin(th), Math.cos(th)];
      const gust = attitude?.gust01 ?? 0;
      // Deterministic buffet: two incommensurate sines per axis, the
      // vertical one stronger — the seat of the pants reads as pitch.
      const bx = (Math.sin(t * 8.7) * 0.014 + Math.sin(t * 13.1 + 1.7) * 0.008) * gust;
      const by = (Math.sin(t * 5.3) * 0.016 + Math.sin(t * 11.3 + 0.6) * 0.01) * gust;
      const eye: [number, number, number] = [
        aircraft[0] - 0.6 * fwd[0] + 1.35 * up[0] + bx,
        aircraft[1] - 0.6 * fwd[1] + 1.35 * up[1] + by,
        aircraft[2] - 0.3,
      ];
      // Look direction: head pitch adds to airframe pitch, then head
      // yaw swings about WORLD up (a prone pilot's head yaw).
      const thL = th + (attitude?.headPitchRad ?? 0);
      const yaw = attitude?.headYawRad ?? 0;
      const dir: [number, number, number] = [
        Math.cos(thL) * Math.cos(yaw),
        Math.sin(thL),
        -Math.cos(thL) * Math.sin(yaw),
      ];
      return {
        pos: eye,
        look: [eye[0] + dir[0] * 30, eye[1] + dir[1] * 30, eye[2] + dir[2] * 30],
      };
    }
    case "binoculars":
      // Orville's field glasses tracking the flight from the ground.
      return {
        pos: [launch[0] + 12, launch[1] + 1.75, launch[2] + 7.4],
        look: [aircraft[0], aircraft[1] + 0.8, aircraft[2]],
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
 * glide feels identical at 30 and 144 fps. Pure in its inputs.
 * Presentation smoothing must NEVER crash the frame: non-finite or
 * huge dt (tab suspend → rAF gap) CLAMPS — a giant dt snaps to the
 * target, which is exactly the right recovery. */
export function easeCameraToward(
  prev: CameraState,
  target: CameraShot,
  dtS: number,
  k = 5.5,
): CameraState {
  const dt = Number.isFinite(dtS) ? Math.min(Math.max(dtS, 0), 1) : 1;
  const t = 1 - Math.exp(-k * dt);
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

/** Base vertical FOV for the game presets. */
export const BASE_FOV_DEG = 50;

/** Base FOV widened by airspeed — a subtle rush cue, clamped to
 * +8.4° at racing speed so presets keep their framing. */
export function speedFov(baseDeg: number, airspeedMps: number): number {
  const extra = Math.max(0, Math.min(24, airspeedMps - 9)) * 0.35;
  return baseDeg + extra;
}

/** Number-key bindings (Digit1..Digit6). */
export const PRESET_KEYS: Record<string, CameraPreset> = {
  Digit1: "free",
  Digit2: "chase",
  Digit3: "wingtip",
  Digit4: "daniels",
  Digit5: "onboard",
  Digit6: "binoculars",
};
