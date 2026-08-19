// E2.2 pose applier (bead wf-root-guzez.3.2): the THIN mutator from a
// computed FlyerPose onto the parametric airframe's articulated groups.
// All rig logic lives in pose.ts (tested headlessly); this file only
// assigns transforms — it must stay too simple to hide a bug.

import type { FlyerAirframe } from "./parametricAirframe.ts";
import type { FlyerPose } from "./pose.ts";
import type { ControlState } from "./pose.ts";
import { computePose } from "./pose.ts";

export function applyPose(a: FlyerAirframe, p: FlyerPose): void {
  // Canard pitches about its hinge axis (+y in the airframe's frame).
  a.canardGroup.rotation.x = p.canardRad;
  // Twin rudders yaw together.
  a.rudderGroup.rotation.y = p.rudderRad;
  // Warp: tip-section twist (the airframe's 0.6 mode-shape factor is
  // applied inside its own update helper; here we drive the groups the
  // interface exposes — upper/lower wing tip rotations mirror L/R).
  a.upperWing.rotation.z = p.warpTipRad * 0.05;
  a.lowerWing.rotation.z = p.warpTipRad * 0.05;
  // Hip cradle slides with the warp command.
  a.cradleGroup.position.x = p.cradleOffsetM;
  // Counter-rotating propellers.
  a.leftPropBlades.rotation.z = p.leftPropRad;
  a.rightPropBlades.rotation.z = p.rightPropRad;
}

/** The scripted demo state (DONE-WHEN driver): a deterministic 12-second
 * control script — canard doublet, warp-with-slaved-rudder reversal, prop
 * spin-up — pure in t so the 60 fps loop and tests share it. */
export function scriptedState(tS: number): ControlState {
  const t = tS % 12;
  const canardDeg = t < 3 ? 12 * Math.sin((Math.PI * t) / 1.5) : t < 4 ? 0 : 0;
  const warpDeg = t >= 4 && t < 9 ? 8.5 * Math.sin((2 * Math.PI * (t - 4)) / 5) : 0;
  const rpm = Math.min(350, 60 * t);
  return {
    canardDeg,
    warpDeg,
    rudderDeg: 0,
    coupled: true,
    propAngleRad: ((rpm / 60) * 2 * Math.PI * t) % (2 * Math.PI),
  };
}

/** One frame of the scripted demo: compute + apply; returns the pose so
 * the HUD can show clamp/preview flags. */
export function driveScripted(a: FlyerAirframe, tS: number): FlyerPose {
  const pose = computePose(scriptedState(tS));
  applyPose(a, pose);
  return pose;
}
