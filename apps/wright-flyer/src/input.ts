// E2.4 tick-quantized input mapping (bead wf-root-guzez.3.4). PURE key
// state -> quantized control commands through the transport quantizer
// (1/4096 grid): the same deterministic values the InputScheduler admits,
// so replays reproduce what the pilot pressed. Transducer families:
//   - keyboard-rate: keys SLEW the command over dt (devices cannot claim
//     to measure force — plan law);
//   - mouse-cradle: the plan §2.3 hip cradle — pointer DRAG OFFSET is a
//     position command (drag down = pull, right = warp right); springs
//     back to neutral on release through decayCradle;
//   - gamepad-stick: left-stick deflection is a position command with a
//     radial deadzone; stick back (y+) pulls, like the keyboard's ↓/S.
// Every family emits onto the SAME 1/4096 grid, so the applied-input
// trace cannot tell devices apart except through the mode label.

import { quantizeControl } from "./transport/inputClock.ts";

export type InputTransducerMode = "keyboard-rate" | "mouse-cradle" | "gamepad-stick";

export interface KeyState {
  canardUp: boolean;
  canardDown: boolean;
  warpLeft: boolean;
  warpRight: boolean;
  recenter: boolean;
}

export interface PilotCommand {
  /** Quantized canard command in [-1, 1] (maps to ±30° at the pose). */
  canard: number;
  /** Quantized warp command in [-1, 1] (maps to ±8.5°). */
  warp: number;
  /** The transducer that produced this (identity ingredient). */
  mode: InputTransducerMode;
}

/** Full-scale slew per second of held key (fraction of full travel). */
export const CANARD_SLEW_PER_S = 1.4;
/** Warp slew per second of held key. */
export const WARP_SLEW_PER_S = 1.1;
/** Recenter spring rate [1/s] toward neutral while `recenter` is held. */
export const RECENTER_PER_S = 3.0;

/** One input tick: slew the previous command by the held keys over dt,
 * clamp to [-1, 1], and QUANTIZE — the returned values are exactly what
 * the input trace records. Deterministic and pure. */
export function stepCommand(prev: PilotCommand, keys: KeyState, dtS: number): PilotCommand {
  if (!Number.isFinite(dtS) || dtS < 0) {
    throw new RangeError(`dt must be finite and non-negative, got ${dtS}`);
  }
  const axis = (value: number, up: boolean, down: boolean, ratePerS: number): number => {
    let v = value;
    if (keys.recenter) {
      const decay = Math.exp(-RECENTER_PER_S * dtS);
      v *= decay;
    }
    if (up !== down) {
      v += (up ? 1 : -1) * ratePerS * dtS;
    }
    return quantizeControl(Math.max(-1, Math.min(1, v)));
  };
  return {
    canard: axis(prev.canard, keys.canardUp, keys.canardDown, CANARD_SLEW_PER_S),
    warp: axis(prev.warp, keys.warpRight, keys.warpLeft, WARP_SLEW_PER_S),
    mode: "keyboard-rate",
  };
}

/** The neutral command. */
export const NEUTRAL: PilotCommand = { canard: 0, warp: 0, mode: "keyboard-rate" };

/** Map key codes to the KeyState (the ONLY place bindings live). */
export function keysFrom(down: ReadonlySet<string>): KeyState {
  return {
    canardUp: down.has("ArrowDown") || down.has("KeyS"), // pull = nose up
    canardDown: down.has("ArrowUp") || down.has("KeyW"),
    warpLeft: down.has("ArrowLeft") || down.has("KeyA"),
    warpRight: down.has("ArrowRight") || down.has("KeyD"),
    recenter: down.has("Space"),
  };
}

/* ------------------- position transducers (plan §2.3) --------------- */

/** Pointer drag radius [px] mapped to FULL travel in each axis. */
export const CRADLE_FULL_TRAVEL_PX = 90;

/** The hip cradle: pointer offset from the grab point [px] -> command.
 * Screen y grows DOWNWARD, and dragging toward the pilot (down) is the
 * pull direction, so +dy maps to +canard. Pure; quantized. */
export function cradleFromPointer(dxPx: number, dyPx: number): PilotCommand {
  if (!Number.isFinite(dxPx) || !Number.isFinite(dyPx)) {
    throw new RangeError(`cradle offsets must be finite, got ${dxPx}, ${dyPx}`);
  }
  return {
    canard: quantizeControl(Math.max(-1, Math.min(1, dyPx / CRADLE_FULL_TRAVEL_PX))),
    warp: quantizeControl(Math.max(-1, Math.min(1, dxPx / CRADLE_FULL_TRAVEL_PX))),
    mode: "mouse-cradle",
  };
}
/** Rest deadband [fraction of full travel]. Quantized exponential decay
 * has sticky fixed points: a quantized value k/4096 stops moving once
 * k·(1 − e^(−rate·dt)) rounds back to k (here every k ≤ 20 sticks).
 * Below this band the spring therefore snaps to EXACT neutral — the
 * resting command becomes bit-neutral regardless of release history,
 * which replay identity cares about. 0.005 travel ≈ 0.15° canard. */
export const CRADLE_REST_BAND = 0.005;

/** Release spring: exponential decay toward neutral at the SAME rate as
 * the keyboard recenter spring, then the declared rest-band snap. */
export function decayCradle(cmd: PilotCommand, dtS: number): PilotCommand {
  if (!Number.isFinite(dtS) || dtS < 0) {
    throw new RangeError(`dt must be finite and non-negative, got ${dtS}`);
  }
  const decay = Math.exp(-RECENTER_PER_S * dtS);
  const rawCanard = cmd.canard * decay;
  const rawWarp = cmd.warp * decay;
  if (Math.abs(rawCanard) <= CRADLE_REST_BAND && Math.abs(rawWarp) <= CRADLE_REST_BAND) {
    return { canard: 0, warp: 0, mode: "mouse-cradle" };
  }
  return {
    canard: quantizeControl(rawCanard),
    warp: quantizeControl(rawWarp),
    mode: "mouse-cradle",
  };
}

/** Structural gamepad view — headless-testable; main.ts adapts the DOM
 * Gamepad object into this shape (never imported here). */
export interface GamepadLike {
  readonly connected: boolean;
  readonly axes: readonly number[];
}

/** Radial deadzone: deflections below this are exactly neutral. */
export const GAMEPAD_DEADZONE = 0.12;

/** Left-stick sample -> command. Radial deadzone with rescaling keeps
 * the stick's DIRECTION across the zone edge (no snap-to-axis) while
 * guaranteeing the resting stick reads exactly neutral. Stick back
 * (screen convention y+ = pulled) is pull/+canard, matching ↓/S. */
export function sampleGamepad(pad: GamepadLike | null | undefined): PilotCommand | null {
  if (pad === null || pad === undefined || !pad.connected || pad.axes.length < 2) {
    return null;
  }
  const x = Number.isFinite(pad.axes[0]) ? pad.axes[0]! : 0;
  const y = Number.isFinite(pad.axes[1]) ? pad.axes[1]! : 0;
  const m = Math.hypot(x, y);
  if (m <= GAMEPAD_DEADZONE) {
    return { canard: 0, warp: 0, mode: "gamepad-stick" };
  }
  const f = Math.min(1, (m - GAMEPAD_DEADZONE) / (1 - GAMEPAD_DEADZONE) / m);
  return {
    canard: quantizeControl(y * f),
    warp: quantizeControl(x * f),
    mode: "gamepad-stick",
  };
}
