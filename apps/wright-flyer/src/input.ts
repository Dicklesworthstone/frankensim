// E2.4 tick-quantized input mapping (bead wf-root-guzez.3.4). PURE key
// state -> quantized control commands through the transport quantizer
// (1/4096 grid): the same deterministic values the InputScheduler admits,
// so replays reproduce what the pilot pressed. InputTransducerMode:
// keyboard is a RATE transducer (keys slew the command; devices cannot
// claim to measure force — plan law).

import { quantizeControl } from "./transport/inputClock.ts";

export type InputTransducerMode = "keyboard-rate";

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
