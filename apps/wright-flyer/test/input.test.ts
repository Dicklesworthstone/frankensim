// E2.4 input/HUD battery (bead wf-root-guzez.3.4): deterministic slew,
// quantization-grid membership of EVERY emitted command, clamp at full
// travel, recenter decay, key-binding map, HUD triad + flag lines.
// Repro: node --test test/input.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";
import {
  CRADLE_FULL_TRAVEL_PX,
  GAMEPAD_DEADZONE,
  NEUTRAL,
  cradleFromPointer,
  decayCradle,
  keysFrom,
  sampleGamepad,
  stepCommand,
} from "../src/input.ts";
import { hudLines } from "../src/hud.ts";
import { computePose } from "../src/airframe/pose.ts";

const K = {
  canardUp: false, canardDown: false, warpLeft: false, warpRight: false, recenter: false,
};

test("slew is deterministic and every command sits on the 1/4096 grid", () => {
  let cmd = NEUTRAL;
  const dt = 1 / 120;
  for (let i = 0; i < 240; i++) {
    cmd = stepCommand(cmd, { ...K, canardUp: true, warpRight: true }, dt);
    assert.ok(Math.abs(cmd.canard * 4096 - Math.round(cmd.canard * 4096)) < 1e-9, "grid");
    assert.ok(Math.abs(cmd.warp * 4096 - Math.round(cmd.warp * 4096)) < 1e-9, "grid");
  }
  // Two seconds at 1.4/s -> clamped at full travel; warp at 1.1/s -> 1.0 too.
  assert.equal(cmd.canard, 1);
  assert.equal(cmd.warp, 1);
  // Determinism: same sequence, same bits.
  let again = NEUTRAL;
  for (let i = 0; i < 240; i++) {
    again = stepCommand(again, { ...K, canardUp: true, warpRight: true }, dt);
  }
  assert.deepEqual(cmd, again);
});

test("opposed keys hold; recenter decays toward neutral", () => {
  let cmd = { ...NEUTRAL, canard: 0.5, warp: -0.5 };
  const held = stepCommand(cmd, { ...K, canardUp: true, canardDown: true }, 1 / 120);
  assert.equal(held.canard, 0.5, "opposed keys cancel");
  for (let i = 0; i < 240; i++) {
    cmd = stepCommand(cmd, { ...K, recenter: true }, 1 / 120);
  }
  assert.ok(Math.abs(cmd.canard) < 0.01 && Math.abs(cmd.warp) < 0.01, "recenter decays");
  assert.throws(() => stepCommand(NEUTRAL, K, Number.NaN), /finite/);
});

test("key bindings map both arrow and WASD, pull = nose up", () => {
  const k = keysFrom(new Set(["ArrowDown", "KeyD", "Space"]));
  assert.ok(k.canardUp && k.warpRight && k.recenter && !k.canardDown && !k.warpLeft);
  const wasd = keysFrom(new Set(["KeyS", "KeyA"]));
  assert.ok(wasd.canardUp && wasd.warpLeft);
});

test("HUD renders the period triad and surfaces the rig flags", () => {
  const pose = computePose({ canardDeg: 31, warpDeg: 0, rudderDeg: 0, coupled: true, propAngleRad: 0 });
  const lines = hudLines({ airspeedMps: 13.86, elapsedS: 12, engineRpm: 1025, camera: "chase", pose });
  assert.match(lines[0]!, /31\.0 mph/);
  assert.match(lines[1]!, /12\.0 s/);
  assert.match(lines[2]!, /1025 rpm/);
  assert.ok(lines.includes("CONTROL AT STOP"), "the clamp flag must surface");
  const clean = hudLines({
    airspeedMps: 0, elapsedS: 0, engineRpm: 0, camera: "arrival",
    pose: computePose({ canardDeg: 0, warpDeg: 0, rudderDeg: 0, coupled: true, propAngleRad: 0 }),
  });
  assert.equal(clean.length, 4, "no flags when clean");
});

test("all five camera presets work against scripted state", async () => {
  const { cameraFor, PRESET_KEYS } = await import("../src/camera.ts");
  const launch: [number, number, number] = [0, 4.1, -625];
  const aircraft: [number, number, number] = [12, 5.3, -625];
  const seen = new Set<string>();
  for (const preset of Object.values(PRESET_KEYS)) {
    const shot = cameraFor(preset, 6.5, launch, aircraft);
    for (const v of [...shot.pos, ...shot.look]) assert.ok(Number.isFinite(v), preset);
    assert.ok(shot.pos[1] > 0, `${preset} camera above ground`);
    seen.add(shot.pos.map((v) => v.toFixed(2)).join(","));
    if (preset !== "daniels" && preset !== "free") {
      assert.ok(Math.hypot(shot.pos[0] - aircraft[0], shot.pos[2] - aircraft[2]) < 25, preset);
    }
  }
  assert.equal(seen.size, 5, "five DISTINCT viewpoints");
  const d1 = cameraFor("daniels", 1, launch, aircraft);
  const d2 = cameraFor("daniels", 9, launch, [50, 8, -600]);
  assert.deepEqual(d1.pos, d2.pos, "the tripod never moves");
});

test("hip cradle maps drag offset to a quantized position command", () => {
  // Center grab = neutral.
  assert.deepEqual(cradleFromPointer(0, 0), { canard: 0, warp: 0, mode: "mouse-cradle" });
  // Drag DOWN = pull (+canard), drag RIGHT = +warp (sign conventions).
  const pull = cradleFromPointer(0, CRADLE_FULL_TRAVEL_PX);
  assert.equal(pull.canard, 1);
  const push = cradleFromPointer(0, -CRADLE_FULL_TRAVEL_PX);
  assert.equal(push.canard, -1);
  const right = cradleFromPointer(CRADLE_FULL_TRAVEL_PX, 0);
  assert.equal(right.warp, 1);
  // Overshoot clamps at full travel; diagonal lands on the grid.
  const over = cradleFromPointer(4000, -4000);
  assert.equal(over.warp, 1);
  assert.equal(over.canard, -1);
  const diag = cradleFromPointer(45, 22.5);
  for (const v of [diag.canard, diag.warp]) {
    assert.ok(Math.abs(v * 4096 - Math.round(v * 4096)) < 1e-9, "grid");
  }
  // Non-finite offsets refuse loudly.
  assert.throws(() => cradleFromPointer(Number.NaN, 0), /finite/);
});

test("cradle release decays on the recenter spring and rests bit-neutral", () => {
  let cmd = cradleFromPointer(CRADLE_FULL_TRAVEL_PX, CRADLE_FULL_TRAVEL_PX);
  assert.throws(() => decayCradle(cmd, Number.NaN), /finite/);
  assert.throws(() => decayCradle(cmd, -1), /finite/);
  for (let i = 0; i < 600; i++) {
    cmd = decayCradle(cmd, 1 / 120);
    for (const v of [cmd.canard, cmd.warp]) {
      assert.ok(Math.abs(v) <= 1, "bounded");
    }
  }
  assert.equal(cmd.canard, 0, "exact neutral rest");
  assert.equal(cmd.warp, 0, "exact neutral rest");
  assert.equal(cmd.mode, "mouse-cradle", "mode label survives decay");
  // Determinism: same replay of steps -> same bits.
  let again = cradleFromPointer(CRADLE_FULL_TRAVEL_PX, CRADLE_FULL_TRAVEL_PX);
  for (let i = 0; i < 600; i++) again = decayCradle(again, 1 / 120);
  assert.deepEqual(cmd, again);
});

test("gamepad radial deadzone kills drift but preserves direction", () => {
  assert.equal(sampleGamepad(null), null, "no pad");
  assert.equal(sampleGamepad({ connected: false, axes: [0.5, 0.5] }), null, "disconnected");
  assert.equal(sampleGamepad({ connected: true, axes: [] }), null, "too few axes");
  const rest = sampleGamepad({ connected: true, axes: [0.05, -0.08] });
  assert.ok(rest !== null && rest.canard === 0 && rest.warp === 0, "inside deadzone = neutral");
  // Pure diagonal just past the zone keeps its direction (no snap).
  const d = sampleGamepad({
    connected: true,
    axes: [GAMEPAD_DEADZONE + 0.01, GAMEPAD_DEADZONE + 0.01],
  });
  assert.ok(d !== null && d.warp > 0 && Math.abs(d.warp - d.canard) < 1e-9, "diagonal stays diagonal");
  // Full deflection is exactly full travel on the grid; stick back pulls.
  const full = sampleGamepad({ connected: true, axes: [0, 1] });
  assert.ok(full !== null && full.canard === 1 && full.warp === 0);
  const grid = sampleGamepad({ connected: true, axes: [0.7, -0.2] });
  assert.ok(grid !== null);
  for (const v of [grid.canard, grid.warp]) {
    assert.ok(Math.abs(v * 4096 - Math.round(v * 4096)) < 1e-9, "grid");
  }
});
