// E5.1 node harness (bead wf-root-guzez.6.2): drives the REAL wasm
// binary (wasm32-unknown-unknown via wasm-bindgen nodejs target)
// through init/step/digest and EXECUTES every documented refusal code
// at the engine surface. JSONL receipts on stdout; exit 1 on any miss.
//
// Build + run:
//   cd crates/fs-flyer-wasm
//   wasm-pack build --target nodejs --release --out-dir <PKG_DIR>
//   node node-harness/engine_harness.mjs <PKG_DIR>/fs_flyer_wasm.js

import { createRequire } from "node:module";
const require = createRequire(import.meta.url);

const pkgPath = process.argv[2];
if (!pkgPath) {
  console.error("usage: node engine_harness.mjs <path-to-fs_flyer_wasm.js>");
  process.exit(2);
}
const wasm = require(pkgPath);

let failures = 0;
function jlog(caseName, payload) {
  console.log(JSON.stringify({ suite: "fs-flyer-wasm-node-e51", case: caseName, ...payload }));
}
function check(name, cond, detail) {
  if (!cond) {
    failures += 1;
    jlog(name, { pass: false, detail });
  } else {
    jlog(name, { pass: true });
  }
}

// --- refusal codes before any init -----------------------------------------
const stepBefore = JSON.parse(wasm.flyer_engine_step(false, 0.0, 0.0));
check("engine-not-initialized:step", stepBefore.refusal?.code === "engine-not-initialized", stepBefore);
const digestBefore = JSON.parse(wasm.flyer_engine_digest());
check("engine-not-initialized:digest", digestBefore.refusal?.code === "engine-not-initialized", digestBefore);

// --- mode-invalid + scenario-invalid (cap AND cap+1) ------------------------
const badMode = JSON.parse(wasm.flyer_engine_init(1n, 1.294, 11.0, 3, 0, 18.3, 40n, false, false));
check("mode-invalid", badMode.refusal?.code === "mode-invalid", badMode);
const capOk = JSON.parse(wasm.flyer_engine_init(1n, 1.294, 20.0, 0, 0, 18.3, 2n, false, false));
check("headwind-at-cap-admits", !!capOk.ok, capOk);
const capPlus = JSON.parse(wasm.flyer_engine_init(1n, 1.294, 20.000000000000004, 0, 0, 18.3, 40n, false, false));
check("scenario-invalid", capPlus.refusal?.code === "scenario-invalid", capPlus);

// --- real Dec-17 init: identity fields --------------------------------------
const init = JSON.parse(wasm.flyer_engine_init(1903n, 1.294, 11.0, 0, 0, 18.3, 40n, false, false));
check(
  "init-identity",
  !!init.ok && typeof init.ok.run_intent_id === "string" && typeof init.ok.tick0_digest === "string" && init.ok.trim_v_mps > 5.0,
  init
);

// --- stepping: on-rail, in-band terminal, run-ended -------------------------
let last = null;
for (let i = 0; i < 40; i++) last = JSON.parse(wasm.flyer_engine_step(false, 0.0, 0.0));
check("terminal-in-band", last.ok?.phase === "ended:max-ticks", last);
const past = JSON.parse(wasm.flyer_engine_step(false, 0.0, 0.0));
check("run-ended", past.refusal?.code === "run-ended", past);
const digest = JSON.parse(wasm.flyer_engine_digest());
check("digest", typeof digest.ok?.digest === "string" && digest.ok.digest.length === 64, digest);

// --- human mode: control-input-missing (absent AND non-finite) -------------
const humanInit = JSON.parse(wasm.flyer_engine_init(7n, 1.294, 11.0, 2, 0, 18.3, 40n, false, false));
check("human-init", !!humanInit.ok, humanInit);
const noInput = JSON.parse(wasm.flyer_engine_step(false, 0.0, 0.0));
check("control-input-missing:absent", noInput.refusal?.code === "control-input-missing", noInput);
const nanInput = JSON.parse(wasm.flyer_engine_step(true, Number.NaN, 0.0));
check("control-input-missing:nan", nanInput.refusal?.code === "control-input-missing", nanInput);
const withInput = JSON.parse(wasm.flyer_engine_step(true, 25.0, 0.01));
check("human-steps", withInput.ok?.phase === "on-rail", withInput);

// --- assist visibility (E5.3c) ----------------------------------------------
const assistInit = JSON.parse(wasm.flyer_engine_init(1903n, 1.294, 11.0, 0, 0, 18.3, 5n, true, false));
check("assist-init", !!assistInit.ok, assistInit);
const assistStep = JSON.parse(wasm.flyer_engine_step(false, 0.0, 0.0));
check("assist-visible", assistStep.ok?.assist_active === true && typeof assistStep.ok?.assist_dc_rad === "number", assistStep);

// --- catapult intent binding (E5.4) -----------------------------------------
const catInit = JSON.parse(wasm.flyer_engine_init(1903n, 1.294, 11.0, 0, 0, 18.3, 5n, false, true));
check("catapult-init", !!catInit.ok, catInit);
const plainInit = JSON.parse(wasm.flyer_engine_init(1903n, 1.294, 11.0, 0, 0, 18.3, 5n, false, false));
check("catapult-binds-intent", catInit.ok?.run_intent_id !== plainInit.ok?.run_intent_id && catInit.ok?.tick0_digest === plainInit.ok?.tick0_digest, {});

// --- determinism across re-init ---------------------------------------------
function transcript() {
  const parts = [wasm.flyer_engine_init(1903n, 1.294, 11.0, 0, 0, 18.3, 12n, false, false)];
  for (let i = 0; i < 12; i++) parts.push(wasm.flyer_engine_step(false, 0.0, 0.0));
  parts.push(wasm.flyer_engine_digest());
  return parts.join("\n");
}
check("determinism", transcript() === transcript());

jlog("summary", { failures });
process.exit(failures === 0 ? 0 : 1);
