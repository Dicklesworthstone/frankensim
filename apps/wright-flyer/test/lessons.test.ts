// E7.5 battery (bead wf-root-guzez.8.7): the three curated lessons
// run END-TO-END (start -> advance through every step -> done, each
// step's overlays named); the undeclared-claim falsifier refuses (a
// curated script never improvises physics); the perception view is
// verbatim cue-state-beside-truth with the gap DISPLAYED and takes
// no renderer input by construction; caps at cap AND cap+1.
// Repro: node --test test/lessons.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  MAX_CUES,
  MAX_STEPS,
  advanceLesson,
  curatedLessons,
  perceptionView,
  startLesson,
  validateLesson,
  type Lesson,
} from "../src/lessons.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-lessons","case":"${kase}",${payload}}`);
}

test("the three curated lessons run end-to-end", () => {
  const lessons = curatedLessons();
  assert.equal(lessons.length, 3, "three v1 lessons");
  for (const lesson of lessons) {
    const started = startLesson(lesson);
    assert.ok(started.ok, `${lesson.id} starts`);
    if (!started.ok) continue;
    let run = started.value;
    let steps = 0;
    while (!run.done) {
      const step = run.lesson.steps[run.stepIndex];
      assert.ok(step !== undefined && step.overlays.length > 0, "every step shows something");
      const adv = advanceLesson(run);
      assert.ok(adv.ok);
      if (!adv.ok) break;
      run = adv.value;
      steps += 1;
      assert.ok(steps <= MAX_STEPS, "terminates");
    }
    assert.equal(steps, lesson.steps.length, `${lesson.id} visits every step`);
    // Advancing past done refuses.
    const past = advanceLesson(run);
    assert.ok(!past.ok && past.refusal.code === "lesson-already-done");
  }
  jlog("end-to-end", `"lessons":3`);
});

test("FALSIFIER: an undeclared claim refuses at validation", () => {
  const improvised: Lesson = {
    id: "rogue",
    title: "rogue",
    declaredClaims: ["a declared claim"],
    steps: [
      {
        title: "improvise",
        overlays: ["glyphs"],
        voicedClaims: ["the flyer was actually stable"], // never declared
      },
    ],
  };
  const v = validateLesson(improvised);
  assert.ok(!v.ok && v.refusal.code === "lesson-undeclared-claim");
  assert.ok(!v.ok && v.refusal.message.includes("the flyer was actually stable"));
  // A blind step (no overlays) refuses too.
  const blind: Lesson = {
    ...improvised,
    steps: [{ title: "nothing", overlays: [], voicedClaims: [] }],
  };
  const b = validateLesson(blind);
  assert.ok(!b.ok && b.refusal.code === "lesson-step-blind");
  jlog("undeclared-claim", `"code":"lesson-undeclared-claim"`);
});

test("lesson step caps at cap AND cap+1", () => {
  const mk = (n: number): Lesson => ({
    id: "caps",
    title: "caps",
    declaredClaims: [],
    steps: Array.from({ length: n }, (_, i) => ({
      title: `s${i}`,
      overlays: ["glyphs" as const],
      voicedClaims: [],
    })),
  });
  assert.ok(validateLesson(mk(MAX_STEPS)).ok, "AT cap");
  const over = validateLesson(mk(MAX_STEPS + 1));
  assert.ok(!over.ok && over.refusal.code === "lesson-steps-invalid");
  const empty = validateLesson(mk(0));
  assert.ok(!empty.ok && empty.refusal.code === "lesson-steps-invalid");
  jlog("caps", `"max_steps":${MAX_STEPS}`);
});

test("perception view: verbatim cue-beside-truth with the gap displayed", () => {
  const rows = perceptionView([
    { cue: "pitch-rate", perceived: 0.11, truth: 0.14, units: "rad/s" },
    { cue: "sink-rate", perceived: -0.8, truth: -0.8, units: "m/s" },
  ]);
  assert.ok(rows.ok);
  if (rows.ok) {
    assert.ok(Math.abs((rows.value[0]?.gap ?? 0) - -0.03) < 1e-15, "gap displayed");
    assert.equal(rows.value[1]?.gap, 0);
    assert.equal(rows.value[0]?.perceived, 0.11, "verbatim");
  }
  // Refusals: unnamed cue, non-finite, caps.
  const unnamed = perceptionView([{ cue: " ", perceived: 0, truth: 0, units: "1" }]);
  assert.ok(!unnamed.ok && unnamed.refusal.code === "perception-cue-invalid");
  const nan = perceptionView([{ cue: "x", perceived: Number.NaN, truth: 0, units: "1" }]);
  assert.ok(!nan.ok && nan.refusal.code === "perception-cue-invalid");
  const mk = (n: number) =>
    Array.from({ length: n }, (_, i) => ({ cue: `c${i}`, perceived: 0, truth: 0, units: "1" }));
  assert.ok(perceptionView(mk(MAX_CUES)).ok, "AT cap");
  const over = perceptionView(mk(MAX_CUES + 1));
  assert.ok(!over.ok && over.refusal.code === "perception-cues-invalid");
  jlog("perception", `"max_cues":${MAX_CUES}`);
});
