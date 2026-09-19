# Simultaneous modal contacts

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/modal-multi-contact.performance \
  /tmp/modal-multi-contact.wav --block 37
```

Three modal bodies are described explicitly. The outside bodies start moving
inward across two different gaps; the center body is initially at rest and is
the **only** pressure observer. Both contacts can engage in the same timestep.
The later force at sample 300 and release at sample 340 act on the left body,
not on the receiver. Setting both contact stiffnesses to zero produces silence.
The example has 1201 samples, retaining a short last callback without padding.
It is an authored reduced mechanical model, not a calibrated instrument preset.

## Version 4 input

The existing `music_render modal` command accepts
`frankensim-modal-performance-v4`. All component, bilateral connection and force
records retain the units and meanings documented in `MODAL_PERFORMANCES.md`,
`COUPLED_MODAL_PERFORMANCES.md` and `MODAL_CONTACT_PERFORMANCES.md`. Between the
bilateral connections and `events`, supply:

```text
multi_contact_limits MAX_CONTACTS MAX_SWEEPS MAX_SETUP_TERMS
contacts COUNT
contact_limits MAX_ROOT_ITERATIONS MAX_FORCE_N MAX_PENETRATION_M FORCE_ATOL_N FORCE_RTOL
contact STIFFNESS EXPONENT CHI GAP_M WEIGHT SOURCE_LABEL
contact_left COMPONENT_INDEX SHAPE_0 ... SHAPE_N_MINUS_1
contact_right COMPONENT_INDEX SHAPE_0 ... SHAPE_N_MINUS_1
```

Repeat the last four rows once per contact. Each attachment must use its
component's existing mass-normalized basis. The same source-labelled `fs-dcontact`
potential and nonadhesive Hunt-Crossley law used by version 3 are used unchanged.
All components require `retain-state`; nonlinear static preload is not inferred.
Earlier file versions keep their existing physics, hashes and provenance.

`MAX_CONTACTS` is in 1..=32 and `MAX_SWEEPS` in 1..=128. The setup screen is
`p*(n*(k+2)+(k+1)^2)+n*p^2`, with `n` total modes, `p` contacts, and `k` bilateral
connections; input caps it at 16777216 terms. Each sweep uses at most the supplied
per-contact root iterations plus endpoint/final residual checks. Exceeding
these limits refuses; no force ceiling or tolerance is enlarged automatically.

## Coupled solution, not sequential impacts

`MultiContactModalSystem` computes the full cross-contact displacement response
of the existing bilateral network. Iterated scalar reaction solves use the same
shared-body endpoint motion. The solver requires every contact residual to pass
jointly, stages one mechanical sample with all reactions, then rechecks every
law on the actual states and closes the combined storage/work/loss balance.
One un-converged sequential pass is not accepted. Redundant compliant maps are
supported without requiring an invertible contact-compliance matrix.

`ScheduledRenderer::from_multi_contact_modal_forces` hosts this system using the
existing physical actuator compiler, sample scheduler and streamed WAV encoder.
Direct system trials publish nothing on refusal/cancellation; callback wrappers
retain the existing poisoned-host rule after a failed partly completed callback.
The new provenance records the actual contact count and v4 input identity.

Attachments are fixed and declared, not found by a collision detector. These are
compliant normal contacts, not rigid impacts, free rigid-body integration or
friction. Nonlinear coordinate convergence is not guaranteed; difficult sets
can exhaust the explicit budget. The contact owner allocates, so this is not a
hard-real-time claim. Time refinement remains necessary and nonlinear contact
can alias even when the original linear modes satisfy their bandwidth guard.

Focused native checks (not executed in the authoring environment):

```bash
cargo test -p fs-couple --test modal_multi_contact --test music_render_multi_contact
cargo test -p fs-couple --test modal_contact --test music_render_contact
cargo test -p fs-couple --bin music_render
```
