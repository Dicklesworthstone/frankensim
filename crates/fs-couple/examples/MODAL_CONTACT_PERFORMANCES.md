# Contact between vibrating modal bodies

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/modal-contact.performance \
  /tmp/modal-contact.wav --block 37
```

The first component starts moving toward the second across a declared gap.
Only the second component contributes observer pressure, and it has no actuator
force. Contact must engage to produce its sound. The reaction acts on both
components and disappears when they separate. The example also applies and
releases a force on the first component at samples 71 and 120. It has 1201
samples, retaining a short final callback, and a declared full scale of 1 Pa.
Changing `contact 300000000` to `contact 0` disables this force law and makes
the receiver silent; it is not a separately triggered recording or oscillator.

## Version 3 input

Versions 1 and 2 keep their existing meaning. Use
`frankensim-modal-performance-v3` as the header, keep the v2 voices and bilateral
connection records, then insert these four records **before `events`**:

```text
contact_limits MAX_ITERATIONS MAX_FORCE_N MAX_PENETRATION_M FORCE_ATOL_N FORCE_RTOL
contact STIFFNESS ALPHA CHI GAP_M WEIGHT SOURCE_TOKEN
contact_left COMPONENT_INDEX SHAPE_0 ... SHAPE_N_MINUS_1
contact_right COMPONENT_INDEX SHAPE_0 ... SHAPE_N_MINUS_1
```

The existing `coupling_limits` and `connections` records are still required.
`connections 0` declares no permanent bilateral links; a contact can also be
added to an already spring/damper-connected network. Exactly one unilateral
normal contact is supported, not a simultaneous multi-contact solver.

The attachment vectors are in each component's admitted mass-normalized modal
basis, in `1/sqrt(kg)`. Relative closure is left displacement minus right
displacement. Positive penetration is `max(closure - gap, 0)`, in metres.
The `fs-dcontact` potential is `weight * K/(alpha+1) * penetration^(alpha+1)`.
`K` therefore has units `N/m^alpha`, `alpha` is at least one, `chi` has units
`s/m`, and the dimensionless quadrature weight is nonnegative. The law uses
the existing nonadhesive Hunt-Crossley unloading rule. Gap may be signed;
initial penetration contributes to initial energy and faces the penetration
limit. Stiffness or weight zero explicitly disables the contact.

`SOURCE_TOKEN` is one whitespace-free label, at most 1024 UTF-8 bytes. It records
caller provenance; it is not a material validation certificate. The native
`ModalContact` API also accepts an `Obstacle` constructed from a contact receipt.
The decoder does not invent stiffness, damping, a restitution coefficient,
geometry, or a material source.

`MAX_ITERATIONS` is in 1..=128. Up to two bracket endpoint evaluations precede
those refinements, followed by a constitutive check on the actual staged state.
Force and penetration limits are positive physical ceilings, not rescaling
targets. Both force tolerances are positive, with relative tolerance below one.
The force gate uses `ATOL + RTOL * max(applied_force, constitutive_force)`.
The network's energy budget and tolerances apply to **network plus contact**,
including contact dissipation and only the actual external actuator work.

Every component uses `retain-state`. Contact-loaded static preload is refused:
settling only the linear network would omit the contact reaction. Supply an
explicit initial state or simulate the loading history. The chosen model keeps
positive-frequency flexible modes; it is not a free rigid-body hammer model.

## Execution and limitations

`ContactModalSystem` reuses the existing bilateral reaction solve, modal
exact-held-force transitions, `fs-dcontact` potential and `fs-phs` discrete
gradient. Its scalar contact solve uses the whole network's displacement
compliance, so a connected body is not treated as a motionless obstacle.
Candidate states publish together only after contact force, penetration,
component limits and full energy checks pass. Direct cancelled or refused
trials remain retryable; the callback wrapper becomes unusable after a physics
error because an earlier part of that callback may have completed.

`ScheduledRenderer::from_contact_modal_forces` uses the existing force compiler
and sample clock. It supports independent actuator release and cancellation at
completed host callback boundaries. The command writes the existing scaled
PCM16 stream and includes the v3 input identity and contact scope in provenance.
No second WAV encoder or event-time conversion is introduced.

This is a finite-step compliant-contact approximation, not exact continuous
impact dynamics or a rigid constraint. Timestep/mode refinement remains
necessary, especially across engagement and separation. Nonlinear contact can
generate frequencies above the linear mode screen: alias-free sound is not
claimed. The existing discrete-gradient contact path allocates during trials;
iteration limits alone do not certify real-time performance. Tangential
friction, multiple simultaneous contacts, moving geometry and experimental
physical validity are outside this contribution.

Focused native checks (not executed in the implementation environment):

```bash
cargo test -p fs-couple --test modal_contact --test music_render_contact
cargo test -p fs-couple --test modal_coupling --test coupled_render --test music_render_coupled
cargo test -p fs-couple --bin music_render
```
