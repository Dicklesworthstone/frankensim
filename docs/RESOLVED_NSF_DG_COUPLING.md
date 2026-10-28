# Resolved compressible flow and the original DG bore

`fs_couple::nsf_dg::NsfDg` connects an original material-owned
`fs_lbm::nsf3::NsfGrid3` to an original `DgPortLoad`. An initial resolved
pulse, vortex or interior obstacle field stays in that fluid owner. Incoming
pressure feedback from the evolving bore reaches the original NSF boundary and
populations during the same accepted coupled interval.

This is a reusable acoustic interface for resolved source zones. It does not add
a flute jet, a reed law, a fitted source coefficient, an output envelope or an
instrument-specific solver.

## Physical interface

The fluid zone must already have its original characteristic axis. Select either
complete face; its outward normal flow is positive into the DG inlet. The other
face retains a caller-supplied incoming characteristic in metres per second.
For a zero-velocity reference, the original boundary uses

`J− = u_normal − (p − p_reference)/(rho_reference c_reference)`.

Only the incoming characteristic on the coupled face is solved. For every trial,
the original NSF owner advances from the same accepted state while the original
DG owner supplies its current step-mean law `p_DG = H + Z Q`. The scalar solve
matches the mean reference-relative NSF face pressure to that law, using the
face's actual integrated continuum volume flow.

The complete reference gas state, including transport properties, must equal the
DG inlet gas. The actual complete face area must match inlet area. The DG interval
must equal a positive integer number of original NSF steps, allowing only
floating-point roundoff in derived clocks. Both characteristic reference
velocities must be zero. These conditions prevent an unowned mean-flow interface
or independently chosen atmosphere from silently entering the acoustic path.
The original thermochemical model, caloric chart, transport support and all
original per-cell material refusals remain active.

Each candidate retains the complete NSF populations, force chart, energy,
material caches and histories. Both fluid and DG candidates remain unpublished
until root convergence, original-source admission, integrated port-work checks,
the caller's projection budget and cancellation checks all pass.

## Work and transported quantities

Every accepted frame exposes both original NSF faces, original fluid and DG
accounts, actual mass and energy changes, and original downstream terminal flows.
The NSF quadrature directly accumulates `integral((p-p_reference) u_normal) dA dt`;
it does not subtract large absolute-pressure work totals afterward.

The frame separately reports:

- **Projection covariance:** resolved acoustic work minus
  `mean(p′) mean(Q) interval`, including transverse and temporal covariance.
- **Acoustic work residual:** resolved acoustic work minus actual integrated DG
  inlet work. Its magnitude must fit the explicit absolute/relative budget.
- **Port work residual:** original DG inlet work minus its own mean-pressure and
  held-flow product.
- **Mass projection difference:** actual lattice mass flux minus
  `rho_reference mean(Q) interval`.

The acoustic duct does not transport the full NSF mass, entropy or enthalpy.
Actual NSF mass and total-energy fluxes remain visible and are not credited as
heat or forced into a false whole-system balance. The original NSF viscous-heat account and DG loss
accounts retain their existing budget-remainder qualifications.

## Pressure playback

Construct a cold `NsfDg` from the retained zone, retained load and explicit
`NsfDgLimits`. Limits bound grid cells, NSF substeps, all root-evaluation
cell-steps, iteration count, finite source duration, the incoming-characteristic
bracket and both pressure/work tolerances. Bracket endpoints must themselves be
admissible source states; an out-of-domain endpoint refuses rather than clipping
a physical amplitude.

`nsf_dg::performance::NsfDgPerformance::new` accepts that system, an
`ApertureObservation`, the matching sample rate, a finite `Vec<f64>` of incoming
characteristics on the other face, a maximum callback length and an
`Arc<CancelGate>`. One input is held per complete coupled interval. Choose
`ApertureObservation::Inlet` for the actual DG mean inlet pressure or
`TubeBaffled(receiver)` for the existing computed-radiation-bound microphone.
The microphone consumes the actual accepted terminal flows and retains its
physical/numerical delay.

Call `with_observer_gas(&propagation)` before playback to use the same physical
atmosphere law as native flute and reed observation. `GasPropagation` selects
the shared law and finite path/filter budgets; the actual DG radiation source
supplies the complete exterior gas, including viscosity, conductivity and heat
capacities. Each original surface-to-microphone path retains its own absorption,
dispersion and causal delay. `observer_gas()` exposes the original material and
path-response report. The fluid and bore reactions are unchanged, and the source's
existing cancellation gate also owns this cold preparation.

The result implements the existing `PressureRenderer`, so the original pressure
decimator, ensemble and PCM writers can consume it. Callback partitioning does
not change either physical owner or source timing. A failed physical step rolls
back both owners; a later observation failure preserves the accepted physical
prefix and poisons further playback.

## Scope and refinement

This first composition is a uniform acoustic-mode projection with a held DG flow
and the original first-order characteristic boundary. Its plane placement,
transverse-mode content, NSF spatial/temporal resolution and projection work
must be refined for each physical source. A loose allowed work residual is not
an accuracy certificate. No full resolved flute/oboe jet qualification,
nonreflecting-boundary target, thermal-transport matching or perceptual accuracy
claim follows from these source-transaction and coupling tests.

Focused tests are in `fs-lbm/tests/nsf3_characteristic_trial.rs` and
`fs-couple/src/nsf_dg/tests.rs`. They cover exact original-step replay, partial
candidate cancellation, reference-relative work, both face signs, real downstream
feedback, gas-pressure substitution, thermochemical support, projection refusal,
cold area/clock/work admission, and physical-pressure/PCM partition identity.
