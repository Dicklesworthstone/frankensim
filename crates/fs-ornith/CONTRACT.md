# fs-ornith CONTRACT

Flagship 1 (plan §15.1, bead mye.2): the ORNITHOID MULTI-INLET
AIRCRAFT, smoke tier — L/D × certified stability × maneuverability,
end to end, with the e-racing payoff measured and every atlas row
carrying its certificates.

## Purpose and layer

Layer **L6 (HELM)**. Composes battle-tested lower crates end-to-end:
fs-bem (panel + Kutta + unsteady wake + the landed cl-adjoint), fs-vpm
(vortex-particle gait metric), fs-race/fs-exec (e-raced generations),
fs-lbm (channel refinement), fs-sos (Lyapunov certificates),
fs-surrogate (conformal bands), fs-dfo (NSGA-II front machinery). No
new numerics live here; the claim is the COMPOSITION with certificates
at every joint.

## Public types and semantics

- `param::OrnithCandidate` — the smoke-tier design: NACA-4 section
  thickness, trim α, inlet chord position, flapping amplitude ×
  reduced frequency; `from_genes` decodes the NSGA gene box.
  `cl_gradient` is the JACOBIAN ACTION: ∂cl/∂α by the fs-bem adjoint
  solve with finite-difference output partials, ∂cl/∂thickness by
  deterministic central difference (documented interim).
  `inlet_mass_flow` reads the upper-surface
  panel speed at the inlet station (suction proxy).
- `screen::{lift_to_drag, flap_metric, screen_generation,
  ScreenReport}` — panel L/D with the DOCUMENTED drag proxy (profile
  cd₀(t) + induced cl²/(π·AR·e); the inviscid panel has no drag of its
  own), the fs-vpm flapping wake metric (BEM-shed vortices advected as
  particles, streamwise vorticity drift), and the E-RACED generation:
  losses normalized against a fixed analytical `LossSpan` per the fs-eproc
  contract (the vessel flagship's measured lesson applied), eliminations and
  the fixed-N-equivalent savings LEDGERED. `screen_generation` returns the
  structured `fs_race::RaceError` rather than manufacturing a winner when the
  support or input contract fails.
- `refine::{refine, RefineReport}` — fs-lbm channel flow around the
  rasterized section; forces by CONTROL-VOLUME momentum balance over
  the PUBLIC cell moments (∮(ρuu + p·I)·n dA, p = ρc_s² — no reliance
  on fs-lbm internals); the report carries its MODEL-FORM HONESTY
  label. Two full flow-throughs (one was MEASURED unsettled:
  steadiness 1.1e-3 → 5.5e-6 after two).
- `certify::{certify, CertifyReport, LdSurrogate}` — the 2-state
  pitch model (stiffness from the BEM lift slope, damping from
  thickness), closed-form Lyapunov P from AᵀP + PA = −I verified by
  fs-sos; the CERTIFIED ROA proxy = P-ellipsoid area under the pitch
  saturation bound, 0.0 WHEN UNCERTIFIED (never pretended); the L/D
  surrogate carries a split-conformal band, coverage GATED.
- `atlas::{build_atlas, Atlas, AtlasRow}` — NSGA-II over (−L/D,
  −ROA, −maneuver, inlet violation); every row carries its stability
  certificate, surrogate prediction, and gene lineage; hypervolume +
  knee attached; the knee design gets an adjoint-direction L/D polish.

## Invariants

1. Certified means CERTIFIED: `roa_volume > 0 ⇔ certificate verified`
   (gated per atlas row).
2. Race losses have a declared maximum paired difference of 1.52: normalized
   base range 1.5 plus total jitter width 0.02. A support breach yields no race
   report (the PairwiseRace contract).
3. Seed replay is bitwise: same seed → identical atlas genes and
   objectives, identical race trajectories (gated).
4. Reports never drop their honesty labels (model-form agreement is
   sign + order-band, not point-matching).

## Error model

Fixture-scale programmer errors still panic (`expect`/`assert`). The e-raced
screen is a `Result` surface: guarantee-voiding loss/span failures propagate as
`fs_race::RaceError`; smoke-tier tests explicitly require success.

## Determinism class

Fully deterministic: BEM/LBM/SOS are deterministic; NSGA-II and the
race jitter derive from fixed seeds through deterministic streams.

## Cancellation behavior

None at fixture scale (seconds per stage). The e-raced screen consumes
an `fs_exec::KillRegistry` and records kills.

## Unsafe boundary

`unsafe_code = "deny"` via workspace lints — none.

## Feature flags

None.

## Conformance tests

`tests/battery.rs`, verdict-JSON rows:

- **orn-001** parameterize: adjoint ∂cl/∂α == central FD to 1.8e-8
  rel; inlet mass-flow responds to the inlet lever (1.35 fore vs 1.10
  aft).
- **orn-002** screen: race winner == deterministic argmax L/D; 23/24
  dominated candidates eliminated early; 1414 evals vs fixed-N 9600
  (6× saved — the P7 payoff measured after checked-span normalization);
  the flap metric responds to
  the gait.
- **orn-003** refine: LBM control-volume lift agrees with the panel
  sign; steadiness 5.5e-6 after two flow-throughs; honesty label
  present.
- **orn-004** certify: Lyapunov certificate verified; certified ROA
  0.66; conformal coverage 0.97 on 60 fresh candidates (target 0.90).
- **orn-005** atlas: 24 certified rows, hypervolume 56.4, knee
  attached, adjoint polish L/D 7.64 → 9.41; roa>0 ⇔ certified on
  every row.
- **orn-006** replay: atlas and race bitwise from seeds.
- **orn-007** what-breaks-first: LBM budget exhausts after one
  refinement → 7/7 candidates degrade gracefully to the
  surrogate+conformal path, 6/7 inside the band — the campaign
  survives its own honesty clause.

## No-claim boundaries

- **F-rep fuselage + IGA-shell wing + manifold harmonics (~200
  coefficients) + local FFD**: the sectional candidate is the smoke
  parameterization; the full chart stack is the recorded successor.
- **Kernel-independent FMM screening at hundreds of candidates per
  generation**: the 3D FMM exists in fs-bem; this flagship screens 2D
  sections — wiring the 3D lifting surface is successor scope.
- **Cumulant LBM/LES on sparse VDB lattices, DWR keyed to L/D, LES
  model cards**: the refinement here is D2Q9 BGK channel flow with a
  control-volume force and an honesty label saying exactly that.
- **Koopman/DMD trim models with conformal e-bands per trim state**:
  the pitch model is a 2-state companion form from the lift slope.
- **NSGA-III reference-point selection**: fs-dfo's landed NSGA-II +
  hypervolume + knee are used; III is named successor scope.
- **LBM adjoints** [M]: deferred per the plan; gradient polish uses
  the BEM adjoint lane only, honestly labeled.
- **FrankenScript study driver + Error Ledger wiring**: the battery
  drives the pipeline directly; fs-script/fs-ledger orchestration is
  the deliverable's production wrapper, tracked by the HELM epics.
