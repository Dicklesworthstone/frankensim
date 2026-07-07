# fs-scenario — CONTRACT

The boundary-condition and load-case ALGEBRA (plan patch Rev D): a
`Scenario` is a typed value answering "what is being done to the domain?"
— with dimensional analysis on every value, provenance (seed + canonical
IR), and admission-time validity checks that catch the class of mistakes
no solver can fix.

Ambition tags: typed BCs/frames/signals/combos [S]; seeded ensembles
(Dryden, Kanai–Tajimi, Carreau bands) [S]; canonical IR [S].

## Purpose and layer

Layer **L3** (FLUX support). Runtime deps: `std`, fs-qty, fs-rand,
fs-cheb, fs-ga, fs-math. The Design Ledger stores scenarios as canonical
IR artifacts; that integration lives ABOVE L3 (fs-ledger is L6) and is
exercised here through a dev-dependency in conformance tests only.
Consumers: fs-solid, fs-flux, fs-lbm, fs-uq, fs-regime, the milestone
flagships.

## Public types and semantics

- `signal::TimeSignal` — `Constant`, `Ramp` (clamped; the vessel tilt
  `(ramp 0deg 65deg 3s)`), `Table` (strictly increasing times + declared
  `Interp` contract: Linear/Hold, clamped ends), `Chebfun` (fs-cheb
  function object). Every signal knows its `Dims`.
- `frame::FrameTree` — frames with `Fixed`, `Rotating`, and `Tilt`
  motions; poses are fs-ga MOTORS composed down the parent chain
  (`world_pose`), cycle/dangling-parent checked. Rotation about an
  off-origin axis is `T(c)·R·T(−c)`.
- `bc` — `BoundaryCondition { region, physics, kind, value, compatibility,
  frame }`; `expectation(physics, kind)` is the dimensional contract
  table (velocity for flow Dirichlet, kg/s for mass-flow inlets, Pa for
  pressure/traction, K / W/m² / W/(m²K) for thermal, m for elastic
  Dirichlet; no-value kinds; everything else structurally Unsupported).
  Flux-carrying inlets MUST declare `Compat::Incompressible`.
- `ensemble::StochasticEnsemble` — seeded generators: Dryden gust PSD,
  Kanai–Tajimi ground-acceleration PSD (spectral representation with
  Gaussian coefficients — a genuine Gaussian process), Carreau parameter
  bands. `realize(member)` is a pure function of `(seed, model-kernel,
  member)` via Philox `StreamKey` — replayable bitwise.
- `scenario::Scenario` — root value: frames, base BCs, `LoadCase`s,
  factored `Combination`s (`1.2D + 1.6L`), ensembles, `ContactLaw`s
  (Frictionless/Coulomb/Tied), explicit `Environment` (gravity, ambient
  temperature/pressure — REQUIRED constructor argument, never silently
  defaulted). `validate()` returns structured `Violation { code, what,
  fix }` values.
- `ir::write_ir`/`ir::parse_ir` — canonical byte-stable s-expression
  encoding; floats in shortest-round-trip form, dims as explicit exponent
  vectors, integers as integers (u64 seeds are exact).

## Invariants

1. **Round-trip losslessness**: `parse_ir(write_ir(s)) == s` for every
   representable scenario; `write_ir` is byte-stable (canonical form).
2. **Dimensional soundness**: `validate()` rejects any BC/frame/ensemble/
   environment value whose SI exponents disagree with the contract table.
3. **Net-flux compatibility**: if any condition declares
   `incompressible`, either declared mass flows balance to 1e−9 relative
   or a pressure outlet exists — otherwise `flux-imbalance` with the
   imbalance quantified in the message and an actionable fix.
4. **Frame chains terminate**: cycles and dangling parents are violations;
   `world_pose` refuses cyclic chains at runtime too (hop budget).
5. **Bitwise ensemble replay (G5)**: identical `(seed, model, member)` →
   identical realization bits; different member or seed → different draw.
6. **Statistical spectrum match**: the ensemble-averaged periodogram of
   Kanai–Tajimi realizations converges to the target PSD (conformance
   holds band-mean relative error < 15% at 48 members with fixed seed).
7. **Nothing defaulted silently**: `Scenario::new` requires an
   `Environment`; `Environment::earth_lab()` exists but must be named at
   the call site.

## Error model

`ScenarioError`: `Dimensions { context, expected, got }`, `Frame`,
`Evaluate`, `Parse { at, what }`. Validation produces `Vec<Violation>`
(code + what + fix) rather than failing fast — agents get the whole
repair list at once.

## Determinism class

**D0**: signal evaluation, frame poses (fs-ga + fs_math::det trig), and
ensemble realizations (Philox + det trig, fixed draw/summation order) are
bit-identical across runs and platforms. IR text is canonical.

## Cancellation behavior

All operations are bounded (validation is linear in scenario size;
realizations are O(samples × harmonics) with caller-chosen sizes). No
long-running loops; P7 satisfied by boundedness.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-scenario/conformance`):

- **sc-001** rich vessel-pour fixture (all constructs) round-trips
  memory ↔ IR ↔ fs-ledger artifact losslessly; canonical text
  byte-stable.
- **sc-002** seeded violations caught with structured fixes: flux
  imbalance (repaired by adding an outlet), cyclic + dangling frames,
  wrong-dimension BC, undeclared inlet compatibility, unknown combo case,
  kinetic > static friction.
- **sc-003** KT/Dryden/Carreau members bitwise-identical across repeated
  realization; members differ from each other; seed matters; Carreau
  draws stay inside declared bands.
- **sc-004** 48-member KT ensemble periodogram vs target PSD: per-bin
  and band-mean tolerances hold (metrics logged as JSON).
- **sc-005** G0 frame laws: chain composition equals the manual motor
  product; the tilt ramp matches a directly-built motor at five times,
  clamps past its end; points on a rotation axis stay fixed.
- **sc-006** G3 unit coherence: deg/rad and mm/m spellings (via
  fs-qty parsing) converge to the same SI values and the same canonical
  IR; validation is spelling-blind.

## No-claim boundaries

- **Physics vocabulary is v0**: IncompressibleFlow / Thermal /
  Elasticity kinds only. New physics extend `expectation` — adding a
  (physics, kind) pair is a table row plus tests, not a redesign.
- **Region names are strings here**: binding to fs-geom `Region` objects
  (existence, patch-measure integration for velocity-inlet flux) happens
  in the consuming solver layer; net-flux checking covers DECLARED
  mass-flow values, not velocity-profile surface integrals.
- **Recorded ground-motion suites (PEER-class) are not bundled**: the
  `Table` signal is the container for imported records; curation of suites is
  data, not code, and lives with fs-uq.
- **No load-combination EVALUATION**: combinations are typed references
  with factors; assembling factored response quantities is solver-side.
- **The ledger `scenarios` integration is a thin artifact row** (canonical
  IR + seed); a dedicated relational table is deferred to the ledger's
  next schema migration if queries demand it.
