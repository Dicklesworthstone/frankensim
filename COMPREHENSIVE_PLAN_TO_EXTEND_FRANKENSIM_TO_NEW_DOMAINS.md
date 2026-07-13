# COMPREHENSIVE PLAN TO EXTEND FRANKENSIM TO NEW DOMAINS

> Mechanisms, Contact, Flexible Bodies, Tribology, Electromagnetics, Power
> Electronics, Control, Thermal/Reactive Flow, Acoustics, Degradation, and
> Coupled Machines on a Sheaf-Cohomological and Port-Thermodynamic Spine.
>
> Companion to `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` and
> `COMPREHENSIVE_ADDENDUM_TO_FRANKENSIM_PLAN.md`. The master plan remains the
> constitution, the addendum remains its refinement, and this document is the
> expansion charter. It answers one question:
> **what is preventing FrankenSim from simulating a Geneva drive, an involute
> gearbox, a Wankel rotor, a PM synchronous motor, a generator, and an internal
> combustion engine — with real material properties and honest evidence — and
> what is the shortest mathematically serious path to doing all of them at a
> world-class level?**
>
> Ambition tags follow the house convention: `[S] Solid`, `[F] Frontier`,
> `[M] Moonshot`. Nothing tagged `[M]` may sit on the critical path. Audited
> **current-state** claims are grounded in live files/contracts; proposed
> mathematics and capabilities are tagged and governed by their proof
> obligations, while external definitions and benchmark families are anchored
> in Appendix D. The implementation pass pins the exact repository commit and
> source artifact for every resulting Bead.
>
> **Amendment rule.** The master plan remains authoritative for mission,
> layering, dependency policy, and the Decalogue. Where this charter explicitly
> adds assumptions or narrows an aspirational master-plan statement—for example,
> that exact incidence alone removes every spurious mode, that IPC is
> unconditionally intersection-free, or that a Dirac interconnection makes an
> arbitrary partitioned time discretization passive—the stricter statement is
> the governing amendment for these new domains until the master text is
> updated. Live crate contracts remain the authority for what is implemented
> now. This prevents document precedence from resurrecting a claim this review
> has already falsified.

### Normative delta and ratification register

This charter is not allowed to become a contradictory third source of truth.
The following deltas are intentional amendments, not accidental drift:

| Delta | Earlier text being refined | Governing rule in this charter |
|---|---|---|
| Coupling/passivity | A Dirac or power-preserving interconnection can read as passive “by construction” | Only the interconnection is lossless by construction. Component storage/dissipation/source laws, time discretization, transfer, iteration, and the closed accounting-window audit determine a coupled passivity claim (§3.7, §6). |
| Contact ownership | Contact appears inside `fs-solid` | Reusable detection/response protocols move to L3 `fs-contact`; `fs-solid` and `fs-mbd` consume adapters. Generic nonlinear/conic algorithms remain in L1 `fs-solver`. |
| Spectral ownership | Some domain eigensolvers appear in `fs-la` | Generic operator spectra, nullity, continuation-health, and multiplier extraction belong to L1 `fs-spectral`; domain crates assemble operators and interpret results. |
| Dimensional algebra | `Qty` and Buckingham-π use five base dimensions | The amount-of-substance exponent is added atomically, with versioned five-to-six wire migration and semantic crosswalks (§7.3). |
| Phase sequencing | The older roadmap did not contain the machine-domain prerequisite graph | The E0–E8 dependency/exit gates in §10 govern this expansion, including a dry-tribology baseline before the E2 Geneva exit. |
| Research governance | The addendum's “one research bet” can be read as one bet for the entire program | One unproven mechanism is admitted **per independently falsifiable proof lane**, under an explicit portfolio WIP/budget cap. Many such lanes may run in parallel; no acceptance gate may hide two unproven mechanisms behind one result. This increases theorem throughput without weakening falsifiability. |

The implementation program creates a ratification Bead that mirrors every
accepted delta into the master plan/addendum and records the exact affected
sections. Until that synchronization lands, this table wins for the named
deltas; it does not authorize unrelated layer or dependency-policy exceptions.
These are present-tense proof boundaries, not prohibitions: a stronger theorem
for a precisely quantified class may supersede any one of them through the
Theorem Foundry, independent reproduction, and an explicit ratification diff.

---

## 0. Executive Summary

FrankenSim today is a broad, working library substrate for certified geometry,
fields, optimization, evidence, and replay, but it is **not yet a packaged
general-purpose machine simulator**. Its strongest paths are still predominantly
static-geometry field and shape workflows: geometry in, elliptic/parabolic
physics on it, adjoints out, evidence attached. The machines
the user wants to simulate — Geneva mechanisms, gear trains, rolling
constant-width mechanisms, Wankel rotors, electric motors, generators, and ICE
engines — are *machines*: assemblies
of bodies in relative motion, exchanging force, flux, heat, charge, and gas
through contacts, air gaps, and ports. Eight pillars are missing or materially
incomplete; one of them, coupling, exists today only as a scalar seed:

| # | Missing pillar | One-line diagnosis |
|---|---|---|
| P1 | **Motion** | Geometry is static per-evaluation. Rigid-motion algebra exists (`fs-ga` PGA motors) and time-parameterized frame trees exist (`fs-scenario::FrameTree`), but **no chart can move**: nothing binds a motor path to a `Chart`, so there are no swept volumes, no envelopes, no moving boundaries. |
| P2 | **Interaction** | No multibody dynamics, no joints, no kinematic constraints, no DAE machinery, no body-to-body contact, no continuous collision detection, no penetration depth, no event location. The only rigid-body code in the tree is a torque-free attitude step (`fs-time/src/lie.rs:53`). |
| P3 | **Fields beyond mechanics** | No electromagnetic, thermal-domain, or acoustic solver stack. `fs-feec` already ships an exact discrete de Rham complex (Whitney 0–3 forms, Nédélec H(curl), Raviart–Thomas H(div), Hodge stars, computable cohomology), and `fs-cutfem` ships a generic constant-coefficient 2-D steady scalar Poisson/diffusion kernel. Those are valuable seeds, but neither is yet a typed Maxwell, transient/nonlinear heat, or acoustic formulation with material laws, sources, BCs, error policy, and validation. |
| P4 | **Gas** | No compressible flow, no equation of state, no species, and no combustion. All production fluid paths are incompressible or low-Mach LBM. An ICE thermofluid cycle cannot yet be expressed at a useful physical fidelity. |
| P5 | **Matter, interfaces, and history** | `fs-material` covers a useful mechanical-law kernel but not a cross-physics database. There are no electrical, magnetic, thermal, transport, tribological, fatigue, fracture, creep, corrosion, process-state, or interface-property records with provenance, covariance, history variables, and validity domains. |
| P6 | **Coupling (stub)** | `fs-couple` has a useful lossless-port/power-balance seed and Aitken, but only a *scalar* added-mass fixture. It lacks explicit dissipative/storage/source relations, stream balances, field transfer, mortar, vector IQN-ILS, monolithic escalation, and multirate stability machinery. |
| P7 | **Machine engineering stack** | No flexible multibody reduction, rotor/bearing dynamics, lubrication/EHL, seal leakage, fatigue/fracture/wear evolution, acoustics/NVH, power-electronic device lane, closed-loop controls, calibration/data assimilation, or tolerance-aware reliability loop. Without these, a mechanism can move in a demo but cannot yet earn a world-class machine-design claim. |
| P8 | **Assurance and adoption** | No machine-level hazard/safety-case, EMC/insulation/containment, standards-conformance, experimental V&V, interoperable engineering workflow, competitive benchmark ledger, or at-scale qualification program. Scientific evidence cannot safely drive engineering decisions without these. |

The central strategic claim of this plan: **all eight pillars are extensions of
the mathematics FrankenSim already committed to, not additions of foreign
engines.** Specifically:

1. **Bulk constitutive operators, interface laws, and internal variables are
   distinct first-class artifacts.** In FEEC/DEC, the exterior derivative `d`
   is exact integer topology. Permittivity, reluctivity, conductivity, heat
   capacity, and similar bulk metrics enter weighted mass/Hodge operators;
   plasticity, hysteresis, chemistry, friction, wetting, contact conductance,
   damage, and wear do **not**. They live in typed local, interface, reaction,
   and history-dependent constitutive maps. This separation preserves `dd = 0`
   without pretending topology proves constitutive stability.
2. **A mechanism has a nonlinear configuration object and a family of tangent
   constraint complexes.** Finite SE(3) loop closure is non-commutative
   holonomy. At a configuration `q`, the constraint Jacobian/rigidity complex
   exposes infinitesimal motions through `ker J(q)` and self-stress through
   `ker(J(q)ᵀ)`; a cellular sheaf can organize compatible local linearizations
   when its restriction maps are actually defined. This precise split retains
   the cohomological leap while avoiding false identifications of nonlinear
   assembly, mobility, and overconstraint.
3. **Motion certificates extend interval machinery with validated dynamics.**
   Swept/envelope and exact-path event enclosures can reuse `fs-ivl` Taylor and
   Krawczyk primitives. Simulated dynamics additionally require a true-flow
   tube and class-qualified complete root accounting; ordinary dense output is insufficient.
   This creates a differentiated no-silent-miss lane without making an
   unsupported universal claim about other tools.
4. **Machines compose through storage, lossless interconnection, dissipation,
   and reaction artifacts.** Motor torque ↔ shaft speed,
   winding voltage ↔ current, gas pressure ↔ volume flow, heat ↔ entropy flow:
   the port-Hamiltonian bus in `fs-couple` is the lossless interconnection
   algebra for an engine or motor-generator set. Resistive, frictional,
   viscous, plastic, thermal, and chemical processes require explicit
   dissipative/reaction relations. Power, energy, entropy, and exergy balances
   are separate runtime audits; a Dirac structure alone does not prove
   passivity of its components or of a time-discrete co-simulation.
5. **Rotating machines exploit the symmetry crate.** A gear with N teeth, a
   motor with p pole pairs, and individual repeated rotor/stator components may
   be Cₙ-symmetric. A three-flank Wankel rotor, two-lobe housing, ports, load,
   defects, and excitation do not magically give the assembled machine one Cₙ.
   `fs-symmetry`'s circulant block-diagonalization is the right foundation when
   the full operator, material state, boundary conditions, and excitation
   preserve the declared group, extended to complex sector coupling and
   (anti)periodic boundary conditions.

The plan is organized as: gap analysis grounded in the current tree (§1),
extension doctrine (§2), the mathematical spine (§3), the new crate atlas (§4),
per-domain technical plans (§5), coupling architecture (§6), solver and
infrastructure upgrades (§7), flagship demonstrators (§8), the Gauntlet program
for new physics (§9), a phase roadmap with exit criteria (§10), an EV-ranked
lever matrix with recommendation cards (§11), proof obligations and failure
modes (§12), galaxy-brain transparency cards (§13), claim boundaries (§14),
and operationalization (§15).

---

## 1. Grounded Gap Analysis — What Exactly Is Blocking, With Citations

This section records what the tree can and cannot do at the audit snapshot, so
the plan's deltas are auditable. Each row was checked against `main`; the
operationalization pass must pin the exact commit again before converting a
row into a Bead because this is an active shared repository.

**Audit snapshot:** `main` at
`51e50ff89e3ef8cf7eb405b8dfd0f463cf50a0cd` (2026-07-13). The shared worktree
also contained unrelated in-flight changes; no claim here treats those dirty
paths as committed architecture, and operationalization must refresh the scan.

### 1.1 What exists and is load-bearing for this plan

| Capability | Where | Why it matters here |
|---|---|---|
| Full SE(3) rigid-motion algebra: PGA motors, screw exp/log, slerp, renormalization | `fs-ga/src/pga.rs:273-457` | A coordinate-robust motor representation for multibody dynamics and mechanism kinematics already exists, with normalization behavior that still belongs in every proof receipt. |
| Time-parameterized rigid frame trees (`FrameMotion::Rotating`, `world_pose(id, t)`) | `fs-scenario/src/frame.rs:26-142` | A rotor-in-stator kinematic tree already types; it just never touches geometry. |
| Discrete de Rham complex: Whitney 0–3 forms, exact incidence `dd=0`, Galerkin + diagonal Hodge stars | `fs-feec/src/{whitney.rs,hodge.rs}`; incidence in `fs-rep-mesh/src/complex.rs:95-129` | A major substrate for Maxwell, heat, and Darcy. Weighted operators, formulations, sources, BCs, stability and solvers remain real work. |
| Nédélec H(curl) and Raviart–Thomas H(div) families, r = 1..4, with discrete grad/curl/div chain maps | `fs-feec/src/highorder/vecfam.rs:362-1050` | The exact element families EM needs. Explicitly deferred to a curl-curl solve: `fs-feec/CONTRACT.md:298-301`. |
| Computable cohomology: exact Betti (Bareiss), Hodge decomposition, harmonic basis, cycle pairing (`circulation`), harmonic deflation | `fs-feec/src/{betti.rs,cohomology.rs}` | Multiply-connected EM domains, flux linkage, and circuit coupling through H¹ generators are *already computable*. |
| Real cellular-sheaf watertightness: δ⁰/δ¹ coboundaries, H¹ obstruction, gauge section solve, ray-parity falsifier, repair/merge | `fs-geom/src/sheaf.rs:41-488`, `sheaf_repair.rs`, `sheaf_merge.rs` | A reusable artifact/falsifier pattern for compatibility problems; it does not prove kinematic or conservative-transfer semantics verbatim. |
| Certified separation lower bound between two bodies (positive gap only) | `fs-query/src/lib.rs:619-733` | Seed of certified contact detection; needs penetration depth and motion. |
| IPC barrier contact vs. fixed SDF obstacle, lagged smoothed Coulomb friction (feature-gated) | `fs-solid/src/contact.rs` | The smooth contact lane exists in embryo; CCD explicitly deferred (`contact.rs:12-13`). |
| Generic constant-coefficient 2-D steady scalar Poisson/diffusion CutFEM on certified SDF cuts | `fs-cutfem/CONTRACT.md:5-9,34-39`; `fs-cutfem/src/fem.rs` | A direct seed for a first steady isotropic heat/current manufactured slice; it is not yet a thermal-domain solver, variable-coefficient constitutive assembly, or transient energy balance. |
| Interval/affine/Taylor-model arithmetic, Krawczyk/Newton certified roots, exact predicates | `fs-ivl` | The engine for certified events, certified sweeps, certified conjugate action. |
| Symplectic Verlet (+discrete adjoint via revolve), SO(3) Lie integrator, generalized-α, IMEX, RK45 | `fs-time/src/*` | Time integration spine; needs SE(3), constraints, events, sparsity. |
| Port-Hamiltonian Dirac ports (mechanical/fluid/thermal), exact interface power audit, Aitken | `fs-couple/src/lib.rs:90-305` | The multiphysics bus; needs electrical/rotational/chemical port kinds and field-level transfer. |
| Cₙ representation theory: circulant DFT block-diagonalization, symmetrization with perturbation bounds | `fs-symmetry/src/lib.rs:108-300` | Gear teeth, motor poles, Wankel flanks. Needs complex sector coupling + Dₙ. |
| Units: `Qty<m,kg,s,K,A>` — **ampere and kelvin are already base dimensions**; Pa·s already aliased | `fs-qty/src/lib.rs:39,156,373-381` | Volts, teslas, henries, W/(m·K), H/m are expressible *today*; only aliases and parser tokens are missing. Mole/candela absent by documented no-claim (`lib.rs:10-12`). |
| Evidence colors with enforced no-laundering, origin-typed claims, receipt re-execution | `fs-evidence/src/color.rs:420-893`, `fs-package/src/lib.rs:469-2181` | New domains can reuse the honesty apparatus once every load-bearing query/operator emits and composes the required receipt. |
| Level-set advection (WENO5 + TVD-RK3), FIM redistancing, velocity extension | `fs-topols/src/{weno,fim,veloext}.rs` | Moving-interface machinery reusable for phase change and moving boundaries. |
| LBM D2Q9/D3Q19 BGK with free-surface prototype and thermal fixture | `fs-lbm/src/*` | Fluid workhorse; needs moving solid boundaries and stronger collision models. |
| Campaign pattern: parameterize → e-race screen → refine → certify → atlas → notebook | `fs-ornith`, `fs-vessel`, `fs-thrust-e2e`, `fs-wasm/flagships.rs` | The delivery vehicle for the new flagships (§8). |

### 1.2 What is absent (the blockers)

| # | Blocker | Evidence |
|---|---|---|
| B1 | No binding of motion to geometry: no `MovingChart`, no swept volume, no envelope, no general Minkowski sum (ball case only), no ALE | `fs-query/src/lib.rs:611` (`minkowski_ball` only); fs-query CONTRACT no-claim on general Minkowski; live production API/contract review |
| B2 | No multibody dynamics: no joints, no constraint DAE, no Lagrange multipliers, no index reduction/GGL/Baumgarte, no SE(3) integrator (SO(3) only), no event location | `fs-time/CONTRACT.md:144-152`; live production API/contract review |
| B3 | No body-to-body contact: no CCD, no penetration depth (`ClearanceField` clamps at 0: `fs-query/src/lib.rs:634`), no Signorini/mortar, contact only deformable-vs-fixed-obstacle | `fs-solid/src/contact.rs:12-13`; live production API/contract review |
| B4 | No production mechanism vocabulary or solver surface: gear, cam, linkage, Geneva, conjugate action, and transmission error are absent from the implemented public architecture even where isolated words may appear in plans/tests | live crate/API sweep |
| B5 | No electromagnetics: no curl-curl assembly, no gauge handling, no EM material fields, no complex/blocked Krylov, no EM forces | `fs-feec/CONTRACT.md:298-301`; `fs-regime/src/scaling.rs:25` (placeholder); fs-solver real-only (`fs-solver/src/op.rs:10-21`) |
| B6 | No thermal-domain solver stack: the generic `fs-cutfem` constant-coefficient steady scalar diffusion kernel, BC types (`fs-scenario/src/bc.rs:105-129`), an adjoint heat fixture (`fs-adjoint/src/timedep.rs:15-44`), and an LBM fixture do not yet provide typed conductivity/storage, transient nonlinear energy transport, radiation/contact/phase change, or thermal validation | cited files and live API/contract review |
| B7 | No production compressible/reactive-gas transport, EOS/species, or combustion API | live crate/API review; Godunov/WENO names that do exist belong to level-set transport, not a gas-dynamics implementation |
| B8 | Materials are mechanics-only; no cross-physics property records; key implemented law parameters use raw `f64` rather than a typed cross-domain property-query path | `fs-material/src/elastic.rs:11`, CONTRACT |
| B9 | Coupling is a scalar fixture: no field transfer, no mortar, no vector IQN-ILS, no multirate | `fs-couple/CONTRACT.md:69-79` |
| B10 | Implicit integrators are dense-only (`galpha.rs`, `stiff.rs` take dense operators); adjoints only for Verlet and BE-heat | `fs-time/CONTRACT.md:153` |
| B11 | `fs-spectral` does not own the production sparse kinematic/sheaf complex this plan needs; its contract places the sparse sheaf-Laplacian path elsewhere | `fs-spectral/CONTRACT.md:73-76` |
| B12 | fs-feec is 3-D simplicial only; no 2-D complexes needed for planar magnetics and other planar field formulations | `fs-feec` fixtures all `TetComplex`; `HexComplex` storage-only |
| B13 | No flexible-machine, rotor/bearing/seal, power-electronics/control, life/degradation, acoustic/NVH, or closed-loop validation stack | live crate/contract atlas and repo-wide symbol sweep |
| B14 | No operational machine V&V, hazard/safety/EMC case, standards-conformance axis, end-to-end engineer/interchange workflow, or scale/competitive qualification program | live crate/contract atlas and workflow audit |

Everything in §3–§10 is the ordered, evidence-gated program for removing B1–B14.

---

## 2. Extension Doctrine

Twelve principles govern the extension. They are additions to the Decalogue, not
replacements, and every new crate CONTRACT must cite the ones it implements.

- **D1 — Motion is a first-class citizen.** A pose is a `Motor`; a trajectory is
  a motor path with a certified enclosure class; a moving shape is a
  `SpacetimeChart` whose immutable snapshots implement the static `Chart`
  contract. Static geometry is the identity-path special case, not an
  assumption baked into every API.
- **D2 — Matter is a typed constitutive graph.** Bulk constitutive weights are
  split into storage/metric laws (ε, μ/ν, ρc) and transport/dissipation laws
  (σ, k and their nonlinear/tensor successors); each enters the appropriate
  weighted mass, Hodge, flux, or dissipation operator. Interface laws
  (friction, wetting, electrical/thermal contact), reactions, and path-dependent
  internal variables (plasticity, hysteresis, damage, chemistry, wear) are
  distinct typed maps. No property is an ad hoc scalar buried in a kernel, and
  no Hodge-star claim is extended beyond the mathematics it actually proves.
- **D3 — Events are certified.** Any discrete transition in continuous time
  (impact, contact switching, valve opening, diode commutation, ignition) must
  carry a class-valid `CertifiedEvent`/`NoEventCertificate`, a typed
  `PossibleEvent`/`Unknown` box, a set-valued simultaneous/grazing/Zeno or
  inclusion outcome, an explicit refusal, or an `Estimated` baseline label.
  Only the certified-event/no-event class may claim complete isolated-event
  accounting. Silent event
  misses are the dynamics analogue of tunneling rays and are treated as Sev-0
  certificate bugs.
- **D4 — Mechanisms use nonlinear descent plus tangent cohomology.** The solid
  foundation represents finite configurations as SE(3) assignments satisfying
  nonlinear constraint and holonomy equations. A groupoid/stack realization is
  promoted only when its restriction/descent data and practical localization or
  composition benefit are explicit. At each configuration a linearized constraint
  complex or rigorously defined cellular sheaf yields infinitesimal mobility,
  compatibility, and dual self-stress artifacts. Unilateral contact remains a
  cone/complementarity problem. Rolling without slip, knife-edge constraints,
  and other velocity constraints live in a nonintegrable Pfaffian distribution,
  not in a fictitious position-level loop closure. The distinctions are
  executable and never collapsed into one decorative `H¹` slogan.
- **D5 — Machines compose through typed ports and explicit constitutive
  relations.** Cross-domain coupling uses power-conjugate ports, lossless Dirac
  interconnections, and separately typed storage, dissipation, reaction, and
  source elements. Per-exchange power is audited; open-system energy, entropy,
  mass/species, charge, and exergy balances include boundary and advective
  terms. No solver reaches into another solver's private state.
- **D6 — Correlation predictions are Validated or Estimated, never promoted by
  code verification alone.** Engineering correlations
  (Woschni heat transfer, Wiebe burn, Chen–Flynn friction, Steinmetz core loss,
  Archard wear, AGMA factors) are admitted eagerly — with model cards, validity
  domains, and `Validated`/`Estimated` colors. Their implementation and algebra
  may be numerically `Verified`; their physical prediction remains tied to a
  named validation domain and cannot inherit that color. This makes honesty a
  testable product feature rather than a slogan about incumbent tools.
- **D7 — Reduced fidelity is a lane, not a lie.** 0-D/1-D system models
  (cycle-averaged engine thermodynamics, lumped thermal networks, MNA circuits)
  are first-class citizens with their own Gauntlet batteries because they are
  indispensable design-space and system-level models. 3-D field solves refine and calibrate the
  reduced lanes through the surrogate certify-or-escalate machinery.
- **D8 — Every new domain ships its adjoint or its excuse.** Design optimization
  is the mission. A new solver merges with a discrete adjoint and a gradient
  verification gate, or with a CONTRACT no-claim naming exactly why (nonsmooth
  event dependence, free-surface non-differentiability) and what generalized,
  smoothed, differentiable-relaxation, or derivative-free optimization lane is
  valid instead.
- **D9 — Verification is not validation.** Discretization proofs, code
  verification, benchmark agreement, material-data credibility, model-form
  validation, and prediction uncertainty are separate evidence axes. A result
  may be numerically Verified and physically only Estimated. Claim receipts
  carry the load-bearing dependency subgraph, so an irrelevant weak datum does
  not demote an unrelated output and a relevant one cannot be laundered away.
- **D10 — Real machines age and vary.** Manufacturing process, lot, texture,
  residual stress, surface finish, lubricant state, temperature/rate/history,
  tolerances, fatigue, creep, corrosion, demagnetization, damage, leakage, and
  wear are state, not footnotes. The flagship definition of done includes
  uncertainty and life-cycle evolution, not only a nominal first cycle.
- **D11 — Closed-loop operation is part of the plant.** Sensors, estimators,
  power electronics, controllers, saturation, quantization, delay, and faults
  are simulated with the same budgets and evidence discipline as the physical
  machine. Stability, reachability, and bifurcation/Floquet artifacts accompany
  time histories.
- **D12 — One research bet per proof lane, many ambitions in the portfolio.**
  The program may remain maximalist, but any acceptance gate isolates one
  unproven mechanism at a time against a boring baseline. Frontier and moonshot
  features have activation criteria, kill criteria, fallback lanes, and
  retained falsifiers; ambition is accumulated by certified ratchets.

---

## 3. The Mathematical Spine

This section is the heart of the plan: eleven constructions that make the new
domains *native* to FrankenSim's sheaf-cohomological architecture rather than
ported-in foreign code. Each construction names its artifacts (types,
certificates, monitors) per the artifact contract — math that does not compile
to a runtime artifact does not appear here.

### 3.1 Constitutive Operator Graph — topology, metric, interfaces, and memory `[S→F]`

**The observation.** FEEC already separates the exact incidence complex
(`dd = 0`) from metric-dependent operators. That is a powerful start, but it is
not a theorem that every material property is a Hodge star. The executable
decomposition is:

| Role | Examples | Artifact |
|---|---|---|
| Topology/balance | incidence, boundary, grad/curl/div complexes | exact integer `d`, trace, boundary and orientation maps |
| Bulk storage/metric | heat capacity, permittivity, permeability/reluctivity, density/inertia | weighted mass/Hodge operator with units and positivity/coercivity bounds in the formulation's admissible regime |
| Bulk transport/dissipation | Fourier, Ohm, Darcy, viscosity | monotone flux operator or dissipation potential |
| Reversible/cross-coupled block | Hall and gyroscopic/Coriolis terms; piezoelectric, magnetoelastic and thermoelectric coupling | block operator/free-energy Hessian plus declared skew/nondissipative part, variable/time-reversal parity, Onsager–Casimir/Kelvin reciprocity and coupled stability |
| Interface law | friction, wetting, contact/gap conductance, radiation exchange, electrical contact, membrane/permeation law | oriented trace relation with action/reaction and power/entropy audit |
| Reaction/source | combustion, electrochemistry, impressed current, remanence/magnetization, sharp-interface phase transfer | stoichiometric/source operator with conservation and admissibility constraints appropriate to that formulation |
| Internal memory | plasticity, damage, hysteresis, viscoelasticity, fatigue, wear, demagnetization | state-transition law with thermodynamic admissibility where applicable, an explicit empirical no-claim otherwise, and a consistent tangent when claimed |

**Artifacts.** `WeightedStar`/`WeightedMass` assemble scalar, SPD-tensor, and
nonlinear bulk weights over existing Whitney/Nédélec/RT spaces; typed block
operators carry reversible skew and cross-coupling without pretending each
block is SPD. A sibling
`ConstitutiveGraph` composes local bulk, interface, reaction, and history maps
without forcing them into a matrix abstraction. Every node declares units,
state, calibration domain, differentiability class, energy storage,
dissipation/entropy production, and a consistent tangent/VJP when claimed.
Piecewise-constant weights land first; quadrature-point anisotropic,
frequency-dependent, and history-dependent maps follow behind `[F]` flags.

**Actual proof obligations.** Weighted operators cannot alter the algebraic
identity `dd = 0`, but that fact alone proves neither stability nor absence of
spurious modes. Each field formulation must also establish the relevant
subcomplex and bounded commuting projection, boundary/gauge treatment,
coercivity or inf-sup condition, discrete compactness/eigenmode behavior where
applicable, material positivity/monotonicity, quadrature consistency, and
convergence of the quantity of interest. These artifacts are the difference
between reusing FEEC and merely borrowing its vocabulary.

**Immediate consequences.** Heat conduction uses a storage mass operator and
a conductivity operator, for example
`M_{ρc}(T) Tdot + d₀ᵀ M_k(T) d₀ T = f + b`; electrostatics and current flow use
the corresponding weighted scalar complexes; 2-D `A_z` magnetics is a
nonlinear scalar curl/Poisson formulation; and 3-D magnetics uses the Nédélec
complex plus gauge and harmonic-space handling. The existing elements make
these paths unusually short, but realistic motor/engine claims still require
boundary/source models, nonlinear solves, force-functional error estimates,
moving-conductor terms, preconditioners, and validation. No calendar promise
is made before an end-to-end baseline is measured.

### 3.2 Nonlinear Kinematic Descent and the Tangent Constraint Complex `[F/M]`

One ordinary linear sheaf cannot simultaneously represent finite SE(3)
configurations, infinitesimal mobility, self-stress, unilateral contact,
nonholonomic rolling, and hybrid mode changes. FrankenSim will implement the
stronger five-level object:

1. **Finite configuration/groupoid layer.** Bodies carry motors `M_i`; each
   joint supplies a nonlinear relation `C_e(M_i,M_j,p_e)=0` or an admissible
   relation for unilateral contact. A loop product is non-commutative holonomy,
   not an additive graph cocycle. `motor_log` gives a local residual only inside
   its chart/reach; global continuation tracks branches, monodromy, and chart
   transitions explicitly. A sheaf/stack of local solution charts is the
   moonshot representation once the plain constraint-groupoid baseline works.
2. **Configuration-dependent tangent complex.** At a specific `q`, assemble
   the geometrically exact Jacobian `J(q)=DC(q)`. After grounding or quotienting
   rigid-body gauge modes, `ker J(q)` is the space of *infinitesimal* motions.
   This is not automatically the dimension of a finite configuration manifold:
   singular mechanisms require second-order rigidity tests and validated
   continuation. Where joint-local restriction maps genuinely define a
   cellular sheaf `K_q`, its `H⁰` is identified with this kernel and is useful
   for localization and distributed assembly.
3. **Nonholonomic distribution layer.** Ideal rolling/no-slip, knife-edge, and
   similar constraints are velocity-level relations `A(q) qdot = 0`. They are
   not rewritten as `Φ(q)=0` unless Frobenius/integrability is actually proved.
   The admissible distribution, its bracket-growth/curvature and holonomy, and
   the Lagrange–d'Alembert reaction/work convention are explicit artifacts.
   An ideal no-slip constraint lane and a resolved finite-friction contact lane
   cross-check each other; creep, microslip, rolling resistance, spin friction,
   and loss of contact are hybrid constitutive modes rather than violations
   silently projected away. This is load-bearing for genuinely rolling
   constant-width mechanisms, rolling bearings, wheels, and traction drives.
4. **Dual equilibrium/self-stress layer.** A constraint multiplier `λ` lives in
   the dual constraint space and induces generalized reaction `J(q)ᵀλ`.
   Self-stress is the special subspace `ker(J(q)ᵀ)`: multiplier patterns in
   internal equilibrium with zero net generalized reaction. Compatibility
   obstructions live in `coker J(q)` and are paired/isomorphic to that dual
   nullspace only after the inner products and finite-dimensional duality are
   declared. Maxwell–Calladine index identities relate mechanisms and
   self-stress, but neither space is called `H¹` unless an explicit cochain
   complex and duality theorem establish that identification. Localized cycle
   representatives are diagnostic artifacts, not automatic proofs of preload.
5. **Unilateral/complementarity layer.** Contact and backlash define tangent
   cones and normal cones. Form closure asks whether
   `{ξ : Aξ=0, Cξ≥0}` contains only gauge motion **at a declared active contact
   mode**; frictional lock requires the appropriate Coulomb cone. Clearance or
   backlash generally creates an inactive-contact neighborhood, so the correct
   claim is then a bounded reachable/play set and non-escape/invariance result,
   not `{0}` form closure. Interval/conic dual witnesses can certify a
   linearized active boundary mode, while robust retention over tolerances,
   loads, impulses, friction and mode uncertainty requires global reachability
   or optimization over the full box—not corner sampling.

**Spectral health, stated correctly.** A bare `σ_{m+1}(J)` is ambiguous for a
rectangular wide Jacobian because a conventional singular-value array has only
`min(rows, cols)` entries and need not contain the domain-nullity zeros. Declare
a domain metric/scaling `G>0` and constraint metric `W>0`, then form
`J̃=W^(1/2) J G^(-1/2)`. Track the full gauge-fixed domain normal operator
`N_V=J̃ᵀJ̃=G^(-1/2)JᵀWJG^(-1/2)`; after `m` certified mechanism modes, the
health margin is `λ_{m+1}(N_V)` with eigenvalues ordered upward. More
invariantly, certify the original-coordinate nullspace `K_x=ker J`, map it to
the whitened nullspace `K_y=G^(1/2)K_x=ker J̃`, and bound
`σ_min(J̃|_{K_y^⊥})`, where `⊥` is Euclidean in whitened coordinates
(equivalently, use the `G`-orthogonal complement in original coordinates),
remembering `λ=σ²`. The dual constraint-space operator
`N_C=J̃J̃ᵀ=W^(1/2)JG^(-1)JᵀW^(1/2)` tracks self-stress multiplicity; it shares
the nonzero spectrum, with kernels mapped back through the stated scalings.
The serialized convention states ordering, weighting,
scaling, gauge, and whether null zeros were explicitly padded. All margins
carry multiplicity and interval eigenvalue bounds. “Smallest nonzero
eigenvalue” is not stable when nullity changes and does not universally
characterize a nonlinear singularity.
`fs-spectral` remains an L1 matrix/operator service; the L3 kinematics crate
assembles its operator and consumes the generic monitor, preserving layer
direction.

**Artifacts.** `ConfigurationConstraintGraph`, `HolonomyReport`,
`TangentRigidityReport { gauge, ker_j, ker_jt, cokernel, index }`,
`SecondOrderMobilityReport`, `NonholonomicDistributionReport`,
`RollingConstraintReceipt`, `FormClosureCertificate`,
`DwellRetentionCertificate`, and `SingularityMargin`. Each report carries
units/scaling, configuration and active-set identity, rank intervals, branch
assumptions, and `Unknown` bands.
Falsifiers include Grübler–Kutzbach where its genericity assumptions hold,
independent SVD/rank-revealing factorizations, Maxwell–Calladine identities,
automatic/finite-difference Jacobian checks, and validated finite continuation.

**Why the sheaf idea still earns its place.** Local joint/body charts and
interface observations can be glued by descent; incompatible local sections
and support of dual residuals can localize bad joints or tolerances; cosheaf-like
pushforward naturally aggregates loads and reactions. An ideal zero-clearance
Geneva boundary mode may earn a cone certificate; a manufactured mechanism
with crescent/slot clearance instead earns a robust dwell-play reachable set,
non-escape/invariance certificate, and boundary-mode tangent margins—not a
misnamed cohomology class or a false claim that the wheel has no local motion.
Gear ratios are relations on a configuration graph; planetary closure is
recovered from the actual constraint quotient and independently checked
against the Willis relation.

### 3.3 Certified Motion — swept volumes, envelopes, and conjugate action `[F]`

**Motor-path enclosures.** A trajectory is `M(t)`, `t ∈ [0, T]`. It is not an
unconstrained vector of Taylor polynomials: that would generally leave `SE(3)`.
`CertifiedMotorTube` stores Taylor models in a local Lie-algebra chart plus a
certified retraction, or component models with unit-motor/Study-constraint
residual bounds. It fixes the motor double-cover sign deterministically,
encloses normalization, composition, inversion, and point/AABB action, and
performs validated chart transitions before a logarithm cut. Every tube is
falsified against pointwise `fs-ga::Motor` evaluation. Composing its certified
action with a chart's Lipschitz evaluation encloses
`φ(M(t)⁻¹ x)` over a time interval—the primitive from which everything below is
built.

- **Swept volume as an implicit enclosure chart.** For a sign-correct implicit
  field, `inf_t φ(M(t)⁻¹x)` represents the union's inside/outside predicate.
  It is **not generally a signed-distance function**; even when every input is
  an exact SDF, distance accuracy must be stated separately by region and may
  fail at medial sets. `SweptChart` evaluates certified lower/upper implicit
  bounds by deterministic branch-and-bound over `t`, returning `Enclosure` or
  `Unknown`, never silently advertising exact distance. It plugs into the Rep
  Router only through conversion edges that preserve that enclosure class.
- **Envelopes and conjugate profiles.** Interior envelope candidates satisfy
  `{x : ∃t, F(x,t) = 0 ∧ ∂F/∂t(x,t) = 0}`, but tangency regularity, parameter-
  endpoint contributions, branch selection, visibility, and equality to the
  relevant swept boundary must be proved rather than assumed. For gears this
  supports the theory of conjugate action: the mating tooth profile is the envelope of the moving
  tooth in the frame of the other gear. For a **Wankel geometry, derive and
  certify separately** the ideal apex-point epitrochoid, the apex-seal center/
  contact locus, and the actual bore as the envelope of the declared finite
  seal-tip geometry and clearance model; the rotor flank is a separate
  conjugate-envelope/design problem. A Wankel rotor is Reuleaux-like
  but is not assumed to be an exact constant-width Reuleaux triangle, and its
  apex seals slide rather than roll. The analytic trochoidal construction,
  with parameters and frame convention derived and symbolically checked,
  becomes the G1 oracle; no mnemonic formula is promoted without that derivation.
  Artifact: `envelope(chart, motor_path) -> EnvelopeChart`
  with interval-certified containment (envelope ⊆ enclosure band).
- **The Law-of-Gearing certificate.** For a candidate gear pair, define the
  transmission function `i(θ₁) = −dθ₂/dθ₁` implied by contact and compare it
  with a specified target `i*(θ₁)`/instantaneous pitch curves; constant ratio is
  one specialization. For planar fixed-axis pairs, the common normal passes
  through the instantaneous pitch point, and classical Euler–Savary governs
  the applicable curvature relation of the planar relative motion/centrodes.
  Spatial/crossed-axis gearing instead certifies `n·(v_b-v_a)=0`, the common
  normal and admitted line/point-contact branch, plus screw-axis/instantaneous-
  center data, and uses family-specific conjugate-surface differential geometry
  and contact-curvature machinery—it never imports planar Euler–Savary by
  analogy. Generator singularities, trimming, addendum/
  root limits and branch visibility separately govern undercut/interference.
  Artifact: `conjugate_certificate(profile_a, profile_b, i_star, ε)`
  — an interval sweep over every active contact branch and handoff in a full
  mesh cycle proving `|i(θ)-i*(θ)|≤ε`,
  plus localized violation intervals when it fails. This is a **law-of-gearing
  ratio certificate**, not transmission error. Transmission error is the
  driven angular-position deviation from the ideal kinematic map. FrankenSim
  reports common-normal/conjugacy defect, integrated unloaded kinematic TE,
  static loaded TE including tooth/shaft/bearing/housing compliance and
  manufacturing error, and dynamic TE/NVH excitation as separate quantities.
  An interval width is uncertainty in one of those quantities, never its
  definition. Properly generated, noninterfering ideal involutes under their
  declared center-distance/contact assumptions supply the constant-ratio
  oracle; modifications and elastic load intentionally create nonzero TE.
- **Clearance and interference over motion.** Extend `fs-query::separation` to
  `separation_over(path_a, path_b, [0,T])`. A lower bound on separation alone
  can prove a positive gap but does **not** enclose the minimum. The result is a
  `ClearanceRange { lower, upper, witness_time, errors }`: branch-and-bound
  supplies the global lower bound, while a feasible time/configuration witness
  supplies the upper bound. If no admissible witness is available, the receipt
  is explicitly one-sided. Chart-conversion, spatial-discretization, motion-
  tube, and optimization errors are all inflated into the range so the claim
  applies to abstract regions, not merely a converted mesh. Backlash analysis,
  cam-follower interference, and assembly feasibility build on this primitive.
- **Overlap witness versus penetration depth.** The minimum of
  `max(φ_A,φ_B)` can expose a deepest common-interior witness/inradius for exact
  SDFs; it is not the minimum translation or SE(3) displacement needed to
  separate bodies. Translational penetration depth for convex bodies has the
  certified Minkowski-difference/GJK-EPA lane; nonconvex translational depth
  adds decomposition or interval global optimization. Minimum separating
  displacement in full SE(3) is a different pose-space optimization requiring
  an explicit translation/rotation metric and symmetry treatment. Contact response
  uses local signed gaps and normals, not a mislabeled global overlap scalar.

### 3.4 Validated Hybrid Events — proof over the true flow, not an interpolant `[F]`

Geneva engagement, gear handoff, valve motion, contact mode changes, and switch
commutation are hybrid events. Dense-output root polishing is useful but does
not prove that a true trajectory has no missed root—especially for grazing,
clustered, simultaneous, or even-multiplicity events.

**Required foundation: `ValidatedStep`.** For prescribed analytic
`MotorPath`s, interval/Taylor evaluation can directly enclose the exact path.
For simulated ODEs/DAEs, an ordinary local polynomial encloses only its
interpolant. A `ValidatedStep` must enclose the *true flow* using an interval
Picard/Lohner or Taylor-model flowpipe, or a defect-plus-Lipschitz a-posteriori
bound. Constrained DAEs additionally carry a manifold-residual tube, index and
regularity assumptions, and consistent-initialization receipt. Nonsmooth
contact time-stepping is a separate inclusion/complementarity lane rather than
being forced through a smooth ODE certificate.

Per validated step, guard range and derivative enclosures first prove exclusion
or partition the interval. Complete finite root accounting is claimed only for
a serialized guard family with an isolability theorem—e.g. polynomial/analytic
or Taylor-model guards with certified nonflatness/finite derivative order, or a
declared definable/o-minimal class on a compact phase box. Interval Newton/
Krawczyk may certify existence and uniqueness of a simple root; class-specific
root-count/exclusion subdivision then proves that all roots have been
accounted for. An arbitrary smooth/black-box guard can have a flat zero interval
or infinitely accumulated roots, so it receives `Unknown`, a set-valued
inclusion/Zeno artifact, or refusal rather than a fake finite certificate. A
zero-containing interval alone does not prove a root. Grazing guards use joint
enclosures of `(g,gdot,...)`; simultaneous events enumerate admissible reset
orders or return a set-valued post-state. Resets are revalidated, and Zeno
accumulation requires a proved limit/inclusion or a budget refusal.

**Artifacts.** `ValidatedStep`, `TrajectoryTube`, `RootCountCertificate`,
`CertifiedEvent { time, guard, pre_set, post_set }`, `NoEventCertificate`, and
`ModeLedger`. `EventLocator` may also return `PossibleEvent` or `Unknown` with
the unresolved spacetime box. A trajectory is event-capture `Verified` only
when every accepted step has true-flow coverage and complete guard accounting;
the classical dense-output lane remains a fast `Estimated` baseline and an
independent falsifier.

### 3.5 Equivariant Sliding Interfaces — sheaf compatibility plus conservative transfer `[F/M]`

**The problem it solves.** Many rotating-machine field models have a moving/
stationary discretization interface: the motor air gap, a Wankel chamber-film
interface, or a shaft/seal fluid interface. Those solvers need to glue moving
and static traces without global remeshing, conservatively, at each coupling
sample.

**The construction.** Model the air-gap interface as a cellular sheaf over the
interface circle's cell complex whose restriction maps from the rotor side are
composed with the rotation `R(θ(t))` — a *rotation-equivariant* version of the
existing chart-overlap sheaf. Then:

- **Compatibility and conservation are different obligations.** An
  equivariant sheaf/descent object can state whether rotor- and stator-side
  traces are compatible and can localize a gluing defect. Power/flux
  conservation additionally requires oriented primal/dual pairings and a
  commuting transfer diagram. If `T_rs:V_r→V_s` is the primal transfer and
  `M_r,M_s` define the two pairings, its unsigned adjoint is
  `T_rs* = M_r⁻¹ T_rsᵀ M_s` on real spaces, or
  `T_rs* = M_r⁻¹ T_rsᴴ M_s` under the declared complex sesquilinear convention;
  traction/flux orientation may add a declared sign.
  Both direction and convention are stored rather than inferred,
  constants/normal fluxes must be preserved, mortar spaces need an inf-sup
  bound, and moving meshes need a discrete geometric-conservation law. These
  identities—not vanishing `H¹` alone—certify zero artificial interface power.
  Existing watertightness code is a design pattern and falsifier, not a
  verbatim proof reuse.
- **Symmetry block-diagonalization.** On a circular interface whose mesh,
  coefficients, boundary data, and excitation preserve the declared cyclic
  action, rotation acts through a Cₙ representation (or S¹ truncated to Fourier modes):
  harmonic-basis interface coupling (the classical "air-gap element") is the
  DFT block-diagonalization that `fs-symmetry` already implements for
  circulant operators — extended with complex inter-sector phase factors
  (Bloch/Floquet boundary conditions) so a p-pole-pair motor solves one pole
  pitch instead of the full circumference **only when sources, winding layout,
  faults, skew, saturation state, and excitation respect the selected periodic
  or anti-periodic representation**. A consuming L3 domain adapter proves that
  symmetry reduction is legal before solving; L1 `fs-symmetry` only supplies
  representation kernels and residual tests. `fs-cheb` supplies the Fourier
  interface basis.
- **Mortar/Nitsche lanes.** General geometry may use dual mortar with Lagrange
  multipliers, or a separately derived multiplier-free consistent/penalty
  Nitsche coupling sharing the machinery style of `fs-cutfem`. Compatibility,
  stability/inf-sup as appropriate, transfer
  adjointness, GCL, and direct global balance are separately audited.

**Artifacts.** `SlidingInterface { static_trace, rotating_trace, θ(t) }` with
`transfer(field, t) -> (field, InterfaceReceipt)`; `SectorSymmetry` (complex
sector-coupled circulant solve, extending `fs-symmetry`); and
`InterfaceReceipt { compatibility, primal_dual_adjointness, flux, power, gcl }`.
The direct balance recomputation is disjoint from the sheaf-localization path.

### 3.6 Cohomology–Circuit Duality — windings, currents, and flux linkage `[F]`

**The observation.** Coupling a field-level machine model to a lumped circuit
(battery, inverter, load resistor) often uses cuts whose admissible independent
classes are organized by relative (co)homology. FrankenSim already computes
harmonic bases and cycle pairings (`fs-feec/src/cohomology.rs:236-284`):

- A closed filament loop can use the de Rham pairing
  `λ = N⟨A,c⟩ = N∫_S B·n dS`, provided the cycle/spanning-surface orientation,
  relative homology class, boundary conditions, and gauge independence are
  certified. A production stranded winding instead uses a distributed,
  divergence-compatible unit-current winding function `J_unit`, with
  `λ = ∫ A·J_unit dV` plus declared end-winding/leakage treatment. Both
  realizations share one neutral `WindingTerminal` protocol and cross-check
  each other where both apply.
- A massive conductor carrying imposed net current uses a formulation-specific
  cut, source field, or relative-(co)homology constraint selected from the
  topology and boundary conditions. A period of `H` is constrained only when
  Ampere's law, orientation, excitation model, and chosen generator make that
  pairing valid. `harmonic_basis` and `deflate_harmonics` provide ingredients,
  not the complete conductor model.
- Kirchhoff's laws themselves are the statement that the circuit graph's
  cycle/cutset spaces are orthogonal complements under a declared graph
  pairing and orientation convention. MNA assembly in `fs-circuit` reuses the
  same integer-incidence discipline as `fs-rep-mesh` complexes, while recording
  the pairing that turns algebraic annihilation into a physical power identity.

**Artifacts.** `WindingTopology`, `FilamentPairing`, `DistributedWinding`, and
`WindingTerminal { turns, orientation, resistance_model, leakage_model }`.
The discrete pairing is linear for a fixed mesh/winding function, so its VJP is
exact; geometry motion, nonlinear materials, circuit coupling, end effects,
and time differentiation still participate in the full adjoint and do not come
“for free.” Gauge-invariance and energy-reciprocity tests are merge gates.

### 3.7 The Port-Thermodynamic Bus — typed composition without false passivity `[S→F]`

Extend the protocol vocabulary to rotational, electrical, magnetic, thermal,
chemical, and stream ports, but distinguish four primitives:

1. `ConservativeJunction` — a Dirac/Stokes–Dirac power-conserving relation;
2. `StorageElement` — Hamiltonian/free-energy state and constitutive gradient;
3. `DissipativeRelation` — resistive/frictional/viscous/conductive/plastic law
   with monotonicity or nonnegative-production evidence;
4. `SourceOrReservoir` — prescribed environment, fuel, voltage, heat, or work
   with explicit accounting boundary.

An open `StreamPort` carries mass flow, species/element flow, momentum flux,
energy and entropy flow together. Each boundary selects exactly one energy
accounting chart. Canonical moving-stream energy is
`ṁ(h+|u|²/2+gz)` plus declared **deviatoric/extra-stress** work, heat and
diffusive terms; pressure-flow work is already inside enthalpy. A full
Cauchy-stress-work chart instead uses internal energy `e`, removes the `p/ρ`
part of `h`, and must prove the exact pressure-work crosswalk before admission.
A conjugate chart such as `(T,Ṡ),(μᵢ,ṅᵢ),(-p,V̇)` is admissible only with an exact Euler/
Legendre identity for the same state, basis and boundary kinematics; definitions
of heat/diffusion must say where partial-enthalpy and chemical-potential terms
live. Admission rejects unproved crosswalks or summing decompositions and
thereby double-counting chemical energy. Species basis, element/charge constraints and
moving-control-surface relative flux/GCL convention are part of the schema.
Reaction affinity `A_r=-Σᵢνᵢr μᵢ` pairs with reaction progress where used. Audit windows
check first-law residual, mass/element/charge balance, and
`ΔS - ∫ Sdot_boundary dt = S_gen ≥ 0`, where the signed total boundary entropy
flux explicitly accounts for advected species, diffusive/chemical-potential
terms, heat and radiation under the selected nonequilibrium convention. Only a
simple convective stream/reservoir specializes this to terms such as
`s mdot + Qdot/T_boundary`; signs and boundary temperatures remain explicit.
An optional exergy ledger reports
`T0 S_gen` and availability destruction relative to an explicit environment.
Individual oriented port contributions need not be nonnegative.

Every machine remains a declarative graph of these typed relations. A graph
admission pass checks dimensions, orientation, causality/DAE index, conserved
quantities, source closure, and audit completeness. Strongly coupled physics
may assemble a monolithic residual/Jacobian through public operators; “no
private state reach-through” remains the rule, while “partitioned Dirac only”
is not.

### 3.8 Entropy-Stable Compressible Flux Ledger `[F]`

For the 1-D/2-D compressible gas lane (§5.8), use an entropy-conservative
spatial flux derived for the **exact admitted EOS/mixture** (classical
Tadmor/Ismail–Roe/Chandrashekar constructions are not universal formulas) plus
explicit dissipation, then
pair them with a fully discrete entropy-stable time/source/boundary update. A
semidiscrete spatial proof alone earns only a semidiscrete claim. The complete
lane targets a **per-step discrete entropy inequality** checked by a dedicated
source/boundary-aware ledger. A shock-capturing scheme whose mathematical-
entropy inequality satisfies its declared signed inequality by construction is
audited with boundary, moving-area, heat, friction, and reaction source terms and
earns a **conditional discrete-stability claim**; it does not by itself validate
the EOS, turbulence, combustion, or shock prediction. Positivity/invariant-
domain preservation covers density, internal energy and species. Reacting
lanes require an EOS-specific convex entropy pair, pressure-equilibrium
preservation where applicable, thermodynamically consistent chemistry, and
source integration satisfying the declared entropy/free-energy law. Moving
meshes add metric identities/GCL and nonlinear-solve tolerance to the per-step
inequality. Phase-equilibrium boundaries can lose smooth convexity and enter a
separate free-energy/inclusion lane rather than inheriting a gas proof. Baseline
MUSCL–HLLC stays a separately labeled engineering lane unless its exact variant
is proved to satisfy both properties.

### 3.9 Paired Sheaf–Cosheaf Machine Balance `[F/M]`

Compatibility and conservation travel in opposite directions. Local state and
trace observations restrict from assemblies to parts—a sheaf-like operation.
Loads, currents, species production, heat and momentum contributions aggregate
from parts to assemblies—a cosheaf-like pushforward. FrankenSim will make that
pairing executable rather than calling every graph a sheaf:

- `CompatibilityComplex` owns local observations, restriction maps, descent
  residuals and obstruction support;
- `BalanceCocomplex` owns oriented extensive quantities and aggregation maps;
- `PowerPairing` proves discrete virtual-work/energy reciprocity between them;
- `RelativeBoundaryReceipt` records what crossed the accounting boundary, so a
  local conservation claim cannot accidentally become a closed-system claim.

The G0 artifact is a commuting diagram between restriction/aggregation and the
actual discrete trace/transfer matrices, plus a pairing identity under the
declared mass metrics. The boring baseline is direct global assembly. This
moonshot activates only if it improves localization, incremental recomputation
or distributed composition while matching that baseline; otherwise it remains
an explanatory layer, not production authority.

### 3.10 Proof-Carrying Adaptive Fidelity and the Living Machine Twin `[F/M]`

A world-class system should not force one fidelity everywhere. Represent each
subsystem implementation—lumped, 1-D, 2-D, reduced-order, 3-D field, empirical
correlation—as an object with typed ports, validity domain, cost/error model,
state projection/lifting maps and falsifiers. A `ModelCrosswalk` records whether
two representations commute on conserved quantities and QoIs; discrepancy is a
measured naturality defect, not hand-waved “calibration.”

`FidelitySheaf` assigns locally admissible models over operating-regime and
spacetime cells; section selection is a constrained planning problem under
accuracy/time/memory budgets. An obstruction means neighboring choices cannot
exchange compatible state or close the error budget, triggering refinement,
model escalation or refusal. Existing surrogate, regime, evidence, plan,
assimilation and as-built crates supply the pieces.

The runtime artifact is a `ProofCarryingMachineTwin`: nominal/as-built geometry,
material/process posteriors, controller version, sensor history, model graph,
receipts and unresolved discrepancies. New measurements update only downstream
claims through self-adjusting provenance; OED selects the next coupon, sensor
or dynamometer test with maximum expected decision value. The baseline is a
fixed-fidelity replay. Activation requires lower cost at equal decision error
on held-out transients; kill criteria forbid a router that thrashes, launders
evidence, or cannot reproduce the fixed-fidelity result.

### 3.11 The FrankenSim Theorem Foundry — invent the missing mathematics `[M]`

The conditional language above is a truth boundary for *today*, not an ambition
boundary for tomorrow. FrankenSim should deliberately attempt theorem artifacts
intended to exceed what any single commercial workflow exposes; roadmap phase
E7 must test that competitive workflow claim with a versioned survey, while
phase E8 tests the
mathematical claims rather than assuming either. The research unit is a
`TheoremCard { statement_id, statement_version, exact_quantifiers, assumptions,
imported_lemmas, statement_state, mathematics_state, formalization_state,
falsification_state, implementation_state, review_state, trusted_computing_base,
kernel_version, proof_hash, unresolved_gaps, formal_kernel, executable_checker,
counterexample_search, downstream_claims }`.
Those states are **orthogonal**: a theorem can be proved but not implemented at
scale, an executable checker can pass without constituting a proof, and a scale
failure does not refute a theorem outside the failed implementation claim.
The axes are statement lifecycle (`Draft|Stable|Superseded|Withdrawn`),
mathematics (`Conjectured|ConditionalProof|Proved|Refuted`), formalization
(`Unformalized|PartiallyFormalized|KernelChecked`), falsification
(`Unchallenged|CampaignRunning|FalsifierSurvived|CounterexampleFound`),
implementation (`Absent|ConformanceChecked|ScaleQualified`), and review
(`Unreviewed|IndependentAuditPassed|Reproduced`). `CounterexampleFound`
forces the exact mathematics statement to `Refuted`; surviving a finite
falsification campaign never forces it to `Proved`. A Rust checker establishes
implementation conformance, not automatically a formal proof. A counterexample marks the exact statement
version `Refuted`, invalidates every downstream claim that depended on it, and
produces a minimal failure artifact plus a separately versioned sharper
successor—it never silently edits the old theorem into something true.

The initial theorem portfolio is:

1. **Relative-cohomology electromechanical power theorem.** Construct a
   representative-independent pairing between field cochains, distributed or
   filamentary winding classes, and circuit terminals such that gauge changes,
   tree/cotree choices, and homologous cut changes leave flux linkage invariant,
   while the discrete field/circuit work identity closes exactly. Extend the
   theorem to moving conductors by transporting relative classes with a discrete
   Lie/ALE map whose motional-EMF term and GCL commute; topology-changing motion
   uses explicit maps, spans, or cobordisms between pre/post relative complexes
   (with zigzag-persistence diagnostics where useful) and becomes a certified
   hybrid event rather than an unspoken identification of non-isomorphic
   classes.
   Artifact: `ElectromechanicalPairingTheorem` plus a checker over arbitrary
   admissible basis/cut representatives.
2. **Equivariant no-spurious-power gluing theorem.** Prove sufficient—and seek
   near-necessary—conditions under which a paired sheaf/cosheaf transfer on a
   moving nonmatching interface simultaneously preserves constants, oriented
   flux, virtual work, discrete energy, and GCL. The theorem should cover cyclic
   harmonic air-gap elements as a corollary and quantify the exact defect when
   symmetry, quadrature, or mortar stability is broken. Artifact:
   `EquivariantPowerGluingTheorem` whose generated runtime receipt is §3.5's
   `InterfaceReceipt`.
3. **Guaranteed electromagnetic force/torque enclosure theorem.** Combine
   compatible FEEC spaces, equilibrated functional majorants, interval geometry
   and constitutive data, and shape-functional duality to enclose Maxwell-stress
   and virtual-work force/torque—including discretization, quadrature, geometry,
   material/stress convention, gauge, boundary conditions, held electrical
   variables, torque origin/frame, and interface-transfer contributions. Where the two methods'
   assumptions overlap, prove that intersecting their enclosures is sound and
   strictly sharper on a characterized class. Artifact: `ForceFunctionalBound`
   with a machine-checkable dependency closure, not merely method agreement.
4. **Nonlinear mechanism descent and singularity theorem.** Build the local
   SE(3) solution groupoid/stack and prove when its tangent cohomology, dual
   complex, and higher obstruction maps recover finite mobility, self-stress,
   branch multiplicity, and second-order rigidity. Extend the construction to
   nonholonomic distributions by identifying what curvature/bracket-growth and
   groupoid holonomy replace position-level closure for ideal rolling. Regular
   mechanisms should emerge as a corollary; Bennett/Bricard, rolling disks and
   constant-width bodies, and unilateral Geneva dwell modes are the adversarial
   cases. Artifact: `NonlinearMobilityTheorem` connecting finite continuation,
   tangent/second-order certificates and nonholonomic accessibility without
   conflating them.
5. **Hybrid event-completeness theorem for constrained machines.** Prove that a
   finite cover of validated regular/index-fixed DAE flow tubes on a compact
   phase domain plus an isolable finite guard family, reset closure, and root-
   count obligations yields a complete mode ledger on a bounded horizon,
   assuming no accumulation or a proved Zeno limit, including a
   set-valued statement for simultaneous/grazing events and a typed alternative
   for unresolved flat/accumulating guards. Moreau inclusions are a separate
   extension. Artifact: `HybridCoverageTheorem` consumed by contact,
   valves, commutation, backlash, and ignition.
6. **Whole-machine discrete thermodynamics theorem.** Find a composition and
   time-discretization class in which Dirac interconnections, storage gradients,
   memoryless maximal-monotone dissipators, and stateful hysteretic/reactive/
   frictional internal-variable laws with explicit free-energy/dissipation
   updates, stream ports, and asynchronous
   multirate windows yield an exact discrete first-law balance and a signed
   second-law inequality for the assembled open machine—not just each component.
   Quantify the defect and restore the theorem under adaptive window refinement
   when exact asynchronous closure is impossible. Artifact:
   `ComposedThermodynamicsTheorem` plus exergy and boundary ledgers.
7. **Proof-preserving fidelity-descent theorem.** Treat model crosswalks as
   morphisms carrying conserved state, QoIs, validity predicates, and evidence.
   Prove a global decision-error bound from local naturality defects and show
   that refinement/escalation is monotone in an evidence preorder: a router may
   spend more or less compute but cannot silently strengthen a claim. Artifact:
   `FidelityDescentTheorem`, the formal spine of the proof-carrying machine twin.
8. **Certified conjugate-machine geometry theorem.** Unify gear generation,
   cams, pumps, Wankel seal/housing families, and general envelope-of-motion
   design in one interval-verified planar and spatial surface-contact envelope
   calculus. Establish when the
   stationary envelope equations plus branch/visibility/contact conditions are
   sufficient for conjugacy, when offsets preserve it, and how manufacturing
   and elastic perturbations bound transmission or chamber-closure defect.
   Artifact: `ConjugateGeometryTheorem` with counterexample-producing failure
   localization.
9. **Thermodynamically closed contact–tribology limit theorem.** Characterize
   when barrier, compliant and nonsmooth contact sequences converge without
   spurious energy while retaining frictional heat, wear/internal state,
   Painlevé/set-valued behavior and generalized sensitivity. Artifact:
   `ContactTribologyLimitTheorem` plus an adversarial multicontact balance
   generator.
10. **Reactive moving-mesh entropy/positivity theorem.** Unify variable
    thermodynamics, detailed-balance chemistry, invariant domains, source/
    boundary entropy, nonlinear-solve defect and ALE GCL in a fully discrete
    inequality, with a separate free-energy inclusion at phase boundaries.
    Artifact: `ReactiveAleEntropyTheorem`.
11. **Multisymplectic spacetime composition theorem.** Seek a discrete
    variational/multisymplectic complex for moving-domain field–structure–
    circuit subsystems whose boundary forms coincide with typed ports and whose
    spacetime gluing yields computable local/global conservation-law defects.
    Artifact: `MultisymplecticCompositionTheorem`; this gives the E8 moonshot a
    real statement, owner and checker.
12. **Closed-loop proof/UQ theorem.** Compose validated hybrid reachability,
    controller implementation error, model/fidelity bounds and anytime-valid
    rare-event evidence into a risk bound without confusing sampling validity
    with truth of the physical population model. Artifact:
    `ClosedLoopRiskCompositionTheorem`.
13. **Acoustic radiation/source-transfer theorem.** Prove a class of
    structure/flow/source-surface transfers that preserves acoustic power,
    radiation, reciprocity where applicable and far-field QoIs under bounded
    dispersion/source defects. Artifact: `AcousticTransferTheorem`.
14. **Failure-mode descent theorem.** Develop dependency-aware composition of
    competing degradation modes, common causes, inspections and repairs without
    false independence, and prove when local failure evidence descends to a
    system risk bound. Artifact: `FailureModeDescentTheorem`.
15. **Thermodynamic constitutive-learning theorem.** Build a guarded learned-law
    class whose architecture and checker enforce objectivity, material symmetry,
    free-energy/dissipation admissibility, state-update consistency and
    extrapolation refusal, with approximation and discretization errors kept
    separate. Artifact: `ThermodynamicConstitutiveLearningTheorem`.

Every card has a boring finite-dimensional baseline, an adversarial
counterexample generator, a paper derivation, a small formal proof kernel where
feasible, an executable Rust checker, and a scale-up experiment. These are not
decorative future claims: cards start beside their enabling work—E0 whole-
machine thermodynamics/identity, E1 geometry/mechanism descent, E2 hybrid and
contact, E3 structural acoustics/NVH, E4 electromechanical/closed-loop power,
and E5 aero-/combustion-acoustic and reactive extensions.
E8 is the independent reproduction, kernel-checking and at-scale integration
summit, not the start date. §12 records their
proof obligations, and a proved card may promote an `[F]` lane only after its
checker survives independent falsification. The goal is not merely to apply
known theory better; it is to turn new mathematics into replayable machine
engineering authority. Any proof-assistant or symbolic-algebra toolchain is an
isolated dev/research oracle; the production graph consumes only versioned proof
artifacts and FrankenSim-native checkers, preserving the dependency policy.

---

## 4. New Crate Atlas and Existing-Crate Upgrades

### 4.1 New crates

Layer assignments follow the existing L0–L6 discipline. `xtask check-layers`
checks allowed layer pairs but does not by itself rule out same-layer cycles;
the explicit DAG below must also resolve under real Cargo metadata. The table
lists domain-significant edges—shared utility/`fs-exec` edges remain mandatory
where their contracts require them. Tags are the *initial* tag; nothing ships
to a default path before its Gauntlet tier (§9).

Within L3, `fs-couple` remains a dependency-light protocol/composition base: it
may depend on the existing dependency-free `fs-iface` static interface checker,
but never on a domain solver. Neutral ports, terminals, clocks, opaque interface-
state handles, transfer *traits*, iteration histories, and balance receipts live
there. Concrete mortar/Nitsche/harmonic operators require function spaces,
traces, meshes, and quadrature, so domain adapters assemble them and implement
the neutral traits; they are not smuggled into the dependency-light crate.
Immutable physical observations, claims, cards, validity predicates, and
initial-state distributions live in `fs-matdb`; executable laws and mutable
runtime state do not. EM, circuit, thermal, control, contact, tribology, and
machine crates depend inward on those protocols. The
same-layer graph is an explicit DAG, not permission for mutual dependencies.
`fs-iface` continues to own function-space/trace descriptors and static pairing
admission; `fs-couple` unambiguously owns power/stream/control/WindingTerminal
protocols, window/clock algorithms, IQN histories, and generic transfer
admission. L3 `fs-material` owns the executable `ConstitutiveGraph`, law-node
interface, and generic typed runtime-state protocols. Immutable `LawId`,
`LawVersion`, parameter blocks, state-schema descriptors/versions and initial-
state policy metadata live in L1 `fs-matdb`; L3 executors consume them rather
than redefining them. L3 `fs-tribo`, `fs-em`,
`fs-thermal`, `fs-gas`, `fs-power`, and other domain adapters depend one-way on
`fs-material` and implement specialized nodes. L1 `fs-thermochem` instead
exports lower-layer thermodynamic/kinetic closures and immutable law data; the
L3 `fs-gas` adapter wraps those closures as graph nodes, so L1 never imports an
L3 protocol.
Likewise, `fs-gas` owns thermochemical/control-volume closures and 0-D/1-D
transport, while existing `fs-flux` owns multi-D spatial/ALE/CutFEM transport
and consumes those closures; `fs-gas` never depends back on `fs-flux`.

| Crate | Layer | Tag | Dependency/ownership rule | Scope (one line) |
|---|---|---|---|---|
| `fs-motion` | L2 MORPH | [S→F] | fs-ga, fs-geom, fs-query, fs-ivl | Certified motor paths, immutable `ChartAtTime`, `SpacetimeChart`, swept/envelope implicit enclosures, separation-over-time; never pretends a timeless `Chart` is moving |
| `fs-kinematics` | L3 FLUX | [S→F] | fs-ga, fs-motion, fs-ivl, fs-solver, generic fs-spectral operators | Joint specifications, nonlinear configuration graph/holonomy, tangent rigidity and dual self-stress reports, solver-backed continuation, gear/cam/Geneva geometry; contains no elasticity or dynamics |
| `fs-contact` | L3 FLUX | [F] | fs-motion, fs-query, fs-ivl, fs-solver, fs-tribo, fs-matdb, fs-couple, fs-iface; **not** fs-solid or fs-mbd | Capability-routed CCD/candidate management over moving bodies, embedding tribological constitutive responses into gap/cone/barrier/complementarity residuals, and response receipts; consumes static geometry from `fs-query` and law evaluation from `fs-tribo`, while generic conic/nonlinear math stays in `fs-solver`; fs-solid and fs-mbd consume adapters |
| `fs-mbd` | L3 FLUX | [F] | fs-kinematics, fs-contact, fs-ga, fs-time, fs-solver, fs-spectral, fs-query, fs-matdb | Rigid and flexible/reduced multibody dynamics over abstract compliance/force operators, constrained SE(3), holonomic and nonholonomic/nonsmooth lanes, reactions and rotordynamics hooks; no hard edge to fs-solid/fs-machine |
| `fs-machine` | L3 FLUX | [F] | fs-kinematics, fs-mbd, fs-contact, fs-tribo, fs-material, fs-symmetry, fs-couple | Machine-element schemas/residuals and typed excitation/load-spectrum artifacts for gears, bearings, shafts, seals, cams and valves; **not** cross-solver orchestration, which belongs to the L6 machine graph. Split a family into its own crate when it needs an independent state model, solver, or two downstream consumers |
| `fs-em` | L3 FLUX | [F] | fs-feec, fs-cutfem, fs-motion, fs-time, fs-solver, fs-matdb, fs-material, fs-symmetry, fs-couple, fs-iface | 2-D/2.5-D/3-D nonlinear magnetics, moving-conductor MQS, gauging/topology, eddy currents, forces/torques, demagnetization and loss maps |
| `fs-circuit` | L3 FLUX | [S→F] | fs-qty, fs-time + fs-couple port/terminal protocol; **not** fs-em | Descriptor MNA, structural index/consistent initialization, admitted impulse-free state continuity plus distributional impulse/refusal receipts for inconsistent topology changes, complementarity and regularized device lanes |
| `fs-power` | L3 FLUX | [F] | fs-circuit, fs-matdb, fs-material, fs-couple | Diode/MOSFET/IGBT/SiC device laws/cards, parasitics, dead time, reverse recovery, switching/conduction loss, averaged and fully switched converters/PWM; concrete thermal/control systems couple above the crate to avoid cycles |
| `fs-control` | L3 FLUX | [S→F] | fs-time, fs-solver, fs-qty + fs-couple clock/signal protocols | Runtime sampled-data/hybrid control, Clarke/Park, FOC/MTPA/field weakening, observers/PLL, anti-windup, delays/quantization/faults, and solver-neutral policy artifacts; existing L4 fs-opt/fs-sos/fs-robust own offline MPC/H-infinity synthesis, SOS/reachability, and robust tuning |
| `fs-thermal` | L3 FLUX | [S→F] | fs-feec, fs-solver, fs-matdb, fs-material, fs-motion, fs-query, fs-xform, fs-couple | Storage + conduction, nonlinear radiation/radiosity, phase change, moving/contact/gap heat, and thermal-side coupling operators; concrete CHT, thermoelastic, and electro/magnetothermal assemblies live above the participating solvers |
| `fs-thermochem` | L1 BEDROCK | [S→F] | fs-qty, fs-matdb, fs-evidence, fs-math, fs-la, fs-ivl | Species/elements/charge, NASA-9 and other standard-state thermo, mixtures, EOS/phase equilibrium, kinetics algebra and conservation checks—no transport solver |
| `fs-gas` | L3 FLUX | [F→M] | fs-thermochem, fs-material, fs-time, fs-ivl, fs-solver + fs-couple stream protocols | Conservative reactive control volumes and 1-D transport, entropy/invariant-domain lanes, valves/manifolds, spray/evaporation and combustion closures consumed by multi-D `fs-flux` adapters |
| `fs-matdb` | L1 BEDROCK | [S] | fs-qty, fs-evidence, fs-blake3 | Immutable cross-physics observations, claims, bulk/interface/process cards, neutral law IDs/versions/parameter blocks/state-schema and initial-state metadata, validity predicates, uncertainty/covariance, provenance/license and query receipts; no mutable per-run state or executable closures |
| `fs-tribo` | L3 FLUX | [F→M] | fs-matdb, fs-material, fs-solver + fs-couple interface-state protocol | Friction and thermal partition, mixed/boundary/hydrodynamic lubrication, EHL, cavitation, wear/scuffing/pitting; does not depend on contact solver internals |
| `fs-acoustics` | L3 FLUX | [F] | fs-feec, fs-flux, fs-solver, fs-spectral, fs-motion, fs-bem, fs-matdb, fs-material, fs-iface, fs-couple | Interior/exterior and moving-medium acoustics, structural/flow coupling adapters, aeroacoustic source lanes, and physical/order/psychoacoustic QoIs; L6 fs-report/flagships render reports |

Machine UQ, tolerance allocation, calibration, assimilation, optimization, and
campaign state reuse existing `fs-uq`, `fs-toleralloc`, `fs-asbuilt`,
`fs-assimilate`, `fs-oed-e2e`, `fs-robust`, `fs-eproc`, and `fs-ascent` rather
than creating an overlapping `fs-mstate` silo. Fatigue, fracture, creep,
viscoelasticity, damage, corrosion, and coupled active materials extend
`fs-material`/`fs-solid` first; a separate life crate is justified only after a
stable shared state-transition contract emerges.

This atlas proposes 15 focused crates and deliberate extensions of existing
ones. Names stay in the flat `fs-*` namespace; each ships
`CONTRACT.md` with the ten required sections before becoming a dependency
target, per AGENTS.md.

E0 also generates and compiles a machine-readable `ProposedManifestFixture`
containing **every** new and changed Cargo edge plus a same-layer topological
order. The fixture covers the root workspace and the standalone `fs-wasm`
workspace; slash notation and “shared dependencies omitted” are forbidden in
that executable artifact. `cargo metadata`/layer checks, duplicate-type checks,
and a minimal compile target fail on any same-layer cycle or undeclared edge.

### 4.2 Upgrades to existing crates

| Crate | Upgrade | Kills |
|---|---|---|
| `fs-rep-mesh` | Add an oriented `TriComplex2`/general cell-complex substrate with boundary/trace maps, stable feature IDs, topological-vs-embedding dimension, planar/axisymmetric metric metadata, and exact `d₁d₀=0` batteries; a surface half-edge mesh alone is not a 2-D FEEC complex | B5, B12 |
| `fs-feec` | `WeightedStar`/`WeightedMass`; explicit primal/dual/twisted-form and trace orientation contracts over the new L2 complex; commuting projection, coefficient-robustness, gauge, compactness/inf-sup, and quadrature proof lanes—not only `dd=0` | B5, B12 |
| `fs-solver` | Own lower-layer Newton–Krylov globalization (then `fs-ascent` may reuse it downward; never depend upward); block operators and Schur preconditioners; SPD HX/AMS-style auxiliary-space lanes, MINRES only for verified symmetric-indefinite systems, and FGMRES for nonsymmetric or variable-preconditioner systems; consume/re-export rather than duplicate L1 spectral eigensolvers | B5, B10 |
| `fs-time` | Sparse/matrix-free integrators; SE(3) Lie-group and variational lanes; true-flow `ValidatedStep`/`TrajectoryTube`, event/root-count ledger, consistent DAE initialization; implement dependency-light `fs-ad` step/reset VJP and saltation/set-valued record protocols beside the primal, with no fs-time↔fs-adjoint edge | B2, B10 |
| `fs-couple` | Version the scalar seed into `PortSchema` (stable ID, effort/flow dimensions, value/field shape, basis/frame/orientation, power pairing, clock/timestamp and conservation roles); add storage/dissipative/source/stream relations, generic transfer traits, vector IQN-ILS, monolithic escalation and multirate waveform relaxation. Domain adapters assemble mortar/Nitsche/harmonic operators; interpolation, iteration, stability and solution errors remain separate | B9 |
| `fs-qty` | Unit aliases/parser tokens plus **amount of substance as a sixth base dimension now**; version and migrate `Dims`, `QtyAny`, canonical IR, scenario/package/material schemas; keep candela out until photometry is real | B8 (units half) |
| `fs-symmetry` | Remain generic L1 group/representation/block-diagonalization machinery; add complex sector phases, dihedral Dₙ, and representation-residual tests. L3 consumers—not this crate—prove their materials, sources, BCs, faults, and excitation admit a reduction | §3.5 |
| `fs-spectral` | Own the L1 generic resumable sparse Lanczos/LOBPCG/eigen-gap service over `fs-la`. L2/L3 consumers assemble operators and push dimensionless samples; no upward type dependency or duplicate L3 eigensolver | B11 |
| `fs-scenario` | Add versioned `Scalar`, `Vector`, `Tensor`, `ComplexPhasor`, `SpeciesBundle`, `CharacteristicState`, `FieldTraceRef`, and `PortRef` payloads with units/basis/frame/orientation, histories/tables/distributions and continuity/reset semantics. Add Physics/expectation rows; retain the L3 `FrameTree` schema/evaluation and add a one-way adapter that lowers its motion specifications into L2 `fs-motion` paths/tubes—`fs-motion` never imports scenario types. Joints, terminals, controllers and resets remain machine-graph relations, not fake `BcKind`s. L6 `fs-ir` lowers one-way into versioned scenario/domain artifacts and preserves the existing canonical identity via explicit crosswalks | B1, B5–B7 |
| `fs-flux` | Promote the live exactly-divergence-free 2-D path to 3-D/high-order curved and CutFEM boundaries; own multi-D incompressible/compressible spatial discretization and consume `fs-gas` closures without a reverse edge; ALE with discrete GCL/conservative remap; monolithic/partitioned FSI; DNS resolution/convergence lanes, LES/RANS closure model cards, cavitation and multiphase ladders | P7 fluid side |
| `fs-lbm` | Curved moving boundaries, Galilean-invariance-corrected momentum exchange with measured residual, fresh-cell mass/momentum receipts, under-resolved lubrication correction, MRT/cumulant stability, thermal/multiphase/cavitation ladders and moving-boundary benchmark suite | B1 fluid side |
| `fs-material` | Consume L1 `fs-matdb` law IDs/versions, canonical parameter blocks, state-schema descriptors and initial-state metadata; own the generic executable `ConstitutiveGraph`, implementation registry, aggregate runtime-state codec/protocol and objectivity/free-energy/dissipation/locality behavior; extend laws to nonlocal/gradient/phase-field behavior, viscoelasticity, mixed hardening, damage/fracture/fatigue/creep/corrosion, magneto/piezo/thermoelectric coupling and guarded learned laws | B8, P7 |
| `fs-ad` | Extend the existing L1 solver-agnostic AD/revolve substrate with minimal local-step/reset/saltation record traits and checkpoint schedules; it never imports time/domain types | D8 |
| `fs-adjoint` | Consume `fs-ad` protocols to own L3 PDE/solver/time-composition adjoints, shape/IFT assembly, and gradient evidence gates without depending on `fs-time` or higher domain crates; local derivatives remain beside primal operators. Nonsmooth lanes use generalized/set-valued sensitivities or honest no-claim | D8 |
| `fs-bem` | Keep the live Laplace screening lane intact; add a separately gated Helmholtz/convected-Helmholtz exterior-acoustics lane with radiation condition, combined-field resonance control, singular/hypersingular quadrature, FMM/operator validation, and frequency-domain QoI receipts | P7 acoustics |
| `fs-query` | Keep `min max(φ_A,φ_B)` as an overlap-inradius witness; add capability-checked `ConvexSupportMap`, `FeatureComplex+BVH`, `ImplicitGapOracle`, `CodimThickness`, `DeformationMap`, conservative support, penetration/gap/manifold queries, and `GeometricMoments` with certified volume/first/second moments, transform covariance, parallel-axis laws and watertightness preconditions. Rep Router conversion/motion errors inflate every contact bound | B3, B2 mass properties |
| `fs-xform` | Receive a canonical L2 `LevelSetGrid` (or hoist `GridSdf`) together with WENO5 advection, FIM redistancing and velocity extension from `fs-topols`; preserve compatibility re-exports, canonical bytes and goldens. Add conservative remap and split/merge lineage receipts for wear, fracture and phase fronts | B6, B1 fluid side |
| `fs-mesh` | Add dynamic hp/anisotropic goal-oriented adaptivity, moving-mesh quality/untangling, conservative remap, topology-change lineage and contact/wear/fracture refinement; error indicators are tied to declared QoIs and conversion error | P7 high fidelity |
| `fs-solid` | Add damage/fracture/cohesive/XFEM-or-phase-field lanes, cyclic plasticity, creep-fatigue, residual stress, bolted/welded/adhesive/interference-joint adapters and remesh/state-transfer receipts; life prediction remains model- and data-bounded | P7 life/manufacturing |
| `fs-uq` + `fs-toleralloc` | Upgrade beyond fixed MLMC/regression PCE and independent first-order tolerances: adaptive multi-index/MLMC, global sensitivity, rare events, correlated/hierarchical/nonlinear and mode-switching GD&T propagation, dependency-aware tolerance allocation and decision-error receipts | P7 reliability |
| `fs-assimilate` + `fs-asbuilt` | Upgrade beyond dense linear-Gaussian and 2-D known-correspondence seeds: sparse nonlinear/non-Gaussian inference, observability/identifiability and discrepancy, authenticated calibration/validation splits, 3-D/CT correspondence-uncertain registration, posterior-predictive checks and lineage | P7 living twin |
| `fs-ir` | Add a versioned L6 machine graph in E0: stable body/surface/contact-feature/port/state-slot IDs; subsystem/model relations; typed terminals/clocks; materials/interfaces; IC/BC, motion/event/reset and tolerance/correlation semantics; sensors/experiments, hazards/faults, ContextOfUse/QoIs/acceptance; accounting and fidelity/escalation policy. Lower one-way into scenario/domain artifacts with stable-ID/hash round trips | P6, HELM |
| `fs-session` | Add capability/budget/idempotency/cancellation admission for long hybrid and theorem-checking jobs, with resumable model/version identity | HELM |
| `fs-ledger` + `fs-package` | Version and migrate material claims/query receipts, mode/event ledgers, interface/balance receipts, model crosswalks, theorem cards, controller versions, V&V dependency closure, semantic lineage and runtime checkpoints. A checkpoint binds `(LawId, LawVersion, StateSchemaVersion, canonical parameters, state bytes, contract/code hash)`; replay refuses unknown semantics instead of silently decoding | HELM |
| `fs-report` | Render independent verification, numerical-uncertainty, data/parameter, model-validation, prediction-domain, performance, and theorem-status axes from receipts rather than hand-authored prose | HELM |

---

## 5. Domain Plans

Each domain plan states: formulation, discretization, solvers, certificates,
Gauntlet battery, evidence-color policy, and no-claim boundaries. The format is
deliberately parallel to how the existing contracts read.

### 5.1 Motion and Swept Geometry — `fs-motion` `[S→F]`

**Formulation.** `MotorPath` is the point-evaluation view of a
`CertifiedMotorTube`: a piecewise validated local-Lie-chart map
`[0,T] → SE(3)` built over `fs-ga` motors, with the constraints and chart
transitions in §3.3. Analytic constructors include constant-twist screws and
the Wankel rotor **pose** (eccentric-center orbit plus rotor phasing); the
epitrochoid is the derived apex locus, not the body motion. A builder trait lets
*higher* layers lower motions into paths (`fs-scenario` frame trees, `fs-mbd`
trajectories)—the dependency direction stays downward. A moving object does
not implement the timeless `Chart` contract. `SpacetimeChart` exposes
`snapshot(t) -> ChartAtTime` and certified `eval_over(x,[t₀,t₁])`; immutable
snapshots may implement `Chart` with their time and path provenance frozen.

**Key operations.**
- `swept(chart, path) -> SweptChart` — certified min-over-time implicit
  enclosure (§3.3), not an unconditional SDF;
  registers Rep Router conversion edges (swept → mesh via existing dual
  contouring; swept → voxel) with `Enclosure` certificates.
- `envelope(chart, path) -> EnvelopeChart` — the characteristic set
  `(F, ∂F/∂t)=0` is traced by parameterized slices or validated implicit-
  manifold continuation, then rank singularities, parameter endpoints,
  trimming, visibility, self-intersection, and swept-boundary membership are
  classified before a branch is admitted; G1 oracle: closed-form Wankel epitrochoid
  and involute rack conjugates.
- `separation_over(a: SpacetimeChart, b: SpacetimeChart, span) -> Evidence<Separation>`
  — branch-and-bound in time over the existing certified static bound.
- Volume-of-motion invariants: `V(θ)` chamber-volume functions for closed
  chambers formed between moving charts (Wankel chambers, cylinder above
  piston), computed by certified quadrature of the region bounded by named
  charts. Closed forms are G1 oracles only for declared ideal nominal
  geometries/seal conventions; finite seals, clearances, port cutouts,
  deformation, and as-built geometry use certified quadrature/enclosures.

**Determinism/cancellation.** All branch-and-bound loops are tile programs
under `Cx` with deterministic work ordering (queue keyed by interval bounds,
ties by lexicographic time interval), resumable via serialized frontier.

**No-claims (initial).** No deformable sweeps (rigid paths only); envelope
`Unknown` on tangency degeneracies; no self-intersecting sweep resolution
(reported, not repaired).

### 5.2 Multibody Dynamics — `fs-mbd` `[F]`

**Formulation.** Bodies carry a normalized PGA motor `M_i` representing a pose
in `SE(3)`, with the double-cover equivalence `M_i∼-M_i` and deterministic sign/
chart convention, plus a twist `ξ_i∈se(3)` (PGA bivector).
`fs-query::GeometricMoments` owns certified volume/first/second
moments for a chart without knowing materials; `fs-mbd` combines those with a
uniform `fs-matdb` density receipt into the 6×6 spatial inertia. Spatially
varying density later uses a neutral weighted-integration operator rather than
making geometry depend on materials. Certified mass, COM, and inertia with
interval bounds are small, high-value artifacts by themselves. Constraints:
`Φ(M) = 0` holonomic (joints, gear ratios, cam profiles as driven constraints),
Pfaffian `A(M)ξ = 0` nonholonomic (ideal rolling/no-slip/knife-edge), and
unilateral (delegated to `fs-contact`). A velocity constraint is never
integrated into a fictitious configuration constraint without a Frobenius proof.

**Integrators (three lanes, all resumable + deterministic).**
- *Structure-preserving lane [F]:* Lie-group variational integrator
  (discrete Euler–Poincaré on SE(3)) with constraints via RATTLE-style
  discrete null-space/projection — symplectic/momentum behavior under stated
  conservative symmetry assumptions, nonlinear-solve-bounded constraint
  satisfaction, and a backward-error/modified-energy bound only for a named
  smooth or analytic conservative class with regular constraints, admissible
  fixed step, solve tolerance, and horizon. Impacts, adaptive steps,
  dissipation, and general chaotic horizons receive measured balance receipts,
  not that theorem. The discrete adjoint is derived from the actual
  discrete residual and verified; variational structure is leverage, not a
  substitute for that derivation.
- *Engineering lane [S]:* Lie-group generalized-α (Brüls–Cardona–Arnold) on
  index-3 DAE with GGL/stabilized-index-2 option and controllable numerical
  dissipation — the industrial workhorse for stiff contact-rich problems.
- *Nonholonomic lane [F]:* a discrete Lagrange–d'Alembert/nonholonomic
  integrator enforces `A(q)ξ=0` without claiming ordinary RATTLE, Hamilton's
  principle, or canonical symplecticity. An ideal reaction zero-work receipt is
  paired with a friction-cone feasibility test; exhausting available traction
  triggers a ledgered no-slip→creep/slip/contact transition rather than
  projection back onto an infeasible distribution.
- Constraint drift is monitored, not hoped away. Interval evaluation of
  `‖Φ(q_h)‖∞` is a `ConstraintResidualReceipt`, not a bound on distance to the
  true constraint manifold. `ConstraintStateEnclosure` additionally requires
  a validated nonlinear solve/Krawczyk enclosure or a certified regularity and
  inverse-Jacobian bound converting residual to state distance. The same
  distinction applies to interval energy/momentum evaluations versus bounds on
  the exact trajectory.
- Any position/velocity projection emits
  `ProjectionReceipt { Δq, Δξ, constraint_impulse, ΔE, Δp, ΔL,
  symmetry_defect, adjoint_defect }` and also monitors `J(q)ξ`, acceleration consistency where
  relevant, and `A(q)ξ` in the nonholonomic lane. Only a work/symmetry-compatible
  projection retains the strong invariant claim; otherwise retry/reject or
  debit the defect and propagate it into adjoints.

**Solvers.** Newton–Krylov on the saddle system per implicit step
(fs-solver upgrade); the KKT blocks reuse the block framework. Joint-space
(reduced-coordinate) fast path for open chains [S]; absolute-coordinate +
multipliers for loops (Geneva, four-bar) [S].

**Flexible-machine ladder.** Rigid bodies are only rung 0. Rung 1 adds
Craig–Bampton/component-mode synthesis with certified interface modes and
truncation indicators; rung 2 adds geometrically exact beams/Cosserat rods for
shafts, belts, blades, valve trains, and seals; rung 3 accepts abstract
nonlinear substructure force/compliance operators implemented by `fs-solid`
adapters and composed one-way in `fs-machine` or L6—`fs-mbd` does not depend
back on `fs-solid`. Rotordynamics includes gyroscopic terms, bearing and
seal coefficients, imbalance/misalignment, Campbell diagrams, critical speeds,
Floquet multipliers for periodic operation, and nonlinear whirl/bifurcation
continuation. Substructuring interfaces carry power reciprocity and mode-error
receipts.

**Certificates.** Energy/momentum evaluation/audit per step, with exact-state
enclosure only when the numerical-state error is also bounded;
`ConstraintResidualReceipt` and optional `ConstraintStateEnclosure`; ideal
nonholonomic reaction-work and traction-feasibility receipts; event certificates
from §3.4; mobility/singularity
reports from the tangent constraint complex. Gradient verification gate: adjoint vs.
FD on regular fixed-mode smooth joint-parameter perturbations (link lengths,
gear ratio)—merge-blocking as house rules require. Hybrid parameters require
event-time/saltation or declared generalized sensitivities; otherwise no claim.

**Gauntlet.** G0: momentum maps only for exact discrete symmetries and
conservative/unforced phases; conditional backward-error/modified-energy bounds
only on their admitted smooth fixed-step fixtures, with attributed balance
defects elsewhere; constraint algebra laws. G1: pendulum (elliptic-integral
period), a manufactured-forced multibody trajectory, torque-free symmetric-top
or declared steady-precession closed form, slider-crank kinematics, and four-bar
Grashof classification. G2: a pinned independent high-precision double-pendulum
short-horizon reference and IFToMM
benchmark problems (the rectangular Bricard mechanism — redundantly constrained
throughout motion with a configuration-dependent dependent-constraint subset,
a showpiece for finite/tangent reasoning; Andrews' squeezer — the standard stiff multibody benchmark).
G3: frame invariance (conjugating all
  bodies, twists, inertias/material frames, gravity, loads, contacts, and BCs by
  one fixed motor preserves coordinate-free invariants — a full-scene PGA
  sandwich test), time reversal only on conservative smooth variational
  fixtures, and gear-ratio metamorphic scaling. The nonholonomic battery adds
  the rolling disk and Chaplygin sleigh, plus a constant-width body whose
  contact branch switches under certified no-slip, finite-friction slip, and
  loss-of-contact variants. G4:
cancellation storms mid-Newton; pause/resume bit-equality. G5: bitwise replay
across thread counts.

**No-claims (initial rung).** The first merge is rigid-body only, but flexible
reduction and rotordynamics are funded successor gates in this charter, not a
permanent boundary. No real-time claim until machine- and controller-level
latency receipts exist; contact-rich scenes remain bounded by `fs-contact`.

### 5.3 Contact and Friction — `fs-contact` `[F]`

**Detection.** Broad phase: spacetime AABB/BVH sweep-and-prune over conservative
supports of `SpacetimeChart`s (deterministic ordering). Narrow phase: certified pair distance from
`fs-query` extensions. `min max(φ_A,φ_B)` is retained only as a common-interior
witness; convex configuration-space distance/penetration uses support mappings
and certified GJK/EPA-style bounds, while nonconvex bodies use decomposition or
interval global optimization. **Certified CCD** uses feature-pair spacetime
inclusion/conservative advancement and validated motion tubes. A global
`separation(t)` is nonsmooth and identically zero during persistent contact, so
it is not treated as an ordinary simple-root guard.

**Response (two explicitly different physical/numerical lanes).**
- *Smooth lane — IPC-family barrier [F extension]:* hoist the reusable barrier,
  candidate, and CCD primitives seeded by the feature-gated `fs-solid` path into
  `fs-contact`, then make `fs-solid` consume its adapter—never the reverse.
  Body-body contact uses a log-barrier on certified pair
  distance, rigorous CCD, and filtered line search. Non-intersection is
  conditional on initially admissible geometry, conservative candidate sets,
  accepted CCD-limited steps, and successful optimization—not on the barrier
  name alone. Differentiability is local to a fixed smooth feature/active
  regime; feature changes, medial axes, contact topology, and frictional mode
  changes need generalized sensitivities or an honest no-claim. The current
  live obstacle implementation is only a seed and explicitly defers rigorous
  CCD.
- *Nonsmooth lane — Moreau–Jean [F]:* measure-differential-inclusion
  time-stepping with configuration admissibility `g(q)≥0`, nonnegative normal
  reaction/impulse measure supported on the active set, and complementarity.
  The separate Moreau–Jean discrete velocity/impulse/restitution law includes
  declared gap stabilization; a velocity-only condition may not excuse an
  already penetrating or drifting configuration. It uses Newton
  restitution embedded in a **globally admissible multicontact impact law**, and
  Coulomb cones solved as a cone-complementarity problem
  (deterministic PGS only as a baseline; semismooth Newton, interior-point, and
  conic solvers are required production candidates). Solver receipts report
  primal/dual feasibility, complementarity, determinacy/nonuniqueness, and
  Painlevé/inconsistent cases; impulse, kinetic-energy, angular-momentum and
  contact-work receipts prevent pairwise restitution coefficients from
  silently injecting global energy. A compliant regularized lane is a
  **different physical model**, so its compliance/contact-duration calibration,
  discrepancy, and refinement/limit study are ledgered rather than substituted
  silently.
  This is the lane for hard impacts and long-dwell stick (Geneva locking,
  ratchets, backlash rattle) where barriers stiffen badly. (The contact laws
  and contact-specific cone residuals live here; generic conic/nonlinear math
  stays in `fs-solver`, while the Moreau–Jean *stepper* itself lives in
  `fs-mbd`, which depends on this crate — keeping the graph acyclic.)
- Lane selection is a Rep-Router-style explainable decision recorded in the
  ledger (physical compliance/contact-duration scale, admissible energy/work,
  impact velocity, stiffness/time-step budget, and differentiability demand).

**Friction/material data.** `fs-tribo` owns constitutive friction/film/wear law
evaluation behind the neutral interface-state protocol; `fs-contact` embeds the
returned traction/state/tangent receipt into contact residuals without reaching
into tribology internals. Friction is queried from an `InterfaceSystemCard`
including both surfaces/coatings/roughness/lay, intervening medium or lubricant,
third body, humidity, temperature, pressure, velocity, direction, and history;
it is not an unordered material-pair constant. Stribeck/mixed/boundary laws;
Hertz closed forms in `fs-tribo` as
G1 oracles for contact stress; Archard wear `V = k·F·s/H` accumulating over
the `ModeLedger` (hardness H enters Archard here; ductility enters plastic-
contact, damage, and fracture criteria elsewhere—both become load-bearing
inputs without pretending they enter the same law).
For an assembled contacting pair, action/reaction exchange cancels; total
contact work closes the change in stored contact energy plus nonnegative
dissipation under the declared sign convention. Dissipated work is partitioned
among declared heat, plastic/damage/defect/surface-energy changes, lubricant
chemistry and exported debris/species; only the heat share feeds thermal, and
the total closes. A passive impact may not increase
assembled kinetic-plus-stored energy absent a source; an active interface must
declare that source.

**Gauntlet.** G0: complementarity KKT residuals; under the smooth lane's stated
admissibility, candidate-set, CCD, line-search, and convergence preconditions,
an interval-checked accepted step remains nonintersecting. G1: Hertz sphere/plane
and cylinder/plane closed forms; bouncing ball with restitution (analytic
impact-time reference); block-on-incline
stick/slip threshold `tan θ = μ_s` for an ideal rigid block with no other
loads. G2: an exactly pinned Painlevé rod deck
whose parameter regime has an independently derived existence/nonexistence/
nonuniqueness classification—the nonsmooth lane must report that classification
rather than silently choose an impulse; feeder/ratchet benchmarks; disk stacking (Chrono/Siconos cross-checks as
dev-only oracles). G3: a mirrored scene transforms every oriented geometry,
frame, normal, friction direction, spin, load, and BC, then agrees within the
declared interval/tolerance (bitwise only for exact signed-permutation paths);
the restitution-to-persistent-contact limit is tested only on named well-posed
fixtures under a specified time-step/refinement path. G4/G5 as usual.

**Initial-rung boundaries and funded successors.** Rigid/deformable and
rigid/rigid land first. Deformable/deformable solids, shells, rods, seals, and
adhesive/cohesive contact are explicit later gates, not excluded from the
charter. Thermal/electrical contact, frictional heat partition, wear-evolving
geometry, and fluid-film traction follow typed interface receipts. The
nonsmooth lane begins without a classical adjoint; smooth fixed-feature and
generalized-sensitivity lanes must state their actual regime.

### 5.4 Kinematic and Machine-Element Libraries — `fs-kinematics` + `fs-machine` `[S→F]`

The user-facing vocabulary surface: this is where "simulate a Geneva drive"
becomes a three-line scenario.

- **Joint catalog (`fs-kinematics`):** revolute, prismatic, cylindrical, universal, spherical,
  planar, screw (pitch-coupled — natural in PGA), rack-and-pinion, gear pair
  (ratio + backlash band + mesh-stiffness hook), cam-follower (profile-driven
  constraint with lift/velocity/acceleration/jerk export), pin-in-slot
  (Geneva). Each joint supplies a finite constraint, verified Jacobian and
  dual reaction map; a sheaf restriction map is added only where the declared
  local observation spaces satisfy the sheaf axioms.
- **Nonlinear/tangent kinematics** (§3.2): holonomy, infinitesimal mobility,
  dual self-stress, second-order rigidity, continuation, singularity margin.
- **Gears (`fs-machine`):** involute generator (module, pressure angle, profile shift,
  addendum modification), cycloidal and circular-arc profiles, plus internal,
  helical/double-helical, bevel/spiral-bevel, hypoid, worm/crossed-axis, face,
  noncircular, and flexible-ring/strain-wave families. Spatial families use
  cutter/tool generation, SE(3) meshing equations/screw theory, contact-line or
  point classification, loaded tooth-contact analysis, and manufacturing flank
  deviations; no cylindrical formula or standard scope is silently reused.
  Conjugate certificates (§3.3) cover the appropriate surface-contact family;
  generator singularities, trimming, branch visibility, addendum/root limits,
  and tool interference complement planar Euler–Savary analysis for planar
  gears and family-specific surface/contact-curvature analysis for spatial
  gears. Report
  rolling/sliding/spin creepage and separated unloaded/static/dynamic TE. Mesh stiffness (Weber–
  Banaschek-style analytic first [S], 2-D FE per-tooth via `fs-solid` +
  `fs-symmetry` sector solve later [F]), AGMA/ISO bending/pitting calculations
  with exact edition/scope and worked-example conformance. Formula
  implementation is numerically `Verified`; the empirical rating model is not
  physically `Validated` by the standard's name, and machine-life prediction
  becomes `Validated` only against independent component/service evidence for
  its named QoI. A separate `StandardConformance` receipt records this axis. Planetary
  train composition (constraint quotient checked against Willis), backlash
  from tolerance stacks, microgeometry and lead/profile modifications,
  misalignment, tooth-root/contact fatigue, scuffing/micropitting, thermal
  distortion, EHL/mixed lubrication, bearing/shaft/housing flexibility, and
  order-domain vibration/noise excitation.
- **Machine-element continuum (`fs-machine`):** rolling-element, journal,
  magnetic, and compliant bearings; bushings; springs/dampers; torsional and
  misalignment couplings; dry/wet clutches and brakes; belts, chains, timing
  drives, capstans, lead/ball screws, ratchets/freewheels, seals, and compliant
  mechanisms. Pure geometry/constraint pieces reuse `fs-kinematics`; loaded,
  flexible, thermal, lubricated, life, and NVH behavior is composed by L6 or
  one-way domain adapters from `fs-mbd`, `fs-solid`, `fs-tribo`, `fs-material`,
  and `fs-acoustics`; `fs-machine` itself owns schemas/residuals and does not
  import the omitted peer solvers. Each element has a
  rigid/reduced/high-fidelity ladder and a failure-mode model card rather than
  one universal empirical coefficient.
- **Geneva mechanism:** parametric generator (n stations, pin radius, crescent
  clearance), engagement/disengagement events, and two honest retention cases.
  An ideal zero-clearance declared active boundary mode may earn a **first-order
  form-closure** cone/duality certificate plus
  `λ_{m+1}(N_V)` or the equivalent squared
  `σ_min(J̃|_{K_y^⊥})` margin in the whitened coordinates of §3.2; degenerate
  modes require second-order or
  global configuration-space proof. A manufactured positive-clearance device
  instead earns `DwellRetentionCertificate { play_set, nonescape, boundary_mode_margins }`
  over declared tolerance, load, impulse, elastic, friction, wear, and mode
  uncertainty. Bounded play/creep is measured rather than denied. Index accuracy,
  peak contact stress and dry-wear baseline use E2 `fs-tribo` receipts.
- **Constant-width/trochoid module:** constant-width curve generators, epitrochoid/
  hypotrochoid paths, the Wankel geometry kit (rotor flank derived and checked
  as a conjugate/envelope design under relative motion; ideal apex-point path,
  finite seal-center/contact loci, and actual bore constructed and certified as
  distinct objects under the declared seal-tip and clearance model — §3.3;
  apex/side/corner-seal sliding kinematics and spring/centrifugal/gas-pressure
  seal dynamics; ideal nominal chamber-volume `V(θ)` closed-form oracle, with
  finite seals/clearance/ports/deformation/as-built geometry handled by certified
  quadrature), plus genuinely rolling constant-width mechanisms as a separate
  nonholonomic family with no-slip feasibility and slip/loss-of-contact
  transitions. The plan never equates those motions.
- **Linkage analysis & synthesis:** four-bar/six-bar position/velocity/
  acceleration analysis with certified loop closure; Grashof classification as
  an interval certificate; dimensional synthesis (3–5 precision points via
  Burmester theory) posed as an `fs-opt` problem with `fs-constraint` typed
  constraints — linkage synthesis was already name-checked in the master plan
  as an SOS-scale problem (plan L410); this is its home.

**Gauntlet.** G0: frame equivariance; analytic/AD/FD agreement for `J=DC`;
primal/dual work identity; Maxwell–Calladine index; exact catalog nullities;
finite-continuation and second-order checks at singular fixtures. `δδ=0` on a
bare joint graph is vacuous and is not an acceptance test. G1:
slider-crank/four-bar closed-form kinematics; a properly generated, compatible-
base-pitch, noninterfering ideal involute pair has constant ratio under declared
rigid center-distance/contact assumptions; Geneva ideal geometry and
positive-clearance dwell-play bounds; epitrochoid
closed form; the exactly parameterized Bennett family supplies a G1 symbolic
closure/nullity oracle (mobility 1 where Grübler says −2). G2 executes a pinned
Bennett/IFToMM-style deck as a nonlinear/tangent acceptance test. A gear-TE G2 case must provide a
reproducible geometry/error/load/support input deck and trace; a definition-only
paper is not promoted into an executable benchmark.
Rolling-disk/Chaplygin-sleigh and constant-width branch-switch decks exercise
Pfaffian constraints, ideal reaction zero work, friction feasibility and slip.
G3: scale/frame invariance of all certificates. G5: mechanism trajectories
bit-replay.

### 5.5 Electromagnetics — `fs-em` `[F]`

**Formulations, in delivery order.**
1. **2-D planar nonlinear magnetostatics [first]:** `A_z` scalar formulation,
   with a weak residual that distinguishes impressed winding current from the
   remanence/magnetization term and records its sign/orientation convention, on
   triangle complexes (new fs-feec 2-D support) or the existing quadtree CutFEM
   machinery. This is the high-value first rung for the PMSM/BLDC flagship.
   Fit/enclose a convex magnetic energy `w(B)` with `H=∂w/∂B`, or its admitted
   Legendre-dual coenergy `w*(H)` with `B=∂w*/∂H`. That potential-based lane
   verifies a symmetric Hessian, convexity/eigenvalue bounds and Legendre
   duality, and may support the declared virtual-work functional. A separate
   strongly monotone, nonintegrable `H(B)` lane verifies positivity of the
   symmetric tangent part/coercivity but earns no stored-energy/coenergy force
   claim unless a potential is separately proved. A monotone scalar spline of
   `ν(B)` alone proves neither an SPD tangent, stored-energy integrability, nor
   Newton convergence. Newton uses globalization; a solution/convergence
   enclosure requires Krawczyk/radii-polynomial or validated inverse-Jacobian
   machinery for the assembled residual.
2. **3-D magnetostatics:** `curl(ν curl A) = J`, Nédélec elements
   (`vecfam.rs`), gauge by tree-cotree (deterministic spanning tree keyed by
   cell IDs) or by Coulomb-gauge Lagrange multiplier; multiply-connected
   domains via `harmonic_basis`/`deflate_harmonics` — the cohomology tools
   exist (`fs-feec/src/cohomology.rs`).
3. **Moving-conductor MQS (A–φ/A–V):** solve the coupled weak field and charge-
   continuity system on explicitly partitioned conducting/insulating domains,
   with `J = σ(-∂A/∂t - ∇φ + v×B)` in the chosen frame, source compatibility,
   interface/boundary conditions, gauge/cohomology constraints, and energy
   reciprocity. ALE/rotating-frame variants must agree on overlap fixtures.
   This rung is mandatory for generators, induction rotors, solid-conductor
   eddy currents, and motional EMF; a sliding interpolation alone cannot
   replace it.
4. **Electrostatics/current flow (0-form/scalar) [cheap]:** separate
   permittivity-weighted electrostatic and conductivity-weighted steady-current
   elliptic solves — compact first consumers of weighted operators,
   with boundary/source/coercivity obligations still explicit—and valuable
   immediately (resistance/capacitance extraction, Joule source for
   `fs-thermal`).

**Fidelity ladder.** 2-D static/rotor-angle sweeps land first; 2.5-D multislice
adds skew and calibrated end corrections; full 3-D moving-conductor MQS resolves
end winding, leakage, axial flux, and broken-symmetry effects within its declared
MQS/material model. The reduced lanes
are labeled and carry an escalation trigger from measured end-effect error.
PMSM/BLDC is the first forcing function, not the taxonomy boundary: induction,
wound-field, synchronous- and switched-reluctance, axial/radial/transverse-flux,
linear machines, transformers, solenoids, magnetic bearings, eddy-current
brakes, and electrostatic actuators each select the formulation/material/loss
ladder their physics requires rather than inheriting a motor template by name.
Lamination homogenization, skin/proximity effects, rotor-bar end rings, common-
mode/cable parasitics, bearing currents, insulation stress/partial discharge,
and conducted/radiated EMC are explicit machine successors. Electroquasistatic
and full-wave rungs are distinct validity regimes with wavelength/displacement-
current escalation tests; MQS is never extended by name beyond its asymptotics.

**Solvers.** HX/AMS-style auxiliary-space preconditioning and harmonic
deflation serve properly gauged SPD magnetostatics, where CG is admissible.
Verified symmetric-indefinite saddle blocks may use MINRES; nonsymmetric,
non-Hermitian real-equivalent, or variable-preconditioner systems use FGMRES,
with block Schur and coefficient/topology-robust preconditioners. Solver admission checks symmetry,
definiteness, nullspaces, and source compatibility rather than selecting by
physics name.

**Sources & materials.** Stranded windings via the neutral winding realizations
(§3.6); solid conductors with topology- and formulation-specific imposed-current
constraints; permanent
magnets as precisely oriented remanence/recoil constitutive offsets and RHS
terms; cards include temperature-dependent demagnetization/recoil limits.
Single-valued nonlinear BH uses the convex-energy/strong-monotonicity cards
above. Hysteresis
(Jiles–Atherton/Preisach) `[M]`, core loss initially via Steinmetz/Bertotti maps
as `Validated` post-processing only against named material/frequency/temperature/
waveform data (D6). Vector/hysteretic laws carry internal state and a discrete
free-energy/dissipation receipt; they are not forced into a path-independent
stored-energy fiction.

**Forces/torques.** Maxwell stress tensor on air/vacuum contours (or an
explicitly declared material stress tensor) with the sliding
interface's separate compatibility and conservation audits, plus virtual-work/
eggshell method cross-check. Method
agreement is an implementation-consistency gate, not physical validation and
not a certificate by itself. Ordinary DWR supplies an error *estimate*;
`Verified` force bounds require a guaranteed functional majorant/equilibrated
or interval-residual bound with stated reliability constants. Disagreement
demotes and localizes. Adjoint consistency requires a domain-specific shape
derivative of the actual discrete EM residual and force functional, verified
against independent perturbations; existing Hadamard infrastructure is reusable
machinery, not a proof that this derivative is already correct.
Virtual work states the held electrical variables and functional—e.g.
`+∂W'/∂x|I`, `-∂W/∂x|λ`, or the coupled-circuit functional—plus deformation/
magnetostriction convention, torque origin and frame. A force theorem is scoped
to that material, boundary, gauge and electrical-control choice.

**Rotation.** Equivariant sliding interface (§3.5), with a sheaf compatibility
layer and separately conservative Fourier/harmonic coupling on
the air-gap circle; sector symmetry with Bloch phases for pole-pair reduction
(fs-symmetry upgrade); rotor angle driven by `fs-mbd` (co-simulated through
the port bus: electromagnetic torque out, angle in).

**Gauntlet.** G0: incidence exactness separately from weighted-map symmetry,
coercivity/monotonicity, units, quadrature and Jacobian consistency; gauge
invariance of B and gauge/source compatibility; reciprocity/adjoint and
field-circuit energy identities. Moving-conductor fixtures close the Poynting /
Lorentz-mechanical / conductor-frame-Ohmic three-way receipt. G1:
MMS with smooth manufactured A (convergence order per element degree, slope
within 0.2 per house rule); coaxial cable, sphere in uniform field, Helmholtz
coil closed forms. G2: **TEAM workshop problems** — the established international
TEAM benchmark family: selected, exactly identified TEAM 10/13/20/24/30a/30b
variants covering nonlinear magnetics, force/torque, eddy-current, and rotating-
machine behavior — their geometry revision, excitation, QoI definition, source,
and acceptance envelopes are stored in-repo rather than inferred from the
problem number. G3: frame rotation equivariance and unit-rescaling invariance
(volts vs millivolts); variational energy monotonicity only on fixtures whose
space nesting, linearity, quadrature, and minimization assumptions imply it,
with ordinary convergence studies elsewhere. G4/G5 house-standard.

**Initial-rung boundaries and successors.** The first production scope is
magnetoquasistatic; full-wave/antenna and superconducting constitutive lanes
remain `[M]`. Moving-conductor EMF is **inside** the funded machine path.
Hysteresis starts off-default but advances through loss-map, dynamic hysteresis,
and vector-hysteresis rungs with calorimetric/loop validation.

### 5.6 Circuits, Power Electronics, and Control — `fs-circuit` + `fs-power` + `fs-control` `[S→F]`

**Circuit kernel.** Descriptor MNA uses exact graph incidence plus structural
rank/index analysis, consistent initialization, topology-aware state mapping,
and, for an admitted **impulse-free** topology change, explicit capacitor-
charge/inductor-flux state continuity. When a
topology change imposes inconsistent post-switch constraints, a distributional
impulse solve may change individual capacitor voltages/charges or inductor
currents/fluxes while enforcing KCL/KVL and whatever charge, energy and source
balances the declared closed boundary actually supports, with an explicit
energy-defect receipt; otherwise the solver refuses or introduces declared
regularizing parasitics. Ideal switches and
diodes form a complementarity/hybrid lane that detects simultaneous
commutation, inconsistent resets, impulses, chatter, and Zeno behavior; a
physically regularized device lane is the production fallback. A certified
event time alone is not a certificate that the post-switch DAE is well posed.

**Power electronics.** `fs-power` adds temperature-dependent diode, MOSFET,
IGBT, and SiC/GaN device cards; reverse recovery, output capacitance, gate
dynamics, dead time, bus/wiring parasitics, saturation, and switching/conduction
loss. Averaged/dq, cycle-averaged loss, and fully switched PWM/SVPWM lanes share
one terminal contract and quantify model discrepancy. Converter heat returns
to the thermal graph; electrical and thermal states iterate or solve
monolithically when feedback is stiff.

**Closed-loop control.** `fs-control` supplies sample/hold, multirate clocks,
Clarke/Park transforms, current/speed/torque loops, anti-windup, MTPA, field
weakening, observers/PLL/Kalman-family filters, sensor noise/bias/faults,
quantization, delay and saturation. L3 executes solver-neutral policy artifacts
and runtime safety monitors; existing L4 `fs-opt`, `fs-sos`, and `fs-robust`
own MPC and H-infinity synthesis, SOS/reachability, robust tuning,
control-by-interconnection design, and hybrid-shield construction, then lower
versioned policies into `fs-control` without a reverse dependency.
Every Clarke/Park artifact declares power- versus amplitude-invariant
normalization, phase order, mechanical/electrical angle, pole-pair count, and
sign convention. Local LTI margins do not certify a saturated, delayed,
quantized, periodic or hybrid loop: those lanes require lifted/Floquet,
IQC/Lyapunov/reachability/region-of-attraction evidence plus observability,
detectability and fault-isolation receipts.
Periodic machines report Poincaré maps and Floquet multipliers; continuation
tracks subharmonics, chatter, pull-in/loss-of-synchronism, and controller-induced
bifurcations. Generated controller code and fixed-point arithmetic replay the
same scenario before any HIL claim.

Reduced dq PMSM/induction models remain valuable `Estimated` or
numerically-verified reduction lanes when compared only with field solves;
they become `Validated` only for named dynamometer/electrical datasets and QoIs.
Field discrepancy still drives calibration/escalation. Gauntlet: KCL/KVL and energy/storage
identities; impulse-free state-continuity maps and distributional impulse/
energy-defect/refusal fixtures; RLC and converter closed forms;
versioned canonical buck/boost, rectifier, two-level-inverter and drive-cycle fixtures;
closed-loop stability margins, anti-windup/saturation metamorphics, fault
injection, and bit-replay.

### 5.7 Thermal and Energy Transport — `fs-thermal` `[S→F]`

**Field lane:** semidiscrete storage and transport are separate. A regular
single-phase/apparent-heat-capacity lane uses
`M_{H_T}(T,z) Tdot + M_{H_z}(T,z) zdot + D₀ᵀ M_k(T,z) D₀T = f + b(T)` with
surface mass,
Neumann/Robin, contact, radiation, and moving-boundary terms oriented and
unit-checked. Implicit first-order steppers, nonlinear enthalpy derivatives,
anisotropic composites, temperature-dependent properties, and direct
energy-balance receipts—with enclosure only where the numerical and data paths
support one—are required. The existing timedependent adjoint is a
pattern, not proof for nonlinear phase/radiation/contact paths.
If internal state `z` is frozen or energetically neutral, the card may omit the
`H_z zdot` term explicitly; otherwise a total-enthalpy residual must include it.
**Network lane:** lumped RC thermal networks (nodes = machine components,
resistances from conduction/convection cards) — storage nodes, lossless
junctions, dissipative conductances, and reservoirs are distinct; engine/motor thermal management is
network-first in practice (D7).
**Sources:** Domain solvers export typed source/feedback operators through
`fs-couple`; `fs-thermal` never reaches into peers. Sources include Joule `J·E`
only for a stationary Ohmic convention. For moving conductors with
`J_c=σ(E+v×B)`, conductor-frame heat is
`q_Ω=J_c·(E+v×B)=J_c·sym(σ⁻¹)J_c≥0` for passive resistivity, while Lorentz mechanical power is
`(J_c×B)·v`; under the declared sign split,
`J_c·E=q_Ω+(J_c×B)·v`. `fs-em` exports distinct field/source, Ohmic,
mechanical and time-harmonic-average receipts, and thermal consumes only the
loss contribution. Hall/thermoelectric cross-effects use the coupled
Onsager–Casimir entropy-production receipt rather than this scalar shortcut.
Eddy/core losses likewise retain their convention and
validation regime. Friction heating from tangential traction
contracted with slip velocity, with sign and heat partition explicit, at
`fs-contact` interfaces (the tribo-thermal loop); combustion
heat release from `fs-gas`. Heat partition and boundary temperatures are
explicit so energy and entropy audits close over accounting windows.
`LossOwnershipReceipt { loss_id, resolved_or_modelled, upstream_debit,
thermal_credit, exclusions }` prevents double counting: Steinmetz/Bertotti
maps are excluded when resolved hysteresis/eddy physics already contains that
loss; winding resistance is reconciled with field-resolved copper loss; and
converter loss maps are excluded from resolved switching events. Thermal credit
must equal the upstream electrical/mechanical debit within the accounting bound.
**Phase change:** a primary-enthalpy/inclusion formulation handles sharp latent-
heat plateaus or jumps, where `dH/dT` may be singular or set-valued; apparent
heat capacity is an explicitly regularized lane. Interface tracking uses the
level-set kernels hoisted from `fs-topols` into `fs-xform`
(§4.2) where sharp fronts matter (Stefan problem as G1 oracle). This is where "latent heat density" from the user's requirements
lives. Level-set transport alone does not conserve latent heat: the phase-front
flux jump and remap receive a direct balance certificate.
**Radiation/cooling/coupled response:** nonlinear `T⁴` surface exchange lands
as a funded `[F]` lane with view-factor/radiosity reciprocity and enclosure;
participating-media radiation follows. Conjugate heat transfer, coolant
networks, boiling/evaporation, thermal contact/gap transfer, thermoelastic
expansion/stress, and bidirectional electrothermal/magnetothermal feedback are
part of the staged charter, assembled by `fs-couple`/L6 adapters rather than
peer-to-peer crate dependencies.
**Gap/contact conductance:** `InterfaceSystemCard`s depend on pressure,
roughness, coatings, gas/lubricant, temperature and history—not merely a pair
key—and couple to contact normal traction.
**Gauntlet.** G1: 1-D/2-D conduction, radiative exchange, laser-flash and
Stefan balances; G2: exact-version NAFEMS cases with licensed/provenanced QoIs,
conjugate and thermoelastic cases, motor/engine calorimetry; G3:
superposition in the linear regime, scaling invariance; G5 replay.

### 5.8 Gas, Thermochemistry, and Combustion — `fs-thermochem` + `fs-gas` `[S→M]`

0-D/1-D system models are the first useful engineering rungs; multi-D chamber
and spray models are indispensable for eventual world-class combustion design.
The roadmap stages them—it does not declare the most ambitious rung unfunded.

- **Tier 0 — typed species and reference thermodynamics [S]:** immutable
  `SpeciesId`, elemental composition and charge; NASA-9 standard-state molar
  `cp°,h°,s°` with reference pressure and temperature regions, while `u°` and
  `g°` are derived under explicit phase, EOS, reference-pressure, and elemental-
  reference conventions; typed molar masses;
  mixture entropy and transport mixing. Mass fractions remain convenient state
  variables, but conversions use the amount-of-substance dimension. Reaction
  matrices satisfy elemental and charge conservation exactly. Ideal,
  calorically/thermally perfect, cubic real-gas, and tabular multiphase EOS form
  an escalation ladder with phase-stability checks.
- **Tier 1 — conservative 0-D control volumes [S→F]:** state includes total
  mass and species. A no-momentum lumped volume may evolve total internal
  energy; if momentum is evolved, the conservative state uses total energy
  including kinetic energy and a declared potential convention. Both use an
  explicit reference/formation convention. The energy derivation chooses
  either total chemical energy or sensible energy plus an explicit reaction
  source—never both, preventing double-counted heat release. Intake/exhaust,
  residual gas/EGR, crevice volumes, blow-by, valve lift/discharge, wall heat,
  fuel evaporation, combustion phasing, friction and leakage close the balance.
  Numerical enclosures, parametric uncertainty, and model-form discrepancy are
  reported as separate typed components of the p–V/indicated-work budget; they
  are combined only under an explicitly declared probability or worst-case
  semantics, rather than implying a correlation has become verified.
- **Tier 2 — 1-D reactive gas networks [F]:** conservative state
  `U=(ρ,ρu,ρE,ρY₁…ρY_{Ns-1})`, area/junction/valve sources, heat, friction and
  chemistry. MUSCL–HLLC is an `Estimated` baseline. A separately proved
  entropy-stable and positivity/invariant-domain-preserving lane includes
  boundary and source entropy flux, well-balanced area terms, and explicit
  mathematical-entropy sign convention. Acoustic waves feed `fs-acoustics`.
- **Tier 3 — chemistry and phase ladder [F→M]:** a calibrated Wiebe lane remains
  `Estimated` until an independent named experimental dataset and QoI validate
  it; calibration data are not reused as validation evidence. Global then
  skeletal/detailed kinetics with stiff integration;
  autoignition, knock/detonation, multicomponent diffusion, liquid films,
  injection/spray breakup and evaporation, real-fuel surrogates, NOx/CO/UHC/
  soot and aftertreatment. Each rung has mechanism/version provenance and
  uncertainty, never a bare “latest” database. Exact element/charge
  conservation is necessary but insufficient: reverse rates/equilibria must be
  thermodynamically consistent with the chosen standard states, and the
  reaction update must satisfy its declared nonnegative entropy-production or
  free-energy-dissipation inequality.
- **Tier 4 — multi-D chamber/port flow [F→M]:** existing `fs-flux` owns the
  compressible FVM/DG spatial/ALE/CutFEM implementation and consumes one-way
  `fs-gas` thermochemical, Riemann, reaction, and boundary closures. Moving
  domains use discrete GCL, conservative remap, valve/piston/
  rotor motion, conjugate heat transfer, wall films, DNS resolution/convergence
  lanes, and LES/RANS closure model cards. The certifiable reduced solver selects escalation regions;
  multi-D results calibrate closures back down.

**Gauntlet.** G0: positivity, `ΣY=1`, element/charge conservation, reference-
state consistency; for frozen-composition ideal gas,
`c̄p-c̄v=R_u` on a molar basis or `cp-cv=R_mix=R_u/M_mix` on a mass basis
(reacting/equilibrium derivative paths do not inherit that shortcut); EOS
Maxwell identities, source and entropy balance. G1: exact Riemann/nozzle and
air-standard cycles with their assumptions; canonical constant-volume and
flame/ignition reactors. G2: versioned/licensed CFR and motored pressure traces,
manifold-wave, spray, combustion, emissions and acoustic datasets with named
QoIs. G3 uses reversible-limit, Galilean/unit, equilibrium and refinement
metamorphics; it does **not** demand monotonic entropy production under mesh
refinement, which is not a general invariant.

### 5.9 Moving-Boundary Fluids and Lubricated Machines — `fs-flux` + `fs-lbm` `[F→M]`

The live `fs-flux` path is already a 2-D H(div)-conforming, exactly
divergence-free incompressible solver and must not be omitted. Its machine
ladder is 3-D/high-order curved geometry → CutFEM/ALE with discrete geometric
conservation and conservative remap → monolithic/partitioned FSI → cavitation,
free-surface/multiphase, DNS resolution/convergence, and LES/RANS closure-card
lanes. This is the
high-accuracy continuum route for pumps, turbines, cooling jackets, bearings,
valves and FSI.

The complementary LBM lane adds interpolated curved moving boundaries,
Galilean-invariance-corrected momentum exchange with a measured residual,
deterministic fresh-cell initialization
with mass/momentum receipts, partial saturation, under-resolved lubrication
correction, and MRT/cumulant collision. Free-surface wetting queries an
`InterfaceSystemCard`: solid, liquid, gas, roughness, contamination and
advancing/receding hysteresis—not a material-pair angle.

`fs-tribo` couples Reynolds/Elrod–Adams cavitating films, thermal viscosity,
elastic deformation, mixed asperity contact and EHL to bearings, gears, piston
rings and apex/side seals. Benchmarks include moving Couette/Taylor–Couette,
oscillating and sedimenting bodies, impellers/gerotors, journal/thrust bearings,
gear-mesh films, cavitating pumps and seal leakage, each with force/flow/mass/
power QoIs and a direct GCL/balance receipt.

Liquid-hydraulic and pneumatic networks are a funded companion lane: compliant
lines, accumulators, pumps/compressors, cylinders, spool/check/relief valves,
orifices, leakage, water hammer, dissolved/aerated gas, cavitation and thermal
state. The 0-D/1-D network uses the same stream/storage/dissipation contracts;
multi-D valve/pump regions escalate into `fs-flux`. Fluid bulk modulus, vapor
pressure, gas dissolution, viscosity and speed of sound are typed material
inputs, never generic constants.

### 5.10 Materials Data and Tribology — `fs-matdb` + `fs-tribo` `[S→M]`

**Layer/ownership boundary.** L1 `fs-matdb` owns immutable observation datasets,
property claims, constitutive-model cards, initial/reference-state
distributions, typed validity predicates, opaque basis/frame identifiers,
evidence links, and query/fusion policy. It owns neither executable closures nor
per-run history. Its neutral law identity/version/schema metadata is consumed
by L3 `fs-material`, `fs-tribo`, `fs-em`, `fs-gas`, and other adapters that own
executable implementations and mutable typed runtime states. L1
`fs-thermochem` owns lower-layer thermodynamic/kinetic closure math and types,
not session state; its L3 consumer owns evolution. L6 ledger/package owns
content-addressed checkpoints. `fs-matdb` never imports L2 `fs-ga` transforms, L3
domain/runtime state types, or L6 persistence types. `fs-thermochem`
and other domains consume cards one-way; basis transforms occur in L2/L3;
`fs-ledger`/`fs-package` persist versioned cards and receipts in L6.

**Schema (see Appendix B).** A `MaterialCard` identifies chemistry, phase,
temper/heat treatment, processing route, lot/spec revision, microstructure,
texture/orientation and aging state. Each property carries a unit- and
frame-aware value/curve/tensor, the single upgraded
`fs-evidence::ValidityDomain` predicate algebra over typed axes (T, p, frequency,
field, strain/rate, composition, history or named dimensionless groups),
interpolation/extrapolation policy, distribution/covariance and correlations,
calibration lineage, source/license/`fs_blake3::ContentHash`, schema revision,
supersedes links, and evidence class. The current box validity type remains a
compatibility projection; a competing second validity-domain type is forbidden.
Conflicting observations remain separate `PropertyClaim`s indexed under one
property key. Fusion/selection is an explicit policy, and joint covariance or
correlation datasets may span properties and claims; no map overwrite invents
a canonical value.

Every query returns `Evidence<PropertySample>` plus a `PropertyUsageReceipt`.
Only load-bearing receipts enter a claim's dependency closure: irrelevant weak
data do not demote an unrelated output, while relevant weak data cannot be
laundered. A source citation alone does not automatically make a handbook
value `Validated`; specimen/process match, uncertainty and validation context
are required.

**Property vocabulary (initial):** density; elastic moduli; yield/ultimate/
elongation (ductility); hardness (HV/HB/HRC with named conversion
correlations); fracture toughness, S–N/ε–N and crack-growth data; creep,
viscoelastic, damage/corrosion and aging laws; thermal conductivity, specific heat,
expansion, latent heats + transition temperatures; electrical conductivity/
resistivity (with α_T); relative permeability + BH curves + remanence/
coercivity/energy product, susceptibility and magnetization/specific-moment
curves for magnets and magnetic media (total magnetic moment is then integrated
over the actual body, not stored as a geometry-free material scalar);
permittivity + loss tangent; viscosity
(Newtonian, Sutherland, Carreau — the last already implemented in scenarios);
dielectric strength/breakdown and partial-discharge inception/insulation aging;
damping/internal friction and strain-rate plasticity; moisture sorption,
diffusion and swelling; thermal decomposition, flammability and outgassing;
fluid bulk modulus, vapor pressure, surface tension, speed of sound, flash/fire
point, chemical compatibility and solubility; roughness spectra, lay and coating
thickness; diffusion/permeability tensors; lamination geometry and vector-loss
data; surface energy; gas permeability/diffusivity; emissivity; phase-diagram/
coexistence data plus activity, fugacity, and EOS parameters from which
`fs-thermochem` computes phase equilibrium;
magnetostriction, piezoelectricity and thermoelectric coefficients. Interface
properties live in `InterfaceSystemCard { surface_a, surface_b,
intervening_medium, lubricant_or_third_body, environment, texture_frame,
constitutive_models }`: directional friction, restitution, thermal/electrical
conductance, wear/scuffing, permeability, adhesion, and advancing/receding
wetting hysteresis. Each constituent model owns its own state schema and initial
policy. When several laws coexist or couple, L3 `ConstitutiveGraph` admission
generates an `AggregateStateSchema` and initialization receipt; no ambiguous
singular interface-wide schema is stored in L1. Current wear, lubricant,
contamination, temperature and memory state are passed separately to each law
query and checkpointed by the run. Wetting is a solid–liquid–gas system, not an
unordered pair.

**Seed dataset (curated, committed, cited):** structural steel 4140 and
AISI 1045; aluminum 6061-T6; gray cast iron and representative aluminum/
wear-resistant-coating systems for explicitly named engine and Wankel
configurations (never treated as universal construction);
copper (OFHC); a pinned bearing steel; a named carburized gear-steel/process;
winding enamel, slot liner and impregnation epoxy; a named bearing grease;
water-glycol coolant and application-specific gear/bearing oils;
electrical steel M-19 (BH curve + Steinmetz coefficients); NdFeB N42 and
ferrite Y30 magnets; PTFE, PEEK, nitrile (seals); air; iso-octane as a reference
fuel, followed by a versioned multicomponent gasoline surrogate; methane;
exhaust-gas mixture; and a fully identified
redistributable 5W-30 reference formulation—viscosity grade alone is not a
material identity. Electrical-steel thickness/anneal/test method and magnet
supplier/process/temperature state are pinned rather than inferred from grade.
Seed data use
redistributable public/licensed sources and preserve conflicts rather than
inventing a canonical number. “Real material” means a named condition and
uncertainty, not a marketing-grade label; flagships may refuse when the needed
process state or interface system is missing.

**`fs-tribo`:** Hertz closed forms (G1 oracles), Coulomb→rate/state/Stribeck
friction, thermal partition and flash temperature, Reynolds/Elrod–Adams films,
mixed asperity contact, EHL, cavitation, scuffing/micropitting/pitting, and
wear-evolving geometry. Archard and film-thickness correlations remain
`Validated` or `Estimated` model cards; a coupled resolved lane is the
certify-or-escalate successor.

An offline L6/`xtask` pack compiler—not L1—performs source parsing, license
admission, unit/basis normalization, covariance preservation and canonical
hashing for material/species/kinetics packs. Runtime L1 consumes only pinned
normalized artifacts. Every usage receipt binds the query point, basis/frame
transform, interpolation/extrapolation decision, selected covariance slice,
law/evaluator version and all source hashes.

### 5.11 Acoustics and NVH — `fs-acoustics` `[F→M]`

**Formulation ladder.** Lumped acoustic impedances and 1-D characteristic duct
networks land first. Interior acoustics then uses first-order pressure/velocity
time-domain systems and frequency-domain Helmholtz formulations with complex
impedance boundaries; structural modal/harmonic solvers export normal velocity
and consume acoustic traction through a power-reciprocal interface. Exterior
radiation uses Helmholtz BEM/FMM with Sommerfeld radiation, singular/
hypersingular quadrature, and Burton–Miller or another proved combined-field
resonance treatment. Convected Helmholtz is admitted only for its mean-flow
class; nonuniform vortical flow escalates to linearized Euler/acoustic-
perturbation equations or direct CAA. Thermoviscous narrow-gap/seal lanes and
thermoacoustic combustion-instability/Rayleigh-index analysis are explicit
successors.

**Source ladder and coupling.** Rotating electromagnetic force harmonics,
gear/bearing orders, imbalance, combustion pressure, valve/flow sources and
structure-borne paths retain source identity and clock/phase convention.
Lighthill and FW–H hybrid lanes record source surface, quadrature, propagation
model and turbulence/source-model discrepancy; agreement between two acoustic
propagators does not validate the source. The discrete structural/acoustic
normal-velocity–traction pairing must close interface power in both time and
frequency domains. Aeroacoustic, vibroacoustic and combustion-acoustic
couplings may escalate monolithically when partitioned stability fails.

**QoIs and evidence.** Report referenced physical pressure, particle velocity,
intensity, sound power, modal participation, order spectra, directivity and
transfer paths. SPL/SWL and psychoacoustic metrics carry reference value,
weighting, window, bandwidth, calibration and observation geometry; logarithmic
units are semantic types, not dimensionless floats. Microphone, anechoic/
reverberant-room, modal-hammer and dynamometer artifacts include clock
synchronization, instrument calibration, repeatability and covariance.
Frequency-domain power uses a declared sesquilinear phasor convention: for peak
phasors, for example, `P̄=(1/2)Re∫Γ(-p n)·v_s* dΓ` (no `1/2` for RMS), with
conjugation, orientation and peak/RMS recorded. Transpose reciprocity is not
confused with the Hermitian energy adjoint in lossy or convected media.

**Gauntlet.** G1 includes duct/cavity modes, pulsating sphere, baffled piston,
plane-wave impedance, scattering and structural–acoustic manufactured
solutions. G2 pins NAFEMS R0083, NASA CAA benchmark decks, measured transfer
functions and combustion-acoustic cases with exact geometry, medium, source,
QoI and acceptance band. G3 checks reciprocity only where its assumptions hold,
energy-flux balance, frame/phase conventions, mesh/time dispersion and source-
surface movement invariance for valid FW–H fixtures. The initial lane makes no
physical-validation claim for turbulence-generated noise, nonlinear acoustics,
psychoacoustic perception or combustion instability without independent data.

### 5.12 Degradation, Life, Manufacturing, and Reliability `[S→M]`

**State evolution.** The coupon-to-machine ladder begins with rainflow plus
Palmgren–Miner and named mean-stress corrections as bounded empirical
baselines; strain-life, multiaxial critical-plane/spectral fatigue,
Paris/NASGRO-family crack growth, creep-fatigue, corrosion-fatigue, fretting,
wear, insulation aging, demagnetization and lubricant degradation follow.
Higher-fidelity fracture uses cohesive, XFEM or phase-field lanes with crack
closure/contact, adaptive remeshing and conservative internal-state transfer.
No damage law becomes “material truth”: every one carries load-spectrum,
temperature/environment, geometry/notch, process, validity and calibration
receipts. Competing failure modes use dependency-aware system reliability and
never multiply independent probabilities unless independence is evidence.

**Manufactured and as-built reality.** Machine IR gains datum systems, GD&T,
surface texture, fits and correlated process tolerances; assembly order/path
feasibility; bolts/preload, welds, adhesives, keys, splines and interference
fits; and process→microstructure→residual-stress→property lineage for casting,
forging, machining, additive manufacture, heat treatment and coatings. Stable
`BodyId`, `SurfacePatchId`, `ContactFeatureId`, `StateSlotId` and lineage
morphisms survive split/merge/remesh/wear where unambiguous. Ambiguity returns a
typed refusal and invalidates caches/contacts/windings/adjoints instead of
silently rebinding them.

**Reliability and maintenance.** Correlated/hierarchical variability,
manufacturing modes, load/environment histories, model discrepancy and
inspection uncertainty feed FORM/SORM only as baselines, then adaptive subset/
importance/multilevel rare-event estimators and sequential inspection/
maintenance policies. Sampling measure, dependence model, likelihood or
importance weights and physical-population evidence are explicit; anytime-valid
statistics do not validate a guessed tolerance distribution. Outputs include
limit-state probability bounds, sensitivity, failure-mode attribution,
remaining-useful-life distributions and an inspect/repair/retire policy—not a
single false-precision “life.” Validation climbs coupon → feature/component →
subsystem → full machine, with blind holdouts at each promoted claim.

**No-claims.** Code-to-code comparison verifies implementation, not service
life. A standards formula earns `StandardConformance`; a prediction earns
`Validated` only for a named population, process, environment, load spectrum,
failure mode and QoI. Regulatory certification remains outside the evidence
color system.

### 5.13 Porous, Capillary, and Active-Material Couplings `[F→M]`

Material properties become load-bearing only when a formulation consumes them.
An `fs-flux`/`fs-gas` adapter ladder therefore adds anisotropic Darcy/Brinkman,
multicomponent Fick/Maxwell–Stefan, solution-diffusion membranes, adsorption/
sorption-swelling, interface resistance, and pressure/scale-dependent
Klinkenberg/Knudsen transitions where valid; poroelastic damage may evolve
permeability. Intrinsic porous `k [m²]`, fluid-dependent hydraulic
conductivity, basis-specific membrane permeability (diffusivity×solubility),
and thickness-normalized permeance are distinct types; a receipt binds
thickness, gas mixture, T/p/humidity and test method. Unsaturated/two-phase
extensions add relative permeability, capillary-pressure/saturation hysteresis,
Richards/two-phase Darcy and Biot poromechanics. Multicomponent diffusion declares barycentric mass or molar
reference frame, enforces the corresponding zero-sum flux constraint, handles
the Gibbs–Duhem/nullspace explicitly, uses chemical-potential/thermodynamic-
factor driving forces, and proves PSD entropy production on the `N_s-1`
independent-flux subspace. Counterdiffusion and any admitted Soret–Dufour
coupling are G0 fixtures. G1/G2 cover Darcy slabs, layered anisotropy, consolidation,
membrane breakthrough/permeation and direct species/mass balance.

Capillary/free-surface lanes consume surface tension and dynamic interface
cards: Young–Laplace/Young–Dupré are limited equilibrium oracles, while
Cox–Voinov, molecular-kinetic and resolved-contact-line rungs account for rate,
roughness, hysteresis, contamination and phase change. Mass, capillary work and
surface-free-energy receipts close separately.

Active-material lanes add piezoelectric electromechanics (open/short-circuit,
reciprocity and energy tests), magnetoelastic/magnetostrictive stress with the
same constitutive energy/dissipation and EM-force accounting convention, and
thermoelectric Seebeck/Peltier/Thomson coupling with Kelvin/Onsager–Casimir
checks under declared magnetic field/time-reversal parity and a
whole energy/entropy balance. Electrochemical battery/fuel-cell/electrolyzer
models form a later expansion portfolio over the same vector chemical,
electrical, thermal and degradation ports.

---

## 6. Coupling and Co-Simulation Architecture — `fs-couple` v2

The port bus becomes the machine-composition runtime:

1. **Typed relations and streams** per §3.7/Appendix A. Every subsystem exports
   public residual/port operators and declares storage, dissipation, sources,
   conserved quantities, causality, differentiability and DAE index.
2. **Windowed physical audits.** Integrate first-law power/enthalpy, mass,
   species/element/charge, momentum, entropy production, and optional exergy
   destruction over a closed accounting window. `fs-couple::EnergyAudit` is a
   useful balance-defect seed, not by itself a proof of subsystem passivity.
3. **Interface transfers.** `fs-couple` owns the neutral transfer schema,
   traits and audits; FEEC/CutFEM/EM/thermal/fluid adapters assemble concrete
   mortar/harmonic/Nitsche operators from their own traces, meshes and
   quadrature. Each separately proves compatibility, primal/dual adjointness,
   constant/flux preservation, inf-sup/GCL, and direct power balance. RBF is a
   fallback with an explicit defect and escalation trigger. Sheaf descent
   localizes compatibility defects but does not replace conservation identities.
4. **Partitioned and monolithic lanes.** Aitken → vector IQN-ILS, block Newton
   and monolithic residual/Jacobian escalation. Convergence of an interface
   iteration is not proof of correct coupled physics; residual, stability and
   solution-error artifacts remain distinct.
5. **Multirate co-simulation.** Waveform relaxation reports interpolation/
   extrapolation remainder, contractivity or energy-stability evidence,
   splitting error and subsystem error separately. If stability cannot be
   certified, shrink the window or assemble the strongly coupled pair.
6. **The machine graph is data.** A FrankenScript graph records subsystem,
   relation, interface, clock, controller, environment, audits and escalation
   policy. Admission rejects dimensional/sign/orientation gaps, algebraic loops
   without a solve policy, missing source closure and unaccounted state.
7. **Hybrid time is superdense and checkpointable.** Simultaneous resets use
   `SuperdenseTime { physical_time, microstep }`, not floating time plus an
   implicit order. `MachineStateSnapshot` binds graph/version, component state
   codecs, active modes, event-subdivision and solver frontiers, coupling/IQN
   history, circuit topology, controller clocks, RNG roots and accounting
   windows. Every component has an explicit state-schema migration/refusal path.
8. **Windows are transactional.** A coupled step is `propose → solve → audit →
   commit`. Cancellation requests drain workers and serialize the last accepted
   state plus resumable frontier; a half-solved or half-audited window is never
   externally committed. Pause/resume and migration replay to the declared
   determinism class.
9. **Stable identity is physical infrastructure.** Body, surface patch,
   feature, terminal, port and state-slot IDs are semantic—not array indices.
   Split/merge/remesh/wear/fracture events emit lineage morphisms. Ambiguous
   lineage invalidates dependent caches and returns a typed ambiguity instead
   of silently reconnecting a winding, contact, constraint or adjoint.

**Worked example — the ICE port graph (§8.5):**
crank `fs-mbd` ⇄ (torque·ω) ⇄ piston/rod/crank mechanism; cylinder `fs-gas`
0-D volumes ⇄ (p·dV/dt) ⇄ piston faces via `V(θ)`; valves = certified-event
orifices driven by `fs-kinematics` cam profiles; manifolds `fs-gas` 1-D ⇄ valve
ports; walls ⇄ (T·entropy flow) ⇄ `fs-thermal` network; bearing/skirt friction
⇄ `fs-tribo` ⇄ crank torque; optional `fs-em` + `fs-circuit` generator on the
crank for the genset capstone.

---

## 7. Solver and Infrastructure Upgrades (cross-cutting)

### 7.1 `fs-solver`
- **Newton–Krylov driver** with line-search/trust-region globalization,
  Eisenstat–Walker forcing, and the house resumability contract — needed by
  §5.2/§5.5/§5.7/§5.8 alike; one implementation, four clients.
- **Block-operator framework:** typed 2×2/3×3 block `LinearOp` composition with
  block preconditioners (extends the existing Stokes block-diag pattern into a
  reusable facility); real-equivalent complex form for time-harmonic EM.
- **Hiptmair–Xu (HX) preconditioner** for H(curl): auxiliary-space correction
  includes an edge-space smoother, scalar-gradient auxiliary solve, vector
  nodal `[H¹(Ω)]³` interpolation/auxiliary solve (or the required coordinate/
  constant-vector data in the admitted lowest-order construction), and a
  harmonic/cohomology coarse space. High order supplies its interpolation
  operator. Boundary, topology and coefficient-contrast assumptions are
  explicit; robustness is measured rather than inferred from the family name.
- **Eigensolver consumption:** L1 `fs-spectral` promotes `fs-la`'s
  Lanczos/LOBPCG behind the one generic resumable service. `fs-solver` consumes
  or re-exports it for tangent-constraint spectral gaps, drivetrain modes, and
  Orr–Sommerfeld work instead of growing a duplicate L3 implementation.

### 7.2 `fs-time`
- Sparse/matrix-free implicit integrators (generalized-α over `LinearOp`).
- SE(3) Lie-group integrators on motor states (the CONTRACT's named successor).
- `ValidatedStep`/`TrajectoryTube`, class-qualified complete root accounting,
  `EventLocator`/`ModeLedger`/Zeno-set handling (§3.4); the fast unvalidated
  dense-output lane remains explicitly separate.
- Implement L1 `fs-ad`'s generic protocols beside each primal step/reset:
  local VJPs and saltation/event-time derivatives for regular isolated
  transitions, generalized/set-valued records or no-claim at simultaneous,
  grazing, active-set and Zeno cases. L1 `fs-ad` supplies solver-agnostic
  schedules/records; L3 `fs-adjoint` composes and certifies them through generic
  interfaces without either crate depending on `fs-time`.

### 7.3 `fs-qty` — amount of substance is load-bearing
Ampere and kelvin already exist; EM/thermal units add aliases and parser tokens.
Chemistry, equilibrium constants, molar NASA reference data, stoichiometric
rates, electrochemistry and Faraday coupling cannot be dimensionally honest if
kilograms per mole are represented as an “ordinary” five-base quantity. Widen
`Dims` from `[m,kg,s,K,A]` to `[m,kg,s,K,A,mol]` atomically across the root and
standalone `fs-wasm` workspaces. This is an intentional source break, not a
source-compatible alias. Migrate every literal dimension vector and generic
arity, especially `fs-regime`'s 5×n Buckingham-π rank/nullspace/named-group
algebra, `fs-opt`, `fs-ir`, `fs-scenario::ir`, `fs-ledger`, `fs-package`,
`fs-crosswalk`, identity digests, parsers, JSON, tests and goldens. Versioned
five- and six-vector decoders map legacy values to `mol=0` under a semantic-
crosswalk receipt; old bytes/content hashes remain immutable and record
`old_hash→new_hash` rather than being overwritten. Candela remains omitted
until photometry is real.

Dimensional equality is still not semantic equality. Add noninterchangeable
quantity kinds/conversions for affine absolute temperature versus temperature
difference; angle/revolution and mechanical versus electrical angle/pole-pair
phase; angular velocity/rpm; torque versus energy; pressure versus stress;
strain and composition bases; mass/molar/concentration bases; instantaneous/
peak/RMS/phasor values; entropy versus heat capacity; and logarithmic acoustic
levels with reference pressure/power. This prevents dimensionally legal `2π`,
pole-pair, `√2`, Celsius-offset and basis errors.

`SpeciesId`, `ElementId`, elemental matrix `A`, stoichiometric matrix `N`, and
charge vector `z` are immutable artifacts; G0 proves `A N = 0` and `zᵀN = 0`.
Solvers may store dimensionless mass fractions `Y_k`, but mass↔amount
conversion always uses typed molar mass and records the reference basis.

### 7.4 Determinism, cancellation, and performance program for the new kernels
Every new hot kernel is a tile program under `Cx` with deterministic reduction
trees, per house rules. Initial roofline target families (stated so they can
be failed, machine-fingerprinted per `fs-roofline`):

| Kernel | Target family |
|---|---|
| Weighted curl-curl apply | Baseline against its matrix-free roofline and an equivalent assembled sparse operator at equal error; ratify a hardware-specific target only after both are measured |
| Contact broad+narrow phase, 10⁴ moving pairs | bounded per-step time with deterministic ordering; measured lane, no absolute claim until baselined |
| 1-D gas FVM | baseline-relative updates/s at fixed model, grid, CFL, error and ISA; target set only after measurement |
| 0-D engine cycle | baseline-relative cycles/s at fixed species/model/tolerance; campaign throughput target set after profiling |
| MBD Newton step, 10² bodies/10² constraints | dominated by sparse KKT solve; measured lane |
| LBM moving-boundary overhead | Provisional E7 target ≤ 15% over the static-boundary lane at matched accuracy and lattice; ratify or replace from measured receipts before making it an acceptance gate |

Baselines are recorded via `fs-roofline` receipts before any optimization
story is told (house rule: performance folklore is a failure mode).

### 7.5 Optimization and UQ tie-ins (why this is all worth it)
The mission is *synthesis*, not just simulation. Every domain plan above ships
classical, generalized, or derivative-free sensitivity lanes with declared
regimes so that, at the end of the roadmap, these are
well-posed FrankenSim studies: maximize gear-train efficiency subject to an
evidence-typed transmission-error band; minimize motor cogging torque over magnet
shape with a force/torque error budget at the claimed rung; maximize ICE brake efficiency over cam
timing + compression ratio under a knock-integral chance constraint with an
anytime-valid statistical bound under declared sampling and model assumptions;
minimize Geneva indexing error CVaR under manufacturing tolerance
distributions. Existing UQ/tolerance/assimilation/OED/robustness crates are
useful seeds but do **not yet close** the digital-thread loop. The explicit
§4.2 upgrades add correlated/nonlinear/hierarchical UQ and GD&T, rare events,
global sensitivity, identifiability/discrepancy, nonlinear/non-Gaussian
assimilation, and 3-D correspondence-uncertain as-built registration. Together
they connect material coupon and dynamometer data to posterior
model discrepancy, controller tuning, manufacturing allocation, health
monitoring and remaining-life decisions. ASCENT consumes these typed artifacts
rather than assuming every nonsmooth machine objective has a smooth adjoint.

### 7.6 Operational V&V, calibration, and prediction assessment

`ContextOfUse`, `ValidationPlan`, `ExperimentArtifact`, `CalibrationSplit`,
`SolutionVerificationReceipt`, and `PredictionAssessment` are machine-readable
artifacts. They require calibration, validation and blind holdout partitions;
instrument calibration, synchronized clocks, repeatability/covariance and data
authenticity; observability/identifiability/confounding diagnostics and inverse-
crime prevention; mesh/time/nonlinear/iterative uncertainty; model-form,
parameter, numerical, data, aleatory and epistemic semantics; validation
metrics that include experimental and numerical uncertainty; posterior-
predictive checks; and an applicability/extrapolation decision. Evidence from a
high-fidelity code, reduced model, or second implementation is verification or
discrepancy evidence—not physical validation without an experiment.

The claim dependency graph is QoI-specific. Numerical/data/model errors may
form a declared bound or probabilistic budget waterfall; categorical evidence
colors never become percentages or an averaged “confidence score.” ASME V&V,
JCGM GUM and NASA-STD-7009-family workflows inform the schemas, but conformance
to a process standard remains a separate receipt from model validation.

### 7.7 Adaptivity, scale, and heterogeneous execution

Dynamic hp/anisotropic refinement, CutFEM cut-cell stabilization, moving-mesh
quality, conservative remap, wear/fracture/topology-change state transfer and
DWR/QoI targeting are a cross-domain program, not a one-time meshing step.
Refinement receipts bind topology lineage, balance defect, projection error and
whether the requested QoI bound actually decreased.

Performance gates use a `ScaleQualification` matrix on Apple Silicon and
many-core EPYC/Threadripper with: DOFs/cells/species/bodies/contact features and
active contacts/coupling windows; memory ceiling/checkpoint size; strong/weak
scaling and NUMA efficiency; deterministic-mode overhead; cancellation latency
and restart cost; nonlinear/Krylov robustness; and throughput at matched QoI
error. HIL adds tail latency, jitter and worst-case execution, never just mean
throughput. Every accepted number pins hardware, scale, partition, accuracy,
budget and fallback.

The first production scale is single-node NUMA with topology-aware partitions.
A deterministic multi-node domain-decomposition and checkpoint/migration track
is funded behind `[M]`, including partition-independence, loss/retry and drain
proofs. An accelerator lane is admitted only by an explicit constitutional
dependency-policy decision and audited safe substrate; production FFI is not
smuggled in under a benchmark. If the policy prevents the required capability,
the plan reports that capability unavailable and funds a Franken-native
substrate program.

### 7.8 Safety, EMC, and assurance `[S→M]`

Scientific certificates are executable evidence, **not regulatory product
certification**. Machine IR nevertheless owns `HazardId`, typed safety
requirements, operating envelopes, fault-containment regions and safety-case
links. The assurance ladder includes FMEA/FMECA, dependency-aware fault trees/
common-cause models, reachability-backed unsafe-state witnesses, fault
injection, runtime monitors, degraded control and emergency shutdown. Program
targets include rotor burst/overspeed and mechanical containment; pressure
vessel/line rupture; runaway combustion, fire/thermal propagation and toxic
emissions; shock, creepage/clearance and dielectric breakdown; partial
discharge/insulation aging and bearing currents; watchdog/timing/sensor-
plausibility faults; and safe cancellation/restart.

EMC is a coupled machine program: differential/common-mode parasitics,
cable/harness and grounding/bonding models, inverter edge spectra, conducted/
radiated emissions and immunity, shaft/bearing currents and sensor/control
upsets. Exact standard editions and lab setups become `StandardConformance` or
experimental artifacts; simulation alone does not claim legal compliance.
Every hazard assumption has an owner, monitor, violation effect and expiry.

### 7.9 Engineer workflow, interoperability, and competitive scoreboard

The end-to-end workflow is explicit:

`import/create geometry → assemble/joint → bind materials/interfaces/as-built
state → declare IC/BC, drive/load cycle, tolerances, faults and hazards → choose
ContextOfUse/QoIs/acceptance/budgets → admit fidelity/solver/escalation policy →
calibrate on training data → verify/validate on holdouts → simulate/UQ/optimize →
export replay, evidence and decision package`.

The canonical interchange is FrankenSim-native and content-addressed. Pure-Rust
FMI 3.0.2/SSP 2.0 parsing/export or quarantined isolated-process adapters may
support adoption; foreign FMU execution is outside the trusted production
dependency graph and its outputs cannot inherit native certificate status.
Round-trip batteries cover FrankenScript → machine graph → scenario/domain
artifacts → package/replay, preserving stable IDs and semantic hashes.

`CompetitiveCapabilityLedger` plus a license-respecting black-box suite compares
named versions of alternatives on QoI accuracy, conservation/event misses,
modeling effort and diagnostics, coupling stability, optimization/UQ throughput,
reproducibility/provenance, evidence transparency, scale and cost. The leapfrog
target is unmatched **certified whole-machine composition** while also winning
selected fidelity/performance decks—not an untestable assertion of universal
peak fidelity on day one.

---

## 8. Flagship Demonstrators — Forcing Functions, Not Demos

The names below are suite/module IDs inside the existing L6
`fs-flagship-e2e` harness, not five additional crates unless measured build or
ownership pressure later justifies a split. Five dependency-ordered core
flagships plus an expansion portfolio each follow the proven campaign pattern
(parameterize → e-race screen → refine → certify → atlas → content-addressed
notebook) with a reduced `fs-wasm` browser companion. Each names "what breaks
first," per house style.

### 8.1 `fs-geneva-e2e` — the certified intermittent-motion machine
**Composition:** `fs-kinematics` + `fs-mbd` + `fs-contact` + `fs-tribo` +
existing UQ/tolerance crates + `fs-matdb`. **Deliverables:** certified mode
sequence for full revolutions (every
engagement/disengagement event enclosed — §3.4); the **dwell-retention
certificate** (an ideal-lock subreceipt with first-order form closure plus
`N_V`/restricted-`J̃` margin only for a declared
ideal zero-clearance active boundary mode; for real clearance, a certified
dwell-play reachable set and non-escape/invariance result over load, impulse,
friction, elasticity, wear, tolerance and mode uncertainty, with boundary-mode
margins—§3.2); indexing-accuracy
distribution under pin/slot tolerance UQ with e-process stopping; peak Hertzian
contact stress against a named material/contact-strength criterion; Archard
wear-life estimate (`Estimated`, honest);
atlas over (station count, pin radius, crescent clearance) illuminating the
accuracy/stress Pareto. **What breaks first:** contact chattering at engagement
→ Moreau–Jean step-size control + restitution law selection; certified event
isolation cost near grazing engagement → `Unknown`-with-interval refusals kept
visible in the report.

### 8.2 `fs-gear-e2e` — the certified gear train
**Composition:** `fs-kinematics` geometry + `fs-machine` + `fs-symmetry` +
`fs-solid` tooth/shaft/housing FE + `fs-tribo` + existing UQ/tolerance crates.
**Deliverables:** ratio/common-normal certificate; unloaded, static-loaded and
dynamic TE; microgeometry, misalignment, mesh stiffness and support
flexibility; bending/contact fatigue, scuffing and micropitting model cards;
mixed/EHL lubrication, thermal distortion, backlash UQ, efficiency, vibration
and acoustic order maps; planetary closure checked against Willis. **What
breaks first:** mesh-stiffness FE cost across roll angles →
sector symmetry + surrogate with conformal bands (certify-or-escalate exists).

### 8.3 `fs-motor-e2e` — PMSM/BLDC motor and generator
**Composition:** `fs-em` 2-D/2.5-D/3-D moving-conductor ladder + audited
sliding interface + `fs-symmetry` + `fs-circuit`/`fs-power` inverter +
`fs-control` FOC + `fs-mbd`/rotordynamics shaft + `fs-thermal` network +
`fs-matdb` (M-19, NdFeB, copper). **Deliverables:** cogging-torque waveform
with a dual-method (Maxwell stress vs. virtual work) consistency receipt;
back-EMF/flux linkage via filament and distributed-winding checks (§3.6);
torque–speed–efficiency and switching/loss maps; demagnetization margin,
thermal derating, rotor critical speeds, closed-loop load steps, sensor and
open-switch faults; generator quadrant; conductor-frame Ohmic/Lorentz/Poynting
closure; `LossOwnershipReceipt`s preventing double counting between resolved
eddy/hysteresis/winding/switching physics and empirical loss maps; insulation,
bearing-current and EMC hazard artifacts at the admitted rung;
atlas over (magnet arc, slot opening, turns). **What breaks first:** air-gap
transfer conservation at coarse interface resolution → direct element/dual-
transfer/GCL balance receipts detect and localize it, while the sheaf layer
localizes trace-compatibility defects; cogging accuracy demands fine discretization → sector symmetry +
p-refinement of the existing high-order spaces.

### 8.4 `fs-wankel-e2e` — the trochoidal rotary-engine capstone
**Composition:** `fs-kinematics` trochoid kit (separate certified ideal apex
path, finite seal/contact locus, actual-bore envelope, and rotor-flank design — §3.3)
+ `fs-mbd` (eccentric shaft + phasing gear as a gear-pair constraint) +
`fs-contact` (apex/side/corner-seal sliding kinematics) + `fs-gas` Tier 1
(three chambers with closed-form `V(θ)` as G1 oracle) + `fs-thermal` +
`fs-tribo` (sealing, friction, wear, leakage, and durability are historically
prominent Wankel bottlenecks, now first-class simulated quantities with
hardness-informed Archard life).
**Deliverables:** certified housing/rotor/seal geometry and phasing; chamber
pressure/species traces and p–V loops per face; port timing; indicated/brake
torque; crevice and inter-chamber blow-by, oil film, seal leakage/friction/
wear, housing/rotor thermal deformation, wall heat, ignition and emissions;
knock-
onset check via Livengood–Wu `[F]`; atlas over (eccentricity ratio K, port
timing, compression ratio). **What breaks first:** apex-seal contact under
combustion pressure is stiff and nonsmooth → IPC lane for design studies,
Moreau–Jean for hard-impact validation runs; leakage between chambers uses
an explicitly calibrated reduced correlation on the first rung, with resolved
gap/CHT escalation in the later machine-fluid path.

### 8.5 `fs-genset-e2e` — the grand synthesis: ICE + generator
**Composition:** the mature core machine stack, explicitly—not an implicit
“everything” dependency. Single-cylinder four-stroke: `fs-kinematics` cam
profiles + valve-train kinematics; `fs-mbd` crank–rod–piston with bearing
friction (`fs-tribo`); `fs-gas` 0-D cylinder (Wiebe + Woschni) + 1-D manifolds
(HLLC baseline plus separately entropy-stable lane); `fs-thermal` component network; crank port coupled to
the `fs-motor-e2e` machine running as a generator into an `fs-circuit`
rectifier + load. **Deliverables:** full-cycle p–V diagram with a numerical
enclosure on indicated work plus separately typed parameter and model-form
uncertainty (Tier-1 interval lane); brake torque/efficiency
maps with every correlation's categorical evidence attached to its exact
QoI-dependency DAG and a separate numerical/model/data budget waterfall; transient
start-up and load-step (multirate co-simulation exercising coupling item 5 in §6); electrical
output waveform with certified commutation events; the report's headline is
the **evidence dependency graph of a complete machine**—which claims are
Verified (eligible geometry/kinematics propositions, discrete balance
identities, and event-capture statements under their assumptions), which
Validated (only correlations/losses with named held-out datasets and QoIs), and
which Estimated (unvalidated Wiebe/Woschni, wear, acoustics). Claim-level honesty
is a competitive target to be demonstrated by a versioned comparative workflow
study, not asserted from aspiration. **What breaks first:** multirate stability between
stiff combustion and slow thermal → waveform-relaxation windowing with separate
interpolation, splitting, subsystem-solution, and contractivity/energy-stability
receipts; shrink the window or escalate to the monolithic pair when stability
cannot be shown. Correlation validity-domain violations at
extreme operating points → automatic `Validated → Estimated` demotion (already
implemented in `fs-evidence::regime_demotion`) rather than silent wrongness.

### 8.6 Expansion flagships — the leapfrog portfolio

- **`fs-induction-drive-e2e`:** squirrel-cage moving-conductor eddy currents,
  switched inverter, FOC, rotor thermal state, broken-bar/eccentricity faults,
  torque ripple, acoustic orders and pull-out/stall/slip-instability atlas
  (converter/grid synchronization is named separately if present).
- **`fs-constant-width-e2e`:** a genuinely rolling Reuleaux/constant-width body
  with Pfaffian no-slip branch, ideal reaction-zero-work receipt,
  finite-friction feasibility and certified slip/loss-of-contact transitions,
  contact-branch event ledger, caging/swept envelope, pressure/contact-work and
  tolerance UQ. This is an explicit falsifier against conflating rolling with
  the Wankel's sliding apex seals.
- **`fs-pump-bearing-e2e`:** motor-driven pump/gerotor with cavitation,
  fluid-film bearings, EHL/mixed lubrication, seal leakage, rotor critical
  speeds, thermal growth, vibration and remaining-life UQ.
- **`fs-turbo-efuel-e2e`:** turbocharged/e-fuel ICE with compressor/turbine
  maps, real-gas/fuel evaporation, knock/emissions/aftertreatment, oil/coolant
  loops and generator load transients.
- **`fs-digital-twin-e2e`:** ingest coupon, metrology, sensor and dynamometer
  data through existing as-built/assimilation/OED crates; identify correlated
  material/interface/model discrepancy; update posterior predictions,
  controller settings and maintenance decisions without laundering validation.

---

## 9. The Gauntlet Program for New Domains

House rule restated: every solver arrives with G0–G5 or documents exactly what
is proven and what is not. Domain-specific batteries are named in §5; this
section adds the cross-domain apparatus.

**New canonical verification/validation registry (G1 analytic and G2 benchmark
entries tagged individually), stored with exact version, source, license, input
deck hash, QoIs and acceptance envelopes:**
G1 entries include Hertz contact closed forms, bouncing-ball impact maps,
block-on-incline stick/slip, Bennett linkage mobility, involute constant-ratio,
and Geneva closed-form
geometry, epitrochoid closed form, coaxial cable/sphere-in-field, Stefan,
Riemann/nozzle, and air-standard-cycle oracles. G2 entries include exactly
versioned IFToMM rectangular-Bricard and Andrews-squeezer problems, TEAM
10/13/20/24/30a/30b, NAFEMS reports, CFR/engine pressure traces, and
oscillating-cylinder and gerotor flow (moving-boundary LBM). A family name such
as TEAM, NAFEMS or CFR is not an executable benchmark until those fields and an
independently reviewed expected result are present.

**Certifying the certifiers (master-plan/addendum falsifier discipline):**
- Mobility certificates vs. randomized Jacobian rank probes and Grübler
  falsifiers (disagreement is Sev-0).
- Event soundness derives from validated true-flow and root-count proofs;
  independent high-resolution/event-time integrations are falsifiers only,
  including grazing/simultaneous families where `Unknown` is correct.
- Swept/envelope soundness derives from interval set inclusion and exhaustive
  subdivision certificates; Monte Carlo/adversarial sampling is a falsifier,
  never the proof.
- Conjugate-action certificates vs. independent numeric kinematics (rigid
  contact simulation of the same pair).
- Entropy/power audits vs. analytic balances and a disjoint implementation;
  double-double replay is a useful precision check but not independence.
- Interface compatibility vs. sheaf/descent localization, and conservation vs.
  a disjoint global flux/power balance plus primal/dual/GCL identities.

**Sensitivity gates:** every classical-adjoint lane merges only with
adjoint-vs-AD-independent FD agreement at regular G1 fixtures; complex-step is
used only for holomorphic, branch-free implementations that preserve complex
perturbations. Hybrid, contact and
nonsmooth paths additionally test saltation/generalized derivatives and refuse
where regularity assumptions fail.

**V&V separation:** every flagship report has independent axes for code
verification, solution verification, numerical uncertainty, parameter/data
uncertainty, model-form validation, prediction-domain relevance and comparison
to experiment. “Validated” always names the experiment/dataset and QoI;
“Verified” always names the mathematical/numerical claim.

---

## 10. Phase Roadmap

Phases are capability ratchets, not dates. Each contains small vertical slices;
no phase gate requires every research branch to succeed simultaneously. Every
exit item expands to an executable gate naming input deck, QoI, acceptance band,
evidence class, problem scale/hardware, time/memory budget, cancellation bound,
determinism class and fallback; broad prose alone cannot close a phase. The
machine-readable prerequisite DAG and critical path live beside the manifest
fixture, while theorem cards start with their enabling phases.

| Phase | Scope | Exit criteria (all Gauntlet-gated) |
|---|---|---|
| **E0a — Units, identity, schema** | six-base/semantic units; stable machine IDs/lineage; scenario/machine IR; ledger/package/crosswalk and standalone-fs-wasm migration | five↔six wire decoding and immutable old/new hashes; all literal/generic consumers migrated; FrankenScript→scenario→package stable-ID/hash round trip; ambiguous lineage refuses |
| **E0b — Operators and solvers** | weighted/coupled constitutive operators; `ConstitutiveGraph`; nonlinear/block solver; neutral versioned ports/transfers | heat/current manufactured order/units plus free-energy/dissipation/cross-coupling checks; manifest fixture compiles without cycles; scalar-port migration and power-pairing battery |
| **E0c — Data and V&V registry** | immutable observation/claim/model cards; offline pack compiler; typed validity/evidence; benchmark and V&V registry | one complete material-query dependency closure with source/license/covariance; blind calibration/validation split; intentionally invalid pack/query refused |
| **E0d — Replayed vertical slice** | one minimal admitted **machine-graph skeleton** using typed toy rotational/electrical/thermal storage, loss and source elements—not a physical motor claim—plus transactional snapshot/cancellation and report | incomplete graph rejected; admitted graph runs, audits, checkpoints, cancels/drains/resumes and bit-replays at the declared scale; the graph is later rebound to real E1–E5 domain artifacts without identity drift; initial thermodynamics/identity theorem cards open |
| **E1 — Prescribed motion and nonlinear kinematics** | oriented 2-D complex, `CertifiedMotorTube`/`SpacetimeChart`, sweeps/envelopes, configuration/holonomy graph, tangent/dual rigidity, `GeometricMoments`, exact-path events | derived trochoid/involute and slider-crank/four-bar oracles; two-sided clearance range; Bennett/Bricard finite/tangent checks; prescribed-path isolable-guard root accounting; conjugate/mechanism theorem cards active |
| **E2 — Dynamic interaction** | rigid `fs-mbd`; holonomic and nonholonomic validated dynamics; capability-routed CCD; smooth/nonsmooth contact; dry `fs-tribo` baseline (Coulomb, Hertz, thermal partition, Archard); robust Geneva | constrained/rolling dynamics decks; conditional nonintersection; global multicontact energy/work and complementarity receipts; ideal active-mode first-order form closure **and** real-clearance dwell-play/non-escape; `fs-geneva-e2e` |
| **E3 — Machine elements, flexibility, life, NVH** | spatial gear families, bearings, shafts, seals/fasteners/joints; component modes/Cosserat; advanced lubrication/EHL/cavitation/flash temperature; fatigue/fracture/wear; acoustics/rotordynamics | loaded TE and mesh/bearing validation; Campbell/Floquet/acoustic fixtures; lubrication/thermal/loss ownership; coupon→component life holdout; `fs-gear-e2e` |
| **E4 — Electromechanical drive stack** | 2-D→2.5-D→3-D moving-conductor EM, audited sliding interface, thermal, descriptor circuits, power electronics and control | exact-version EM/thermal/converter benchmarks; gauge/energy/GCL checks; closed-loop inverter-fed motor with thermal derating and faults; `fs-motor-e2e` |
| **E5 — Reactive-flow machine stack** | typed thermochemistry, conservative 0-D/1-D gas, EOS-specific entropy/invariant-domain lane, porous/capillary transport, moving-boundary `fs-flux`/LBM, CHT, acoustics and seal leakage | thermo/species/detailed-balance laws; Riemann/nozzle/engine/acoustic validation QoIs; fluid/body impulse-torque-work and GCL; Wankel geometry/volume/seal/blow-by/thermal receipts; `fs-wankel-e2e` |
| **E6 — High-fidelity escalation** | dynamic hp/adaptivity; real-gas/multiphase/sprays, multi-D reacting ALE/CutFEM, emissions/aftertreatment, advanced hysteresis/radiation/FSI; single-node scale and multi-node decision gate | reduced↔high-fidelity discrepancy/escalation validated on at least two machines; ScaleQualification at matched QoI error; at least one reproducible 3-D reacting moving-domain deck closes mass/element/energy/GCL and its named QoI band, so 3-D combustion is an admitted measured rung rather than a permanent exclusion |
| **E7 — Whole-machine synthesis** | genset, induction, constant-width, pump-bearing and turbo-e-fuel flagships; safety/EMC cases; mature campaigns; upgraded correlated UQ/GD&T, nonlinear assimilation/as-built, rare-event reliability, OED/optimization/interoperability | closed accounting/stability/error/loss-ownership windows; blind experimental validation; one safety/fault-containment case; FMI/SSP quarantine conformance; a versioned `CompetitiveCapabilityLedger` on named whole-machine decks; at least two robust design improvements and one posterior digital-twin update |
| **E8 — Theorem Foundry integration summit** | independent reproduction and formal adjudication/promotion attempts for all §3.11 cards; nonlinear stacks, multisymplectic composition, set-valued hybrid reachability and adaptive model-category routing | every card has exact quantifiers, TCB, baseline, counterexample generator, checker, retained falsifier and reproduction pack; at least three independent cards spanning (i) topology/coupling, (ii) mechanisms/hybrids/contact and (iii) thermodynamics/fidelity/life are both `KernelChecked` and `ScaleQualified`; failed moonshots remain versioned/refuted and do not block solid lanes or become hidden prerequisites |

Parallel tracks may start after their actual prerequisites, but one unproven
research mechanism is isolated per proof lane. E6 high-fidelity promotions and
the E8 Foundry summit are not hidden prerequisites for a solid reduced-fidelity
flagship; only claims that consume those rungs depend on their gates, while an
unsuccessful moonshot remains visible and unavailable. In particular the motor
flagship cannot close before its thermal/power/control dependencies, and the
Wankel cannot close before contact, thermal, leakage and thermochemistry.

---

## 11. EV Opportunity Matrix and Recommendation Cards

Per the alien-graveyard selection discipline, score
`(Impact × Confidence × Reuse)/(Effort × AdoptionFriction)` on 1–5 inputs.
This exposes “brilliant but unadoptable” ideas instead of hiding adoption cost.

| Lever | I | C | R | E | A | Score | Phase |
|---|---:|---:|---:|---:|---:|---:|---|
| R1 Weighted constitutive operators | 5 | 5 | 5 | 2 | 1 | **62.5** | E0 |
| R2 Typed material/interface data + receipts | 5 | 5 | 5 | 3 | 2 | **20.8** | E0 |
| R3 HX/AMS + block preconditioners | 4 | 4 | 4 | 3 | 2 | 10.7 | E4 |
| R4 Certified spacetime/swept enclosures | 5 | 4 | 4 | 3 | 2 | 13.3 | E1 |
| R5 Validated true-flow event capture | 5 | 3 | 5 | 4 | 2 | 9.4 | E2 |
| R6 Nonlinear/tangent kinematic complex | 4 | 4 | 4 | 3 | 2 | 10.7 | E1 |
| R7 SE(3) rigid→flexible MBD ladder | 5 | 4 | 5 | 5 | 2 | 10.0 | E2–E3 |
| R8 Conservative 0-D/1-D reactive stack | 5 | 4 | 4 | 5 | 2 | 8.0 | E5 |
| R9 Body-body contact + rigorous CCD | 5 | 3 | 5 | 5 | 2 | 7.5 | E2 |
| R10 Nonsmooth complementarity lane | 4 | 3 | 3 | 4 | 2 | 4.5 | E2 |
| R11 Equivariant conservative interface | 5 | 3 | 5 | 4 | 2 | 9.4 | E4 |
| R12 2-D→3-D moving-conductor EM | 5 | 4 | 4 | 5 | 2 | 8.0 | E4 |
| R13 Topological/distributed winding coupling | 4 | 4 | 4 | 3 | 2 | 10.7 | E4 |
| R14 Entropy/invariant-domain flux ledger | 4 | 4 | 4 | 4 | 2 | 8.0 | E5 |
| R15 Port-thermodynamic + multirate/monolithic bus | 5 | 4 | 5 | 4 | 2 | 12.5 | E0–E7 |
| R16 Add amount-of-substance dimension now | 5 | 5 | 5 | 2 | 1 | **62.5** | E0 |
| R17 Power-electronics/control/bifurcation stack | 5 | 4 | 4 | 4 | 2 | 10.0 | E4 |

Before conversion to Beads, **every** row gets a full recommendation card:
failure signature, change, adoption wedge, baseline comparator, primary-paper
status, budgeted mode, cost p50/p95/p99, assumptions, calibration trigger,
fallback/rollback, interference tests, verification and validation artifacts,
reproduction pack, demo linkage, activation criterion and kill criterion.
Three representative cards follow:

**Card R1 — Weighted constitutive operators.** *Failure/change:* B5/B6 leave
complete FE spaces without material-aware field assembly; add typed
`WeightedMass/Star` plus constitutive-map contracts. *Adoption wedge:* linear
heat/current flow on existing meshes. *Baseline:* existing scalar assembly and
an isolated disjoint FEM oracle. *Primary-paper status:* FEEC stability and HX
primary sources are seeded in the registry and their relevant results identified;
each formulation-specific derivation and proof must still be reproduced. *Budgeted
mode:* scalar cell weights; tensor/quadrature/nonlinear disabled. *Risk and
trigger:* loss of coercivity or Newton stall activates load stepping/Picard and
refuses the nonlinear claim. *Cost:* p50/p95/p99 measured only after the heat
slice; no invented number. *Verification:* units, coupled free-energy/block
stability, PSD dissipative part, declared reversible skew/Onsager–Casimir
structure, consistent Jacobian, quadrature, commuting-space convergence,
gauge/BC and MMS;
`dd=0` is a separate topology test. *Validation/demo:* exact-version thermal/EM
benchmarks, E0 heat slice then motor. *Interference/rollback:* run unweighted
goldens and coefficient-jump cases; feature-gated assembler can be disabled
without schema loss. *Activation/kill:* activate tensor/nonlinear by G1 order
and robustness; kill a formulation that cannot obtain its stability artifact.

**Card R5 — Validated true-flow event capture.** *Failure/change:* the live
time crate has no event location, and roots of an interpolant are not roots of
the true flow; add `ValidatedStep`, class-qualified complete root accounting, resets and
`ModeLedger`. *Adoption wedge/budget:* exact prescribed paths first, then small
ODEs; cap subdivision and return `Unknown`. *Baseline:* classical dense-output
polish, plus disjoint fine integration as falsifier. *Risk/trigger/fallback:*
wrapping blow-up or DAE irregularity shrinks steps, changes enclosure method,
or falls back to `Estimated`; simultaneous/grazing modes return sets. *Cost:*
p50/p95/p99 subdivision counts and certificate rate from the adversarial suite.
*Verification:* true-flow inclusion, root existence/uniqueness/count over finite
isolable guard families, reset coverage and honest flat/grazing/simultaneous/
Zeno alternatives. *Demo:* Geneva before valves and
switches. *Activation/kill:* “Verified capture” activates only with complete
coverage; kill any shortcut that certifies only dense output.

**Card R6 — Nonlinear/tangent kinematic complex.** *Failure/change:* B4/B11
lack mechanism reasoning; add finite constraint holonomy, `J`, `ker J`,
`ker(Jᵀ)`, second-order/continuation and cone certificates, using a sheaf only
when constructed. *Adoption wedge:* planar joint catalog and four-bar.
*Baseline:* Grübler, symbolic kinematics and numerical rank/SVD. *Risk and
fallback:* near-singular rank or branch ambiguity returns interval `Unknown`
and numerical SVD `Estimated`; no false `H¹` label. *Budgeted mode:* tangent
report without global continuation. *Cost:* p50/p95/p99 factorization and
subdivision receipts. *Verification:* frame equivariance, `J=DC`, primal/dual
work, Maxwell–Calladine, Bennett/Bricard finite continuation, ideal active-mode
Geneva first-order form closure, and real-clearance dwell-play/non-escape.
*Interference:* gauge/scaling/multiplicity tests; L1 `fs-spectral` stays
generic. *Activation/kill:* activate nonlinear sheaf/stack only when it beats
the plain constraint graph on localization or composition; otherwise retain
the simpler artifact.

---

## 12. Proof Obligations, Assumptions Ledger, Failure Modes

### 12.1 Proof obligations (executable, per the house "certificates over vibes")

| ID | Obligation | Mechanism |
|---|---|---|
| PO-1 | Incidence exactness is independent of material. Constitutive blocks prove units/objectivity/material symmetry as applicable, coupled free-energy Hessian/block stability, PSD symmetric dissipative part, reversible skew/gyroscopic structure, Onsager–Casimir/Kelvin relations with declared time-reversal parity, tangent consistency, quadrature and coefficient robustness; not every map is individually symmetric or monotone | algebra + spectrum/majorant + MMS |
| PO-2 | Conservative junctions have zero power defect; storage/dissipation/source/stream elements close energy, momentum/angular momentum where applicable, mass, element, charge and open-system entropy over one accounting chart. Loss IDs are disjoint and each thermal credit equals its upstream debit | runtime receipts + `LossOwnershipReceipt` + disjoint recomputation |
| PO-3 | Smooth variational MBD earns only its conditional backward-error/symmetry statement. Position, velocity, acceleration and Pfaffian residuals plus state-distance bounds are distinguished; every projection reports impulse, energy/momentum/angular-momentum/symmetry/adjoint defects; nonholonomic reactions satisfy declared virtual-power and friction-feasibility laws | conditional theorem + residual/state/projection/rolling receipts |
| PO-4 | Complete event coverage requires a compact phase domain, regular/index-fixed modes, true-flow tubes, a finite isolable guard family, reset closure, and no accumulation or a proved Zeno limit; flat/black-box/inclusion cases return the typed alternative | `ValidatedStep` + class-specific root-count/coverage proof; scans only falsify |
| PO-5 | Swept-union inner/outer enclosures are sound under declared sign semantics; relevant envelope-boundary branches are sound **and complete** across endpoints/rank loss/trimming/visibility/self-intersection; conjugacy additionally covers common-normal/contact/no-interference conditions | interval inclusion/subdivision/continuation proof; adversarial samples only falsify |
| PO-6 | The fully discrete reactive moving-mesh lane uses an EOS/mixture-specific convex entropy pair, preserves density/internal-energy/species admissibility and pressure equilibrium where required, satisfies chemistry detailed-balance/entropy production, GCL/source/boundary terms and nonlinear-solve defect; phase boundaries use a separate free-energy inclusion | per-step/source ledger + refinement battery |
| PO-7 | Evidence closure contains every load-bearing material/model/correlation receipt and excludes irrelevant ingredients; verification and validation axes cannot launder one another | receipt re-execution + dependency mutation tests |
| PO-8 | `⟨Av,w⟩_W=⟨v,A*w⟩_V` for new linear transfers, with `A*=M_V⁻¹AᵀM_W` in real spaces or `M_V⁻¹AᴴM_W` in complex sesquilinear spaces and any orientation/peak-RMS sign/factor declared, plus verified discrete sensitivities for regular nonlinear paths; nonsmooth exceptions refuse | weighted-adjoint/VJP gates + saltation/generalized tests |
| PO-9 | Sliding/nonmatching interfaces satisfy compatibility, signed primal/dual adjointness, constant/flux preservation, inf-sup where relevant, direct power balance and moving-interface GCL | separate interface receipt fields |
| PO-10 | Tangent mobility/self-stress/index reports agree with independent rank/duality checks or return `Unknown`; finite mobility needs second-order/continuation evidence | kinematic falsifier suite |
| PO-11 | Reactions conserve elements/charge and satisfy reverse-rate/equilibrium thermodynamics. Diffusion declares reference frame, zero-sum constraint/nullspace, thermodynamic driving force and PSD entropy production on the independent-flux subspace. Mass/amount bases and ideal-gas `cp-cv` identity use the correct frozen-composition molar or mass-specific convention | exact `AN=0`, `zᵀN=0`, detailed-balance/Onsager and thermo identities |
| PO-12 | Impulse-free circuit transitions preserve the admitted state relation; inconsistent switches solve distributional MNA under declared boundary balances with an energy-defect receipt, consistently initialize the new descriptor state, regularize explicitly, or refuse—individual charge/flux continuity is not universal | structural MNA + impulse/switch batteries |
| PO-13 | EM passes gauge/source/winding/moving-frame checks; forces name stress/material convention, held electrical variables, origin/frame and numerical/geometry/material bounds; moving conductors close field/source, Lorentz-mechanical and conductor-frame-Ohmic power | field/circuit/force/Poynting batteries |
| PO-14 | Moving meshes/fresh cells satisfy GCL and mass/momentum/energy remap as applicable, plus equal-and-opposite fluid/body impulse, torque and boundary-work/power receipts | ALE/LBM direct audits |
| PO-15 | Material/interface queries are unit/frame/domain/definition correct, covariance preserving and provenance/license complete. Intrinsic porous permeability, hydraulic conductivity, basis-specific membrane permeability and thickness-normalized permeance cannot crosswalk without definition/test receipts | schema/property mutation battery |
| PO-16 | Performance gates bind model, resolution, error/QoI, machine fingerprint and baseline; no naked throughput number is a claim | `fs-roofline` receipt schema |
| PO-17 | A §3.11 theorem can strengthen a production claim only when its exact statement/assumptions are machine-readable, its formal/executable checker reproduces, and adversarial representative changes/counterexamples fail to break it | `TheoremCard` lifecycle + independent reproduction pack |
| PO-18 | The electromechanical and gluing theorem cards prove representative/gauge/cut invariance and whole-interface work/GCL identities, while the force card encloses every named numerical/geometry/material contribution | commuting-diagram checker + guaranteed functional bounds |
| PO-19 | The whole-machine thermodynamics card composes open-system first/second-law statements across multirate windows with explicit boundary entropy/exergy; unresolved defect cannot be hidden in a component audit | compositional theorem kernel + runtime boundary ledger |
| PO-20 | Fidelity descent and conjugate-geometry cards preserve evidence monotonicity and bound model/envelope naturality defects; a counterexample revises the theorem rather than being discarded | crosswalk/envelope proof kernels + counterexample corpus |
| PO-21 | Multicontact action/reaction, stored contact energy, dissipation, passive impact work, frictional heat and wear-state updates close globally; compliant/barrier/nonsmooth substitutions expose model discrepancy and limit path | contact/tribology work and state-transition ledger |
| PO-22 | Acoustic lanes satisfy radiation, complex-power/energy-flux, dispersion and interface reciprocity where applicable; source transfer/surface equivalence and turbulence/combustion source validity are separate claims | analytic radiation/dispersion + structural/source-transfer batteries |
| PO-23 | Sampled/periodic/hybrid controls declare transform conventions and pass stability/region-of-attraction, observability/detectability, implementation-error, timing/fault-isolation and closed-loop risk obligations appropriate to their nonlinearities | lifted/Floquet/IQC/Lyapunov/reachability + fault batteries |
| PO-24 | Life/reliability claims bind load/environment/process/failure mode, dependence model, sampling measure and importance weights; competing risks never assume independence without evidence and validation uses held-out coupon→machine data | limit-state/rare-event reproducibility + blind validation |
| PO-25 | Stable IDs/lineage, transactional snapshots and safety cases reject ambiguous rebinding or half-committed state; regulatory/standard conformance remains orthogonal to scientific evidence | lineage mutation, cancel/resume, hazard/fault and conformance batteries |

### 12.2 Assumptions ledger (top-level; each crate CONTRACT refines)

Every artifact refines a row using the schema
`AssumptionId | predicate | scope | evidence | runtime monitor | violation
effect | owner | expiry/review gate`:

| ID | Predicate/scope | Evidence and monitor | Violation effect / owner / review |
|---|---|---|---|
| A-001 | rigid/reduced body rung is adequate for named QoIs | mode/truncation and interface-compliance indicators | escalate flexible rung or refuse / fs-mbd / E3 |
| A-002 | MQS wavelength/displacement-current regime is valid | dimensionless regime receipt and frequency monitor | electroquasistatic/full-wave escalation / fs-em / E4–E6 |
| A-003 | 0-D cylinder is spatially mixed; 1-D ducts are section-averaged | stratification, knock, wave/secondary-flow indicators | multi-zone/multi-D escalation / fs-gas / E5 |
| A-004 | smooth/contact capability and continuum scale are adequate | curvature/thickness/roughness/Knudsen/contact-duration receipts | change representation/law or refuse / fs-contact + Rep Router / E2–E6 |
| A-005 | symmetry group preserves geometry, state, BCs and excitation | representation residual at every update | expand sectors/full model / domain owner / each solve |
| A-006 | material/process/query lies inside validity predicate | `PropertyUsageReceipt` monitor | demote, recalibrate, extrapolation refusal / fs-matdb consumer / each query |
| A-007 | closure/correlation/turbulence/lubrication regime applies | named dimensionless groups and held-out discrepancy | resolved/high-fidelity escalation / law owner / phase gate |
| A-008 | probability/dependence/population model represents manufacture/use | data lineage, posterior predictive and drift monitor | decision-risk demotion or robust worst-case lane / fs-uq / E7 |

Leaving any domain demotes or refuses according to policy—it never silently
extrapolates.

**Program risk register.** This is distinct from runtime physics fallbacks and
uses `RiskId | owner | likelihood | impact | leading indicator | trigger |
mitigation | contingency | residual risk | review gate` in the ledger. Seed
risks are: flowpipe wrapping/certificate cost; contact candidate/certificate
explosion; coefficient/topology-robust EM preconditioning; chemistry stiffness;
material-data licensing and experimental scarcity; hybrid-adjoint invalidity;
same-layer Cargo cycles; schema/hash migration; theorem-tool trusted computing
base; single-node scale/accelerator-policy mismatch; interoperability adoption;
and misuse of scientific receipts as safety/regulatory certificates. Each has a
named owner, quantitative trigger and kill/escalation action before E0 planning
can close.

### 12.3 Failure modes and fallback policy (selected)

| Failure mode | Detection | Fallback |
|---|---|---|
| DAE/nonholonomic constraint drift | position/velocity/acceleration/Pfaffian residuals and optional state-distance bound | symmetry/work-compatible projection only with full `ProjectionReceipt`; otherwise reject/retry or debit impulse/energy/momentum/adjoint defect |
| Contact chattering / Zeno | event-rate budget | Moreau–Jean lane switch + restitution model review; typed refusal on budget exhaustion |
| Newton divergence on saturating iron | stall diagnosis (existing pattern) | Picard lane; load stepping |
| Curl-curl null-space pollution | residual stagnation + gauge check | tree-cotree re-gauge; harmonic deflation (`deflate_harmonics`) |
| Sliding-interface non-conservation | direct power/flux/GCL audit, with sheaf compatibility localization | refine basis/interface; mortar or monolithic lane; refuse exchange |
| Stiff combustion blow-up | invariant-domain/positivity trip | IMEX substepping; calibrated Wiebe fallback with its card-specific `Validated` or `Estimated` status |
| Multirate interface instability | residual + energy/contractivity monitor | shrink window; monolithic strong-coupled fallback |
| Grazing/simultaneous events | true-flow root-count or reset-order ambiguity | set-valued `PossibleEvent`/`Unknown`; inclusion/time-stepping lane |
| Material/process mismatch | typed domain and usage-receipt failure | calibrated substitute with model discrepancy or refuse claim |
| Contact indeterminacy/Painlevé | cone primal/dual/determinacy receipt | set-valued response or explicitly different calibrated compliant model with discrepancy/limit study, never arbitrary impulse |
| Reduced model leaves domain | discrepancy/escalation monitor | 2-D→2.5-D→3-D or network→field escalation |

---

## 13. Galaxy-Brain Transparency Cards

**Card 1 — Tangent mobility, dual self-stress, and finite motion.**
*Equation:* at configuration `q`, infinitesimal motions are `ker J(q)` after
gauge removal and self-stresses are `ker(J(q)ᵀ)`; a proven sheaf realization may
identify the former with `H⁰(K_q)`.
*Substitution:* Bennett linkage — 4 bodies, 4 revolutes with special skew
geometry: Grübler gives `6(4−1) − 4·5 = −2` (predicts an over-constrained,
immobile structure); the free-floating tangent kernel built from the *actual*
joint axes is 7-dimensional (6 global rigid motions + 1 internal motion)
→ one infinitesimal internal mode. Validated continuation/second-order analysis
must still prove a finite one-dimensional branch. *Plain English:* counting
assumes generic independent constraints; the actual Jacobian measures this
geometry, and continuation asks whether the infinitesimal motion really moves.
*Would change:* an interval-rank or continuation `Unknown` remains visible.

**Card 2 — Topology is fixed; constitutive evidence is not.**
*Equation:* heat uses `M_{ρc}Tdot + d₀ᵀM_kd₀T=f+b`; magnetics uses a gauged
`d₁ᵀM_νd₁A=j` formulation; steady current uses
`d₀ᵀM_σd₀φ=i`. *Substitution:* a named copper process-state card returns
`σ(T,history)` with units, covariance and validity receipt; that sample weights
the dissipative conductivity form `M_σ`. *Plain English:* connectivity stays exact,
but stability, boundary conditions, quadrature and material credibility still
need proof. *Would change:* anisotropy, saturation, hysteresis or extrapolation
selects a richer constitutive node and may demote the physical claim.

**Card 3 — Certified event capture.**
*Equation:* a `TrajectoryTube Q(I)` encloses the true flow; exclusion proves
`0∉g(Q(I))`, while existence/uniqueness and root-count certificates establish
and enumerate events. *Substitution:* the Geneva benchmark records the actual
subdivision count and time enclosure—no unmeasured nanosecond claim appears in
the plan. *Plain English:* prove over the true motion, not sampled dense output.
*Would change:* grazing/simultaneous contact returns a set or `Unknown`.

**Card 4 — Flux linkage as a relative-topology/material pairing.**
*Equation:* a certified filament loop uses
`λ = N⟨A,c⟩ = N ∫_S B·n dS`; a production stranded coil uses
`λ = ∫ A·J_unit dV`; `v = Ri + dλ/dt` under the declared sign.
*Substitution:* the same slot winding is represented both ways on a fixture and
must agree within discretization/end-effect bounds. *Assumptions:* closed
cycle/spanning surface, gauge compatibility and winding orientation.
*Would change:* solid conductors, motion, skin/proximity and end windings select
MQS/distributed-current lanes.

**Card 5 — The dwell-retention proof.**
*Equation:* for an **ideal zero-clearance declared active boundary mode**, the
first-order cone `{ξ:Aξ=0,Cξ≥0}` may reduce to gauge motion, with a Farkas/
conic-duality witness and `λ_{m+1}(N_V)≥γ²` (equivalently
`σ_min(J̃|K_y^⊥)≥γ` in the whitened coordinates defined in §3.2). Degenerate
modes additionally require second-order or global
configuration-space proof. For manufactured positive clearance, the artifact
is instead `DwellRetentionCertificate { play_set, nonescape,
boundary_mode_margins }`.
*Substitution:* the 4-station Geneva globally bounds the dwell reachable set
over the full tolerance, active/inactive-mode, load, impulse, friction,
elasticity and wear box; margins and maximum play/creep are outputs, not
invented constants or corner samples.
*Plain English:* the proof says how far the wheel may move and why it cannot
escape the dwell region under declared disturbances—not falsely that a wheel
with clearance has no local motion. *Would change:* failed invariance, friction
feasibility or elastic bounds demote/refuse retention and localize the escape
mode.

---

## 14. Claim Boundaries and Ambition Horizon

Boundaries are rung-specific, not ceilings. The first default path makes no
full-wave EM, vector-hysteresis, turbulent reacting-flow, emissions,
multi-D chamber, EHL, deformable/deformable contact, HIL or real-time claim.
Those capabilities remain explicit `[F]/[M]` rungs in §§5, 8 and 10 and may
advance only with their Gauntlets. Reduced 0-D/1-D/network/dq/2-D lanes are
labeled on every downstream artifact and carry escalation criteria.

Wear/fatigue/fracture/corrosion life begins `Estimated`; it can become
`Validated` for a named material/process/load spectrum and experimental
domain, but numerical verification never turns model-form life prediction into
universal truth. Smooth-contact derivatives are local-regime claims;
nonsmooth/generalized sensitivities state their exact regularity. Any
unresolved event, material extrapolation, symmetry violation, nonconservative
interface, gauge failure, solver indeterminacy or unclosed balance prevents the
strong claim while preserving the partial artifacts.

---

## 15. Operationalization

1. **Beads conversion.** This plan should be decomposed into `br` issues per
   phase with dependency edges mirroring §10 (E0 items unblock E1–E5;
   flagships depend on their §8 composition lists), including explicit unit +
   e2e test beads per feature cluster, per the house beads workflow. That
   conversion is a follow-up work session with 4–5 plan-space refinement
   passes; this document is written to be self-contained input for it.
2. **Contracts first.** Each new crate lands `CONTRACT.md` (all ten required
   sections, honest no-claims from §14) *before* it becomes a dependency
   target — `xtask check-contracts` enforces.
3. **First program increment (E0 + thin E1 slices):** six-base unit/schema
   migration; typed material/interface query receipts; weighted storage/
   transport heat G1; neutral port relation protocol; minimal versioned machine
   IR/admission plus ledger/package migrations; nonlinear/block solver;
   exact prescribed `SpacetimeChart` with a derived trochoid oracle; classical
   event baseline followed by validated exact-path root accounting. Each slice
   is independently mergeable and useful; no calendar estimate precedes task
   decomposition and a measured baseline.
4. **Dev-only oracles.** External references (FEniCS/FEMM for EM fixtures,
   Chrono/Siconos for contact fixtures, Cantera for thermochemistry tables)
   are used strictly as isolated dev-time comparison oracles per the
   dependency policy — never in the production graph.
5. **DSR lanes.** New crates join `dsr quality --tool frankensim`; flagship
   smoke stages get provenance-pinned metric goldens with explicit tolerances.
   Hash equality is reserved for canonical artifacts whose determinism contract
   actually requires it; roofline receipts precede every performance claim.
6. **Layer/cycle admission.** Before a new crate is generated, add its proposed
   manifest edges to an architecture fixture, run the `xtask` dependency rules,
   and require `cargo metadata` to resolve the real same-layer DAG without a
   package cycle. Contact never depends on its solid/mbd consumers;
   circuits never depend on EM for shared types; L1 spectral services never
   consume L2/L3 backends; L3 solvers never depend upward on ASCENT.
7. **Research-bet ledger.** Each `[F]/[M]` Bead names one baseline, activation
   metric, kill criterion, cost distribution, fallback and retained falsifier.
   Parallel portfolio work is encouraged; a single proof lane never entangles
   multiple unvalidated inventions.

---

## Appendix A — Extended Port Table

| PortKind | Effort | Flow | Effort·Flow | First users |
|---|---|---|---|---|
| Mechanical (translation) | force [N] | velocity [m/s] | power | existing |
| **Rotational** | torque [N·m] | angular velocity [rad/s] | power | fs-mbd ↔ fs-em, fs-gas |
| Fluid | pressure [Pa] | volume flow [m³/s] | power | existing |
| Thermal | temperature [K] | entropy flow [W/K] | power | heat boundary relation; audit includes boundary temperature |
| **Electrical** | voltage [V] | current [A] | power | fs-circuit ↔ fs-em |
| **Magnetic** | mmf [A] | flux rate [Wb/s] | power | fs-em internal decompositions |
| **Chemical** | species/reaction electrochemical-potential vector `μ⃗` [J/mol] | amount-flow/reaction-progress vector [mol/s] | `μ⃗·ṅ⃗` or affinity·progress power | basis, element/charge constraints, electric-potential convention and mass/molar crosswalk are explicit |
| **Stream bundle** | pressure/velocity/temperature/composition state | mass, species, momentum, total-enthalpy and entropy fluxes | multi-balance | gas/liquid network interfaces |

The effort/flow table defines conjugate views; it does not mean a flowing
stream can be reduced to one scalar port. Every relation records orientation
and sign, and admission proves dimensional power before execution.

## Appendix B — Material Card Schema (sketch)

```rust
/// Immutable measured evidence; parsing/normalization happens offline.
pub struct ObservationDataset {
    pub id: ObservationDatasetId,
    pub specimen_and_process: SpecimenRecord,
    pub method_and_instrument: MeasurementMethod,
    pub observations: ObservationArtifactRef,
    pub covariance_and_censoring: ObservationUncertainty,
    pub provenance: SourceRecord,
}

/// One non-destructively retained interpretation of observations.
pub struct PropertyClaim<K> {
    pub id: PropertyClaimId,
    pub key: K,
    pub value: PropertyValue,        // scalar/curve/tensor/distribution/model parameters
    pub validity: fs_evidence::ValidityDomain,
    pub uncertainty: UncertaintyModel,
    pub joint_uncertainty: Vec<CorrelationDatasetId>,
    pub interpolation: InterpolationPolicy,
    pub evidence: EvidenceDescriptor,
    pub observations: Vec<ObservationDatasetId>,
    pub provenance: SourceRecord,    // citation, license, specimen/process, ContentHash
}

/// Immutable L1 declaration of a law identity and state metadata.
/// LawId/LawVersion/parameter/schema/policy types live in fs-matdb; the
/// executable implementation and runtime state behavior live in L3.
pub struct ConstitutiveModelCard {
    pub law_id: LawId,
    pub law_version: LawVersion,
    pub parameters: CanonicalParameterBlock,
    pub state_schema: StateSchemaVersion,
    pub initial_state_policy: InitialStatePolicy,
    pub validity: fs_evidence::ValidityDomain,
    pub source_hashes: Vec<fs_blake3::ContentHash>,
}

pub struct MaterialCard {
    pub id: MaterialStateId,          // chemistry + phase + process/temper + revision
    pub schema_version: MaterialSchemaVersion,
    pub revision: RevisionId,
    pub supersedes: Vec<MaterialStateId>,
    pub claims: BTreeMap<PropertyClaimId, PropertyClaim<PropertyKey>>,
    pub constitutive_models: BTreeMap<LawId, ConstitutiveModelCard>,
    pub by_key: BTreeMap<PropertyKey, Vec<PropertyClaimId>>,
    pub content_hash: fs_blake3::ContentHash,
}

pub struct InterfaceSystemCard {
    pub surface_a: SurfaceState,
    pub surface_b: SurfaceState,
    pub intervening_medium: MediumState,
    pub lubricant_or_third_body: Option<MediumState>,
    pub environment: EnvironmentState,
    pub texture_frame: BasisFrameRef, // opaque L1 ID; transforms are owned above L1
    pub claims: BTreeMap<PropertyClaimId, PropertyClaim<InterfaceKey>>,
    pub constitutive_models: BTreeMap<LawId, ConstitutiveModelCard>,
    pub by_key: BTreeMap<InterfaceKey, Vec<PropertyClaimId>>,
    pub content_hash: fs_blake3::ContentHash,
}

/// Owned and evolved by the consuming L3 domain law, never by fs-matdb;
/// its neutral IDs/schema versions are imported from fs-matdb.
pub struct ConstitutiveRuntimeState {
    pub law_id: LawId,
    pub law_version: LawVersion,
    pub state_schema: StateSchemaVersion,
    pub canonical_state_bytes: Vec<u8>,
}

/// Persisted by L6 fs-ledger/fs-package, not declared as an L1 dependency.
pub struct StateCheckpointReceipt {
    pub state_slot: StateSlotId,
    pub law_id: LawId,
    pub law_version: LawVersion,
    pub state_schema: StateSchemaVersion,
    pub runtime_state: fs_blake3::ContentHash,
    pub canonical_parameters: fs_blake3::ContentHash,
    pub contract_and_code_hash: fs_blake3::ContentHash,
}

pub struct PropertyUsageReceipt {
    pub query: PropertyQuery,
    pub selected_claims: Vec<PropertyClaimId>,
    pub selection_or_fusion: ClaimSelectionPolicy,
    pub correlation_datasets: Vec<CorrelationDatasetId>,
    pub sample: Evidence<PropertySample>,
    pub evaluator_and_law_version: EvaluatorVersion,
    pub frame_transform: Option<FrameTransformReceipt>,
    pub interpolation_or_extrapolation: EvaluationDecision,
    pub source_artifacts: Vec<fs_blake3::ContentHash>,
}
```

## Appendix C — Glossary (selected)

- **Conjugate action:** the kinematic condition under which two contacting
  profiles transmit a constant (or specified) velocity ratio. Under the
  declared regularity, active-branch, visibility, endpoint/trimming, contact,
  and noninterference hypotheses, this is equivalent to the relevant envelope/
  common-normal condition (the fundamental law of gearing —
  classical Euler–Savary governs the associated curvature analysis for planar
  relative motion/centrodes; spatial families use their own conjugate-surface
  differential geometry and contact-curvature machinery).
  Properly generated, noninterfering ideal involutes satisfy it under their
  declared center-distance/contact assumptions; the certificate in §3.3
  measures deviation outside that oracle.
- **Trochoidal Wankel geometry:** the ideal apex-point path, finite apex-seal
  center/contact locus, and actual housing-bore envelope under a declared
  seal-tip/clearance model are distinct artifacts. The rotor is Reuleaux-like,
  not assumed to be an exact constant-width Reuleaux triangle.
- **IPC family:** incremental-potential barrier contact. Nonintersection is
  conditional on admissible initial state, conservative candidate/CCD and an
  accepted solve; differentiability is local to regular feature/active regimes.
- **Moreau–Jean:** measure-differential-inclusion time-stepping for nonsmooth
  dynamics, with configuration admissibility/active-set support plus a discrete
  velocity–impulse/restitution law; it is not velocity-only penetration repair.
- **Hiptmair–Xu:** auxiliary-space preconditioner for H(curl)/H(div) systems
  combining edge smoothing, scalar-gradient and vector `[H¹(Ω)]³` auxiliary
  solves/interpolation, and topology-aware harmonic coarse treatment under
  declared boundary/coefficient assumptions.
- **Dirac structure / port-Hamiltonian:** a lossless power-conserving
  interconnection structure. Passivity additionally requires admissible storage
  and dissipative relations plus a compatible time discretization.
- **Tangent kinematic complex (or proved sheaf realization):** a configuration-
  dependent linear object whose kernel can encode infinitesimal compatibility
  when its maps are explicitly defined; finite motion is nonlinear holonomy/
  continuation and self-stress belongs to the dual kernel.
- **Equivariant sliding-interface sheaf:** organizes time-dependent trace compatibility
  and defect localization. Conservative transfer is separately proved by
  primal/dual, flux, inf-sup and GCL identities.
- **StandardConformance:** evidence that a named edition's formulas, procedures
  or test protocol were implemented and reproduced inside their published
  scope. It is orthogonal to numerical `Verified` and empirical `Validated`.

## Appendix D — Primary Reference and Benchmark Registry Seed

These references anchor definitions and benchmark provenance; they are not
authority-by-citation. Each implementation Bead records `unread`, `read`,
`derived`, `reproduced`, or `independently_falsified` status and pins the exact
artifact/version used.

| Topic | Seed source | Plan consequence |
|---|---|---|
| FEEC stability | Arnold, Falk & Winther, [Finite element exterior calculus: from Hodge theory to numerical stability](https://arxiv.org/abs/0906.4325) | `dd=0` is insufficient; subcomplex + bounded cochain projection and formulation stability are explicit gates |
| Cellular sheaf spectra | Hansen & Ghrist, [Toward a Spectral Theory of Cellular Sheaves](https://arxiv.org/abs/1808.01513) | sheaf Laplacians/cohomology become relevant only after explicit stalks/restrictions are constructed; L1 spectral machinery stays generic |
| Sheaf–cosheaf foundations | Curry, [Sheaves, Cosheaves and Applications](https://arxiv.org/abs/1303.3255) | restriction and aggregation directions are distinct; §3.9's power pairing and commuting diagram are new theorem obligations, not borrowed authority |
| Species thermodynamics | McBride, Zehe & Gordon, [NASA/TP-2002-211556](https://ntrs.nasa.gov/citations/20020085330) | NASA-9 uses standard-state molar thermodynamics; amount-of-substance, phase/EOS derivations, and reference metadata are load-bearing |
| Entropy-stable discretization | Tadmor, [Entropy Stable Schemes](https://doi.org/10.1016/bs.hna.2016.09.006) | entropy-conservative spatial flux is not the whole fully discrete/source/boundary theorem; the exact scheme and sign convention are pinned |
| Port-Hamiltonian dissipation | Cervera, van der Schaft & Baños, [Interconnection of port-Hamiltonian systems and composition of Dirac structures](https://doi.org/10.1016/j.automatica.2006.08.014) | Dirac relations are power-conserving; dissipation is a separate resistive structure |
| Rigidity/self-stress | Rocks et al., [Integrating local energetics into Maxwell–Calladine constraint counting](https://arxiv.org/abs/2208.07419) | mechanisms use `ker J`; self-stress uses `ker(Jᵀ)`; index and energy structure are distinct |
| Contact | Li et al., [Incremental Potential Contact](https://doi.org/10.1145/3386569.3392425) and [Convergent IPC](https://arxiv.org/abs/2307.15908) | nonintersection/convergence claims list candidate, CCD, solve and refinement assumptions; the Bead pins the exact evolving preprint version/hash; smoothness is regime-specific |
| Codimensional contact | Li et al., [Codimensional Incremental Potential Contact](https://arxiv.org/abs/2012.04457) | shell/rod/seal thickness and nonintersection claims remain conditional on the declared codimensional model, candidate generation, CCD, accepted solve, and thickness assumptions |
| Nonsmooth multicontact | Acary, [Energy conservation and dissipation properties of time-integration methods for nonsmooth elastodynamics with contact](https://arxiv.org/abs/1410.2499), plus the exact global impact-law source pinned by the Bead | Moreau–Jean energy behavior and multicontact admissibility are law/parameter dependent; pairwise restitution is not a global energy theorem |
| Validated flowpipes/events | Walawska & Wilczak, [An implicit algorithm for validated enclosures of the solutions to variational equations for ODEs](https://arxiv.org/abs/1509.07388), plus a pinned validated-DAE/root-isolation source selected by the Bead | true-flow enclosure, DAE regularity and finite-guard isolability are separate obligations; arbitrary smooth guards do not get complete finite root claims |
| Gear transmission error | Athavale, Krishnaswami & Kuo, [SAE 2001-01-1507](https://saemobilus.sae.org/papers/estimation-statistical-distribution-composite-manufactured-transmission-error-a-precursor-gear-whine-a-helical-planetary-gear-system-2001-01-1507) | TE is angular-position deviation from ideal conjugate motion, not interval width; an executable G2 case still needs a complete reproducible deck |
| Wankel seals/kinematics | Handschuh & Owen, [NASA/TM-2010-216353](https://ntrs.nasa.gov/citations/20100036253) | reduced seal-load/friction reference with explicit assumptions; it does not prove universal finite-radius housing geometry |
| Multibody benchmarks | IFToMM, [Library of Computational Benchmark Problems](https://www.iftomm-multibody.org/benchmark/) | exact input/result artifacts for Andrews, Bricard, slider-crank and other cases |
| TEAM electromagnetic benchmarks | COMPUMAG official definitions: [10](https://www.compumag.org/jsite/images/stories/TEAM/problem10.pdf), [13](https://www.compumag.org/jsite/images/stories/TEAM/problem13.pdf), [20](https://www.compumag.org/jsite/images/stories/TEAM/problem20.pdf), [24](https://www.compumag.org/jsite/images/stories/TEAM/problem24.pdf), [30a](https://www.compumag.org/jsite/images/stories/TEAM/problem30a.pdf), [30b](https://www.compumag.org/jsite/images/stories/TEAM/problem30b.pdf) | pin the exact revision, geometry, excitation, material law, circuit, QoI, and acceptance data; the number alone is not an input deck |
| Thermal benchmarks | NAFEMS, [Thermal analysis benchmark index/guide](https://www.nafems.org/publications/glossaryofbenchmarks/thermalanalysis/) | exact case ID/report/license/QoI required; “NAFEMS set” alone is not executable |
| NASA-9 dev oracle | Cantera 3.2, [NASA9 parameterization documentation](https://cantera.org/3.2/reference/thermo/species-thermo.html) | oracle receipt pins Cantera release, species/mechanism hash, temperature region, reference pressure, and units; never a production dependency |
| Nonholonomic mechanics | Modin & Verdier, [What makes nonholonomic integrators work?](https://doi.org/10.1007/s00211-020-01126-y), plus the exact discrete Lagrange–d'Alembert source selected by the implementation Bead | rolling/no-slip is Pfaffian and generally nonintegrable; it does not inherit ordinary Hamilton/RATTLE claims |
| H(curl) auxiliary space | Hiptmair & Xu, [Nodal auxiliary space preconditioning in H(curl) and H(div) spaces](https://doi.org/10.1137/060660588), and official [hypre AMS documentation](https://hypre.readthedocs.io/en/latest/solvers-ams.html) | edge smoothing, scalar/vector auxiliary solves, coordinate/interpolation and topology/boundary assumptions are explicit |
| Switched descriptor impulses | Yildiz, [A MNA-Based Unified Ideal Switch Model for Analysis of Switching Circuits](https://doi.org/10.1142/S0218126613500461) | inconsistent ideal switching can cause impulses and discontinuous device states; impulse-free continuity is conditional |
| Nonlinear magnetic energy | Mandlmayr & Egger, [On implicit interpolation models for nonlinear anisotropic magnetic material behavior](https://arxiv.org/abs/2311.02380) | convex energy/coenergy and differential-tangent claims are distinguished from merely strongly monotone nonintegrable laws |
| Reacting entropy/positivity | Ching, Johnson & Kercher, [Part I](https://arxiv.org/abs/2211.16254) and [Part II](https://arxiv.org/abs/2211.16297) | entropy/positivity constructions are tied to the exact thermodynamics, chemistry, source and discretization; no universal flux formula is assumed |
| Acoustics radiation and sources | Burton & Miller, [combined-field exterior acoustics](https://doi.org/10.1016/0022-247X%2884%2990146-X); Ffowcs Williams & Hawkings, [moving-surface aeroacoustics](https://doi.org/10.1098/rsta.1969.0031) | exterior resonance treatment and hybrid source transfer have distinct assumptions and receipts |
| Acoustic benchmarks | NASA, [CAA benchmark proceedings](https://ntrs.nasa.gov/citations/20050217396); NAFEMS, [R0083 acoustic benchmark](https://www.nafems.org/publications/resource_center/r0083/) | exact deck, medium/source/QoI, radiation convention and acceptance uncertainty are pinned |
| Tribology/EHL | Hamrock & Dowson, [NASA-TP-1342](https://ntrs.nasa.gov/citations/19780025504), followed by independent optical-film and traction datasets selected by the Bead | film-thickness correlations are regime-bounded baselines; full EHL additionally verifies pressure/film/temperature/mass and force balance and validates against held-out measurements |
| Fatigue and fracture data | ASTM [E466-21](https://store.astm.org/standards/e466) and [E647-24](https://store.astm.org/standards/e647), plus a pinned [NASA NASGRO reference artifact](https://ntrs.nasa.gov/citations/20250011200) as a dev oracle | coupon method, environment, load ratio, crack scale/closure and transferability are explicit; test-method conformance or oracle agreement never validates a component/service-life prediction |
| Machinery and motor assurance | [ISO 12100:2010](https://www.iso.org/standard/51528.html), [IEC 60034-1:2026](https://webstore.iec.ch/en/publication/89961), and [IEC 60034-2-1:2024](https://webstore.iec.ch/en/publication/67756) | exact scope, rating/loss test method and hazard process become conformance artifacts; simulation evidence remains distinct from regulatory approval or a laboratory test |
| V&V and uncertainty | ASME, [VVUQ standards registry](https://www.asme.org/codes-standards/publications-information/verification-validation-uncertainty); BIPM/JCGM, [GUM publications](https://www.bipm.org/en/committees/jc/jcgm/publications); NASA, [NASA-STD-7009](https://standards.nasa.gov/standard/NASA/NASA-STD-7009) | context of use, solution verification, experimental uncertainty, calibration/validation split and prediction assessment are separate artifacts |
| Combustion validation | Sandia, [Engine Combustion Network](https://ecn.sandia.gov/) | exact ECN configuration, diagnostic uncertainty, fuel/injector/ambient state and blind/held-out QoIs are versioned; a family name is not a validation result |
| Interoperability | Modelica Association, [FMI 3.0.2](https://fmi-standard.org/docs/3.0.2/) and [SSP 2.0](https://ssp-standard.org/docs/2.0/) | conformance and quarantined adapter behavior are tested without admitting foreign runtimes to the trusted solver graph |
| Gear standards | ISO, [ISO 6336-1:2019 scope](https://www.iso.org/standard/63819.html), followed by exact applicable parts/editions/worked examples | formula implementation earns verification/standard conformance only inside scope; physical drive life needs independent validation |

Before accepting TEAM/CFR/AGMA/ISO or any experimental family named elsewhere,
the corresponding registry row must resolve the exact problem/report edition,
redistribution terms, geometry/material/excitation, QoIs, uncertainty and
acceptance band. A standards calculation earns numerical `Verified` and/or
`StandardConformance` for the exact edition and scope. Its physical prediction
is `Estimated` unless independent held-out evidence for the named QoI and
population earns `Validated`; no color is inherited from the publisher's name.

## Appendix E — Requirement-to-Evidence Traceability Seed

This table is a seed for a generated ledger, not a hand-maintained status
dashboard. Each Bead expands its row to
`RequirementId | capability/property | blocker | owner/artifact | prerequisite
phase | milestone | flagship | benchmark/data | proof obligation | claim
boundary | status`; catalog generation fails if a requirement has no owner,
gate, evidence route or honest boundary.

| ID | Capability/property and blocker | Owner/artifact | Phase / forcing evidence | Obligation and boundary | Status |
|---|---|---|---|---|---|
| B1 | moving geometry, sweeps, ALE absent | `fs-motion::CertifiedMotorTube`, Rep Router, fs-flux | E1/E5; trochoid, moving-wall/GCL decks | PO-5/14; rigid first, deformable successor | proposed |
| B2 | MBD/joints/DAE/nonholonomic dynamics absent | fs-kinematics, fs-mbd, fs-time | E1–E2; IFToMM, rolling disk/sleigh | PO-3/4/10; regular/inclusion classes explicit | proposed |
| B3 | body contact/CCD/penetration absent | capability-routed fs-query/fs-contact/fs-tribo | E2; Hertz, Painlevé, multicontact | PO-21; nonintersection conditional | proposed |
| B4 | gears/cams/Geneva/linkages absent | fs-kinematics/fs-machine | E1–E3; Geneva/gear flagships | PO-5/10/20; family-specific standard/contact scope | proposed |
| B5 | no EM formulation | fs-rep-mesh/fs-feec/fs-em | E0–E4; TEAM + motor/generator | PO-1/8/9/13; MQS first, full-wave separate | proposed |
| B6 | no thermal-domain stack | fs-thermal/fs-xform/fs-couple | E0–E5; Stefan/NAFEMS/calorimetry | PO-1/2; data/model/radiation rungs separate | proposed |
| B7 | no compressible/species/combustion | fs-thermochem/fs-gas/fs-flux | E5–E6; Riemann/nozzle/ECN/engine | PO-6/11; closure and experimental validity named | proposed |
| B8 | no cross-physics data/units/interface history discipline | fs-qty/fs-evidence/fs-matdb/fs-material | E0; pack/query/cross-version batteries | PO-7/11/15; immutable evidence vs runtime state | proposed |
| B9 | scalar coupling seed only | fs-couple `PortSchema`, domain transfer adapters | E0–E7; coupled vertical slices | PO-2/8/9/19; no passivity from topology alone | proposed |
| B10 | solver/time/adjoint gaps | fs-solver/fs-time/fs-ad/fs-adjoint | E0–E6; MMS/event/cancel/gradient gates | PO-3/4/6/8; nonsmooth no-claim where needed | proposed |
| B11 | mechanism spectral/continuation gaps | fs-spectral/fs-kinematics/fs-time | E1–E3; Bennett/Bricard/Campbell | PO-10; scaling/gauge/multiplicity serialized | proposed |
| B12 | no oriented 2-D complex/FEEC path for planar fields | fs-rep-mesh `TriComplex2` + fs-feec 2-D families | E0–E4; exactness/orientation/MMS and planar EM | PO-1/13; a surface half-edge mesh is not sufficient | proposed |
| B13 | machine flexibility/life/NVH/control missing | fs-mbd/fs-solid/fs-acoustics/fs-control/UQ upgrades | E3–E7; gear/motor/Wankel/genset | PO-22–24; held-out population/QoI required | proposed |
| B14 | assurance/V&V/workflow/scale absent | Machine IR, V&V artifacts, safety case, Scale/Competitive ledgers | E0–E8; blind tests, faults, FMI/SSP, scale suites | PO-16/17/25; scientific receipt is not certification | proposed |
| RQ-ROLL | truly rolling constant-width/Reuleaux mechanism | nonholonomic fs-mbd + fs-contact; `fs-constant-width-e2e` | E2/E7; support/centroid, branch-switch, tolerance decks | PO-3/4/21; distinct from sliding Wankel seals | proposed |
| RQ-GEAR | broad real gear mechanisms | conjugate/spatial fs-kinematics + fs-machine | E1–E3; law-of-gearing, loaded TE, ISO scope | PO-5/10/20; each family has its own contact/standard class | proposed |
| RQ-FRICTION | real friction/contact/wear | `InterfaceSystemCard`, fs-tribo runtime state | E2–E3; incline/Hertz/EHL/wear holdouts | PO-15/21; system/history-specific, not pair constant | proposed |
| RQ-CONSTITUTIVE | weighted, nonlinear, coupled and history-dependent matter laws | fs-material `ConstitutiveGraph` + domain nodes | E0–E8; coupled-law MMS and held-out material tests | PO-1/7; reversible skew, dissipation and empirical validity remain distinct | proposed |
| RQ-DENSITY | density→mass/COM/inertia | fs-matdb receipt + `fs-query::GeometricMoments` | E1; analytic solids/cross-representation | PO-7/15; nonuniform density escalates integration | proposed |
| RQ-MECHMAT | ductility, hardness, toughness, fatigue/creep/corrosion | fs-material/fs-solid life ladder | E3–E7; coupon→component→machine holdouts | PO-15/24; no universal service-life claim | proposed |
| RQ-ELEC | electrical/thermal conductivity, dielectric/insulation | weighted operators, fs-em/fs-thermal/fs-power | E0/E4; resistance/heat/insulation fixtures | PO-1/2/13; named frequency/T/process regime | proposed |
| RQ-MAG | permeability/BH/remanence/magnetic response/moment | energy/state cards + fs-em; moment integrated over body | E4; TEAM and motor dynamometer | PO-1/13; no geometry-free total moment or scalar-spline theorem | proposed |
| RQ-PHASE | latent heat/phase density | fs-thermal total enthalpy + state/remap | E0/E5; Stefan/laser-flash | PO-1/2; apparent-heat regularization labeled | proposed |
| RQ-FLUID | viscosity, bulk modulus, vapor pressure, hydraulics/EHL | fs-matdb/fs-flux/fs-gas/fs-tribo | E3/E5; Couette/bearing/water-hammer | PO-2/14/15/21; rheology/phase regime explicit | proposed |
| RQ-PERMEATE | intrinsic permeability, gas transport and membranes | §5.13 transport adapters | E5–E6; Darcy/layer/membrane decks | PO-11/15; permeability/conductivity/permeance definitions distinct | proposed |
| RQ-WET | hydrophilic/phobic/capillary behavior | dynamic `InterfaceSystemCard` + free-surface lane | E5–E6; Young–Laplace/advancing-receding tests | PO-2/15; rate/roughness/contamination bound | proposed |
| RQ-MOTORGEN | motors and generators with electronics/control/thermal | `fs-motor-e2e` | E4/E7; TEAM, converter, dynamometer, EMC | PO-2/8/9/13/23/25; fidelity and held variables explicit | proposed |
| RQ-ICE | piston/Wankel ICE and genset | `fs-wankel-e2e`, `fs-genset-e2e` | E5–E7; pressure/ECN/calorimetry/emissions | PO-2/4/6/11/14/24; correlation colors stay QoI-specific | proposed |
| RQ-ACOUSTIC | NVH/aero/combustion acoustics | fs-acoustics | E3–E7; NAFEMS R0083/NASA CAA/measurements | PO-22; source validity separate from propagation | proposed |
| RQ-ACTIVE | piezo/magnetostrictive/thermoelectric/electrochemical | §5.13 coupled-law adapters | E4–E8; reciprocity/energy/entropy decks | PO-1/2; guarded `[F/M]` until validated | proposed |

---

*End of plan. The constitution remains `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`
for mission, layers, dependency policy and the Decalogue. The explicit
Normative Delta and Ratification Register at the beginning of this charter
governs its named refinements until they are mirrored upstream; live crate
contracts govern implemented reality. Every other conflict is resolved loudly,
with a reviewable diff and ratification Bead—not by silent precedence.*
