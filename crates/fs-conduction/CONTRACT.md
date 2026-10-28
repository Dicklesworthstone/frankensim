# CONTRACT: fs-conduction

## Source-paired spatial sensible-enthalpy transport

`transient::enthalpy::material::MaterialEnthalpyTransport::from_isotropic_card`
binds Cp(T), rho(T) and scalar k(T) from one immutable material card at fixed
typed non-temperature conditions. The constructor explicitly declares an
isotropic transport model; source quantity dimensions alone do not establish
isotropy. All three properties must have explicit temperature validity over
the requested path. Missing-temperature producer scalars therefore cannot
silently become heating curves.

The adapter retains the original Cp/density chart, conductivity receipts and
a separately queried reference-density answer. Reference temperature must lie
inside the admitted span with identical other conditions, and its selected
density claim must match the caloric chart. Reference mass stays frozen as
temperature evolves. The caller explicitly selects phase, enthalpy datum and
caloric interpolation-error budget; source uncertainty is not propagated or
replaced by this numerical budget.

`prepare_uniform` uses the existing `EnthalpyBackwardEuler`.
`reference_material` supplies the existing heterogeneous owner with the same
chart and frozen density. Distinct material interfaces retain separate nodes
and explicit contact resistance; no mixed-node enthalpy law is invented.

The `material_transient` example accepts compiled NASA 316, first-heating
AA2024-T3 or NASA PEEK model packs, one/two layer thicknesses, mesh/time bounds
and optional contact resistance. Its declared experiment uses 300 K initial
and reference states, 360/300 K convective reservoirs, 1000 W/(m² K) end films
and a 20×20 mm cross-section with insulated sides. It emits accepted temperature profiles,
separate reservoir/contact heat rates and step/cumulative energy balances.
Geometry, mass and source conditions stay fixed; expansion, phase change and
microstructure evolution remain outside this experiment.

Four focused regressions cover source/reference pairing and refusal, the
independent variable-Cp heating balance, frozen mass, and nonuniform
equilibration with doubled k(T) without net heat creation.

```sh
cargo test -p fs-conduction --test material_enthalpy
```

## Sourced planar layered steady conduction

`layered::LayeredWall` solves an ordered equal-area planar stack from the existing
`ConductivityTable` owner. Each layer retains its source tables and usage
receipts; `PlanarThermalLayer::from_card` additionally records the actual
material-card identity. All contact and film resistances are explicit scenario
inputs in m2 K/W. Explicit zero means perfect contact or a prescribed surface
temperature; no boundary correlation is inferred.

The constant-conductivity solve returns signed heat flux, total heat rate,
U-value, and separate temperatures on both sides of every material/contact.
It checks each layer's actual face temperatures against that source table's
span and reports independently recomputed Fourier-law and boundary residuals.
The immutable input is unchanged on failure. Layer/knot caps bound work, and
the solver checks the existing executor cancellation gate between layers.

Varying conductivity is refused by the constant solve entry point. The separate
`LayeredWall::solve_temperature_dependent` entry point integrates every positive
piecewise-linear source segment with a scaled analytic inverse and bounded
monotone common-flux shooting. Each source table must support the entire reservoir
interval at its original fixed nonthermal context. Caller-declared iteration,
temperature and constitutive-flux limits gate the returned candidate. Cancellation
is checked between layers and every 128 source segments. Effective resistance and
transmittance belong to this boundary pair; they are not tangent conductances.
Seven cases in `tests/layered_nonlinear.rs` cover affine and interior-knot laws,
contacts/films, reversal, equal temperatures, original source receipts and bounded
refusal. No density/Cp, transient history, lateral spreading,
radiation, moisture transport, propagated material uncertainty or experimental
qualification is inferred. Geometric/process matching beyond the explicit
temperature span remains the caller's duty, including thickness-dependent
product declarations. G0/G1/G3/G4 cases in `tests/layered.rs` exercise the
analytical series-resistance law, contact jumps, exact source context,
material/thickness substitution, reversal, subdivision and transactional refusal.

## Growing deposition thermal network (issue #1)

`transient::deposition` accepts a resolved absolute-time schedule with stable
cell/material/body/segment IDs, incoming kg and J/kg, world-frame metre geometry,
and explicit thermal interfaces. Birth adds `m*h` once. Contact covers ambient
surface area; separation and reservoir release restore it. Old enthalpy and
source correspondence are retained. Interface thermal conductance is not a
mechanical bond strength. `enthalpy::deposition` binds the existing equilibrium
phase charts and executor Cx; it does not change the stationary P1 solver.

The executable `examples/deposited_strip.rs` deposits four passes of four cells
on a declared heated reservoir, cools, and applies a prescribed bed release.
It prints located temperature/exposure/energy samples and asserts supplied
mass, total energy closure, and first-pass reheating by the second pass. All
coefficients are illustrative. It predicts no release displacement or adhesive
detachment time. The example imports the exact std-only production kernel so
it can also be compiled directly when workspace dependencies are unavailable.

The sparse control-volume model subcycles forward Euler using a global chart
temperature-slope CFL bound and lands exactly on every event. A caller law must
supply a valid bound. Each advance returns an immutable candidate plus incoming
mass/energy, external heat, actual storage change, signed joule defect and work
counts. An absolute energy gate, chart domain, cell/link/event/substep budgets,
or cancellation refusal returns no candidate. Polling occurs every 256 entries
in thermal loops and between events/steps. Clone checkpoints preserve the event
cursor and state; bitwise restart equivalence is tested on the same timestep
grid and ISA. Adding query stops changes the integration grid and can change
the numerical solution. Schedule construction and cloning are budget-bounded
but not internally cancellable.

G0/G1/G3 tests cover birth/no preactivation, repeated endpoint idempotence,
successive-layer reheating and exposure, two-body exchange, analytic lumped
cooling/time refinement, uniform spatial partition, a nonuniform insulated
strip cosine-mode diffusion reference with independent spatial refinement
(and a separately halved time step), latent-plateau energy storage, sources, heterogeneous
laws, bed release, invalid topology/domain, and failed/cancelled trials.
The ignored cost test exercises 1,000 timed cells and 999 links; observations
are hardware-specific, not a performance guarantee. The core is std-only and
its exact source can be tested with `rustc --edition 2024 --test
crates/fs-conduction/src/transient/deposition.rs` when sibling repositories
prevent workspace resolution. This does not compile or validate the Cx/chart
adapter or the entire workspace.

No claim: automatic mesh/contact generation, conservative remeshing, geometric
deformation, thermal radiation, crystallization kinetics, viscoelastic history,
healing, cohesive failure, printer control, or calibrated FDM prediction. Each
cell is uniform; the caller must justify that approximation and conductances.
Geometry grows by adding supplied cells and changing interfaces, but does not
move. This is a useful reduced R1 foundation, not completion of umbrella #1.

> Beads `frankensim-extreal-program-f85xj.5.1`, `.5.3`, `.5.4`, and `.5.6`. The
> cooling vertical's steady heat-conduction solve, matching-P1 contact rung,
> first card-backed surface-radiation rung, and receipt-preserving PCB
> homogenization handoff.

## Purpose and layer

`radiation::PointwiseAmbientRadiation` supplies spatially resolved ambient
radiation on named natural P1 traces. It integrates
`epsilon sigma N_i (T_h(x)^4 - T_ambient^4)` and the exact temperature tangent
`4 epsilon sigma N_i N_j T_h(x)^3` with a positive degree-five triangle rule.
It retains the actual nodal temperature field; a hot spot is not replaced by
the fourth power of a patch mean. Every surface node must lie in the retained
emissivity claim's support. The immutable prepared mesh/partition is reused;
3-by-3 face tangents apply directly with linear storage and no dense matrix.

`EnthalpyBackwardEuler::advance_with_pointwise_ambient_radiation`, including
the heterogeneous wrapper, incorporates this operator in the original
enthalpy Newton residual, analytic tangent and Jacobi preconditioner.
Reference masses, sharp latent plateaus, conductivity, contact, original
convection and prescribed-temperature reactions remain in that single solve.
The result's `pointwise_radiation` report is a separate outward heat account
already included in `energy_residual_j`; only the original Robin heat goes
to a convection/coolant consumer. Physical history remains immutable on
refusal. Existing mean-patch APIs retain their explicitly reduced model.
The old adjoint binder refuses these endpoints, including zero-net-radiation
equilibrium, until its boundary ownership includes this operator.

Focused G1/G3 tests independently expand barycentric simplex moments, check
the analytic tangent and passive quadratic form, compare uniform cooling to
a scalar implicit balance, preserve latent-plateau storage, and solve a
nonuniform irradiated surface with prescribed reservoir reactions. This is
fixed-emissivity, fixed-geometry exchange with a black ambient reservoir and
view factor one; it does not model enclosure reflection, shadowing, participating
media, a temperature-dependent emissivity derivative, or supply an error bound
for the spatially discretized temperature field.

`PointwiseEnclosureRadiation` extends local P1 emission to closed diffuse-gray
exchange: the radiosity source is the exact area mean of `sigma T_h^4`, while
the nodal flux is the integral of `epsilon N_i (sigma T_h^4 - G_patch)`.
Irradiation is uniform within each declared view-factor patch; local emission
is resolved, with no inferred spatial visibility kernel or temperature-varying
emissivity. Geometry, emissivity and the bounded patch factorization are
prepared once. Sparse face/trace storage and a small patch response matrix
apply the full reflected tangent without a dense nodal matrix.

`advance_with_pointwise_enclosure_radiation` in both homogeneous and
heterogeneous enthalpy owners solves these terms with conduction, original
natural boundary heat, latent storage and prescribed-node reactions in one
Newton-FGMRES system. The reflected tangent is generally nonsymmetric; the
Picard preconditioner uses the positive local-emission diagonal, while the
reported exact diagonal retains the reflected term even when it is negative.
Each candidate independently gates the radiosity equations and closed heat
sum in watts. The final `pointwise_enclosure_radiation` report participates
exactly once in the joule balance; original Robin and Neumann accounts retain
only their authored heat. The base adjoint refuses these endpoints until its
owner includes the reflected spatial tangent. Tests in
`tests/enthalpy_enclosure/pointwise.rs` compare independent simplex moments,
reflected derivatives, nonuniform nodal balances, finite latent heat transfer,
and a strong-self-view case with a negative exact diagonal.

`transient::MaterialHeatCapacity` resolves constant density and specific heat
from one immutable material card at an explicit reference state, then pins both
selected sources at the endpoints of a declared positive temperature interval.
Only validity-wide scalar static laws are admitted: temperature-varying
density/Cp is not silently frozen. Those laws belong to the existing
reference-mass enthalpy path. The `isotropic_conductivity` pairing method uses
the same card, external conditions and interval with the original conductivity
table owner; an unrelated card refuses. Legacy hand-composed conductivity
models remain independent because they carry claim receipts rather than card
identity.

`MaterialHeatCapacitySource` additionally binds exact catalog density/Cp keys
and an explicitly named temperature axis, without aliases or re-entered
coefficients. Caller-listed source selectors admit only typed dimensionless
Boolean applicability flags constrained to matching singleton values; they
cannot drive property values or replace physical field/rate coordinates.
Their original values and quantity descriptors remain in every endpoint and
the paired `isotropic_conductivity_with_property` source. Both caloric and
conductivity sources must explicitly bound temperature; point-only or unknown
temperature data cannot authorize a finite heating trajectory. Bounded frozen
engineering references retain that model's domain and authority limitations.
`BackwardEuler::uniform_material` and `material_assigned`, and the history
runner's `MaterialUniform`/`MaterialById` capacity variants, retain these source
bundles while using the original lumped P1 capacity matrix. Every initial,
trial and returned temperature lies in every incident capacity source's domain,
including prescribed nodes and coupled endpoint evaluation. Out-of-domain
steps publish no new history. Frozen reference geometry and constant scalar
density imply fixed mass; this is a single-phase temperature formulation, not
expansion, latent heat or a material-error bound. G1/G3/G4 regressions in
`tests/material_heat_capacity.rs` check an independently derived tetrahedron
diffusion rate under a one-card swap, heterogeneous caloric heating, selected
source retention, and transactional refusals.

The original theta/Crank–Nicolson owner now uses the SAME prepared storage
through `march_with_capacity`. `transient::HeatCapacity` is the shared
declaration; `ThermalHistoryCapacity` remains a re-export of that exact type.
Uniform, per-element, material-ID and card-derived assignments therefore use
one preparation and source-domain policy. The theta path validates every
incident conductivity and capacity interval before and after prescribed
boundary lifting and at every accepted endpoint. Borrowed material source
objects retain the original receipts; a returned temperature field does not
itself own a new material card.

`march_with_capacity_and_source_scale_gradient` computes the existing exact
linear discrete adjoint using that same prepared capacity and an accepted
physical temperature for conductivity assembly. It does not query a
finite-domain material at an invented zero kelvin. The legacy gradient-only
API retains its historical zero-valued assembly reference; use the combined
entry point when source support requires an actual admitted trajectory.
Only the volumetric source is differentiated: Robin reference heat, Neumann
flux and Dirichlet values stay fixed. The source-only load is assembled by the
original P1 rule with an adiabatic boundary, avoiding both a second integration
rule and cancellation against large environmental loads. The material tests
include independent theta eigenmode/heating/gradient identities, heterogeneous
storage, Robin/Neumann boundary invariance, finite differences and endpoint
refusals. No phase-transition, nonlinear-material or support-boundary derivative
is claimed.

`incandescent` also exposes a deliberately lumped steady balance for a
declared filament operating point: `V²/R = εσA(T⁴ - T_ambient⁴)`. It composes
the crate-owned Stefan-Boltzmann constant with explicit electrical and geometry
inputs. It does not infer carbon resistance, emissivity, surface area, lead
conduction, gas conduction, lifetime, or visible-light efficacy.

Assemble and solve the STEADY heat-conduction equation on `fs-feec` tet
complexes, with anisotropic and temperature-dependent conductivity from
`fs-matdb` cards whose query receipts travel with the solve. Layer: **L3**
(FLUX).

Surface radiation has three explicit entry points. A linearized
surface-to-ambient model supplies
`h_rad = 4 ε σ T_mean³` over a caller-declared temperature-departure budget and
reports its pointwise discrepancy against `ε σ (T_s⁴ - T_a⁴)`; the departure
budget is not a derived accuracy bound. A
gray-diffuse enclosure solves
`J_i - (1 - ε_i) Σ_j F_ij J_j = ε_i σ T_i⁴`, then applies
`q_i = J_i - Σ_j F_ij J_j` as a frozen outward face flux in a declared outer
fixed point. Emissivity comes from an `fs-matdb` material card; this crate
generates the closed-form two-surface concentric-sphere geometry and supports
finite-mesh diffuse view-factor sampling through `radiation::view_factors`.
Caller-supplied analytic or external-QMC matrices remain admitted with their
own evidence. These radiation rungs
retain and enforce the selected emissivity claim's temperature-validity
interval.

Gray-diffuse radiosity solves use departures from one surface's blackbody
emission. Factored fourth-power differences retain small thermal contrasts;
the net heat uses the original emissive balance `epsilon * (sigma T^4 - G)`
in those same coordinates. This avoids erasing near-equilibrium exchange by
subtracting two large absolute radiosities, especially on reflective surfaces.
Absolute J/G remain diagnostic outputs and can round to the same value while
their physical exchange remains representable. View-factor row defects remain
in the shifted load and heat law; neither geometry nor enclosure heat is
normalized. The original matrix, pivot checks, material domains and prepared
factorization remain unchanged. Analytic plate tests cover adjacent floating
temperatures, reflective surfaces, exact isothermal equilibrium, and retained
imperfect-row leakage. This numerical correction adds no view-factor accuracy,
material calibration or continuum radiation claim.

`solve_with_convective_gray_diffuse_enclosure` explicitly composes that same
closed-enclosure exchange with original convection on each exact same-named
uniform Robin trace. A frozen radiative flux is represented by the algebraic
reference `T_air - q_rad/h`; it is never an air temperature. The original
radiosity/solid fixed point and final nonlinear residual gate remain the solver
owners. The result retains the exact combined boundary and separately reports
original convective heat for air coupling. Only that separated heat goes to
air; closed-enclosure radiation remains internal exchange. Partial/nodal
convective traces, non-Robin laws, nonpositive reservoir temperatures and
unrepresentable composition refuse. No coupled radiation adjoint follows from
the shifted boundary. `tests/radiation/convection.rs` checks an independent
two-slab heat balance, channel separation, exact replay and admission/refusal.

`estimate_diffuse_view_factors` samples area and cosine-weighted outward
directions using independent Owen-scrambled Sobol replicates. Every boundary
triangle remains an opaque first-hit obstacle, even when it is not a named
patch. The raw result retains hit counts, escaped/unassigned/back-facing
counts, empirical standard errors, actual intersection work, and an identity
binding the physical mesh, patch partition, seed and replicate counts.
`project_closed` is separate and explicit: it refuses any sampled loss,
symmetrizes area-weighted exchange, preserves the symmetric sampled zero
pattern, and uses bounded symmetric scaling to close rows. Excessive factor
adjustments, infeasible support and exhausted work refuse. Discrete conservation
does not establish geometric closure or a confidence bound. The current direct
triangle scan admits 2..=64 patches, 2..=32 replicates, up to one million
triangles, power-of-two samples and an explicit total intersection budget.

`solve_with_ambient_radiation` adds the nonlinear area-mean patch law
`epsilon sigma A (T_mean^4 - T_ambient^4)` to existing named uniform Robin
convection. It combines positive radiative secants with the original Robin
rows and solves through the existing material/contact producer until both
unrelaxed temperature updates and applied/nonlinear watt differences close.
Its conduction report retains the combined boundary actually solved; separate
original-convection rows and totals are the only heat carried to air. The
exact final boundary/configuration, card receipts, signed radiative heat and
cumulative solid/Krylov work remain available. Surface temperatures obey the
emissivity claim's validity; the external reservoir is not a surface-material
query. `SurfaceEmissivity::from_card_pinned` preserves exact claim selection.

An optional contact operator couples two duplicated, exactly coincident P1
traces without identifying their temperature degrees of freedom. For an
area-specific resistance `R'' > 0`, side A to side B heat flux is
`q_ab = (T_a - T_b)/R''`. The ordered `fs-matdb::InterfaceSystemCard` identity
and selected-property receipt travel with every reported interface flux.

Strong form:

```text
  −∇·(k(T, x) ∇T) = f             in Ω
                T = T_D           on Γ_D
        (−k∇T)·n  = q_n           on Γ_N
        (−k∇T)·n  = h (T − T_ref) on Γ_R
```

Weak form solved (find `T ∈ H¹` with `T = T_D` on `Γ_D`, for all
`v ∈ H¹_{0,Γ_D}`):

```text
  ∫_Ω ∇v·k(T)∇T dV + ∫_{Γ_R} h v T dA
      + ∫_{Γ_C} (1/R'') (T_a - T_b) (v_a - v_b) dA
      = ∫_Ω v f dV − ∫_{Γ_N} v q_n dA + ∫_{Γ_R} h v T_ref dA
```

Discretization: continuous P₁ Lagrange — the FEEC 0-form (Whitney vertex-hat)
space. Element geometry (signed volumes, the constant barycentric gradients
`∇λ_a`, the gradient Gram) comes from `fs_feec::element_geometry`; nothing here
re-derives it. Assembly stages triplets into `fs_sparse::Coo`.

### Sign conventions (declared once, applied everywhere)

- `q = −k∇T` is the heat-flux vector, W/m².
- `q_n = q·n` with `n` the OUTWARD unit normal: **positive `q_n` means heat
  LEAVING the domain**. `fs-scenario`'s `Physics::Thermal / BcKind::Neumann`
  row fixes the dimensions (W/m²) and not the sign, so this is this crate's
  declaration.
- `f` is a volumetric SOURCE, W/m³: positive adds heat.
- `EnergyBalance::dirichlet_in_w` is heat ENTERING through Dirichlet rows.

### Quadrature (the exactness claim, stated)

| integral | rule | exact for |
| --- | --- | --- |
| `∫_e ∇λ_a·K∇λ_b dV` | `V · gᵀ_a K g_b` | any constant `K` on the element — the P₁ gradients are constant, so this is not a quadrature at all |
| `∫_e λ_a f dV` | `Σ_b V(1+δ_ab)/20 · f_b` | `f` linear on the element (second-order otherwise) |
| `∫_F λ_a q_n dA` | `Σ_b A(1+δ_ab)/12 · q_b` | `q_n` linear on the face |
| `∫_F h λ_a λ_b dA` | `h̄_F · A(1+δ_ab)/12` | `h` constant on the face (second-order otherwise); `h̄_F` is the face mean |
| `∫_F h λ_a T_ref dA` | `h̄_F · Σ_b A(1+δ_ab)/12 · T_ref,b` | `h` face-constant AND `T_ref` linear on the face |
| `∫_F (1/R'') λ_a λ_b dA` | `(1/R'') · A(1+δ_ab)/12`, staged as signed A/A, B/B, A/B, and B/A blocks | face-constant `R''` on exactly matching P1 traces |

`K` is evaluated once per element at `T̄_e`, the element mean of the P₁
iterate. For P₁ that mean is the EXACT element average of `T_h`, so for `k`
linear in `T` the element integral `V·k(T̄_e)` equals `∫_e k(T_h) dV` exactly;
for a curved `k(T)` it is the midpoint rule, second-order in the element
diameter.

## Public types and semantics

| type | meaning |
| --- | --- |
| `ConductionMesh` | a tet complex with element geometry and its EXTRACTED boundary (faces incident to exactly one tet), each with area, outward normal, centroid; `new_region_owned` deterministically clones each used source vertex per incident region so an audited conformal labeled volume can expose independent, exactly coincident P1 traces while preserving each tet's ordered geometric coordinates and element order |
| `ScalarField` | uniform or per-vertex nodal data; the single carrier for every boundary value and source. There is no closure-valued field — a closure cannot be snapshotted or content-addressed |
| `ThermalBc` | `Dirichlet{temperature}`, `Neumann{outward_flux}`, `Robin{htc, t_ref}` |
| `ThermalBoundaryBuilder` / `ThermalBoundary` | a boundary PARTITION: regions may not overlap, leftovers refuse unless an adiabatic remainder is declared out loud |
| `InterfaceResistance` | a positive area-specific resistance in m² K/W selected from an ordered `InterfaceSystemCard`; retains both the card identity and exact property-use receipt. `from_card_pinned` preserves an explicit project claim pin without weakening property, dimensional, validity, finiteness, or positivity checks |
| `InterfaceSurface` / `ThermalInterfaces` | named, explicitly oriented matching-face pairs; every coincident boundary pair must be bound exactly once and may not also carry an external boundary condition |
| `ThermalResistanceTerm` / `SeriesThermalResistance` | named K/W terms and a deterministic series sum; stated half-widths combine conservatively by addition, while any unstated term keeps the complete band unknown |
| `ThermalResistanceTerm::slab_from_card` | resolves an exact static scalar conductivity key at the complete reference query point and computes `L/(kA)`. `ResistanceOrigin::BulkMaterialCard` retains card identity, usage receipt, original conductivity uncertainty and caller geometry. Resistance uncertainty stays `Unstated`: geometry and reciprocal uncertainty are not propagated. This is a reference-state series resistance, not a nonlinear temperature-span solve, transient heating law or inferred contact resistance. Two focused API regressions cover source retention, geometry scaling and property/state/geometry refusals |
| `ThermalResistanceTerm::slab_over_temperature_span_from_card` | steady uniform one-dimensional slab with fixed face temperatures, no generation and no lateral losses: integrates every crossed segment of one continuously supported scalar conductivity source and computes `R=L/(A*mean_T(k))`. `source_span_heat_rate_w` returns the corresponding hotter-to-colder heat rate using only those retained temperatures. Original source answers, uncertainties, receipts, card and geometry remain in `BulkMaterialTemperatureSpan`. The integral is exact for the admitted piecewise-affine law up to floating-point rounding; source/interpolation/geometry uncertainty remains unstated for resistance. No point-only extrapolation, nonlinear multilayer interface solve, transient storage, isotropic inference or experimental qualification is claimed |
| `ComponentPower` / `PowerMap` | per-component dissipation in W bound to a mesh vertex set, plus the declared system total. Names are unique and sorted; bound vertices are sorted and deduplicated so a repeated binding cannot take a double share while the totals still balance |
| `PowerMap::volumetric_source` | projects the map onto the nodal `W/m³` source AND returns the audit computed from that same field, so the reported delivered power is what the solve receives rather than a parallel estimate of it |
| `PowerMap::regional_volumetric_source` | the s93ej.2 generation companion of `from_region_ids`: one finite W/m³ per region id, lumped to nodes so delivered power is exactly `Σ_e f_e V_e`. Unmapped regions and length mismatches refuse. This is generation, not conductivity |
| `SurfaceComponentPower` / `SurfacePowerMap` | the boundary analogue: dissipation bound to boundary-face SLOTS, projected onto a nodal OUTWARD-flux field for a `ThermalBc::Neumann` region. `region_selector(mesh)` hands back the matching face predicate so the flux and the region cannot describe different face sets |
| `SurfacePowerAudit` / `SurfacePowerAuditRow` | the surface balance, reporting both the geometric `face_area_m2` and the map-wide `lumped_area_m2` the distribution actually divides by — they differ when footprints share vertices, and the lumped one is what keeps the total exact |
| `PowerAudit` / `PowerAuditRow` | the retained three-number balance — declared total, component sum, delivered total — with a per-component table of bound volume, applied density, and delivered watts. `total_half_width_w()` sums stated bands conservatively and returns `None` if any component is `Unstated` |
| `SurfaceEmissivity` | one dimensionless hemispherical-total emissivity selected from an `fs-matdb` material card; retains material state, card identity, exact property-use receipt, selected claim's optional closed temperature-validity interval, and the source uncertainty statement |
| `LinearizedSurfaceRadiation` / `LinearizedRadiationPoint` | a declared `4 ε σ T_mean³` Robin coefficient with a caller-chosen departure budget, enforced emissivity-card validity, full `T⁴` comparison, measured point discrepancy, and first-order emissivity-width propagation |
| `ViewFactorMatrix` / `ViewFactorEvidence` | an admitted closed-enclosure matrix whose row sums, area-weighted reciprocity, tolerances, evidence tag, and content identity are checked and retained; `concentric_spheres` supplies one closed-form two-surface geometry |
| `DiffuseViewFactorConfig` / `DiffuseViewFactorEstimate` | bounded mesh visibility sampling in `radiation::view_factors`; raw loss-aware factors and empirical replicate errors |
| `ViewFactorProjection` / `ProjectedViewFactors` | explicit sampled-matrix conservation correction with actual adjustment, work and geometry-QMC evidence retained |
| `RadiationSurface` / `GrayDiffuseEnclosure` / `RadiosityReport` | named non-overlapping P1 exterior traces, deterministic diffuse-gray radiosity, net heat rates, system residual, and enclosure energy closure |
| `CoupledRadiationConfig` / `CoupledRadiationSolution` | a budgeted, under-relaxed outer fixed point that applies one uniform frozen radiation flux per named trace to explicit adiabatic-remainder faces |
| `ConvectiveEnclosureSolution` | the same steady enclosure solver on exact original uniform Robin traces, with the actual combined boundary and separately identified convective heat for air coupling |
| `AmbientRadiationPatch` / `AmbientRadiationConfig` / `AmbientRadiationSolution` | nonlinear ambient patch radiation sharing the original uniform Robin trace, with separate convection/radiation accounting, heterogeneous material and finite contact preserved, independent temperature/watt gates and exact final binding |
| `ConductivityTable` | one scalar `k(T)` as sampled knots plus the `fs-matdb` receipts that produced them. `from_claims_pinned` samples one explicit claim pin at every knot and still refuses invalid or out-of-domain use |
| `LumpedThermalTransport` | reduced-body transport: explicitly provenance-free declared constants, sourced conductivity/emissivity, or sourced conductivity with radiation explicitly disabled; every active property receipt is retained, sampled models refuse extrapolation and bind to the same card as the phase curve |
| `ConductivityModel` | constant tensor, isotropic `k(T)`, or orthotropic `Σ_i k_i(T) e_i e_iᵀ`; every construction is checked symmetric and positive definite. `from_pcb_homogenization` consumes fs-matdb's immutable laminate result and retains one property-use receipt per copper/matrix material use |
| `MaterialId` / `MaterialTable` / `ElementMaterials` | a checked constitutive table plus one id per tet. `from_region_ids` maps labeled mesh regions onto that table. `bind_labeled_volume` is the mesher-agnostic s93ej.1 handoff: an independently audited `TetComplex` + region ids + region→material map become a `ConductionMesh` and assignment together. Assembly, Newton, `assemble_operator_with_element_materials` / `assemble_jacobian_with_element_materials`, and `element_heat_flux_assigned` use the named model at the element mean temperature. This is constitutive selection, not `element_scale`. Unknown/duplicate/empty/unmapped-region/length-mismatched assignments refuse. `solve_with_element_materials` is the product entry; `solve_with_element_materials_and_interfaces` is the joint contact path. Region materials cannot erase undeclared coincident faces. A bound assignment's table receipts and provenance travel on `ConductionReport` in place of the unused fallback. The uniform `ConductionProblem::material` path is unchanged. Two- and three-layer series, assigned `k(T)` Newton, assigned Jacobian-vs-FD, assigned IFT vs FD, one- and two-material contact coexistence, assigned matdb receipts, unused-table/mixed-provenance honesty, rotated orthotropic x-series, contrast/refinement/coalescing/map-swap metamorphics, and a labeled-volume lowering through `bind_labeled_volume` are covered; that last case is not yet the parent solve-stage E2E |
| `AssembledSystem`, `DofMap` | the full `n×n` operator and load, and the free/prescribed bookkeeping the Dirichlet elimination uses |
| `ConductionSolver`, `ConductionState` | the resumable nonlinear iteration and its snapshot payload |
| `ConductionSolution`, `ConductionReport`, `EnergyBalance`, `LinearSolveEvidence` | the field and everything established about how it was produced |
| `ResidualClaim`, `StallDiagnosis` (re-exported from `fs-solver`) | the types of `LinearSolveEvidence::reported` and `LinearSolveEvidence::stall`; re-exported so every type reachable through this crate's public fields is nameable from this crate (bead qvmwp) |
| `RobinFlux` | `EnergyBalance::robin_out_w` restricted to ONE declared Robin region: face count, area, area-weighted mean `h` / wall temperature / `T_ref`, and the integrated heat rate. This is what lets a conjugate driver attribute heat to a named convective trace instead of to the whole domain |
| `ThermalBoundary::region_for` | the region index owning a boundary-face slot, or `None` in the declared adiabatic remainder — the attribution the `RobinFlux` decomposition is built on |
| `ConductivityDesign` | the IFT adjoint hook for the LINEAR case: `dJ/dρ` over per-element conductivity multipliers. `new_with_element_materials` uses each tet's own linear `K_e`; any `k(T)` assignment is refused rather than linearized |
| `march_with_element_materials` | the linear theta-method march using each tet's assigned constant `K_e`. One-material assignment matches the uniform march; a two-layer assignment relaxes to the series interface. Any `k(T)` assignment is refused rather than frozen across a step |
| `source_scale_gradient_with_element_materials` | the discrete backward adjoint of `J = w · T^N` w.r.t. a uniform source scale, using each tet's assigned linear `K_e`. One-material assignment matches the uniform adjoint; a two-layer assignment matches independent finite differences. Any `k(T)` assignment is refused |
| `march_duty_cycle_with_element_materials` | the declared duty-cycle march using each tet's assigned linear `K_e`. One-material assignment matches the uniform duty march; a two-layer assignment relaxes to the series interface. Any `k(T)` assignment is refused |
| `ConductionError` | the total typed refusal set; `rule()` gives a stable slug |

`ConductivityTable::from_claims_at_query_points` and its pinned counterpart
accept an exact conductivity property key, an explicit temperature-axis name
and complete source query points. They preserve the declared property kind
through typed selection, retain typed coordinates in every receipt and require positive absolute
temperature in kelvin, an increasing sample grid, and fixed non-temperature
coordinates and axis descriptors. They infer neither axis aliases nor missing
source conditions. The single-axis convenience constructors use `T` only.
All four source-backed constructors use the canonical typed material source
line and require one continuously supported claim throughout the requested
temperature span. They retain every source curve knot and crossed surface-cell
edge, with the original fixed context, exact property key and replayable
receipts. Caller-added samples pin that same admitted claim. A coarse caller
grid cannot replace a source curve or surface slice with its endpoint chord;
support holes, source handoffs and tabulated-only temperature support refuse
even when both endpoints exist. Source coordinates named as temperature must
also carry temperature dimensions.
Declared curves remain explicitly authored models without database authority.
The real NIST 304/316/6061 source-compiler/store/conduction test in
`xtask/tests/matdb_pack_cli.rs` checks a manufactured steady temperature field,
heat flux, material substitution and context/domain refusals. Its IAPWS-water
counterpart checks 293.15–298.15 K and 333.15–338.15 K at exactly 100000 Pa,
with the liquid-phase context retained. Both use eight free interior degrees
of freedom, direct temperature/flux references and an energy-balance check.
This is numerical integration evidence for the declared interpolants;
water is conduction-only, with no buoyancy, fluid-motion, experimental,
transient-heat or phase-change qualification.

### `ThermalBc::from_scenario_row` — the Robin seam

`fs-scenario`'s Robin row is `expectation(Physics::Thermal, BcKind::Robin) ==
Value(HTC)`: the transfer coefficient ONLY. `BoundaryCondition` carries no
companion `T_ref`. This crate therefore requires `T_ref` to be named at the
lowering call — explicitly, or from `Environment::ambient_temperature`. That is
the seam the E05 correlation rung plugs into: a correlation computes `h`, the
reference temperature stays a declared property of the row, and the solve never
guesses which ambient a coefficient was fitted against.

Time signals, Chebyshev spatial profiles, and typed payloads are REFUSED at the
lowering boundary — the scenario lowering path is steady-only. The staged
`march_*`/duty-cycle drivers documented below DO step over declared capacity
data, but each step re-assembles and marches the SAME steady operator; they are
not a general transient frontend (see the no-claim section).

### The nonlinear driver: globalization and stop rule, both declared

Stop rule (`StopRule`), accepted when either clause holds:

```text
  ‖R(T)‖₂ ≤ residual_rtol · ‖b_f(T)‖₂ + residual_atol      (residual clause)
  ‖ΔT‖∞  ≤ step_atol on the last accepted step             (step clause)
```

Exhausting `max_iterations` is `ConductionError::NotConverged`, never a quiet
answer. The reason a run stopped is in `ConductionReport::stop_reason`.

Globalization (`Nonlinearity`):

- `FixedPoint { relaxation ω, max_backtracks }` — damped Picard,
  `T ← T + ω(y − T)` with `A(T)y = b(T)`. `A(T)` is SPD, so this path runs
  preconditioned CG.
- `Newton { line_search }` — Armijo sufficient decrease on `‖R‖₂`:
  accept when `‖R(T + αδ)‖ ≤ (1 − c α)‖R(T)‖`, shrink by `shrink` per
  backtrack, at most `max_backtracks`. The Jacobian carries the
  `(V/4) gᵀ_a K'(T̄_e) ∇T_h` term, which is NOT symmetric, so this path runs
  FGMRES.

Both share the **admissibility guard**: a trial iterate whose element mean
leaves the material's sampled temperature span has no merit value at all, so it
is rejected and the step shrunk. Extrapolating material data is never a fallback.

Preconditioner: ILU(0) where the factorization exists, Jacobi where it breaks
down. The choice is a pure function of the matrix, so a replay picks the same one.

### Residual claims are re-measured, not relayed

`fs-solver` reports CG's residual as a `ResidualClaim::RecursiveEstimate` — an
estimate that drifts in the dangerous direction. This crate therefore RECOMPUTES
`‖b − Ax‖₂/‖b‖₂` after every linear solve and gates on that recomputed number
(`LinearSolveEvidence::true_relative_residual` / `converged_true`). The
producing solver's own typed claim is carried verbatim in
`LinearSolveEvidence::reported` so the two are never confused.

### Goal-controlled regional maximum solves (`q61wp.73`)

`LinearGoalAnalyzer::solve_maximum_to_goal` uses the cached linear operator
and checked inverse to improve a field until its regional maximum's full
algebraic error meets an absolute Kelvin tolerance. The analyzer's zero-weight
preparation goal cannot admit this stopping decision: every candidate is
checked with `analyze_maximum`, including changes of the hottest vertex.
The method shares the existing linear-goal defect-correction loop, its one
primal iteration budget, bounded check cadence, best-field retention and
post-observer cancellation gate. It returns `LinearMaximumSolve`, with an
explicit missing-bound, exhausted-budget or no-progress outcome when needed.
No dual or stability solve is repeated during correction. This is a
stored-system goal solve, not a continuum accuracy or energy-balance claim;
the returned field is not a `ConductionSolution`. The existing
`tests/linear_maximum.rs` checks the zero-weight counterexample against a
dense solve, budget limits, prescribed regions and cancellation/retry.

`polish_linear_maximum` connects this correction to an existing physical
`ConductionSolution`. It independently reassembles the baseline and requires
its actual residual to meet the original report's finite threshold. A changed
best candidate is adopted only if independent reassembly meets that same
threshold; otherwise the baseline field and its maximum analysis survive with
a physical-gate refusal. Even a budget-limited improvement may be accepted,
while `goal_met` refers strictly to the returned field. Energy, Robin/contact
fluxes and final residual are recomputed by the same report helper as the
ordinary solve. Historical nonlinear iterations and linear evidence are
preserved; correction work has separate counters. The threshold is declared
caller policy, not authenticated provenance. Nonlinear conductivity refuses;
no frozen-radiation or air-feedback guarantee is inferred. Focused tests in
`tests/algebraic_goal/polish.rs` cover real residual-accepted correction, an
independent contact/Robin heat-balance oracle, missing inverse evidence,
zero-work retention, false baseline summaries and cancellation.

`LinearGoalAnalyzer::new_for_maximum` additionally supports a bounded inverse
proposal when positive diagonal scaling cannot verify the stored operator.
This covers consistent contact and Robin operators whose comparison matrices
are not strictly dominant. At most 256 free unknowns are admitted, and both
`n*n` stored proposal entries and `n*(n+nnz)` verification visits must fit the
existing residual nonzero limit. Every column solve shares the remaining
`max_stability_iterations` allowance with the original scaling proposal;
`stability_iterations` reports their combined actual work. Each analysis
rechecks `I - A R` against the same retained matrix before using the inverse.
The proposal is never treated as exact, and exhausted work or insufficient
verification leaves the bound unavailable. Existing finite inverse bounds and
the ordinary linear-goal constructor keep their prior behavior. The current
affine-feedback checker retains its separate solid-inverse verification route;
this fallback establishes no additional coupled-error authority.
`tests/algebraic_goal/inverse.rs` checks a strong consistent-contact operator
against an independent dense solution, real maximum-goal correction, shared
iteration limits and structural admission.

When neither dominance nor the dense proposal proves the inverse (typically
above 256 free unknowns on obtuse-dihedral P1 meshes), `new_for_maximum`
prepares `fs_solver`'s `VerifiedSparseInverse`: bounded interval LDL
elimination with positive pivot intervals and an outward `A-H` asymmetry
check. It is charged to the same `max_stability_iterations` allowance (one
iteration = one matrix pass of interval updates; zero allowance does no
inverse work), and its stored entries are capped by the residual nonzero
limit. Only an established bound is retained. Every later analysis, including
the coupled affine-feedback maximum (which uses it before the shifted-Gram
spectral fallback), re-checks the exact stored CSR before using its norm.
On the 43,840-tet heatsink ladder rung this proves `||A^-1||_inf <= 9.1e6` over
9,295 rows with about 1.2M peak entries (measured 2026-09-27).

## Invariants

1. **Operator symmetry and definiteness.** The Dirichlet-reduced conduction +
   Robin operator is symmetric to `1e-12` relative and has a positive Rayleigh
   quotient on the probed set. Evidence:
   `tests/conformance.rs::assembled_operator_is_symmetric_and_positive_definite`.
2. **Element nullspace.** Every element matrix annihilates the constant vector
   (`Σ_b K^e_ab = 0`) — the discrete `∇·(K∇c) = 0`. Evidence:
   `tests/conformance.rs::element_stiffness_rows_sum_to_zero`.
3. **Elimination consistency.** The free rows of `A·T_full − b` are the reduced
   system's residual, and prescribed vertices survive the solve bit-exactly.
   Evidence: `tests/conformance.rs::dirichlet_elimination_matches_the_full_residual`.
4. **Boundary partition.** Regions never overlap; leftover faces refuse unless
   an adiabatic remainder is declared. At a vertex shared by two Dirichlet
   regions the LOWEST declared region index wins — a deterministic tie-break,
   not a traversal artefact.
5. **Jacobian correctness.** The assembled Newton Jacobian agrees with central
   differences of the residual to `1e-5` relative on a non-uniform iterate.
   Evidence: `tests/conformance.rs::newton_jacobian_matches_central_differences`.
6. **Material provenance is never invented.** `ConductivityModel::provenance()`
   returns `MatdbReceipts` only when EVERY component came from an `fs-matdb`
   query with a retained receipt; a declared constant says `Declared`.
   Evidence: `tests/conformance.rs::matdb_receipts_travel_with_the_solve`,
   `::declared_conductivity_never_claims_matdb_provenance`, and
   `tests/pcb_homogenization.rs`.
7. **No implicit extrapolation.** Sampling conductivity outside a claim's
   validity or evaluating outside its sampled grid refuses. A selected
   emissivity claim's temperature interval is retained after query:
   linearization mean/ambient/evaluation temperatures, radiosity driving
   temperatures, and newly solved coupled surface means all refuse outside it.
   Evidence:
   `tests/conformance.rs::material_queries_outside_validity_refuse` and
   `tests/radiation.rs::radiosity_refuses_temperatures_outside_retained_emissivity_validity`.
8. **Energy balance.** `closure_w = source + dirichlet_in − neumann_out −
   robin_out` equals `−Σ_{free} r_i` identically; on a converged solve it is
   below `1e-9` relative. Evidence: `tests/conformance.rs::energy_balance_closes`.
9. **Contact is complete and explicit.** Exact coincident boundary triangles
   are never treated as perfect contact and never silently left as an
   adiabatic gap. Every pair requires one named surface, a positive in-domain
   `area-specific-thermal-contact-resistance` claim, and an explicit
   `ThermalInterfaces` solve/assembly path. Reusing a bound object with a
   different face set, vertex mapping, area, or boundary condition refuses.
   Evidence: `tests/contact.rs::missing_interface_card_and_missing_binding_refuse_typed`,
   `::surface_declaration_and_binding_guards_each_refuse_their_own_defect`
   (blank/duplicate names, empty surfaces, duplicate/self/out-of-range/
   noncoincident pairs, and external-boundary overlap),
   `::coincident_geometry_refuses_nonopposing_and_multi_owner_traces`, and the
   four `a_bound_interface_refuses_*` tests that isolate face-set, boundary,
   exact-vertex-correspondence, and area-bit drift.
10. **Radiation authority is explicit.** Every emissivity retains its material
    card, property-use receipt, and selected temperature-validity interval;
    every view-factor matrix passes finite range, row-closure, and
    area-weighted reciprocity admission before a radiosity solve, and its
    declared areas must match the bound mesh traces. Enclosure exchange closes
    globally on an unequal-area three-surface fixture whose residual cannot
    collapse to exact pairwise negation. The outer conduction coupling stops
    only on its declared surface-temperature rule, refuses an exhausted budget,
    and cannot overwrite a Robin or other occupied base row. Evidence:
    `tests/radiation.rs`, including the named row-closure, area-mismatch,
    singular-limit, exhausted-budget, and occupied-Robin refusal tests.
11. **Component power is delivered EXACTLY, not to discretization order.**
    Summing the exact P1 element load `∫ λ_a λ_b dV = V(1+δ_ab)/20` over `a`
    gives `V/4 · Σ_b f_b`, so the power the assembled operator injects is
    `Σ_a f_a w_a` with the lumped nodal volume `w_a = (1/4) Σ_{e ∋ a} V_e`.
    Distributing a component's watts by those same weights makes its delivered
    power equal its declared power identically, at every mesh resolution.
    Distributing a per-element density onto nodes instead would lose power
    wherever an element is shared at a component boundary — and the loss would
    be invisible, because it looks exactly like discretization error.
    Overlapping components are ADMITTED and still balance, since the delivered
    functional is linear in the nodal field. Evidence: `tests/power.rs`.
12. **A power map that does not add up refuses, and names the rows.** The
    declared system total is a checked input, not a comment: a mismatch beyond
    the caller's relative tolerance returns the per-component table inside the
    refusal, so the reader sees which row is wrong. A forgotten component is
    the commonest data-entry error in a thermal model and the cheapest to
    catch. Evidence: `tests/power.rs`.
13. **Surface power is delivered exactly too, by the same argument one
    dimension down.** The exact P1 face load `∫λ_aλ_b dA = A(1+δ_ab)/12` has
    row sum `A/3`, so the delivered surface power is `Σ_a q_a s_a` with the
    lumped nodal AREA `s_a = (1/3) Σ_{f ∋ a} A_f`. The weight is computed over
    the WHOLE map rather than per component, which is what keeps adjacent
    footprints sharing a vertex exact: the assembled Neumann load integrates
    over the whole region, so the distribution weight must too.
    Evidence: `tests/power.rs`.
14. **Boundary-face SLOTS and canonical FACE indices are different numbering
    spaces.** A component binds slots into `ConductionMesh::boundary`; a
    region predicate is handed a `BoundaryFace` whose `face` field indexes the
    complex's canonical face table. Both are `usize`, so confusing them
    compiles and silently selects the wrong faces.
    `SurfacePowerMap::region_selector` takes the mesh precisely so it can do
    that translation once, in one place. Evidence: `tests/power.rs` asserts
    the selector reproduces exactly the claimed slots — it caught this
    confusion during development.
15. **The Robin decomposition partitions the convective total without
    redefining it.** `report.robin_fluxes` carries one row per declared Robin
    region owning at least one boundary face, in declaration order, each
    integrating `∫ h (T_h − T_ref) dA` over that region's faces with the SAME
    per-face rule the balance uses. `EnergyBalance::robin_out_w` keeps its
    historical running-accumulator value bit-for-bit — it is NOT re-derived as
    a sum of these parts — so the whole and the sum of parts agree only up to
    floating-point summation order, and the test pins that gap rather than
    asserting it away. A Robin region matching no face emits no row: it moves
    no heat and has no area to average over, so inventing zero-valued means
    for it would be a fabricated measurement. Evidence:
    `tests/robin_flux.rs::each_robin_region_matches_its_own_closed_form_heat_rate`
    (per-region closed forms on a doubly-convected generating slab whose two
    ends differ by ~20x, so a swapped attribution fails by magnitude),
    `::the_two_end_heat_rates_account_for_the_whole_generated_power`,
    `::the_decomposition_sums_to_the_whole_domain_total_up_to_summation_order`,
    and `::the_rows_cover_exactly_the_faces_whose_condition_is_robin`.

## Error model

Every fallible path returns `ConductionError`, a total typed enum with a stable
`rule()` slug, the offending quantity, and — where one exists — the fix. No
panics cross the crate boundary: `fs-feec`'s degeneracy ASSERTION is pre-empted
by a `DegenerateElement` refusal, and `ConductionMesh::new` validates sizes and
finiteness before any kernel runs.

Structurally impossible internal states are `expect`/`debug_assert` (e.g. the
converged step caching its own assembled system); they are not reachable from
any public input.

Three refusals exist specifically to stop a plausible wrong answer:

- `SingularPureNeumann` — with neither a Dirichlet nor a Robin row the steady
  operator is singular up to a constant. A Krylov method would return something
  that looks like a temperature field. This refuses instead.
- `LinearSolveFailed` — raised on the RECOMPUTED Euclidean residual, not the
  solver's recurrence estimate.
- `Interface` — raised when a resistance card is missing, ambiguous,
  dimensionally wrong, non-positive, when exact face ownership is incomplete,
  or when a mesh-bound interface object no longer matches the supplied mesh
  and boundary.
- `Radiation` — raised for a missing or inadmissible emissivity card, an
  emissivity-use temperature outside the selected claim's retained validity, an
  out-of-domain linearization, an invalid or nonreciprocal view-factor matrix,
  overlapping/wrong-mesh surface bindings, an occupied boundary row, a singular
  radiosity operator, or an exhausted outer fixed-point budget.

## Determinism class

**`DeterministicPerIsa`** — bit-identical across runs, worker counts, and build
modes on ONE ISA. Cross-ISA bit-identity is NOT claimed.

Why it is achievable, mechanism by mechanism:

- Assembly stages into `fs_sparse::Coo`, whose canonicalization is a pure
  function of the staged `(row, col, sequence, value)` tuples; element and face
  loops traverse in canonical index order regardless of tile size, so tiling
  cannot move a bit.
- Every inner product and norm goes through `fs-solver`'s fixed-shape chunked
  reduction (shape depends on length only, never on threads).
- Preconditioner selection is a pure function of the matrix. Every Krylov method
  starts from `x₀ = 0`. There is no RNG anywhere in this crate.
- Conductivity tables are piecewise-linear with a declared tie-break (a
  temperature landing exactly on an interior knot uses the segment to its
  RIGHT), so `eval` and `derivative` agree on which piece is in force.
- Contact surfaces, face pairs, and series-resistance terms are sorted by
  stable names/slot keys before assembly or summation. Input permutation cannot
  change their floating-point accumulation order.
- View-factor rows and radiation surfaces retain admitted order. The radiosity
  solve uses deterministic partial pivoting with the lowest row as the exact
  tie-break, and the outer fixed point traverses surfaces and faces in retained
  order. External QMC metadata may be admitted, but this crate runs no RNG.
- The library calls NO platform transcendental. Only `sqrt` (IEEE-754
  correctly rounded, hence exempt), `fs_math::det::{sin, cos}` in
  `fixtures::annulus_sector`, and `fs_math::det::exp` in
  `lumped::response_at`.

Why the class is not stronger: the cross-ISA claim requires a G5 audit on both
reference ISAs per `docs/DETERMINISM_CLASSES.md`, and none has been run for this
crate. The code is already libm-free on its default path, so the promotion is a
matter of running the audit, not of changing the arithmetic.

Evidence at the claimed class: `tests/conformance.rs::replay_is_bitwise_identical`
(same solve twice, bitwise) and `::snapshot_resume_reproduces_the_uninterrupted_trajectory`
(pause → seal → restore → resume reproduces the straight run bitwise).

## Cancellation behavior

Element and boundary-face loops run in tiles of `ASSEMBLY_TILE` (512) and poll
`Cx` at every tile boundary; contact-face assembly uses the same tile bound,
the nonlinear driver polls once per iteration, and radiation polls during
radiosity factorization and every outer fixed-point iteration. A
cancelled run returns `ConductionError::Cancelled { stage, at }` naming the
stage (`"assemble-elements"`, `"assemble-boundary"`, `"assemble-jacobian"`,
`"assemble-interfaces"`, `"nonlinear-iteration"`,
`"radiation-radiosity-assemble"`, `"radiation-radiosity-factor"`, or
`"radiation-coupling"`) and the tile or iteration index that was about to run.

Cancellation DRAINS: the partially staged `Coo` is dropped, and a cancelled
iteration leaves `ConductionState` untouched — no half-applied update can
survive. Evidence: `tests/conformance.rs::cancelled_assembly_is_a_structured_refusal`,
`::cancellation_mid_iteration_drains_with_a_structured_refusal` (which asserts
the state is unchanged).

Legacy-v1 resumability: `ConductionState` implements
`fs_exec::solver::LegacySolverStateV1`
(`TYPE_ID_V1 = 0xf5c0_4d75_c710_0001`, `SCHEMA_VERSION_V1 = 1`), and
`snapshot()`/`restore()` name `LegacySnapshotV1Adapter` explicitly. Restore is
the adapter's unbounded compatibility path: it validates magic, type id, schema
version, length, and checksum BEFORE the payload decoder runs, but accepts no
caller-pinned exact root and has no cancellation probe. Admission/migration
owners must bound and authenticate bytes before calling it. This does not claim
v2 identity or migration authority. A flipped byte, a truncation, and an append
are all refused. Evidence:
`tests/conformance.rs::snapshot_envelope_refuses_tampered_bytes`.

## Unsafe boundary

None. `unsafe_code` is denied workspace-wide and this crate registers no capsule.

## Feature flags

None. Everything here is `[S]` solid work on the default path.

## Conformance tests

- `tests/conformance.rs` — G0 algebraic laws, G4 cancellation drills, G5
  determinism and snapshot round-trips, and the typed refusal surface (22 tests).
- `tests/mms.rs` — the G1 battery: seven primal manufactured-solution ladders gated by
  `fs_mms::OrderGate`, plus the declared `MmsMatrix` and a complete
  `fs-vvreg` Level-A target crosswalk whose gaps carry reasons.
- `tests/analytic.rs` — slab (Dirichlet–Dirichlet), slab with a uniform source,
  slab (Dirichlet–Robin), a rectangular affine patch, cylindrical shell,
  spherical-shell patch, and an adiabatic-tip straight fin, each with a
  declared envelope and a stated reason for its size; a complete Level-A
  analytic crosswalk keeps unsupported rows visible.
- `tests/contact.rs` — G1 two-slab series-resistance comparison against the
  Level-A `0.3 K/W` row. The registry tolerance is applied to the
  SOLVE-DERIVED end-to-end resistance `deltaT/Q`, and a deliberate value beyond
  that envelope is rejected; the separately assembled three-term budget is
  retained as arithmetic evidence, not mislabeled as the solve result. G0
  covers receipt replay, uncertainty, overflow, missing-card/binding, every
  surface-declaration and geometric-ownership refusal, and isolated bound-state
  drift in face set, boundary condition, exact vertex correspondence, and area
  bits. A three-slab/two-interface fixture declares the two named surfaces in
  both orders, requires canonical name order, and compares jump-bearing flux
  records bitwise while distinct resistances make a swap observable.
  Mapped-resistance coverage (`f85xj.12.9`): per-face application; the
  SORT-PERMUTATION test — a reversed declaration, pairs and values reversed
  together, must bind to the same per-face values as the forward one, which is
  what proves the deterministic face sort cannot separate a face from its
  resistance; the metamorphic all-equal map runs BOTH nonzero-jump solves and
  reproduces the uniform temperature field and interface report BITWISE;
  doubling `R''` runs both solves and moves heat rate from the independent
  `100 W` series value to `75 W` while increasing the interface jump; and
  refusals cover a count mismatch, mixed cards, a blank measured rationale,
  and non-finite/non-positive measured values.
- `tests/power.rs` — component power maps (`f85xj.5.9`): delivered power equal
  to declared power at `1e-12` and INVARIANT ACROSS MESH REFINEMENT (n = 2..5),
  which is what separates an arithmetic identity from a discretization
  property; overlapping components superposing without loss; a duplicated
  vertex binding collapsing instead of taking a double share; the analytic
  cross-check that a whole-domain map reproduces the equivalent
  `ScalarField::Uniform` SOLVE to `1e-9`, pinning the power path to the
  crate's existing source path; audit refusal carrying the per-component
  table; binary64 rounding admitted inside tolerance while a forgotten
  component refuses; one `Unstated` component making the total band `None`;
  the full admission-refusal surface; and declaration-order independence.
  Surface coverage (`f85xj.5.14`): delivered surface power exact at `1e-12`
  and across refinement; adjacent footprints sharing vertices still balancing,
  with an assertion that the fixture genuinely exercises the shared-vertex
  case; the analytic cross-check that a full-face map reproduces the
  equivalent uniform `ThermalBc::Neumann` SOLVE, which also pins the sign
  (heating must raise the body above its Dirichlet wall); the region selector
  reproducing exactly the claimed slots; and refusals for a face claimed
  twice, an out-of-range slot, and a total that does not add up.
- `tests/power_table_e2e.rs` — file-to-solve ingestion (`f85xj.11.7`), the only
  battery in this crate that crosses a crate boundary: a tracked power-budget
  export is read FROM DISK, parsed by `fs_io::power_table` (a DEV-dependency,
  so the production graph stays free of the importer), lowered through an
  explicit name-to-region catalog, and solved. It pins delivered power equal to
  the file's declared total at `1e-12`, the per-row wattage surviving the parse
  bit for bit, the deterministic by-name component order regardless of
  spreadsheet row order, one `Unstated` row still making the aggregate band
  `None`, the assembled `EnergyBalance` agreeing with the audit including the
  SIGN of the Dirichlet rail, and the hot spot landing inside the largest
  dissipator's footprint — which is what a mis-bound name-to-region seam would
  break and a balance check alone cannot catch. G3: identical digits declared
  `mW` instead of `W` move every nodal temperature RISE by exactly 1000, which
  measures what the ingestion-boundary unit directive prevents, since this
  crate reads no dimension (`COMPONENT_POWER_DIMS` is consumed nowhere).
  Refusals: a declared total disagreeing with its rows must surface `PowerMap`'s
  balance refusal rather than being corrected at parse time, and a component
  with no footprint refuses rather than being skipped — skipping would remove
  real dissipation and still balance.
- `tests/transient.rs` — theta-method marching (`f85xj.5.11`): the EXACT
  uniform-adiabatic-heating identity `dT/dt = f/(rho c_p)` at machine
  precision across four (theta, step) combinations, which pins the
  capacitance assembly and the stepping algebra together — a lumped or
  mis-scaled mass matrix would drift; the capacitance row-sum equalling
  `rho c_p V` directly; marching to large time recovering the independently
  tested STEADY solution; a boundary-inconsistent initial condition being
  lifted rather than carried; observed temporal order through
  `fs_mms::OrderGate` (`f85xj.5.15`) on a ladder against a much finer
  reference on the SAME mesh, so the spatial error cancels and what is fitted
  is purely temporal — backward Euler observed 0.994 against a theoretical 1,
  Crank-Nicolson 2.008 against 2, plus a negative control gating the
  first-order ladder as second order and requiring the refusal; the 1-D slab
  step response against its series solution (below); refusals for `k(T)`,
  `theta < 0.5`, non-positive step or capacity, and a mismatched initial
  vector; determinism across independent runs; and a cancelled march
  publishing no partial solution.
  Adjoint coverage: the source-scale gradient passing the crate-wide
  `fs_adjoint::verify_gradient` gate under BOTH backward Euler and
  Crank-Nicolson — with `informative_directions` asserted so the pass cannot
  be vacuous — plus a closed-form check on the uniform adiabatic case where
  `d(mean T)/ds = f t / (rho c_p)` exactly, so any drift is adjoint algebra
  rather than discretization; and refusals for `k(T)` and mismatched
  functional weights. G1 corpus binding: the `fs-vvreg` Level-A
  `thermal-a-lumped-transient` row (`theta/theta0 = exp(-1)` at one time
  constant) resolved at the row's DECLARED Biot ceiling, plus the check that
  matters more — the discrepancy SHRINKS with Biot, so the agreement is
  attributable to the regime the row declares rather than to a tolerance that
  happened to fit.
- `tests/lumped.rs` — the reduced rung (`f85xj.5.13`): the closed-form
  response checked against the ODE IT CLAIMS TO SOLVE by numerical
  differentiation, rather than against the exponential formula it already is;
  endpoint and time-constant properties; the crate's Biot ceiling asserted
  equal to the corpus row so the duplication cannot drift; the gate admitting
  inside and refusing outside, with the refusal required to point at
  escalation; a loosened gate reporting its own ceiling; extraction
  round-tripping the conductance that produced a rise; and THE PAIR THAT
  JUSTIFIES THE GATE — the rungs agreeing within 5% inside the regime, and
  DISAGREEING by more than 10% outside it (measured through a deliberately
  loosened gate), with the physical direction asserted: a non-isothermal body
  retains more heat than lumping predicts. Comparison uses a VOLUME-WEIGHTED
  mean; a plain vertex average is dominated by boundary nodes on a coarse
  cube (56 of 64 at n = 3) and would have made the comparison meaningless
  while looking plausible.
- `tests/duty.rs` — duty-cycle drivers (`f85xj.5.12`): the analytic per-segment
  energy integral; declared interpolation across constant and ramp segments;
  the past-the-window REFUSAL (with the refusal required to say why);
  Crank-Nicolson with aligned steps delivering the declared energy EXACTLY;
  backward Euler showing a real first-order residual that halves with the
  step; the deliberate counter-case where a SYMMETRIC cycle cancels that error
  to zero, pinned so a zero residual is never read as proof of exactness; a
  ragged window refusing; constant load integrating exactly under either
  scheme; the sampled peak lagging the burst (thermal mass doing its job);
  step-quantized time-above-limit monotone in the limit; and refusals,
  determinism, and cancellation.
- `tests/radiation.rs` — G0 card, retained-validity, dimensional, view-factor,
  reciprocity, area-binding, cancellation, and refusal checks, including the
  singular radiosity limit, occupied Robin overlay, and exhausted outer budget;
  an unequal-area three-surface algebra fixture exercises pivoting, multi-term
  back-substitution, and non-tautological enclosure closure. The two-surface
  radiosity and two-solid outer fixed point are algebra/coupling harnesses
  checked against closed-form radiosity and an independent scalar Newton
  equilibrium. Their literal `F₁₂ = F₂₁ = 1` matrices are caller-declared and
  geometrically false for the cube-face and separated-finite-slab traces, so
  these tests are not analytic view-factor kernels or G2 geometry benchmarks.
  `concentric_spheres_emit_the_closed_form_reciprocal_enclosure` executes the
  closed two-surface concentric-sphere formula and pins its reciprocal areas,
  row closure, self-view, and analytic evidence tag; its companion refusal test
  rejects nonphysical radii, nonrepresentable (underflowing or overflowing)
  areas, and a positive factor that underflows.
- `tests/geometry_view_factors.rs` — finite-square analytic integral in a
  rectangular cavity made from real tetrahedral panels, unequal-area
  conservation projection and radiosity closure, an opaque unassigned blocker,
  outward escape, named-stream replay, and cancellation/work/refusal cases.
  These target geometric sampling and discrete algebra, not physical calibration
  or a certified error interval.
- `tests/adjoint.rs` — the linear IFT gradient against central differences
  through `fs_adjoint::verify_gradient`, plus G1 manufactured P1 and P2 dual L2
  ladders with per-rung discrete primal/dual identity checks (4 tests).
- `tests/heterogeneous.rs` — regionwise constitutive assignment: one-material
  parity, two- and three-layer series slabs, anisotropic and rotated
  orthotropic x-series, assigned `k(T)` Newton, assigned Jacobian vs central
  differences, assigned linear IFT vs FD, labeled-mesh handoff, contact
  coexistence, matdb receipts, unused-table/mixed provenance, contrast
  monotonicity, identical-k coalescing, swapped region maps, refinement,
  mapping/span/length refusals, and the claim that a material map does not
  erase undeclared coincident faces.
- `tests/pcb_homogenization.rs` — G3 model/QoI handoff: fs-matdb's
  receipt-backed laminate tensor enters the production steady solve; the
  bounded lower/nominal/upper in-plane conductivities bracket the heat-flux
  QoI; the reference copper stack exhibits the intended order-of-magnitude
  difference from an isotropic-FR4 placeholder; material receipts remain on
  the nominal solve report. This is model integration evidence, not physical
  validation.
- `tests/level_b_crosscode.rs` — the executing Level-B cross-code lane
  (bead f85xj.4.3): all four `fs_vvreg::thermal_level_b` cases are rebuilt
  from their typed catalog definitions, solved with this crate, and compared
  probe-by-probe against frozen scikit-fem/SuperLU references within the
  declared envelopes; the canonical BLAKE3 mesh-identity hash must match the
  external run's record (bit-identical Kuhn mesh across languages); an
  interpolation utility is proven exact on a synthetic affine field and
  consistent with nodal extraction; the nonlinear case must take more than
  one Picard step; each solve's energy closure must stay citable. Agreement
  is same-discretization independent-implementation evidence with Estimated
  colour — not discretization error, not physical validation.

Stable evidence locator (historical wording): `Every test prints a JSON-lines verdict`.
The precise current scope is narrower: every analytic comparison and MMS gate
prints a JSON-lines verdict carrying the numbers it asserted, so a failure in
those batteries reproduces from its log line alone.

### Level-A registry binding is executable, not retained authority

The analytic and MMS tests dev-depend on `fs-vvreg` and resolve canonical
Level-A case rows at runtime. Eight analytic cases use the catalog
parameters and reference values directly: the two slab fluxes, uniform-source
center rise, rectangular affine probe, cylindrical- and spherical-shell
conductances, `mL=1` fin efficiency, a two-slab matching-P1 contact network
whose three `0.1 K/W` terms reproduce the `0.3 K/W` series reference. The
infinite-parallel-plate row is intentionally REFERENCE-ONLY here:
`ViewFactorMatrix::infinite_parallel_plates` stores `F12 = F21 = 1` literally,
so admitting that matrix and consuming it in radiosity is not an
executing-kernel reproduction of the geometry formula. Four primal P1 L2
ladders take their
theoretical order and two-sided gate from the catalog: isotropic Dirichlet,
combined anisotropic temperature-dependent conductivity, mixed Neumann, and
Robin. A fifth P1 L2 ladder resolves the heat-adjoint target through
`fs_adjoint` and checks the discrete dual identity on every mesh. The combined
primal case uses a rotated orthotropic principal frame with three distinct
linear `k_i(T)` laws and a chain-rule manufactured source, so it exercises
off-diagonal `K(T)`, nonzero `dK/dT`, and nonlinear iteration in one executing
ladder. The P2 isotropic Dirichlet target executes through the existing
`fs_feec::highorder::simplex::SimplexSpace` order-2 stiffness, load, boundary
mask, and L2 machinery. That is high-order FEEC-kernel evidence; it does not
change this crate's P1 material/boundary frontend into a general P2 thermal
frontend. The P2 heat-adjoint target uses the same high-order operator, solves
the transposed system through `fs_adjoint::ift_gradient_matfree`, and checks the
discrete primal/dual identity on every mesh.
The analytic and MMS crosswalks are complete over their respective catalog
partitions, so adding, removing, or silently renaming a Level-A row fails the
battery and every absent binding retains a reason.

This is a **test-time execution binding**, not a corpus receipt. The catalog's
`1e-12` analytic tolerance establishes reproduction of the closed-form
definition; it is not substituted for this solver's geometry, discretization,
or model envelopes. The test verdicts label that distinction explicitly, and
no ladder or machine fingerprint is persisted into `fs-vvreg`. Consequently
the registry query remains numerical `NoClaim`, all Level-A physical caps
remain `Estimated`, and the four rows not bound in this crate are still
reference/target-only here. Two Nusselt rows execute separately in
`fs-convection`, and the normalized lumped-transient row executes through the
first-order generalized-alpha paths in `fs-time`. Across those owning crates,
18 of 19 Level-A rows have executing-kernel bindings. The infinite-parallel-
plate view-factor row is the explicit remaining reference-only gap; catalog
completeness does not promote it, or any other row, beyond its registry
authority.

## No-claim boundaries

**Physics scope.**

- THE CORE SOLVE IS STEADY. `ConductionState` is an ITERATION state, not a
  time state, and the assembled operator carries no time derivative. The
  transient surfaces that DO exist are the separately documented staged
  drivers above (`tests/transient.rs` theta-method marching, `tests/duty.rs`
  duty-cycle drivers, and `march_duty_cycle_with_element_materials`): they
  march the steady operator with declared capacity data and are scoped by
  their own rows — this bullet no longer denies them (that denial was a live
  README/CONTRACT contradiction, bead f85xj.16.9). No claim beyond those
  staged drivers: no general transient frontend, no adaptive time stepping,
  no k(T) marching.
- RADIATION IS A SURFACE-ONLY FIRST RUNG. Linearized surface-to-ambient exchange
  and diffuse-gray closed-enclosure exchange are implemented. This crate
  generates a closed-form concentric-sphere matrix and finite-triangle
  geometry-QMC estimates with explicit opaque first-hit visibility. Sampled
  losses cannot be redistributed into a closed matrix. Neither zero observed
  losses, empirical replicate errors nor conservation projection prove closure
  or bound geometric error. No BVH acceleration, robust-predicate guarantee,
  external-generator validation, participating medium,
  spectral/specular/semitransparent model or spatially varying emissivity is
  provided. The geometry and view factors remain fixed during a solve; no
  geometric sensitivity follows from the radiosity tangent. In particular, the
  infinite-parallel-plate constructor materializes the analytic limiting value
  `1`; it is an input convenience and admission fixture, not an executing
  view-factor kernel. The linearized departure budget is caller-declared and
  carries no accuracy guarantee; the point discrepancy is a measurement, not a
  bound. The enclosure replaces each named P1 trace by its exact area-mean
  temperature and uses that mean's fourth power, so it neither computes nor
  bounds `area_mean(T⁴) - area_mean(T)⁴`.
- THE ENCLOSURE COUPLING IS PARTITIONED. It applies one face-uniform frozen
  flux per named P1 trace in an under-relaxed outer fixed point. The original
  radiation-only entry point requires explicit adiabatic-remainder faces;
  the convective entry point explicitly composes exact uniform Robin traces
  with separate air/radiation accounting. Neither overwrites an undeclared
  physical row or adds radiation to Neumann/Dirichlet rows, nor represents a
  spatially varying trace flux, and is not a monolithic nonlinear operator.
  The report retains input emissivity uncertainty and measures the linearized
  point discrepancy, but it does not propagate a nonlinear radiosity or
  coupled-solution uncertainty certificate. No radiation adjoint is exposed.
- AMBIENT PATCH RADIATION uses a constant emissivity sampled under an explicit
  card query, and the fourth power of each trace's area mean. The inner Robin
  flux varies linearly with local P1 temperature; it is neither a uniform
  applied flux nor an integration of pointwise T^4. Only existing uniform
  Robin regions are admitted, with a black isothermal reservoir and unit view
  factor. No enclosure reflection, occlusion, spatial emissivity variation,
  radiative tangent/adjoint or continuum/uncertainty bound is inferred.
  `tests/ambient_radiation.rs` checks the independent scalar slab heat balance,
  hot-reservoir heat signs, changed emissivity, exact reconstruction, joint
  heterogeneous/contact resistance, claim pins, cancellation and both outer
  convergence gates.
- NO CONVECTION PHYSICS. The Robin row is a convective BOUNDARY COUPLING: `h`
  is an input. This crate computes no correlation, solves no boundary layer, and
  has no fluid side. Conjugate coupling is a separate bead.
- CONTACT IS MATCHING P1 ONLY. The implemented jump law requires exact
  coordinate-bit equality between duplicated boundary triangles and a positive
  `R''` that is constant ON EACH FACE. `InterfaceSurface::new` declares one
  value for the whole surface; `InterfaceSurface::new_mapped` declares one per
  face pair (bead `frankensim-extreal-program-f85xj.12.9`), which is what lets
  a measured as-built bond-line thickness or roughness map drive the operator.
  A map is still a GEOMETRY-DRIVEN SPATIAL VARIATION, not a constitutive law:
  there is still no nonmatching or mortar projection, no
  pressure/gap/temperature-dependent resistance, no changing contact area, no
  imperfect geometric match tolerance, and no automatic lowering from an
  `fs-scenario` interface object. Perfect contact is never inferred, and a
  mapped surface still resolves nothing about where its numbers came from —
  see the value-origin boundary below.
- A RESISTANCE MAP VARIES THE VALUE, NOT THE AUTHORITY. Every face of a mapped
  surface must cite the same material card; a surface whose faces disagree
  about their card refuses rather than averaging them or attributing the
  result to whichever card came first. Per-face values are bound to their own
  face pair before the deterministic face sort, so the reorder cannot separate
  a face from its resistance.
- `ResistanceValueOrigin` SEPARATES CITING A CARD FROM LAUNDERING A
  MEASUREMENT THROUGH ONE. A value produced by
  `InterfaceResistance::with_measured_value` retains the card identity and
  receipt for MATERIAL authority only and records
  `CallerSupplied { rationale }`; the card is not evidence for that number.
  This crate does not validate the supplied value, does not know how it was
  measured, and does not compose its metrology uncertainty with the card's —
  the caller supplies an uncertainty describing the supplied value, and
  `Unstated` stays an honest unknown. A measured map entering a solve does not
  make the solve validated.
- THE TRANSIENT PATH IS LINEAR, FIXED-STEP, AND ADJOINT-FREE. `transient`
  (bead `f85xj.5.11`) advances `C dT/dt + K T = b` by the theta method over
  the SAME `K` and `b` the steady path assembles, so a transient run cannot
  drift from the steady physics — it is the steady physics with one extra
  term, and `tests/transient.rs` pins that by marching to the steady solution.
  What is NOT here: a temperature-dependent `k(T)`, which is REFUSED rather
  than frozen across a step (freezing is a different scheme with its own error
  behaviour, and adopting it silently would make the observed time order
  depend on how strongly `k` varies); adaptive step control; and duty-cycle
  drivers, which remain staged as `f85xj.5.12` (E17 duty cycles) and
  `f85xj.5.13` (lumped reduced rung).
- THE TRANSIENT ADJOINT COVERS ONE PARAMETER AND ONE FUNCTIONAL FAMILY.
  `source_scale_gradient` differentiates `J = w · T^N` with respect to a
  UNIFORM SCALE on the volumetric source, and nothing else: there is no
  gradient with respect to conductivity, capacity, boundary values, per-
  component power, or step size, and no path-integrated (time-windowed)
  objective. A max-over-region functional is deliberately absent because it is
  not differentiable, and smoothing it silently would return the gradient of a
  different objective than the caller asked for.
- THE LUMPED RUNG IS CHEAPER *AND* WEAKER, AND REFUSES OUTSIDE ITS REGIME.
  `lumped` solves a first-order network ANALYTICALLY — no time stepping, so it
  contributes no integration error of its own — and is gated on the Biot
  number, the quantity that decides whether a body may be treated as
  isothermal at all. `solve_gated` REFUSES above the ceiling rather than
  returning a cheap wrong number, and the refusal points at escalating to the
  full transient rather than at loosening the gate. A caller may loosen the
  ceiling, and the verdict then reports the ceiling actually applied so a
  loosened run is never mistaken for a corpus-gated one.
- `LUMPED_BIOT_CEILING` DUPLICATES THE CORPUS ROW, ON PURPOSE. This crate must
  not depend on `fs-vvreg`, so the `thermal-a-lumped-transient` Biot ceiling is
  restated as a constant. `tests/lumped.rs` asserts the two agree, so the
  duplication cannot drift silently — the test is the mechanism that makes the
  duplication safe, not a nicety.
- NODE CONDUCTANCE EXTRACTED FROM A STEADY RISE IS MODEL-FORM EVIDENCE, NOT A
  MEASURED PROPERTY. `extract_node_from_steady_rise` computes `hA = P / dT`,
  which assumes the observed rise is dominated by the surface path the reduced
  model represents and inherits every assumption of the run it came from.
- THE REDUCED RUNG HAS ONE AMBIENT AND NO NODE-TO-NODE COUPLING. Each node has
  its own surface path to a single ambient, so the network is diagonal and its
  nodes do not exchange heat with each other. A genuine multi-node network
  with inter-node conductances is not implemented.
- A DUTY CYCLE IS A SCALE ON AN EXISTING POWER MAP, AND ITS WINDOW IS CLOSED.
  `duty::DutyCycle` declares a piecewise multiplier on the volumetric source,
  not a second description of the hardware, so the cycle cannot disagree with
  the map about which components exist. Segments TILE the window by
  construction, and asking for a scale past the window REFUSES rather than
  holding the last value — extending a finished schedule invents load history
  nobody declared. A window that is not a whole number of steps refuses for
  the same reason: a step straddling the end injects undeclared load.
- THE WINDOWED ENERGY AUDIT IS EXACT ONLY UNDER A STATED CONDITION. The theta
  method's load weighting `θ s^{n+1} + (1−θ) s^n` is the TRAPEZOID rule at
  `theta = 0.5`, which reproduces a piecewise-linear profile exactly when no
  step straddles a segment boundary — so Crank-Nicolson with aligned steps
  delivers the declared cycle energy to floating point, and
  `WindowedEnergyAudit::steps_aligned` reports whether that condition held. At
  `theta = 1` it is the right-endpoint rule and the residual is real. The
  audit reports declared, delivered, and residual WITHOUT normalizing one to
  the other: the gap is the integration error the caller chose.
- A ZERO ENERGY RESIDUAL IS NOT EVIDENCE THAT THE INTEGRATION IS EXACT. On a
  cycle whose ramps mirror each other, backward Euler over-integrates the
  rising ramp by exactly what it under-integrates the falling one, so the
  windowed residual vanishes while every individual step is still first-order
  wrong. Symmetric duty cycles are an extremely common shape, so this
  cancellation is likely rather than exotic. `tests/duty.rs` pins both sides
  of it deliberately.
- CYCLE SUMMARY QUANTITIES ARE STEP-SAMPLED. `CycleSummary::peak_temperature_k`
  is the maximum over recorded STEP ENDS; the true continuous peak can fall
  between steps and is not resolved. `time_above_limit_s` is counted in whole
  steps, so it is quantized and can be off by up to one step at each crossing.
  `excursion_k` is an amplitude, NOT a fatigue quantity: cycle counting,
  damage models, and acceleration factors are out of scope and must not be
  inferred from it. These are plain reported values; binding them to the
  evidence discipline of the thermal QoIs waits on `f85xj.5.10`.
- THE E17 SEAM IS NOT WIRED. `DutyCycle` is declared directly in this crate's
  own types. Consuming an `fs-scenario` operating-envelope or duty-cycle
  scenario object — the actual intent of `f85xj.5.12`'s (a) — waits on
  `f85xj.17.5`, which defines those objects. What exists here is the driver
  and the audit, not the scenario adapter.
- THE LUMPED BINDING IS A REGIME CHECK, NOT A TRANSIENT VALIDATION. Resolving
  the Level-A `thermal-a-lumped-transient` row shows this solver agrees with a
  closed-form APPROXIMATION inside that approximation's declared small-Biot
  context, and that the residual gap is controlled by Biot. It is not evidence
  against measured data, it says nothing about the solver outside small Biot,
  and it does not make any capability L4. A spatially resolved transient is
  now bound separately — see the slab step response below — but that is
  VERIFICATION, not validation, and neither binding is evidence against
  measured data.
- CHECKPOINTING IS VACUOUS IN THE ADMITTED REGIME, AND IS THEREFORE NOT BUILT.
  The staged plan assumed the checkpointed pattern from `fs-adjoint::timedep`,
  which exists to avoid retaining the forward trajectory. Here `A = C/dt + θK`
  and `B = C/dt − (1−θ)K` are STATE INDEPENDENT, because `k(T)` is refused, so
  the backward march never consults a forward state and there is nothing to
  checkpoint. Adding the machinery would install a mechanism that does nothing
  and imply a generality the solver does not have. Checkpointing becomes
  necessary exactly when `k(T)` is admitted, and that is the point at which it
  should be built.
- ONLY UNCONDITIONALLY STABLE SCHEMES ARE ADMITTED. `theta` is restricted to
  `[0.5, 1]`. Below one half the scheme is only conditionally stable and the
  limiting step depends on the mesh, so an explicit-leaning theta would hand
  callers a configuration that looks fine on a coarse mesh and oscillates on a
  refined one. That is a refusal, not a default.
- VOLUMETRIC HEAT CAPACITY IS DECLARED, NOT CARD-BACKED. `rho·c_p` carries no
  material card, no property-use receipt, and no provenance class — unlike
  conductivity, which does. A transient result therefore has strictly weaker
  material provenance than the steady result it settles onto, and that gap is
  real rather than cosmetic.
- POWER-MAP INGESTION IS THE CALLER'S. `PowerMap` admits already-parsed typed
  rows. Reading a component power table off disk — CSV/tabular ingestion
  through the `fs-io` quarantine and catalog path — is NOT implemented here,
  deliberately: this crate takes no `fs-io` dependency, and file parsing is a
  quarantine concern, not a physics one. The crate validates what it is handed
  and refuses what does not balance; it asserts nothing about where the
  numbers came from or whether the file they came from was well-formed.
- A SURFACE POWER MAP IS A NEUMANN REGION, AND THE REGION MUST BE THE MAP'S
  OWN. `SurfacePowerMap::neumann_flux` returns a nodal outward-flux field
  whose lumping assumes the assembled load integrates over EXACTLY the faces
  the map claims. Declaring it over any other face set delivers a power the
  audit does not report, silently, because both halves still look
  individually reasonable — so derive the region from
  `SurfacePowerMap::region_selector`. Two components may not claim the same
  face: a face carries one surface load. Adjacent footprints that share
  VERTICES are fine and still balance, because the lumped nodal area is
  map-wide rather than per-component.
- SURFACE FLUX IS SIGNED POSITIVE-LEAVING, so a surface heat SOURCE produces a
  NEGATIVE `outward_flux`. `SurfacePowerAuditRow::outward_flux_w_per_m2` is
  therefore negative for a heater, while `declared_w` and `delivered_w` are
  positive powers into the solid. The two conventions coexist deliberately;
  the sign is pinned by a solve test rather than by comment.
- A BALANCED POWER MAP IS NOT A VALIDATED ONE. The audit proves the declared,
  summed, and delivered totals agree — an internal consistency and correct
  discrete distribution claim. It says nothing about whether the declared
  watts are the component's real dissipation, whether the footprint binding
  matches the physical die, or whether the resulting temperatures are right.
  Component uncertainty is a conservative arithmetic sum of stated
  half-widths, not a probabilistic convolution, and no correlation between
  components is modelled.
- `BackwardEuler::linearize_prescribed_endpoint` binds a retained capacity-storage
  endpoint without another primal solve. It reassembles the physical free-node
  residual `C(T-T_old)+dt(K(T)T-b)` at the supplied actual Robin values and
  enforces the caller's fixed accepted target in joules and unchanged energy
  budget. Prescribed temperatures, nonlinear conductivity, regional assignment
  and matching contact remain in that physical assembly. The accepted field
  is retained bit-for-bit; no replacement state or hidden tolerance widening
  is allowed. Its zero primal iterations/inner residual mean no primal Krylov
  solve ran, not that the physical residual is exactly zero.
- `StepLinearization::pullback_with_robin_feedback` uses the same consistent
  P1 low-rank Robin factors, bounded nonsymmetric FGMRES and rounded-multiplier
  residual check as the steady mean-feedback owner, applied to `C/dt+J(T)`.
  Local coefficient/reference slopes and a row-major cross-reference matrix
  retain their exact selected-port order; diagonal cross entries add to local
  reference slopes and must not duplicate them. The complete watt-residual
  multiplier feeds existing history, P1 source and capacity-scale pullbacks
  directly, without another factor of dt. Combined convection/radiation rows
  require their physical owner to map Robin reference/log-h partials back to
  original controls. At most 64 ports and four nodal factors per port fit the
  explicit entry cap; this is not an RSS bound or derivative of a correlation,
  geometry, material law or hydraulic model. The existing radiative transient
  test checks a nonuniform trace with k(T), a prescribed-temperature jump,
  independent fourth-power scalar sensitivities, retained-field identity and
  refusal of invalid physical endpoints and undersized derivative budgets.
- `transient::enthalpy::EnthalpyBackwardEuler` evolves a spatial specific-enthalpy
  field on the existing P1 tetrahedral mesh. It reuses the existing conductivity,
  source, Neumann, Robin and finite-contact assembly; temperature-dependent k(T)
  is evaluated at each actual nonlinear trial. Fixed reference masses are
  `m_i = sum_e rho0 V_e/4`; equilibrium density changes never create or remove
  specimen mass. The nodal equilibrium chart returns temperature and liquid
  fraction, preserving a true isothermal latent plateau without apparent-Cp
  smoothing. The discrete residual is `m*(h-h_old)+dt*(A(T(h))*T(h)-b)` in joules.
  Its exact generalized Jacobian is `M+dt*J_T*diag(dT/dh)`, with column scaling;
  a zero temperature slope leaves a positive mass column. The shared bounded
  Newton/FGMRES driver uses a physical Picard-Jacobi inverse, and each current
  tangent is assembled once per Newton attempt. Sparse actions allocate no
  scratch. Typed chart, spatial, solver, work, cancellation and energy failures
  return no new state and never mutate the supplied history. The accepted
  `sum m_i*(h_new_i-h_old_i)` must match independently integrated external heat
  within the caller's absolute joule tolerance; internal contact cancels from
  that account. Spatial caps precede mass allocation and checked products bound
  Newton workspace/work arithmetic. Constitutive, assembly and operator tiles
  poll Cx, with final polling before publication.
  This initial lane binds one spatially uniform equilibrium chart and a declared
  reference density. Prescribed endpoint temperatures select enthalpy by
  inverting each node's original piecewise-linear chart. Temperatures equal
  to a latent plateau (including either endpoint) or outside source support
  refuse; a prescribed bath never chooses an interior liquid fraction.
  Fixed rows impose `m_i*(h_i-h_prescribed_i)=0`, while free rows retain the
  full nonlinear transport balance. The inward reservoir reaction includes
  both transport and `m_i*(h_i-h_old_i)/dt`, so changing boundary temperature
  charges its storage once. Fully prescribed bodies are admitted for this
  transient balance even though steady spatial solvers still require a free
  degree of freedom. `dirichlet_in_w` remains separate from volumetric and
  natural-boundary heat. Ambient/enclosure acceptance recomputes the reaction
  under the physical endpoint radiation law; `physical_dirichlet_in_w` belongs
  to that audit while the inner solution retains its applied-boundary reaction.
  The mixed-bath/radiation regression checks both against an independent scalar
  balance and the shared physical endpoint evaluator. The enthalpy adjoint
  retains prescribed constraints and their independent bath-temperature controls.
  Fixed residual rows impose `m_i*(h_i-h_prescribed_i)=0`; their history and
  physical-load derivatives are zero. Free residual rows retain their coupling
  to prescribed columns. `prescribed_temperature` returns
  `m_i*lambda_i/(dT/dh)_i` per prescribed vertex, while `nodal_load_pullback`
  returns `dt*lambda_i` on free residual rows and zero on fixed rows. P1 source
  basis columns at prescribed vertices remain active wherever they load
  neighboring free residual rows.

  `EnthalpyHistoryRun` admits unique prescribed temperatures without changing
  initial physical enthalpy. Event-aligned, adaptive and ambient-radiation
  histories charge prescribed storage jumps to the bath once. Held boundary
  knots retain the exact left limit for the interval ending at a switch; the
  new bath acts on the next finite interval. Ambiguous latent temperatures
  and out-of-domain charts refuse without advancing state. Scanning and
  state-dependent-contact histories retain their explicit Dirichlet refusal.
  The natural-boundary history runner remains unsupported.
  Geometry and energetic internal
  variables are frozen; there is no material motion, expansion, pressure work,
  remapping, phase kinetics or vaporization. The basic `advance` entry point
  has no implicit radiation law; the separate radiative entry point follows.
  `tests/enthalpy_transport.rs` compares constant-Cp stepping with the existing
  temperature solver, latent heating with independent energy formulas and the
  exact uniform lumped limit, a nonlinear-conductivity endpoint with independent
  P1 algebra, mixed sensible/latent exchange with an exact one-Newton solution,
  immutable-history/cancellation refusals, prescribed-wall P1 balances,
  boundary-jump/repeated-step heat, full-bath cooling, nonlinear conductivity
  and interior latent storage driven by a hot wall. These are numerical checks,
  not experimental material validation or interval error bounds.
  `cargo run -p fs-conduction --example enthalpy_stefan -- 40 120` evolves a
  melting front in a fixed tetrahedral slab from an initial similarity field
  and time-dependent boundary flux. The future front is never imposed. The
  example reports the computed 50% liquid contour, equivalent molten length,
  nodal temperature RMS discrepancy and independently integrated whole-run
  heat balance. Its material values are synthetic numerical-reference inputs.
  `tests/enthalpy_stefan.rs` holds 480 time steps fixed and refines from 40 to
  80 axial cells, checking improving front, molten-length and temperature
  discrepancies plus energy closure. This finite comparison does not establish
  a formal spatial or temporal order, or validate liquid-flow physics.
- `EnthalpyBackwardEuler::evaluate_endpoint` and its heterogeneous facade
  evaluate an arbitrary supplied enthalpy candidate under the current physical
  boundary law without advancing history or taking a solver iteration. They
  reuse the primal residual and source/face integration, including declared
  phase conductivity, contact and optional nonlinear ambient radiation.
  Convection and radiation remain separate signed exchanges. The returned
  nodal residual and energy defect are measurements: the coupling caller must
  enforce its own fixed acceptance gates. Equal temperatures or zero total
  energy defect do not hide opposite latent-storage defects at separate nodes.
  No tangent or differentiability claim is needed at an enthalpy-chart corner.
  `tests/enthalpy_radiation/endpoint.rs` checks those latent defects, exact
  phase-scaled P1 balances, signed reservoir heat and invalid/cancelled calls.
- `EnthalpyBackwardEuler::advance_with_ambient_radiation` and its heterogeneous
  phase-assignment facade add implicit gray ambient radiation to this same
  reference-mass step. The existing steady ambient-patch iterator is shared:
  unique uniform Robin traces use `epsilon sigma A (T_mean^4-T_ambient^4)` with
  the secant applied through the existing pointwise P1 Robin operator. Initial
  and trial area-mean surface temperatures must lie in the emissivity claim's
  support. Convection and radiation retain distinct heat reports, including
  negative radiative heat from hotter reservoirs. Every inner enthalpy solve
  starts from the same immutable history; an outer radiation iteration never
  advances physical time or adds latent heat twice. Heterogeneous conductivity,
  declared phase charts and constant finite contact remain in the inner solve.
  Raw patch-temperature and per-patch watt gates are followed by a reassembled
  endpoint residual in joules and a direct storage-minus-net-input audit using
  the nonlinear radiative heat. The residual target is the original Newton
  `max(atol,rtol*||R(h_old)||)` for the complete radiative initial residual;
  inner solves use one quarter of this target. The original absolute energy
  tolerance is never widened to accommodate the patch-watt tolerance.
  Cancellation, chart-domain failure, exhausted radiative/Newton work, or
  failed physical residual/energy closure publishes no endpoint. Spatial
  Newton iterations now also poll through the shared cancellable Krylov and
  globalization callbacks; a callback panic has a typed refusal and discards
  the incomplete attempt. `tests/enthalpy_radiation.rs` independently checks
  the mean and three spatial P1 modes on a regular tetrahedron against scalar
  backward-Euler algebra, both reservoir heat directions, unsmoothed latent
  phase exchange, deterministic replay and admission/closure refusals.
  This is a fixed-step estimated discretization, not pointwise `T(x)^4`
  quadrature, an enclosure/view-factor model, an adaptive transient frontend,
  or experimentally validated heat transfer.
- `EnthalpyBackwardEuler::advance_with_enclosure_radiation` and its heterogeneous
  facade couple the same immutable h history to the existing gray-diffuse
  radiosity owner. One to 64 nonoverlapping radiative patches must exactly match
  positive uniform Robin or uniform Neumann regions. Each patch emits at its P1
  area mean and applies a uniform outward flux. On convective traces this uses
  the algebraic reference shift `T_air - q_rad/htc`; on Neumann traces exchange
  adds to the original authored flux. An insulated trace needs no artificial
  convection, allowing finite bodies to exchange reflected radiation in vacuum
  while their independent sensible and latent enthalpies evolve. Original
  convection alone is exposed to an air consumer. `prescribed_neumann_out_w`
  independently integrates the original Neumann loads, excluding exchange;
  this distinction is essential for mixed Robin/Neumann enclosures, whose
  internal radiative heat generally does not cancel within either subset.
  Supplied view factors retain their existing closure/reciprocity admission.
  Every accepted endpoint passes fresh radiosity equations, applied patch heat,
  the complete nonlinear h residual and the original joule energy gate. The
  initial physical residual sets the Newton target, with inner solves tightened
  to one quarter of that target. Closed exchange changes local storage and its
  signed sum cancels from whole-system heat; it is not an external reservoir.
  The matching `linearize_accepted_with_enclosure_radiation` rechecks accepted h
  and adds the exact tangent `dt*G*(dq/dT_patch)*W*diag(dT/dh)` to the existing
  sparse h Jacobian. The patch derivative uses the same pivoted radiosity solve
  and checks its original equations. Two nodal vectors per patch are admitted
  before allocation; these factors participate in forward/transpose actions
  and the preconditioner. Latent zero temperature slopes remain zero, while
  direct h history, nonlinear/phase conductivity and contact remain active.
  Original convection ports are retained for inlet and HTC derivatives.
  Source/history derivatives also retain reflected exchange on vacuum traces;
  requesting a convective control port on a Neumann region refuses.
  Emissivities, view factors and geometry are fixed; emissivity validity
  endpoints refuse classical derivatives. `tests/enthalpy_enclosure.rs` adds
  independent reflected-heat/energy references, mixed latent transfer,
  source/history finite differences, transpose pairing, deterministic replay
  and admission/cancellation checks. Its vacuum tests independently solve
  finite-body equilibration, latent transfer and mixed prescribed-flux balances,
  check relabeling and source derivatives, and retain failed-step history.
  This is the existing patch-mean radiation
  discretization, without inferred geometric visibility or pointwise T^4 flux.
- `transient::enthalpy::adjoint` supplies discrete endpoint derivatives for
  this same spatial enthalpy balance. `linearize_step` solves a real endpoint;
  `linearize_accepted` rechecks a supplied endpoint's original Newton target,
  chart fields and integrated energy balance before attaching its sparse
  tangent. The linearization owns the endpoint and tangent and borrows only the
  fixed mesh, so step-local sources and boundaries need not remain alive.
  Temperature observation seeds become enthalpy seeds by `diag(dT/dh)`;
  direct enthalpy/history seeds are added before the shared bounded FGMRES
  transpose solve. For `J_h^T lambda = h_bar`, history receives
  `M_ref*lambda` and P1 nodal source density receives
  `dt*M_source^T*lambda`. The source uses the existing consistent tetrahedral
  integration, while reference storage remains lumped. This explicit `dt`
  follows the enthalpy residual's joule units. The transpose implements
  `M_ref+dt*diag(dT/dh)*J_T^T`, sharing the primal conductivity Jacobian.
  A true transpose residual, exact inner-column cap, allocation-size checks
  and cancellation gates precede gradient publication. Smooth chart knots
  and plateau interiors are admitted; temperature slope corners, chart
  validity endpoints and conductivity slope/validity boundaries refuse an
  ordinary two-sided gradient. No chart, density, individual tensor-entry or
  table-knot, geometry, boundary-control or phase-fraction derivative is inferred.
  `tests/enthalpy_adjoint.rs` checks nonuniform-slope transpose identities,
  nonlinear history/source derivatives against five-point forward differences,
  zero temperature sensitivity on a latent plateau, differentiability policy,
  and refused/cancelled/tampered inputs.
  `cargo run -p fs-ascent --example enthalpy_calibration` composes these
  derivatives through a bounded full-storage spatial trajectory to fit two
  localized heater pulses from temperature histories using the existing SQP
  optimizer. The example's synthetic observations and numerical recovery are
  separate from experimental material validation.
- `EnthalpyStepLinearization::conductivity_scale_pullback` contracts the same
  complete joule-residual multiplier into independent element conductivity
  scales. The accepted binder retains four local conductive loads per element
  from the production element temperature, assigned tensor, stiffness and
  typed phase multiplier at mean nodal liquid mass fraction. The result is
  `-dt * lambda_e^T * (K_eff,e*T_e)`, already including dt, in mesh element
  order. It differentiates a scale of the full effective tensor/curve, with
  phase-law settings and material identities fixed; the existing transpose
  already carries k(T), phase, contact, natural/air and radiation feedback.
  Summing element entries gives a regional scale derivative. No second
  state derivative, material query or solve is performed by the pullback.
  The owner retains 32 bytes per element; temporary typed phase scales during
  binding and the bounded output during pullback each use 8 bytes per element
  at separate times, plus vector headers. Allocation-size, element-cap,
  finite-full-dual and cancellation checks precede publication. The mixed-phase
  nonlinear test in `tests/enthalpy_adjoint.rs` compares this control with fresh
  perturbed solves and checks uniform-temperature, cap and refusal behavior.
- `linearize_step_with_ambient_radiation` and
  `linearize_accepted_with_ambient_radiation` bind total physical derivatives
  of that same radiative enthalpy endpoint, including through the heterogeneous
  phase-assignment facade. Original convection and reference densities remain
  fixed. The accepted binder reconstructs chart-resolved temperatures, the
  nonlinear radiative boundary and the complete original initial-residual
  target; it rechecks both nodal joule residual and independently integrated
  physical energy before publishing a linearization. Public primal telemetry
  is not evidence of acceptance. For patch surface mass `B`, mean weights `w`,
  reservoir `a` and secant `s`, the tangent includes
  `dt * (ds/dT_mean) * B(T-a) * (w^T diag(dT/dh))` in addition to the existing
  conductivity/Robin tangent. This includes conductivity's `k'(T)`, the
  correct enthalpy column scaling and both the radiative coefficient and
  reference-temperature feedback. The genuine transpose retains this rank-one
  term; neither forward iterations nor a frozen secant are differentiated.
  `EnthalpyRadiationStepGradient` returns history and P1 source-density
  pullbacks plus dt, per-patch black-reservoir temperature and absolute
  emissivity controls. The dt control holds endpoint source/boundary data
  fixed; time-dependent schedules must supply their own data derivatives.
  `into_transport` retains all radiative state feedback for existing history
  tapes. At most 64 patches and an explicit `2 * vertices * patches` scalar
  factor-entry cap bound retained feedback; no dense global tangent is formed.
  Cancellable FGMRES, true transpose-residual verification, original primal
  gates, phase/material smoothness and emissivity-support neighborhoods precede
  any gradient publication. Plateau interiors are admitted; chart corners and
  unsupported two-sided neighborhoods refuse. No phase-chart, density,
  geometry, convection or conductivity-parameter derivative is inferred.
  `tests/enthalpy_radiation.rs` compares all supported controls against actual
  perturbed nonlinear steps with mixed sensible/latent slopes and nonconstant
  conductivity for both hot and cold reservoirs. It also checks transpose
  identities, true returned dual residuals, exact latent-interior controls,
  facade/tape parity and tampered-primal/budget/cancellation refusals.
- `HeterogeneousEnthalpyBackwardEuler::linearize_accepted_storage_with_ambient_radiation`
  binds that same physical radiative endpoint to its declared reference-storage
  assignments. It returns `HeterogeneousEnthalpyStepLinearization`, retaining
  the complete sparse-plus-low-rank radiation tangent through `into_transport`
  and the existing history, P1 source and element/material reference-density
  pullbacks. The density derivative holds initial specific enthalpy, charts,
  geometry, conductivity and radiation parameters fixed; it changes only
  `rho0 * V / 4` storage masses. Unused material records receive zero. The
  old-state difference and material assignments are bound by the storage owner,
  so an unrelated tangent or history cannot be attached to these controls.
  No additional nonlinear solve is performed. The existing physical residual,
  energy, smoothness, cancellation and feedback-entry gates remain in force.
  This entry point does not return dt, reservoir or emissivity controls.
  `tests/enthalpy_radiation/adjoint.rs` checks nonuniform mixed-phase storage
  with k(T), hot/cold reservoirs, history/source parity, the genuine transpose,
  unused-record aggregation and density derivatives against fresh perturbed
  nonlinear solves.
- `HeterogeneousEnthalpyBackwardEuler::linearize_accepted_parameters_with_ambient_radiation`
  retains both reference-storage binding and physical radiation metadata in
  `HeterogeneousEnthalpyRadiationStepLinearization`. Its pullback runs exactly
  one complete radiative transpose, then moves that multiplier into the shared
  private reference-density contraction. The result contains `storage`
  (history, P1 source, element/material density and one Krylov report), plus dt,
  radiating-region names, reservoir-temperature and absolute-emissivity
  derivatives. Original convection, charts, geometry and material parameters
  remain fixed; dt holds endpoint data fixed and does not imply a schedule or
  horizon derivative. Radiative controls include the coefficient/reference
  chain rule, while density changes only the declared reference storage.
  `into_storage` drops the extra control metadata while retaining all state
  feedback. The existing fixed-radiation-parameter storage binder uses this
  same owner and remains source-compatible. The existing radiative storage
  regression checks exact parity of all combined gradients and the single
  reported transpose with the separate physical owners; independent nonlinear
  perturbations remain the controls' numerical reference.
- `EnthalpyStepLinearization::robin_response` binds original uniform convection
  ports to that same checked enthalpy endpoint. The producer retains the exact
  convection boundary; a radiative binder restores the original boundary after
  forming the complete physical tangent. Neither `transport()` nor
  `into_transport()` loses its radiative feedback. A response cannot be built
  from caller-supplied temperatures, heat reports or replacement port data.
  Its pullback accepts direct specific-enthalpy carry separately from nodal
  temperature, area-mean wall and outward convection-heat seeds. Temperature
  seeds are multiplied by `dT/dh` before adding carry: no inverse slope or
  apparent heat capacity is used on a latent plateau. The same true-residual
  FGMRES solve returns history/source gradients, convection-reference and
  `ln(HTC)` controls, and assembled nodal watt-load gradients `dt*lambda`.
  Consistent triangular Robin mass quadrature and explicit heat-functional
  terms enter the controls. All other material/chart/geometry parameters stay
  fixed. Only convection rows are exposed to the air coupling owner; radiation
  is retained inside the solid Jacobian and never counted as air heat.
- `with_robin_feedback` consumes any of the uniform, heterogeneous, radiation,
  or combined enthalpy linearization owners to attach declared convection-law
  derivatives at that accepted endpoint. A bounded matrix gives each original
  convective reference's derivative with respect to the selected area-mean
  wall temperatures; one slope per port gives its own wall-dependent HTC.
  Exact consistent P1 rank-one updates complete both tangent actions and
  Jacobi, retaining all conductivity, phase and radiation feedback. Columns
  carry `dT/dh`, including exact latent-plateau zeros. Physical convection-heat
  objective seeds also include the attached laws. History, P1 source,
  reference-density and radiation-parameter controls use the same completed
  transpose; conversion between owners preserves the feedback. Returned
  original-reference and log-HTC controls remain independent parameter
  partials at fixed wall means. The caller owns physical law evaluation and
  its derivatives; attaching a finite derivative is not validation of its
  constitutive accuracy or acceptance of a different primal endpoint. One
  attachment admits at most 64 unique original ports and an explicit total
  factor-entry cap: `4 * vertices * convection_ports` plus existing radiation
  factors. Shape, ownership, nonfinite arithmetic, work and cancellation
  refusals precede publication. The focused Robin-feedback regressions check
  full nonlinear re-solves, physical heat objectives, tangent transposes,
  latent carry, storage-density and shared radiation-reservoir derivatives.
- The Biot-gated `LumpedEnthalpyBody` admits equilibrium solid-liquid phase
  change and latent heat on a caller-supplied, bounded specific-enthalpy curve.
  It couples constant internal power with convection and surface radiation and
  returns temperature, density, and phase fraction without material-name
  branches. `LumpedThermalTransport` may retain temperature-sampled
  conductivity and emissivity receipts from the exact same material card as
  that phase curve; each accepted state evaluates those properties at its own
  temperature. Each property requires one continuously supported source claim
  throughout the requested temperature span. Every source curve knot in that
  span is retained; the typed conductivity path also retains every crossed
  surface-cell edge at its declared fixed pressure or other source context.
  Coarse caller grids cannot erase interior extrema or bridge tabulated-only
  data. The Biot gate conservatively uses the resulting
  minimum conductivity and maximum emissivity. The explicit
  `from_material_card_without_radiation` path accepts complete typed query points,
  retains conductivity source knots and receipts, and requires no emissivity
  observation. Zero radiation is a declared boundary assumption; it does not
  assert that the material has zero emissivity. Fixed source context, same-card
  identity, temperature coverage and the convection Biot gate still apply.
  Declared constants remain available
  but explicitly carry no material provenance. This is an isothermal reduced
  rung only: it refuses when the
  convection-plus-radiation Biot gate fails or when a step leaves the supplied
  phase-data domain. It does NOT solve spatial temperature gradients,
  nonequilibrium phase kinetics, thermal stress, finite-strain deformation,
  moving boundaries, free surfaces, fluid flow, ablation, or vaporization;
  those regimes require escalation to the corresponding spatial multiphysics
  solvers rather than reinterpretation of this trajectory.
  `LumpedThermalEnvironment::advance` consumes the shared
  `fs-material::phase::UniformEnthalpyStepInput` to propose one such step for
  a coupled owner. Mass, exposed area and volume come from that owner; the
  characteristic length is their `V/A`. Internal heat is converted once to
  mean power. Boundary heat is integrated from the actual endpoint fluxes,
  with a separately checked caller energy-residual budget and the complete
  transport report retained. The environment is a prescribed infinite
  uniform fluid and enclosure with independent positive absolute temperatures;
  it publishes no state and debits no finite source. Convection uses the fluid
  temperature; diffuse-gray whole-surface radiation uses the enclosure
  temperature. The conservative radiative Biot bound uses the latter, and
  the march identity includes both. Setting both temperatures equal recovers
  the shared-temperature model. Directional illumination and view-factor
  networks are outside this rung.
  `solve_lumped_contact_step` adds finite-body exchange with signed supplied heat through
  the existing `SeriesThermalResistance`. Backward Euler solves one signed
  heat transfer jointly with both supplied energies and returns both enthalpy/temperature/phase states together;
  latent heat is already included in each body's chart. It enforces both
  source domains, a contact Biot limit using each declared exchange area,
  and explicit contact-law/total-energy residual budgets. Cancellation and
  refusal publish no state. G1 sensible equilibration/refinement, G3 latent
  heat and side exchange, and G4 cancellation/domain refusals are covered by
  `tests/lumped.rs`. Resistance is held fixed for the step. Callers own
  interface-law refresh, operator composition and atomic mechanical-state
  publication. Spatial spreading, geometry evolution and dynamic contact
  pressure are not resolved by this reduced operator.
  `LumpedFiniteThermalContact::advance` connects this operator to the shared
  `UniformEnthalpyStepInput` callback. It accounts for internal specimen energy
  once in the joint contact balance, returning the actual signed boundary
  transfer and the support candidate alongside both energy residuals.
  Only accepted initial and final states must lie in-domain; an isolated
  source-only intermediate need not exist. The scheme is first-order in time;
  its residual budget does not bound temporal discretization error. Analytic
  source/contact tests cover simultaneous heating and cooling, both-body
  sources, signed cooling and side exchange.
  The caller commits the support only after its coupled specimen transaction
  accepts. The existing string-runtime integration test verifies temperature,
  damping and pressure feedback, support energy debit and late-refusal retry.
  The Biot gate alone does not prove initially uniform temperature or the
  validity of homogenizing a rapid local source. Cross-crate G1/G4 tests in
  `fs-couple/tests/string_specimen.rs` consume this production path.
- PCB HOMOGENIZATION IS CONSUMED, NOT REINTERPRETED. This crate uses the
  nominal tensor and retains its constituent property receipts. Coverage
  provenance, Reuss/Voigt structural bounds, propagated principal bounds,
  scale-separation evidence, and the explicit via-array no-claim remain on the
  caller-owned `PcbHomogenizedConductivity`. The PDE comparison runs its
  lower/nominal/upper tensors; it does not certify probability, physical board
  accuracy, trace-scale temperatures, or a nonlinear QoI uncertainty band.

**Numerical scope.**

- The production material/boundary frontend is P1 (FEEC 0-form) only, on
  body-fitted TET complexes. Level-A P2 fixtures directly exercise the
  high-order `fs-feec` simplex operator; there is no general higher-order
  thermal frontend, no hex frontend, and no CutFEM frontend.
- Flux recovery (`element_heat_flux`) is the natural element-wise-constant
  `−K∇T_h`. It is one order lower than the temperature, it is NOT
  flux-conserving, and NO superconvergence is claimed for it.
- The energy balance checks the assembled operator against an independent
  post-integration of the SAME boundary data. It catches a sign or scale error
  in one of those two paths, not an error in the data itself.
- Internal contact contributions cancel from the global energy balance. The
  report exposes each interface's integrated A-to-B heat rate, but that is not
  a local flux-conservation certificate or an interface-error estimator.
- `RobinFlux` is an ATTRIBUTION of the convective total, not a local error
  bound. Restricting the balance's own post-integration to one trace inherits
  exactly that rule's exactness — face-constant `h`, linear `T_ref` — and
  nothing stronger. In particular, `mean_htc · area · (mean_wall −
  mean_reference)` is NOT `heat_rate_w` when `h` varies over the region: the
  heat rate is the integral, the means are area-weighted diagnostics, and the
  two coincide only in the uniform-`h` case. A downstream consumer that
  reconstructs heat from the means is computing a different quantity.
- A per-surface heat rate is not a conjugate coupling. Nothing here solves a
  fluid side, raises a downstream `T_ref` by the heat a trace removed, or
  closes a solid-to-air balance; `T_ref` remains declared data. Those are
  bead `f85xj.5.7`.
- The series uncertainty budget is a conservative sum of stated symmetric
  half-widths, not a probabilistic convolution and not a confidence
  certificate. Correlation among terms is not modeled; one `Unstated` term
  makes the complete half-width unavailable. Contact-resistance adjoints are
  not exposed.

**Evidence scope — read this before quoting an order.**

- The G1 orders are OBSERVED, not proven. Each gate says: on these meshes, for
  these coefficients and data, the fitted log-log slope sits within
  `fs_mms::ORDER_GATE_TOLERANCE` of the P₁ rate. Nothing is established about
  other meshes, other coefficients, or non-smooth data.
- **P₁ on these structured Kuhn meshes reproduces exact solutions up to CUBIC
  at the nodes** (measured:
  `tests/conformance.rs::polynomial_reproduction_at_the_nodes`). A ladder run on
  a quadratic or cubic solution therefore measures the INTERPOLATION error
  `‖T − I_h T‖`, which does not depend on `K` at all. The two Dirichlet ladders
  use a QUARTIC solution for exactly this reason; the Neumann and Robin ladders
  stay on lower-order solutions because those keep the BOUNDARY DATA exactly
  representable in the P₁ trace space, which is what those cases test. Anyone
  reading those two ladders as evidence about the interior operator is reading
  more than is there.
- A residual is not an error bound. This crate reports residuals, an energy
  balance, and observed orders. It ships NO certified bound on `‖T − T_h‖` and
  no continuum goal-oriented (DWR) estimate. The `ConductionReport` deliberately
  has no field that could be mistaken for one.
  `ConductivityDesign::discrete_goal_error` evaluates the discrete identity
  `J(T_discrete) - J(T_approximate) = z^T (b - A T_approximate)` for a linear
  free-dof goal, using the same operator and Dirichlet lift as the primal.
  It recomputes and checks the dual residual, retains signed free-dof
  contributions, and refuses cancellation, nonfinite data/results or a stalled
  dual. Goal normalization protects Krylov norms under unit rescaling.
  These contributions are not elementwise DWR indicators. A mesh-error caller
  still needs an enriched space, coarse-field prolongation and consistent goal
  assembly; evaluating the solved mesh alone cannot establish mesh convergence.
  The focused `tests/adjoint.rs::discrete_goal_error_*` tests compare against an
  independent primal goal difference and check scaling, resolved/zero goals,
  invalid data, exhausted iterations and cancellation. They do not establish
  continuum effectivity, a certified bound, or an adaptive CLI solve loop.
  `adjoint::compare_discrete_goal` additionally compares two full fields under
  the production linear or nonlinear operator, including heterogeneous k(T)
  and declared matching-P1 contact. It checks the reference primal residual
  and prescribed temperatures, then reuses the PCG or actual nonsymmetric
  K'(T) transpose solve. Its signed contributions are
  `lambda_i [R(T_reference)-R(T_approximate)]_i`; its separately retained
  `linearization_remainder` is the measured goal change minus their sum.
  This includes finite nonlinear truncation and remaining dual arithmetic
  error, so a consumer must account for its magnitude. It is an observation
  between supplied discrete fields, never a bound on the continuum solution.
  A kink at the reference, extrapolation, unverified reference, changed fixed
  data, failed Krylov solve or cancellation refuses. The existing adjoint
  target checks linear equivalence/unit scaling, quadratic remainder scaling
  for a genuinely nonlinear material, and these refusal boundaries.
  `RobinGoalLinearization` extends this two-field seam to consumer-owned
  reference feedback. It retains the same checked solid residuals and actual
  material/contact tangent without another primal solve, exposes the existing
  Robin pullback, and verifies `J_s^T lambda - W^T wall_adjoint = weights`
  against the ORIGINAL goal norm. Supplied reference changes correct both
  the reference primal and the two-field residual by the exact integrated
  Robin loads. Fixed adjoints remain zero and both complete residuals must
  pass the declared tolerance. Inner pullbacks use one tenth of that tolerance;
  their total Krylov count is supplied by the consumer. The feedback data is
  a numerical input, not proof of an external derivative: `fs-airflow` owns
  the physical exponential-law implementation and its independent checks.
- The fin case is a MODEL comparison against the 1-D fin equation, not a
  discretization check. Its 2% envelope carries the fin model's own error; the
  Biot number that bounds it is computed and printed by the test.
- The cylindrical case's CONDUCTANCE envelope (0.5%) is dominated by GEOMETRY:
  the annulus is meshed as a polygon, so the solved surface sits a chord sagitta
  inside the true cylinder. Its L2 envelope (0.2% of the radial drop) is the
  discretization claim, and it is checked to shrink like `h²`.
- The spherical case avoids coordinate singularities by solving a pole-free
  patch with exactly adiabatic conical and azimuthal faces. Its heat rate is
  normalized by the patch's analytic solid angle. The 1% conductance envelope
  carries faceted-surface geometry error; the separate 0.1% L2 envelope is the
  discretization claim and must shrink like `h²`.
- The conductivity-design parameter adjoint covers the LINEAR case only, and refuses a
  temperature-dependent material rather than silently linearizing. The
  manufactured P1 dual uses `z(s) = s - s^4/4`, `-Delta z = 3s^2`, homogeneous
  Dirichlet data at `s=0`, and a natural boundary elsewhere. Its L2 ladder and
  discrete primal/dual identity establish observed consistency for that
  fixture. The separate verified design gradient establishes consistency
  between the assembled `∂R/∂T`, the analytic pullback, and the primal solve.
  Neither result establishes shape derivatives, mesh sensitivity, or a
  goal-oriented error bound.
- The P2 primal ladder covers constant isotropic conductivity, homogeneous
  Dirichlet temperature-excess data, and the existing high-order FEEC
  tetrahedral operator. It does not establish P2 temperature-dependent or
  anisotropic materials, Neumann/Robin/contact/radiation boundary lowering,
  nonlinear solves, cancellation-bounded assembly, or a retained registry
  receipt.
- The P2 adjoint ladder covers a constant-isotropic, homogeneous-Dirichlet
  manufactured dual on that same high-order operator. It establishes the
  observed L2 dual order and the discrete primal/dual identity for this
  fixture. It does not add a P2 parameter pullback or design gradient,
  nonlinear or anisotropic adjoints, shape derivatives, DWR authority,
  cancellation-bounded assembly, or a retained registry receipt.

**Deferred, and why.**

- **G5 golden replay.** This crate pins no golden hash. Every golden/replay hash
  must have a row in `golden-couplings.json` (see below), which is outside this
  bead's file boundary. The determinism evidence here is self-consistency
  (replay and snapshot-resume bitwise), which is weaker than a registered
  golden and is labelled as such.
- **Cross-ISA determinism.** No G5 audit has been run on a second ISA, so the
  class stays `DeterministicPerIsa`.
- **Nonlinear conductivity-design adjoint.** Still refused. The separate
  Robin/load response and two-field goal comparison use the actual k(T)
  Jacobian; they do not supply conductivity-parameter or shape derivatives.
- **Parallel assembly.** The assembly is single-threaded. The determinism
  argument above states what a parallel implementation would have to preserve
  (staging in canonical element order), but no parallel path exists to test.

## Golden couplings (bead y4pt)

This crate pins NO golden or replay hash, so it has no row in
`golden-couplings.json`. Adding one is the precondition for a G5 golden-replay
claim; until then the determinism evidence is the self-consistency lane
described above.

## The E17 duty-cycle seam (`duty::lower_envelope_duty_cycle`, bead f85xj.5.12)

`fs_scenario::envelope` declares operating envelopes and duty cycles;
`lower_envelope_duty_cycle` turns one into the power-versus-time schedule the
transient march consumes. The scenario object stays the single declaration of
the operating history, and this crate keeps the single declaration of the
spatial power distribution. Neither can restate the other.

### Only the power axis is applied, and the rest must not vary

A `DutyCycle` is a scale on a volumetric source, so the only envelope
coordinate it can carry is power. Ambient temperature, pressure and discrete
states act on BOUNDARY CONDITIONS, which are assembled once and are not
re-assembled per segment.

Applying the power axis and silently ignoring the others would march a
varying-ambient cycle at a fixed ambient and report a confidently wrong peak
junction temperature, with nothing downstream able to detect it. The lowering
therefore REFUSES when any non-power axis takes more than one value across the
cycle, naming the axis and both disagreeing dwells, and records the held values
in `LoweredDutyCycle::held_axes` when it does not. A caller can check those
against the boundary they built.

The power axis is resolved BEFORE the other axes are scanned. Otherwise a
misspelled axis name makes the power axis look like just another held axis, and
since power is the one axis that is SUPPOSED to vary, the caller would be told
"total-power varies across the cycle" — true, and a diagnosis of entirely the
wrong problem.

Fraction-weighted scenario cycles are refused: a transient march needs absolute
time, and dimensionless fractions carry none.

### The lowered schedule is DISCONTINUOUS, and that reverses an exactness claim

A scenario duty cycle declares dwells AT points, not transitions BETWEEN them,
so it lowers to piecewise-CONSTANT segments. No implicit smoothing: a step load
is real and common, and rounding it off changes the peak.

The consequence must be stated because it reverses a property this crate
otherwise advertises. The theta load weighting is the trapezoid rule at
`theta = 0.5`, which is exact for a piecewise-LINEAR schedule when no step
straddles a boundary — but a jump is not piecewise-linear. Across a boundary
where the scale jumps by `ds`, the endpoint weighting misses exactly
`(ds / 2) * dt` scale-seconds, first order in the step and independent of where
the jump sits.

So a scenario-declared cycle does NOT deliver its declared energy exactly, and
`WindowedEnergyAudit::residual_j` reports that gap rather than zero. The
existing exactness result was demonstrated on a CONTINUOUS fixture (constant ->
ramp -> constant -> ramp) and does not carry across the seam.

The residual is analytic, not numerical noise, and the tests pin it as such:
one asserts it equals the predicted `(ds / 2) * dt * base_power` to 1e-6 rather
than merely being small, and another asserts that halving the step halves it
EXACTLY, which is what makes it a step-resolution artifact of representing a
discontinuity rather than a defect in the lowering. A third shows the remedy —
declaring the transition as a ramp restores exact balance at the same total
energy — so the limitation comes with its own escape.

### No-claim boundaries (E17 seam)

- `held_axes` records that an axis was CONSTANT across the cycle. It does not
  check that the caller's boundary conditions correspond to those values;
  nothing at this layer can, because the boundary is built independently.

- A scale above one is admitted and means the declared envelope draws more than
  the reference map's total. The SPATIAL DISTRIBUTION is assumed unchanged
  under that scaling, which is a modelling assumption the caller owns: real
  hardware at higher load often redistributes as well as increases.

- The lowering is a pure re-expression of a declared schedule. It performs no
  solve, validates no envelope (that is
  `EnvelopeDutyCycle::validate_against`), and adds no evidence of any kind to
  the numbers it carries.

## The slab step response and the temporal order gate (bead f85xj.5.15)

Two verification items named by `f85xj.5.11` and left as prose caveats there.

### A spatially resolved analytic transient

The Level-A lumped binding is spatially TRIVIAL by construction: small Biot
means no internal gradient, so nothing in it exercises the coupling between the
spatial operator and the time integrator. The 1-D slab step response does. A
unit box with both `x` faces held at the step temperature and the other four
adiabatic is one-dimensional in `x`, and its closed form is the odd-mode series
`theta/theta_i = sum (4/(n pi)) sin(n pi x / W) exp(-n^2 pi^2 Fo)`.

Refining the MESH at a fixed Fourier number, the max nodal error against that
series fits an observed order of **1.929** against the P1 theoretical **2**,
inside `fs_mms::ORDER_GATE_TOLERANCE`. That is the first evidence that the
transient path inherits the SPATIAL order rather than only its own temporal
one.

Three conditions make the measurement mean what it says, and each is asserted
rather than assumed:

- **The series truncation is dominated, not merely small.** Forty-one odd modes
  at these Fourier numbers leave a tail bounded by ~1e-160 against measured
  errors of order 1e-1. A ladder whose reference carries a comparable error
  measures the reference.

- **The temporal error is negligible at every rung.** Doubling the step count
  on the finest-but-one mesh moves the error by about one part in 1e4. Without
  that check the fitted slope would be a blend of two discretizations and the
  reported order would mean nothing.

- **The march starts at a POSITIVE Fourier number.** The textbook step response
  begins from a uniform interior against walls already at the step temperature
  — a discontinuity at `t = 0` which is not in the P1 finite element space at
  ANY mesh size. Its interpolation error would neither vanish nor converge
  cleanly, and it would dominate exactly the quantity being measured. The march
  therefore starts from the analytic solution at `Fo = 0.02`, where the series
  has already smoothed and is consistent with the walls. It is the same step
  response; only the window excludes the instant the model cannot represent.
  This is a general point about convergence studies, not a quirk of this
  fixture: you cannot measure spatial order against an initial condition that
  is not in the finite element space.

### One order gate, not two

Temporal order now routes through `fs_mms::OrderGate`, the same gate the
spatial MMS ladders use. The ladder holds the MESH fixed and refines the step
against a much finer reference on that same mesh, so the spatial error is
identical in every run and cancels in the difference; what is fitted is purely
temporal.

The Richardson triple that previously measured this was retired. It computed
the same quantity with hand-written bands and produced no structured record, so
it amounted to a second, looser notion of "order was checked" living beside the
crate's real gate. The gate is STRICTER — `ORDER_GATE_TOLERANCE` is 0.2, giving
[0.8, 1.2] and [1.8, 2.2] against the old (0.8..1.3) and (1.7..2.3) — and emits
the same JSON-lines verdict a spatial ladder does, so a temporal regression is
now reported the same way a spatial one is. The scheme-separation assertion
moved across unchanged.

A negative control gates the backward-Euler ladder as if it were second order
and requires the refusal, because routing through a gate proves nothing unless
the gate can fail.

### No-claim boundaries (f85xj.5.15)

- This is VERIFICATION — the solver solving its own equations correctly —
  and not validation against measured data. The Level-4 story still needs the
  E04 corpus, and `f85xj.4.5`'s instrumented rig remains the blocker there.

- The observed spatial order is measured on ONE analytic transient with
  Dirichlet walls, no source, and constant properties. It is not evidence about
  order under Robin boundaries, a spatially varying source, or anisotropic
  conductivity.

- The temporal ladder is SELF-convergence against a finer step of the same
  scheme on the same mesh. It shows the scheme converges at its nominal rate;
  it does not show the semi-discrete system it converges TO is the right one.
  That is what the slab response and the steady-recovery test cover.

The existing `material::LiquidMassFractionConductivity` facade delegates its
constitutive arithmetic to `fs-material::phase::conductivity`, retaining the
original constructors, derivatives, `ConductionError` variants and messages.
`material_law()` and `From` transfer that exact declared law between the FEM and
flowing consumers; the existing FEM residual/Newton/adjoint paths are unchanged.

### Environmental response of finite-cell thermal contacts

`transient::deposition::contact::apply_material_contacts` evaluates the original
interface source at both actual endpoint enthalpy states. Scalar, curve and
sampled-surface contact laws share `ThermalInterfacePoint::is_response_axis`:
declared temperature/phase coordinates evolve, while explicitly typed static
pressure, composition, moisture, concentration or external temperature remain
at their supplied values. Every response axis must be admitted; matching stress
dimensions, untyped external axes, field/rate/history coordinates and unsupported
source triangles do not acquire a frozen physical model.

The existing source evaluator and series resistance owner retain their domains,
claim pins and receipts. Fresh bulk conductance becomes `1/(1/G_bulk + R''/A)`;
the two body temperatures stay independent and publication remains atomic.
This fixed-geometry, zero-interface-storage reduction supplies no evolving
contact pressure or composition transport. The composed current/heat regressions
in `fs-couple/tests/electrical_contacts/conditions.rs` include actual two-body
heat exchange and a source-support-hole refusal with unchanged coefficients.

### Material-dependent exchange between finite thermal reservoirs

`interface::lumped_pair::MaterialContactLaw` composes the existing ordered
`ThermalInterfacePoint` with `advance_contact_pair` and its finite-body enthalpy
solver. It selects the interface claim once, then refreshes that exact source at
both actual side temperatures and declared liquid mass fractions. An optional
mean temperature is an explicit source coordinate; other supplied conditions,
including pressure and composition, remain fixed.

Every trial starts from the original states and supplied heats. Relaxation
changes only the trial enthalpies. Acceptance requires the endpoint source,
the physical contact heat-law residual, both reservoir energy balances and the
endpoint Biot bound to pass. One signed transfer leaves A and enters B; zero
area opens the link while retaining each body's internal heating. Returned
proposals include the applied inner solve and the refreshed endpoint receipt.
Cancellation, exhausted iterations or unsupported source conditions return no
partially accepted pair.

This is a fixed-geometry, two-reservoir, zero-interface-storage reduction.
Source checks cover the visited query points, not continuous paths through
sparse multiaxis support. Arbitrary nonlinear laws have no uniqueness or
convergence guarantee. The reduced Biot check does not qualify local spreading,
flash temperatures, evolving pressure/composition or ambient/radiative exchange.
The coupled impact regression uses the actual mechanical contact area and
publishes both reservoir proposals through the existing coupled thermal owner.

## Pointwise radiation in coupled endpoint evaluation

`EnthalpyBackwardEuler` and its heterogeneous wrapper expose
`evaluate_endpoint_with_pointwise_ambient_radiation` and
`evaluate_endpoint_with_pointwise_enclosure_radiation`. They evaluate supplied
enthalpy candidates through the same original transport, radiation loads and
prescribed-temperature reactions as the implicit primal. They do not advance
history, require a constitutive tangent, or assert convergence. Latent plateaus
and their corners remain admissible physical endpoints.

The complete nodal residual includes local radiation; the heat balance retains
the selected operator's actual net outward heat, including a closed enclosure's
signed closure defect. Original Robin heat alone is exposed to convection or
coolant consumers. Each prepared radiation owner must match the exact transport
mesh and current physical boundary. Outer couplers must check both the original
nodal-residual and joule tolerances: zero total enclosure heat can conceal
opposite local defects. Existing cancellation and source-domain refusals remain
active. Focused regressions compare independent surface moments and tetrahedral
balance rows, latent-corner storage, heterogeneous zero-sum/nonzero-local
enclosure balance, and the accepted prescribed-temperature reaction.

### Finite-body contact entropy

Every accepted LumpedContactStep now carries LumpedContactEntropy from the
same two material charts and conservative heat transfer. Storage is the exact
relative integral of m*dh/T, including latent plateaus once. Contact production
uses Q*(1/T_B - 1/T_A) at the accepted backward-Euler endpoints, evaluated
without subtracting nearly equal reciprocals. It is a flux quadrature, not an
exact unresolved time integral. Supplied heat divided by each receiving body's
temperature, the temporal entropy defect, and energy-residual/T terms remain
separate. External source entropy cannot be inferred from supplied joules.
Entropy is evaluated only after acceptance, outside contact bisection. G1/G3
checks cover analytic storage, isolated balance, timestep refinement, latent
transfer, side reversal, supplied heat, equilibrium and zero duration.
