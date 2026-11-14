# CONTRACT: fs-conduction

## Centered physical heat in spatial enthalpy transport

Every original assembled system retains a temperature reference and the source,
Neumann and Robin load expressed relative to it. Physical residuals evaluate
`A*(T-T0) - b_relative`; Robin quadrature integrates `h*(T_bath-T0)` directly.
This preserves small deposited heat beside large absolute-temperature Robin
loads and removes spurious conductive heat from a uniform assembled state.
Conductivity and its derivatives still use the original absolute material
temperature. The temperature shift neither changes the enthalpy datum nor
replaces reference mass, sharp latent storage or boundary physics.

Native backward Euler and both SDIRK2 stages, their prescribed-temperature
reactions, entropy accounts and compensated interval heat reuse these same
physical rows. Sensible-temperature transients and duty cycles use them too.
The original CSR and absolute load remain available for elimination and
linearization. Assembly adds one nodal load buffer and a scalar, with no extra
material query or element pass; no measured speed or general roundoff bound
is asserted.

Focused native regressions cover exact sensible and latent equilibria on a
skew tetrahedron, a weak latent heat pulse inside a much larger Robin load,
and representable temperature translations with independently integrated
nonuniform Robin heat. These are numerical conservation checks, not material
calibration or time/spatial error certificates.

## Constant-temperature invariance of native conductive loads

Bulk P1 gradients use the three original nodal differences before multiplying
by FEEC shape gradients. Uniform temperature therefore gives exactly zero
local gradient, recovered heat flux and material-control load, including
source-bound phase and melting fields. Material tensors and their derivatives
are still evaluated at the original absolute temperature/enthalpy coordinates;
subtracting a reference does not change k(T) or its source support.

The element stiffness computes six independent contractions in the symmetric
three-edge basis and completes the constant null mode. Local matrix-free loads
and nonlinear coefficient columns use the same independent gradients and
balanced fourth nodal load. Steady, nonlinear Robin, first/second-order design
and enthalpy-stage conductivity pullbacks reuse these loads. This removes
absolute-temperature cancellation from these local operations. The centered
physical rows above also protect the assembled residual; solve and global
energy roundoff remain. No mesh-accuracy or measured-speed
claim follows from the reduced arithmetic.

Focused G0/G1/G3 tests check exact uniform-state zeros, skew affine Fourier
flux/dissipation, symmetric weak-form stiffness and representable temperature
translations. Native source-law regressions also check that shifting a
nonuniform absolute temperature still changes a temperature-dependent k(T).

## Material-phase transport on the original lumped enthalpy owner

`lumped::LumpedPhaseThermalTransport` binds the original prepared conductivity
and gray-surface law to the same complete caloric object. Conductivity, emissivity,
typed temperature/fraction axes and fixed query conditions must agree. Preparation
retains exact pinned emissivity paths and the full-chart conductivity minimum;
the existing Biot gate therefore includes the least-conductive admitted phase.

Every native march, spectral-heating residual and coupled endpoint evaluates
radiation at its actual trial enthalpy. Positive radiation feedback includes
changing emissivity on isothermal latent plateaus and uses original minimum Cp
on sensible spans. The cached interval bound combines with electrical/source
feedback in the existing monotone implicit solver. Native heat storage,
convection, radiation and fixed mass/geometry retain their existing owners.
`LumpedPhaseThermalEnvironment` exposes this same path to uniform coupled steps.

Three tests of the exact production radiation-feedback helper passed using the
pinned Rust toolchain. Native integration regressions cover latent emissivity,
source-path support and full-chart admission; their execution remains pending
the missing constellation dependencies. This bound controls implicit-root
admission, not temporal accuracy or measured-material calibration.

## Conductivity design through both physical SDIRK2 stages

`EnthalpySdirk2Linearization::conductivity_scale_pullback` retains each
accepted stage's original effective conductive element loads and returns
independent first/second-stage full-law scale derivatives plus their shared
material-control sum. Contractions use the complete eliminated-block stage
adjoints, `gamma*dt`, and each stage's own prescribed-row mask; no extra
off-diagonal tableau or state-derivative term is added. Accepted phase
multipliers and complete source-resolved directional tensors use the existing
constitutive integration and native element stiffness. Contact/radiation
feedback stays in the complete state transpose; bulk controls do not scale
storage, interfaces, boundaries or heat sources.

The pullback needs no additional solve, material query or sparse assembly.
Its two retained loads use 64 bytes per element; the bounded stage and shared
result vectors use 24 bytes per element. Element/address/allocation limits,
finite arithmetic and cancellation refuse without a partial gradient.
Explicit objective dependence remains caller-owned.

The exact production local contraction passed four native Rust checks against
the independent diffusion amplification factor, the uneliminated Butcher
stage equations, changing prescribed masks and interstage carry, and
nonfinite refusal. Full nonlinear material/phase finite-difference and latent
equilibrium integration regressions are defined but have not been executed
in this workspace's incomplete native build.

## Material-orientation controls for anisotropic conduction

`ConductivityModel::rotated` applies a declared proper world-frame rotation
to the complete retained conductivity law while preserving its source receipts,
temperature support and principal conductivity tables. The geometry, thermal
contacts and boundary laws remain independently owned.

`RobinResponse::conductivity_rotation_pullback` returns three derivatives per
element for infinitesimal right-handed world-axis rotations, in objective units
per radian. It contracts the complete thermal nodal-load adjoint with the native
stiffness assembled from `Omega*K - K*Omega`; rotation commutes with the original
P1 material integral, including crossed conductivity knots. It retains the
Dirichlet lift and adds no second K'(T) term. Grouped material orientations sum
their element vectors. Isotropic conductivity has zero rotation sensitivity.

This is a bulk material orientation control at fixed geometry, storage,
interface laws and source context. A caller's adjoint must already include all
state feedback, and explicit objective dependence needs its own direct term.
No new global matrix or primal/transpose solve is created by the contraction.
The production rotation differential passed a direct native Rust check against
independent finite rotations and isotropic/repeated-eigenvalue null directions.
Full heterogeneous nonlinear re-solve regressions are defined but require a
complete native build for execution.

## Conductivity design across fixed thermal contacts

`ConductivityDesign::with_interfaces` retains a checked, immutable matching
or nonmatching thermal-contact binding in the same operator used by primal,
adjoint, discrete-goal and retained Hessian solves. Bulk multipliers scale
only their element conductivity tensors; contact and boundary conductances
remain fixed, and the complete Dirichlet lift is retained. This design hook
does not differentiate temperature-, pressure- or design-dependent contact-law
updates. Its contact regressions use an independent series-resistance
temperature/gradient/curvature reference, interface side relabeling, and
per-element Hessian finite differences. These native integration regressions
are defined but still require execution in a complete native build.

## Material-supported and residual-checked conductivity design

`ConductivityDesign` checks every tetrahedron vertex against that element's
actual conductivity source support, including prescribed vertices and regional
material assignments. A supported element mean cannot hide an unsupported hot
or cold corner. Constant table-backed bonded films retain their independent
trace-temperature domains as well. The primal, objective, gradient, retained
Hessian goal and discrete-goal comparison all apply this same physical check.
No source values are extrapolated to rescue an otherwise converged design.

First- and second-order design paths share the existing normalized, preconditioned
PCG solve for the symmetric linear conduction operator. Both primal and dual
require finite solutions and independently recomputed residuals, poll cancellation
between iterations and obey the actual iteration allowance. A stalled or failed
adjoint returns a refusal instead of a gradient with a nonconverged report.
The first-order gradient now assembles its operator once and uses the same
conductivity pullback as the retained Hessian goal; no solver-iteration derivative
or replacement heat model is introduced.

Focused G1/G3 regression definitions in `tests/adjoint.rs` use the independently
derived single-tetrahedron balance `T_free = 300 + f/(4 k rho)`. They cover
unsupported extrema hidden by valid averages, regional overrides, bonded-film
support and signed goal scaling from `1e-200` to `1e200` under a one-iteration
budget. Pinned Rust syntax checks passed. These native integration regressions
remain unexecuted here because the full external dependency build is unavailable;
no whole-solver, continuum-error or experimental-validation claim is made.

## Mean-constrained equilibrium with live contact sources

`solve::solve_mean_temperature_newton_with_contacts` joins the existing
mean-constrained Newton solve to `TemperatureContacts`. Each residual and
line-search trial queries the original pinned contact sources at its actual
ordered trace temperatures. Matching face means, pointwise matching and retained
nonmatching quadrature keep their existing geometry and source owners. Fixed
pressure, finish and medium coordinates are retained; their evolution is not
inferred. Material contacts needing latent fractions still require original
enthalpy coordinates and refuse this temperature-only entry point.

The full tangent adds only `dK_contact/dT * T` to bulk, sheet and frozen-contact
terms. Both reduced solves retain that feedback and the reference-column
response, so enforcing the P1 volume mean does not change the solved law.
Source-domain misses permit bounded Armijo backtracking without switching
claims or extrapolating; other source failures abort. Initial and accepted heat
operators must remain connected. Every original heat equation and global heat
compatibility are checked, with no artificial reservoir reaction. The returned
`MeanTemperatureContactSolution` retains the exact accepted contact evaluation
alongside its physical conduction solution; no replacement endpoint solve runs.

Focused G1 definitions in `tests/nonmatching_contact/newton.rs` cover analytic
matching and unequal-trace equilibria, original source/pressure receipts, and
independent inadmissible-full-step/admissible-half-step source-domain trials.
These are discrete implementation checks, not experimental validation or an
inference of the mean from conserved caloric energy.

## Mean-constrained contact assemblies

`solve::solve_mean_temperature_with_interfaces` and
`solve_mean_temperature_newton_with_interfaces` use the existing sparse
matching/nonmatching contact operator for a flux-driven assembly with one
prescribed P1 volume-mean temperature. Independent traces retain their finite
temperature jump and original resistance receipts. The complete admitted heat
operator determines connectivity; geometrically separate bodies joined by a
positive contact need one mean, while thermally disconnected components refuse.
The connectivity pass reuses the assembled sparse operator and needs only one
vertex-sized disjoint-set array, without dense constraints or a second contact
assembly.

The mean is a datum, not an inferred conserved-energy temperature. Every
original heat equation is rechecked, including the temporarily eliminated
algebraic reference row. Internal contact heat is reported once and contributes
no external heat or artificial Dirichlet reaction. The Newton variant retains
the original temperature-dependent bulk/sheet laws and includes the fixed
contact in both constrained tangent solves and line-search residuals. Live
temperature-dependent contact resistance requires a separate constitutive
owner and is not inferred by these fixed-interface entry points.

Focused G1 tests in `tests/contact.rs` and `tests/nonmatching_contact.rs` check
piecewise analytic temperature fields, a manufactured nonlinear field, finite
contact jumps, original source receipts, unequal trace meshes, physical heat
compatibility and zero thermostat reaction. They establish the declared
discrete behavior, not a continuum or experimental validation claim.

## Original Cp(T) in enthalpy storage and conductivity moments

Material enthalpy storage and prescribed-temperature/equilibrium inverses use
the retained source heat-capacity integral. Source boundaries may be coarse
without replacing the physical Cp law. Nodal reference mass and the original
stored-energy and boundary-reaction accounts remain unchanged. Classical
adjoint admission compares the actual one-sided caloric slopes; a continuous
Cp source is not a kink merely because adjacent h/T chords differ.

`ExactMaterialEnthalpySource` and the `from_isotropic_card_exact` /
`from_orthotropic_card_exact` constructors prepare this same native solver
from the original Cp/density knot union, without an artificial subdivision
allowance. Reference mass, source context, directional frames and receipts
share the existing admission. The material-enthalpy regressions compare two
original boundaries with thousands of bounded chords, exercise nonlinear
heating and directional equilibration, and check the physical energy balance.

`MaterialPhaseTransport` retains the original sensible Cp endpoints beside
each prepared scalar or principal conductivity interval. P1 specific enthalpy
still defines the conductivity sampling field. For affine source Cp(T), h(T)
is quadratic and k(T) is affine on a source interval. Transforming the existing
positive simplex moment integral to temperature produces degree-six
polynomials for both mean conductivity and its nodal derivatives: four-point
Gauss integrates them exactly in real arithmetic. The derivative transform
cancels Cp analytically, retaining nonzero tangents for adjacent representable
nodal enthalpies. Only original source and nodal knots split the integral.

The original affine moment path remains in use when no nonlinear conductivity
interval is crossed. Scalar/tensor assembly, source-context binding, sparse
residuals, Newton, adaptive stepping, radiation, flux output and adjoints use
the same moment owner. Original source k(T) slopes and shared Cp continuity
determine sourced conductivity corners. There is no per-element allocation,
new source query, numerical differentiation, or extra conductivity table
approximation inside the numerical leaf.

The native P1 temperature interpolant, fixed reference geometry and mass,
explicit isotropy/orthotropy declaration, source domains and existing
semismooth/adjoint restrictions remain part of this model. Arithmetic rounding
and source uncertainty are separate from polynomial integration exactness.
Six isolated production moment tests execute successfully; whole-consumer
regressions also cover quadratic storage, Dirichlet heat, source knots, physical
flux and assembled Jacobians. Those native consumer tests were not executed in
this change because the real asupersync dependency build was killed by the
environment's memory limit.

## Original caloric extrema in lumped radiative feedback

The shared implicit-step bound for `solve_lumped_enthalpy` and prescribed
spectral heating uses the original minimum Cp on each sensible interval to
bound `dT/dh`. Affine source Cp retains its endpoint extrema even when the
public caloric chart has only two boundaries. Their average h/T chord is not
a bound on the inverse caloric derivative. Authored sensible laws keep their
constant capacity, and sharp latent plateaus contribute zero temperature slope.

The existing source-emissivity feedback bound, Biot admission, native enthalpy
balance, step/work budgets and cancellation remain the owners. The bound
requires one allocation-free pass over the original intervals; it adds no
source sampling, capacity regularization or thermal solver. It is sufficient
for the existing monotone implicit balance, not a temporal-accuracy certificate.
The focused regression defines a variable-Cp three-root counterexample, checks
the physically reachable equilibrium branch and exact storage/heat balance
through native and spectral heating, compares coarse and subdivided charts,
and compares time refinement with an independent temperature-coordinate ODE.
The new regression is source/syntax checked but unexecuted in the current
workspace; no native runtime or full-workspace verdict is claimed here.

## Original Cp continuity in contact and radiation adjoints

Contact and material pointwise-radiation derivative admission compare the
original one-sided `dT/dh = 1/Cp(T)` at a sensible caloric boundary. Different
average h/T chords do not create a derivative discontinuity when the original
Cp is continuous. The numerical tangent still uses the original pointwise
caloric derivative. Authored affine laws, true temperature/fraction corners,
latent boundaries, source knots and source-domain limits retain their existing
strict admission; no material smoothing or source extrapolation is introduced.

Focused regressions place a nonconstant-Cp material exactly at a redundant
source knot. Contact actions and transposes are compared with finite
differences of independently evaluated P1 heat transfer. A radiating spatial
endpoint with nonzero heat exchange is manufactured at that same knot; its
complete tangent and transpose are compared with an independent caloric
inverse, tetrahedron gradients and exact simplex moments. The history
pullback is also compared with a forward-perturbation reference. These regression
definitions are syntax checked but unexecuted in the current workspace;
native execution remains required and no runtime verification is claimed.

## Source-surface contact with evolving material radiation

`TemperatureContacts::advance_enthalpy_with_material_ambient_radiation` and
`advance_enthalpy_with_material_enclosure_radiation` compose the bounded
temperature-contact fixed point with the original material pointwise owners.
Triangulated temperature/pressure contact laws participate without a contact
derivative. Each native enthalpy response evaluates emissivity at the actual
nodal temperature and phase; enclosure reflection follows that same state.
The immutable old enthalpy, original caloric objects, material identities and
selected source claims remain bound throughout every discarded response.

Acceptance independently re-queries actual contact flux and material radiation
under the original residual and energy targets. Existing Newton/Krylov limits,
cancellation, and refusal behavior remain active. Ambient and open-reservoir
heat remain external; closed reflected exchange and interbody contact cancel
internally. Radiation does not become convective heat or acquire solid storage.
The accepted endpoint retains its source/phase radiation reports and remains
usable by the existing complete material-radiation/contact adjoint binder.
This path iterates contact temperature and still refuses phase-dependent
contact resistance, whose evolution requires the enthalpy-coordinate owner.

Focused regressions manufacture nonuniform latent states from independent
degree-six triangle moments and analytic two-patch radiosity, then compare the
coupled ambient, closed and open solutions with those states and their energy
accounts. They also retain source receipts, exact chart ownership, failed-step
history, repeatability and cancelled-work checks. These native regressions
require execution on a build host with the complete physics dependencies.

## Nonuniform original nonlinear film inventory

`SurfaceHeatStorage::with_inventory_profile` represents tapered or unevenly
distributed inventory while borrowing the original face-total C(T) law and
retaining its exact amounts, work constraint and source identity. The exact
P1 row moments differ from normalized corner amounts: profile `[0,1,3]` gives
fractions `[1/4,5/16,7/16]`. A zero corner still receives the positive moment
from neighboring basis functions. Scaling before summation admits finite large
and subnormal profiles without overflowing their total. Original uniform
bindings and explicit uniform profiles retain their prior arithmetic.

The same fraction multiplies original capacity, relative energy, entropy and
double-double stage storage. Existing BE/SDIRK2, adaptive, boundary-reaction,
heterogeneous and history-adjoint owners therefore consume one consistent
spatial inventory. The prepared sparse incidence encodes the original local
corner in its existing index array; no extra incidence array, table clone,
per-query allocation or source re-resolution is introduced. Nonlinear nodal
lumping, perfect bonding and the original sensible-storage limits remain.

Two native tests of the unchanged production coefficient helper passed against
independent degree-two triangle quadrature, extreme scales and invalid input.
Three integration regressions were added to the existing surface-capacity
enthalpy seam for original-source identity, nonuniform BE/SDIRK2 energy and
entropy, and manufactured transient/history-adjoint responses. Those three
integration tests remain unexecuted because the required external constellation
dependencies are absent; helper execution does not establish whole-solver proof.

## Sensible bonded-film storage in primary enthalpy

`EnthalpyBackwardEuler::with_surface_heat_capacity` and its heterogeneous counterpart bind the existing `SurfaceHeatCapacity` laws to the original host enthalpy. The shared P1 rule gives `C_i = integral(C_s N_i dA)` in J/K, including nonuniform thickness/capacity and its original material source. The native row storage is `E_i(h_i) = m_i h_i + C_i T_i(h_i)` up to a fixed datum, with increment `m_i Delta_h_i + C_i Delta_T_i`. Reference mass, host caloric charts and phase fractions keep their original owners. A sharp host latent plateau has zero sensible film increment at fixed temperature; no latent energy or apparent heat-capacity regularization is added.

`with_surface_heat_storage` adds nonlinear sensible films on the same uniform or heterogeneous host. Each `SurfaceHeatStorage::prepared` or `::tabulated` binds the original borrowed `PreparedHeatStorage` or `TabulatedHeatStorage` to one boundary triangle. Its law describes that face's **total fixed inventory** in J/K, never an areal capacity. By default each P1 row receives one third of that inventory. `with_inventory_profile([a0,a1,a2])` instead distributes the same total by the exact P1 row fractions `w_i = (2*a_i+a_j+a_k)/(4*sum(a))`, in the original boundary-face vertex order. The finite nonnegative values describe relative areal inventory, with at least one positive value; they infer no new total amount. Uniform profiles retain the original division-by-three arithmetic. At nonuniform temperatures this remains nodal lumping, not exact nonlinear surface quadrature. Reusing one law on another face declares another equal inventory. The original law integrates `C(T)` across every crossed breakpoint; it is never sampled into a constant capacity or a new independent enthalpy chart. Source-qualified amounts, constraints and evidence remain inspectable through that exact object. Declared tables acquire no material-source authority.

Backward Euler, both SDIRK2 stages, adaptive whole-window audits, retained endpoint evaluations, spectral heating, and native contact/radiation paths use this same storage. Prescribed-temperature reactions include the film's actual heat once, while prescribed constraint rows retain their original reference-mass scaling. Native endpoint/history adjoints use the independent storage slopes `m_i + [C_i + sum(w_face,i C_face(T_i))] T_i'(h_i)` at the new and old states; the eliminated SDIRK2 carry uses the first-stage slope. Active-film history temperature kinks, host-chart endpoints and active film temperature-support edges refuse classical two-sided derivatives. Interior host latent plateaus remain admissible even at a film temperature-support edge because their temperature stays fixed. Continuous capacity-table knots require no artificial corner refusal: the integrated storage derivative is the continuous original capacity. Film parameter, geometry and source-selection derivatives are not inferred.

Material-bound constant films retain their original density/Cp bundles; nonlinear films retain their original prepared or declared caloric laws. Both history and candidate temperatures must stay inside all participating source supports. Empty and zero constant additions retain no film storage or inactive sources. A face may occur once per addition; further additions explicitly represent further layers. Nonlinear face incidence is compiled once into contiguous sparse rows; caloric queries allocate nothing and touch only incident films, without copying tables or resolving claims during the solve. A separately meshed coating must not be counted again. The film is thermally thin, perfectly bonded, and purely sensible; its original capacity work constraint must be applicable to the supplied physical inventory. Independent film phase change requires another state model. Tangential film conductivity remains a separate boundary law.

G1 regressions in `tests/surface_heat_capacity/enthalpy.rs` cover constant and affine-temperature film capacities through sensible/latent/sensible host branches, independently integrated energy, BE/SDIRK2 history and source sensitivities with distinct old/stage/new capacities, finite-difference history checks, and actual prescribed-bath reactions. Constant-table equivalence exercises nonuniform heat spreading through both uniform and heterogeneous owners. Source-identity, original-domain, duplicate/unknown-face, cancellation and two-sided-derivative tests cover transactional refusals and valid latent storage at a support edge. These are numerical implementation references, not material-validation or continuum-error certificates.

## Tangential heat conduction in bonded thin coatings

`SurfaceHeatConduction` declares finite, nonnegative isotropic sheet conductance
`k*t` or a symmetric positive-semidefinite in-plane sheet tensor in W/K on
original boundary triangles. Uniform and P1 nodal coefficients
use the original face vertex order; a declared thickness or conductivity is not
inferred from the bulk material. `ThermalBoundary::with_surface_heat_conduction`
binds those laws to the original mesh and thermal partition. Perfect bonding
and negligible through-thickness gradients give the film the solid's trace
temperature, without additional temperature degrees of freedom.

The shared boundary assembly adds the exact P1 tangential Fourier block
`A * grad_s(N_i).mean(K_s).grad_s(N_j)`. Its symmetric internal transfers
preserve a constant-temperature null mode and add no reservoir load. Existing
Dirichlet reactions therefore include coating heat spreading; Neumann, Robin,
contact and radiation accounts keep their original physical meanings.
`mean_heat_flux_w_per_m` exposes the exact face-average tangential flux, in W/m
rather than the bulk flux density in W/m². Empty and zero laws leave the original
operator unchanged. Duplicate active faces refuse, including repeated additions;
coincident films require an explicit sum of their sheet conductances.

Steady, ordinary-capacity, enthalpy and their existing assembled adjoints share
this operator. Enthalpy tangents column-scale the same film block by `dT/dh`,
including zero slopes on a latent plateau. Both endpoint and event-left-limit
boundary histories preserve the declared coating. Each active triangle adds
nine COO entries between vertices already connected by its owning tetrahedron,
so the CSR graph is unchanged; staging-size arithmetic is checked. This is not
a peak-memory guarantee. The roundoff estimate counts the extra contributions.

Film storage remains an independent `SurfaceHeatCapacity` input to its existing
owner, and a separately meshed layer must not be counted again.

Tensor constructors require an explicit pair of orthonormal world-space
tangent axes, checked to absolute tolerance `1e-12`; boundary binding checks
their original face plane. Each nodal 2-D coefficient must be finite, symmetric
and positive semidefinite. Rank-one sheets retain an insulated in-plane
direction; zero sheets add nothing. P1 interpolation and exact face averaging
preserve the declared dissipation law. Tensor laws do not silently become
scalar means: the isotropic coefficient accessor returns `None` for them.
The frame is fixed on each original triangle, and no projection or reduction of
a general bulk tensor with normal/tangential coupling is inferred.

`uniform_from_bulk_zero_normal_flux` makes one bulk reduction explicit. For a
finite symmetric bulk tensor in world coordinates, positive supplied thickness,
and the declared tangent frame, it eliminates the normal gradient under the
constraint `q_n=0`: `K_s=t*(K_tt-K_tn*K_nt/K_nn)`. Normal conductivity must be
positive and the Schur complement positive semidefinite. A zero normal heat
flux can require a nonzero normal temperature gradient. The caller must justify
the thin-film trace approximation, including the resulting temperature drop;
perfect bonding and small thickness alone do not establish it. A nonzero normal
flux would additionally drive `t*K_tn*q_n/K_nn` tangential heat, which this
reduced law does not include. No externally heated or general through-thickness
anisotropic layer is inferred from this constructor.

The reduction scales the original tensor before projection, computes normal
coupling through `K_tn/sqrt(K_nn)`, and multiplies each reduced coefficient,
bulk scale and thickness with balanced magnitudes. Lost normalization terms,
unrepresentable nonzero outputs and nonpassive reductions refuse; negative
Schur entries are not clamped. The resulting law uses the same tensor assembly,
heat flux, storage composition and boundary-history paths as an explicitly
declared intrinsic sheet. G1/G3 regressions independently reconstruct 3-D
Fourier flux with its relaxed normal gradient, compare the assembled operator,
and check rigid rotations, normal reversal, isotropic/rank-one limits and
finite-output scaling.

The constant constructors add no material-source authority, temperature
dependence, finite through-thickness resistance or film-property derivative.
Continuum mean verification explicitly refuses active coatings until its flux
model represents tangential transport. Regressions in
`tests/surface_heat_capacity.rs` check exact uniform/varying coated-slab heat and
flux, the regular-tetrahedron decay mode with combined storage, enthalpy scaling,
nonlinear bulk tangents, event histories and zero/input-refusal behavior.
The `surface_heat_capacity/tensor.rs` G1/G3 regressions add rotated directional
P1 flux, dissipation and heat closure, isotropic degeneration, a backward-Euler
analytic decay mode, rank-one null directions and input admission. These are
discrete numerical references, not experimental material validation or a
continuum-error certificate. Tensor transport uses fixed-size stack arithmetic
and the same nine triplets per face; no throughput improvement is claimed.

### Temperature-dependent tangential film conductivity

`SurfaceHeatConduction::from_isotropic_conductivity` and
`from_orthotropic_conductivity` retain the original positive `ConductivityTable`
sources with a fixed positive thickness. Orthotropic principal axes remain fixed
in the face plane; no coupled bulk-tensor reduction is inferred. Source knots,
support and receipts remain accessible through `conductivity_tables`; declared
inputs gain no sourced authority. The immutable tables are shared by `Arc`,
so boundary histories clone their identity without copying the tables. The
surface carrier is now `Clone`, with borrowed read-only accessors.
`on_boundary_face` shares the same source allocation across a coating mesh;
it preserves caller-declared nodal ordering and frame, with target face/frame
admission deferred to the original boundary binder.

The original B-spline moment kernel integrates the full P1 face temperature,
including every crossed source knot, using the triangle's three-node simplex
density. Its analytic columns are `<lambda_c dk/dT>`, not `k'(mean(T))/3`.
The shared boundary owner assembles `A*g_a.<K_s>g_b` or the full tangent
`J_ac=A*g_a.<K_s>g_c+A*g_a.<lambda_c K_s'>grad_s(T)`. Coefficient derivative
columns sum to zero; internal film transport adds no reservoir power. The same
nine entries use the existing sparse graph, and the original enthalpy owner
column-scales the complete tangent by `dT/dh`. Surface storage stays separate.

Steady initialization intersects incident bulk and film temperature support.
Linear-only time/adjoint/publication paths refuse temperature-dependent films;
native nonlinear steady, backward-Euler, enthalpy and Robin response paths use
the full current law. Classical derivative admission refuses source endpoints
and a uniform true material kink. Nonuniform crossed-knot integration remains
smooth. Film thickness, source-selection and frame derivatives are not inferred.

G1/G3 definitions in `surface_heat_capacity/nonlinear.rs` compare independently
integrated directional flux/Jacobian values, a manufactured nonlinear endpoint
with ordinary and enthalpy storage, flat-table equivalence and adjoint admission.
Three std-only production moment tests cover triangle/subtriangle integrals,
coincident or adjacent temperatures and the original tetrahedral density.

## Retained pointwise-radiation and contact enthalpy endpoints

`linearize_temperature_contact_endpoint_with_radiation` binds the supplied
pointwise `EnthalpyRadiation` owner and original live contacts to a real accepted
`EnthalpyStepSolution`. Ambient, closed/open enclosure and native material-
emissivity variants delegate to their existing typed endpoint binders before
the complete contact-source feedback attaches. Mean-patch and absent-radiation
variants require their own entry points and are refused here.

The result retains the original physical target, chart/material smoothness,
workspace caps and source receipts. A fixed-emissivity owner cannot consume an
endpoint produced by state-dependent material emissivity; the material owners
retain their actual source and phase derivatives. Both contact traces, original
prescribed columns, local emission and reflected irradiation remain in one
physical transpose. Existing nodal-load and radiation parameter contractions
account for the timestep once, and the black environment remains outside solid
storage. The adapter creates no synthetic source-state receipt or new solve.

## Retained live-contact total-enthalpy endpoints

`HeterogeneousEnthalpyBackwardEuler::linearize_temperature_contact_endpoint`
consumes the actual accepted `EnthalpyStepSolution`, original previous specific
enthalpy and original absolute physical residual target. It re-queries the
retained contact laws and revalidates the unchanged endpoint through the
existing heterogeneous energy owner. The complete source-contact tangent then
attaches to that same endpoint without a new nonlinear solve or synthetic
forward receipt.

Both trace-temperature partials act through the regional chart slopes. Native
material-bound contacts additionally retain independent phase-source feedback.
Prescribed columns remain in the full tangent, and bath, source, old-enthalpy,
reference-density and phase-aware conductivity controls share its transpose.
History carry stays in specific-enthalpy coordinates across latent plateaus.
The existing natural/air Robin feedback can attach to the returned owner.
This entry point does not own radiation; accepted radiation fields require
their distinct complete boundary adapter.

Focused regressions exercise nonuniform latent histories, ordered two-sided
source laws, both prescribed-bath orientations, direct fixed-node objectives,
original residual and energy gates, source/chart corners and cancellation.
Syntax and source review do not establish numerical execution. This derivative
supplies no continuous-time or constitutive-validation certificate.

## Capacity adjoints with live contacts and ambient radiation

The complete `ContactCapacityLinearization` composes with the existing ambient
coefficient/reference feedback and remapping on the unchanged physical endpoint.
The watt-coordinate dual supplies source, whole contact-law, capacity and old-
temperature contractions without another timestep factor. Both independent
contact-temperature partials and nonlinear bulk conductivity remain in the
same transpose. The original endpoint residual and energy gates remain in force.

`tests/temperature_contact/radiation/capacity_adjoint.rs` compares six independent
controls against complete perturbed histories: whole resistance-law scale,
source scale, emissivity, reservoir temperature, total capacity and uniform
initial temperature. It covers vacuum curve/Newton and Robin source-surface
fixed-point histories, with a final and an earlier selected state. No primal
field or accepted old state is replaced. Native runtime execution remains
required; this discrete derivative supplies no temporal-error certificate.

## Complete steady contact and ambient-radiation adjoints

`AmbientRadiationLinearization` binds the existing area-mean secant at an
unchanged field and exposes the physical boundary, coefficient/reference slopes,
and weighted convection-reference feedback separately from the transpose solve.
The original ambient adjoint entries now use this same binding and remapping.

Rebinding that physical boundary through `linearize_steady_endpoint` retains
both contact-temperature partials and heterogeneous bulk K'(T). The existing
contact/Robin feedback owner supplies one complete dual; the ambient owner
then recovers original convection, absolute emissivity and independent bath
controls. Radiation-only traces have exactly zero convective derivatives.
No new primal solve, contact model, radiation solver or field mutation occurs.

Existing source smoothness, domain, cancellation, trace and work allowances
remain in force. Remapping checks nodal/trace layout, finite controls and the
dual residual; it does not independently authenticate an arbitrary supplied
dual. `tests/temperature_contact/radiation/adjoint.rs` compares four physical
controls against fresh nonlinear solves with and without convection, and checks
stale endpoints, insufficient work and malformed dual refusals. Native runtime
execution remains required. This is a local discrete steady derivative, without
a continuum, unique-maximum or gradient-error certificate.

## Prescribed temperatures through live contact laws

`ContactSteadyLinearization<'a>` retains its exact original conduction problem.
Its `prescribed_temperature_pullback_with_flux_feedback` contracts the same
complete nodal dual as the selected goal with all prescribed-temperature
columns. The result is `w_D - J_FD^T lambda_F`, including the direct fixed-node
objective, nonlinear bulk conductivity, frozen contact once, both source-law
temperature partials, and the original Robin/reference and mean-flux feedback.
Fixed vertices remain in every wall observation. No second primal or dual solve
or field mutation occurs, and another problem cannot replace the retained one.

The existing static Robin lift shares the checked contraction and keeps its
public signature. The work cap includes both bulk/contact CSR entries and the
live-law correction; it is not a peak-memory bound. Source smoothness, physical
endpoint, shape, fixed-node dual and cancellation gates remain authoritative.
`tests/temperature_contact_tangent/prescribed.rs` compares fresh bath solves on
both sides of an asymmetric temperature law with nonlinear conductivity and
Robin/flux feedback. It checks the direct prescribed-node term and explicitly
omitted contact/feedback columns. Native execution remains required; the lift
does not certify the provenance of an arbitrary caller-supplied dual.

The transient contact linearization, `ContactTemperatureLinearization`, gives
the same lift through `prescribed_temperature_lift`. The caller supplies the
problem; the temperature must equal the linearization's own field bit for bit,
and any other field refuses. As in the steady lift, the frozen contact
interfaces enter the bulk/Robin assembly, which therefore declares the
contact's coincident faces, and the retained source-law correction is added
once, at the fixed-node trace temperatures as well. A material-state contact
refuses: its transpose needs the original enthalpy coordinates. No primal or
dual solve occurs.

## Fully prescribed steady fields

The steady spatial solver and its initial-field preparation accept meshes with
no free temperature degrees of freedom. The solve still assembles and validates
the original physical operator at the prescribed temperatures, including bulk
material laws, sources and boundary/interface exchange. It returns the prescribed
field and its reservoir heat reactions with zero nonlinear/Krylov iterations.
The zero free-node residual says that no temperature equations remain; it does
not certify the continuum accuracy of a coarse, fully constrained mesh. Separate
adjoint, verification and coupled-owner admission limits remain unchanged.

## Heat-imbalance solves for linear transient conduction

The theta march solves `(C + theta*dt*K) delta_T = dt*(b - K*T_old)` with zero
increments on the already imposed Dirichlet vertices. This is the original
theta discretization, including heterogeneous capacity, conductivity and Robin
terms. Its relative Krylov tolerance is measured against the actual heat
imbalance, so an absolute-temperature storage term cannot hide a weak heat pulse.
The constant reduced operator and preconditioner are retained across steps.
Both the primal and its unchanged discrete adjoint use the shared scaled,
cancellable CG driver with a recomputed residual and at most two defect repairs
inside the original iteration budget. The reported residual concerns the
returned physical increment; it is not a spatial/time error bound or a bound
on rounding when that increment is added to an absolute temperature.

## Complete steady contact and reflected-radiation adjoints

`ContactSteadyLinearization::pullback_with_robin_and_flux_feedback_and_resistance_controls`
composes the existing complete bulk/contact tangent with the boundary owner's
mean-Robin and independent mean-flux factors. The original mean-only entry
delegates to the same implementation. No new primal, radiation solver, time
factor or frozen-contact response is introduced.

`EnclosureRadiationLinearization` supplies those factors for closed/open gray
enclosures on original uniform Robin and prescribed-flux patches. Rebinding its
fresh physical boundary to `linearize_steady_endpoint` rechecks the unchanged
field against the original residual and heat tolerances. Both contact-temperature
partials, bulk K'(T), reflected exchange and any original natural/air feedback
participate in one true transposed solve. Whole resistance-law controls contract
that same dual; the radiosity owner remaps original convection, emissivity and
the shared open-bath temperature controls.

The existing source-slope, cancellation, trace-element and feedback-memory gates
remain in force. `tests/temperature_contact/radiation/enclosure/adjoint.rs`
compares contact, source, emissivity and bath derivatives with independent
physical re-solves on mixed Robin/flux surfaces with nonlinear conductivity.
This is a local discrete steady derivative. The ambient composition above and
capacity trajectory composition below have separate owners; pointwise contact
radiation and continuum error certificates are not supplied by this entry.

## Capacity storage with open reflected radiation and live contacts

`ContactCapacityLinearization::pullback_with_robin_and_flux_feedback_and_resistance_controls`
now composes the complete stored-energy/contact response with the existing
enclosure mean-Robin and mean-flux factors. One original endpoint transpose
supplies source, history, total-capacity and whole-contact-law controls;
`EnclosureRadiationLinearization` remaps independent physical emissivities and
the single open-bath temperature. The unchanged accepted field is revalidated
against its original joule/energy targets before any derivative is used.
The retained capacity matrix includes declared volume and bonded surface
storage, without a second time factor or reconstructed volume-only operator.

The existing enclosure/contact adjoint fixture reverses three nonuniform steps
with nonlinear bulk transport, live source contacts, mixed Robin/flux patches
and nonuniform surface capacity. Independent complete trajectory perturbations
check contact, source, emissivity, bath, capacity and initial-temperature
controls. A stale endpoint or insufficient factor allowance refuses. Native
execution remains required; syntax checks do not establish a numerical pass.

The existing backward-Euler enclosure owner now accepts `OpenDiffuseEnclosure`
alongside closed exchange. `advance_prescribed_with_open_enclosure_radiation`
retains ordinary capacity responses; the matching
`TemperatureContacts::advance_capacity_with_open_enclosure_radiation` and
`advance_capacity_newton_with_open_enclosure_radiation` entries retain the
original source-surface fixed point and complete bulk/contact Newton responses.
Every trial uses the same old field and admitted step. Geometry and the prepared
radiosity operator stay fixed; no additional time integrator or material law is
introduced.

Fresh contact and radiative loads enter the complete nodal joule residual,
including prescribed-temperature storage reactions. Only physical-surface watts
enter the solid's outward radiation account; the independently evaluated black
reservoir closes the radiosity system without cancelling the solid's external
heat. Authored Neumann heat and original convection remain distinct, so only
convection can heat transported air. Applied patch loads contain only physical
surfaces, while fresh radiosity retains the separate reservoir row.

The original surface, residual, energy and work gates remain in force. All
nested responses, Newton updates, backtracks and Krylov iterations are counted.
Cancellation, source refusal or exhausted work publishes no accepted endpoint.
Focused regressions in `tests/temperature_contact/radiation.rs` use independent
two-surface gray exchange and original nodal assembly to check both heat signs,
source-curve/surface agreement, prescribed reactions, work and rollback. Native
execution remains pending; syntax checks establish no runtime pass. The complete
discrete adjoint above preserves this forward owner; no pointwise model or
temporal-error certificate is added.

## Live source contacts with open radiative cavities

`TemperatureContacts::solve_steady_newton_with_open_gray_diffuse_enclosure`
and `solve_steady_with_open_gray_diffuse_enclosure` compose the original complete
contact Newton or source-surface Picard response with `OpenDiffuseEnclosure`.
Closed and open exchange share one contact/radiosity iteration and unchanged-field
endpoint audit. The pinned source laws, heterogeneous conductivity, physical
traces, original thermal anchor and bounded response counts remain authoritative.

The black environment has no solid degree of freedom. Fresh physical-surface
loads enter free residuals and prescribed-temperature reactions; only those
surface heat rows enter the solid energy balance. The independently evaluated
reservoir row closes the radiation system and never cancels the solid's external
heat or enters transported air. Original convection stays separate. Acceptance
requires the original contact residual, heat closure, patch mismatch and radiosity
watt gates; cancellation and source-domain refusal publish no candidate. This
steady patch-mean producer adds no pointwise model or continuum-error certificate.
Capacity open-cavity storage and complete steady contact/radiosity derivatives
are composed by the shared owners described above.

## Nonlinear generated-layer temperature derivatives

`NonlinearGeneratedLayeredWallSolution::temperature_profile` borrows one
original solved planar layer. Its gradient and jet queries retain the original
Kirchhoff segment inversion, admitted conductivity source, stored face values
and cancellation checks. Successful queries allocate nothing and do not solve
another thermal problem; work remains bounded by the retained source knots.

For the original signed volumetric source, `q'=s`, `T'=-q/k(T)` and
`T''=-(s+k'(T)*(T')²)/k(T)`. The conductivity derivative contributes curvature
even without internal generation. A continuous conductivity knot retains
first-order temperature/gradient capability. At a genuine interior slope jump,
nonzero flux refuses unique second-order curvature; exact zero flux removes
the ambiguous term. Physical layer faces use the source derivative approached
from within that layer, including reversed heat flow. Profiles preserve source
and physical-layer support and never smooth contacts or substitute a fit to
rounded face temperatures. The original solver and its tolerances are unchanged.

## Original generated-layer profiles for thermal lensing

`GeneratedLayeredWallSolution::temperature_profile(layer_index)` borrows the
original layer, solved face records and declared volumetric source. It admits
constant conductivity once; `GeneratedLayerTemperatureProfile::temperature_jet_at`
then returns temperature, heat flux and both spatial derivatives in O(1) work
without allocation or source-knot rescans. The convenience solution-level
`temperature_jet_at(layer_index, offset_m)` retains the same semantics.

Interior values use the original solved inlet flux and Fourier law:
`q(x) = q_left + s*x`, `T(x) = T_left - x*(q_left + s*x/2)/k`,
`T'(x) = -q(x)/k`, and `T''(x) = -s/k`.
Stored face temperatures and fluxes are retained exactly at endpoints. Queries
refuse outside the physical layer or original conductivity temperature span,
and malformed required result rows refuse. Profiles never cross contacts,
resample a mean state, fit rounded faces or perform another thermal solve.
The renderer can therefore consume the same accepted thermal state for
refraction and focusing while retaining explicit world placement and separate
material regions. Existing conservation, source receipts, source-domain and
cancellation controls remain in the generated-conduction test suite alongside
independent profile/derivative and endpoint controls.

## Original material phase states in spatial thermal contacts

`TemperatureContact::with_material_states(MaterialContactState { point, curves })`
binds the original ordered thermal-interface source and exact native caloric
objects to matching face means, matching pointwise traces, or retained
nonmatching quadrature. Each enthalpy trial resolves the original nodal states
and interpolates temperature and liquid mass fraction independently on each
side. It never averages material enthalpies or invents a trace caloric state.
Fixed pressure, finish, medium, source selection, and geometry remain owned by
the original interface.

The existing enthalpy Newton operator now includes both
`dR/dT * dT/dh` and `dR/df * df/dh`. Phase-dependent resistance therefore
continues to evolve on a latent plateau. Signed contact exchange remains
internal and conservative. Original node states, queried fraction traces, and
source receipts remain available in the physical report.

Forward phase-bound surface laws use the original selected source-cell tangent,
including inward motion from zero or unit fraction; every residual still checks
the pinned source support. Accepted adjoints and public `linearize_enthalpy`
require unique two-sided participating chart and source derivatives. Accepted
contact and radiation adjoints rebind the exact enthalpy endpoint, so distinct
latent states with identical temperature cannot be substituted. Fixed and
temperature-only contacts also receive strict participating-vertex T(h) checks.
Temperature-only evolution refuses material-state bindings explicitly.

## Original source-bound phase conductivity on spatial enthalpy fields

`transient::enthalpy::material::phase::MaterialPhaseTransport` binds an original
`PreparedPhaseConductivity` to the existing 3D enthalpy FEM owner. The caller
declares either one isotropic source or three ordered principal-direction
sources in a fixed orthonormal frame, and supplies a positive fixed reference
density. The adapter retains every exact prepared source and their shared
original caloric chart, including source paths, selected claims and fixed
coordinates. Orthotropic sources must share that exact chart and compatible
nonthermal state; their admitted enthalpy interval is the intersection of the
three original supports. The ordered sources, frame and reference density are
part of the transport identity.

Each element uses the volume-normalized integral of `k(h_h)` over its P1
enthalpy field, with consistent nodal derivatives given by the integral of
`N_j dk/dh`. The shared scalar moment leaf splits at crossed source and nodal
knots and integrates the retained piecewise-linear law without geometric
subdivision or reinterpreting enthalpy as temperature. Residual-only evaluation
skips derivative moments; a Newton stage computes values and derivatives
together. In the orthotropic case each principal path retains its own source
knots. The physical tensor is the sum of the three integrated principal values
times their axis dyads; all four nodal derivative tensors use the same rotation.
The spatial tangent includes both `K grad(N_j) dT_j/dh_j` and
`(dK/dh_j) grad(T)`. The second term remains active during a latent plateau;
it is not multiplied by `dT/dh` a second time. A uniform isotropic source or an
entirely isotropic material table retains compact scalar buffers: one value
per element for residual evaluation, plus four nodal derivatives for Newton.
Mixed or orthotropic fields retain complete tensor values and derivatives.
Newton reuses its combined constitutive moment pass in both cases. Physical
flux recovery streams each element's tensor without retaining a whole-mesh
conductivity buffer.
Exact source endpoints and interpolation from the nearer endpoint
preserve positive conductivity across large source contrasts.

Temperature remains the native P1 interpolation of nodal chart temperatures.
Storage remains the original fixed-reference nodal mass times specific-enthalpy
change; latent heat is not counted again. Existing sparse assembly, nonlinear
acceptance, endpoint energy accounting, boundary/contact owners, adaptive
stepping and heat-source/history adjoints remain authoritative. Physical
element-average heat flux is exposed through the adapter and uses the same
original conductivity moments; it is not a face-conservative recovery.

The source requires strictly positive conductivity throughout its admitted
common support. Unsupported initial and direct endpoint enthalpies, replacement
numerical bases, element material overrides and additional phase multipliers
are refused. A finite Newton trial leaving the source enthalpy interval is a
dedicated recoverable domain refusal in the existing line search; nonfinite
values remain fatal. Classical adjoint admission conservatively rejects source
domain endpoints and a uniform element exactly at a true slope corner in any
participating principal source;
nonuniform elements may cross internal source knots. Existing caloric-chart
differentiability guards remain active.

The uniform adapter retains its original isotropic or orthotropic sources and
one caloric chart.
`HeterogeneousEnthalpyBackwardEuler::from_phase_materials` binds one explicitly
declared `MaterialPhaseTransport` per admitted material identity, retaining
each original chart, enthalpy datum, source support and fixed reference density.
Storage and source assembly share immutable element and vertex assignment maps.
The shared uniform/assigned source-moment resolver feeds the existing residual,
Newton, adaptive, contact and adjoint owners; material tables are not copied
into per-element source objects.

The physical element-average flux and the accepted conductivity-scale pullback
use the same complete tensors as the residual. A full-law scale multiplies all
principal responses together; it does not differentiate frame orientation,
source selection or individual source knots. Neither directional conductivity
nor changing anisotropy adds another mass or latent-energy account. Geometry
and principal frames remain fixed; deformation-driven rotation and an inferred
phase-mixture rule are outside this binding.

Distinct material IDs require separate traces and declared contact. Neither
enthalpy datums nor phase charts are averaged across materials. Material-local
source support and classical-adjoint smoothness remain explicit.
`phase_transport_material`, `phase_transport_problem` and
`phase_transport_heat_flux` expose the retained source, exact mesh/numerical
basis and physical element-average flux. Existing history, heat-source and
reference-density adjoints retain their original meanings; differentiating
source claims or material assignments is not established by this binding.
Tensor, sharp-melting and declared-multiplier paths remain available through
their original owners.

`from_phase_materials_with_nodal_reference_density` binds those same original
conductivity sources and caloric charts together with an explicitly authored,
positive P1 reference-density field. It assembles the immutable nodal masses
directly as `sum_e V_e*(sum_a rho_a + rho_i)/20`; no intermediate scalar-density
inventory is required. Density and enthalpy correlations therefore survive in
storage, latent evolution, endpoint energy and the existing density adjoints.
The original source still owns `k(T,f)` and its domain. The profile is a separate
reference-inventory declaration, not a density inferred from conductivity or
updated from the evolving temperature.

For this profile, `reference_material` and `phase_transport_material` return
`None` rather than expose a superseded scalar inventory. The source-only
`phase_conductivity_source` accessor retains the exact prepared conductivity
and borrowed caloric owner; `phase_transport_problem` and
`phase_transport_heat_flux` retain their original geometry and flux meanings.
Late profile replacement on an already source-bound problem remains refused.
Focused regressions use independent P1 storage and conductivity integrals,
graded sensible/latent heating, two-material conservative contact, density
JVP/adjoint pairing and finite differences, uniform scalar parity and explicit
source, budget and cancellation refusals. Geometry and reference mass remain
fixed; this constructor establishes no motion or mechanical coupling.

## Enthalpy-dependent spatial enclosure emissivity

`MaterialPointwiseEnclosureRadiation` binds original material/card/claim
emissivity to the exact caloric chart of every existing gray trace. It wraps
the original pointwise closed or open enclosure: geometry, ordered surfaces,
view factors, explicit black bath and watt gates retain their owners.
At every h candidate it integrates P1 epsilon and temperature separately,
retaining the full correlated mean of epsilon*sigma*T^4. The numerical
area-mean epsilon updates the reflection matrix; it is not a material claim.
Degree-six face quadrature is exact for the declared P1 fields.

Every candidate factors its actual radiosity system. The full h action and
transpose include the derivative of reflection as well as local emission and
material emissivity. The existing small factorization is reused per Krylov
action, with quadratic patch work and no dense nodal matrix. Original centered
radiance differences preserve small near-equilibrium exchanges. The positive
emission-only diagonal supplies preconditioning; it does not replace signed
material feedback. Actual nodal source receipts remain in patch reports.

Uniform and heterogeneous primary-h solvers, independent endpoint audits,
source contact Newton and adaptive time stepping all retain this same owner.
Adaptive comparisons still gate h, temperature and liquid fraction separately:
emissivity can evolve and redistribute latent energy while temperature stays
constant. Closed exchange remains internal; the physical surface sum of an
open enclosure supplies its signed external heat once. No phase kinetics,
moving geometry, spectral material emissivity or material-enclosure parameter
adjoint is inferred.

Focused G1/G3 fixtures use exact simplex moments, explicit two-surface
radiosity, full derivative/transpose pairing, two distinct latent materials,
energy-datum shifts and an independent continuous-time latent exchange
integral. Source/domain and cancellation refusals retain immutable history.

## Original-material ambient radiation adjoints

Explicit accepted-state material-ambient binders retain the complete h
action/transpose, including liquid-fraction emissivity feedback on isothermal
latent plateaus. Binding rechecks original source-state receipts, exact
storage/boundary ownership, the physical nodal residual and total energy.
Within explicit numerical-entry allowances, one existing transpose supplies
history/source, heterogeneous reference-density, optional contact-resistance,
original scalar, curve-ordinate or surface-sample emissivity, and independent bath-temperature
pullbacks. Shared claim controls are summed by original card/claim identity.

Geometry, caloric knots, source coordinates/cells, source selection and nonthermal
conditions stay fixed. Unsupported two-sided source/phase neighborhoods
return typed smoothness errors. Forward one-sided Newton slopes do not admit
an adjoint. Material-dependent reports cannot be rebound as frozen-emissivity
adjoints. Focused G1/G3 regressions compare analytic latent action/gradients
and independent five-point source-resolved differences; G4 covers budget,
altered-receipt and cancellation refusal.

Native compilation/execution of these spatial increments remains pending the
shared thermal dependency build. Rust parsing and independent numerical
references are not reported as native integration-test passes.


## Spatial resistance maps on independently meshed interfaces

`NonmatchingSurface::new_mapped` binds one positive, source-backed resistance
to every explicitly named A-side face. The original common-refinement geometry
is prepared once; each overlap stencil retains its owning face's coefficient.
The B-side tessellation does not redistribute or average the resistance map.
Both nonmatching overlaps and exactly matching delegated pieces use the same
assigned values for assembly, signed heat and whole-law resistance controls.

A map must cover the A-side faces exactly once and retain one original card
identity. Values and their authority stay bound through deterministic sorting.
`nonmatching_face_resistance_sources` exposes every keyed original source;
the aggregate flux's first canonical receipt is representative, not a source
claim for an averaged resistance. Fixed geometry, declared plane/coverage
tolerances and original work budgets still apply. This is spatially varying
static resistance, not a nonlinear nonmatching temperature or pressure law.
Focused G1/G3 fixtures cover a heterogeneous bondline on independent meshes,
uniform-map parity, conservation, delegation and resistance-scale contraction.

### Spatially mapped original contexts on nonlinear thermal contacts

`TemperatureContact::new_nonmatching_mapped` assigns one complete original
reference query to each explicitly named A-side face of a planar nonmatching
contact. Entries cover A exactly once, in any input order. The existing
`NonmatchingSurface::new_mapped` retains the one original card identity and
physical trace ownership; the opposite tessellation does not resample the map.
`reference_face_states` exposes every canonical original query and selected
source answer. The older singular reference accessors describe the first
canonical mapped face, not an averaged source.

Each retained overlap point updates its own ordered trace temperatures while
keeping that A face's pressure and other query coordinates and individually
pinned claim. Full temperature tangents and transposes use the same local
context and selected chart. Each distinct selected source chart is prepared
once per interface linearization, shared by all mapped faces that use it.
Frozen assembly, signed heat, contact Newton, sensible/enthalpy evolution and
endpoint adjoints continue through their original owners with no geometry
reintersection or second heat transfer.

The original pair-test allowance bounds map preparation before local queries;
Cx is polled between source admissions. Missing, duplicate, foreign, unsupported
or cancelled contexts refuse without publishing a partial law. Actual trial
source failures retain the existing transactional thermal behavior.
Pressure, geometry and all nonthermal contexts remain fixed during a thermal
solve. This supplies no mechanical pressure dynamics, gap/contact-area
evolution, source extrapolation or unsampled source-domain certificate.

Five regressions in `tests/nonmatching_contact/mapped_temperature.rs` cover
uniform-map operator/receipt parity, independent per-face heat under opposite
remeshing, distinct pinned pressure charts, fresh-query tangents/transposes and
ordered reversal, native conservative thermal evolution, and admission/support
refusals.

## Steady source-backed thermal contact

`TemperatureContact::new_pinned` admits a face-mean temperature law using an
already selected claim. Native project bindings can retain their exact source
selection even when the card contains competing claims. Missing pins, source
domain exits, incompatible units and nonpositive answers retain the original
source refusals; later queries cannot switch claims.

`TemperatureContacts::solve_steady` couples the original steady spatial solver
to the same pinned matching or nonmatching contact laws. The initial full
Kelvin field is explicit; prescribed temperatures replace its fixed entries.
Every counted inner solve uses the relaxed query field as a warm start and
retains heterogeneous and nonlinear bulk conductivity through the original
conduction owner. Constitutive trials do not create physical time steps.

Acceptance freshly queries the source at the unchanged candidate, assembles
its complete physical operator and checks the Euclidean free-node residual
in watts together with the independent whole-domain heat balance. Opposing
local contact defects cannot pass through cancellation of total heat. The
returned physical source states, free-node residuals, reservoir reactions and
Robin reports describe that accepted field; the retained inner report remains
explicitly tied to its applied contact coefficients. No uncounted replacement
solve, tolerance widening or caller-state mutation occurs on refusal.

`SteadyContactConfig` declares contact iterations, relaxation, a watt residual
target and a separate watt heat-closure tolerance. Counted nonlinear/Krylov
work is bounded by the existing inner solver limits and the contact allowance.
This is a bounded fixed point, with no claim of arbitrary-law convergence,
or uniqueness. Pure-Neumann and fully prescribed limitations remain those
of the existing spatial solver. Analytical slab
regressions cover matching/nonmatching traces, Dirichlet/Robin reservoirs,
source support and a frozen-contact candidate with cancelling global heat
but an unacceptable local residual.

## Coupled Newton solve for strong thermal-contact feedback

`TemperatureContacts::solve_steady_newton` solves the full steady physical
residual through the original `fs_solver::NewtonKrylovState`. The original
bulk-plus-frozen-contact Jacobian receives only the missing contact-temperature
feedback; its sparse assembly and physical preconditioner are reused throughout
each Krylov solve. Both matching face means and nonmatching overlap quadrature
retain their independent side temperatures and exact pinned source context.

The existing Armijo driver reduces trials that leave the positive Kelvin or
admitted bulk/source temperature domains. A crate-private callback preserves
the original typed material query error before its public diagnostic wrapper;
only the specified domain errors are retryable. No source extrapolation, claim
switch, finite-difference Jacobian, dense solve or repeated intersection is
introduced. Structural/source-arithmetic errors and cancellation abort.

The numerical absolute tolerance cannot exceed the separate fixed physical
residual target; relative stopping is disabled for this entry point. Iteration,
physical residual evaluation, Krylov and contact work have explicit limits.
Acceptance reuses the original physical evaluation only when its free field
matches the accepted Newton state bit for bit, then checks the free-node watt
residual and independent reservoir heat closure. The result retains actual
contact receipts, Dirichlet/Robin heat and the unmodified Newton report. It
performs no hidden final PDE solve and creates no frozen-solver report.
A run that stops because its final line search is rejected reaches those
gates only when the smallest trial, which the line search evaluates last,
equals the accepted field bit for bit: the Newton step then moves no free
temperature, and a numerical target below the residual's rounding floor
cannot be met at that field. Its report keeps `converged == false` and
`LineSearchRejected`. Every other stall refuses.

The stiff slab regressions have the independent admitted root of 100 W and
check the actual first Newton direction, source-domain damping, local residuals
and both reservoir reactions on matching and nonmatching meshes. A frozen
Picard step leaves source support for the same fixture. Continuous source
curves use the original generalized primal tangent at knots and inward source
endpoints; surfaces retain their unique-gradient requirement. Regressions
start at both endpoints and a genuine slope corner while independently
checking that strict sensitivity queries still refuse there. Arbitrary-law
convergence, root uniqueness, evolving pressure/contact area and interface
storage are outside this stationary contact model.

`TemperatureContacts::solve_steady_newton_with_ambient_radiation` composes
that complete bulk/contact Newton solve with the original area-mean gray
radiation driver. Each radiative response reuses its last accepted temperature
as an initial field. Source support, prescribed values and all contact queries
remain active. The original initial-field policy is also available without
claiming a conduction-only boundary anchor, so a declared radiating zero-flux
surface can supply the real thermal reservoir.

After the original radiation gates pass, the unchanged endpoint is independently
reassembled with newly queried contact resistance and its exact radiative
coefficient. Free-node residual, whole-domain heat closure and the separate
convection/radiation decomposition must meet the original watt allowances.
Original convection alone is available to transported-air consumers. The result
retains both the applied secant boundary and the physical endpoint boundary,
as well as cumulative Newton/Krylov work. Its last Newton history describes
the applied secant; this is a bounded radiation fixed point and does not supply
a monolithic radiation tangent or a contact/radiation adjoint.

## Coupled transient capacity contact feedback

`TemperatureContacts::advance_capacity_newton` solves constant-capacity
backward Euler with the source-backed contact response inside the same Newton
system as nonlinear bulk conductivity. Every residual trial re-queries the
pinned contact law at its actual trace temperatures. The tangent adds only the
missing contact-temperature feedback to the ordinary bulk-plus-frozen-contact
Jacobian, so the conservative contact block is counted once. Constant bulk
conductivity uses this same owner when the contact response is nonlinear.
Continuous source curves use their generalized primal tangent at knots and
inward source endpoints; strict endpoint sensitivities keep their independent
two-sided derivative requirement.

The existing Armijo globalization backtracks positive-Kelvin, bulk-material and
typed contact-source domain misses without switching claims or extrapolating
source data. The original frozen-contact operator supplies the physical
preconditioner. A caller-declared face cap bounds contact tangent assembly, and
the original linear/Newton work budgets remain authoritative. Focused capacity
regressions advance inward from both source endpoints and a slope corner,
recheck the physical endpoint and refuse outward source-domain trials without
changing old history.

Publication requires both the nonlinear solver's endpoint gate and a fresh
complete contact endpoint audit against the caller's absolute joule residual and
energy tolerances. The returned contact receipts therefore describe the exact
published temperature field; no outer Picard replacement solve is hidden.
A Newton run whose final line search is rejected skips the solver's numerical
residual gate only when the smallest trial, which the line search evaluates
last, equals the field bit for bit. The Newton step then moves no free
temperature, and a numerical target below the residual's rounding floor cannot
be met at that field. That field still meets the step's energy gate and the
fresh contact audit above, and its report keeps `residual_j` above
`threshold_j`. Every other stall refuses.
Geometry, pressure, contact area and all nonthermal source coordinates remain
fixed. This is backward Euler with constant prepared heat capacity; variable
caloric storage and latent heat remain owned by the primary-enthalpy path.

`TemperatureContacts::advance_capacity_with_ambient_radiation` and
`advance_capacity_newton_with_ambient_radiation` compose those same capacity
responses with the existing area-mean ambient gray-patch iterator. Source
surfaces retain the bounded contact fixed point; continuous curves can use the
complete bulk/contact Newton tangent. Each contact and radiation trial receives
the original full old temperature field. Radiation changes only its applied
Robin secant and never advances capacity history or reapplies source energy.

The original patch-temperature and watt gates remain active. Acceptance also
re-queries the pinned contact sources and evaluates the complete endpoint with
its actual radiation coefficient. Both the free-node joule residual and the
separate physical and reported storage/heat accounts must pass the original
absolute contact and energy gates. Inner contact and numerical tolerances only
tighten. Prescribed reactions include fixed-node storage once; internal contact
heat is not an external input. Radiation-only zero-flux traces use their actual
combined-boundary area and temperature without inventing convection rows.

`ContactCapacityRadiationStep` keeps the final applied-secant numerical response
separate from its fresh physical contact states and `EndpointEvaluation`.
`radiation.iterations` counts radiation callbacks; `thermal_responses` counts
all nested contact responses. Newton updates, backtracks and Krylov work include
discarded trials and are bounded by checked products of the existing controls.
This is a bounded nested iteration, not a monolithic radiation tangent or a
contact/radiation derivative. Enclosure and pointwise exchange are outside this
capacity entry point.

Focused G1/G4 cases in `tests/temperature_contact/radiation.rs` compare curve
Newton with original source-surface fixed points, independently reconstruct
P1 storage/residual and reservoir heat, cover radiation-only and convective
traces, charge both prescribed-bath jump directions, and check bounded refusal,
cancellation and exact retry. Native execution of these added cases remains
pending the authoring environment's unchanged dependency/tooling blocker;
syntax validation does not establish numerical execution.

`TemperatureContacts::advance_capacity_with_enclosure_radiation` and
`advance_capacity_newton_with_enclosure_radiation` extend the same live
capacity/contact responses to the existing closed area-mean gray enclosure.
They supply complete response and fresh-contact audit callbacks to the original
capacity enclosure iterator. The reflected radiosity solve, exact trace binding,
shifted boundary loads and physical heat decomposition retain one shared owner.
Source surfaces use bounded contact fixed points; continuous curves use the
complete contact/bulk Newton tangent. All trials receive the unchanged original
old temperatures, including prescribed nodes.

The original patch-temperature and applied-versus-fresh heat gates remain
active. Acceptance additionally requires the freshly re-queried contact and
reflected-radiation endpoint to meet the original absolute joule contact gate
and independent storage/heat gate. Inner contact and numerical tolerances only
tighten. Authored Neumann input and original convection remain separate from
signed internal enclosure exchange; prescribed bath reactions charge fixed-node
storage exactly once. `ContactCapacityEnclosureStep` retains the last applied
numerical response separately from the fresh physical endpoint and contact
states. Enclosure callbacks and all thermal responses are counted separately;
Newton updates, backtracks and Krylov work cover every response under checked
combined caps. The ordinary enclosure entry point still preserves its original
predictor-based target, including the all-prescribed zero free-residual target.

Focused G1/G4 closed-contact cases in `tests/temperature_contact/radiation.rs`
compare curve Newton and source-surface fixed points, reconstruct nodal storage
and residuals, check independent two-gray-patch heat and external balances,
exercise both bath-jump directions, and preserve history through refusal and
retry. These new cases have syntax validation only under the same native-build
blocker. Pointwise/open exchange and contact/radiation sensitivities remain
outside these capacity entry points.

## Material-state spatial ambient radiation

`MaterialPointwiseAmbientRadiation` resolves original surface emissivity at each
participating node's enthalpy-derived temperature and liquid mass fraction. Each
named trace retains its original material card, caloric object, fixed nonthermal
conditions and pinned source claim. Responses include static scalars,
piecewise-linear curves in the declared temperature or liquid-fraction axis,
and original triangulated `Surface2D(T, liquid fraction)` sources in either axis
order. Surface tangents and original-sample control weights reuse matdb's
prepared source response; the original cell topology is never replaced by a
grid. A different material chart requires a distinct trace; conflicting shared
vertex charts and mismatched storage/mesh bindings refuse before nodal lookup.

The declared P1 discretization interpolates both temperature and emissivity.
Positive 16-point Duffy/Gauss quadrature integrates the resulting degree-six
surface load and tangent polynomials up to floating-point roundoff. Prepared
face/node bindings and one source query per participating patch node avoid
repeated material queries at quadrature points. Factored fourth-power differences
retain the sign and precision of emission or absorption near a reservoir.

The complete enthalpy tangent includes both emission times dT/dh and the direct
material term (T^4 - T_bath^4) times d epsilon/dh. Liquid-fraction feedback survives
isothermal latent plateaus. The original positive frozen-emission diagonal remains
the preconditioner; signed material feedback stays in the true matrix-free action.
The complete h action and transpose are available for physical operator checks.
Reports retain actual nodal phase, emissivity and material receipts; fixed-snapshot
adjoints explicitly refuse those reports before reconstructing an endpoint.

For a joint source the chain rule is `d epsilon/dh = epsilon_T*dT/dh +
epsilon_f*df/dh`. A strict derivative requires source support only along locally
varying chart coordinates: a held-temperature latent seam or held-fraction
solid/liquid support edge can remain admissible. Changing-coordinate corners,
support endpoints, chart corners and active applicability endpoints refuse.
Primal trials use the same active-axis gradient wherever strict admission is
possible and otherwise retain an explicit one-sided source-cell slope.

Ambient and closed/open enclosure adjoints retain original surface axes, sample
coordinates, triangle order and sample-value derivatives. Each node contributes
to at most three source samples; enclosure controls include the current
emissivity-dependent reflection solve. Numerical control workspace charges
eight words per surface sample and eight per source triangle, covering retained
metadata and contraction output, while the existing 128-word face allowance
covers the enlarged basis. Immutable prepared source caches remain owner memory,
separate from this adjoint allowance. Focused source and thermal fixtures cover
canonical basis weights, latent energy/action/control oracles, both source axis
orders and simultaneous temperature/fraction tangents; native execution remains
pending an available Rust build environment.

Uniform, heterogeneous, generic `EnthalpyRadiation::MaterialAmbient`, independent
endpoint-audit and coupled temperature-contact methods use the same original
enthalpy/Newton owner. Adaptive material-ambient methods, with and without
source-backed contacts, reuse the original h/T/fraction controller and accepted-only
heat integration. Source-domain trials backtrack without clipping, extrapolation
or source switching. Tests cover hot/cold latent evolution, enthalpy datum shifts,
nonuniform exact simplex moments, full tangent/transpose checks, source/mesh refusal
and phase-controlled trajectories against an independent exponential solution.
These source laws also feed the existing material enclosure reflection and
original-source adjoint owners; spectral or directional emissivity is not added.

### Original material emissivity on native finite-body surfaces

`MaterialGraySurfaceLaw` exposes the same pinned `MaterialEmissivitySource`
through `fs-thermal`'s existing `ImplicitGraySurfaceLaw`. Native gray ambient
boundaries and two-finite-body radiative links can therefore use original scalar,
temperature, phase-fraction and joint source-surface claims. Every trial queries
the actual bulk enthalpy and retains the analytic chain rule, including changing
emissivity on an isothermal latent plateau. No second interpolation or caloric
chart is created. The native conductivity leg and thermal store must retain the
same original caloric object; matching card contents alone are insufficient.

Construction checks fixed caloric-source conditions, pins the emissivity claim
and verifies both endpoints of explicitly declared finite enthalpy bounds. Each
runtime sample still checks complete source support. `evaluate` returns the
original phase, emissivity and resolved material receipts with the native value
and tangent; `reference_emissivity` is only the preparation snapshot. The source
owner's original emissivity range `(0,1]` remains unchanged. Convection is zero,
leaving nonradiative gap or boundary exchange to its separately declared owner.
This is an explicit diffuse-gray bulk-state closure: a total hemispherical
measurement does not itself establish spectral or directional independence.
The existing finite-face owner alone computes heat, massless surface balance,
enthalpy derivatives and entropy exchange.

## Coupled primary-enthalpy contact Newton

`TemperatureContacts::advance_enthalpy_newton` closes source-backed R(T)
inside the original spatial enthalpy residual. Reference masses, distinct
material charts, nonlinear bulk conductivity, phase-conductivity response and
sharp latent plateaus remain owned by the original enthalpy step. Matching
and nonmatching contacts retain their existing physical quadrature and pinned
source claims. The exact additional tangent is
`dt * (dK_contact/dT * T) * dT/dh`; the original tangent already contains
frozen contact conduction. Latent storage remains `m * delta_h`, including
where `dT/dh` is zero. There is no apparent-capacity approximation.

The original fs-solver Armijo/FGMRES driver bounds Newton and linear work.
One complete sparse tangent is reused at each Newton field, with the original
positive Picard preconditioner. An explicit residual-callback budget includes
initialization and rejected trials. Only original chart, bulk-temperature or
typed contact-source domain misses are retryable. One-dimensional resistance
curves admit generalized Newton slopes at their knots and inward source ends;
ambiguous surface gradients and other constitutive failures still refuse. Absolute numerical
stopping is no looser than the independent physical nodal joule gate, and the
unchanged energy-balance gate checks the accepted state before publication.
The result's applied and physical contacts describe that same endpoint.

G1 tests compare spatial melting with independent eight-node enthalpy roots
and a stiff R(T) case with an independent 60-digit scalar/matrix root inside
source support that the frozen-contact response leaves. G4 tests exercise
bounded residual work, chart-domain refusal and cancellation without changing
history. The original entry point is fixed-geometry, fixed-pressure,
non-radiative thermal contact; it asserts neither arbitrary-law convergence
nor uniqueness or experimental
validation of the synthetic constitutive fixtures.

## Adaptive coupled enthalpy transport

`AdaptiveEnthalpyState::advance_to_with_temperature_contacts` and its pointwise
ambient/enclosure variants reuse the original full-step/two-half-step controller.
Every comparison and fine response solves the actual source-backed contact and
selected radiation law with the adaptive state's unchanged Newton/energy policy.
The three local error gates remain specific enthalpy, temperature and liquid mass
fraction, with reference-invariant relative scales. Only complete accepted fine
pairs change history or contribute external heat; source support, binding and
cancellation refusals preserve the last accepted field. Observed methods expose
only accepted fine pairs and retain exact event clipping and stop/resume behavior.

Per-endpoint residual and contact work bounds combine with checked three-endpoint
trial and linear-work products; reports count actual attempted endpoints, including
failed and rejected comparisons. Finite independent local residual failures use
`EnthalpyError::EndpointResidual`, distinct from Newton stalls and whole-domain
energy defects, and trigger the original bounded timestep reduction. Tests cover
an independent time-continuous P1 melting reference, accepted-only radiative heat,
closed internal exchange, bounded failed work and transactional numerical retries.

## Coupled enthalpy contact and pointwise radiation

`TemperatureContacts::advance_enthalpy_newton_with_pointwise_ambient_radiation`
and `advance_enthalpy_newton_with_pointwise_enclosure_radiation` compose the
original local-radiation owners with source-backed contact inside one primary-
enthalpy Newton solve. Each residual uses its actual chart-resolved P1
surface temperatures. The complete Jacobian adds the original analytic
radiation action, column-scaled by `dT/dh` and multiplied by `dt`, to the full
contact/transport tangent. Reflected patch coupling remains matrix-free;
the original positive local-emission diagonal remains the preconditioner.
There is no outer secant-radiation response or patch-mean replacement for
local emission.

The prepared radiation owner must retain the exact mesh and current boundary.
Fixed emissivity snapshots retain their original nodal temperature support;
enclosure view factors and uniform patch irradiation retain the owner's
existing model. Contact and closed-enclosure exchange are internal transfers;
ambient/open-enclosure heat, convection, prescribed flux and fixed-temperature
reactions keep distinct signed accounts. The unchanged endpoint is independently
reassembled under both physical laws before local residual and energy acceptance.
Physical pointwise reports remain in `ContactEnthalpyStep.conduction`; legacy
outer-radiation counters and secant-boundary fields remain empty.

The explicit `linearize_accepted_enthalpy_with_pointwise_ambient_radiation` and
enclosure counterpart rebind the original accepted state, source contacts and
radiation model before retaining one complete physical Jacobian and transpose.
`pullback_with_pointwise_radiation_controls` uses one transpose solve for
history, heat-source, reference-density, whole-resistance-law, emissivity and
reservoir-temperature controls. Outward radiation contracts the negative
incoming-load dual. Frozen contact conduction and contact-temperature feedback
are each included once; no derivative of solver iterations is substituted.
Geometry, view factors and source emissivity snapshots stay fixed. Accepted
adjoints retain strict source/chart smoothness admission, including refusal
of ambiguous source corners. The nonradiative binder rejects radiative results.

G1 references manufacture immutable histories from exact triangle polynomial
moments and analytic two-patch radiosity, then recover nonuniform latent
endpoints with unequal reference masses, heating/cooling baths, independent
material enthalpy offsets and separate Robin convection. A closed-enclosure
case checks local exchange despite zero total radiative heat. Complete
ambient/enclosure adjoints are compared with independently perturbed primal
solves for resistance, emissivity, bath, history and source controls. G4 cases
check exact owner binding, unsupported local source temperatures and
cancellation without advancing history. These are numerical-model checks,
not geometry-derived view-factor or experimental material validation.

## Primal contact tangents at tabulated source knots

`TemperatureContacts::linearize_temperature_for_newton` selects a generalized
tangent of the original continuous piecewise-linear resistance law: the right
segment at an interior knot and the inward segment at either source endpoint.
Scalar laws retain zero temperature derivatives. The complete contact flux
still includes the original frozen operator and its resistance feedback once.
The steady, constant-capacity and primary-enthalpy Newton owners use this tangent
without smoothing the source curve or replacing latent storage. Every physical
trial retains its pinned source-domain check; unsupported steps are not extrapolated.

The existing `linearize_temperature` and all accepted-endpoint adjoints retain
their unique two-sided derivative requirements. The new primal rule does not
assert a unique physical sensitivity at a constitutive corner. Two-dimensional
source surfaces retain their existing strict gradient admission.

The focused G1/G3 regressions use the exact consistent-triangle contact flux at
an equilibrium resistance kink, both perturbation directions, known source
segment slopes, and analytic uniform enthalpy heating/cooling from knots and
source endpoints. Outward source-domain steps preserve the old history on refusal.

## Complete steady contact tangents and adjoints

`TemperatureContacts::linearize_steady_endpoint` binds the unchanged physical
temperature field after rechecking exact prescribed values, the original
free-node residual in watts and whole-domain heat closure. Insulated singular
steady systems refuse even if a constant field has zero residual. The existing
Robin response owner retains bulk `K + K'(T)T`, frozen interface conductance
and both contact-temperature feedback terms exactly once. Its bounded FGMRES
tangent and true transpose support heterogeneous bulk laws, matching face
means and nonmatching common-refinement quadrature without another primal
solve or differentiation through iteration history.

Nodal temperature, selected Robin mean-temperature and outward-heat goals
share the same complete adjoint. Controls include assembled nodal source loads,
uniform Robin reference temperatures, log heat-transfer coefficients and
each named interface's whole-law resistance multiplier. Resistance controls
are at fixed geometry, source claim and nonthermal context, with no physical
time factor. The local derivative does not evolve pressure or contact area,
select a derivative at source/material corners, or quantify source uncertainty.

Focused regressions compare independent slab-series derivatives and fresh
coupled solves on both trace meshes; mixed reservoir/source controls exercise
the true transpose. Endpoint tampering, singular insulated equilibria and
exhausted trace budgets refuse. These analytical fixtures establish no
measured-data validation.

## Mesh-goal comparisons with live source contact laws

`TemperatureContacts::compare_steady_goal` compares an unchanged, physically
admitted steady reference with another full nodal field on the identical mesh.
The original selected contact source is evaluated independently at both fields,
with both A/B temperatures and the declared face-mean or pointwise quadrature.
Bulk material assignments and temperature-dependent conductivity remain with
that same physical problem. Geometry, nonthermal source context, loads and
prescribed Dirichlet/Neumann/Robin data stay fixed. Boundary laws with feedback
need their own complete residual owner; this method does not add that feedback.

The reference must pass the original free-node residual and whole-domain heat
checks in watts. The reported relative primal residual is diagnostic and never
substitutes the derivative solver's tolerance for those physical gates. One
bounded transpose solve uses the complete bulk and contact Jacobian, including
both resistance temperature partials. No extra primal solve, claim selection,
contact freezing across fields or derivative through iteration history occurs.

Full nodal contributions are `lambda_i * [R(T_ref)-R(T_approx)]_i`, with zero
fixed-node entries. The signed direct goal difference and its observed finite
nonlinear remainder are retained separately; their identity includes remaining
linear-solve and arithmetic error. The remainder must remain in a refinement
estimate. These are discrete refinement indicators, not continuum error bounds.

Malformed or nonfinite fields/weights, altered prescribed values, invalid
linear/watt/trace budgets, extrapolated source queries, a dirty or nonsmooth
reference, failed transpose solve and cancellation refuse without publishing a
partial comparison. A finite approximate field may cross source segments or
land at a source corner because it requires only the physical law, not a
unique derivative there. Phase-state contact keeps its temperature-only refusal.

Independent two-slab regressions cover matching and nonmatching trace meshes,
explicit two-state residuals, the analytic complete transpose and a nonzero
quadratic finite remainder. Frozen-law comparisons fail those same oracles.
A constant source recovers the linear goal identity, and simultaneous nonlinear
bulk/contact feedback retains a quadratic remainder. Source corners, physical
endpoint tampering, source-domain exits, numerical budgets and cancellation
have direct regression cases. Their synthetic sources imply no measured-data
validation and the observed mesh difference implies no continuum guarantee.

`compare_steady_goal_with_robin_feedback` additionally binds each field's own
uniform Robin values on selected original regions. The boundary-law owner
supplies `[h, reference, dh/dmean, dreference/dmean]` at the reference, actual
`[h, reference]` at the approximation, and an additional row-major reference
feedback map between region means. Diagonal map entries add to local slopes;
the same derivative must never be supplied twice. Neither contact queries nor
the approximate boundary reuse the reference field's constitutive values.

The same complete contact/bulk transpose receives consistent P1 Robin
coefficient, local-reference and cross-region reference corrections. Both
physical watt gates remain in force before the dual, including for a zero
goal. Up to 64 distinct existing uniform Robin regions and four nodal factor
vectors per region must fit the separate feedback-entry allowance. Malformed
law/map dimensions, nonfinite values, nonpositive coefficients, duplicate or
unsupported regions, and exhausted storage allowances refuse before boundary
cloning. The original boundary partition, source claims, prescribed values,
input fields and contact ownership remain unchanged. Law evaluations and
slopes retain the consumer's responsibility; this supplies neither pointwise
boundary feedback nor an unrequested radiation model.

An independent four-temperature slab circuit exercises nonlinear Robin laws,
asymmetric affine reference coupling and live two-sided contact on matching
and nonmatching traces. Its directly solved transpose and finite nonlinear
remainder distinguish frozen approximate boundaries and omitted reference
feedback. Stale-reference, malformed/budget and cancellation regressions retain
the same refusal contract.

`compare_steady_goal_with_ambient_radiation` composes the same live-contact
comparison with the original area-mean gray-surface ambient law. It takes the
original convection/zero-flux partition, actual convection values at both
fields, reference convection slopes and additional reference feedback, and
the retained emissivity patches. Each field gets its own physical radiative
secant at its actual area-mean wall temperature. The reference's complete
contact/bulk transpose receives the radiative coefficient and weighted-reference
slopes; original convective reference feedback is scaled by
`h_convection/(h_convection+h_radiation)` exactly once. Local convection-reference
slopes and diagonal map terms are additive, with no duplicated derivative.

Original uniform zero-flux traces may carry radiation without convection.
Their convective coefficients, local slopes and feedback rows must remain
zero; their wall means may still drive another region's declared reference
through a feedback column. A positive placeholder convection reference cannot
replace the patch's retained radiation bath or create an air heat contribution.
All radiation patches must occur once in the selected regions. The existing
64-region, contact-trace and four-factor feedback budgets apply, alongside the
unchanged physical reference residual and whole-domain watt gates.

Only the reference requires an emissivity/source tangent. The approximate
field is evaluated physically and may land on a supported emissivity interval
endpoint; unsupported temperatures still refuse, including for a zero goal.
The mean observation uses the same face-sum and fused area accumulation as
the primal owner. No replacement primal, frozen approximate radiation, relaxed
reference gate or continuum accuracy claim is introduced. This composes the
mean-temperature secant model; it does not substitute pointwise fourth-power
quadrature or reflected enclosure exchange.

An independent four-temperature circuit with two radiation baths, nonlinear
convection, asymmetric reference coupling and two-sided live contact checks
the complete transpose and residual on matching and nonmatching traces.
The direct fourth-power oracle and quadratic finite remainder distinguish
frozen approximate radiation, frozen convection and omitted cross feedback.
A radiation-only slab checks physical emissivity-endpoint admission, source
domain exits, retained bath ownership, dirty references, exhausted budgets,
invalid zero-convection feedback and cancellation.

`AmbientRadiationGradient::compare_goal_at` uses the same reference binding and
physical-only approximate binding for continuous solids and fixed thermal
interfaces. Original uniform zero-flux radiating regions are admitted with the
same zero-convection ownership rules; ordinary Robin values are rebound from
the supplied law at both fields. Complete weighted-reference feedback and
physical source-domain checks remain active, including at approximate validity
endpoints. The existing mean-Robin goal owner retains the reference residual
check and complete transpose. No artificial convection is needed for vacuum.

`EnclosureRadiationGradient::compare_goal_at` and `compare_open_goal_at` supply
the corresponding complete comparisons for fixed or absent contacts and closed
or open reflected gray-diffuse enclosures. The live-source entry points are
`TemperatureContacts::compare_steady_goal_with_gray_diffuse_enclosure` and
`compare_steady_goal_with_open_gray_diffuse_enclosure`; all four share one
physical comparison owner. Original Robin rows contain `[h,r,dh,dr]` at the
reference and `[h,r]` at the approximation, with an additional row-major
reference-feedback map. The selected Robin list may be empty or omit fixed
convection. All physical radiation patches contribute their complete
`G (dq/dmean) W` operator independently of that list, including Robin/Robin
reflection and genuine authored nonzero Neumann fluxes. No placeholder
convection or shifted `r-q/h` reference enters either residual or tangent.

Both fields are reassembled by the original enclosure equilibrium owner with
their own fresh physical radiosity and original convection. Live contact laws
are queried anew at each field; only the reference takes complete independent
A/B source derivatives. A private endpoint binder admits the material/contact
tangent only after the actual direct-load free-node residual and global energy
closure pass explicit absolute watt gates. Radiosity equations and closed
radiative energy separately pass the explicit heat allowance at both fields.
An open bath is retained at its original temperature, including every reflected
path, and its algebraic heat is excluded from the solid's external account.
Original conduction must retain its thermal anchor.

Approximate material, contact and emissivity queries need physical support,
including supported source endpoints, but no unused tangent. Reference source
corners and validity endpoints still refuse, even for a zero objective. The
same source, exact mesh/region coverage, prescribed values and cancellation
checks remain authoritative. At most 64 Robin controls and 64 physical patches
are admitted; `4*n*controls + 2*n*patches` factor entries must fit the declared
budget, separately from the original contact trace allowance. Reported nodal
terms use the difference of the two physical residuals, and the finite
linearization remainder closes their sum to the unchanged goal difference.
The primal relative diagnostic uses the physical eliminated load norm; its
admission is still absolute W. The dual relative diagnostic is the residual of
the complete thermal transpose; the iteration count additionally includes
the independently checked radiosity transposes needed for its patch Jacobian.
Reflected radiation is a nonlinear tangent even with constant bulk/contact
laws. These are discrete refinement diagnostics, not continuum bounds.

## Temperature-dependent source laws on nonmatching thermal interfaces

`TemperatureContact::new_nonmatching` assigns an ordered source law to two
explicit planar P1 trace meshes. `TemperatureContacts::with_nonmatching`
admits complete face ownership and coverage once, then retains the original
common-refinement geometry. Every trial interpolates the independent A/B
temperatures at the same physical quadrature point and queries the pinned
source there. A mean-temperature law must be selected explicitly; it uses
the two local trace values. Fixed pressure, finish, medium, history and claim
selection retain their original context.

The existing positive three-point rule on each overlap triangle integrates
one shared flux into equal-and-opposite weak loads. This is a quadrature
approximation to a nonlinear contact law, not an exact continuum integral.
Geometry and quadrature-count budgets are explicit. Exactly coincident pieces
of a variable nonmatching law use this pointwise rule exactly once; unrelated
matching contacts retain their existing face-mean semantics and arithmetic.
No geometry intersection is repeated by a constitutive trial.

The projection reference and clipping operands are canonical across both
traces. Reversing A/B together with the ordered source coordinates retains
the same overlap fan and quadrature points, including nonlinear contact
tangents. This symmetry leaves the original quadrature order unchanged.

`TemperatureContactEvaluation::quadrature_states` retains each ordered face
pair, canonical quadrature index, area weight, actual trace temperatures and
original source receipt. The complete contact Jacobian differentiates both
temperature coordinates through the local P1 traces and includes the
resistance feedback. Its transpose and whole-law resistance controls feed
the existing capacity and primary-enthalpy endpoint owners. A frozen
resistance snapshot alone still omits that feedback.

Source refusal, derivative corners, exhausted work limits or cancellation
publish no partial contact state or thermal history. Mechanical pressure,
gap, surface-area evolution and automatic quadrature error control remain
outside this fixed-geometry extension. Focused regressions exercise an exact
constant heat flux with spatially varying source resistance, mixed overlap
and coincident patches, a manufactured transient energy balance, source
reversal, the full Jacobian/transpose and fresh-solve resistance sensitivities.

## Pointwise source laws on matching thermal interfaces

`TemperatureContact::new_pointwise` binds an explicitly ordered mean-temperature
or two-sided source law to exactly paired P1 triangles.
`TemperatureContacts::with_matching_quadrature` requires separate face and
quadrature allowances; mixed independently tessellated contacts continue
through `with_nonmatching`. The original matching owner admits every pair and
exclusive boundary ownership. Three positive degree-two points per triangle
are built directly from its exact vertex correspondence, so faceted curved
surfaces need no global planar approximation or clipping. The existing
pointwise source, weak-form, flux and derivative owner is reused, with each
contact assembled once. `surface_is_nonmatching` remains false for exact pairs.

Every integration point queries the original pinned card at both local traces
under the declared coordinate convention. Paired vertices also pass the source
query, preventing an unsupported hot/cold corner from hiding behind a supported
mean or interior rule. These finite guards do not certify that an arbitrary
nonconvex source domain has no unsampled holes. Strict derivatives additionally
require the source-coordinate neighborhoods at guards; forward scalar/curve
Newton keeps its existing semismooth policy, with every actual residual still
source-checked. All nonthermal coordinates, geometry and source selection stay
fixed.

A frozen snapshot retains all three actual source receipts per face in
`quadrature_states`; `states` remains the face-mean API's output. Residual,
signed A-to-B heat, complete contact tangent/transpose and whole-law resistance
pullbacks use the same samples. Existing steady, capacity, primary-enthalpy and
radiation-coupled consumers inherit this rule. Three-point integration is exact
for constant resistance and the corresponding polynomial integrands; general
nonlinear reciprocal laws remain quadrature approximations. Existing face-mean
constructors preserve their declared semantics. Seven regressions in
`tests/pointwise_matching_contact.rs` and its transient child cover unresolved
hotspots, analytic flux, bent surfaces, source support, budgets, derivatives,
closed capacity energy and latent enthalpy.

## Two-sided source laws on matching thermal interfaces

`interface::temperature::TemperatureContact::new_two_sided` binds separate
absolute-temperature source axes to the independently resolved A/B P1 trace
means. An original interface card may therefore supply R(T_A,T_B) without
replacing unequal body temperatures by their arithmetic mean. `new` retains
the explicitly declared mean-temperature closure. Ordered face assignments
and source axes are never silently sorted or exchanged; pressure, medium,
finish, history and the selected claim stay fixed at their original context.

The original source evaluator owns interpolation, holes and support limits.
Each snapshot retains the actual local source receipt. The reported arithmetic
mean is diagnostic for a two-sided law. Resistance is constant on each paired
triangle at those trace means, while the original consistent face mass
integrates its nonuniform temperature jump conservatively. This is not a
quadrature-point material law or an evolving mechanical contact-area model.

The same snapshots feed existing capacity, primary-enthalpy and radiation
endpoint owners; a failed source query publishes no candidate history.
The complete tangent and its true transpose include both source partials,
using three P1 mean weights on each side and the original consistent face
mass for the full nodal jump. Scalar, curve and original triangular surface
responses share this assembly; each surface derivative is prepared once per
law per linearization. Fixed pressure and other external axes do not become
thermal derivative coordinates. Actual pinned source queries still admit
both neighbours of every declared thermal coordinate, including coordinates
that affect applicability alone. Equal-gradient cell edges are admitted;
source corners, unsupported wedges and thermal domain ends refuse a classical
tangent even when the primal point value is valid.

Focused G0/G1/G3 tests cover ordered source queries, surface/axis reversal,
actual implicit energy exchange, per-side finite differences, transpose and
conservation identities, fixed-pressure support boundaries, and invalid
coordinate or source-support refusals. A primary-enthalpy endpoint adjoint
also checks resistance-scale, prior-enthalpy and source-density directions
against fresh solves while a non-contact node remains inside a latent plateau.

## Source-bound implicit finite-cell contact

`transient::deposition::contact::implicit::PreparedMaterialContact` binds the
original ordered material-contact source and caloric objects to the native
implicit chain/graph. It pins the initial physical interface claim and prepares
its original source line or triangulated response once, retaining source
coordinates, ordinates and the original physical query receipt.
Constant scalars, linear curves and declared triangulated linear surfaces in
either side's temperature or liquid mass fraction, the declared mean temperature,
or a fixed external coordinate are admitted. Surfaces use the original material
evaluator and its cached analytical cell gradients; each response coordinate
contributes to both enthalpy derivatives as required by the caloric chain rule.
Source holes remain active inside declared enthalpy bounds. No finite-difference
source law, filled support box or temperature/fraction mixing rule is inferred.
`source_line()` and `source_surface()` expose the applicable prepared response.

Every trial uses original numeric caloric values and the complete analytic
enthalpy chain rule. Mean-temperature response includes both half-weighted
temperature derivatives; phase response remains active on a latent plateau.
All original applicability bounds follow their actual two-sided coordinates,
including coordinates absent from the value formula. Scalar/curve runtime
evaluations allocate nothing; surfaces retain the canonical robust predicate
and expansion allocation behavior. Neither path performs source queries or hashing.
Typed preparation/source failures remain available; the native trait preserves
its `MaterialDomain` contract. Original reference receipts remain labelled as
reference queries, while trial values are compiled numerical responses.

The native thermal owner applies resistance once in series with actual bulk
conductivity and owns the residual, Jacobian, energy, entropy and time history.
Focused cross-owner regressions cover unequal-mass latent/sensible exchange,
material-only resistance changes, simultaneous nonlinear bulk conductivity,
constant-source compatibility, simultaneous phase/temperature surface response,
both-side and mean-temperature/phase chain rules, missing source cells, source
domain refusal, cancellation and clean retry. No pressure/gap/area evolution or
interfacial heat storage is added by this stationary-contact reduction.

## Radiation rate bounds on the shared thermal timeline

`RadiationCellSources::conductance_bounds_w_k` supplies conservative per-cell
rate scales from the same resolved gray-diffuse surfaces, material emissivities
and caloric temperatures as the actual radiation powers. It adds no heat links,
source powers, radiation solve or thermal state. For exactly closed view-factor
rows, each surface contributes at most `4 sigma A epsilon Tmax^3`; disjoint
surfaces on one cell add. Self-view and same-cell exchange do not increase that
bound. The derivation uses the original nonnegative net-radiation operator.

Admitted view-factor row defects are retained. A positive reflected-row margin
and an explicit amplification factor bound their contribution; a missing or
unrepresentable bound refuses. Geometry is never renormalized, and the original
radiosity/power and closed-energy gates remain authoritative. The estimate is
for frozen emissivity and does not certify passivity of defective geometry,
temperature-dependent emissivity derivatives, or temporal accuracy.

`fs-couple::MaterialHeatObserver` supplies these rates to the original thermal
additional-conductance limiter. Conduction, caloric mass/slope and radiation
therefore restrict the SAME accepted explicit substep. Source/event epochs and
rollback are preserved. Focused finite-plate regressions compare the native
trajectory with an independently evaluated explicit recurrence, vary only
emissivity, and check work-budget refusal. Source/time refinement remains needed.

## Material contact sources on the shared thermal timeline

`transient::deposition::contact::MaterialContactObserver<O>` joins the original
ordered interface-source adapter to any compatible bulk-only deposition
observer. Each selected active link must appear explicitly in the inner owner's
fresh conductance map. The existing `apply_material_contacts` composes
`1/G_bulk + R''/A` once; untargeted coefficients and all other physical hooks
remain with the inner observer. Do not wrap an owner that already included the
same interface resistance. The same corrected G reaches native pairwise heat
and the native CFL bound, with no new thermal solve or heat-source injection.

Attachment resolves and pins the original contact claim, checking each side's
exact borrowed caloric owner, ordered material card and physical cell. Probing
uses a cloned inner observer in its original begin/sources/conductances order;
the retained inner physical history does not advance. Bound reference mass,
body/segment/material identity, cell geometry, link endpoints, nominal area and
separate surface coverage cannot silently change. Reprepare for changes to
these physical bindings. Authored link conductance may change because the
fresh bulk/contact law replaces that coefficient explicitly.

Every actual starting state reevaluates contact resistance from the two distinct
caloric temperatures and phase fractions. Pressure and other admitted external
conditions remain explicitly held; this does not infer mechanical contact load,
area, interfacial storage, spreading resistance or an arbitrary moving-contact
law. A zero area/bulk conductance remains insulating but still validates source
data. Original selected claims cannot switch when a support boundary is reached.
The endpoint and post-event checks use the same source owner before the candidate
is accepted, including when no later substep would query that endpoint.

`last_samples` retains the attachment or most recent starting-state bulk/contact
series evaluation. `endpoint_sources` separately reports actual final interface
sources, without calling frozen bulk coefficients final-state evaluations.
All observer hooks—including interval-dependent heat, reservoir conductances,
internal transfers, events and timestep restrictions—are forwarded. The original
bulk/source/checkpoint errors remain typed. The existing deep-clone transaction
rolls back physical observers and enthalpy together on any refusal.

This is the existing first-order explicit coefficient-freezing scheme; frozen
CFL admission does not establish arbitrary nonlinear time accuracy. Tests in
`deposition_contact/observer.rs` exercise complete hook forwarding, fresh series
composition, actual latent-phase/contact feedback, conserved internal exchange,
original bulk errors, source pinning and failed-step retry. The cross-physics
fixture `fs-couple/tests/material_contact_motion.rs` composes the same observer
with continuing material-driven friction, heat and articulated motion. Native
execution must be reported separately from parsing and source review; these
synthetic regressions do not establish measured contact accuracy.

## Source-paired spatial sensible-enthalpy transport

`transient::enthalpy::material::MaterialEnthalpyTransport::from_isotropic_card`
binds Cp(T), rho(T) and scalar k(T) from one immutable material card at fixed
typed non-temperature conditions. The constructor explicitly declares an
isotropic transport model; source quantity dimensions alone do not establish
isotropy. All three properties must have explicit temperature validity over
the requested path. Missing-temperature producer scalars therefore cannot
silently become heating curves.

`from_orthotropic_card` uses that SAME storage/source path with three
conductivity property keys and a declared finite orthonormal principal frame.
All directions must have supported source data over the full caloric interval;
no missing principal value is copied from an unrelated scalar. The existing
tensor owner rotates both k(T) and dk/dT into the physical mesh frame, so
material orientation changes the actual nonlinear transient operator. Source
identities and the fixed physical frame participate in the paired material
identity; orientation cannot silently reuse an isotropic or differently
oriented material identity. This stationary model does not evolve the frame
with strain, phase or microstructure.

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

Focused G0/G1/G3 regressions cover source/reference pairing and refusal, the
independent variable-Cp heating balance, frozen mass, and nonuniform
equilibration with doubled k(T) without net heat creation. Directional tests
check the independently rotated tensor and nonlinear tangent, missing and
nonfinite frame/data refusals, orientation-dependent heat spreading in the
original PDE, and material/field rotation covariance with conserved energy.

```sh
cargo test -p fs-conduction --test material_enthalpy
```

## Centre-regular generated solid cylinders and spheres

`layered::SolidRadialLayeredBody` binds an original `LayeredWall` to a solid
core followed by concentric shells. The first physical layer thickness is the
core radius. Symmetry requires Q(0)=0 and zero inner film; the only reservoir is
outside the body. Cylinders have insulated endcaps and spheres cover the full
solid angle. Original material tables, cards, receipts, physical contacts and
the outer film remain retained.

`solve_temperature_dependent_with_generation` admits one signed uniform
volumetric source per layer. Enclosed source determines every radial heat rate.
For a core, q(r)=s r/m and K(T_c)-K(T(r))=s r^2/(2m), with m=2 for cylinders and
m=3 for spheres. This explicit regular limit never evaluates a hollow-shell
logarithm or inverse radius at zero. The outer film fixes its face temperature;
one inward pass uses the original source-segment Kirchhoff inverse, shared shell
moments and actual-area contact jumps. No inlet-flux shooting or second material
law is introduced. Work is linear in layer count after the existing bounded
source preparation; cancellation retains the original layer/segment checkpoints.

Every face and possible internal extremum must remain inside the declared
common original source interval. The result reports actual areas, signed local
Q/A flux, total face heat, layer generation, hot/cold spots and the one-reservoir
balance against an independently accumulated source-times-volume inventory.
The centre has exactly zero area and regular zero local flux.
`temperature_at` and `heat_flux_at` recover the retained nonlinear profiles.
Constitutive and film defects are independently recovered from the rounded
temperatures; absolute flux tolerances use the core's outer area as reference.

This is fixed radially symmetric steady conduction with generation independent
of temperature. It supplies no transient heat capacity, phase-front motion,
reaction/Joule feedback or material uncertainty enclosure.
Five integration regressions cover signed nonlinear exact profiles, real shell
contacts/films, independently derived internal cold spots and source refusal,
homogeneous subdivision, length scaling, zero source and cancellation. Two
native unit tests exercise the actual regular core geometry and dimensional
generation integrals across signed sources and physical scales.

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

`LayeredWall::solve_with_generation` extends that same constant-conductivity
owner with a declared uniform volumetric source in every layer, in W/m3.
It integrates `dq/dx = s` and `dT/dx = -q/k` analytically. The face fluxes
therefore differ by the layer's generated heat per area; every contact and
surface film uses its own local flux. One scalar balance closes the two
reservoir temperatures, with work linear in the layers and retained source
knots and no mesh or iterative solve. Original material tables, cards,
receipts, contact assumptions and source densities remain in the result.

The result reports both outward reservoir heat rates, total generation and
energy residual, separate layer faces, spatial mean temperatures and exact
quadratic extrema with their positions. Both faces and every interior
extremum must be positive Kelvin and within that material's source span; a
supported pair of faces cannot conceal an unsupported hot spot or cold spot.
Signed source densities declare production/removal without inventing a Joule,
reaction or refrigeration mechanism. Fourier/film/contact residuals recomputed
from the rounded temperatures and the final boundary defect are numerical
diagnostics, not material uncertainty or an experimental validation claim.
The nonlinear generation entry point below handles k(T). Temperature-dependent
generation and lateral heat spreading still require the spatial/coupled owners.

Seven focused G0/G1/G3/G4 regressions in `tests/layered/generation.rs` encode
the independently derived heated-slab parabola, unequal reservoir heat split,
heterogeneous contacts, signed and cancelling sources, reversal, area scaling,
subdivision, zero-generation parity, original receipts and interior support
refusals. Native execution was blocked in the authoring environment by shared
disk exhaustion and the unchanged asupersync dependency exceeding the shared
8 GiB memory limit; no passing Cargo result is claimed for this addition.

Varying conductivity is refused by the constant solve entry point. The separate
`LayeredWall::solve_temperature_dependent` entry point integrates every positive
piecewise-linear source segment with a scaled analytic inverse and bounded
monotone common-flux shooting. Each source table must support its own actual
layer-face temperatures at its original fixed nonthermal context. Layer-local
intervals may be disjoint across finite contact resistance, and reservoir
temperatures may lie outside a material's support behind a surface film.
Positive conductivity makes every downstream temperature decrease with the
shooting flux; an unsupported hot/cold trial therefore narrows that same flux
bracket without evaluating an out-of-domain constitutive value. Caller-declared iteration,
temperature and constitutive-flux limits gate the returned candidate. Cancellation
is checked between layers and every 128 source segments. Effective resistance and
transmittance belong to this boundary pair; they are not tangent conductances.

`LayeredWall::solve_temperature_dependent_with_contacts` composes that same
bulk solve with indexed `LayeredTemperatureContact` source bindings. The shared
contact-coordinate owner preserves explicit mean or ordered left/right source
temperatures; pressure and other supplied context remain fixed. Each selected
interface claim stays pinned. A bounded Picard loop reuses the working stack,
changes only contact coefficients, and relaxes only source temperature
coordinates. Acceptance queries the original laws at the unchanged returned
interface temperatures and checks physical contact flux against the common
stack flux using the original flux tolerances. Reports retain applied and fresh
source answers, both interface temperatures, physical residuals and all counted
wall/shooting responses. Unlisted contacts and surface films stay unchanged.
Equal reservoir temperatures resolve the actual equilibrium resistance before
the bulk solve, so zero heat cannot mask a wrong local resistance or receipt.
Focused G1/G3/G4 tests in `tests/layered_nonlinear/contact.rs` cover nonlinear bulk
plus varying and fixed contacts, equilibrium, ordered-source heat reversal,
source support, bounded failure and cancellation. No extrapolation, claim
switching, uniqueness, arbitrary-feedback convergence, continuous multiaxis
source-path support or propagated uncertainty is inferred. This entry point is
planar and source-free; the generated variant below shares its contact-feedback
kernel. Radial and transient contact feedback remain with their corresponding
owners.

`LayeredWall::solve_temperature_dependent_with_generation_and_contacts` composes
that same bounded feedback with the existing nonlinear generated-wall owner.
Every interface uses the original outgoing/incoming local fluxes, and each side
must match the fresh pinned contact law under its own flux tolerance. Both
native fluxes and their continuity residual are retained. Source-coordinate
relaxation does not blend resistance values, returned temperatures or receipts.

Prescribed signed generation is passed once to the original bulk solve. Its
material-domain and interior-extremum checks, films, both reservoir conditions,
heat split and complete source profiles remain authoritative. The original
reservoir-minus-generation residual must also pass the native maximum-face-flux
tolerance multiplied by area. Equal reservoir temperatures or zero NET generation
do not imply zero internal flux; only a complete all-zero layer source vector
permits equilibrium contact initialization. Independently supplied volumetric
Joule loss is supported as prescribed heat; no electrical or temperature-dependent
generation feedback is inferred.

Focused analytical tests in `tests/layered_nonlinear/contact/generation.rs`
cover equal-bath nonlinear heating/cooling and reflection with contact support
excluding the baths, opposite local contact fluxes under zero net generation,
source-free and fixed-source parity, unchanged films, original receipts,
interior extrema, bounded refusal and cancellation. These definitions and
independent numerical equation checks do not establish native test execution.

`LayeredWall::solve_temperature_dependent_with_generation` combines the same
source-segment integral/inverse with `dq/dx=s` for finite, signed, uniform
per-layer production/removal. The Kirchhoff potential is quadratic, while
temperature is recovered from the original positive piecewise-linear k(T).
Its possible interior extremum is recovered analytically and admitted before
a candidate can be returned. The inlet flux is determined by bounded monotone
shooting; face fluxes, contact jumps and surface films all use local flux.
`NonlinearGenerationOptions` declares a positive temperature envelope plus the
existing work/residual limits. Each layer intersects that envelope with its own
original source interval and checks both faces and interior extrema against
this local support. Reservoirs may lie outside the envelope behind explicit
nonzero films. No rejected trial evaluates
an unsupported conductivity. A single O(layers) trial buffer is reused, and
source inversion begins at the active segment found by binary search.
The result retains material-card identities and receipts, both outward reservoir
heat rates, signed generation and energy balance, local face fluxes, extrema,
the numerical inlet-flux bracket and independently recovered constitutive
residuals. `temperature_at` recovers the nonlinear interior profile from the
retained segment integrals without a mesh or new material queries. Five G1/G3/G4
tests in `tests/layered_nonlinear/generation.rs` add affine/kinked source oracles,
internal removal, heterogeneous contact/film balance, reversal, subdivision,
area scaling, constant/zero-source parity and interior-domain/budget refusal;
existing source/card and film-support tests exercise the new path too.
Generation is prescribed independently of temperature; neither this path nor
its roundoff residuals claim Joule/reaction feedback, material-error enclosure,
transient accuracy or experimental qualification.

`RadialLayeredWall` binds the same original source stack and generated nonlinear
solve to positive-radius hollow cylinders (insulated ends) and full spherical
shells. It retains physical radial thicknesses, the actual area of each contact
and each reservoir film, and the exact generated volume of every layer. A small
cached geometry leaf supplies integrated reciprocal area and distributed-source
moments to the existing Kirchhoff/shooting owner. Cylindrical logarithms use
ln_1p, thin-shell source moments avoid subtractive cancellation, and stationary
radii invert swept volume without subtracting nearly equal radii. Conductivity
integration, layer-local source admission, interior extremum checks, work
limits and cancellation remain shared with the planar path.

Reports distinguish total signed radial heat rate from the local Q/A flux;
outward reservoir heat closes the signed volume generation. Profile recovery
uses the retained original material integrals. The inner face area is only a
numerical normalization: the absolute flux tolerance refers to that area, and
reported constitutive residuals are converted to watts. Cylinder length scales
powers without changing its radial temperatures. Five focused tests in
`tests/layered_nonlinear/radial.rs` cover nonlinear logarithmic/inverse-radius
conduction, internal hot/cold spots, multilayer contacts and films, subdivision,
length/thin-shell limits and source/geometry/budget/cancellation refusal.
A zero-radius solid centre, axial/lateral spreading, moving geometry, nonlinear
generation feedback and transient radial storage are outside this shell model.
Seven cases in `tests/layered_nonlinear.rs` cover affine and interior-knot laws,
contacts/films, reversal, equal temperatures, original source receipts and bounded
refusal. No density/Cp, transient history, lateral spreading,
radiation, moisture transport, propagated material uncertainty or experimental
qualification is inferred. Geometric/process matching beyond the explicit
temperature span remains the caller's duty, including thickness-dependent
product declarations. G0/G1/G3/G4 cases in `tests/layered.rs` exercise the
analytical series-resistance law, contact jumps, exact source context,
material/thickness substitution, reversal, subdivision and transactional refusal.

Local-source regressions additionally recover independent Kirchhoff/contact
balances for disjoint hot and cold materials, films outside material support,
signed generation, cylindrical and spherical shells, and heated solid cores.
Prescribed faces, actual interfaces and interior extrema outside their original
source support still refuse. The source integrals, inverse profiles, physical
areas, heat accounts and caller tolerances retain their original owners.

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

`enthalpy::deposition::fracture::FractureCooling::synchronize_with_heat`
optionally converts an explicit per-interface fraction of NEW cohesive
dissipation into heat at the original two finite cells. One authoritative
cohesive trial supplies the work; elastic storage, initial damage, unloading,
reclosure and numerical defects never become heat. The original caloric law
resolves temperature and latent storage, and contact/history/enthalpy changes
publish together. Nonthermal fracture energy remains separately reported.
Partition, aggregation, cell-storage and report rounding share one declared
joule budget; compensated totals retain small terms beside large work values.
The legacy contact-only `synchronize` retains its arithmetic. This is fixed-area,
staggered endpoint heating with immutable toughness, not resolved crack-tip
temperatures, within-step thermomechanical equilibrium or a full entropy model.
Seven focused heat-conversion regressions extend the existing contact tests.

`synchronize_over_time` accepts the original constant-rate straight cohesive
segment, including a physical surviving-bond viscosity. It requires the explicit
duration to equal current thermal time minus the retained last mechanical
synchronization time. Intermediate heat-transport stops do not reset that clock.
The caller derives duration by subtracting those accepted absolute times and
uses that represented interval in the mechanical update, rather than a separately
rounded nominal timestep.
The original `trial_over_time` law replays and checks the complete supplied point;
foreign laws, inconsistent irreversible history, unsupported rates, stale/zero durations and
cancellation refuse atomically. Untimed calls still reject viscous laws and
retain the rate-independent zero-time replay behavior.

`synchronize_over_time_with_heat` uses independent fracture and viscous heat
fractions and independent receiving-side splits through `CohesiveHeatPartition`.
Both physical losses and both nonthermal remainders stay separate. The same
compensated partition and native caloric update account for all converted work,
including viscosity accumulated before complete failure when endpoint dashpot
traction is already zero. A second synchronization at the accepted time refuses;
a later stationary hold produces no repeated loss. Six focused G1/G3/G4 groups
cover independent intact/bilinear dashpot integrals, unloading, mixed shear
failure, separate heat channels, actual mechanical time, source/history/rate
refusals, latent-domain failure and cancellation after staged heat deposition.
These endpoint APIs do not reconstruct unobserved mechanical paths, localize
fracture events, evolve temperature-dependent coefficients or perform simultaneous
thermomechanical equilibrium. The prescribed event path retains its existing
rate-independent admission. Every physical reversal requires its own supplied
segment: an undamaged closed loop has no retained point-history signature, so
neither its existence nor its viscous work can be recovered from endpoints.

`FractureCooling::transfer_temperature` explicitly connects the actual caloric
cell temperatures to source-defined common- or mixed-mode bond weakening. Each
complete update retains the original ordered interface card, selected claims,
nonthermal coordinates, area, jump and irreversible history. The successor's
temperature must equal the explicitly weighted two-cell temperature returned by
`interface_temperature_k`; a stale or substituted source refuses. Difference-first
weighted interpolation is exact for equal temperatures and
receiving-side endpoint weights, so a uniform field at a source-domain boundary
cannot drift outside it through weighted-sum roundoff.
The original cohesive temperature law supplies fixed-jump parameter work, storage change and
new damage loss, integrated on original area with a declared joule gate.
Thermal failure immediately opens the existing thermal contact or selects its
separately declared closed-failed conductance. Cooling cannot heal or rebond it.
Both mechanical and thermal clocks must already describe the same endpoint, so
an unfinished timed mechanical segment cannot be replayed under changed laws.
No heat is inferred from parameter work or damage release, and no dashpot work
occurs during the zero-duration transfer. Four focused G1/G3/G4 test groups cover
heating/cooling through common, mixed and viscous sources, partial mode history,
actual postfailure heat exchange, source/temperature refusals and cancellation.
This is staggered endpoint coupling, not within-step fracture localization,
resolved interfacial temperature, a cure law or a complete thermodynamic closure.

> Beads `frankensim-extreal-program-f85xj.5.1`, `.5.3`, `.5.4`, and `.5.6`. The
> cooling vertical's steady heat-conduction solve, matching-P1 contact rung,
> first card-backed surface-radiation rung, and receipt-preserving PCB
> homogenization handoff.

## Frequency resolution on the native thermal history clock

The shared event endpoint selector applies `LinearChirp::resolved_step_end`
to authored boundary values, Robin reference temperatures, volumetric sources
and additional supplied histories before selecting any solver stages. It uses
the scenario default of at most 1/16 cycle per interval, including both tails
of a sweep. Existing knots and exact held-value left limits remain intact.
Fixed event steps and adaptive full/two-half comparisons therefore cannot all
skip complete source cycles and mistake oscillating heat for a constant input.
The native nonlinear solve, heat balance, local error controller, lifetime
budget and cancellation continue to decide acceptance and atomic publication.

Direct caller-prescribed endpoint APIs retain their existing semantics; this
automatic guard belongs to event-aware history stepping. The phase allowance
is not an antialias filter, a nonlinear bandwidth bound or a global time-error
certificate. `thermal_history_adaptive` includes an independent adiabatic
cosine-heating regression whose exact integrated heat differs from the
deliberately aliased raw endpoint result.

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
Explicit pointwise accepted-state binders attach the same spatial operator to
the enthalpy transpose. The base binder retains its ownership refusal for a
pointwise endpoint, including zero-net-radiation equilibrium.

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
only their authored heat. Explicit pointwise accepted-state binders retain the
complete reflected spatial tangent and its genuine transpose. Tests in
`tests/enthalpy_enclosure/pointwise.rs` compare independent simplex moments,
reflected derivatives, nonuniform nodal balances, finite latent heat transfer,
and a strong-self-view case with a negative exact diagonal.

`PointwiseEnclosureRadiation::new_open` uses that same local P1 emission,
reflected irradiation and complete nonsymmetric tangent with the explicit
black reservoir of `OpenDiffuseEnclosure`. The reservoir's reciprocal row
contributes radiation and an independently evaluated heat account, never a
solid node or artificial storage. `PointwiseEnclosureReport::outward_heat_w`
remains the physical-surface sum used exactly once by the existing homogeneous,
heterogeneous and adaptive enthalpy energy accounts. `ambient` retains bath
temperature, effective area, radiosity, irradiation and signed outward heat;
`system_energy_closure_w` gates the sum of surface and bath powers. Original
Robin/Neumann heat remains unchanged. Finite storage permits an otherwise
unanchored, insulated body to cool radiatively or absorb real latent heat.

The open spatial owner retains an immutable snapshot of the prepared model;
its lifetime only borrows mesh/boundary. `with_ambient_temperature_k` reuses
the same factors and sparse geometry for a new prescribed endpoint bath
temperature without advancing history. The existing accepted-state adjoint
retains the exact spatial tangent/transpose and surface-emissivity controls;
its reservoir-temperature gradient contains one shared black-bath control.
`tests/enthalpy_enclosure/pointwise/open.rs` checks independent simplex-moment
loads and reflection, tangents/transpose/parameter derivatives, every nodal
thermal balance, radiation-only transient cooling, latent-plateau heating and
cooling, and the implicit bath sensitivity. These remain fixed-emissivity,
patch-uniform-irradiation models: no spatial visibility kernel, spectral bath,
emissivity-state derivative, or automatic bath-temperature time schedule is
inferred.

`PointwiseEnclosureRadiation::new_spectral` and `new_open_spectral` extend
that same spatial/enthalpy owner to 1..=64 explicit contiguous wavelength
bands. Each band supplies independent opaque diffuse surface data with its
source and temperature support. Emissivity equals absorptivity in that band;
reflectivity is their complement. The admitted enclosure supplies traces,
areas and reciprocal geometry only: its total-gray material cards are not
used or reported as spectral measurements. The returned patch source is
explicitly `PointwisePatchEmissivity::Gray` or `Spectral`.

Each band prepares its radiosity factorization once and shares sparse trace
topology. Local Planck-band emission uses the actual P1 temperature field;
the complete reflected tangent and true transpose sum the original band
operators without a dense nodal radiation matrix. The existing homogeneous,
heterogeneous, adaptive, endpoint and accepted-adjoint paths consume this same
owner. Spectral emissivity cotangents retain their wavelength rows separately
from gray controls, and the shared open-bath derivative sums every band once.
Accepted-adjoint workspace admission includes the retained band count.

The shared `radiation::spectral::BlackbodyBand` primitive supplies finite-band
power, its analytic derivative with wavelength endpoints fixed, stable
near-equilibrium differences and spectral coverage fractions. It uses NIST
Planck conventions and defining SI constants while retaining this crate's
Stefan--Boltzmann normalization. Positive adaptive spatial quadrature estimates
both emitted loads and temperature tangents with explicit depth/evaluation
limits; nonpolynomial spectral emission has no degree-five exactness claim.
The finite window must meet its explicit omitted-blackbody-fraction budget
over the complete P1 temperature envelope and the black bath. This coverage
gate is not a certified net-heat error bound, and omitted wavelengths carry
no invented heat. Per-band reports retain original sources, powers, algebraic
closure and spatial work/error estimates.

Focused spectral regressions exercise selective plates against independent
Planck/resistance equations, nonuniform open loads and transposes, every
original nodal enthalpy equation, actual latent storage, accepted bath
sensitivities, source/coverage refusal and bounded cancellation. Five actual
Planck leaf tests passed through an isolated Rust wrapper around the original
source; full-crate spectral integration tests remain pending native execution.
This is fixed spectral-band emissivity with patch-uniform irradiation. It
does not introduce a spectral material database, temperature-dependent
emissivity, directional/specular optics, participating media in the solid
owner, automatic visibility, a native spectral project declaration or an
experimental-validation claim.

`transient::enthalpy::heterogeneous::EnthalpyRadiation` selects the existing
prescribed, mean ambient, mean enclosure, pointwise ambient or pointwise
enclosure owner through `advance_with_radiation`. The selection borrows the
original prepared operators and retains the same reference masses, phase
charts, conductivity, contacts and immutable history. It introduces no new
iteration or heat law. `EnthalpyHeatStep` preserves the complete original
result: mean-patch physical residual/energy and separated convection accounts
remain available beside their applied inner boundary; pointwise reports remain
inside the native step. `conduction()` is a state/inner-solve view, while
`physical_energy_residual_j()` selects the authoritative complete heat audit.
Consumers must not send a mean-patch combined Robin heat into a coolant port.

Both pointwise assemblies expose `tangent_apply_transpose` and
`parameter_pullback`. The latter contracts a caller-supplied full nodal
outward-load weight with absolute surface emissivity and, for ambient exchange,
reservoir-temperature controls in declared surface order. The enclosure
contraction keeps each spatially weighted emission integral and differentiates
reflected irradiation through the already prepared patch response. It adds no
parameter linear solve. Sharing an emissivity source does not merge distinct
surface controls. Geometry, view factors and the immutable emissivity query
remain fixed; every surface node must be strictly inside the source temperature
support for a classical derivative.

`linearize_accepted_with_pointwise_ambient_radiation` and
`linearize_accepted_with_pointwise_enclosure_radiation` are available on both
homogeneous and heterogeneous enthalpy engines. Callers supply the original
accepted physical residual target and an explicit radiation workspace cap.
The producer freshly evaluates the complete physical residual and energy,
checks constitutive smoothness and retains one active face/patch assembly.
Its transpose applies the physical free-row mask before the transposed
radiation action, and the timestep and temperature-to-enthalpy slopes on the
correct side. Fixed columns remain present for prescribed-temperature controls.
Natural/air Robin and contact feedback attachments retain that operator when
refreshing their diagonal.

The existing history, source, fixed-bath, conductivity and reference-density
pullbacks consume the same complete dual. `pointwise_parameter_pullback` on
the checked transport linearization accepts outward-load weights; an objective
uses minus `nodal_load_pullback`, whose timestep and prescribed-row mask are
already included. The shared `adjoint::pointwise::workspace_entries` admits
bounded nodal scratch, sparse face storage and the small patch matrices with
checked arithmetic. This is an entry allowance, not a peak allocator/RSS proof.
A latent plateau still has zero direct temperature seed while energetic
enthalpy history and its boundary parameter sensitivities can remain nonzero.

`BackwardEuler::with_surface_heat_capacity` adds thermally thin, perfectly
bonded coating storage on original boundary-face slots. `ArealHeatCapacity`
declares finite nonnegative J/(m²·K), and `SurfaceHeatCapacity` supplies uniform
or P1 face coefficients in the original vertex order. The exact row-sum moment
is `A (2 C_i + C_j + C_k) / 12`, preserving the correlation of a P1 capacity
(or thickness) with the temperature trace. It is added to the existing volume
capacitance; linear/nonlinear endpoints, contacts, radiation, prescribed-node
reactions and discrete adjoints therefore share one combined storage operator.
The prescribed reaction includes the surface energy change divided by the
step duration. Empty and zero additions preserve the original operator exactly.

The surface shares the solid boundary temperature: perfect thermal bonding
and negligible through-thickness gradients are explicit model assumptions.
This is frozen linear storage, with no separate film temperature, added
contact resistance, moving geometry or latent heat. Tangential transport is
declared independently through `ThermalBoundary::with_surface_heat_conduction`.
Storage in a separately meshed coating must not be added again.
`SurfaceHeatCapacity::from_material` binds three declared nonnegative reference
thicknesses to an already admitted constant `MaterialHeatCapacity`; it retains
the original density/Cp source and intersects temperature support at every face
vertex with a nonzero storage moment. It never samples a varying caloric law
and freezes it. Duplicate faces within one addition, disjoint active source
domains, non-finite/unrepresentable moments and cancellation refuse preparation.
Successive additions intentionally accumulate further layer storage.

`transient::MaterialHeatCapacity` resolves constant density and specific heat
from one immutable material card at an explicit reference state, then pins both
selected sources at the endpoints of a declared positive temperature interval.
Validity-wide scalar static laws and curves/surfaces depending only on typed
held external conditions, including pressure, moisture and composition, are
admitted. Every source coordinate and quantity meaning remains in the reference
and endpoint receipts. Any
temperature dependence still belongs to the existing reference-mass enthalpy
path, even when sampled endpoint values agree. Held conditions set an initial
inventory on the fixed reference volume; stepping does not evolve density,
mass, composition or pressure-volume work. The `isotropic_conductivity` pairing
method uses the same card, external conditions and interval with the original conductivity
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
Uniform, per-element, nodal, material-ID and card-derived assignments therefore use
one preparation and source-domain policy. The theta path validates every
incident conductivity and capacity interval before and after prescribed
boundary lifting and at every accepted endpoint. Borrowed material source
objects retain the original receipts; a returned temperature field does not
itself own a new material card.

`HeatCapacity::Nodal`, `BackwardEuler::nodal` and
`assemble_nodal_capacitance` admit a positive P1 volumetric heat-capacity field
on the fixed reference mesh. They reuse the enthalpy owner's exact P1 row
moments: each tet contributes `V (sum(c) + c_a) / 20` to vertex `a`.
Consequently, stored heat integrates the correlation between a graded capacity
and a P1 temperature, which cell-averaged capacity discards. The diagonal is
assembled once, remains positive, and is reused by theta, backward-Euler and
history owners, including their nonlinear and radiative endpoint operations.
For constant specific heat this represents P1 reference density directly. If
both density and specific heat vary, the declared field is their P1 volumetric
capacity approximation; it does not exactly multiply two independent P1
fields. Values refer to the same supplied reference volume throughout the run;
current-state density cannot change closed-body mass. This explicit field
carries no material-source receipt and introduces no expansion or phase change.

`HeatCapacity::MaterialVolumeFractions` and
`BackwardEuler::material_volume_fractions` bind fixed P1 constituent volume
fractions to original `MaterialHeatCapacity` sources. This explicitly declares
a stationary, nonreacting composite at local thermal equilibrium with additive
reference volumes and no heat of mixing. Its capacity is
`C_vol = sum(phi_i rho_i Cp_i)`, consistent with constituent masses; separately
averaging density and Cp and multiplying them is not this physical model.
Fractions lie in `[0,1]` and sum to one at each vertex within `64*f64::EPSILON`,
without clipping or normalization. The resulting P1 capacity uses the same
exact nodal row moments, assembled once and reused by theta, backward-Euler,
history and duty-cycle owners. All constituent source bundles survive; their
temperature domains constrain every vertex of a tet where that constituent's
P1 fraction is nonzero. Entirely absent constituents add no domain restriction.
This conservative support check covers mixed interiors even when one corner's
fraction is zero. Conductivity remains independently declared, and no evolving
composition, expansion, reaction or latent heat is inferred. Regressions in
`tests/material_volume_fraction_capacity.rs` check mass-correct composite
heating, independent constituent-energy conservation, immutable preparation,
source retention and active-support admission.

`duty::march_duty_cycle_with_capacity` uses that same prepared storage for
scheduled volumetric power. Uniform, per-element, material-ID and source-bound
density/Cp capacities retain their distinct thermal inertia throughout the
cycle. The original conductivity and capacity domains admit every initial,
prescribed and accepted nodal temperature; an unsupported endpoint returns no
partial cycle. Existing uniform-capacity duty APIs delegate to this owner.
Boundary heat remains independent of the power schedule, while capacity,
conductivity and the linear preconditioner are prepared once. The existing
windowed audit continues to report scheduled source energy. Independent
regressions in `tests/material_heat_capacity/duty.rs` cover material-specific
diffusion eigenmodes under a ramp, conservative heat exchange across a
two-material interface and transactional source-domain refusals.

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
The steady solve scales by the fourth root of `εσA`, so representable results
do not require the coefficient or either `T⁴` to be representable separately.
Net radiation retains a factored temperature difference near equilibrium.
Zero drive preserves any admitted bath temperature exactly. Nonfinite derived
outputs and positive electrical power that underflows to zero return
`IncandescentRadiativeError::Unrepresentable`; browser bindings emit a typed
refusal rather than invalid JSON numbers. G0/G1/G3 tests cover these numerical
limits, a small net heat flow, and area/emissivity rescaling.

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

`radiation::OpenDiffuseEnclosure` admits reciprocal surface factors with row
sums below one. Each row deficit sees one explicitly black, isotropic ambient
reservoir. The completed reciprocal system adds a bath of effective exchange
area `sum(A_i F_i,bath)`, not a fabricated mesh surface; no escaped view is
redistributed among solid patches. The existing shifted radiosity kernel and
factorization evaluate both surface powers and ambient power independently.
Their signed sum remains a diagnostic, never a forced zero. Source emissivity
domains constrain solid temperatures, not the ambient temperature. The model
supports 1..=255 physical surfaces and refuses an absent escape, overfull rows,
overlapping/foreign traces, unsupported source temperatures and nonphysical or
unrepresentable arithmetic.

`radiation::solve_with_open_gray_diffuse_enclosure` drives the actual steady
solid solver on 1..=64 exact named uniform Robin/Neumann traces. It shares the
closed enclosure's bounded outer iteration and fresh full nonlinear residual
gate, retaining heterogeneous conductivity and fixed contact. Original
convection, original prescribed heat, net solid radiation and independently
computed ambient absorption remain separate. The original conduction problem
must have a thermal anchor before radiation is applied. Frozen emissivity,
isothermal patch means and a black diffuse bath are explicit approximations;
this adds no spatial visibility kernel, directional illumination,
participating-medium model or radiation-only nullspace solver.
`EnclosureRadiationLinearization::bind_open` and
`EnclosureRadiationGradient::pullback_open_with_reference_feedback` reuse the
closed thermal transpose with the exact completed open radiosity equations.
One absolute bath-temperature control includes all reflected paths; only
physical surfaces supply FEM traces, emissivity controls and state feedback.
Original Robin/Neumann laws, smooth heterogeneous conductivity, matching
contact and prescribed-temperature lifts retain their complete tangents.
The unchanged primal must pass its physical residual gate before derivatives
are returned. These are local discrete sensitivities, not error certificates.
`tests/radiation/open.rs` compares partial views and reflection with
independent scalar equations, checks equilibrium/near-equilibrium and cold-bath
limits, and solves an actual conduction/convection/radiation slab against an
independent heat-balance oracle. It also compares bath, emissivity and source
derivatives with complete nonlinear re-solves on Robin and prescribed-flux
traces, and checks that the open radiosity derivative tends to the closed
derivative as the escape vanishes.

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
triangle scan admits 1..=64 patches, 2..=32 replicates, up to one million
triangles, power-of-two samples and an explicit total intersection budget.
`project_open` instead retains sampled escape as exchange with an explicit
black surrounding reservoir and runs the same bounded reciprocal projection
on the augmented system. Its adjustment budget includes the ambient fractions.
Unassigned opaque and back-facing hits still refuse: they cannot be relabeled
as escape. The returned physical-surface factors and evidence feed
`OpenDiffuseEnclosure` directly. Geometry tests exercise an actual open cavity,
a single exposed patch, preserved blocking, replay and solid-plus-bath closure.

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
| `∫_e ∇λ_a·K∇λ_b dV` | `V · gᵀ_a K_bar_e g_b` | retained piecewise-linear `K(T_h)` in a fixed principal frame; its material moments integrate the entire P₁ temperature field, including crossed source knots |
| `∫_e λ_a f dV` | `Σ_b V(1+δ_ab)/20 · f_b` | `f` linear on the element (second-order otherwise) |
| `∫_F λ_a q_n dA` | `Σ_b A(1+δ_ab)/12 · q_b` | `q_n` linear on the face |
| `∫_F h λ_a λ_b dA` | `Σ_c h_c ∫_F λ_a λ_b λ_c dA` | `h` linear on the face; cubic simplex moments are integrated exactly |
| `∫_F h λ_a T_ref dA` | `Σ_b R_ab T_ref,b` with the same exact Robin block `R` | independently linear `h` and `T_ref` on the face |
| `∫_F (1/R'') λ_a λ_b dA` | `(1/R'') · A(1+δ_ab)/12`, staged as signed A/A, B/B, A/B, and B/A blocks | face-constant `R''` on exactly matching P1 traces |

The spatial Robin block is shared by steady, temperature-transient and
enthalpy residuals/tangents. `tests/robin_flux.rs` checks its matrix and load
against independent Duffy/Gauss integration, then verifies steady generation,
transient cooling and latent-plateau heat removal against that same independent
surface integral. This removes coefficient-averaging error for P₁ trace data;
it does not certify the supplied convection law or the spatial solution error.

`K_bar_e = (1/V_e) integral_e K(T_h) dV` integrates each retained
piecewise-linear conductivity component across all crossed material and nodal
temperature knots. The tetrahedral temperature distribution is a quadratic
B-spline; weighting it by a local hat function lambda_c gives a cubic B-spline
with that node's temperature repeated. Two-point Gauss integration on each
interval therefore integrates both `K_bar_e` and
`D_c = <lambda_c K'(T_h)>_e` exactly as polynomials, up to roundoff. Positive
scalar averaging preserves the fixed principal frame and constitutive positivity.
Constant tensors, uniform fields and single-segment tables retain fast paths.
No per-element heap allocation or new material query occurs; primal/flux calls
skip unused derivative moments. Cancellable integrations poll every 128 intervals.

The shared temperature Jacobian is
`J_ac^e = V g_a^T K_bar_e g_c + V g_a^T D_c grad(T_h)`. The enthalpy owner
applies its existing temperature-chart and declared mean-liquid-mass-fraction
scale derivatives to these same moments. Uniform/assigned flux recovery returns
`q_bar_e = -K_bar_e grad(T_h)`, the element-average flux of the discrete P1
field. Both steady and enthalpy conductivity-scale pullbacks reuse the integrated
stiffness. A nonuniform field crossing an interior material knot has a
differentiable material integral, even when its mean is exactly on that knot.
Uniform fields at a true slope kink retain the right-segment semismooth Newton
choice; the total-adjoint API refuses a unique two-sided material tangent there
and at sampled validity endpoints. Full nodal source-span admission remains
required. This is exact integration of the admitted P1/table model, not a claim
about unresolved continuum variation, interelement flux conservation or roundoff
bounds.

Seven focused crate regressions in `material/p1/tests.rs` cover independent
conductivity/derivative moments, constant and anisotropic limits, full source
support, consistent operator/flux/pullback results, and temperature/enthalpy
Jacobians against perturbed production residuals. Eight additional independent
native checks executed the actual scalar numerical leaf, including analytical
Beta moments, a multiknot oracle, all vertex permutations, adjacent floating-point
temperatures and interval cancellation. Full crate compilation and integration
execution remain separate validation requirements.

Material validity is checked over the whole P₁ field, independently of this
quadrature choice: every tetrahedron's four nodal temperatures must be finite
and inside its assigned conductivity model's span. These are the extrema of a
linear field on the tetrahedron. A valid element mean cannot authorize an
unsupported hot or cold corner. Primal assembly, temperature/enthalpy
Jacobians, and uniform/assigned heat-flux recovery use the same check. The
orthotropic span is the intersection of its principal conductivity components;
shared vertices must satisfy every incident material. Tests in
`tests/heterogeneous.rs` exercise both source endpoints and unsupported nodal
extrema whose means remain supported. Linear backward-Euler publication also
checks the new field against the conductivity spans; assembling a bounded
constant law at an admitted old temperature does not admit arbitrary heating.

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
| `ThermalResistanceTerm::slab_from_model` | reduces the original homogeneous conductivity tensor at an explicit reference temperature, retaining the model, every source receipt, normalized world-space normal and geometry. `NormalTemperatureGradient` gives `R=L/(A nᵀKn)` and permits tangential heat flux; `NormalHeatFlux` gives `R=(L/A)nᵀK⁻¹n` and permits a tangential temperature gradient. A diagonally equilibrated Cholesky solve avoids forming an inverse. The resulting term enters the existing series and finite-body consumers with uncertainty unstated. These are uniform-field reductions: rotated anisotropy with simultaneous isothermal ends and adiabatic sides, end/spreading effects and temperature-span evolution still require spatial or nonlinear treatment. Regressions distinguish the 25.5025-fold resistance difference for a rotated 100:1 tensor, compare both reductions with native spatial Fourier solves under their corresponding boundary conditions, and retain source/state/geometry refusals |
| `ThermalResistanceTerm::slab_over_temperature_span_from_card` | steady uniform one-dimensional slab with fixed face temperatures, no generation and no lateral losses: integrates every crossed segment of one continuously supported scalar conductivity source and computes `R=L/(A*mean_T(k))`. `source_span_heat_rate_w` returns the corresponding hotter-to-colder heat rate using only those retained temperatures. Original source answers, uncertainties, receipts, card and geometry remain in `BulkMaterialTemperatureSpan`. The integral is exact for the admitted piecewise-affine law up to floating-point rounding; source/interpolation/geometry uncertainty remains unstated for resistance. No point-only extrapolation, nonlinear multilayer interface solve, transient storage, isotropic inference or experimental qualification is claimed |
| `ComponentPower` / `PowerMap` | per-component dissipation in W bound to a mesh vertex set, plus the declared system total. Names are unique and sorted; bound vertices are sorted and deduplicated so a repeated binding cannot take a double share while the totals still balance |
| `PowerMap::volumetric_source` | projects the map onto the nodal `W/m³` source AND returns the audit computed from that same field, so the reported delivered power is what the solve receives rather than a parallel estimate of it |
| `PowerMap::regional_volumetric_source` | the s93ej.2 generation companion of `from_region_ids`: one finite W/m³ per region id, lumped to nodes so delivered power is exactly `Σ_e f_e V_e`. Unmapped regions and length mismatches refuse. This is generation, not conductivity |
| `SurfaceComponentPower` / `SurfacePowerMap` | dissipation bound to boundary-face SLOTS; `bind_boundary` attaches a separate uniform Neumann row `-P/A` to each physical footprint, preserving each heater's own power even at shared vertices. `neumann_flux` admits one nodal field only when shared vertices require the same flux; use `region_selector(mesh)` for its exact covered face set |
| `SurfacePowerAudit` / `SurfacePowerAuditRow` | the surface balance integrates each component's actual face flux; `face_area_m2` and the normalization accessor `lumped_area_m2` both refer to that component's own geometric area, never a neighbouring footprint's area |
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
| `MaterialId` / `MaterialTable` / `ElementMaterials` | a checked constitutive table plus one id per tet. `from_region_ids` maps labeled mesh regions onto that table. `bind_labeled_volume` is the mesher-agnostic s93ej.1 handoff: an independently audited `TetComplex` + region ids + region→material map become a `ConductionMesh` and assignment together. Assembly, Newton, `assemble_operator_with_element_materials` / `assemble_jacobian_with_element_materials`, and `element_heat_flux_assigned` use the named model integrated over the element P1 temperature field. This is constitutive selection, not `element_scale`. Unknown/duplicate/empty/unmapped-region/length-mismatched assignments refuse. `solve_with_element_materials` is the product entry; `solve_with_element_materials_and_interfaces` is the joint contact path. Region materials cannot erase undeclared coincident faces. A bound assignment's table receipts and provenance travel on `ConductionReport` in place of the unused fallback. The uniform `ConductionProblem::material` path is unchanged. Two- and three-layer series, assigned `k(T)` Newton, assigned Jacobian-vs-FD, assigned IFT vs FD, one- and two-material contact coexistence, assigned matdb receipts, unused-table/mixed-provenance honesty, rotated orthotropic x-series, contrast/refinement/coalescing/map-swap metamorphics, and a labeled-volume lowering through `bind_labeled_volume` are covered; that last case is not yet the parent solve-stage E2E |
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

Stop rule (`StopRule`), accepted only when the heat-residual criterion holds:

```text
  heat_target = residual_rtol · ‖A(T) T_anchor − b‖_free,2 + residual_atol
  ‖R(T)‖₂ ≤ max(heat_target, residual_roundoff_floor)
```

`T_anchor` takes the midpoint of the original prescribed temperature range on
free vertices, retaining the original fixed temperatures. Robin reservoirs do
not move that reference when Dirichlet anchors exist. Without fixed vertices,
the reference is `integral(h*T_ref)/integral(h)` over the Robin boundary, using
the same exact P1 face moments as assembly. Weak remote reservoirs therefore
cannot create an artificial temperature difference against stronger anchors.
The current iterate still owns every conductivity query. The reference therefore
measures source/flux loads and imposed temperature differences without making
the relative tolerance depend on a common absolute-temperature baseline or on
the caller's initial guess. Source-free unequal Dirichlet values retain a
nonzero physical driving scale. The reference is derived from the problem and
does not change the resumable state or legacy snapshot format.

The public residual and physical reports retain the assembled `A*T-b`. To avoid
chasing cancellation and temperature-representation roundoff near isothermal
states, each free row with `m` entries contributes
`gamma_(2m+2) * (|b_i| + sum_j |A_ij*T_j|)`, where
`gamma_n = n*EPSILON/(1-n*EPSILON)`; their Euclidean norm is the numerical floor.
`ConductionReport::residual_threshold` records the actual maximum used for
acceptance. This floating-point allowance is not an interval enclosure, an
assembly-error certificate, or a bound on spatial discretization error.

An accepted step with `‖ΔT‖∞ ≤ step_atol` and an unresolved heat residual returns
`ConductionError::NotConverged`, as does exhaustion of `max_iterations`. A small
temperature update alone never publishes a converged field. The legacy
`StopReason::StepTolerance` tag remains readable but is no longer emitted.
`tests/steady_heat_imbalance.rs` checks weak generated and prescribed heat at
large temperature baselines, strongly unequal Robin conductances with/without
Dirichlet anchoring, anisotropic isothermal equilibrium, and the small-step
refusal against independent slab solutions.

Globalization (`Nonlinearity`):

- `FixedPoint { relaxation ω, max_backtracks }` — damped Picard,
  `T ← T + ωδ` with `A_ff(T)δ = −R(T)`. The Krylov tolerance applies to the
  physical heat imbalance and the solved temperature increment. `A(T)` is SPD,
  so this path runs preconditioned CG.
- `Newton { line_search }` — Armijo sufficient decrease on `‖R‖₂`:
  accept when `‖R(T + αδ)‖ ≤ (1 − c α)‖R(T)‖`, shrink by `shrink` per
  backtrack, at most `max_backtracks`. The Jacobian carries the
  `V gᵀ_a <λ_c K'(T_h)>_e ∇T_h` term, which is NOT symmetric, so this path runs
  FGMRES.

Both share the **admissibility guard**: a trial iterate with any element vertex
outside that material's sampled temperature span has no merit value at all, so it
is rejected and the step shrunk. Extrapolating material data is never a fallback.

`InitialGuess::DirichletMean` projects the prescribed-temperature mean into each
free vertex's incident-material support. With no Dirichlet data it uses the
midpoint of that support, or the existing zero initial guess for unbounded
models. Per-element assignments completely replace the unused uniform fallback
for initialization as well as assembly. Empty support intersections at shared
vertices refuse by name. This is numerical initialization only: prescribed
temperatures and converged values are never clipped. The all-Robin two-layer
regression checks the resulting temperatures and heat flux against the exact
slab-plus-film series solution for two different unused fallback models.

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
13. **Surface power stays on its declared footprint.** Each component owns
    a separate uniform-flux Neumann row with `q_n = -P/A` over its own area.
    The native P1 face load integrates the constant exactly: each face adds
    `-q_n A_face/3` to its three vertices. Temperature remains continuous;
    prescribed flux may jump across a heater edge. A zero-watt neighbouring
    footprint receives zero direct heat, regardless of shared vertices.
    The single-nodal-field entry point refuses conflicting vertex fluxes
    rather than redistributing power while passing a total-only audit.
    Evidence: `tests/power.rs`, including per-footprint integrals and native
    assembly, steady solves, and a transient `P * elapsed` calorimeter.
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
    per-face rule the balance uses. The P₁ product is integrated exactly,
    including spatial correlation between the transfer coefficient and the
    temperature difference. The reported area means remain separate diagnostics;
    their product need not equal the heat rate. `EnergyBalance::robin_out_w`
    uses its own running accumulator, not a sum of region totals, so the whole
    and the sum of parts agree only up to
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
  and across refinement; adjacent footprints sharing vertices retain their
  individual watt inputs, including an off component; independent constant-
  flux P1 load moments and compatible equal-flux nodal fields; a native
  steady comparison and an insulated transient `P * elapsed` balance;
  the analytic cross-check that a full-face map reproduces the
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
  Crank-Nicolson and backward Euler delivering the declared ramp and dwell
  energy across aligned and unaligned boundaries; insulated-body calorimetry
  for a pulse missed by endpoint samples; mismatched full-load watts refusing; a
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
- RADIATION IS A SURFACE-ONLY FIRST RUNG. Linearized surface-to-ambient exchange,
  diffuse-gray closed enclosures and explicitly black-bath open enclosures are implemented. This crate
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
- CONTACT RETAINS INDEPENDENT P1 TRACES. The matching path requires exact
  coordinate-bit equality between duplicated boundary triangles. The explicit
  planar nonmatching path integrates their common refinement under declared
  geometry/coverage tolerances. Positive static resistance may vary per
  matching face pair or per nonmatching A-side face. The temperature-contact
  owner resolves declared mean or two-sided source laws at matching face means
  or nonmatching overlap quadrature points. A spatial map alone implies no
  evolving pressure, gap, temperature or contact-area law. Perfect contact and automatic scenario lowering are not
  inferred. A mapped surface still distinguishes material authority from
  caller-supplied values, as specified below.
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
  `LumpedNode::temperature_deviation_half_life_s` is `ln(2) C/(hA)` for
  deviation from the constant-power equilibrium, not for absolute temperature.
  The analytic response uses shared `fs-math::decay` arithmetic to retain
  small heating increments and scaled late tails. Focused G0/G1/G3 cases in
  `tests/lumped.rs` cover these limits, half-life ratios and interval composition.
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
- PRESCRIBED DUTY POWER IS INTEGRATED ANALYTICALLY OVER EVERY STEP.
  `DutyCycle::energy_scale_between` integrates each crossed constant or linear
  segment, including discontinuities and pulses shorter than a timestep. The
  resulting heat enters the thermal increment directly; theta still controls
  the discretization of heat transport, without changing input joules. The
  local integral avoids subtracting large cumulative energies, uses bounded
  segment lookup, and requires no per-step allocation or matrix rebuild.
  Unrepresentable intervals refuse instead of silently losing a pulse.
- THE SOURCE AUDIT IS TIED TO THE ASSEMBLED LOAD. The supplied full-load watts
  must agree with the spatial integral of the volumetric source (relative
  tolerance `1e-10` plus `256 * EPSILON` times its absolute nodal load sum).
  Declared energy uses the authored segment durations and declared watts;
  delivery uses the interval integrals and assembled watts. Neither is
  normalized to the other. The audit covers prescribed volumetric input,
  including heat at constrained nodes; boundary heat remains separate.
- EXACT INPUT ENERGY DOES NOT MAKE THE TEMPERATURE HISTORY EXACT.
  `WindowedEnergyAudit::steps_aligned` describes whether the solver resolved
  schedule boundaries, not whether prescribed source energy is conserved.
  Coarse implicit steps can smear a short burst's temperature response even
  when every joule is retained. `tests/duty.rs` verifies ramp, dwell and
  sub-step pulse energy with insulated-body calorimetry and boundary heat.
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
- SURFACE POWER BELONGS TO ITS PHYSICAL FACE SET. `SurfacePowerMap::bind_boundary`
  attaches separate named Neumann rows to faces of the original mesh left in
  the boundary's explicit adiabatic remainder. It refuses existing face/name
  ownership and partition-size mismatch; the caller retains the original
  mesh/boundary identity. Two components may not claim the same face. The
  optional `neumann_flux` representation requires equal prescribed fluxes at
  every shared vertex and refuses otherwise. Its region must still cover
  EXACTLY the map's faces; derive it from `region_selector(mesh)`. A map-wide
  nodal normalization that preserves only the summed watts is insufficient:
  it can heat a switched-off adjacent footprint.
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
  `m_i = sum_e integral(rho0 N_i)`, reducing to `sum_e rho0 V_e/4` for uniform
  density; equilibrium density changes never create or remove
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
  This lane binds one spatially uniform equilibrium chart and declared uniform
  or positive P1 nodal reference density. Prescribed endpoint temperatures
  select enthalpy by
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

  `from_nodal_reference_density` assembles each graded tetrahedron's exact
  row masses `V*(sum(rho)+rho_a)/20`. These preserve total mass and the integral
  of a P1 specific-enthalpy field, including density/enthalpy correlations that
  a cell mean discards. Normalized density accumulation avoids overflowing the
  raw density sum when the final reference mass is representable. Globally
  uniform input delegates to the original constructor with identical bits.
  The resulting immutable masses enter the existing nonlinear heat equation,
  latent storage, constrained-node reservoir heat, endpoint energy and
  history/source adjoints. Heterogeneous material-density pullbacks retain
  their original element/material controls; this constructor adds no nodal
  density-control derivative, material mixture law, evolving geometry or
  current-state density substitution. Native quadrature regressions check
  independent degree-two integration, graded total enthalpy, exact geometric
  partitioning, uniform bits and extreme representable density/volume scales.
  The full public regressions in `tests/enthalpy_transport/nodal_density.rs`
  cover nonuniform sensible heating, latent-phase evolution and adjoints,
  prescribed-bath storage, replay and admission. Their full spatial execution
  remains unverified in the current dependency-constrained environment.

  `HeterogeneousEnthalpyBackwardEuler::with_nodal_reference_density` binds the
  same exact P1 reference-mass integration to the existing heterogeneous owner.
  Original material assignments and borrowed caloric charts remain fixed;
  separate material traces retain their explicit contact convention. The
  positive nodal profile defines immutable storage for all existing primal,
  endpoint, radiation and history paths. It is a pre-evolution reference
  inventory declaration, not a density-at-temperature update that changes a
  closed specimen's mass. `reference_material` returns `None` for such a profile
  so scalar-only coupled source owners refuse an unsupported inventory binding.
  An already bound scalar phase-source inventory refuses profile replacement;
  the explicit `from_phase_materials_with_nodal_reference_density` constructor
  instead declares the profile when the original sources are first bound.
  The original element/material density controls remain additive constant
  perturbations about the frozen profile. Explicit nodal controls additionally
  provide a residual JVP and an optional nodal pullback using the same symmetric
  degree-two P1 mass action. Prescribed residual rows retain their original
  zero density derivative; all controls use the accepted endpoint's existing
  adjoint multiplier. Native kernel tests check independent quadrature, signed
  actions and transpose identities. Public integration regressions in
  `tests/heterogeneous_enthalpy/nodal_density.rs` cover graded sensible/latent
  storage, source-binding refusal, uniform parity and density derivatives;
  full thermal integration execution remains unverified in this environment.

  `EnthalpyHistoryRun` admits unique prescribed temperatures without changing
  initial physical enthalpy. Event-aligned, adaptive and ambient-radiation
  histories charge prescribed storage jumps to the bath once. Held boundary
  knots retain the exact left limit for the interval ending at a switch; the
  new bath acts on the next finite interval. Ambiguous latent temperatures
  and out-of-domain charts refuse without advancing state. Both iterated and
  coupled-Newton state-dependent-contact histories use that same admission and
  reservoir accounting. Contact transfer remains internal heat even when a
  prescribed bath heats or cools one of the bodies. The scheduled-contact
  regression in `tests/temperature_contact_enthalpy/surface_dirichlet.rs`
  compares event-aligned and explicit endpoints against an independent P1
  contact root, including bath storage jumps, replay and transactional refusal.

  `ScanningEnthalpyRun` admits the same unique prescribed temperatures for
  moving-source endpoints, event-aligned histories and adaptive half steps,
  including nonlinear conductivity and ambient radiation. The physical bath
  reaction accounts for both prescribed-node storage and absorbed source power.
  Prescribed equation rows are eliminated without dropping their source-basis
  columns from the consistent P1 mass action. Forward and reverse-history
  sensitivities retain free-state propagation, all seven heater controls and
  prescribed-temperature pullbacks; fixed caloric history has no endpoint-state
  influence. `tests/scanning_enthalpy.rs` checks independent latent-interior
  balances, held bath/source events, transactional refusal and adaptive reservoir
  sums. `tests/scanning_sensitivity.rs` checks those controls against finite
  differences with nonlinear conductivity, radiation and reverse propagation.

  Whole-scan objectives, point responses and smooth peak responses retain the
  original radiative reservoir and emissivity pullbacks. The shared reverse
  sweep contracts each endpoint reservoir seed into its authored constant,
  ramp, table or Chebyshev values on the same side of every event. One
  emissivity control per patch accumulates across the whole history, in
  original patch order. Full-tape and checkpointed evaluation share this path
  and its transpose work; neither solves radiation a second time. Objective
  and response scalar allowances include these returned controls. Value-only
  evaluation keeps the original physical signal admission, including valid
  chirps and constitutive corners, without demanding a derivative. Knot times,
  material claims and interpolation modes remain fixed. The focused
  `tests/scan_objective/radiation.rs` regressions compare whole-history and
  response derivatives with physical finite differences and check checkpoint
  identity, held-event sides and exact scalar allowances.

  The same shared reverse sweep returns authored prescribed-temperature row
  gradients through the original constrained-row multipliers. Boundary
  construction retains the actual winning region at each prescribed vertex;
  adjacent faces never multiply one bath contribution, and a lower-index static
  region still takes precedence over a scheduled region. Each selected region
  seed is contracted into its original uniform/constant, ramp, table or
  Chebyshev values with the original event side. Static partition data remain
  fixed inputs. Full-tape, point, peak and checkpointed results retain these
  gradients without another physical or transpose solve; scalar allowances
  include every returned bath control. Value-only paths retain valid physical
  signals independently of derivative support. `tests/scan_objective/baths.rs`
  checks coupled nonlinear-conductivity/radiation finite differences, replay
  identity, shared-vertex precedence, held-event causality and exact budgets.

  A scan can also retain the original `TemperatureContacts` owner through
  `with_contacts` or `with_contacts_newton`. Each moving-source trial updates
  live contact temperatures and, for coupled Newton, independently bound liquid
  fractions on both material traces. The same owner handles explicit endpoints,
  held-event sides, adaptive half steps and replay, retaining original contact
  states and the physical heat account in `ScanTransport::Contact`. Temperature
  contact iteration supports the existing area-mean ambient radiation. Explicit
  `with_pointwise_ambient_radiation` attachment instead composes the original
  local P1 emission with either contact policy, including phase-aware Newton.
  Natural convection/flux and shared prescribed-node reactions keep their own
  accounts; local radiation is neither shifted into a Robin surrogate nor
  counted twice. Coupled Newton still refuses area-mean radiation because that
  is a distinct physical convention. Fixed geometry, source claims and prescribed pressure remain
  unchanged. Endpoint adjoints bind the accepted complete contact residual and
  use one transpose for source, history, radiation and whole-law resistance
  multipliers. Whole-scan objectives, point responses and smooth peak responses
  also retain these complete-law controls in canonical interface order. The
  common reverse sweep sums the original endpoint contractions without another
  timestep factor or transpose; source identity, mapped status and trace count
  must remain unchanged. Full-tape and checkpointed paths use the same ordered
  accumulation. Scalar allowances count each original named interface once,
  including fixed interfaces and delegated pieces within an attached live model.
  Separate fixed-only interfaces remain fixed design inputs on this scan API.
  `tests/temperature_contact_enthalpy/scanning.rs` and
  `tests/phase_contact/scanning.rs` cover moving deposition, independent total
  energy and latent-contact balances, supported radiation, control derivatives,
  reverse propagation through melting and atomic refusal. Native execution of
  these additions remains unverified while the pinned dependency build fails.

  `tests/phase_contact/scanning/history_controls.rs` adds full-history contact-law
  finite differences for phase-aware matching, pointwise and delegated
  nonmatching traces; it checks point and peak limits, exact checkpoint identity,
  both held-event conventions and exact scalar admission. These tests are
  authored regressions, not claimed native passes in the blocked build lane.

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
- `EnthalpyBackwardEuler::adaptive_state` and its heterogeneous facade provide
  a full-window controller for directly assembled physical laws, alongside the
  existing authored-history adaptive stepper. Each trial
  compares one full backward-Euler step with two half steps from the same
  accepted history. Acceptance requires componentwise enthalpy, temperature
  and liquid-fraction discrepancies to meet their explicit tolerances; an
  isothermal latent plateau cannot conceal different stored energy. Relative
  scales use field changes, so a large absolute enthalpy reference does not
  loosen accuracy. The two-half-step state is accepted without extrapolation.
  Reference masses, nodal phase charts, temperature/phase conductivity,
  contacts and prescribed-temperature reservoir reactions retain the endpoint
  solver's physical meaning. Separate entry points retain the existing
  pointwise ambient and enclosure radiation laws.
  Only the two accepted fine endpoints contribute external heat; comparison
  solves and rejected trials never add heat. Whole-run stored energy is
  recomputed directly from reference masses and the original enthalpy field.
  Its reported balance defect is measured, not bounded by one step's energy
  tolerance. Each call lands exactly on its declared event time under fixed
  source and boundary laws, allowing explicit workload changes between calls.
  A bisectable trial that would leave exactly one floating-point clock tick
  before an event includes that tick before its coarse and fine solves when
  the extended duration exceeds the request by at most `8*EPSILON*request`.
  Coarsely resolved clocks retain the minimum-step refusal. Actual half-step
  times determine heat input and all acceptance checks.
  Rejection backoff uses the unextended requested duration, preserving the
  minimum-step stop. An already unbisectable requested interval still refuses
  without changing accepted time, state, heat or observer history.
  A trial budget or minimum-step stop returns a distinguishable partial
  trajectory. Cancellation and endpoint refusals preserve accepted physical
  history; continuing a bounded call retains the same controller state.
  Observed entry points publish a borrowed `EnthalpyAcceptedTrial` only after
  accepting a complete pair of fine half steps. The view includes both physical
  endpoints and their actual times, so consumers can track intermediate peaks
  without cloning fields or retaining a whole trajectory. Rejected comparison
  work never reaches the observer. `ControlFlow::Break(())` returns the distinct
  `ObserverStop` status after committing the accepted pair; continuation resumes
  after that pair without duplicate observations. Direct, pointwise ambient and
  pointwise enclosure radiation paths share this observation mechanism.
  This is an in-memory resumable trajectory with local time-error estimates,
  not a persisted checkpoint, accumulated-error certificate or timestep
  adjoint. `tests/adaptive_enthalpy.rs` checks independent spatial diffusion
  and latent-plateau solutions, analytic quartic radiative cooling,
  tighter-tolerance convergence, phase-fraction error control,
  event-aligned heating/cooling, accepted fine-step heat,
  rejected work, bitwise budget/resume replay and cancellation.
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
  from the production integrated assigned tensor, stiffness and
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
  Its shared bounded implicit solve requires both enthalpy interval refinement
  and the actual discrete energy residual. The constant-power allowance is
  `mass * enthalpy_tolerance_j_kg` joules; spectrum callers retain their explicit
  joule allowance. Exact roots, including adiabatic latent states, are preserved.
  At adjacent representable enthalpies only an endpoint within that allowance
  may be accepted; otherwise the step refuses rather than publishing an
  unclosed balance. These numerical tolerances do not bound time-discretization
  error. Nonfinite boundary/source power and residual arithmetic are refused.
  Disabled radiation does not evaluate an irrelevant fourth power, and active
  radiation factors the fourth-power difference to retain near-equilibrium flux.
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
  `LumpedThermalEnvironment::prepare_endpoint` exposes the same native
  backward-Euler storage and boundary powers to a simultaneous coupled solve.
  Preparation retains the original body, phase curve, transport and Biot
  admission once; `LumpedEnthalpyEndpoint::evaluate` queries actual represented
  wall enthalpy and endpoint fluid temperature without an inner thermal root.
  Internal heat remains an explicit complete source supplied by the caller.
  Sharp phase plateaus remain in the original enthalpy relation. The samples
  are unaccepted residuals, and the enclosing owner must close the opposite
  finite-store accounts. The original `advance` still supplies the accepted
  wall trajectory, followed by the enclosing owner's actual-endpoint recheck.
  This seam adds no wall heat-capacity approximation or independent source law.
  `LumpedEnthalpyEndpoint::solve_with_internal_power` now delegates an
  additional signed endpoint source to the original bounded enthalpy bisection.
  Its callback receives original h and temperature, preserving latent phase
  coordinates; prepared baseline heat remains included exactly once. The
  caller supplies a whole-chart positive-feedback bound in W/J of total body
  energy. `maximum_step_s_for_internal_power` combines that bound with native
  emissivity feedback before a solve, without replacing source data or choosing
  a temporal-accuracy tolerance. The returned native endpoint closes both the
  environment's joule budget and mass times its enthalpy tolerance. Cancellation,
  unrepresentable heat and chart escape refuse without publishing state.
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
  The …1720 tokens truncated… pullback, and verifies `J_s^T lambda - W^T wall_adjoint = weights`
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

### The lowered schedule preserves discontinuous switching and dwell energy

A scenario duty cycle declares dwells AT points, not transitions BETWEEN them,
so it lowers to piecewise-CONSTANT segments. No implicit smoothing: a step load
is real and common, and rounding it off changes the peak.

The source integral splits at every dwell boundary, including those inside a
thermal timestep. Each dwell contributes its duration times its power, for
both backward Euler and Crank-Nicolson transport. No invented ramp is needed
to recover the authored input energy. Tests cover aligned and unaligned
switches and compare a ramped profile with the same total energy.

This fixes the endpoint-quadrature deficit previously reported at switching
boundaries. It does not resolve a burst's continuous temperature peak inside a
coarse step: callers still choose timestep resolution for thermal dynamics.

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

`advance_contact_pair_with_environments` and
`MaterialContactLaw::advance_with_environments` additionally compose two
independent `LumpedThermalEnvironment` owners with that same finite contact.
Contact-covered area is removed from each body's convection/radiation area;
separation restores it. Bounded boundary-feedback trials reuse the original
contact solver and exact-enthalpy endpoint power evaluator. Fresh local and
combined heat balances, the environmental owners' residual allowances and a
conservative combined contact/environment Biot bound precede acceptance.
The original whole-surface environmental admission remains conservative with
partial coverage. Neither state is advanced in isolation or published during
iteration. Reports separate the internal contact transfer from actual ambient
and radiation heat, so the collection's external heat excludes contact Q.
This remains a fixed-geometry, uniform-body backward-Euler reduction with
explicit iteration limits and independent timestep-refinement requirements.

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

## Original heterogeneous endpoint residual for coupled feedback

`EnthalpyEndpointEvaluation::residual_j` exposes the original assembled heterogeneous endpoint residual in mesh-vertex order. It is moved from the existing assembly into the report, without a second operator or a second assembly. Prescribed-temperature rows are zero after their exact enthalpy checks; their physical reaction remains in `dirichlet_in_w`. Coupled owners can therefore close every free nodal energy balance instead of substituting a global scalar energy condition.

`HeterogeneousEnthalpyBackwardEuler::evaluate_endpoint_with_radiation` dispatches the existing `EnthalpyRadiation` choice to the original endpoint evaluation for no radiation, mean ambient/enclosure exchange or pointwise ambient/enclosure exchange. Existing radiation controls and source admission remain authoritative; this read-only dispatch performs no nonlinear solve or history advancement. The material-graph spatial integration tests exercise consistent P1 heating, heterogeneous charts and prescribed bath reactions. New native test execution is pending the authoring environment's unchanged runtime dependency/resource blocker.

## Source-bound spatial sharp melting and phase-specific conductivity

transient::enthalpy::material::melting::MaterialMeltingTransport resolves one
original ResolvedMeltingCurve, independent solid/liquid conductivity tables
and a single-phase reference-density answer from the same material card.
Pressure, composition, phase selectors and other fixed source coordinates
remain those admitted by the original caloric owner. Each branch's Cp, density
and conductivity requires explicit temperature support. Conductivity queries
never extend the solid source into the liquid interval or vice versa. All
original caloric, conductivity and reference-density receipts remain available;
the combined identity includes the explicit coexistence closure.

`MaterialMeltingTransport::from_isotropic_card_exact` consumes
`ExactMaterialMeltingSource` to retain only the original Cp/density boundaries
and fusion plateau. Each sensible temperature, inverse and tangent still uses
the original affine Cp integral. It does not interpolate T linearly between
those compact public knots. The separate absolute fusion-rounding allowance
admits only the represented plateau width at the selected enthalpy datum;
it is not a sensible interpolation or source-uncertainty bound.

Both ingresses share the same conductivity preparation, source-domain checks,
phase-volume closure and frozen reference inventory. The resulting owner
feeds the original uniform/regional, scheduled, adaptive, BE/SDIRK2 and adjoint
paths without a separate field operator. Original conductivity knots and their
classical-adjoint refusals remain present even when the caloric chart is compact.
This removes artificial caloric subdivisions, not physical source transitions,
spatial discretization error or material uncertainty.

Four focused regression definitions in `tests/material_melting/exact.rs`
check independent affine-Cp inverses and conductivity tangents, compact versus
refined source representations, mixed-phase spatial residual/Jacobian columns,
native BE/SDIRK2 heating and freezing with fixed reference mass, and retained
source/corner refusals. These new definitions passed syntax parsing but were
not natively executed in the resource-exhausted implementation environment.

At coexistence the caller chooses a parallel arithmetic or series harmonic
conductivity mixture in **volume** fraction. The liquid volume fraction is
computed from the original liquid mass fraction and the two sourced densities
at the melting point. Neither mixture is an inferred microstructure or a
material-validation claim. Outside coexistence, the corresponding source k(T)
and the admitted caloric T(h) determine k(h). The exact one-sided derivative
includes phase-conductivity variation on an unsmoothed latent plateau.

The adapter prepares the existing EnthalpyBackwardEuler owner with its
original immutable reference masses. It interpolates nodal k(h) with P1 basis
functions and supplies the exact element mean and dk/dh to the existing sparse
assembly. A retained unit tensor is only the numerical basis for this assembly,
not an additional conductivity claim. Replacing that tensor, adding element
material overrides or stacking a second phase multiplier refuses; physical k
cannot be multiplied twice silently. No duplicate solver or apparent-Cp
approximation is introduced. Existing boundary reactions, contact, radiation,
adaptive stepping, energy gates and cancellation keep their original owners.
Existing endpoint adjoints use the complete conductivity chain rule and refuse
additional source k(T) corners as well as the original caloric corners.

The six original regressions in tests/material_melting.rs cover independent
mixed-phase P1 balances and Jacobian columns; history-adjoint finite
differences; original lumped-owner uniform heating/freezing with unchanged
reference mass; a boundary-driven computed phase contour whose independent
semi-discrete equation is dh/dt=800+0.5h, with heat closure and temporal
refinement; validity/budget/cancellation refusals; and a conductivity-source
corner refusal. The contour case is a single-tetrahedron temporal check, not
a continuum Stefan or mesh-convergence claim. The existing Stefan refinement
suite remains the spatial-front reference. Focused actual-source native
verification passes all 12 uniform/regional/scheduled material-melting cases and
15 adjacent enthalpy-adjoint regressions. This is not a full workspace/DSR verdict.

These source bindings use isotropic materials over a stationary
mesh. Evolving geometry, density-driven mass changes, pressure work, convection,
phase kinetics, anisotropic phase mixtures and uncertainty propagation are
not supplied by this path.

### Regional original-source sharp melting

`HeterogeneousEnthalpyBackwardEuler::from_melting_materials` binds a retained
`MaterialMeltingTransport` per material and an explicit material ID per element.
Each region keeps its own fusion temperature, latent heat, enthalpy datum,
solid/liquid Cp and conductivity sources, phase densities and coexistence
closure. Immutable vertex assignments share the existing heterogeneous storage
owner; neither sources nor charts are copied into per-element models.
Distinct material IDs require disjoint nodal traces and declared interfaces.
Conflicting shared vertices, unknown IDs, orphan vertices and exhausted spatial
budgets use the existing preparation refusals rather than averaging charts.

Reference masses use each material's original single-phase density answer.
`from_melting_materials_with_nodal_reference_density` instead admits an explicit
positive P1 reference-density profile through the existing exact mass assembly.
That declared initial inventory supersedes scalar source density; phase density
still converts mass fraction to volume fraction and never changes stored mass
during heating or melting. `melting_transport_source` exposes each retained
source independently of this inventory choice.

`melting_transport_problem` retains the original unit numerical basis and
refuses extra conductivity tensors or phase multipliers. The existing nodal
k(h) interpolation, exact P1 element mean and enthalpy derivative assemble the
same physical operator for BE, SDIRK2, adaptive stepping and strict discrete adjoints.
Each nodal source sees only its assigned material's h. Uniform
`MaterialMeltingTransport::element_heat_flux` and assigned
`melting_transport_heat_flux` use the same source-specific nodal k(h) and T(h).
These are constitutive element averages, not a conservative face-flux recovery;
contact heat remains in the existing interface report. Genuine source
or caloric corners retain the classical-adjoint refusal.

The regional coupon tests use distinct fusion temperatures, latent heats,
reference densities, coexistence closures and disjoint enthalpy datums. They
independently reconstruct conduction/contact balance, element heat flux and
Jacobian columns, compare source-adjoint perturbations, and check both time
integrators. Latent exchange with graded reference inventory, uniform-source
parity and original-assignment/source-domain refusals use the same test seam.

### Original melting sources on scheduled thermal histories

`EnthalpyHistoryRun::from_melting_materials` and its
`from_melting_materials_with_nodal_reference_density` variant retain the same
original source assignment through the existing thermal timeline. They consume
`PhaseEnthalpyHistoryProblem` boundary, contact and volumetric-source schedules;
a one-entry material table also covers a uniform melting specimen. Initial
enthalpy is admitted against each original caloric and conductivity law before
the trajectory is returned. Reference inventory is prepared once and remains
unchanged across heat pulses, cooling, phase change and continuation.

The existing history problem dispatch distinguishes declared conductivity,
prepared phase-source transport and original sharp-melting transport. Sharp
melting keeps its P1 interpolation of nodal k(h), while prepared phase sources
keep their integral of k(h_h); neither discretization silently substitutes for
the other. All existing fixed/adaptive BE, SDIRK2, contact and radiation history
wrappers use that same selected storage and physical problem.

Event times, boundary/source quadrature, phase-aware time-error controls,
accepted-prefix state, cancellation and step budgets remain with the existing
history owners. A refused future source or cancelled advance cannot reset
accepted enthalpy, material identity or the physical clock. This ingress adds
no integrator, latent-energy account, phase-kinetic law or moving-domain model;
each wrapper retains its existing supported combinations and limitations.

G1 tests integrate held heating/cooling through sensible and latent states,
independently check two-material contact balances at source events for both
integrators, and reconstruct nodal radiative heat from the original boundary
areas. G4 checks retain accepted enthalpy and clocks through refused future
heating/bath values and cancellation, then compare resumed evolution with
replay. Six existing prepared-phase history cases also pass, guarding the
other original source family through the same dispatch.

## Interference-driven linear lumped thermal response

The native lumped rung now accepts an explicitly bound prepared quadratic
power spectrum through `lumped::interference::PreparedLumpedInterference`.
It borrows the actual `LumpedNetwork`, original nodes and shared
`fs-fft` spectra. Preparation applies the original declared `BiotGate`,
then requires exactly one named binding per node, nonblank source units and
a common declared source clock, and a finite positive watts-per-spectrum-unit
scale. For a spatially uniform W/m³ source, that scale can be an explicitly
declared volume; no spatial integral or absorption fraction is inferred.

Each independent node retains its original constant C, G and common ambient.
For an absolute source-clock interval, the shared causal spectral integral
solves C dT/dt = P(t) - G(T - T_ambient) analytically, preserving all authored
difference and sum frequencies without a time grid, carrier averaging or
network re-solve. A general quadratic observable is not proof of dissipation:
the caller explicitly binds signed applied thermal power. Reversible field
storage, acoustic intensity and mechanical work do not become heat merely
because they have a spectrum. An ensemble-mean power spectrum produces the
linear model's mean temperature, with no invented fluctuation statistics.

The atomic result retains original node/network/gate/clock provenance,
absolute endpoint and DD temperature change, unweighted signed source input
joules, stored energy change and outward signed ambient exchange. The source
integral, causal convolution and initial ambient relaxation are used once
each; observation does not consume energy or advance external thermal state.
DD initial-temperature queries retain weak accumulated heating across successive
intervals, while the shared one-minus-decay kernel avoids subtracting a rounded
exponential from unity. The balance residual is an algebraic diagnostic, not
an independent error certificate.

Preparation bounds nodes, combined prepared/result vector payload and total
per-query spectral work (counting reused spectra again for distinct thermal
rates), with a hard ceiling of 1,048,576 terms. Queries allocate one result
vector; their per-node causal mathematics allocates nothing. Cancellation and
any invalid unit/clock/name/scale/domain/finite-range/budget admission discard
partial results. Ambient, initial and endpoint Kelvin temperatures must be
finite and positive. Endpoint admission does not prove positivity or a material
validity envelope throughout the intervening trajectory. Nonlinear capacities,
temperature-dependent conductances, radiation and exchange between nodes are
outside this constant-property rung.

The focused tests in `src/lumped/interference/tests.rs` compare constant
forcing with the original native response; full real carrier/beat forcing and
proper optical-carrier beats with an independent RK4 ODE; interval composition
and energy signs; weak DD temperature accumulation and tiny relaxation; and
original Biot, bindings, physical domains, work/storage budgets and cancellation.

The exact shared implementation, the original native LumpedNode/LumpedNetwork/
BiotGate source and the thirteen new thermal tests compiled, linked and ran
with real fs_math: all 48 combined tests passed, with no compiler warnings.
Only the outer ConductionError variants were nominal. Full Cargo workspace,
Clippy and the repository quality gate were not run in this scoped executable.

## Interference-driven native nonlinear enthalpy evolution

The child module `lumped::interference::enthalpy` advances one actual
`LumpedEnthalpyBody` under an explicitly bound prescribed thermal-power
spectrum. Original mass, equilibrium phase curve, latent fractions,
temperature-dependent transport receipts, fluid/enclosure temperatures and
convection-plus-radiation Biot gate remain authoritative. Its backward-Euler
residual calls the original phase and boundary-power helpers and the same
bounded native bisection implementation. Both callers retain their original
physical residual; the shared solver refines until each caller's actual joule
budget passes instead of treating a narrow enthalpy interval as energy closure.

Every accepted substep consumes the full shared spectrum's analytic interval
integral in DD joules, with an explicit positive W-per-input-unit scale and
declared source clock. Unit conversion acts on source amplitudes BEFORE time
integration, so a tiny authored source in large physical units is not first
erased by an unscaled integral. The consumer takes the converted integral
directly as source joules; it does not divide through a rounded mean watt value.
A spectrum already in W uses scale one. The absolute clock adds original start
and accepted elapsed time in DD, so late short substeps are not reset or erased.
All DC power belongs in the bound spectrum; the original config's separate
internal-power field must be zero.

Source resolution and thermal numerical accuracy are separate controls. The
shared outward-rounded primitive-ripple bound limits each step's departure
from interval-mean forcing, including cycles invisible to endpoint or
step-doubling checks. Positive unit scaling rounds that bound outward again.
One admitted maximum-step bound covers all shorter steps without repeating a
spectral traversal. Fast carriers can remain analytically integrated when
their primitive ripple is below the declared joule budget. The original
maximum thermal step still controls a first-order nonlinear discretization;
root/energy residual and source-ripple tolerances do not certify its global
temperature error.

The atomic trajectory retains every actual accepted native phase state and
its source, convective, radiative, storage and residual energy. Boundary
energies use original endpoint fluxes times accepted duration, with positive
sign inward. Total storage is also computed directly from original and final
enthalpy coordinates. Source energy is counted once; observations do not debit
a field, mutate a material card or publish external thermal state. Exact dark
equilibria are retained without injecting bisection rounding. Input below the
original f64 phase-state resolution stays visible in source joules and the
declared numerical residual, rather than being relabeled as darkness.
Attempts, accepted steps, spectral visits, output samples and vector storage
are bounded; cancellation, unsupported phase/transport domains or an unclosed
energy step discard the complete partial march.

A spectrum explicitly represents applied thermal power, which may include
signed extraction; reversible field storage is never reclassified as heat.
Nonlinear evolution driven by ensemble-mean power is a trajectory under that
mean forcing, not an ensemble-mean temperature or fluctuation prediction.
Electromagnetic material feedback requires an explicit rebuild by the original
field owner. The existing restrictions on uniform-body validity, source
homogenization and equilibrium phase kinetics remain.

Focused tests cover the original latent plateau under beating, an alias case
with zero full/half-cycle input but a real melting excursion, admissible fast
carrier primitives, original constant-power bisection, sourced varying
emissivity/radiation with thermal-step refinement, signed heat accounting,
low-component source energy, destructive darkness, a large absolute clock,
physical-unit restoration and fidelity/domain/metadata/resource/residual/
cancellation refusals.

Native Rust compilation, linking and execution passed all 13 new thermal
regressions and all 49 shared temporal/spectral/causal regressions: 62 passed,
0 failed. The focused executable uses the exact current equilibrium phase
chart, body/transport, conductivity table, emissivity and implicit helper
source, with real fs_math, fs_matdb, fs_evidence, fs_qty and fs_blake3 libraries.
Only the unrelated crate ConductionError envelope is nominal. No constitutive
substitute or alternate production thermal solver is used. This is not a full
Cargo workspace, Clippy or DSR quality-gate run. Native rustfmt passed.

## Scheduled material pointwise radiation on primary enthalpy

`MaterialRadiativeEnthalpyHistoryRun` retains the original
`EnthalpyHistoryRun` and adds original material ambient emissivity with an
independent black-reservoir temperature schedule per trace. Construction pins
the selected claims, checks exact caloric-object identity against storage and
admits initial nodal source support. Every nonlinear endpoint delegates to the
existing heterogeneous primary-h material radiation owner. Evolving epsilon(T)
or epsilon(f), pointwise epsilon*T^4 integration, source domains and physical
energy gates remain owned by that solver; no second phase or heat solver is
introduced.

Direct endpoints retain right-continuous schedule sampling. Event-aligned
steps merge reservoir, boundary and source knots and use interval-end left
limits so a short held pulse is not skipped or backdated. Adaptive proposals
reuse the original h/T/liquid-fraction step-doubling controller. The second
half starts from the first half's complete h state, with its own refreshed
emissivity. Both accepted half-step reports survive; all accepted fields and
the clock publish atomically through the original run. Refusal or cancellation
preserves the complete accepted state and source bindings.

Each sampled boundary needs one prepared material radiation owner; Newton
reuses it within the interval. Preparation-time emissivity answers are exposed
as references only; accepted reports contain the actual nodal source states.
The model is stationary P1 conduction with black ambient view factor one. It
does not add enclosure histories, spectral radiation, a finite bath, melt
flow, remeshing, source uncertainty propagation or a global time-error bound.

Four focused regression cases were added: closed-form latent heating/cooling
under epsilon(f), enthalpy-reference invariance, held-pulse event timing and
signed energy closure, phase-only adaptive rejection, atomic rollback/resume
and exact chart identity. All four changed Rust files pass parser checks and
tracked diffs pass whitespace checks. Native conduction compilation/execution
remains unverified under the existing memory/disk constraint; these cases are
included by the existing thermal-history CI target, not claimed as passed.

## Phase-continuous source clocks across enthalpy feedback intervals

The additive `solve_lumped_enthalpy_spectrum_from_clock_dd[_cancellable]`
entrypoints accept a canonical DD absolute origin. Finite components are
normalized without rounding the origin through f64; unrepresentable origins
refuse before source work. That explicit origin is authoritative and replaces
the supplied controls' legacy `start_s`. The returned controls retain its
rounded compatibility view, while `source_start_s_dd()` retains the complete
clock. The original f64 entrypoints delegate with a zero low component.

A consumer that refreshes the original field/material owner can continue from
the returned origin plus the final accepted DD elapsed time. Its next interval
therefore retains carrier/beat phase even when a short interval lies below an
f64 epoch's spacing. The consuming `into_samples()` moves the actual accepted
thermal history out of a march without a clone or a borrow of a temporary
source spectrum. Energy, work and other report getters remain available before
consumption. These APIs do not add or imply a constitutive feedback law,
transient polarization-history transfer or a nonlinear temperature-error bound.

Three focused regressions cover repeated accepted intervals at a large epoch
with a significant fractional-clock contribution, exact accepted-state/history
composition and source-energy closure, zero-copy history extraction, canonical
origin admission/cancellation, and unchanged native results through the
existing f64 API. All three new regressions and the thirteen existing nonlinear-enthalpy tests
compiled and ran successfully after recovery: 16 passed, 0 failed. The focused
executable uses the exact native body, phase chart, conductivity/emissivity,
transport and bisection helpers, with the actual fs_math, fs_matdb, fs_evidence,
fs_qty and fs_blake3 libraries. Only the unrelated ConductionError envelope is
nominal. The 49 shared temporal/spectral/causal tests also passed. Pinned native
rustfmt and independent source review passed. This is not a full Cargo/DSR gate.

## Original material enclosure controls on the complete enthalpy adjoint

Material pointwise closed/open enclosures now have explicit accepted-state
adjoint binders for uniform storage, heterogeneous reference masses and
material temperature contacts. The retained physical h tangent includes source
emissivity feedback and the current reflected field. Original scalar claims
and curve ordinates keep their material/claim identity, axis and fixed knot
abscissae. One shared bath-temperature control exists only for open enclosures.
History, source, reference-density, contact-resistance and radiation controls
compose through the original complete thermal transpose.

For K=I-diag(1-mean(epsilon))*F, a parameter contraction uses one K-transpose
solve with the retained current factorization. The right-hand side is F^T times
weighted patch absorption. Local source derivatives integrate
Na*Nb*(sigma*T_h^4-G_patch) and include the reflected correction divided by patch
area; the open bath contributes its original black-row derivative. Centered
radiance is retained to avoid cancellation near equilibrium. No nodal dense
Jacobian or separate reflection solver is introduced.

Binding rechecks the unchanged accepted physical report and original source
receipts, admits declared workspace before parameter retention, and polls
cancellation. The old fixed-emissivity pullback refuses material controls.
Geometry, view factors, caloric charts and source knot abscissae remain fixed.
Moving source corners and caloric transitions refuse typed two-sided admission.
A source corner at a coordinate constant on a smooth open h interval is admitted:
latent temperature and pure-phase fraction have zero composed h slope there.
Ambient and enclosure controls share that same source-composition rule.

Three native regression cases were added for closed/open latent source, bath,
history, density and heat controls; exact nonuniform triangle moments; datum
invariance; budgets, receipts and cancellation; and typed kink refusals. The
existing ambient constant-fraction case was updated. All twelve touched Rust
files pass parser checks; native conduction compilation/execution is not claimed
under the current resource constraint. An independent exact-polynomial-moment
calculation on unequal-area open patches agrees with finite differences within
7.1e-8 for nodal emissivity controls and 1.1e-10 W/K for bath controls. This checks
the reflection formula, not execution of the native spatial owner. Both material
radiation test targets are included in the existing thermal CI lane.


## Interference heating in the original spatial enthalpy residual

`transient::enthalpy::interference::PreparedSpatialInterference` binds unique
named spectra to explicit vertices of one original `EnthalpyBackwardEuler`.
The authored field is P1 volumetric thermal-power density. Each binding declares
its input units and a positive conversion to W/m3; unbound vertices are exactly
zero. The prepared owner retains the actual storage and mesh objects, so equal
labels, topology or node counts cannot authorize a foreign geometry.

For each requested interval, the shared DD-clock causal integral applies the
unit conversion before time integration. The existing P1 rule assembles
`V*(sum(E)+E_i)/20` as dual nodal joules. The physical source account separately
integrates `V*sum(E)/4`; it is not inferred from a Newton residual or from a
point density times an invented nodal volume. The source enters the original
nonlinear balance directly as `m*delta_h + dt*(A*T-b) - source_j`, before the
solver's final f64 residual interface. There is no J-to-average-W-to-J round trip.
The original Newton/Krylov driver, exact transport tangent, phase chart, reference
mass and material-domain gates remain authoritative.

The owned endpoint reports actual h, temperature and liquid fraction, its DD
source load, local residuals and independent DD source/storage/boundary accounts.
Prescribed-temperature reaction is the complete constrained-node physical
balance, including the negative spectral load and any initial boundary jump.
Internal contact heat is reported by the original owner and excluded from
external input. Latent heat is already in h. All volumetric input, including DC,
must belong to the explicit spectra; a nonzero static source refuses to prevent
a second credit. The result exposes no fictitious instantaneous or mean source
power. Both free-node and total physical joule defects must pass the declared
endpoint energy tolerance before any result is returned.

The maximum thermal step and source-ripple allowance are separate controls.
Ripple admission uses an outward-rounded positive spatial sum of the shared
local primitive-deviation bounds. Opposite nodal sources and complete carrier
cycles therefore cannot hide unresolved spatial forcing behind a zero total
source integral. A fast carrier can be admitted when its physical primitive
ripple is already below the requested allowance. This is not a certified bound
on nonlinear thermal time error, continuous-time peaks or spatial interpolation
error. The caller refines time and mesh independently and may chain accepted h
and DD end clocks without resetting source phase.

Original Neumann/Robin boundaries, fixed thermal interfaces and all four original
pointwise radiation owners retain their mesh, boundary and caloric identities.
Temperature-dependent contact-feedback and automatic field/material feedback are
not introduced by this ingress. A prescribed expected power does not establish
an ensemble-mean nonlinear temperature. Source spectra may be signed; a general
quadratic observable is not proof of dissipative heat.

Preparation bounds bindings, retained spectral terms and additional source
workspace. Each interval makes two bounded spectral passes and linear mesh
quadratures; native spatial/Newton budgets still apply separately. Cancellation,
domain failure, ripple refusal or nonconvergence leaves input history and source
owners unchanged. Nonzero unrepresentable products refuse instead of becoming
false darkness. This is a bounded MR19/MR06 capability, not closure of their
broader coupled-physics acceptance.

Eight scoped arithmetic tests compile and pass against the exact production DD
residual, P1 source quadrature and outward ripple helpers with actual fs_math
and shared fs_fft. The probe supplies a fixed prepared-mesh/Cx interface; it
does not exercise or substitute a Newton or thermal solver. Twelve separate
original-owner integration regressions cover P1 reactions, nonlinear spatial
conduction, latent heating, cancelling local carriers, DD continuation, tiny
unit conversions, original radiation/Robin and refusal paths. Those full FEM
tests were added but not executed in this recovered environment. Production
and test sources pass pinned native rustfmt parsing, and independent source
review covered the native source and radiation lifetime paths. Full Cargo,
Clippy and DSR verification remains unperformed.

## Original-source conductivity throughout enthalpy histories

`EnthalpyHistoryRun::from_phase_materials` and
`from_phase_materials_with_nodal_reference_density` retain the original
`MaterialPhaseTransport` sources, exact caloric borrows and admitted material
assignments. `PhaseEnthalpyHistoryProblem` supplies boundary, volumetric-source
and interface histories without an unrelated conductivity tensor or additional
phase multiplier. Initial h must satisfy both the original caloric chart and
the admitted conductivity support.

Plain, event-aligned and adaptive endpoints, fixed/material ambient radiation,
material contact, scanning and accepted scan adjoints obtain their physical
problem from the same retained storage owner. Original k(h_h) is integrated on
the actual P1 enthalpy field at every Newton candidate, including interior
source-knot crossings. Frozen k, averaged phase and averaged nodal k are not
substituted for that source. Reference masses remain immutable, whether their
declared inventory is scalar or P1. Refusals preserve accepted h/T/f, time and
step allowance.

Three focused G0/G1 regressions compare an independent nonlinear P1 latent
moment equation, two enthalpy datums, scalar/P1 reference inventories, held
source events, phase-only adaptive rejection and source/cancellation rollback.
The actual updated library compiled natively against the matching cached
dependency rlibs; all 16 thermal-history tests passed, including those three.
The original material ambient and enclosure targets also passed 9 and 7 tests.
A broken nested test-module path was corrected so the four earlier material
ambient-history regressions actually execute. DSR and the required sibling
checkouts were unavailable: these are focused native fallback results, not a
full workspace/DSR verdict, a continuum Stefan bound or experimental validation.
The stationary source interpretation (explicitly isotropic or fixed-frame
orthotropic) and contact/radiation domain limits remain authoritative.

## Material enclosure radiation throughout enthalpy histories

`MaterialEnclosureEnthalpyHistoryRun` composes the original primary-h history
with source-dependent gray reflection in a closed or open enclosure. The
original pointwise owner supplies mesh, ordered named traces, reciprocal view
factors and watt tolerances. The wrapper pins original emissivity claims and
checks exact caloric-object identity against every participating storage node.
A closed enclosure takes no reservoir schedule; an open enclosure requires
one explicitly dimensioned black-bath temperature history.

Each sampled boundary must preserve the complete original named face
partition and natural radiation traces. Immutable face geometry and optical
factors are shared across candidates instead of rebuilding surface quadrature
or refactoring the fixed snapshot. The material owner still refreshes actual
nodal epsilon(T/f) and the current reflected operator at every Newton state.
Original-source k(h_h) histories compose through the same retained transport
dispatcher, including nonuniform bulk temperature gradients.

Direct endpoints use right-continuous values. Event steps merge source,
boundary and bath knots and use interval-end left limits. Adaptive proposals
reuse the h/T/f full-step/two-half-step controller; the second half retains
the first half's complete caloric state. Both half-step reports survive and
publish atomically with the accepted endpoint count and clock. Invalid source,
bath, budget or cancellation leaves all accepted fields unchanged.

Closed reflection remains internal exchange. Open-bath radiation, Robin
convection and volumetric generation retain separate signed heat accounts and
enter the existing reference-mass energy check exactly once. This is a fixed
geometry gray-patch irradiation model with material pointwise emission, not
spectral transport, participating gas, moving visibility or a global time-error
bound.

Five new G0/G1/G3/G4 history regressions cover independent two-patch latent
endpoints, large enthalpy-datum changes, distinct material temperatures, merged
bath/source/Robin events, phase-only adaptive rejection, source identity,
failure/cancellation/resume and simultaneous original-source conductivity with
active spatial transport. The actual updated library compiled natively against
matching cached dependencies. The enclosure target passed 11 tests, then its
new twelfth source-composition test passed through a focused invocation.
The source-history and material ambient targets also passed 16 and 9 tests.
Two earlier enclosure-test import paths were corrected; all three earlier
material enclosure adjoint regressions now execute successfully. DSR and
required sibling checkouts were unavailable, so this is focused native
fallback evidence, not full workspace or experimental validation.

## Phase-aware contact throughout primary-enthalpy histories

`ContactEnthalpyHistoryRun::new_newton` retains original A/B caloric objects
and phase-aware interface sources through direct, event-aligned and adaptive
histories. Every candidate evaluates independent P1 T(h)/f(h) traces and the
existing complete contact Jacobian. Latent fractions remain active where
dT/dh is zero. Matching, pointwise matching and nonmatching quadrature,
original conductivity and immutable reference masses keep their owners.

Prescribed-temperature chart admission and the independent physical nodal
residual/global joule gates remain authoritative. Coupled controls are retained
and must match each advance; a mismatched policy refuses. The existing
temperature-only outer-iteration constructor and its mean-radiation attachment
remain available. Policy inspection distinguishes those methods explicitly:
`coupling()` returns the optional outer policy and `newton_config()` returns
the optional full-h Newton policy. Initial attachment and failed candidates
publish no physical state.

Three new G0/G1 regressions execute independent latent quadratic exchange
across all three contact discretizations, complete melting/freezing through
energy-predicted sensible equilibration, held heater pulses, phase-only
adaptive rejection, exact chart/source identity, source support, cancellation,
policy refusal and deterministic resume. A malformed partial-support fixture
was corrected so its authored surface samples lie inside their declared
validity, preserving the same affine resistance law.

The actual updated library and test targets compiled against matching cached
native dependencies. All three new history tests and all three existing
thermal-contact history tests passed. The larger phase-contact target retains
an unrelated all-fixed-contact fixture that the current temperature-contact
constructor refuses; this is not a claim of an all-green broad target.
DSR and required sibling checkouts were unavailable. These are focused native
implementation checks, not measured-interface validation or a continuum
error bound.

## Phase-dependent contact with material enclosure histories

`ContactEnthalpyHistoryRun::with_material_enclosure` composes the retained
enthalpy Newton contact policy with the original closed/open pointwise material
radiation owner. The existing coupled Newton system simultaneously resolves
reference-mass storage, source conductivity, both bodies' contact T/f states,
phase-dependent emissivity and gray reflection. Contact resistance and
emissivity are reevaluated from the same candidate h field. The original
combined endpoint audit checks the full nodal joule residual and global heat
balance before any accepted state or clock is changed.

A shared immutable enclosure binding retains the original mesh, named traces,
view factors, source contexts, pinned emissivity claims and exact caloric
objects. It also serves the standalone material enclosure history, preserving
that public API. Open enclosures require one typed black-bath temperature
history; closed enclosures have no external radiative reservoir. Its knots
join source and boundary events, including coincident held-signal changes.
Direct endpoints sample right limits; event-aligned and adaptive candidates
use interval-end left limits. Rejected candidates add no heat or accepted
steps. Adaptive phase/enthalpy control commits only the two actual half-step
states and their original physical contact/radiation reports.

Internal contact and closed reflection cancel in the global balance. Open-bath
exchange, authored heat sources and original boundary heat retain separate,
once-only accounts. Attachment requires the Newton policy, an unadvanced
history and exclusive radiation ownership. Source-domain, bath, identity,
work-budget and cancellation refusals preserve the accepted h/T/f and clock.

Three new native regressions compare eight-node implicit enthalpy trajectories
with independently solved two-body contact and gray reflection equations,
including exact triangular emissivity moments, unequal reference densities,
isothermal phase evolution, closed energy conservation, scheduled bath/heater
changes, phase-only time-error rejection and continuation after a future bath
refusal. All six phase-contact history tests and the twelve existing material
enclosure radiation tests pass against the composed production library.

## Second-order spatial primary-enthalpy stages

`transient::enthalpy::sdirk` adds `advance_sdirk2` to the original uniform and
heterogeneous enthalpy owners. Two `EnthalpySdirk2Input` records carry physical
laws sampled at t+gamma*dt and t+dt, with gamma=1-1/sqrt(2). This is Alexander's
two-stage stiffly accurate SDIRK2 method (SIAM J. Numer. Anal.14(6),1977,
1006–1021; https://doi.org/10.1137/0714068). It reuses the existing reference
mass, spatial assembly, exact material/contact/radiation tangents and bounded
Newton/FGMRES solver. Source k(h_h), heterogeneous caloric charts, fixed and
source-dependent contact, prescribed temperatures and original pointwise or
material ambient/enclosure radiation participate at both actual stage states.
Mean-temperature outer radiation iteration is explicitly outside this API.

Both stages retain the real accepted old h. The first stage's authentic net
nodal heat creates a dual-joule carry for the second residual. That carry is
never queried as a fictitious caloric state or credited as external heating.
The second solve starts from the real first-stage h. A valid stiff step may
therefore succeed even when the usual algebraic effective-old-h expression
lies outside the original source support. All actual stages and endpoints
must still satisfy their original material domains and physical joule gates.

`EnthalpySdirk2Stage` retains instantaneous physical powers, original material
contact/radiation reports, prescribed-bath reactions and nonlinear work.
`EnthalpySdirk2Step` separately integrates the tableau-weighted external heat
and reconstructs original-to-final reference-mass storage. Independent global
and nodal balances include the actual constrained-node reactions, with internal
contact/conduction and closed reflection excluded from external heat. A failed
second stage returns no accepted interval or fabricated backward-Euler receipt.

Accuracy is second order in smooth admitted regimes. Phase/source corners and
forcing jumps need event resolution/refinement. L-stability alone supplies no
unconditional positivity, phase monotonicity or entropy theorem; the original
BE adjoint is not a second-order trajectory adjoint. Focused native regressions
check an independent spatial eigenmode and its second-order refinement, the
exact SDIRK stability function, valid real stages with an out-of-domain
algebraic carry, enthalpy datum shifts, cancellation/refusal/replay, and changing
prescribed baths with independently integrated source and reservoir heat.

## Scheduled second-order enthalpy with contact and material radiation

`into_sdirk2` converts an existing `EnthalpyHistoryRun`, Newton
`ContactEnthalpyHistoryRun`, `MaterialRadiativeEnthalpyHistoryRun`, or
`MaterialEnclosureEnthalpyHistoryRun` into `Sdirk2EnthalpyHistoryRun` without resetting its accepted h/T/f, clock, reference
mass, source conductivity, pinned contact/emissivity claims or lifetime budget.
Each interval samples the original boundary/source/bath schedules at both
SDIRK2 tableau times. The original spatial stage kernel closes source k,
phase-dependent contact and local material ambient emission or reflection
at their actual stage states. Only the stiffly accurate endpoint publishes; stages are not separate
accepted intervals and their carry is never a physical heat source.

Direct advancement uses right-continuous stage samples. Event advancement
ends at the next source, boundary or open-bath knot and uses its exact final-
stage left limit. Coincident heater/bath switches participate once. Both stage
owners retain their own sampled boundary, source and reservoir; one endpoint
sample cannot silently stand in for both stages. A nonrepresentable interior
stage time refuses without altering physical history.

Adaptive advancement reuses the existing phase-aware interval controller:
one full SDIRK2 interval is compared with two actual half intervals in h, T and
liquid fraction. The full discrepancy and conservative unequal-half inflation
are retained rather than assuming smooth second-order asymptotics at material
corners. No Richardson-extrapolated state is accepted. The actual fine pair's
integrated heat reports remain separate, two interval slots are reserved, and
publication is atomic after the final physical/error/cancellation checks.
`endpoint_solves` counts SDIRK2 interval attempts, each containing two implicit
stage solves. Native nonlinear nonconvergence and endpoint residual/energy
refusals now bisect within the same attempt/step limits. Source-domain,
material, budget and cancellation refusals propagate immediately; exhausted
physical recovery retains its original typed refusal.

`Sdirk2EnthalpyHistoryRun::with_moving_sources` adds immutable normalized
absorbed-power footprints alongside the original background heat, phase
contacts and material ambient/enclosure radiation. Both actual stages and
adaptive half intervals sample their own source position and power. The native
`source_w` and `source_j` accounts include the complete field exactly once.
Backward-Euler scanning and SDIRK2 share source admission, superposition and
footprint-motion limiting through the original moving-source owner; the scan
continues to retain its per-source deposition observations. Path and power
knots join the original boundary/reservoir events, and optional motion limits
include zero-power travel. Direct overlong intervals refuse; event/adaptive
intervals shorten while retaining the original phase, energy and entropy gates.

An accepted `ScanningEnthalpyRun` can transfer directly through `into_sdirk2`
with its source declarations, motion policy, caloric state, clock and lifetime
budget intact. Original Newton contact laws transfer with their policy.
Explicitly selected fixed pointwise ambient radiation transfers its original
emissivity declarations and reservoir histories as well. Both SDIRK stages
prepare their local radiation against their own sampled boundary and bath;
direct `with_pointwise_ambient` attachment uses the same original admission.
Area-mean scan radiation and iterated-only contact policies refuse because the
native stage owner cannot preserve those conventions; no physics is silently
dropped or replaced. Moving geometry, optical absorption and automatic SDIRK2
objective/checkpoint tapes are outside this integration. Focused regressions in
`tests/phase_contact/history/enclosure/scanning_sdirk.rs` use independent P1
source loads, latent-contact equations and gray-reflection equations to check
stage sampling, held events, motion-limited adaptive entropy, continuation and
transactional refusal. Their native execution remains unverified: the pinned
`asupersync` build was killed before `fs-conduction` compiled.

The history `sdirk::sensitivity` module binds explicit or event-limited intervals
through this same stage sampler and one original physical solve. The accepted
h/T/f, clock and count publish only after both native stage tangents, moving
source bindings and observation metadata admit. Original source declarations
are shared immutably across retained interval bindings. Each source and authored
Dirichlet history contracts its two actual stage cotangents at the recorded
time and side; common authored controls count once, with no extra tableau weight.
Temperature and liquid-fraction observations use the original per-node charts,
and the complete previous-h seed remains available for reverse history.

Optional parameter-aware history pullbacks reuse the same two native transposes
for contact resistance and fixed/material radiation controls alongside sampled
or authored source/bath controls. The full native parameter representation and
both stage reports survive. Radiation parameters remain explicit per-stage
controls; authored radiative schedule maps are not inferred. Selected time grids,
event/adaptive decisions, geometry and caloric charts remain fixed. Five focused
regressions in `tests/phase_contact/history/enclosure/scanning_sdirk/sensitivity.rs`
cover independent source/contact/reflection equations, exact held-source and bath
coefficients, prior-state derivatives, wrapper identity and atomic refusal.
Pinned parsing and bounded source review passed; native execution is unverified.

Five new native history tests compare original P1 conductivity moments and
source-clock quadrature with independent scalar stages; preserve graded
reference masses and enthalpy datum invariance; resolve held sources and
closed/open contact-reflection events; reject temperature-invisible phase
error; and resume after second-stage source, future-bath or cancellation
failure. The full thermal-history target passes 19 tests and the focused
contact-history group passes 8, alongside 16 spatial enthalpy transport tests.

## Discrete adjoints of the actual SDIRK2 stages

`linearize_sdirk2` binds the actual two-stage enthalpy solve for uniform or
heterogeneous reference-mass storage. The primal retains its distinct SDIRK2
stage and whole-interval heat accounts. Sparse stage Jacobians and bounded
transpose solves share the existing enthalpy adjoint operator; no fictitious
previous enthalpy or backward-Euler acceptance receipt is constructed.

Let a=gamma*dt and r=(1-gamma)/gamma. At the accepted first stage, including
its prescribed-node reaction, a*Q1=M*(h1-h0). The second free-row equation is
therefore M*h2-(1-r)*M*h0-r*M*h1-a*Q2=0. The reverse sweep first solves the
second physical transpose and carries r*M times its free-row multiplier into
the first-stage seed. The old-state pullback includes both the first-stage
free-row storage multiplier and (1-r) times the second-stage one. Distinct
prescribed masks remain stage-local, including nodes released between stages.

`pullback_stages` accepts objective cotangents on both physical stage states
and returns the previous enthalpy gradient plus separate stage nodal-source
and prescribed-temperature gradients. Temperature objectives use each stage's
original dT/dh, so a latent plateau does not acquire a fictitious temperature
slope. The retained stage linearizations can be swept backward over a bounded
caller-owned trajectory without retaining Newton or Krylov iteration history.
Each transpose retains its independently recomputed residual and work count.

`pullback_with_parameters` and `pullback_stages_with_parameters` retain those
same two thermal transposes while exposing original physical controls at each
stage. Whole-law contact resistance multipliers use the retained native fixed
or evolving interface operator. Fixed gray/spectral radiation and material
ambient/enclosure emissivity controls keep their original scalar, curve or
surface coordinates and source identities. Each stage contracts a*lambda on
its own free rows, reversing that sign for outward radiation; prescribed-row
multipliers never enter physical heat controls. The native reflected-radiation
contraction is reused without another thermal solve or extra tableau factor.
Contact trace work and combined radiation workspace have separate explicit
allowances. Independent stage owners stay separate; shared declarations may
sum their corresponding stage contributions without another time weight.

`tests/material_ambient_radiation/adjoint/sdirk_parameters.rs` and
`tests/phase_contact/history/enclosure/sdirk_parameters.rs` exercise independent
latent-stage equations, original emissivity ordinates and reservoirs, changing
prescribed masks, coupled contact/reflection, energy invariance, retained fixed
interfaces and bounded refusal. These regressions are authored but remain
unexecuted under the unchanged native dependency build blocker.

These are derivatives of the admitted discrete solve on a fixed time grid,
fixed reference geometry and retained material/source selection. Chart or material
corners require the existing smoothness admission; the method does not invent
a derivative at a refused kink. Automatic history tapes, differentiation of
adaptive/event decisions, geometry and caloric-chart parameter gradients are
outside this interface.

## Bind original lowered boundaries to scenario histories

`ThermalBoundaryHistory::from_boundary` accepts the original mesh and its
already-lowered thermal boundary, retaining face assignments, row kinds and
shared-vertex Dirichlet precedence. It reconstructs the partition through the
existing boundary builder and then uses the original history admission and
sampling owner. Geometry selectors are not rerun. Callers must keep the
original mesh and boundary together; the face-count check is not a geometric
identity proof.

The native cooling consumer uses this seam for prescribed temperature,
signed heat flux and independent Robin reference-temperature histories.
Sampling rebuilds the original prescribed-node lift. Existing capacity and
total-enthalpy solvers retain their old physical state and account prescribed
jumps through their existing boundary-node storage reactions. This interface
adds no boundary-aware trajectory adjoint or temporal accuracy claim.


### Prescribed baths with source-backed enthalpy contact fixed points

`TemperatureContacts::advance_enthalpy` now admits prescribed-temperature
boundaries through the original heterogeneous enthalpy owner. Its initial
contact query uses that owner's boundary-compatible endpoint candidate; every
inner transport response and fresh physical audit still receives the original
accepted `old_h`. The prescribed-node storage jump therefore contributes to the
reservoir reaction exactly once. Ambiguous latent-plateau bath temperatures,
source-domain errors, work exhaustion and cancellation still refuse without
publishing a changed caloric history. The existing ambient-radiation fixed-point
entry shares the same boundary initialization; no new radiation law is introduced.

The G1/G3 surface-contact regression compares positive and negative bath jumps
against an independent eight-node temperature system with analytic reference
mass, tetrahedral stiffness and consistent triangular contact blocks. Separate
regressions retain the independent spatial-melting roots and verify latent-bath,
work-limit and cancellation refusals. These prove the bounded discrete model;
they do not establish arbitrary source-surface convergence or continuum error.

## Pointwise radiation inside source-backed enthalpy contact fixed points

`TemperatureContacts::advance_enthalpy_with_pointwise_ambient_radiation` and
`advance_enthalpy_with_pointwise_enclosure_radiation` compose the existing
source-contact fixed point with the original local ambient and closed/open
enclosure owners. Each contact response passes the same accepted old specific
enthalpy to heterogeneous enthalpy Newton. Radiation is solved at that response's
actual local field, without replacing simplex emission by mean temperature.

A fresh source query at the candidate produces the physical contact operator.
Acceptance requires both its complete nodal joule residual and the independent
storage/heat balance. Prescribed-temperature storage enters the reservoir
reaction exactly once. Internal contact and closed reflection cancel globally;
ambient or explicit open-bath radiation remains separately accounted external
heat. Global cancellation cannot excuse a failed local residual.

Response, Newton and Krylov counts retain all discarded contact trials.
Pointwise radiation contributes to the inner Newton work, with zero separate
outer radiation iterations; retained Newton history belongs to the final inner
solve. Owner mismatch, late source-domain exit, cancellation and exhausted work
return no accepted endpoint and preserve the caller's old history.

Four native regressions use independent simplex emission/radiosity and
manufactured transport equations to check local latent states, datum shifts,
closed/open heat, both prescribed-bath jump directions and refusal/retry.
The interface retains fixed reference mass/geometry and the prepared owner's
emissivities/view factors. It adds neither a surface-contact tangent nor a
convergence guarantee for arbitrary resistance surfaces.


## Entropy admission of the original nonlinear P1 enthalpy step

Uniform and assigned heterogeneous storage expose `advance_with_entropy`
and `evaluate_endpoint_with_entropy`. The first performs the original
backward-Euler Newton solve, physical energy checks and an additional caller-
toleranced entropy admission before returning a state. The second measures a
supplied endpoint without accepting it. Both retain the original full physical
operator before prescribed-row elimination; the Newton tangent is not a heat
transport operator.

The account integrates each original caloric chart's dh/T against frozen
reference nodal mass. Sensible films contribute their original caloric entropy
at the same host temperature; latent heat is already contained in the host
chart and is never deposited again. Dirichlet reactions include fixed-node
storage at the actual time-quadrature bath temperatures. Each active Robin
triangle retains its positive P1 bath field B(x), transfer coefficient h(x)
and body trace T(x). Bath entropy is the actual integral
-integ[h(x)*(T(x)-B(x))/B(x) dA]. The same native consistent face heat rows
q_i give signed production integ[h*(T-B)/B dA]-sum_i(q_i/T_i).
Heat-flux and volumetric inputs remain explicitly unclassified; dividing their
watts by receiving nodal temperature does not invent a reservoir or convert
work into heat.

Bulk, sheet and contact production use the same sparse transport entries as
energy. The original consistent Robin matrix contributes its actual heat and
signed entropy separately. Positive constitutive tensors alone do not imply
nonnegative nodal finite-element entropy production. Admission refuses negative
bulk/contact, Robin or admitted radiation production, negative thermal storage minus inputs, or
unclosed entropy beyond the caller's absolute J/K allowance. It never clips
matrix entries, lumps Robin faces, or repairs heat from the discrepancy.

Exact state entropy, physical transport, time-discretization defect, free-row
algebraic residual and remaining reduction closure are separate quantities.
The account adds one linear pass over vertices, sparse entries and boundary
faces, plus the original caloric queries; Newton iterations perform no entropy
work. Spatial bath integration uses at most two median-isotherm subtriangles,
exact polynomial transverse moments and closed logarithmic reciprocal moments.
Near equal baths use a positive series with a bounded omitted analytic tail;
all loops are bounded and cancellation-aware. No mean temperature, face
refinement, altered heat matrix or extra quadrature state is introduced.
An algebraic reciprocal-reference correction retains the near-constant limit
without subtracting two first-order entropy totals. These are f64 evaluations
of the analytic integrals, not certified floating-point enclosures.
Photon entropy requires its own owner. Acceptance is a check of this admitted
discrete state, not a whole-system second-law certificate or a continuum,
mesh-monotonicity or time-accuracy theorem.

G1/G3 regressions compare spatial baths and nonuniform transfer coefficients
against independent positive Duffy/Gauss integration, all face-vertex
permutations, repeated bath values, small near-equilibrium production and an
independent logarithmic reference over a 1e12 bath contrast. Accepted native
BE and SDIRK2 steps retain their original energy trajectories and actual
stage-weighted bath entropy; nonpositive bath temperatures,
original material limits and cancellation retain transactional refusal.


## Entropy admission on actual second-order spatial stages

`advance_sdirk2_with_entropy` is available on uniform and assigned heterogeneous
enthalpy storage. It reuses the original two nonlinear SDIRK2 stages, their
final physical operators, source-dependent contacts and original film calorics.
The prior-stage dual-joule carry stays inside the stage equation; it is never
added as a heat source or evaluated as an extrapolated material history.

Every physical source, bath, bulk/contact, Robin and admitted radiation entropy contribution uses
its actual tableau temperature and the same positive weight as energy. Exact
material and film entropy storage is counted once, from original to final state.
The native double-double stage loop supplies prescribed-bath reactions and free
residuals, including subtraction of prior heat on constrained nodes.

With tau=gamma*dt and r=(1-gamma)*dt/tau, the signed temporal defect is
deltaS - r*sum(deltaE1/T1) - sum((deltaE2-carry)/T2). The separate algebraic
account is r*sum(residual1/T1)+sum(residual2/T2), over the actual free rows.
These are tableau-stage terms, not final residual/Tfinal. The temporal defect
can be negative even during ordinary physical cooling and is never clipped or
converted to material heating. The combined identity is independently closed.

The optional gate checks each weighted stage's signed bulk/contact, Robin and
admitted radiation production before checking the whole interval's
storage-minus-input and closure. An energy-conserving coarse
prescribed-temperature history can fail this gate; no state or partial interval
is returned on refusal. Pointwise ambient and reflected-enclosure radiation
admit the closed material/field/black-bath accounts below; mean-patch entropy
remains unsupported. Photon entropy requires its separate original owner. The existing
energy-only API takes the same implementation with entropy disabled. Additional
work consists of scalar accumulation in the existing stage loop and one entropy
pass per converged stage, with no additional spatial assembly or Newton solve.
Five focused regression definitions cover an independent uniform Robin ODE,
latent and film storage, enthalpy-datum shifts, heterogeneous reference masses,
negative discrete spatial production, and carry-dependent prescribed-bath refusal.


## Closed-reservoir entropy of native spatial ambient radiation

Uniform and assigned heterogeneous storage expose
`advance_with_radiation_and_entropy` and
`evaluate_endpoint_with_radiation_and_entropy`, using the existing
`EnthalpyRadiation` selection. `None`, `PointwiseAmbient` and
`MaterialAmbient` retain the ordinary native residual and nonlinear solver.
Source-dependent emissivity keeps the exact original material and caloric
objects, source support, local temperature and latent-phase response.
The existing SDIRK2 entropy method accepts these same pointwise ambient
owners at each authentic tableau time.

For native outward nodal heat q_i, original patch heat Q_p and positive
black-reservoir temperature B_p, radiative entropy input is
-dt*sum_p(Q_p/B_p). Signed discrete material/field/bath production is
dt*(sum_p(Q_p/B_p)-sum_i(q_i/T_i)). Every reservoir retains its own B_p.
The actual final native nodal load and per-patch report are borrowed directly
from the energy owner; emissivity queries and face quadrature are not repeated
for entropy. Radiative q_i enters each physical row and prescribed-temperature
reaction once. Stage weights, free-row residuals and prior-stage carry remain
those of the native energy equation.

This model encloses the stationary radiation field and the material/black-bath
exchange, with no net radiation-field storage. Its reservoir transfer and total
exchange generation are not entropy carried by photons. Total emissivity does
not determine a photon spectrum or gray occupation law, and no spectral
entropy is inferred here. The separate photon owner must retain its explicit
physical assumptions; its emitted entropy must not be added again to this
closed heat-reservoir balance. Mean-temperature outer radiation iterations
remain unsupported. Native reflected enclosures use their own independent
external bath row in the account described below.

Positive radiative material coefficients do not make the P1 receiving-node
entropy test positive for every nonuniform face. Signed production is retained
and checked for each positively weighted SDIRK stage and the whole interval.
Double-double bath and nodal reductions reduce near-equilibrium cancellation,
but the native heat integrations retain their own floating-point roundoff.
There is no unconditional relative-accuracy or positivity claim at vanishing
production. Caloric state entropy, radiative and conductive transport,
temporal defect, algebraic residual and final closure stay separate.

`thermal_entropy_change_j_k(cx, old_h, new_h)` provides the same exact original
bulk/film state difference without a solve or transport assembly. It retains
reference masses, original caloric support and latent integrals, validates both
history shapes and polls cancellation. It supports accepted interval owners
that must measure original-to-final storage rather than sum rounded subtotals.

Focused regression definitions cover radiative heating/cooling without changing
the native energy trajectory; independent simplex moments for multiple baths,
nonuniform faces and fixed-node reactions; authentic SDIRK quadrature; signed
negative finite-element production; near-equilibrium heat exchange; original
material emissivity feedback on a latent plateau and caloric datum shifts; and
unchanged source/chart/cancellation refusal. They do not assert a photon-flux
or continuum entropy certificate.

## Entropy across the accepted adaptive enthalpy history

`adaptive_state_with_entropy` prepares a resumable
`AdaptiveEnthalpyEntropyState` using the shared backward-Euler and SDIRK2
adaptive controller. The physical history is exposed read-only through `history()`;
`entropy_balance()` retains the complete accepted window. The original
energy-only methods continue through the same controller with the optional
account disabled. No separate Newton solve, spatial discretization, time-error
estimator or controller is introduced.

Each comparison trial keeps its coarse solve energy-only and entropy-admits
the two actual fine intervals through their original physical operators.
Backward Euler uses its endpoint account; `advance_to_sdirk2` and its observed
counterpart retain both native stages and their actual heat/reservoir
quadrature. The methods may alternate between event calls without resetting
accepted enthalpy, heat or entropy.
Only a pair accepted by the native time-error controller enters the retained
account. Rejected comparison work contributes to work counters, never to
storage or heat/entropy totals. Ten double-double accumulators retain signed
reservoir transfers, conductive/Robin/radiative production, unclassified
source and flux terms, temporal defects and algebraic residuals across
resumption. Material and bonded-film storage is remeasured directly from the
original initial enthalpies to the accepted final field. The unclassified-input
flag survives cancellation of signed source totals.

One caller-supplied absolute J/K tolerance applies to each fine interval and
to the whole accepted window; it is not multiplied by the number of steps.
The next account is checked before the physical field, clock and entropy are
published together. A refused pair leaves earlier accepted pairs intact.
Observer cancellation after publication retains that complete pair, matching
the existing controller. Original recoverable physical endpoint refusals keep
their bounded step backoff; entropy refusals propagate without a new retry
policy. Source laws may change between calls while the original caloric owners
and fixed reference masses remain bound to the history.

`advance_to_with_radiation` and its observed counterpart admit the existing
fixed and material-state pointwise ambient or enclosure radiation owners with
their actual internal-reflection and external black-bath accounts. Mean-patch
radiation remains unsupported. No photon entropy is inferred. Regression definitions compare
native accepted trajectory bits, budget-resumed totals, rejected-work
exclusion, latent storage and signed source cancellation, observer
cancellation, and radiative cooling with the original native energy path.
The SDIRK2 definitions additionally compare original stage bath entropy,
method switching through latent storage, whole-window resumption and native
negative-production refusal. L-stability is not an entropy-positivity proof.


### Entropy admission for scheduled adaptive enthalpy histories

`EnthalpyHistoryRun::advance_adaptive_step_with_entropy` and
`Sdirk2EnthalpyHistoryRun::advance_adaptive_step_with_entropy` use the original
event-limited, phase-aware adaptive controller. The coarse solve remains an
energy-only error comparator. Each of the two actual fine intervals uses its
native entropy finish, and only the accepted fine pair enters the returned
`EntropyAdaptiveHistoryStep`. Its physical report preserves the original
accepted trajectory and work counters; `half_entropy` retains both native
accounts; `entropy` combines their exchanges and signed defects with storage
queried directly between the original and final caloric states.

The declared absolute J/K allowance applies to each fine interval and once to
the complete accepted pair. It is not doubled for the two intervals. The pair
must pass before fields, time or accepted-step counters are published.
Cancellation, material-domain, nonlinear, energy and entropy refusals preserve
the previously accepted history. Native nonlinear nonconvergence and endpoint
residual/energy failures share the original controller's bounded bisection;
material/source admission, cancellation, explicit work-budget and entropy
refusals remain immediate. All retries retain event boundaries, phase
coordinates, step floors and attempt limits. Discarded trials contribute no
entropy to the result.

Scheduled SDIRK2 retains its original two-stage quadrature, including genuine
stage temperatures and interval-end left limits for held sources and boundary
schedules. Its internal stages and stage carry are not independently accepted
physical intervals or external heat. The shared stage sampler supplies the
same native source, boundary, material and contact declarations to either
finish mode. The scheduled material reflected-enclosure binding retains its
original radiation at every sampled stage, including actual open-bath
temperatures. Unsupported mean-patch entropy models retain their refusal.

These scheduled reports cover one accepted pair per call; they do not create a
persistent lifetime entropy account or multiply a lifetime allowance by the
number of calls. Entropy admission is separate from the local temporal-error
estimator and does not provide a global discretization or phase-transition-time
certificate.

## Entropy admission of native reflected spatial radiation

The original pointwise enclosure variants, `PointwiseEnclosure` and
`MaterialEnclosure`, now participate in the native backward-Euler,
SDIRK2 and supplied-endpoint entropy paths on uniform and heterogeneous
enthalpy storage. The same shared finish admits the existing scheduled SDIRK2
material-enclosure binding. Resumable fixed-law adaptive histories retain the
original energy-only coarse comparison and entropy-admitted fine endpoints.
The original caloric objects, fixed reference masses, latent heat, films,
material-emissivity queries, view factors and reflected heat operator remain
unchanged.

For the native nodal heat q_i oriented out of material, a closed enclosure has
zero external radiative entropy input and signed generation
-dt*sum_i(q_i/T_i). Reflection is internal to this material/field account.
An open enclosure retains the independently evaluated black-bath row:
Q_b is positive OUT OF THE BATH, so entropy input is dt*Q_b/B and generation
is -dt*(sum_i(q_i/T_i)+Q_b/B). B is the actual positive bath temperature
in that step or tableau stage. Physical surface heat is never counted again
as an independent reservoir transfer.

The bath row is not replaced by minus the physical surface sum. Their signed
`system_energy_closure_w` remains separately available in the original
radiation report under the original watt tolerance. No entropy correction,
heat redistribution, coefficient clipping or fictitious material-temperature
bath is created to force closure or production positivity. The measured P1
nodal entropy production remains signed and uses the ordinary admission gate,
including each positively weighted SDIRK2 stage.

The entropy finish borrows the converged nodal heat and complete native report.
It performs no extra radiosity solve, emissivity query or radiation quadrature.
Source-dependent local emissivity retains its original temperature and phase
state, and finite-window spectral heat retains its original model/truncation
controls. This is a stationary radiation-field and thermal-reservoir account;
it supplies neither spectral photon entropy nor an omitted-band entropy bound.
Mean-temperature radiation approximations remain unsupported.

Focused regression definitions compare closed exchange with independent
parallel-patch equations, open exchange with actual stage bath temperatures,
and material/latent endpoint accounts with independent simplex and radiosity
moments. Native energy trajectories, entropy signs, original storage, adaptive
accepted-pair parity and independent bath heat are retained.

## Reciprocal transmitting thermal sheets

`radiation::sheet::ReciprocalSheetOptics` admits source-bearing, passive
two-sided diffuse optical partitions. Each side may have different absorptance
and reflectance; one reciprocal transmittance couples the two sides. Both
faces must share one temperature and the same exterior radiometric / refractive
index convention. Consumers bind equal physical face areas. Kirchhoff emission
uses each side's own absorptance; visible alpha or gray emissivity does not
silently become spectral data.

Gray and finite-band evaluations report incident, absorbed, reflected,
transmitted, emitted and outgoing power per area. Material heat is explicitly
absorbed minus emitted. A single face's net incoming-minus-outgoing optical
power includes transmission and must not be deposited into storage. The two
thermal face powers belong to one physical sheet inventory.

The finite-band path uses the original `BlackbodyBand` emission, fixed-endpoint
analytic tangent and omitted-spectrum fractions without renormalization.
Consumers retain their original spectral coverage admission. Invalid optical
partitions and temperatures outside the declared source interval refuse.
These scalar laws describe diffuse, isothermal sheets; coherent films,
specular/refraction geometry and through-thickness temperature gradients
remain separate models.

## Layer-resolved absorbing and emitting material

`radiation::slab::AbsorbingSlab` prepares the hemispherical response of a
front-to-back stack of sourced `AbsorbingLayer` inputs. Each layer retains
its physical thickness, absorption coefficient in inverse metres, optical
depth, positive temperature support and original source. A layer has uniform
temperature and absorption; different layers may have different temperatures.
The caller chooses a spatial partition adequate for its material gradients.

Preparation retains one direction cosine through the entire stack before
angular integration. The exact ray factors use exponential transmission and
expm1 absorption. Positive adaptive GL8/GL16 integration prepares transmission,
reciprocal face escape/absorption and symmetric interlayer exchange once.
Splitting an isothermal layer preserves the full path instead of introducing
artificial angular redistribution. The homogeneous diffuse transmittance is
2 E_3(kappa d); thin diffuse absorptance tends to 2 kappa d.

Warm emission and layer-heat applications reuse immutable shared coefficients
without angular quadrature. Each layer receives its own external absorption
minus emission plus reciprocal heat exchange with all other layers.
Interlayer exchange is retained even when no light traverses the full stack.
The stable common-background form accepts signed emission deviations, provided
every input uses the same background and gray/band/per-wavelength scope.
Only complete finite results reach caller-owned heat output.

At most 64 layers and bounded angular work are admitted; cancellation and
exhausted accuracy/work controls refuse preparation. The retained embedded
coefficient discrepancy and incident-energy partition residual are numerical
diagnostics, not certified material or continuum error bounds. A finite
positive thickness/absorption product that overflows takes the opaque limit.

An undeclared layer uses refractive index one. `new_refractive` accepts a
sourced real index n >= 1, common to the entire stack, in an index-one exterior.
It preserves Snell refraction and the separate s/p Fresnel reflection series
before averaging unpolarized power. All incoherent round trips are summed
analytically, including boundary-return paths between nonisothermal layers.
Total-internal-reflection modes remain in the internal exchange, with the n²
radiance/etendue factor; adjacent opaque layers therefore exchange n² times
their vacuum blackbody-exitance difference. The common-index stack has equal
front/back reflectance, while individual layer escape fractions can differ.
Both reflection and transmission enter the incident partition diagnostic.

`AbsorbingLayer::with_refractive_index` supplies a distinct sourced index per
layer to the same `AbsorbingSlab::new` preparation and all existing heat consumers.
The common Snell coordinate p = n sin(theta) retains angle memory through every
interface. Integration splits at every critical p=n; the common measure d(p²)
includes local n² thermal radiance and channels trapped in high-index subregions.
Prefix/suffix reflection sweeps sum all repeated passes with positive independently
retained reflection deficits. The full layer exchange costs O(N²) per ordinate;
there is no growing transfer matrix or iterative bounce truncation. Lossless
trapped channels exchange no heat, while small positive absorption survives
reflection rounded to one. A real low-index gap blocks its nonpropagating channels.
`reflectances()` retains distinct front/back values; `reflectance()` denotes the
front. The common `refractive_index()` is now optional, with original individual
indices and source strings retained by `layers()` for heterogeneous materials.
Uniform indices reuse the original preparation path. G0/G1/G3 regressions compare
independent dense two-direction boundary equations, lossless stacks, low-index
gaps, weak trapped absorption, reversal, subdivision and equilibrium.

`SlabBeamResponse::new` prepares collimated polarized illumination of the same
ordered layer model without preparing diffuse angular coefficients.
`AbsorbingSlab::beam_response` shares the existing stack and its sources.
The admitted exterior cosine is in (0,1]; the supplied s fraction is incident
power, with the remainder in p. Conserved Snell direction, independent positive
Fresnel R/T and Beer factors feed the original prefix/suffix elimination. Beam
preparation and application cost O(N); no thermal pair matrix, angular
quadrature, iterative bounce truncation or additional optical solve is needed.
Matched-index grazing cosines and small positive absorption/escape are retained.

`exchange` accepts nonnegative front/back irradiation already projected onto
physical slab area, in W/m². It applies no second cosine and returns directional
reflected/transmitted power separately from layer heat and the signed arithmetic
energy residual. Every layer's absorption comes from positive propagation loss,
never subtraction from rounded unit reflection/transmission. Only finite complete
results reach the heat output. `layer_net_fluxes_with_beam` adds that heating to
the original diffuse/emitting heat balance after verifying the same ordered
source-bearing layers; outgoing beam power stays directional. Opposing beams
are mutually incoherent. The caller owns each source's gray/band scope and any
angular/spectral ensemble weights; no diffuse-port redistribution is inferred.
Focused regressions compare independent dense boundary equations, Beer and
Brewster limits, directional integration against the original diffuse response,
reversal/subdivision, thermal superposition, weak loss and transactional refusal.

`layer_fluxes_at_depth` resolves forward/backward beam flux inside an original
material layer; interface queries select that layer's side. `layer_absorbed_flux`
integrates positive heating over a caller-owned depth interval, giving W/m² of
slab area or cell-average W/m³ after division by interval thickness. Preparation
retains the Snell cosine and both endpoint fluxes for each incident side, with
separate s/p histories combined only by their supplied power weights. Each query
costs O(1), without rebuilding rays, integrating angles or subdividing material.
Interval heat uses positive Beer attenuation and expm1 loss, not subtraction of
nearly equal fluxes. Original finite absorption and queried distance retain shallow
penetration even when the whole layer's optical depth overflows. These beam-only
profiles add no inferred diffuse thermal depth distribution. Regressions compare
the original independent dense two-stream fields, whole-layer heat, additive
intervals, subdivision/reversal, weak narrow cells and shallow opaque profiles.

`AbsorbingSlab::diffuse_depth_response` prepares `SlabDiffuseDepthResponse` for
one caller-owned depth interval in an original physical layer. Uniform,
unpolarized external radiance is integrated with 2 mu dmu weighting of the
original beam fields, including all Snell/Fresnel return paths. Positive exact
Beer interval losses retain weak and narrow-cell heating; the physical optical
stack and source identities remain shared and unchanged. Preparation has bounded
adaptive GL8/GL16 work and retains accepted panels, actual angular evaluations
and front/back coefficient discrepancies. Subsequent `absorbed_flux` applications
cost O(1), without allocation, quadrature or additional optical solves. Inputs
are the hemispherically integrated front/back irradiances in W/m²; the result is
interval heating in W/m², or cell-average W/m³ after division by positive cell
thickness. Zero-width and nonabsorbing intervals need no angular evaluations.
This response assigns only external incident absorption, not a depth profile for
thermal emission or interlayer exchange. Regressions use independent E3 interval
oracles, original whole-layer heat, additive spatial cells, refractive reversal,
weak narrow and shallow opaque cells, cancellation and exhausted work controls.

`AbsorbingSlab::depth_response` prepares net heating of an original depth
interval from external diffuse irradiation and all original thermal layers.
`SlabDepthResponse` retains positive face and interlayer couplings, integrating
the existing conserved-Snell elimination across every relevant critical-index
band. The common d(p²) measure includes n² thermal radiance and trapped channels.
Original finite Beer interval factors preserve weak losses without inserting
physical layers or interfaces.

`net_flux` applies each coupling to the incident/source exitance minus target
exitance. Detailed balance cancels self-emission and reabsorption analytically;
uniform equilibrium is exactly zero, including large common backgrounds.
Finite signed deviations from one common background are supported. Each
original layer keeps its uniform temperature and material scope. Applications
cost O(N), without allocation, angular integration or repeated ray solves.
Bounded preparation retains work counts and embedded coefficient discrepancies,
not certified error bounds. Tests cover E3 thermal-cell integrals, full-layer
correspondence, interval sums, reciprocal cell couplings, port balance, opaque
n² exchange, weak trapped absorption and refusal paths.

G0/G1/G3 tests compare transparent two-interface reflection, opaque n² exchange,
independent polarized boundary-return ray sweeps, isothermal subdivision,
reversal, and equilibrium. The index source must support the layer temperature
intervals and declared gray/band scope. Real-index Fresnel boundaries assume
negligible interfacial extinction: this is not a conductor model. There is no
scattering, within-layer index gradient, rough interface, coherent interference or
evanescent tunnelling. Diffuse thermal incident radiance has emergent hemispherical
power is generally not Lambertian.
A scalar view-factor consumer explicitly owns its exterior diffuse-port
approximation. Gray inputs do not become spectral measurements. Consumers
admit actual layer temperatures against every source, supply their own Planck
scope, and separately assess angular, spectral and spatial model errors.

## Source contacts for original Cp, Cv and assembled storage

`transient::deposition::contact::implicit::PreparedCaloricContact` binds an
original ordered interface source to two `CaloricConductionLayer` objects.
It accepts source-exact Cp enthalpy, fixed-density Cv internal energy,
fixed-volume sensible storage and declared mass/volume assemblies. The actual
active link determines A/B order and area; the original caloric objects must
be registered on its cells and validate their mass/volume. Exposed material
cards must belong to those original bodies. For assemblies, the caller
explicitly chooses a constituent or effective-conductivity surface card;
membership is not an inferred mixture geometry.

The adapter reuses `PreparedMaterialContact`'s compiled scalar/curve/surface
engine. One physical query at the actual initial body temperatures selects
the original source and retains its `ResolvedInterfaceStatePoint` receipt.
Trial evaluations preserve every applicability coordinate, source knot,
triangulated cell and support hole without source queries or reselection.
Temperature and optional declared mean-temperature response use the original
inverse-storage tangent: `1/Cp`, `1/Cv`, `rho_ref/C_vol`, or
`m_ref/C_total(T)` as appropriate. Declared phase-fraction axes refuse;
the existing phase contact retains its original phase/fraction semantics.

`energy_bounds_j_kg` and `state_at_specific_energies` retain each body's native
energy coordinate. The historical h-labelled native trait forwards those
same values and derivatives without promoting internal energy or assembled
storage to enthalpy. `IntegratedCaloricLink::new_with_contact` inserts this
resistance exactly once and applies its two energy derivatives in the existing
integrated-conductivity residual/Jacobian. The original thermal owner remains
responsible for conservation, entropy and accepted time evolution. Contact
temperature coordinates refer to the declared bulk endpoint states; fixed
pressure, gap, finish and history are not mechanically evolved. No contact
heat inventory or additional thermal solver is introduced.

The public `deposition_contact` regressions add a mixed Cv/Cp native step with
an independent quadratic-energy, one-scalar transferred-heat oracle and exact
flux tangents, plus mass/volume assembly scaling, original identity, support
hole, axis and cancellation refusals. These native tests were authored but
not executed in the implementation environment. Syntax parsing and an
independent scalar arithmetic audit are not compilation or native execution.

## Exact fixed-reservoir convection through sensible and latent enthalpy

`lumped::exchange::solve_lumped_linear_exchange` evolves an original
`LumpedEnthalpyBody` against a fixed-temperature bath with its held convection
coefficient, area and mass. It delegates to the original phase curve's
`linear_heat_exchange_with_checkpoint`, using exposure `G * duration / mass`.
The material operation integrates `dh/dX = T_bath - T(h)` through retained
source-affine Cp branches, authored affine T(h) intervals (including mushy
spans), and sharp fusion plateaus. A plateau consumes its latent inventory
at the actual constant driving temperature. Heating and cooling use the
same supplied equilibrium chart; latent heat enters the storage account once.

The phase operation visits each crossed interval once, allocates no trajectory
or temporary chart, and reuses the original bounded sensible-exchange kernel.
No temporal subdivision or apparent heat capacity is introduced. It reports
the integrated relative heat, exposure-weighted temperature, represented
endpoint enthalpy change and their actual residual separately. Original
material/curve identities, state density and source support remain intact.
The bath may lie outside source support; the material trajectory may not.
Completely lost absolute-enthalpy increments refuse. Partial quantization
remains visible rather than being hidden by the arbitrary enthalpy datum.

The native body adapter retains the conservative full-span conductivity/Biot
gate and requires radiation to be explicitly disabled. It admits the actual
stored-energy residual against the caller's absolute joule tolerance before
returning a complete immutable proposal. Zero conductance or duration retains
the initial state. Cancellation before work, during chart traversal or before
publication returns no accepted proposal. This is a fixed-geometry, fixed-mass,
fixed-conductance, infinite-bath reduction. The caller retains the original
chart's pressure/work convention; internal power, finite baths, radiation,
geometry evolution and nonequilibrium phase kinetics require their existing
physical owners rather than reinterpretation of this exchange.

Six material and four native lumped regression definitions exercise independent
analytical sensible/latent heating and cooling, latent residence, equilibrium,
source subdivision, exchange composition, datum shifts, authored mushy laws,
weak heat, original identity/support and cancellation, integrated convection
means, Biot/radiation admission and the absolute state-quantization budget.
These definitions were syntax-checked but not natively executed in the
resource-exhausted implementation environment.

## Recoverable nonlinear trials in scheduled enthalpy adaptation

The native scheduled `EnthalpyHistoryRun` and `Sdirk2EnthalpyHistoryRun`
adaptive methods, including their entropy variants, opt into the same physical
recovery classes as the standalone adaptive enthalpy owner:
`NotConverged`, `EndpointResidual` and `EnergyBalance`. The existing selector
halves a failed interval within its original limits. It never changes a
physical tolerance, accepts a failed trial or treats a refusal as a time-error
estimate. A floor, attempt cap or nonrepresentable smaller clock returns the
actual last physical error.

Every retry resamples the original scheduled endpoint/stage laws and starts
from the unchanged accepted enthalpy. Only a complete admitted fine pair
publishes fields, clock and heat/entropy. Attempt counts include exactly the
interval solves entered before a failure: one failed coarse interval followed
by a successful comparison reports four solves. Existing selector callers
that have not opted into recovery retain immediate refusal. Scheduled SDIRK2
contact/reflection retain their native owners; separate legacy BE wrappers
keep their existing policies.

Four focused `thermal_phase_adaptive/recovery.rs` definitions use an independent
one-free-node sharp-phase spatial equation to exercise melting/freezing,
authentic SDIRK2 stages, accepted-only heat/entropy, continued history,
typed floor/cap refusal, cancellation and a scheduled bath outside source
support that must not be hidden by a smaller interval. These new sources
passed syntax/API review but native compilation and execution remained blocked
by the shared implementation environment's exhausted resources.

## Analytic enthalpy response for exact phase exchange

The original phase curve's `linear_heat_exchange_response` and checkpointed
variant return the same fixed-bath exchange with its analytic enthalpy/heat
Jacobian. One traversal composes the original sensible-kernel inlet gain,
stable transferred-heat loss, bath response and positive mean-temperature
moments. A latent segment contributes unit inlet gain and bath response equal
to its exposure. Continuous T(h) preserves the variational flow through
interval transitions without perturbed marches or finite-difference timing.

The exposure derivative retains the attenuated driving temperature even when
the rounded endpoint equals the bath. Normalized positive moments preserve
weak mean-temperature response without subtracting nearly equal quantities.
At positive-time equilibrium, an interior chart knot with unequal one-sided
caloric slopes has no unique two-sided Jacobian and refuses this response
entry. The value-only entry still accepts that equilibrium. A latent interior
retains its actual phase-inventory response; support endpoints have admitted-side
derivatives.

`lumped::exchange::solve_lumped_linear_exchange_response` applies the original
body's mass and conductance to return integrated-heat derivatives with respect
to initial specific enthalpy, bath temperature and duration, in kg, J/K and W.
It shares the same Biot, radiation-disabled, source-domain, cancellation and
actual stored-energy residual gate as the value-only operator. These derivatives
describe the continuous original caloric map, not floating-point quantization
of absolute h. They do not change geometry, material context or reservoir
ownership and do not themselves accept a coupled finite-reservoir history.

Four added material and two added native regressions, with existing mushy and
energy-budget cases extended, check analytical cross-fusion heating/cooling,
source subdivision/composition, weak and rounded-to-bath response, latent
equilibrium/corner admission, native chain units and independent three-phase
finite differences. They passed syntax and bounded arithmetic/API review;
native compilation and execution remained unverified in the resource-exhausted
implementation environment.

## Scheduled material ambient radiation in native SDIRK2 enthalpy

`MaterialRadiativeEnthalpyHistoryRun::into_sdirk2` retains the accepted h/T/f
fields, physical clock, step budget, original chart objects, pinned emissivity
claims and every independently scheduled black reservoir. An unadvanced native
or material-contact SDIRK history may attach those same bindings through
`with_material_ambient`. Ambient and reflected-enclosure radiation have
exclusive ownership; attaching a second radiation owner refuses.

Both stages sample each reservoir at its actual physical stage time. Original
pointwise material emissivity and its enthalpy/phase feedback are evaluated by
the existing radiation residual and tangent, including latent plateaus.
Source, boundary and all bath knots share the existing event clipping; the
event endpoint uses the original left-limit semantics. Only the stiffly accurate
endpoint, or the admitted adaptive fine pair, becomes accepted history.

Native contact, source and convection accounts remain distinct from radiation.
The existing stage energy, physical Newton, source-domain, cancellation and
adaptive h/T/phase gates are unchanged. Entropy uses the native stage-weighted
material ambient radiation account and original-to-final caloric storage; no
new photon, spectral or participating-medium model is inferred.

Five focused regressions cover retained advanced history and source claims,
independent bath event timing and entropy, ramp-driven adaptive phase rejection,
second-stage/cancellation/identity refusals, and spatial material contact with
independent radiating faces and baths. Their scalar and eight-coordinate
references derive stage heat and entropy independently. All six Rust files
passed the pinned parser; native compilation and test execution remained
unavailable in this resource-exhausted session.

## Direct finite-body phase exchange on the lumped contact owner

`lumped::exchange::pair::solve_lumped_pair_exchange` integrates one isolated
fixed-resistance contact between two admitted `LumpedEnthalpyBody` values.
The existing `SeriesThermalResistance` supplies G=1/R, and each body retains
its original fixed mass, identity and equilibrium enthalpy/phase chart. The
original contact Biot check Lc/(R A k_min) uses the body's effective exchange
area and conservative full-chart conductivity minimum.

The shared material finite-pair kernel follows the original two-body
trajectory through sensible branches and latent plateaus. Affine T(h)
intervals advance analytically; varying original Cp uses bounded source-split
integration and inversion with separate heat/exposure error gates. This
avoids an endpoint-power temporal approximation for held isolated contact.
The existing backward-Euler and entropy-conjugate contact APIs continue to
own their externally supplied-heat models.

`LumpedPairExchangeConfig::energy_tolerance_j` independently gates EACH
represented stored-energy residual and their sum. A small combined residual
cannot hide offsetting local errors. The result retains both identities,
duration, resistance, maximum Biot number, original-chart entropy and all
material trajectory/representation diagnostics. Cancellation, source
exhaustion, equilibrium overshoot, unrepresentable heat and either physical
energy-budget violation refuse the complete proposal without mutating a body
or accepted history.

The contact has no external bath, radiation law, applied heat or mass flow.
Callers retain composition with other physical owners and the charts'
pressure/work convention. Density and volume changes are not added again as
latent heat or pressure work; spatial spreading is outside the lumped model.
Focused material/native regressions cover melting and freezing, reciprocal
transfer, energy and entropy, source limits, original Cp variation, contact
Biot refusal, representation errors and cancellation. Pinned Rust parsing
was available; native Cargo compilation and test execution were unavailable
in this resource-exhausted session.

## Immersed dielectric slabs

`AbsorbingSlab::new_in_medium` prepares diffuse incidence and thermal emission
with one explicit, sourced, lossless exterior dielectric shared by both faces.
`new_between_media` instead retains independent [front, back] exterior indices
and sources. Diffuse `transmittances()` stores power fractions FROM each face,
with n_front² T_front = n_back² T_back; singular index/source/transmittance
accessors designate the front face. `outgoing_fluxes` applies both directional
transmissions together with reflection and the original layer emission.

`SlabBeamResponse::new_between_media` and `beam_response_from_side` author
an angle in either exterior medium. The returned exterior cosines identify the
conjugate Snell channel, where directional power transmission is reciprocal.
A nonpropagating exterior channel has zero cosine and refuses nonzero incident
power; it never acquires invented evanescent heat. The original `new_in_medium`
and `beam_response` entries retain the same common-medium/front-angle behavior.
All layer indices remain absolute. Conserved p=n*sin(theta), split at every
critical index, carries polarized Fresnel reflection and all propagating
internal returns without angular mixing or exponentially growing matrices.

Incident irradiation is physical W/m² through the face. `face_coupling` stores
fractions of that incident power. Layer blackbody inputs remain vacuum
Planck/Stefan--Boltzmann exitances; escaping thermal emission uses exterior
n² times those face coefficients. Interlayer coefficients retain the absolute
d(p²) mode measure, including trapped rays. Thus an equal-temperature bath
uses irradiance n²*B and gives exact zero layer heat. Scaling every index by
s preserves R/T/incident absorption and scales thermal emission/exchange by s².

Depth-resolved diffuse and beam heating share the same original endpoint
fields. Critical-angle quadrature uses positive Snell discriminants and a
sine-squared band map; it does not reconstruct nearly tangent internal angles
from rounded exterior cosines. A blocked layer has zero propagating field and
zero beam heat. A low-index gap can separate independently illuminated regions.
Repeated prepared heat applications require no angular solve.

This is smooth, incoherent geometric optics with bulk Beer absorption. It does
not infer evanescent absorption, frustrated reflection, coherent films, lateral
scattering or interfacial complex extinction. Sources must cover the original
temperature and gray/band scope. Beam/thermal superposition requires the same
pair of exterior indices; the air-enclosure consumer explicitly requires
index one on BOTH faces. Native gray thermal cell sources and timestep bounds
use each face's own n² emission factor. Local diffuse and net-depth queries
retain the same exterior pair and critical-angle bands.

Six additional G0/G1/G3 groups cover independently evaluated asymmetric
Fresnel return series and diffuse transmission, each face's thermal equilibrium,
whole-layer/local heat balance, reversal, subdivision, index scaling and closed
channel admission. Native gray-cell integration has a separate authored
regression; its full dependency lane is not implied by the slab source probe.

Regressions cover independent Fresnel escape-cone integration, matched-medium
Beer exchange, equilibrium, energy balance, index rescaling, reversal, local
heating, beam superposition, cancellation and bounded work. The optical slab,
beam, local-depth and original air-enclosure admission checks passed in a
55-test native source probe with the actual fs-math implementation. Native
gray/spectral thermal-cell integration remains unexecuted in this environment;
this probe is not a full workspace or material-admission dependency test.

The Fresnel, radiance and incoherent-layer assumptions follow
[PBRT 4e §9.3](https://pbr-book.org/4ed/Reflection_Models/Specular_Reflection_and_Transmission)
and [Byrnes §6](https://arxiv.org/abs/1603.02720).

### Radiation admission for distinct conducting surfaces

`LumpedThermalTransport::radiation_is_disabled` reports whether the retained
maximum emissivity is exactly zero over the admitted transport span. This
read-only query lets a composing heat owner require zero body-node radiation
before eliminating a distinct conducting surface. It does not change the
original emissivity model, source support, Biot admission or radiation
temperature; relocating radiation requires its own physical surface balance.

## Typed scalar thermal scenarios on the original boundary and history owners

`ThermalBc::from_scenario_row_at_in_frames` consumes canonical typed scalar
Dirichlet temperature, outward Neumann heat-flux density and Robin heat-transfer
coefficient through `BoundaryCondition::thermal_at`. Exact dimensions and
physical semantics are checked before lowering: semantic temperatures must be
absolute, even-parity quantities, never temperature differences. The original
native boundary owner retains flux signs, coefficient positivity and the
separately supplied Robin reference temperature. Legacy scalar/time-signal
rows retain their existing behavior. The steady entry still refuses to invent
a time for a typed history, including a fixed scalar with a time-dependent frame.

`ThermalBoundaryHistory::new_in_frames` and `from_boundary_in_frames` borrow
the original frame tree for the trajectory. The existing constructors retain
the implicit world tree. At most 256 frame rows are admitted. Source and frame
tables are not cloned at endpoints; the same immutable exterior-face partition,
shared Dirichlet precedence, surface conduction and original material/thermal
steppers remain authoritative. Coordinate frames validate at actual physical
time even when source values clamp or wrap, and do not move the bound mesh.

The existing event clock now visits held/linear typed-table knots and used
ancestor translation/tilt knots; fixed scalars add no source knots. Refuse and Clamp tails retain their
policies; exact held left values are taken only after quantity, domain and
physical-time frame checks. Periodic typed tables remain usable for explicit
endpoint sampling, but event stepping requires a finite unrolled schedule and
otherwise refuses without advancing state. Existing accepted-step budgets,
cancellation, energy closure and atomic publication are unchanged. The shared
history boundary feeds the original linear/nonlinear, radiation and enthalpy
paths; no new time integrator or constitutive law is introduced.

Focused regressions extend the original boundary-history and thermal-event
targets. They specify frame-bound partition preservation, native-operator
equivalence for all three typed thermal boundary families, a 120 W/m² inward
pulse on a 1 m² exterior face delivering 15 J, ancestor-frame event clipping,
and cancellation/tail/budget refusal with accepted-state preservation. Scenario
tests specify canonical replay and semantic admission. These new thermal tests
have not executed, and the new source has not been typechecked or formatted by
Rust tools: the execution environment disconnected during dependency
compilation before edits. Source review is the available evidence; no passing
runtime or physical-validation claim is made.


## Physical vector flux on the original thermal faces

`ThermalBc::NeumannVector` prescribes a spatially uniform WORLD heat-flux
vector in W/m². Its checked constructor retains all three components. Each
original exterior face receives `q · n_out`; shared vertices never average
normals, and one named region may cover differently oriented faces. The same
projection feeds original P1 surface assembly, total prescribed heat, entropy
exchange, isolated-equilibrium admission and optional flux verification.
Scalar nodal flux retains its previous face integration. The native boundary
builder refuses a nonfinite projected normal flux before publishing a boundary.

Typed scenario vectors lower through `heat_flux_world_at` into this native
row. The original history owner preserves its exact partition, Dirichlet
precedence and static surface conduction. The shared event clock visits vector
source knots and used ancestor-motion knots. At held jumps, the previous three
components are rotated at the actual endpoint; neither physical time nor the
frame is shifted by an epsilon. Periodic event histories still require finite
unrolling. Native heat/cancellation/budget checks precede accepted-state publication.

The existing radiation helpers continue requiring their declared uniform
normal-flux or Robin traces on radiating regions; they do not silently replace
a vector row. Unrelated vector boundaries remain in total enclosure heat.
This adds no spatially varying vector field, moving mesh, continuous-time work
integral for rotating forcing, or vector-control adjoint.

Regressions compare a multi-normal native vector solve with independently
assigned scalar face rows and check analytic heat balance. History regressions
retain the partition under moving axes and compare the original backward-Euler
pulse account with scalar forcing, including a rotated held left value, source
and frame knots, tails, cancellation and exhausted budgets. Pinned Rust syntax
parsing and source review completed; compilation and runtime tests remain unrun.

## Original scalar and vector flux controls in complete thermal adjoints

`PrescribedFluxControls` binds the original P1 Neumann faces in explicit control
order. Scalar controls are uniform normal-flux offsets; vector controls retain
all three WORLD components and each face's own outward normal. Their tangent
loads enter the existing `RobinResponse::apply`; their pullback consumes the
complete masked nodal-load multiplier from steady or transient owners. The
outward-load minus is applied once. No additional solve, state feedback,
Dirichlet lift or timestep factor is introduced by the contraction.

The retained basis is sparse in selected faces. It replaces the dense
zero-feedback trace workaround in native CLI capacity/enthalpy adjoints and
also supplies the CLI steady scalar flux derivative. Source amplitudes, history,
physical face partition and accepted fields remain independently owned.
Original scalar/vector coordinate, fixed-node, finite-arithmetic, face-budget
and cancellation checks apply. Geometry, frame/schedule chain rules and direct
objective dependence remain outside this fixed native-flux derivative.

Focused regressions define native-load differences, scalar reduction, skew-face
vector chain rules, and all-component physical finite differences plus transpose
identities through nonlinear contact steady and backward-Euler solves. Pinned
Rust syntax parsing and scoped diff checks passed; Cargo typechecking and these
runtime regressions were not executed in this session.
