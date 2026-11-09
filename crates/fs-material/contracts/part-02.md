<!-- fs-material contract, part 2 of 4; a higher number is newer. Index: ../CONTRACT.md. Keep this file under 900 KB: the GitHub contents API returns empty content above 1 MiB, and a writer that sees an empty file overwrites it. -->
## Implicit spatial freezing with original source tangents

`phase::water::freezing::transport::FreezingConductivityState` now retains
analytic local dT/dh and dk/dh with its original material/transport response.
Sensible branches use the reciprocal original Ice Ih or IF97 Cp and the
original conductivity temperature derivative. Fusion retains exactly zero
dT/dh. Its conductivity derivative follows the physical liquid VOLUME
fraction, dphi/dh = vi*vl/(v_mix^2*L), through the caller's series or parallel
closure. Both exact fusion endpoints select the plateau-side generalized
derivative. No apparent heat capacity, smoothed transition, finite-difference
runtime tangent, extra latent source or shifted energy reference is used.

`FreezingWaterImplicitLaw` borrows the exact original caloric and transport
objects, verifies chart equality, and provides one fused source query to the
shared thermal chain adapter. Bounds and all source-domain refusals remain
those of the original immutable material chart. The extra scalar tangents are
computed during the existing query and add no allocation or inverse query to
its getters.

`freezing::spatial::FreezingSlab::advance_implicit_to` uses the original
backward-Euler residual and the reusable `ImplicitConductionWorkspace`.
Its endpoint conductances include the actual phase closure and face half-cell
resistance. The same slab can continue with explicit or implicit transport.
Both paths share independent original-source Delta H=Q and Delta U+p*Delta V=Q
checks and retain original thermal/reservoir entropy. The nested implicit
report repeats the physical energy account for comparison; callers must not
add it a second time. Cancellation or failed nonlinear/work/energy admission
leaves the original slab intact, and numerical scratch remains reusable.

The existing `water_freezing_front` example accepts `--implicit`,
`--cells N` and `--dt SECONDS`. For example,
`cargo run -p fs-material --example water_freezing_front -- --implicit --cells 96 --dt 0.5`
computes the 120-second 10 mm slab with 244 accepted intervals and a maximum
independent energy residual about 1.1e-10 J in the native run. This is algebraic
closure, not physical or time-discretization accuracy. At 120 seconds the
frozen-reference equivalent length is 5.456292577 mm; halving the maximum
step gives 5.449258361 mm. These are phase-inventory lengths on the fixed
reference mesh, not mechanically resolved ice thicknesses.

All 24 focused freezing-namespace tests pass against the actual updated
thermochemistry and thermal libraries. New tests compare analytic source and
phase-volume tangents with independent differences, verify exact borrowed
bindings, resolve an actual source-backed front beyond the initial explicit
CFL limit, compare three implicit timesteps with a fine explicit reference,
reverse the bath at an exact authored time, reproduce chunked continuation
bit-for-bit, and reuse scratch after cancellation/work refusal. Refinement
errors remain separate from the much smaller joule acceptance tolerance.

All previous physical limits remain: atmospheric source/transport overlap,
declared phase morphology, fixed reference geometry and external pressure.
This addition supplies no convection, stress, nucleation, supercooling,
rigid-vessel coupling or moving surface. Implicit interval retries do not
provide automatic temporal error control or unconditional nonlinear
convergence. Full workspace/DSR qualification remains unavailable in the
current checkout; the native validation compiles the original module sources
against actual production dependency artifacts, without numerical stand-ins.

## Source-defined expanding thermal surfaces

`thermal::transport::expansion::surface::PreparedExpansionSurface` binds explicit
reference areas grouped by principal normal to one original
`PreparedFreeExpansion`. At each admitted temperature it derives
`A_i = A_i,ref * s_j * s_k` and the exact temperature tangent. Another consumer
can reuse an already evaluated original tangent only when the full prepared
path identity matches. Area evaluation performs no runtime allocation, hashing
or material query. Original card/claim references, fixed query conditions,
reference normalization and the declared principal/isotropic basis remain
authoritative. Zero selected faces remain zero; electrical L/A does not infer
a thermal surface or boundary coverage.

The `expansion_surface` regressions cover normalized anisotropy and one-sided
source knots, exact tangents and selected/zero patches, legitimate derivative
cancellation, foreign paths and domain/area refusal. `with_oriented_patches`
also admits bounded planar reference-area vectors in the same principal frame.
The exact Nanson map `a = cof(F) a_ref` gives each current area and its temperature
tangent; oblique faces are not counted as independent projected principal faces.
Totals include both portions and retain allocation-free runtime evaluation.
The triangle-cross-product, normal-sign/axis-permutation, isotropic-limit and
extreme-area tests exercise this geometry; the circuit expansion boundary test
passes its total area into the original implicit heat/entropy solve.
Curved surfaces require an explicit facet/quadrature approximation. Shear,
constrained stress and moving interbody contact coverage require their own
geometry owners. Thermal boundary laws and disjoint coverage admission remain
with the consumer; surface geometry alone adds no heat or storage.

## Sealed dense-water heat exchange on the shared thermal timeline

`phase::water::dense::thermal::DenseIsochoricThermalLaw` retains the original
prepared IF97 Region 3 isochore. Its declared storage basis is internal energy,
with original T(u), whole-interval positive-Cv slope bound, local analytic
dT/du = 1/Cv, and original-source entropy increments. Birth admission requires
m = rho V within a 64-epsilon arithmetic allowance and refuses unsupported u.
Mass advection and the existing sublimation owners refuse this storage basis.

`DenseWaterThermalNetwork` couples sealed rigid volumes through the original
thermal contacts, reservoir boundaries, events, timestep owner and immutable
candidate publication. All volumes are born at time zero. Stable source
differences and represented endpoints independently check Delta U = Q and
Delta H - V Delta p = Q. Moving-boundary work is exactly zero; V Delta p is
an enthalpy contribution, not wall work. Thermal-coordinate entropy and original
source entropy are reported separately, exposing finite inverse representation
effects without discarding tiny representable heat.

Four new integration tests, together with four existing dense-source tests, cover
source Cv/tangents and tiny-heat entropy, mass/volume admission and actual
advection refusal, finite-body equilibration, exact bath events, continuation,
cancellation, budgets and source-domain failure. The bath transient converges
against independent quadrature of t = (m/G) integral Cv(T)/(T_bath - T) dT:
at 20.695034947 s, temperature errors for maximum steps 0.4/0.2/0.1 s are
0.241386644/0.120461412/0.060267601 K.

`water_dense_vessel --bath` computes finite-conductance heating and cooling,
pressure and reservoir entropy; its original prescribed-heat mode remains
available. The example's 0.5 kg, 1-litre vessel with authored 25 W/K conductance
reaches 732.426872 K and 68.727035 MPa at 100 s, with Q = 120810.286679 J,
Delta H = 163953.619467 J and V Delta p = 43143.332788 J. Across 800 steps its
maximum history energy residual is 2.6663e-7 J.

This is a spatially uniform network on one fixed-density chart. Conductances are
authored; dense-water conductivity, container strength, wall storage/compliance,
internal convection, multiphase behavior and regional stitching are not supplied.
Native validation compiles the actual dense namespace and complete thermochemistry
and thermal owners; it does not constitute a full workspace gate.

## Adaptive spatial freezing with actual phase-inventory accuracy

`FreezingWaterImplicitLaw::time_observables` obtains temperature and liquid
mass fraction from one original caloric inverse. `FreezingSlab::advance_adaptive_to`
uses the shared implicit accuracy controller and retains the existing independent
source Delta H, Delta U + p Delta V and entropy accounts. No apparent heat
capacity, duplicate fusion source, extrapolated accepted state or synthetic
temperature-to-phase relation is added.

`FreezingSlabAdaptiveAudit` separates the physical audit from the controller's
local h/T/phase estimates and all attempted nonlinear work. Source and controller
failures remain transactional, including cancellation after the final original-
source audit. A returned `suggested_dt_s` explicitly seeds continuation.
Scheduled face-temperature changes retain their exact event times and actual
reservoir entropy temperatures.

Two additional source tests bring the actual freezing namespace to 26 passing
tests. A 24-cell, 12-second comparison against independent explicit transport at
a 0.002-second cap reduces mass-weighted absolute enthalpy-field error from
1.637753637 to 0.226923967 J when local h/T/liquid-fraction allowances tighten
by 100. Maximum temperature discrepancy falls from 0.046597616 to 0.006387667 K,
maximum liquid-fraction discrepancy from 0.011288234 to 0.001564972, and bath
heat discrepancy from 1.498449183 to 0.207880577 J. Accepted comparison pairs
increase from 141 to 1382: tighter accuracy costs additional original solves.
The second case tests an exact cold-to-warm bath event, continuation at that
event with the returned proposal, two-bath entropy and final-audit cancellation.

`water_freezing_front --adaptive` exposes dimensioned enthalpy and temperature
allowances, a liquid-mass-fraction allowance, a relative change allowance, and
minimum/maximum comparison spans. The 24-cell, 120-second default adaptive run
accepts 762 pairs (1524 physical half steps) from 817 comparisons. Its maximum
accepted error ratio is 0.925807011 and maximum independent calorimetric residual
is 1.1102e-9 J. Computed frozen reference-equivalent depth is 5.225752655 mm.

These are local temporal estimates and computed reference-geometry results,
not global or spatial accuracy certificates. Spatial refinement, source ranges,
the explicit mixed-cell conductivity closure, and the fixed-pressure/no-flow
limits of the original slab remain independent requirements.

## Complete concentric rod thermoelastic potential

`ConcentricRodThermoelasticLaw` extends the original axial thermoelastic material law with the native annular bending, shear and torsional potential. It exposes the same potential, entropy and mechanical heat-capacity derivatives for all six small-strain coordinates. The derivatives come from the pinned material sources, with integrated thermal expansion and the original section moments; there is no independent coefficient table or material-name dispatch.

The admitted chart has fixed reference geometry, scalar reference density and expansion, and scalar or globally affine Young and shear moduli. Source domains, finite results and the original strain envelope are checked. More general temperature laws and phase-dependent geometry are not inferred. The `fs-couple` rod enthalpy regressions compare all components against independent formulas and temperature differentiation, including parity with the native rod energy. These new native tests await execution after the unchanged runtime dependency exceeded the authoring environment's resources.

## Shared tensor dynamics material binding

`visco_tensor::ResolvedTensorMaterial` resolves density, equilibrium elasticity,
every selected Maxwell tensor/time and strain limit from one material card/query.
It retains the source and intersects density and constitutive frequency domains;
dispersive density refuses. `fs-solid::TetViscoMaterial::from_source`, sourced
`fs-couple::TensorMembrane` facets and the production 3D valve parser share this
binding. Membranes rotate the complete hereditary law from declared material axes
into each original facet before the existing plane-stress elimination. Consistent
and lumped inertia, constitutive memory and physical heat retain their native
owners. Three `fs-couple/tests/tensor_material_source.rs` tests cover analytic
inertia/frequency, source-only heat changes, rotation/reindexing invariance and
source/domain/cancellation refusal; native execution remains unverified in the
memory-limited authoring worker. This is frozen-state, fixed-geometry, small-strain
dynamics; no thermal evolution, prestress, bending or acoustic radiation is inferred.

## Shared fluid acoustic response

`fluid_acoustics::thermoviscous::absorption::ClassicalBulkAbsorption` prepares
classical longitudinal absorption from one original `ThermoviscousFluid`
evaluation and an explicitly supplied nonnegative bulk viscosity. It uses
`alpha = omega^2 [4 mu/3 + zeta + kappa(1/Cv - 1/Cp)]/(2 rho c^3)` for pressure
amplitude, with the same density, adiabatic sound speed, transport and calorics
as the original liquid or real gas. No ideal-gas state, atmospheric default,
bulk viscosity from shear viscosity, or molecular-relaxation data is inferred.
The frozen owner performs no per-frequency material queries and exposes
`omega*tau` for weak-loss admission. DC is exactly zero; invalid frequencies
and unrepresentable positive losses refuse. G1/G3 tests cover the independent
specific-heat formula, separate bulk/shear/thermal effects, frequency scaling
and the existing gas classical limit. This is a uniform Newtonian weak-loss
model with constant bulk viscosity over the caller's supported band, not a
dispersion, phase-change, wall-loss or experimental-validation claim. Adding
resolved relaxation already included in the supplied bulk viscosity is invalid.

The same prepared owner exposes `thermoviscous_wavenumber(omega)` for the
acoustic branch of the full homogeneous linear Navier--Stokes--Fourier spatial
dispersion relation. It retains the original real-fluid density, adiabatic sound
speed, Cp/Cv, viscosity and conductivity. For `exp(i k r - i omega t)`, `Im(k)`
is amplitude attenuation and `Re(k)` supplies the associated phase; the separate
thermal diffusion mode is not returned. A rationalized quadratic root avoids
cancellation at weak loss and Cp/Cv close to one, with allocation-free constant
work per frequency. Zero frequency is exact; nonpassive, unrepresentable or
unresolved acoustic/thermal branches refuse, including roundoff near the
discriminant cut. The leading attenuation recovers Stokes--Kirchhoff, while
finite-frequency corrections satisfy the original continuity, momentum and
energy equations. This does not supply molecular spectra, boundary layers,
nonlinear response or experimental validity outside the original coefficients'
supported band. The parabolic transport equations do not establish a finite
wavefront; an observer's `r/c` reference delay is a numerical convention.

`ClassicalBulkAbsorption::thermoviscous_wave` returns the typed
`BulkThermoviscousWave` from that same acoustic root, with distinct total
compressive traction/velocity and thermodynamic pressure/velocity impedances
[Pa s/m]. With `eta = 4 mu/3 + zeta`, `t = p - eta du/dx`, hence
`Z_t = rho omega/k` and `Z_p = Z_t + i eta k` in `exp(i k x - i omega t)`.
The root owner retains the caloric contribution before subtracting viscous
stress; pressure remains resolved even when those latter terms would cancel,
including coincident viscous/thermal diffusion scales. DC returns the continuous
limit `k=0`, `Z_t=Z_p=rho c`. Nonfinite/unrepresentable impedances refuse.
Mechanical power uses `Re(t conj(u))/2`; complete thermal acoustic-energy flux
also includes `Re(T' conj(q'))/(2 T0)`. This acoustic eigenmode does not supply
the independent thermal mode or thermal end conditions. Native regressions
close the original momentum/heat equations and full energy flux, the exact
Cp=Cv viscous solution, separate shear/bulk heating, and weak/strong-loss limits.
The existing wavenumber evaluation order and branch/refusal policy are retained.

`BulkThermoviscousEquilibrium` binds the original `ThermoviscousFluid` and
`ThermalFluidTransport` views for consumers that retain both acoustic and
thermal diffusion. Both views must supply the identical transport tuple and
known signed isobaric expansion. The snapshot retains the actual absolute
temperature, signed expansion and explicit independent bulk viscosity, and
checks `beta_S = beta_T - alpha_p² T0/(rho Cp)` with
`beta_T = (Cp/Cv)/(rho c²)`. Unknown expansion, mismatched source views and
inconsistent caloric/acoustic derivatives refuse; Cp/Cv never supplies an
invented expansion sign. Positive, negative and zero expansion are admitted
when the original thermodynamic identity closes. Propagation and thermal
terminal conditions remain the responsibility of the consuming operator.

`fluid_acoustics::FluidAcoustics` admits one immutable acoustic snapshot with
positive finite density, adiabatic sound speed, bulk modulus `K_s = rho c²`,
compressibility `1/K_s` and plane-wave impedance `rho c`. Gas adapters reuse the
original complete gas coherence check. Original IF97 Gibbs, dense and unified
regional states supply their own density and sound speed; no ideal-gas cp/cv
substitution is made for a liquid or dense fluid.

The checked Peng–Robinson thermochemical owner supplies its actual stable-root
density and adiabatic sound speed through the same view. Its retained ideal
reference is provenance only and never supplies these coefficients. The
existing dense-gas caloric regression checks acoustic/transport parity and
preserves the actual departure from the ideal reference.

`ResolvedFluidAcousticState` retains the complete material query, original source
claims and every usage receipt. The caller explicitly selects a sound-speed or
adiabatic-modulus role and exact scalar property keys. Pressure semantics,
isothermal-modulus inference, tensor components and missing source properties
are not silently converted. Refresh pins the same original claims and refuses
unsupported state changes atomically. A numerical snapshot alone confers no
material-source authority.

`ResolvedAcousticLiquidState` resolves density, heat capacity, conductivity,
viscosity and acoustics at the same query. `LiquidPropertyRequirements` assigns
exact original source names to transport roles without aliases or claim mutation.
The dynamic/kinematic basis is explicit; kinematic viscosity converts with the
current density. Transport and acoustics must agree on the exact density bits.
Existing optional liquid expansion remains an independent source derivative.

The combined acoustic liquid retains its original `ResolvedLiquidState`, rather
than discarding that owner after extracting a transport tuple. `try_from_liquid`
adopts an already complete bundle with no query or claim reselection, including
an admitted density-expansion path. `ThermalFluidTransport` delegates to that
same owner: typed temperature and known, unknown, zero or signed expansion retain
their original meaning. Refresh evaluates the complete bundle once, retains the
original density path and nonthermal context, and rebuilds every numerical view
atomically. A transport-only bundle cannot silently acquire an acoustic claim.

The material regressions exercise both viscosity bases, exact source keys,
shared density bits, source/tangent refresh and refusal without modifying the
original state. Downstream regressions use this same object for natural
convection, accepted finite-cell enthalpy and cavity modes. This supplies no
compressible-flow, phase-change, attenuation or nonlinear-acoustic model.

Seven focused material regressions cover gas pressure scaling and stale-state
refusal, original IF97 derivative identities, source refresh/domain boundaries,
adiabatic-modulus semantics, both liquid viscosity bases, exact source vocabulary
and unrepresentable coefficients. This is local small-signal nondispersive bulk
acoustics, without attenuation, shear waves, phase selection or cavitation claims.


## Source-bound hysteretic magnetic cores

`field_storage::hysteresis::source::ResolvedScalarHysteresis` assigns explicit
original property keys to background permeability, fitted field limit and each
play branch's reversible slope, pinning field and saturation polarization.
Branch count/order are data; shared coefficients retain shared keys. The native
1024-branch limit and shared distinct-property budget apply. The complete card,
selected claims, quantities, evidence and original query survive lowering into
`UniformHystereticCore::from_source`. The existing play update, monotone inverse,
remanence, terminal work and physical pinning-heat owner remain unchanged.

This rung admits constant fitted parameters inside their declared domains.
Temperature-dependent fits and unowned dynamic field/frequency/stress/rate
applicability refuse. Heating additionally requires explicit temperature support
for every selected coefficient; a caller interval cannot widen it. Refresh
requeries the same claims, preserves geometry and constitutive memory, and
returns a candidate. It neither changes coefficients nor resets stored energy.
This supplies no thermal magnetic entropy law, vector/core spatial field model,
eddy/excess loss, Curie extrapolation or measured authority for synthetic fits.

Three focused G0/G1/G3 regressions cover independent two-branch work/remanence,
geometry reduction, data-defined branch sharing, source-domain/quantity refusal
and memory-preserving temperature continuation.


### Source tensor thermal potentials

`visco_tensor::source::thermal::TensorThermalLaw` admits source-derived
stored tensors over every original selected temperature knot.
The original affine isotropic E(T) path retains fixed Poisson ratio.
Independent affine bulk and shear sources use exact Mandel stiffness
derivatives: normal diagonal K'+4G'/3,
normal cross terms K'-2G'/3, and shear diagonal 2G'. Their derived Young modulus
and Poisson ratio need not be affine. Independent affine Young/shear sources
retain K(T)=E(T)G(T)/(3(3G(T)-E(T))) and its exact first and second derivatives.
Positive E, G and 3G-E across the complete interval remain required.
Independent affine isotropic E/nu and orthotropic E_i/G_ij/major nu_ij sources
use exact compliance derivatives and inverse derivatives:
C'=-C S' C and C''=2 C S' C S' C-C S'' C. Normal compliance uses dimensionless
scaling D=sqrt(E_reference) and fixed-size Cholesky solves. Isotropic interval
admission uses exact volumetric/deviatoric compliance eigenvalue bounds,
including auxetic sources. Orthotropic admission uses a sufficient scaled
Gershgorin bound over every temperature; it can refuse a valid non-dominant
path. For these inverse-compliance paths, positive endpoint tensors alone do not
establish interior stability.
Full `Mandel21` equilibrium sources instead admit an independent affine slope
for each original stiffness coefficient, including normal/shear couplings that
change without a common modulus scale. The existing source-path checks visit all
original knots. SPD endpoint tensors then establish SPD throughout their affine
interval by convexity; no diagonal-dominance restriction is imposed on this
direct stiffness family. The same Mandel rotation carries the complete derivative
into the native projected free energy, entropy and internal-energy response.
Each of the 21 coefficients of every symmetric Maxwell branch may have an independent
affine slope, including signed cross terms. Positive piecewise-linear
relaxation times remain supported. Reference density and strain validity
remain fixed.

`AffineTensorThermal` is the strict affine subset of that same source admission;
consumers retaining one constant derivative still refuse independent Young/shear
paths. A general derivative snapshot is tied to its resolved temperature and
must be refreshed before evaluating another endpoint. Rotated branch slopes
are shared immutably. Nonlinear equilibrium refresh uses fixed-size ratio
and matrix arithmetic without source queries or branch allocations. Anisotropic
paths cache their Mandel rotation once; rotated first/second derivatives and
frame-invariant operator bounds reach the same projected thermal potential.
The strict affine source adapter retains its original common-scale admission
and numerical limits.

`ProjectedTensorProny::thermal_state` evaluates Helmholtz free energy psi,
entropy S=-psi_T, internal energy U=psi+T*S and mechanical heat capacity
dU/dT=-T*psi_TT using the original strain map
and complete native factor memory. Observable Maxwell strain recovery shares
the native checked transfer inverse. Varying rank-deficient factors are admitted
when their original nonzero rows are independently recoverable and the thermal
derivative acts entirely within that physical strain range. A derivative requiring
hidden strain refuses even at zero actual memory. Unchanged tensors need no
strain reconstruction. The expansion evaluator includes equilibrium and every Maxwell spring.
The fixed-size TensorThermalExpansion descriptor admits affine alpha(T) and
carries exact-temperature alpha/alpha' snapshots with whole-interval bounds.
Its source authority comes from the existing integrated expansion owner;
it does not replace the free-strain integral. Curvature includes
-alpha':C_eq:e and -sum alpha':C_j:q_j even where alpha itself is zero.
Constant singular branches evaluate this work from their original factors,
without inversion, under an observable |alpha':C_j:q_j| <= 2L ||C_j alpha'||
envelope. Varying branches retain physical spring-strain reconstruction and
the 2L envelope; the advertised bound includes the same roundoff allowance.

The shared whole-interval capacity bound includes nonlinear stiffness and
expansion curvature. The receiving caloric owner must establish positive
total capacity. Mechanical capacity is not an extra caloric or dissipative
heat source. Unsupported source paths, expansion derivative jumps and
finite-strain thermodynamics remain outside this admission.

`ProjectedTensorProny::thermal_stiffness_derivatives` contracts the same admitted
source derivatives through the original strain columns and physical volume:
`V S^T (dC/dT) S`. It returns every equilibrium and Maxwell cross term in
original branch order under explicit coordinate/work budgets and cancellation.
Signed or indefinite derivatives are not positive-semidefinite storage and are
never Gram-factored. Kelvin viscosity and relaxation-time derivatives carry no
stored free energy.

A consumer retaining physical equilibrium/Maxwell spring strains can evaluate
entropy directly from these complete derivatives, including a changing
rank-deficient Maxwell tensor. This does not relax the existing
`thermal_state` refusal when a consumer retains only factor memory and cannot
recover the required physical strain. Source identity, frame, affine interval
admission and thermodynamic accounting remain the consumer's responsibility.

`transfer_factor_memory` supports changing rank-deficient Maxwell factors with
the same observable strain space, including bulk-only and shear-only arms.
It normalizes every original nonzero Gram row independently and constructs a
checked row-space inverse; only exactly zero original rows are absent. Thus
even a very weak positive physical direction survives without an eigenvalue
cutoff or a fabricated elastic floor. Old and new row counts and physical
range projectors must agree. Changing rank/nullspace, unresolved dependent
rows, nonrepresentable recovery and memory outside the retained factor range
refuse before publication. This is a bounded extension for independent nonzero
factor rows, not a general redundant-factor pseudoinverse. The same recovered
observable strain supplies affine branch entropy and expansion curvature;
the existing 2L envelope and receiving caloric capacity gate remain binding.
Four focused G0/G1/G3 projection regressions cover rank-five thermal round trips,
independent entropy/capacity contractions, weak positive-mode retention, and
transactional refusal of lost/hidden directions and inconsistent memory.

## Source density retained through bounded phase transport

`phase::transport::resolve_phase_transport_state_with_density_departure` keeps
the original source density even at exactly solid states, while the original
enthalpy/phase chart retains its constant reference density and latent energy.
The caller supplies an explicit finite relative departure bound in [0,1).
Exactly solid classification precedes liquid-only Cp/viscosity queries. Mobile
kinematic viscosity converts to dynamic viscosity using RAW source density;
projection never changes that physical viscosity or conductivity.

`PhaseTransport::transport_density_kg_m3` makes the proposed inertial/transport
density explicit even without a liquid transport tuple. The ordinary phase-flow
owner refuses a raw density departure; a bounded Boussinesq owner must retain
reference mass separately. Unequal phase densities, shrinkage and phase-volume
sources remain unsupported by this fixed-density path. Existing exact source,
context, phase identity and chart-domain checks remain active.

The new material regression covers raw density, exact-solid property admission,
kinematic-to-dynamic conversion and strict departure bounds. It is authored but
unrun in this implementation environment; no source or experimental validation
is inferred from the synthetic fixture.

## Phase transport on changing chart density

`phase::transport::resolve_phase_transport_state_for_volume_change` resolves
the original phase/source bundle when its consumer accounts for physical
volume change. The original enthalpy chart and `MushyPhaseLaw` retain their
temperature, phase, constituent-volume and mobility rules, including admitted
unequal phase densities and sensible expansion. At each queried h, source
density must match the chart's current bulk density exactly, including bits.
No density tolerance, reference projection or replacement claim is introduced.

The returned `ResolvedPhaseTransportState` preserves the original material
receipts and existing `FluidTransport`/`PhaseTransport` interfaces. Kinematic
viscosity converts using that same current density. Exact-solid classification
still precedes liquid-only Cp/viscosity requirements; a solid has no fabricated
fluid tuple. Original chart, mobility, source-domain and identity refusals remain
active. Agreement at one sampled state does not establish source/chart agreement
throughout an interval.

`resolve_phase_transport_state_with_shear_for_volume_change` additionally
owns the native physical `sqrt(2 D:D)` rate in 1/s on an existing typed
coordinate. Original h-derived temperature, fraction and current density share
one source snapshot with rate-dependent viscosity, conductivity and admitted
extra responses; pinned selections and fixed pressure/composition remain intact.
Kinematic viscosity converts with that same actual density. Exact solids require
zero rate and no liquid-only property claims. Shear-dependent density or heat
capacity still refuses even at a matching sampled value, because the original
caloric chart does not supply that storage law. The entry introduces no
constitutive iteration, independent phase state or volume-flow solve; original
fixed-density shear entries retain their admission.

The consumer must conserve mass and extensive enthalpy and account for occupied
volume `M/rho` in geometry or continuity. Pressure work and phase support remain
with their physical owners; this resolver supplies no EOS pressure update or
momentum solve. Existing fixed-density, shear and bounded density-departure
entry points retain their previous admission, and fixed-volume phase flow still
requires its original constant-density condition.

## Affine eigenstrain in source tensor memory and thermodynamics

`ProjectedTensorProny` retains explicitly admitted stress-free Mandel strain in
the original quadrature frame. Equilibrium affine load/constant terms represent
one physical elastic potential; original branch Gram factors lower free strain
to the same native memory coordinates. Temperature continuation subtracts
sqrt(V)*B_new*delta(e_th) after original factor transfer, preserving viscous
strain instead of resetting spring memory. Unchanged singular factors remain
supported without inversion. Kelvin eigenstrain-rate coupling is not inferred.

For constant stored tensors and expansion derivative alpha, native thermal
state evaluation includes equilibrium V*alpha:sigma and the corresponding
volume-scaled branch-factor entropy. Internal energy remains U=psi+T*S. The
complete equilibrium plus Maxwell tensor capacity coefficient includes even
constrained directions absent from a reduced displacement basis. Source-driven
coupling admits positive total fixed-strain caloric capacity on every supported
chart segment; material uncertainty and source coverage remain separate.

The affine-stiffness/constant-alpha combination uses
`affine_expansion_capacity_bound_j_k2` on the original projected material
point and source derivative. At fixed total and viscous strains its curvature
is V*[alpha:C_inst:alpha - 2 alpha:C_eq':e - 2 sum alpha:C_j':q_j].
Absolute endpoint matrix bounds cover the complete affine temperature interval,
with the original Mandel limit L for equilibrium elastic strain and 2L for
each varying Maxwell spring. The shared thermal-state evaluator enforces that
memory envelope whenever expansion consumes it; accepted, rebound and final
states use the same check. Constant tensors retain their original coefficient,
and zero alpha retains the original affine-C state admission.

Focused coupling regressions retain nonzero and singular Maxwell memory,
rotated expansion, zero-alpha compatibility and source/domain refusals. They
are authored but unrun here; syntax review and independent physical-coordinate
arithmetic do not establish native runtime or experimental validation.

### Shared intrinsic permeability in continuum hydraulic consumers

`DarcyGasBarrier::continuum_geometry` borrows the original ordered intrinsic tensors, source receipts, compiled sum(L/k_effective) and intersected T/p support for native liquid or other continuum hydraulic owners. It constructs no gas state and introduces no second Darcy or transport law. Consumers use their actual viscosity, density, superficial area and solved flow, and admit both absolute face pressures and local temperature against the same original source limits. The original barrier's explicit transverse constraint is retained. Zero intrinsic permeability remains exactly sealed.

Every nonzero gas-slip calibration and every prescribed gas-only wet-mobility state refuses this view, including unit-mobility endpoints. Full single-fluid pore occupation is the admitted case; liquid relative permeability, capillarity, changing saturation, skeleton deformation and storage are not inferred. Native tests cover original rotated/layered geometry, unchanged source receipts and domains, signed speed limits, exact seals, and gas-calibration refusal.

### Temperature-dependent intrinsic permeability

`porous::permeability::temperature::PermeabilityTemperatureLaw` binds three exact
principal permeability keys in square metres to the original material card,
complete selected source bundle, and explicit orthonormal material frame.
Temperature refresh changes only the admitted absolute-temperature coordinate;
pressure, other context, quantity meanings, and every selected claim stay pinned.
It refuses unsupported source continuation instead of silently selecting another
claim. Repeated keys explicitly declare equal principal values, and zero remains
an exact sealed direction. This source owner retains the existing
`IntrinsicPermeability` tensor and its normal-gradient convention.

Consumers supply their actual fluid viscosity, momentum discretization and heat
account. The law does not infer gas slip, relative permeability, evolving pore
geometry, saturation or skeleton heat storage. A consumer whose native operator
cannot represent oblique or sealed directions must refuse them explicitly.
Focused synthetic regressions exercise rotated temperature-dependent tensors,
fixed pressure and quantity metadata, whole-bundle source pinning, exact seals,
and immutable refusals; native test execution remains pending in this environment.

### Source-bound separable Forchheimer coefficients

`PermeabilityTemperatureLaw::with_forchheimer` attaches three exact beta keys
already present in the SAME original resolved material bundle. The static
dimensional quantity is inverse metres. Beta is independently supplied,
finite and nonnegative; neither intrinsic permeability nor an Ergun relation
invents it. The convention in each declared principal direction is
`f_i = -rho beta_i |u_i| u_i`, with superficial velocity and NO factor one half.
This is an explicitly calibrated separable directional closure. Even equal
beta components do not make it the rotationally invariant vector law
`-rho beta |u| u`. The consumer must represent the declared material axes.

Temperature continuation preserves all selected claims, auxiliary conditions,
source requirements and the explicit positive finite maximum magnitude of each
superficial velocity component. The caller must supply a calibration applicable
to the chosen fluid and specimen. Consumers enforce the speed support on their
actual velocity degrees of freedom; they must not admit an over-limit face by
averaging cancelling velocities. Zero beta remains exactly zero. This source
owner introduces no momentum solver, pore-speed conversion, thermal store,
turbulence, saturation or compressibility closure. The native fs-lbm coupling
uses the original face force and its separate cubic physical heat account.

## Conservative frozen-gas endpoint coordinates for native moving walls

`ClosedGasParcel::at_volume_internal_energy_increment` reconstructs a complete
source state from a positive represented volume and an actual internal-energy
increment. The original temperature-anchored caloric inversion and EOS own the
result. Parcel mass, composition, source model and auxiliary conditions remain
unchanged, and original-law increments cancel common formation datums before
rounding. Unsupported source states, unrepresented coordinates, cancellation
and energy mismatch return no changed parcel.

`ClosedGasParcel::stored_energy_change_to` compares two admitted endpoints of
that same inventory directly through the original species caloric increments.
It cancels common formation references before summation and permits volume and
temperature changes while requiring exact mass, complete source/composition
model and auxiliary-query identity. The original caloric owner admits the full
temperature interval. Changed sources, inventory or auxiliary conditions refuse;
cancellation publishes no changed parcel. This is a state difference, not a heat
input, an adiabatic-path assertion or a new thermodynamic closure. Focused
regressions cover temperature-dependent calorics beside a large formation datum,
reversal, volume changes and source/inventory/context mismatch.

The volume/energy method is a coordinate adapter, not an extra heat/work integrator. Its
enclosing mechanical/thermal owner must close actual gas energy against its
native boundary work and heat. The deformable-cavity owner in fs-couple does so
using the original tetrahedral pressure footprint and retained Newmark state.
The source regression checks isothermal volume/pressure, exact mass and context,
invalid energy/volume and immutable refusal; native execution remains pending.

## Source-bound Kelvin tensor viscosity

`visco_tensor::source::ResolvedTensorMaterial::with_kelvin_viscosity` binds the full
symmetric Mandel rate tensor from 21 explicit Pa s property names on the same
card and physical query. Missing entries are not inferred as zero. Original
claims, external-coordinate meanings and intersected frequency domains remain
attached through refresh and rotation. The existing tensor law supplies PSD
admission; projection retains its actual Gram rows for physical rate forces
and nonnegative heat without refactoring a rounded reduced matrix.

Thermal feedback admits validity-wide constants or piecewise-linear temperature
curves, checking every knot of their common refinement with a bounded total of
4096 original knots. The existing material owner's PSD arithmetic policy applies
at every union knot; convexity of the PSD cone admits the intervening linear
segments. Endpoint-only sampling, extrapolation and nonlinear coefficient
interpolation do not establish this admission. Both constant-storage thermal
relaxation and affine stored-tensor thermal owners retain this rate law. Kelvin
viscosity contributes no storage or mechanical entropy; its consumer owns the
actual strain rate and heat. The rigid retained-layer consumer applies it to
recoverable shear before series plastic flow, at the same trial temperature as
the other material kinetics.

## Oriented reciprocal thermoelastic material storage

`elastic::thermal::AffineThermoelastic` carries independently selected
isotropic or oriented orthotropic stiffness and expansion through the existing
resolved-response owner. Source-pinned affine Young's moduli, Poisson ratios,
shear moduli and instantaneous expansion coefficients retain separate fixed
proper frames. The native compliance inverse supplies stiffness, with exact
`C'=-C S' C` and `C''=2 C S' C S' C-C S'' C`; affine engineering constants
are not treated as affine stiffness. Original source integration supplies the
full free-strain tensor.

`ThermoelasticSymmetry::IsotropicBulkShear` also accepts independently affine
bulk and shear modulus sources directly. In Mandel coordinates the exact law
is `C=3 K P_vol+2 G P_dev`, with constant `C'` and zero `C''`; it does not
synthesize Young/Poisson claims or treat their rational transforms as affine.
Expansion remains independently isotropic or oriented orthotropic. Original
positive source bounds and conservative whole-interval capacity admission
still apply. Bulk/shear is an elastic source schema and refuses as an expansion
schema. Independent invariant-energy/entropy/capacity and rotation regressions
cover this path.

`ThermoelasticSymmetry::IsotropicYoungShear` retains independently affine E/G
claims with the exact rational bulk response `K=E G/[3(3G-E)]`. Its first and
second derivatives, including nonzero `C''`, enter the same thermodynamic
potential. Scaled ratios avoid cubing dimensional moduli. Whole-interval
admission requires positive E, G and `3G-E`; outward-rounded derivative bounds
feed the original combined heat-capacity check. Expansion remains independent,
and no derived bulk or Poisson claims are inserted. The same regression seam
checks both source schemas, with the E/G oracle differentiating volumetric
compliance independently and refusing an interval that crosses the elastic pole.

The same potential supplies stress, entropy, internal energy, fixed-total-strain
heat capacity and stress-temperature derivatives to the original spatial solid,
finite-heat, inertial and adiabatic consumers. Full-interval positive-compliance
and capacity admission uses outward-rounded cofactor and entrywise matrix
bounds. Conservative refusal is possible for wide or poorly conditioned source
intervals; sampled endpoints alone do not establish admission. The original
isotropic arithmetic is retained. This is smooth affine-source, fixed-frame,
fully solid small-strain thermomechanics, without phase evolution or plasticity.

## Molecular Rayleigh scattering on shared gas states

`gas_scattering::MolecularRayleigh` admits an explicitly sourced reference
refractivity (`n - 1`), molecular number density, visible wavelength and
unpolarized 90-degree depolarization ratio. The Lorentz–Lorenz cross section
uses the matching King correction. No species name supplies missing optics.
`MolecularRayleighMixture` retains the original `GasSpec` rows, shared mass/mole
composition and explicit water label, then delegates thermodynamics and
transport to the original `GasState`/Wilke/WMS owners. Molecular cross sections
add by number fraction; normalized phases add by scattering weight. The
prepared coefficient is `N_A * sum(x_i sigma_i) / M_mix` in m²/kg, so the actual
gas density supplies m^-1 attenuation and responds to pressure and temperature.

The same references also supply coherent refractivity through
`RayleighGasState::refractive_index[_at_density]`. The prepared Lorentz–Lorenz
coefficient is `K = (N_A/M_mix) sum_i x_i L_ref,i/N_ref,i` [m³/kg], with
`L_ref = (n_ref²-1)/(n_ref²+2)`. The local phase index is
`sqrt((1+2 rho K)/(1-rho K))`; its exact density derivative is returned for
smooth-field ray bending. Polarizabilities add before this coherent-index
evaluation, independently of the cross-section sum used for scattering, and
no King factor changes the phase law. `n-1` is separately retained without
subtractive cancellation. Vacuum is admitted; negative/nonfinite density and
local `n-1 > 0.01` refuse before the Lorentz–Lorenz pole. This boundary is an
explicit dilute-use policy, not measured accuracy. Phase and group indices
coincide only under the existing constant-polarizability approximation.
The v2 state identity includes this model's semantics while retaining the
original sources, composition and accepted gas state.

`refractive_index_jet_at_density` additionally supplies the analytic second
density derivative for spatial ray-bundle focusing. With `L = rho K`, it uses
`n'' = (n')² (1 + 8 L)/(3 n)` and retains the original first-order result.
Unrepresentable nonzero curvature refuses this optional query without narrowing
the first-order phase law. Spatial Hessians remain the reconstruction owner's
responsibility. G1 regressions compare an independent second derivative of the
closed-form index, finite differences of its first derivative, vacuum, null
polarizability and the dilute-domain limits.

`gas_scattering::composition::MolecularRefractionModel` prepares component
polarizabilities and original molar masses from an accepted
`MolecularRayleighMixture`/`GasState` pair or `ThermochemicalRayleighState`.
Thermochemical species retain canonical order; ideal-gas rows retain original
order. Nonnegative partial mass densities [kg/m³] or molecule densities
[molecules/m³] enter `L = sum_i rho_i K_i = sum_i N_i v_i` before the shared
nonlinear Lorentz–Lorenz phase law. The result retains `n-1`, `dn/dL` and
`d²n/dL²`; no fraction normalization or division by total density occurs.
All-zero rows are vacuum. Empty/mismatched counts, invalid densities,
unrepresentable positive contributions and out-of-band/dilute-domain queries
refuse. Preparation retains the original 128-component bound; evaluation uses
bounded linear work and allocates no memory. The identity binds original
model sources, component optics, masses and order. Varying concentrations do
not establish new thermodynamic acceptance, absorption or new species.
Constant polarizability remains the default; explicit dispersion sources may
be attached as described below. G1/G3 fixtures
cover original ideal/NASA state recovery, mass/number basis equivalence,
composition sensitivity at fixed total density and analytic nonlinear
derivatives.

`composition::MolecularCauchyDispersion` supplies an explicitly sourced Cauchy
refractivity law anchored at each original molecular reference number density:
`r(lambda)=r_ref+B(lambda_um^-2-lambda_ref_um^-2)+C(lambda_um^-4-lambda_ref_um^-4)`.
Finite nonnegative shape coefficients, nonzero source identity and an explicit
source interval inside 380–780 nm containing the reference are required.
Monotone endpoint admission covers the dilute domain throughout that interval.
`with_dispersion` requires original-reference/order/count agreement, intersects
source support and binds sources and coefficients independently of attachment
history. `at_wavelength` prepares additive coefficients and analytic
`lambda dK_i/dlambda` once, sharing reciprocal powers across species and keeping
original molar masses. Phase and spatial derivatives use the full Lorentz–Lorenz
sum; group index is `n_g=n-n_L sum_i density_i lambda dK_i/dlambda`, with a
separately retained cancellation-safe group refractivity. Vacuum and zero-shape
fits preserve constant-law results. Extrapolation and unrepresentable values
refuse. Three material regressions cover reference anchoring, frozen constant
bits, density-basis equivalence, independent spectral differentiation and
source/order/domain boundaries. This declared normal-dispersion fit adds no gas
presets, measured accuracy, anomalous/resonant dispersion, dense-gas correction
or absorption.

`MolecularRefractionWavelength::scattering` derives elastic molecular scattering
from the same source-bound molecular volumes prepared for phase/group refraction.
Original component order, density basis, molar masses, source support and full
dilute/representability admission remain binding. Independent King-corrected
species cross sections add as `beta_s = sum_i N_i sigma_i`; the coherent mixture
polarizability is never squared to obtain scattering. The result retains physical
mass density and the scattering-weighted depolarization matching the molecular
angular law. Relative species dispersion therefore changes the effective phase
with wavelength. Evaluation uses bounded stack storage and preserves original
molecular arithmetic for absent or zero-shape dispersion. Nonnegative componentwise
density bounds conservatively bound the collision rate at each prepared wavelength.
Vacuum has zero rate/depolarization; optically inert components still contribute
mass. Positive underflow refuses. Four G1/G2 regressions cover an independent
electric-dipole/phase-mixture oracle, density-basis and original-law parity,
vacuum/inert mass/refusals and componentwise majorants. Species anisotropy
dispersion, correlated scattering, absorption and emission remain outside scope.

Preparation is bounded by the original 128-component gas limit. The original
`MolecularRayleighMixture`/`RayleighGasState` spectral and phase queries retain
allocation-free constant work and their exact lambda^-4 law under constant
molecular polarizability on 380–780 nm. The composition-prepared extension above
uses bounded linear work for explicit dispersive sources. Neither law adds
resonance, Raman shift, dense-gas correction, Stokes transport, absorption or
thermal emission. The original gas owner still owns phase and
thermodynamic validity. Source hashes retain caller input identity, not
measurement authority. `bind_gas_state` checks complete numerical coherence;
the caller must retain composition identity because a bare `GasState` has none.

G0/G1/G3 regressions cover the independent dipole cross-section reference,
King correction, exact wavelength law, molecular addition, phase normalization,
mass/mole conversion, pressure/temperature scaling, original-state binding and
invalid-source/domain refusal. `fs-render::material_state::gas_scattering`
consumes the same coefficients and phase without a sampled spectral fit.

`gas_scattering::thermochemical::ThermochemicalRayleighState` attaches one
explicit optical source to each active canonical species of an original
`ResolvedThermochemicalGasState`, including direct `AxialGasField` samples.
The complete accepted owner is retained. Molecular cross sections and
Lorentz--Lorenz polarizabilities use its original mole fractions and evaluated
molar mass without replaying a calorically perfect `GasSpec` or another gas law.
Missing, duplicate and foreign optical rows refuse; input ordering is immaterial.
The shared molecular reduction supplies both scattering and refractivity.

`RayleighGasState::model_identity` separates the frozen composition and optical
model from its local thermodynamic state. Thermochemical source packs and
transport cards belong to that model; selected NASA regions and evaluated
caloric/EOS/transport receipts remain in the complete state identity. A field
can therefore cross a NASA temperature-region boundary while retaining the
same declared optics. Native cell-grid binding borrows these accepted states,
checks exact density coherence and rejects composition/source changes.
The declared constant-polarizability approximation and optical omissions above
remain. Focused tests cover original-state retention, source-order matching,
NASA region transitions, independent molecular sums and gas-field attachment.
Refractivity regressions cover reference recovery, the independent density
derivative, coherent mixture weighting, pressure/temperature dependence,
dilute precision, vacuum and domain/pole refusal. The renderer's explicit
density-and-gradient adapter keeps the original molecular law and snapshot
identity; it never treats piecewise-constant cell averages as smooth gradients.

### Finite mixed-gas condensation and re-evaporation

`phase::condensation::{WaterCondensationModel, WetGasParcel}` retains one
source-backed NASA-9 gas, its canonical species inventories and a finite
liquid-water store. At prescribed temperature and fixed gas volume, the
original IAPWS saturation pressure determines the vapor capacity through
water's **partial** pressure. Condensation retains the removed water as liquid;
re-evaporation cannot exceed that liquid inventory. Noncondensing species,
gas volume, auxiliary material coordinates and original source packs survive.
The original gas owner re-resolves density, pressure, sound speed and transport.
Positive phase inventories and their transferred increments must remain
representable; a large film cannot hide a lost transfer in its mass tolerance.

The caller explicitly identifies water and chooses a reference temperature
aligning NASA water-vapor enthalpy with the IAPWS vapor datum once. H2O molar
mass must agree within 10 ppm. Liquid enthalpy and internal energy carry the
same constant shift. The combined internal-energy increment cancels formation
references analytically before evaluating sensible and latent contributions;
the reported latent enthalpy is a diagnostic, never an additional heat source.
A heat-transfer consumer must solve this combined increment against actual
net heat before publishing gas, liquid and thermal partners together.

This approximation admits 273.16--323.15 K and 20--200 kPa total pressure.
Liquid properties use the saturated-liquid approximation, with no gas-volume
displacement, humid-air fugacity/enhancement/Poynting/Henry correction, ice,
nucleation kinetics or interface transport. No measured error bound or full
humid-air thermodynamic-equilibrium certification is asserted. Tests in
`phase/condensation/tests.rs` cover G0 water/carrier/energy conservation and
finite-film exhaustion, G3 invariance to large formation-reference changes,
G4 cancellation/source/domain and unrepresentable-inventory refusal, and G5
repeated deterministic cycles.

`ThermochemicalPhase::LiquidWaterIapwsLimit` carries that same saturation
reference through thermochemical records and `AxialGasField`. The legacy
`LiquidWaterLimit` retains Buck unchanged. The IAPWS policy accepts only its
declared near-atmospheric domain and 256 machine epsilons of partial-pressure
reconstruction roundoff; real excess water and unsupported liquid phases
refuse. It is admission, not phase conversion. Interpolation between admitted
stations is rechecked and can still refuse supersaturation; no implicit clamp
discards water or latent energy. `gas_field/condensation_tests.rs` exercises a
real wet-parcel endpoint through the shared field and those refusal boundaries.

## Dilute aqueous osmosis and selective membrane transport

`osmosis` owns one shared isothermal constitutive law for independent dissolved
particle concentrations in mol/m³. `IdealDiluteSolution` retains finite solute
amounts and solvent volume separately; pure-solvent exchange preserves every
solute amount and recomputes dilution, osmotic pressure, and water activity.
Its pressure is `RT sum(c_i)` and its activity obeys
`ln(a_w) = -Vbar_w sum(c_i)`, with an explicit incompressible solvent partial
molar volume. This is the declared dilute model, not a non-dilute ideal
mole-fraction model. The separate activity-to-pressure primitive accepts
externally supplied nonideal activity through `-RT ln(a_w)/Vbar_w`; it does not
invent a species decomposition for selective transport.

`MembraneTransport` uses distinct hydraulic permeance [m/(Pa s)],
concentration-based solute permeance [m/s], and per-species reflection [0,1].
Left-to-right water flux is `Lp [delta_p - RT sum(sigma_i delta_c_i)]`.
The constant-coefficient Spiegler–Kedem convection-diffusion equation is
integrated across each membrane species path with a stable exponential fit.
Both donor coefficients are nonnegative; zero flow reduces to Fick diffusion,
zero solute permeance reduces to upwind solvent drag, and full reflection
excludes drag. A nonselective interface creates no osmotic water pumping.
Geometry is supplied by consumers and is not counted twice. The explicit
Darcy-slab crosswalk consumes thickness, intrinsic permeability and the original
fluid owner's current dynamic viscosity rather than confusing permeability,
conductivity and permeance.

The prepared 0–64-species kernel and finite-inventory replacement allocate no
memory; fixed-order compensated sums and the existing deterministic elementary
functions define repeatability. Invalid shapes, negative inventories and
unrepresentable outputs refuse atomically without publishing partial results.
G0/G1/G3 regressions cover independent-particle and SI-unit oracles, activity,
forward/reverse osmosis, selective mixtures, exact membrane integration,
orientation, finite dilution, actual viscosity, and positive-concentration
dissipation. Zero concentration has a finite transport limit but no finite
ideal chemical-affinity claim. Consumers own integration, conservation,
cancellation, pressure/storage work and any thermal balance.

There is no ionic chemistry, Donnan potential, solute partial-volume, phase
stability, transient concentration-layer storage, or experimental-validation
claim. Raw coefficients carry no material-card authority; source-backed users
must retain the original coefficients' full solution and operating support.


### Finite saturation and retained solid solute

`osmosis::solubility::FiniteSolubilitySolution` composes the original L1
`IdealSolubilityV1` pair laws with the shared dissolved-particle thermodynamics.
It retains 1–64 explicit neutral components sharing one solvent, each with a
current total and separate dissolved/solid amounts. Only dissolved particles
contribute to osmotic pressure and water activity. Pure-solvent removal can
precipitate excess solute; dilution redissolves the same finite solid stock,
including complete solid exhaustion. The dissolved-only `AqueousSolution` API
keeps its existing meaning and is a view, not a substitute for the phase owner.

Fixed-temperature pair evaluations are cached in shared immutable components.
Each solvent-volume trial performs direct finite partitions without repeated
temperature exponentials or chemical root solves. Explicit retemperature
re-evaluates every original law within its domain. Original laws, identities,
totals, capacities and floating-point phase balances remain inspectable.
Optional distinct solid/dissolved heat capacities remain part of that same
law: their integrated temperature correction reaches membrane and moisture
consumers through the original cached saturation, without another phase model.
G0/G3 cases check two-component precipitation/redissolution, original ideal
unsaturated and zero-stock limits, prescribed heating/cooling, phase balances,
species-layout/domain refusal and immutable failed trials.

Transport consumers can prepare an `independent_clone` once and then atomically
replace solvent volume and TOTAL component inventories in the same species
layout. The bounded replacement reuses all prepared buffers and cached pair
coefficients, admits every candidate phase and dissolved thermodynamic state
before mutation, and refuses shared stocks instead of allocating a hidden copy.
A borrowed `partition_at_inventory` evaluates transport trial stocks without
altering the accepted solution. Subsequent dilution or retemperature retains
the updated totals. Heap accounting includes phase storage, original pair
identifiers and their retained allocation capacity. Focused regressions cover
phase crossing, continuation with transported stock, stable buffer addresses,
and unchanged state after invalid, shared-storage or water-activity refusal.

The model assumes independent nonreacting neutral solutes, instantaneous
equilibrium, negligible solute/solid volume and an incompressible solvent.
It does not supply electrolyte dissociation, common-ion or pH chemistry,
nonideal activities, mixed crystals, hydrate water or nucleation/rate laws.
Energy-balanced updates require the explicit caloric adapter below.
Caller-declared pair parameters retain their original numerical
authority; the constructor does not authenticate measured solubility data.


### Original material sources for neutral saturation and film kinetics

`osmosis::solubility::source::NeutralSolubilitySource` admits an explicitly
pinned, full-content model identity from an immutable material card. Version-one
`neutral-solubility-vant-hoff` supplies exactly reference absolute temperature,
saturation in mol/m³ solvent and dissolution enthalpy in J/mol. Original typed
temperature support and singleton exact solvent/solute qualifiers are required;
foreign pairs, unsupported conditions, missing sources, wrong units and
extra/missing parameters refuse. Both the complete material card and original
model remain inspectable. This binds an authored/calibrated source convention;
it does not create measured evidence or supply phase heat capacities.

Its independently sourced `NeutralSolubilityKineticsSource` admits prescribed
fixed-film conductance with an explicit precipitation policy, or conductance at
a declared fixed-count particle reference stock. kA retains its apparatus and
geometry basis and is never inferred from equilibrium, density, bath volume or
a material name. Native kinetics use the intersection of both temperature
domains; if needed, the original constant-enthalpy saturation curve is
numerically rebased inside that intersection without extending source support.
Existing material, caloric and transport owners consume the same native laws.
Callers retain the source binding alongside numeric owners to retain the
original receipts. Focused cases cover drying, heating, finite-bath particle
dissolution/growth, source-declared growth surfaces, exact pair/units and refusal
of unsupported source domains.


### Energy-balanced finite material solutions

`FiniteSolubilitySolution::add_heat` connects the current material reservoir to
the original L1 `SolubilityCaloricModelV1`. It requires exact agreement on solvent
identity, incompressible molar volume, ordered components and every original
solubility-law parameter. It prepares the live dissolved and TOTAL inventories,
including stocks already changed by membrane transport, and preserves the
original solvent volume exactly. Positive or negative heat is closed by the
shared thermal owner; the returned material solution is repartitioned at that
solved temperature using the same original laws and stocks.

`FiniteSolubilityHeatStep` carries both the accepted material endpoint and the
original first-law balance. Aqueous concentrations, osmotic pressure and water
activity therefore follow the new temperature and dissolved amounts. The
operation borrows its input and publishes nothing on incompatible sources,
unsupported heat, original material refusal or late cancellation. Ownership
can be moved out of the result for subsequent material evolution without
copying the accepted phase inventories.

This adapter adds no caloric constitutive law or transport integrator. It uses
either the L1 owner's explicit common phase Cp or the exact distinct phase
curves already bound to its temperature-dependent solubility law, instantaneous
independent neutral-solute equilibrium and fixed-pressure, fixed-volume
incompressible-solvent convention. It does not claim simultaneous thermal
membrane transport or variable-pressure work. Focused regressions check
independent sensible-plus-dissolution heat across exhaustion and recrystallization,
unequal-phase caloric feedback into osmotic pressure and solvent activity,
updated transported stocks, exact solvent volume, consistency of both accepted
views, source incompatibility, unsupported heat and transactional cancellation.

### Finite-rate soluble material inventories

`osmosis::solubility::kinetics::KineticSolubilitySolution` retains 1–64 original
ordered kinetic laws and actual dissolved/solid stocks. Its aqueous view uses
only the current dissolved amounts for osmotic pressure and water activity;
the instantaneous-equilibrium target is a separate diagnostic. Construction
accepts explicitly unequilibrated stocks or the complete original equilibrium
owner through `from_equilibrium`, with caller-supplied conductances and growth
policies. Original phase heat-capacity bindings remain attached to each law.
The `new` constructor also retains each original law's opt-in particle surface:
fixed particle count and self-similar shape give `G(s)=G_ref*(s/s_ref)^(2/3)`.
The original reference solid stock remains unchanged during drying, dilution,
temperature changes and transport; actual crystals shrink or grow with their
local inventory. This option requires existing solid and supplies no nuclei.

Prescribed solvent-volume and temperature changes retain both phase stocks.
Only subsequent kinetic steps cause dissolution or precipitation. Repeated
steps reuse prepared phase/aqueous buffers and bounded stack storage, with no
heap allocation. Every original component and the resulting solvent activity
are admitted before any publication; early, intermediate and final cancellation
leave all accepted inventories unchanged. Clones share immutable laws and own
independent writable state.

`replace_dissolved_inventories` admits transported solvent/dissolved amounts
while retaining every original crystal stock. `copy_state_from` transfers
complete compatible receipts into prepared buffers without allocation or
repartitioning. Explicit two-phase construction retains donor tails below one
ulp of the dominant phase; immutable source payloads and Arc control storage
remain in the preparation memory estimate. Native membrane transport uses this
boundary and evolves the original kinetic owner on its adaptive timeline.

The basic `advance` operation is isothermal and infers no thermostat heat,
particle breakup/agglomeration or spontaneous nucleation. Nonideal phase
exchange requires the explicit original thermodynamic binding below. The
coupled thermal operation separately accounts for heat and temperature.
Focused cases verify opposed phase transfers, actual osmotic response,
finite-rate dilution/heating lag, retained crystals through transport,
original source sharing, buffer reuse, layout refusal and cancellation.
Particle regressions check drying-induced growth followed by analytical
redissolution after dilution, evolving conductance, unchanged reference geometry
through dissolved transport and refusal of incompatible geometry checkpoints.

### Activity-consistent neutral crystallization

`KineticSolubilitySolution::with_nonideal_thermodynamics` binds the original
`VirialSolutionModel` to the SAME retained kinetic owner. Temperature, solvent
molar volume and ordered species count must match, and every original pure-pair
saturation must lie inside that model's declared convex concentration box.
The original kinetic laws, particle reference stocks and precipitation policies
remain unchanged. A different source cannot replace an existing binding.

The aqueous-minus-pure-solid chemical potential is
`RT*[ln(c_i/csat_i) + (B*c)_i - B_ii*csat_i]`. Thus each original binary
saturation remains exact, while other dissolved species shift conditional
saturation through their explicitly supplied cross interactions. No measured
mixture calibration, association or concentrated-solution applicability is
inferred from the pure-pair source. This neutral model retains separate pure
crystals; it does not supply ionic or mixed-solid equilibria.

Each declared-order coordinate stage holds the other current concentrations
fixed, solves the monotone conditional saturation inside the original box and
passes that concentration to the original thermochemical `advance_at_saturation`.
The same finite-bath and changing-particle-area trajectories perform all phase
movement. The root is analytic for zero diagonal interaction; otherwise a
safeguarded Newton/bisection solve permits at most 128 iterations. Unsupported
or unrepresentable roots refuse without a concentration floor or extrapolation.
After preparation, advancement uses existing buffers and bounded stack storage
without heap allocation. This sequential phase split is first order; an outer
full/two-half controller must control cross-component temporal error. Receipts
retain each coordinate stage's actual target, transfer and exhaustion event,
not a simultaneous equilibrium claim after later coordinates move.

`solid_relative_chemical_potentials_into` and
`solid_relative_inventory_free_energy_change_j` use that same calibration. The
work combines the original stable ideal primitive, the symmetric virial excess
difference and the pure-solid reference term, including changing solvent volume.
The calibration is linear in component moles and contributes no pressure at
fixed amounts, so the original virial osmotic pressure remains consistent.
Complete dissolved-plus-solid component and solvent conservation must be checked
separately before cancelling common reference terms across a closed network.

Cloning, compatible copying, transported inventories and dilution retain the
bound model and actual phase stocks. A fixed-temperature-only virial source
supplies neither temperature-dependent interactions nor excess calorics, so
temperature changes and `advance_with_heat` refuse unless the explicit thermal
family below has been admitted. Numerical source applicability remains separate
from measured material authority.

Focused native tests check unchanged binary particle trajectories, positive and
negative cosolute shifts, two-solute refinement against independent continuous
RK4 trajectories for both surface laws, chemical-work and pressure derivatives,
reference invariance, trace roots, domain refusal and atomic rollback.
`cargo run -p fs-material --example cosolute_crystallization` demonstrates
cosolute-induced growth followed by complete redissolution after prescribed
dilution. Its coefficients are explicitly authored demonstration data.

### Temperature-dependent activities and crystallization heat

`kinetics::VirialSolubilityThermalModel` explicitly binds the original caloric
model to a symmetric interaction family
`B(T)=B_ref+C*(1/T-1/T_ref)`, with `C` in K m³/mol, one declared temperature
interval and the original concentration box. The supplied reference slice is
retained exactly. Convexity checks on both endpoint matrices qualify the affine
inverse-temperature family; they do not establish empirical mixture validity.
All pure-pair saturation curves must stay inside the original concentration
support throughout the thermal interval.

The original solid caloric datum is retained. The calibrated AQUEOUS excess
potential is `RT*[n^T B(T) n/(2V)-sum(n_i B_ii(T) csat_i(T))]`. Its consistent
additional enthalpy is
`R*n^T C n/(2V)+sum(n_i csat_i(T)*(B_ii(T)*delta_H_i(T)-R*C_ii))`.
The linear reference correction is required even for constant nonzero B when
the original saturation changes with temperature. Baseline phase Cp curves
therefore remain reference calorics of this explicit potential; the effective
mixture Cp and dissolution heat include the new correction. Already measured
effective mixture heat capacities cannot be added again as reference Cp.
No temperature-dependent common energy datum is silently cancelled as heat.

Source preparation uses exact extrema of original piecewise-linear Cp curves
and bounded interval arithmetic across 32 temperature subintervals to require
positive effective dissolved Cp throughout the interval. The corrected
derivative includes `csat_i*[B_ii*(delta_H_i^2/(R*T^2)+delta_Cp_i)
-2*C_ii*delta_H_i/T^2]`. Original solid and solvent Cp remain positive.
A source whose thermal branch cannot be qualified refuses; sampled endpoint
Cp alone is insufficient. Direct saturation-log and Kirchhoff increments
retain neighboring-temperature changes in the reference heat, while symmetric
bilinear increments retain small mixture work without subtracting large
absolute enthalpies.

`KineticSolubilitySolution::with_thermal_nonideal_thermodynamics` attaches this
family to the SAME actual phase owner. Original pair laws, reference calorics,
solvent convention and current fixed-temperature slice must match. Cloning,
dilution, transported dissolved inventories and explicitly supported temperature
changes retain the original immutable family. A different thermal slope cannot
replace it merely because its current B matrix happens to agree. Preparation
memory includes the retained caloric source and its actual vector capacities.

The existing `advance_with_heat` adapter passes that SAME family Arc to the
original thermochemical thermal owner and reconstructs the accepted aqueous
view at the solved temperature. The native controller advances conditional
phases at the stage's current temperature and then closes heat at frozen final
stocks. Its existing full/two-half comparison controls the first-order split;
the positive frozen-stock capacity qualifies the scalar heat root. No separate
phase solver, adaptive timeline or replacement kinetic state is introduced.
Whole-window heat is checked directly from original/final phases and temperature.
Fixed-temperature-only fits still refuse heat evolution.

Regression sources cover Gibbs--Helmholtz consistency, constant-B reference
heat, interior thermal instability, adjacent-temperature and trace increments,
independent coupled ODE refinement, particle geometry, retained source identity
and transactional refusal. `cargo run -p fs-material --example
cosolute_crystallization -- --heat` adds adiabatic growth and prescribed
heating/cooling using explicitly authored reference and interaction data.
Thermal windows retain fixed solvent volume and do not simultaneously solve
spatial heat/mass transfer, pressure-dependent solubility or ionic precipitation.

### Finite-rate thermal material evolution

`KineticSolubilitySolution::advance_with_heat` binds the exact original kinetic
models and solvent convention to `SolubilityThermalKineticsV1`. It prepares the
current transported stocks, advances the shared thermal/phase equations and
rebuilds the aqueous view from the accepted actual dissolved amounts at the
solved temperature. No equilibrium projection or second caloric equation is
introduced. The returned material and original energy/temporal report can be
moved together for subsequent evolution. Refusal or cancellation leaves the
input material untouched.

Zero-time phase receipts refresh the accepted endpoint targets; whole-window
heat, residual and adaptive statistics remain on the thermal report. The
underlying method is first order, keeps volume and declared reference conductances
fixed over the heat window, and does not simultaneously transport solvent or solute.
Selected particle geometry changes the actual surface conductance through the
same kinetic owner; exact source matching includes its reference stock.
Adapter regressions cover live transported inventories with an independent
heat/pressure/activity calculation, rate/policy/solvent-source mismatch,
final-publication cancellation/retry and sensible-only heating with trace
crystals. An independent constant-Cp energy balance also covers particle
dissolution and its slower cooling relative to a fixed-area bath.
`cargo run -p fs-material --example solubility_dissolution` demonstrates
fixed-area and shrinking-particle dissolution, concentration-induced crystal
growth, coupled adiabatic finite-rate cooling and a separate equilibrium heat calculation.
All demonstration coefficients are explicitly authored.

### Jacket-driven dissolution and crystallization

`KineticSolubilitySolution::advance_with_heat_transfer` sends a prescribed
reservoir temperature and complete contact conductance in W/K to the SAME
`SolubilityThermalKineticsV1` owner. The native heat root includes the actual
solution-temperature feedback; a frozen starting-temperature power is not
substituted for the boundary law. The accepted thermal report retains the
independently integrated supplied heat and direct original-to-final enthalpy.

The material handoff is shared with `advance_with_heat`: exact original phase
laws, particle reference geometry, caloric sources and optional nonideal
thermal family remain bound. Actual transported dissolved/solid stocks enter
the solve, and the accepted aqueous view is rebuilt at the solved temperature.
Volume and component totals stay fixed over the window. The jacket temperature
is prescribed; finite coolant storage and simultaneous solute transport are
not inferred. Zero contact, source mismatch and late cancellation preserve the
original transactional behavior.

`cargo run -p fs-material --example cosolute_crystallization -- --jacket`
cycles one nonideal particle bath through cooling, a cold hold, heating and
recooling. Its CSV output includes jacket temperature/conductance, actual phase
stocks, conditional saturation, pressure, supplied heat and energy residual.
All heat-transfer and interaction coefficients are authored demonstration data.

### Finite coolant storage during soluble-phase evolution

`KineticSolubilitySolution::advance_with_finite_heat_reservoir` couples its
actual phase owner to an explicit `SensibleHeatReservoirV1` and contact
conductance. It retains the SAME source checks and endpoint reconstruction as
the prescribed-heat and prescribed-temperature methods. The new
`KineticSolubilityReservoirStep` returns the accepted material together with
the original core reservoir outcome, including equal/opposite supplied heat,
both actual endpoint temperatures and the closed-system energy residual.

The finite reservoir has a declared constant complete heat capacity, not an
inferred coolant mass, composition or phase transition. Its compensated energy
inventory and original reference travel with the returned reservoir across
windows. Solvent and solute remain in the solution; heat crosses the contact.
Nonideal excess/reference calorics and evolving crystal surface area remain
on the original material/thermochemical owners.

`cargo run -p fs-material --example cosolute_crystallization -- --finite-jacket`
first connects a cold finite coolant charge, then explicitly replaces it with
a hot one. Each contact continues through four windows with the same reservoir
energy inventory. Output includes both temperatures, phase stocks, each body's
enthalpy change, cumulative coolant energy and the closed energy residual.

## Shared cavitation coefficients and frozen mixture acoustics

`cavitation` owns a checked immutable liquid state at one explicit absolute
temperature: liquid density, dynamic viscosity, surface tension, equilibrium
vapor partial pressure and liquid sound speed. Optional saturated vapor density
and latent vaporization enthalpy must be supplied together; no gas EOS, gas
content or missing coefficient is inferred. Numeric zero-viscosity, zero-surface-
tension and zero-vapor-pressure limits are explicit. Source wrappers retain
exact property roles and the original material/query receipts, convert kinematic
viscosity using the same current density, and refresh all claims atomically.
Temperature-dependent values remain frozen during a consumer advance.
`SaturatedWaterCavitation` retains the original IF97 saturation pair and industrial
transport; its surface tension is caller-supplied for the actual interface.
It uses the saturated-liquid incompressible approximation and supplies neither
metastable EOS continuation nor a critical-region validity claim.

`fluid_acoustics::FrozenHomogeneousMixture` supplies Wood's compressibility through
the existing `FluidAcoustics` trait, preserving exact pure-phase endpoints and
converting mass quality with additive phase volumes. It assumes homogeneous
motion, frozen phase masses, phase-adiabatic response, common pressure, negligible
capillary pressure and negligible bubble inertia. It does not model resonance,
attenuation, heat or mass exchange, slip, equilibrium boiling sound speed, bubble
inception or collapse. G0/G3 regressions cover source retention/refresh, quantity
refusal, original water coefficients, pure limits, phase exchange, mass/volume
conversion and independent compression of the separate phase volumes.
Mass-basis admission retains its authored quality and both constituent volume
fractions independently. Density and compliance are computed directly from those
volumes, avoiding a lossy mass-to-volume-to-mass round trip at high density
contrast. `liquid_volume_fraction()` preserves a small liquid volume even when
`void_fraction()` rounds to one; neither rounded fraction alone identifies a
pure phase for a mass-basis state. Unrepresentable positive phase volumes refuse.

## Real-gas pressure forces and mechanical tangents

Closed ideal and Peng--Robinson parcels share the typed, borrowed
`GasMechanicalLoad` for work-conjugate generalized forces, mechanical power and
matrix-free tangent actions. The existing real-gas compression/heat owner supplies
the accepted material state. Its adiabatic bulk modulus is rho*c^2 from that
state's acoustic law, not the ideal-gas substitution (Cp/Cv)*p.

For caller-supplied V(q), gradient and Hessian actions, F=(p-pe)*grad(V) and
DF[dq]=(p-pe)*Hess(V)[dq]-(rho*c^2/V)*grad(V)*(grad(V).dq).
Exterior pressure is explicitly constant during the local linearization.
The default small-chart entry retains its 256-coordinate limit. Both closed
ideal and real-gas parcels also expose `mechanical_load_with_limit`, borrowing
a spatial chart under the geometry owner's explicit coordinate budget.
Cancellable force, power and tangent actions reuse the same pressure and
acoustic response, poll in blocks of at most 256 coordinates, and publish no
partial result on refusal. The geometry limit is not a constitutive restriction;
the pullback introduces neither mesh ownership nor another volume law.

Borrowing the parcel preserves the state epoch across force and tangent use.
The adapter owns no geometry, piston integrator or phase-equilibrium closure;
the existing homogeneous-branch restrictions continue to apply.

Focused regressions compare dense-gas stiffness with adiabatic pressure
differences, the full geometric tangent with displaced-state force differences,
and refined generalized-force quadrature with caloric compression work plus
the declared exterior-pressure work.

## Original-law real-gas property-temperature queries

`ThermochemicalPengRobinsonModel::evaluate_temperature_from_snapshot` retains
the accepted bulk pressure, active composition and declared basis, source
coordinates, transport claims and exact original PR law while querying a new
property temperature. Snapshots share immutable critical data and interactions;
normal repeated queries avoid rebuilding the PR vectors and pair matrix. A
different source bank must reproduce the actual active critical data and both
interaction-matrix axes. Inactive species remain inactive. Source admission
also runs when the original mechanical law is shared.

The caller explicitly selects Vapor or Liquid, which must reproduce the
accepted bulk root. Stable-root selection refuses. The returned homogeneous
property state supplies real density, Cp, acoustic and expansion derivatives;
it has no new inventory or storage. This does not establish equilibrium or a
continuous single-phase path through the intervening temperatures. Focused
regressions cover active-subset mapping, source substitution, branch mismatch
and original property-domain refusal.

## Conservative real-gas mixture transfer

`IsochoricRealGasParcel::transfer_mass_to` transports each upstream species and
its enthalpy between fixed volumes. Compatible mixtures may differ in fractions
or active species. The original ideal source join retains NASA, transport,
water and pack bindings; original PR critical properties and directed pair
coefficients must agree wherever shared. Every pair in the resulting catalog
must be present in an original model. Separate zero-interaction pure-species
catalogs do not supply a missing cross-species coefficient.

The receiver rebuilds its native mixture from actual extensive inventories and
reconstructs density, temperature, pressure, calorics, acoustics and transport.
At the new mixture's final density and old receiver temperature, conserved ideal
species reference energies cancel analytically. Donor sensible-energy change,
upstream p/rho and all original PR departure energies remain. The existing
density/energy increment solver owns inversion; no second EOS or root is added.
The identical-model path preserves its original numerical arithmetic.

The canonical `GasMassTransferStep`, also named `RealGasMassTransferStep`, reports
each species flux and its maximum actual inventory residual. Total mass, donor
withdrawal, receiver gain, pair species conservation and energy each close to
their declared tolerances. Existing nozzle and piston consumers inherit this
composition support. Source/domain refusal, exhausted work and cancellation
publish neither endpoint. This remains a homogeneous first-order transfer,
without equilibrium flash, jet momentum or an inferred mixing-entropy account.

Focused regressions compare an independent 80-digit dense-mixture reference,
continue native pressure work after mixing, shift separate species formation
datums by order 1e24 J/kg, introduce a species with known pair data, and refuse
missing/conflicting sources and unsupported or unfinished endpoint solves.

## Conservative real-gas heating at fixed mass and volume

`thermochemical_gas::nonideal::energy` reconstructs a complete source-bound
Peng--Robinson state from conserved density and internal energy. Absolute
energy, a fixed-reference energy change, and a temperature-anchored increment
are explicit entry points. The last integrates the original NASA ideal-reference
energy change and adds the PR departure change at the same density, avoiding
subtraction of large common formation energies. Total Cv supplies the tangent
to the existing bounded safeguarded root solver. Nonlinear iterations use
direct density/temperature EOS evaluation without cubic root solves or transport
card resolution. Initial and accepted material states retain source pressure
and temperature limits; the final PT inventory remains available separately
from the retained density root.

`nonideal::energy::parcel::IsochoricRealGasParcel` is a continuing rigid-vessel
consumer. Initial real-gas density and explicit volume determine conserved mass;
signed heat advances `Q=m Delta u`. Density and volume stay fixed while pressure,
Cp/Cv, thermal expansion, sound speed, impedance and source transport refresh
together. The original ideal-reference source remains provenance, not a second
mechanical density. A whole update commits only after its joule balance and final
cancellation check. Zero heat preserves the accepted state. Source refusal,
instability, work exhaustion and cancellation cannot publish partial state.

The caller declares a continuous, homogeneous, frozen-composition branch.
Positive pressure, negative isothermal pressure/volume derivative and positive
total Cv are checked at evaluated states. Mechanically stable metastable states
are allowed; this is not an equilibrium flash, phase-change energy model,
reaction model, vessel-wall model or spatial compressible-flow solver.

G1/G3/G4 regressions cover an independent high-precision cubic/departure heating
reference, pressure differing from the ideal fixed-volume shortcut, shared
material views, heat subdivision and cooling, large formation-energy offsets,
source-domain refusal, and late cancellation rollback. The `real_gas_vessel`
example reuses existing synthetic source-card input and emits pressure, density,
temperature, conserved mass, heat capacity, sound speed and joule residuals as
CSV. Its manufactured mixture supplies no empirical or vessel-safety claim.

### Conserved mass with real-gas volume work

`nonideal::energy::parcel::ClosedRealGasParcel` retains the same inventory and
source model while physical volume changes. Density follows conserved mass
divided by the accepted volume. Its reversible adiabatic operation reconstructs
temperature from the original NASA ideal-reference entropy plus the PR residual
entropy, with total `Cv/T` as the safeguarded root tangent. Both initial and final
densities enter caloric and entropy increments; the reference departure is never
silently recomputed at the candidate density. Species energy/entropy references
and fixed-composition mixing constants cancel before subtraction.

Accepted work by the gas is `-m Delta u` on the admitted isentrope. The separate
volume/internal-energy endpoint supports a consumer-owned heat/work balance,
including a zero-energy-change expansion whose real-gas temperature changes.
The endpoint API supplies no heat/work split, piston dynamics or time integration.
Fixed-volume heat reuses `IsochoricRealGasParcel`; conversions retain its exact
accepted state and inventory. Complete volume updates publish only after source
admission, entropy/work consistency and the final cancellation check. The same
homogeneous-branch limitations and original transport-card bounds apply.

`ClosedRealGasParcel::stored_energy_change_to` compares two already admitted
endpoints using the original NASA caloric integral and PR departure increment
at both endpoint densities. Exact conserved mass, the complete source model,
auxiliary query coordinates and quantity conventions must agree. Common
formation references cancel before evaluation. A real-fluid expansion can cool
at unchanged internal energy; this API reports the state difference without
assigning a heat/work split or imposing an adiabatic path. Cancellation and
refusal leave both endpoints unchanged. `at_volume_isentropic` exposes the
existing immutable proposal directly, avoiding a preliminary parcel clone in
mechanical consumers. Focused tests check independent PR compression work,
isoenergetic cooling, small work under large formation datums and source identity.

Independent high-precision binary-mixture regressions cover isentropic
compression pressure, temperature and work, constant-internal-energy expansion
cooling and entropy increase, partitioned work and reversal, shared material
views, post-compression heating, source refusal and final-audit cancellation.
The optional fourth `real_gas_vessel` argument is a final/initial volume ratio;
after heating it emits the volume path, cumulative gas work and entropy residual.


## Native axial and plane-stress plasticity

`plastic::reduction::UniaxialJ2` is the exact free-sided axial stress restriction
of the original rate-independent J2 law, with separate isotropic/Prager hardening,
signed axial plastic strain, accumulated equivalent strain and retained lateral
plastic contraction. Its checked work/heat entry delegates to the original 3-D
energetic owner. No additional heat fraction or cyclic calibration is inferred.

`PlaneStressReduction<L>` enforces zero zz/yz/zx stress using the original law
and frozen accepted history throughout local Newton and line search. It returns
full tensor strain/history and the Schur-condensed consistent tangent in tensor-
shear order [xx, yy, xy]. Existing Perzyna laws advance one physical interval;
local iterations do not accumulate material time. `PlaneStressWorkspace` reuses
the original dense factorization storage. Failed admission or iteration changes
no accepted state. These are small-strain material reductions; nonlinear shell
assembly, finite-strain plasticity, damage and experimental calibration are
separate capabilities. Focused native material tests cover axial/cyclic/energy
limits, nonproportional plane stress, retained transverse history and rate-law
semantics; no full-workspace or runtime qualification follows from those tests.

### Finite plastic sheet forming and thickness recovery

`plastic::finite::plane_stress` reduces the original multiplicative finite J2
return to traction-free sheet kinematics. Given the in-plane deformation and a
positive transverse seed, bounded Newton/backtracking solves all three Piola
tractions `P[:,2]=0`. The shared positive transverse factorization supplies the
Schur-complement tangent of that same discrete return. Every trial borrows the
same accepted plastic metric and accumulation; trial iterations cannot harden
the specimen or consume additional physical time.

`FinitePlaneStressResponse` retains the full three-dimensional deformation,
candidate metric, plastic/creep accumulations, native stresses and separate
incremental-potential, stored-energy and physical-dissipation accounts. Normal
thickness is recovered as `h_reference * F22`. Plastic volume is preserved by
the native metric update, while elastic and thermal volume changes remain
possible; an isochoric initialization is never a final total-volume constraint.
The identical reduction supports the original rate-independent, Perzyna,
Norton, combined creep/plastic and frozen-temperature families. Kinetic thermal
inspection without an actual physical interval still refuses. Source adapters
continue to own the original law, model card and validity domain.

The fixed-temperature `FinitePlaneStressMaterial` interface retains the complete
native constitutive behavior, not merely its elastic/plastic backbone. Original
J2/source bindings, Perzyna, Norton and simultaneous Norton/J2 laws share the
same spatial consumer while keeping their own history validation, accepted-state
conditions, independent storage inspection and physical interval. Positive
viscous overstress is admissible when the original law admits it; inspection
never performs another flow step. The kinetic reduced potential is not physical
dissipation. Rate-independent laws report no physical clock, and thermal
evolution still requires the native two-temperature work account.

The reduction uses fixed-size arrays and the original native spectral workspace;
spatial consumers poll cancellation between material points. Rejection
or exhausted local work returns no candidate and changes no accepted history.
G1/G3 tests in `finite_plastic_plane_stress` compare equibiaxial compression and
tension against independent logarithmic-strain equations, check residual shape
and thickness after unloading, differentiate frozen-history stress/potential,
and exercise transverse-shear memory and thermal/rate ownership. This is a
locally stable transverse branch of an isotropic yield law. It does not model
shell bending, wrinkling, necking, anisotropic forming limits, tool contact or
fracture, and synthetic constitutive tests do not establish measured metal
formability.

## Finite-rate liquid-water evaporation and condensation

`evaporation::WaterEvaporation` binds a strictly liquid interval of the original
IF97 `IsobaricWater` chart to `fs_thermal::evaporation`'s signed Stefan-film
transport. The same source owns liquid enthalpy, density, temperature inversion
and entropy; IF97 Region 4 owns saturation pressure and its derivative, and
Region 2 owns saturated-vapor enthalpy. The carried vapor-minus-liquid enthalpy
gap includes the original liquid pressure contribution. Latent energy is a
subaccount of carried vapor energy, never an additional thermal source.

The water adapter admits 273.16–373.15 K and liquid pressure at most 0.2 MPa,
restricted further by its caller-selected pure-liquid branch below boiling.
The shared original heat timeline evolves actual liquid mass and enthalpy;
reference geometry stays fixed and recession is reported separately. Explicit
gas-film temperature, volatile partial pressure, total pressure and binary
diffusion conductance are required. This does not supply activity/fugacity,
Poynting or humidity-enhancement corrections, freezing, nucleation, boiling,
droplet transport or native free-surface mass/topology evolution. Tests cover
source reuse, independent Stefan fluxes, signed heat/mass balance, analytic
temperature tangents and unsupported-regime refusal.

`WaterEvaporation::pure_water_vapor_surface_tension` exposes the original
IAPWS R1-76(2014) coefficient and analytic temperature tangent at an admitted
interface temperature. This explicitly pure-water-vapor reference query
neither qualifies an arbitrary carrier-gas pairing nor changes the original
evaporation rate, saturation pressure, calorics or heat account. No Kelvin or
Poynting correction is introduced. Literal source values and the original
adapter's narrower temperature interval are checked through this entry point.

The shared `fs_thermal::evaporation::KineticEvaporation` adapter may borrow this
same water source and select finite Hertz--Knudsen or small-drift Schrage
interface kinetics. The explicit accommodation probability changes only the
interfacial resistance, solved in series with the nonlinear Stefan film. All
water properties, source limits, caloric references and thermal hooks continue
to delegate to the original source. Finite-gas regressions exercise both phase
directions with original IF97/NASA-9 calorics, unchanged nitrogen inventory,
joint water mass/energy closure, and upstream vapor sensible heat. Their
authored accommodation is illustrative, not fitted water data. The local
equal-interface-temperature and continuum-film limitations remain explicit.

### Finite atmosphere and shared thermal history

`evaporation::gas::advance_surface_with_gas` couples the signed Stefan-film
transfer to the original finite gas species inventory and caloric inverse.
The explicit volatile binding retains one fixed liquid-to-gas enthalpy datum;
source disagreement is bounded at both surface endpoints and its admitted
remainder is carried in the energy transfer. Evaporating vapor carries its
surface enthalpy, while condensing vapor carries its upstream gas enthalpy,
including sensible heat once. Carrier inventories, source identities, gas
volume and auxiliary coordinates survive the original-owner reconstruction.
Finite humidity, pressure and temperature feed the following interval. Donor
inventory, temperature departure and equilibrium-overshoot limits refuse an
oversized step without publishing either phase. This frozen-film update is
first order and requires physical timestep refinement.

`evaporation::gas::equilibrium::advance_equilibrium_surface` explicitly binds
that same finite-gas film owner to an original sharp solid/liquid enthalpy chart.
Pure-liquid constituent enthalpy supplies vaporization compatibility and density;
bulk condensed enthalpy supplies sensible and fusion storage. Their difference
is retained in the original extensive mass/energy update, so evaporative cooling
can freeze liquid at melting temperature and condensation can melt old solid.
Upstream vapor sensible enthalpy remains included exactly once. Only initial
liquid inventory can fund evaporation; solid inventory cannot supply a liquid
interface. Final solid/liquid masses and signed fusion release are reconstructed
from the original chart, and fusion release is a diagnostic, never a second heat
source. Complete melting continues onto the original liquid sensible branch.

Exposed liquid area is explicitly authored; phase mass fraction does not infer
surface coverage. A positive liquid seed must remain at both endpoints. Complete
freezing requires an interface event owner before sublimation or nucleation can
take over. The original frozen-rate step bounds and fixed gas/reference geometry
remain in force. This reduction does not infer mushy alloys, pressure-dependent
phase diagrams, separate constituent temperatures or moving-interface geometry.

`evaporation::gas::equilibrium::boundary::advance_to_freezing` consumes at most
the requested phase-only interval and returns the conservative all-solid fusion
endpoint when reached, together with accepted and unconsumed physical time.
The event uses the original bulk/constituent energy difference and frozen film
rate, including cooling by sufficiently cold condensing vapor. It locates an
early endpoint before full-horizon exposure or donor admission, then retains
the original accepted-duration inventory, temperature and equilibrium gates.
The physical-rate mass request and represented terminal mass remain distinct
and share the original mass allowance. Exact source-knot selection must pass
the sum of absolute condensed and original-gas energy defects under one joule
budget. There is one gas reconstruction and no residual liquid or latent energy
is discarded. The ordinary scalar API continues to require a liquid seed.

The terminal state does not select a new exposed interface or switch from
continuum film transport to a molecular sublimation law. That handoff remains
explicit. If network or reciprocal gas heat is split over the same interval,
the complete heat/phase proposal must be retried at the shortened event time;
a full-duration heat state cannot be relabeled with the earlier phase clock.
Frozen-rate event times still require physical timestep refinement.

`evaporation::gas::equilibrium::heat::advance_equilibrium_surface_with_heat`
adds reciprocal gas-film heat through the existing implicit heat owner, evaluated
on that same bulk enthalpy chart. Fusion plateaus require neither apparent heat
capacity nor a unique condensed-body inverse h(T). Explicit controls name a bulk-enthalpy
interval; accepted original, heated and final endpoints retain a liquid seed.
The exact post-heat bulk state supplies the following vapor substep, including
its changed liquid inventory. Heat acts first and vapor transfer follows over
the same physical interval: this is a first-order Lie split that requires
complete-step refinement, not a simultaneous heat/phase solve.

The whole-interval phase partition includes heat-induced and vapor-induced
fusion once, with separate heat and vapor subreceipts. Signed storage defects
are reported, and the sum of both absolute substep energy defects must fit the
single original joule allowance. Original-to-final temperature limits and late
cancellation cover the combined proposal. Zero heat coefficient delegates
exactly to the phase-only path before inspecting unused heat controls. The
sensible coefficient must exclude vapor-carried sensible and latent energy.
The original fixed-volume, authored-area and interface-event limitations apply.

`evaporation::gas::equilibrium::heat::boundary::advance_to_freezing_with_heat`
locates conservative complete freezing with reciprocal gas heat and vapor on
one duration. Each bounded trial replays both original owners from unchanged
inputs. An early vapor endpoint only brackets time; it is accepted only when
its elapsed time equals the complete heat interval. The returned whole-body
phase account includes both stages, and unconsumed time belongs to an explicit
next interface owner. Zero sensible coefficient delegates exactly to the
analytic phase boundary without inspecting unused heat or event controls.

The heat bracket retains the exact original solid knot as its terminal lower
bound. Original nonzero-heat bisection endpoints may consequently differ within
their original temperature and energy allowances from the strict partial-phase
facade. With zero original vapor flux, including exact saturation, heat itself
may finish freezing. Any exact-knot correction must fit both the original
heat-equation tolerance and the single sum-of-absolute energy-defect allowance;
it is reported explicitly. Active vapor is never suppressed because heat has
exhausted the liquid seed: those trials only bracket the event or refuse with
a smaller-step request. Source, donor, endpoint and cancellation gates remain
shared. The terminal event belongs to the heat-first split; physical timestep
refinement remains necessary, and neither continuous-time uniqueness nor an
automatic solid-interface handoff is implied.

`advance_surface_with_gas_mass_selection` accepts a native inventory selector
before the existing caloric solve. The selected signed transfer must preserve
direction and nonzero mass within the original absolute mass budget. Both
requested and selected donor withdrawals obey the inventory limit, and all
subsequent temperature/equilibrium gates remain unchanged. The report retains
the original film request and the selected debit. The original entry point
delegates with the identity selector and retains its numerical behavior.

`advance_surface_with_gas_and_heat_mass_selection` preserves that native mass
selection when reciprocal sensible film heat precedes phase transfer. The
original implicit heat solver and source calorics remain unchanged. Initial
zero-duration phase admission uses the identity selector; the caller's selector
runs exactly once on the actual post-heating phase request. With zero sensible
coefficient it delegates directly to the selected phase path. The original
heat API delegates with identity selection and retains its existing behavior.

`evaporation::gas::timeline::advance_deposition_with_gas` advances one admitted
pure-liquid cell, its original heat/contact/source network and one finite gas
on the same event-clipped clock. Births and every heat substep retain original
caloric/storage admission. Catalog and control admission also run for no-op
and before-birth requests. Initial exposure survives endpoint covering; no
phase transfer precedes birth. Selected birth-datum offsets enter incoming
and stored energy once. Source-relative defects accumulate before diagnostic
formation references, with unchanged whole-window mass/energy tolerances and
one original thermal work allowance. Explicit inventory, temperature and
saturation step limits discard the complete proposed interval and retry from
the same accepted state at a strictly shorter representable endpoint. Rejected
thermal work spends that allowance; rejected events, births and transfers do
not enter physical accounts. The original anchored time grid resumes when its
endpoint is reached. Source/domain errors still refuse directly.

`evaporation::gas::timeline::equilibrium` admits the original sharp fusion chart
on that same event/retry loop. Native conduction, contact and source heat,
optional reciprocal gas-film heat, and signed vapor exchange evolve one bulk
enthalpy coordinate. The registered thermal caloric identity must match the
original phase chart, including scheduled births and zero-duration requests.
Numeric phase admission at heat substeps avoids rebuilding hashed phase reports;
the complete original-or-birth-to-final partition is reconstructed once. It
includes heat-driven fusion after the liquid interface becomes covered, which
a sum of vapor-only reports would miss. The original fixed geometry, birth and
exposure timing, attempted-work limit and whole-window conservation accounts
remain shared with the pure-liquid lane. A positive liquid seed is required at
each accepted heat and vapor endpoint; terminal freezing needs an explicit
interface event. The full heat/phase split requires timestep refinement, and
neither failure nor cancellation publishes either original owner's candidate.

`evaporation::gas::timeline::equilibrium::boundary::advance_deposition_to_freezing`
adds conservative completion to that same network/phase loop. Every shortened
phase event is localized by replaying the complete native heat and vapor proposal
from its unchanged accepted start. A clipped vapor endpoint is only a search
bound; acceptance requires its duration to equal the complete thermal interval.
All neighboring cells, native events and gas therefore share the accepted clock.
The observer retains the first all-solid native substep even if later heat would
remelt it, and rejected probes spend the original thermal substep allowance.
Event localization also has a fixed 64-trial bound; constitutive step limits
still request finer physical intervals through the original retry owner.

Covered or zero-transfer surfaces, including an exactly saturated film, can
complete freezing by native heat alone.
An exact source-knot selection is admitted only after a detected crossing and
within the original thermal joule allowance, including preceding thermal defects.
Its represented correction enters stored energy and the thermal residual, never
external heat. Active vapor retains its original donor, source and gas gates.
The whole-window phase partition includes the terminal zero liquid inventory;
the receipt returns unconsumed time for an explicit subsequent interface owner.
Initial/scheduled liquid-seed admission and source identity remain strict.
This endpoint belongs to the first-order split and requires timestep refinement;
it does not infer coverage, nucleation, moving geometry or a sublimation law.

`advance_deposition_to_freezing_with_heat` in that same boundary module adds
reciprocal gas heat to the complete event clock. Each native thermal trial
uses one shared gas-heat/phase trial; scalar and timeline facades consume that
same physical owner and candidate assembly. A clipped gas event exposes only
an event bracket, never a publishable shorter-clock state. The timeline owns
all event retries, avoiding nested searches and charging discarded native work
to the original allowance. Native-only terminal selection is permitted only
when any following gas sensible heat is inactive at the fusion endpoint.
Whole-window sensible, vapor, fusion and neighbor-cell accounts retain their
original owners. Zero gas heat coefficient delegates exactly to the native
boundary API before inspecting unused heat controls. The same split-domain,
fixed-geometry, source and explicit interface-handoff limitations apply.

### Several retained films sharing one finite atmosphere

`evaporation::gas::chamber::equilibrium` extends the original chamber reduction
to several sharp solid/liquid inventories with explicitly exposed liquid area.
Each original chart supplies pure-liquid constituent compatibility, density and
the actual liquid donor inventory; its bulk enthalpy retains sensible and fusion
storage. Only caloric preparation differs from the existing pure-liquid lane.
One canonical work vector, compensated species sums, original old-vapor withdrawal
gates and one original gas reconstruction remain shared. Per-body phase reports
retain source identity and distinguish freezing from signed vapor transfer.
Their constituent mass defects spend the same chamber-wide mass allowance;
the original sum of absolute caloric defects shares one energy allowance.
No concurrent evaporation finances excessive condensation elsewhere. Accepted
bodies retain a liquid seed; terminal freezing and reciprocal gas heat require
their explicit owners. This frozen-rate step requires timestep refinement and
does not infer surface coverage, nucleation or moving geometry.

`evaporation::gas::chamber` advances multiple positive liquid inventories against
one original gas snapshot. It reuses the single-surface Stefan-film preparation,
liquid caloric update and endpoint checks. Stable surface IDs canonicalize input
and native mass-selection order. Gross condensation of each species is bounded
by its old vapor inventory; simultaneous evaporation cannot finance an excessive
withdrawal. Canonical compensated species sums retain sensible energy when
opposing transfers have zero net mass. One original gas reconstruction supplies
the common final pressure, temperature and composition for every film's
equilibrium check. Native mass selectors run before caloric updates and preserve
both requested and represented donor limits.

The complete interval has one mass and energy allowance. Actual local liquid,
whole-gas and per-species inventory changes are checked independently. The sum
of absolute local caloric defects and the sole gas defect bounds acceptance,
so increasing surface count does not renew tolerance or hide opposing errors.
An empty active roster admits and preserves the original gas. This remains a
first-order phase-only frozen-film step with retained positive liquid mass and
fixed gas volume. Dry-wall nucleation, film extinction and volume-displacement
work remain outside this owner. Regression sources exercise shared
inventory, opposite-signed exchange, input-order invariance and source-relative
energy accounting; native execution remains a separate required validation.

The opt-in `evaporation::gas::chamber::heat` path solves simultaneous reciprocal
sensible heat before the same joint phase interval. Each film retains its own
explicit coefficient and original increasing liquid enthalpy branch. The common
gas temperature is a scalar energy root; independent bounded liquid responses
use their original caloric inverses. Safeguarded secants, retained brackets and
compensated heat sums require no frozen heat capacities, numerical derivatives,
dense matrix or per-trial gas-source reconstruction. Local response uncertainty
is bounded in joules before choosing an outer-root sign, including stiff films.
The original gas caloric branch is prepared once, and the actual represented
liquid increments determine the one heat-stage gas debit and reconstruction.
Every heat equation is rechecked at that reconstructed endpoint.

Films initially at the gas temperature remain active when other films move it;
intermediate films may reverse their original heat-flow direction. Original
source domains and global initial temperature extrema remain binding. All-zero
coefficients delegate to the existing phase path, and a single active film
reuses the original implicit heat solver. Exact heated liquid states are retained
for native mass-carry accounting. Combined storage and absolute energy defects
include sensible heat once; vapor-carried sensible and latent energy remain in
the phase stream. This heat/phase split is first order and requires full-step
refinement. Separate per-film and outer iteration bounds cap work, with
cancellation and transactional refusal throughout.

`evaporation::gas::timeline::multiple::advance_deposition_with_gas` connects
that same chamber to several original thermal cells. Heat/contact/source
transport advances once per event-clipped interval, then all exposed films
exchange with one common gas. Each cell retains its source chart, temperature,
liquid mass and fixed datum. Births receive no earlier phase exchange, and
endpoint covering cannot rewrite the preceding exposure. Step-limit refusals
retry the complete heat/phase proposal while spending the original attempted
thermal-work allowance. Whole-window total/species mass and absolute caloric
defects are checked without renewing the budget per film or interval. Per-cell
accounts retain both gross directions and recession; geometry remains fixed.
Regression sources cover cold/warm counterflow, original heat and birth energy,
covering events, finite-gas retries, cancellation, singleton compatibility and
refinement of the complete shared clock.

`timeline::multiple::advance_deposition_with_gas_and_heat` inserts that same
joint sensible-heat owner between original wall transport and joint phase
exchange. Every enabled film retains exactly its original thermal enthalpy
branch and explicit coefficient, including admission before birth. All three
stages share each event-clipped interval and discard together on a retry; only
accepted sensible heat enters per-cell and whole-window physical accounts.
The existing window-wide absolute energy budget includes wall, joint heat and
phase defects once. All-zero coefficients preserve the phase-only endpoints
and audit exactly. Five regression sources cover an independent nonlinear
thermal-network solution, opposing phase transfers, birth/covering timing,
retry/cancellation/zero-heat behavior, and complete-clock time refinement.

### Simultaneous thermal-network and finite-atmosphere phase endpoints

`evaporation::gas::implicit::solve_deposition_with_gas` solves every original
thermal cell, each selected volatile inventory and the original finite-gas energy
at the same backward-Euler endpoint. Selected positive liquid cells retain
their individual phase laws, caloric offsets, pure-liquid domains and exposed
reference areas. Distinct pure-liquid interfaces may name distinct original
volatile species in the same finite atmosphere. At each interface the other
gas species are nonabsorbing carriers under the supplied Stefan-film closure;
this does not add mixed-liquid activities or multicomponent cross-diffusion.
Unselected carrier inventories, source identities, gas volume and query context
remain preserved.
Endpoint total/vapor pressure and gas temperature feed the original Stefan
films; native heat sources, contacts, baths and conservative heat transfers
remain inside the same thermal solve. Each supplied sensible film coefficient
exchanges `dt*h*A*(Tgas-Tliquid)` with the opposite gas debit exactly once.

Evaporating mass carries the original vapor enthalpy at its liquid surface;
condensing mass carries upstream enthalpy at the endpoint gas temperature.
The shared corrected thermal ingress includes source alignment and incoming
vapor sensible enthalpy in complete stream energy while retaining latent heat
as `dm*L`. Opposing surfaces retain their individual receipts, gross withdrawal
limits and energy transfers even when their net vapor mass is zero. Condensation
spends each species' own original gas inventory; simultaneous evaporation and
another species' inventory cannot finance an oversized withdrawal. Every
original liquid and source-domain gate remains.

The original `fs_solver` finite-difference Newton--Krylov owner solves native
cell-energy coordinates plus one numerically scaled mass per selected volatile
and one gas temperature. Species coordinates use canonical identity order.
Initial gas Cv supplies only a coordinate scale. Each gas residual evaluates
the original final-composition NASA-9 internal-energy increment directly:
`dU_sensible - Tgas0*sum(dm_species*R_species)
- sum(dm_surface*dh_v_species(Tgas0,Tupstream_surface)) + Q_sensible`.
Opposing transfers of different species retain their distinct gas constants,
calorics and latent heats even when total gas mass is unchanged. Each initially
absent vapor receives its own source-owned mass predictor; no trace inventory
is manufactured. A species with zero net mass can still carry nonzero energy.
Formation constants cancel before arithmetic; no frozen heat capacity or
per-trial nested caloric root is used. One final original species-boundary
reconstruction supplies the accepted gas, followed by a fresh thermal/film
evaluation at that actual gas state. Native cell/global energy and mass gates,
carrier/volatile/total gas mass, source-relative gas energy, and the combined
mass/energy account must all pass before either owner is returned. The gas
reconstruction receipt and final coupled consistency residuals remain distinct.

The residual cache holds two numerical points and never accepted payloads.
The source-owned initial mass predictor, every uncached stencil/line-search
trial and both final physical checks spend one evaluation budget; the sole
gas inversion retains its original bounded work. Original source refusals and
cancellation propagate atomically. Invalid numerical probes refuse rather than
clipping humidity or manufacturing trace vapor. Callers must choose source-
supported difference scales and durations. This event-free, first-order owner
requires physical timestep refinement: start events must already be processed,
interior events refuse and endpoint events remain pending. It adds no boiling,
dry-film onset, moving interface, liquid displacement work or bulk melting model.

The opt-in `implicit::adaptive::advance_deposition_with_gas_adaptive` advances
an event-free window with physical step-doubling control. The original
first-order controller compares a full interval with two actual half intervals
and accepts only the conservative fine trajectory. It measures every native
cell's enthalpy and temperature, each selected liquid's mass, gas temperature,
and every gas species inventory. Relative tolerances scale changes from the
accepted interval start, preserving small latent/inventory changes beside large
reference enthalpies or stationary inventories. The estimate controls local
time error and requires independent validation near constitutive corners.

Only accepted fine intervals enter per-cell mass, both gross phase directions,
latent heat, complete native-datum vapor energy, vapor-enthalpy correction,
sensible heat and recession totals. One whole-window energy budget bounds
absolute thermal and gas defects; the thermal contribution also bounds the sum
of absolute selected-surface defects. Actual final total/species/carrier mass
and cumulative absolute mass defects are independently gated. Source-relative
accounts preserve caloric-reference cancellation without subtracting absolute
formation energies. Coarse and accuracy-rejected trajectories spend the same
Newton, residual-evaluation and gas-iteration budgets as accepted work. Final
payloads come only from the last accepted fine interval, and late refusal or
cancellation discards the entire candidate. Source-domain/numerical failures
propagate; only an admitted physical time-error estimate requests reduction.
Events inside the requested window refuse upfront; exact-target events remain
pending for the original timeline owner.

Both endpoint and adaptive reports retain net transferred mass by original
species identity. Whole-window closure independently gates the actual change
of each volatile inventory against that species' accepted transfers; the signed
aggregate vapor residual cannot conceal opposing species errors. Declared
future species retain zero transfer until activation, and each new active
species joins the same coupled endpoint without restarting request budgets.

An empty phase roster uses only the original thermal energy coordinates. Native
nonlinear heat and contact transport still advance, while the complete gas state
is preserved without a volatile identity, film or gas iteration. Original gas
catalog identity, heat-source freshness, numerical work and cancellation gates
remain. The endpoint still requires active heat geometry.

The opt-in `implicit::timeline::advance_deposition_with_gas_timeline` carries the
same adaptive coupled solve through native deposition births, contact/covering
changes, separations and boundary changes. Comparisons stop exactly at each
event. The original event owner applies same-time changes in authored order,
including start and target events, with zero permitted heat substeps. Empty
geometry advances directly to the next event or target; existing heat cells can
advance before any selected liquid exists. Future selected births are admitted
structurally, but their phase/caloric source is queried only after activation.
Each phase interval uses its own original exposure and positive donor inventory.
The canonical active roster is refreshed at event boundaries, not every trial.

One comparison loop retains numerical, thermal-substep and conservation budgets
across every event. Only accepted fine phase intervals enter transfer totals;
birth mass and original-datum birth energy are reported separately, together
with selected birth mass times its explicit gas/liquid enthalpy offset. Whole-
window mass closure compares each final cell with its own initial or actual
birth mass before reduction. Native event defects and accepted source-relative
thermal/gas defects enter the same absolute energy allowance without subtracting
large aggregate formation energies. A birth exactly at target is admitted but
receives no earlier phase exchange. The last payload describes the preceding
accepted fine endpoint; an all-empty/event-only advance has none. Source,
topology, domain, budget or cancellation refusal discards the complete candidate.
Callback-only source discontinuities require an authored event or explicit
advance boundary, retaining their left-limit value on the preceding interval's
endpoint. Native boundary events enforce this ordering directly. This adds
scheduled integration, not boiling, dry-out, moving
interfaces or a new phase constitutive law.

### Live molecular transport inside the coupled evaporation endpoint

`implicit::solve_deposition_with_gas_and_films` permits a source-owned resolver
to update representative film temperature and binary mass-transfer coefficient
from the actual thermal and finite-gas trial. The original state, every uncached
physical residual and both fresh final checks use that resolver, including the
check after the original gas caloric reconstruction. Pressure and vapor partial
pressure remain the exact original gas values, and the callback must preserve
the exact active cell roster. Neither a manufactured humidity nor a stale
atmosphere can substitute for the finite inventory in this ingress.

The matching `adaptive::advance_deposition_with_gas_adaptive_and_films` and
`timeline::advance_deposition_with_gas_timeline_and_films` resolve transport
independently on each coarse and fine trial. Future cells are queried only
after their native birth. Source history remains immutable during speculative
evaluation; source discontinuities retain the existing event/left-limit rules.
Only accepted fine states contribute transfer accounts. Original callbacks,
static-film entry points, finite-gas sources, nonlinear budgets and whole-window
conservation checks share the same implementation. A no-op resolver preserves
the original endpoint and work report exactly.

The resolver supplies transport, not an independent heat or phase update.
Sensible film heat and carried vapor enthalpy retain their existing owners.
The source-backed laminar plate/strip adapter in `fs-airflow` uses this ingress
to refresh original binary diffusivity and gas viscosity/density with evolving
temperature, pressure and composition, under that correlation's admitted
laminar dilute-vapor regime. The generic callback adds no turbulence, droplet,
nonideal gas, multicomponent cross-diffusion or strong-blowing model.

G0/G1 regressions compare signed endpoint rates against the returned liquid
and gas, preserve carrier inventory and the original energy accounts, exercise
adaptive and birth timelines, and reject altered gas pressures, changed rosters,
invalid coefficients and a late source refusal without publishing either state.

### Spherical liquid inventory in the original finite gas

`evaporation::gas::droplet::{advance_droplet_with_gas,
advance_droplet_with_gas_and_heat}` connects the shared spherical geometry to
the existing finite-gas evaporation and reciprocal sensible-heat owners. The
liquid stores only mass and original enthalpy. Each interval derives its radius
and exposed area from original liquid density, applies the stagnant-sphere
`k_m=D_AB/r` coefficient, then reconstructs accepted geometry from final mass
and final temperature-dependent density. Molecular binary diffusivity is
explicitly supplied; no default, heat/mass analogy or forced-sphere correlation
is substituted. The existing Stefan film and selected interface kinetics remain
responsible for phase transfer.

Actual gas pressure, vapor inventory, composition, gas caloric inversion,
carrier preservation and original liquid/vapor energy alignment retain their
existing owners and refusal gates. The optional heat-first split uses the same
initial area over its physical interval; latent/carried vapor energy and
reciprocal sensible heat are accounted once. The immutable candidate reports
both geometries, initial `D/r`, the phase substep, total exchange and optional
sensible heat, and returns the original gas catalog. Failed final geometry or
cancellation publishes neither liquid nor gas. Zero diffusivity seals mass
transfer while admitted sensible heat can still change density and radius.

Initial geometry and supplied diffusivity are frozen within a step: the whole
step is first order and must be refined, including changes in the transport
source. The original finite gas remains a fixed external control volume; liquid
volume is not subtracted and no sealed-vessel displacement/pressure-work model
is implied. Positive spherical seeds, uniform liquid temperature, quasi-steady
stagnant continuum radial diffusion and negligible forced/buoyant convection
are explicit scope limits. Complete dry-out, nucleation, moving/deforming
sprays, breakup, coalescence, Kelvin effects and rarefied transport remain
outside this owner.

Added regressions exercise independent spherical signed rates, sequential
geometry refresh, reciprocal sensible heat at zero diffusivity, thermal density
changes, original species/energy conservation, source identity, fixed gas
volume, inventory refusal and cancellation. A source-backed case uses the
existing IF97 water and NASA gas models for evaporation and condensation.
These regressions were added but could not be executed in this session because
the pinned build was killed in the unchanged upstream dependency `asupersync`.

### Condensate onset and complete dry-out at a prescribed wall

`evaporation::gas::wall_film` admits an explicitly fully wetting, locally planar
film in thermal equilibrium with a prescribed wall temperature. Its
`GasWallFilmInventory` distinguishes `Dry` from a positive `Wet` inventory:
a dry surface carries no liquid temperature, heat capacity or hidden seed mass.
The interface temperature comes from the supplied liquid enthalpy on its original
chart. A supersaturated gas can create a film from exactly zero mass; an
undersaturated dry surface cannot supply vapor. This is a declared wetting
closure, not a prediction of heterogeneous nucleation or a contact-angle model.

The original Stefan law, species identity, caloric alignment and finite-gas
reconstruction remain authoritative. Evaporation stops at the exact inventory
exhaustion time of the frozen rate, records that time, and returns `Dry` without
an artificial residual film. The complete requested interval remains the shared
clock; phase flux is zero after that event. Vapor-donor, source-domain,
temperature-departure and final equilibrium checks still refuse oversized
intervals. The existing positive-film owner retains its stricter withdrawal
limit and unchanged semantics.

The wall supplies both film sensible storage and the difference between carried
vapor enthalpy and retained liquid enthalpy. Incoming vapor retains its upstream
gas enthalpy; outgoing vapor uses the interface temperature. Original NASA
increments and the admitted phase-source remainder determine this heat before
large formation references can cancel it. The reported wall heat closes gas
internal energy plus film enthalpy, including film birth and complete exhaustion;
latent heat is a diagnostic part of that transfer, never a second source. A
finite wall must debit this same heat through its own thermal owner. Gas volume
and surface area remain fixed. The frozen-rate/time-equilibrated closure requires
physical timestep refinement and does not include film drainage, liquid motion,
boiling, ice, gas momentum or liquid displacement work.

### Simultaneous zero-capable wall-film endpoints

`evaporation::gas::wall_film::endpoint::GasWallFilmEndpoint` exposes the original
wall-film phase law and finite-gas residuals to a coupled thermal solver. Its
trials resolve the actual proposed gas inventory, pressure, temperature and
transport through the original source catalog without an inner caloric root.
The same wall-film preparation serves the existing frozen-initial-gas advance
and this endpoint evaluation; no second saturation, latent-heat or gas law is
introduced. The caller supplies the original liquid-chart wall enthalpy and an
explicit nonadvective gas heat port, separate from vapor-carried enthalpy.

At the endpoint, condensation uses the upstream FINAL gas enthalpy. Evaporation
uses the interface enthalpy and is limited by the actual initial liquid mass.
Dry surfaces cannot evaporate, and supersaturation can create a positive film
without seeding. The wall debit includes original film equilibration plus the
complete vapor/liquid enthalpy gap, including at complete dry-out. Condensation
is bounded against the original vapor donor; neither newly transferred mass nor
another species can enlarge that allowance. The gas balance retains species
pressure-flow energy, with formation references canceled before evaluating the
small thermal defect. Volume, carrier inventories, source identities and
auxiliary gas query coordinates remain fixed.

Mass and energy residuals have independent physical gates. `finalize` performs
one original species/energy reconstruction and one fresh endpoint-rate and
energy evaluation against the returned gas. Its absolute audit additionally
charges the original reconstruction defect. Failure or cancellation returns no
joint candidate. A finite wall must close the opposite wall heat port and
accept its own storage together with this gas/film candidate. The endpoint
owner supplies no nonlinear algorithm, wall heat capacity or external heat law.

The endpoint complementarity permits exactly zero liquid storage. Its dry-out
time is a diagnostic obtained by freezing the ENDPOINT rate over the interval,
not an exact event time for the continuous trajectory. Physical timestep
refinement remains necessary. Fully-wetting, fixed-area, fixed-gas-volume and
pure-liquid source limits remain unchanged; this does not add nucleation,
boiling, drainage, gas momentum or liquid displacement work. Focused regressions
manufacture independent nonlinear-caloric endpoints for dry condensation,
complete evaporation, separate gas/liquid datum shifts and atomic refusal.

### Reciprocal sensible heat through the gas film

`evaporation::gas::advance_surface_with_gas_and_heat` adds a caller-supplied
heat-transfer coefficient on the same exposed area. Implicit backward-Euler
heat exchange uses the original liquid enthalpy inverse and a prepared original
gas caloric branch. Scalar trials evaluate the original gas internal-energy
increment; a single original-owner reconstruction publishes the heat endpoint.
The same represented heat enters the liquid and leaves the gas. The gas keeps
its mass, composition, volume, sources and auxiliary coordinates, while its
temperature, pressure, transport and acoustic properties are reevaluated.

Implicit sensible heat precedes the existing phase transfer over the same
interval. `GasFilmHeatSolve` supplies the admitted liquid enthalpy range,
temperature residual tolerance and iteration bound. Original-state datum
compatibility, whole-step temperature limits and the same absolute combined
energy budget apply. A zero heat coefficient delegates to the original phase
path. Film heat excludes the sensible energy already advected by vapor.
The composed heat/phase method is first order and requires timestep refinement;
implicit heat stability does not certify split accuracy.

`timeline::advance_deposition_with_gas_and_heat` uses this same exchange on the
original event clock and bounded retry driver. Heat controls are admitted even
for no-op or before-birth requests, and the complete heat branch lies within
the original caller-admitted liquid range. Signed sensible heat is reported
separately from vapor-stream energy and included once in both storage accounts.
The `water_evaporation` example loads the existing N2/H2O species, NASA-9 and
transport packs and exposes evaporation/condensation modes, both film
coefficients and timestep refinement. Its CSV retains liquid mass/temperature,
gas humidity/pressure/temperature, transport and acoustic response together.
The supplied film coefficients and geometry are explicit example choices,
not experimental evaporation-rate evidence.

The conserved thermal account is liquid enthalpy plus gas internal energy in
fixed reference geometry. It is not the total internal energy of a sealed
moving two-phase vessel. The regression suite covers independently computed
storage balances, finite backpressure, heat/event coupling, datum shifts,
cancellation/refusals and timestep refinement. Two source tests retain the
checked-in NASA-9 water/nitrogen coefficients and original IF97 water chart;
their transport cards and film coefficient are explicit numerical choices,
not experimental evaporation-rate evidence. Nonlinear sensible-heat regressions
cover reciprocal storage, zero-heat parity, stiff approach to common
temperature, complete split refinement and transactional source/cancellation
refusal. Moving interfaces, dryout, boiling, gas momentum and spatial two-phase
flow remain outside this bounded owner.

### Anisotropic permeability with an explicit transverse constraint

`IntrinsicPermeability::confined_permeability` evaluates `1/(n.K^-1.n)` for
fully developed flow constrained to the prescribed normal. The existing
`normal_permeability` evaluates `n.K.n` for a prescribed normal pressure
gradient with unconstrained transverse flow. `DarcyFlowConstraint` makes this
physical distinction explicit; no tensor is replaced by an axis-aligned
average. An inverse square-root norm retains representable small positive
permeability without overflowing `1/k`. A normal with a nonzero component in a
sealed material direction gives exactly zero confined permeability.

`DarcyGasBarrier::with_flow_constraint` compiles the selected reduction once
through the original ordered layers. Gas flow, interface pressures, slip
calibration, repeated wet-mobility updates and continuum hydraulic consumers
all retain it, together with the original sources and support. This is a local
one-dimensional fully developed constitutive reduction. It does not solve
finite-plug entrance flow or transverse pressure matching across heterogeneous
interfaces. Those require a spatial tensor flow solve.

G0/G1/G3 tests in `porous::permeability::tests::confined` independently invert a
rotated tensor, verify zero transverse velocity, exact seals, extreme positive
coefficients, shared gas/continuum layer balances and source/wet-state retention.
`fs-airflow`'s `liquid_pressure::porous` test exercises the changed operating
point through the original material liquid pressure and heat owners.

### Heterogeneous dense membrane layers and interlayer permeance

`MembranePath::Series` composes up to 64 ordered dense layers and explicitly
sourced interlayer molar permeance on one common area. Each species has one
series flux: `G_i = 1 / (sum(L_j/P_ji) + sum(1/H_ji))` and
`J_i = G_i (p_i,left - p_i,right)`. Resistance ratios are scaled before
summation, avoiding direct reciprocals of tiny permeances. Effective laminate
permeance is reported separately from every original bulk permeability and
interface coefficient; no effective bulk material claim is invented.

Each element retains its original selected material or interface claims.
Ordered `InterfaceSystemCard` surfaces must match adjacent material identities.
The solver reconstructs Henry-equivalent species partial pressures, sums them
at each interface and checks each element's same source claims at its actual
local temperature and pressure. Selective counterpermeation may place interior
total pressure outside both reservoir pressures; reservoir-only admission
would miss that failure. These are equivalent gas potentials, not polymer
stress or additional gas inventories.

An explicit zero permeability or permeance seals the species exactly. Multiple
exact seals isolating an undetermined interior species potential refuse this
steady no-storage model. The single-film API delegates to the same mechanism,
and `MembraneNetworkLink::path` selects a film or series path. All rates enter
the original shared species and donor-partial-enthalpy network update, including
its gross withdrawal, conservation, final-state entropy and cancellation gates.
No second energy update or additional dissipative heat is introduced.

Focused G0/G1/G3 regressions cover analytic layer/interface resistance, retained
identity and local source-domain refusal, exact seals, reversal, small-permeance
arithmetic and shared mass/energy conservation. This is a near-isothermal,
linear Henry-sorption, negligible-storage closure. Porous supports, bypass
defects, nonlinear competitive sorption, polymer stress, thermal storage and
transient breakthrough remain outside this model.

## Source-bound dielectric polarization memory

`field_storage::dielectric::ResolvedDebyeCapacitor` retains one complete material
bundle, explicitly named high-frequency absolute permittivity, positive absolute
permittivity increments and positive relaxation times, plus uniform facing area
and gap. The shared source gate enforces SI dimensions, original claims and fixed
external conditions; frequency/field/strain/stress-dependent parameters cannot be
frozen into this time-domain law. At most 16 processes are admitted, with distinct
source properties and no fitted defaults. Empty process lists recover ordinary
lossless capacitance.

For the `exp(i omega t)` convention the law is
`epsilon = epsilon_inf + sum(delta_epsilon_k/(1+i omega tau_k))`. Its native
realization is C_inf in parallel with series R_k/C_k arms, where
`C_k=delta_epsilon_k A/d` and `R_k=tau_k/C_k`. Each process retains a separate
polarization charge. Storage is `q_inf^2/(2 C_inf)+sum(q_k^2/(2 C_k))`, rather
than applying a frequency-dependent capacitance to total charge. The bounded,
allocation-free admittance evaluator retains the passive high-frequency loss
limit even when omega*tau overflows.

This is an explicitly declared isotropic linear Debye model at frozen material
state and geometry. It does not infer measured pole data, parameter-temperature
work/entropy, leakage, electrode ESR/ESL, fringing, hysteresis or breakdown.
`fs-couple::material_circuit::dielectric` supplies its native circuit consumer.
G0/G1 source tests cover geometry and SI scaling, source retention, independent
single-pole frequency limits, multipole additivity, the zero-memory limit and
passivity/source-domain refusals. These synthetic fixtures verify the model's
implementation and do not establish experimental material accuracy.

### Bulk polarization memory on native field cells

`field_storage::dielectric::bulk::ResolvedDebyeMaterial` admits the same original
parameter roles and complete material bundle without assigning capacitor
geometry. Each fixed isotropic pole retains its own 1D, 2D or 3D polarization
vector and obeys `tau Pdot + P = delta_epsilon E`. Total displacement is
`D=epsilon_inf E+sum P`; physical storage density is
`epsilon_inf |E|^2/2+sum |P|^2/(2 delta_epsilon)`. Polarization loss is separate
from Ohmic conduction and must not be counted again as conductivity.

Allocation-free midpoint elimination supplies a positive effective permittivity
and vector history load for the native field solve. The matching bounded local
recovery returns every polarization endpoint, total displacement, storage change,
nonnegative polarization heat and midpoint electric work. Storage change uses
actual component increments instead of subtracting biased total energies. The
floating-point balance residual remains explicit for the field owner's physical
tolerance. This midpoint law is second order and passive; it does not promise
monotonic decay or an exact exponential response at unresolved relaxation times.

Coefficients, their source query, and material orientation stay fixed. This
opt-in constitutive owner adds no thermal refresh, anisotropic pole, ferroelectric
hysteresis, breakdown, or measured-fit claim. Existing zero-pole storage and
lumped uniform-gap reductions remain available. G0/G1/G3 tests in
`tests/dielectric/bulk.rs` compare source retention and circuit storage, the
independent pole ODE and work identity, second-order exponential-decay convergence,
discrete harmonic response, vector rotation, lossless limits, and invalid inputs.

## Source-backed optical tables on finite vacuum-wavelength bands

`state_point::SpectralConductorStatePoint` and
`SpectralTabulatedDielectricStatePoint` resolve nine original vacuum wavelengths
in SI metres together with their absolute complex-index `(eta, k)` or real-index
and intensity-absorption `(n, alpha)` values. Their explicit resolvers each
require 27 properties from the same immutable card and complete physical state.
Wavelength coordinates use `SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES`; the
separate conductor and dielectric property arrays select the constitutive
family explicitly. The existing visible schemas and their selection rules are
unchanged. There is no material-name inference, spectral resampling or RGB fit.

The conductor ingress delegates its complete n/k projection to the existing
`SpectralComplexIndexStatePoint` owner through `complex_index()`. It reuses that
owner's samples, nominal interpolation, support and identity instead of defining
another complex-index law. The SI wavelength properties and their original
receipts remain in the complete resolved bundle; the wrapper rebuilds the
projection's wavelength/key bindings on refresh, so temperature-dependent source
coordinates cannot leave stale wavelengths behind. Particle, EM and surface
consumers can therefore share the exact same projection. The separate dielectric
schema preserves original alpha interpolation; it does not interpolate a
converted k table and thereby change the source absorption law.

Every wavelength is positive, finite and strictly increasing. The first and
last coordinates define the closed source band. Index conventions are absolute
relative to vacuum. The dielectric retains the native `[1,3]` real-index,
weak-absorption model; its nonnegative `alpha` uses the exact static
power-neper-per-metre quantity in `I(d) = I(0) exp(-alpha d)`. Amplitude-neper,
decibel, dimension-only attenuation and external slab losses cannot substitute.
Conductor `k` uses the explicit static `OpticalExtinctionCoefficient` quantity,
with nonnegative dimensionless values distinct from attenuation per metre; the
opaque rendering consumer additionally admits its own strictly positive `k` domain.
These snapshots add no extrapolated tails, complex dielectric interface,
birefringence, thermal-emissivity law or measurement authority.

`try_from_material` can adopt a larger original snapshot while retaining its
other properties and receipts. Optical quantities must remain canonical scalars
at fixed external conditions; an additional frequency, field or deformation
coordinate cannot be silently frozen into a wavelength knot. `at_query_point`
pins every original claim and coordinate meaning, then readmits all wavelengths
and coefficients together. A missing source, invalid grid or passivity refusal
publishes no partial candidate. Schema-specific identities prevent different
families adopted from the same larger snapshot from substituting on refresh.
SI coordinates convert once to the shared complex-index owner's vacuum nm;
dielectric coordinates convert once at rendering ingress. The renderer separately
admits the caller's retained interpolation-gap limit. That limit is not an interpolation
error certificate, and endpoint refresh does not certify the intervening path.

Focused regression sources exercise a nonuniform 8--14 micrometre grid,
temperature-dependent coefficients and narrow absorption, source receipts,
competing-source refresh, grid/family identities, missing quantities,
intensity/amplitude conventions, rebuilt source-wavelength bindings, domain
refusals and transactional continuation.

## Source-backed visible dielectric dispersion

`state_point::VisibleSellmeierStatePoint` retains three nonnegative,
dimensionless oscillator strengths and three squared vacuum resonance
wavelengths in SI m², together with the same Beer-Lambert source schemas
used by Cauchy dielectrics: either four RGB transmittance/reference-distance
properties or nine typed spectral intensity-absorption knots. Every required
property resolves atomically from one immutable card at one complete state
and selected source policy. Spectral absorption retains its power-neper
quantity kind and wavelength grid; it is not converted through RGB. Partial
spectral declarations do not fall back to RGB, and ambiguous simultaneous
absorption schemas require the existing explicit property-selection policy.
The source index convention is absolute relative to vacuum; no ambient-air
correction or material-name default is inferred.

The generic visible optical resolver selects exactly one complete conductor,
Cauchy or Sellmeier schema and refuses ambiguous complete alternatives.
`fs-render::material_state` converts squared resonance wavelengths to µm²
once, admits the renderer's transparent band through its existing dielectric
owner, and checks that optics and mechanics share the original card, query
and coordinate meanings. Both generic scenes and the Euler scene adapter
consume this binding. Surface finish remains independent of bulk dispersion.

Focused G0/G1/G3 regressions cover source/temperature retention, SI optical
evaluation, source-only coefficient changes, both surface finishes, ambiguous
families, missing or wrongly dimensioned coefficients, gain, extrapolation,
cross-state substitutions and a resonance inside the render band. These are
synthetic source and algebra checks, not experimental glass validation.

## Native Voce hardening and exact axial reduction

`plastic::voce::VoceJ2` adds an independently selected isotropic saturation
`Q (1 - exp(-b alpha))` to the existing linear isotropic/Prager J2 law. The
bounded scalar return shares the original tensor flow and consistent tangent.
`J2State` retains all plastic strain and accumulated history; `Q = 0` delegates
exactly to the original law. The energetic interpretation stores the integral
of hardening, dissipates only the initial-yield offset, and keeps the convex
backward-Euler remainder out of physical heat.

`plastic::voce::source` requires independent saturation magnitude/rate and
plastic-strain hardening data, retaining the complete material bundle. A single
post-yield tangent cannot identify these parameters. Coefficients are fixed over
the load path; changing source temperature alone does not establish thermal
history or energy transfer. `plastic::reduction::UniaxialVoce` uses the same
scalar return with the exact free-sided stress constraint, lateral plastic
contraction, and native energy account. Borrowed `SmallStrainLaw` forwarding
lets reductions retain source ownership without copying material tables.

G0/G1/G3 tests in `plastic_voce`, `plastic_voce_source`, and `plastic_reduction`
check nonlinear reference curves, return consistency, tensor/condensed tangents,
cyclic energy, source refusals, and the exact zero-saturation limit. These are
numerical implementation checks, not empirical calibration. Finite strain,
nonlinear backstress recovery, damage, and temperature-dependent Voce energy
coupling remain outside this law.

### Quasi-steady concentration polarization

`osmosis::polarization::PolarizedMembraneTransport` retains the original
selective membrane and optional per-species liquid films on both faces.
`FilmPermeance` is a strictly positive finite mass-transfer coefficient [m/s];
the explicit diffusion-layer conversion is D/delta. An absent film is an exact
identity. Film mass transfer is distinct from hydraulic permeance and does not
add an invented mechanical pressure loss or a second membrane area.

At each proposed solvent flux, the films use the original exponential-fit
convection/diffusion law with solvent convection velocity Jv; the selective
membrane keeps its original (1-sigma_i)Jv velocity. Analytic series elimination
gives nonnegative bulk-to-bulk donor coefficients and positive reconstruction
of both membrane-face concentrations. The original water law consumes those
face concentrations. Leaking species retain their actual shared solute flux;
the perfect-rejection exponential formula is not applied to permeating species.

Prepared boxed arrays are bounded by the original 64-species limit. A cached
coefficient value reuses each fit for the finite inventory update and face
reconstruction without allocation or a nested water solve. Scaled arithmetic
and logarithmic recovery preserve representable extreme-Peclet limits; there
is no epsilon leak or silent enrichment cap. Unrepresentable physical faces
refuse. `heap_bytes` counts the owned arrays for consumer memory admission.

G0/G1/G3 checks use independent layer flux equations, harmonic diffusion
resistance, full-rejection enrichment/depletion, orientation reversal, absent
species, nonselective uniform transport, passive dissipation, and the limiting
water flux at large hydraulic permeance. These are prescribed, negligible-
storage film models, not resolved crossflow turbulence, evolving fouling,
finite film storage, electrochemical polarization or calibrated membrane data.



### Calibrated inertial resistance of porous layers

`ForchheimerCalibration` retains explicit per-layer beta [1/m], specimen normal,
transverse constraint, provenance and absolute T/p/superficial-speed support.
The convention is `-dp/dx = mu*q/k + rho*beta*|q|q`, without an implicit factor
of one half. Original intrinsic tensors and material sources remain unchanged;
compiled support intersects every calibration. Confined flow uses the original
harmonic tensor reduction. A prescribed normal gradient requires a permeability
eigen-direction for each inertial layer, so unresolved transverse velocity is
not silently discarded. No beta is inferred from k, porosity or an Ergun preset.

Isothermal no-slip ideal-gas transport solves
`mu*sum(L/k)*j + sum(beta*L)*|j|j = (p_left-p_right)*p_mean/(R*T)` analytically.
Scaled dominant-term branches avoid squaring dimensional coefficients and
subtracting nearly equal roots. The original profile owner reconstructs each
layer's own viscous and inertial pressure-squared loss. Both face volume rates
use their actual densities while one mass rate is conserved. The original drag
work diagnostic remains a diagnostic, never a second heat source. Exact zero
beta retains the existing Darcy/Klinkenberg arithmetic.

The borrowed continuum geometry and frozen-gas coefficient pair expose both
viscous and inertial resistance. Linear-only consumers and the molecular/
Knudsen dusty-gas closure explicitly refuse positive inertia. Slip plus inertia,
wet-pore inertia, momentum acceleration/choking, turbulence and extrapolation
beyond the supplied directional/fluid calibration are not supplied by this rung.
Focused kernel tests cover independent quadratic balances, extreme numerical
scales and tensor-direction admission. Original-owner regression sources cover
layer pressures, mass/drag work, zero-beta behavior, actual speed admission,
source support, confinement, seals and unsupported combined-law refusal.


## Formation-reference-independent real-gas nozzle energy

`ThermochemicalPengRobinsonModel::specific_enthalpy_change_between_densities`
combines the original NASA enthalpy increment with a finite PR enthalpy departure.
The prepared source reference retains both density coordinates and factors small
attraction, volume and logarithmic changes before combining them; no absolute
residual enthalpy or rounded pressure difference is subtracted. The nozzle uses
`velocity^2 = -2 Delta h`; absolute upstream enthalpy remains the stream datum.
The existing frozen-composition homogeneous source domain and bounded
entropy/sonic roots remain authoritative. Independent high-precision weak-change
references, endpoint reversal and formation-reference shifts test this coordinate.
Existing near-shut pressure-drive and variable-Cp PR throat references remain.

Real-gas energy/enthalpy trial samples use the original PR scalar density
evaluator and NASA property-only calorics. The repeated reconstruction path
therefore needs no component-vector allocation or fugacity calculation; it
retains the same residual U, complete Cv, pressure and temperature derivative.
The accepted material state still resolves the full EOS and original transport
sources, preserving source-domain, fugacity and complete caloric/acoustic checks.

Relative energy and entropy reconstruction, including simultaneous midpoint
heat/volume stages, use the same prepared finite owner and actual source Cv
tangents. Entropy combines source NASA changes at unchanged pressure, the ideal
density correction and the PR density departure, avoiding cancellation of
real-pressure and ln Z terms. Requested tolerances, cancellation, exact density
and full accepted material/transport checks remain unchanged. Cross-composition
receiving energy retains each mixture's own departure; the frozen-composition
finite owner is used only after that receiving composition is fixed.

## Reversible volume work on source-bound real gas

`IsochoricRealGasParcel::at_volume_isentropic` constructs an immutable reversible
adiabatic endpoint on the existing frozen-composition NASA/Peng--Robinson storage
owner. Exact mass, source/composition identity and extra query axes are retained.
The candidate density is `mass / requested_volume`; later `add_heat` remains
isochoric at that new volume.

The original safeguarded entropy root uses complete `Cv/T`. Source NASA relative
entropy cancels reference constants, and PR residual entropy follows
`S^R = (H^R - G^R)/T`. Returned `GasVolumeStep` work is minus accurate source-law
`m Delta u`, including departure changes without subtracting large absolute
formation energies. `work/Delta V` is the work-conjugate secant pressure; the
resolved material pressure is the endpoint pressure.

Every evaluated density/temperature retains the original positive-pressure,
covolume, mechanical-stability and complete-Cv admission. Final material resolution
also enforces original transport support. The caller owns applicability of a
continuous reversible homogeneous single-phase path: local stability does not
prove global phase equilibrium or exclude coexistence. Latent heat, phase
conversion, piston dynamics and an independent mechanical-work certificate are
outside this operation.

Invalid or unrepresented geometry, unbracketed or exhausted roots, source refusals
and cancellation return no candidate. Zero volume change is an exact identity
after control validation. Focused regression sources include an independent
75-digit dense-gas entropy oracle, pressure quadrature, the dilute gamma limit,
large reference shifts, heating at the accepted density and atomic refusal.


## Coherent porous and fluid-gap reflection

`porous::layers::RigidBackedLayers` retains the original porous laws and fluid
source, evaluates that fluid exactly once per prepared frequency, and composes
up to 64 ordered rigid-backed layers. Every porous response retains the original
frequency support and fluid-dependent JCA(L) density/bulk law. Lossless gaps use
that same fluid's density and sound speed; no new attenuation preset is supplied.

`PorousFluidLaw::effective_relaxing` and `Layer::RelaxingPorous` retain an original
`RelaxationState` and couple its complex Cp to the JCA(L) pore thermal function.
The underlying gas energy equation is
`(s rho Cp*(s) - kappa Laplacian) T' = s p'`, with `s=i omega`.
For JCAL, write `A=kappa phi/(rho k0')`,
`B=4 k0'^2 rho/(kappa Lambda'^2 phi^2)` and
`W=Cp* + A/s sqrt(1+s B Cp*)`; then
`K_eff=(p/phi)/(1-R/W)`. JCA uses `k0'=phi Lambda'^2/8`.
The original molecular calorics are conjugated once from their opposite phasor
convention. Complex Cp controls both pressure/temperature coupling and thermal
penetration, so free-gas attenuation is not added to a separately cooled pore.
The isothermal limit removes molecular compression loss; negligible wall
exchange recovers the original relaxing gas modulus divided by porosity.
Zero reservoir strength preserves the previous response bit for bit.

The relaxing gas must match the stack's one coherent equilibrium pore-fluid
snapshot. Porous density and drag remain the original geometry and mu/k0 law.
Frequency preparation caches the response for all reflection, transmission,
interior-field and heat consumers. Their bulk-compliance/thermal heat channel
includes both molecular relaxation and pore heat diffusion, without assigning
a gas/frame temperature partition. Exterior fluids and gaps remain explicitly
lossless. This ideal-gas, rigid-frame constitutive extension is not a new
measured gas/foam qualification. `porous_relaxation` tests use declared synthetic
reservoirs to check the isothermal, wide-pore and frozen-gas limits, passivity,
coherent field energy and layer-resolved R+T+absorption balance.

Tangential wavenumber and pressure/normal velocity continuity give exact planar
extended-reaction reflection within the admitted equivalent-fluid model.
Attenuating-wave transfer and homogeneous state normalization preserve thick
absorbers, thin loss, grazing phase and impedance poles without growing matrices.
Angle queries allocate nothing. Absorbed power is the actual normalized surface
work, not a cancellation-prone subtraction from unit reflected power; values are
not clipped. All-gap stacks retain exactly zero dissipative work.

Finite normalized admittance is a separate fallible view, with exp(+i omega t)
convention and no area factor. A local Helmholtz boundary explicitly chooses the
normal-incidence view and converts the phasor convention once. This does not
supply a nonlocal finite-panel operator, Biot frame motion or gap wall losses.
The focused `porous_layers` regression source covers original-owner equivalence,
independent transfer matrices, coherent source sampling, layer order/subdivision,
lossless resonances, thick/grazing limits and sub-roundoff absorption.

## Native Chaboche cyclic hardening and source binding

`plastic::chaboche::ChabocheJ2<N>` retains up to eight independent tensor
backstresses over the original Voce and optional linear Prager law. Each
explicit branch follows `dX_i = (2/3) C_i d(epsilon_p) - gamma_i X_i d(alpha)`
with finite nonnegative C_i and gamma_i. Backward Euler eliminates each branch
as `X_i,new = (X_i,old + (2/3) C_i delta_epsilon_p)/(1+gamma_i delta_alpha)`.
The bounded scalar consistency solve updates its flow direction as the plastic
increment changes; the exact tensor tangent includes that direction derivative
and is generally nonsymmetric. A symmetric linear solver cannot assume otherwise
from the isotropic elasticity. Trial responses keep incoming history unchanged.

For each active branch, recoverable storage is `3 X_i:X_i/(4 C_i)`. Its physical
dynamic-recovery dissipation is `3 gamma_i delta_alpha X_i,new:X_i,new/(2 C_i)`;
this is already included once in the physical total alongside
`sigma_y delta_alpha`. The separate backward-Euler term
`3 delta_X_i:delta_X_i/(4 C_i)`, elastic/linear hardening remainders and the
convex Voce remainder remain numerical loss. Neither stored hardening energy
nor numerical loss is relabeled as heat. The energy API checks endpoint yield
and work balance with explicit physical-unit tolerances. A heat-conversion
fraction or thermal evolution requires a separate declared coupling.

`plastic::chaboche::axial::UniaxialChaboche` restricts this same consistency
equation to free-sided axial stress by replacing the elastic correction 3 mu
with E. It preserves every backstress, lateral elastic contraction and plastic
incompressibility; the reconstructed full tensor must satisfy the declared
small-strain ceiling, never exceeding 0.05. Zero recovery gives an additional
linear kinematic branch; absent branches recover the original native Voce law.

`plastic::chaboche::source` binds one independently sourced branch C/gamma and
retains the complete original Voce material bundle, quantities and claim-use
receipts, including unrelated physics roles. Both coefficients are mandatory
even when explicitly zero. A monotonic total-strain tangent cannot identify
nonlinear recovery. Wrong units, missing or ambiguous data, extrapolation and
evolving strain/stress/history/rate/frequency coefficient coordinates refuse.
Supported temperature curves may resolve at a fixed external state; every
coefficient remains fixed throughout the load path. Re-resolution does not
transfer history or establish a temperature-dependent energy law. Static
recovery, time-dependent creep, finite strain, damage and experimentally
calibrated ratcheting or fatigue life remain outside this implementation.



## Nonideal aqueous thermodynamics and finite-bath Donnan equilibrium

`osmosis::nonideal::VirialSolutionModel` prepares an explicit, bounded,
isothermal free energy on solvent-volume concentrations. The symmetric input
matrix B [m³/mol] enters `F_mix/(RT V) = sum(c_i*(ln(c_i/c_ref)-1)) + c^T B c/2`;
the conventional osmotic second-virial matrix is B/2. Pressure, water activity,
solute chemical potentials, their Jacobian and dilution stiffness derive from
that same potential. Native Cholesky certifies a numerically convex declared
concentration box at preparation. Coefficients, concentration support, species
ordering and fixed temperature are caller supplied; none establish experimental
material authority. Exact absent species have finite solvent thermodynamics;
singular solute logarithms/derivatives return named refusals, never floors.

`osmosis::AqueousSolution` retains the original finite inventory with either
ideal thermodynamics or an immutable shared virial model. Trial pressure and
atomic inventory replacement allocate nothing. Existing ideal inventories
convert directly. Owned volume replacement copies the amount array. Consumers
must use the selected model's activity/pressure consistently: a fully retaining
membrane or vapor-equilibrium boundary is supported, while the existing ideal
leaky membrane and Fick film laws are not silently made nonideal.

`osmosis::donnan::IdealDonnan` solves ideal, fast-ion electrochemical equilibrium
in two finite baths. Every mobile ion's total moles are conserved; fixed charge
is supplied separately in mol of elementary-charge equivalents. Independent
trapped osmotic particles are distinct from fixed charge. Both baths must admit
electroneutrality without fabricated counterions. Equal valences share a
monotone scalar solve, with stable minority partitions, bounded work and
cancellation. Counterion-only endpoints return typed infinite-potential limits
with finite physical inventories and pressure. Both output arrays publish
atomically. The report includes full ionic-plus-trapped pressure, charge
residuals and a nonnegative water-transfer stiffness when representable.

These are separate admissible constitutive closures. Ideal Donnan partition
must not be mixed with a nonideal pressure correction alone. Neither owner
supplies association, precipitation, double-layer storage, imposed currents,
transient Nernst--Planck transport or solute partial-volume corrections. Focused
native regressions exercise thermodynamic derivatives, convexity including the
64-species blocked factorization, analytical finite-bath partitions, multivalent
charge balance, salt-free limits, free-energy minimization and atomic refusals.

### Finite-inventory dense membrane breakthrough

`membrane::transient::ResolvedMembraneDiffusion` adds isothermal transient
permeation for one pure permeant in a planar laminate. Original material cards
separately supply `gas_membrane_diffusivity` [m²/s] and
`gas_molar_henry_solubility` [mol/(m³ Pa)]; `P = D S` uses the same typed
permeability/permeance owner as steady transport. No intrinsic porous k, mass
basis, STP volume or measured diffusion coefficient is inferred from permeability.

Cell-centered finite volumes retain actual `S dx` storage and harmonic
half-cell/interlayer resistance. The existing tridiagonal backward-Euler solver
owns elimination; equal-duration steps on the same immutable grid reuse its
factorization. Each immutable candidate reports actual face fluxes, reservoir
transfers, stored amount/change, reconstructed Henry pressures and local/global
amount-balance residuals. Stability does not imply resolved breakthrough:
first-order time stepping and spatial resolution require refinement. Geometry
and work are bounded to 64 layers and 4096 cells, with bounded cancellation polls.

Original material and ordered-interface claims stay pinned and are requalified
at actual old/new cell and face pressure extrema. Unsupported source states,
unrepresentable arithmetic, failed amount budgets and cancellation publish no
physical state. Exact zero D or interface permeance closes that path while
retaining authored loading, including inventory trapped between multiple seals.
Positive Henry solubility is required for finite-storage cells. An undetermined
zero-volume trace between coincident seals cannot qualify an explicit interface
source and is never assigned an invented pressure.

Reservoir pressures and uniform temperature are prescribed externally. This
rung supplies no finite-gas energy coupling, heat of sorption, competitive mixed-
gas uptake, nonlinear sorption, swelling or concentration-dependent diffusivity.
Analytical regression sources compare initial delay, Fourier-series breakthrough
and storage, temporal/spatial convergence and the homogeneous `L²/(6D)` lag;
other cases cover heterogeneous/interface steady limits, Henry partitioning,
unloading, seals, amount conservation and transactional cancellation. Source
regressions check original receipts, local domains, ordered identities and
unsupported coefficient bases without claiming measured membrane accuracy.

### Finite immobile trapping during membrane breakthrough

`membrane::transient::slab::trapping` couples mobile Henry storage S*q to an
explicit finite Langmuir site population in each selected layer. Only mobile
solute transports through the original D*S face law. The immobile pool obeys
finite adsorption/desorption kinetics with explicit capacity per membrane
volume, pressure affinity and desorption rate. Initial trapped loading is an
independent declared inventory. Saturation, delayed breakthrough, preloading
history and unloading emerge from the coupled conservation equations.

The allocation-free scalar limit of the original competitive-site owner
analytically eliminates each backward-Euler site inventory and supplies its
positive rate derivative. The original multicomponent site entry point reuses
this same one-species limit. Remaining mobile pressures and optional finite
pressure-capacity reservoirs share one nonlinear tridiagonal solve using native
fs-la LU, bounded Newton iterations and positivity-preserving backtracking.
There is no diffusion/uptake splitting or transport of the trapped pool.

Reusable O(cells) scratch, bounded cancellation and immutable candidates retain
transactional failure. Original local, membrane, reservoir and closed-system
amount checks include mobile and trapped changes. An independent scalar-site
kinetic residual must also fit the same increment/flux-based budget; stored
conservation alone cannot conceal uptake rounded away beside a large preloaded
site. Reports distinguish both inventories, gross adsorption/desorption and
boundary/interface transport.
The coefficients and temperature are fixed; competitive trap populations,
swelling, sorption heat and a total energy closure are not inferred. A single
implicit step is first-order and requires time/space refinement. Authored
membrane and scalar-site regressions cover independent rational/quadratic
solutions, the native Henry limit, finite nonlinear equilibrium, saturation-
delayed breakthrough, preload/desorption history, seals, sub-ULP kinetic uptake
refusal and atomic nonlinear/cancellation failure; these are not measured
polymer validation.

### Prepared shared-site kinetic response

`adsorption::competitive::sites::backward_euler_sites_into` evaluates the
original common-vacancy backward-Euler equations into a reusable bounded
`CompetitiveSiteBuffer`. All species consume one finite capacity. The existing
allocating entry point delegates to this same owner; its scalar limit retains
the original scalar arithmetic and derivative.

Occupied amounts, gross streams, vacancy and independent per-species kinetic
residuals remain distinct. Analytic adsorption-rate derivatives use diagonal-
minus-rank-one factors with O(species) storage. A dominant species retains its
small diagonal response through a complementary sum. Evaluation stages on the
stack and publishes only after admission and cancellation checks. Capacity-
scaled closure remains separate from transport's increment/flux-scaled kinetic
accuracy, so capacity rounding cannot obscure a physical uptake discrepancy
or force refinement of an otherwise representable tiny uptake.

### Competitive multicomponent membrane diffusion

`membrane::transient::slab::competitive::CompetitiveMembraneDiffusionState`
couples 1..32 explicit species on the same ordered membrane geometry and clock.
Each keeps its original Henry storage and D*S face law. Adsorbing species share
ONE exactly agreed site-capacity density in each layer; capacities are not
summed. Original species IDs determine canonical channel and boundary order.
Finite shared vacancy, saturation and displacement evolve with mobile storage
in the same implicit solve. Only mobile material crosses the original faces.

The prepared competitive-site response eliminates local occupied endpoints and
supplies their cross-species derivatives. Bounded Newton/backtracking uses
native per-species tridiagonal factors and block Gauss-Seidel corrections.
Independent reconstruction of the complete Jacobian residual admits each
correction. Each sweep costs O(cells*species); no dense mesh Jacobian or per-cell
kinetic allocation is introduced. An explicit byte ceiling covers retained
workspace, with bounded nonlinear iterations, linear sweeps and backtracks.

Prescribed partial pressures or `FiniteCompetitiveMembraneDiffusionState`
with explicit whole-bath molar pressure capacities are supported. Every species
independently audits local mobile/trapped storage, original face transfers,
kinetic residual and finite-bath closure. Aggregate mass cannot hide conversion
between species. Shared occupancy cannot exceed capacity, and any work,
representability or cancellation refusal publishes no partial species state.

Fluxes are relative to a stationary rigid membrane matrix; this is not a
zero-sum Maxwell-Stefan gas-mixture frame. Temperature, geometry, capacities,
Henry coefficients and kinetics remain fixed. Original-source mixed-gas
qualification, gas EOS/calorics, swelling, sorption heat, cross diffusion and
adaptive competitive windows are not inferred. First-order time and spatial
refinement remain caller responsibilities.

Focused independent checks cover a common-vacancy finite-bath root, nonlinear
competitive equilibrium, invading-species displacement, Riccati time refinement,
canonical species permutation, original scalar/nonadsorbing limits, sub-ULP
uptake refusal, bounded work and late cancellation. These are numerical model
checks, not measured multicomponent membrane qualification.

### Adaptive finite-trap membrane windows

`membrane::transient::slab::trapping::adaptive` reuses the original membrane
full-step/two-half-step controller for simultaneous Henry diffusion and finite
Langmuir kinetics. Mobile pressure, mobile amount, trapped amount, finite
reservoir inventories and integrated face transfers enter separate local
accuracy tests. A constant total loading cannot hide unresolved exchange
between mobile and trapped stores. Only accepted native halves contribute to
the physical timeline and signed/gross face totals.

The shared native adapter exposes optional immobile inventories and may tighten
iterative algebraic controls; original controls are never relaxed. Per-probe
and original-to-final balances include both stores. Reports retain separate
mobile/trapped amounts and changes, plus finite-reservoir and closed balances.
Reusable native workspaces, explicit evaluation/iteration/backtrack limits and
bounded cancellation retain transactional failure. The empty, zero-adsorption
scalar limit keeps exactly zero loading and exactly full vacancy without
discarding its nonzero adsorption-rate derivative.

Focused G0/G1/G3 regressions compare sealed exchange against its continuous
Riccati solution, refine finite-reservoir transport against native steps, check
work/native/cancellation refusal and exact target times, and exercise the
scalar empty-site limit. These local temporal estimates do not certify a
whole-window or spatial error bound, measured membrane accuracy, changing
temperature, sorption heat, competitive species or swelling.

### Original-source binding for finite membrane traps

`membrane::transient::trapping::ResolvedTrappingMembraneDiffusion` attaches
original `LangmuirAdsorbent` sources to the original prescribed-pressure Henry
membrane. Explicit dry-skeleton density converts the original mol/kg capacity
to trap mol/m3 of the complete layer. Distinct material cards explicitly declare
composite constituents; no porosity or mobile transport coefficient is inferred.

`ResolvedFiniteTrappingMembraneDiffusion` binds the same original film to two
explicit whole-reservoir molar pressure capacities. Both reservoir inventories,
mobile Henry loading and immobile Langmuir loading share the original native
implicit solve. The wrapper owns one native finite state; immutable source
metadata and receipts do not duplicate its physical inventories. Construction
requalifies the actual finite initial pressures, including baths different from
the previous prescribed boundaries.

Each candidate requalifies original mobile D/S, ordered interfaces and trap
capacity, affinity and kinetics across actual old/new cell and layer-face
pressure envelopes. Complete typed physical contexts must match. Unowned
composition/concentration coordinates refuse in this pure-permeant bridge.
Original source identities and both families of receipts remain retained.
Source failure or cancellation publishes neither inventory, clock nor receipts.
These fixed-temperature owners do not infer sorption heat, gas-mixture
composition, an EOS, caloric storage or thermostat work.

Both source owners expose adaptive windows through the SAME kinetic time
controller. Every coarse and fine native probe passes original source admission;
an unsupported fine release transient cannot hide behind a supported coarse
endpoint. Source refusals are hard errors, not time-step retry signals.
Only accepted half steps advance physical inventories and integrated transfers.
The returned state retains the final accepted step's support and all original
bindings. Local temporal estimates are not whole-window or spatial certificates.

Langmuir source-property evaluation admits pressure zero only when its original
selected material domain covers it; actual gas-parcel pressure admission is
unchanged. Focused G0/G1/G3 source-stack regressions exercise source-qualified
vacuum, native/source trajectory equality, held humidity/mixture context refusal,
finite-bath rebinding and conservation, preloaded release beyond source support,
adaptive raw/source equality, fine-probe source rejection and transactional late
cancellation. These source-owner tests do not establish measured material
qualification or full workspace/application validation.

### Original gas parcels around finite trapping membranes

`membrane::transient::trapping::reservoir::IsothermalTrappingGasMembrane`
couples the original source-qualified mobile/trapped film to two original pure
`ClosedGasParcel` reservoirs at one exact held temperature. Shared gas admission
supplies finite capacities, species/caloric/transport/pack identity and gas
enthalpy. Each original adsorbent must independently admit that same catalog;
a shared species identifier alone cannot authorize a substituted caloric law.

Each step uses the original source-relative gas boundary owner and its exact
isothermal flow-work offset. Actual gas mass and pressure closure must pass
before the native reservoirs are reseeded from actual EOS pressures; membrane
inventories are moved without copying or changing them. Original material
support covers old, numerical and actual-EOS cell/face pressure envelopes, and
the selected trap coefficients must remain unchanged.

Its adaptive adapter performs this complete gas/source reconstruction on every
coarse and fine probe. Unsupported gas/material sources remain hard errors.
Accepted face totals determine net fixed-temperature gas enthalpy and thermostat
streams, with an independent audit of actual gas, mobile and trapped mass.
Fixed steps retain the original gas caloric reports; a window does not fabricate
a summed gas-energy residual or whole-film energy closure. Preparation and
failure retain original gas sources, volumes, temperature and accepted state.

Focused original-source regressions compare nontrivial no-trap flow against the
existing Henry gas owner, exercise trapped uptake and adaptive conservation,
preserve trajectories under a 1e12 J/kg caloric-reference shift, and refuse
unsupported gas domains, foreign trap caloric laws and late cancellation.
No mixed-gas competition, mobile/trap caloric storage, sorption heat or
nonisothermal membrane is inferred.

### Accuracy-controlled transient membrane windows

`membrane::transient::slab::adaptive` supplies bounded backward-Euler
step-doubling windows for layered membranes with fixed prescribed pressures
or finite molar pressure capacities. Acceptance compares every cell/reservoir
pressure, stored amount and integrated face transfer. Accepted half steps alone
advance inventory. Compensated signed/gross face totals are audited against
original-to-final cell, reservoir and whole-window balances.

Every native probe retains its original admission and additionally reserves
1% algebraic residual allowances below pressure/storage accuracy and
duration-weighted window balance shares. Separate reusable full/half caches
retain the native O(n) factorizations. Work includes rejected probes, and the
exact representable requested endpoint is reached or the complete window
refuses without changing accepted state. An immediately adjacent final clock
remainder is included before solving, never silently discarded.

`ResolvedMembraneDiffusion::advance_candidate` uses that same driver and the
original source-bound step on every probe, preserving selected D/S/interface
laws and treating source failures as hard refusals. Temperature and constitutive
coefficients remain fixed; callers split windows at prescribed boundary events.
Local temporal estimates do not certify spatial or accumulated global error.
Authored tests compare Fourier breakthrough/cumulative permeation and exact
finite-reservoir eigenmodes, amount sensitivity, sealed equilibria,
loading/unloading, cache reuse and atomic work/accuracy/algebraic/cancellation
refusal, with separate source-bound receipt and support regressions.

## Debye dielectric temperature and thermal potentials

`field_storage::dielectric::thermal::ThermalDebyeCapacitor` binds the original
complete Debye source and fixed geometry to one typed absolute-temperature
interval. Instantaneous and polarization permittivities must be positive affine
source paths over that whole interval, including every interior knot. Relaxation
times may follow positive continuous piecewise-linear paths; their temperature
dependence changes friction, not stored energy or entropy. Every refresh retains
the original card, selected claims, pole roles, quantity conventions and fixed
nonthermal coordinates, and re-evaluates all other properties in the bundle.

The shared field thermodynamic primitive evaluates each independent charge's
`psi=q²/(2 C)`, `S=psi C'/C`, `U=psi+T S` and
`C_q=-2 T psi (C'/C)²`. Bounded evaluations allocate no scratch vectors. The
analytic charge-envelope bound on the complete negative field heat capacity
must be combined with the original zero-field caloric chart; the coupled owner
requires positive total capacity and enforces every charge bound. This layer
does not advance electrical or thermal history.

`ResolvedDebyeCapacitor::at_query_point` exposes the same complete source refresh
as an endpoint law. `fixed_charge_free_energy_change_j` admits temperature-only
parameter changes at fixed original charges and returns signed Helmholtz work
using relative capacitance differences, retaining small work beneath a large
unchanged energy reservoir. `thermal_offset_change_j` independently resolves
the change of `T S` through stable charge and temperature increments and
compensated pole sums. Both remain separate from native irreversible resistor
heat. Unsupported storage kinks, changed source/geometry/pole correspondence,
missing histories and nonthermal refreshes refuse; measured material accuracy,
moving geometry, switching impulses and general non-affine caloric laws are
outside this model. `tests/dielectric_thermal.rs` supplies independent analytic
potential/derivative and capacity checks, source/limit refusals, and large-bias
parameter-work and entropy-offset regressions with synthetic material data.


## Real-gas simultaneous heat and volume work

`IsochoricRealGasParcel::midpoint_heat_volume` returns the common
`GasHeatVolumeCandidate<_, ResolvedThermochemicalRealGasState>` for signed heat
and volume change at conserved exact mass and frozen composition. It solves
midpoint in internal energy and volume at `Vm=(V0+V1)/2`, uses
`Wgas=pm*(V1-V0)`, and reconstructs the endpoint through the shared
between-densities energy-increment solver with the original initial density
retained. Both accepted stages require full material and transport admission.

NASA reference constants cancel in source increments. Original PR departures
account for both temperature and density; entropy uses the native PR departure
owner. Analytic derivatives of the same discrete midpoint pressure use complete
Cv and original pressure derivatives. With `D=Cv+DeltaV*p_T/(2m)>0`, they are
`dp_mid/dV1=(p_nu*Cv-T*p_T^2)/(2m*D)` and
`dp_mid/dQ=p_T/(2m*D)`, using `du/dnu|T=T*p_T-p`.
The heat tangent retains the source-derived sign. Gas and stage-energy residuals
must meet absolute tolerances before any candidate is returned.

The caller owns the continuous homogeneous single-phase path and timestep
refinement. Local stability does not establish phase equilibrium or latent heat;
actual entropy change is reported without asserting exact finite-step isentropy
or nonnegative entropy production. Identity preserves the original state, and
cancellation/source/domain/iteration refusal returns no candidate. The common
candidate's default generic specialization preserves the existing ideal API.
Focused regression sources cover independent dense PR oracles, dilute limits,
discrete derivatives, reference invariance, fixed-volume heat, adiabatic
refinement, both stage domains and cancellation.

## Source-bound pure-fluid coexistence and fixed-volume calorics

`thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel` now exposes
`pure_saturation_at_temperature` through the original thermochemical EOS owner.
The accepted coexistence roots retain the original NASA-9 evaluation and receipt
at their common temperature and pressure. An explicit molar volume then fixes
the vapor amount by the original lever-rule closure. Full mass-specific h, u
and s combine that one ideal reference with the amount-weighted PR departures;
no second latent-heat source or entropy of mixing identical phases is added.

The existing source temperature domain remains authoritative. Gas transport
cards are neither required nor averaged for this caloric evaluation, and the
result does not implement a two-phase transport, heat-capacity or acoustic law.
Volumes outside the resolved dome and multicomponent saturation requests refuse.
The numerical coexistence controls, work cap, unresolved critical-region
refusals and original phase roots remain owned by `fs-thermochem`.

Focused G0/G3/G4 regressions check original source receipts, independently
computed latent h/u, common phase Gibbs energy, h-u=p/rho, conservative phase
amount, energy-reference shifts counted exactly once, and source/domain
refusals. The synthetic caloric reference tests implementation consistency;
they do not assert experimental methane caloric accuracy.


## Coherent porous-layer reflection and transmission

`porous::layers::TransmittingLayers` reuses the existing prepared layer laws and
one original incident/saturating-fluid snapshot. A borrowed `FluidAcoustics`
source supplies one receiving snapshot, explicitly modeled as a lossless
homogeneous half-space. The original rigid-backed API and source domains remain.

Pressure R is referenced to the entrance, pressure T to the exit. Conserved
tangential wavenumber, outgoing causal normal wavenumber and pressure/normal
velocity continuity retain oblique refraction, critical incidence and evanescent
receiving fields. An evanescent field has complex pressure T and zero outgoing
normal power. Transmitted power uses its actual normal-admittance weight.

Bounded transfer retains the common complex amplitude with compensated log
magnitude and phase; it neither infers T from R nor treats an unspecified backing
as opaque. Weighted power is evaluated before squaring can underflow. Queries
allocate nothing and preserve the 64-layer limit. All-lossless-fluid stacks have
exactly zero absorption. For absorbing transmitting stacks, absorption is entrance work
minus outgoing power: sub-roundoff loss is unresolved, and values are not clipped.
The focused regression source covers independent matrices, reciprocity,
subdivision, critical/TIR phase, matched delays, thick absorption, extreme
contrast, original-source sampling and admission.

The normal-incidence local-boundary view uses
`TransmittingLayers::normalized_surface_admittance` directly from the same
scaled entrance pressure/velocity pair; it does not reconstruct a weak load
from rounded reflection. Native beta is `Z0*v_n/p` in exp(+i omega t), and a
true entrance pressure node refuses a finite Robin representation. Ordinary
reflection, transmission and admittance queries perform no heat integration.

The explicit `normal_surface_power` query returns `LayerSurfacePower` with
independent irreversible-loss and outgoing-receiver conductances. Each
coefficient multiplies total entrance `|p|^2/(2*Z0)` to obtain W/m². The original
positive viscous/thermal integrals are normalized to that entrance pressure
before any squared-field underflow; heat is not a reflected/transmitted power
difference. Compensated log amplitudes retain receiver flux independently.
The residual against `Re(beta)` remains visible and is never redistributed.
Matched lossless backing has zero heat while accepting transmitted power.
This query retains the native explicit layer allowance, fixed storage and
cancellation in both propagation and loss passes. Focused regressions cover
weak receiving conductance, sub-roundoff porous heat, matched lossless loading,
native incident-power equivalence, pressure-node refusal and bounded work.
This is a planar local normal-incidence power account, not a finite receiving
geometry, transient thermal solve or elastic-backed extension.

`Layer::FluidLayer` additionally retains a borrowed original `FluidAcoustics`
source for each finite homogeneous layer, so differing density and sound speed
participate in the same oblique reflection/transmission solve. Repeated uses of
one source share its prepared snapshot, including reuse as the incident or
receiving half-space source. Distinct sources are each sampled independently.
`FluidGap` still means the original
incident fluid, and porous constitutive responses still use the original pore
fluid. No viscosity, thermal loss or source calibration is inferred from the
lossless finite-layer view.

`Layer::RelaxingFluid` carries the original `RelaxationState` through these
same finite layers without a porous skeleton. Its equilibrium gas supplies
the real density; the original molecular reservoir sum supplies the complex
bulk modulus, conjugated once into the layer owner's `exp(+i omega t)`
convention. Frequency-dependent phase, refraction, attenuation and positive
compression heat therefore remain one constitutive response. The numerical
`PorousFluidResponse` view has exactly zero static resistivity and real density
in this source role. No shear viscosity, Fourier heat diffusion, interface
thermal mode, wall loss or artificial attenuation is inferred. Lossless
exterior fluids retain their declared approximation.

Repeated references share one prepared dynamic response. With no active
reservoirs the existing lossless `FluidLayer` path is used exactly, retaining
critical and evanescent limits. Coherent exterior/interior beams, elastic
terminations, layer-resolved local heating and angular power all consume the
existing bounded propagator and original source identity. Free-layer tests
reconstruct the molecular modulus and conventional transfer matrix
independently, integrate local heat against intensity loss, check differing
fluid states, zero-reservoir/zero-thickness limits, splitting invariance,
opaque attenuation and cancellation. This remains a linear stationary
infinite-plane molecular model, without a measured-spectrum qualification or
transient source feedback claim.

Critical-angle classification retains the original sound-speed pair. Factored
source contrast preserves nearly matched grazing propagation, and only
cancelling critical evaluations use the existing double-double precision ladder
with exact power-of-two speed scaling. A rounded speed ratio cannot decide
whether a channel carries outgoing power. Equal-medium and normal-incidence
fast paths remain exact; no angular band is clipped to artificial criticality.

`porous::layers::ElasticTransmittingLayers` terminates the original finite
porous/fluid stack in an explicitly isotropic, homogeneous, unstressed, lossless
`ElasticAcousticState` half-space. Shared tangential slowness, continuous normal
traction/velocity and zero shear traction determine outgoing P/S fields and
separate normal powers. A normalized boundary pressure/velocity pair retains
critical and evanescent modes and Rayleigh surface-pressure zeros without
impedance division. Original density, tensor, identity and response regime are
retained; Christoffel speeds are prepared once and angle evaluation allocates
nothing. Existing bounded layer attenuation supplies only finite-layer loss.
This infinitesimal-wave model does not project anisotropy, supply finite solid
slabs, move a porous frame or infer damping from static source moduli.

Normal-wave matching retains finite transmission through evanescent layers
(frustrated total internal reflection). Such decay is not heat. At critical
incidence, the transfer uses the finite pressure/normal-velocity limit instead
of dividing by zero normal wavenumber. Bounded decaying matrices preserve thick
evanescent layers without growing exponentials. Focused G1/G3 regressions cover
independent heterogeneous fluid/porous matrices, analytic evanescent tunneling,
critical-angle continuity, quarter-wave impedance matching and original-source
reuse/refusal. `layered_acoustic_refraction` is a runnable warm-air-layer example
using the existing gas-state owner. Infinite planar, linear, stationary-fluid
and rigid-frame scope remains; finite-aperture boundary conditions, shear waves,
flow and continuous sound-speed gradients are not supplied by this layer model.

### Coherent angular acoustic fields

`TransmittingLayers::prepare_angular_spectrum` binds an explicit finite set of
incident plane-wave directions and complex peak pressures to the SAME prepared
layer scattering response. Each component retains its conserved tangential
wavenumber and original complex pressure reflection/transmission; no coefficient
is inferred from a power fraction. Original incident, layer and receiving fluid
snapshots are reused without additional material queries.

`PreparedAcousticSpectrum::incident_field` observes the incident plus reflected
field at `z <= 0`, measured from the entrance face. `transmitted_field` observes
the receiving half-space at `z >= 0`, measured from the exit face. Coordinates
are metres and phasors use `exp(+i omega t)`. Pressure and all three particle
velocity components add coherently BEFORE evaluating the complex intensity
`p * conj(v) / 2`; its real and imaginary parts retain active and reactive
intensity in W/m2. Thus destructive interference and evanescent reactive fields
are not replaced by an incoherent power average. The original outgoing normal
wavenumber retains critical incidence and exponential evanescent decay.

Preparation has explicit component/work limits and O(component count) storage;
observation is O(component count), allocates nothing and checks cancellation.
This is an explicitly supplied finite angular sum, not an automatically sampled
aperture, a certified angular quadrature or a finite-specimen diffraction model.
The existing incoherent angular-power API
retains its separate meaning.

`PreparedAcousticSpectrum::prepare_layer_field` optionally prepares one original
finite layer's coherent pressure, three-component particle velocity and complex
intensity. Its `field` query uses depth from that layer's front face, includes
both faces, and keeps the original transverse phase origin. One-sided material
traces retain pressure and normal-velocity continuity without incorrectly
equating the tangential velocities across a density jump. A zero-thickness
identity admits only depth zero. Porous velocity retains the equivalent-fluid
volume-flow convention and the original complex dynamic density.

Optional preparation costs O(components × layers) through the original bounded
scattering owner and stores O(components) for the selected layer; exterior-only
spectra incur no layer cache. Queries allocate nothing, re-query no source,
perform O(components) local matrix work and poll cancellation. Full/remaining
normalization ratios cancel exponential growth analytically, preserving finite
entrance fields even when transmission through the whole layer underflows.
The original finite critical-angle matrix limit is reused without dividing by
normal wavenumber. Focused controls in `porous_angular_spectrum` compare
independent propagating/critical/evanescent interior fields, coherent interface
work, opaque penetration, and flux divergence against the original local
viscous/thermal loss law. This observation does not evolve thermal state or
extend the planar equivalent-fluid model's physical validity.

The `layered_acoustic_refraction` example accepts a fourth argument `field` to
emit a coherent pressure, velocity and active/reactive intensity CSV map through
the original warm-air layer. Its nine directions and binomial pressure weights
are explicit; the original three-argument angle/power table remains available.
The fourth argument `interior` instead resolves the same beam inside that warm
layer, including its two boundary faces and evanescent/standing fields.

### Finite acoustic sensor footprints

`RectangularAcousticWindow` adds an exact uniform rectangular observation
footprint parallel to the original planar interfaces. Exterior windows use
the corresponding entrance/exit point-field origin; interior windows use depth
from the selected layer's front face. Analytic sinc overlaps integrate the
original coherent tangential phases, yielding spatially averaged complex
pressure and velocity, separately integrated time/spatial mean squares, and
active/reactive intensity. Squaring an averaged pressure does not substitute
for the integrated acoustic energy.
Double-double wavevector differences, sinc decrements and quadratic weights
retain positive narrow-window energy and heat at a dark interference fringe.

`PreparedLayerAcousticSpectrum::window_with_heat` integrates the original
coherent viscous and thermal constitutive heating, including all velocity
directions; heat is not reconstructed from spatially averaged fields.
Zero-thickness layers admit interface fields but refuse heat density. The
observation uses O(N) scratch storage and explicitly budgeted N(N+1)/2 component
pairs, polls cancellation during component preparation and pair integration,
and performs no source queries or new boundary solve. It remains a passive
sensor footprint on the original infinite planar solution; it adds no aperture
clipping, edge diffraction, temporal exposure or thermal evolution. Nonfinite
or negative quadratic results refuse rather than being presented as physical
power.

### Local coherent acoustic heating

`PreparedLayerAcousticSpectrum::field_with_heat` observes the same original
interior pressure/velocity field and separately reports viscous and thermal
heat density in W/m3. The original dynamic density and bulk modulus supply
`-omega Im(rho) |v|^2 / 2` and `-omega Im(1/K) |p|^2 / 2`, respectively.
Fields are summed coherently before evaluating these terms; transverse
velocity can therefore dissipate at a pressure node, while identical
opposite-phase sources cancel both heat mechanisms. Lossless layers have
exactly zero heat. Reactive intensity and evanescent decay are not assigned heat.

Two cached constitutive weights add constant work after the original field
observation, with no allocation, new solve or material query. Storage inherits
the original spectrum component budget. The existing `prepare_layer_field`
and `field` APIs, layer ordering, zero-thickness interface observation and
critical/opaque propagation remain unchanged. A zero-thickness identity has
no material volume and refuses the heat-density query. Cancellation publishes
no partial result and leaves the prepared field reusable.

These are equivalent-fluid dissipation mechanisms, not a fluid/frame thermal
partition, temperature evolution or experimental calibration. Focused G0/G1/G3
regression source compares independently derived complex boundary equations,
integrates sampled heat against the original analytic layer-absorption owner,
preserves coherent nodes/cancellation, and checks lossless/source/cancellation
semantics.

The `porous_acoustic_heat` example emits a complete interior CSV map for two
coherent rays in the original Basotect Foam-02 layer with explicitly authored
characterization viscosity and the USSA-1976 dry-air state. Frequency, thickness,
pressure scale and relative phase are arguments. It resolves one transverse
interference period, reports all pressure/velocity/intensity components and
both heat mechanisms, and buffers the bounded map before output. The sampling
is a visualization grid, not an asserted spatial integration accuracy or a
thermal-feedback model.

### Diffuse and authored angular incident power

`RigidBackedLayers::diffuse_power` and `TransmittingLayers::diffuse_power`
integrate their SAME prepared extended-reaction responses over an isotropic
incident hemisphere or an explicitly truncated cone. With incidence cosine
`mu`, the incident normal-power measure is `mu dmu`; polar midpoint panels
retain this projection exactly once. Original internal refraction, porous loss
and receiving normal-admittance weighting are unchanged. Independent incoming
directions add powers, never complex reflection amplitudes. Evanescent receiving
fields carry zero transmitted normal power and do not become heat.

`angular_power` instead consumes explicit incoherent directions and their
already projected normal incident watts, including any authored solid-angle
weights. Both paths independently accumulate reflected, transmitted and absorbed
work with compensated sums. Absorption uses the original positive analytic local
viscous/thermal integrals, including tangential drag, rather than a rounded
entrance-minus-transmitted power difference. No clipping or energy-balance
renormalization is applied. Queries allocate nothing or re-query material
sources, poll cancellation per sample and per layer in both native passes, and
obey an explicit sample allowance plus a 65,536-sample hard cap.
All-zero authored weights refuse because no normalized distribution exists;
a declared diffuse distribution at zero watts retains coefficients and zero work.

`angular_layer_absorption` and `diffuse_layer_absorption` retain these same
independent powers plus each original layer's viscous/thermal absorbed fractions
and watts in `AngularLayerAbsorption`. Fractions use total incident angular power;
lossless and zero-thickness rows stay exactly zero. An explicit layer allowance
bounds each angle's propagation/loss work. Fixed storage and O(angles × layers)
work preserve sub-roundoff heating without a second material or field solver.
This is a frozen-material heat account, not transient temperature feedback or
separate pore-fluid/frame calorics. Focused G1/G3/G4 controls cover independent
thin-layer angular and exact diffuse limits, spatially integrated directional
heat, subdivision, zero power, retained sources and cancellation/layer-budget
refusal.

Polar resolution is a convergence axis, not a certified error tolerance. The
`porous_layers` regression includes a three-rung comparison against independent
transfer matrices integrated in polar angle, a truncated cone, separate
directional power controls, exact lossless/matched limits, critical-refraction
power bounds, original-source retention and cancellation/budget refusal. Infinite
planar stationary rigid-frame assumptions remain; this does not predict finite
specimen diffraction or certify a reverberant room's diffuse field.

### Layer-resolved viscous and thermal acoustic absorption

Rigid-backed, fluid-transmitting, elastic-half-space and finite-elastic-plate
prepared layers expose `layer_absorption` for explicit incident normal watts at
one angle. Each original entrance layer retains separate viscous and thermal
absorbed fractions and their combined watts. The original R/T API
is unchanged. Its SAME bounded propagation supplies the layer-face states;
no additional material query, fitted attenuation or second field solver enters.

Work follows `omega/2*(-Im(rho)*|v|^2 - Im(1/K)*|p|^2)` integrated through
each retained equivalent fluid. Oblique viscous loss includes tangential
velocity. Positive analytic wave integrals avoid subtracting nearly equal
incident/outgoing powers. Thin-layer midpoint parity retains the `d^3` normal
viscous loss near a rigid velocity node. Bounded inward-decaying waves and the
existing compensated log-amplitude owner handle thicker layers. Lossless fluid,
critical/evanescent fields and zero-thickness identities have exactly zero heat.

Reflection, transmission and local heat remain independent accounts with an
explicit floating-point balance residual; no residual is redistributed or
coefficient clipped. Queries use fixed storage, linear layer work, an explicit
caller allowance up to 64 layers and cancellation in both passes. Zero incident
watts retain coefficients and zero work. Existing `porous_layers` G0/G1/G3/G4
regressions compare interface fluxes and independent spatial quadrature, weak
transmitting loss, thin standing waves, subdivision, opaque/refraction limits,
source retention and atomic work/cancellation refusal.

For `ElasticTransmittingLayers`, the original longitudinal/SV boundary solution
supplies the layer states. Transmitted power is the sum of the two outgoing mode
powers; evanescent modes retain their loading but carry no normal power, and the
lossless elastic half-space produces no heat. `porous_elastic_layers` G0/G1/G3/G4
regressions compare an independent three-equation elastic boundary plus spatial
loss quadrature across both mode cutoffs, the matched thin-layer limit below
roundoff, normal-incidence fluid equivalence, subdivision, opaque layers and
cancelled/budget-refused queries.

`ElasticPlateLayers` likewise retains its original finite P/S plate solution,
including resonances and evanescent tunneling into the receiving fluid. Only the
declared entrance layers count toward report rows and the layer budget; the one
lossless plate adds bounded work and no heat. `porous_elastic_plate` regressions
compare independent six-wave boundary equations plus local spatial quadrature,
the exact zero-thickness and normal longitudinal-fluid limits, quarter/half-wave
changes in porous heating, sub-roundoff heat with evanescent plate transmission,
and cancellation/budget refusals.

`ElasticPlateStack::layer_absorption` extends this same positive work account
to porous layers before, between and after multiple fluid-coupled elastic
plates. Rows and the allowance count every original section, including lossless
plates and zero-thickness identities; those rows have exactly zero heat. The
original coherent P/S propagation retains each section's bounded face pair and
logarithmic scale. Existing analytic pore-viscous and thermal integrals then
resolve local heat without a second field solve, source query or spatial
quadrature. Solid reactance, evanescent tunneling and outgoing transmission are
never booked as heat. Allocation-free queries use fixed storage for at most 64
sections, with cancellation in both linear passes and before publication.
Regressions independently solve all twelve boundary conditions of a double-plate
system and integrate local viscous/thermal work through the gap, check normal
longitudinal limits before/between/after plates, retain heat below reflection
roundoff, shield downstream layers behind opaque plates, and exercise budgets,
identities, zero watts and atomic cancellation.

`ElasticPlateStack::prepare_angular_spectrum` also accepts the shared
`IncidentPlaneWave` and `AngularSpectrumBudget` interface for coherent
reflection/transmission. `PreparedElasticAcousticSpectrum::prepare_layer_field`
selects an original fluid/porous section before, between or after plates, with
indices counting every section. Original P/S propagation retains each plate's
signed transmission, bounded face states and logarithmic scale, preserving
phase for downstream interference. Fluid-only and elastic stacks share the
exterior field, local fluid propagation, coherent sum and positive heat owners.

Selected-layer queries return pressure, all three fluid velocities,
active/reactive intensity and separate local viscous/thermal heat after coherent
summation. Opposite-phase sources can cancel both losses; pressure nodes can
retain transverse viscous heat. Lossless fluid sections have zero heat and
opaque plates shield downstream pores. Plate indices refuse fluid-interior
queries. Zero-thickness fluid identities admit interface fields but no heat
density. Preparation uses the existing component/section caps and retains
O(components) storage per selected section. Queries allocate nothing, perform
O(components) work, re-query no source and poll cancellation. Original linear,
stationary, planar, rigid-frame and lossless-solid limits remain; no solid
interior field, intrinsic solid damping or temperature evolution is inferred.
Independent twelve-boundary-equation complex fields across both mode cutoffs,
spatially integrated versus analytic heat, normal-fluid equivalence and mixed
three-dimensional interference exercise this extension, together with nodes,
opaque shielding, identities and atomic cancellation.

The mechanism split does not assign mean heat separately to the pore fluid and
solid frame or evolve their temperatures. Original linear, rigid-frame,
infinite-planar and frequency/source restrictions remain; no additional
experimental validity or error certificate is supplied by the local account.

### Finite reservoirs in transient dense membrane transport

`membrane::transient::slab::finite` adds two explicit whole-reservoir molar
pressure capacities [mol/Pa] to the original sorbed-state finite-volume owner.
Both boundary pressures and all film cells advance in one tridiagonal
backward-Euler solve; no prescribed-pressure trial can spend a finite reservoir
before its depletion is solved. At connected equilibrium, total capacity is
`C_left + C_right + area*sum(S_i L_i)`, including the actual sorbed inventory.

The fixed-pressure and finite-pressure paths share normalized storage rows,
the original tridiagonal solver and face/storage reporting. Cached factors bind
the original immutable grid, both capacities and duration. Each reservoir face
transfer cancels exactly once against the existing membrane face owner; original
cell, reservoir and combined amount budgets gate every immutable candidate.
Exact seals retain separate conserved inventories. Empty numerical reservoirs
are allowed, while any physical gas adapter must enforce its own EOS domain.

Capacities and temperature are fixed external inputs. This numerical owner
supplies no EOS, caloric source, gas heat budget, or whole-system energy claim.
The same spatial/time-refinement obligations as the original transient owner
apply. Six additional G1/G3 regressions cover independent three-inventory
transients, first-order time refinement, heterogeneous capacity equilibrium,
large-reservoir limits, disconnected seals, capacity-cache identity and
transactional cancellation; all fifteen original/new numerical tests passed
against the actual coefficient/slab/finite sources and original tridiagonal
solver in the focused native lane.

### Material-source calibration of porous inertial resistance

`ForchheimerCalibration::from_material` selects the actual constant beta [1/m]
from an existing resolved material bundle; `PorousGasLayer::with_material_forchheimer`
selects it directly from the same bundle as intrinsic k. The complete original
card, selected claim, observations, uncertainty, query coordinates and quantity
metadata remain attached through an immutable shared source binding. This adds
no Ergun fit, dimensionless-to-dimensional conversion or inferred fluid identity.

Only exact inverse-metre, nonnegative validity-wide constant beta is admitted.
Its original source temperature/absolute-pressure domain intersects the existing
calibration support. An explicitly dimensional `superficial_speed` [m/s] axis
can tighten the maximum speed consumed by the existing gas, frozen-graph and
liquid owners. Pore velocity, unowned evolving coordinates and a source speed
range excluding zero refuse. Direction and transverse constraint still require
explicit calibration applicability; wrong specimen card, query context or
quantity metadata cannot be attached to the original permeability layer.

The existing scalar-authored API, anisotropic geometry, exact Darcy limit and
wet/slip exclusions remain in force. Observation-backed synthetic regressions
check original source retention, actual beta-dependent gas/continuum pressure
loss, wrong quantity/curve/domain refusal and specimen/context/direction
substitution. Three exact original source-domain tests passed in the focused
native lane; complete material-card integration execution is not implied by
those narrower domain checks.

## Primary creep with retained strain hardening

`NortonCreep::with_strain_hardening` explicitly selects
`StrainHardeningCreep { reference_strain, exponent, residual_rate_fraction }`.
The native accumulated equivalent creep strain supplies the history coordinate:
`h(alpha) = f + (1-f) (1+alpha/a)^(-m)`, with `a > 0`, `m >= 0` and
`0 <= f <= 1`. The original Norton reference rate is the **initial** rate at
zero accumulated creep. The residual fraction specifies the asymptotic
secondary/initial rate ratio at the same stress and temperature; zero gives
a regularized primary power law. No parameters are inferred from yield,
melting temperature or material names. Absent hardening, exponent zero and
unit residual fraction retain the original numerical Norton path exactly.

Tensor, free-sided axial and simultaneous isotropic J2 returns solve the
endpoint history equation together with stress balance. A bounded scalar
solve in `log(K delta_alpha / q)` retains small creep increments even when
their stress decrement rounds away, as well as very small relaxed stresses.
The exact algorithmic tangent includes the derivative of the implicit
hardening multiplier. The coupled yield-branch decision uses the accumulated
creep at the yield crossing; perfect plasticity still solves the creep
history equation at fixed yield stress. Effective stiffness remains logarithmic
through very small positive plastic hardening.

The same signed creep strain, nondecreasing accumulated path, trial/commit
separation, tensor shear convention and small-strain bounds remain authoritative.
Rate resistance is a kinetic memory, with no inferred recoverable hardening
energy: physical creep heat remains `q delta_alpha`, separately from plastic
heat, elastic/plastic storage and backward-Euler numerical loss.
Frame interval rebinding preserves the optional primary parameters.

`tests/creep_hardening.rs` supplies independent quadratic returns,
constant-stress continuum refinement, tensor covariance, nonproportional
consistent tangents, simultaneous yielding, exact disabled reductions, and
sub-stress-roundoff history/very-small-hardening regressions. This is an
explicitly authored isotropic small-strain model with first-order time
integration. It adds no cyclic recovery, tertiary damage, rupture, finite-strain
primary creep or experimentally qualified service-life prediction.

### Original source claims and native consumers

`NortonCreepSource::new_with_strain_hardening` and
`NortonJ2Source::new_with_strain_hardening` select the extended family explicitly.
They require three additional dimensionless original claims:
`creep_hardening_reference_strain`, `creep_hardening_exponent` and
`creep_residual_rate_fraction`. Their complete admitted temperature lines must
be supported and constant, including interior knots; all ten Norton claims
or twelve joint creep/plastic claims stay pinned on refresh. The original
constructors remain secondary-only. A card containing extra properties does
not silently change the selected constitutive model.

The existing `CreepLaw` and `CreepSourceFamily` paths carry the same history
and kinetics through tensor and free-sided axial materials, thermal expansion,
frame fibers, tetrahedral solids and solid/enthalpy coupling. The frame
continuation check refuses a changed primary selection. Source temperature
refresh changes Arrhenius kinetics without resetting the accumulated path.
Coupled thermal trials still begin from the same accepted mechanical history,
and only physical creep/plastic dissipation reaches the heat owner.

Native tests exercise sustained-load frame deceleration, permanent unloading
and reloading, simultaneous frame yielding with a residual secondary floor,
tetrahedral constant-stress time refinement and source-bound self-heating.
The self-heating reference couples the Arrhenius and primary laws to finite
heat capacity. These are manufactured analytical/consistency cases, not new
experimental calibration data or service-life qualification.

## Finite primary creep with shared kinetic history

`FiniteNortonCreep::with_strain_hardening` and
`FiniteNortonJ2::with_strain_hardening` select the same normalized
`StrainHardeningCreep` multiplier as the additive owner. The finite stress
measure is explicitly equivalent Kirchhoff stress; a Cauchy-stress fit is not
silently converted. Reference viscosity describes the virgin rate. Standalone
Norton retains accumulated creep in its original `alpha` field; simultaneous
Norton/J2 retains it in `creep_alpha`, with plastic `alpha` still controlling
energetic linear or Voce hardening. No new history schema or clock is created.

The original bounded logarithmic primary root is shared across kinematics.
The active finite joint return eliminates the plastic increment as
`p = (q_trial-q(c))/(3 mu)-c`, leaving one monotone scalar solve, including
perfect, linear and Voce plastic hardening. The original finite owner retains
isochoric metric reconstruction, objectivity, full analytic `dP/dF`, history
bounds and representability refusals. Absent or disabled primary selection
keeps the original secondary numerical path.

The reduced kinetic potential is the integral of the inverse endpoint rate
`q(c) = s0 [eta c/(dt s0 h(alpha_old+c))]^(1/n)`. Physical heat remains
`q delta_c`, with separate plastic yield heat and stored hardening. The common
`n=m=1`, zero-floor case uses its exact cubic primitive. Other selections use
an allocation-free transformed Gauss 7 / Kronrod 15 integral with a scaled
absolute/relative error estimate, explicit sharp-endpoint detection, a depth
cap of 28 and a 512-panel work cap. This is a bounded numerical estimate,
not interval certification. Unresolved quadrature returns an error.
`FiniteJ2Response::flow_potential_j_m3` retains the native integrated value
without subtracting large stored energies. Thermal tetrahedra and hexahedra
consume this field directly; plane-stress sheets already use native potentials.

`AffineThermalFiniteJ2::with_norton_strain_hardening` preserves the selection
through every physical-duration/temperature freeze and rejects Perzyna or
rate-independent selections. Kinetic hardening adds no fixed-history free
energy, entropy or heat capacity. Exact source selectors
`finite-primary-norton-hencky-kirchhoff-thermal`,
`finite-primary-norton-j2-hencky-kirchhoff-thermal`, and
`finite-primary-norton-j2-voce-hencky-kirchhoff-thermal` require all original
coefficients plus dimensionless `creep_hardening_reference_strain`,
`creep_hardening_exponent` and `creep_residual_rate_fraction`. They retain
20/25/28 parameters respectively and the original standalone/combined state
versions. Old selectors continue rejecting these extra parameters.

Focused cases cover finite quadratic/cubic references, nonproportional
potential and tangent derivatives, separate plastic/creep histories, primary
resistance through heating/cooling, P1/Q1 plane-stress continuation, loaded
thermal tetrahedra and source-bound adiabatic enthalpy feedback. These
manufactured cases do not supply experimental calibration, cyclic recovery,
tertiary damage, rupture or service-life qualification. Backward Euler remains
first order in physical time.

## Simultaneous permanent creep and isotropic yielding

`creep::plastic::NortonJ2` combines the original Norton kinetics and native
energetic linear-isotropic J2 model at one endpoint stress. Separate original
`CreepState` and `J2State` histories retain permanent creep, plastic strain and
both accumulated paths. A monotone logarithmic branch test selects creep-only
or active plastic flow. The active branch reduces to the original bounded
Norton scalar root with the effective elastic/hardening stiffness; perfect
plasticity has a closed endpoint. The tensor predictor and consistent tangent
remain shared with the Norton owner. No sequential full-duration return maps
or separately added elastic numerical losses are used.

Stored energy includes the native elastic and energetic hardening terms.
`CreepEnergy` exposes creep and plastic physical dissipation separately and
retains their total; numerical elastic/hardening loss is never thermal input.
Scalar balance, yield complementarity, reconstructed stress, separate history
admission and the complete work identity gate each checked candidate. A changed
yield offset may induce plasticity at fixed strain without changing storage.
Nonzero kinematic hardening is refused because its shifted flow need not be
coaxial with creep. Elastic and hardening storage coefficients stay fixed.

`NortonJ2Source` selects all nine original kinetic/elastic/yield/hardening
claims together. Continuation pins every selected identity and the original
source context. Complete temperature lines require constant elastic, kinetic
and hardening coefficients; yield may vary over its positive supported curve.
The original J2 owner performs the explicitly selected H or uniaxial-tangent
conversion. Immutable prepared source tables are shared across clones; the
composed resolver refreshes the full bundle once per endpoint. Missing data,
unsupported temperature support, evolving mechanical source axes and changing
hardening storage refuse. Source binding is not measured material calibration.

`CreepLaw` and the default-generic `CreepIncrement` carry this complete history
and work account to existing solvers. `CreepSourceFamily` supplies the same
original Norton/caloric admission to fixed-reference thermal consumers.
`tests/creep_plastic.rs` checks simultaneous analytical solutions, consistent
tangents, rotation/reversal, energy, exact activation and refusal;
`tests/creep_plastic_source.rs` checks heated yielding, source continuation,
hardening-basis conversion, whole-path admission and cancellation. This remains
small-strain secondary creep with energetic isotropic hardening; it adds no
finite-deformation, primary/tertiary creep, damage, rupture or service-life claim.

## Source-backed reversible nonlinear dielectric storage

`field_storage::dielectric::reversible::ReversibleDielectricMaterial` compiles
the explicit absolute `reversible_electric_displacement` D(|E|) source [C/m2]
on a caller-named electric-field axis [V/m]. The existing typed line-query
owner pins one original claim and retains every canonical curve or fixed-context
surface breakpoint, with its existing 4096-knot bound. The immutable card,
zero-field material snapshot, receipts, slopes and direct integrals are shared
across clones. Construction requires D(0)=0, strictly positive representable
segment slopes and inverse slopes, continuous support, and typed fixed external
conditions for every other query/source/applicability coordinate. Relative
permittivity, polarization alone, source switching, extrapolation and unowned
frequency, strain or history conditions are not alternative ingress paths.

`evaluate` takes a nonnegative field magnitude; `evaluate_displacement` and
`field_at_displacement` invert the same source at a nonnegative displacement
magnitude. The allocation-free O(log knots) response supplies D/E and dD/dE,
the transverse and longitudinal eigenvalues of the isotropic electric tangent.
At zero the first slope is used; interior knots use their right derivative and
the final knot its left derivative. Physical charge-form energy is
`w(D)=integral_0^D E(s) ds`; voltage-form coenergy is
`w*(E)=integral_0^E D(s) ds`. Each is integrated directly, preserving small
positive energy that subtraction from E.D would lose. Their sum is the Legendre
pairing E.D, not a license to replace either potential by E.D/2.

Signed scalar `energy_change_density_j_m3` uses even w(D), integrating only the
interval between displacement magnitudes. `energy_change_remainder_density_j_m3`
integrates `E(s)-E(reference_D)` over signed displacement with the original odd
inverse law. `coenergy_change_density_j_m3` integrates D dE between nonnegative
field magnitudes. These allocation-free, compensated crossed-segment operations
retain small circuit/mechanical work and line-search changes beneath large
stored biases; the remainder is not physical heat. Work is bounded by the
original path, and nonfinite, out-of-range or unrepresentable operations refuse.

The owner is a reversible, isotropic, zero-remanence fixed-context snapshot.
It supplies no leakage, hysteresis, polarization memory, thermal entropy,
electrostriction or moving geometry. A static temperature-dependent source
slice does not establish a closed conservative internal-energy model; such
consumers must either retain a source-independent potential or separately own
its thermal corrections. Original source axes remain inspectable for this
admission. `tests/dielectric_reversible.rs` checks independent piecewise energy
and coenergy, canonical source/inverse agreement, analytic tangents and the
low-field limit, signed crossed-knot work, large-bias remainders, surface-slice
receipts, and source/domain refusals. These authored fixtures establish numeric
model behavior, not measured dielectric accuracy or floating-point enclosures.

### Finite original gas reservoirs around a transient membrane

`transient::reservoir::IsothermalGasMembrane` couples the original source-bound
dense laminate to two original `ClosedGasParcel` reservoirs. Both gases must
contain exactly the membrane's single permeant at the exact held membrane
temperature. The existing gas catalog join checks caloric, transport, water and
source-pack identity; matching species names alone do not establish source
compatibility. The rigid ideal-gas volumes supply whole-reservoir molar pressure
capacities `V/(R T)` to the same finite-reservoir tridiagonal owner. Reservoir
pressures and all sorbed cell inventories therefore advance together.

Each accepted end-face transfer enters the original gas species-boundary owner
in its source-relative enthalpy coordinate. The explicit gas-side thermostat
heat is `-dm*(R*T/M)`, using the original owner's exact flow-work grouping.
Formation-reference constants cancel before the sensible-energy reconstruction.
The original gas owner reconstructs EOS pressure, transport, acoustics and
source-qualified state; the adapter refuses a changed held temperature. Stream
enthalpy, external thermostat heat and original gas energy residuals are
reported separately.

The proposed gases' EOS pressures must agree with the joint diffusion solution
within the declared equivalent-mass budgets. Only reservoir pressure coordinates
are then reseeded from those actual gas states; the same membrane inventory is
moved without duplicating its cell vector. Before returning an immutable
candidate, the adapter audits actual gas-plus-sorbed mass changes and actual
before/after total mass, and checks pressure-capacity/mass consistency. Original
material and ordered-interface claims are requalified across initial, numerical
final and reconstructed final cell/face pressure ranges. Failure or cancellation
publishes neither gas nor membrane changes.

This rung is pure-species, fixed-volume, fixed-temperature, linear Henry/Fick
transport. It reuses the existing factor cache and original gas/source owners.
The numerical reservoir owner permits empty reservoirs, but this adapter retains
the original gas owner's positive-density and positive-pressure domain rather
than inventing trace gas. D and S do not specify sorbed caloric storage or heat
of sorption: gas-side energy checks are not a whole-membrane energy certificate.
Mixed-gas competition, adiabatic evolution, swelling and nonlinear sorption need
their own constitutive and energy closures.

### Nonideal solution diffusion relative to finite stationary solvent

`osmosis::diffusion::SolutionMaxwellStefan` binds one original
`VirialSolutionModel`, explicit solvent and thermodynamic-order solute identities
and molar masses, and a complete positive symmetric matrix of constant binary
Maxwell–Stefan diffusivities [m²/s]. The original virial owner supplies
reference-cancelled chemical-potential differences. No ideal-gas EOS, inferred
tracer diffusivity, concentration floor or material authority is introduced.
Every endpoint must retain the same original model, fixed temperature, solvent
convention and nonnegative solute concentrations within its declared convex domain.

The admitted model neglects solute partial volumes and holds solvent stationary:
`c0=1/Vbar0`; logarithmic-mean face concentrations define x and total c. The
same symmetric mean supplies source sampling and the shared friction matrix.
The discrete ideal chain rule `c_log Delta(ln c)=Delta(c)` forms the solute
forces directly as `d_i=[Delta(c_i)+c_log,i sum_j B_ij Delta(c_j)]/(c L)`.
This equals `x_i Delta(mu_i)/(RT L)` at positive endpoints and keeps the
finite limiting force at an exact absent endpoint without a concentration floor.
The solvent reaction
`d0=-sum_solutes d_i` closes the constrained force system. The existing shared
molar Maxwell–Stefan owner performs the sole friction solve;
`N_i=J_i-(x_i/x0)J0` gives the published solvent-frame flux and exactly `N0=0`.
Every original friction row is checked again after conversion. Original
pair-friction entropy and `-sum_solutes N_i grad(mu_i)/T` must be representable,
agree within explicit budgets and remain nonnegative for regular endpoints.
A nonzero concentration jump into an exactly absent species has a singular
instantaneous endpoint diagnostic, reported as `None` for both entropy fields.
An identically absent spectator contributes exactly zero and does not disable
the regular entropy checks. These optional diagnostics do not waive the
evolution owner's finite-step free-energy and species-conservation checks.
Solute fluxes retain the caller's original order even when the shared solver
canonicalizes species differently.
This default is neutral fixed-solvent diffusion, not unrestricted liquid mechanics,
electrochemical migration, Soret/Dufour transport or an absolute caloric balance.

Successful face evaluations reuse the prepared original workspace and allocate
nothing; failures and cancellation preserve outputs. Source-owned fixed-volume
free-energy differences cancel the common linear concentration-reference term
before arithmetic and are limited to closed species-conserving accounting.
Independent public tests exercise binary solvent-frame reduction, ternary uphill
diffusion against separate two-unknown elimination, positive dissipation, species
permutation and concentration-reference invariance, adjacent representable
concentrations, pure-solvent uptake, absent spectators, disjoint multicomponent
fronts, bounded means across extreme finite scales, and atomic cancellation.
The binary ideal flux has the dilute Fick limit `-D Delta(c)/L` for arbitrary
concentration ratios. Its explicit finite-solvent factor `1/x0` remains; this
is not exact integration of a finite-concentration steady profile. The
logarithmic-mean chain rule follows the finite-volume construction of Cancès,
Ehrlacher and Monasse (https://doi.org/10.1093/imanum/drad032); no theorem for
their full time discretization is claimed for this constitutive face.

### Floating and prescribed-potential ionic Maxwell--Stefan diffusion

`osmosis::diffusion::electrical::ZeroCurrentDiffusion` assigns explicit original
solute valences and absolute inventory-charge/current tolerances. Both the
constant and material-source owners expose `flux_zero_current_into`. Every
finite endpoint must already satisfy its charge budget; no ion, dissociation,
activity model, concentration floor or compensating inventory is invented.
The original virial free energy and pair data must be qualified for those ions.

Charged face concentrations use ONE convex interpolation of the neutral
endpoints, weighted by reciprocal total ionic concentration along their straight
path. This retains neutrality and the logarithmic mean for proportional salt
profiles. Componentwise logarithmic ionic means would generally violate
neutrality for independent salts. Neutral co-solutes keep the original
logarithmic mean; both source resolution and friction call the SAME face helper.
Ionic force uses the actual face concentration times the original chemical
potential jump. A one-sided-zero charged species refuses its singular finite
face potential; exact species absent at both endpoints remain absent.

The native MS owner evaluates chemical and electric responses in the original
stationary-solvent frame. Eliminating the scalar electric gradient enforces
`sum z_i N_i=0`; a final combined-force solve rechecks all original rows, actual
charge flux and conjugate chemical/friction entropy. Electrical power vanishes
under this same current constraint, so it is not added as spurious heat. At most
three original friction solves use the retained workspace per face evaluation.
Reports return the right-minus-left floating voltage, actual current residual
and native solve count. No absolute voltage datum enters. A carrier-free equal
state has an undetermined voltage, explicitly `None`. Publication is atomic.

Both owners also expose `current_response_into`, returning an immutable
`IonicDiffusionResponse` with the chemical species flux and affine charge law
`q(dpsi)=q_chemical-G*dpsi`, where `dpsi=F(phi_right-phi_left)/(RT)` and q is in
mol/(m² s) of elementary-charge equivalents. Chemical and normalized electric
responses use the SAME full pair matrix. G must be positive when carriers
exist; an exactly absent ionic pool has exactly zero q and G. A response
retains the original neutral face, chemical forces, declared solve tolerances
and ionic policy. Constant data borrow the original matrix without allocation;
source data own the actual n-by-n resolved matrix and original left/face/right
receipts. Three bounded face vectors are retained, not a maximum-sized binary
matrix or cloned endpoint solutions. Source receipt/query allocations retain
the existing source path's explicit payload scope.

A network can assemble its global potential solve from these affine responses,
then call each response's `flux_at_potential_into` for ONE native combined-force
replay. That replay uses the immutable original source snapshot without another
material callback and checks every solvent-frame friction row and actual affine
current before publication. Preparation uses two native solves, or one for an
exactly absent carrier pool. Direct `flux_at_potential_into` on either owner
also accepts a prescribed finite voltage in one native friction solve; the
material owner first qualifies the actual sources as usual.

Finite current is permitted in the prescribed-field path. Chemical dissipation
alone may be negative, representing electric work driving uphill transport.
`PrescribedPotentialDiffusionReport` separates chemical entropy from signed
electrical contribution `-R*(dpsi/L)*sum_i z_i N_i`; their sum must be
nonnegative and match original pair-friction dissipation within explicit
budgets. The report returns actual signed charge flux, dimensionless voltage,
volts, combined electrochemical entropy and solve count. Exactly singular
neutral concentration fronts keep the existing absent chemical/combined entropy
diagnostic; the finite electrical term remains explicit. Nonfinite fields,
unsupported charged fronts, source failures and late cancellation preserve
caller outputs.

This is an isothermal, local electroneutral limit; it supplies no Poisson/Debye
layer dynamics, electrode reactions, ion association or caloric closure.
Consumers must enforce charge continuity and globally compatible potentials;
independent zero-current constraints on arbitrary cyclic graphs are not
justified. The
force/friction and current/entropy structure follows Van-Brunt, Farrell and
Monroe (https://arxiv.org/abs/2211.13505), equations 1.1--1.2, 6.5 and 16.13;
the dilute ambipolar limit and its scale assumptions follow Ghosal and Chen
(https://doi.org/10.1098/rspa.2010.0028). Public network-module regressions compare
1:1 salt flux/voltage, divalent nonideal salt with a neutral co-solute, and
independent Na/K/Cl salt fluxes against separate reduced equations and a
quadrature-constructed neutral face. Direct material regressions additionally
compare finite-current flux and conductance against independent two-row
elimination, admit electrically driven negative chemical entropy while checking
total dissipation, recover the floating-field limit, preserve orientation and
absent carriers, retain the actual source matrix across replay, and check late
cancellation and invalid-field output atomicity.

### Finite gas–adsorbent sensible heat on the original pore timeline

`adsorption::competitive::network::step_candidate_with_heat` and
`AdsorbingPoreDomain::{step_candidate_with_heat,advance_with_heat}` add explicit
sensible heat between each finite gas and its bound finite adsorbent. A pinned
nonnegative W/K material source owns exactly typed `gas_temperature`,
`solid_temperature`, and actual gas `pressure`. Unowned composition or coverage,
legacy aliases, dynamic response, extrapolation and claim substitution refuse.
Both actual source endpoints retain their original claim receipts.

Frozen-start `Q=G(Tg−Ts)Δt` enters the original gas inversion through a separate
heat record and the same dry-plus-bound inversion with equal opposite sign.
Transported partial enthalpy and adsorption binding energy remain separate.
Closed finite-system entropy, original absolute species/mass/energy budgets and
transactional cancellation remain mandatory. The existing bounded pore timeline
retries typed entropy decrease and accumulates actual directional heat over the
complete window; accepted temperatures feed subsequent adsorption and porous
transport. This supplies first-order feedback, not unconditional stability,
monotonic relaxation or a timestep-error certificate.

The public `adsorbent_gas_heat` integration target checks independent gas Cv and
dry-plus-bound storage in both directions, actual source refresh and ordering,
analytic two-body relaxation under temporal refinement, later adsorption and
porous-flow feedback, bounded entropy retries and atomic source/cancel refusals.
Its manufactured caloric fixtures establish numerical behavior, not empirical
conductance or adsorption accuracy.


## Source-bound thermal Voce hardening

`ThermoplasticSource::with_voce_hardening` and `CaloricPoint::with_voce_hardening` explicitly select source-pinned energetic Voce saturation on the existing fully solid caloric/expansion owner. Q(T) >= 0 must be globally affine and b > 0 fixed over the complete source line, including internal knots and surface crossings. Existing affine E/H_iso/H_k and native J2State remain; endpoint stress, consistent tangent, retained history and physical/numerical work delegate to VoceJ2.

With g(alpha) = alpha - (1 - exp(-b alpha))/b, additional free energy is Q g, entropy is -Q' g, and frozen-history coefficient work is (Q1 - Q0) g(alpha0). These enter the existing reversible heat and delta(T s) accounts. Neither saturation storage nor backward-Euler loss becomes physical heat. An exactly zero-Q endpoint retains its thermal entropy through a stable unit-Q shape; unselected and identically zero-Q paths preserve original arithmetic. Affine Q and fixed b add no frozen-history heat-capacity derivative.

Kinked/negative Q, varying b and tangent-only inference refuse. Feedback proposals retain the same old history and caloric chart; only balanced endpoint candidates publish. Focused regressions cover the full tangent, loaded heating, coefficient work and entropy, simultaneous thermal expansion/moduli changes, source support/pins, zero-Q limits, and retry/cancellation boundaries. These are small-strain, affine-coefficient verification results, not empirical thermal calibration or a Chaboche thermal-history extension.

## Stable vector work for reversible dielectric response

`ReversibleDielectricMaterial::coenergy_change_vector_density_j_m3` retains
physical work from component changes even when rounded field magnitudes agree.
It combines the original crossed-segment radial integral with compensated,
permittivity-weighted component products. A small transverse field change under
a large bias therefore retains its quadratic work, including in the linear
material limit. This is the source-owned work operation used by spatial
nonlinear backtracking; endpoint stored energies remain separate diagnostics.


## Independently sourced Chaboche backstress branches

`plastic::chaboche::source::ResolvedChabocheBranches<N>` binds one through eight independently named C_i/gamma_i pairs to the existing native ChabocheJ2<N>. `BackstressProperties` order is branch/history order, including explicit zero-C branches. There is no sorting, averaging, padding or inference from total stiffness or monotonic tangent. Each pair has exact Pa/dimensionless quantities, finite nonnegative values, fixed external conditions and independent claim selection.

Duplicate or aliased branch/backbone/reference-density roles, missing or ambiguous pins, evolving mechanical/rate/frequency coordinates, unsupported counts and total-strain tangent inference refuse. The complete original ResolvedVoceJ2 bundle, extra density/physics roles and source-use receipts survive adoption; ordered getters expose each original branch claim. SmallStrainLaw delegates stress, state, analytic tangent and symmetry declaration to the unchanged native constitutive owner, so existing generic continuum paths can use the resolved law directly.

Synthetic cyclic/source regressions cover two and three distinct recovery scales through reversal, explicit inactive/linear branch limits, equal total modulus with differing cyclic response, full source identity and receipts, complete quantity/condition admission, and the eight-branch bound. Canonical N=1 agrees exactly with the original single-branch source API. These tests verify integration rather than empirical calibration. Re-resolution at another temperature does not supply thermal Chaboche history transport or its energy account.

## Conserved-volume pure-fluid phase-changing energy closure

`thermochemical_gas::nonideal::equilibrium` recovers temperature, pressure and
phase amount from conserved density and original-reference internal energy.
Every subcritical trial resolves the original PR coexistence: volumes inside
the dome use the lever rule; outside volumes use the admitted liquid or vapor
density EOS. Unresolved coexistence refuses instead of silently retaining a
mechanically stable metastable state. The caller supplies a temperature interval
and pressure bracket that supports saturation throughout that interval.

The safeguarded solve uses the actual equilibrium isochoric tangent, physical
energy residual and temperature correction/bracket tests. Temperature and EOS
budgets cover the entire call, including nested phase searches and the initial
reference state; cancellation is polled inside those searches and before
publication. Increment solves integrate the original NASA energy change before
combining PR departures, so even a large formation datum is counted once and
does not erase small heat input. Source domains remain authoritative; no
two-phase transport, interfacial or reaction closure is inferred. The base
caloric state selects no acoustic relaxation policy; the explicit view below
admits the homogeneous-equilibrium limit.
Absolute-energy inversion refuses if adjacent represented energies cannot
resolve the requested energy/temperature tolerances or the entire caloric span
collapses; the source-anchored increment API remains available for large datums.

`equilibrium::parcel::IsochoricEquilibriumParcel` conserves fixed mass and
volume through signed heating/cooling and atomically commits the complete
equilibrium state only after independently evaluated source increments close
the joule balance. Failure or cancellation preserves the previous state.
Public integration targets `nonideal_equilibrium` and
`nonideal_equilibrium_parcel` exercise independent Maxwell-area/caloric targets,
heating and cooling across both liquid and vapor dome boundaries, large-datum
invariance, retained density, total-work limits and atomic refusal. Synthetic
NASA fixtures check constitutive consistency rather than experimental accuracy.


### Simultaneous creep and yielding in free-sided axial fibers

`creep::plastic::uniaxial::UniaxialNortonJ2` condenses the original joint
Norton/J2 law by eliminating transverse stress. The same prepared logarithmic
scalar return map uses E for an axial fiber and 3G for a three-dimensional
point. It solves creep and isotropic plastic flow together, including perfect
plasticity, and supplies the exact condensed tangent. It retains separate signed
creep/plastic strains and separate accumulated equivalent paths through unloading
and reversal. Each history reconstructs its original traceless tensor, including
lateral inelastic contraction; the transverse elastic strain uses the original
Poisson ratio.

The optional source-bound law retains all original kinetic, elastic, yield and
hardening claims through a shared immutable source bundle. Checked energy delegates
to the native full-tensor joint law and independently compares stress, free-side
conditions and both reconstructed histories. Hardening remains stored energy;
creep heat, plastic heat and numerical loss remain distinct. The generic energy
carrier also reserves explicit reversible thermal-storage and transfer accounts,
both exactly zero for these unexpanded laws.

The `creep_plastic_uniaxial` regressions compare linear closed-form solutions,
nonlinear tangents, the native three-dimensional law, zero-stress unloading,
reversal, original source selection, and the yield/creep activation boundary.
This remains small-strain, isotropic secondary creep with isotropic J2 hardening.
It does not imply finite-strain validity, thermal expansion, kinematic hardening,
primary/tertiary creep, damage, rupture or measured service-life qualification.

Validation for this increment: 22 actual material tests passed (8 joint axial,
8 three-dimensional joint, 6 original Norton axial) using the pinned Cargo
runner and genuine cached dependencies. This focused exact-source material
assembly is supplementary evidence, not the complete workspace/DSR gate.

### Wetting-liquid retention and relative permeability

`porous::wet::liquid` supplies an explicitly calibrated van Genuchten--Mualem
wetting-liquid branch. Residual liquid and gas saturations, shape m, nonnegative
connectivity exponent, pressure scale in Pa, saturation/temperature/absolute
liquid-pressure/superficial-speed support and provenance are explicit.
`WetPoreLiquidLaw::sample` evaluates effective saturation, liquid relative
permeability and Pc=pg-pl; `at_capillary_pressure` inverts the same retention
branch. Stable logarithmic/expm1 expressions retain small positive dry/wet
increments. The exact dry endpoint reports zero liquid mobility and
`LiquidCapillaryPressure::DryLimit`, never a fabricated finite pressure.

`ContinuumDarcyGeometry::with_liquid_states` binds one retained liquid state
per original ordered layer. It compiles sum(L/(k_effective*kr_l)) while retaining
the original intrinsic tensors, frames, transverse-flow constraint and source
bundles. Rebinding always starts from original k and cannot compound old
relative permeabilities. Original and liquid calibration T/p/speed limits are
intersected; an intrinsic seal or zero liquid mobility stays exactly sealed.
Positive Forchheimer inertia refuses this view, including a unit-mobility
endpoint, because a saturated inertial calibration does not supply a wet one.
Gas-slip and gas-only saturation calibrations retain their existing refusals.

The constitutive and Darcy rebinding interfaces prescribe instantaneous liquid
mobility on one fixed branch. They do not solve gas pressure, liquid storage,
a wetting front, Richards/two-phase dynamics, capillary work or hysteresis. Capillary pressure is
a phase-pressure relation and must not be added as another hydraulic drag or
heat source. Synthetic analytic and source-preservation fixtures check the
equations and coupling; they do not establish a measured material calibration.

### Implicit conservative unsaturated-liquid columns

`porous::wet::liquid::transient::LiquidSaturationColumn::step_candidate_implicit`
advances the existing rigid one-dimensional finite-volume column by backward
Euler. Each cell retains its original retention calibration, porosity,
thickness and explicit intrinsic normal permeability. Prescribed gas pressure,
temperature, viscosity, boundary liquid pressures and signed body-force
pressure gradient define one physical interval. The explicit and implicit
owners share the same source-qualified Darcy face flux. Half-cell resistances
add in series, including unequal widths and permeabilities; an exact dry or
intrinsic seal has exactly zero connected flux and never acquires an epsilon
mobility or a fabricated finite dry pressure.

Analytic capillary-pressure and mobility derivatives assemble an O(n)
tridiagonal Newton direction using the native adjacent-pivot solver. If the
full mobility tangent loses diffusion signs, a frozen-mobility pressure
direction supplies the nonlinear correction. Every trial still samples the
complete original endpoint law and flux. Bounded backtracking admits only
original saturation, actual liquid-pressure, temperature and superficial-speed
support; it neither clips states nor widens source domains. Unsupported initial
coordinates, exhausted supported iteration, entirely unsupported trial steps
and cancellation remain distinct refusals.

Acceptance independently checks each dimensionless storage/flux equation and
one fixed original-to-final column liquid-volume budget per transverse area.
Compensated actual cell changes avoid subtracting large absolute inventories;
rejected trial fluxes never enter the accepted report. Counts include all
evaluated and rejected states, linear solves and safeguarded directions.
Failure preserves the original saturation, material states and clock. The
endpoint method has no explicit stability timestep restriction, but still
requires supported convergent states and representable time advancement; it
does not certify temporal accuracy, introduce hidden substeps, solve gas or
heat storage, infer fluid density, or add hysteresis or capillary work.

`advance_candidate_implicit` and its boundary-conditioned counterpart advance
one complete requested interval with bounded adaptive substeps. A coarse
backward-Euler endpoint and two fine endpoints estimate local saturation error;
only the two original fine endpoints and their signed face transfers are
accepted. Absolute tolerance and relative tolerance on the actual saturation
change control refinement independently of the nonlinear residual tolerances.
The comparison accounts for unequal representable half intervals. Minimum and
maximum comparison steps and an explicit implicit-solve budget bound work.
Unsupported original source coordinates remain fatal; supported nonlinear
failure or excess temporal error can shorten a trial. Rejected probes never
contribute to the compensated integrated face and net boundary transfers.
Acceptance recomputes the original-to-final water-volume balance against one
fixed whole-window budget, without multiplying it by the accepted step count.
No partial state escapes a refusal or cancellation. The local first-order error
estimate is not a certificate of global, spatial or constitutive accuracy.
Independent nonlinear ODE refinement, closed unequal-cell equilibration and
prescribed-flux inventory fixtures exercise the public whole-window path.

`LiquidSaturationBoundary` supplies sealed, signed superficial-flux and
pressure/resistive-pressure conditions to that same explicit/implicit flux
owner. Positive flux points left to right on both ends. Surface hydraulic
resistance [Pa s/m] acts in series with the actual porous half-cell resistance,
with no additional storage or body-force span. Its exact pressure and mobility
sensitivities enter the existing tridiagonal nonlinear solve. The original
pressure APIs select zero surface resistance.

Initial and accepted states qualify actual speed and reconstructed porous-side
boundary pressure. An imposed flux cannot force an unsupported hydraulic head;
external pressure behind a resistive surface is distinct from material-side
pressure. Nonzero flux through exact dry/intrinsic seals refuses. The typed
explicit entrypoint requires a fixed volume budget; implicit local residual and
whole-column budgets retain their original meaning. Public fixtures check
inventory-determined closed equilibration, infiltration/drainage and throughflow,
manufactured nonlinear Robin balance and its mirror, actual source domains,
budget refusal and late cancellation. No finite-reservoir compliance, coupled
gas/heat storage or empirical hydraulic calibration is inferred.

### Passive dissolved transport in evolving porous liquid

`LiquidSaturationCandidate::transport_dilute_solutes` carries original dilute
neutral-solute inventories over the SAME accepted explicit or implicit water
interval. The candidate retains exact original saturation/time binding and
checks the original cell/calibration identities. Both original and final liquid
storage must be positive; drying deposition or reopening requires another
phase/storage owner. Every cell's water-volume change and the whole-column
balance must close against the original hydraulic face integrals.

One backward-Euler donor-cell/diffusion matrix per solute uses actual final
water storage and `theta*D/tau`, with theta=porosity*saturation. D is explicit
aqueous diffusivity; tau is the transport divisor in this formula. Geometric
or acoustic tortuosity conventions require an explicit source-supported
conversion. Native half-cell resistances and explicit surface film resistance add
in series. A no-diffusion boundary can prescribe hydraulic inflow concentration;
a reservoir can drive molecular exchange at zero water flow. Zero D blocks
diffusion without suppressing the original hydraulic transport.

The existing native adjacent-pivot tridiagonal solver gives O(cells*species)
work and bounded reusable scratch. Every original cell equation and one fixed
constituent balance is checked independently after solving; signed advective
and diffusive face transfers are separate. No concentration clipping, fictitious
water flux or explicit diffusion CFL is introduced. This is a first-order
passive dilute approximation with caller-qualified D/tau and applicability,
without inferred activity, cross-friction, ionic migration, reaction, adsorption
or hydrodynamic dispersion. Local time error requires physical interval refinement.

Eight tests of the exact production transport operator and original native
tridiagonal solver pass: independent two-volume BE and eigenmode limits, unequal
wet resistance, changing-storage uniform concentration and reversal, film-driven
exchange, stiff advection, positive trace response, budgets and cancellation.
Two additional native hydraulic binding/interval integration tests are authored
but unexecuted; the full material dependency build remains resource-blocked.

### Reversible solute sorption in porous liquid

`LiquidSaturationCandidate::transport_sorbing_dilute_solutes` couples the same
accepted hydraulic interval to an explicit linear immobile pool per cell and
species. The declared capacity B [m] is dry-solid mass per area times Kd [m³/kg],
so equilibrium loading is b=B*c. Equilibrium mode requires an already consistent
original loading within the absolute amount tolerance; it does not project the
initial state. Finite-rate mode admits disequilibrium and evolves
`db/dt=k*(B*c-b)` with a separately declared nonnegative rate [1/s]. Zero rate
isolates the original pool, while zero B still permits release of an old pool.

Backward Euler eliminates each bound variable analytically. With
`alpha=k*dt/(1+k*dt)`, only `alpha*B` and `alpha*b0` enter the concentration-system
storage diagonal and RHS. The original water volume, hydraulic balance,
advective concentration and `theta*D/tau` half-cell resistances remain unchanged.
The same native tridiagonal solver needs one O(cells) solve per species. Bound
output and report payload are admitted before allocation. Cancellation publishes
no partial hydraulic, aqueous or sorbed state.

Aqueous and sorbed endpoint amounts and actual net bound exchange remain separate.
The new report's `total` account explicitly counts both pools; its signed
boundary advection/diffusion integrals contain aqueous transport only. Original
cell/column total-amount balances and represented bound-exchange defects have
independent checks. No substep adsorption/desorption history is inferred. This
is a first-order passive dilute linear model with caller-qualified fixed B/k;
competitive or saturating adsorption, ionic effects, reactions, sorption heat,
drying deposition and hydraulic feedback require other constitutive owners.

Sixteen focused native tests pass using the exact transport/sorption modules and
original fs-la tridiagonal source: eight existing aqueous cases plus independent
two-pool kinetics, retarded diffusion, stiff/equilibrium and zero-rate limits,
changing water, time refinement, trace release, budgets and cancellation. One
hydraulic adapter/clock-binding regression is authored but unexecuted because
the full material dependency graph remains resource-blocked.

### Source-bound finite-capacity sorption during porous solute transport

`LiquidSaturationCandidate::transport_langmuir_dilute_solutes` couples the
original aqueous finite-volume owner to `AqueousLangmuir` equilibrium storage
or its independently sourced `LocalLangmuirTransfer` LDF law. Retained dry
adsorbent mass per bulk transverse area converts original mol/kg loading into
mol/m². The same analytical element is required throughout each constituent;
independent slots do not compete for sites. Inert slots retain zero site stocks.
Actual hydraulic water capacity times an explicit physical transverse area
supplies the local source's solution-volume coordinate. Prescribed temperature,
additional chemical conditions, concentration and loading retain original typed
source admission at the initial state and every nonlinear trial/end point.

Backward Euler eliminates local bound loading as a convex combination of its
original stock and the original isotherm target. The occupied and vacant stocks
remain independent so a representable vacancy tail survives rounded saturation;
vacant sites never enter the transported analytical-element balance. Both site
balances use cancellation-sensitive compensated sums. Equilibrium initialization
must already agree with the original isotherm within the declared amount budget.
Finite-rate mode permits disequilibrium, adsorption and desorption; zero local
rates and inert columns retain the original free-transport solve arithmetic.

Damped, concentration-box-preserving Newton reuses native adjacent-pivot
tridiagonal LU; only the local storage tangent changes. Water geometry, accepted
hydraulic face integrals, donor-cell advection and wet diffusion/film resistances
remain unchanged. The original source owners expose shared borrowed admission
and stable loading/tangent seams, avoiding source-model clones and duplicate
isotherm evaluation in nonlinear trials. Every trial consumes the cumulative
source-point budget; linear solves, iterations and backtracks are independently
bounded. Source refusal, loss of representability, cancellation and unclosed
equations discard the whole private proposal. The numerical workspace allowance
includes all Newton/phase arrays; original source query/receipt allocations and
borrowed source carriers remain outside that numerical-payload allowance.

Nine new authored G1/G3/G4 groups passed natively with the exact production
numerical and original source-owner modules: independent kinetic quadratics,
heterogeneous diffusion face roots, saturating breakthrough/outlet integrals,
zero-rate/inert parity, changing-water uniformity, first-order time refinement,
independent vacancy tails, source/site refusals, work limits and cancellation.
All 16 existing porous transport/sorption tests and 27 existing aqueous
isotherm/batch/flow tests also passed (52 total); compiler stderr was empty.
This bounded native run does not type-check the full hydraulic adapter or
replace the resource-blocked fs-material/workspace gate. Changed Rust files
passed the pinned formatter check. This remains first-order passive dilute
analytical-element transport: no competitive uptake, ionic migration, reaction,
sorption heat, evolving site capacity or hydraulic feedback is inferred.
Newton convergence and residual budgets do not certify global convergence,
temporal error, or experimental validity of the caller's material sources.

### Multiple source-bound sorbent populations during porous solute transport

`transport_multi_langmuir_dilute_solutes` extends both the original hydraulic
candidate and its prepared dispersion interval with separately identified
independent sorbent populations. Each retains its original `AqueousLangmuir`
or `LocalLangmuirTransfer`, dry mass per physical area, occupied stock and
independent vacant stock. Different finite LDF rates retain distinct uptake and
desorption histories even when aggregate bound inventories coincide. These are
disjoint site populations, not competing species on a shared capacity.

The one-population owner and this extension share the same source-admitted
local response. Each population's analytic backward-Euler storage change and
nonnegative tangent are summed with compensated arithmetic into ONE original
damped Newton/native tridiagonal aqueous solve. Darcy advection and molecular
or source-qualified mechanical dispersion keep their original physical face
laws. The accepted report retains each population in caller order as well as
aggregate bound and combined aqueous/bound balances; continuation uses the
separate occupied/vacant stocks, without repartitioning an aggregate.

Sparse cell/constituent bindings are grouped once with CSR offsets. Up to the
original batch owner's 128 populations may occupy a slot, under an additional
explicit total-population cap. Newton trials allocate no population arrays.
A shared source-evaluation cap counts initial, trial and accepted points;
linear-solve and numerical-payload caps cover the complete solve, including
CSR, retained per-population outputs and prepared dispersion storage. Every
source checks its actual concentration/loading, physical cell solution volume
and prescribed temperature. Source boxes must intersect, all populations in a
physical cell must share one liquid temperature across constituents, and
equilibrium populations require consistent original loading. Refusals identify
the original population index. Inputs remain unchanged on refusal/cancellation.

Eight focused native groups plus the existing 60 passed (68 total), using
exact production owners and genuine unchanged dependencies. Independent
two-rate closed-cell quadratics check opposing uptake/desorption and the stiff
limit; frozen/slow limits and one-population arithmetic are checked separately.
Equal-total/different-partition states retain distinct desorption memory.
A heterogeneous advection/dispersion face oracle, sparse constituent ordering,
original multi-population `AdsorptionBatch` parity at the same two BE half steps,
source/temperature refusals and combined resource limits also pass. This
bounded proof does not type-check the full hydraulic adapter or full workspace,
which remain resource-blocked. The model is first-order passive transport with
fixed independent site populations; it adds no shared-site competition,
sorbent motion, thermal feedback, evolving capacity or global time-error bound.

### Source-qualified longitudinal mechanical dispersion in porous liquid

`LiquidSaturationCandidate::with_solute_dispersion` prepares one exact hydraulic
interval for the original free, linear-sorbing or source-Langmuir solute owner.
A distinct original constitutive card supplies nonnegative longitudinal
dispersivity alpha_L [m], with typed speed, saturation, porosity and physical
transport-scale support. Calibration scale is independent of numerical cell
width. Additional original conditions remain explicit caller-owned context;
consistency with hydraulic temperature/chemistry not retained in the candidate
is not inferred. No material name, molecular diffusivity or tortuosity supplies
a hidden dispersivity.

Each adjacent cell uses its actual endpoint water fraction theta and the same
accepted face volume divided by interval duration: v_face=|q_face|/theta.
The bulk half-cell coefficient is theta*D_molecular/tau + alpha_L*|q_face|.
Molecular and mechanical contributions add before original half-cell series
and boundary-film assembly. Mechanical dispersion is never divided by tau.
Opposing faces remain separate, so converging/diverging water does not erase
physical mixing through a signed cell-average cancellation. The original
molecular-only calls retain their original resistance/arithmetic/allocation
path; zero-alpha and stagnant prepared intervals retain molecular flux and
inventory arithmetic exactly.

The returned nested `diffusive_face_amount_mol_m2` is explicitly the total
gradient-driven molecular-plus-mechanical transfer. No separately conserved
mechanism split is invented at heterogeneous interfaces. Advection and every
aqueous/bound conservation gate remain original. Preparation retains both
face speeds, both mechanical coefficients and actual cell coordinates with
borrowed original source bindings. It consumes the existing source-point and
memory budgets before allocation; its retained arrays also count inside the
subsequent whole transport memory cap. The consuming interval moves its arrays
into the result without a deep clone; source qualification clones no cards or
query points. Failure/cancellation publishes no partial constituent state.

Eight new authored G1/G3/G4 groups passed with exact production numerical/source
owners and genuine native dependencies: zero-alpha/stagnation parity, opposed
faces and tortuosity independence, heterogeneous series/film/reflection oracles,
a three-level advecting Gaussian, linear and Langmuir combined balances,
source/unit/scale refusals, and cumulative payload/work/cancellation limits.
The Gaussian checks both the continuous profile's first-order refinement and
the independent discrete variance, separating physical dispersion from donor-cell
and backward-Euler broadening. All prior 52 porous/aqueous-owner tests also
passed (60 total), with empty compiler stderr. Four changed Rust files passed
the pinned formatter check. This bounded run does not type-check the full
hydraulic adapter or replace the blocked fs-material/workspace gate.

This is one-dimensional longitudinal dispersion for passive dilute constituents
on the candidate's time-frozen water interval, using endpoint storage. Explicit
and implicit hydraulic candidates retain their own chosen flux time level.
Unresolved within-interval changes/reversals require smaller water intervals.
No transverse/tensor dispersion, turbulent mixing, evolving calibration,
charged-species migration, experimental accuracy or global error bound follows.

### Adaptive coupled gas, adsorption and sensible-heat transport

`AdsorbingPoreDomain::{step_candidate_adaptive_to,advance_adaptive_to}` control
the complete adsorbing-pore transient by comparing one original native interval
with two original half intervals. Componentwise actual gas-species masses,
bound-species moles, gas temperatures and solid temperatures each use explicit
absolute tolerances plus an optional relative allowance on changes. Large
background inventories and absolute temperatures cannot hide a small transient.

Only the fine path can publish or contribute directional gas-species,
adsorption, skeleton-heat or gas/solid-heat transfers, original caloric defects
or entropy. Every probe is a single existing unsplit candidate. Its original
source, material, geometry and physical refusal owners remain active; source
and domain errors remain fatal. Existing pore work limits count all native
probes and accepted physical half steps. Representable midpoint/final clocks,
late cancellation and no-prefix publication are mandatory.

Complete original-to-final gas-plus-bound species inventories and compensated
accepted energy/entropy must fit the same single physical budgets, independent
of subdivision. Reports retain the last native source state and accepted gross
directional transfers. The first-order local step-doubling estimate is
asymptotic; it is neither a global error certificate nor a replacement
thermodynamic model. The public heat integration tests check independent
two-body thermal decay, tighter error convergence, independent gas/bound/thermal
accuracy controls, coupled porous transport, exact replay, and atomic work,
source, entropy and cancellation refusals.

### Material-sourced, composition-dependent liquid diffusion

`osmosis::diffusion::source::MaterialSolutionMaxwellStefan` binds every original
solvent/solute and solute/solute Maxwell–Stefan pair to an explicit immutable
material card and exact property key. Claims are selected once at an admitted
reference solution and remain pinned. Each evaluation resolves actual left,
logarithmic-mean concentration face and right states against those same claims;
endpoint support cannot be inferred from an admissible midpoint. Curves and
joint sampled surfaces retain the native resolver's interpolation support,
units, evidence and refusal behavior. The face matrix enters the existing
solvent-frame solver through one shared internal path.

The source adapter owns only the original thermodynamic model's absolute
temperature and explicit amount-concentration coordinates, including solvent
concentration from its molar volume. Unowned pressure, missing pairs, ambiguous
or incompatible quantity declarations and unsupported physical states refuse.
Binary kinetics are never inferred from virial free energy or substituted with
tracer coefficients. Original thermodynamic species order and canonical solver
order are mapped explicitly. Zero concentrations require actual original source
support and are never replaced by a seeded trace. This remains the original
isothermal incompressible-solvent description, without pressure evolution.

Successful face evaluations return complete left/face/right source receipts
and the current symmetric binary matrix. Source queries allocate transient
material evidence; the constant-coefficient API retains its allocation-free
advancement. Prepared source storage retains original cards and claim IDs,
without a dynamic receipt cache. Its estimated payload charge includes all
retained claims, observations, models, source geometry and indexes, including
unselected data, and may conservatively count a shared card per pair. It is
not an allocator-overhead, transient-memory or resident-memory bound.

Public source regressions check an independent concentration-dependent
binary flux, complete original receipts, canonical permutations, pinned
endpoint and joint-surface refusals, required kinetics/coordinates and retained
unselected source data, plus source-qualified pure-solvent uptake and explicit
refusal when zero lies outside the pinned source domain. Liquid-core regressions
cover subnormal positive concentration ratios: the chemical log ratio forms its
dimensionless relative change before the atanh argument, preserving dilute
forces without concentration floors or invented species.

### Finite multicomponent membrane separation

`transient::mixture::IsothermalGasMembraneMixture` couples two original finite
ideal-gas mixtures to independently sourced linear Henry/Fick channels on one
shared rigid laminate. Every joined catalog species needs an explicit channel,
with the same layer thicknesses, grid and held temperature. A source-selected
D=0 channel blocks transport and retains its authored sorbed inventory; missing
species data never imply impermeability.

Each channel reuses the original conservative finite-reservoir tridiagonal solve
with C=V/(R T). All species transfers enter one original source-relative gas
boundary reconstruction per endpoint. The complete original catalog remains
available when a species is absent from one gas. Fluxes are relative to the
stationary membrane and driven by partial pressure, with no gas-diffusion
zero-sum constraint. Changing composition updates the original gas EOS,
thermochemistry, transport and acoustics together.

Original layer and ordered-interface source support is checked against SUMMED
gas-equivalent pressure at cells and both local face traces, including original
EOS total pressures at exterior boundaries. Old, solved and actual reconstructed
states are admitted; interior totals may exceed both reservoir totals during
selective counterpermeation. Per-species and combined actual gas-plus-sorbed
mass budgets are checked after EOS pressure reseeding. Canonical species order
and each channel's existing model/capacity/duration factor cache are retained.
Failure or cancellation returns no partial gas, membrane or clock update.

Independent sorption is an explicit model assumption, and original claims must
already qualify the intended mixture. Selected composition, humidity or
concentration applicability restrictions refuse rather than staying silently
fixed as composition evolves. This does not promote pure-gas measurements,
infer competitive sorption, cross-diffusion, swelling or sorbed calorics.
Original gas enthalpy streams and exact isothermal thermostat flow-work offsets
are separate reports, without a whole-membrane energy certificate.

`IsothermalGasMembraneMixture::advance_candidate` completes one requested window
with the original native implicit steps and separate full/half-step factor
caches. One full step estimates the temporal error of two accepted half steps;
every sorbed-cell Henry pressure, gas species partial pressure and actual gas
species amount receives its own comparison. The relative allowance scales the
change from the last accepted state, so large equilibrium inventories or
opposing species streams cannot hide separation error. This local first-order
estimate is not a global temporal error certificate. The original complete
catalog also admits pure gas through the same one-species API.

Accepted half steps alone contribute signed and gross directional species
streams, original gas enthalpy streams and explicit thermostat heat. Rejected
probes spend the explicit work budget and contribute no physical account.
Actual original-to-final gas inventories and cellwise sorbed changes must obey
fixed per-species boundary, film and closed-system amount budgets, plus both
direct-total and stable-increment mass budgets, independent of substep count.
Whole-window gas energy compares original-source sensible storage with the
accepted source-relative input, cancelling partial-enthalpy reference transport
before arithmetic; the separately bounded inventory residual is not multiplied
by an arbitrary formation datum. Reported absolute enthalpies retain that datum.
No sorbed caloric store or heat of sorption is inferred from D and S.

Total retained profile cells and all original complete-mixture probes are
bounded. Representable shared clocks, explicit step limits and cancellation
remain admission gates. Source, caloric, linear and gas reconstruction refusals
are hard failures rather than guessed timestep-retry signals. Any late failure
returns no partial gas, film or clock candidate. Original material selection
and actual cell/face/bath source-domain checks remain native at every probe;
the final candidate retains current source support and the final native step
receipt. Source queries and immutable native candidates still allocate; the
prepared numerical factor caches do not imply an allocation-free advance.

`new_with_surfaces` admits independently authored exterior molar interface
permeance [mol/(m2 s Pa)] through the existing `InterfaceSystemCard`. Left order
is gas to polymer; right is polymer to gas. Finite-gas attachment checks the
exterior against the original species transport material. Surface and adjacent
half-cell resistances share the original scaled series calculation and native
tridiagonal owner. `None` selects perfect Henry contact; zero surface permeance
seals exactly. No extra gas-film volume, sorbed storage or caloric law is inferred.

Gas and polymer traces remain distinct. Original polymer D/S support follows its
actual side; original interface support covers both sides. Mixture polymer traces
sum species Henry pressures while exterior gas traces use actual EOS totals.
Coincident surface and bulk seals refuse source qualification when their
intervening trace is undetermined. Original source selections, fixed species/mass
budgets, isothermal thermostat accounts and transactional adaptive advancement
remain native. Pure-gas evidence does not acquire mixed-gas kinetic authority.
Public surface fixtures check an independent series-resistance steady flux and
finite relaxation mode, selective mixture evolution, actual two-sided source
domains, exact seals, ordered identities and cancellation.

### Original-material liquid retention calibration

`WetPoreLiquidLaw::from_material` selects residual liquid/gas pore-volume
fractions, van Genuchten m, the connectivity exponent and the capillary pressure
scale from an original `ResolvedMaterialStatePoint`. The first four parameters
require exact dimensional-one quantities; the pressure scale requires Pa.
These are validity-wide constants. Sampled curves, missing parameters, pressure
head and inverse-length/pressure alpha coefficients cannot silently supply this
calibration. The complete original claims, observations, uncertainty, selection
and query context remain accessible with the selected property names.

Source applicability intersects the declared saturation, temperature, absolute
LIQUID pressure and superficial-speed support. The explicit `liquid_saturation`
coordinate denotes pore-volume fraction; relative humidity and dry/wet mass
moisture bases are not convertible substitutes. Superficial speed uses m/s and
a domain beginning at zero. Unsupported evolving axes and pressure-dimensional
aliases such as capillary pressure refuse. The constitutive output Pc=pg-pl
remains separate from absolute liquid-pressure applicability.

A source-backed liquid law can attach only to intrinsic permeability retained
from the same immutable material card, complete original query and quantity
metadata. Separate property selections may share that context. New prescribed
saturation samples retain their original calibration context while obeying its
supported interval. Authored calibrations remain available as explicitly
authored inputs, with no measured authority inferred from their provenance.
The shared continuum and native pressure/heat receipts retain this source
through the exact liquid states used to compile mobility.


### Sourced thermal expansion of permanent creep and joint yielding

`creep::thermal::expansion::ExpandingCreepSource` composes an original Norton
or joint Norton/J2 source with `AffineIsotropicThermoelastic`. It requires an
explicit expansion claim, stress-free reference, original card and fixed source
context on the canonical absolute T coordinate. The native thermoelastic owner
integrates affine instantaneous expansion and supplies entropy and a conservative
full-interval mechanical heat-capacity bound. E, nu and energetic hardening remain
constant through the original creep source admission; Arrhenius kinetics and
sourced yield still follow the actual temperature.

`ThermalCreepState` retains the accepted temperature with every original creep
and plastic memory. Constitutive trials subtract sourced free strain from total
strain and call the original return map and tangent. A checked increment resolves
both actual thermal endpoints and delegates the mechanical energy audit. Explicit
tensor-Voigt/Mandel conversion preserves shear work. The potential receives total
strain minus all original irreversible strains, while each component history and
the total/elastic strain remain independently bounded.

The energy carrier keeps free-energy storage, physical creep/plastic dissipation,
backward-Euler numerical loss, change of T times mechanical entropy and signed
reciprocal caloric transfer separate. With f the integrated expansion, reciprocal
heat is trace(sigma_new)*delta(f)-delta(T*s). Its sum with physical dissipation is
`caloric_heat_j_m3()`; this sum can be negative during thermoelastic cooling even
though physical dissipation is nonnegative. Every accepted increment closes total
stress work against all these independently evaluated terms. Numerical loss is
never deposited as material heat.

The source exposes a conservative mechanical capacity floor over the complete
thermal/history envelope and checks a virgin zero-displacement birth against its
actual stress-free reference. A caloric consumer must add its own original chart
capacity and admit a positive combined capacity. `creep_expansion` contains five
manufactured regressions for free and restrained affine expansion, mixed-shear
tangents/work, simultaneous thermally driven creep/yield, retained temperature
through cooling and reversal, capacity bounds and original-source refusals.
This is a bounded small-strain model with fixed storage coefficients, reference
mass and geometry. Variable elastic/hardening storage, non-affine thermal source
paths, phase change, finite deformation and measured life qualification remain
outside this composition.

Validation for this increment: all 56 focused creep/material regressions,
including the five expansion cases, compiled, linked and executed successfully
with the pinned compiler and genuine dependency libraries. This was an exact-source
material-owner assembly, not a Cargo workspace acceptance gate.

### Variable-volume nonideal solution free-energy changes

`VirialSolutionModel::reference_cancelled_inventory_free_energy_change_j`
evaluates the same original mixing potential between finite solute inventories
whose solvent volumes may also change. It admits both endpoint concentration
states and removes only the common linear concentration-reference term before
arithmetic. The caller must close actual species inventories independently;
an open compartment still owes its chemical-reference transfer.

Relative amount and volume logarithms preserve dilution even when the rounded
endpoint concentrations coincide. Symmetric bilinear differences retain the
virial excess contribution without subtracting absolute background energies.
Exact absent-species limits remain continuous, with no concentration floor.
The bounded O(N²) evaluation uses stack scratch and supplies no caloric energy
or mechanical work. The coupled finite solution-network owner composes this
change with its existing pressure law's exact potential primitive and checks
one complete-window mixing-plus-mechanical account.

Four public inventory-energy regressions check independent finite potential
differences, concentration-reference invariance, the fixed-volume limit,
next-representable-volume dilution and its osmotic-pressure derivative,
chemical-potential derivatives, zero-species limits and original support.

### Joint native and externally sampled gas transport

`species::network::transport::step_candidate_with_paths_and_fluxes` joins
integrated species and additional heat records to the same immutable snapshot
as all native `GasTransportPaths`. Their physical owner must qualify the
original constitutive sources, frame, geometry and interval and declare a
distinct exposure. Added heat excludes the original donor partial enthalpy
already carried by each species. Existing entry points delegate with empty
additional records.

The original combined work bound, gross original-species withdrawal, source
catalog join, one caloric reconstruction per endpoint, fixed absolute mass and
energy budgets and closed entropy gate apply jointly. No incoming species can
fund outgoing allowance. Cancellation or any refusal publishes no endpoint.
The seam supplies no new constitutive law or temporal controller.

The finite continuum adapter in `fs-gas::material::transport::parcel` consumes
this seam for the native mass-frame Maxwell–Stefan/Soret/Fourier/Dufour face.
Its additional heat is reduced heat times area/time plus each signed species
mass times the original analytic donor-to-face partial-enthalpy increment.
Formation references cancel before arithmetic, and donor enthalpy is charged
once. Actual initial and final participating gas states must also satisfy the
original transport sources. Continuum exposures remain distinct from dusty-gas
pore transport; this composition supplies no porous thermal-transpiration law.

`transport::audit_closed_window` closes a rigid gas network directly between
its original and final inventories. It reuses the original source-qualified
entropy traversal and, only when requested, accumulates analytic NASA internal
energy increments with one common reference removed per conserved species.
The identity is `n1*du(T0,T1) + (n1-n0)*du(Tref,T0)`; absent-species limits query
only temperatures at which that species actually exists. Original source cards
and explicit caloric-reconciliation metadata must match. Actual total and
individual species mass residuals, internal energy and entropy obey one fixed
budget independent of accepted-step count. Energy/entropy datum roundoff is not
multiplied into the separate mass allowance. The existing closed `PoreDomain`
window uses this same audit instead of subtracting absolute mixture energies.
The existing molecular-diffusion adaptive window also checks actual endpoint
energy independently of its compensated accepted-step residual ledger.
`DiffusionNetworkError::EntropyDecrease` distinguishes finite-step entropy
refusal from hard source/domain errors. Public gas integration tests exercise
large formation references, sub-resolution absolute-energy defects, species
appearance/disappearance, irreversible mixing, source/cancellation refusal and
porous counterdiffusion. This is a closed endpoint account, not a temporal error
estimate or an audit for externally heated, reacting or moving-volume systems.

`transport::audit_reservoir_window` supplies the separate rigid, nonreacting
held-bath chart: Delta U minus bath species partial enthalpy and reduced heat,
and Delta S minus bath species partial entropy and reduced heat divided by its
actual bath temperature. Each explicit bath model must reproduce its retained
resolved state. Original caloric, transport, pack and water-species bindings
remain identified even when a species vanishes at one parcel and appears at
another; a final endpoint cannot authorize a new source through a mass tolerance.

Analytic species increments cancel common formation and entropy references
before rounding. Pressure-flow work appears once as the ideal-gas R*T_bath term;
no absolute large-reference subtraction or mass-tolerance-scaled datum enters
the energy/entropy result. Source-qualified and aggregate per-species mass,
direct total mass, energy and total entropy retain fixed absolute budgets.
An absent finite species may enter from an authored bath; nonzero exchange at
zero bath activity refuses. This is an endpoint account for explicitly supplied
bath-chart streams, not a rate law, arbitrary open-stream chart, trajectory
certificate or finite-reservoir storage model. Public gas fixtures independently
check ideal mixing entropy, gas cooling, hidden energy defects, original source
identity across absent-species transitions, fixed budgets and late cancellation.

### Source-bound conductivity through latent phase change

`thermal::conduction::phase::PreparedPhaseConductivity` binds one selected
scalar static conductivity claim to an explicit enthalpy interval of the original
borrowed `EquilibriumEnthalpyPhaseCurve`. Pure-phase, latent and solidus/liquidus
paths use the canonical pinned source-path evaluator. Temperature and liquid
mass fraction can evolve together, or temperature can change at a fixed mixed
fraction. Every original crossed source cell/curve knot and caloric knot is
retained. The claim selected at the initial state is pinned throughout; this
does not select a new source at later states. Caloric and conductivity sources share the same material
card and fixed external conditions. An evolving fraction requires a declared
fraction axis, even for a phase-independent conductivity source.

The compiled law implements the native `ConductivityLaw` and
`ImplicitConductionLaw`, including nonzero `dk/dh` on a plateau where `dT/dh=0`.
Interior breakpoints use right derivatives and the admitted upper bound uses
left derivatives. The adapter also supplies the original curve's analytic
`d(1/rho)/dh` for fixed-mass implicit planar expansion/contraction, preserving
the same one-sided convention at a clipped phase boundary. Sensible intervals
differentiate the retained rho(T) through original Cp(T); latent intervals
retain their specific-volume law even when `dT/dh=0`. This volume tangent
changes transport resistance without changing inventory or adding pressure
work. The caloric owner's public numeric sample supplies unchanged
temperature, fraction and density without constructing a state identity. Runtime
evaluation allocates nothing and performs no source query or hashing. Original
source paths retain their receipts through `source_paths`; samples identify
their original path with `source_path_index`. Compiled samples are numerical responses,
not new measurements. A combined knot cap bounds preparation, with cancellation
around source queries and during copying. Typed source/caloric failures remain
available; native thermal traits retain their `MaterialDomain` refusal contract.

The existing explicit observer and implicit chain/graph own storage, heat flux,
time stepping and acceptance. No phase mixing rule, liquid fluid properties,
reference-mass update or second thermal solver is introduced. The same law's
piecewise-affine derivative includes both T and fraction dependence through
the original caloric segment. Negative conductivity and unsupported source
paths refuse. Explicit zero insulation is
admitted only with the zero tangent required by the native implicit owner.
Constitutive regressions cover source crossings, one-sided tangents, original
caloric values, source-only material swaps, insulation, source/domain refusals
and cancelled preparation. Cross-owner latent transport, finite hot-body
solidus/liquidus transport and contact tests live in `fs-couple`; independent
backward-Euler energy elimination supplies their finite-body heat oracle.
The compiled-law identity domain is v2 for the generalized original source paths.

## Joint Voce hardening and Perzyna overstress

`plastic::rate::VoceOverstressJ2` jointly integrates independently selected Voce saturation and linear or power-law associative overstress at an explicit physical duration. Its single consistency equation includes the elastic/linear-hardening correction, the exact incremental Voce strength and the inverse kinetic stress. Summed secants feed the original J2 tensor return; summed endpoint derivatives feed its consistent tangent. No intermediate rate-independent return commits or changes incoming history. The linear case reuses the native Voce solve; nonlinear kinetics use a safeguarded logarithmic bracket with a 96-iteration ceiling and no per-response heap allocation.

Physical dissipation is `(sigma_y + f_visc) delta_alpha`. Recoverable saturation energy and its convex backward-Euler remainder retain their original storage and numerical accounts. Q=0 delegates exactly to the original selected rate law; n=1 recovers linear kinetics. Explicit positive viscosity/duration/stress normalization, n>=1, finite native histories and the declared small-strain ceiling remain mandatory. Nonlinear Chaboche recovery, sub-yield creep, finite deformation and empirical calibration are outside this law.

`tests/plastic_voce_overstress.rs` verifies independently prescribed joint rates, all tangent components under reversal, cyclic work partitions, exact reductions, a closed-form constant-stress Voce creep solution with first-order refinement, tiny representable rate increments and invalid input refusal. Focused native execution does not qualify the whole Cargo workspace or the live cancellation runtime.

## Finite dissolved gas in cavitating liquids

`cavitation::dissolved_gas::HenryDissolution` exchanges one selected
noncondensable species between a finite liquid inventory and an original
`ClosedGasParcel`. Source-resolved molar Henry solubility (C*=H p_i) and
molecular diffusivity combine with authored film area/thickness. The exact
affine finite-pair relaxation includes changing gas partial pressure and uses
the native deterministic expm1 owner; it neither assumes an infinite reservoir
nor clips an explicit Euler update. Other gas species retain their inventories,
and exsolution can introduce an initially absent selected species into an
existing carrier without a numerical seed.

The original species-boundary/caloric owner reconstructs gas composition,
pressure, transport and acoustics. Original gas partial enthalpy and an
explicitly aligned dissolved partial-enthalpy chart determine separate gas
and solution thermostat heat; pressure-flow work is counted once. Gas and
liquid state/clock publish together only after source, mass, energy and final
cancellation checks. Both bulk and interfacial concentrations must fit the
original material support and authored dilute limit.

`DissolutionFilm::gas_side_transfer` adds an explicitly fixed, pressure-based
gas coefficient in mol/(m² s Pa), using the existing `MolarPermeance` quantity.
The original scalar and simultaneous owners use the series law
`1/K_l = delta_l/D_l + H/k_g`. The finite equilibrium, gas reconstruction and
thermal accounts are unchanged; only the physical rate and interface traces
change. `None` retains the original one-film arithmetic and `Some(0)` seals
the gas path. Two simultaneously sealed films refuse an undetermined massless
interface. Coefficient values carry no inferred material or experimental
authority, and must apply over the complete pressure/composition interval.

Bulk liquid and the ACTUAL liquid-side interface are queried against the
original material support and dilution bound. With gas resistance, `H*p_bulk`
is only a driving potential. A conducting interface cannot exceed the gas total
pressure, and the mixture path bounds the combined interfacial partial pressures
throughout the interval. Opposing transfers can require a shorter step to satisfy
that conservative envelope. The liquid timescale check uses the reduced actual
relaxation rate `lambda*delta_l²/D_l`. It does not certify the separate gas film.
The two-film derivation follows the interfacial equilibrium and flux continuity
in [ETH Zurich's film model](https://www.hyper-tvt.ethz.ch/fundamentals-masstransfer-filmodel.html).

This is a fixed-temperature, fixed-volume dilute-solution model with negligible
solute partial volume and already quasi-steady films. The authored
film-relaxation bound is an applicability check, not a
transient diffusion-layer solution or accuracy theorem. Nucleation, bubble
radius motion, evaporation, Stefan flow, composition-dependent gas transport and
adiabatic thermal evolution require separate
physical owners; no generic entropy-production claim is made. Four native
tests cover independent finite equilibrium/relaxation and thermostat energy,
carrier-preserving exsolution, small-step/partition identities, and atomic
source/dilution/film/cancellation refusals.

`dissolved_gas::mixture::advance_mixture` composes up to 64 independently
sourced dilute species into the same finite liquid and gas. The unchanged
scalar film solve prepares all transfers from one original gas snapshot; one
native species-boundary solve reconstructs the mixture. The liquid volume and
clock must agree exactly across species, and the volume is counted once.
Absorption and exsolution may occur simultaneously. The global energy account
includes the gas storage once, all dissolved partial-enthalpy changes, and both
explicit thermostat inputs. Species ordering is canonical, and all gas/liquid
histories publish together after final cancellation admission.

The summed solute concentrations must fit the declared dilute bound. Because
opposite species transfers can make total pressure nonmonotone, the complete
interval uses a conservative envelope of the monotone partial pressures for
both Henry and original gas-source admission. A source may therefore require
a shorter interval even when its final state fits. Six G0/G3 regressions cover
independent finite-inventory oracles, opposite transfers and heat, input-order
invariance, time subdivision, exact scalar reduction, combined dilution and
layout refusal, final cancellation/retry, and an unsupported interior pressure
excursion with supported endpoints. These are synthetic equation checks;
cross-solute interactions, reaction, moving bubbles and adiabatic temperature
feedback remain outside this fixed-temperature dilute model.

Additional scalar and mixture regressions check two-film finite-inventory
oracles, both interface fluxes, unchanged equilibrium and heat, actual source
support, sealed-film limits, reduced liquid relaxation rate, time subdivision,
shared-interface refusal, canonical ordering and transactional cancellation.


### Independently sourced gas-film resistance

`HenryDissolution::with_gas_film` optionally adds a separately selected
`GasFilmSource` for the solute's diffusivity in the declared carrier gas.
The gas coefficient `k_g = D_G/(R T delta_G)` enters the same prepared series
resistance as `DissolutionFilm::gas_side_transfer`; selecting both gas-film
descriptions refuses ambiguity. Independent thicknesses and diffusivities give
`1/K_L = delta_L/D_L + H R T delta_G/D_G`. The original exact finite-pair
relaxation, native gas reconstruction and thermostat accounts use that shared
conductance. The default liquid-only path retains its original coefficient
and interface admission.

Actual interfacial partial pressure and concentration satisfy Henry equilibrium
and the two films' common molar flux. Both endpoint traces retain these values,
the gas-film material answers and each film's quasi-steady applicability ratio.
Source support and authored dilute bounds apply to bulk and interface states.
Simultaneous exchange retains one original carrier inventory and conservatively
bounds the combined gas dilution and full-interval pressure support. Zero
diffusivity seals a single film; two sealed films refuse because the interface
state is undetermined.

Four synthetic native regressions check independent series resistance and
finite equilibrium, both film fluxes, absorption/exsolution heat, time
subdivision, liquid-only and sealed limits, source/timescale/cancellation
refusals, combined species dilution and permutation replay. Equivalent authored
and sourced coefficients share the same trajectory; dual selection refuses
without changing state. This extension assumes dilute Fick diffusion with
admitted constant diffusivities and negligible film storage at prescribed
temperature, fixed volumes and fixed area. Stefan flow, finite diffusion-layer
transients, nucleation, moving bubbles and adiabatic heat feedback remain
outside these isothermal entry points; a selected numerical fixture is not measured validation.

### Adiabatic finite-solvent absorption

`HenryDissolution::advance_adiabatic` couples the existing finite-inventory
Henry/two-film exchange to explicitly supplied solvent mass and original liquid
Cp(T). Gas and solution share one evolving temperature under a fast thermal
equilibration assumption. Original NASA-9 gas Cv and source-exact solvent/solute
caloric integrals close the combined energy balance, with no external thermostat
heat. Original Henry, diffusivity and optional gas-film sources continue at the
actual new temperature before further uptake. The thermal balance remains
anchored to the beginning of the requested interval to avoid accumulating
independent substep residuals.

The kinetic integration is first-order splitting with bounded full/two-half-step
refinement in both dissolved moles and temperature. Reported absolute error
estimates are not certified global errors or source uncertainty. Source, mass,
energy and final cancellation checks precede publication of law, gas and liquid
together. Any unsupported source or thermal-bracket trial refuses atomically.
Synthetic tests cover absorption heating, exsolution cooling, nonlinear gas
calorics, temperature-dependent solubility against an independent closed-energy
ODE, and thermal/work-budget/final-cancellation refusal.

Rigid volumes, dilute solution, quasi-steady films, negligible liquid expansion
and negligible dissolved partial volume remain explicit assumptions. Solvent
mass is authored independently of chart density. Finite thermal contact
resistance, solvent evaporation, phase change, reaction and simultaneous
multisolute adiabatic exchange are outside this entry point. Authors must align
caloric datums and qualify Henry/caloric laws together; no entropy-production
theorem is inferred.


## Thermal Voce with physical overstress kinetics

The original ThermoplasticSource and CaloricPoint now compose selected Voce saturation with linear or power-law Perzyna kinetics in either selection order. Source coefficients remain pinned over the full original caloric interval. Endpoint stress, state and consistent tangent use the joint VoceOverstressJ2 return; actual physical dt is mandatory. The shared thermal first-law account retains Q(T) storage, its coefficient work and entropy alongside positive yield/viscous dissipation. Kinetic coefficients add no recoverable energy or entropy. Numerical integration loss remains separate from heat.

Thermal regressions check independent flow, endpoint free energy and entropy during heated stress relaxation, all tangent components, preserved source receipts, duration refusal and exact zero-Q recovery of both original kinetic laws. The common thermal-account completion also admits separately owned internal-variable entropy and coefficient storage without duplicating expansion or caloric physics.


## Thermal Chaboche history and reversible branch energy

ThermoplasticSource::with_chaboche_branches binds one through eight ordered independently selected C/gamma source pairs over the original caloric span and explicit Voce backbone. Each active C(T) is globally affine and strictly positive; identically zero branches are admitted, and gamma remains fixed. Branch names cannot alias existing thermal or mechanical roles. Full source pins and branch order persist.

ThermalChaboche keeps z=3X/(2C) fixed under coefficient change, transports X*=C1 X0/C0 before the unchanged native Chaboche return, and reuses its generally nonsymmetric tangent, physical recovery and numerical accounts. Shared thermal energy adds branch entropy -C' z:z/3 and coefficient storage (C1-C0)z0:z0/3. Thermal rescaling and numerical loss are never physical heat. One resolved old-to-new law represents one accepted physical increment; Newton proposals reuse its frozen old history. Resolve a new law for the next physical increment. Retain accepted stress: another timed return evolves flow even at unchanged strain. With rate-independent coefficients only, a new-to-new law can inspect accepted state without repeating thermal transport.

Seven synthetic material regressions check complete source support and pins, reversible heating/cooling, all tangent components, combined affine elastic/Prager/Voce/branch coefficients and expansion, cyclic reversal, zero-branch reductions and invalid domains. Native tetrahedral component checks live in fs-solid. Changing gamma, active zero-crossing C, nonlinear-temperature branch free energy, finite deformation and empirical thermal-cyclic calibration are outside this scope.

## Prepared numeric gas caloric samples

`GasCaloricBranch::sample_density_internal_energy_increment` exposes
temperature, pressure, Cv and the original-law energy increment from the SAME
bounded density/Delta-u inversion used by full material reconstruction.
The original source/composition borrow, continuous caloric branch, positive
capacity bounds, residual checks and cancellation are retained. The existing
temperature solve is factored once; no second EOS or caloric root is supplied.

This numeric stage view avoids rebuilding transport/source receipts during
each stage of a coupled ODE trajectory. It is not a complete resolved material
state. Accepted physical endpoints must still use the original full
reconstruction, which resolves transport and all retained query axes at the
new temperature and pressure. Three native prepared-branch tests pass,
including exact numeric/full-state agreement for heating, cooling and zero
increment, source reuse, invalid domains and cancellation.


## Joint multi-branch Chaboche and Perzyna kinetics

plastic::chaboche::rate::ChabocheOverstressJ2<N> integrates the original Chaboche/Voce/Prager history with explicitly timed linear or power-law Perzyna flow. It solves native consistency F(delta_alpha) minus the inverse kinetic stress, retaining the increment-dependent nonproportional direction and every implicit backstress. The added kinetic derivative enters the original generally nonsymmetric tensor tangent. Linear viscosity reuses the native consistency solve with eta/dt added to its scalar correction; nonlinear kinetics use a safeguarded logarithmic bracket with a 96-iteration ceiling. Inverse kinetics are shared with VoceOverstressJ2. Responses use fixed branch arrays without heap allocation.

Physical dissipation contains yield work, viscous overstress work and native dynamic-recovery loss exactly once. Voce, linear and branch storage and all backward-Euler numerical losses retain their original owners. Absent branches reduce exactly to VoceOverstressJ2; exponent one recovers linear kinetics. Unselected native Chaboche and its axial restriction retain their original path. Positive physical duration, admitted kinetics and finite reachable history are required. Static recovery, sub-yield creep, finite deformation and empirical cyclic-rate calibration remain outside this law.

Seven focused regressions check independently prescribed increments and each branch, all 36 tangent entries under nonproportional reversal, complete cyclic energy partitions, exact reductions, analytical nonlinear stress relaxation with time refinement, tiny rate increments and invalid history/kinetics. Existing Chaboche, axial and branch-source regressions protect the unchanged untimed paths.


## Timed thermal Chaboche and Perzyna source integration

Thermal Chaboche supports explicitly selected linear and power-law Perzyna kinetics through the joint native Chaboche return. resolve_enthalpy_for_step and advance_with_energy_for_step require the actual positive physical interval; endpoint iterations retain the same original history and duration. Endpoint viscosity, stress scale and exponent resolve from the original pinned source bundle. Thermal transport of each independent backstress precedes the joint return. Yield, viscosity and dynamic recovery contribute physical dissipation once; viscosity adds no free energy or entropy. Duration-free rate calls refuse.

Three additional thermal regressions independently check nonproportional tangent/flow, complete two-temperature energy and entropy, positive viscous/recovery dissipation, source kinetic refresh, physical-duration relaxation, frozen-history retries and exact absent-branch/n=1 reductions. Static recovery, sub-yield creep, varying gamma, finite-deformation plasticity and experimental thermal-cyclic calibration remain outside this implementation.

## Explicit equilibrium acoustic view and conservative pressure Jacobian

`PureFluidEquilibriumState::equilibrium_acoustics` borrows the original admitted
state and explicitly selects instantaneous phase and thermal equilibrium with
common phase velocity. Its `HomogeneousEquilibriumAcoustics` view implements
the existing `FluidAcoustics` contract, so original acoustic consumers receive
coherent density, sound speed, bulk modulus, compressibility and impedance
without reconstructing a substitute EOS. Homogeneous liquid/vapor branches
reuse the original PR caloric mechanics. Strict two-phase interiors use the
original coexistence Clapeyron slope and moving-phase Cv with the source NASA
ideal Cv counted once. Exact or unresolved saturated endpoints refuse a
bidirectional acoustic law rather than silently choosing a one-sided response.

The same view exposes pressure and entropy derivatives in physical `(T,rho)`
coordinates, the isentropic temperature slope, and the conservative derivatives
`p_u|rho = p_T/Cv` and
`p_rho|u = p_rho|T + (p_T/Cv)(T p_T-p)/rho^2` for mass-specific internal energy.
They close `c^2 = p_rho|u + p p_u|rho/rho^2` without subtracting formation
energies. The original state, source receipt and admitted density remain
available through the borrowed view. Gas transport cards are not required;
no mixed Cp, finite-rate relaxation, phase slip, surface storage, attenuation
or frequency-validity prediction is supplied.

The public `nonideal_equilibrium_acoustics` target checks independent
Maxwell-area/isentropic methane EOS oracles, full-state pressure/entropy and
conserved-energy finite differences, lower equilibrium sound speed than the
existing frozen Wood closure for matching phases, original homogeneous
mechanics, endpoint refusal and bitwise invariance under a `1e20 J/kg` source
energy shift. The NASA caloric fixture is synthetic and checks constitutive
consistency, not experimental methane sound-speed accuracy. Generic acoustic
trait use retains the original thermodynamic state instead of flattening it
into an unrelated material model.

## Source-backed nonlinear dielectric temperature scaling

`field_storage::dielectric::reversible::thermal` explicitly admits
`D(E,T)=s(T) D_ref(E)` from the original reversible dielectric claim and a
separate original positive dimensionless affine scale claim with `s(T_ref)=1`.
`ThermalReversibleDielectricMaterial` retains the complete material bundle,
selected source identities, original field path, nonthermal coordinates and
supported absolute-temperature interval. The reference must have no temperature
value axis. This owner does not infer separability or thermodynamics from a
general `D(E,T)` surface. Original source receipts remain reference evidence;
resolved snapshots share the immutable prepared table and evaluate actual scaled
displacement, both tangents, inverse response, charge-form potential, coenergy
and scalar/vector work. `is_temperature_scaled()` makes this additional physical
role explicit to coupled owners.

For actual displacement magnitude `D`, `z=D/s` and volume `V`, the field terms
are `F=V*s*w_ref(z)`, `S=V*s'*g_ref(z)`, `U=F+T*S`, and
`C_D=-V*T*(s')^2*z^2/(s*epsilon_ref,diff)`. Here `w_ref` and `g_ref` are the original
charge-form potential and coenergy, respectively. The charge-form response is
Helmholtz free energy on a thermal snapshot, not its complete internal energy.
The zero-field caloric chart supplies the thermal datum; signed field entropy
and internal energy are reversible contributions. The full-strip capacity bound
includes both sides of original tangent jumps and conservative outward rounding.
A finite thermal owner must keep the combined caloric and field capacity positive.

`admit_source_refresh` permits only the same original source family and owned
thermal coordinate. Fixed-history parameter work remains the responsibility of
the native electrical owner. `thermal_offset_change_j` retains `Delta(T*S)` by
compensated original inverse-segment coenergy integration and explicit affine
scale increments, including work below endpoint spacing. It never labels that
exchange as irreversible dielectric heat. Unsupported source strips, malformed
scale laws, changed bundles, negative-magnitude state queries and unrepresentable
nonzero field contributions refuse. Focused evidence is in
`tests/dielectric_reversible_thermal.rs`; existing dielectric kernel and Debye
thermal tests continue to exercise the shared material boundaries.

### Source-bound multicomponent caloric and enthalpy-constrained flash

`thermochemical_gas::nonideal::flash` adds complete PT calorics and a bounded
fixed-pressure enthalpy inversion over the original Peng–Robinson two-phase
owner. The immutable result retains the original critical-property/mixing law,
canonical feed, full NASA receipt, phase roots and original material/fugacity/
Gibbs residuals. Phase mole, mass and volume fractions remain distinct; bulk
density is conserved feed mass over additive phase volume.

The exact feed carries each ideal species enthalpy/internal-energy reference
once. Amount-weighted original phase departures provide the nonideal part.
Entropy also includes the phase redistribution of ideal mixing:
`R [sum(z ln z) - (1-beta) sum(x ln x) - beta sum(y ln y)]`.
Thus a rounded phase material residual is never multiplied by a large species
formation datum, and feed ideal mixing is not incorrectly retained unchanged.
Each phase must have stable original mechanical and caloric response using its
own composition; a feed Cv is not substituted for separated-phase Cv.

The enthalpy solve uses safeguarded secants and periodic bisection, physical
enthalpy residual and temperature-step/bracket acceptance. Its caller-declared
temperature interval must retain original source support and explicit caloric
continuity. Both temperature endpoints must resolve a two-phase state when a
search is needed. At identical pressure and zero enthalpy increment, an admitted
reference within the declared bracket is returned exactly after one PT call.
The increment API integrates original NASA changes and then adds departure
differences, preserving small heating/cooling and pressure-change calorics under
large independent species reference shifts. It also returns the independently
evaluated entropy change. Absolute-energy spacing that cannot resolve the
requested tolerance refuses with an increment-API alternative.

The work limit counts all original PT calls, including the reference and both
bracket endpoints. The reported cubic-EOS count is explicitly an upper bound,
`PT_calls * (6 * maximum_iterations + 2)`, not a measured iteration count.
Cancellation is retained through nested phase searches, bounded source
admission and final publication. Single-phase or unresolved search endpoints
refuse without a stability fallback. This remains the original bounded
two-phase solution rather than a global Gibbs-minimum certificate; no mixed
Cp, acoustic relaxation, gas transport, interface kinetics or flow-rate law is
invented. A process consumer must impose its own entropy/work convention.

The public `nonideal_flash` target compares independent 50-digit simultaneous
cubic/fugacity PT and cubic/fugacity/enthalpy PH solutions, independently
reconstructs full phase NASA/mixing/departure calorics, and checks species/mass/
volume closure, isenthalpic pressure-drop entropy, source-datum invariance,
small heat cycles, bounded work, exact fixed points and cancellation. The
equimolar methane/ethane fixture uses synthetic Cp/R values 3.5 and 4.5; these
checks verify constitutive consistency, not experimental mixture accuracy.
The phase enthalpy convention is documented by the original-law cubic EOS
reference at
<https://idaes-pse.readthedocs.io/en/stable/explanations/components/property_package/general/eos/cubic.html>.

### Multicomponent enthalpy paths across screened phase boundaries

`thermochemical_gas::nonideal::flash::equilibrium` adds source-bound screened homogeneous/separated PT states and fixed-pressure absolute/increment enthalpy recovery. A temperature bracket may span the mixture's liquid, coexistence and vapor regions. Every sampled state uses the original PR two-search phase screen; a witnessed instability must produce the original admitted split. Homogeneous states additionally require the original complete mechanical, Cp and Cv admission. The retained screen is finite evidence, not a global-minimum certificate; unresolved searches or amount-indeterminate equal-Gibbs roots refuse.

The strict two-phase public APIs retain their existing contract. Both APIs share one bounded scalar enthalpy solver and the original two-phase caloric owner. Exact feed ideal enthalpy/entropy references occur once; phase redistribution supplies the original EOS departures and phase-resolved ideal-mixing entropy correction. Source-integrated enthalpy and entropy increments cancel species reference datums before rounding. The increment path can heat or cool across both bubble/dew boundaries. At unchanged pressure and zero increment, an admitted reference inside the bracket returns exactly after one PT selection.

The shared owner admits the complete original caloric source interval, demands a representable increasing enthalpy bracket, checks sampled monotonicity, and accepts only the specified enthalpy residual together with its temperature safeguard. Absolute targets whose floating-point spacing cannot resolve tolerance refuse with guidance to use increments. Publication retains the original source receipts and polls cancellation after all closure calculations. Work receipts distinguish observed PT selections from the sum of original screening work and conservative nested-flash EOS ceilings.

The existing `nonideal_flash` integration target checks independent 50-digit binary liquid/vapor and coexistence oracles; heating and cooling across both phase boundaries; the original EOS identity h-u=p/rho; exact fixed points and large species datum invariance; mandatory refusal of unresolved/negative-witness splits and pure-fluid amount ambiguity; total work, source support and final cancellation. The feature infers no mixed Cp, acoustic relaxation, transport, reactor, flow-rate or process entropy law.

### Original phase-chart identity in native enthalpy consumers

The `ThermalLaw` adapter for `EquilibriumEnthalpyPhaseCurve` exposes the original
curve identity through `caloric_identity`. It is the same identity already
binding the material card, enthalpy datum, phase knots and retained caloric source
receipts. Optical laws indexed by native enthalpy can therefore verify that their
declared caloric reference is the actual registered thermal owner, including on
isothermal latent plateaus. This adds no duplicate phase inversion or latent heat.
It does not authenticate caller-authored optical data or infer a phase mixing law.

The `fs-render` example `phase_change_optical_heating` uses this adapter directly
for illumination, infrared exchange and phase observations during heating and
convective refreezing of one synthetic specimen. The original chart supplies
temperature and liquid mass fraction throughout both directions of the cycle.

### Conserved multicomponent internal-energy and volume equilibrium

`thermochemical_gas::nonideal::flash::equilibrium::volume` recovers pressure,
temperature and phase amounts from the original feed's mass-specific U,V.
Each pressure trial reuses the existing caloric PH owner with
`h_target = u_target + p*v_target`; a safeguarded outer pressure solve closes
the original additive phase volume. The declared pressure endpoints must
bracket its volume residual from positive to negative. The off-root PH curve
is not assumed globally monotone. All sampled states retain the original
two-search stability and optional split evidence; failure or indeterminate
phase amounts refuse without a homogeneous fallback. This bounded bracketed
solution is not a global Gibbs minimum, entropy maximum or unrestricted
uniqueness certificate.

The increment path samples its reference once and uses
`Delta h_target = Delta u_target + p*(v_target-v_initial)
                  + (p-p_initial)*v_initial`.
Original NASA internal-energy integration plus the difference of original
phase U departures independently checks the final Delta-u; large species
formation datums cancel before rounding. Original entropy increments retain
the phase-mixing correction. No mixed heat capacity, extra latent source,
transport, acoustic relaxation or process heat/work law is introduced.
At identical represented volume and zero increment, an in-bracket reference
is preserved exactly after one PT selection with no pressure inversion.

Acceptance requires both direct U and V residuals and the pressure safeguard.
The inner PH tolerance and outer volume tolerance reserve parts of the energy
budget. A nonzero small energy increment also tightens that budget so it cannot
vanish into a large configured tolerance. Absolute energy or intermediate
enthalpy coordinates whose spacing cannot resolve tolerance refuse.
Temperature-source continuity is admitted once. Final state publication
retains the complete original source receipt and a final cancellation poll.

Separate limits bound pressure inversions, per-pressure PH selections and all
original PT selections in the entire call, including the reference and all
endpoints. Public work receipts distinguish attempted pressure/PT counts from
the sum of observed stability-screen EOS calls plus bounded split work.
The existing PH root and original caloric/phase owners are shared; intermediate
pressure trials do not rebuild complete material receipts.

The existing `nonideal_flash` target adds G0/G1/G3/G4 cases for independent
homogeneous and coexistence U,V oracles; heating/cooling cross-boundary caloric
increments; bitwise species-datum invariance; fixed-volume heat, pressure and
phase-amount response without latent double counting; exact one-call rest;
nested work/source/phase refusals; and final-cancellation replay. Synthetic
methane/ethane Cp/R fixtures prove constitutive consistency, not empirical
mixture-property accuracy. A consumer must still select its own mass, volume,
heat, work, phase-relaxation and entropy assumptions.

## Continuing heat-driven mixture phase equilibrium

`nonideal::flash::parcel::IsochoricMixtureEquilibriumParcel` retains one immutable
source model, conserved feed mass and fixed vessel volume. Initial conditions
may be supplied as PT (which defines the matching vessel volume) or as mass,
volume and original-reference specific internal energy through the U,V owner.
The retained phase state includes the original species/reference receipts and
the finite phase-stability evidence; no source inventory or phase label is
inferred from a homogeneous result.

`add_heat` supplies signed joules to the original internal-energy/volume solve.
Species formation and entropy datums cancel inside the source-integrated
increments, before rounding. Sensible and latent contributions remain parts of
one thermodynamic potential; no additional latent heat source is applied.
Phase compositions and amounts may change while the total feed, vessel mass
and volume remain unchanged. Original source volume remains visible within its
explicit residual allowance rather than being overwritten by the target.

Each accepted step retains actual energy and entropy increments, separate
energy/volume residuals, pressure and PT counts, and the original conservative
EOS evaluation bound. Supplied heat does not identify the reservoir temperature
or establish a process entropy-production theorem. Zero heat leaves every
stored field unchanged and performs no phase solve. Invalid heat/inventory,
exhausted work, source/phase refusal, failed physical balances and final
cancellation publish no mutation. Cloning the complete parcel retains the
source and accepted state for deterministic checkpoint/retry.

This is a closed, rigid, spatially uniform caloric equilibrium consumer.
Instantaneous phase equilibration does not supply evaporation kinetics,
wall transport, moving-boundary work, slip, chemical reactions, acoustic
response, global entropy-maximum certification or measured-fluid validation.
G0/G1/G3/G4 cases exercise independent binary cubic/fugacity isochoric heating
and cooling through phase change, partitioned heat/extensive scaling,
independent species-datum invariance, conserved initialization, unrepresentable
inputs, exact no-work rest, work refusal and atomic final-cancellation replay.

## Constant-pressure mixture phase change with mechanical work

`nonideal::flash::parcel::isobaric::IsobaricMixtureEquilibriumParcel` evolves a
closed original mixture against a prescribed constant quasistatic pressure
load. PT or original-reference PH initialization retains the complete model,
mass and accepted phase state. Volume is the original additive phase volume;
it may change while total feed composition and mass remain fixed.

Heat updates reuse the existing phase-aware PH increment owner. The original
source enthalpy increment equals supplied heat. Original NASA internal-energy
integration and phase U departures independently supply Delta U; boundary
work is the fixed pressure times the actual published volume change. Both
Delta H-Q and Delta U+p Delta V-Q must pass the explicit joule tolerance
before state and volume are updated together. Sensible and latent energy
remain within the same source potential, with no additional latent source.
Species formation and entropy datums cancel before increment rounding. A
small nonzero heat input tightens the inner enthalpy tolerance so a coarse
caller tolerance cannot erase it.

Reports retain signed heat, source enthalpy/internal-energy increments,
expansion work, volume and entropy changes, both residuals, actual PT work and
the original conservative EOS bound. Positive work leaves the mixture.
Zero heat preserves all accepted bits and performs no phase solve. Invalid or
unrepresentable inventory/heat, source/phase failures, exhausted work and
final cancellation leave the entire parcel unchanged. Cloning retains the
original immutable source and accepted state for deterministic retry.

This is a spatially uniform, closed, constant-pressure mechanical-equilibrium
consumer. It does not integrate piston inertia, varying external pressure,
heat-transfer rates, phase-relaxation kinetics, chemical reactions, slip or
acoustics. Entropy change alone does not determine entropy production without
the heat reservoir. Original finite phase-screen evidence is retained without
promotion to a global Gibbs-minimum or measured-mixture accuracy certificate.

The existing `nonideal_flash` target adds G0/G1/G3/G4 cases for independent
binary cubic/fugacity states across both phase boundaries; original Delta U,
Delta H and expansion-work closure; heating/cooling return paths;
partitioned heat and extensive scaling; large independent species datums;
small heat under coarse tolerance; initialization/refusal; exact rest; work
exhaustion; and atomic final-cancellation replay. Synthetic NASA calorics are
constitutive numerical tests, not experimental fluid validation.

## Pressure–entropy equilibrium through the original mixture phases

`thermochemical_gas::nonideal::flash::equilibrium::entropy` provides absolute
pressure/entropy and original-source entropy-increment constraints. The original
screened PT owner resolves homogeneous states or separated compositions and
amounts at every trial. Both stationary searches must converge; a witnessed
instability must resolve its split. This finite evidence is retained without
claiming a global phase minimum or resolving amount-indeterminate coexistence.

The existing safeguarded caloric temperature root now serves enthalpy and entropy
coordinates. It retains explicit source temperature bounds, monotonic branch
admission, residual and temperature safeguards, a bounded number of original PT
queries, conservative nested EOS ceilings, and cancellation before publication.
Existing PH and U,V consumers continue through the same enthalpy wrapper.
Neither inversion substitutes a frozen or invented mixed heat capacity.

Entropy increments integrate the original NASA temperature and pressure terms
before adding the original phase mixing/departure difference. Independent
species entropy and formation datums cancel before rounding. Accepted increment
solutions also report separately source-integrated Delta h and Delta u. Absolute
entropy refuses when the requested tolerance is below its floating-point spacing.
Nonzero entropy increments tighten the root tolerance to at most a quarter of the
requested magnitude. At unchanged pressure, zero increment, and a reference
inside the declared bracket, the exact reference returns after one PT query;
changing pressure at zero Delta s still solves the isentropic endpoint.

This material-level constraint supplies thermodynamic endpoints. Process work,
heat transfer, flow rate, phase kinetics, reactions, and acoustic closure remain
with their physical consumers. Source or work refusal and cancellation publish
no partial result.

Eight focused regression cases cover independent homogeneous and split calorics,
an independently solved binary expansion and reversible compression, both phase
boundaries, source datum invariance, small signed entropy changes, exact rest,
physical/source/work refusals, and final cancellation/replay. Source and tests
passed Rustfmt parsing. Native compilation and execution remain pending because
the shared build filesystem is full; the new suite is not claimed green.

## Entropy/volume closure through the original mixture equilibrium

`thermochemical_gas::nonideal::flash::equilibrium::entropy::volume` recovers
pressure, temperature and phase amounts from specific entropy and specific
volume at conserved complete feed composition. The increment API integrates
original NASA entropy before phase mixing/departure differences, preserving
independent species reference-datum invariance. It reports separately integrated
Delta h and Delta u; these are thermodynamic changes, not an independent process
work quadrature. Zero Delta s supplies a reversible adiabatic endpoint.

The original safeguarded pressure bracket now serves both U,V and S,V closures.
Each trial pressure uses the original bounded caloric temperature root and
screened PT phase law. The original U,V coordinate and energy tolerances remain
owned by its existing adapter. No mixed heat capacity, gamma, kinetic phase
rule or new generic nonlinear solver is introduced. Every sampled phase must
pass both original stationary searches and any required composition split;
finite phase evidence is not a global equilibrium or unrestricted uniqueness
certificate.

One original reference sample is retained across all pressure trials. Full
material receipts are rebuilt only for the accepted endpoint. Explicit pressure
and temperature brackets, entropy/volume residuals, a pressure correction or
bracket safeguard, per-pressure PT limits and whole-call pressure/PT ceilings
bound the solve. Work accounting includes reference and endpoint selections,
and per-pressure refusals retain the whole-call PT count. Small requested entropy
or volume changes reserve residual resolution against their magnitude. At zero
entropy increment and identical volume, an in-bracket reference returns exactly
after one PT selection and no pressure inversion. Cancellation and refusal
publish no partial endpoint.

Six focused regressions recover independent liquid/coexistence/vapor states,
the independently solved binary isentropic expansion and reverse compression,
phase-boundary changes under large independent species datums, exact rest and
whole-call budget boundaries, invalid/source/bracket refusals, and final
cancellation/replay. They also check Delta h - Delta u against the actual change
in p*v. Rustfmt parsing passed. Native type checking passed with
`cargo check --offline --locked -p fs-material --test nonideal_flash`, including
the six new regressions. Runtime execution remains unverified: the native test
build was killed by the 8 GiB memory limit while compiling upstream asupersync,
before any test ran. No executed-test pass is claimed.

## Prescribed-volume reversible mixture equilibrium

`nonideal::flash::parcel::isentropic::IsentropicMixtureEquilibriumParcel` retains
one original mixture model, conserved feed mass, prescribed physical volume and
accepted phase state. PT initialization defines its matching volume; S,V
initialization retains the caller's volume and original entropy reference.
`change_volume` uses the original source-increment S,V owner with zero entropy
change. Temperature, pressure, phase compositions and amounts can change while
the complete species inventory remains fixed.

The process is explicitly closed, spatially uniform, reversible, adiabatic and
quasistatic. Boundary work BY the mixture is inferred as minus the original
source-integrated Delta U. This is a thermodynamic endpoint work convention,
not independent pressure-work quadrature or a time-discrete mechanical energy
audit. The step also reports independently integrated Delta H, actual Delta S,
the original EOS volume minus the prescribed published volume, and all original
pressure/PT work counts and conservative EOS ceilings. Species formation and
entropy datums cancel before increment rounding; latent energy remains within
the original thermodynamic potential.

Small nonzero motion reserves entropy and volume resolution against its declared
pressure-work scale. These numerical safeguards do not certify a relative work
error. A motion below the retained initial volume resolution, nonrepresentable
inventory or tolerances, zero or wrong-sign inferred work, failed entropy or
physical-volume closure, source refusal, exhausted work and final cancellation
publish no mutation. State and physical volume change atomically after the final
poll. Nonzero changes require the accepted original phase policy, including its
flash and stability controls, so reference reevaluation cannot silently rebase
stored caloric energy or entropy. Source rectangles and outer caloric tolerances
may change without changing that phase policy. Exact unchanged volume preserves
every accepted bit with no phase query
and reports its retained source-volume residual without refining that state.
Its valid controls are admitted but a different phase policy is not applied.

This consumer does not infer piston force, inertia, elapsed time, heat transfer,
phase kinetics, slip, chemical reactions, acoustic response or a global phase
minimum. Its finite phase-selection evidence retains the original source scope.
Six focused G0/G1/G3/G4 cases cover an independent binary isentropic expansion
and compression, three-rung pressure-volume quadrature against source work,
partitioned motion, extensive scaling, small motion under loose controls,
large independent species datums, entropy initialization, exact rest,
unrepresentable inputs, bounded work and final-cancellation replay.

Native validation: `cargo check --offline --locked -p fs-material --test
nonideal_flash` passed with both the S,V closure and all six parcel regressions
registered. Rustfmt and whitespace checks passed. Runtime test execution remains
unverified following the upstream compiler memory failure described above.

## Frozen finite-strain viscoelastic potentials for spatial equilibrium

`visco_hyper::ViscoHyperelastic::freeze(accepted, dt_s)` borrows the original
native law, complete accepted `ViscoState` and one nonnegative physical interval.
The resulting immutable `FrozenViscoHyperelastic` implements
`hyper::cards::FiniteStrainPotential` using the original branch minimization,
Piola stress and symmetric algorithmic tangent. Newton and backtracking trials
cannot advance the accepted viscous metrics or elapsed interval. Stress and
potential evaluation share one branch preparation.

`step_response` proposes the same native update as `step` and exposes its
reduced incremental potential, physical stored energy, frozen-state loading-work
increment and nonnegative relaxation-loss increment separately. The increments
are retained before cumulative-ledger rounding. The incremental potential is
the function whose gradient gives the trial Piola stress; it is not physical
stored energy or heat. `dt_s = 0` retains the native elastic-jump semantics.
Nonfinite work/storage, negative relaxation loss, invalid histories and original
constitutive refusals publish no candidate. No alternate material equation,
fitted relaxation parameter or mechanical integration convention is introduced.

Neo-Hookean and Ogden branch storage remove the analytically vanishing linear
terms of the principal logarithmic strains before summation. Nearly relaxed
branches therefore retain their positive quadratic elastic energy, consistent
with the cancellation-free stress and loss calculations, in both full and
remainder endpoints. This does not recover strain information already rounded
out of the supplied deformation or stored inverse viscous metric.

Near isotropic trial metrics, principal logarithmic stretches use a centered
metric increment and the shared small-strain `log1p` evaluation. Representable
tiny shear is resolved before forming unit-size stretch eigenvalues. Relaxed
metrics likewise sum small exponential increments before adding the identity,
preserving viscous memory through repeated small-strain holds. Larger
distortions retain the original spectral path; no linear constitutive surrogate
replaces the finite-strain energy, stress or algorithmic tangent.

`tests/visco_incremental.rs` checks directional energy/stress/tangent consistency
after nontrivial history, unchanged accepted states under repeated trial queries,
and the distinct potential/storage/heat accounts during held-strain relaxation.
The spatial consumer is `fs-solid::hyper3::visco`; its physical endpoint-work
quadrature discrepancy belongs to that consumer and must not be booked as heat.

## Original-energy mixed displacement/pressure materials

`hyper::mixed::MixedFiniteStrainPotential` supplies the original bulk-free
remainder plus `p*(J-1)-U_star(p)`, with tensile-positive `p=U'(J)`.
`Hyperelastic` and an immutable `FrozenViscoHyperelastic` implement the same
state-only and full-Hessian queries. Pressure coupling is exactly `cof(F)`;
the displacement block retains `p*d²J/dF²`, and pressure compliance is the
negative inverse original bulk curvature. Preparation admits parameters/history;
actual point queries retain the native deformation and stretch limits.

`OriginalMaterial` preserves the supplied quadratic or log-squared bulk law.
The latter uses only its locally convex branch `ln(J)<1`, including an explicit
bounded pressure inverse; neither its other branch nor a global Legendre
certificate is claimed. `ExactIncompressible` supplies `J-1`, zero compliance
and zero bulk storage without penalty arithmetic. Coupled Neo-Hookean refuses.
Full-invariant HGO fibres remain in their original remainder, so this multiplier
represents only the separable bulk contribution, not necessarily mean stress.

The response distinguishes the saddle potential, remainder potential, original
bulk storage at pressure-reconstructed volume, and complementary potential.
`mixed_remainder_free_energy` evaluates physical equilibrium/branch storage
directly, avoiding subtraction of large native bulk energies. A frozen viscous
incremental potential is never relabeled as stored energy or heat. Spatial
owners must use their weak pressure discretization's bulk storage and establish
their own stability, quadrature and convergence evidence.

`tests/hyper_mixed.rs` verifies independent pressure-dual derivatives, recovery
of original energy/stress/tangent by pressure elimination, every mixed Hessian
block, exact-constraint bulk independence, physical versus incremental storage,
original isochoric histories, rotation/unit scaling and domain refusals.

## Retained mixed-viscoelastic remainder histories

`visco_hyper::{ViscoRemainderState, FrozenViscoRemainder,
ViscoRemainderStepResponse}` give the mixed spatial owner a distinct physical
remainder ledger. `initial_remainder_state`, `freeze_remainder`,
`step_remainder_response` and `remainder_free_energy` retain the native
inverse viscous metrics, actual deformation gradient and branch relaxation law.
The existing `ViscoState` continues to mean full material energy.

Both history types share the original branch preparation and return map.
The remainder endpoint evaluates the equilibrium remainder and original branch
energies directly: it never evaluates or subtracts the unused bulk energy,
constructs a surrogate full-energy state, or projects the actual deformation to
unit determinant. Full-invariant HGO fibres therefore retain their actual-volume
dependence. The spatial mixed owner supplies pressure stress and reconstructed
bulk storage separately.

A frozen remainder borrows one accepted history and one nonnegative physical
interval throughout Newton trials. A proposed step leaves that history unchanged.
Its physical storage, frozen-metric loading work, incremental potential and
nonnegative branch loss have separate fields. Malformed histories, nonfinite
ledgers, negative loss, invalid intervals and unsupported coupled laws refuse.

`tests/visco_remainder.rs` checks native-history/loss equivalence, the physical
bulk difference where representable, immutable trials, tangent consistency,
full-invariant fibres, refusals, and a discriminating `kappa=1e308, F=3I` case:
the full bulk energy overflows while the remainder endpoint stays finite and
matches `kappa=1e3` bit for bit. These are isothermal material accounts; neither
the incremental potential nor a signed numerical work defect is heat.

## Original elastic-source derivatives in finite thermoelastic storage

`hyper::thermoelastic::source::ResolvedNeoHookeanThermalCoefficients` binds the
selected native compressible Neo-Hookean model to the original same-card,
same-context Young's-modulus and Poisson-ratio sources at the actual absolute
temperature. Pinned source lines retain every answer and receipt. Both paths
must be affine over an explicit smooth two-sided interval, with E positive
and nu in (-1,1/2). The derived current mu and lambda must be finite and match
the selected model to roundoff. Missing evidence, source kinks, different
contexts, unsupported model families and inconsistent coefficients refuse.
The caller qualifies these sources as thermally relaxed infinitesimal moduli;
a prestressed tangent measurement is not silently reinterpreted.

Exact rational first and second derivatives of mu(E,nu) and lambda(E,nu)
supply the original native energy's partial W_T, W_TT and P_T at fixed elastic
deformation. The existing `FiniteThermoelasticState` natural-map chain rule
then includes both mixed expansion/coefficient terms in the original-volume
Helmholtz derivative. It derives mechanical entropy, internal energy,
fixed-total-F capacity `-T Psi_TT` and the complete stress-temperature tangent
from that one potential. Signed derivative coefficients are used only in the
derivative of the original energy; they are not substitute mechanical laws.
The existing fixed-coefficient entry and exact constant-source limit remain.

The retained affine sources also support explicit
`thermoelastic_state_at_temperature` evaluation throughout their original
interval. Each query rebuilds the actual native coefficients, their exact
derivatives, stress and tangent at the requested temperature with the original
stretch envelope. The original model metadata remains anchored at its
preparation temperature. The original caloric owner and reference density
still supply base storage; total capacity must remain positive before the
existing adiabatic acoustic response is admitted.

`FiniteThermoelasticState` retains its native law and exposes direct free-energy
and internal-energy changes between compatible finite states. These reuse the
original native finite-pair potential and reference-volume increment owners;
a symmetric compensated change of T times mechanical entropy preserves weak
caloric coupling against large stored-energy baselines. No difference of two
rounded absolute internal energies is required. The finite-solid consumer
owns the heat transaction and its branch limits. Focused native checks compare
independent potential/entropy/internal-energy and stress derivatives,
noncommuting natural maps, constant limits, coefficient-curvature Cp and actual
wave responses; finite heating tests additionally retain a 5,000 J/m3 change
when both absolute internal energies round to the same value near 1.125e20.

## Objective discrete-gradient elastic work

`hyper::discrete_gradient::DiscreteGradientPoint` freezes an accepted deformation
and evaluates a two-state elastic response using the original hyperelastic
energy. `EnergyPart::Full` retains the complete Neo-Hookean, Mooney-Rivlin,
Ogden, or HGO energy; `Remainder` uses the existing mixed split without evaluating
an unused bulk modulus. Full-invariant fibres retain their actual volume
dependence.

The Gonzalez correction is applied in Mandel coordinates of
`C = F^T F`, around `(C_before + C_after)/2`. A Cholesky representative feeds
the original scalar energy and native automatic differentiation. The returned
algorithmic Piola stress satisfies `P_alg : (F_after - F_before) = Delta W`
within the caller's explicit local work allowance. Its material moment
`P_alg F_mid^T` is symmetric, which supports angular-momentum preservation
when paired with centered reference inertia. The endpoint derivative is the
true, generally nonsymmetric Jacobian of this two-state stress; it is not a
symmetric elastic Hessian or an endpoint potential. The full C-chain derivative
also applies when equal endpoint strains are related by a finite rotation.

`DiscreteGradientState::energy_change_pa` retains the original finite-pair
change independently of the two absolute storage values. The correction and
work-defect gates use that change, including the original-volume reference
pullback. Symmetric endpoint products form `Delta C` without subtracting rounded
absolute strain tensors; compensated midpoint and work contractions retain
the corresponding weak-work identity under large bulk stress. The actual
reported f64 stress must still pass the caller's work allowance.

`DiscreteGradientSettings` explicitly supplies the small-increment switch,
absolute work tolerance in Pa, and relative work tolerance. The midpoint limit
reports its actual work defect and refuses an insufficient allowance.
`volume_state` and `volume_response` separately expose the objective discrete
cofactor, its endpoint Jacobian, and the endpoint cofactor. These local volume
identities do not establish a global pressure gauge or supply a complete mixed
time integrator.

G0/G3 evidence in `tests/hyper_discrete_gradient.rs` covers original-energy
work, all nine endpoint tangent columns, volume work and derivatives, moment
symmetry, frame covariance, equal-strain finite rotations, small-increment
refusal, invalid domains, and unused-bulk independence. The focused native-source
lane passed all five tests; its worst normalized material finite-difference
tangent discrepancy was `2.74e-10`. This is local constitutive evidence, not
a long-horizon or frequency-matching claim for a solid integrator.

Reference: O. Gonzalez, *Exact energy and momentum conserving algorithms for
general models in nonlinear elasticity*, CMAME 190 (2000), 1763-1783,
doi:10.1016/S0045-7825(00)00189-4.

## Original finite-bulk discrete gradients

`hyper::discrete_gradient::bulk::DiscreteBulkPoint` freezes an original
quadratic or log-squared bulk law and one accepted quadrature volume. Interval
pressure is the original energy secant
`p_alg = (U(J_after)-U(J_before))/(J_after-J_before)`. It is distinct from
endpoint `U'(J_after)`. Quadratic inversion is analytic; the logarithmic law
uses bounded safeguarded Newton/bisection on an explicit positive bracket
strictly below `J=exp(1)`. No bulk coefficient or energy is replaced.

`DiscreteBulkSettings` supplies the volume bracket, scalar iteration cap
(1 through 128), and separate absolute/relative pressure and work allowances.
Cancellation is the spatial owner's quadrature-boundary responsibility, with
each scalar query bounded by that cap. Cancellation-safe nearby logarithmic
secants and their derivatives retain the coincident limit. The response
separates original endpoint storage, stable energy change, endpoint pressure,
interval work, positive inverse derivative, pressure residual, actual work
defect and iteration count. Both local pressure and work gates must pass.

`from_endpoint_pressure` initializes through the existing original
`pressure_potential`; the spatial owner must still admit the weak geometric
volume constraint. Compatible initialization matters: arbitrary geometric
quadrature volumes can contain unresolved pressure-space components that the
quadratic secant reflects on alternating steps. Accepted retained volumes are
never silently replaced by geometric determinants. `rebind` creates a new
immutable point on the admitted branch. Coupled Neo-Hookean and invalid original
parameters refuse.

G0/G1 tests in `tests/hyper_discrete_bulk.rs` cover original energy/work,
independent pressure integration, quadratic and logarithmic compression and
expansion, inverse derivatives, equal/near-equal volumes, heterogeneous moduli,
reversal, compatible initialization, branch/fold/domain refusal, bounded work,
and independent pressure/work tolerances. The focused native-source lane passed
all five new tests and all five existing objective discrete-gradient tests.
The largest relative inverse-compliance finite-difference discrepancy was
`4.43e-8` for the fixed probe including a near-fold fixture. This is a local
constitutive closure; spatial stability, pressure-space consistency and temporal
convergence belong to its mixed consumer.

### Adaptive competitive membrane breakthrough

`membrane::transient::slab::competitive::adaptive` advances all competing
channels and optional finite molar-pressure reservoirs on one accepted
timeline through the original shared membrane controller. Every native probe
retains the coupled shared-site solve. Error estimates compare each species'
mobile pressure and amount, trapped amount, bath pressure and amount, and
integrated face transfer. An additional shared-vacancy comparison uses
compensated occupancy changes, preserving small changes near saturation.

Algebraic residuals, including independent site kinetics, must remain below
temporal pressure/storage and duration-weighted balance allowances. Native
tolerances are tightened to the strictest channel capacity, while relative-only
requests retain their original controls and must pass independent absolute
audits. Accepted half steps alone contribute signed/gross face totals.
Original-to-final cell, membrane, bath and closed-species accounts are audited
separately; a quiet first species cannot hide another species' transient.

`adaptive_workspace(maximum_bytes)` applies its byte ceiling to both retained
native workspaces together. Work limits, failed native/source checks and
cancellation publish no partial physical state. The private
`CompetitiveNative` seam lets source-bound owners retain their physical
admission on every joint probe. These are local temporal estimates for fixed
coefficients, geometry and temperature, not spatial or whole-window error
certificates, mixed-gas source qualification, swelling or sorption calorics.

Focused G1/G3 tests compare competitive loading with an independent Riccati
solution, finite-bath evolution with continuous-time eigenmodes, nonfirst
species error control, competitive displacement and original face balances.
Late work/cancellation refusal preserves the original state and clock.

## Competitive membrane diffusion against actual finite gas mixtures

`transient::mixture::competitive` attaches one explicit original
`CompetitiveAdsorbent` lattice per layer to an existing source-qualified
`IsothermalGasMembraneMixture`. `with_competitive_sites` requires an explicit
trapped profile for every canonical gas species, including zeros for species
without a site law. The original competitive owner checks common capacity,
dry-skeleton calorics and reference temperature; independent capacities are
never summed. Every lattice species must exist in the original gas/mobile
roster and retain its original gas-source identity.

`IsothermalCompetitiveGasMembrane::step_candidate` advances mobile Henry/Fick
storage, competing finite sites and both finite reservoirs simultaneously.
The prepared numerical owner retains its bounded nonlinear/linear work and
workspace byte ceiling. One original source-relative gas reconstruction per
endpoint accepts all species streams together. Actual EOS partial pressures
become the next reservoir coordinates only after pressure-to-mass closure.
Separate species accounts include mobile and trapped storage, both actual gas
inventories, signed stream enthalpy and explicit gas thermostats; aggregate
mass conservation cannot conceal species conversion.

Original mobile/interface and adsorption source bindings are re-admitted
across old, numerical and actual endpoint envelopes. Partial pressures drive
transport and adsorption, while summed polymer pressures and actual gas EOS
totals qualify source applicability. Finite-resistance polymer traces remain
distinct from gas traces. Unowned composition, humidity and concentration
dependence refuses; this adapter does not promote pure-species measurements
to mixed-gas evidence. Failure or cancellation publishes no physical endpoint.

The model is isothermal in a rigid, stationary membrane frame. Gas enthalpy
and thermostat reports exclude an unrepresented mobile/trapped caloric store.
No heat of sorption, swelling feedback, Maxwell–Stefan gas-frame diffusion,
cross diffusion or unrestricted mixed-gas validity is inferred.

### Adaptive competitive breakthrough with actual gas reconstruction

`transient::mixture::competitive::adaptive` gives
`IsothermalCompetitiveGasMembrane` a complete-duration adaptive candidate and
two retained native workspaces. It reuses `MembraneMixtureAdvanceControls`
and the original competitive time controller. Every coarse and half-step
probe runs the full source-qualified gas/membrane candidate, including the
original gas reconstruction, pressure-capacity closure, common-site solve
and actual source envelopes. Actual gas species inventories receive their
own temporal comparisons, in addition to all native pressure, mobile,
trapped, vacancy and face-transfer comparisons.

Only accepted half-step histories contribute signed and gross species
transfers, original stream enthalpies and explicit gas thermostats. Fixed
whole-window species, total-mass and gas-energy budgets constrain native
probes through duration-weighted algebraic shares. Stored native controls
remain hard ceilings and are restored after each proposal. Final audits
independently compare actual original-to-final gas species amounts,
mobile/trapped storage, accepted gas-boundary transfers, whole mass and the
original source-relative gas thermal increment. A large formation-enthalpy
datum therefore does not contaminate the thermal residual.

The maximum total cell count and both native workspaces have explicit bounds.
Rejected proposals, source/native refusals, exhausted work and late
cancellation leave all original gas and membrane inventories, time and
source receipts unchanged. Local temporal estimates remain distinct from
global or spatial error certificates. Temperature and geometry stay fixed;
gas-side caloric accounts do not supply sorption heat or whole-membrane
energy closure.


## Steady rougher-cleaner recycle with explicit purge

`slurry::recycle::RougherCleanerCircuit` composes two original selective separators,
a source-bound fresh slurry, and a declared nonselective cleaner-tail return in
[0,1]. It solves each exact dry-component identity and the separate carrier
analytically from its positive exit paths, then replays the existing separator,
splitter and mixer operations. Results retain every density query/receipt, both
separator prescriptions, the original circuit provenance, all internal streams,
final product, purge, combined external tailings, per-component circulating
loads/recovery, and internal/external conservation residuals.

The external residual is checked relative to fresh feed, never inflated recycle;
the mixer residual is checked relative to solved mixed flow (relative tolerance
2e-12). No corrective material stream or coefficient normalization is introduced.
Missing/extra partitions, blank circuit provenance, invalid return fractions,
nonfinite or erased positive flows, and unacceptable closure refuse. A component
with no external exit refuses even at zero fresh flow because its circulating
state is nonunique. Zero fresh flow otherwise has no defined circulating-load
or overall-recovery ratio. Inputs are immutable and work uses the existing
128-solid bound and deterministic component order.

This is a steady prescribed-partition balance. The supplied partitions must be
applicable to the recycled blend; density evidence does not establish separator
calibration, capacity, liberation, particle breakage, reagent demand or kinetic
recovery. `mineral_recycle` demonstrates a closed circuit with real mineral/water
density packs and explicitly authored stage partitions. Focused tests cover the
independent flow solution, both balance boundaries, evidence retention, zero/full
return, small positive exits, and trapped or numerically ill-conditioned loops.


## Grinding-circuit underflow purge

`comminution::circuit::ClosedGrindingCircuit::solve_with_underflow_recycle`
extends the existing full-return solve with a nonselective split of classifier
underflow. The caller supplies the existing `Partition`; product returns to the
mill and reject is an explicit external purge. Both original fractions are
retained independently, including a small purge beside rounded unit recycle.
The circuit still uses the original size-redistribution and classifier owners.

The triangular equations use effective return s*r and the positive exit sum
`t + s*(q + r*fragmented_fraction)`, with separate water escape
`t + s*water_overflow`. A positive purge can discharge a previously driven
trapped loop. Zero purge preserves the original full-return semantics and trap
refusal. No coefficient normalization, new breakage law or hidden makeup is
introduced. The solved actual mill/discharge conditions must remain inside
both original source domains; the circulating-load ceiling applies to returned
underflow, not the raw classifier underflow before purge.

Results retain raw underflow, recycle, purge, overflow, and the caller's split.
Local split and mixer accounts and the combined external dry/element/water
boundary are independently checked with the original numerical settings.
`product_minus_fresh_*` remains the literal overflow-minus-fresh difference;
`external_minus_fresh_*` reports overflow plus purge minus fresh and is the
conservation residual when purge is nonzero. Purge is not fine product.

The existing `mineral_circuit` consumer accepts an optional purge fraction,
prints every outlet, and forms any circuit-scale Bond estimate from fresh
throughput and the combined external PSD. Zero dry flow has no P80. This
prescribed split supplies no calibrated bleed-control policy, metal recovery
optimization, selective rejection or purge-treatment prediction.

## Continuous population-balance grinding

`comminution::population_balance::ContinuousPopulationBalanceMill` solves the
steady, perfectly mixed size-mass balance from original per-class selection
frequencies and conditional daughter mass fractions. `at_residence_time` uses
an admitted common residence; `at_dry_holdup` derives it from positive dry
hold-up and actual feed. The triangular coarse-to-fine solve is bounded to
128 classes. Source context covers both rates and daughters, with explicit
typed support for residence and dry flow; actual solids fraction is also owned.

Every measured solid-phase element follows its parent mass congruently.
Carrier water is unchanged at the outlet and shares the same residence only
under the explicit perfect-mixing assumption. Hold-ups equal outlet rates times
residence. Original daughter fractions close within eight binary64 epsilons
and are never renormalized; floating-point results are not interval certificates.
Bad layouts, nonconservative/coarsening daughters, unsupported source points,
nonpositive operating conditions and unrepresentable arithmetic refuse without
mutating the feed. No ore-rate inference, liberation,
size-dependent discharge, media dynamics, power or wear prediction is supplied.

The six focused population-balance tests cover independent two-class and
multiclass balances, every dry/element/water account, throughput/residence
response, source admission and invalid kinetic declarations. The executable
`continuous_ball_mill` example uses explicitly synthetic rates and distributions.

`advance_inventory` extends the same source-bound local law to finite time.
It retains class dry/element masses and carrier water in `transient::MillInventory`.
Constant-feed intervals use bounded positive backward-Euler substeps; integrated
discharge uses those same accepted substeps. Final inventory, instantaneous
outlet, integrated input/output and actual dry/element/water closure residuals
are returned together. The caller can carry that inventory through startup,
feed changes and washout. Actual feed and residence retain the original source
admission; an empty feed cannot provide a required solids-fraction coordinate.
Cancellation, insufficient work budget or failed conservation publishes no
partial state. This is a first-order numerical trajectory, with an explicit
step limit; a small balance residual is not a trajectory-error bound.
The closed-circuit transient adapter below supplies instantaneous classifier
recycle. Rate-power calibration and plant transfer remain outside this local model. The `transient_grinding` example and seven independent
analytical/refinement/conservation regressions exercise the complete path.

## Simultaneous shared-reagent batch kinetics

`reaction_batch::network::BilinearBatchNetwork` joins 1..=32 original bilinear
models sharing one explicitly named consumed reagent. Every other reactant is
unique; no product can feed a reactant in this scope. Rates use the same staged
inventory and fixed volume/temperature. Original model cards, mappings and
conserved equations remain attached to the completed result.

RK4 step doubling advances reaction extents simultaneously in canonical
reaction-ID order. Every trial inventory is reconstructed by `ReactionNetwork`
from the original stocks and total extents. Inputs and all source conditions
are checked at actual concentrations and absolute elapsed times, including the
full requested duration. Positive initial reactants remain positive for finite
time; invalid trials shrink the step, with no recovery clipping or priority.

Callers provide a maximum step, absolute/relative accumulated local-error
estimate budget and at most 65536 attempts. Reported error estimates are not
certified global bounds. Source departure, exhausted work, unrepresentable
progress or invalid topology returns no partial batch and does not mutate input.
Scaled RK weights preserve positive subnormal rates; unrepresentable positive
error estimates refuse. This is not a stiff solver, equilibrium/speciation
model, calibrated ore response, temperature feedback or reaction feedback loop.

Focused tests compare to independent closed forms, exercise unequal proton
coefficients, preserve permutation invariance, and cover finite shared acid,
source-time support, numerical range and atomic refusal. The executable
`competing_mineral_leaching` uses explicitly synthetic CuO/calcite kinetic cards
with the source-defined neutral carbonic endpoint for the gangue reaction.


## Complete isothermal competitive-membrane calorics

`transient::mixture::competitive::calorics` closes the mobile, trapped and actual-gas heat accounts through `with_caloric_closure`. Every original species/layer material explicitly supplies `gas_membrane_mobile_internal_energy_relative_to_gas_enthalpy` [J/mol]: condensed mobile internal energy minus the original gas molar enthalpy at the held temperature. Original Langmuir calorics supply the bound-phase difference through their adsorption heat. Solubility and gas heat capacity do not substitute for this condensed-phase datum.

The extra property is pinned into the existing D/S source bundle. Every fixed or adaptive native probe therefore enforces its actual total-pressure domain using the original source. Missing data and unsupported source trajectories refuse. Summed mobile and bound storage determines membrane thermostat heat. Original gas thermostat accounts separately remove pressure-flow work, while actual original/final gas species inventories close a whole-system first-law residual under an explicit absolute joule allowance. Conserved species formation references cancel analytically; original independent mass/species budgets still apply.

Layer/species dissolution and trapping accounts decompose storage and count binding once; they are not interlayer thermal-flux predictions. The dry skeleton has no sensible change at held temperature. The original transport solver, adaptive controller, accepted-only histories, finite work budgets and immutable failure behavior are retained.

Three focused G0/G3 groups passed on the native material test executable: fixed/adaptive mobile heat is unchanged by a 1e24 gas-reference shift; layered competitive displacement agrees with an independent internal-energy reconstruction; and missing/late-invalid caloric sources and late cancellation publish no state. These are numerical implementation checks, not experimental validation or a DSR verdict. The existing Langmuir source's representable reference calibration remains its own obligation.

Scope: stationary rigid geometry, additive per-species calorics and external thermostats. Temperature evolution, sorption heat capacity, thermal conduction, excess-mixture energy and entropy production are not claimed by this held-temperature account.

## Binary liquid miscibility and shared phase transport

`liquid::solution::BinaryLiquidSolutionLaw` binds the L1
`BinaryRegularSolutionV1` to one original liquid-property bundle. The first
species' composition coordinate must explicitly carry the mole-fraction
semantic kind. Temperature and composition refresh together; pressure and every
other coordinate remain fixed, and every selected property claim stays pinned.
The caller supplies compatible species-pair parameters, component molar masses
and mixture-property data. Neither a material label nor a density determines
miscibility.

`equilibrate` returns one homogeneous liquid or two coexistence phases with
their actual compositions. The molar lever rule is converted to separate
component amounts, mass and source-derived volume for each phase. Small
component allocations are evaluated directly, with the larger allocation
retained as the complement of the original inventory. Conservation holds to
floating-point roundoff. The smaller feed mole fraction is normalized directly
and retained by the thermodynamic owner instead of being recovered by subtracting
from one. A nonzero component that vanishes in the source's scalar first-species
coordinate refuses. Invalid or unrepresentable inventories, a source
failure in either phase, and cancellation publish no partial equilibrium.

Every `ResolvedBinaryLiquidPhase` implements `FluidTransport` and
`ThermalFluidTransport` using its own density, heat capacity, conductivity and
viscosity source evaluations. Kinematic viscosity is converted with that same
phase's density. The actual molar concentration is available to the existing
Maxwell-Stefan driving-force path. No pure-fluid property average, additive
volume assumption or invented expansion coefficient is introduced. A
two-phase equilibrium intentionally has no homogeneous transport tuple.

G0/G3 tests in `tests/liquid_miscibility.rs` cover temperature-driven splitting
and remixing, unequal molar masses, independent species/mass/volume balances,
source claim retention, nonlinear composition tables, pure limits, wrong
composition basis, unsupported sources and atomic cancellation. Fixtures are
synthetic numerical checks, not measured pair fits. The prescribed-temperature
entry point is bulk binary liquid equilibrium at fixed pressure. The separately
bound caloric path below closes thermal changes. Both paths exclude vapor/solid
phases, finite-rate mixing, nucleation, capillarity and spatial phase separation.

### Finite isothermal mixing and its heat duty

`BinaryLiquidSolutionLaw::mix_isothermal` combines two finite equilibrium
inventories at the same temperature and fixed pressure. Both must retain the
same original pair law, ordered species/molar masses, material claims, query
context and source roles. Immutable binding data are shared across results;
identically reconstructed bindings are accepted. Independent component sums
enter the original equilibrium/source resolver, so the result carries actual
phase amounts, masses, source-dependent volumes and transport properties.

The required heat **into the liquid** is the change in equilibrium mixing
enthalpy, with the actual coexistence phases included. Pure-component caloric
datums cancel because each species and temperature are unchanged. The source
mixture Cp is not a second heat contribution. Gibbs and liquid entropy changes
are retained separately; a reservoir at the prescribed temperature receives
`-Q` and `-Q/T`. Combined entropy generation is `-delta G/T`, equal to
`delta S_liquid-Q/T` to roundoff. Common molar contributions are subtracted
before inventory scaling, so identical intensive states incur exactly zero
heat. Raw near-zero signed roundoff is not clamped into a physical source.

G0/G1/G3 tests in `tests/liquid_miscibility/mixing.rs` cover analytic ideal,
endothermic and exothermic mixing, phase-split heat duty, source-resolved volume
changes, both component balances, input interchange, immutable source/law
ownership, final-source refusal, overflow and cancellation atomicity. This is
an equilibrium isothermal heat duty, not an adiabatic temperature solve,
kinetics, stirring work or an interfacial-energy model. The external heat owner
must accept the opposite account; no transport Cp consistency is inferred.

### Source-bound caloric mixing and thermal storage

`liquid::solution::caloric::BinaryLiquidCaloricLaw` binds the original mixture
law to two explicit pure-component `IsobaricHeatCapacityLaw` sources. Species
order, common datum and the original typed static pressure must match. The
adapter preserves each source law, identity and usage receipt, converts every
original Cp breakpoint from mass to molar units using the declared molar mass,
and delegates all thermal balance to `fs-thermochem::solution::caloric`.
The caller qualifies the sources as the named pure liquids; a label or the
dimensions of Cp cannot establish that chemical qualification. Transport Cp
is neither used as a second excess heat term nor certified as consistent with
the independently supplied pair free energy.

`equilibrate`, `mix` and `add_heat` preserve both component inventories and
immutable source ownership. Different-temperature mixing and heating/cooling
resolve the temperature from sensible plus nonideal mixing enthalpy. Numerical
temperature trials perform no material-source queries. The final material
refresh consumes the EXACT accepted thermodynamic equilibrium, avoiding another
coexistence solve or normalization. Its phase-specific density, viscosity,
conductivity, volumes and masses therefore belong to the accepted caloric
temperature and compositions. A failed final source refresh, mismatched owner,
unrepresentable sum or last-poll cancellation publishes no result.

`caloric::homogeneous::ResolvedHomogeneousBinaryLiquidFeed` explicitly retains
an initially homogeneous parcel, including metastable and spinodally unstable
compositions. `BinaryLiquidCaloricLaw::homogeneous` queries its actual initial
transport sources and preserves its original homogeneous excess enthalpy; it
does not replace that energy with an isothermal equilibrium. The distinct feed
type exposes its local thermodynamic response, conserved component amounts,
source-derived mass/volume and transport state without claiming phase
equilibrium or entering equilibrium-only tank consumers.

`relax_homogeneous` delegates signed heat and that unchanged initial caloric
state to the existing L1 energy solver. Zero heat therefore permits demixing
with its accompanying temperature change, while supplied heat can cross the
miscibility boundary. Only the accepted equilibrium refreshes the original
phase transport sources, and the returned balance remains compatible with
ordinary heating, mixing and thermal storage. Source mismatch, incomplete
phase-domain coverage or cancellation leaves the feed unchanged. This is a
closed isobaric equilibrium endpoint, without separation kinetics, nucleation,
spinodal dynamics or interfacial energy. Synthetic regressions cover unstable
and metastable feeds, analytic variable-Cp heat balances, remixed endpoints,
existing storage/consumer compatibility and atomic source refusals; native
integration execution remains required to validate those cases.

`caloric::storage::BinaryLiquidThermalStorage` makes the same fixed-composition
equilibrium available to the existing `fs-thermal` enthalpy timeline. The
`thermal_storage` constructor uses the cold equilibrium endpoint as zero
specific enthalpy; heat increments and the inverse temperature use the original
caloric owner. Its conservative slope
bound comes from the original caloric owner's source knots and includes a
signed heat capacity of mixing. Changing equilibrium compositions adds
nonnegative heat capacity for the symmetric pair law; negative fixed-phase
excess Cp must reduce the lower bound instead of retaining a pure-only CFL.
No second integrator or approximate latent-heat chart is introduced. This is
bulk, fixed-composition isobaric storage; phase-resolved transport must still
be refreshed through the original material adapter. The entropy capability
uses the SAME pure Cp/T integrals and accepted phase mixing entropy, divided
by the retained reference mass. Both enthalpy endpoints must be admitted;
a nonzero interval unresolved by the caloric inversion refuses. Existing
entropy-aware thermal advances can therefore audit changing liquid phases.
External source watts still do not define a reservoir entropy. This closed
storage constructor has no advection or caloric-identity binding; the explicit
stream constructor below supplies the common reference required for flow.

The same storage exposes equilibrium specific heat in J/(kg K), including
phase redistribution from the original regular-solution owner and both pure
Cp sources. It remains independent of the normalization inventory. A strict
Cp query refuses a phase-onset corner or unresolved phase branch instead of
claiming a unique derivative there.

`caloric::storage::conduction::BinaryLiquidConductionLaw` binds that exact
storage object to the existing implicit thermal chain/graph. One original
enthalpy inversion supplies T and the analytic inverse slope `mass/(dH/dT)`;
conductivity and its enthalpy derivative use the same temperature. At phase
corners the original caloric owner explicitly selects a homogeneous one-sided
Newton tangent; inside an unresolved numerical binodal bracket it selects that
same approximate tangent. Acceptance still uses complete equilibrium energy.
The bulk conductivity curve is caller-declared for the fixed composition,
pressure and geometry, with complete support over the storage domain. It is
not inferred from phase conductivities or treated as a resolved two-phase
fluid. No interphase kinetics, morphology law or phase-fraction error observable
is supplied. Synthetic regressions compare Cp
against independent high-precision values and enthalpy differences, and check
native implicit heat exchange across miscibility against an independent
energy balance with unequal cell masses and an interval beyond explicit CFL.

Focused source tests cover ordered molar conversion, original source retention,
unequal-temperature mixing with variable Cp, closed heating/cooling, exact
accepted phase identity, pressure/datum/source mismatches and atomic refusal.
Thermal timeline tests check prescribed heat through a phase-count transition;
the hydrostatic consumer checks independently calculated interfaces, buoyancy
and native motion from the resulting thermal endpoints. These are synthetic
numerical fixtures. They do not establish measured pair accuracy or a combined
thermal, gravitational and mechanical energy balance while interfaces move.

### Nonideal liquid streams, mixers and recirculating thermal flow

`BinaryLiquidCaloricLaw::stream_storage_from_mass_fractions` prepares the SAME
storage and conduction consumers on a common original pure-species enthalpy
datum. Each inlet carries its own actual sensible-plus-mixing enthalpy;
receiving cells recover temperature using their own retained composition and
original caloric owner. Composition-dependent heat of mixing, variable pure
Cp, signed excess Cp and equilibrium redistribution enter once. Independent
cold-endpoint zeroes cannot be used for energy exchanged between compositions.

Exact authored mass fractions bind directly to the existing stationary-species
solution. Conversion uses each original molar mass and refuses lost positive
traces. The common reference binds ordered species, masses, pure-source
identities, pressure and the complete original A/B/C/T0 interaction law and
domain. It excludes composition and numerical budgets; the separate storage
identity retains composition and temperature support. Inversion subtracts the
same cold-endpoint offset, and admitted endpoint temperatures are exact.
Unrepresentable common-datum ranges and unresolved interior offsets refuse.

Existing `Advected::with_stationary_species` and implicit chain/graph consumers
own mass and energy transport, including inlet-fed recirculation and simultaneous
conduction. This remains prescribed, stationary, mass-balanced bulk flow at
one fixed pressure. It does not evolve species inventories, select individual
liquid phases for transport, solve momentum or infer geometry. Original
fixed-composition storage entropy remains available.

Each stream also retains one `CaloricEntropyAnchor`, prepared from the original
pure Cp/T integrals and its actual cold equilibrium's mixing entropy. The
shared entropy-reference identity is independent of composition and prepared
temperature bracket. It removes only common linear pure-species entropy
constants, which cancel in the admitted stationary species balance.
`with_implicit_advection_entropy` uses the existing implicit chain/graph owner
to retain boundary entropy, signed composition/thermal mixing and numerical
temporal defect separately on the same accepted BE or SDIRK stages. This is
source-relative fixed-pressure bulk-flow accounting, without inferred absolute
species entropy, selective-phase work or a complete pressure/chemical-work
closure. Enabling the audit preserves the energy trajectory.

Native integration tests run an open implicit mixer at ten cell residence
times and an inlet-fed recirculating conduction graph against independent
coupled energy equations and analytic variable-Cp enthalpy inverses. They
cover unequal component masses, both signs of excess Cp, equal-temperature
heat of mixing, unequal-temperature feeds, conservative external energy,
exact composition binding, source/datum rejection and unchanged closed storage.
Additional native cases compare original entropy integrals with actual mixer
stage quadrature and recirculation/conduction balances, including distinct
cold brackets, signed excess Cp, equal-temperature mixing, source mismatch
and cancellation. Audited and energy-only cell states and energy reports are
bit-for-bit equal.


### Qualified molecular momentum moments

`thermochemical_gas::energy::parcel::transfer::species::effusion::momentum`
adds an explicit angle-independent, nondeflecting transmission qualification
bound to the original effusive interface. Scalar mass transmission alone does
not qualify the angular moment. The capability borrows already sampled beams,
retaining original source identity and the separately stored molecular energy
correction; no additional gas or material query is made.

An outward species has mass-flux-weighted normal speed
`sqrt(pi*R*T/(2*M_i))` and force moment `A*tau_i*p_i/2`. Accepted finite impulse
uses the original species mass actually withdrawn by the conservation owner.
For reciprocal streams the transmitted normal stress adds both directional
moments, including at zero net mass flow. It is not a wall-force difference or
an additional continuum pressure-thrust term. Collimating or accommodating
channels require a different angular law. Moving consumers must carry mass,
momentum and the Galilean energy transform together.

G0/G1/G3/G4 fixtures cover half-pressure moments, species selectivity, arbitrary
formation datums, reciprocal equilibrium stress, accepted finite withdrawals,
source identity, closed paths, cancellation and binary scale extremes. Native
execution remains unverified with the current missing constellation siblings.

### Prescribed molecular reservoirs with simultaneous vacuum escape

`thermochemical_gas::energy::parcel::transfer::species::network::transport::effusion::reservoir`
advances original internal paths, vacuum outlets and prescribed molecular
reservoir ports in one finite species/energy transaction. The external source
retains its original thermochemical and aperture qualification; its nominal
volume is not a finite donor budget. Its owner must accept the opposite
exported species and energy accounts. This scalar owner assumes a common rest
frame; a moving composition must supply the corresponding momentum and
inertial energy receipts.

Each reservoir aperture retains both sampled molecular beams, including at
zero net mass flow. Incoming species may extend the finite gas's original
catalog through the existing source-identity join. All gross outgoing species
spend the original finite inventory across every internal and external path;
simultaneous supply cannot fund an otherwise excessive withdrawal. All accepted
rows enter one original caloric reconstruction. Vacuum sampling and endpoint
qualification are shared with the existing vacuum-only owner.

The boundary energy is evaluated relative to each original species enthalpy
at the finite gas's initial temperature. A source-owned enthalpy increment and
the molecular `-R*T/(2*M_i)` correction preserve counterflow heat beside large
formation-energy datums. Incoming molecular energy includes its complete
relative kinetic and internal energy exactly once. No additional
thermalization heat, continuum pressure work, source depletion or closed-gas
entropy condition is invented. Actual endpoints are re-admitted against every
active original source and transmission; refusal and cancellation are atomic.

Five focused fixture groups cover reciprocal equilibrium, hot counterflow and
formation-datum independence, simultaneous supply/vacuum with shared gross
withdrawal, introduction of an original species, vacuum-only parity and
source/endpoint/identity/cancellation refusal. These are numerical checks with
synthetic source coefficients, not empirical validation. Native compilation
and execution remain unverified while required constellation siblings are
absent.

### Drifting Maxwellian incidence through qualified molecular apertures

`thermochemical_gas::energy::parcel::transfer::species::effusion::drift`
binds the existing shifted-Maxwellian half-space owner to original gas calorics,
species transmission and aperture area. `EffusiveDriftInterface` adds an authored
relative-speed support and provenance to the original angular qualification.
That declaration covers an undisturbed incident Maxwellian, exposed collisionless
transmission and the original temperature/pressure/catalog range, including
finite-thickness, grazing and counterstreaming-collision assumptions. Rest-frame
thinness alone does not qualify moving incidence.

`ClosedGasParcel::drifting_effusive_flux` receives source velocity minus aperture
velocity and an outward unit normal. Each species retains original source
enthalpy, its separate molecular crossing-energy correction, the complete
conditional vector velocity and the original zero-drift beam. Accepted finite
mass supplies both impulse and molecular energy. Area and original transmission
are applied once. Incoming relative kinetic energy is already included; no
additional thermalization heat or continuum pressure force is introduced.

Exactly zero drift preserves the original beam. Tangential-only drift preserves
its mass flux and adds tangential kinetic energy to the original thermal
coordinate. The direct positive gross incidence preserves rare outward-drift
thermal tails; a separately retained mass-rate change preserves small drift
beside rounded equilibrium beams. This source adapter reuses the existing
thermal kinetic kernel and does not implement a second Gaussian integrator.
Zero area is inactive; a calibrated zero transmission still admits its source
and drift support. Moving consumers must refresh the sampled relative velocity.

Six focused genuine-source groups cover rest parity, tangential/internal-mode
energy, opposing-halfspace conservation, positive tails and small increments,
qualification/cancellation, and rotation/area/caloric-reference invariance.
Rustfmt parsing and source/API checks passed. Native compilation and execution
remain unavailable because required constellation siblings are missing.

### Finite gas exchange with drifting molecular reservoirs

`reservoir::step_candidate_with_drift` composes vacuum outlets, comoving sources,
qualified drifting sources and existing internal gas paths on one original
species budget and one caloric reconstruction. The original comoving entry
delegates with an empty drift slice. A `EffusionDriftingReservoirStep` retains
distinct source-qualified incoming moments and the finite gas's original
outgoing vacuum beam. It never labels a shifted distribution as a rest-frame
reciprocal beam. The prescribed source's volume is not a finite donor budget.

All original accepted masses carry their own molecular energy. The common
caloric helper evaluates original species enthalpy increments, the rest
molecular correction and the drift selection change separately, preserving
heat beside large formation-energy references. Incoming drift energy is
included once. Gross outgoing streams spend the original finite inventory
across every path; simultaneous supply cannot fund excessive withdrawal.

Near zero drift, the authoritative signed mass retains the stationary net plus
the separately resolved drift increment. Far outward tails use direct gross
incoming mass, so baseline cancellation cannot close a real molecular path.
The existing boundary owner checks agreement of signed and gross accounts and
the resulting mass/energy residuals against its independent physical budgets.
Species may enter from zero using their original source catalog.

Initial samples and final gas endpoints retain original catalog, caloric,
transport and transmission qualification, including calibrated zero
transmissions. Closed geometric paths remain inactive. The enclosing mechanical
owner supplies the inertial transform and re-admits/resamples relative drift
as its aperture velocity changes; this material owner evolves no bulk momentum.

Five genuine-source groups cover exact comoving parity, tangential heating
under large reference shifts, shared three-port withdrawal, tiny signed drift
and nonzero thermal tails, and endpoint/source/cancellation refusal. Source and
fixture files passed Rustfmt parsing and bounded API review; native execution
remains blocked by the missing constellation checkouts.

### Reciprocal drifting streams on one shared molecular interface

`species::effusion::drift::exchange::sample` samples both original gases once
through the SAME qualified aperture in one declared translating control frame.
Its normal points left to right; the original incidence owner uses the opposite
normal for the left donor and that normal for the right donor. Separate gross
species masses retain their complete conditional vector impulses and original
crossing energies. No finite donor is replaced with a prescribed reservoir.

The sampler also retains one signed species transfer and the two caloric rows
relative to each gas's original partial enthalpy. A shared formation-independent
energy coordinate and original enthalpy increment construct both rows before
rounding; gross reference-energy reports are not subtracted for thermal closure.
Near-rest signed drift increments remain separate from rounded gross masses;
far-tail beams use their direct positive sampled incidence. Zero net mass does
not discard counterflow momentum or energy.

The consuming finite owner must admit BOTH original gross species withdrawals,
transform the same stream impulse into its common inertial frame, account for
actual bulk kinetic changes and publish both endpoints atomically. This sampler
does not evolve gas, add thermalization heat, infer remaining-wall forces or
supply flight/collection geometry between separated vessels. Its retained
`admit_endpoints` operation checks both final source/drift states and reciprocal
selected transmission support without sampling or spending a beam again.
Zero area is inactive; finite frozen rates require time refinement.

Focused original-source cases specify pressure moments at zero net mass,
thermal/tangential energy and datum crosswalks, tiny signed drift and positive
tails, orientation/scaling invariance and atomic source/cancellation refusal.
Rustfmt parsing and direct API review are available; native compilation and
execution remain blocked before Rust compilation by missing sibling manifests.

### Reusable original gas species inventory

`ClosedGasParcel::species_inventory_kg` exposes the SAME canonical active-species
inventory and representability checks used by the original finite material
owners. `species_mass_kg` delegates to it and retains exact zero for absent
species. Coupled molecular withdrawal checks and accepted-state balances query
each original/final gas inventory once; adaptive comparisons query each trial
state once. Per-species lookups no longer repeat full composition conversion
and allocation, and no independent mass-fraction model or approximate cache is
introduced. Existing finite-species and conservation cases exercise the shared
values; no runtime performance measurement is claimed.

## Original liquid sources in a spatial miscibility field

`liquid::solution::gradient` exposes `homogeneous_cell` and
`BinaryLiquidSquareGradientLaw` for finite-rate spatial liquid compositions.
Each cell retains the original homogeneous regular-solution state, including
metastable and spinodal states, and resolves the actual liquid properties
through the original card, claim selection and query context. A bulk
coexistence solve is not performed independently in each spatial cell.

The adapter binds the exact original pair to an explicit
`BinarySquareGradientV1`. Every resolved cell compares its source-derived
volume against `(n_A+n_B)/c`; the caller supplies a numerical relative
allowance in `[0,1e-6]`. An incompatible density source refuses instead of
being replaced by an invented mixing rule. Agreement at one composition does
not qualify equal partial molar volumes for an entire mixture. The gradient
coefficient and the physical range of this volume model remain explicit
constitutive inputs.

The existing equilibrium-phase resolver and homogeneous caloric feed now
share this same source-refresh path. Original component inventories, molar
masses, phase mass/volume and source receipts remain intact. A positive trace
that cannot survive the source's single first-species composition coordinate
refuses explicitly, including in a coexisting phase. The spatial cell wrapper
does not claim local equilibrium or change the homogeneous caloric feed's
relaxation model.

Focused integration tests cover homogeneous spinodal cells versus the
original bulk split, source identity and volume/mass consistency with unequal
molar masses, incompatible density/pair refusal, source-domain limits, traces,
empty inventory and cancellation. These material integration tests require
the normal material/runtime build; the dependency-light thermochemistry tests
do not substitute for that gate.

### Original-source binary liquid diffusivity

`liquid::solution::diffusion::MaterialBinaryDiffusivity` binds one explicit
ordered species pair to the original scalar Maxwell–Stefan diffusivity claim
in m²/s. Absolute temperature and first-species mole fraction require their
actual static quantity declarations. Typed external conditions, including
pressure, remain fixed at the original query; stress, strain, rate and field
coordinates cannot acquire a frozen physical owner through this adapter.

Each evaluation resolves the left endpoint, arithmetic-mean composition face
and right endpoint through the shared material resolver and the same pinned
claim. It evaluates nonlinear source data at the face rather than averaging
endpoint diffusivities. Joint interpolation support, numerical nonnegativity
and source validity apply to all three samples. Explicit zero diffusivity is
admitted as a closed kinetic path; negative values, unsupported points and
unrepresentable traces refuse. The owner neither extrapolates nor switches
claims when temperature or composition changes.

The adapter supplies stage-lagged kinetics to a transport owner; it does not
solve transport, infer a thermodynamic factor from diffusivity, or certify an
entire path from endpoint samples. Query and complete receipt allocations are
explicit; no mutable source cache is retained. Focused source tests exercise
temperature/composition changes, fixed pressure, original claim identity,
endpoint and joint-face support, zero/negative coefficients, typed-coordinate
refusals and cancellation without partial publication.

## Coherent reflection and mode conversion between elastic solids

`elastic::acoustic::interface::IsotropicElasticInterface` borrows two original explicitly isotropic `ElasticAcousticState` values. Their Christoffel speeds, density, Mandel stiffness, identities and original strain ceilings remain authoritative. `prepare_incidence` fixes conserved tangential slowness from one propagating P/SV/SH channel and an incidence cosine; one row-scaled native 4×4 complex LU supplies P/SV scattering from both sides, while SH uses its independent traction match. Preparation is cancellable and retains fixed-size wave and scattering arrays.

The prepared operator applies six coherent incoming modal velocity amplitudes without allocation or another factorization. It returns each outgoing field, modal normal powers, shared-interface velocity/traction and independently reconstructed mechanical flux. Original tensor contraction preserves weak grazing traction; the near-critical normal-wave calculation retains the original speed pair and selectively uses double-double arithmetic. Evanescent outgoing waves retain complex phase and stress but carry zero individual normal power. Nonzero evanescent/grazing incoming excitation refuses. An exactly identical mechanical medium transmits without a fictitious interface, including undriven grazing modes.

Peak phasors use `exp(i omega(s.x-t))`; +z points first to second and traction is `sigma.e_z` on both sides. Finite drives must satisfy a conservative sum of modal Mandel-strain norms within each original componentwise source ceiling. This is a floating-point model check, not an interval certificate. No power clipping, attenuation, anisotropic projection, contact compliance, finite-layer resonance, head-wave or finite-source approximation is introduced. The scope is an infinite, stationary, perfectly bonded interface between frozen lossless isotropic half-spaces; static moduli do not acquire measured ultrasonic authority. Singular or unrepresentable boundary systems refuse.

Ten focused G0/G1/G3/G4 regressions check independent normal-incidence impedance formulas, original-tensor field/traction continuity, coherent two-sided power, weighted unitarity and reciprocal conversion, supercritical phase and evanescence, critical branches, exact transparency, auxetic/high-contrast media, common material scaling, zero-Lame grazing traction, original domain refusal and cancellation. Native actual-source execution passed with the original elastic constructors, Jacobi, Cholesky, complex LU and scalar GEMM; the narrow wiring omits unrelated source-selection/constitutive modules. Full Cargo/workspace verification was unavailable.

## Source-bound planar binary liquid interfaces

`BinaryLiquidSquareGradientLaw::equilibrium_interface` joins the original
finite-inventory equilibrium, its actual coexisting material phases and the
same-pair planar interface calculation. The private-construction
`ResolvedBinaryLiquidPlanarInterface` exposes both retained owners to
capillary consumers. A one-phase inventory returns None even when a different
overall composition of that pair could split.

Each coexisting phase must pass the same source-density/fixed-site-volume
check as a spatial cell. The check reuses resolved source properties rather
than evaluating or replacing them again. Original species inventories,
temperature, card/claim selection and query context remain intact. The
operation has at most two bounded coexistence solves, two source refreshes
and the explicit quadrature work allowance; failure or cancellation publishes
no partially bound interface. It supplies neither a substrate wetting law nor
a claim that a diffuse interface is thin enough for a particular flow geometry.

Integration regressions retain actual mass/volume and source identity, compare
the bound interface with its original core calculation, reject source-domain
and density incompatibility, and preserve one-phase/cancellation behavior.
These material integrations require the normal runtime build.

## Source-backed finite-band analytic Sellmeier

`SpectralSellmeierStatePoint` retains six original SI oscillator coefficients,
two explicitly sourced fit limits in vacuum metres, and nine independently
sourced intensity-absorption coordinates and values. All 26 claims, typed
quantities, receipts and complete material state remain attached. The fit band
and absorption support remain separate; consumers admit the complete fit
before using their common domain. Refresh pins the original claims. A card
with either explicit fit-limit property must use this resolver; the legacy
visible resolver refuses rather than discarding that domain, even with pins.

Explicit spectral tabulated indices admit ordinary real `n >= 1` with
`n² <= f64::MAX/64` arithmetic headroom. Legacy visible source bounds remain.
These source schemas provide neither missing transparency nor a measured
material preset, and do not certify extrapolation beyond original support.

## Finite elastic plate refraction and transmission

`ElasticPlateLayers` prepares the original incident thermoviscous state,
existing equivalent-fluid entrance layers, one explicitly isotropic lossless
elastic state and a borrowed receiving fluid. Density and original Christoffel
speeds are retained. `scattering(cosine)` returns coherent pressure R/T,
plate-front pressure, outgoing P/S and receiving normal wavenumbers, and normal
reflected/transmitted/porous-absorbed power fractions. Incidence queries allocate
nothing; elastic plate work has fixed size, plus the original entrance stack.

Both plate faces enforce continuous normal traction/velocity and zero shear
traction. Even/odd P/S potentials retain finite-thickness mode conversion,
interference, Lamb resonances and evanescent tunneling. Bounded sin/cos/sinc
pairs and a factored logarithmic transmission determinant avoid growing
evanescent transfer entries. Critical angles, exact zero thickness and the
normal P-only limit remain finite; the latter reuses the original finite-fluid
transfer owner. Receiving evanescence carries zero normal power, and elastic
evanescence is never heat. Only original porous entrance loss is absorbed.

The scope is one homogeneous, frozen, unstressed, lossless isotropic plate.
Original state provenance, response regime and strain admission remain intact;
no anisotropic projection, solid damping, lateral finite size or reconstructed
interior elastic field is supplied. Nine new native controls and 19 existing
controls pass, including an independent six-boundary-equation oracle, normal
and zero-thickness limits, critical neighborhoods, reciprocity, porous loss,
thick evanescence and a narrow air-loaded flexural tunneling resonance.

## Explicit broadband molecular dispersion

`MolecularCauchyDispersion::try_new_in_band` admits an explicitly supplied
normal-dispersion molecular fit on its original positive finite vacuum-wavelength
interval, including ultraviolet and infrared support. The original reference
wavelength and number density remain attached, and the reference must lie inside
the fit interval. Monotone endpoint checks retain the dilute refractivity and
representable-derivative requirements. The existing `try_new` and unmodified
constant-polarizability sources retain their visible-domain behavior.

After every component has an explicit fit, `MolecularRefractionModel` uses the
intersection of those original source bands. Analytic phase/group indices,
spatial composition derivatives and independent molecular cross sections share
one wavelength preparation. A zero-shape explicit broadband source still uses
its molecular volume for nonvisible scattering; visible constant-source arithmetic
is retained. The original wavelength-independent depolarization and nonabsorbing,
dilute independent-molecule assumptions remain. No measured dataset, resonance
model, gas absorption or unbounded extrapolation is inferred.

`MolecularRefractionWavelength::refractive_index_jet` also supplies analytic
first and second wavelength derivatives at fixed original partial number
densities. It differentiates each original Cauchy law through the full nonlinear
Lorentz–Lorenz mixture, retaining weak dispersion even when the phase index
rounds to unity. The scalar and spatial refraction paths do not pay for this
additional jet. Constant-polarizability and vacuum limits have zero curvature.

## Atmospheric-air phase, group and state refraction

`air_refraction` owns the historical NIST Toolbox Ciddor equations A21–A41.
Inputs explicitly retain water-vapor mole fraction, dry-air CO2, absolute
pressure, temperature and vacuum wavelength. Admitted support is 300–1700 nm,
233.15–373.15 K, 60–120 kPa, water fraction 0–0.2 and CO2 0–2000 µmol/mol.
This numerical envelope is not phase-stability or uncertainty evidence; the
state owner must exclude condensation, droplets and incompatible constituents.

Stack-owned wavelength preparation supplies phase and group indices, spectral
curvature, and complete pressure/temperature first and second derivatives at
fixed composition. It retains the source compressibility law and never changes
the fluid equation of state. Queries allocate nothing and never extrapolate.
`AirRefractiveSample::absolute_from_relative` applies the phase/group product
rule when the relative source uses the same vacuum-wavelength coordinate;
missing relative dispersion remains unavailable. Absorption, scattering and
moving-medium effects require separate sources.

Focused controls compare 14 published NIST dry/moist-air table entries,
independent spectral/state differences, source-domain corners and the relative
index conversion. These validate this implementation of the named formulation,
not a new experimental accuracy claim.

`PreparedCiddorAirSpectrum` prepares the original dry-air/vapor spectral poles
once independently of composition. Its `sample_jet` supplies the complete phase
gradient and Hessian in [temperature K, absolute pressure Pa, water mole fraction,
dry-air CO2 micromol/mol], including all six mixed partials. The original
compressibility closure and fixed-composition arithmetic remain shared.
`ciddor_model_identity` identifies the formulation independently of a particular
atmospheric snapshot. Production derivatives are analytic and allocate nothing;
the scalar/frozen owner does not evaluate the added composition derivatives.

Focused controls differentiate independent NIST density-ratio equations across
three states and three wavelengths, check all state-jet entries and source-domain
corners, and retain all ten frozen output bit patterns at three captured states.
Spatial providers still own their gradients, Hessians and phase validity;
humidity is a mole fraction, and no fluid EOS or phase-equilibrium closure is
introduced.

## Coherent fluid-coupled elastic plate stacks

`ElasticPlateStack` composes up to 64 ordered original isotropic elastic plates
and existing equivalent-fluid sections at one frequency. Every nonzero plate
has independent shear-free faces; positive fluid or porous-fluid separation is
required between nonzero solids. Original density, Christoffel speeds, material
identity, response regime and strain ceiling remain attached. One original
incident/saturating fluid snapshot is shared, with distinct fluid sources reused
across gaps and receiver.

The same tangential wavenumber and bounded pressure/normal-velocity pair cross
every section. Shared P/S parity propagation, original normalized fluid transfer
and compensated logarithmic transmission retain mode conversion, coherent
plate/cavity resonances, critical angles and evanescent tunneling. Pressure or
velocity nodes and exact transmission zero do not become admittance poles.
Preparation allocates bounded storage; incidence queries allocate nothing.

Solid plates and receiving fluid are lossless; absorption belongs only to original
porous sections. Bonded laminates, anisotropy, solid damping, lateral boundaries,
finite-drive admission and reconstructed interior fields remain outside this
owner. Thirty-eight native controls pass, including an independent 12-boundary-
equation double-plate oracle over 128 configurations, coherent cavity tunneling,
reciprocity, critical/zero/thick limits and original single-plate regressions.

## Coherent reflection through a finite bonded elastic slab

`IsotropicElasticLayer` binds three original explicitly isotropic, unstressed, lossless elastic states and a finite nonnegative layer thickness. `PreparedElasticLayer` retains the six exterior P/SV/SH scattering channels at one positive frequency and conserved tangential slowness, with each incoming wave referenced to its physical face. Both boundaries enforce all velocity and traction components. P/SV share one row-scaled 8 by 8 factorization and SH one 4 by 4 factorization; at most two residual corrections reuse those factors when a critical traction row is small.

Opposite-face anchoring keeps evanescent exponentials decaying. A regular even/odd cosine and sinc basis near internal criticality retains the second independent solution at zero normal wavenumber. Interior velocity, original-tensor traction and Mandel strain are available at depth. Total mechanical flux retains evanescent interference power; evanescence is not heat. A conservative whole-layer strain envelope admits the original material domain even between interface nodes. Zero thickness reduces to the original interface without an interior field. Identical mechanical media retain transparent propagation. Repeated drives and depth queries use fixed storage without refactorization or allocation.

Eight independent slab regressions and ten existing interface controls pass in the focused original-source numerical build. They cover independent Fabry–Perot coefficients, homogeneous and absent-layer limits, original-tensor traction and displacement gradients, coherent power, weighted unitarity and reciprocity, an analytic evanescent barrier across more than fifty decades of transmitted power, adjacent critical-angle branches, hidden interior strain antinodes and cancellation. The build uses the existing narrow validation facade and is not full-workspace or experimental qualification. The model is one infinite planar perfectly bonded slab; attenuation, prestress, anisotropic states, compliant contact, finite edges, arbitrary layer stacks and initial free guided modes are not supplied.

### Finite elastic coatings terminated by physical supports

`elastic::acoustic::interface::layer::supported::IsotropicElasticCoating` binds
the original incident solid and coating to a free, clamped, scalar-impedance or
reciprocal tensor-impedance rear support. The coating retains all P/SV/SH modes
and the same original-tensor critical/decaying basis as the two-sided slab.
One row-scaled native 12 by 12 solve retains reflected amplitudes, six coating
coordinates and the actual rear velocity; cached maps make coherent drives and
depth sampling allocation-free. No fictitious exterior incident evanescence or
scalar effective coating substitutes for the full boundary problem.

Both original strain domains remain authoritative. Coherent front, interior
and rear fluxes are checked independently against positive support work. The
coating itself is lossless: evanescence and standing waves do not create heat.
Explicit rear support velocity/traction retain weak motion beside stiff loads;
the unmodified modal field reconstruction can have lower relative accuracy at
such an endpoint. Exactly zero thickness delegates to the original surface
and has no coating interior. Source isotropy, small strain, stationary linear
support, infinite lateral extent and perfect front bonding remain the domain.

Six focused regression definitions cover exact absence, analytic impedance
transformation and quarter-wave matching, homogeneous tensor phase, an
independent 70-digit oblique critical-mode solution, critical/evanescent SH
transport, stiff support work, original domains and cancellation. Native
compilation and execution status are reported separately.

## Coherent two-sided fluid and elastic interface

`FluidElasticInterface` borrows the original `FluidAcoustics` owner and explicitly isotropic elastic state. The fluid's coherent acoustic coefficients are evaluated once and retained with the original owner; no material-name or ambient-air fallback is selected. `PreparedFluidElasticInterface` admits fluid P and solid P/SV/SH excitation from either side, including simultaneous coherent illumination at one tangential slowness. A row-scaled native 3 by 3 complex factorization enforces normal velocity, normal traction and zero solid shear traction. SH reflects independently. Tangential fluid and solid velocities may slip.

Outgoing evanescence retains its complex pressure and mechanical fields with zero individual radiated normal power. Incoming evanescent or grazing channels cannot be driven. Repeated drives reuse fixed-size coefficients without allocation. Original tensor strain bounds apply throughout the solid half-space; explicit caller-declared `FluidDriveLimits` bound fluid relative-density and particle-Mach amplitudes, including evanescent velocity amplification. These declared limits do not establish experimental validity of the fluid closure. Each response retains the original solid identity, fluid snapshot, coherent channel amplitudes, boundary fields, and independently reconstructed fluid pressure-work and solid traction-work fluxes.

Six independent regressions pass alongside all eighteen solid interface/slab controls in the focused original-source build: impedance formulas, an independent oblique Lamé/stress oracle, coherent unitarity and reciprocity, original-tensor continuity, receiving-fluid evanescence, adjacent critical branches, density scaling, original source retention, amplitude limits and cancellation. This uses the existing narrow material validation facade and original numerical kernels, not full-workspace or experimental qualification. The interface is infinite, planar, lossless, quiescent and small-signal with an unstressed isotropic solid; it does not add viscosity, thermoviscous boundary layers, anisotropy, finite apertures, radiation-reaction source evolution or cavitation. Singular or unrepresentable contact systems refuse explicitly.

## Original IAPWS ordinary-water refraction

`water_refraction::IapwsWaterRefraction` evaluates IAPWS R9-97 equation (1) and its released coefficients from original ITS-90 temperature, mass density and vacuum wavelength. It admits only 261.15–773.15 K, 0–1060 kg/m³ and 200–1100 nm. Fixed-size wavelength preparation and allocation-free state queries retain absolute phase index, stable dilute refractivity, analytic group index, and first/second temperature and density partials, including the mixed derivative. `sample_if97_state` consumes the original homogeneous IF97 state's T/rho without replacing its EOS or claiming IAPWS-95 pressure precision.

The caller retains ordinary-water purity, homogeneous phase validity and accepted state identity. The law supplies no absorption, scattering, ice/brine or two-phase mixture closure, and does not extend the published optical uncertainty. Four focused actual-source controls cover 24 published refractive-index values using independently printed IAPWS-95 density tables with rounding allowances, 20 derivative states, four original IF97 states including the complete isobaric chain rule, source-domain refusals, and vacuum/dilute limits.

## Declared WLF and Arrhenius viscous temperature shifts

`visco::temperature::TemperatureShift` retains caller-supplied WLF or molar
Arrhenius coefficients, the reference temperature and an inclusive calibrated
temperature interval. Its convention is `a_T = tau(T)/tau(T_ref)`: WLF uses
`log10(a_T) = -C1 (T-T_ref)/(C2+T-T_ref)`; Arrhenius uses
`ln(a_T) = (E/R) (1/T-1/T_ref)`. Temperatures are absolute kelvin and E is
J/mol. The reference shift is exactly one. The complete WLF interval must lie
strictly above its pole. Out-of-band temperatures and unrepresentable shifts or
logarithmic derivatives refuse. No universal material coefficients are inferred.

The existing tensor, scalar Maxwell, fractional Zener, modified Zener and lowered
Prony owners consume the same shift. Maxwell times and complete Kelvin viscosity
tensors scale by `a_T`; Kelvin factors scale by `sqrt(a_T)`. Elastic moduli,
anisotropy, branch order, elastic factors and strain limits remain unchanged.
No eigensolve, NNLS refit or second time integrator is introduced. Existing
harmonic and history operators therefore evaluate
`C_T*(omega)=C_ref*(a_T omega)` and reduced time `t/a_T`.
A lowered model retains its original fit measurements on the correspondingly
mapped frequency grid; it gains no experimental or inter-sample authority.
Bounded source frequency support shifts by `1/a_T`, and nonzero coefficients
may not silently vanish by underflow.

The returned `ShiftedViscoelastic<M>` retains the actual temperature and complete
shift declaration alongside the original constitutive implementation.
`temperature_shift.rs` and `temperature_shift_scalar.rs` exercise the formulas,
reference fixed point, refusals, harmonic response, plane stress, reduced-time
histories, positive viscous heat and transformed source support.
These are authored numerical fixtures, not measured material calibration.
The caller must qualify the same material as thermorheologically simple.
Moisture shifts, vertical modulus changes, thermal expansion, integration along
a changing temperature history, and caloric/entropy coupling remain separate.
This is a temperature-dependent constitutive slice of IRP M0.2-W4
(`frankensim-irp-v5-instrument-realism-pexfp.24`), not closure of that whole bead.

## Temperature shifts of finite-strain viscous kinetics

`TemperatureShift::shift_hyperelastic` applies the shared calibrated WLF or
Arrhenius factor to the original `ViscoHyperelastic` branch relaxation times.
It retains equilibrium and branch energies, ordering, deformation limits and
local Newton controls. Original public parameters are checked by the existing
constructors. Shift and original material errors remain separately typed.

The original full and bulk-free histories, mixed pressure formulation, return
map and consistent tangents remain the numerical owners. Accepted viscous
metrics and work/storage/dissipation ledgers carry between temperature intervals.
At frozen temperature, the same update corresponds to reference duration
`dt / a_T`; temperature alone does not change stored energy or reset memory.
No eigensolve, fit or evolution is performed when preparing the shifted law.

This is thermorheologically simple kinetic dependence with unchanged elastic
energies. It does not supply expansion, entropy, caloric coupling or exact
integration through a changing-temperature interval. Such trajectories require
explicit temporal refinement, and the caller retains the original calibration
and source qualifications. Shift the original reference law on each preparation.

## Calibrated temperature kinetics in the original thermovisco and heat owner

`ThermoMaxwell::with_temperature_shift` binds the shared WLF/Arrhenius
declaration to an original constant relaxation-time multiplier. The multiplier
remains the reference-time normalization; varying tabulated time scaling and
second shifts refuse. Calibration must cover the complete mechanical temperature
interval, contain the original reference temperature there, and leave every
scaled branch time finite and positive. The analytic shift is evaluated at each
actual mechanical/thermal trial temperature, without resampling a nonlinear
curve or re-evaluating exponentials for elastic storage.

`ResolvedThermoMaxwell::with_temperature_shift_from_source` binds explicit
reference-temperature/C1/C2 or molar-activation-energy keys from the same immutable
material card and fixed query context. Original claims are pinned, their full
temperature-knot union and the reference-temperature sample are retained, and
varying calibration coefficients or incomplete support refuse. Absolute kelvin
and temperature-difference roles remain distinct when semantic quantities are
provided; Arrhenius requires J/mol. The original mechanical samples and role
keys remain intact; `temperature_shift_source()` exposes the calibration
declaration and receipts.

Existing backward-Euler and exact strain-ramp history/tangent paths consume
these kinetics. Caloric feedback recomputes the same kinetics at the returned
thermal temperature and deposits only physical dashpot loss through the existing
entropy/internal-energy owner. Relaxation-time changes introduce no extra
thermal work, branch resets, or second heat source. Exact strain ramps retain
endpoint-frozen temperature coefficients and do not claim exact continuously
varying-temperature integration; thermal refinement is still required. This is
the executable temperature slice of M0.2-W4, without inferred moisture response
or empirical material qualification.

Seven focused regressions cover analytic off-grid kinetics, time-rescaling
equivalence, mixed heating/cooling history and energy, source/quantity/support
ownership, refusals, and implicit self-heating for both integrators.

## Source-owned nonlinear convection on integrated material layers

`ReservoirLayerGeometry::evaluate_convective_material_to_bath` composes an
original nonlinear film callback with the same supported conductivity chart
and analytic Kirchhoff inverse. The callback supplies signed surface-to-bulk
heat, an explicit film conductance and its retained source payload at each
actual surface temperature. One bounded root uses nonnegative solid power
along the physical heat direction. The integral-coordinate solid power is
retained directly, avoiding amplification of rounded surface temperature in
thin, highly conducting layers.

Only the actual material profile requires conductivity support; a bulk state
outside that chart is admitted when a supported passive surface balance exists.
The exact source-limit face closes the bracket. Acceptance requires the
declared absolute/relative power residual, not a small temperature or power
interval. Safeguarded secant steps reuse bracket residuals and fall back to
bisection. No derivative, nested root, extra thermal inventory or global
uniqueness certificate is introduced.

Actual callback heat is authoritative. Its cell-to-bulk secant lowers once to
the original reciprocal heat owner. At exact equilibrium the local original
conductivity slope is in series with the declared film conductance; an explicit
zero film conductance remains zero. Driven zero/active films, unsupported roots
and unrepresentable inverses refuse. Original callback errors, thermal
failures, cancellation and exhausted surface balances remain typed. Each solve
uses at most the authored number of interior trials plus two bracket film
calls, with cancellation checks around every call.

Six focused G0/G1/G4 regressions cover constant-film parity, manufactured
nonlinear film/k(T) profiles, both heat directions and equilibrium,
source-limited roots, thin-metal flux accuracy, and typed
source/cancellation/budget refusals.
All six and eight existing reservoir tests passed in a bounded build of the
original material modules with the native thermal error and series helpers.
This does not establish full Cargo or graph/cavity integration execution.

### Two passive baths on one conducting surface

`ReservoirLayerGeometry::evaluate_material_to_two_baths` closes one supported
conductivity integral against two distinct boundary laws on the same surface.
It reuses the original safeguarded power iteration with a signed material-power
bracket. Each callback supplies actual surface-to-bath heat, its nonnegative
secant and retained source payload. Zero coefficient and power explicitly
insulate a bath. Opposing powers and exact zero net material heat remain
distinct, including heat relayed between baths when material and one bath have
equal temperatures. There is no substitute weighted ambient temperature.

The physical material/bath temperature hull is intersected with the conductivity
chart and an explicit original boundary-source surface interval. Only the
actual material path and boundary trials require source support; a bath never
becomes a substitute material or gas query. Acceptance uses the declared power
residual, scaled by the largest individual material/bath power. Original
callback errors, source-domain refusal, cancellation and exhausted work remain
typed. At most two bracket callbacks supplement the bounded interior trials.

Returned powers remain authoritative for separate energy accounts. Three
nonnegative pair conductances describe the frozen surface elimination and
capacity rows, including the bath-to-bath path. They use the material integral
mean conductivity and normalized products, preserving thin-wall power and
avoiding avoidable overflow. They are not derivatives of changing source
coefficients. This adds no heat storage, second root, uniqueness certificate or
new physical model-error allowance.

Five direct regressions cover independent three-reservoir rates, opposing
drives and exact neutral relay, insulation and one-bath parity, nonlinear
convection/radiation with affine k(T), source-limited surfaces, thin-metal
accuracy, source errors, cancellation and work-budget refusal. All five plus
the previous eighteen material/rectangular controls passed in the existing
bounded original-source build (23 total). Full source ownership and native
graph/LBM compilation are outside that numerical lane.

## Source-integrated implicit reservoir and ambient fluxes

`thermal::conduction::network::ImplicitLayerBoundary` borrows an original
`ResolvedConductionLayer` and prepares one planar centre-to-face geometry and
areal film resistance. It uses the existing material-to-bath Kirchhoff integral
and bounded face solve at every native implicit trial. The secant conductance
and complete outward dQ/dh share the same source and caloric branch, including
the admitted left derivative at an upper source endpoint. No source query,
allocation or separate heat-capacity account occurs during evaluation.

The original chain or graph supplies each actual reservoir or ambient bath
temperature, owns accepted enthalpy, heat and entropy, and splits authored bath
events at their exact times. A cell may have both separately covered channels.
Original area, distance, resistance and caloric object must match; unselected
faces retain their native coefficients. A bath outside conductivity support is
admitted only through a finite film whose solved material face stays supported.
Zero film requires the prescribed material face itself to remain supported.

Focused source-network tests cover affine-k analytic backward Euler, both bath
channels, simultaneous distinct baths, finite-film temperatures and tangents,
source boundaries, exact events, cancellation and deterministic retry. These
are numerical implementation checks. The half-cell reconstruction is
quasi-steady and pure-phase; it does not resolve internal face storage, fluid
boundary layers, multidimensional spreading or a moving phase interface.

## Calibrated thermal shear-thinning viscosity

`fluid_transport::rheology` owns a checked SI Carreau-Yasuda master curve
(`0 < n <= 1`, `a > 0`) with finite declared reference-rate support. It retains
the existing WLF/Arrhenius `TemperatureShift` and explicitly distinguishes
`mu(T,g) = a_T * mu_ref(g)` from `mu(T,g) = a_T * mu_ref(a_T*g)`; the latter
checks support at the reduced reference rate. Stable log evaluation preserves
representable extreme responses, analytic shear/temperature tangents and
positive differential stress. Singular or unrepresentable scalar shear
derivatives are explicit `None`; stress remains regular at rest.

Binding queries the supplied `FluidTransport` once, retains the exact original
source/snapshot, and replaces dynamic viscosity once while preserving current
density, heat capacity and conductivity. Since this trait has no temperature
coordinate, supplying a source at the evaluated temperature remains the caller's
obligation. Source refusals propagate unchanged. No calibration authority,
elastic memory, yield surface, phase law, extensional-flow fit or extrapolation
is inferred. Ordinary constitutive queries allocate nothing and parameter logs
are prepared once.

`ThermalRheologyState::bind_thermal` additionally admits the source's original
`ThermalFluidTransport` snapshot. Exact source temperature and agreement of its
thermal/momentum tuples are checked before returning `ThermalRheologicalTransport`.
Both public views then share the same frozen, viscosity-replaced coefficients
without repeated source queries. Signed, zero and unknown source expansion
remain distinct; temperature is not reconstructed from the viscosity law.
This enables source-owned thermal buoyancy and conservative material rebinding
without losing the calibrated rheology. Flow evolution must still resolve its
actual shear/temperature; the adapter does not qualify Newtonian empirical
correlations for non-Newtonian fluids. G0/G3 controls check mixed-temperature
and incoherent-source refusals, retained owners, one-time queries, and the
coherent change of Prandtl/Rayleigh groups under viscosity replacement.

`fluid_transport::rheology::pipe::LaminarCarreauPipe` resolves fully developed,
no-slip circular-pipe flow of that original law at one frozen temperature.
Both rest and the explicit physical wall-rate ceiling must lie inside source
support, including the original reduced-rate shift. A bounded bracketed Newton
solve inverts the positive stress law `tau=mu(g)*g`. Total volume flow uses a
positive shear-rate-coordinate integral of the canonical stress tangent, so it
needs one wall inversion rather than an inverse at every integration node.
Constant-viscosity declarations recover Poiseuille flow directly. Nested
positive Simpson quadrature splits at the original transition rate and requires
two consecutive refinements to agree within the declared tolerance; exhausted
inverse or quadrature budgets refuse. Observed quadrature changes are numerical
estimates, not certified error bounds or material uncertainty.

The prepared radial profile retains the original pipe/law and controls. It
resolves axial velocity, outward-face radial shear stress, physical velocity
gradient, local viscosity and nonnegative `tau*g` heating. Factored wall distance
preserves small positive velocities near the no-slip wall. Accepted pressure is
derived from the actual wall constitutive response; the requested pressure and
signed inversion residual remain explicit. Consequently the integrated flow and
`accepted_pressure*flow` dissipation describe the same accepted state. Heat is
reported for a thermal owner, not deposited or used to change the frozen law.
No density, Reynolds-number criterion, turbulence qualification, entrance loss,
wall slip, elastic memory or pressure-dependent source response is inferred from
the caller's laminar rate ceiling. Focused regressions cover independent radial
velocity/heat integration, both Newtonian plateaus and constant-law limits,
pressure reversal/rest, singular scalar zero-rate slopes, shifted source limits,
near-wall velocities, finite work budgets and extreme intermediate products.

`LaminarCarreauPipe::profile_at_flow` additionally resolves the required
pressure for a prescribed signed volume flow through that same radial law.
A bounded Newton solve acts directly on wall shear rate, with the original
source ceiling, positive quadrature and analytic Rabinowitsch derivative.
It does not nest constitutive inversions, finite-difference pressure trials,
or replace the radial viscosity by a single apparent value. The returned
solution retains requested flow, accepted flow residual and actual iteration
count; pressure, profile and irreversible power use the accepted state.
The existing temperature calibration and both shift conventions are unchanged.

Each accepted profile exposes the positive differential pressure conductance,
computed from its retained radial integral without a new solve or quadrature.
This derivative is even under reversal, differs from the secant in thinning
flow and recovers Poiseuille conductance at rest. Flow residual tolerance,
quadrature estimate and source support remain separate; neither Newton
convergence nor a declared wall-rate ceiling establishes laminar validity or
material accuracy. Focused G0/G1 controls cover analytic Newtonian pressure
and heat, an independently inverted nonlinear stress law, pressure-derivative
finite differences, signed flow, source ceilings and bounded-work refusal.

`LaminarCarreauPipe::profile_at_wall_rate` admits an authored signed wall rate
directly through the original prepared source. The exact declared ceiling can
therefore supply an operating-envelope profile without an external pressure
round trip rounding beyond its supported stress. Original radial quadrature,
temperature, flow, pressure and positive heat remain unchanged.

`fluid_transport::rheology::pipe::transient::RadialCarreauPipe` resolves actual
nonlinear radial viscous startup, relaxation and reversal. It evolves the
axially developed equation `rho*u_t = delta_p/L + (r*tau)_r/r`, with
`tau=mu(|u_r|)*u_r` evaluated by the same immutable original Carreau preparation.
Uniform radial finite volumes retain every cell's axial velocity and exact
annular mass. No frozen apparent viscosity, constant mean-flow inertance or
fitted Newtonian memory spectrum replaces the changing radial stress field.

Implicit midpoint and shared signed face stresses give a discrete pressure-work
balance against full profile kinetic-energy change and nonnegative viscous heat.
The reported time-average flow is midpoint quadrature, not an exact continuous
trajectory integral. Nonlinear momentum residual work and the observed balance
defect remain separate from physical heat. Both endpoint shear fields must stay
inside the original source/rate support; their linear time interpolation cannot
have a larger face-rate magnitude. Temperature is held during each interval;
density and geometry stay fixed.
Density is an explicit physical input, not new material calibration authority.

Analytic original stress tangents form a tridiagonal Newton system, solved by
the existing native `fs-la::tridiagonal::TridiagonalWorkspace`. The accepted
same tangent system also returns the positive derivative of midpoint flow with
respect to pressure at fixed initial history and duration. Bounded backtracking
does not relax source ceilings. Cold construction allocates O(radial cells)
state and scratch; warm steps reuse them without allocation. Local material
queries and fixed-size radial tiles have bounded cancellation checks, and the
native structured solve retains its own checks. A final checked swap publishes
history and time together; source refusal, exhausted work or late cancellation
leaves the accepted state unchanged and the workspace exactly retryable.

`RadialCarreauPipe::step_from` evaluates an immutable initial radial history
into a separate reusable output state using the same candidate solve as
`step`. Pressure trials therefore retain their actual initial velocity field
without allocating or advancing accepted history. Model identity, source,
nonlinear work and final cancellation checks precede the output swap; refusal
leaves that output unchanged. The `fs_airflow::graph::rheology::transient`
consumer uses this operation for coupled nodal pressure and linear storage
solves, while this material owner remains a single original radial pipe.

`LaminarCarreauPipe::at_temperature` prepares the same original calibration at
a prescribed temperature without changing geometry or the physical rate ceiling.
Both rest and that unchanged ceiling must remain inside the original reduced-rate
support. `RadialCarreauPipe::at_temperature` carries this preparation into the
same radial discretization, preserving every annular mass and its compatibility
identity. Existing states and scratch/output buffers remain usable without
copying or resetting velocities or time; kinetic storage is unchanged by the
viscosity switch. The diffusion timescale is recomputed from the actual shifted
zero-rate viscosity. This no-allocation radial rebind supplies a uniform held
temperature, not heat evolution, expansion or density feedback. The caller must
retain the original fixed-density approximation and separately account for
thermal energy. G1/G3 tests cover the exact discrete Arrhenius diffusion mode,
both nonlinear shift conventions, retained energy/history, warm workspace reuse
and refusal of unsupported temperature or reduced-rate ceilings.

This is a spatially and temporally convergent discretization, not an exact
propagator. Midpoint is second order but under-resolved stiff modes can oscillate;
energy balance alone does not establish time accuracy or monotonicity. The
explicit radial mesh, diffusion timescale, nonlinear tolerances and iteration,
backtracking and material-evaluation budgets remain visible. Spatial and time
refinement must resolve the actual forcing. The model assumes incompressible,
laminar, axially developed flow with stationary no-slip walls. It adds no axial
entrance/convection, compressibility, cavitation, thermal evolution, transition
prediction or network/junction dynamics, and does not deposit its heat. Focused
regressions independently compare Newtonian Bessel startup/decay across three
radial resolutions, measure temporal order on an exact discrete eigenmode,
recover the original nonlinear steady profile, differentiate the transient
pressure response, retain unequal radial histories at equal mean flow, and
check thermal shifts, signed reversal, work and atomic refusal/retry.

`fluid_transport::rheology::slit::LaminarCarreauSlit` resolves fully developed,
no-slip planar Couette--Poiseuille flow using the original thermal Carreau-Yasuda
source law. The circular-pipe and slit reductions share one immutable prepared
constitutive state and the original bounded stress inverse; neither duplicates
the material fit or reinterprets dynamic viscosity as a kinematic coefficient.
Geometry, uniform temperature and an explicit supported shear-rate ceiling are
fixed during a solve. The zero-shear core and both actual wall states must remain
inside the original temperature and reduced-rate support.

`LaminarCarreauSlit::with_gap` validates and rebinds a positive physical aperture
while retaining the original prepared temperature, law and source-rate ceiling;
it does not repeat thermal preparation. `zero_shear_viscosity_pa_s` exposes the
original prepared rest viscosity, and the solved profile's `stress_offset_pa`
exposes its midpoint stress without subtracting large wall stresses. The
`fs_tribo::resistive_film::carreau::CarreauFilm` consumer uses these source-owned
quantities for coupled pressure, squeeze and sliding Reynolds graphs while the
slit itself remains a single-passage constitutive owner.

Momentum gives `tau(y)=c+delta_p*(h/2-y)/L`. A bounded, bracketed offset solve
enforces the upper-minus-lower wall speed. Pressure-driven and sliding shear
therefore belong to one velocity profile. The resulting flux cannot be obtained
by adding separately solved nonlinear Couette and Poiseuille responses or by
substituting one apparent viscosity into a Newtonian conductance. Positive
quadrature resolves the inverse stress law and splits at the original
transition and zero stress. Near pure pressure flow, an interval identity avoids
subtracting nearly equal opposite-shear primitives.

The response retains pressure flow separately from mean-wall entrainment, both
work-conjugate resisting wall tractions, pressure work, requested wall work and
actual-profile dissipation. Heat is the nonnegative sum
`delta_p*Q_pressure + width*length*c*delta_U_accepted`. The accepted relative
wall-speed residual and its signed work defect are explicit; neither is hidden
as extra heat. Common translation transfers equal opposite pressure/wall work
and adds volume transport without changing viscous heat. Profile sampling
returns actual velocity, signed velocity gradient and shear stress, source
viscosity and local `tau*gamma` heating at the same frozen state.

Pressure mobility, relative wall resistance and the cross derivative use the
same accepted inverse stress tangent. Pressure mobility is a positive weighted
central variance, accumulated without subtracting nearly equal second moments;
the reciprocal wall/pressure cross response has the corresponding opposite
sign. Positive nested Simpson weights reuse previous samples and require two
successive refinement agreements. Exact constant-law and pure-Couette limits
avoid numerical quadrature. Evaluation uses bounded scalar scratch and allocates
no per-query storage. Checked solve and sampling entry points poll between
capped constitutive inverses and before publication; all refusals leave the
immutable owner available for an identical retry.

`LaminarCarreauSlit::profile_at_flow_and_walls` and its checked counterpart
resolve the pressure required by a prescribed signed net volume flow and both
wall velocities. Blocked flow can require adverse pressure and retain positive
wall-powered heat. The inverse uses the same profile's analytic pressure
conductance and a bounded physical stress bracket; the combined sliding and
pressure wall states must remain inside the original source envelope. A
constant-viscosity source uses the exact Couette--Poiseuille inverse. Pure wall
entrainment needs exactly zero pressure. Flow residuals and their relative
tolerance are evaluated after subtracting mean-wall entrainment, preserving
small pressure-flow corrections under common translation. The accepted profile,
requested flow, residual and trial count remain explicit; requested flow never
replaces accepted flow in work or heat. At most one pure-Couette evaluation and
the caller's capped number of nonzero-pressure trials occur, using scalar
scratch and the existing bounded local checks. Analytical blockage, independent
nonlinear pressure/heat, combined source ceilings, translation and late
cancellation regressions cover this inverse.

This is a source-owned planar passage, not a Reynolds graph, free-surface,
cavitation, entrance-flow or thermal-time solver. Width must be sufficient for
side walls to be neglected. The caller owns thin-gap/laminar validity, pressure
limits, source provenance and heat deposition; the rate ceiling alone certifies
none of them. Observed quadrature changes are numerical estimates, not rigorous
error bounds or calibration uncertainty. Focused regressions recover Newtonian
profiles and all work/tangent limits, compare combined nonlinear flow and heat
against an independent closed-form stress inverse and physical-gap quadrature,
differentiate pressure/sliding response, retain both thermal shift conventions,
check translation/reversal and source limits, and exercise bounded-work and
late-cancellation refusal with exact retry. All original pipe regressions also
remain applicable to the shared constitutive preparation.

Existing LBM callbacks consume the bound state through their native nonlinear
collision and viscous-heating path; the broader legacy lattice Carreau domain,
units and prepared relaxation-time policy remain intact. Constitutive regressions
cover independent formula/finite-difference tangents, limits, both thermal
conventions, support and extreme products. `material_thermal_rheology`
independently prescribes shear stress, checks actual relaxation/rate and
`mu*g^2` heating, source replacement invariance, temperature/shear changes to
momentum/thermal populations, common `a=2` SI/lattice parity, and transactional
calibration refusal.

## Concentric annular Carreau transport

`fluid_transport::rheology::annular::LaminarCarreauAnnulus` resolves one steady,
fully developed axial profile between two concentric no-slip cylinders. Signed
pressure drop and both axial wall velocities share the original thermal
Carreau-Yasuda source and constitutive inverse. Cylindrical momentum fixes
`tau(r)=C/r-delta_p*r/(2L)`; a bounded monotone stress solve enforces the relative
wall speed. Both actual wall shear states must remain inside the original
temperature/rate support. The inner wall is not replaced by an axis or a planar
slit, and independently solved nonlinear pressure/sliding fields are not added.

Logarithmic radius resolves the inner-cylinder shear. Gap and normalized annular
area retain narrow-gap geometry without subtracting squared radii. Positive
quadrature reuses the original scalar and central-variance integrators, splitting
at stress zeros, extrema and original source transitions. Constant-viscosity
limits use the analytic logarithmic Couette and annular Poiseuille response,
including a cancellation-free small-gap mobility. All per-query scratch is
bounded and allocation-free; checked profile and fractional-gap sampling poll
between capped constitutive inverses and before publication.

The accepted profile retains signed wall tractions and actual cylindrical
forces, pressure work, requested wall work, positive integrated `tau*du/dr`
heat, wall-speed residual and work defect. Mean-wall entrainment is separate
from relative flow; curvature-driven Couette transport belongs to relative
flow even at zero pressure. Common translation changes pressure/wall work by
equal opposite amounts without changing heat or the reported work defect.
Pressure mobility is a positive weighted central variance of normalized
annular area; the same tangent supplies sliding resistance and reciprocal
pressure/wall cross response. Numerical refinement changes remain estimates.

Geometry and uniform temperature are frozen. The caller owns laminar/entrance
admission, pressure limits and heat deposition; there is no rotation,
eccentricity, cavitation, slip, material transport or transient claim. Tests
recover analytical Newtonian profiles/forces/heat, compare nonlinear cylindrical
momentum against an independent exact stress inverse, differentiate both ports,
recover the planar narrow-gap limit, preserve thermal shifts and translation,
and exercise curved-wall source refusal and cancellation with exact retry.

