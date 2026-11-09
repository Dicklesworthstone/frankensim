# fs-em contract

## Original SCHOTT temperature response in Maxwell transport

`ThermoOpticOpticalResponse::try_lossless` composes the shared material-owned
TIE-19 law with an original transparent spectral glass reference. Native
interfaces, coherent films, angular spectra and bulk/interface frequency
dispersion use the same temperature-corrected absolute index. The material
owner differentiates the reference-index prefactor and thermal pole analytically
through ordinary wavelength curvature; no heated Sellmeier refit or frequency
differencing occurs. Fixed-size wavelength queries allocate no storage.

Actual temperature and SI permeability belong to a separate complete snapshot
on the original immutable material card. Reference-state glass coefficients
and catalog reference/coordinate atmospheres stay fixed. Rebinding always
starts from that reference. Unsupported coordinates, a wrong reference
temperature, disjoint bands and unsupported derivative neighbourhoods refuse.
Admitted phase-only table knots remain usable.

`try_lossless` continues to refuse nonzero original absorption.
`try_with_reference_absorption` explicitly holds the original intensity
absorption spectrum fixed while correcting real phase index. The complete
imaginary index and its analytic wavelength derivatives reach complex
permittivity, interfaces, film interiors/heating and angular spectra. This
retains zero attenuation dispersion for constant alpha and the original loss
dispersion for sloped spectra. The selected absorption policy participates in
source continuity and survives temperature rebinding. It is an authored
temperature-independent-loss approximation, not a newly supplied thermal
absorption law. No expansion, stress or causal transient reservoir is inferred.

Three material curvature controls and seven native Maxwell regression groups
cover independent TIE-19 derivatives, source/domain ownership, changed Fresnel
and film response, interface frequency derivatives, fixed reference air,
complex loss jets, local/interior heat and multifrequency source continuity.
The actual shared material's 16 tests passed in an isolated Rust compilation;
the native EM integration tests remain unrun.

## Original temperature-resolved Corning 7980 refraction

With `plane-wave`, `SilicaOpticalResponse::try_lossless` and
`DielectricResponse::SilicaOptical` carry the original `fs-material`
`Corning7980Silica` law into native interfaces, coherent films, angular spectra
and bulk propagation delays. The NASA CHARMS Table-3 fit retains its absolute
vacuum index and 30--300 K, 400--2600 nm domain. Preparation freezes its six
temperature polynomials once; wavelength queries allocate nothing. The same
complete scalar snapshot supplies the unique absolute temperature and explicit
SI permeability. Unowned coordinates, unsupported states and out-of-band
frequencies refuse. Phase index forms Maxwell permittivity; original group
index and analytic Sellmeier curvature separately drive group delay and GVD.

The caller explicitly selects the named silica source and zero absorption.
This is a homogeneous isotropic lossless approximation, without inferred
extinction, scattering, pressure/stress response, thermal expansion or a
causal transient energy reservoir. Native multifrequency source continuity
retains the exact original snapshot; finite-window reservoir requests refuse.
Six focused G0/G1/G3 controls in `plane_wave/optical/silica/tests.rs` cover
independent high-precision phase/group/thermal/GVD values, original quantity
and domain admission, oblique Snell/Fresnel/Brewster/TIR behavior, coherent
thermal film phase and interior fields, and multifrequency source identity.

## Fundamental Lorentz forcing in an original static magnetic bias

`SpatialHarmonicInductionResponse::static_bias_lorentz_loads` transfers the
peak carrier-frequency term `J_ac cross B_dc` from an accepted original
magnetostatic field and the retained AC current. Its cosine and sine loads
have no half-factor. Add them to the existing AC mean and twice-frequency
loads; both coefficients retain consistent P1 nodal force, native wrench
and virtual work. The original static gauge and trace values may differ,
while geometry, material identity and prescribed-trace mask must agree.
No electromagnetic solve, temporal sampling or state advancement occurs.

This path admits stationary conductors and selected cells without impressed
DC current. Motion or DC-driven selected cells require the additional
DC-current/AC-field term and are refused. Linear original magnetic laws,
source-bound fields, bounded selection/storage, canonical ordering and
transactional cancellation remain explicit. Nonlinear saturation, recoil
operating envelopes, DC terminals, impressed-current force and mechanical
feedback are outside this owner. Tests reconstruct full real J and static+AC
B before independent tetrahedral quadrature, exercise gauge/phase/scaling,
and reject incompatible source, trace, DC current, motion and work budgets.


## Coherent finite-window heat and loads from native 3D harmonic induction

With `harmonic-exposure`, `SpatialInduction::prepare_harmonic_exposure` retains
the already solved original A–phi fields at their actual frequencies. Every
response must share the exact source owner and prescribed steady velocity.
Equal frequencies combine in double-double arithmetic before all products;
distinct frequencies retain both their beat and carrier terms through the
existing analytic `fs-fft` temporal Gram kernel. Scalar and oriented tensor
conductor sources remain unchanged, including their resolved skin/proximity
current distributions. No additional field solve or temporal sampling occurs.

Uniform physical intervals return source-bound per-cell Joule heat and separate
signed Lorentz work. Original tetrahedral P1 mass weights integrate `J·E` and
`J×B` without vertex lumping. The same existing consistent load integrator
exports selected cells' mean nodal forces, native wrench and P1 virtual power;
duration times those means gives the corresponding impulse. Pair/workspace
budgets precede allocation, cancellation publishes no partial result, and
observations never advance accepted histories. Tests in
`tests/spatial_induction/harmonic/coherent_exposure.rs` compare real space/time
quadrature, moving and rotated anisotropic conductors, full-cycle recovery,
coherent cancellation, interval additivity and source/motion refusals.

`observe_with_static_bias` also retains the finite-window carrier impulse
`mean(J_ac) cross B_dc` from an admitted original static magnetic field with
the same problem and prescribed-edge support. It reuses the native static-bias
admission: conductors must be stationary, and all original conducting cells
must have zero impressed DC current. Each distinct carrier uses one shared
analytic zero-carrier moment;
DD currents sum before multiplication by the static field. The original limits
must admit an additional `12*C*F` temporal contractions and `F` cached kernels,
using the authored line count even when coherent lines coalesce. The native
bias field is borrowed without a field-sized copy. Tests in
`tests/spatial_induction/harmonic/bias.rs` independently reconstruct multitone
force, torque, virtual work and impulse; they check interval additivity,
full-period and short-window limits, linear bias-force scaling, exact added
work/storage bounds and transactional refusals.

These remain settled linear observations, excluding DC current forces/heat,
switch-on transients, nonlinear/hysteretic frequency mixing and feedback.
The bias extension transfers force on the resolved AC conduction current.
Static sources in insulating cells may provide the bias, but their own
forces/heat and AC interactions are outside the observation. DC-driven
conducting cells refuse, preserving the existing native force admission.
AC Joule heat is unchanged by the static magnetic bias in this admitted model.
Only the reported physical Joule heat may be deposited; adding the original
cycle-mean losses again would double count it. No complete finite-window
source/storage balance, continuum error bound or empirical validation is
implied by this observation path.

## Prepared source-temperature trials for coupled induction motion

`SpatialInduction::prepare_materials_trial` retains one source-pinned endpoint
operator set over the unchanged accepted A–phi history. Its reusable view
serves multiple prescribed-velocity trials without repeating material-weighted
assembly or replacing the accepted magnetic source, impressed-current history,
gauge, trace or clock. Every native solve still reports original magnetic
parameter work, conductor-frame Joule heat, Lorentz work and source supply.
The ordinary material-refresh entry point delegates through this same view.
This supplies operator reuse for finite-inertia thermal feedback; it does not
advance mechanics, choose temperature or deposit heat on its own.

## Native circuits on nonlinear and hysteretic solid electrodes

The `nonlinear::terminals::circuit` and feature-gated
`hysteretic::terminals::circuit` modules close the original 3D solid-conductor
electrodes against native RLC, mutual-inductance and ideal-transformer graphs.
The accepted owner retains original graph fluxes/charges, field cochains,
physical time and, for hysteresis, every vector polarization. Actual electrode
currents come from the same spatial charge reactions used by voltage control.

The nonlinear interface uses the existing bounded Newton–Krylov solver. Every
residual, Jacobian and final replay restarts the identical accepted circuit and
field histories. The constitutive field is evaluated at each proposed drive;
no affine field admittance is reused across saturation or magnetic memory.
Endpoint material/conductor operators are prepared once and shared across the
terminal trials. The original linear circuit retains its affine response cache.

All field variants share original graph preparation and physical accounting.
The numerical reference resistor and compensating current source cancel from
the physical terminal law and heat. Actual graph resistor losses, graph storage,
transformer work, electrode current balance and external-source work are
reconstructed before publication. Original source-bound graph resistance may be
refreshed without changing graph storage coordinates or field history.

Native circuit midpoint efforts/flows pair with the field's endpoint drive and
reaction through the same `dt * V * I`. The coupled nonlinear/hysteretic method
remains first order. Its stable balance includes graph and magnetic storage,
bulk Joule and pinning heat, authored graph heat, numerical field dissipation,
external and impressed supply, and the original parameter/mechanical terms.
Private interface work is internal. Trial exhaustion, cancellation, material
refusal or failed KCL/energy replay publishes neither history.

This remains fixed-geometry MQS with the underlying field's admitted material
model. Hysteresis requires the existing `magnetic-hysteresis` feature and fixed
fits; moving/rotating magnetic memory and temperature-dependent fit free energy
are not supplied by electrical graph coupling.

## Physical 3D electrodes with nonlinear and hysteretic induction

`spatial::vector_potential::induction::nonlinear::terminals` supplies voltage and
total-current electrodes on the original reversible A–phi field. With
`magnetic-hysteresis`, the parallel `hysteretic::terminals` owners retain the
original vector play memories as well. These are solid-conductor boundary
conditions: distributed current, skin/proximity response and local Ohmic loss
are determined by the same tetrahedral field and charge equations.

Voltage drive prescribes endpoint equipotentials and recovers incoming current
from charge reactions. Current drive prescribes signed incoming totals, groups
each electrode's charge unknowns and solves the equipotential voltage in the
same bounded Newton/FGMRES system. Totals must balance separately on every
connected conductor. Both forms reuse the original exterior-face, conductor,
disjoint-vertex and connected fixed-A-trace admission. Relative voltage
normalization makes the result invariant under common-mode shifts. Conversion
between controls preserves accepted edge potential, impressed load, physical
time and every hysteretic cell memory.

The linear, nonlinear and hysteretic owners share electrode layout, charge
regrouping, source RHS and transfer recovery. Source-bound conductor refresh
retains the electrode maps; reversible material refresh retains its original
parameter-work account. Each receipt exposes signed incoming currents, voltage
relative to its conductor reference and endpoint terminal work. Physical charge
checks include free nodes, imposed totals and whole-conductor balance. Field,
source, work and cancellation gates precede immutable candidate publication.
There is no isolated-field stepping accessor on an electrode owner.

The reversible work balance adds actual terminal work to impressed supply and
retains magnetic parameter work and mechanical export. The stationary hysteretic
balance is `delta W + Q_Joule + Q_pinning + D_numerical = W_impressed + W_terminal`.
Numerical endpoint dissipation remains separate from physical heat. Electrode
coupling preserves the underlying first-order MQS, fixed-geometry scope and does
not add moving hysteretic memory or temperature-dependent magnetic fits.

## Native three-dimensional hysteretic winding circuits

With the default-off `magnetic-hysteresis` feature,
`spatial::vector_potential::induction::hysteretic::winding` connects the original
tetrahedral vector-memory field to native Kirchhoff circuits. Its
`SpatialHystereticWindingCircuit` retains every accepted Whitney edge, vector
branch polarization, graph storage coordinate and winding current. All nested
terminal trials restart those same accepted histories; the native A–phi solve
continues to enforce conductor charge balance at each endpoint.

Reversible and hysteretic fields share one private spatial terminal kernel.
Original closed, signed winding profiles provide `J = F I` and reciprocal
`V = F^T delta(A)/dt + L_end delta(I)/dt + R I`. Current-input and voltage-input
ports retain their native graph meaning. Explicit passive end modes preserve
self/mutual storage and numerical endpoint loss; they do not replace the
distributed core. Source-bound strand resistance and ordinary graph-resistor
refresh preserve their original material choices, geometry and accepted states.
Resolved core conductors and homogenized stranded cells remain disjoint.

The combined receipt separates core eddy heat, vector pinning heat, strand copper
heat and graph-resistor heat. Private port work cancels once against actual
field/winding work. The stable total storage increment and independent numerical
loss close against external graph supply. Retained field source preparation is
shared across terminal trials, while each constitutive evaluation freezes the
accepted cell memory. Gauge, source closure, terminal, field and energy admission
all precede candidate publication; failures and cancellation advance no state.

The admitted scope remains stationary, fixed-geometry, first-order MQS with
source-selected isotropic vector hysteresis and constant fitted coefficients.
Moving magnetic memory, temperature-dependent fit thermodynamics and empirical
device qualification are not supplied by terminal coupling.

## Retained three-dimensional vector hysteretic induction

The default-off `magnetic-hysteresis` feature exposes
`spatial::vector_potential::induction::hysteretic`. Its
`SpatialHystereticProblem` retains the original tetrahedral geometry and one
source-selected isotropic vector play law per cell.
`SpatialHystereticInduction` retains the Whitney edge potential, impressed
current, physical time and every branch polarization. It starts demagnetized;
actual loading steps establish remanence and minor-loop history. A branchless
source supplies the ordinary linear background without a second material law.

The original A–phi coupling owns relative tree gauge, floating charge
components, topological admission, source compatibility, scalar or oriented
tensor conductivity and exact local P1 current/Ohmic recovery. Both reversible
saturation and hysteresis use the same bounded Newton/FGMRES kernel. Every
constitutive trial starts from identical accepted cell memories; the original
three-component return map and its differential reluctivity supply H and the
curl-curl tangent. Failed trials cannot advance magnetic memory.

The field receipt separates actual conductor Joule heat, irreversible pinning
heat on every magnetic tetrahedron, stable magnetic storage change and numerical
endpoint dissipation. The accepted work identity is
`delta W + Q_Joule + Q_pinning + D_numerical = W_impressed`.
Only the two physical loss channels are available for thermal deposition;
there is no added empirical core-loss map. Original Ampere, charge, flux and
energy gates all precede immutable candidate publication. Source-bound
conductivity refresh retains material occupancy and accepted vector history.

This is fixed-geometry, stationary-material, first-order MQS induction with
constant fitted magnetic coefficients and explicitly selected isotropy.
Prescribed traces remain fixed. Moving hysteretic material requires transport
and rotation of its stored polarization; a velocity argument cannot supply
that physics and is intentionally absent. The fit's original source evidence
is retained without claiming measured rotational validity, anisotropic
hysteresis, temperature-dependent magnetic entropy or continuum error bounds.

## Intrinsic switching kinetics on native three-dimensional induction

An explicitly sourced intrinsic viscosity profile now uses the accepted
`SpatialHystereticInduction` duration in every native B-driven constitutive
trial. The same accepted vector memory and physical `dt` remain frozen across
Newton, backtracking, physical electrode and graph-terminal trials. The native
branch potential contains `eta*|delta polarization|^2/(2dt)`; its consistent
effort and inverse permeability enter the original A–phi Newton/FGMRES solve.
That shared solver globalizes on the complete residual norm. Stored magnetic
energy remains the original free energy and includes no dissipative term.

Local and complete-field receipts distinguish dry pinning, intrinsic viscous
switching and resolved Ohmic heat. The complete physical identity is
`delta W + Q_pinning + Q_viscous + Q_Joule + D_numerical = W_impressed + W_terminal`.
Existing solid electrodes, their native circuits, stranded winding circuits
and finite thermal delegates inherit this behavior through the same field
owner. Their existing physical totals deposit the full constitutive heat once;
numerical endpoint dissipation is never thermalized. Original kinetic source
validity and temperature support participate in existing material admission.

Zero viscosity preserves the original rate-independent kernel. Positive
viscosity is an independently fitted intrinsic loss, not an estimate from
conductivity or total device core loss; resolved eddy loss must not be included
again in its fit. Fixed geometry, fixed fitted coefficients and the existing
default-off feature remain the scope. The independent six-tetrahedron tests in
`tests/spatial_induction/hysteretic/rate.rs` cover timed holds and relaxation,
separate local/global work accounts, late-cancellation replay, exact zero-viscosity
recovery and equivalence of physical voltage and total-current electrodes.

## Axisymmetric magnetic loads on native mechanical bodies

The original linear and nonlinear axisymmetric fields expose
`maxwell_wrench` and `impressed_current_wrench`, returning native
`fs_ga::Wrench` values in SI units. Their actual accepted geometry, material
snapshots, impressed currents and affine per-cell B remain authoritative.
A complete `2*pi` rotational sweep is required; a sector is not multiplied
into a full-body result and its missing azimuthal cut faces are not inferred.

The original `TriComplex2::trace_for_faces` owns the closed Maxwell surface.
Every edge of positive swept area must lie inside the mesh, with zero impressed
current and the same source-resolved induced-linear material on both sides.
All contour edges share that one exterior snapshot, rather than merely matching
numeric permeability. Nonlinear cores and bounded recoil magnets may lie inside
this collar. The exterior is assumed mechanically inert, without striction or
a field-dependent pressure law. Edges at r=0 have zero swept area and need no
fictitious negative-radius neighbour.

The axial traction is
`Bz*Hr*n_r + (Bz*Hz-Br*Hr)*n_z/2`, integrated under `dS=2*pi*r*ds`.
Two-point Gaussian edge integration resolves its cubic radial/stress product
exactly for the original affine discrete field, including sloping edges.
The separate impressed-current method integrates only
`Fz=-integral J_phi*Br dV`, using exact radially weighted triangle moments.
It is the free-current Lorentz contribution, not a substitute for the load on
a magnetic core or permanent magnet. It must not be added to an overlapping
Maxwell load.

The supplied `PlanarFrame` places the meridional [r,z] plane: the first
direction is radial at phi=0, the second direction is the magnetic axis, and
the frame origin is an axis point at z=0. The moment origin is expressed in
that same receiving frame. Azimuthal integration cancels radial net force and
intrinsic torque; the returned torque is the moment of the axial force about
the requested origin. Radial tractions can still produce hoop stress, which
this net wrench does not distribute to structural nodes. Native mechanical
dynamics can consume the result when body frame and moment origin agree.

The existing `MaxwellForceLimits` bounds original cells and edges, with
bounded polling, canonical mesh-order compensated reduction and no partial
publication on cancellation or refusal. This is an instantaneous static load;
moving geometry requires another field solve and a coupled mechanical/source
work account. It provides neither a force enclosure nor a qualification claim.

`tests/axisymmetric_force.rs` specifies source-free quadrupoles on axis,
sloping and toroidal surfaces; coil-force refinement against independent
Lorentz and fixed-current interaction-coenergy references; permanent-magnet
force refinement without an impressed-current surrogate; saturated-core loads;
polarity, frame/origin and native body acceleration; and source, contour,
sector, budget and cancellation refusals. Method agreement is numerical
consistency evidence; mesh and contour refinement remain necessary.

## Axisymmetric reversible saturation and permanent-magnet recoil

`planar::axisymmetric::nonlinear::AxisymmetricNonlinearMagnetostatics`
reuses the original regular `u=A_phi/r` geometry, exact azimuthal-current
moments and positive-swept-area exterior trace admission. Each triangle retains
an original `IsotropicMagneticMaterial`: source-bound linear permeability,
reversible H(B), or bounded oriented recoil. Recoil directions are meridional
[r,z]; radial remanence on an axis-touching cell refuses because it cannot have
a regular azimuthally invariant direction there. No material fit or source
amplitude continuation is synthesized.

Positive native Gauss quadrature through a Duffy triangle map defines one
consistent discrete material potential, residual and differential tangent.
The caller selects order 3 through 16, with an explicit total sample budget;
order three integrates the original linear and fixed-recoil pairings exactly.
Reversible H(B) supplies the original differential and transverse secant
reluctivities, assembled as a positive directional Gram form. The native
Newton--Krylov owner supplies bounded line search and FGMRES; sparse tangents
and positive Jacobi diagonals are cached per Newton point. Explicit cell,
coordinate, triplet, quadrature, Krylov-basis and iteration limits bound work.
Cancellation or failed final physical residual admission publishes no field.

Every trial checks all three vertices of the affine B field against the
original reversible source range. Recoil checks the affine H field's parallel
interval and transverse-norm bounds at those same vertices. Convexity then
covers the entire triangle, including extrema missed by interior quadrature.
Inadmissible source trials contract through native line search, without
clamping, extrapolating or advancing material state. Initial fields must be
admissible and must match the original prescribed traces.

The accepted field exposes actual A_phi/B/H samples, per-cell original recoil
operating margins, and quadrature-integrated storage W and coenergy W*.
Storage is zero at B=0 for induced cells and at B=Br for fixed recoil cells;
recoil coenergy includes Br.H and can be negative. The static diagnostic is
`W + W* = integral J_phi*A_phi dV + prescribed-boundary reaction pairing`.
It does not substitute B.H/2 for nonlinear storage or claim a transient source
work history. H is evaluated from the original law at the sample position,
without interpolation across material interfaces.

`tests/axisymmetric_nonlinear.rs` specifies exact saturated-current fields
and source energy for both current signs and sweep angles, agreement with the
original linear owner in both B components, shifted recoil storage and release,
whole-cell source admission, permeable saturated-solenoid refinement against an
independent Ampere integral, and bounded-work/cancellation refusals.
Native compilation remains unverified: the available environment killed the
existing asupersync dependency before reaching project code. Independent
quadrature, energy-gradient and tangent checks support the formulation without
claiming execution of those Rust tests.

This is a fixed-geometry, fixed-state isotropic magnetostatic estimate.
Nonlinear quadrature and the physical mesh require refinement near constitutive
segment changes; original source margins are not continuum error bounds.
Anisotropy, hysteresis, irreversible demagnetization, temperature evolution,
eddy currents, force qualification and automatic infinite-domain treatment
remain outside this owner.

## Axisymmetric magnetic fields from azimuthal currents

`planar::axisymmetric::AxisymmetricMagnetostatics` retains the original
meridional `TriComplex2`, complete per-cell material state points and admitted
linear isotropic permeability. Its continuous P1 unknown is
`u = A_phi/r` [T]. The physical cylindrical curl is
`B_r = -r u_z`, `B_z = 2u + r u_r`, and `A_phi = r u`.
Finite axis values therefore allow nonzero axial B while A_phi and radial B
vanish at r=0. No inverse radius or scalar-voltage gauge shift is used.

The weak equation pairs `curl(r u e_phi)/mu` with
`curl(r v e_phi)` under `dV = angular_span*r*dr*dz`; its source is
`integral J_phi*r*v dV`. Exact cubic triangle moments integrate this
operator, azimuthal-current load and recovered energy. Immutable preparation
reuses the original cancellable deterministic sparse assembly, mesh trace
admission and checked CG with a positive Jacobi diagonal. Explicit cell,
coordinate and staged-entry limits bound preparation. Cancellation or failed
material, trace or solver admission returns no field.

Every connected component requires a complete prescribed exterior edge of
positive swept area, with prescribed values of u. Axis-only and isolated-point
traces refuse. Unprescribed physical edges carry homogeneous tangential H;
the axis carries regularity. The angular span represents an azimuthally
invariant sweep and scales extensive quantities without changing local fields.
The caller supplies the finite outer boundary and its physical approximation.

Returned fields borrow the original source owner and retain actual impressed
currents. A_phi is quadratic, while B and H are affine within each triangle;
one-sided material-interface values remain separate. Per-cell and total
`W = integral B.H/2 dV`, the source pairing `integral J_phi*A_phi dV`,
prescribed-boundary reaction work and the actual defect in
`2W = source_pairing + boundary_work` are exposed with the native algebraic
solver report. These static pairings do not claim a transient work history.

`tests/axisymmetric_magnetic.rs` specifies uniform axial and source-free
quadrupole fields, an exact azimuthal-current solution with energy/work and
angular scaling, thick-solenoid refinement against independent Ampere-law
fields with a permeable core, and trace/source/budget/cancellation refusals.
Native execution has not yet completed in the available environment.

This is a fixed-geometry, fixed-state, linear isotropic magnetostatic estimate.
It supplies no nonlinear, anisotropic, recoil or hysteretic law, eddy currents,
force qualification, automatic infinite-domain treatment, independent angular
face excitation, experimental material validation or continuum error enclosure.

## Native winding circuits with spatial vector hysteresis

The default-off `magnetic-hysteresis` feature exposes
`winding::hysteretic::HystereticWindingCircuit`. It connects the original
`PlanarHystereticEddy` history to native Kirchhoff circuits through the shared
nonlinear terminal solver. Every trial retains the SAME accepted potential and
all original vector branch polarizations. Remanence, saturation and minor loops
are solved from that history; no effective inductance or scalar memory replaces
the spatial constitutive state.

The original signed winding profiles impose `J = F I` and receive the reciprocal
voltage `V = F^T delta(A)/dt + L_end delta(I)/dt + R I`. Native current-input and
voltage-input private ports support voltage and ideal-current excitation without
adding physical impedance. Existing end-winding self/mutual modes retain their
separate storage and endpoint numerical loss. Original source-bound stranded
resistance and graph-resistor reassembly reuse their existing owners; reassembly
preserves native fluxes, charges, field history, current and accepted time.
Stranded copper and resolved core eddy cells must remain disjoint.

Private circuit-port work cancels actual field/winding work once. The complete
receipt includes core Joule heat, vector pinning heat on every magnetic cell,
stranded copper heat and graph heat, each exactly once. Physical heat excludes
backward-Euler numerical dissipation. The combined storage increment uses the
field's stable per-cell magnetic increment, native circuit increment and original
end-mode increment. Independent voltage/current, work, heat and storage gates
must all close before one immutable candidate can publish circuit and vector
history together. Failed trials, exhausted budgets and cancellation leave the
accepted state intact.

Geometry and fitted magnetic coefficients remain fixed. Temperature-dependent
conductivity may refresh only through the original material sources and the
magnetic fit's admitted constant-coefficient interval. Motion, solid terminal
current injection, changing magnetic fits, anisotropic hysteresis and empirical
qualification remain outside this estimated first-order planar MQS composition.
A finite thermal receiver must accept the corresponding complete heat candidate
atomically.

The focused cases in `tests/hysteretic_eddy/winding.rs` compare independent
scalar-play voltage-driven minor loops with both winding polarities, ideal-current
drive with core eddy diffusion and end storage, the branchless backward-Euler RL
limit, and source/reassembly, budget and cancellation rollback.

## Fixed-recoil midpoint and Eulerian moving-conductor induction

`PlanarEddy` accepts actual fields from the original `PlanarMagnetostatics`
owner with source-admitted recoil magnets. Its midpoint residual uses the full
weak remanence moment, `K A_mid - load(Br/mu)`, including interface and exterior
terms. The original projected conductor mass, current-return constraints and
checked sparse solve remain authoritative. Recovered endpoint and per-cell
storage use the same fixed-remanence energy `V*|B-Br|^2/(2*mu)`; stationary work
closes `delta W + Q_Joule = W_impressed` without inventing electrical work or
heat for the retained remanence.

The existing Eulerian moving-conductor path reuses that affine gradient and
the actual midpoint B for Ohm's law and Lorentz work. Its account remains
`delta W + Q_Joule + W_mechanical = W_impressed`. The magnet's remanence,
orientation and occupancy stay fixed in the original frame. Every recovered
field retains the original recoil operating-envelope checks. The immutable
magnetic owner prepares its optional nodal remanence load once and reuses it
in static and transient solves; induced-only problems allocate no such cache.

Initialization rejects incremental winding derivative fields. Conductivity
refresh on a recoil cell requires the original fixed-magnetic-source admission;
replacing magnetic materials through the induced-only parameter-work path
refuses. This supports fixed recoil bias, not changing remanence, irreversible
demagnetization, or motion of a finite magnet or material interface.

Manufactured-source regressions cover interface and natural-boundary remanence,
analytic midpoint decay and second-order refinement, per-cell shifted storage,
source replacement refusal, a stationary magnet with a separate moving
conductor, zero-velocity agreement and complete-candidate recoil-envelope
refusal. The moving reference independently checks both heat and Lorentz work.

## Shared end flux in nonlinear differential and harmonic winding response

The planar and spatial incremental bases accept the same original passive
`EndWindingInductance` through `new_with_end_modes`. Operating linkage and
current-rate voltage evaluate its signed factors; the local differential matrix
contains every end self/mutual term. The resolved nonlinear tangent, fixed bias,
field source/gauge checks and singular-coordinate admission remain authoritative.
These end parameters are fixed linear storage, not a new saturation or loss law.

Planar harmonic winding impedance includes the complete end matrix. Its retained
field/copper block and original end factors independently recover terminal voltage
and phase-dependent cycle-mean storage, including native circuit replay. Large
end terms cannot mask the separately checked active-field power pairing. Copper,
solid-conductor and floating-conductor heat retain their original spatial owners.

With `harmonic-exposure`, finite windows contract the original shared-mode
currents before forming temporal products. Beat and carrier cross terms enter
both end storage and its rate; modal payload and pair work count against the
existing explicit budgets. Regression cases in `winding_incremental`,
`spatial_incremental_winding`, `winding_harmonic` and `winding_harmonic_circuit`
compare nonlinear finite differences, signed rank/cancellation, analytical lossy
impedance, loaded circuits and independent real-time two-tone quadrature.

## Shared end flux in linear planar and three-dimensional winding bases

`PlanarWindingBasis::new_with_end_modes`, its fixed-remanence counterpart and
`SpatialWindingBasis::new_with_end_modes` retain the original
`EndWindingInductance`. The resolved field Gram matrix and the complete passive
end matrix supply the native circuit's self/mutual block. End modes represent
only flux outside the resolved geometry; they cannot change the mesh fields,
permanent-flux offset, source admission or recoil operating envelope.

Physical linkage evaluates the original signed end factors separately from the
active field matrix, retaining cancellation of shared flux. Three-dimensional
field reconstruction also evaluates external stored energy through those factors.
Legacy constructors retain independent diagonal ends; their diagonal accessors
exclude shared-mode self terms. The complete matrix and original mode model are
separately inspectable. Singular terminal blocks still obey the original circuit
admission; no numerical leakage or regularization is introduced.

Planar Ohmic material replacement requires the complete original end model to
match. Equal independent diagonals or equal matrix diagonals cannot authorize a
change in mutual-flux orientation. G1/G3/G4 regression cases in `winding_circuit`,
`spatial_winding`, `winding_permanent` and `winding_ohmic` cover independent
two-coil evolution and power balance, polarity reversal, field/terminal energy,
unchanged permanent bias, full-model replacement and cancelled preparation.

## Spatial hysteresis with resolved eddy-current diffusion

The default-off `magnetic-hysteresis` feature also exposes
`planar::eddy::hysteretic::PlanarHystereticEddy`. It retains the original
`PlanarHystereticState`, all vector branch polarizations, fixed boundary
traces, electrical regions and clock. The native isolated-conductor
`Conduction` projection, algebraic air/floating-region admission, exact P1
Joule integral and backward-Euler quadratic operator remain their original
owners. Impressed current needs an explicit balanced return on every mesh
component; disconnected or relabeled electrical regions refuse.

The shared variational field solve adds `C_sigma delta A/dt` to both residual
and tangent and half its increment-work pairing to the objective. Every Newton
trial still starts from the SAME accepted magnetic history. The complete account
is `delta W + Q_pinning + Q_viscous + Q_Joule + D_numerical = W_supplied`.
Per-triangle dry pinning, intrinsic viscous switching, and exact nodal Ohmic E/J
and heat remain distinct. The original material endpoint remainder is numerical
loss, never deposited as heat.
The assembled quadratic work must independently agree with the spatial Joule
integral within the explicit magnetic work budget.

Electrical and magnetic sources retain the same original material card and
complete coordinate schema, while selecting their own property requirements.
Conductivity may refresh with the owned temperature; the original constant
magnetic fit is rechecked across that explicitly source-supported interval.
Changed fitted coefficients, nonthermal coordinates or source identities refuse.
Initial algebraic constraints use the accepted constitutive H directly, without
equilibrating away remanence or introducing an initialization impulse.

An explicit `ResolvedVectorHysteresis::from_isotropic_source_with_viscosity`
fit selects one original nonnegative coefficient `eta` [m s/H] per branch.
The field passes its owned positive `dt` into the native B-driven return map:
each branch adds `eta*|delta polarization|^2/(2dt)` to the incremental objective,
while physical switching heat is the FULL `eta*|delta polarization|^2/dt`.
The original magnetic storage and consistent inverse permeability remain the
field's energy and tangent owners. An explicit all-zero profile follows the
existing rate-independent kernel exactly; clockless statics and its authored
path continuation refuse active kinetics instead of inventing a timescale.

The native winding circuit already delegates every frozen-history terminal
trial to this induction owner, so its accepted interval also controls intrinsic
switching without a second circuit or field solver. Existing finite thermal
consumers deposit `irreversible_heat_j = pinning_heat_j + viscous_heat_j` once,
alongside independently resolved Joule heat. Kinetic coefficients retain the
same original source card, query, receipt continuity and explicit temperature
support as the other magnetic coefficients. They must describe intrinsic
switching independently of the eddy-current dissipation already resolved here;
conductivity or a total measured core-loss curve does not supply that fit.

`tests/hysteretic_eddy/rate.rs` adds independent scalar-kinetic/P1 conducting
relaxation accounts, intrinsic switching without a resolved conductor, exact
zero-viscosity recovery, and clockless/failure/cancellation replay cases. These
manufactured fits verify the model's numerical contract, not experimental
material authority.

Four G1/G3/G4 regressions cover independent scalar-play/P1 relaxation with
remanence, the exact branchless backward-Euler limit, balanced driven nested
loops with rate-independent pinning heat and rate-dependent Joule heat, and
whole-candidate rollback. This remains an estimated, first-order, fixed-geometry
isotropic MQS model. Reversals and rotation must be resolved/refined in the
physical loading history; motion, terminal circuits, variable magnetic
coefficients and experimental material validation are separate extensions.
Pinned-nightly formatting and bounded source/API/energy review completed.
Native execution remains blocked by missing constellation sibling manifests
and shared storage pressure; these regression claims are not execution receipts.

## Spatial vector hysteresis on the planar field solver

The default-off `magnetic-hysteresis` feature exposes
`planar::nonlinear::hysteretic::PlanarHystereticMagnetostatics`. Every original
triangle retains a `ResolvedVectorHysteresis` source and independent accepted
polarization memory. Branchless source laws supply linear background cells;
the caller explicitly selects the isotropic extension of the fitted coefficients.
Material names, permeability alone and scalar coupon evidence do not establish
a vector hysteresis law.

The original planar P1/FEEC Newton solver is shared with reversible magnetics.
For every global trial it evaluates the native B-driven vector return map from
the SAME accepted history, uses its consistent inverse permeability rotated
into the gradient frame, and minimizes the incremental potential
`delta W + Q_pinning + Q_viscous/2`; the static rate-independent limit has no
viscous term. Minimizing stored energy alone would omit pinning, while using
full viscous heat would double its constitutive resistance. Shape gradients,
physical cell volumes and immutable source
preparation are retained outside the Newton loop.

Immutable candidates contain actual `curl(A)`, all native material states,
per-cell stored-energy changes, separate physical pinning/viscous heat, and a separately
identified nonnegative endpoint numerical remainder. Impressed-current work
and prescribed-boundary reaction work are both included. Acceptance independently
checks the full field residual, local constitutive flux residuals, the sum of
absolute cell energy defects, assembled work balance and cancellation. A failed
trial or discarded candidate cannot alter accepted magnetization.

Gauge lifting and work pairings use a separate prescribed reference for each
connected component. `candidate_from_initial` accepts a distinct numerical
initial field while retaining the original physical history, so a feasible
equilibrium can be reached even when a stale-interior boundary lift would leave
the fitted field domain. Previously magnetized specimens are prepared through
their actual load history; zero-flux initialization never guesses remanence.

Seven G1/G3/G4 regressions cover independent one-free-coordinate play cycles,
remanence and nested minor loops, exact saturated-loop physical heat with
decreasing numerical loss under history refinement, heterogeneous linear-limit
parity, rotating prescribed fields and boundary work, transactional refusals,
admissible initialization near source limits, and disconnected gauge changes.
The implementation is a fixed-geometry, fixed-condition, isotropic estimated
rung. It does not add eddy currents, mechanical motion, thermal feedback,
anisotropic fits, path-independent coenergy, or experimental material authority.
Continuous-path claims require reversal/rotation resolution and refinement.

`planar::nonlinear::hysteretic::path` resolves one explicitly authored affine
current-density/exterior-potential segment through the SAME native field and
return-map owners. Starting impressed current is supplied explicitly and checked
against the retained constitutive effort on free coordinates, without solving
away an inconsistent initial state. Exterior differences are interpolated using
the original deterministic reference for each connected component. Every native
solve installs the authored target component gauge; accepted potentials are
never shifted afterward independently of their retained curl and memory.

A full increment is compared with two sequential halves. Only the two-half
history is retained, without extrapolated or averaged polarization. Independent
indicators control branch polarization, L1 per-cell physical pinning-heat error,
and the maximum of supplied-work defect and stable native numerical remainder.
Absolute allowances are weighted by interval fraction; relative allowances use
local changes. Compensated accumulation retains accepted cell heat, storage,
numerical loss, current work and boundary work, with an independent final local
and assembled energy gate. Numerical loss remains excluded from physical heat.

The complete segment has explicit field-solve and Newton-update caps including
failed and discarded comparisons; native per-solve Krylov, backtracking and
material caps remain unchanged. A typed flux-range refusal or bounded nonlinear
failure subdivides the loading interval: an invalid coarse boundary lift need
not prevent reaching a physically admissible endpoint. Failed Newton work still
counts against the complete-path cap, and depth exhaustion preserves its native
cause. Invalid input, cancellation and other source/constitutive refusals remain
failures. The bounded stack holds interval descriptors, not a history tree, and
reuses a rejected first half only from its identical starting memory.
Any refusal or cancellation leaves the complete accepted field unchanged.
`tests/hysteretic_magnetic/path.rs` covers analytical linear work, exact saturated
loop heat and minor-loop memory, an independently uniform vector-history
reference, disconnected gauges, valid continuation near the source limit after
a coarse failure, strict budgets and late-cancellation replay.
Indicators are empirical, not certified global bounds or event detection.
Waveform turns and source events require separate authored segments. This static
continuation does not integrate time, eddy currents or thermal chronology.

Verification status for this addition: source/energy review and the available
format/diff checks completed; native regression execution remains pending.
The execution environment ran out of shared storage and then disconnected
before any Cargo test could run. The fallback Rust 1.91.1 compiler smoke test
is not a substitute for the repository's pinned nightly/DSR gates.


## Frequency-resolved collisional plasma half-spaces

`plane_wave::IonizedGasWaveResponse` binds the original prepared
`IonizedGasHarmonic` to explicit lossless background epsilon/mu at one positive
angular frequency. It reuses the original neutral and mutual Coulomb drag and
forms `epsilon_eff = epsilon_background - i*sigma(omega)/omega` for
`exp(+i omega t)`. The background excludes those free carriers; the caller
declares its applicability to the frozen mixture. Supplied T/T_h, T_e and typed
pressure coordinates must match the retained gas; unmapped condition axes
refuse. No vacuum coefficients or composition-to-permittivity fit are inserted.

`DielectricResponse::IonizedGas` uses the existing interface kernel, including
negative real epsilon and the passive outgoing branch. The frequency-specific
cache preserves allocation-free normalized scattering. A different frequency
refuses and requires explicit preparation from the original transport.
`InterfaceReflection::transmitted_conjugate_h` retains the independent native
Fresnel numerator, avoiding subtraction from a reflection rounded to one.

`PlanarDielectricInterface::prepare_ionized_gas_wave` applies an actual complex
peak incident tangential field and recovers complete TE/TM vectors, including
normal TM E. The original material validates surface carrier and pair-relative
peak drift and supplies physical current, friction partitions and reversible
kinetic storage. The outgoing homogeneous wave only decays, so this surface
admission bounds all deeper points. Normal flux must match the native
transmission power, and `-2 Im(k_z) S_z` must equal original friction heat.
There is no conversion of reactive power into heat.

The retained wave evaluates arbitrary finite nonnegative depths with no
allocation or repeated carrier solve. Currents and carrier velocities acquire
the same decaying complex phase; heat and stored energies acquire its squared
amplitude. Depth-integrated heat uses `-expm1(2 Im(k_z) z)` rather than subtracting
nearly equal fluxes. Background electric/magnetic energy and carrier kinetic
energy remain separate. Thermal consumers apply area/volume/time and deposit
each physical heat source once.

Four G1/G3/G0 regressions cover independent Drude dispersion and plasma cutoff,
native vector Maxwell/Poynting identities through depth, original coupled
Coulomb thermal partitions, condition/frequency mismatch, excessive drift,
cancellation and zero/opaque fields. This is the inherited frozen, local,
unmagnetized weak-drift model, without sheath/ionization evolution or a measured
broadband qualification. Finite films and their beam/exposure consumers
explicitly refuse this law until whole-field drift admission covers their
standing-wave enhancement.

## Magnetized plasma Faraday rotation and circular absorption

`plane_wave::plasma::magnetized::MagnetizedIonizedGasInterface` binds the original
prepared Hall/cyclotron transport, explicit bound-charge epsilon/mu and a
lossless incident source. Normal incidence must be parallel or antiparallel to
the prescribed static B; zero B admits any normal. The original gas/background
condition checks are shared with the unmagnetized owner. Circular labels refer
to the outgoing normal, so reversing B swaps the two signed responses.

Each circular eigenmode uses its original conductivity and
`epsilon_background - i*sigma/omega`, with the existing passive Fresnel kernel.
Arbitrary transverse complex incident E produces actual reflected fields and
two transmitted waves with distinct phase/attenuation: Faraday rotation and
circular dichroism. No gyrotropic tensor is cast into a scalar material law.

`prepare_wave` validates the original transport at the complete surface field.
For normalized circular amplitudes, a carrier's instantaneous peak speed is
`(|u_plus|+|u_minus|)/sqrt(2)`, independent of relative phase. Both magnitudes
decay, so this bounds every deeper carrier and Coulomb pair-relative drift.
Admitting two separate modes would not establish the combined-field bound.

Retained modal responses provide allocation-free depth fields, currents,
carrier velocities, source-owned signed thermal partitions and kinetic storage.
Their circular orthogonality removes cycle-mean quadratic cross terms. Native
Maxwell residuals, interface power and attenuation/friction balance are checked;
integrated physical heat uses expm1. Background field energy stays separate
from carrier kinetic energy. Signed gyroscopic reactive power is retained,
and static magnetic no-work is never deposited as heat.

G1/G3/G0 regressions cover independent signed Drude/Fresnel modes and resonance,
zero-B native parity, spatial Maxwell/Poynting closure, Coulomb heat partitions,
Faraday rotation/dichroism, B reversal, frame rotation, combined drift refusal,
source conditions and cancellation. These are analytical implementation checks,
not experimental qualification. The model remains frozen, local and classical,
with normal-incidence homogeneous half-spaces only: no oblique/longitudinal
modes, finite films, kinetic spatial response, pulses or nonlinear wave force.

### Finite-volume magnetized carrier absorption

`PreparedMagnetizedIonizedGasWave::slab_absorption` integrates a positive-area,
positive-depth control volume within the original homogeneous half-space.
The two outgoing circular modes retain their separate attenuation, signed
electron/heavy pool powers and original drift-dependent exchange correction.
Their analytic exponential integrals use expm1, with a logarithmic range
fallback for otherwise representable products; no depth quadrature, port-power
subtraction, transport reevaluation or new Maxwell solve occurs. Disjoint
slabs add to the original total absorption, including weak, deep and opaque
limits. Signed pool closure is checked against total friction.

The report borrows its exact wave/source and returns the actual geometric
volume, original friction shares and already-included signed exchange
correction. Coherent carrier kinetic and lossless background field energies
are separate storage, never extra heat. A finite-energy consumer must bind
this volume to its actual parcel and declare the held-field interval. This
control volume has no second reflecting boundary and does not add a finite
film, temperature evolution or a continuously refreshed wave solve. Focused
tests compare independent local-field quadrature, Poynting loss, subslab
additivity, signed Coulomb partitions and weak/opaque/zero-source limits.

## Conduction at independent common-voltage references

The planar static, material-response and prescribed-motion current paths retain
the existing scalar owner's relative potentials through E/J recovery, contact
jumps, boundary reactions and terminal work. Each electrically disconnected
component uses its own prescribed reference; finite conducting contacts join
those references through the original operator. A large common voltage cannot
erase a physically solved resistive or motional drop before Joule heat or
Lorentz work is evaluated.

`ConductionField`, `MovingConductionField` and `ConductionResponseReport`
retain both relative potentials and their component references. Absolute
`voltage_v` is a reporting view, reconstructed only after physical admission,
with authored voltage traces preserved exactly. It can round away differences
smaller than its own floating-point spacing. Two-terminal circuit recovery
scales the retained relative solution and reference separately, then rechecks
the original field/current/heat laws at the actual drive.

G1/G3 tests cover analytical series current and heat, independent disconnected
references, finite-contact jumps and losses, actual nonlinear material queries,
open motional screening, and the signed electrical/mechanical/Joule split under
common biases larger than the physical drops. This changes numerical gauge
handling, not the admitted fixed-geometry Ohmic/material/motion models.

## Stationary resistive interfaces beside moving conductors

Planar and spatial `solve_with_motion` now admit original finite-resistance
contacts when BOTH adjacent original cells have exactly zero prescribed velocity.
Moving bulk elsewhere retains `J = Sigma(E + u cross B)`; a moving or sliding
contact still refuses even when its magnetic field is zero. The admitted
zero-thickness stationary interface retains its original source, paired traces,
exact P1 current and `integral([V]^2/R'')` heat. Its mechanical work is zero.

Both moving-field owners expose `contact_response` and `contact_power_w`.
Cell Joule arrays contain only bulk heat; total `joule_power_w` includes bulk
and surface heat once. `electrical_power_w` remains the signed BULK integral
of laboratory `E.J`. Boundary and Hodge work compare that integral PLUS
surface work with the terminal supply, and the complete balance remains
`terminal = bulk heat + contact heat + mechanical output`.
The spatial owned snapshot retains these original surface results for explicit
thermal partition by the coupling owner.

Independent planar and tetrahedral series references exercise stationary parity,
shorted braking, loaded generation, motoring, exact contact jumps/heat and
refusal of motion on either contact side. A current-driven planar case also
crosses the contact with one component reference under a common voltage larger
than its physical drops. These are synthetic quasistatic model checks, with
no sliding-contact, frictional heating, moving occupancy or induction claim.

## Resolved electrostatic forces, motion and floating electrodes

The native planar and spatial electric fields now expose `virtual_work` for a
complete original nodal displacement. One shared allocation-free cell
contraction differentiates the actual P1 dielectric energy using original
fields, gradients and physical volumes. Local displacement differences make
common translation cancel exactly. With held electrode voltages the returned
mechanical work is the energy derivative; with held charges it is the negative
re-equilibrated energy derivative. This path admits fixed explicitly isotropic
coefficients, with zero free volume/sheet charge, impressed polarization and
imposed displacement. Tensor transport, electrostriction and unowned source
work are refused. Virtual mesh work alone is not a finite trajectory or an
automatic solid nodal-traction distribution.

`SpatialElectrostatics::capacitance_matrix` extracts the full tetrahedral
Maxwell matrix with one native solve per nonreference terminal. Electrical
components follow physical volume connectivity and shared equipotential
electrodes, with independent voltage references. The original stiffness Gram
form is checked against native charge reactions; scalar and tensor source laws
stay attached. Actual drives superpose retained relative potentials and recover
the complete original field, checking actual-amplitude charge and energy.
This is the homogeneous source-free response, not an implicit offset for fixed
polarization, sheets or bulk space charge.

`SpatialCapacitanceMatrix::prepare_floating` adds prescribed TOTAL floating
free charge through the existing equilibrated dense Cholesky owner. Each
electrical component must retain a voltage reference before factorization.
The Schur matrix, affine induced retained charge and trapped-charge energy all
remain available; trapped energy is already part of the native field store.
Small charge-induced fields use component-relative voltages with separately
retained reference values, so large independent common biases do not erase
them. Recovery needs no further FE solve and retains original field/material,
charge, energy and force-validity checks.

`fs_couple::electrostatic::motion` connects a two-electrode planar field to
a source-derived translating body. An authored affine mesh path retains vertex
lineage, topology, cells, dielectric laws and extrusion thickness; all cell
signed-area ratios must remain positive over its complete coordinate interval.
This does not establish global nonintersection or contact. New coordinate
trials rebuild the original Hodge and solve its field; accepted and reference
responses are reused. The same analytic virtual work supplies instantaneous
force. The existing native articulated impulse and bounded scalar feedback
owners determine a common mechanical step, with fixed isolated charge or an
ideal held voltage. Independent field-storage, signed electrical supply,
native kinetic and external mechanical work checks gate publication. No
discrepancy is assigned to heat or repaired by changing state.

That motion is one guided electrostatic-equilibrium coordinate with frozen
strain-independent scalar permittivity. The caller owns its physical mesh
motion and complete moving inertia. Its reduced axial support wrench is not a
recovered full electrode traction or torque. Finite source circuits,
three-dimensional moving interfaces, electrostriction, polarization evolution,
charge leakage/transport and electromagnetic propagation remain separate
extensions. Energy admission does not certify continuum force or temporal
error.

Four multi-electrode, four floating-charge, seven native virtual-work and five
guided-motion integration regressions are included with synthetic materials.
The two exact production contraction-kernel tests compiled and passed with the
pinned Rust toolchain. Full native field/coupling compilation and those twenty
integration tests remain unexecuted in this environment because the shared
host exhausted memory and writable disk before its required dependency build
completed. Formatting/static source checks are not a substitute for that run.

## Three-dimensional terminal power under prescribed conductor motion

`SpatialConduction::solve_with_motion` solves the original tetrahedral current
problem with `J = Sigma(-grad(V) + u cross B)`. The actual motional weak load
changes voltage and terminal current while reusing the original SPD Hodge,
source-bearing scalar/rotated-tensor conductivity, current traces and gauges.
An open component screens a conservative motional field without an invented
short-circuit current. The stationary public API and field semantics remain
unchanged; its recovery allocates no motion buffers.

The accepted result retains the original field owner and prescribed cellwise
velocity/B snapshots. It separates nonnegative conductor-frame Joule heat,
`J cross B` force, signed mechanical output and signed boundary electrical input.
Independent local, Hodge, terminal and current gates enforce
`P_terminal = Q_Joule_rate + P_mechanical`; generation reverses electrical power
without making heat negative. Only Joule power belongs in the thermal receiver.
Additional assembly/recovery work and retained storage are linear in cells and
vertices; the native sparse solve and its iteration limits remain authoritative.

This is a nonrelativistic fixed-domain DC snapshot with prescribed piecewise
constant B and velocity. It does not evolve induction, magnetic sources,
mechanics, finite moving interfaces or thermal advection. Finite-resistance
contacts require stationary adjacent cells as described above; motion at either
surface still needs an owned moving-surface work/heat law. G0/G1/G3/G4
fixtures in `tests/spatial_conduction/motion.rs` cover open-circuit screening,
shorted braking, loaded generation, motoring, layered current drive, anisotropic
scene rotation/volume scaling, stationary parity and transactional refusals.
These synthetic-law fixtures do not qualify an empirical machine or material.

`SpatialMovingConductionField::lorentz_loads` transfers its SAME retained
`J cross B` through the native P1 load kernel into original-vertex nodal forces
and canonical `fs_ga` wrenches. Constant cell density integrates as `V*f/4` per
vertex; selected-cell scattering preserves total force, torque and P1 virtual
work. `lorentz_loads_for_velocity` additionally admits the receiver's actual
P1 velocity only when each selected cell mean matches the prescribed motion.
This local check prevents cancelling velocity/work defects and does not assume
that arbitrary discontinuous cell motions admit a shared nodal trace. Original
cell/source ownership, limits, canonical selection and cancellation are retained.
These volume loads must not be added to overlapping induction or Maxwell loads.
Consuming `into_snapshot()` transfers accepted arrays without copying them to
an owned `SpatialMovingConductionSnapshot`; a coupled candidate retains the
exact refreshed conductor alongside this diagnostic record. The snapshot has
no topology-dependent force extraction or reattachment operation.

Solved-field tests in `tests/spatial_conduction/motion/mechanical.rs` cover
braking/generating/motoring force and torque, mesh refinement and shared nodes,
selected-body physical power, affine P1 and transformed rigid virtual work,
local velocity mismatch refusal, and native articulated-body linear/angular
impulses. The finite held-load impulse uses its actual average-velocity kinetic
work, including the quadratic impulse term; extraction does not re-solve the
DC field for the subsequent mechanical velocity or close dynamic feedback.

## Consistent three-dimensional Lorentz loads for mechanical momentum

`SpatialInductionCandidate::lorentz_loads` and its nonlinear counterpart
transfer a selected set of original conducting cells from the same accepted
induction update into `SpatialLorentzLoads`. For each tetrahedron the native
P1 load is `F_a = volume * (sum_b f_b + f_a) / 20`, where `f = J cross B`.
This preserves net force, torque and every P1 virtual-work pairing. Vertex
sampling with lumped volume weights or replacing the density with its mean
does not generally preserve these mechanical quantities.

Original vertex order is retained for nodal-force consumers. A canonical
`fs_ga::Wrench` transforms both force and torque to an explicitly supplied
receiving frame and moment origin, ready for native `fs-mbd` inputs. Selected
cell IDs are canonicalized; duplicate, nonconducting, over-budget and cancelled
selections refuse without changing state. Near-source moment reduction and
compensated accumulation limit avoidable translation and cancellation error.
Receiving-frame origin subtraction precedes addition of the rotated local
reference, preserving its lever arm under a large common frame translation.
The extraction is linear in mesh vertices plus selected element work, with
sorting and binary source lookup only for selected cells; wrench queries are O(1).

Linear midpoint loads and nonlinear endpoint loads retain their own update's
time stage. Pairing with the admitted prescribed velocity recovers the selected
Lorentz work. These are free-conduction-current volume forces, excluding
magnetization/interface forces; an overlapping enclosing Maxwell load must not
be added. This transfer neither moves geometry nor closes mechanical feedback.

Production-kernel regressions cover analytic nonuniform
force/torque, general P1 and rigid virtual work, arbitrary frames and origins,
refinement/shared-node assembly, translated geometry, transactional refusals,
and actual native articulated-body linear/angular impulse response. Additional
linear and saturated solved-field regressions compare the transferred loads
with their original induction work.

## Three-dimensional dielectric fields, circuit storage and force

`spatial::electrostatic::SpatialElectrostatics` solves the fixed-geometry,
linear dielectric problem on the original tetrahedral P1 complex. It uses the
shared spatial geometry, checked Jacobi-CG and weighted Whitney Hodge
`K = d0^T M_epsilon d0`. The complete material state at each cell supplies
absolute-SI scalar or explicitly oriented principal permittivities through
`SpatialDielectricLaw`; geometry, gradients and the immutable operator are
prepared once for repeated boundary and charge excitations.

Typed exterior faces carry P1 voltage or uniform outward displacement D.n.
Unlisted faces have D.n=0. `solve_with_sources` additionally accepts uniform
free-charge sheets on original internal faces and fixed mesh-frame impressed
polarization per cell. Constant volume charge and face loads use exact P1
source moments: rho*volume/4 and sigma*area/3 respectively. The signed Neumann
load is minus the outward displacement. Polarization contributes
`b_P = integral grad(N).P` to `K V = q_volume + q_sheet + b_P - integral N D.n`.
It must not duplicate induced polarization already contained in epsilon.
Compatibility includes volume, sheet and boundary load on each disconnected
component before admitting a numerical gauge. Physical prescribed potentials
and numerical gauges remain distinct, so a gauge cannot hide an incompatible
charge. Existing source-free spatial magnetic/current consumers keep the same
boundary wrapper.

Solved potentials remain relative to each component through E, D, electrode
reactions, source work and stored energy. Absolute voltage is reconstructed only
for reporting, retaining authored voltage traces exactly. Reactions use voltage
differences and compensated accumulation. Admission checks total and component
charge balance, free-node/gauge residuals, and dielectric energy against both
the weak-source pairing and Hodge quadratic. Scalar and tensor material jumps
are retained cellwise; raw interior normal-D jumps are diagnostics, not a
smoothed field, equilibrated flux or continuum error certificate. Corrected
interface residuals account for the prescribed sheet's normal-D jump.
Cell displacement is `D=epsilon E+P`; storage remains
`integral E.epsilon.E/2` at fixed impressed P, with the full signed
`integral E.P` reported separately. Reaction charges subtract all authored
loads, including b_P, so bound polarization is never relabeled free charge.

`SpatialCapacitance` admits two disjoint unions of original exterior faces,
solves their zero-volume-charge unit excitation, and checks electrode charge
against C=2U. Its capacitor branch belongs to the existing native circuit
owner. Actual terminal voltages or charge reconstruct the original spatial
field and recheck amplitude-scaled charge/energy without another Krylov solve.
Do not count circuit capacitor storage and the same field energy twice.
This owner does not infer floating/multiple electrodes or assign impressed
volume charge to a passive two-terminal capacitor.

`SpatialElectricField::maxwell_wrench` uses the existing closed-surface
Maxwell integral and SE(3) wrench convention with D and E. Every integration
face is internal, and both neighboring tetrahedra must retain one identical
complete isotropic source, zero prescribed volume charge and impressed
polarization, and no charged-sheet/electrode/imposed-displacement support.
The homogeneous exterior is assumed
nondeforming, without electrostriction or fluid-pressure forces. This is a
fixed-geometry instantaneous load; advancing motion requires explicit
electrical, stored and mechanical work accounting and field/surface refinement.

`fs_couple::electrostatic::spatial_coulomb_volume_load` separately transfers
the accepted field's free-volume-charge force density `rho_free E` into the
original solid P1 volume-load integrator on exactly matching positions and
ordered tetrahedra. It retains constrained-node forces for support reactions
and returns the native displacement/stress load case. Explicit selected
charges must transfer their momentum to the mechanical material. This term
does not supply electrode/sheet/polarization/interface forces, charge transport,
electrostriction or moving-field feedback; it must not overlap an enclosing
Maxwell load on the same body. The four direct integration regressions in
`fs-couple/tests/electrostatic_loads.rs` cover SI resultant/torque/virtual work,
an actual supported solid solve, selection and geometry/cancellation refusals.
They are included but remain unexecuted under the same build resource limit.

The independent `electrostatic::tetrahedral::TetDielectricOperator` remains
the explicit-tensor/polarization weak-form primitive for constitutive coupling.
Both paths use the native FEEC kernels. The prescribed polarization does not
infer measured material authority, intrinsic polarization energy, source work
for evolving P, hysteresis or electrostriction. Nonlinear/dispersive 3D
dielectrics, infinite exterior boundaries and moving-field work remain outside
this boundary-value owner. Source preparation reuses the original operator;
it adds no field solve or independent dielectric discretization.

`tests/spatial_electrostatic.rs` and its force submodule cover analytic layered
capacitance and D continuity, retained high-bias fields, a 3D Poisson refinement
sequence, pure-Neumann/disconnected compatibility, rotated anisotropy, native
RC charge/heat/storage, actual-voltage/charge recovery, Q E force refinement,
wrench frames and the original source/cancellation/work-budget refusals.
Its source submodule adds analytic sheet/polarization fields, shorted and
open-circuit limits, reversed-cell rotated tensors, mixed/pure-Neumann charge
compatibility, source ordering, gauge invariance and force-collar refusals.

## Gauge-stable dielectric charge and energy

Planar electrostatic elimination and physical recovery remove the independent
voltage reference of each connected component before multiplying field or
charge operators. Electrode-basis superposition uses each connected terminal
group's reference, including ideal electrodes that join disconnected mesh
pieces. Large unrelated voltage offsets cannot become fictitious charge or
inflate the Krylov load scale of another component. Charge and energy retain
their original absolute acceptance gates; no tolerance is relaxed to hide
gauge cancellation. Focused charged-domain and independent-capacitor tests
exercise these physical invariances.

The shared planar scalar solve uses an admitted Jacobi inverse in the existing
checked-CG recurrence. Diagonal setup shares Dirichlet elimination's row pass
and retains one scalar per free coordinate; matrix and RHS remain unscaled.
Acceptance still uses the independently recomputed residual of the physical
equation. Three native exact-source scalar tests pass for heterogeneous
capacitor response, strict iteration reduction against identity preconditioning,
unrepresentable diagonals, gauge invariance and cancellation/replay. This is a
fixture iteration-count result, not a wall-clock throughput claim or a full
crate/repository test pass. New planar/spatial source integration tests are
present; their full Cargo execution remains unverified in the implementation
environment because shared memory/disk pressure blocked dependency compilation.

## Source-bound load resistance on continuing three-dimensional windings

`SpatialNonlinearWindingCircuit::with_graph_resistors` attaches the existing
`GraphResistorMaterial` sources to every ordinary resistor on the native graph.
Its source-resolved geometry and electrical basis determine resistance; later
updates return a separate owner at the same physical time with exactly retained
A, winding current/linkage, native fluxes, capacitor charges and algebraic receipt.
Resistance changes add no magnetic-storage parameter work. The next original
native descriptor step computes the actual stage efforts without an artificial
initial equilibrium projection.

Reassembly retains original private terminal order, source indices and all native
storage coordinates. The unchanged Hamiltonian is checked explicitly. Winding
profiles and exact Whitney load caches are shared through `Rc`, so an ordinary
resistor update does not copy or rebuild those field payloads; identical source
bindings also reuse the original native assembly. This is an allocation reduction,
not a measured runtime claim. Every accepted state and step retains complete graph
material bindings beside its distinct graph-local heat records. A finite thermal
consumer separately owns original source pinning and temperature-refresh admission.

The focused case in `tests/spatial_induction_winding/resistors.rs` changes source
resistivity after nonzero-time evolution and independently checks backward-Euler
winding response, native trapezoidal RLC response, heat, work and storage. It checks
retained source metadata and typed cancellation/support refusals. Pinned Rustfmt
and source/accounting review passed; native execution and DSR remain unrun here.

## Native winding circuits on nonlinear three-dimensional induction

`spatial::vector_potential::induction::nonlinear::winding` connects the original
3D Whitney A–phi history directly to native Kirchhoff circuit terminals.
`SpatialNonlinearWindingCircuit` replaces selected original inductor placeholders
with private power ports, retaining unrelated inductors, capacitor charges,
transformers and external source indices. The same interface Newton owner admits
current-input and voltage-input causality; ideal-current-driven windings need no
invented series resistance or magnetic storage. Unsupported algebraic topologies
remain native refusals.

Each `SpatialNonlinearWindings` profile is the complete signed cellwise `J/I`
including turns and closed return paths. Internal face continuity and exterior
closure are checked even where tangential A is prescribed. Its original exact
Whitney load `F` is cached once. Every terminal trial starts from the same accepted
field and circuit histories and uses
`V = F^T delta(A)/dt + L_end delta(I)/dt + R I` at endpoint winding current.
Mutual response, reversible saturation, bulk eddy diffusion and admitted prescribed
bulk-conductor motion therefore enter through the actual field rather than a
fitted inductance. The pure insulating-core/air-core limit uses this same owner
with empty bulk conducting support and no electric-potential unknowns.

`SpatialWindingConductor` supplies the original scalar electrical source and
explicit conductor volume fraction. Resistance is the physical integral
`int rho |J/I|^2 / filling dV`; recovered intrinsic strand current is `I*(J/I)/filling`.
Distinct insulated bundles may share a tetrahedron only within its total volume
fraction. Every active profile cell needs its own source, and stranded cells are
disjoint from resolved bulk eddy conductivity. Spatial `E dot J` heat independently
checks terminal `R I^2`; it is not added again as bulk-conductor heat.

The combined account includes native graph storage, original magnetic free energy,
explicit omitted end self-inductance, external supply and signed mechanical work.
Private terminal work cancels internally. Magnetic parameter work is separate;
backward-Euler magnetic and end-inductor numerical dissipation never becomes
constitutive heat. Stable field-energy increments are retained instead of
subtracting aggregate magnetic stores. Independent native-work, field-pairing,
terminal-current/voltage, local heat and complete energy gates must all pass.

Original-source conductivity and linear or source-separable magnetic temperature
refresh preserve A, physical current history, geometry, ports and circuit storage.
The endpoint material and conductivity operators are prepared once per circuit
candidate and shared across its terminal Newton/Jacobian trials. Unchanged sources
reuse the prepared operators. Bounded profile storage, native graph dimensions,
field work, nested trial counts and cancellation remain explicit. This is an
allocating first-order coupling; no timing or error-certification claim is made.

This extension admits stationary stranded profiles, scalar stranded resistivity,
isotropic reversible/fixed-recoil magnetism and the existing scalar or tensor bulk
conductor owner. Prescribed Eulerian motion inherits its fixed-occupancy and material
transport restrictions; it does not move winding geometry or evolve mechanics.
Solid terminal current injection, moving interfaces, strand skin/proximity loss,
magnetic hysteresis, composite caloric homogenization and engineering validation
remain separate work. The standard source/flux-linkage convention is also described
in the [3D A–phi formulation reference](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v261/en/Subsystems/Maxwell/Content/Maxwell/APhiFormulationinMaxwell3DTransient.htm).

Four G1/G3/G4 tests in `tests/spatial_induction_winding.rs` independently integrate
Whitney linkage and source B–H energy, compare coupled two-loop backward Euler
with unrelated native RLC trapezoidal evolution, continue nonlinear eddy history,
and check actual trial exhaustion and late cancellation/retry. Rustfmt parsing and
source/accounting review were performed; these newly authored native tests and the
required `dsr quality --tool frankensim` gate have not run in this environment.

## Three-dimensional distributed windings and native mutual circuits

`spatial::vector_potential::winding::SpatialWindingBasis` connects complete
closed 3D unit-current profiles to the original source-bound Whitney field.
Profiles include physical turns, orientation and return paths in the original
mesh frame. Fixed induced-linear materials and homogeneous tangential magnetic
traces admit one native unit-current solve per winding; original local current
continuity, complete relative curl-kernel compatibility, gauge and source
checks remain authoritative. Current cannot disappear through a prescribed
boundary. Preparation explicitly bounds winding count and retained field data.

Every self/mutual entry compares `integral(B_i.H_j dV)` against the original
Whitney pairing `integral(A_j.J_unit_i dV)` in both directions. The complete
reciprocal matrix includes source geometry and all meshed end effects; only
explicitly declared omitted series self-inductance is added. There is no
empirical coupling coefficient or invented leakage for dependent modes.

`SpatialWindingCircuit` replaces selected original inductors through
`assemble_circuit_with_mutual_inductance`. It preserves caller winding order,
original branch/source indices, unrelated circuit stores and transformers.
Consistent initialization, actual coupled currents and audited stepping use
the native DAE. Full positive-definiteness admission remains with that owner;
dependent modes need physical leakage or an explicit ideal-transformer model.

Ordinary circuit steps reuse the small native system. Optional field
reconstruction performs no new 3D solve: original curl, stiffness, source and
constitutive operators recover potential, face flux, B/H, current and cell
energy at the actual currents. Current compatibility, Ampere residual,
divergence/conforming-flux and work allowances are rechecked at that amplitude.
Declared external inductive energy is separate from meshed field energy;
both already belong to the native winding storage and are never added twice.

Geometry, material sources, winding profiles and boundary traces stay fixed.
This linear reduction does not imply saturation, remanence, skin/proximity
diffusion, thermal evolution, moving winding transport or an infinite exterior.
G0/G1/G3 fixtures in `tests/spatial_winding.rs` independently normalize a compact
closed current loop, compare direct and reconstructed fields, check gauge,
polarity and material response, and solve a native mutual RL circuit against
independent coupled equations with unrelated LC storage. Return, source, budget,
dependent-mode and cancellation refusals preserve the native physical gates.
These authored fixtures do not constitute measured-device qualification.



## Driven solid-winding motional EMF and reciprocal Lorentz work

The nonlinear moving field and native winding circuit now admit driven solid
sections through their original current-constrained regions. With
`q = delta(A)/dt + v.grad(A)` and `u = n I/integral(sigma dA)`, the same
projected conductivity operator solves the magnetic residual and reconstructs
`E_conductor = u - P_sigma(q)`. Actual total `J = sigma E_conductor` supplies
both exact P1 Joule heat and `J cross B` work. The conductor-frame heat,
laboratory electric field and motional electric field remain distinct receipts.

The native terminal law adds the source-derived back-EMF
`V_motion = depth*n*mean_sigma(v.grad(A))` to direct
`F_new^T delta(A)/dt`, explicit end-inductance voltage and original DC loss.
This voltage is computed directly even at zero current. Separate DC and signed
motional drive work reconcile the complete field source work with the circuit;
neither creates an extra heat deposit. Negative mechanical work is supplied
by the imposed-motion owner, while physical Joule heat stays nonnegative.

`step_with_eulerian_velocity_and_sources_candidate` refreshes only original
conductor/resistor sources, retaining field/circuit history, magnetic laws,
end storage and topology. Solid profiles follow the refreshed conductivity;
an unchanged magnetic potential does not acquire artificial transformer EMF.
The sparse nonlinear operator, nonsymmetric FGMRES solve, true residual,
Peclet limits and original energy/cancellation gates remain authoritative.

This is prescribed Eulerian motion with invariant material occupancy and
coefficient fields along the motion. It does not translate finite interfaces,
move permanent magnetization, or account for thermal advection. Stranded cells
and fixed-recoil magnets remain stationary. The conductor-frame convention and
transport restrictions follow the original
[quasistatic Lorentz formulation](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.09.html)
and [Eulerian velocity scope](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.040.html).

Focused G0/G1/G3 fixtures in `tests/nonlinear_winding_solid/moving.rs` cover
zero-motion equivalence, open-current motional voltage, motor/generator work,
total spatial heat and force, original conductivity refresh, and unchanged
histories on source/budget/cancellation refusal. They are authored numerical
fixtures; native execution and empirical device qualification remain separate.



## Causal solid-winding diffusion in native nonlinear circuits

`winding::nonlinear::NonlinearWindingCircuit::new_with_solid_regions` and its
explicit-port variant drive whole connected solid regions through the original
native terminal solve and `PlanarNonlinearEddy` history. The accepted spatial
potential retains causal skin/proximity memory under finite-amplitude drive
and source-defined magnetic saturation. The existing sparse magnetic operator,
projected consistent P1 conductivity mass, winding/end storage and native
`CircuitDae` remain authoritative; no terminal-impedance realization is fitted.

Harmonic and transient paths share `planar::eddy::solid::SolidWindingRegion`
and the original source normalization. With `S = integral(sigma dA)`, each
whole signed section enforces `integral(J dA) = n I` using
`J_unit = n sigma/S`, `R_DC = depth*n^2/S`, and
`E = n I/S - (delta(A)/dt - mean_sigma(delta(A)/dt))`. One region belongs to
one winding coordinate; explicit return sections, floating conductors and
disjoint stranded copper retain their existing physical roles.

Conductivity refresh updates the source functional and DC drive from the same
original cells while preserving accepted potential, current, circuit and clock.
Actual induction is `F_new^T delta(A)/dt`, plus explicit end-inductance voltage.
An endpoint-linkage change from reweighting a stationary field creates no EMF.
The old potential paired with new profiles is prepared once per terminal solve.
`inductive_voltage_v` reports physical induction; endpoint linkage remains a
diagnostic that can change when the conductor profile changes.

The original concrete field receipt carries actual total local E, J and
physical Joule heat, including uniform/circulating-current cross terms.
`uniform_drive_work_j` is already included in supplied field work and physical
heat; it is never an additional deposit. Stranded heat remains disjoint, and
backward-Euler magnetic Bregman loss remains numerical dissipation. Independent
terminal, field, end-storage and global work gates precede publication. Source,
topology, budget and cancellation refusals preserve accepted histories.
The stationary entry points retain fixed geometry; the Eulerian moving
extension above retains explicit material-transport restrictions.

G0/G1/G3 tests in `tests/nonlinear_winding_solid.rs` and its `slabs` module
compare bipolar diffusion/saturation pulses, native loaded-circuit recurrences,
stationary biased-field conductivity refresh, and two-grid skin/proximity
profiles against independent scalar and continuous-slab oracles. These authored
fixtures are not experimental material qualification; native execution requires
the workspace dependency repositories.



## Bias-dependent winding impedance with resolved eddy currents

`planar::eddy::harmonic::PlanarHarmonicWindings` prepares frequency-domain
small-signal responses from an admitted `PlanarIncrementalWindingBasis`.
The original DC operating point, signed distributed winding profiles, complete
nonlinear differential magnetic Hessian, source temperature, geometry, remanence
and boundary trace remain authoritative. Both longitudinal `dH/dB` and transverse
`H/B` enter the field tangent. Phasors are peaks with `exp(-i omega t)`, matching
`fs-phs::harmonic`; inductive impedance therefore has negative imaginary part.

Every unit winding solves `(K_t - i omega C_sigma) a_j = f_j`. `C_sigma` is the
existing consistent P1 conductivity mass minus each isolated region's weighted
mean projection. It enforces zero integrated axial current independently in
each floating conductor, with original topology/material admission. The sparse
Hessian and matrix-free regional projection feed the existing real-equivalent
block and cancellable FGMRES owners. The original diffusion/storage Jacobi
preconditioner is reused on both real components, with explicit restart, cycle,
field-payload and Krylov-storage limits. Full-mesh operator buffers are reused;
no dense field matrix is formed. Conductor-free frequencies reuse DC derivative
fields only after their true residual passes the current AC allowances.

Original source-paired linkage gives the complete reciprocal terminal matrix
`Z = R_copper - i omega (F^T a + L_end)` and effective inductance `-Im(Z)/omega`.
Every active stranded triangle binds its complete original conductor source;
its resistance integrates `rho J_unit^2` over the actual extruded volume.
Stranded and resolved conductor cells must be disjoint. Explicit end/leakage
storage is retained; end resistance is not inferred. Resolved core/shield eddy
currents produce actual frequency-dependent loss and inductance. Stranded
profiles remain fixed; driven solid sections use the additive path below.
Fitted hysteresis loss and full-wave effects remain outside this reduction.

Acceptance audits complex reciprocity and the COMPLETE mutual power identity
`Z_ij = R_ij + omega^2 a_i* C_sigma a_j - i omega a_i* K_t a_j`, including
declared end storage. Arbitrary phase-shifted terminal currents reuse the unit
fields, then independently reconstruct original-cell B, projected E, J and exact
P1 cycle-mean Joule power. Real terminal power equals spatial copper plus eddy
heat; reactive power equals `-2 omega` times mean quadratic perturbation storage.
DC copper heat is separate. A caller-declared cell-diameter/skin-depth maximum
uses the larger local incremental permeability, including the transverse
direction. This spatial admission bound is not an error certificate.

`tests/winding_harmonic.rs` compares native responses with independently soluble
field modes across DC bias, frequency, conductivity and extrusion depth. It
covers signed mutual coupling, phase-dependent power, floating-current
neutrality, original heat sources, end storage, the no-eddy limit, affine bias
and gauge invariance, transverse tangent/skin depth, double-counting/topology
refusal, work budgets, nonconvergence and cancellation. Manufactured sources
exercise the implementation; they do not validate an empirical material/device.
Native test execution remains a separate requirement where workspace dependency
repositories are available. Finite-amplitude nonlinear cycles, material/geometry
evolution and causal transient realization of the terminal matrix are not claimed.

### Resolved skin and proximity currents in solid windings

`new_with_solid_regions` adds `SolidWindingRegion { region, winding, signed_turns }`
bindings to the same prepared owner. Each whole connected conductive section
has the signed constraint `integral(J dA) = n I_w`. With
`S = integral(sigma dA)`, eliminating its axial potential gradient gives
`J_unit = n sigma/S`, `R_DC = depth*n^2/S`, and actual peak field
`E = n I_w/S + i omega (A - mean_sigma(A))`. The existing sparse magnetic
pencil, projected conductivity operator, current-closure topology, original
source/flux pairing and complex mutual power identity are unchanged.

The retained DC winding profile must already match this source-normalized
distribution. Its volume-weighted absolute discrepancy, scaled by
`max(1, abs(I_DC))`, must satisfy both the caller's weak-load and integrated
current allowances. No operating point or profile is silently replaced.
One conductive region has one winding coordinate and a declared polarity;
every return section is explicit. Other regions remain floating. A winding
may combine solid and stranded sections without sharing source cells.

`Conduction::heat_with_region_drives` adds the uniform drive after the original
conductivity projection and reuses its exact P1 total-field Joule quadrature.
Every driven region checks its actual complex integrated current against `n I`.
Solid-cell heat includes the local uniform/circulating-current cross terms;
separate DC and eddy densities would give the wrong spatial thermal load.
Returned solid, stranded and floating source-cell partitions are disjoint.
`mean_dissipated_power_w()` sums their actual heat once for circuit/thermal
consumers. An unpowered solid winding can still develop resolved proximity
currents while its integrated current remains zero.

`tests/winding_solid.rs` exercises four separated slab sections with two
independent currents against `J'' = -i omega mu sigma J`, with interface data
from Ampere's law. It checks spatial-current refinement, independently
integrated heat, changing neighboring-current phase at fixed winding current,
the uniform low-frequency limit, local Joule cross terms, mixed solid/stranded
sections plus a floating shield, explicit returns, source-profile refusals and
unchanged accepted bias after cancellation. These are manufactured planar
MQS tests, not empirical winding validation. Axial end redistribution, individual
insulated strands inside an authored homogeneous region, full 3-D end turns,
and displacement-current regimes are not inferred.

### Finite-window winding heat and magnetic work

The opt-in `harmonic-exposure` feature adds
`PlanarHarmonicWindings::prepare_exposure` over original frequency owners and
authored peak-current phases. Every tone must retain the exact same prepared
source, geometry, boundary and nonlinear DC linearization. Equal frequencies
combine currents with compensated arithmetic before one native field recovery.
Different frequencies retain their original resolved skin/proximity response.
The compact bank stores numerical fields without copying material bundles;
further observations do not solve a field or circuit again.

Uniform physical windows use the existing `fs-fft::harmonic_window` product
kernel, retaining both beat and sum-frequency terms. Kernels are prepared once
per frequency pair, and double-double accumulators retain interference through
the complete contraction. Original-cell P1 Joule integrals use the actual total
solid/floating conductor E. Stranded-current moments are computed once per
winding and reused for its material cells. Each cell retains its original
source and returns mean watts and interval joules; no cycle-mean subtotal is
deposited a second time.

The same fixed Hessian and explicit end inductances give quadratic perturbation
storage `U = a^T K_t a/2 + sum L_end i^2/2`. Independently contracted terminal
work, spatial Joule heat and `dU/dt` must close within the caller's watt allowance.
Storage change is the integral of its rate, preserving short intervals without
subtracting two rounded endpoint energies. Signed returned source work may be
negative during energy return. Line count, pair contractions and additional
numerical workspace are bounded before allocation; native response budgets
still apply separately. Queries cost
`O((F^2+1)*(vertices+eddy_cells+windings)+stranded_cells)` for F distinct frequencies
and cancel without publishing a partial result. This is an operation bound,
not measured throughput.

Five G1/G3/G4 regression groups in `tests/winding_harmonic.rs`, enabled with
`--features harmonic-exposure`, cover independent real-time field/heat/work
integration, common-period cycle averages, large cancelling equal-frequency
drives, solid winding heat, darkness at nonzero bias, source identity, budgets
and cancellation. The source-backed tests remain unexecuted in this environment
because the workspace's seven sibling repositories are unavailable. These are
settled infinitesimal AC perturbations; total DC heat and DC/AC cross work,
switch-on history, finite saturation cycles, thermal evolution and hysteresis
are not inferred.


### Native circuits drive and load the resolved winding fields

`HarmonicWindingFrequency::circuit_response` connects the same full terminal
impedance to original native R/L/C, mutual-inductance, ideal-transformer and
voltage/current-source constraints. Each field-owned inductor branch is replaced
temporarily by a positive reference resistor and a compensating native current
port. One `HarmonicModel` factorization supplies the open-source and unit-port
responses; a bounded terminal Schur solve and final native replay cancel the
reference resistance exactly. Its value is numerical conditioning, not physical
damping, storage or heat. Existing sparse unit field solutions are reused.

The result retains original branch/node/transformer identities, individual
source outputs, authored control-index gaps and complex phases. Repeated controls
expand to private native ports so their individual currents are not lost; one
control cannot mix voltage and current units. Original external mutual links
must not touch field-owned windings, whose complete magnetic coupling already
belongs to the field terminal matrix. Actual source constraints, transformer
work, Kirchhoff currents, native residuals and winding terminal consistency
are checked independently of the small interface solve.

Spatial solid, stranded and floating-conductor losses are recovered once through
the same field response. Original circuit resistor heat and unrelated native
storage remain separate. Supplied complex power closes against their complete
real and reactive account, using the distinct conjugation conventions for
voltage-driven and current-driven native ports. Reference resistors and private
compensating sources never appear as physical losses or external supplies.
Input graph, native dense payload, interface storage, solve work and cancellation
are bounded before publication; true poles and inadmissible constraints refuse.

Five authored-source regression groups cover loaded coupled windings and
reference-resistance invariance, mixed ideal current/voltage clamps, repeated
source controls, ideal transformers with external mutual RLC at series resonance,
and binding/budget/cancellation/true-pole refusal. The same fixed-bias planar MQS
and source-domain limits apply; native test execution is a separate requirement
where the workspace dependency repositories are available.

## Nonlinear mixed-terminal conduction

`conduction::response::terminals::ConductionTerminalDrive` closes complete
original scalar bulk and explicit interface laws under simultaneous terminal
voltage, conserved total-current or mixed controls. Both planar and axisymmetric
response owners use the original constitutive iteration with one native grouped
field solve at the actual operating state per trial. Floating equipotential
electrodes retain zero net current while transferring current between faces.

Accepted candidates pair the continuing source owner, actual native field and
terminal report. Terminal voltages are physical differences from terminal zero;
source work uses authored current controls and native voltage-source reactions.
Swept geometry, local bulk/contact heat and the native algebraic, constitutive,
current and work gates remain binding. No artificial basis probe, scalable
nonlinear resistance or separate nonlinear solver is introduced.

Every field solve, including discarded constitutive trials, consumes the
explicit total field-solve cap. `ConductionResponseError::FieldBudget` reports
completed work and the next required solve before starting it. Source/domain
refusal, budget exhaustion and cancellation leave the input owner reusable.
The existing scalar memoryless, fixed-geometry and uniform electrically
nonlinear contact boundaries remain.

Focused G1/G3/G4 regression definitions cover simultaneous control modes,
floating-face transfer, actual E/J source samples, new-amplitude continuation,
terminal permutation, common references, swept sector/local heat, nonlinear
surface losses, and exact budget/cancellation retry. A narrow-domain fixture
admits the physical solution while excluding artificial unit-voltage probes.
Pinned syntax checks and independent analytical/bisection references were
performed; these new full native integration regressions remain unexecuted
in the current dependency environment.

## Direct sparse multi-terminal conduction

`PlanarConduction::solve_terminals` and its axisymmetric counterpart solve
simultaneous electrode voltages, conserved total currents or mixed controls
through the original sparse conductivity operator. Complete exterior vertices
of each ideal electrode share one unknown. A floating zero-current electrode
can transfer current between separate physical faces. Unlisted surfaces
insulate; all terminals require one conducting connection through bulk,
original finite contacts or shared electrodes, and each swept electrode needs
its own positive physical area.

Sparse coordinate identification uses the original difference-form couplings.
Internal ideal-electrode terms vanish before summation, preserving small
external connections beside large internal conductances. Signed off-diagonal
entries remain intact. The original checked Jacobi-CG solve runs once at actual
drive; no artificial material-domain probe or dense terminal factor is needed.
The existing reusable conductance-matrix API remains available for frozen
operator response extraction. No wall-clock speedup is claimed.

`ConductionTerminalLimits` bounds terminal count and conservatively admits
extra scalar/index workspace before allocation. Native source-domain,
free-current, total-current and heat/work gates remain active. The returned
`ConductionTerminalReport` retains ordered physical potential differences and
reactions. Independent source work uses authored current controls and native
voltage-source reactions, once each, with double-double physical differences
before rounding. A large common reference remains separate from local fields
and contact jumps. Refusal and cancellation leave the original owner reusable.

Five executed G0/G1/G3/G4 sparse-identification kernel tests cover floating
transfer and energy, signed couplings, strong internal versus small external
conductance, malformed arithmetic and cancellation/retry. Full native
regression definitions additionally cover all control modes, common references,
terminal reordering, exact swept disk reactions, local heat, source work and
admission. Those new native integration regressions have not been executed in
the current dependency environment.

## Axisymmetric nonlinear scalar conduction

`conduction::axisymmetric::response::AxisymmetricConductionResponse` closes
complete original scalar conductivity/resistivity laws at actual meridional
|E| and |J|. Its private shared response engine retains the native swept-volume
operator, exact radial boundary/contact moments, source selections, and separate
algebraic/current/work and constitutive residual gates. Prescribed voltage and
current-density traces, signed total-current electrodes, and regulated-power
drives all solve at actual amplitude within the existing bounded iteration.
Every current/power terminal must have its own positive-area physical trace.

Static/phase contacts retain the native nonuniform P1 surface response.
Electrically varying interface laws bind their original temperature, phase and
electrical coordinates; their endpoint voltage jumps and normal densities must
be uniform within declared gates. A varying trace refuses rather than using a
surface-average constitutive law. Exact surface losses remain separate from
bulk cell powers and enter total power once.

Refreshed bulk/contact seed APIs let an outer owner supply its actual thermal
context while preserving geometry, source identities, coordinate roles and
phase ownership. Fixed or bounded adaptive source-coordinate relaxation carries
through accepted candidates. No frozen linear operator, scalable nonlinear
resistance or planar motion API escapes this facade; continuation retains the
axisymmetric model and its actual-field report together. Rejected domains,
budgets, nonconvergence and cancellation leave the input owner unchanged.

Focused G1/G3/G4 regression definitions cover heterogeneous axial transport,
radial nonlinear refinement against an independent coaxial law, angular-span
scaling, exact axis endpoint loads, local heat, voltage/current/power controls,
strong-law relaxation, source-dependent interface heat and context refresh,
nonuniform-contact and zero-area-terminal refusal, and cancellation/retry.
The new full native integration regressions have not been executed in the
current dependency environment. The same fixed-geometry, scalar memoryless
model boundaries apply; positivity does not establish convergence, uniqueness,
dynamic stability or measured-material qualification.

## Axisymmetric Ohmic conduction and radial current moments

`conduction::axisymmetric::AxisymmetricConduction` reuses the native scalar
conductivity Hodge, checked solve and current/Joule recovery with the exact
`angular_span * r dr dz` measure. Source-bound scalar and oriented meridional
conductivities retain their original material queries and physical units.
Voltage traces and inward current-density boundaries cover actual surfaces of
revolution; the latter integrate `N_i*r` exactly instead of dividing their
total equally. The zero-area axis supplies natural regularity, never an
electrode or prescribed current surface.

`AxisymmetricResistance` retains one native probe solve and recovers actual
voltage, signed total-current or compliant power drives without another
Krylov solve. Local Joule heat, terminal current and source work are checked
at every recovered amplitude. This is the same dissipative circuit component;
the spatial heat must not be credited a second time. Temperature or material
changes require a new source-bound preparation.

`AxisymmetricConductanceMatrix` prepares the reciprocal terminal map with
exactly N-1 native field solves and a bounded, equilibrated dense factor.
Disjoint positive-area surfaces of revolution define ideal equipotential
electrodes; one electrode may join multiple original mesh pieces. Recovered
voltage, conserved total-current and mixed drives use the retained native
operator and source-domain/current/work checks. A zero-total-current electrode
floats while its local surface currents can transfer charge between faces.
`AxisymmetricMixedConductance` factors an arbitrary current-terminal subset
once; repeated controls require no new field solve. Common voltage bias is
removed before arithmetic, and simultaneous local Joule heat includes all
quadratic cross terms. Tests in `tests/conduction/axisymmetric_matrix.rs` use
independent cylinder/coaxial laws, floating-electrode current transfer,
geometry scaling, actual local heat and bounded/cancelled rejection cases.

Fields and anisotropic axes are local radial/axial components with a decoupled
azimuthal direction. A partial angular span is a symmetry sector. This owner
does not add azimuthal current, induction, Hall response,
motion or thermal evolution. Its private core exposes no planar motion API.
The G1/G3/G4 tests in `tests/conduction/axisymmetric.rs` use analytical axial
cylinders, layered logarithmic coaxial resistance, exact disk-current moments,
rotated meridional tensors, metric/axis refusals and cancellation/retry. These
are implementation checks, not measured-material qualification or continuum
error enclosures; executing them remains a separate validation obligation.

## Axisymmetric electrical material contacts

`AxisymmetricConduction::new_with_contacts` and
`new_with_phase_contacts` retain the original ordered interface sources and
bulk material bindings on matching surfaces of revolution. Finite R'' in
ohm m2 produces a resolved voltage jump, local normal current and separate
interface Joule heat. Perfect contact uses shared vertices. The zero-area
symmetry axis cannot supply a physical contact surface.

The shared native surface operator caches its exact radial P1 moments with
the conductivity stiffness. For g=A/R'' and q_i=r_i/(r_0+r_1), its entries are
g/6 times [[1+2q_0,1],[1,1+2q_1]]. Assembly, current and the positive-square
heat integral use these same moments. Radius normalization avoids overflowing
the radius sum. Planar weights retain the previous planar arithmetic.
Resistance and multi-terminal recovery include the contact topology and heat
without a fictitious thin volume or an extra field solve at each drive.

`cell_power_w()` reports bulk loss; `contact_response()` reports each
interface loss; `power_w()` includes both. Thermal coupling must refresh the
original interface state and explicitly partition surface heat between its
receivers. No partition, thermal contact conductance, electrical surface
storage or moving-interface law is inferred.

The focused axisymmetric contact tests cover exact axial disk/annular
resistance, nonuniform radial current and heat, coaxial refinement, source
pressure, orientation, large common bias, matrix recovery, geometry/budget
refusal and cancellation/retry. They use manufactured source laws.

## Three-dimensional Ohmic transport and local Joule heat

`spatial::conduction::SpatialConduction<C>` consumes the original tetrahedral
geometry and source-bearing `SpatialConductorMaterial`. Scalar and explicitly
oriented principal conductivity/resistivity share `d0^T M_sigma d0`, the existing
Whitney weighted Hodge, component-relative voltage gauges and checked CG.
The original geometry is shared across source refreshes; only material weights
and the electrical operator are rebuilt. No magnetic constitutive law is used.

Complete exterior P1 voltage faces and uniform inward-current-density faces
define the boundary. Unlisted faces are insulating; voltage/current faces are
disjoint. Each current-only component must have zero net imposed current before
one deterministic voltage gauge is applied. A current-density trace is not an
equipotential total-current electrode. Source work includes authored trace loads
and solved voltage reactions exactly once, using component-relative potentials.

`resistance[_at_reference_voltage]` and `solve_current[_with_reference_voltage]`
add ideal equipotential electrodes declared as disjoint unions of complete
exterior faces. The total current is signed; the native field determines its
local distribution, including heterogeneous parallel paths and source-backed
contacts. Repeated faces, shared vertices between opposite terminals, internal
contact faces and a pair with no bulk/contact-connected path refuse. An explicit
zero drive retains the finite positive resistance established by a nonzero
reference solve; it does not make disconnected terminals admissible.

`SpatialResistance` retains that original source-bound operator, reference
voltage, relative potential basis and field. Repeated `field_at_voltage` or
`field_at_current` calls perform native field/contact recovery and physical
acceptance without another Krylov solve. Actual signed currents, free-node
residuals, bulk/contact Joule heat and terminal/Hodge work are recomputed;
current control also checks both terminal reactions and authored I*V. The
returned scaled field has no new linear-solve report. Material or contact
refresh requires a new response; no nonlinear, inductive or moving-geometry
extension is inferred. G1/G3/G4 regressions cover series resistance, unequal
parallel current injection, rotated conductivity and SI scaling, separate
contact heating, current reversal/zero, canonical face order and cancellation.

`SpatialConduction::conductance_matrix` prepares a reusable response for N
disjoint equipotential electrodes using exactly N-1 native field solves. The
full reciprocal matrix is the original bulk/contact stiffness energy Gram
matrix, independently checked against terminal reactions and gauge closure.
A grounded equilibrated factor from fs-la handles balanced signed terminal
currents, including a zero-net-current floating electrode that can transfer
current between its faces. Shared electrode admission retains actual bulk,
finite-contact and equipotential connectivity; untouched bodies keep the
original native gauges. Caller controls bound terminal count, basis solves,
retained values, dense work and coefficient residuals.

`SpatialConductanceMatrix` superposes relative potentials and then recovers the
combined E, J and bulk/contact powers through the original field owner. It never
adds isolated basis powers, so simultaneous-drive cross terms and cancellation
are retained. Voltage/current excitation changes reuse the admitted factor;
material, contact or geometry changes require a new response. Native current,
power and cancellation gates remain authoritative. The focused three-terminal
regression checks an independently derived heterogeneous/contact matrix,
floating electrodes, polarity, gauge invariance, local heat and work bounds.
This remains a frozen-source stationary Ohmic response, without induction or
field-dependent constitutive feedback.

`SpatialConductanceMatrix::prepare_current_control` closes mixed terminal
drives on that same response: an arbitrary proper subset prescribes signed
total current while the remaining electrodes prescribe voltage. Its bounded,
diagonally equilibrated current block is factored once and reused. A floating
electrode is the zero-total-current case, with potentially nonzero currents
between its faces. No repeated field solve, local insulation approximation or
independent resistor fit is introduced. The accepted original field retains
bulk/contact Joule heating and all simultaneous-drive cross terms. Separate
actual-state checks enforce the authored currents and terminal work, including
current-source work. Common bias is removed before solving and heat recovery;
returned nodal voltages use the existing terminal-zero reference convention.
The heterogeneous/contact series oracle covers voltage-driven floating
electrodes, nonzero signed injection, several current controls, a controlled
reference terminal, source-work/local-heat closure and cancellation/retry.

Actual cell E=-grad(V), J=sigma E and volume-integrated E.J remain unsmoothed.
Acceptance separately checks free-node/current balance, recovered versus Hodge
power and supplied terminal work. These dimensional residual gates are not
spatial-error certificates. Native synthetic G1/G3/G4 tests cover layered 3D
series conduction under voltage/current/mixed traces, a rotated tensor with
transverse current, SI geometry scaling, source swaps, gauges and cancellation.

The explicit liquid-mass-fraction role admits source-supported confined phase
snapshots in this stationary Ohmic owner. Three principal sources retain their
own exact keys, independent conductivity/resistivity basis and fixed axes.
Existing induction keeps its stricter static-source admission. This adds no
phase mixing rule, Hall law, dielectric storage, induction, contact layers,
moving geometry, field-dependent conduction or measured-material validation.

## Nonlinear scalar source response on original tetrahedra

`spatial::conduction::response::SpatialConductionResponse` closes the original
`ScalarConductorResponse` laws at actual three-component E/J magnitudes. The
shared material query retains all source roles, conditions, pins and support;
ordinary frozen 3D snapshots still refuse unowned electrical coordinates.
Each bounded Picard iteration reuses the original geometry and rebuilds only
the native bulk/contact operator and its terminal response. Boundary traces
and simultaneous voltage/current/floating electrodes use the same field owner,
including actual source work and unsmoothed bulk/contact Joule heat. Nonlinear
operating points are never scaled as fixed resistances.

`with_relaxation(weight)` optionally damps owned E/J source coordinates between
iterations, using the same source-query helper as the planar owner. The default
weight is one; a weight in `(0,1)` requeries the original source at intermediate
coordinates without blending coefficients or moving external context. Actual
source support and all final physical gates remain mandatory. The original
iteration and field-solve budgets still bound work; source evaluation counts
include extra damping queries. Strong positive current-dependent resistivity
and its current-driven conductivity dual are checked against independent
quadratic roots with the entire field in the third spatial component.

Acceptance separately gates conductivity, source-evaluated current density and
the sum of absolute cellwise constitutive power defects, alongside native
current/Hodge/work balance. Every N-1 terminal basis solve is charged across
discarded iterations before the work begins. The immutable candidate retains
applied and evaluated complete sources. Unsupported states, budget exhaustion,
nonconvergence and cancellation leave the original owner unchanged. This is a
stationary scalar, memoryless closure with prescribed temperature/phase context.
Contacts remain frozen unless the explicit local contact-response entry point
below owns their electrical coordinates. Positive conductivity establishes
neither global convergence nor uniqueness. No nonlinear tensor or nonuniform
contact-quadrature law is inferred. Manufactured 3D regressions use independent quadratic series/contact
and current-dependent bar solutions, including third-component fields, floating
and injecting electrodes, polarity, common bias, local heat and source-work
closure, constitutive refusal and accumulated basis-solve budgets.


The explicit `solve_coupled_terminal_candidate` seam rebuilds each frozen
N-electrode response inside this same bounded source closure. An external
driver returns actual terminal voltages and its private proposal; the original
matrix itself recovers E/J, reactions and work for constitutive admission.
No caller-constructed field can substitute for that native recovery. Every
N-1 basis solve is charged before the external callback runs. The continuing
driver must separately charge its descriptor steps and keep all source trials
on one accepted storage history. `into_frozen_parts`/`from_frozen_parts`
transfer and readmit the same complete bulk, contact and relaxation bundle
without exposing the retained linearization as a nonlinear material law.

## Source-dependent resistance on uniform electrical contact traces

Both `PlanarConductionResponse::new_with_contact_response` and
`SpatialConductionResponse::new_with_contact_response` bind the complete
original `ScalarElectricalInterfaceResponse` to each matched native trace.
Every iteration queries each interface at its actual A-to-B voltage jump and
normal current density, rebuilding the same positive P1 contact mass matrix
beside the same native bulk operator. Actual source values are queried at every
trace vertex. Resistance, voltage and current-density defects have separate
dimensional gates; contact constitutive work defects join the absolute bulk
work-defect gate. Surface Joule heat is still the original exact quadratic
trace integral, credited once, without a constitutive-iteration heat source.

Electrically varying contacts are restricted to traces whose voltage jump and
normal density are uniform within the declared contact tolerances. A varying
trace returns `NonuniformTrace`; no face average or surrogate resistor stands
in for nonlinear quadrature. Fixed contact laws retain their original exact
nonuniform P1 response, including in a mixed specimen. General spatially varying
nonlinear contact laws, opening, arcing, surface storage and measured-material
qualification remain outside this rung. Signed and magnitude coordinates retain
their distinct A/B conventions. Temperature and explicitly declared side phase
coordinates may be refreshed by an enclosing owner; electrical relaxation
changes neither of them.

`solve_with_materials_and_contacts_candidate` on the spatial owner retains
geometry, pins, coordinate roles and contact pairing while closing explicitly
refreshed thermal contexts. Planar `ConductionResponseSources` transfers bulk
laws, contact laws, their controls and relaxation together through a frozen
operator handoff. Missing contact bindings cannot turn an electrical response
axis into a static source. Analytical layered-joint regressions check both
voltage/current drives and signs, temperature refresh, current continuity,
bulk/surface work, policy retention, nonuniform refusal and cancellation/retry.

## Ionized-gas current and Joule heating

`PlanarConduction` accepts the sealed `PlanarOhmicResponse` material response,
including `fs-material::ionized_gas::IonizedGasConductivity`. The same weighted
FEEC assembly, boundary current solve, source-work audit, cell current and
`E.J` heat integration are used for source-property conductors and ionized gas
with actual frozen populations. Gas populations, separate electron/heavy
temperatures and collisional rates remain in the original
material snapshot. Contact construction and material refresh keep their
stronger source-material requirements; no gas/solid contact law is inferred.

Field recovery checks the material's admission against the actual electric
field before publishing current or heat, including circuit voltage recovery.
This prevents a weak-drift gas law from silently extending into breakdown.
Resistance and total-current preparation accept an explicit positive reference
voltage, so a valid small physical drive need not admit a one-volt computational
basis. The retained reference field is admitted at its actual voltage;
`G=P/Vref²=I/Vref` supplies resistance, and circuit recovery scales once by
`V/Vref` before independently checking the actual field, current and power.
The gas response is static, unmagnetized and local, with a prescribed massive
neutral collision background. It does not add a discharge, Hall, radiation,
space-charge or magnetohydrodynamic solver. Thermal, flow or finite-rate
chemistry feedback must carry its actual new populations and temperatures into
a new accepted snapshot and resolve collision rates there. A kinetic endpoint
is never projected back onto LTE merely to obtain conductivity.

The material's checked Joule partition consumes the actual recovered FEEC
electric field and supplies separate electron/heavy thermal powers. Each
power density uses that field cell's physical volume before entering the
conservative thermal-pool update. Their sum is the original field's `E.J`;
the native resistor heat is not added a second time. Updated temperatures feed
the same source-owned material and field owners.

G1/G3 tests in `tests/ionized_gas.rs` use analytic H/H+/electron thermodynamics
and explicit synthetic collision sources, then solve the actual native mesh.
Affine voltage, current density, Joule power, native circuit resistance and
voltage reversal agree with independent slab formulas at two temperatures.
G0 tests refuse actual PDE and circuit fields above the carrier drift domain.
Further G1/G3 tests carry an accepted finite-rate ionization endpoint through
actual neutral-density momentum drag, field current, Joule power and native
circuit resistance, and vary electron temperature with populations held fixed.

`tests/ionization_feedback.rs` carries separate-temperature monatomic atomic
ionization through that same path. Original NIST event rates and source-derived
atomic detailed balance change carrier amounts; native cell-volume Joule powers
are held explicitly during the chemistry interval. Independent monatomic
calorics audit electrical work, changing electron heat capacity and chemical
storage once. Refreshed populations and collision rates then change the native
current and circuit resistance. This is an ionization substep with a subsequent
field feedback evaluation, not a continuously fixed-voltage transient or an
empirically validated discharge; elastic exchange remains its separate owner.

`tests/coulomb_feedback.rs` uses the same native field/circuit path with
reciprocal Coulomb momentum transfer. The material resolves correlated carrier
velocities, retaining both its original background drag and every charged pair.
Original FEEC cell volumes convert its friction heating and signed
drift-dependent thermal correction to two held powers. The existing gas energy
owner consumes the identical population-bound Coulomb source together with the
original neutral momentum source, refreshes their zero-drift exchange at actual
temperatures, and closes both pool budgets against one electrical-work input.
Refreshed two-temperature material changes native current and resistance.
This remains an explicitly held-power energy interval followed by field
feedback, without a second PDE, a free plasma common-velocity gauge, or an
additional `E.J` charge to either thermal pool.

`tests/collisional_ionization_feedback.rs` combines the two paths in one
reacting energy evolution. Both an ionizing parcel and a recombining parcel
consume the original atomic impact/reverse law, neutral momentum source and
Coulomb screening policy while carrier amounts and both temperatures change.
Independent monatomic budgets include changing electron heat capacity,
original chemical storage, both internal exchange channels and the one held
native electrical-work input. The final populations invalidate the original
Coulomb transport binding; an explicit refreshed model retains its original
screening and domain limits. Tests check actual neutral-density drag, compare
the final coupled conductivity with an independent two-fluid inverse and a
counterfactual holding the old populations, then verify native field current,
Joule power and circuit resistance. No equilibrium projection or additional
material/PDE owner is introduced.

## Prescribed Eulerian moving-conductor induction

`PlanarEddy::step_with_eulerian_velocity_candidate` continues the original
accepted vector potential, impressed-current history, floating electrical
regions and clock with a prescribed continuous P1 planar velocity. It solves
`K A_mid + P_sigma(Adot + v.grad(A_mid)) = load(J_impressed_mid)` on the same
fixed mesh and source-bound linear material snapshots. The existing
conductivity-weighted P1 mass/mean projection enforces zero integrated axial
current separately in every isolated conducting region. Uniform motional EMF
is cancelled by the floating scalar-potential gradient, without fictitious
eddy current or Joule heating.

Velocity feedback makes the operator nonsymmetric. The existing cancellable
FGMRES owner supplies the recurrence, with a positive diffusion/storage Jacobi
preconditioner and a correct matrix-free transpose. Full-mesh work buffers and
cell geometry are reused during Arnoldi actions; neither the advection matrix
nor regional projection is densified. Explicit restart/cycle/Krylov-storage
budgets bound the solve. Acceptance recomputes the actual retained increment's
true Euclidean residual against the larger of caller-supplied absolute [A m]
and relative targets. A caller-declared cell magnetic Peclet limit
`mu*sigma*max_vertex_speed*diameter/2` bounds admission of this unstabilized
Galerkin spatial discretization; it is not an error certificate.

Returned original-cell fields distinguish laboratory E from conductor-frame
`E'=E+v cross B`, with `J=sigma E'`. The SAME stage B and P1 current supply
`J cross B`, exact quadratic Joule heat, and exact P1 Lorentz-mechanical work.
The accepted account is `delta W + Q_Joule + W_mechanical = W_impressed`;
positive mechanical work is delivered by the field to prescribed motion.
Thermal and mechanical receivers must close their own balances with these
reported transfers. These are free-conduction-current forces, not an additional
Maxwell-stress body load or a magnetization/interface-force model.

This is a nonrelativistic, fixed-geometry, fixed-coefficient Eulerian MQS lane.
Motion must preserve the represented material occupancy/coefficient field,
including homogeneous through-flow with separately supplied boundary data.
Translating finite bodies, moving material interfaces or magnetic sources,
ALE/sliding meshes, electrical terminals, nonlinear magnetic laws, remanence,
hysteresis and displacement current remain outside it. Prescribing velocity
does not advance mechanical state or infer its energy supply. The stationary
and temperature-refresh operations remain separate. Physical convention and
fixed-domain applicability follow the official COMSOL AC/DC guide,
[Velocity (Lorentz Term)](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.040.html).

Focused G1/G3/G4 tests compare the spatially varying-velocity field with an
independently assembled nonsymmetric two-mode system and separate degree-two
quadrature for heat and mechanical work. They cover source/material/depth
substitutions, genuine field feedback, uniform-field charge cancellation,
zero-velocity stationary continuation and work/source/velocity/cancellation
refusals. Synthetic coefficients do not qualify an empirical material/device.

## Purpose and layer

L3 material-resolved electromagnetic fields. The current owner is a linear,
isotropic, fixed-state planar magnetic vector-potential problem:
`-div(mu^-1 grad(A_z)) = J_z`, `B = (d_y A_z, -d_x A_z)`,
`H = B/mu`. Coordinates and extrusion thickness are SI metres, potential is
Wb/m, impressed current is A/m2, and permeability is absolute H/m.

## Owners, public API and invariants

`PlanarMagnetostatics` owns an immutable `TriComplex2`, one complete resolved
material state per triangle, and the existing FEEC stiffness. Permeability
uses `fs-material::field_storage::linear_magnetic_permeability`, the exact
admission shared with source-derived circuit inductance. Unsupported
field/frequency dependence or operating-field bounds refuse. A different
material or temperature snapshot requires a new prepared problem.

The existing weighted Whitney one-form Hodge and incidence assemble
`d0^T M(1/mu) d0`. Impressed cell current enters exact P1 load moments.
Prescribed potentials must cover complete exterior edges and anchor every
connected component. Remaining exterior edges have zero natural tangential H.
The existing checked CG owner supplies the recurrence and a recomputed true
Euclidean residual; exhausted convergence never publishes a field.

`PlanarField` borrows its exact prepared owner and retains nodal potential,
cell B/H, volume-integrated constitutive energy and the solver report. Induced
cells use B.H/2; recoil cells retain the fixed-remanence reference below. Fields remain
separate on each side of a material interface. Constant local gauge offsets
are subtracted before differentiation.

`PlanarFrame` admits a finite right-handed orthonormal 3-D placement.
`CellSample` names an exact source triangle and consumer position in metres.
`sample_cells` validates barycentric and extrusion support before expressing
B in the consumer's global frame. It never guesses a nearest cell. The output
feeds the existing LBM moving-conductor/current/heat/force owner.

## Errors, work and execution

`MagneticError` distinguishes shared material/assembly refusals, unsupported
inputs, cancellation and an actual nonconvergence report. Assembly limits,
explicit sample limits and solver iteration controls bound work. Polls occur
through source admission, assembly, load construction, component traversal,
Krylov boundaries, field recovery and transfer. Failure leaves immutable input
owners unchanged and publishes no partial output.

Traversal and arithmetic use fixed mesh/source storage order and the existing
numerical owners. No cross-platform bitwise identity or certified roundoff
enclosure is claimed. There is no unsafe code or feature flag.

## Conformance and boundaries

Synthetic regressions cover affine fields and energy, layered permeability and
tangential-H continuity, impressed-current scaling, manufactured-field
refinement, temperature-source refresh, component/boundary admission, source
model refusal, cancellation and explicit geometry/frame transfer. The actual
two-step LBM regression changes only one region's permeability and checks
computed B, current/Joule heating, Lorentz force and evolved momentum/thermal
populations. Sources parse with the pinned Rustfmt; native compilation and
regression execution remain separate verification requirements.

This is a discrete P1 estimate on a finite planar domain with fixed linear
material state. It does not model axisymmetry, open/infinite boundaries,
nonlinear or hysteretic ferromagnetism, irreversible demagnetization, magnetic loss,
spatial eddy-current/flow induction feedback, distributed magnetostriction or
moving-geometry field-work coupling. The bounded static net-force connection
is described below; it is not a dynamic electromechanical energy integrator.
The existing flow connection treats the solved field as imposed in its
low-induction regime. Point transfer is not a conservative remap and does not
certify an arbitrary target grid's discrete magnetic divergence. Numerical
fixtures are not empirical material or device validation.

## Static magnetic load to native mechanics

`PlanarField::maxwell_wrench` integrates the computed Maxwell stress
`T = B tensor H - (B.H) I/2` over the closed boundary of explicitly selected
mesh triangles. Those triangles include the target and an exterior collar.
`TriComplex2::trace_for_faces` remains the topology owner. Both cells incident
on every contour edge must share one complete resolved exterior material
snapshot, use an induced-linear law without remanence, and have zero original
impressed current. Enclosed recoil magnets remain supported by the original
field solve and its source-envelope operating checks. The contour cannot touch
the outer mesh boundary or run along a material interface; incident triangles
must also lie on opposite geometric sides rather than form a folded embedding.
The original static solve retains its actual cell-current field for this admission;
induction-history recovery alone cannot satisfy it. Incremental responses on
magnetized problems also refuse because they omit the fixed bias and carry no
actual recoil operating-point admission. Unbiased current responses remain usable.

The model assumes a homogeneous linear exterior medium without mechanical
pressure response or magnetostriction. Its permeability comes from the same
material query as the field. No material name selects a constant or force law.
Contrasting enclosed permeability changes the actual field and thus the load.
The contour's exterior-side P1 B/H values remain piecewise constant: edge
force and first moment are exact for those discrete values, while the physical
force remains a mesh- and contour-dependent numerical estimate. No gradient
recovery, interfacial averaging or error certificate is implied.

Actual extrusion thickness multiplies the surface load. Geometric triangle
orientation determines outward normals for either clockwise or
counterclockwise complexes. Force and torque use the caller's `PlanarFrame`;
the torque origin is explicit in that same frame. Force and first moment are
reduced about a nearby mesh reference before one final frame/origin transform.
Common frame/origin translation is subtracted before adding the local reference,
so it cannot erase the contour or half-depth lever arms. The returned `fs_ga::Wrench`
is the canonical native articulated-dynamics input, not another force type or
mechanical integrator. End-face contributions cancel in the invariant planar
model. End fringing and distributed magnetostriction remain unsupported; the
electric owner below uses the same geometric integration. The fictitious contour is not a physical body
surface, so its traction must not be used as an inferred local solid load.

This is an instantaneous static connection. Motion requires updated geometry
and field solves with explicit source, field and mechanical work accounting;
it cannot reuse a frozen field as a complete energy-conserving actuator. An
overlapping J-cross-B load must not also be applied to the enclosed body.
Original mesh-cell and edge limits bound contour traversal/allocation. Fixed
edge order and compensated summation preserve deterministic selection-order
behavior. Cancellation and admission failures leave the original field intact.

Focused G1/G3/G4 regressions exercise zero force in a uniform field, refinement
toward an independent Lorentz-force integral, thickness and source scaling,
source-only permeability changes propagated through native free-body
acceleration, frame/origin transformations, and unsupported contour/history
refusals. These synthetic calculations do not validate an empirical material
or device. The stress convention follows J. L. Kirtley, MIT 6.685, Class Notes 1,
section 4.2: https://web.mit.edu/6.685/www/chapter1.pdf.

`PlanarField::virtual_work` and `PlanarNonlinearField::virtual_work` also
differentiate the original static discrete magnetic potential along a complete
P1 vertex displacement. Each cell's integrated impressed current, original
Dirichlet A_z values, material state and extrusion thickness remain fixed.
The returned mechanical work is `-d(Pi)/ds`, with
`Pi = sum(W_cell) - sum(J_z*volume*mean(A_z))`. A unit rotation tangent of a
body, extended through the surrounding mesh, extracts its axial torque. At
zero Dirichlet A_z the work is the fixed-current coenergy derivative; without
impressed currents it is minus the fixed-flux energy derivative. Mixed controls
retain the full stationary potential rather than substituting either shortcut.

The native A_z pullback gives `W'_cell div(u) - volume H.(grad(u) B)`.
Reversible nonlinear materials use their original `W'=B.H-W` and H, not a
secant-linear replacement. The contraction shares the allocation-free
geometric kernel with electric virtual work; local displacement differences
annul rigid translation before compensated accumulation. Work costs linear
time in original vertices/cells and requires no field re-solve. Recoil magnets,
hysteresis and induction-recovered fields refuse because no static material
motion/source convention is inferred. This is an instantaneous derivative,
not finite-step electromechanical energy accounting or magnetostriction.
Analytic scaling/rigid-motion tests and independently re-solved geometry
perturbations exercise current-only, flux-only and mixed controls, nonlinear
coenergy, nonzero rotor torque, cancellation and original-state retention.

The separate `rigid_magnet_virtual_work` methods admit source-backed recoil
magnets under an explicit rigid material-frame convention: Br rotates with
the original magnet cell while its magnitude, recoil slope, material state
and energy reference stay fixed. Every recoil cell must have a skew P1
displacement gradient, checked against the local contraction's roundoff scale
without an absolute strain deadband. Connected magnets move rigidly; separate
magnets may have different motions, including zero. The complete displacement
declares torque origin, angle and exterior extension. To extract one rotor's
torque, other physical source regions stay stationary and the moving extension
lies in their current-free surroundings.

The original recoil storage `|B-Br|^2/(2*mu)` contributes
`-volume H.Omega(B-Br)=0` because `delta Br=Omega Br`. The implementation retains
this exact shape/source cancellation; freezing Br in the spatial frame would
produce a different, incorrect load for the declared material motion. Native
nonlinear surroundings still contribute their actual coenergy. No alternate
recoil energy reference or formation energy is introduced, and negative recoil
cell coenergy is allowed. Shear/dilation of any magnet, incremental responses
and induction histories refuse. This is static virtual work: finite movement
requires rotated geometry AND source directions followed by original field
solves and source-envelope admission. Independent finite-rotation tests cover
alignment torque, saturated surroundings, fixed stator currents, a frozen-Br
negative control, rigid-motion objectivity and atomic refusals. They establish
implementation consistency, not empirical magnet/device validation.

## Distributed winding and circuit connection

`winding::PlanarWindingBasis` admits signed trianglewise `J_z / terminal_current`
profiles [1/m2], including turn count and outgoing/return orientation. Profile
integrals use cross-sectional area; flux linkage and energy use the mesh's
actual extruded volume. Each profile must close separately on every connected
mesh component to the caller's current-balance tolerance. Closure residuals
remain available. Unit-current fields use homogeneous prescribed potential
and the existing planar boundary admission and checked solver.

The active-section inductance is the symmetric field-energy Gram matrix
`L_ij = integral(B_i.B_j/mu dV)`. Every entry must agree with the original P1
`integral(A_i J_unit,j dV)` pairing within an explicit absolute tolerance.
Original fields, material snapshots and solver reports remain reachable.
Callers separately declare nonnegative series end/leakage self-inductances;
zero explicitly neglects those ends. These entries contribute to total flux
and stored energy, but not to the represented cross-section field. End mutual
coupling is unsupported. Singular total inductance refuses; no numerical
leakage is invented. The existing Cholesky owner constructs the inverse, whose
measured `L L_inverse - I` row-sum residual must pass the declared threshold.

`winding::circuit::WindingCircuit` binds every original graph inductor exactly
once through explicit branch indices. Positive winding current follows that
branch's first-to-second terminal orientation. Native `delta_lambda = L i` replaces
placeholder branch inductance; the full inverse replaces only the original
flux storage block. Capacitor charges, algebraic rows, Kirchhoff stamps and
source ports retain their original correspondence. The existing `CircuitDae`
owns consistent initialization and audited evolution. Circuit-state currents
reconstruct the same field through the existing triangle/position/frame
transfer; resistor branches continue to own dissipation.

Extraction admits at most 16 windings with explicit basis-storage and numerical
budgets. The dense circuit adapter admits at most 32 nodes, 64 branches,
16 transformers, 32 source ports and 128 extended states. Cancellation or
refusal publishes no partial owner. Regressions cover actual mutually coupled
circuit evolution, an interleaved capacitor and reversed branch mapping,
source-work/heat/storage accounting, material-only substitution, winding sign
and turn scaling, extrusion depth, end storage, gauge closure and refusals.
Native test execution remains a separate verification requirement. This lane
keeps geometry and material state fixed during evolution; it does not infer
nonlinear magnetics, motion, skin/proximity losses, changing-material work or
three-dimensional winding-end fields.


## Differential winding inductance at nonlinear bias

`winding::incremental::PlanarIncrementalWindingBasis` first solves the actual
nonlinear DC field from signed spatial winding profiles and terminal currents,
with the original material laws, fixed recoil sources and prescribed potential
trace. The operating field must pass the native nonlinear residual gate.
Its source-derived constitutive Hessian is then assembled once and reused by
checked CG for every unit-current derivative; perturbations are zero on exactly
the original prescribed vertices. The mesh and tangent remain sparse.

The resulting symmetric terminal matrix is `L_diff = d(lambda)/di`, including
declared linear end self-inductances. It is not `lambda/i` or `2W/i^2` in a
saturated or biased core. Both radial `dH/dB` and transverse `H/B` tangent
eigenvalues are retained. Every entry of the Hessian-energy Gram matrix must
agree with the original `integral(dA/di_j * J_unit,i dV)` pairing within the
explicit henry allowance. The existing inverse-residual gate admits the small
terminal inverse. The original operating field, actual absolute flux, signed
profiles, boundary, source identities and true-residual reports remain reachable.
Voltage evaluates `L_diff di/dt` with all mutual terms and no voltage from a
fixed remanent flux offset.

Source H(B) corners refuse within a caller-declared tesla exclusion band;
maximum-flux and recoil-envelope edges refuse because a two-sided derivative
is not supported there. The band is not a field-error certificate. This is a
local, reversible, fixed-temperature and fixed-geometry DC derivative, not an
arbitrary-amplitude circuit law, changing-temperature work law, hysteresis
model, skin/proximity impedance or three-dimensional end-field calculation.

G1/G3/G4 regressions compare saturation with an independent two-mode analytic
solution, including the full transverse tangent under nonzero boundary bias;
check signed mutual/turn/depth scaling, explicit end storage and inverse action;
compare independent nonlinear current perturbations; retain permanent flux
without spurious induced voltage; and exercise source corners, field bounds,
work limits, singular modes, nonconvergence and cancellation. The formulation
follows [ANSYS Theory Reference §5.5](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_thry/thy_emg6.html).
These manufactured sources do not qualify an empirical core or device.

## Spatial material winding resistance and native-circuit heat

`winding::ohmic::OhmicWindingCircuit` binds every active original magnetic
triangle to its complete `ResolvedConductorState`. The same signed unit-current
profile and actual extruded mesh volume give each winding's intrinsic series
resistance, `R_w = sum(rho_e * J_unit,w^2 * volume)`. Conductivity and resistivity
use the canonical scalar conductor owner; no metal names or preset coefficients
enter the solver. Current-carrying windings must occupy disjoint triangles;
unresolved mixtures/fill factors and skin, proximity, eddy and end-turn losses
are not inferred.

Resistance is stamped on the mapped native flux coordinates, conjugate to the
full mutually coupled `L^-1 delta_lambda` currents. All original magnetic/capacitive
storage, Dirac structure, ports, ideal constraints and numeric graph resistors
are retained. The existing descriptor midpoint step remains the only electrical
time integrator. Per-triangle `E.J` heat, every graph-resistor heat, integrated
series `I^2 R` work, native dissipation and total stored-energy/source-work balance
must agree within explicitly supplied absolute power/energy allowances. These
checks do not certify temporal or spatial approximation error.

Each candidate begins with the original accepted fluxes/charges and time; source
refresh, nonlinear thermal trials and cancellation cannot debit accepted history.
A refresh retains the complete original card, coordinate schema, requirements,
selected source claims and independent electrical source basis. Only the owned
absolute-temperature coordinate may change. Fixed external conditions and the
prepared magnetic permeability stay fixed. A temperature-dependent permeability
source refuses changed temperature even if its endpoint values coincide, since
fixed-L evolution cannot account for magnetic parameter work along that path.
Data-only material replacements are prepared as a new magnetic/electrical model.
`OhmicWindingCircuit::rebind_material_candidate` transfers accepted history to
that model at the same clock. It retains the exact mesh/metric, original
homogeneous magnetic boundary, signed current profiles, end-inductance model,
branch correspondence and all nonmagnetic graph data. Every capacitor charge is
preserved bit-for-bit. Total physical linkage is conserved by translating native
incremental linkage as `delta_lambda_new=delta_lambda_old+(lambda0_old-lambda0_new)`;
identical bias retains native linkage bit-for-bit without adding and subtracting
the larger permanent flux. The native consistency solver re-equilibrates
algebraic voltages/currents. Incompatible ideal constraints refuse instead of
creating an unaccounted flux or charge impulse. Both actual endpoint fields must
satisfy their source recoil envelopes before a replacement is returned.

Replacement adopts the new full mutual-inductance block and spatial conductor
sources, including a different admitted resistivity/conductivity source basis.
The signed external work is the change in `W0 + delta_lambda^T L^-1 delta_lambda / 2`
at fixed total linkage. Its incremental part is independently compared with the
native field/capacitor storage change before adding the bias reference, so large
W0 cannot mask electrical-accounting errors. The full signed work is checked
against the caller's absolute work bound. The report retains before/after bias
and total field/capacitor energy, current, resistance and source-state identities.
This is work in the declared ZeroFlux/FixedRemanence field potential; it does not
predict material formation energy, an intermediate replacement trajectory,
hysteresis or irreversible demagnetization. No time, source-port work or Joule heat is
advanced. Subsequent original midpoint steps and field sampling use the new
material basis; resistance-only replacements have zero instantaneous storage work.
Source, geometry, ideal-constraint, work and late-cancellation refusals publish
nothing. G1/G3 regressions in `tests/winding_ohmic/rebind.rs` cover analytic mutual
current/work, nonzero-time continuation with actual new spatial heat, reverse and
zero-work replacements, and incompatible geometry/constraints.

This is externally supplied material reconfiguration, not continuous thermal
evolution, motion, hysteresis, mass transfer or a caloric heat debit. Changing
geometry or automatically closing magnetic parameter work into thermal storage
still requires its own coupled owner.

## Temperature-refreshed spatial induction

`PlanarEddy::step_with_materials_candidate` retains the original mesh, A history, fixed boundary traces, electrical region identities, source choices and clock while assembling source-pinned endpoint permeability. It reuses the original checked CG with endpoint magnetic effort (K_old A_old+K_next A_next)/2; insulating/floating algebraic constraints are checked at the actual endpoint. The accepted candidate owns its rebuilt operator through a shared immutable owner; a field recovered from that candidate borrows the candidate lifetime.

Each triangle reports actual old/new curl A, magnetic free energy, volume and signed parameter work P=V B_old dot B_next (1/mu_next−1/mu_old)/2. The native account is delta psi + Joule heat = impressed-source work + P. Parameter work is not automatically heat: a thermal or external-work owner must close it. Fixed-coefficient steps retain P=0. No constitutive loss is added to resolved sigma E². Fixed geometry and linear isotropic nondispersive mu remain required; moving media, hysteresis and unconditional coupled entropy/passivity are not claimed. Synthetic tests cover analytical one-mode evolution, repeated operator/history continuation, algebraic insulating refresh, unchanged-source refusal, assembly budget and cancellation.

## Source-resolved reversible permanent magnets

`PlanarMagnetostatics::new_with_material_models` accepts induced-linear and `PlanarRecoilMaterial` cells on the same mesh. Each magnet supplies a complete `ResolvedRecoilMaterial` and an explicit unit direction in the original planar frame. The same FEEC stiffness and checked CG solve impressed current together with the exact P1 weak load `integral (Br/mu_recoil).rot(grad(phi)) dV`. Retaining this whole moment includes material-interface and natural-boundary source contributions even when Br is constant inside a cell. Unprescribed exterior edges continue to impose zero tangential physical H.

Recovery reports `H=(B-Br)/mu_recoil`, source-envelope operating margins, and cell storage `V*(B-Br).H/2`, with the `FixedRemanence` reference explicit. Induced cells retain `V*B.H/2`. A source-envelope exit refuses the complete field, including fully prescribed solutions. This is fixed-state constitutive storage, not a force, hysteresis-loss, demagnetization-history or live parameter-work result. Existing explicit frame/triangle sampling consumes the resulting B. Fixed-remanence winding response requires the explicit affine constructor described below. Stationary midpoint and Eulerian moving-conductor histories retain the original fixed Br source and its shifted storage reference as described above.

The recoil convention follows [COMSOL's remanent-flux constitutive relation](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.032.html); the shifted fixed-state energy reference follows the linear-recoil discussion in [FEMM's permanent-magnet energy note](https://femm.info/doku/doku.php?id=pmenergy). `material_permanent` adds manufactured-source regressions for constant-Br natural boundaries, current-plus-remanence loading, magnet/air interface fields, direction reversal, material/temperature swaps, source/domain refusals and consumer guards. These are analytic model checks, not measured catalog-material validation.

## Fixed-remanence affine winding and circuit response

`PlanarWindingBasis::new_with_fixed_remanence` uses the same linear operator for a source-admitted zero-current field and explicit incremental unit-current responses. Its additional stationarity allowance is in Wb and checks `integral H0.B_unit dV`; the inductance-pairing allowance remains in H. Retained-basis capacity includes the bias field. The generic entry also accepts entirely induced-linear meshes as zero bias.

Total linkage is `lambda=lambda0+L*i`. Native `CircuitDae` coordinates remain `delta_lambda=L*i`, including declared end inductance; initialization and current recovery do not interpret permanent flux as current. Native magnetic/capacitor energy is relative to the zero-current state. `bias_energy_j` retains its constant field reference W0, and `stored_energy_with_bias_j` adds it for reporting without subtracting a large constant during step accounting. Changing fixed Br data changes bias flux/field/storage without inventing electrical power.

`field_for_currents` reconstructs `A0+sum(i_w*A_unit,w)` and uses original physical recovery to admit every recoil cell. It also superposes the original static current snapshots, allowing this actual field to use the existing instantaneous Maxwell-force owner; arbitrary induction recovery still lacks that provenance. Sampling, checked circuit initialization/`step_audited`, and Ohmic candidate endpoints share the recoil admission. Incremental responses are explicitly marked and carry no actual magnet operating-point claim. Envelope exits/cancellation return no accepted replacement. Endpoint checks do not bound the continuous trajectory between steps. The raw DAE exposes incremental coordinates and carries no material-domain policy itself; fixed-remanence callers use the checked wrapper.

An active conducting recoil triangle must reuse its exact original complete material bundle, including its electrical claim. This admits unchanged recoil/conductivity sources without generic permeability and refuses temperature, remanence, slope, limit or source refresh that would leave stale fixed bias. Ordinary induced-material electrical refresh retains its previous policy. Geometry, material conditions and remanence stay fixed; no motional EMF, hysteresis or thermal parameter work is supplied. Manufactured-card `winding_permanent` regressions cover affine linkage/storage and layered inductance oracles, native RL/Ohmic work, mutual branch mapping, data swaps, source-envelope refusal and fixed conducting-magnet source binding.

## Material-driven reversible nonlinear planar fields

`planar::nonlinear::PlanarNonlinearMagnetostatics` admits a linear or reversible
H(|B|) source per original triangle. It shares mesh/metric, exterior-edge and
component-gauge admission, shape gradients, weighted Whitney Hodge, incidence,
sparse assembly and checked CG with the existing planar owner. It assembles
the actual weak residual `integral H(curl A).curl(v)-J v`, and the symmetric
tangent with radial eigenvalue dH/dB and transverse eigenvalue H/B. The tensor
is rotated into the gradient frame before the existing degree-1 Hodge assembly.

The complete initial field must satisfy prescribed potentials and source
envelopes. Newton uses explicit iteration/linear/backtracking limits and an
energy-decreasing Armijo search; trial fields outside the original H(B) range
are rejected, never extrapolated. Acceptance recomputes the free-coordinate
nonlinear Euclidean residual [A m], independent of Krylov convergence. The
reported target is the maximum of the absolute allowance and the relative
allowance times max(initial residual, impressed-current load). Exhaustion or
cancellation returns no partial field. Material, geometry and source owners
are immutable throughout; traversal and tie-breaking use fixed storage order.

Returned B/H, per-triangle energy/coenergy, total storage and actual work counts
remain tied to the prepared owner. Nonlinear energy integrates the source H(B),
not B.H/2. `sample_cells` reuses the original explicit triangle/frame transfer,
so nonlinear material changes reach existing imposed-field flow consumers.
This is a finite-domain fixed-geometry P1 estimate, not a discretization error
certificate, nonlinear winding inductance, hysteresis, eddy-current trajectory,
force/coenergy motion derivative, thermal parameter-work or remanence model.
The same reference conventions are described in the
[FEMM manual](https://www.femm.info/doku/lib/exe/fetch.php?media=upload%3Afiles%3Amanual.pdf),
magnetic energy/coenergy discussion. Synthetic tests do not validate a catalog
material. Focused G1/G3/G4 tests cover exact discrete saturation, source-only
substitutions, oblique mixed-material interface continuity, energy derivatives,
the linear limit, source envelopes, bounded nonconvergence and cancellation.

## Spatial dielectric electrostatics and field-derived circuit capacitance

`electrostatic::PlanarElectrostatics` retains one complete material state per
original triangle. Its absolute SI permittivity is admitted by the same
`linear_electric_permittivity` law used by the uniform `ResolvedCapacitor`.
The existing weighted Whitney Hodge, incidence, sparse assembly and checked CG
owners solve `-div(epsilon grad(V)) = rho_free`; no magnetic material or fictitious
conductivity substitutes for a dielectric. Fixed external temperature, pressure
and composition remain attached to their original source queries. Frequency,
field, strain and stress-dependent or applicability-bounded laws require an
appropriate constitutive model and refuse here.

Exterior voltage traces cover complete mesh edges and anchor every connected
component. `solve_with_sources` admits explicit uniform outward D.n on other
exterior edges, uniform free charge on internal edges, and complete fixed
impressed-polarization vectors in the original mesh frame. Unlisted exterior
edges carry zero normal displacement. Surface moments use the actual extrusion
area and distribute sigma*area/2 to each endpoint; outward D contributes with
the opposite weak sign. Internal sheets and exterior displacement traces are
distinct and cannot be interchanged or duplicated; displacement cannot occupy
a voltage edge. Endpoint/source permutations preserve canonical integration.
Recovery reports cell `E=-grad(V)`, `D=epsilon E+P` and
`integral E.epsilon.E/2`, without smoothing material
interfaces. Prescribed-node weak reactions `K*V-f` are electrode free charges;
their sum plus volume and sheet charge must equal imposed outward displacement.
The RHS includes `integral grad(N).P`; its voltage pairing is
`-integral E.P`, exposed as a separate full signed physical pairing. Fixed P
is additional to the induced response already present in the material law.
Free-node residuals and `U - V.(all_source_loads + electrode_charge)/2` have independent caller-supplied
absolute charge/energy gates in addition to the checked linear residual.
These finite-element balance checks do not certify spatial approximation error.

`capacitance` solves two disjoint complete terminal traces at one and zero volts
with no volume charge. The source-bound unit-voltage field gives `C=2U`, checked
against its independently summed positive-terminal charge. Its `circuit_branch`
feeds the existing Kirchhoff/pHS capacitor with that geometry-derived value;
the original field object retains all source and geometry identities. Changing
only dielectric data or a supported fixed temperature changes both the spatial
field and native circuit storage/decay. Live coefficient changes still require
parameter-work accounting, and are not implemented by overwriting a branch.

The model is fixed-geometry planar, linear, nondispersive and lossless, with
scalar or explicitly oriented SPD dielectric laws.
Fringing within the specified domain is solved; missing exterior/end geometry,
leakage/conduction, dielectric loss, breakdown, electrostriction, intrinsic
polarization energy, changing-P source work and evolving thermal work are not
inferred. Maxwell contours exclude source-bearing collars while allowing
enclosed polarized/charged bodies. Source-free capacitance and retained-charge
owners keep their original material law and storage semantics.
`tests/electrostatic.rs` covers analytic layered
dielectrics, volume-charge/electrode Gauss balance, voltage-offset and energy
identities, native RC response, source-state refresh, and explicit refusals using
synthetic claims rather than measured material qualification. Its source
submodule covers internal charged interfaces, open/shorted polarization,
combined Neumann/sheet/polarization fields, units, high reference voltage,
invalid/duplicate sources and cancellation/retry.

## Multiple dielectric electrodes and native circuit storage

`PlanarElectrostatics::capacitance_matrix` extends the same source-resolved
field owner to a complete ordered set of disjoint electrode traces. The scalar
and multiple-electrode APIs share trace admission; every terminal covers whole
exterior edges, and unassigned exterior edges retain the natural insulated
boundary. Exactly one original checked field solve is performed per terminal,
with explicit terminal, total-solve and retained-field capacity bounds.
Every unit field retains the original mesh, material states and charge reactions.

The full Maxwell matrix is the symmetric electric-energy Gram matrix in farads.
Its coefficients must independently agree with every terminal charge reaction,
reciprocity and the common-voltage nullspace within the caller's coefficient
tolerance. No reference terminal is silently removed. An explicitly grounded
submatrix preserves the remaining original terminal order. A floating electrode
is not eliminated by deleting its row: it still needs its physical charge
constraint in the receiving model.

`response` forms voltage differences before superposition and integrates the
resulting original cell fields for stored energy. `circuit_at_voltages` retains
the source matrix and returns native Kirchhoff capacitor branches with
`C_pair=-C_ij`, together with consistent signed initial branch charges.
Its explicit one-to-one electrode/node map and i-before-j orientation determine
the exact native capacitor charge ordering. Actual terminal charge and network
storage must match the field at the supplied voltages within independent C/J
tolerances. Every terminal needs a capacitive connection. Zero mutual entries
produce no branch; positive off-diagonal entries refuse passive circuit lowering
without clipping, even though the original field matrix can remain useful.

The native circuit owner supplies subsequent KCL/KVL, charge dynamics, source
work and resistor losses. The extracted network is the dielectric store itself,
not an additional parallel energy account. Geometry and dielectric state remain
fixed; parameter work, dielectric loss and leakage remain outside this
reduction. The static field-force connection below is not a motion/work closure.
A balance gate is not a spatial-error certificate.
Focused synthetic tests compare a three-electrode mesh with independent P1
coefficients and direct field charge/energy, verify gauge and terminal ordering,
propagate a regional dielectric-data substitution through actual native coupled
RC dynamics, and refuse a positive-mutual coarse discretization without altering
its field matrix. Full native integration execution is still required in an
environment able to compile the pinned dependency graph.

### Charge-constrained floating electrodes

`electrostatic::network::floating` prepares mixed voltage/total-charge control
on the original Maxwell matrix. Each floating conductor remains equipotential
while its total free charge is held at the supplied value, including nonzero
trapped charge. The retained terminals keep original order; floating inputs use
their explicitly selected order. Topological gauge admission includes common
electrodes joining disconnected dielectric pieces and refuses an entirely
unreferenced component before numerical factorization.

The original `fs-la` Cholesky owner factors the diagonally equilibrated floating
block once. Subsequent queries reuse that factor and the retained-to-floating
voltage map, without new FE solves. Explicit terminal and dense-work allowances
bound preparation; cancellation is checked around the bounded dense factor and
through assembly and field recovery. The Schur complement retains the actual
floating path rather than deleting its matrix row. It is a differential
capacitance at fixed floating charge; nonzero charge also produces an affine
charge offset and stored energy, recovered by the mixed-control solve.

Every accepted solution recovers the original spatial field and independent
weak electrode reactions. Actual-amplitude total-charge and terminal/field
energy gates apply separately from reduction residuals. Voltage differences
remain separate from their common reference to preserve small physical fields.
The full electrode state initializes the existing native capacitor network and
its charge/energy history without discarding the floating nodes. Positive
mutual entries on a coarse field discretization do not prevent this SPD block
solve, although they still refuse pair-capacitor circuit lowering.

Focused `electrostatic::network::floating` regressions cover analytic neutral
and charged three-electrode fields, several simultaneous charge constraints,
capacitance scaling, common-potential and terminal permutation invariance,
native RC decay at conserved floating charge, physical resistor heat, gauge
refusal, independent actual-state gates and cancellation/retry. This is frozen
linear planar electrostatics: no changing geometry/material work, dielectric
dispersion, spatial-error certificate or empirical device validation is added.

## Spatial material DC current and native resistor heat

`conduction::PlanarConduction` binds each original planar triangle to the
canonical `ResolvedConductorState`, preserving its complete source bundle and
explicit independent resistivity/conductivity basis. The weighted Whitney Hodge
and incidence assemble `d0^T M(sigma) d0` in actual conductivity dimensions.
Dielectric and conductor fields share scalar Dirichlet elimination, checked CG,
and complete-electrode trace admission through one internal owner. They retain
their distinct physical laws: conductor reactions are current [A], and
`integral E.J dV` is irreversible power [W], with no dielectric `1/2` factor.

The stationary equation is `-div(sigma grad(V))=0`. Prescribed voltage traces
cover complete exterior edges and anchor each connected component; remaining
edges are insulated. Recovery keeps original cell `E=-grad(V)` and `J=sigma E`
without averaging across interfaces. Weak reactions are current entering the
conductor at prescribed nodes. Net terminal current, free-node current residual,
and the difference between integrated Joule and terminal power must satisfy
independent absolute physical allowances. These discrete balance gates do not
certify spatial approximation error.

`resistance` solves the actual heterogeneous geometry at one volt and checks
its integrated power against the independently summed terminal current before
returning a native `Branch::Resistor`. Its retained source/geometry field is the
same dissipative component as that branch. `PlanarResistance::field_at_voltage`
recovers current and per-original-triangle `cell_power_w` at the native circuit's
actual stage voltage without re-solving. For a midpoint circuit step, the
midpoint branch voltage and `dt*cell_power_w` partition that same native resistor
heat; the spatial report is not an additional loss. Complete material snapshots
and actual extruded volumes remain accessible for an explicit receiving thermal
cell map. This electrical owner neither advances thermal state nor invents a
heat partition between unrelated bodies.

`PlanarResistance::field_at_current` reuses that same prepared linear field for
signed total-current sweeps, including zero. Preparation retains the original
electrode pair and its reference reaction; each recovery independently checks
both terminal currents and authored I*V against the recovered bulk/contact heat.
No new Krylov solve or uniform injection-density assumption is introduced.
`solve_current_with_reference_voltage` uses this same recovery, retaining its
reaction-based voltage scaling. The basis is valid only for its fixed source,
geometry and contact state; material feedback must prepare a replacement.

`conduction::matrix::PlanarConductanceMatrix` extends the same source-bound
operator to arbitrary disjoint complete exterior-edge electrodes. Exactly N-1
native solves prepare a symmetric energy Gram matrix, independently checked
against terminal reactions and the common-voltage nullspace. Every physical
body needs an electrode, and all terminals must connect through bulk, finite
contacts or shared ideal electrodes. Scalar and admitted tensor sources use
the same path. No pair-resistor fit is imposed on the full terminal matrix.

Voltage and balanced total-current drives reuse the original potential basis;
mixed voltage/current partitions factor their equilibrated current block once.
A zero-total-current electrode floats while its local surface currents remain
free to transfer charge. Terminal zero can itself be current-controlled.
Common voltage bias is removed before field/current/heat arithmetic. Native
recovery computes bulk and contact heat together with all simultaneous-drive
cross terms, then independently checks actual terminal reactions and I.V.
Preparation bounds terminal count, native solves, dense work and response
storage; factorization, triangular solves and recovery poll cancellation.
Source/material/contact changes require a new prepared response.

Five native regression groups cover analytic three-terminal bulk/contact
networks, contact heat cross terms under reordered terminals, floating/mixed
controls, disconnected/work refusals and late cancellation/reuse. Three
isolated checks executed the actual matrix module and small native Cholesky
against an analytical field scaffold; complete FEEC integration remains
unexecuted where dependency compilation is blocked by shared resource limits.

Temperature, composition and fixed pressure act through original resolved
source queries. New supported material snapshots require a newly prepared
operator. Both independent electrical source bases use the common scalar
conductor law; no material-name branch or inferred unit-length wire enters.
Frequency-, field-, strain- or rate-bounded sources refuse the stationary scalar
model. The mesh and state remain fixed during each solve. Insulating regions
must be represented by the conductor domain's actual boundary, not artificial
small positive conductivities. Floating unknown-potential electrodes, nonlinear
conduction, displacement current, induction,
Hall effects and geometry/thermal evolution are outside this admission.

`PlanarConduction::new_with_contacts` adds source-backed zero-thickness Ohmic
interfaces to that SAME native stiffness. Each `ConductionContact` binds two
distinct original exterior edges with exactly coincident paired endpoints and
opposite adjacent cell interiors to a complete `ResolvedElectricalInterface`.
The ordered source material states must match those cells, including revision.
The source must support finite positive R'' at static external conditions;
pressure or material changes act through source data and a new preparation.
Explicit zero resistance requires shared mesh vertices rather than a penalty.

The contact law is `J_n = (V_a - V_b)/R''`. Its exact P1 trace matrix is
`A/(6 R'') * [[2,1],[1,2]]`, assembled with equal/opposite side blocks using
actual extruded edge area. Contacts connect component gauges; a voltage trace
can therefore anchor another body through the interface. Internal contact
edges cannot also be external voltage/current electrodes. Duplicate or
nonmatching interfaces refuse before publication.

`ConductionField::contact_response` retains endpoint voltage jumps/current
densities, signed integrated current and exact nonnegative surface Joule power.
`cell_power_w` remains bulk heat; `power_w` is bulk PLUS surface heat, checked
against the same independent terminal work. Native current drive, two-terminal
resistance and circuit-voltage reconstruction all use this operator and total
loss. Surface heat is not silently assigned to either neighboring bulk cell.
Moving-contact, dielectric charge evolution and the bulk-only enthalpy/circuit
adapters explicitly refuse these contacts until their surface work, charge and
heat-partition evolution is owned. Static field/circuit resistance extraction
is supported. `tests/conduction/contacts.rs` checks analytic series resistance,
pressure and bulk substitutions, current-controlled contact connectivity,
nonuniform surface quadrature, total Joule balance and unsupported interfaces.

`new_with_phase_contacts` retains that SAME operator and admits only explicitly
owned liquid-mass-fraction axes in addition to its original fixed-external
conditions. Each declared fraction is an existing dimensional scalar in [0,1];
aliases, missing/type-mismatched axes and unowned response/applicability refuse.
Axis rows follow their original ordered traces through canonical sorting and
remain inspectable on the returned problem. The unchanged `new_with_contacts`
entry does not acquire permission to freeze arbitrary phase or field axes.
The composing thermal owner supplies the exact phase from its original enthalpy
chart and refreshes pinned source snapshots. `fs-couple` now provides that
explicit spatial-contact/enthalpy composition, including two-sided surface heat;
the bulk-only adapters above still refuse contact losses they do not own.
No phase fraction is inferred from temperature, and no phase-interface law,
motion or automatic electrical mixing rule is invented by this constructor.

Existing mesh/assembly and Krylov budgets bound work; cancellation is polled
during material admission, assembly, elimination, field/heat recovery and branch
reconstruction. Refusal never modifies an accepted source or field. G1/G3/G4
regressions in `tests/conduction.rs` compare independent series and parallel
material geometries, source-basis and temperature changes, thickness, voltage
reversal and gauge shifts. Actual repeated native RC steps must lose precisely
the recovered spatial Joule heat. Unsupported source models, trace/assembly/
iteration limits and late cancellation refuse. Synthetic fixtures are not
experimental material validation. The continuum current and dissipation laws
follow MIT *Electromagnetic Fields and Energy*, sections
[7.2](https://web.mit.edu/6.013_book/www/chapter7/7.2.html) and
[11.3](https://web.mit.edu/6.013_book/www/chapter11/11.3.html).

## Electric field, circuit and native mechanical load

`PlanarCapacitanceMatrix::field_at_voltages` reconstructs an actual
`ElectricField` from its retained unit-electrode potentials. This invokes the
same native gradient, weak reaction and physical-energy recovery as a direct
solve, with independent actual-voltage charge and energy tolerances. It does
not perform another Krylov solve. Returned potentials use the first electrode
as their explicit zero reference; a common offset cannot create a field.
Circuit consumers supply their instantaneous endpoint voltages in the original
terminal ordering, preserving the exact dielectric mesh and source snapshots.

`ElectricField::maxwell_wrench` integrates `D tensor E - (D.E) I/2` through the
same contour, outward normals, extrusion thickness, frame/moment transform and
compensated reduction used by magnetic forces. It returns the existing
`fs_ga::Wrench` accepted by native articulated dynamics. Both incident cells of
every contour edge must have the same complete dielectric snapshot, zero
original volume-charge source, and no prescribed electrode vertex on the edge.
The field retains per-cell charge sources even when their global sum vanishes.
Material interfaces, mesh boundaries, folded neighbours and invalid budgets
refuse before any load can be published. Enclosed material contrast and free
charge are supported through their actual computed field.

The assumed homogeneous exterior has no electrostriction or mechanical fluid
pressure response. The result is a net instantaneous load on the enclosed body,
not a distributed traction on the artificial contour. Motion still requires
updated geometry and a source/field/mechanical work closure; no frozen-field
actuator or overlapping free-charge force may be counted as that closure.
Synthetic G1/G3/G4 tests compare against an independent volume-charge force,
check material-only changes through native acceleration, reconstruct the force
after a native multi-electrode RC discharge, and exercise gauge/frame,
source-domain, budget and cancellation behavior. These are discrete numerical
connections, not measured dielectric-device qualification. Force-density and
exterior-medium conventions follow the energy discussion in MIT 6.013,
chapter 11: https://web.mit.edu/6.013_book/www/chapter11/11.9.html.

## Nonlinear reversible induction and spatial Ohmic heat

`planar::eddy::nonlinear::PlanarNonlinearEddy` retains the actual source-bound
P1 field, electrical regions, fixed boundary trace and physical clock. Its
backward-Euler endpoint equation is `F_H(A_next) + C_sigma delta A/dt = load(J_next)`.
The same material H(B) response, FEEC tangent, bounded Newton/Armijo search and
checked CG used by the nonlinear static owner solve this positive variational
problem. The induction term contributes its quadratic objective, residual and
tangent together. No second magnetic solver or inferred permeability is used.

Linear midpoint and nonlinear backward-Euler histories share conductor source
admission, connected-region topology, isolated-current projection, complete
electrical roles, balanced impressed-current loads and algebraic constraints.
They also share the exact consistent-P1 `sigma E_z^2` integration. Every output
`EddyCellHeat` retains its original source, triangle, nodal E/J, volume, power
and heat. A conductivity refresh pins the complete original card, requirements,
coordinate schema, source selection and electrical basis. Nonlinear conductors
also retain the H(B) claim. A temperature change is admitted only when the full
original H(B) range is still supported and its constitutive law is unchanged;
the explicit linear-temperature refresh described below owns changes in mu(T).

For a fixed law, the source-integrated endpoint energies satisfy
`delta psi + Joule heat + numerical dissipation = endpoint impressed work`.
The separately reported numerical term is the nonnegative magnetic Bregman
loss `V*(H_next dot (B_next-B_old) - (psi_next-psi_old))`. It uses the original
positive source slopes and angular convexity to avoid cancellation; the linear
limit is `V*|B_next-B_old|^2/(2*mu)`. This is time-discretization loss and MUST NOT
heat receiving material. Only the existing spatial Joule account is deposited.
Per-cell and total numerical losses, actual source work and the full nonlinear
residual remain visible. Failure of source support, nonlinear budget, current
balance, energy accounting or cancellation returns no candidate replacement.

This is a first-order stationary, isotropic, reversible magnetoquasistatic
section with fixed geometry. It does not add hysteresis, remanence evolution or
conductor motion. A recovered transient field uses the existing explicit
position/frame transfer without re-solving static equilibrium. Synthetic
`nonlinear_eddy` tests cover independent P1 decay, material-data substitution,
Joule/numerical separation, time refinement, conductivity refresh, balanced
excitation through conducting/insulating cells, source refusal and rollback.
The std-only production conduction leaf also has direct analytical projection,
power, floating-current and gauge/reversal tests; those leaf results do not
establish execution of the complete material/mesh/Newton integration.

## Circuit-driven nonlinear material windings

`winding::nonlinear::NonlinearWindingCircuit` couples actual stranded terminal
currents to `PlanarNonlinearEddy` and the original `fs-phs` Kirchhoff descriptor.
The authored signed `J_z/I` profiles reuse the original winding component
balance and exact P1 linkage pairing. A bound original inductor branch becomes
a private native power port. The default current-input causality has input `-I`
and output `V_a-V_b`.
Its placeholder inductance contributes no storage. Unbound inductors,
capacitors, transformers, graph resistors and external ports retain their
original circuit owners and coordinate ordering.

The existing bounded Newton--Krylov and fallible directional-difference owners
solve `V = delta(lambda)/dt + R_winding*I`, with endpoint winding currents and
`lambda = integral(A_z*J_unit dV) + L_end*I`. Every trial restarts both accepted
histories. Declared nonnegative end/leakage self-inductances are optional;
unrepresented end geometry, resistance and mutual coupling are not invented.
Disjoint stranded triangles use the same source/resistance and spatial Joule
owners as linear windings. They cannot also own resolved eddy conductivity.
Complete material/card/claim/basis continuity is retained across conductivity
temperature refresh. The original full nonlinear H(B) law remains frozen;
linear mu(T) may refresh through the explicit material-work path below.

The private port work `-dt*I*V` cancels the field/winding terminal work.
`NonlinearWindingStep.supplied_work_j` contains only original external source
work; its native descriptor receipt still visibly includes private ports.
The combined audit checks field, native circuit and end-inductance storage;
spatial core and stranded Joule heat; explicit graph heat; and separately
reported magnetic/end-inductance backward-Euler loss. Numerical loss is never
physical heating. Terminal voltage residual, source/linkage work pairing and
each heat decomposition must pass physical allowances before publication.
The candidate publishes the circuit and field together, preserving original
histories on source, budget, native constraint, nonlinear or cancellation refusal.

`NonlinearWindingCircuit::new_with_ports` also admits explicit voltage-input
causality for current-driven windings, including ideal current sources. That
private native voltage-source stamp receives the actual
`delta(lambda)/dt + R_winding*I` and returns `-I_native`; the same Newton owner
solves `I_native-I_field = 0`. Current-input and voltage-input windings can share
one nonlinear field, including signed mutual linkage. A declared positive
`residual_scale_ohm` only converts the current residual to the shared numerical
voltage scale. It is not a physical resistor and adds no heat, voltage drop,
stored energy or parasitic regularization. The independent declared
`current_tolerance_a` checks the actual current mismatch; the existing voltage
and combined energy allowances remain mandatory. Step receipts expose both
physical voltage and current residuals, separate from scaled Newton residuals.

Voltage-input initialization requires an explicit initial terminal voltage,
compatible original native fluxes/charges and matching admitted field current.
The original rank-aware native initializer owns that consistency check. Its
zero-time algebraic midpoint receipt determines the current; ambiguous initial
current splits may refuse instead of silently replacing the field history.
Finite voltage-input trials retain the accepted native state and let the
original descriptor step evolve its storage. They do not reinitialize or
project capacitor charge against a trial step-average voltage. Source-bound
graph-resistance refresh retains the same terminal causality and port ordering.
Original external inputs are held over each finite interval, and endpoint field
currents pair with native step-average voltages in this first-order method.
The receipt describes finite-step work; it is not an instantaneous-source
impulse model or an automatic resolution of incompatible ideal-source constraints.

This first-order coupling admits at most 16 windings within the existing
32-node/64-branch/32-port/128-state allocating circuit bounds. The native graph
must admit its explicitly selected private-port trials; singular ideal-source
topologies can refuse this partitioned lane. Geometry, proposed magnetic law and
native graph coefficients stay fixed during each native trial. No hysteresis
or experimental material-qualification claim is introduced.

`nonlinear_winding_circuit::current` supplies independent manufactured oracles
for signed ideal current drive through linear/saturating fields, declared end
storage and continuation; smooth-ramp work and numerical-loss time refinement;
native capacitor/inductor history; mixed causal and reversed mutual windings
with source-bound graph resistance; core eddy shielding
and heat; and initialization/current-budget/cancellation refusals. These cases
exercise the actual field and native circuit owners, not a replacement lumped
inductance model; their material laws do not claim measured qualification.

`NonlinearWindingCircuit::step_with_eulerian_velocity_candidate` resolves
moving bulk conductors inside this same terminal Newton solve. Its actual
nonlinear moving-field candidate changes the winding flux linkage and solved
load current; no fitted speed voltage or independent force current is added.
The public controls, step and candidate retain stationary default field types
and carry the complete moving-field step when this operation is selected.
The combined balance is `delta(storage) + Q + D_BE + W_mechanical = W_external`;
signed mechanical work is positive out of the field, while stranded, graph and
bulk conductor-frame Joule heat remain separate nonnegative physical losses.
Every trial preserves the accepted circuit/field history, and terminal voltage,
source/linkage reciprocity and the combined energy account are checked before
publication. Complete finite P1 velocity is required. Stranded winding cells
and fixed-recoil magnet cells must be stationary; source laws, conductivity,
mesh, material occupancy and profiles stay fixed. The underlying moving-field
divergence, connected-current, source-envelope, magnetic Peclet and nonlinear
residual bounds apply. Moving winding geometry, transported magnetization and
sliding interfaces remain outside this operation.

### Oriented permanent magnets with saturating cores and induction history

`IsotropicMagneticMaterial::Recoil(PlanarRecoilMaterial)` admits the same
source-backed, specimen-oriented recoil law as the original linear owner.
Linear regions, single-valued reversible H(B) cores and permanent magnets share
the original nonlinear P1/FEEC solve. Magnet cells evaluate
`H=(B-Br)/mu_recoil`, shifted storage `psi=|B-Br|^2/(2 mu_recoil)`, and
Legendre coenergy `psi*=mu_recoil |H|^2/2+Br.H`. This coenergy can be negative.
Remanence enters the constitutive weak residual independently of impressed
current, including on natural boundaries; it is not a fitted current source.

Each initial and trial field must satisfy the original parallel and transverse
recoil envelope. Newton backtracking handles both nonlinear H(B) range and
recoil-range refusals. Accepted fields expose per-cell operating points and
remaining source margins beside their declared energy reference. These are
numerical estimates, not certified bounds or an irreversible-demagnetization
model. Construction checks the positive recoil tangent without requiring an
inadmissible zero-field initial state; callers can supply a valid biased seed.

The same fixed law participates in `PlanarNonlinearEddy` and
`NonlinearWindingCircuit`. Physical storage uses the shifted energy, while
the recoil backward-Euler Bregman loss is computed directly as
`|B_next-B_old|^2/(2 mu_recoil)`. Constant Br cancels from this increment.
Winding linkage retains the actual magnetic bias; its time difference supplies
Faraday voltage. Spatial Ohmic loss remains the only magnetic material heat
in this lane. Accepted induction fields can transfer an exterior Maxwell
wrench while retaining enclosed magnets and saturating material.

Conductivity refresh requires the original complete material source and
unchanged recoil quantities. Magnetic-law refresh preserves the exact recoil
law and orientation. Temperature-dependent remanence, permeability or envelope
changes are refused; the caloric nonlinear-recoil path is also refused until
remanence work and entropy data have an owner. Geometry and remanence remain
fixed. No hysteresis, irreversible loss, moving-magnet induction or 3D end-field
claim follows from this extension.

### Source-pinned magnetic temperature refresh on continuing circuits

`PlanarNonlinearEddy::step_with_materials_candidate` and
`NonlinearWindingCircuit::step_with_materials_candidate` reassemble the same
original mesh at proposed source-pinned linear mu(T) states. Original material
cards, complete claims, coordinate schema, fixed context and constitutive kinds
remain pinned. Nonlinear H(B) and oriented recoil entries must remain exactly
unchanged. A refreshed
candidate owns its assembled problem and retains original field, circuit and
linkage histories through clone and continuation.

Each magnetic cell reports signed `parameter_work_j` equal to
`P=psi_next(B_old)-psi_old(B_old)`. The backward-Euler Bregman loss uses the
endpoint law for both field values. Native energy closure is consequently
`delta psi + Joule heat + numerical loss = external work + P`, including native
circuit and declared end storage for windings. Source work and parameter work
remain distinct. A caloric owner must separately resolve field entropy and
reversible heat; neither P nor numerical loss is automatically material heating.
Original fixed-law stepping uses the same implementation with zero P.

`nonlinear_winding_circuit::thermal` supplies independent RL and conducting-P1
checks for data-only mu(T), signed local/global work, physical/numerical heat
separation, owned restart, source correspondence and late cancellation. Native
execution must be reported separately from parser and source review.

`NonlinearWindingCircuit::with_graph_resistors` binds every ordinary graph
resistor to a complete `GraphResistorMaterial`. It reassembles only the original
native resistor coefficients from the common `ResolvedConductor` owner. Original
branch endpoints, transformers, private winding ports, state-coordinate maps,
charges, fluxes, field history and winding linkage stay intact. Exact native
storage equality is admitted before returning an independent candidate at the
same clock. A static scalar Ohmic source and complete unique original-resistor
coverage are required; winding branches cannot be replaced through this API.
The full source bundles and geometry remain with the accepted circuit and its
step report beside the original graph-local heat records. Explicit resistance
reconfiguration stores no parameter energy. Thermal consumers must admit the
accepted sources and refresh them at their actual endpoint through the shared
material/enthalpy owner; this API supplies no temperature or contact law.
Synthetic `nonlinear_winding_circuit` regressions exercise independent nonlinear
RL and coupled-winding/capacitor equations, source-only H(B)/resistivity swaps,
core/strand loss separation, temperature refresh and transactional failure.

## Accepted induction and nonlinear circuit fields drive mechanical loads

`PlanarEddy::maxwell_wrench` and `PlanarNonlinearEddy::maxwell_wrench` recover
the actual accepted B/H from their retained potential and share one source
admission with the original geometric Maxwell-stress integral. They do not
re-equilibrate the field. Every contour edge needs the same complete linear
material snapshot on both sides, zero original impressed current, and no
resolved eddy conductor. Enclosed source-defined nonlinear magnetic material, fixed recoil magnets,
stranded current and eddy currents remain admissible. The history supplies the
source support that a detached recovered static field intentionally lacks.
Mesh work bounds apply before field recovery. Original edge orientation,
extrusion, exterior sampling, compensated reduction, mechanical frame and
moment origin remain owned by the common contour implementation.

`NonlinearWindingCircuit::field()` exposes this same accepted induction state,
so source saturation/resistivity and actual voltage-driven current determine
its native mechanical wrench. The result can feed the existing articulated
body dynamics without a fitted force coefficient or material-name branch.
This is an instantaneous net magnetoquasistatic load on a fixed planar section,
not a local nodal traction map or a moving-field work closure. No energy is
debited by this observation, and an overlapping Lorentz body load must not be
added. Nonlinear/conducting exterior stress, displacement-current electric
stress, end fields, striction and exterior-medium pressure response remain
outside this admission.

G1/G3/G4 tests in `tests/magnetic_force/induction.rs` cover the linear history
limit and independent volume-Lorentz refinement, source-only saturation and
voltage changes through the native winding circuit and body acceleration,
retained decaying-eddy force versus an incorrect static reset, source/contour
refusals and late cancellation. Full native integration execution remains
required; source parsing and metadata compilation are not physical validation.

## Distributed material-conductor current boundaries

`PlanarConduction::solve_with_current_traces` accepts signed uniform inward
current density on actual exterior edges beside prescribed voltage traces.
The original planar metric supplies extruded edge area; exact P1 integration
supplies half the edge current to each endpoint of the existing FEEC solve.
Current traces are not assumed equipotential. Duplicate/interior edges,
nonfinite loads and simultaneous voltage/current traces on the same edge refuse;
every connected component retains a real voltage trace and return-current path.

Free-node current residuals subtract the authored load. Terminal power uses
the prescribed current at current-only vertices and the actual combined weak
reaction at voltage vertices. `nodal_applied_current_a` separately preserves
the applied trace contribution at mixed corners, allowing the voltage-source
reaction to be recovered without double counting. The original cell E.J and
terminal work must agree before acceptance. Zero-load voltage and native
two-terminal reductions retain their existing owner. G1/G3/G4 regressions check
layered material fields, extrusion, reversal, gauge, mixed corners, input order
and transactional cancellation. This is fixed-geometry isotropic DC conduction,
not an ideal equipotential current electrode, contact law or electrical storage.

## Oriented material conductivity in the native spatial solve

`PlanarConduction<C>`, its field and its resistance reduction accept the shared
source-bearing `PlanarConductorMaterial`. The default remains the original
scalar conductor. Oriented principal data and mixed scalar/tensor cells use the
same FEEC assembly, voltage/current boundaries, solve and power acceptance.
The existing weighted owner admits tensor positivity; current is `J=sigma E`
and the original per-cell integral of `E.J` must match terminal work. Tensor
fields expose their full matrix instead of an effective scalar surrogate.
Source-fixed frames must name the same geometry coordinate frame throughout
one problem. Fixed geometry, symmetric DC and static external-condition limits
remain; no Hall response or out-of-plane closure is inferred. Four focused
tests cover independent rotated fields, current traces, heterogeneous resistance,
isotropic degeneration, missing principal sources and inconsistent frames.

The same static field owner also admits an explicitly retained liquid-mass-
fraction coordinate from the shared conductor material type. Its physical
owner must supply that actual coordinate; no arbitrary source axis becomes
fixed merely by appearing in a query. All other source conditions retain the
existing static admission. The field remains a stationary snapshot; latent
enthalpy evolution, fixed-volume phase admission and circuit feedback belong
to the existing `fs-couple` thermal owner. No phase geometry, phase-interface
resistance or automatic electrical mixture law is inferred here.

## Native spatial leaky-dielectric charge and energy history

### Material replacement at retained native charge

`PlanarLeakyDielectric::rebind_material_candidate` substitutes new complete
dielectric/conductor source bundles on the original mesh at the accepted clock.
The native scalar-potential solver recomputes voltage at conserved free-node
charge; those charge coordinates are retained bit-for-bit. Original electrode
voltages remain held. Their charge transfer and signed electrical work are
explicit, separately from signed material-parameter work. The shared native
coefficient-work calculation includes all tensor cross terms and independently
closes `delta_U = W_electrodes + W_material`. Neither work account is Joule heat.

New material cards, electrical basis and material orientations are data; declared
mesh frames must agree. New dielectric laws are validity-wide constant and
subsequent ordinary temperature/phase updates remain pinned to the replacement
sources. Conductivity-only changes retain voltage and charge exactly and alter
subsequent native current/heat. Work bounds, assembly budgets, field/charge
gates and cancellation refuse transactionally. This explicit instantaneous
constitutive replacement does not predict an insertion trajectory, material
formation energy, thermal entropy, mass transport or radiative transients.

G1/G3/G4 tests in `electrostatic/rebind.rs` cover conserved-charge and held-voltage
capacitor limits, heterogeneous redistribution and reversal, rotated dielectric
data, electrical-basis substitution through continued relaxation, gauge
invariance, work refusal and late-cancellation retry. These are synthetic
analytical model checks, not measured material validation.

### Charge relaxation

`electrostatic::relaxation::PlanarLeakyDielectric<C>` retains physical P1
potentials, weak nodal charges, cell dielectric storage and time on the original
planar mesh. It uses the existing dielectric and scalar/oriented Ohmic Hodge
operators and checked potential solve. The midpoint interval satisfies
`(K_epsilon + dt K_sigma/2) V1 = (K_epsilon - dt K_sigma/2) V0`
at free vertices. No lumped fitted capacitance or duplicate cell constitutive
law replaces the native operators.

Original complete exterior voltage traces remain fixed in topology and may
ramp in value. Other boundaries isolate total conduction plus displacement
current. Electrode reactions include both terms and pair with midpoint voltage
for source work. Shared native cell recovery computes physical local E.J heat,
E.D/2 storage and displacement current. The candidate checks free-node current,
net electrode current, weak-charge/cell-storage agreement and
`delta U_field + Q_Joule - W_electrodes` before publication.

The dielectric source must declare validity-wide constant absolute scalar
permittivity. Conductivity may refresh at the actual state, retaining original
material identity, full source claims, basis, frame, phase axis and coordinate
meanings. The original accepted charge is the starting point of every trial;
failure or cancellation publishes no state. Variable permittivity, dispersion,
magnetic induction, breakdown, carrier chemistry and moving interfaces are not
represented.

Held-coefficient midpoint accuracy is second order, but charge sign is not
preserved for arbitrary intervals; the physical relaxation time and spatial
field still require independent resolution. Native regression cases live in
the fs-couple `conduction_phase_enthalpy::dielectric` lane and compare analytic
epsilon/sigma decay, heterogeneous interface charge and voltage-ramp work.
Source/parser review completed for this change. These regression cases have not
been executed in the current environment because pinned sibling dependencies
and configured remote runners are unavailable.


### Directional dielectric storage in native fields and charge relaxation

`PlanarElectrostatics::new_with_laws` assembles scalar or explicitly oriented
per-cell source permittivities through the original weighted FEEC Hodge
operator. Field recovery uses `D = epsilon E` and storage uses the physical
integral of `E.D/2`, including cross terms. The native multi-electrode
Gram/reaction/circuit extraction and superposed response retain that same
tensor storage. Existing pair-capacitor reduction still refuses nonpassive
pair coefficients. No scalar average or conductivity-derived frame is used.
Scalar construction remains available. The scalar accessor
`permittivity_f_m()` now returns `Option<&[f64]>`, with `None` for any
tensor/mixed field; the complete coefficient and source-law views are explicit.

`PlanarLeakyDielectric::new_with_dielectric_laws` retains these laws through
the original continuing charge owner and complete conductivity/temperature/
phase material refresh. Dielectric and conductor material axes may differ;
their geometry frames must identify the same physical mesh. Every selected
principal permittivity must declare validity-wide constancy, and every
rotated matrix entry must remain unchanged during a trajectory. The original
charge, current, storage, source-work and rollback gates remain active.

Maxwell force contours retain their explicitly isotropic exterior convention.
An anisotropic interior may be enclosed by a supported scalar collar.
Declared tensor exterior laws refuse even when their numerical eigenvalues
coincide. This does not supply universal anisotropic mechanical stress,
normal-direction closure, dispersion, induction or moving geometry.

The focused tests specify oblique displacement/storage, specimen/frame
rotation, native capacitance extraction and independent RC decay, and
anisotropic interiors with scalar force collars. Reconstructed runtime
sources match the previously parsed source lengths; final reconstructed
sources/tests received source/API review but could not be parsed, compiled
or executed after the executor became unavailable.


### Retained charge under source-backed epsilon(T)

`PlanarLeakyDielectric::new_with_thermal_dielectric_laws` retains the original
material derivative laws in the SAME native spatial charge owner. Its step uses
`(Kepsilon1+dt Ksigma/2)V1=q0-dt Ksigma V0/2`: q0 is the original accepted
charge, never reconstructed with replacement epsilon. Original conductivity
recovery supplies midpoint E.J; the original dielectric recovery supplies both
endpoint D and Helmholtz storage. Displacement current uses actual D1-D0.

Each original triangle reports the exact two-endpoint coefficient work
`P=-V E0.(epsilon1-epsilon0).E1/2`, closing
`delta_psi+Q_Joule-W_electrode-P=0`. This signed transfer is separate from
Ohmic loss. A thermal owner must also close field entropy/internal energy;
merely applying epsilon(T) and relabeling the storage residual as heat is not
admitted. The ordinary constant-epsilon construction and phase-refresh gates
remain intact. Fixed mesh, voltage traces, independent conductivity/dielectric
frames, original source bundle and selected claims retain their native checks.

Focused coupled regressions cover scalar and oriented source laws, both slope
signs and zero slope, continued charge/heat, electrode work and the isolated-charge
limit. Changed native source passed Rust syntax and bounded independent
accounting/source review; native compilation/test execution was unavailable.

### Saturating and fixed-recoil moving induction

`PlanarNonlinearEddy::step_with_eulerian_velocity_candidate` advances
source-bound reversible H(B) and stationary recoil bias with prescribed P1
conductor velocity. Backward Euler uses the original analytic FEEC tangent plus
the nonsymmetric motion derivative, bounded Newton/FGMRES and residual
backtracking. The linear and nonlinear paths share their matrix-free motion
action/transpose, floating-conductor projection and exact P1 current/heat/work
recovery; Arnoldi actions reuse full-mesh scratch.

One projected endpoint current supplies conductor-frame Joule heat and signed
Lorentz work, closing delta psi + Q + W_mechanical + D_BE = W_impressed. The
original convex Bregman D_BE remains numerical loss, never a thermal credit.
Complete-source slope bounds govern the caller-admitted magnetic Peclet limit.
Recoil cells require zero velocity, whether conducting or insulating. Stationary
magnets may bias adjacent moving metal, but permanent magnetization is not
transported. Geometry, material sources, remanence axes and electrical regions
remain fixed; every failed candidate preserves accepted history.

The scheme is first order in time. The centered P1 Peclet allowance is a
formulation admission, not an error certificate. Moving interfaces, thermal
parameter evolution within the operation, hysteresis, irreversible
demagnetization and terminal-circuit motion coupling require separate owners.
Five physical regressions cover an independently assembled nonlinear P1
residual and heat/work quadrature, an analytic stationary-magnet/moving-metal
transient with time refinement, stationary continuation, saturated uniform-field
charge cancellation and source/domain/work/cancellation refusals.

### Anisotropic moving-conductor charge, heat and force

`PlanarConduction::solve_with_motion` reuses the original source-resolved
scalar/rotated SPD conductivity Hodge and checked charge solve for per-original-cell
planar velocity and an imposed normal magnetic field. The weak motional load,
total-current boundary reactions, `J = Sigma(E + u cross B)`, conductor-frame
Joule heat and `J cross B` force share one current/material/volume owner.
`MovingConductionField` retains the authored motion and field snapshots and
separately reports heat, signed mechanical work, volume `E dot J` and independently
recovered terminal work. Admission checks current continuity, the sum of absolute
local work-pair defects, and `heat + mechanical = terminal power`. Unloaded
exterior boundaries insulate the total current. Velocity or field changes reuse
the prepared conductivity matrix. Original stationary and zero-motional-field
recovery remain intact.

Scope is a prescribed fixed-geometry quasistatic snapshot; no self-consistent
induction, Hall law, moving mesh or fluid/solid time integrator is inferred.
Six analytical tests cover open-circuit charge separation without spurious drag,
shorted braking, loaded generation and motor signs, rotated anisotropic currents,
imposed-current EMF, heterogeneous counter-moving layers and stationary/refusal
behavior. Synthetic sources verify these equations, not an empirical material fit.

### Three-dimensional current-free magnetic fields

`spatial::SpatialMagnetostatics` solves `H = -grad(phi)`,
`B = mu H + Br` and weak `div(B) = 0` on actual tetrahedra. Each original
cell retains its resolved positive permeability or source-backed oriented recoil
law. The original Whitney 1-form material Hodge, incidence/stiffness assembly
and checked CG own the numerical operator. Prepared geometry and stiffness are
reused across boundary solves. `try_element_geometry` supplies a checked
companion to the original tetrahedral geometry owner and refuses malformed,
singular or unrepresentable cells before field publication.

Explicit exterior faces accept continuous nodal magnetic potential [A] or
outward normal B [T]; unspecified faces impose zero normal flux. Each disconnected
pure-flux component must satisfy its own integrated flux compatibility and
receives one deterministic gauge. The independent full weak residual includes
the artificial gauge node, so the gauge cannot supply missing flux. Relative
potentials preserve field and energy under representable common gauge shifts.

Recovery returns actual cell H/B, original recoil operating-point admission,
weak boundary fluxes and interior normal-flux jumps. Recovered cell energy is
checked against the assembled Hodge energy. Fixed remanence uses the existing
shifted reference `mu |H|^2 / 2`; this is not an independent force or variable-
magnetization energy law. Failed domain, solve, work or cancellation checks
return no partial field.

This is a single-valued scalar-potential model on current-free domains with
fixed linear/recoil sources and explicit finite boundaries. Linked free-current
circulation, general 3D current-driven curl-curl fields, induction, saturation,
hysteresis, motion and infinite-domain accuracy require other owners. The mesh
must be a conforming nonoverlapping tetrahedral complex; local adjacency checks
do not certify absence of remote cell intersections. Raw P1 B is cellwise
constant and its interface jumps are reported, not claimed to be an exactly
equilibrated H(div) flux.

Seven analytical tests cover vector fields with potential/flux boundaries,
large gauge offsets, oriented recoil/air demagnetizing layers and energy,
harmonic-field mesh convergence, disconnected compatibility and malformed
geometry/boundaries, original recoil limits, true solver-budget exhaustion and
late cancellation. These are synthetic physical oracles, not measured device
validation or a general 3D electromagnetic capability claim.

### Saturating three-dimensional current-free fields

`spatial::nonlinear` uses the original current-free H = -grad(phi) weak flux
problem and source-resolved linear, reversible H(B), or oriented fixed-recoil
materials on the same tetrahedral mesh. It inverts the original monotone H(B)
segments without fitting a second law. The exact radial differential
permeability supplies a tensor Whitney tangent; bounded Newton/checked-CG
steps use decreasing integrated coenergy with source-range-preserving
backtracking. Geometry, scalar boundary admission, component gauges and the
reduced checked-CG owner are shared with the original linear spatial solver.

Results retain distinct source-integrated energy and coenergy, a Legendre
identity audit, full weak flux residuals (including gauge vertices), raw
interface flux jumps and recoil operating points. Physical boundary work uses
unknown reactions only on physically prescribed potential traces; elsewhere
it uses the authored normal flux. A loose residual allowance cannot turn an
interior defect into boundary supply. Invalid source ranges, unbalanced
Neumann flux, exhaustion and cancellation return no accepted field.

The focused spatial_nonlinear target covers reduction to the original linear
solver, rotated saturating/air layers, recoil/core/air flux and exact energy,
a refined nonlinear radial manufactured field, differential material response,
and source, energy, budget and cancellation refusals. This is a fixed-state
single-valued scalar-potential model; currents, hysteresis, eddy currents,
moving work and thermal evolution require their respective owners. P1 flux
jumps and finite exterior truncation remain discretization errors; these tests
do not establish an experimental material calibration or a performance claim.

### Three-dimensional Maxwell force and torque

`SpatialMagneticField` and `SpatialNonlinearField` expose `maxwell_wrench`
through one private 3D surface integrator. Original tetrahedron selection
defines a closed outward surface. Exact triangle area/centroid quadrature
integrates outer-cell Maxwell stress and torque. Every contour neighbor must
share one complete source-identical induced-linear exterior snapshot;
nonlinear/recoil interfaces and mesh boundaries refuse. Saturated and remanent
targets may be enclosed inside a linear collar. Canonical Se3 placement and an
explicit receiving-frame torque origin return the native Wrench with a
canonical compensated reduction, work budgets and cancellation.
The final origin shift subtracts common frame/origin translation before adding
the rotated nearby mesh reference, retaining small physical moment arms.

The five focused force tests cover independent nonzero cuboid refinement on
two surfaces, saturated-core symmetry and interface refusal, frame/origin
behavior, and source/region admission. The cuboid result has first-order raw
field/surface convergence. Test-only extrapolation checks its independent
continuum limit; the public API returns the raw discrete stress result, with
no extrapolation or inferred error bound. Field and surface refinement are
required for useful force accuracy. This is fixed-geometry static loading;
moving-geometry work exchange, material striction, equilibrated stress and a
virtual-work cross-check require separate physical work.

### Current-carrying three-dimensional vector-potential magnetostatics

`spatial::vector_potential` prepares the actual lowest-order 3D Whitney curl-curl operator d1^T M_(nu,2) d1 for source-bound isotropic linear and oriented fixed-recoil media. It admits exact piecewise-constant impressed-current loads, exterior tangential-potential edge traces, and natural n cross H=0 elsewhere. Current continuity and the entire relative curl-kernel compatibility are checked, including harmonic modes on handled domains, without silently changing the source. The existing checked CG solves the explicitly compatible invariant ranges of semidefinite operators; neither the full curl-curl operator nor the auxiliary face operator is claimed globally SPD. The bounded original cohomology-metric projection removes all curl-free gauge components and must preserve every magnetic face flux. Accepted fields retain truthful free-edge equations, divergence-free face flux, reconstructed B/H, original recoil-envelope admission, shifted storage, and independent current/trace constitutive pairing. The shared Maxwell integrator additionally requires one source-identical induced-linear current-free exterior collar. This is static lowest-order magnetostatics; higher-order spaces, Coulomb saddle/HX solvers, 3D eddy or nonlinear vector evolution, and TEAM qualification remain open.

The six focused spatial_vector_potential tests cover an independent
current-driven quadratic-potential manufactured solution under refinement,
exact uniform vector fields with gradient-shifted traces, a handled-domain
harmonic source that passes local continuity but must be refused globally,
recoil energy and limits, actual budgets/cancellation, and the current-free
Maxwell collar. The stored static source/trace pairings are not time-integrated
work or Ohmic heat.

### Current-driven reversible saturation in three dimensions

`spatial::vector_potential::nonlinear` extends the original Whitney edge
magnetostatics to source-bound reversible H(B), linear and oriented fixed-recoil
cells. Exact original source energy supplies H and its full symmetric radial
reluctivity tangent; only the changing material Hodge is rebuilt. Geometry,
current loading, boundary admission and complete relative-kernel compatibility
remain shared with the linear vector owner. Newton minimizes integrated source
energy minus impressed-current pairing at fixed tangential A, with bounded
source-preserving backtracking and the existing checked CG/coexact owners.

Newton's assembled dual residual is projected into the admitted coexact tangent
space before its relative linear solve. The measured compatibility defect must
stay inside the original physical current allowance. Impressed current is never
changed: every accepted field must satisfy the original full free-edge equations,
source Legendre identity, div-curl/reconstructed face flux and independent
current/physical-trace energy pairing. Energy and coenergy remain distinct.
Maxwell loads use the existing current-free, source-identical linear collar.

The six focused nonlinear-vector tests cover reduction to the original linear
current solver, gradient-shifted traces, first-order saturated current-driven
manufactured-field convergence, rotated source energy derivatives, recoil/core/
air tangential-field and dual-energy agreement, actual source/iteration/work/
cancellation refusals, and force admission. The six original vector regressions
also pass after the shared extraction. This is fixed-state static reversible
magnetism; hysteresis, nonlinear induction, thermomagnetic evolution and
validated engineering force bounds are separate extensions.


## Prepared distributed Maxwell force and torque

`spatial::force::eggshell::PreparedSpatialEggshell` caches an original-mesh
P1 virtual-displacement functional. Nodal weights are one over a complete
target plateau and zero on every exterior mesh vertex. All changing-weight
cells and intermediate constant plateaus require one identical original
induced-linear, current-free exterior material when the functional is applied.
Source/current-carrying targets cannot hide inside fractional plateaus.
Preparation binds the actual immutable mesh/positions, not merely matching
dimensions or topology, and is reusable across accepted excitations.

Translation integrates `-T grad(w)` using original cell B/H and volumes.
Rotation uses the original P1 interpolant
`sum_i N_i w_i (omega cross (x_i - origin))`, including its complete gradient
moment. It is exactly rigid over the target plateau and vanishes outside.
Both functionals integrate the discrete constant cell fields exactly without
smoothing interfaces, inferring a stress inside magnetic matter, or requiring
another field solve. Deterministic compensated reduction returns the canonical
`fs_ga::Wrench`, transformed through the caller's SE(3) frame and moment origin.
As in the surface integral, the final origin shift preserves the nearby mesh
reference under common frame/origin translation.

Scalar linear/recoil, scalar reversible nonlinear, current-driven linear/recoil
vector, and current-driven nonlinear vector owners share this same sealed
geometry/stress implementation. Original source identity, current-free exterior,
field provenance, geometry budgets and cancellation gate every application.
This is instantaneous mechanical load transfer. It does not supply a moving
domain energy exchange, material striction, a source/control shape adjoint,
equilibrated stress or a verified physical error bound. Collar-profile and
mesh refinement remain necessary.

The six added cases in `spatial_maxwell_force/eggshell` check raw cuboid
force/torque against analytic targets and two collar profiles, independent
re-solved coenergy translation finite differences, reuse/frame/origin/geometry
identity, complete closed weights and source/fractional-plateau refusal,
a saturated enclosed core, and vector/current-collar admission. The original
surface-force regressions remain in the same target.


## Moving-conductor 3D A–phi induction on the native Whitney owner

`spatial::vector_potential::induction` continues an admitted
`SpatialVectorField` with fixed tangential-A traces and original linear or
fixed-recoil magnetic materials. Prepared sparse curl-curl and conductivity
blocks couple the magnetic increment to electric scalar potential:
`J = sigma * (-A_dot - grad(phi) + v cross B_mid)`. The same P1 stage current
enters Ampere's law, charge conservation, Ohmic heat and Lorentz mechanical work.
Exact tetrahedral P1-product integration supplies the quadratic transfers.

A deterministic relative tree gauge removes gradient freedom while retaining
harmonic circulation. Dynamic A is not independently projected into a static
coexact representative. Face-connected conducting components own separate
floating potentials, including distinct potential DOFs at incidental point or
edge contacts. Conductor boundaries/interfaces are naturally insulating.
The sparse nonsymmetric motion operator uses the original bounded, cancellable
FGMRES owner with a recomputed true residual of the actually retained increment. Rounded
zero-physics A increments are admitted only when this residual and all physical
gates still pass; solver tolerances are not relaxed. Full physical Ampere equations,
including gauge-eliminated rows, every component's charge equations, endpoint
insulating equations, discrete div-curl, recoil operating bounds and complete
field/Ohmic/Lorentz/impressed-work balance independently gate acceptance.
Unknown/assembly/Krylov budgets and a caller-admitted cell magnetic Peclet
allowance bound the prepared operation.

Conductivity must retain the identical original magnetic material bundle and
fixed external conditions. Impressed cell current is piecewise linear in time;
nodal velocity is held over a step. Coefficients and Eulerian occupancy stay
fixed. This is a lowest-order, linear/fixed-recoil MQS rung: moving interfaces,
ALE, displacement current, saturation, hysteresis and thermal feedback are not
inferred. The terminal extension below separately admits physical electrodes.
Heat and mechanical transfers require their own receivers; returned candidates
do not silently advance either receiver.

`spatial_induction` checks real floating charge cancellation of uniform motion,
gradient-gauge invariance, an independently integrated charge-eliminated cube
mode and second-order time refinement, analytic three-dimensional motional
feedback/Ohmic/Lorentz work, separated floating conductors, fixed-recoil energy
and source limits, and bounded/cancelled deterministic rollback. Synthetic
materials qualify equations and discretization, not empirical devices.

### Prescribed magnetic boundary excitation

`SpatialInduction::step_with_boundary_candidate` advances a piecewise-linear
tangential-A trace on exactly the original exterior edge support. The known
`-delta_A/dt` enters the same Whitney magnetic and floating-charge equations;
the midpoint boundary curl also enters `v cross B`. The original gauge,
conducting components and sparse operators are retained. Accepted endpoints
become the held trace of subsequent ordinary steps, preserving eddy-current
memory through a ramp followed by a hold. Support changes and conflicting or
nonfinite face traces refuse transactionally.

`magnetic_boundary_work_j` is the signed original magnetic reaction paired
with the boundary increment. It is a separate external supply in
`delta W + Q_Ohmic + W_mechanical = W_impressed + W_boundary + W_terminal + W_parameter`.
It is neither impressed volume-current work nor a thermal credit. The material
trial and `step_with_materials_and_boundary_candidate` retain the same split
under source-pinned conductivity/permeability refresh. Physical electrodes,
nonlinear/hysteretic evolution and moving geometry do not acquire a new
boundary-drive API by implication.

G0/G1 fixtures check uniform-field storage and signed boundary supply, an
independently integrated one-tetrahedron P1 Faraday/charge/Ohmic oracle,
time-dependent gradient cancellation, ramp-to-hold diffusion, fixed-trace
equivalence and rejected-candidate rollback. These are numerical-model checks,
not an external-coil model or an empirical induction-heating qualification.

### Physical electrodes on three-dimensional inducing conductors

`induction::terminals::SpatialTerminalInduction` attaches disjoint equipotential
electrodes to an accepted native A–phi history. Each electrode is a union of
complete exterior faces of one face-connected conductor, with fixed tangential
A. All electrodes on a conductor share a connected fixed-A reference path;
disconnected trace pieces cannot silently define a gauge-dependent voltage.
The original relative magnetic gauge, cell materials, anisotropic/scalar
conductivity, physical state and clock are retained. Conductors without terminals
keep their floating charge gauge. Prescribed midpoint potentials enter the same
Whitney/charge form; the scalar field and induction current are solved together.

The native nodal weak reactions supply signed total electrode currents rather
than a presumed uniform current profile. Interior/floating charge residuals and
the sum of terminal currents on each conductor independently gate acceptance.
Voltages are shifted to each conductor's first electrode before assembly and
power accounting. `terminal_work_j` is the incoming electrical supply and enters
the existing magnetic-storage/Ohmic/Lorentz/impressed-work balance separately.
Candidates reuse the original consistent nodal Lorentz-load transfer, including
its exact P1 virtual-work pairing. No terminal impulse changes accepted A.

The voltage-driven rung holds geometry, occupancy and material sources fixed;
it does not supply contact impedance, external circuit/lead fields, ALE,
nonlinear or hysteretic magnetic evolution, or displacement current. Existing
floating nonlinear and thermal histories retain their separate contracts.
`spatial_induction::terminals` checks an independently integrated conducting
cube's RL transient, second-order time refinement, the DC resistance limit,
free-charge/electrode conservation, exact voltage common-mode invariance,
motional open-voltage and short-circuit response, mechanical transfer, and
budget/geometry/cancellation refusal without partial state publication.

### Total-current electrodes on three-dimensional inducing conductors

`SpatialTerminalInduction::into_current_driven` prepares
`terminals::current::SpatialCurrentInduction` without changing accepted A,
material sources or time. All charge nodes on a physical electrode share one
scalar unknown; its equation prescribes the total incoming reaction current.
Interior charge nodes remain independent. Each face-connected conductor keeps
its own reference electrode, and prescribed current must balance separately on
that conductor. The grouped mass is the original `P^T M P`: distinct local basis
cross terms are counted twice even when they share a global diagonal entry.

Preparation caches the electrode space and sparse operator. Each stage uses one
original coupled A–phi solve; no pointwise current profile, DC resistance or
sampled dense electrode response substitutes for induction. The solved scalar
potential supplies the actual electrode voltage. All electrode currents,
including the reference reaction, are independently checked against their
inputs. The same native current supplies Ohmic heat, consistent mechanical
loads and signed terminal work in the original midpoint energy balance.
Conversion back to voltage control retains the identical accepted field/clock.

G0/G1 regressions exercise an independently integrated current-driven cube
transient and DC limit, voltage/current duality with interior charge nodes,
motional open-circuit voltage and generator work, refusal of globally balanced
but separately incompatible conductor currents, and budget/cancellation
rollback. The parent fixed-geometry/material MQS scope remains in force; this
adds no contact impedance, displacement current or moving material interface.

### Native external circuits on physical 3D induction electrodes

`spatial::vector_potential::induction::terminals::circuit::SpatialTerminalCircuit`
joins the original physical electrode layout to an authored `fs_phs::circuit`
RLC graph, reciprocal mutual inductance, ideal transformers and original source
ports. Native fluxes and charges retain their original branch ordering and
storage; `native_seed` does not project them. Every finite-stage field and circuit
trial starts from the same accepted histories. Only the complete jointly checked
candidate advances either history.

For fixed linear/recoil magnetics, conductivity, geometry and prescribed
velocity, the actual midpoint field step supplies the affine electrode relation
`I = I0 + Y V`. This retains diffusion, magnetic history and motional generation;
`Y` is not a static resistor model and its power is never labeled heat. One
positive numerical reference resistor and compensating native current port per
unshorted coordinate give `j = V/Rref - I`. The native circuit response
`V = V0 + H j` closes through a small scaled LU solve of
`[Id - H (Gref - Y)] V = V0 - H I0`. Final native and full-field replays must
recover the same terminal voltages and currents. An electrode wired to the same
circuit node as its conductor reference has exactly zero relative voltage and
creates no artificial self-loop source; a wholly shorted layout needs one field
step and no response factorization.

Explicit voltage/current probe scales, a field-trial budget, native state and
interface-memory bounds apply before response construction. Probe fields retain
the original constitutive operating-domain checks. At most 16 unshorted terminal
coordinates, 32 original source indices, 48 total native ports, 32 circuit nodes,
128 total branches and 128 native states enter this dense interface. A cached
response requires identical duration, velocity, fixed configuration, field solve
settings and probe scales. History offsets and final physical replays are always
fresh; the current call's voltage, current, Schur and energy allowances still
gate acceptance. A cold nontrivial interval uses `coordinates + 2` field trials;
a response reuse uses two.

Positive electrode current enters the conducting body. Its terminal supply is
`dt * sum((Vt - Vref) * It)` and cancels native interface work after removing the
numerical reference loss. Independent original branch currents reconstruct KCL,
resistor heat, capacitor charge work, inductor flux work including mutual terms,
ideal-transformer work and original external-source supply. Combined field and
circuit storage, conductor-frame Ohmic heat, physical graph-resistor heat,
Lorentz work, impressed-current work, signed magnetic parameter work and external
supply must close. The native
diagnostic record retains the canceled reference contribution; the physical heat
fields exclude it. Consistent Lorentz loads use the same accepted field current.

`SpatialTerminalCircuit::step_with_materials_candidate` admits the complete
original source-pinned induced-linear material refresh used by voltage/current
electrodes. Endpoint operators are prepared once and shared by the baseline,
each affine probe and final replay; each starts with the original accepted A,
old magnetic source and native circuit flux/charge. A material refresh invalidates
the previous field response. The new response may be reused after acceptance
under its exact endpoint operators, duration, velocity and numerical controls;
history offsets and physical final replays remain fresh. `field_trials` counts
actual full-field solves, including the single solve of a wholly shorted layout.
Independent changing-material RC/RLC regressions retain before-state magnetic
effort, capacitor charge and inductor flux, check original source work and
reference-impedance invariance, and verify cache/budget/support/cancellation
behavior. This entry reports parameter work without choosing temperature or
depositing reversible magnetic heat; the finite-caloric consumer owns that law.

This owner adds no switching, source-dependent circuit storage,
geometry refresh, contact impedance, thermal deposition, nonlinear field step,
moving interface or displacement-current law. The native circuit retains its
own topology and initial-constraint refusals. Focused electrode-circuit tests
compare independent reduced FE/RC/RL midpoint equations, numerical-reference
invariance, motional generation, ideal source and shorted-node topologies,
recoil-domain refusal and preservation of both histories on rejection.

### Source-backed resistors in linear 3D electrode circuits

`SpatialTerminalCircuit::with_graph_resistors` binds every physical authored
resistor branch to its complete `GraphResistorMaterial` source and uniform path
geometry. The shared native resistance reassembly retains graph topology, mutual
inductance, terminal/private-port ordering and exact accepted flux/charge state;
it checks unchanged storage before returning. Field A, magnetic source and clock
also remain unchanged. This explicit source reconfiguration can attach the first
bindings; a finite-heat consumer separately pins accepted source selections and
admits their temperature successors. Resistance changes create no stored-energy
parameter work and never reseed native algebraic efforts.

`graph_resistor_materials` exposes the retained original-branch ordering.
Identical sorted bindings retain the response cache; changed bindings invalidate
it before another affine circuit closure. Every accepted candidate, including a
wholly shorted electrode layout and a simultaneous field-material refresh,
retains the actual source bindings alongside per-branch physical heat receipts.
Shared native accounting excludes artificial reference resistors from those
receipts. Focused regressions check independent R(T)/field/RLC midpoint dynamics,
old flux/charge preservation, original source work, heat partition, cache reuse
and invalidation, and complete/unique branch-map and cancellation refusal.

## Source-backed orthotropic conductivity in native 3D induction

`fs_material::conductor::spatial` supplies a sealed shared conductor boundary
for the original scalar source, three explicit principal sources and heterogeneous
scalar/orthotropic cell maps. `ResolvedOrthotropicConductor3` retains the complete
resolved material bundle, exact conductivity/resistivity basis and property keys,
and finite orthonormal principal directions in physical mesh coordinates.
Resistivity is inverted in the principal frame; the original tensor rotation
owner forms the symmetric positive coefficient. Missing directions, incompatible
units, nonstatic sources and unrepresentable positive tensors are refused.

Both native induction histories, their candidates, temperature refresh and
Maxwell/eggshell load adapters retain the same conductor type. Existing scalar
`new` constructors preserve all-insulator type inference; `new_with_conductors`
admits the shared source owner. One coefficient supplies Whitney/charge mass,
both motional operators, J=sigma E, exact P1 E.J heat, charge conservation and
J-cross-B mechanical work. The actual matrix infinity norm bounds the magnetic
Peclet gate. Scalar arithmetic keeps its original multiplication order.

Temperature refresh keeps every principal key, source basis and exact fixed
frame; it cannot replace orientation or fall back to an unrelated scalar source.
This does not add Hall response, frame evolution, ALE, electrical memory or a
phase-mixing law. Magnetic constitutive and finite-caloric ownership remain in
their existing subsystems.

The focused `spatial_induction::orthotropic` cases specify independently
integrated cube-mode decay and radial-motion heat/Lorentz work for aligned and
oppositely rotated tensors, explicit principal resistivity inversion, source and
frame refusal, scalar compatibility, heterogeneous cells and backward-Euler
numerical-loss separation. The analytical charge-eliminated mass is
tr(sigma)/40-(sigma_xy+sigma_xz+sigma_yz)/120, so a rotation changes the native
trajectory and heat without a change to the physics code.

## Reversible saturation in three-dimensional A–phi histories

`spatial::vector_potential::induction::nonlinear::SpatialNonlinearInduction`
continues an admitted original nonlinear vector field with coupled conductivity,
floating charge and prescribed nodal motion. Backward Euler uses H(B_next) and
J=sigma*(-delta_A/dt-grad(phi)+v cross B_next). The original source's consistent
differential reluctivity enters the original FEEC tensor curl-curl operator and
bounded Newton/FGMRES. The shared `induction::coupling` implementation owns the
relative tree gauge, conducting components, sparse conductivity/motion blocks,
source compatibility and exact P1 current/Ohmic/Lorentz recovery for both
transient discretizations. No independent copy of those equations is introduced.

Original linear, reversible isotropic saturation and fixed-recoil cells may
coexist. Fixed tangential traces and actual prescribed-edge masks are retained.
Every nonlinear trial starts from the same accepted A and source history.
The retained full Ampere and charge equations, div-curl, source/recoil envelopes
and total work independently gate acceptance. Iteration, assembly and Krylov
work bounds and cancellation return no partial state.

The source-exact 3D Bregman loss D=H_next.delta_B-delta_psi is computed from
positive original segment slopes and angular convexity. Reports separately
retain signed magnetic storage change, physical Joule heat, Lorentz work,
impressed work and nonnegative numerical dissipation. D is never thermal heat,
hysteresis or a fitted core-loss term. This is first-order fixed-geometry,
fixed-magnetic-law MQS evolution; vector hysteresis, ALE, terminals, resolved
displacement current and mechanical state advancement remain separate work.

Seven new tests check an independently charge-eliminated BE mode, analytic
saturated decay and first-order refinement, exact multi-knot work/heat/loss,
three-dimensional motion/charge/Lorentz feedback, rotation/gauge covariance,
fixed-recoil energy/envelopes and transactional source/work/cancellation
refusal. Six original linear midpoint induction regressions remain passing.

## Source-pinned temperature refresh of three-dimensional induction

`SpatialInduction::step_with_materials_candidate` evolves the retained linear
A–phi history with endpoint source-backed conductivity and permeability.
Original geometry, Whitney bases and incidence allocations are shared; only
material-weighted operators are rebuilt. The complete original magnetic and
electrical source selections, conductivity basis, temperature coordinate,
fixed traces, conducting support and accepted history remain authoritative.
Per-cell receipts separate signed permeability parameter work from native
Ohmic heat and prescribed Lorentz work. Magnetic free-energy change plus
Ohmic heat plus mechanical work equals impressed work plus parameter work.
The endpoint field and source bundle are retained together; a static re-solve
cannot discard induction memory. This lower owner does not choose temperature
or deposit reversible heat: the finite-caloric consumer must account for
magnetic entropy/internal energy and positive combined heat capacity.

### Source-pinned temperature refresh of physical 3D electrodes

`SpatialTerminalInduction::step_voltage_with_materials_candidate` and
`SpatialCurrentInduction::step_current_with_materials_candidate` apply the same
complete original source admission to physical equipotential electrodes.
Only material-weighted magnetic and conductivity operators are rebuilt; the
original trace, relative gauge, electrode layout and grouped scalar rows remain
intact. Every solve pairs the accepted old magnetic effort with the new endpoint
effort. Updating the material cannot replace or re-equilibrate the old A history.
The total-current path still requires one sparse field solve and recovers the
electrode voltage, spatial current distribution and terminal work together.
Original charge, Ampere and work checks include the signed permeability
parameter transfer separately from physical Joule heat. These entries accept
original source-pinned induced-linear cells, not temperature-varying recoil or
nonlinear magnetics; they neither choose temperature nor mutate a heat receiver.
Focused six-tetrahedron transient tests change both permeability and conductivity,
compare independent endpoint-coefficient equations and voltage/current duality,
and exercise source/support/cancellation refusal without accepted-state mutation.

`SpatialNonlinearInduction::step_with_conductors_candidate` provides the same
source-pinned conductivity refresh for retained saturating, linear and fixed-
recoil magnetic histories. The complete original conductor bundle, electrical
basis and selected claims stay fixed; only the source-owned absolute temperature
may change. Original magnetic admission checks the whole H(B) or recoil range,
so equal endpoint coefficients cannot conceal changed magnetic storage. These
magnetic laws must remain temperature independent. Original relative gauges,
conducting components, geometry and A history are retained while the shared
conductivity operator is rebuilt. Local Joule heat, signed Lorentz work and
positive backward-Euler numerical loss remain separate native transfers.
The finite-enthalpy consumer in `fs-couple` determines temperature and deposits
only physical heat; this lower owner never mutates a thermal receiver.

## Maxwell loads from retained three-dimensional induction

`SpatialInduction` and `SpatialNonlinearInduction` now expose
`maxwell_wrench` and `eggshell_wrench` through the same private load bridge.
Original constitutive B/H is recovered directly from accepted A, with no
static re-solve that could discard eddy-current memory after source removal.
The existing closed-surface integration and prepared P1 virtual-displacement
functional retain canonical SE(3), moment-origin, deterministic reduction and
geometry-identity behavior. No duplicate magnetic stress or force model is added.

Every surface neighbor or eggshell transition/fractional plateau must be one
complete source-identical induced-linear exterior, without resolved conductivity
or impressed current. A momentarily idle conductor is not accepted as insulating
air. Conducting, recoil or saturating targets may lie on the enclosed weight-one
plateau. Mesh budgets bound complete surface-field recovery; cancellation and
refusal mutate no accepted electromagnetic state. Returned instantaneous
wrenches do not advance mechanics or debit mechanical work, and overlapping
Lorentz body loads must not be added again.

Two substantive tests use an independently integrated compact P1 current loop.
For manufactured m_z=0.125 and imposed B_x=0.4, raw initial torque_y is
0.0499977666885 N m versus the independent 0.05 N m reference. Retained torque
after impressed-current removal is 0.0383014007176 N m and decreases to
0.0116541343558 N m over four subsequent steps; the fresh static source-off solve
instead has zero torque. Both transient owners and both force functionals,
analytically equivalent linear-limit endpoints, frames/origins, late cancellation
and strict conductor/current/source collar refusals are exercised.
This finite-mesh manufactured result is not a universal force-error bound or
experimental device validation; field/collar refinement remains required.

## Nonlinear scalar material response in the original spatial conduction solve

`conduction::response::PlanarConductionResponse` closes original pinned
conductivity/resistivity sources at the actual recovered |E| and |J| of each
triangle. Its bounded iterations use the SAME FEEC assembly, voltage/current
traces, field recovery and charge/terminal-work admission as
`PlanarConduction`. The ordinary constructor keeps its static-source refusal.
The frozen linearization is private outside the explicit coupled-solve callback;
nonlinear response is never obtained by scaling a converged unit-voltage result.

Accepted candidates retain both the applied material coefficients that produced
the native fields and their complete source evaluations at those fields.
Independent conductivity and current-density defects, plus the sum of absolute
cell power defects, must pass alongside native algebraic/current/work gates.
Numerical iteration defects are not Joule heat. Rejected sources, budgets,
nonconvergence and cancellation leave the accepted owner unchanged.

`with_relaxation(weight)` selects optional fixed under-relaxation of the owned
E/J operating coordinates. Each damped trial requeries the complete original
source, including its actual support, pins and other material roles; temperature,
phase and all other context coordinates are unchanged. No coefficient averaging
or relaxed physical acceptance is allowed. The default weight one retains the
undamped arithmetic and query count. Smaller weights can close strongly
alternating positive-resistance laws; they do not establish global convergence.
The existing iteration budget bounds trials and `source_evaluations` counts
both actual-field and damped-coordinate queries. Accepted candidates retain the
selected weight. `into_frozen_parts` transfers the physical operator together
with a complete `ConductionResponseSources` bundle, retaining numerical
relaxation and any original interface response bindings on rejoin.
Independent regressions demonstrate unique passive resistance/conductivity
roots that fail the undamped budget but converge with weight one-half, and
preserve signed drive, local heat, source-domain refusal and cancellation/retry.

Both planar and tetrahedral owners also expose `with_adaptive_relaxation`.
`AdaptiveConductionRelaxation` bounds the initial/minimum/maximum weights in
`(0,1]`; the shared `fs-ports::BoundedAitkenRelaxation` kernel selects subsequent
weights from the original bulk and contact operating-coordinate defects. Each
coordinate has its own scale, shared by both residuals in the secant calculation.
Trials still requery the complete original sources; no resistance/conductivity
blending, support clipping or numerical heat is introduced. Adaptive history is
local to one solve. Accepted candidates and planar frozen handoffs retain the
configured policy, while `with_relaxation` explicitly restores fixed relaxation.
Added analytical regressions specify steep positive table laws whose fixed
half-relaxation is unstable, both drive polarities, 2D/3D current and heat
balances, original contact-source refresh, budget/domain refusal and replay.
These new native regressions have not been executed in the current environment.

This fixed-geometry scalar memoryless rung accepts prescribed native voltage and
current-density traces. A nonlinear equipotential total-current terminal needs
its own constraint; tensor, Hall, moving-conductor, field-storage and hysteretic
paths retain their own owners. Positive source transport does not prove Picard
convergence, constitutive uniqueness or dynamic stability.

Independent regressions exercise two-region nonlinear series continuity,
current-driven resistivity, simultaneous E/J dependence, data-only source swaps,
local E.J power, rotation/reversal/gauge invariance, constant-source identity,
separate power admission, continuation and refusal rollback.


## Total-current electrodes with original source response

`PlanarConduction::solve_current` imposes signed total current on complete equipotential electrodes using the original frozen FEEC operator. Actual-amplitude recovery checks both terminal currents and authored I*V independently; source-dependent callers repeat this constraint inside the existing E/J closure. `PlanarConductionResponse` supplies ordinary and thermally refreshed current candidates without unit-field material queries or scaling accepted nonlinear fields. Actual local E.J remains heat; constitutive defects remain numerical residuals. Six manufactured regression groups cover series and parallel heterogeneous specimens, nonuniform electrode reaction, data-only source changes, both signs/zero, constant-law parity, source/work gates and cancellation/retry. Bounds remain fixed geometry and memoryless scalar response; convergence and uniqueness are not guaranteed.

## Original source response under coupled electrical drives

`PlanarConductionResponse::solve_coupled_with_materials_candidate` admits an
external owner's actual native field through the same source iteration used by
voltage/current boundaries. The returned field must borrow that exact frozen
operator; another operator's field refuses. Each trial retains its proposed
driver state, but only the state paired with the admitted constitutive field
is returned. Circuit history, source work and external work limits remain with
the caller; response iteration does not authorize advancing physical time.

The complete original laws and applied operator can be transferred together to
a continuing owner and rejoined without matrix assembly. Original source
admission rejects mismatched law/operator pairs. The original conductivity,
current-density and local-power defect gates and cancellation behavior remain
unchanged. The conduction regressions compare direct and coupled execution and
check source continuity and cancellation after a native trial.

## Original surface contacts inside nonlinear bulk response

`PlanarConductionResponse` now assembles its original scalar E/J-dependent bulk
sources and source-backed surface contacts through the same native FEEC operator.
Each actual-amplitude voltage, total-current or external-circuit trial preserves
the contact traces, voltage jumps, current continuity and exact surface heat.
Contact losses remain separate from bulk E.J in the report and are included
once in total terminal power. Ordinary static bulk admission is unchanged.

The coupled material/contact-context method retains original surface cards,
selected claims, ordered material assignments and coordinate meanings. A thermal
owner may refresh its explicitly owned temperatures and phase fractions; no
source replacement or electrical-field-dependent interface law is inferred.
The same original sources and contact geometry survive continuation. Focused
manufactured regressions compare independent nonlinear series equations under
voltage and current control, surface-source changes, both polarities, local and
total heat, and missing/replaced-source refusal. Interfaces remain stationary,
Ohmic at each supplied context and without storage; nonlinear or hysteretic
contact laws require their own physical state and closure.

## Radiation pressure from the original reflected and absorbed waves

`PlanarDielectricLayers::radiation_pressure` and
`PreparedLayerIllumination::radiation_pressure` return the cycle-mean net
electromagnetic stack/interface stress resultant per area in Pa, including
oblique incidence and coherent illumination from both faces. Each exterior
half-space retains its own original positive lossless constitutive law at the
query frequency; the two materials may differ. The exterior stress is
`<T_ij> = Re(epsilon E_i E_j* + mu H_i H_j*)/2
- delta_ij (epsilon |E|² + mu |H|²)/4` with the existing peak phasors.
The complete finite stack may contain the original lossy/metallic films.
Each face uses its own epsilon and mu in that declared Maxwell stress.
The resultant includes the interface of an empty stack and can be nonzero
without reflection or absorption. Material striction, ambient mechanical
pressure and a microscopic field/material momentum partition are not inferred.
Lossy exterior half-spaces refuse; loss inside the original finite films remains
admitted. Finite-window stress retains its separate homogeneous-source gate.

One-sided normal momentum uses
`[k_z,left*(2 P_reflected + P_absorbed)
+ (k_z,left-Re(k_z,right))*P_transmitted]/omega`;
tangential momentum uses `k_parallel*P_absorbed/omega`. Physical absorption
comes from the original positive loss modes, so weak force is retained when
subtracting incident and transmitted powers would erase it. Two-sided normal
force contracts the positive one-sided diagonal responses, each exterior's
own admittance and normal wavevector, and the actual complex scattering cross
terms; changing relative source phase can reverse its sign. Transmission's
refractive contribution has the same signed coefficient for either incident
face. A rationalized squared-index difference with a fused product retains
weak refraction when subtracting rounded normal wavevectors would erase it.
Prepared two-sided operators cache these momentum coefficients once per angle.
This force is not another thermal credit.

`PlanarAngularSpectrum::radiation_force` returns net newtons over the declared
periodic cell. Each cached Fourier mode retains its original reflected,
transmitted and absorbed watts and direction; Parseval orthogonality removes
cross-mode terms only in this full-cell integral. The original workspace budget
includes these three f64 values per mode. Queries are allocation-free O(points),
with no new FFT or constitutive solve. `AngularSpectrumCoherence` sums the
original factor forces as second moments, allocation-free O(factors*points).
Both poll cancellation in bounded batches and publish no partial result.
These are operation/storage bounds, not measured throughput.

Focused regressions cover slab resonances, oblique absorption, coherent force
reversal, independent six-field exterior-stress integration, propagation-plane
invariance, ensemble/darkness limits and source/cancellation refusal. Refracting
fixtures additionally cover impedance-matched transmission, unequal exterior
Maxwell stresses, coherent inputs, Brewster incidence, critical incidence and
TIR. The production momentum-kernel tests executed natively with the original
fs-math dependency; the full source-backed integration regressions have not
executed in this environment. This does not provide local force density,
periodic-cell torque, finite-edge scattering, transient momentum storage or
moving-body power closure. A mechanical consumer must apply this load once,
without adding a second Maxwell/Lorentz force on the same body.

## Exact planar dielectric reflection

The opt-in `plane-wave` feature exposes
`plane_wave::PlanarDielectricInterface`, an exact monochromatic boundary solve
for two homogeneous isotropic half-spaces. It retains both original material
states. Incident epsilon and both permeabilities use the existing lossless
field-storage admission. Transmitted epsilon is either the explicitly chosen
nondispersive law or the original `ResolvedDebyeCapacitor` response, converted
from terminal admittance to bulk permittivity without a second Debye law.
`DielectricResponse::Ohmic` additionally consumes an original
`ResolvedConductorState`, with absolute epsilon and mu in that same complete
bundle. Its explicit stationary isotropic law retains the selected conductivity
or resistivity basis and uses `epsilon_eff = epsilon - i sigma/omega`.
The shared outgoing-wave calculation therefore retains conductor reflection
phase, skin-depth attenuation and finite Joule absorption instead of imposing
an ideal mirror. Missing displacement/magnetic data and unsupported dynamic
source coordinates refuse; no vacuum constant or optical carrier law is inferred.
`DebyeWithOhmic` combines the original polarization and conduction responses
when both retain the identical complete resolved material/source snapshot.
Their imaginary permittivities add in parallel. Independent adoption of the
same snapshot is allowed; a foreign source bundle refuses before scattering.

`OpticalConductorResponse` adapts the original visible spectral complex-index
state to this same wave solver. Its complete retained bundle includes all n/k
knots and an explicitly sourced positive absolute permeability; no electric
constant, nonmagnetic assumption or DC-conductivity fit is inserted. Under this
owner's convention `N=n-i*k` and `epsilon=N*N/(mu*c0*c0)`, with exact SI vacuum
light speed. Original n/k are interpolated linearly within their source band;
out-of-band frequencies refuse, including zero-thickness films. Negative real
epsilon is admitted while finite nonpositive imaginary epsilon and the decaying
outgoing branch remain mandatory. The original loss is used once by half-space,
film, beam and interior-heat owners. There is no added Ohmic term or inferred
causal time-domain model. G1/G3 cases compare independent Fresnel phase and
integrated heat, thin/opaque metal films, subdivision, magnetic contrast and
vector-beam power. These synthetic fixtures do not qualify a measured specimen.

The convention is exp(+i omega t), with a forward transmitted wave
exp(-i k_z z) and Im(k_z) <= 0. TE and TM amplitudes refer to the same
tangential electric-field axis on both sides. The TM coefficient therefore
differs in sign from a rotating optical p-polarization basis. Returned power
fractions are normal Poynting fluxes; transmitted power is evaluated at the
interface, before subsequent absorption in the half-space. Total internal
reflection retains its complex phase and evanescent field. The exact critical
angle avoids division by vanishing normal wavenumber.

Near critical incidence, the shared isotropic medium retains the constitutive
product and incidence-square residuals through the existing double-double
precision rung. Only a cancellation-sized discriminant takes this path; ordinary
angles keep the scalar calculation. This distinguishes exact critical incidence
from adjacent propagating or evanescent waves, preserving weak transmitted flux,
TIR phase and finite-film tunnelling. A product rounded to unit index cannot
erase its nonzero residual. These guarantees concern the admitted binary material
ratios and angle, without increasing the accuracy of their source data.

Synthetic source-bound tests cover normal incidence, magnetic impedance,
Brewster cancellation, independent Poynting balance, total internal reflection,
critical incidence, grazing transparency, Debye phase and capacitor-geometry
invariance, retained source identity and explicit refusals. The Ohmic tests
independently integrate half-space Joule heat (including the TM longitudinal
field), compare both electrical source bases and recover the good-conductor
skin-depth limit. This is not a mesh
Maxwell solver, transient wave model, or experimental
qualification of an otherwise static material at an optical frequency.

### Visible optical sources in the original Maxwell owners

`OpticalConductorResponse` retains the original visible metal n/k behavior.
`OpticalMaterialResponse` additionally borrows a `VisibleOpticalStatePoint`
containing the original conductor, Cauchy, tabulated dielectric or Sellmeier
family. `DielectricResponse::OpticalMaterial` feeds these sources through the
same interface, coherent film, internal-field and vector-beam owners. The
incident port requires positive real permittivity and sourced permeability. The
original family can be extended with explicit permeability using
`VisibleOpticalStatePoint::with_additional_requirements`; original claims and
receipts are retained, newly required properties come from the same immutable
card and exact state, and all properties refresh together. Real positive
absolute mu is required; no nonmagnetic default or DC loss is inserted.

Dielectric sources must supply the typed spectral intensity absorption alpha
[power nepers/metre]. The bridge retains their phase-index law, interpolates
alpha (not transmission or k), and uses `k = alpha*lambda_vac/(4*pi)` with vacuum
wavelength in metres. With exp(+i omega t), `N = n-i*k` and
`epsilon = N*N/(mu*c0*c0)`. This reproduces power attenuation exp(-alpha*d);
the field amplitude decays at half that exponent. Interface reflection and
scattering cannot be relabeled as bulk absorption. RGB transmittance, missing
spectral knots and unowned frequency/field/deformation source coordinates
refuse rather than supplying a substitute spectrum or zero absorption.

Every family stays within its original 380--780 vacuum-nm spectral convention.
Cauchy B/C retain SI m²/m⁴; Sellmeier resonance parameters retain SI m². Active
Sellmeier poles anywhere in the complete band refuse, including between source
knots, and the phase index must remain finite and positive. Source refresh
readmits these restrictions. A passive homogeneous monochromatic response does
not establish causal time-domain dispersion, broadband emissivity, scattering,
anisotropy or experimental material qualification. Focused source-bound tests
compare all three dielectric families with independent Fresnel and film sums,
preserve original receipts and same-state magnetic refresh, and check RGB,
source-card, pole and adjacent-frequency refusals. The pure constitutive leaf
has separate native tests for SI dispersion, attenuation and band admission.

`PlanarDielectricInterface::try_new_with_incident_response` and
`PlanarDielectricLayers::try_new_with_incident_response` retain an original
frequency-dependent incident law as well. Their existing angle-based R/T and
film/beam paths require strictly positive real permittivity and exactly zero
incident loss, so incoming/reflected normal powers keep their additive
interpretation. Cauchy, tabulated and Sellmeier incident glass uses its
original vacuum-wavelength parameterization, without a second static epsilon
property. The incident phase wavelength is `lambda_vac/n`; its conserved
transverse wave number controls Snell refraction, total internal reflection,
film phase and the vector beam's light-cone admission.

The prepared stack caches the actual incident epsilon/mu once per frequency.
Interior fields, simultaneous illumination, cycle-mean forces and beam
reconstruction use that same cache. Frequency/thermal rebuilds retain the
original incident response rather than replacing its law with a material
snapshot. Exposure source-identity checks include the incident law; finite-time
stress/momentum observations still require an explicitly nondispersive law in
the observed half-space, even if a spectral law happens to be lossless at all
queried frequencies. Synthetic G1/G3 regressions cover the three dielectric
families, independent Snell/Fresnel and Fourier-vector fields, homogeneous
film phase and internal E/H, both illumination ports, wavelength-dependent
light-cone admission, absorbing/out-of-band refusals and the original static
constructor equivalence. These remain monochromatic-source numerical checks,
not a causal pulse-propagation or measured-material claim.

### Reflection launched inside an absorbing medium

`PlanarDielectricInterface::response_at_transverse_wavenumber` separately
admits passive absorbing incidence at a single stationary isotropic interface.
It accepts signed real conserved kx [rad/m], ky=0, frequency and actual complex
peak [Ex,Ey] [V/m] at z=0. A real angle is not substituted for the complex
incident wave direction. Both original material laws and sourced permeabilities
remain attached, including Debye, Ohmic, Lorentz and same-snapshot combined
responses; original spectral-band and material-domain refusals remain active.

The result retains full Cartesian incident/reflected/transmitted E and H,
complex normal wave numbers, constitutive values and tangential reflection and
transmission gains. Four unnormalized normal powers [W/m²] are explicit:
incident, signed reflected, incident/reflected interference, and transmitted.
Their boundary balance is `P_incident + P_reflected + P_interference =
P_transmitted`, with reflected power nonpositive. The cross term can have
either sign. Transmitted power can exceed the separately evaluated incident
power without gain; neither that result nor the interference is clipped or
relabelled as interface heat. Subsequent bulk absorption remains separate.
These are local boundary fields from an authored finite-source amplitude;
no upstream launch distance or attenuation correction is inferred.

The shared outgoing branch and bounded Fresnel pair are reused. Binary
constitutive references and compensated dispersion preserve the original real
transverse wave number without forming an overflowing omega² or epsilon*mu.
Power-conjugate magnetic traces keep finite longitudinal TM fields at a
critical exit; evanescent exit fields retain their phase and zero normal flux.
Incoming pure evanescent/grazing ports, gain and unrepresentable fields/powers
refuse. Zero authored illumination retains zero fields after material admission.
Actual plasma fields refuse because the two-wave incident field has no original
drift-amplitude admission owner. Existing angle, film and beam ports retain
their lossless-incidence restrictions.

G1/G3 native kernel tests cover independent complex Fresnel equations, both
Maxwell curls, full boundary traces, signed Poynting flux, a closed-form
absorbing-incidence interference case, exact critical/evanescent exits,
lossless-owner reduction, material-product residuals and field/frequency scaling.
Source-backed integration cases specify the original Debye, both Ohmic source
bases and Lorentz/combined-law APIs; their Cargo execution remains pending in
this resource-limited session. The power convention follows the Poynting
derivation in [Byrnes, sections 4.3–4.4.1 and Appendix B](https://arxiv.org/pdf/1603.02720),
with this owner's conjugate time convention and explicitly sourced permeability.

### Coherent films illuminated from an absorbing medium

`plane_wave::layers::absorbing_incidence::PlanarAbsorbingDielectricLayers`
extends the explicit-real-kx owner to up to 64 coherent finite films. It retains
all original incident, film and exit laws and thicknesses; resolves their actual
permittivities and sourced permeabilities once at the requested frequency; and
prepares both TM/TE fields together. `illuminate([Ex,Ey])` binds actual complex
peak fields at the first face. Repeated illumination and `layer_field` depth
queries allocate nothing and do not rebuild the stack or reevaluate a material.
Empty stacks reduce to the absorbing single interface; zero-thickness films
remain identities after their original sources are admitted.

The shared bounded E/H slab recursion runs in a real numerical reference port,
then matches its entrance against the original complex incident admittance.
That reference is only a coordinate scale, never a replacement constitutive
law. Decaying propagation, compensated absolute-q phase, exact-critical entire
fields, and opaque limits reuse the existing slab and depth-field owners.
Positive even/odd dissipative modes give each film's heat independently of
exterior subtraction, including weak loss and coherent standing waves. Full
Cartesian E/H and local `-omega Im(epsilon) |E|²/2` heating include longitudinal
TM electric fields. The result satisfies
`P_incident + P_reflected + P_interference = P_exit + sum(P_film_heat)`.
Entrance interference remains signed and is never assigned to film heat;
exit flux excludes subsequent half-space absorption. Film heat is W/m² and
local heat is W/m³. The consuming thermal model owns area and time conversion.

Native G1/G3 regressions compare against independent dimensional transfer
matrices and six-component fields, continuity of tangential E/H and normal
D/B, independently integrated volume loss, empty-stack interference, film
subdivision, exact-critical fields, opaque loss and weak positive heat.
Existing single-interface and depth-field regressions also pass with the shared
helpers. The public adapter was compiled and exercised with source-resolution
plumbing stubbed; a source-card Debye/Ohmic integration regression is present,
but its Cargo execution remains pending under this session's resource limits.

This owner covers one-sided monochromatic stationary isotropic planar films.
It supplies no finite-coherence, roughness, anisotropic, magnetic-loss,
transient, absorbing angular-spectrum or upstream launch model.
Two-sided absorbing incidence uses the separate opt-in preparation below. Pure evanescent/grazing incidence, active media, unrepresentable
contrasts/fields/powers and plasma without original whole-field drift admission
refuse. The older angle, beam and simultaneous-illumination owners retain their
lossless incident-port admission.

### Coherent illumination from both absorbing half-spaces

`PlanarAbsorbingDielectricLayers::prepare_two_sided` adds simultaneous coherent
sources at both physical faces of the original stack. Both incoming modes must
carry positive normal power. The existing one-sided owner continues to admit
critical/evanescent exits, which cannot supply the second incoming wave and
therefore refuse this opt-in preparation. Sources share frequency, signed real
kx and clock; each actual complex peak [Ex,Ey] is specified locally at its own
face on the same fixed Cartesian electric axes.

One additional bounded sweep reverses the retained original film order using
the SAME real numerical reference. It does not reevaluate any material law.
The odd dissipative-mode phase reverses with the film orientation, so coherent
source amplitudes are combined before their positive squared loss is taken.
The right-source interior fields are reoriented from their original exact face
traces; queries retain the original global +z direction and depth from each
film's first face, without a rounded thickness-minus-depth phase conversion.
Preparation and evaluation are O(films), with allocation-free repeated source
phase/amplitude and interior queries and the existing 64-film cap.

The prepared owner exposes the complex physical port admittances and the full
fixed-axis tangential-E scattering matrix. Actual illumination returns incoming,
outgoing and total Cartesian E/H at both faces, per-film heat, and local depth
fields/heating. Each face reports separate nonnegative incoming/outgoing
powers, a signed interference contribution into the stack, and signed net
inward flux. With inward coordinates `E=a+b`, `H=Y(a-b)`, the cross flux is
`Im(Y) Im(b conj(a))`; it includes transmission from the opposite source and
does not require an undefined reflection ratio when a source is zero. The
physical balance is `net_inward_first + net_inward_last = sum(film_heat)`.
Neither cross term is relabelled as film heat or clipped. Total boundary E/H
retain independent stable transfer numerators, including near cancellation.

Native regressions cover an absorbing matched film's closed-form phase-controlled
heating and both cross fluxes, independent dimensional two-source transfer
matrices, six-component fields and interface continuity, integrated volume
loss, either one-sided reduction, zero lossless-port cross flux, exact-critical
and opaque films, weak odd-mode heat, dark fields, empty stacks and exterior
admission. The public adapter is compiled and exercised in the native probe
with only source-resolution/error plumbing stubbed; original Debye source-card
and lossless-owner parity tests are present but their Cargo execution remains
pending under the session's resource constraint. These are G1/G3 numerical
claims for stationary coherent isotropic planar films, with no finite-coherence,
roughness, anisotropy, magnetic loss, transient or measured-material claim.

### Coherent finite dielectric films

`plane_wave::layers::PlanarDielectricLayers` prepares up to 64 ordered finite
films between an explicit lossless incident source and an explicit lossless,
Debye or Ohmic exit source. Every film retains its original `DielectricResponse`,
thickness and sourced permeability. Constitutive data are cached once at the
selected frequency; angle/polarization queries are allocation-free. Zero
thickness is an exact identity after source admission. The half-space owner
and this stack share material-contrast and outgoing-normal-wave calculations.

Composition retains bounded homogeneous E/H states in common real reference
ports, then reconstructs forward scattering gains. It preserves resonant
coupling without subtracting products of rounded nearly unit reflections.
Propagation uses only the decaying/oscillatory exponential; an analytic
zero-normal-wave limit includes exact internal critical angles. Returned
tangential-E phases refer to the first and last film faces. Total internal reflection, finite-gap
evanescent tunnelling, Fabry-Perot resonance and coherent antireflection are
retained. Exit power is the normal Poynting flux entering the last half-space;
finite-film absorption is reported separately. Dissipative even/odd slab modes
retain weak Debye loss that subtracting reflected/transmitted powers would
round away. Wholly lossless films have exactly zero dissipative work.

The existing plane-wave tests add half-space/zero-thickness reductions,
quarter-wave and resonant phases, grazing subdivision and distinct-material
resonances, independent E/H transfer matrices with magnetic media, original
Debye/geometry invariance, weak-film loss, critical films/exits and thick absorbing/evanescent limits. This remains a coherent
homogeneous isotropic planar model: finite coherence length, roughness,
anisotropy and mesh Maxwell evolution are outside its scope.

Ohmic films use this same prepared scattering, internal-field and vector-beam
path. The local cycle-mean heat is `sigma |E|^2/2`, already included in physical
film absorption; it is never deposited a second time. Focused synthetic cases
cover the resistive-sheet limit (including 50% absorption for sheet conductance
`2/Z`), opaque-metal reflection and vector-beam volume heat versus boundary
Poynting flux. These numerical fixtures do not qualify DC electrical sources
at optical frequencies or include carrier dispersion, magnetic loss or roughness.
The combined-law case also checks reflection phase and integrated Debye plus
Joule heat across frequency, using independent source equations and rejecting
a conductivity attached from another material snapshot.

### Incoherent propagation between native coherent coatings

`plane_wave::layers::IncoherentLayerCascade` composes original
`PreparedLayerIllumination` blocks across explicitly incoherent, lossless
connections. Each block keeps its coherent thin-film interference and original
material/geometry owners. Adjacent blocks must share the exact connecting
source law, frequency, polarization and conserved tangential wavevector.
The native two-port admission requires positive propagating lossless exterior
modes. Evanescent and critical tunnelling belongs inside a coherent block;
an absorbing connecting medium is not silently assigned additive power ports.

Preparation caches the normal-power scattering matrix and left/right absorption
fractions of every original film. Native positive dissipative modes supply
film loss, including longitudinal TM contributions. Multiple reflected passes
are summed through positive transmission-plus-absorption escape factors rather
than `1 - R_left * R_right`, preserving weak transmission when both reflected
fractions round to one. Native opaque limits can isolate an undriven cavity;
no arbitrary trapped initial energy is inserted. Work is O(blocks + films),
bounded to 65 blocks and 64 total films. Repeated independent-face illumination
is allocation-free O(films), with physical incoming/outgoing W/m² and film heat
in original block/layer order. It does not repeat a field solve, material
lookup, phase sampling or trigonometry. Refused illumination leaves the heat
output unchanged; any thermal consumer still owns its area/time conversion.

Complete decorrelation between distinct reflected traversals is an authored
model, not an inference from unknown phase, physical thickness or linewidth.
The API does not supply a partial temporal-coherence model or complex outgoing
fields. A single connection agrees with uniform coherent round-trip-phase
averaging; several connections represent fully incoherent intensity transport,
not arbitrary averaging of correlated static thicknesses. No roughness,
magnetic loss, thermal evolution or measured-material validation is inferred.
G1/G3 regressions in `tests/plane_wave/incoherent.rs` compare the analytic thick
glass plate, lossy coherent coatings against the original phase-averaged fields
and individual heat, lossless multi-interface transport resistance, weak-escape
and opaque limits. G0 cases retain original source/frequency/direction domains,
bounded preparation and transactional illumination refusal.

### Coherent vector beam diffraction through the same films

`plane_wave::angular_spectrum::PlanarAngularSpectrum` accepts complex peak Ex/Ey
samples at the incident face of an original `PlanarDielectricLayers`. The input
is one explicitly periodic, band-limited rectangular cell in physical metres;
x varies fastest. Each axis admits 1..=4096 samples at its original period.
The existing `GeneralFft` uses cached Bluestein convolution for arbitrary
lengths and preserves power-of-two transform bits. Odd axes retain both edge
orders; even axes retain negative Nyquist. A pencil uses at most two 8192-point
transforms between cancellation polls. Every nonzero incident Fourier mode
must lie strictly inside the incident light cone. Exactly zero unsupported
bins may oversample that spectrum; they retain their physical wavevectors and
FFT positions without evaluating an invalid film response or dividing by an
absent longitudinal wave number. No tolerance removes small nonzero modes,
including analysis-FFT residuals. No modes are silently removed, windowed or
renormalized. A caller approximating an isolated beam must resolve
its bandwidth and converge the cell/padding independently; an authored hard
aperture's unresolved high spatial frequencies are not synthesized.

The original `fs-fft` plans analyze both complex sampled components once.
`prepare_modes` instead accepts authored peak Ex/Ey Fourier amplitudes in V/m,
with the same signed FFT order and negative Nyquist convention, without an
extra 1/N factor or analysis FFT. This permits explicit band-limited
interference fields on fine observation grids without inferring a source
filter. Every supported mode uses its conserved transverse wave number and
TE/TM response
from the original source/material/film owner. Complex amplitudes combine before
forming intensity. Source epsilon and mu remain distinct; scalar refractive
index or normal-incidence transmission cannot replace their vector response.
For the existing exp(+i omega t) convention, normalized inverse FFTs analyze
the field and unnormalized forward FFTs reconstruct exp(-i k.r) modes.

`PlanarAngularSpectrum::prepare_aplanatic_pupil` adds the vectorial
Richards–Wolf sine-condition source. `AplanaticPupil` declares a strictly
sub-grazing sine aperture `sin(theta_max)=NA/n`, the geometric focus in the
homogeneous incident medium, and positive cycle-mean incident watts. The
once-per-admitted-point callback receives direction cosines and radius-one
pupil coordinates, and returns a dimensionless complex x/y Jones shape.
Spatial polarization, apodization, stops and authored aberration phases are
retained. The requested watts refer to the beam after that authored shape;
they do not infer the upstream power or transmission of a physical stop.

The lens map rotates radial polarization into the ray's meridional plane and
keeps its azimuthal component. The reference-sphere `sqrt(cos(theta))` factor
and Cartesian solid-angle Jacobian `1/cos(theta)` combine to
`1/sqrt(cos(theta))` in each Fourier coefficient. Modal longitudinal E/H then
come from the original Maxwell owner; no scalar or paraxial substitution is
made. Shape magnitudes are scaled before squaring, with a double-double power
reduction and exact modal Poynting normalization. A compensated phase aligns
the original wavevectors at the requested focus. The common lens-to-focus
phase is explicitly zero; authored Jones phase remains coherent with other
beams. Subsequent film scattering can move or aberrate the transmitted focus.

The complete circular angular support must fit strictly inside both grid
Nyquist limits, and grazing, dark positive-power pupils and unrepresentable
amplitudes refuse. Preparation takes three bounded linear passes, with at
most 32 pupil samples/rotations or 256 final scalings between cancellation
polls. The original workspace budget additionally includes N complex-pair
source samples during construction. The returned original beam needs no new
field-query workspace or lens evaluation. This is a Cartesian quadrature of
an ideal aplanatic focal model: physical-cell/angular-resolution and output
sample-spacing convergence remain caller-owned. No material lens, arbitrary
finite-screen Maxwell solution or subwavelength metasurface is claimed.
Reference: Liu et al., <https://arxiv.org/html/2502.03170v1>, equations (14),
(24). G1/G2/G3 kernel tests compare the high-NA on-axis peak against an
independent continuous Richards–Wolf integral over three angular-grid rungs,
and check radial/azimuthal longitudinal fields, translation, coherent phase,
six-component power/transversality and extreme Jones scaling. Native
`tests/plane_wave/angular_spectrum/aplanatic.rs` exercises the source-bound
constructor, original magnetic impedance, films, field queries and budgets.

Each subsequent plane reuses the cached modes and plans, applies exact kz
propagation, enforces k.E=0, and obtains H=(k cross E)/(omega mu). It returns
all six complex field components and physical cycle-mean Poynting flux.
Incident, reflected, their coherent sum, and transmitted fields have explicit
face/distance conventions. Film-induced evanescent fields decay without a
growing transfer matrix. Deep transmission applies the physical E/H amplitudes
before a unit transfer can underflow, with one phase and two split attenuation
factors per mode. Fields, probes, coherence and exposure preserve representable
evanescent/absorbing tails; ordinary-path arithmetic is unchanged. G1/G3 tests
cover tail values, depth composition and ordinary-path parity.
The original film recursion retains transmitted
power-conjugate H alongside E, relative to the original incident admittance.
Transmitted fields use these two tangential traces without dividing by kz.
At exact critical refraction, the finite longitudinal TM electric field and
TE magnetic field remain present even though normal transmitted power is zero.
No perturbation of the incidence angle or artificial transmitted-power floor
is introduced. The modal storage budget includes the additional two complex
magnetic traces; further-plane reconstruction still requires no film solve.

`PlanarAngularSpectrum::probe_fields` evaluates arbitrary exterior observations
directly, including tilted/irregular receiver layouts and different sides or
depths in one query. All six peak phasors and the full complex cycle-mean
Poynting vector use the original native Maxwell modes. Its real part is signed
energy transport; its imaginary part retains reactive evanescent flow. Probe
coordinates observe the same periodic envelope and authored carrier as the FFT
planes, with no interpolation, paraxial approximation or implicit detector area.
`PreparedAngularSpectrumLayer::probe_fields` supplies the same nonuniform field
query inside one original finite film, plus local constitutive heat from the
complete coherent electric field. Film depth includes either boundary; normal
D/B and tangential E/H retain their respective material continuity laws.
Neither query reevaluates material laws or rebuilds the prepared film response.
For P probes and A exactly nonzero incident modes, direct synthesis costs A*P
native mode/phase evaluations, two linear source-array scans and O(P) carrier
and observable work. No tolerance
removes a source coefficient. Point, modal-work and numerical-storage limits
precede allocation; storage includes the original owner's conservative workspace,
A indices, P outputs and one six-component double-double accumulator. Cancellation
returns no partial result. Uniform planes retain the faster existing FFT path.

Modewise power integrates over the original cell: reflected outward power,
exit-face transmitted power and physical film absorption are separate. Loss
after propagation into an absorbing exit is computed from the actual field,
without relabeling it as film absorption. Local flux retains interference and
its sign. Numerical-storage and FFT-work budgets precede allocation; source
response work is bounded by the admitted grid and original 64-layer maximum.
Storage/work bounds include both axis plans, Bluestein chirps/kernels,
convolution transforms and reusable scratch. Original power-of-two allowances
are unchanged. Dormant
bins retain empty entries in the shared film and two-sided response caches;
field, heat, force, coherence and exposure queries preserve their zero fields.
Cancellation polls between rows, pencils and bounded mode batches, publishing
no partial field. Further planes require six separable FFTs and linear modal
work, with no repeated constitutive or film solve. This is an operation-count
description, not a measured throughput claim.

Native leaf tests exercise both Maxwell curls/divergences, polarization and
signed power, and evanescent decay. Source-backed integration regressions cover
original-film reduction, an independent direct Fourier diffraction series,
TIR, exact critical TE/TM fields and finite tangential flux, independent
critical-load film matrices, both one-sided critical limits, Debye film/exit
loss, budgets and cancellation. Oversampling regressions additionally cover
uniform sampled/modal equivalence, analytic opposite-mode fringes and finite
apertures, coherent versus independent-factor observations, interior electric
nodes under counterillumination, smallest-subnormal unsupported-mode refusal,
budget boundaries, cancellation and replay. Their Cargo execution
requires the workspace constellation dependencies. This feature does not add
transients, incident near-field excitation, rough/finite film edges, arbitrary
three-dimensional scattering, or empirical material-frequency qualification.

### Arbitrary-angle carriers for refracting vector beams

`PlanarAngularSpectrum::prepare_with_carrier` accepts an explicit transverse
wavevector in rad/m and samples of a periodic tangential-electric envelope.
The physical field is that envelope times `exp(-i*k_carrier·x)`. Every
Fourier wavevector receives the carrier before its original TE/TM film solve,
Snell refraction, vector Maxwell reconstruction and propagation. Incidence
angles are therefore continuous rather than quantized to the envelope's FFT
bins. Every occupied shifted mode must fit inside the incident light cone;
exactly empty unsupported bins retain their positions without a film solve.
Critical and evanescent exit waves retain their existing physical treatment.
The original `prepare` is exactly the zero-carrier case.

`prepare_modes_with_carrier` accepts the same authored peak Fourier amplitudes
as `prepare_modes`, with an explicit carrier before physical mode admission.
This permits exactly band-limited oblique beams on fine observation grids
without analysis-FFT leakage into unsupported bins. It neither discards nor
renormalizes any nonzero mode; `prepare_modes` is its exact zero-carrier case.
One-sided field queries skip propagation of identically empty modes, while
retaining supported port geometry so subsequent right-face illumination can
populate those bins. Focused controls compare sparse shifted modes with
analytic vector Fresnel fields, whole-cell power and finite-detector fringes,
and check exact-support refusals, cancellation and right-only illumination.

Carrier phases retain double-double products during argument reduction and
are cached once, adding one complex value per point to the admitted workspace.
Field queries add O(points) complex multiplications without new trigonometric
work, constitutive queries or film solves. Internal fields, cycle-mean
absorption, finite apertures, momentum, partial coherence and simultaneous
illumination use the same shifted modes; right-face inputs use the same
envelope convention. Shared response caches require the same carrier.
Cycle-mean aperture kernels cancel the common carrier analytically and retain
the original envelope-bin differences, including separations below the
carrier's floating-point spacing. Direct mutual-coherence observations retain
the envelope phase separately and apply the carrier once per observation.
Finite-window exposures admit a shared carrier across their original lines
and covariance banks, retaining both real-field spatial sum and difference
terms. No paraxial approximation, discarded modes, near-field incident
excitation or finite-edge diffraction is introduced. Focused controls cover
arbitrary oblique Snell/Fresnel fields, Brewster cancellation, TIR decay,
interior and two-sided fields, source-cache identity and resource refusals.

### Finite-window observations of refracting carrier beams

The original exposure owner now admits nonzero transverse carriers shared by
all lines and covariance factors. The same grouped Maxwell fields, source
clock, material responses and uniform/triangular shutter kernels feed exterior
observations, simultaneous illumination, film interiors, finite apertures,
polarization/control observations, energy flux, Maxwell stress, field momentum
and its weighted rate. Real
sources retain sum and difference terms; only an authored proper-complex
factor removes the real-field sum term.

Spatial kernels retain double-double shifted wavenumbers and the original
envelope offsets. Continuous cell power, stress and momentum now share a
separable centered cosine/sine Gram operator for nonzero real carriers.
Its x and y contractions cost O(F*N*(width+height)), followed by O(N*F²)
frequency products, replacing the former O(N²*F²) full-cell contraction.
The source-bank pair-work getter and cap count N*F² + F*N*(width+height)
when the carrier is nonzero, and retain exactly N*F² at zero carrier.
Axis matrices, center phases, retained DD coefficients, vector headers and
reusable dual-Gram scratch are admitted before allocation. Finite DD results
and cancellation are checked before publication. Proper-bank cell products
cancel the common carrier analytically and keep sparse O(N*F²) work.
Ensemble query admission counts N*F² for each proper bank and adds carrier
axis work only for the mean and real factors. Workspace remains conservatively
bounded by the generic prepared bank; preparing that reusable owner still
admits its separate deterministic carrier cache.
No extra material query, film solve or synthesis FFT is introduced.

Continuous film heat uses original envelope offsets for conjugate products
and DD carrier-shifted arguments for real sum products; the common carrier
center phase rotates completed sum coefficients. Nonzero real banks on
`ContinuousPeriodicCell` use the exact rectangle over the declared cell,
charging D(N²+1) times their frequency-pair count. Proper banks retain D*N
work per frequency pair, independently within a mixed bank. Original sampled
area and depth quadrature semantics remain; irreversible Debye/Joule heat is
still separated from stored/reversible material energy.

Individual stress planes and finite patches retain all six independent tensor
components. Field-momentum volume-rate and joint stress/momentum observations
admit the same carrier source. The joint query shares the spatially contracted
DD fields between both observables, retaining carrier-sum and beat derivatives
under the original uniform or triangular shutter. The periodic two-plane
boundary-resultant shortcut still refuses nonzero carriers because lateral stresses cannot
be assumed to cancel for finite-window real fields. These observations do not
infer a closed-surface mechanical force or thermal feedback.
Focused controls compare independent spatial/time Fresnel quadrature, two-face
standing waves, analytic tensor stresses and actual material-reservoir heat.
Additional regressions use independent positive-time sinc/phase products to
check continuous power, all stress/momentum components and weighted rates;
they retain mixed real/proper factors and verify exact new work admission,
cancellation and reuse. A separate public control test reconstructs actual
complex-weighted beams and checks energy gradients/Hessian actions against
finite-aperture energy. An original-film Debye/Joule test checks real versus
proper carrier heat and exact film work. These new consumer regressions await
native execution; the shared shifted harmonic-kernel tests passed in a bounded
build using the actual mathematical kernel.

### Finite detector apertures preserve spatial interference

`PlanarAngularSpectrum::measure_aperture` observes the original coherent vector
field through an explicitly positioned, axis-aligned rectangular area in the
observation plane. The result separates area-mean vector Poynting flux (W/m²)
from signed integrated normal power (W). It retains every cross-mode term in
E cross H* before integration. Detector width and position therefore resolve
fringe visibility and washout; neither an average amplitude squared nor a
whole-cell sum substitutes for the actual detector.

The original modal-field reconstruction supplies incident, reflected, total
incident-side or transmitted E/H, including critical and evanescent exit
fields. Analytic rectangular Fourier integrals reuse the shared
`fs-fft::harmonic_window` difference kernel, with coordinates in metres and
wavenumbers in radians/metre. Actual admitted axis-wavevector differences,
source products and reductions retain double-double precision. There is one
cycle-mean factor of one half. Reflected power and local energy backflow keep
their physical negative signs; no detector absorption law is inferred.

A query factors the exact Cartesian rectangle kernel into two double-double
axis actions. For an N=W*H grid it performs N*(W+H) three-component magnetic-
vector contractions and N final vector cross products, retaining every
original cross-mode term algebraically. The existing `maximum_mode_pairs`
field bounds those actual axis contraction terms. The query retains
O(W²+H²) cached spatial kernels and one O(N) complex double-double intermediate
grid, with no FFT or new material/film solve. Explicit work and logical-storage
limits admit the computation before reconstruction; storage includes the
original beam workspace, both axis caches and the intermediate vectors.
Cancellation polls during cache construction, modal reconstruction and bounded
batches of both axis actions, returning no partial measurement. These are
operation/storage bounds, not a throughput claim. The source remains the
authored periodic field;
apertures outside its cell observe its periodic continuation. Isolated-edge
convergence, detector back-action, quantum efficiency and temporal exposure
remain separate questions.

An O(N) preflight conditionally rescales extreme magnetic amplitudes by an
exact power of two before the axis actions, restoring both DD components only
after the electric contraction and cycle-mean factor. Inexact source scaling,
lost DD halving components and unrepresentable restored flux refuse explicitly.
This bounds intermediate magnetic-sum overflow; it is not a general extension
of the DD exponent range for every possible reassociation.

Focused native kernel tests cover analytic fringe position/width response,
positive subwavelength dark-fringe power, independent direct vector-field
area quadrature including evanescent and backward waves, whole-cell
orthogonality, detector subdivision, darkness, resource limits and cancellation.
They also compare a nonsquare grid with the complete all-pairs expression,
including distant-detector phases, positive two-axis dark-fringe power, exact
axis-work admission, and cancellation and replay during both separable actions.
Additional range regressions compare huge-H/small-E transport with a finite
all-pairs reference and check explicit scaling, halving and output refusals.
Source-backed integration tests exercise the public original-film/beam adapter,
TIR, reflected power and refusal replay; their full Cargo execution requires
the constellation dependencies.

### Diffraction-resolved fields and dielectric heating inside finite films

`PlanarDielectricLayers::prepare_layer_field` retains one original film and
polarization under unit tangential-E incidence. It reuses the existing bounded
slab/state recursion, including source epsilon/mu, lossy exits, evanescent
tunnelling, zero-thickness identities and opaque limits. Both physical E/H
boundary traces come from that same stack. Interior waves decay from their
nearest face; near a critical normal wave number, an entire cosine/sinc E/H
propagator gives the exact linear kz=0 limit. No growing exponential, division
by a vanishing kz, independent material law or second film solver is used.
Depth is measured from the selected film's first face. The public tangential
H axis is power-conjugate to E, with positive z directed into the stack.

`PlanarAngularSpectrum::prepare_layer` caches those responses for the actual
coherent beam, with the original incidence cosines preserved exactly. Only
the selected layer is retained, in O(grid points) storage. Further depth
planes reuse the original six FFTs and reconstruct longitudinal E/H from
tangential E AND H, so internal critical-angle fields remain representable.
At an interface, longitudinal E belongs to the selected material; normal D/B
and tangential E/H satisfy their distinct continuity conditions.

The returned `AngularSpectrumLayerField` contains all complex vector fields,
signed normal Poynting flux and local cycle-mean dielectric heat in W/m³:
`Q = -omega*Im(epsilon)*|E|²/2` for this exp(+i omega t), peak-phasor convention.
All spatial modes and polarizations interfere before Q is computed. The
transverse heat integral has units W/m of layer depth. Separately, the
prepared layer reports total absorbed watts from the original positive
even/odd loss factors, retaining weak absorption without flux subtraction or
depth quadrature. Lossless films return exactly zero heat. Source mu is real
and lossless; no magnetic dissipation is inferred.

Preparation explicitly budgets both polarization/layer/mode work and logical
numerical storage, including the borrowed beam workspace, selected response
cache, fixed stack scratch and one E/H/heat plane. Query FFT work has its own
cap. Cancellation polls between bounded source solves, modal batches and FFT
pencils; refusal publishes no partial result and leaves prepared owners
unchanged. Reusing a layer for many depths requires no new material or stack
solve. These are storage/operation bounds, not measured throughput claims.

Native production-kernel regressions cover critical/near-critical/opaque
limits, independently constructed counterpropagating waves, both Maxwell
curls and divergences, and volume dielectric heating versus signed Poynting
loss. Source-backed integration regressions additionally exercise beam/film
interface continuity, coherent spatial heating fringes, independent depth
integration, original positive film absorption, budgets and cancellation;
their full Cargo execution requires the workspace constellation dependencies.
This extension remains a monochromatic periodic planar-field model. It does
not evolve a temperature, claim isolated finite-edge convergence, or add a
nonlinear, anisotropic or three-dimensional Maxwell solver. Thermal consumers
must debit the reported absorption only once.

### Coherent vector illumination from both physical faces

`PlanarAngularSpectrum::prepare_illumination` adds an explicitly coherent
right-face Ex/Ey grid to an original left-face beam. Both fields share the
original frequency, transverse coordinates, periodic cell and phasor clock.
Every supported spatial bin uses the existing `PreparedLayerIllumination` for
both TE and TM, with fixed tangential axes on both faces.
`prepare_illumination_modes` accepts the right-face peak Fourier amplitudes
directly, with the same ordering and normalization as `prepare_modes`.
Both exterior media must remain lossless, including for a dark source. Every
bin occupied from either face must propagate in both exterior ports; an
occupied critical or evanescent port refuses even when the other face is dark.
Only bins exactly empty on both faces may bypass an unsupported two-port solve,
retaining their original grid positions. No nonzero coefficient is discarded.

The prepared owner reports positive incoming/outgoing watts per face and
original positive per-film absorption. Each outgoing field is the coherent sum
of reflection and opposite-face transmission. Exterior queries reconstruct
complete E/H with the correct signed normal wavevector, including total
standing-wave fields at either face. Interior queries contract the cached
native source columns and synthesize every spatial/polarization contribution
before applying the same local dielectric loss used by the one-sided owner.
Counterbeam phase therefore changes spatial heat without deriving weak loss
from an incident-minus-outgoing remainder. Integrated radiation force reuses
each original lossless exterior's Maxwell-stress law, including refractive
momentum between unlike materials; force is not another heat or mechanical-work
credit.

Preparation costs O(points*layers) work/storage, explicitly bounded before
allocation. It borrows the original beam and FFT plans and retains both native
source columns per polarization/bin for subsequent O(points) field assembly
and six original FFTs. No material law or stack scattering solve is repeated
for a depth or exterior query. Force queries allocate nothing and contract
the native loss modes in O(points*layers). Cancellation publishes no partial
owner or field. These are operation bounds, not measured throughput claims.

Four focused native integration groups compare independent direct Fourier and
interface equations, characteristic-matrix interior fields and depth-integrated
heat, exterior Maxwell stress, phase-dependent absorption, one-sided reduction,
exact budget boundaries and source/cancellation refusal. They remain unexecuted
in the implementation environment because dependency compilation exceeds its
memory limit. The model remains monochromatic, planar and periodic: finite
aperture edges, evanescent illumination sources, transient momentum storage,
thermal state evolution and arbitrary coherence between different frequencies
are not inferred by this addition.

### Exact depth-integrated heat for monochromatic beams

`PreparedAngularSpectrumLayer::absorbed_power_in_depth_interval_w` and
`PreparedAngularSpectrumIllumination::layer_absorbed_power_in_depth_interval_w`
return physical watts over one film depth interval and the complete original
periodic transverse cell. Each polarization/bin reuses the native analytic
interval loss, retaining coherent counterbeam and internal standing-wave
interference before its power is formed. Transverse Fourier orthogonality and
TE/TM electric-field orthogonality then permit the full-cell power sum.
These queries allocate nothing, require O(grid points) work, poll cancellation,
and use neither depth quadrature nor FFTs. Endpoint bounds and original source
admission remain attached. Refusal leaves the prepared owners unchanged.

The returned watts already include area and depth integration; a thermal
consumer supplies time and deposits each disjoint interval once. Finite
transverse apertures, spatial heat maps and finite-window multifrequency
heating retain their existing field/exposure paths. Focused regressions compare
independent local Poynting work and refined spatial/depth heat quadrature,
interval additivity, one-face reduction, invalid cells and cancellation.

### Spatial partial coherence through the same diffracting films

`angular_spectrum::AngularSpectrumCoherence` borrows original prepared beams
as explicit positive cross-spectral factors: `W(r,r') = sum B_m(r) B_m(r')*`.
Every factor retains the exact same source operator address and periodic grid;
its field amplitude includes its square-root statistical weight. An explicit
anchor fixes that operator/grid without adding power, so an empty factor bank
is exact darkness. No independence between Fourier bins is assumed: the
existing six-component Maxwell reconstruction interferes each factor first,
then the mutually uncorrelated factors add as second moments. The output
retains the complete Hermitian electric coherency tensor, componentwise
magnetic RMS squared and all three complex Poynting components. Peak-phasor
time averaging contributes one factor of one half. Electric coherency is in
V²/m²; physical energy flux comes from E cross H*, not an impedance guess.

The selected-layer consumer shares ONE original TE/TM response cache across
the entire factor bank. Local positive dielectric heat is evaluated after
within-factor field interference, and original positive film-loss modes
supply the integrated watts without subtracting fluxes or quadrature. Local
and integrated heat remain the same physical loss, never two thermal inputs.
Repeated depth queries do not repeat constitutive evaluation or layer solves.
No new source or film kernel is introduced.

Caller limits admit at most 4096 factors, total factor/point pairs, total FFT
work and conservative numerical workspace before allocation. Queries stream
one reconstructed factor at a time: additional storage is O(grid points),
and FFT work is O(factors × points × log(points)). The workspace bound includes
one original beam's complete workspace, result/accumulation planes and the
shared layer cache; other caller-owned factor beams, earlier retained results
and allocator overhead are excluded. Cancellation polls through the original
FFTs and bounded moment-reduction batches, returning no partial statistics.

G1/G3 tests check analytic fringe visibility and signed Maxwell flux,
independent Airy standing-wave heat, depth-integrated absorption, complex
unitary source-factor changes, polarization, zero-rank and lossless limits,
exact source binding, work limits and cancellation. This is a monochromatic
spatial ensemble; finite bandwidth, temporal beats, source-coherence inference,
temperature evolution and isolated-aperture convergence remain outside it.

### Finite exposures of prescribed spectral lines

`angular_spectrum::AngularSpectrumExposure` borrows the original prepared
vector beams at their authored frequencies. Every line must retain identical
original incident/exit/film material owners, layer order, thicknesses and
transverse grid. Each line includes its physical peak amplitude and phase at
time zero; neither a bandwidth nor a spectral-density normalization is inferred.
The material/film responses remain individually resolved at each frequency.

The shared `fs-fft::harmonic_window` integrates real E/H products over an
explicit rectangular time window, retaining both beat and sum frequencies.
Short shutters preserve carrier phase; long shutters suppress resolved beats
according to the actual window. Equal-frequency lines combine with compensated
field accumulation before powers and share one six-component FFT synthesis.
An in-place O(L log L) preparation groups source indices by frequency while
preserving their original order within each group; subsequent queries visit
each source exactly once. All reserved index/offset storage is budgeted.
Temporal kernels are cached per distinct frequency pair. For L lines, F distinct
frequencies and N spatial samples, zero-carrier query work is
O(L*N + F*N*log N + F²*N); nonzero real carriers add the separable pass above,
with explicit line, modal, pair, FFT and numerical-storage budgets. This is an
operation count, not measured throughput. Cancellation returns no partial field.

Results contain real componentwise mean-square E/H, signed mean Poynting flux
and integrated mean normal watts over the original periodic cell. The integrated
watts use exact products of the represented signed Fourier modes before their
quadratic products are sampled. At zero carrier, sparse spatial orthogonality
applies. In particular, a negative Nyquist
mode is not paired with itself as though it had a positive-frequency spatial
partner. Opposite modes are combined into real sine/cosine quadratures and use
the original stable temporal Gram kernel, preserving subcycle standing-wave
nulls. This adds O(F²*N) modal contractions without another FFT or retained bank
of frequency fields. Local samples are unchanged; summing them on the original
grid need not equal the continuous cell integral when products exceed that
grid's spatial bandwidth. Exactness here concerns the represented periodic
field, not an isolated aperture or unrepresented source bandwidth. They do not
assign a monochromatic reactive flux or infer dielectric heating from exterior
field moments alone. The selected-layer reservoir owner below supplies explicit
settled-harmonic heat and storage rates. Initial-value dispersive transients,
source switch-on, stochastic coherence inference and nonlinear Maxwell dynamics
remain outside this prescribed linear spectral-line model. Source-bound G1/G3 regressions compare
against independent real-time quadrature after a Debye film, and cover beats,
subcycle shutters, polarization, signed reflection, exact cancellation and
small residuals among large equal-frequency fields. G4 checks source identity,
budget boundaries, cancellation and replay. G1/G3 cell-power regressions cover
x/y Nyquist traveling waves, mixed real/proper ensembles, spatial refinement
of the same physical wave and a short-window opposite-mode interference null.

### Prepared broadband interference between analyzed observation arms

The original `AngularSpectrumExposure` and `AngularSpectrumExposureEnsemble`
owners expose `prepare_interferometer` for 1–4096 exterior observation arms,
each with an explicit Cartesian analyzer and one of the original six
illumination sides. Frequency-grouped Maxwell fields and the original
`HarmonicWindow` kernels are prepared once for the authored shutter. Distinct
points on a common side/depth share modal reconstruction; direct spatial
synthesis reuses the monochromatic coherence owner's DD envelope/carrier path.

Readings apply an explicit complex amplitude gain and time delay to each arm.
With the original positive-time convention the delay multiplies each actual
frequency by `exp(-i*omega*delay)`. Delays are relative source-clock offsets;
no optical path, splitter transmission or dispersion is inferred. Equal-
frequency arm fields combine before products with their DD low components
retained. The deterministic mean and real covariance factors retain carrier
sum and difference tones; proper factors retain only their original difference
tones. Independent factors never acquire cross products from frequency equality.
Uniform and triangular shutters keep the shared temporal owner's normalization.

Preparation admits conservative A*N*sum(L_r) mode work, triangular temporal
pair counts sum(F_r*(F_r+1)/2), retained amplitudes/kernels, one modal plane and
reusable frequency scratch before allocation. Repeated gain/delay readings
allocate nothing and cost O(A*sum(F_r)+sum(F_r²)), without FFTs or material
solves. Geometry, gains, delays, finite arithmetic, budgets and cancellation
remain checked for dark banks. These are electric mean squares in V²/m²,
not detector watts, source PSD inference or detector back-action. The shared
DD temporal arithmetic gives no relative guarantee at an arbitrary beat null.

Focused regressions compare independent analytic multi-arm delay fringes under
both shutters, mixed real/proper source moments, original two-face fields,
sub-ulp modal residuals at a dark port, exact budgets and transactional refusals.

### Prepared broadband interference through a finite physical detector

`AngularSpectrumExposure::prepare_aperture_interferometer` and its ensemble
counterpart cache the original six-component E/H fields for 1–4096 parallel
arms in one explicitly selected observation half-space. Each arm supplies an
outward face distance and transverse translation; it represents
`F(x+dx, y+dy, depth, t-delay)` in the same medium and Cartesian orientation.
Every read applies one authored complex amplitude gain and actual-frequency
delay to E and H together. No splitter geometry, impedance conversion, passive
throughput, ray rotation or cross-factor coherence is inferred.

The returned existing aperture-moment owner retains all 21 real E/H second
moments, vector Poynting flux and signed Cartesian z-directed watts through the
rectangle. Counterpropagating and evanescent fields retain their original
Maxwell relation; electric mean square is never substituted for transported
power. Deterministic means and real source factors retain carrier sums and
differences, while proper factors retain only differences and their original
one-half normalization. Each independent bank recombines all arms before
products and adds its moments afterwards.

Preparation uses the original modal reconstructor, aperture kernels and shared
harmonic-window owner. Its work/storage limits admit original mode-line work,
six-component cached source entries, separable Gram products and peak logical
workspace before allocation. Repeated reads allocate nothing, perform no FFT
or material reconstruction, and reuse a bounded workspace. For N=Nx*Ny and
bank frequency counts F_r, the Gram application costs
`sum(N*F_r*(Nx+Ny+F_r))`, rather than forming a dense mode-frequency pair
table. The real parity calculation is shared with the existing aperture energy
control; proper products use the corresponding conjugate complex Gram.
Compensated arm sums preserve weak residuals before spatial/time integration.
Cancellation leaves the preparation reusable.

G1/G3 regressions in `tests/plane_wave/exposure/aperture_interferometer.rs`
exercise independent analytic 2D, two-frequency E/H moments, both shutters,
mixed real/proper factors, signed backward flux, standing-wave cancellation,
sub-ulp dark-port fields, all preparation budgets and cancellable reuse. This
remains an ideal observation of the source fields; detector absorption,
back-action and an optical network implementing the authored arms require
their own physical owners.

### Finite exposures of proper complex spectral covariance

`angular_spectrum::AngularSpectrumExposureEnsemble` composes the existing
frequency-resolved exposure banks into an explicitly authored partially
coherent source. One bank is one covariance factor across its frequencies,
with square-root statistical weights included in its field amplitudes. For
factor coefficients `z_r`, the source declares `E[z_r]=0`,
`E[z_r z_s*]=delta_rs` and `E[z_r z_s]=0`. Its covariance is therefore
`W(omega,x;omega',x') = sum_r B_r(omega,x) B_r(omega',x')*`, positive
semidefinite without an inferred eigenbasis or sampled random phases. Factors
need not be orthogonal. Properness alone asserts neither Gaussianity nor
stationarity; nonzero cross-frequency covariance can retain time-dependent
beats in the ensemble averages.

The shared proper-product temporal kernel evaluates one half of the real
cross-frequency covariance integrated over the actual rectangular shutter.
Sum-frequency covariance terms vanish because the authored pseudocovariance
is zero, even for subcycle shutters. Each factor retains the existing
compensated equal-frequency field sums, original frequency-specific material
and film operators, complete vector Maxwell reconstruction and signed flux.
An optional deterministic mean bank uses the original sum-and-difference
kernel; its second moments add without a mean/fluctuation cross term.
All banks retain the same original material owners, geometry and grid.

Preparation visits R prepared factor banks without copying their sources or
frequency groups. Queries stream one bank at a time and retain DD arithmetic
through pair contraction and ensemble accumulation. For N spatial samples,
L_r authored lines and F_r distinct frequencies per bank, work is
`O(N + R + sum_r(L_r*N + F_r*N*log N + F_r²*N))`, including the optional
mean as another bank. Additional query storage is
`O(N*max_r F_r + max_r F_r² + N)`; complete borrowed owners are not summed
into that bound. Public limits cap factor count, aggregate line/pair/FFT work
and conservative peak logical numerical storage before a query. Empty factors
and no mean are darkness. Cancellation publishes no partial result.

Results are the existing real mean-square E/H and signed mean Poynting
observables on exterior planes, including integrated normal watts. They are
ensemble moments, not an individual noisy detector realization, a PSD, an
inferred coherence length or transient dispersive heating. Nonzero
pseudocovariance requires a different explicitly authored source model and is
not silently discarded. Focused regressions cover direct four-state ensemble
and real-time quadrature through a Debye film, monochromatic agreement with
the existing spatial coherence owner, partial beat visibility, optional mean
illumination, unitary factor invariance, weak residuals near a distinct-frequency
null, source identity, budget boundaries, cancellation and replay. Extended
precision accumulation does not certify arbitrary relative accuracy near nulls.

### Finite-window interference-dependent film heat and storage

`AngularSpectrumExposure::prepare_layer` and
`AngularSpectrumExposureEnsemble::prepare_layer` return one
`PreparedAngularSpectrumExposureLayer` for a selected original film. One TE/TM
interior response cache per distinct frequency needing a one-sided response is
shared across the deterministic mean and proper covariance factors. Two-sided
lines reuse their already prepared native two-port columns. Cache reuse requires
the same original material owners, layer order, thicknesses, frequency and grid. Compensated
equal-frequency modal sums precede each six-component FFT. Repeated depth
queries perform no new film solve or constitutive evaluation.

The selected law supplies its original instantaneous and dissipative reservoirs:
lossless, Ohmic, Debye, or Debye with Ohmic conduction. The high-frequency
permittivity is retained before capacitor geometry reduction. For each original
Debye pole, `tau*P_dot + P = Delta_epsilon*E` gives irreversible heat
`(tau/Delta_epsilon)*|P_dot|²` and reversible storage
`|P|²/(2*Delta_epsilon)`. Every spectral line keeps that pole's complex
transfer phase. Mixed-frequency interference is contracted within each
reservoir before its nonnegative heat is added to the others. Ohmic heat is
`sigma*|E_real|²`. Unrelated reservoir currents are not added before squaring.

The shared real or proper temporal kernels average these products over the
authored window. The deterministic mean retains sum-frequency terms; proper
covariance removes them through the explicitly declared zero pseudocovariance.
The proper storage derivative forms the frequency difference in DD before
multiplying pair moments, preserving slow beat changes against a large carrier.
Local mean Ohmic heat, polarization heat, total heat and signed storage rate
are W/m³. Storage includes `epsilon_infinity*|E|²/2 + mu*|H|²/2` and every
original polarization reservoir. Its rate is `(U_end-U_start)/duration`;
multiplication by duration yields storage change in J/m³, not deposited heat.
Transverse heat/storage-rate integrals are W/m of depth. Integrating their sum
through the film balances the difference of signed normal boundary powers in W
in the continuous model. The reported heat/storage depth rates remain original-
grid spatial quadratures and need refinement when their products exceed the
field grid's bandwidth. Interior normal watts use the same continuous modal
cell integral as exterior exposures; no aliased point sum replaces that power.
Positive heating can coexist with net outward power when stored energy decreases.

Preparation caps authored frequency groups, distinct-frequency layer/mode
work, aggregate line work, reservoir pair work, FFT work and logical numerical
storage. At most 4096 authored frequency groups are sorted. With R banks,
G authored groups, F distinct frequencies, N samples, S source layers and P
poles, preparation costs `O(R*S + G*log(1+G) + F*N*S + F*P)`.
A depth query costs
`O(N + R + sum_r(L_r*(N+S) + F_r*N*log(N) + (P+1)*F_r²*N + F_r*log(1+F)))`,
plus sum_real F_r*N*(width+height) for a nonzero common carrier. The layer's
reservoir-pair cap counts N*sum(F_r²)*(3+P) plus that spatial pass once;
proper factors add no spatial pass and poles do not multiply its work.
The S term covers structural source checks when original stack pointers differ.
Additional storage is
`O(F*N + F*P + N*max_r(F_r) + max_r(F_r²) + N + G)`, with conditional
O(width²+height²) cached carrier axes and the admitted DD field representation.
Other complete borrowed owners and allocator overhead are excluded.
Cancellation publishes no partial plane.

These are settled prescribed harmonics at fixed parameters, without a switch-on
solution, inferred initial polarization history, temperature update or thermal
ledger debit. A selected empirical optical n/k law without a causal reservoir
is explicitly refused; the existing monochromatic optical heat path retains
its original semantics. Authored regressions cover real-time reservoir
quadrature, positive heat with net outward energy release, heat-plus-storage
Poynting balance, proper covariance with mean illumination, cache reuse,
large-carrier slow-beat storage, original single-tone absorption, source-law
refusal, budgets, depth and cancellation. The nine new Rust regressions await
native execution. Independent JavaScript references are mathematical checks,
not compilation or execution of this owner. No arbitrary-null relative
accuracy or interval certificate is claimed.

### Finite-window joint illumination from both faces

`AngularSpectrumExposure::try_new_illumination` borrows complete
`PreparedAngularSpectrumIllumination` lines. Every line retains both original
face amplitudes, the same Cartesian tangential axes, exp(+i omega t) clock,
material owners and exact physical stack. Lines may differ in frequency and
spatial/polarization content. Equal-frequency modal E/H fields are coalesced
before one six-component spatial synthesis; temporal products are formed only
after all counterpropagating and spatial interference is present.

The existing deterministic and proper-complex ensemble owners handle these
lines through the same frequency grouping and moment contractions.
An optional one-sided or two-sided mean can coexist with one-sided or two-sided
factors. One factor declares a JOINT amplitude across faces and frequencies;
placing sources on opposite faces does not imply statistical independence.
The previously stated zero mean, identity factor covariance and zero
pseudocovariance remain explicit assumptions, without Gaussianity, stationarity,
source bandwidth or an inferred individual realization.

`illumination_field` selects incoming, outgoing or total fields on either
physical side. Cartesian normal Poynting flux retains its sign; a right-incident
wave flows in the negative z direction. On a two-sided bank, the original
`field` selectors map to incoming-left, outgoing-left, total-left and
outgoing-right, respectively. An outgoing field includes reflection and
opposite-face transmission coherently. For a one-sided bank, incoming-right
is darkness and total-right is its original transmitted field.

The SAME `prepare_layer` and reservoir contraction compute finite-window
Ohmic/Debye heating and signed field/polarization storage. Already prepared
two-port columns are reused without another film-response solve. Only distinct
frequencies requiring a one-sided response allocate shared one-sided layer
caches; mixed mean/factor banks share those caches as before. The native
two-sided owners retain their original O(N*S) response storage and prior
admission, while additional film-preparation work is proportional to
F_one_sided*N*S. Temporal pair/FFT work, frequency-group and numerical-storage
budgets still cover all authored banks. Cancellation publishes no partial
exposure or heat plane.

Both exterior channels inherit the native propagating, lossless-medium
requirement. Original internal films may be lossy, critical or evanescent.
The selected causal material reservoir and positive film thickness remain required
for finite-window heating; no switch-on transient, temperature evolution or
thermal energy debit is inferred. Five authored regressions cover analytic
counterpropagating Maxwell fields at all six exterior selections, an independent
two-port Debye/Ohmic heat-plus-storage balance, joint face covariance with a
one-sided mean, coherent cancellation versus independent factors, original
source identity, resource limits, darkness and cancellation/replay. Rustfmt
parsed the source and test modules; Rust compilation and test execution remain
unverified. Independent two-interface and real-time quadrature gave a
1.77e-11 relative energy-balance residual for the two-frequency fixture.

### Two-sided coherent illumination and layer heating

`PlanarDielectricLayers::prepare_illumination` retains the same source-bound
film operator and prepares simultaneous illumination from its two physical
faces. Inputs are peak tangential electric fields in V/m, with one frequency,
transverse wavevector, polarization and exp(+i omega t) clock. Both exterior
channels must propagate in lossless media. Internal films may still be lossy,
evanescent, exactly critical or opaque; the existing one-sided API retains
its separate admission of lossy and total-internal-reflection exit media.

Two bounded sweeps prepare the complete physical E-field scattering matrix
and each film's positive even/odd loss transfers. Evaluation combines complex
source fields before computing outgoing powers or heat, including phase
cancellation and enhancement. Returned boundary E and power-conjugate H are
in V/m and A/m; port powers and each authored film's heat are in W/m². The
different exterior modal admittances enter transmission and reciprocity.
Layer heating comes directly from dissipative modes, so weak absorption is
retained without subtracting nearly equal incident and outgoing fluxes.
Wholly lossless films have exactly zero heat.

Coherent two-source contractions retain double-double multiplication residuals
until both complex paths have combined. This also serves prepared interior
fields and depth-interval loss modes, so a representable near-node field is
not erased by separately rounded products before its positive heat is formed.
The fixed-size contraction remains allocation-free; transfer coefficients are
still binary64 approximations, not certified exact continuum amplitudes.
G0/G1/G3 tests cover exact dyadic residuals, weak bare-interface dark-port
leakage, interior lossy-film nodes, power-of-two scaling and exact cancellation.

Prepared coherent evaluations allocate nothing and cost O(layers). Explicit
source factors describe a positive second moment C=B B*, with one peak-phasor
time-average factor 1/2; an independent-source helper supplies only diagonal
factors. Statistical evaluations preserve interference within each factor,
sum mutually uncorrelated factors, and obey caller mode/layer-work budgets
and cancellation. These are monochromatic ensemble averages, not spectral
densities, an inferred bandwidth, or a spatial coherence model. The report
does not advance temperatures or debit heat into a thermal owner; a coupled
consumer must book the returned finite-film absorption exactly once.

`PreparedLayerIllumination::prepare_point` and `prepare_point_from_right`
retain the same two-source solution inside an authored positive-thickness film.
A prepared point has complete peak Cartesian E [V/m] and H [A/m] transfers;
repeated coherent or factorized-coherence queries allocate nothing and reuse
those columns. The x/z plane is the incidence plane; z points from the left
stack face toward the right. The right-face coordinate preserves shallow exit
fields even when subtracting a depth from an enormous thickness would round
back to the face. Zero-thickness identities have no interior volume and refuse
point queries; both one-sided face traces of a positive film are admitted.

Local fields use only decaying exponentials from physical film faces and the
analytic critical-wave limit. There is no growing transfer matrix, unstable
forward/backward division at zero normal wavenumber, or second global solve.
The complete longitudinal TM electric field contributes to dielectric heat;
an electric tangential node alone does not imply zero absorption. Positive
local heat is omega*(-Im epsilon)*<E squared> in W/m³, using the original
frequency-resolved Debye constitutive law, not a difference of rounded fluxes.
The report also returns componentwise time-averaged squared fields and complex
Poynting flux E cross conj(H)/2, including active and reactive components.
Explicit coherence factors add these observables after within-factor field
interference; mode/work budgets and cancellation apply. Local and integrated
film heat are two descriptions of the same physical loss and must never be
debited independently into a thermal state.

`PreparedLayerField::absorbed_power_in_depth_interval_w_m2` returns unit-source
heat between two original film depths; `prepare_depth_interval` on the
two-source illumination prepares the corresponding coherent/partial-coherent
transfer. Bounds obey `0 <= start <= end <= thickness`. These W/m² results
integrate homogeneous-film loss analytically through the existing positive
slab modes and endpoint fields, retaining skin-depth heating, standing-wave
interference and longitudinal TM loss without depth quadrature or a flux
subtraction. Preparation and coherent evaluation are allocation-free and O(1);
factor evaluation has explicit work/cancellation bounds. Empty intervals and
lossless cells return exact zero, and a whole-film interval reuses the original
loss transfer. Optical phases retain low parts of endpoint differences as well
as `q*k*d`; the existing rounded-distance attenuation safeguards remain.
The caller converts area/time and deposits this heat once. These are stationary
spectral-line averages, not finite-window heat for a transient reservoir.

Focused regressions compare simultaneous-source interface continuity,
admittance-weighted reciprocity, independent E/H boundary equations and
per-layer Poynting work; phase-controlled thin-film heating; full complex
source covariance; and opaque, critical and evanescent internal limits.
Material frequency qualification and non-planar Maxwell evolution remain
outside this analytic boundary model.

## Source-backed nonlinear electrostatics and differential capacitance

Source-driven Gauss fields may have a single grounded physical conductor;
every electrical component still needs its own voltage reference. General
electrode admission shares the original complete-edge and positive-area checks,
while capacitance/network extraction retains its two-terminal requirement.

Both planar and axisymmetric nonlinear electrodes accept prescribed volume
charge, internal free-charge sheets, outward displacement traces and independent
impressed polarization through `solve_with_sources`. The original linear source
assembler supplies exact P1 moments and geometric boundary admission; a normal
flux cannot also prescribe a voltage-controlled or floating electrode. Each
electrical component independently closes
`Q_electrode + Q_volume + Q_sheet - Q_outward = 0`.
The original material `D(E)`, energy, coenergy and analytic tangent are retained;
impressed P is added only to total displacement, with its `E.P` pairing reported
separately. It does not invent a source storage or heat law. Differential
capacitance holds these sources fixed. Maxwell contours require a source-free
exterior, and shape derivatives with nonzero surface/P sources refuse until a
source transport law is supplied. The common bounded Gauss/Newton iteration
retains original-law backtracking, final residual recomputation and checked-CG
work limits, reusing residual scratch for each linear right-hand side.
Synthetic regressions cover nonlinear interface jumps, depolarization, prescribed
flux signs, mixed electrodes and exact radial source-work moments.

`electrostatic::nonlinear::spatial::SpatialNonlinearElectrostatics` extends
the same Gauss/Newton owner to the original native tetrahedral geometry.
Each cell retains a source-resolved scalar/oriented linear dielectric or the
original isotropic reversible D(|E|) curve. Geometry and gradients are prepared
once; the analytic three-dimensional material tangent is applied matrix-free
through native checked CG with diagonal equilibration and compensated charge
reductions. Its positivity uses the same interval gate as the original weighted
Whitney Hodge, including accepted fields requiring no Newton update.

Complete exterior-face electrodes accept mixed prescribed voltage and total
charge, including one grounded conductor with volume sources, floating metals,
and conductors joining dielectric components. Exact constant-cell charge
moments are rho*Volume/4. Every electrical component retains its own physical
reference and independently closes Gauss balance; physical storage and coenergy
remain separate and their sum is checked against the original charge/potential
pairing. The accepted field exposes original E, D, full dD/dE, per-cell storage,
terminal reactions and independent physical residuals.

`spatial::tangent::SpatialDifferentialCapacitance` under this nonlinear owner
reuses one prepared bias tangent across terminal columns. Floating total
charges and volume sources stay fixed; reciprocity, gauge invariance and
component charge are checked without replacing the nonlinear field by a global
Q=C V model. Source endpoints and genuine constitutive kinks refuse an
unspecified two-sided derivative. The native 3D snapshot has no polarization
history, finite charge transport or inferred open-domain boundary. Six synthetic regression groups
specify rotated nonlinear plates, anisotropic native-Hodge parity with space
charge, trapped middle-conductor charge, independent extreme gauges, a
single-conductor charged tetrahedron and source/work/cancellation bounds.

The nonlinear 3D electrode owner also accepts
`solve_with_sources(SpatialElectrostaticSources, displacement_faces, ...)`.
The existing native source assembler is factored through geometry and shared
with linear 3D electrostatics; it integrates the same exact volume/face moments
and polarization weak loads without preparing a throwaway linear Hodge.
Original interior free-charge sheets and outward total-D traces retain their
physical signs and geometry admission, including conflicts with voltage or
floating metal faces. Per-component admission closes
Q_electrode + Q_volume + Q_sheet - Q_outward independently.

Reported total displacement is D_material(E)+P_impressed; the original
constitutive displacement, energy, coenergy and bias tangent remain separately
available. The signed integral E.P is source work, with an independent check
against its weak nodal pairing. The energy identity uses the original
electrode/free-source pairing minus this polarization term. Differential
capacitance holds all prescribed sources fixed. Nonzero impressed P, sheets
and outward flux refuse shape motion until their transport is defined;
explicitly zero sources preserve the existing motion path. Five additional
synthetic regressions specify nonlinear sheet jumps, open-circuit polarization
under an extreme common gauge, flux signs, combined native tensor/source
parity, physical refusals and deterministic retry.

The nonlinear 3D field's `shape_forces` and `virtual_work` differentiate
that same original mixed electrical functional with voltage controls and
floating total charges held fixed. Each tetrahedron transports its original
total volume free charge, so its exact q/4 source moments remain fixed while
rho changes inversely with volume. Frozen isotropic material laws retain the
independently integrated longitudinal storage and transverse coenergy stresses;
no E.D/2 approximation or cancellation-prone ED-W* subtraction supplies them.
The direct work query contracts native cell gradients without assembling nodal
forces, and local displacement differences preserve exact rigid-translation
invariance. Returned shape forces retain their originating field. Tensor-frame
transport, electrostriction, finite trajectories and a mechanical interpretation
for exterior mesh vertices remain the caller's separate physics. Synthetic
regressions compare analytic nonlinear pressures and re-solved mixed-functional
geometry differences, including conserved space charge, rigid motion,
cancellation and unsupported tensor laws.

`electrostatic::nonlinear::axisymmetric::AxisymmetricNonlinearElectrostatics`
uses this same prepared nonlinear Gauss/Newton and differential-capacitance
owner on the original revolved metric. Constant cell fields integrate against
the exact radial volume; volume-charge loads retain the exact radial P1 moments
shared with linear electrostatics. Physical electrodes accept voltage or total
charge, including floating conductors, and each needs its own positive-area
surface. Fields, storage, charge and capacitance include the authored angular
span. The typed result keeps planar force/history APIs inaccessible.
For frozen isotropic laws without volume charge, meridional `virtual_work`
differentiates the actual discrete functional, including the variation of the
radial Jacobian as well as the field gradient and planar area. Axis vertices
remain on the axis. Tensor-axis transport, moving charge, leakage and finite
mechanical trajectories are not inferred.
Synthetic native tests cover nonlinear disk and coaxial responses, charge
control, differential capacitance, exact source moments and re-solved
voltage/charge shape derivatives; they do not qualify empirical materials.

`electrostatic::nonlinear::axisymmetric::hysteretic::AxisymmetricHystereticElectrostatics`
retains cylindrical polarization behind `dielectric-hysteresis`. The isotropic
vector-play law acts in local `[r,z,theta]` coordinates with zero azimuthal
component. Native endpoint solves and adaptive affine paths share the original
Gauss iteration, electrical-component references, source support and work
budgets. Cell storage, physical switching heat and numerical endpoint loss use
the actual revolved volumes and stay separately accounted. Changing volume
charge pairs voltage with the exact radial P1 charge moments used by Gauss;
an arithmetic potential mean is invalid on a general radial triangle.
Typed endpoint/path candidates retain cylindrical geometry and do not expose
planar mechanical or thermal adapters. The public planar history constructor
continues to reject radial metrics. Source coefficients and geometry stay fixed;
rate effects, leakage, moving charge/geometry and coefficient-changing thermal
work are not inferred. Synthetic regressions cover saturation and reversals,
angular-span scaling, changing space charge on an off-center radial mesh,
independent extreme common biases, control changes and bounded cancellation.

Nonlinear fields and retained polarization histories use one independent
voltage reference per electrical component, including components joined by a
shared conductor. References stay separate from the small physical potentials
through charge control, differential response, changing control modes and
adaptive paths. Volume/electrode Gauss balance is gated independently for each
component; opposite residuals cannot cancel globally. Per-vertex and
per-terminal reference accessors retain every authored bias. The scalar
reference accessors describe only the first voltage-controlled terminal overall.

`electrostatic::nonlinear::PlanarNonlinearElectrostatics` retains the original
planar P1 mesh, physical extrusion depth and material sources. Linear scalar or
explicitly framed tensor cells can coexist with the shared reversible radial
`D(E)` material law. The voltage functional integrates `D.dE`; reported physical
storage integrates `E.dD`. Their sum, rather than twice either quantity, is
checked against the actual volume-charge and electrode work pairing.

Prepared complete exterior-edge electrodes accept voltage or total free charge.
A charged floating conductor has one unknown potential and its own total-charge
equation. Shared traces can join disconnected mesh pieces; each resulting
component must have an actual voltage reference. Relative potentials preserve
representable small fields under a large common voltage bias. Volume charge,
individual floating-charge constraints, free-node residuals and total charge
balance are checked independently before a field is returned.

The analytic radial material tangent distinguishes longitudinal differential
permittivity from transverse secant permittivity. Bounded Newton steps use the
original checked-CG owner, native interval coefficient admission and diagonal
equilibration of the equipotential quotient. Prepared cell geometry applies the
tangent without assembling a global matrix. Line search uses direct material
coenergy increments, including source work, rather than subtracting two large
stored objectives. Material support, cancellation and finite nonlinear/linear
work budgets remain hard refusal boundaries; rejected trials publish no field.

`NonlinearElectricField::differential_capacitance` differentiates that accepted
operating point with volume charge and every floating total charge held fixed.
One prepared tangent is reused across voltage-terminal columns. Recovered
charge, reciprocity and gauge residuals must pass in farads. The result retains
its field and source identity; it is a local `dQ/dV`, not a global `Q=C*V` law.
A genuine source-law kink or exhausted field support refuses a unique two-sided
differential response. These fixed-geometry, memoryless, reversible laws do not
infer dielectric relaxation, hysteresis, thermal evolution or electrostriction.

## Nonlinear electric forces and original discrete virtual work

`NonlinearElectricField::maxwell_wrench` uses the original shared geometric
contour owner with an explicitly frozen isotropic constitutive potential.
Its stress is `E tensor D - w* I`. Parallel and transverse principal stresses
use independently integrated physical energy `w` and coenergy `w*`, retaining
small longitudinal pressure without subtracting two large biased quantities.
The original linear magnetic/electric callers retain their existing dyadic
pressure arithmetic through the same contour geometry and first moments.

Every contour edge requires one identical complete isotropic material source
on both sides, zero original free-charge density, and no electrode vertex.
Actual extrusion depth, winding, frame, origin, work bounds and cancellation
remain with the original owner. The enclosed region may contain charge and
contrasting materials. This gives a net instantaneous native `Wrench`, not an
invented solid-surface traction map or an additional overlapping rho E load.

`shape_forces` differentiates the original accepted mixed electrical P1
functional with respect to vertex positions and physical extrusion thickness.
Voltage controls and floating total charges stay fixed. Each cell transports
its original total free charge, so its density changes inversely with volume;
the consistent q/3 source remains unchanged. Isotropic material contexts stay
fixed. Anisotropy refuses until its deformation/frame convention is supplied.
The result retains the field and its actual controls. `virtual_work` contracts
the same derivative directly over cells without allocating a nodal force
vector; local displacement differences preserve rigid-translation invariance.
Its virtual displacement is in metres and its first-order result is in joules.

These envelope derivatives require no extra field solve or tangent assembly,
but inherit the accepted primal residual and spatial error. Independent
fixed-voltage/coenergy and fixed-charge/energy perturbations test the declared
held variables. A source-free nonlinear radial field tests the pressure term
against a distinct half-pairing falsifier. Refinement and original linear force
parity are implementation evidence, not a force enclosure or measured material
validation. Electrostriction, thermal/material evolution, ambient mechanical
pressure, finite moving-mesh work and a mechanical time step remain with their
corresponding physical owners.

## Full field-derived coupled capacitance in native circuits

`PlanarCapacitanceMatrix::coupled_circuit_at_voltages` transfers the original
grounded Maxwell matrix into a prepared native `CapacitanceBlock`. Each
nonreference electrode supplies one charge coordinate and one branch toward the
explicit reference; its charge is the electrode's total free charge. The full
matrix owns `q=Cv` and `H=q^T C^-1 q/2`. Either off-diagonal sign is supported when
the actual grounded block is positive definite. No mutual coefficient is clipped
or fitted, and the existing restricted pair-capacitor export is unchanged.

The export retains the original field/material owner, reference and coordinate
ordering, explicit graph-node map, native prepared matrix and initial charges.
It independently compares all native electrode charges, including reconstructed
reference charge, with original weak electrode reactions, and native H with the
superposed original cell-field energy at the actual voltages. The existing
mesh/electrode connectivity owner refuses unreferenced components before
factorization; an equipotential electrode joining mesh pieces counts as a real
connection. Native preparation owns SPD/numerical admission. Cancellation returns
no export and leaves the original field reusable.

Floating electrodes remain native nodes and charge coordinates. With no external
current their KCL conserves total charge while mutual coupling changes voltage;
charged floating energy is retained. This API performs no Schur elimination that
would discard affine charge offsets or a trapped-energy constant. Append the
returned branches together and use `block_at_branch_offset` to reindex the SAME
shared prepared block without another factorization. Those branches and block
jointly represent one dielectric store; branch diagonal values alone do not
represent its coupled energy. Material and geometry remain fixed; thermal or
mechanical coefficient changes require their existing explicit work owners.

`tests/electrostatic/network/circuit.rs` checks a positive-mutual field matrix
through actual native RC dynamics at microfarad and picofarad scales, independent
terminal equations and recovered field energy/heat, reference and terminal
permutation, common-voltage invariance, charged floating history, append offsets,
physical source identity, additional-gauge refusal and cancellation replay.

### Finite-window exterior Maxwell stress and surface traction

The existing deterministic and proper-complex exposure owners now expose
`stress_field` and `illumination_stress_field`. They contract the complete
symmetric tensor
`T_ij=epsilon<E_i E_j>+mu<H_i H_j>-delta_ij*(epsilon<E²>+mu<H²>)/2`
from the original real Maxwell fields. Off-diagonal components retain
polarization, oblique shear and spatial interference; counterpropagating fields
are combined before any product. Deterministic means retain carrier and beat
terms, while proper factors remove sum terms only through their authored zero
pseudocovariance. Every average uses the supplied normalized shutter.

Each observed exterior supplies its original homogeneous, lossless,
nondispersive material and its own epsilon and mu. The left and right
materials may differ, including distinct electric and magnetic responses.
Left-side observations remain valid when the exit is lossy or dispersive;
right-side observations then refuse. A zero loss at one frequency is not
treated as an instantaneous constitutive law for a finite window.
This electromagnetic stress excludes mechanical medium pressure,
electrostriction and magnetostriction. Results expose local pascals, traction
on validated unit normals, and continuous integrals over the original
transverse periodic cell in newtons. Signed Fourier modes and exact spatial
products integrate the represented quadratic fields before sampling: sparse
orthogonality applies at zero carrier and to proper factors, while real carrier
observations use the shared separable Gram pass. An unpaired negative Nyquist
mode is not its own opposite.
Integrated stress can therefore differ from summing the returned local
samples. This does not add an isolated finite-aperture edge model or certify
physical spatial convergence.

`boundary_stress_resultant_n` evaluates both TOTAL exterior fields with
outward normals, retaining one plane at a time. The signed resultant is not
automatically a finite-window mechanical force. Even an empty vacuum region
can have a nonzero short-window resultant as its field momentum changes.
For a uniform vacuum shutter, mechanical force additionally requires
`-(P_field,end-P_field,start)/duration`; a shaped shutter requires the
corresponding weighted momentum rate. Material-filled interiors need their
own consistent field/material momentum law. None is inferred or set to zero.

Original exposure admission bounds N*F² temporal frequency-pair units plus
the conditional separable carrier work described above. Each stress
query performs twelve local and at most twenty-four additional modal scalar
temporal contractions per admitted unit, at most thirty-six total. DD
arithmetic is retained through factor accumulation, continuous-cell
integration, trace subtraction and subtraction of boundary resultants.
Equal-frequency grouping and synthesis reuse the original exposure owner.
One bank's modal fields are integrated and then synthesized once in place:
no extra FFT, constitutive solve or film solve is added. A nonzero real carrier
also retains its admitted centered and contracted DD coefficient planes.
Work is O(sum_r(L_r*N + F_r*N*log N + F_r²*N)), plus the admitted
sum_real F_r*N*(width+height) carrier pass, with O(N) output/accumulation.
Original work caps remain authoritative. The explicit
stress-workspace cap includes two spatial-quadrature work arrays and three
six-component modal accumulators and admits the complete query before
allocation. The two-face resultant performs two separately admitted plane
queries and therefore uses twice the per-plane work, retaining only one
plane at a time. A general normal contracts the continuous transverse-cell
tensor; it does not rotate or resample that plane. Cancellation publishes
no partial field or resultant.

Seven original G1/G3/G4 regressions cover full tensor real-time quadrature,
proper ensembles, Cartesian covariance, ultrashort carrier phase, joint-face
standing waves, source-bound shear/traction, momentum-storage distinction,
source refusal, budgets, darkness, cancellation and replay. Before the
continuous-cell extension, the temporal/stress contraction and three test
definitions passed native no-codegen type analysis; that earlier result
does not establish compilation of the current extension.

Three additional source-bound regressions cover x/y Nyquist alias
counterexamples under uniform/triangular shutters, deterministic means and
real/proper factors, true opposite-mode interference, and original
work/workspace admission, cancellation and retry. Independent analytic and
dense JavaScript quadrature references were evaluated. In the unit-medium
Nyquist stress fixture, sampled quadrature overestimates the continuous
tensor by factors 1.9352 and 1.9643 for uniform and triangular shutters.
The added Rust code and regressions remain uncompiled and unexecuted because
no native execution environment is available. Independent mathematical
references do not constitute execution of the Rust owner.

### Real and nonproper spectral source fluctuations

`AngularSpectrumExposureEnsemble::try_new_with_statistics` accepts an
explicit `AngularSpectrumFactorStatistics::{Real, ProperComplex}` for each
original factor bank. The existing `try_new` remains an all-proper declaration.
Statistics are immutable borrowed inputs, with exactly one entry per factor,
including empty banks. The accessor `factor_statistics` reports that declaration.

Every factor has zero mean and unit second moment. A real coefficient has
E[z²]=1 and retains phase-locked carrier sum-frequency terms; a proper complex
coefficient has E[z²]=0 and retains only difference-frequency terms.
Real-real covariance and proper-proper conjugate covariance are identity;
cross-factor unconjugated moments and mixed real/proper moments vanish.
No Gaussian law or statistical independence is required. Unknown source
phase does not justify choosing either law.

Complex fields driven by real factors can represent general nonproper
second moments: a real square root of the stacked real/imaginary source
covariance supplies the factor columns, including their square-root weights.
The deterministic mean contributes separately. Cross-frequency and cross-face
coherence remains inside each factor. No source spectrum, PSD normalization,
sampling bandwidth, stationarity or physical source calibration is inferred.

One shared bank traversal selects the temporal law once per factor across
exterior E/H/Poynting observations, two-sided illumination, original-layer
Ohmic/Debye heat and reversible storage, and full Maxwell stress/traction.
Real factors use the deterministic REAL-PRODUCT kernel while remaining
zero-mean fluctuations. Proper factors use the difference-only product and
storage-frequency law. The distinction is therefore retained through the
entire heat/storage/stress observation, not just through an exterior intensity.

No field, covariance matrix or constitutive solve is added by the statistics
mapping. Original source identities, frequency caches, scalar work limits,
workspace limits and cancellation paths are retained. Queries still stream
one factor at a time and have no statistical-mode branch in their spatial or
frequency-pair contraction loops. Uniform and triangular detector windows
remain supported. Storage rate is a weighted average; its endpoint-change
interpretation applies only to a uniform window.

Four authored native-owner regressions cover phase-locked subcycle energy
and stress, real/proper mixed joint-face fields and dispersive heat/storage,
equivalence of a proper factor to two real quadratures, and explicit shape,
identity, budget, darkness, cancellation and replay behavior. Source formatting
passed. The exact new statistics iterator/accessors and shared temporal owner
passed focused native no-codegen type analysis with nominal source-owner
declarations; this is not a complete EM crate build or executed owner tests.
Independent eight-state real/proper ensemble quadrature matched the analytical
18-component field/stress/heat/storage fixture to maximum absolute errors
5.0e-15 (uniform) and 1.0e-14 (triangular). These are mathematical reference
checks, not execution of the Rust propagation owner.

### Beat-resolved local Ohmic and Debye heat spectra

The original `PreparedAngularSpectrumExposureLayer::heat_spectrum` prepares
an `AngularSpectrumExposureLayerHeatSpectrum` for one row-major transverse
sample at a depth in its selected original film. The result borrows that
prepared layer and exposes its original grid, sample index, depth and reservoir
selection. `spectrum()` returns the shared analytic product spectrum for
downstream causal time evolution; scalar-time and Uniform/Triangular detector
queries are also available, with cancellation.

`AngularSpectrumExposureHeatChannel` selects Ohmic, one DebyePole by original
material index, Polarization (all original poles), or Total. Ohmic heat is
sigma |E_real|^2. Each Debye contribution is
(tau/Delta_epsilon) |Pdot_real|^2, using the original cached
`PoleTransfer::heat`. No new constitutive evaluation or film solve is added.
Reversible epsilon-infinity, magnetic and polarization storage rates are
excluded even when net boundary electrical power is negative. The output is
local W/m^3, requiring an explicit physical volume/quadrature to become watts;
it must not be debited again with the same integrated film absorption.

Every actual frequency plane is reconstructed once through the existing
coherent grouping and internal-field owner; only the requested point is
retained. Original real/proper bank traversal and joint moments remain in use.
Reservoir filtering runs once per retained line before lazy Cartesian
frequency-pair products enter the shared fs-fft engine. Original EM phasors use
exp(+i omega t): the complete filtered response is conjugated into the shared
exp(-i omega t) convention. Native angular frequencies remain original DD
sum/difference frequencies rather than rounded f64 Hz carriers. DD filtering,
contraction and tone accumulation preserve additional arithmetic precision;
they add no accuracy to the original field reconstruction.

Square-root Ohmic scaling keeps both operands balanced, with a DD correction
retaining the authored sigma rather than its rounded square root squared.
The three Cartesian squared components and separate irreversible reservoirs
add; no unphysical cross-pole or cross-component square is introduced.
Proper factors retain difference tones only; phase-locked real factors and
deterministic means retain carrier sum tones. These are ensemble-expected heat
waveforms under the original authored joint second moments, without inferred
Gaussian statistics, independence, stationarity or detector bandwidth.

`AngularSpectrumExposureHeatSpectrumLimits` admits the actual scalar product
count, generated tone count and conservative peak logical workspace before
new allocations. Workspace includes the original layer's complete bound plus
retained point/reservoir fields and shared unsimplified tones. The shared hard
tone ceiling and original source reconstruction work limits still apply.
Additional work is original reconstruction plus O(G*R + K log K), with
O(G*R + K) additional storage for G bank-frequency groups, R selected reservoirs
and K generated tones. Queries allocate nothing and retain no partial result
after cancellation or refusal. Invalid depths, points and original pole
indices refuse. Negative rounded heat observations refuse instead of being
clipped. No relative error certificate at arbitrary destructive nulls is given.

Two focused G1/G3 tests compare two-pole, three-component mixed-bank heat with
independent original-ODE currents over the exact eight real/proper states and
both detector profiles. They also cover proper-as-real quadrature invariance,
an original high-carrier angular beat and resource/cancellation refusal.
Three additional native-owner integration tests compare actual internal film
fields and original finite-window dissipation, including joint illumination
faces, positive heat during negative reversible storage rate, coherent weak
residuals, original source/position identity and atomic refusal/replay.
No switch-on polarization state, temperature-dependent material feedback,
bulk spatial spectrum API or full thermal evolution is inferred by this owner.

The two focused reservoir/spectrum regressions compiled, linked and executed
natively with the real fs_math library: 2 passed, 0 failed. Exact new heat
code, original PoleTransfer and shared spectrum code were used; upstream
Maxwell/material admission was nominal and was not executed. The three
actual film-owner integration tests are authored and formatted but were not
executed or fully typechecked with the complete EM dependency graph here.
This evidence is scoped to reservoir/filter/product mathematics, not a full
Maxwell solve, temperature-feedback integration or repository quality gate.

### Physically integrated film heat-load spectra

`PreparedAngularSpectrumExposureLayer::integrated_heat_spectrum` integrates
actual coherent internal fields into watts. The returned
`AngularSpectrumExposureLayerHeatLoadSpectrum` borrows the original prepared
layer and the complete authored volume quadrature. It exposes the selected
Ohmic/Debye channel, original grid, declared support, actual DD sums of cell
areas/depth weights and represented volume, together with the shared analytic
spectrum. This spectrum is ALREADY in W; a thermal consumer must not multiply
it by volume again or debit overlapping local/integrated absorption.

`AngularSpectrumHeatVolumeQuadrature` names an explicit depth span within the
original film and positive `AngularSpectrumHeatDepthSample` weights in metres.
Every node lies in the declared span. The weight sum must match its length
within 32*f64::EPSILON relative allowance. Actual authored weights are retained
and used unchanged; no normalization forces agreement.
`AngularSpectrumHeatAreaQuadrature::PeriodicCell` uses original dx*dy for each
row-major field sample. `CellAreasM2` supplies one finite area in [0,dx*dy] per original
cell; zero masks a cell and positive fractions describe partial cells without
relocating field samples. An all-zero mask is an explicitly empty support.
`ContinuousPeriodicCell` instead integrates the represented transverse Fourier
fields continuously across the full cell. Equal signed modes contribute to
difference products; genuinely opposite modes contribute to sum products.
The unpaired negative Nyquist mode is never its own opposite. Real carrier
terms remain for the zero mode and true counterpropagating pairs; only authored
proper factors remove all sum terms. Original source coherence, pole phases,
depth weights and native angular clocks remain unchanged.

For every source bank, each original frequency plane is reconstructed once
per depth through the existing coherent field owner and retained only for that
depth. The spatial DD contractions are the volume quadratures of
E_f dot conj(E_g) and E_f dot E_g. Continuous full-cell support contracts the
original modal vectors before synthesis, using the physical complete cell area
once per admitted mode pair. Sampled supports synthesize the original fields.
Original cached reservoir-transfer products
are summed once per frequency pair, outside the spatial/depth loops:
sigma plus sum H_rf*conj(H_rg), or sigma plus sum H_rf*H_rg.
This avoids repeating spatial work for every selected Debye pole.
The original PoleTransfer remains the sole constitutive response owner.

Spatial and reservoir real/imaginary accumulators each carry a DD mantissa
and separate binary exponent. Exact power-of-two decomposition keeps tiny
field squares or tiny reservoir-response squares from underflowing before
large balancing physical coefficients are applied; the inverse large-square
cases also avoid premature overflow. Area/depth weights participate in this
same arithmetic. Physical watt amplitudes are restored only after spatial,
reservoir and diagonal scales meet. Nonzero final tones outside the DD
amplitude range explicitly refuse instead of becoming a dark spectrum.
This extends dynamic range, not DD precision or arbitrary-null accuracy.
The enlarged actual pair-cache type is included in workspace admission.

The heat is a symmetric sum of squared real fields. Upper frequency pairs
combine their symmetric counterparts: diagonal terms retain factor 1/2 and
off-diagonal terms receive factor 1. Real factors and deterministic means
emit sum and difference tones; proper factors emit differences only.
Complete positive-time EM coefficients are conjugated into the shared
negative-time convention, retaining DD native angular frequencies.
Spatial/depth contraction precedes emission: adding samples does not multiply
the number of retained tones. Reversible electromagnetic/polarization storage,
boundary momentum work and parameter work are excluded from deposited heat.

`AngularSpectrumIntegratedHeatLimits` separately admits depth sample count,
repeated original point/source-line and FFT work, spatial point/frequency-pair
work, reservoir/frequency-pair work, generated tones and conservative peak
workspace. The shared hard term ceiling still applies. With Z depth samples,
N original grid points, F_r frequencies per bank, R selected reservoirs and K
generated terms, added work is O(Z*N*sum F_r^2 + R*sum F_r^2 + K log K), besides
Z original field reconstructions. Additional retained storage is
O(max_r(N*F_r + F_r^2) + K), including both contracted raw terms and shared
tone preparation. The conservative bound also includes the original complete
prepared-layer workspace. Continuous full-cell support needs zero additional
synthesis FFT work and accepts `maximum_fft_work=0` for this query; its spatial
pair budget counts modes instead of samples. Original beam analysis and layer
preparation keep their own work bounds. Borrowed inputs and allocator overhead are excluded.
Cancellation/refusal publishes no partial load and leaves the original owner
reusable. Scalar/window observations allocate nothing; negative rounded heat
refuses rather than being clipped.

This retains explicit physical-volume support. Continuous full-cell integration
removes transverse product aliasing for the represented Fourier field; it does
not certify that the physical field is spatially converged. Sampled area options
can still alias field PRODUCTS even when fields are resolved. All options use
numerical depth quadrature, whose samples must resolve standing waves and
attenuation. The caller owns support and refinement.
No quadrature error enclosure, isolated-aperture convergence, switch-on state
or material-temperature feedback is inferred.

Five focused G1/G3/G4 regressions compiled, linked and executed natively with
real fs_math: 5 passed, 0 failed. They exercise exact new spatial contraction
and volume admission, original PoleTransfer and shared spectral preparation.
Independent positive-time original-ODE currents and the exact eight mixed
real/proper states agree across two poles, vector fields, partial/full cells
and both shutter profiles. Tests verify one frequency-plane callback per
bank/frequency/depth independent of reservoir count, fixed spectral term
count when depth sampling doubles, preserved authored weights, invalid support,
cancellation and malformed-plane refusal. Extreme-scale regressions retain
finite heat for sigma=1e200 with E=1e-200, tiny original Debye responses with
large fields, large original Debye responses with tiny fields, and large
fields integrated over tiny physical areas, for real and proper laws.
Independent references filter or balance physical operands before their
products, matching the existing local-owner convention. Truly unrepresentable
nonzero final watt tones refuse; genuinely zero fields remain dark.
Upstream Maxwell propagation and material source admission are explicit nominal boundaries in this focused
executable; those boundary methods are not executed.

Two additional original-owner integration tests are authored and formatted:
a partial-cell mask retains a proper-source beat that cancels over the full
periodic cell, with original finite-window heat comparisons; depth refinement
converges to original monochromatic absorbed watts, with exact budget and
cancellation/replay checks. Those full-owner tests were not executed or fully
typechecked against the complete EM dependency graph here. An independent
normal-incidence film reference confirms the latter fixture's 19.7507815275 W
absorption and Simpson errors 2.89e-3, 1.02e-5 and 6.34e-7 W for 8,32,64
intervals; this separate reference is not execution of the Rust film owner.

Continuous-cell regressions add signed Nyquist falsifiers, true standing-wave
and zero-mode controls, mixed real/proper two-frequency original-reservoir
comparisons against refined spatial integration, and zero-synthesis-budget
admission. The native thermal feedback test uses continuous support while
checking enthalpy deposition and its late DD clock. Three bounded native tests
execute the verbatim spatial-pair kernel, signed-mode mapping and scaled DD
arithmetic with real fs_math, including refined vector-field products and
extreme balanced scales. These checks exclude material/FFT orchestration;
the new full-owner EM and thermal integration tests are authored, not executed
against the unavailable complete dependency graph in this environment.

## Tetrahedral dielectric and polarization weak forms

`electrostatic::tetrahedral::TetDielectricOperator` prepares the full P1
dielectric matrix on the original `TetComplex`: `K = d0^T M_epsilon d0`.
It uses checked Whitney geometry, the native weighted-Hodge interval SPD
admission and native deterministic sparse products. Every cell has an explicit
absolute symmetric permittivity tensor in the mesh frame. These numerical
inputs do not themselves establish measured-material authority; a coupled
material owner retains its source, frame transform and applicability domain.

The same cell gradients and positive reference volumes compute polarization
load `b_P = sum V grad(N).P`, electric field `E = -grad(phi)` and displacement
`D = epsilon E + P`. The nodal equation is `K phi = q_free + b_P`, and the
free-charge reaction is `K phi - b_P`. Voltage differences are formed before
field multiplication, preserving constant-potential gauge invariance.
Coefficients and reconstructed fields remain separate across material jumps.
Preparation is reused across load cases; field/load evaluations take linear
work and retain no evolving state.

Node/cell and geometric-quality bounds are admitted before assembly. Invalid
connectivity, isolated nodes, degenerate cells, nonfinite arithmetic and
unproved SPD weights refuse. Cancellation polls surround native batched
geometry/sparse products and occur during cell admission, weighted assembly
and field recovery; no interruption inside those native batch calls is claimed.
No partial operator or response is returned on refusal. Electrode conditions,
mechanical coupling, convergence, source authority and actual field solving
belong to consumers. This primitive adds no conduction, nonlinear dielectric,
dispersion, moving geometry, mesh-overlap certification or continuum error bound.
Focused synthetic G1/G3 tests are in `tests/tetrahedral_electrostatic.rs`.

### Finite-window exterior field momentum and storage rate

The original exposure and real/proper ensemble owners expose momentum_field
and illumination_momentum_field beside their Maxwell-stress queries. They
return the explicitly named Minkowski density <D cross B> in N s/m³ and
<d(D cross B)/dt> in N/m³, using the same original homogeneous lossless
observation-side exterior epsilon and mu. Each half-space retains its own
original lossless nondispersive law; the two materials need not be identical.
Lossy/dispersive exit laws refuse right-side queries while leaving valid
left-side observations available. In vacuum this is the electromagnetic
momentum S/c². A material exterior does not thereby acquire an inferred mechanical,
striction or moving-body momentum law.

The complete vector fields retain both physical faces, polarization, spatial
interference and equal-frequency coherent grouping. Deterministic real
products keep carrier and beat terms; only the authored proper factor law
removes sum terms. Real derivatives separate carrier and difference terms
before multiplying their frequencies. Proper derivatives form the frequency
difference in DD first, so a large common carrier cannot replace a slow beat.
A proper monochromatic second moment has exactly zero time derivative.
Every rate uses the original normalized shutter; only uniform weighting
reduces to an endpoint momentum difference divided by duration.

Results retain original point samples plus compensated exact continuous-cell
momentum and rate integrals in N s/m and N/m. The native signed modal pairing
and cosine/sine quadratures integrate the represented band-limited fields
before their quadratic products are sampled, so a Nyquist carrier cannot
alias into a false spatially uniform momentum or storage rate. DD modal
sums and differences remain DD through temporal contraction and frequency
differentiation. Spatial quadratures combine before frequency factors; the
already conjugated modal parts rotate by -i while original point amplitudes
rotate by +i. A consumer still owns the enclosing depth/volume integration
and interior material accounting. One plane is not a body force, and this
query does not advance mechanics or debit energy.

The original fields/temporal kernels are reused without new film solves.
One real/proper bank is streamed at a time; its exact modal integral precedes
in-place synthesis of the original point samples, without a second retained
bank or additional FFT at zero carrier. Nonzero real carriers retain the
admitted DD cell representation. The original N*F² temporal units plus the
conditional F*N*(width+height) spatial units remain authoritative.
The fixed modal pass contracts Q <= 2*N spatial quadratures as well as N point
samples: each real scalar-pair vector density/rate uses eighteen products,
and each proper pair uses at most twelve. The additional complete-workspace
limit includes the frequency coefficient cache, DD point accumulation,
returned points and fixed modal scratch. Bounded cancellation applies before
publication.

Two new contraction regressions cover independent real-time density/rate
quadrature, proper ensembles, uniform/triangular shutters and large-carrier
slow beats. Three original source-backed stress fixtures additionally cover
joint-face momentum, its relation to Poynting flux, vacuum volume/surface
momentum balance, budgets, cancellation, darkness and replay. Two further
source-backed regressions exercise actual x/y Nyquist waves, real/proper
factors and a signed cross-frequency beat. They distinguish aliased point
sums from continuous-cell density and rate: a single traveling spatial mode
has zero integrated rate while its original point rates remain nonzero.

Independent JavaScript mathematical checks agreed to 4.22e-15 absolute on
the unit-scale vector/time fixture and 3.49e-9 relative on the 128-interval
vacuum depth quadrature. Additional modal checks agreed to 7.11e-15 absolute
with the analytical Nyquist two-frequency beat and 3.93e-12 with direct
space/time quadrature, including paired modes, a constant and Nyquist bins.
In the one-mode fixture, sampled density was 1.974 times its exact cell
integral and sampled rate was -37.15 N/m although the exact integrated rate
was zero. The Rust regressions are authored but unexecuted: no native
execution environment is available. These checks do not certify the full
Maxwell owner, a material force law or a repository-wide quality gate.


### Finite-window momentum across distinct exterior materials

The common exposure-side selector binds every stress and Minkowski query to
the actual observed material, including all incoming/outgoing/total variants
and the real/proper ensemble owners. Source identity across a frequency bank
is still enforced by original exposure preparation. No constitutive source,
permittivity or permeability is reconstructed from a rounded refractive index.

Three source-bound regression groups exercise electrically and magnetically
distinct interfaces with independent Fresnel/tensor/momentum formulas,
both-face real/proper illumination and observation-side refusal for Debye
and ohmic exits. Their finite-window rate references differentiate the actual
real fields. Right-side refusals leave the valid left observations unchanged.
The changed Rust was parsed with the pinned rustfmt, and the new regression
module was formatted; independent analytic checks were evaluated. Native compilation and tests remain unexecuted
because the shared build disk is full and memory is under pressure. This extends
field observations, not a material mechanical-force or striction law.

### Joint finite-window stress and momentum observation

The original exposure and real/proper ensemble owners expose
`stress_momentum_field` and `illumination_stress_momentum_field`.
One query returns the existing stress and momentum result owners together:
local Maxwell stress/traction, explicitly named Minkowski D cross B density,
its shutter-weighted instantaneous time derivative, and their continuous
periodic-cell integrals. Each observed exterior retains its own original
scalar nondispersive lossless material. The original field, surface/volume
integration, material-momentum and striction no-claims still apply.

A private native frequency-bank owner now serves the separate and joint
queries. The joint path reconstructs each coherent frequency field once,
prepares one temporal kernel cache and synthesizes each retained field once.
Stress projects the original product from the momentum kernel cache.
Both exact modal contractions run before in-place synthesis; original
point samples follow. The original dyad and momentum contractions and
publication/scaling functions remain shared numerical owners, preserving
real/proper factors, coherent sum/beat terms, both-face interference, DD
frequency differentiation and signed Nyquist integration.

The joint workspace cap includes both DD accumulation planes and outputs,
one original reconstruction/FFT workspace, the momentum frequency coefficients
and conservative modal scratch. It equals the two individual workspace
bounds less one original workspace bound. The nonzero-carrier DD cell
representation is shared by both observables. Preparation's N*F² temporal
pair work and conditional F*N*(width+height) spatial work remain authoritative;
at most ninety real or seventy-two proper
scalar temporal products occur per admitted pair unit across point and
modal observations. Work is O(sum_r(L_r*N + F_r*N*log N + F_r²*N)), plus
conditional separable carrier work once per real bank. Avoiding a duplicate
field reconstruction/cache/FFT is a source-level work
reduction, not a measured wall-clock speedup claim. Cancellation or a failure
in either observable returns no partial pair.

Three source-bound G1/G3/G4 regressions cover equivalence to the individual
owners on all physical faces and one-sided entry points, mixed mean/real/proper
statistics, unequal exterior media, analytical Nyquist cell stress and
momentum/storage, complete-workspace refusal, mid-query/final cancellation,
darkness and deterministic retry. The new inputs pass pinned-nightly
rustfmt parsing/formatting. The Rust regressions remain unexecuted; no native
typecheck, runtime performance benchmark or broader quality gate is claimed.

## Planar moving-induction loads on mechanical receivers

Both `MovingEddyCandidate::lorentz_loads` and
`NonlinearMovingEddyCandidate::lorentz_loads` resolve the original conducting
triangle selection into nodal force and arbitrary receiving-frame wrenches.
They use the linear midpoint or nonlinear endpoint J cross B stage already
used for native Lorentz work; no field re-solve, source substitution or second
Maxwell-stress force is introduced. Saturation, admitted fixed recoil and
explicit solid-current drive retain their native current and stage.

The same compensated load-scatter and wrench owner as spatial induction now
admits three or four simplex nodes. Extruded planar P1 weights are
V*(1+delta_ij)/12; tetrahedral weights remain V*(1+delta_ij)/20. Uniform
extrusions act at their physical mid-depth, preserving all three torque
components about arbitrary origins. Planar coordinates are mapped lazily, so
no additional full-mesh coordinate array is allocated. Selected triangles are
canonicalized and checked for uniqueness, conductor membership and budgets;
refusal/cancellation publishes no partial load.

The exported `PlanarLorentzLoads` uses original vertex order, with local
x/y forces and zero z force. Pairing with [vx,vy,0] recovers the selected
native mechanical power; `wrench` applies the existing rigid placement and
moment-origin transformation. A mechanical impulse uses the same duration.
The receiving dynamics must close its actual work on the same time stage.
This transfer does not evolve material geometry, close mechanical feedback,
or add magnetization/interface forces. The original fixed Eulerian occupancy
and source assumptions remain in force.

Nine focused native load/mechanics tests pass: six existing tetrahedral
regressions plus three extruded-triangle checks for analytic nonuniform torque,
P1 virtual work, subdivision/shared-node scatter, and actual articulated
linear/angular momentum and kinetic work. Linear and nonlinear field-integration
tests also gain full/subset native-work and selection/cancellation checks;
those two full field test binaries were not executed in this resource-limited
session. No full Cargo, DSR or cross-platform result is claimed.

### Joint finite-aperture and finite-shutter optical observations

`AngularSpectrumExposure` and `AngularSpectrumExposureEnsemble` now expose
`measure_aperture` and `measure_illumination_aperture`. They integrate the
original continuous modal fields over a uniform axis-aligned rectangle and
the authored normalized uniform or triangular shutter. The reading contains
six electric/magnetic mean squares, the full signed vector Poynting flux,
and signed Cartesian z-directed power in watts through the rectangle.
Frequency grouping, source clocks, polarization, propagating/evanescent
fields and both original illumination faces retain their existing owners.

Deterministic means and real factors retain carrier sum and difference
terms. Proper-complex factors retain difference interference with their
one-half normalization applied once. Independent factor banks add moments;
within-bank spatial, temporal and inter-face coherence remains intact.
Exact separable harmonic kernels require no spatial quadrature, time stepping,
new FFT, material solve or dense source covariance. Real-field contractions
use the original shared `HarmonicProduct::centered_gram_dd` weights, retaining
grouped double-double amplitudes and applying the combined spatial/temporal
center phase once per mode and frequency. Start coordinates and half-width
phases remain separate, preserving tiny translated interference nulls without
rounding the detector width into its position.

Checked admission bounds ordered mode/frequency contractions by
`N² * sum(F_bank²)` before query allocation. Spatial pair caches require
`O(width² + height²)` storage, center rotations require `O(N * F_bank)` work,
and banks stream one at a time. The conservative peak workspace includes the
largest original bank query, retained double-double fields and spatial caches;
overflow, finite-value and bounded cancellation checks remain explicit.
The original periodic source continuation still applies across cell edges.
This observation does not infer absorptivity, detector edge diffraction,
back-action or unprovided source statistics.

Seven focused joint-observation tests and five original monochromatic
aperture tests pass natively. They cover independent space/time/phase
quadrature of all nine moments, carrier/beat averaging, proper coherence,
partitions, full-cell/Nyquist behavior, resource admission and cancellation,
and translated weak-null regressions for uniform/triangular shutters with
unequal widths down to `1e-50`. Native probes exercise the production
contraction and shared harmonic/math owners; they do not execute material
or Maxwell reconstruction. Three additional original-source integration
tests cover film/evanescent fields, multi-frequency coherence and mixed
two-face ensembles, but remain unexecuted without the sibling constellation.
No full Cargo, DSR, cross-platform or experimental-validation result is claimed.

## Original finite-band material optics in Maxwell propagation

`plane_wave::SpectralOpticalResponse` consumes the original
`SpectralComplexIndexStatePoint`, the SI source-wavelength conductor wrapper,
or the typed phase-index/intensity-absorption dielectric snapshot through
`SpectralOpticalSource`. Arbitrary canonical complex-index knots remain on
their original grid. The caller supplies a positive maximum knot span in
vacuum metres, and scalar absolute magnetic permeability must belong to the
same complete material snapshot. No nonmagnetic default, visible clamp,
spectral extrapolation or extra Ohmic loss is supplied.

`DielectricResponse::SpectralOptical` reaches the existing native interface,
layer and angular-spectrum beam owners, including transparent dispersive
incidence. Original n/k interpolation delegates to the shared material owner.
For the dielectric family, original n and power-absorption alpha interpolate
before k=alpha*lambda_vac/(4*pi); interpolating a converted k table would change
the declared absorption law. The existing n-i*k convention and
epsilon=(n-i*k)^2/(mu*c^2) drive phase, attenuation, impedance and Fresnel work.
Frequency admission precedes any roundoff correction at source-band endpoints.
Monochromatic evaluation allocates nothing.

Source refresh returns an owned, fully readmitted snapshot while retaining the
same interpolation-gap policy and original claims. Generic complex-index
bindings retain their authored wavelengths; the SI conductor wrapper refreshes
its source wavelength properties as well as coefficients. The same-state
permeability refreshes with them. Failure leaves the original response intact.
Coherent source compatibility retains the original source owner and gap policy.
Finite-band monochromatic tables do not establish causal time-domain heat
reservoirs or a two-sided group derivative at piecewise-linear knots; the
existing reservoir preparation explicitly refuses them.

Four native dispersion-module tests pass against real fs_math, including SI
intensity attenuation, the complex Maxwell dispersion relation and strict
finite-band frequency edges. Four source-owner integration regressions cover
infrared complex Fresnel, same-state permeability refresh, transparent incident
Snell refraction, coherent film/beam phase, source-owner compatibility, original
alpha interpolation, evolving SI wavelength support and retained gap limits.
Those complete material/FFT integration tests remain unexecuted without the
constellation dependency checkouts.

### Polarization correlations through finite optical apertures and shutters

The same exposure/ensemble owners expose `measure_aperture_moments` and
`measure_illumination_aperture_moments`: all 21 unique real uncentered field
second moments, comprising symmetric EE/HH and full EH in the original
Cartesian frame. The antisymmetric EH terms recover signed Poynting flux.
These are not analytic coherency, Stokes parameters or absorbed detector power.

`AngularSpectrumExposureApertureProjection` authors finite dimensionless real
electric/magnetic analyzer weights without implicit normalization. Direct
projection queries apply them to the original grouped double-double fields
before center rotation and quadratic products. This retains weak crossed-
analyzer signals that may vanish in post-hoc moment contractions; no general
relative-error guarantee near every null is implied. Longitudinal/evanescent
components, real/proper statistics, both illumination faces, source clocks,
spatial/shutter kernels and original admission/cancellation remain shared.

The same direct projection queries also accept
`AngularSpectrumExposureApertureJonesProjection`. Its finite dimensionless
complex weights define the real observed signals
`Re(sum_j a† E_j exp(i omega_j t))` and
`Re(sum_j b† H_j exp(i omega_j t))`, with conjugated analyzer weights.
They select circular and elliptical polarization using the original phasors;
real moment matrices alone cannot recover that information. Projection keeps
the grouped double-double low parts through a crossed-circular null, before
any field products. Real sources retain carrier sums over subcycle shutters;
proper factors keep their original phase statistics. These explicitly authored
frequency-independent observations over the supplied spectral lines are not a
causal broadband polarizer, detector absorption model or source back-action.
The real-weight API and post-hoc real moment projection remain available.

The original reading, full moments and direct projection require respectively
12, 21 and 3 scalar products per ordered mode/frequency pair. Direct field
projection costs O(NF), with no new material solve or FFT. Five new numerical
tests pass alongside the original 12 aperture regressions, including independent
space/time/phase quadrature, linear-analyzer scaling, basis rotation and a
crossed-analyzer mean square near `5e-41`. Two additional original-source
integration tests are committed but unrun without the sibling constellation;
the native numerical probes do not execute Maxwell/material reconstruction.

## Source-backed analytic finite-band Maxwell dispersion

`SpectralOpticalSource::DielectricSellmeier` borrows the same original SI
coefficients, explicit fit limits and independent intensity-absorption spectrum
as the material and native-render owners. `SpectralOpticalResponse` admits the
complete fitted interval and evaluates the original Sellmeier equation directly;
only absorption is interpolated before `k = alpha*lambda_vac/(4*pi)`. The caller's
maximum knot span applies to the original absorption samples, not a fabricated
phase-index table. Frequency evaluation stays inside the common source domain.

Scalar absolute permeability must still belong to the same complete state.
Refresh retains exact source selection and rechecks the full domain. Original
interface/layer consumers receive the same monochromatic response; this does
not introduce a causal time-domain reservoir, inferred permeability, or
out-of-band material behavior.

### Continuous rectangular support for original film heat

`AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(AngularSpectrumAperture)`
integrates the represented transverse electric field continuously over the
authored x/y rectangle, before applying the original Ohmic and Debye reservoir
products. It retains distinct-mode interference, real-bank sum and difference
tones, proper-factor difference tones and their original normalization, and the
physical source clock. The resulting spectrum is already in watts and uses the
original positive depth weights. The rectangle observes the original periodic
continuation, including across cell edges; it introduces no patch-edge field or
material boundary. Depth convergence remains the caller's responsibility.

The implementation reuses original grouped DD modal reconstruction, centered
harmonic Gram weights, reservoir transfer caches, exponent-scaled contractions
and tone emission. It sums center fields before products and contracts spatial
transform decrements separately to retain weak destructive-fringe variation.
This adds no general relative-error guarantee at an arbitrary cancellation.
Work admission charges `(N*N+1) * depths * sum(F_r*(F_r+1)/2)` spatial products;
the axis caches contain `width²+height²` kernels. All caches, retained DD/scaled
fields and output storage are admitted before allocation, cancellation is
polled throughout, and empty banks prepare no rectangle cache. Spatial synthesis
FFT work is zero; source analysis and original film reconstruction retain their
existing bounds.

Only original causal material reservoirs are admissible: sampled finite-band
optical n+k remains insufficient for this heat law. Reversible storage is not
deposited heat. The full-uniform-body coupled heating guard continues to reject
finite patches until a physical spatial mapping to thermal cells is supplied.

## Finite-aperture Maxwell stress and field energy

`AngularSpectrumExposure{,Ensemble}::measure_aperture_stress` and the two-face
`measure_illumination_aperture_stress` use the original continuous rectangular
aperture/shutter engine to retain the complete symmetric Maxwell stress,
physical field-energy density and signed vector Poynting flux. Real fields and
real factors retain carrier-sum and beat terms; proper factors retain their
original law. Actual total incoming/outgoing fields, longitudinal components,
polarization and spatial interference enter before products. The new accumulator
forms eighteen scalar field products per admitted modal/frequency pair, without
an additional spatial FFT or a full correlation matrix. Original whole-query
source, pair-work, storage and cancellation admission remains authoritative.

Only the original observed lossless nondispersive half-space supplies epsilon
and mu, through the existing material admission. A lossless left observation
remains admissible with a lossy exit; stress in that lossy exit refuses. Dyads,
trace subtraction, traction and area conversion retain double-double arithmetic.
Energy receives no second cycle-average factor; signed stresses are not clamped.
Derived integrated tensors, tractions and normal watts refuse unrepresentable
unit conversions. Observation normals are finite unit vectors in the original
frame and do not rotate or resample the transverse aperture.

The result is an electromagnetic surface observation. Mechanical body force
requires a complete closed-surface balance and the field-momentum rate, together
with the applicable material momentum law. Lateral walls do not cancel for an
arbitrary finite patch; opposite patches alone do not replace that balance.
The represented field remains the original periodic continuation, without a
new diffraction, back-action or finite-beam source model.

Six native tests execute the exact aperture kernels, original stress helpers
and new accumulator with real FFT/math: independent space/time quadrature,
plane-wave flux/traction, a standing electric node with magnetic stress,
real/proper statistics, a narrow dark fringe, and work/cancellation/overflow
refusals. Public original-Maxwell/material tests cover cell/patch partition and
two-face node observation but remain unrun without the constellation checkouts.
Native material/source declarations are signature-only links; this does not
establish full material admission, Maxwell assembly or spatial convergence.

### Finite-aperture traction moments and angular-momentum transport

`AngularSpectrumExposure{,Ensemble}::measure_aperture_wrench` and its
two-face `measure_illumination_aperture_wrench` counterpart integrate the
original Maxwell traction and its first spatial moments over the declared
rectangle and shutter. `AngularSpectrumExposureApertureWrench` reports
`integral(T*ez)dA` [N] and
`integral((r-centre) cross (T*ez))dA` [N m]. The positive-z normal is
explicit. `traction_moment_about_nm` translates the origin using an authored
centre-to-origin displacement and DD force/moment contraction before rounding.
A displaced fringe contributes its actual lever arm; multiplying mean traction
by the rectangle centre cannot replace this observation.

The observer reuses the original lossless nondispersive half-space admission,
two-face modal reconstruction, spatial wavevectors, positive-time phasors and
normalized shutter. Deterministic means and real factors retain carrier sums
and differences; proper-complex factors retain their original difference law.
Independent factors contribute moments separately. Both longitudinal fields
and coherent incoming/outgoing cross terms enter before traction is formed.

Centered first-moment kernels integrate `xi*sin(k*xi)` analytically. A DD
small-width series preserves the derivative near a uniform phase, and split
large-phase divisions avoid squaring an unrepresentable argument. Actual
carrier-plus-envelope wavevectors remain in DD. Separate start/half-width
rotations preserve a tiny aperture phase at a distant coordinate.
No spatial sampling, numerical phase derivative, modal-pair matrix or new
FFT is introduced. The original stress-only query keeps its existing cost.
This query performs thirty scalar field products per admitted mode/frequency
pair and retains normal plus first-moment x/y axis caches. Its work remains
`N^2 * sum(F_bank^2)`; original reconstruction bounds and the additional
first-moment caches/headers are charged before allocation. Banks stream one at
a time, cancellation is polled through preparation and contraction, and invalid
budgets, unrepresentable moments or derived units return errors.

This is a surface transport observation, not a net mechanical body torque.
A complete angular-momentum balance still requires the entire closed surface,
the field-angular-momentum rate and the applicable material law. Arbitrary
finite patches have side walls. The original periodic source continuation,
material restrictions and Maxwell model remain unchanged.

Six new native tests execute this exact production module with the original
aperture/stress kernels and real harmonic-window/math owners; all six existing
aperture-stress regressions also pass. They check an analytic coherent fringe
with nonzero moment despite zero individual-mode moments, independently
integrated two-dimensional mixed-polarization fields under uniform/triangular
shutters, proper versus independent real quadratures, patch additivity and
origin translation, narrow-width derivative and remote-coordinate phase
survival, and work/storage,
cancellation/replay and overflow boundaries. Public original-field tests also compare integrated traction, split-patch
moments about a common origin and a two-face standing node; those tests remain
unrun without the constellation checkouts. Cold source/material declarations
in the focused native build are signature-only links; these results do not claim a
full workspace build or public Maxwell/material assembly validation.

### Coherent source control through finite optical apertures and shutters

`PreparedAngularSpectrumApertureEnergyControl` borrows one fixed original
`AngularSpectrumExposure` and up to 64 original coherent control banks. Its
dimensionless complex coefficients multiply complete SI peak field banks on
their existing common phase clock. The fixed illumination is immutable.
Source materials, layer geometry and transverse grid must match exactly;
equal numerical material values do not replace original source identity.
Both one-sided and two-face observations use the original modal reconstruction.

The objective is the aperture/shutter mean of
`(epsilon E.E + mu H.H)/2`, in J/m³, in the original observed lossless
nondispersive half-space. It is field energy density, not detector power,
absorbed heat or a body's integrated energy. The shutter's original normalized
uniform/triangular profile remains authoritative. A valid left observation
continues to work with a lossy exit; observation in a medium lacking the
required storage law refuses through the original material admission.

`value_gradient` returns that objective and Cartesian source derivatives:
`gradient.re = dJ/dRe(a)` and `gradient.im = dJ/dIm(a)`. `hessian_vector`
applies the real Cartesian Hessian to a complex direction. Finite shutters
retain carrier-sum as well as difference-frequency terms, so this Hessian is
not generally complex-linear. The shared aperture kernels own all temporal
and spatial normalization; no additional optical-cycle half is inserted.

Preparation freezes sparse original bank/frequency modal fields and existing
x/y/time kernels. Each query first combines all coherent contributions in DD,
applies separable even/odd x/y/time contractions, then pulls the covector back
to the original control order. Small odd-window variances stay separate from
near-unit even terms. The scalar energy is contracted against the combined
field, rather than reconstructed from cancelling individual source terms.
No dense source Hessian, source covariance or mode-pair table is formed.
Repeated value/gradient and Hessian actions allocate nothing and perform no
Maxwell/material solve, FFT or trigonometric preparation.

For `N = width*height`, union frequency count `F`, and original bank frequency
counts `F_s`, the contraction work is
`O(N*F*(width+height+F) + N*sum_s F_s)`. Explicit limits admit the original
line reconstruction, sparse source-field entries, total separable kernel
pairs and complete peak logical workspace before the corresponding allocation.
Frequency-vector capacity includes duplicates before deduplication. Source
groups, four union-field work arrays, x/y/time kernels, gradient and original
reconstruction scratch are included; borrowed owners and allocator overhead
are excluded. Cancellation polls source/axis/frequency work and buffer tiles.
An error returns no observation or gradient; the prepared owner can be reused.

This differentiates authored coherent amplitudes with geometry and constitutive
operators fixed. It does not differentiate materials, geometry or thermal
feedback, infer coherence between independent stochastic factors, or construct
a causal broadband controller. The original periodic continuation remains the
field model. DD arithmetic is not a relative-error certificate at every null.

Validation: six passing native tests execute the production action and original
harmonic/aperture kernels. Independent space/time integration checks energy and
both complex source derivatives under uniform and triangular shutters; further
tests check Hessian symmetry/positivity and finite differences, the genuinely
real-linear short-shutter Hessian, exact fixed-source cancellation, empty source
slots, a width-1e-8 node's positive energy and derivative, cancellation/replay,
invalid coefficients and overflow. A 2-D transverse, two-frequency case matches
the existing full mode-pair integrator and its directional derivative. Native
fields are manufactured retained coefficients; original material/Maxwell
preparation is signature-checked, not executed. Two public original-assembly
regressions are authored but unrun, covering source perturbations through
films/evanescent exits, two-face magnetic standing nodes, source identity and
work/storage admission. Full Cargo/DSR remains blocked by absent constellation
sibling repositories.

## Collisional-plasma wavelength and skin-depth scales

`PreparedIonizedGasWave::propagation_scales` now exposes the distinct physical
normal phase wavelength, phase velocity, amplitude skin depth, intensity skin
depth, and phase cycles per amplitude e-fold from the already-admitted complex
`k_z`. The query is O(1), allocation-free, and does not repeat carrier transport
or replace the original frequency-dependent conductivity with a vacuum or DC
surrogate. Undefined finite positive scales are `None`; a monochromatic sample
does not invent group velocity.

A focused regression checks every reported scale directly against the original
complex wave number across three decades of frequency and oblique incidence,
then verifies one phase wavelength against the propagated complex field and one
intensity skin depth against the Poynting decay. This is a homogeneous local
plasma-wave diagnostic, not a pulse/group-velocity, sheath, kinetic-spatial, or
broadband-material claim.

## Original IAPWS water dispersion in native interfaces, films and beams

`plane_wave::optical::WaterOpticalResponse::try_lossless` and `DielectricResponse::WaterOptical` bind the original IAPWS R9-97 phase-index law to temperature, bulk density and absolute SI magnetic permeability from one complete `ResolvedMaterialStatePoint`. Absolute-temperature semantics, original stationary property sources, native state bounds and strict vacuum-frequency/wavelength support remain attached. The caller explicitly selects ordinary-water applicability and a zero-absorption approximation; material names do not infer purity, an EOS, or measured optical loss.

The existing constitutive/interface/film/angular-spectrum owners consume epsilon = n_phase²/(mu*c0²), with the same retained permeability. The original group index and analytic material partials remain separately available. Evaluation allocates nothing, and existing prepared film/beam values retain their normal caching. Multifrequency identity checks retain the original source; the finite-window heat owner explicitly refuses this finite-band real-index law because it supplies no original causal material reservoir.

Six focused actual-source controls pass through the complete native plane-wave subtree: published index and group/phase/mu distinctions; source/state/band refusals; oblique TE/TM, Snell, reciprocity, Brewster and TIR; independent Fabry–Perot and interior Maxwell fields across two water states; coherent three-mode/two-polarization beam fields and Poynting flux; and source-continuity/causal-reservoir behavior. This is homogeneous isotropic finite-band optical dispersion under an explicitly lossless approximation, not measured transparent-water absorption, scattering, phase-mixture optics, or a causal time-domain energy law.


## Original molecular-gas dispersion in native Maxwell waves

`MolecularGasOpticalResponse` and `DielectricResponse::MolecularGasOptical`
connect the original `RayleighGasState` or `MolecularRefractionModel` plus
ordered partial densities to existing interfaces, finite films and coherent
angular spectra. The caller explicitly selects a lossless coherent-medium
approximation and supplies a separate original scalar magnetic source. No
vacuum permeability, material-name mapping, replacement EOS or optical fit is
inserted. `material()` denotes that magnetic source only; `source()` retains
the molecular owner and actual immutable state/densities.

For a fixed gas state, every magnetic condition must match the original
absolute temperature, pressure or mass density, with its static quantity
meaning. The partial-density model owns no temperature/pressure snapshot and
therefore admits only an unconditioned magnetic source. Missing, unsupported
or mismatched conditions refuse instead of being inferred. Ordered density
basis, component cap, dilute limits and the original wavelength band remain
those of the molecular owner; declared composition dispersion reaches infrared
frequencies without extending a constant visible polarizability model.

Preparation evaluates the molecular law once, allocating nothing. The cached
phase permittivity is `epsilon = n_phase²/(mu*c0²)`; group index remains a
separate derivative. All subsequent queries enforce the exact prepared angular
frequency. Separately prepared frequencies share an exposure only when the
original model/state, ordered partial densities and basis, magnetic owner and
stack geometry agree; wavelength-dependent caches do not change source identity.

Rayleigh scattering is not deposited heat or an imaginary permittivity in this
explicit approximation. Native monochromatic heat is exactly zero, while the
finite-window reservoir owner refuses to promote finite-band phase dispersion
to an unsourced causal storage law. Molecular attenuation, changing gas state,
spatial GRIN Maxwell propagation and measured-gas qualification are not supplied.

Focused G1/G3 regressions use independent Lorentz--Lorenz/Cauchy phase/group
equations, oblique TE/TM refraction and TIR, Fabry--Perot transmission and
interior fields, coherent multimode beam phase/power, zero true heat despite
nonzero Rayleigh cross section, and source/domain/frequency admission. These
regressions are authored; native execution is not claimed by this extension.

### Absorbing exit media: coherent local heat and analytic substrate watts

`PlanarAngularSpectrum::transmitted_absorption_field` reconstructs the original
transmitted E/H plane and reports local dielectric/Joule heat
`q = -omega Im(epsilon) |E|²/2` for peak exp(+i omega t) phasors. The complete
vector E includes longitudinal TM fields and coherent spatial interference.
The original exit permittivity already owns its Ohmic, dielectric or optical
loss; no separate conductivity or fitted attenuation is added. A lossless
evanescent field has exactly zero heat despite its spatial decay.

`transmitted_absorbed_power_w(z0, z1)` integrates that loss over the entire
periodic cell and a finite depth interval or infinite substrate tail. Modal
orthogonality applies only to the whole-cell integral. Each positive mode
contributes its entering watts times exact exponential transmission and
`-expm1(-2 alpha dz)`. A logarithmic fallback preserves representable watts
when intermediate transmission or opacity underflows. The O(modes) query
allocates no workspace and needs neither a depth mesh nor additional FFTs.

Local heat adds one explicitly admitted scalar plane and O(grid points) work
to the existing six-FFT field reconstruction. Scaling E by the square root of
the loss coefficient before squaring preserves weak-loss heat at large finite
amplitudes. Finite films reuse this same owner through the existing
`AngularSpectrumLayerField` alias. Depth, work, storage and cancellation
refusals publish no partial observation. Finite-film absorption remains
separate; these queries do not evolve temperature or alter constitutive state.

Five integrated regressions cover volume/Poynting/film closure, oblique Ohmic
skin heating, weak absorption, lossless TIR and admission/cancellation.
Four native arithmetic tests passed with actual fs-math before the final
deep-tail/opacity fallback; its additional decimal-reference regression and
the full integration tests remain unrun after execution-host failure.

### Partially coherent illumination heats the same exit substrate

`AngularSpectrumCoherence::transmitted_absorption_field` streams the original
coherent factor fields through that same substrate-loss owner. Full vector and
longitudinal interference remain inside each factor; mutually uncorrelated
factors add only as heat and second moments. The returned
`AngularSpectrumCoherenceLayerField` describes the queried exit half-space as
well as finite-film observations. Its local W/m³ and transverse W/m values are
not additional heat beside a volume's integrated watts.

`transmitted_absorbed_power_w` adds the original factors' analytic depth
integrals in double-double precision, including infinite absorbing tails and
exactly zero lossless-evanescent heat. This whole-periodic-cell query allocates
nothing and uses O(factors * grid points) work without FFT or depth quadrature.
Local observations use the bank's already admitted whole-factor FFT/work and
workspace bounds, which include the streamed and output heat planes. Material
identity, finite-film heat, source coherence and temperature remain unchanged.
Invalid depths refuse even for an empty ensemble; cancellation publishes no
partial result and the original bank remains reusable.

Five focused G0/G1/G3/G4 regression groups use original source/material/Maxwell
owners: independent lossy-interface fringe formulas, one-factor parity and
coherent versus independent cancellation, complex unitary factor mixing with
volume/Poynting/film power balance, lossless TIR, and exact workspace admission,
darkness and transactional refusal. These synthetic constitutive checks do not
qualify a measured substrate or infer thermal evolution.

### Two-point electric coherence and explicitly analyzed interference

`AngularSpectrumCoherence::prepare_interferometer` extends those same Maxwell
observations to 1–4096 explicitly analyzed arms. Preparation retains the DD
complex analyzer amplitude for every original factor/arm pair. Readings apply
the supplied complex amplitude gains and sum each factor's arm fields before
its positive cycle-averaged square. All coherent arm cross terms survive;
independent source factors still add only as second moments. No splitter gain,
free-space path, energy normalization or source statistics are inferred.

For F factors, A arms and N Fourier modes, preparation admits F*A*N mode
reconstructions and retained O(F*A) storage before allocation. Repeated gain or
phase scans cost O(F*A), without allocation, FFT, material solves or new field
reconstruction. Limits include the retained amplitudes and the existing
conservative reconstruction workspace. Every arm is validated for dark banks,
and cancellation checks bound arm/mode work and precede publication. This is
monochromatic electric mean square, not detector watts or back-action; changed
sources, observation geometry or analyzers require preparation again.

Focused G1/G3 tests compare three-arm fringe scans with the analytic geometric
sum for independent uniform and oblique waves, verify one-/two-arm parity,
retain a 2^-59 dark-port residual and a common gain-phase invariance. G0/G4
checks exercise exact work/storage boundaries, invalid dark-source arms and
gains, cancellation and source reuse.

The original monochromatic source-factor owner now exposes
`AngularSpectrumCoherence::mutual_coherence` between two exterior observations.
Each observation supplies a side, a finite nonnegative distance from that side's
film face, and finite physical x/y coordinates in the original periodic cell or
its continuation. The same per-mode Maxwell reconstruction serves both existing
FFT planes and these direct samples: material response, longitudinal fields,
critical transmitted modes and incident/reflected propagation signs are shared.
Prepared inverse-FFT coefficients already include 1/N; direct synthesis uses
exp(-i(kx*x+ky*y)) without another normalization. DD arithmetic retains the
authored wave-number/coordinate products through phase reduction and field sums.

The returned local blocks and cross block retain
W_ab[i,j] = sum_m E_a,m[i] E_b,m[j]* / 2 in V²/m². Swapping observations takes
the conjugate transpose. An explicit complex Cartesian analyzer is w†E, including
conjugation of its dimensionless weights. Projected moments therefore report
Gamma_ab = <a b*>, and combining a + exp(+i phi)b gives
Ia + Ib + 2 Re(exp(-i phi)Gamma_ab). Visibility is 2|Gamma_ab|/(Ia+Ib);
it is undefined for two dark arms and zero for exactly one dark arm. Projection
checks nonnegative finite intensities and the Cauchy-Schwarz bound, with a bounded
final-rounding allowance. The direct recombination query combines each factor's
fields before its positive square, preserving weak dark-port residuals that
subtraction of rounded local/cross moments can erase. Uncorrelated source factors
still add only as second moments. These are electric mean-square observables,
not detector watts, an inferred polarizer output state, or detector back-action.

Each query admits exactly 2*M*N mode reconstructions for M factors and N original
Fourier modes. It uses no FFT, heap allocation, fresh material solve or dense
spatial covariance. The fixed logical numerical workspace includes result
tensors, two synthesized vector fields and a mode's E/H buffers; scalar compiler
temporaries, borrowed source/factor storage and earlier retained results are
excluded. The direct query admits the same conservative workspace bound.
Cancellation is checked on entry, between factors, at most every 32 mode pairs,
and before publication. Rank-zero darkness still validates geometry, analyzers,
relative phase, budgets and cancellation. No isolated-aperture, continuum-
convergence or relative-error guarantee at arbitrary cancellation is inferred.

The focused correlation tests supply independent Young-style visibility with
identical arm intensities, a positive 2^-59 dark-port residual, lossy Airy
reflection/transmission and standing-wave phases, evanescent phase/decay,
two-dimensional FFT normalization and longitudinal-field parity, conjugate
symmetry, unitary factor rotation, complex analyzer handedness and admission
boundaries. Native execution of this addition was blocked by an unavailable
execution host; these cases are authored regressions, not a claimed passing run.

## Finite detectors for the original partially coherent beam

`AngularSpectrumCoherence::measure_aperture` observes the original ensemble
through a finite rectangle using the same vector diffraction and analytic
rectangle operator as a coherent beam. Every factor retains its spatial
interference; independent factors add signed Poynting components and original
integrated watts in DD, without summing their fields or reapplying a half factor.
This retains partial fringe visibility and complex unitary factor-basis
invariance. Reflected transport keeps its physical negative sign.

For F factors on N=W*H modes, the query admits aggregate F*N*(W+H) axis
contractions before reconstruction. The common W²+H² spatial kernels are
prepared once from the exact original source/grid. Each factor streams through
one shared operator; peak numerical buffers comprise one coherent-beam
workspace, the axis caches and one DD contraction plane. No additional FFT,
material solve, factor-sized scratch bank or sampled-field interpolation is
introduced. Empty factor banks validate the depth and check cancellation, then
return darkness with zero query budgets and no detector allocation.

The three G1/G3/G0/G4 regression groups in `tests/plane_wave/coherence/aperture.rs`
cover independent analytic fringe and vector-flux references, whole-cell power,
complex factor-basis changes, original coherent-owner parity, signed reflection,
aggregate work and peak-storage admission, rank-zero depths and cancellation
with replay. Detector heat, quantum efficiency, finite-edge convergence and
source temporal coherence remain outside this monochromatic periodic-field
observation model.

## Axisymmetric dielectric fields and capacitance

`electrostatic::axisymmetric::AxisymmetricElectrostatics` solves linear
electrostatics on the original `TriComplex2` meridian with coordinates `(r,z)`
and an explicit `Metric2::Axisymmetric` angular span. It reuses the same
source-resolved material laws, native weighted Whitney Hodge, scalar boundary
elimination, checked Jacobi-CG and charge/energy recovery as the planar owner.
The public owners reject each other's metrics and return distinct field types;
no planar force or moving-geometry API is exposed for a revolved field.

The native Hodge exactly integrates the affine radius in
`dV = angular_span * r dr dz`. Constant free-volume-charge P1 loads use
`Q_i = Q * (1 + r_i / sum(r))/4`; surface-charge and normal-displacement loads
use `Q_i = Q * (1 + r_i / sum(r))/3` on the two endpoint radii. Equal nodal
shares would be incorrect. Radial ratios are scaled before summation. Fixed
impressed polarization uses the same exact weak gradient load and contributes
to D, while storage remains `integral E.epsilon.E/2`; its separate E.P pairing
does not invent source energy. The zero-area symmetry axis cannot anchor an
electrode or carry an authored nonzero surface displacement source.

E and D are local radial/axial components, independent of azimuth. Directional
material frames are explicitly meridional and rotate with azimuth; the omitted
azimuthal response must be decoupled. Charge, capacitance and energy already
include the supplied angular span. TAU gives a complete body of revolution;
a smaller span describes its symmetry sector, not a separately driven wedge.
The two-terminal capacitor checks Q/V against 2U/V² and lowers into the original
native circuit branch, with one storage account for the same field.

`axisymmetric::network` extends that same radial operator to complete Maxwell
capacitance matrices, repeated voltage-field recovery and mixed voltage/total-
charge floating conductors. Each conductor must have its own positive-area
exterior trace; an axis-only conductor cannot borrow support from edges joining
different terminals, including in the two-terminal API. Preparation reuses the
existing bounded unit solves and floating-block factorization; queries retain
trapped-charge offsets, physical storage and independent component references
without another FE solve. Recovered fields keep the axisymmetric type. The full
native circuit block retains total electrode charges and the original radial
storage, including floating-node charge conservation and resistor heat. Branches
and block are one store, not extra field energy. G0/G1/G3/G4 regressions in
`tests/axisymmetric_electrostatic/network.rs` cover analytic disk and coaxial
conductors, direct-solve parity, angular span, trapped charge, independent large
biases, reference refusal, native RC power and cancelled preparation. This is
the frozen, source-free linear capacitance response; it adds no changing-material
work, distributed leakage, fringing outside the mesh or mechanical force law.

G0/G1/G3/G4 cases in `tests/axisymmetric_electrostatic.rs` cover coaxial and
layered-coaxial capacitance, charged-cylinder refinement, exact axial fields
including the axis, sector scaling, exact radial volume/surface moments, free
sheet charge, impressed polarization, large common voltage offsets, an oblique
meridional tensor, original source retention, cancellation and admission limits.
The model is fixed-geometry linear electrostatics with explicit algebraic
acceptance, not a continuum error bound or measured-device qualification.
External fields absent from the mesh, azimuthal variation, dielectric memory,
conduction, breakdown, electrostriction and electromechanical force/motion
remain outside this owner.

## Distributed Debye dielectric memory

`electrostatic::dispersion::PlanarDebyeDielectric` advances source-backed,
isotropic multipole polarization on the original extruded P1 triangle mesh.
Every cell retains its complete `ResolvedDebyeMaterial` and a vector state per
original pole. `D = epsilon_infinity E + sum(P)`; polarization is not absorbed
into a fictitious fixed permittivity or counted again as Ohmic conduction.

`prepare_step` eliminates the local midpoint pole equations into the existing
weighted FEEC operator. Reusing that prepared object avoids Hodge assembly on
subsequent equal-sized steps; no global polarization unknowns or alternate
linear solver are introduced. Native checked CG solves only the potential
coordinates. Geometry, coefficients and external material conditions are fixed.

Unprescribed nodes retain their initial weak free charges. The original complete
voltage traces exchange physical charge with linearly ramped voltage sources.
Each candidate reports cell polarization heat, original field and pole storage,
boundary charge exchange, and source work. Independent source-work/energy and
free-charge budgets gate publication. Gauge-separated fields preserve physical
response under independent reference changes on disconnected components.
Cancellation, invalid trace changes and foreign prepared operators leave the
accepted field, polarization and clock untouched.

`tests/dielectric_dispersion.rs` checks held-voltage and isolated-charge analytic
relaxation, second-order convergence, heterogeneous field redistribution,
multiple poles, the instantaneous zero-pole limit, source power, independent
component gauges, deterministic replay and transactional refusal (G1/G3/G4).
This is synthetic numerical verification, not measured dielectric qualification.
Midpoint does not promise monotone relaxation at unresolved time steps; resolve
the pole times and refine time and mesh separately. Changing material/temperature,
anisotropic memory, magnetic induction, carrier transport, breakdown,
thermal evolution and mechanical work require their own coupled owners.
Ohmic leakage is available through the opt-in total-current lane below. Local heat is available to a thermal consumer but is not silently
deposited into an unspecified reservoir.


## Distributed total-current dielectrics in planar, axisymmetric and 3D geometry

`electrostatic::dispersion` retains the planar API and adds
`AxisymmetricDebyeDielectric`; `spatial::electrostatic::dispersion` adds
`SpatialDebyeDielectric` on conforming P1 tetrahedra. All three use one
constitutive/time-evolution kernel and their original native weighted Whitney
operators. Axisymmetric volumes include the authored angular span. Physical
voltage traces must cover complete exterior edges/faces; an axis-only,
point-only, edge-only 3D or unanchored component cannot act as an electrode.

The default constructors explicitly model insulation. The default-off
`resolved-leaky-dielectric` feature exposes `new_leaky` with one original
`ResolvedLeakyDielectric` source bundle per cell. It admits independently
sourced nonnegative isotropic DC conductivity, including exact zero. Free
nodal charge follows total-current continuity: its stored-charge increment
plus integrated Ohmic outflow must vanish. This includes Maxwell charge
relaxation and charge buildup at heterogeneous dielectric interfaces.
Prescribed traces supply both stored electrode charge and transported charge.

Local midpoint pole updates and conductivity are condensed into a positive
effective-permittivity operator, without global polarization unknowns or a
fictitious conductivity floor. A prepared operator is reusable at one fixed
time step for the same immutable geometry/material owner. Original local
polarization loss and Joule heat are recovered independently, retained
separately and summed once. Stable storage increments and separate source-work
contributions avoid subtracting large biased total energies; only their
combined source work is gauge independent. Charge, work and energy acceptance
and cancellation remain transactional. Native candidates expose physical
cell volumes and source identities for receiving coupled owners.

`tests/dielectric_leakage.rs` adds six analytic/G1/G3/G4 groups covering DC
through-current, Maxwell relaxation, Maxwell-Wagner interface charge, combined
Debye/Ohmic eigenmodes, sourced-zero parity and axisymmetric scaling.
`tests/spatial_electrostatic/dispersion.rs` adds eight groups covering oblique
3D fields, physical volume/orientation, heterogeneous and isolated relaxation,
independent component gauges, trace admission, cancellation and Ohmic loss.
These are synthetic numerical checks, not measured material qualification.
Fixed geometry and source coefficients are required. Induction, carrier
transport, anisotropic dispersion/leakage, breakdown and moving geometry are
outside this owner. Passive second-order midpoint does not guarantee monotone
relaxation for unresolved time steps; refine against the physical pole and
charge-relaxation times as well as the spatial mesh.

## Whole-film admitted original ionized-gas illumination

`PreparedIonizedGasLayerIllumination::try_new` binds the original finite dielectric
stack to fixed simultaneous complex peak tangential electric inputs at its
two faces. Its exterior ports must be lossless, propagating and non-ionized;
finite films may retain original ionized-gas, lossless, Ohmic, Debye or optical
material responses. The existing unrestricted normalized layer constructors
continue to refuse ionized gas. Their private construction seam is accessible
only to the admitted owner and never escapes in its public results.

The original coherent Maxwell solution supplies both boundary traces and
positive integrated film loss. Each positive-thickness gas film then receives
an analytic conservative envelope for its complete vector electric field.
Counterpropagating waves decay from their respective faces; near the critical
normal wavenumber the existing entire E/H propagator provides a bounded
continuation. TM admission includes its longitudinal electric component.
This admits the original weak carrier drift and Coulomb pair-relative drift
for the whole film, including interior standing-wave antinodes. Separate
admission of the incident beams, endpoint-only checks and sampled-depth checks
do not replace this condition. The envelope is conservative and may refuse an
actually admissible field; it is not an interval certificate or a claim that
its bound is attained.

Each admitted film retains its original source descriptor, physical two-face
field kernel and original isotropic harmonic carrier response. Local queries
return peak E/H, RMS squares, signed complex Poynting flux, positive absorbed
power, free-carrier current, every original carrier velocity and the original
electron/heavy thermal partition. Original signed exchange corrections are
included once. Drift kinetic energy, reactive carrier power and lossless
background field storage remain separate; a complex effective permittivity
is not interpreted as stored energy. Local field access is read-only so the
validated field cannot be replaced before cached carrier rescaling.

Integrated carrier heating uses the original positive slab-loss result and
the original material's thermal fractions, including coherent source
interference. It does not subtract nearly equal boundary fluxes or introduce
depth quadrature. Repeated point and heat queries allocate nothing, perform
no stack solve and solve no carrier modes. The prepared wave cannot accept a
new source amplitude without fresh whole-film admission. Zero-thickness films
retain source identity and zero integrated heat but have no interior volume.
Cancellation and any source/material/domain failure publish no partial wave.

The physical domain remains frozen populations with local unmagnetized
collisions and the original source applicability bounds. This is a coherent
plane-wave film capability, without evolving ionization, sheath dynamics,
magnetized response, transient reservoirs, arbitrary angular beams or lossy
and evanescent exterior input ports.

Focused regressions specify independent bidirectional Fabry-Perot amplitudes
and closed-form volume heating for TE/TM below, near and above the original
plasma frequency; pointwise original carrier responses and energy storage;
multiple films with Coulomb heat partitions and zero-thickness identity; and
coherent inputs whose individual amplitudes pass but combined internal field
exceeds the original drift domain. Kernel regressions cover analytic field
bounds through standing-wave, critical, evanescent and branch-transition
limits. Pinned formatter/parser checks passed; native Rust regressions were
added but were not executed in the shared editing environment.

## One-sided ionized films with passive exit half-spaces

`PreparedIonizedGasLayerIllumination::try_new_one_sided` binds one actual
finite left-face drive and requires the right input to be exactly zero.
The incident medium remains original and lossless; the original non-ionized
exit may be absorbing, evanescent or exactly critical. This adds no fictitious
incoming channel in a lossy medium and preserves the simultaneous-source
constructor's lossless propagating port rule.

The selected-film and whole-stack owners now share one native bounded field
preparation. Its backward and forward sweeps cost O(layers), followed by O(1)
field-kernel selection per film. Each actual amplitude-scaled film passes the
same whole-vector analytic carrier envelope and original thermal partitions.
No normalized ionized layer solver or replacement source amplitude escapes.
Native positive slab losses supply film heat, and physical boundary powers
must satisfy the original Poynting balance.

The exit cache retains independent electric and power-conjugate magnetic
traces. At exact TM critical refraction, vanishing tangential E and normal
wave number do not erase finite tangential H or longitudinal E. `exit_field`
continues this actual outgoing field with the original passive complex phase.
`absorbed_exit_power_to_depth_w_m2` integrates the same attenuation stably
with expm1; finite-film heat remains separate and a lossless evanescent or
critical exit has exactly zero absorbed power. Repeated field and depth-heat
queries allocate nothing and perform no material, layer or carrier solve.

The original frozen-population, unmagnetized weak-drift, homogeneous finite
film and prescribed-plane-wave scope remains. Ionized exterior half-spaces
are not admitted by this finite-film owner. All query depths must be finite
and nonnegative; cancellation or admission failure returns no partial wave.

Three focused regression groups specify original multi-film/two-port parity;
independent dimensional Maxwell/Fabry-Perot amplitudes for absorbing,
evanescent and critical TE/TM exits, local Joule loss and shallow/opaque
substrate heat; and standing-wave drift rejection, identity/silence,
cancellation and source/depth/port restrictions. Pinned Rust formatting/parser,
bounded API/physics review and independent fixture-based power-balance
calculations passed. Native Rust regressions were authored but not executed
in the shared editing environment.

## Conservative heat cells inside admitted ionized-gas films

The admitted film owner now exposes
`layer_absorbed_power_in_depth_interval_w_m2` and
`layer_carrier_heating_in_depth_interval` for its original physical wave.
Both simultaneous two-face sources and one-sided illumination with a passive
exit retain the same amplitude and whole-film drift admission.

Each query restricts the existing positive slab loss operator to two actual
cached E/H traces. Original wave coefficients and reference admittance are
retained at preparation; DD endpoint differences preserve the interval phase.
The two positive loss modes are summed in DD. This O(1), allocation-free
calculation preserves standing-wave interference and longitudinal TM loss
without depth quadrature, boundary-flux subtraction, new material evaluation
or another stack/carrier solve. A complete-film interval returns the cached
original result exactly; disjoint cells recover that integrated film power.

Gas-cell heat retains the original electron/heavy fractions and signed
exchange correction. That correction is included in the partition and must
not be deposited a second time. Empty gas cells and zero-thickness identity
films return an explicit zero partition; ordinary films return no carrier
partition. Invalid indices or depth ranges refuse before either outcome.
A thermal consumer owns area/time conversion and one deposition per cell;
these spectral cycle means do not evolve a reservoir or population.

Three focused regression groups specify independent coherent two-wave volume
integrals for TE/TM, positive cell heat below boundary-flux subtraction
resolution, cell additivity and exact whole-film recovery; original Coulomb
partitions under absorbing/critical exits; and opaque skin-depth heat,
zero-volume/lossless cases and depth/index restrictions. Pinned Rust
formatting/parser checks and bounded source/API/math review passed. The
native Rust regressions were authored but not executed in the shared editing
environment.

### Differential winding inductance at a nonlinear 3D bias

`spatial::vector_potential::nonlinear::winding::SpatialIncrementalWindingBasis`
extracts local self/mutual inductance from the original source-bound Whitney
magnetic problem. It solves the physical operating currents, permanent
remanence and prescribed tangential traces, then reuses one sparse constitutive
Hessian and one coexact perturbation topology across at most sixteen windings.
Perturbation traces vanish on the same prescribed edges. The reciprocal block
`L_ij = delta_A_i^T K_t delta_A_j` must agree with every ordered original
Whitney current pairing. Both longitudinal and transverse differential
reluctivity participate; secant permeability, `lambda/i` and `2W/i^2` cannot
substitute for this derivative.

Profiles retain turns, polarity and complete current closure, including
prescribed exterior faces. Explicit omitted end self-inductance is added
without invented leakage. Operating linkage remains separate from local flux
increments and voltages. Source corners and range endpoints refuse a unique
two-sided response. Recoil laws require an interior physical bias; derivative
fields are infinitesimal responses, not finite one-ampere magnet states.
Allocation bounds, original field/gauge allowances, complete-kernel
compatibility, post-gauge residuals, div-curl checks and cancellation remain
active. Dependent observation coordinates are allowed without claiming an
invertible terminal matrix.

`spatial_incremental_winding` defines regressions for the linear limit,
nonlinear flux finite differences, reciprocal tensor energy, polarity/gauge
invariance, fixed bias, source-domain refusal, dependent modes, bounded work
and cancellation. This is a fixed-geometry, fixed-material-state local
derivative. It does not supply finite-amplitude circuit evolution,
frequency-dependent impedance, hysteresis, propagated source uncertainty or
a spatial-discretization certificate.

## Photoionization state in existing electrical consumers

The `ionization_feedback` integration regression also specifies illumination
of a gas with exactly zero initial free electrons, using the published H I
bound-free source and the photon-conserving `fs-gas` column update. The returned
original thermochemical state supplies the existing momentum-relaxation,
DC conductivity, harmonic carrier response and FEEC conduction owners.
Independent Drude current and electrical-work/heat accounts use the actual
photo-created carrier populations. A subsequent electron-impact interval uses
the same atom/ion/electron channel and retains the absorbed photon energy in
the original gas energy account.

These are sequentially declared illumination, field evaluation and collision
intervals. The regression does not assert a simultaneous radiation/field/
chemistry integrator or attach an arbitrary conductivity to illumination.

## Coupled unresolved end-winding magnetic storage

`winding::end::EndWindingInductance` admits at most 16 windings and 16
caller-authored shared end-flux modes:
L_end = diag(L_self) + sum(ell*n*n^T). Finite nonnegative coefficients
and signed participation define passive storage directly. Dependent and zero
modes remain valid; no inverse, leakage floor or regularization is invented.
These terms describe only magnetic flux outside the resolved field domain.

The complete symmetric matrix is inspectable. Linkage, storage, finite
increments and work use compensated projections of the original modes,
retaining cancellation between opposing currents. Stable energy changes pair
modal midpoint current with its increment instead of subtracting large endpoint
energies. Explicit current/mode bounds keep this terminal work independent of
mesh size. Unrepresentable arithmetic refuses.

`SpatialNonlinearWindingCircuit::new_with_end_modes` retains this same model
for initial linkage, terminal voltage, accepted storage and backward-Euler
numerical dissipation. Source/conductor/material refresh and graph-resistor
reassembly preserve the model. Existing constructors retain independent
diagonal end inductances. End numerical dissipation is never conductor heat,
and the resolved field and graph retain their original work/heat owners.

Focused integration regressions cover independent coupled-loop evolution,
saturated open-secondary EMF, reversed winding polarity, source refresh,
resistor reassembly and late cancellation. Six native kernel regressions
plus five independent actual-source tests passed, including 1,000 expanded
matrix/work comparisons. The enclosing field/circuit integration regressions
are authored but could not be executed in the available dependency environment.

### Shared end flux in planar, solid and moving winding histories

`NonlinearWindingCircuit::new_with_end_modes` and
`new_with_solid_regions_and_end_modes` use the same passive end model for
stranded and source-normalized solid windings. Initial linkage, finite terminal
induction, storage change and backward-Euler numerical dissipation all use the
original signed flux factors. The independent diagonal, complete matrix and
original modes remain inspectable.

Graph-resistor reassembly, conductor/material refresh and prescribed-motion
candidates retain the model and original field/circuit histories. Conductivity
reweighting still pairs F_new with delta(A); changing source weights alone
cannot invent voltage in a stationary biased field. Motional voltage and work
remain the original solid-field transfers. The new `stored_change_j` receipt
uses directly evaluated end-mode work instead of differencing large end-energy
endpoints. End numerical dissipation remains separate from physical Joule heat.
Legacy constructors retain diagonal-only end storage; end parameters stay fixed.

Three focused new regressions cover coupled nonlinear loops with native
capacitor storage, reversed mixed-causality windings through live hot-source
refresh, and continued solid diffusion, conductivity reweighting and nonzero
motion. Independent loop equations confirm both induced-current polarity and
whole-system physical work, rather than merely matching two implementations.

All changed source/test files passed native rustfmt. The shared end kernel's
11 native tests remain the executed numerical evidence. The enclosing 3D and
planar field/circuit tests could not run because the available fs-em dependency
graph failed to load its production mesh/material/FEEC/solver crates.

## Escaping recombination radiation and carrier feedback

The source-bound ionization feedback regression now continues the original
photo-created and collision-evolved gas through spontaneous radiative capture.
The shared frozen gas publishes the depleted electron population to the existing
collision material, lowering electrical conductivity. Gas energy retained since
photoionization plus escaping recombination photons closes the original absorbed
photon account. This composes the existing owners for successive prescribed
intervals; it does not claim simultaneous radiation/field/chemistry integration.

The spectral-ionization feedback regression uses two bands of the original
Verner H I source, checks their distinct transmitted fractions and chemical
plus electron-sensible energy, and passes the same accepted populations and
temperatures into the existing collision-derived material. It specifies the
resulting harmonic current/heat/storage and planar current/Joule-heat balances
without assigning an illumination-dependent conductivity or resetting to Saha
equilibrium. The electrical field is evaluated after the spectral exposure;
the test does not claim a continuously coupled radiation/electrical transient.

## Oriented magnetic permeability in current-free three-dimensional fields

`SpatialMagnetostatics::with_magnetic_laws` retains an explicit scalar or
orthotropic magnetic law per original tetrahedron. Original principal material
sources and their fixed physical frame produce the complete Whitney Hodge
`d0^T M_mu d0`; recovered flux is `B=mu H+Br`, and induced storage is
`H.(mu H)/2` on the original cell volumes. The existing `new` entry point keeps
the scalar path, including its original recoil arithmetic. Recoil sources keep
their own isotropic slope and operating envelope; an unrelated principal law
cannot be attached to a magnet.

Material and orientation preparation is reused across boundary excitations.
Existing disconnected gauges, dimensional residual/energy budgets, bounded
assembly, cancellation and raw interface-flux diagnostics remain unchanged.
Surface and eggshell Maxwell loads require their original homogeneous isotropic
exterior material; an anisotropic target may be enclosed, but the target's
magnetic response supplies no universal material stress law. G1/G3 regressions
cover off-diagonal affine fields under potential and flux traces, layered
normal-flux/tangential-field continuity, source-only substitutions, physical
co-rotation and energy, scalar parity, and stress-domain refusal. They qualify
the discrete fixed-state linear law, not experimental anisotropic material data,
nonlinear anisotropy or moving-frame thermal/mechanical exchange.

## Oriented reluctivity in current-carrying fields, windings and induction

`SpatialVectorMagnetostatics::with_magnetic_laws` carries the same original
principal sources and fixed magnetic axes into the edge curl-curl operator.
The native face Hodge uses the complete inverse tensor, formed from reciprocal
principal values. Field recovery applies `H=nu(B-Br)` and reversible energy
uses `(B-Br).H/2`. Winding self/mutual terms use `B_i.nu B_j` with the original
Whitney source pairing. Explicit shared end-flux factors remain independent of
the resolved field. No scalar average or componentwise Cartesian inverse is
introduced; the original isotropic/recoil constructor remains available.

The retained A-phi induction owner uses those same coefficients in evolution,
endpoint fields, stored energy and its conservative magnetic Peclet bound.
Temperature continuation keeps the original material identity, selected claims,
nonthermal conditions and magnetic axes. Its endpoint-symmetric parameter work
is `V B0.(nu1-nu0)B1/2`, computed directly from principal coefficient changes;
held sources give exactly zero. This signed work remains separate from Ohmic
heat and impressed/motional work. Recoil refresh and moving material axes refuse.
Mechanical surface/eggshell sampling retains its isotropic exterior collar.

Focused manufactured regressions cover rotated uniform fields and gauge
invariance, a six-tetrahedron current/induction reduction with independently
derived stiffness and charge-eliminated mass, signed thermal parameter work,
and reciprocal winding energy/material scaling. These extend the existing
lowest-order bounded-domain model; anisotropic saturation, hysteresis,
frequency-dependent permeability and experimental qualification remain outside
this linear fixed-frame law.

## Principal-plane uniaxial refraction

`plane_wave::uniaxial::PrincipalPlaneUniaxialInterface` binds the original
lossless isotropic incident response to a transparent uniaxial half-space.
`UniaxialOpticalResponse` retains two explicitly ordered original
`SpectralComplexIndexStatePoint` projections, their same complete material
snapshot, source wavelength/property roles, interpolation-gap limit and
same-state scalar permeability. Either axis with nonzero extinction refuses
at the queried frequency. The authored optic axis lies in the incidence plane;
TE excites the ordinary mode and TM the extraordinary mode.

The extraordinary wave satisfies the rotated index ellipsoid. Positive normal
Poynting flux selects its outgoing branch; a negative normal phase component
is retained when physically appropriate. Results distinguish phase-normal and
energy-ray directions, preserving optical walk-off. Signed complex tangential
E and conjugate H amplitudes, SI wavenumbers and independently evaluated normal
power fractions retain the original exp(+i omega t), exp(-i k.z) convention.
Critical TE/TM fields remain finite without evaluating singular admittances;
evanescent modes decay into the crystal and publish no fictitious bulk ray.
Power-of-two scaling and double-double discriminants protect critical-side
classification and weak positive transmission. Queries allocate nothing.

The scope is one homogeneous principal-plane half-space interface, including
an exact isotropic constitutive limit. It does not extend the existing isotropic
film, angular-spectrum or rendering owners to arbitrary anisotropic media.
Biaxial crystals, absorbing anisotropy, mode-coupled films and broadband
pulse evolution require separate capabilities.

### Local dispersive crystal group transport

`UniaxialOpticalResponse::group_indices` differentiates both original
principal spectral laws in their authored wavelength coordinate.
`PrincipalPlaneUniaxialInterface::group_velocity_response` and its signed
tangential-index counterpart return the unchanged native carrier together
with both principal group indices and local `grad_k(omega)` in [x,z] m/s.
The result borrows its complete interface owner. No second Snell/root solve,
fitted spectrum, angular finite difference or allocation is introduced.

The extraordinary velocity differentiates the full dispersion relation
`k^T M(omega) k=omega^2/c0^2`, retaining both ordinary and extraordinary
spectral slopes and native crystal walk-off. The ordinary mode reduces to
`c0/ng_o` along its phase direction; equal phase indices do not erase unequal
principal slopes. Positive dispersive energy aligns group and Poynting
directions in this transparent reciprocal model. Anomalous signs remain
signed; they are not a signal-front or positive energy-speed claim.

Original source endpoints, interpolation cusps, gaps and extinction remain
authoritative. Critical modes may have finite tangential group velocity with
zero normal power. Evanescent modes and exactly vanishing dispersive
denominators have no finite bulk group vector; unrepresentable nonzero
velocities refuse. This is local homogeneous principal-plane transport,
not a normal-transfer phase delay or broadband propagation solver.

### Local crystal normal-transfer dispersion

`PrincipalPlaneUniaxialInterface::frequency_dispersion` and
`frequency_dispersion_at_tangential_index` return a borrowed
`UniaxialFrequencyDispersion` with the unchanged native interface carrier and
complex first and second derivatives of the transmitted normal wavenumber.
The caller explicitly chooses `FixedAngle` or
`FixedTransverseWavenumber`, using the same constraint type as isotropic
dispersion. Fixed angle retains the original incident spectral derivatives;
fixed dimensional transverse wavenumber needs no incident-law derivative.
Both original crystal spectra retain their own slopes, curvatures, source
support, interpolation gaps, endpoint/cusp and loss refusals.

The analytic kernel differentiates the original outgoing eigenbranch using
fixed-size double-double contractions and one common binary index scale.
It preserves the original carrier and performs no neighboring-frequency
boundary solve or frequency differencing. Factored discriminant roots and the
original incidence cosine retain matched grazing limits even when the normal
energy contribution is smaller than the rounded tilted normal phase.
Principal limits and common spectral scaling use original bulk wavenumber jets
directly, preserving exact zero slope or curvature. Exactly critical or
unresolved roots refuse the derivative query without removing the available
native carrier.

`transmitted_propagation` uses the existing signed `BulkPropagation` owner
for the native exp(-i*kz*z) transfer, normal phase/group delay, group-delay
dispersion and evanescent logarithmic attenuation derivatives through a
nonnegative normal distance. A tilted branch may have negative or zero normal
phase without being critical. The outgoing evanescent sign is retained.
These local normal-transfer slopes are distinct from travel along the
Poynting ray; fixed-angle queries include the changing transverse wavevector.
Fresnel boundary phase, finite-crystal interference and broadband pulse
evolution are outside this propagation-only result.

Nine new native numerical regressions plus the 55 original principal/oriented
carrier tests passed in a bounded build against the original `fs_math`.
They cover constraint-dependent geometric dispersion, unequal principal
spectral slopes, evanescent branches, tilted zero phase, matched grazing and
exact stationary/zero-GVD bulk identities. Source-bound integration fixtures
were added and parsed/reviewed; full `fs-em` material integration was not
executed in this editing environment.

### Principal-crystal Fresnel boundary dispersion

`PrincipalPlaneUniaxialInterface::boundary_frequency_dispersion` and its signed
tangential-index counterpart retain the unchanged native carrier and normal
propagation response, and add raw complex first/second frequency derivatives
of reflected tangential E, transmitted tangential E and normalized conjugate H.
The existing `FresnelAmplitudeDispersion` owner supplies phase delay and GDD
where an amplitude is nonzero. Raw jets remain available at a Brewster zero;
an undefined logarithmic phase is never replaced by zero.

Both frequency constraints require the original incident-law jet for boundary
admittance, including fixed dimensional transverse wavenumber. The existing
normal-only query retains its broader admission. An already prepared
`UniaxialFrequencyDispersion` can add the boundary response without another
native Maxwell/root solve, using its exact retained incidence coordinate.

The analytic fixed-size DD contractions use the original independent electric
and magnetic transmission numerators, the outgoing gradient root, and the
source-bound principal spectra. This retains tilted negative/zero normal phase,
matched grazing and representable unmatched grazing derivatives. Reflected and
transmitted normal-power slopes retain both field jets; TIR keeps exact constant
normal powers while preserving its complex phase dispersion. Critical roots,
unsupported source jets and unrepresentable derivatives refuse. These local
boundary slopes do not provide finite-crystal interference, pulse reshaping or
a single delay for an interfering pair of modes.

Seven new numerical fixtures and three source-owner fixtures cover original
isotropic parity, independent derivatives, both frequency constraints, magnetic
contrast, zeros, TIR, source endpoints and grazing limits. They were parsed and
independently source/math reviewed. They were not executed because the bounded
native-build capacity check failed; the earlier 66-test native execution covers
the previously published carrier and normal-propagation kernels only.

### Arbitrarily oriented transparent crystals

`uniaxial::oriented::OrientedUniaxialInterface` retains the same original
incident law and two-axis `UniaxialOpticalResponse`, and admits an authored
Cartesian optic axis with any azimuth. Each frequency query resolves both
original spectral laws and permeability before solving. Caller-supplied complex
tangential [Ex,Ey] fields retain absolute V/m, elliptical polarization and phase.
The ordinary and extraordinary eigenmodes are solved together with both
Cartesian tangential E and H boundary conditions. An out-of-plane axis can
therefore produce reflected cross-polarization and split the transmitted
energy into two differently directed rays.

The response contains full SI electric and magnetic boundary vectors, normal
wavenumbers, independently calculated Poynting power, and separate phase and
energy directions for each transmitted eigenmode. The conserved signed
`kx/k0` query admits exact critical inputs and their adjacent representable
neighbors. Either mode may become evanescent; its decaying field is retained
without publishing a fictitious bulk ray. Exact isotropic and optic-axis
alignment limits use a deterministic nonsingular polarization basis.

This fixed-size, allocation-free solve eliminates reflection analytically to
leave a complex 2x2 transmission system. Power-of-two index scaling and
factored double-double dispersion discriminants protect critical classification.
Normalized double-double axis weights precede the positive index-squared sums;
this retains exact symmetric critical coordinates without a tolerance window
and preserves small principal values at large contrast. Both signed adjacent
representable tangential inputs retain their propagating/critical/evanescent
classification.
Tests reconstruct both Maxwell curls from the dimensional constitutive tensor,
verify Cartesian boundary continuity and independent normal power, and compare
principal sections against the existing decoupled owner. This extension remains
one lossless planar half-space; it does not claim general anisotropic film,
beam, absorbing-medium, geometric traversal or broadband transient support.

`OrientedUniaxialInterface::group_velocity_response` and its signed
tangential-index counterpart pair that unchanged boundary solution with
`OrientedCrystalGroupResponse`: original principal group indices and separate
ordinary/extraordinary group-velocity vectors in Cartesian [x,y,z] m/s.
`PreparedOrientedCrystalField::group_velocity_response` adds the same original
spectral derivatives to an already solved field without another boundary/root
solve. Both paths allocate nothing and retain their original source and axis.

The full anisotropic dispersion supplies the modal derivatives, including
nonzero y walk-off despite conserved ky=0. The two crystal paths share a
double-double signed spectral contraction; transverse cross-product weights
avoid subtracting nearly equal squared norms. Exactly unexcited and evanescent
modes return no bulk group vector. Critical modes retain exact zero normal
group velocity, including tilted modes with negative normal phase. Anomalous
signs remain signed, and a zero or unrepresentable denominator retains the
same explicit boundary as the principal-plane owner.

At exact equal principal phase indices with unequal principal group slopes,
an out-of-plane optic axis requires a spectral eigenbasis that the original
deterministic TE/TM carrier does not supply. The opt-in group query refuses
that interpretation while monochromatic fields remain usable. Principal-plane
crossings and equal-group isotropic limits remain admitted. Separate modal
velocities are not averaged into a velocity for a coherent superposition.
Focused regressions compare independent dispersive Maxwell energy and flux,
principal-section reduction, axis reversal, original-source/Prepared parity,
near-cancelling slopes, spectral crossings and critical/evanescent mode limits.

`OrientedUniaxialInterface::prepare_field` and
`prepare_field_at_tangential_index` retain the original interface, admitted
frequency and solved boundary in `PreparedOrientedCrystalField`. Repeated
Cartesian point queries allocate nothing and repeat neither source resolution
nor the boundary solve. Incident, reflected, total incident-side, ordinary,
extraordinary and total transmitted selections enforce their physical
half-spaces; both boundary traces remain available at z=0.

Complex peak E/H retain the original exp(+i omega t - i k.r) convention.
The infinite plane wave has ky=0 even when energy walks off in y. Both
transmitted modes interfere before the full signed cycle-mean Poynting vector
is evaluated. Separate compensated carrier and relative phases retain optical
path products and birefringent retardance; the shared decay owner preserves
representable late evanescent tails. Decaying fields create no heat in this
lossless law. Spatial Maxwell curls, Fresnel standing waves, polarization
retardance, cross-mode energy transport and critical/evanescent fields are
covered by focused tests. The prepared point evaluator and the existing
`spectrum::PreparedUniaxialAngularSpectrum` share one compensated field
accumulator, retaining relative phase between both eigenmodes and authored
incident directions before coherent power contraction. Sparse spectra keep
their original source resolution, global-axis rotations, coefficients,
mode/workspace limits and cancellation. Finite-slab observers pass their solved
waves, signed normal wavevectors and physical face anchors into this same
propagation owner. The original face-distance subtraction retains its DD
residual before phase and attenuation multiplication; relative carriers are
combined before phasors round. Equal wavevectors anchored at distinct faces
therefore retain their relative phase even at very large observation
coordinates. Workspace admission accounts for the actual shared accumulator.
FFT-grid attachment remains outside this change.

`stack::UniaxialStackResponse::field_at` and
`stack::PreparedUniaxialStackSpectrum::field_at` use this same spatial
propagator for the original coupled multilayer solution. Internal observations
share a compensated layer-local depth across every incident direction, with
forward waves anchored at zero and backward waves at the authored thickness.
Both physical faces retain that exact authored thickness even when cumulative
global face coordinates have rounded. Exterior observations retain global
entrance and final-exit anchors. Relative phase, weak interference, signed
Poynting flux and scaled evanescent tails are therefore evaluated before
individual propagated fields round. The original layer selection, empty-mode
region admission, source preparation, boundary solve, cancellation and
workspace limits remain with their existing owners. Point queries perform no
allocation or material/boundary reevaluation.


### Arbitrary-axis normal-transfer frequency dispersion

`OrientedUniaxialInterface::frequency_dispersion` and its signed
tangential-index counterpart return `OrientedCrystalFrequencyDispersion`.
It retains the original interface and unchanged vector boundary solution,
with separate ordinary and extraordinary complex first/second normal
wavenumber derivatives. Both explicit frequency constraints use the same
source-bound material jets and analytic kernel as the principal-plane owner.
Mode-selected `transmitted_propagation` supplies the existing signed
`BulkPropagation` transfer, normal delay/GDD and evanescent attenuation
derivatives. Modal delays are not combined into one delay for interfering fields.

`PreparedOrientedCrystalField::frequency_dispersion` supplies the same
derivatives without solving the boundary or its roots again. Both original
prepared-field constructors retain their exact private incidence input;
a grazing cosine is never reconstructed from a rounded transverse wavevector.
The existing prepared field and monochromatic carrier remain unchanged.
Exact local isotropy through both derivative orders reduces to the original
scalar-root identity, preventing normalized-axis residuals from shifting a
grazing cutoff. Unequal original spectral jets retain separate branches.

Zero excitation does not remove a mathematically defined eigenmode derivative.
Evanescent modes retain their complex outgoing decay branch. Either exactly
critical branch refuses the complete two-mode derivative query. At equal
principal phase indices with an out-of-plane axis, both first and second
original principal jets must coincide: otherwise the native deterministic
TE/TM basis is not a spectral eigenbasis. Principal-plane crossings and
equal-local-law isotropic limits remain admitted.

The two additional native arbitrary-axis tests bring the bounded numerical
build to 66 passing tests, using the original `fs_math` and native vector
carrier. They compare independent inverse-metric derivatives and native
frequency perturbations, both frequency constraints, axis reversal, distinct
principal spectral slopes/curvatures and matched grazing incidence.
Original-source fixtures retain carrier/source identity, mode-specific
propagation, prepared/direct parity, source-domain and spectral-crossing
refusals. These source-bound fixtures passed parser/API review but were not
executed in the full material dependency build. This remains local modal
propagation. The separate vector-boundary owner below supplies interface
amplitude derivatives; finite slab/stack transfer derivatives remain separate.

### Arbitrary-axis vector boundary frequency dispersion

`OrientedUniaxialInterface::boundary_frequency_dispersion` and its signed
tangential-index counterpart provide source-bound first and second angular
frequency derivatives of the complete incident, reflected, ordinary and
extraordinary Cartesian E/H fields. `PreparedOrientedCrystalField` adds the same
response without another native boundary/root solve. The unchanged carrier,
exact original incidence, normalized eigenbasis and solved coefficients are
retained in a fixed-size preparation witness; the ordinary steady solve captures
no witness.

Both explicit incidence-frequency constraints require original two-sided
lossless incident and principal material jets for boundary admittance. The
normal-only owner keeps its broader fixed-transverse-wavenumber admission.
Critical branches, unsupported source endpoints/cusps, unrepresentable results
and unresolved out-of-plane spectral degeneracy refuse. No neighboring
frequency, fitted law, frequency-dependent normalization switch or per-query
allocation enters the derivative calculation.

The native 2-by-2 Maxwell boundary system and vector eigenfields are
differentiated analytically in a frozen local gauge. Double-double normalized
wave-coordinate derivatives retain weak material slopes that would disappear
when subtracting bulk wavenumber terms. Exact local index-scaling, matched
isotropic and decoupled principal-plane identities protect admitted grazing
limits without a tolerance cutoff. Shared axial cosine factors cancel before
transverse material derivatives, including mixed incident polarization. The original steady carrier arithmetic remains authoritative.

`OrientedCrystalVectorFrequencyJet` exposes raw E/H derivatives and optional
per-Cartesian-component logarithmic phase/delay/GDD. Exact component zeros keep
their raw jets, including nonzero derivatives, and have no invented phase.
The transmitted boundary query coherently sums both modal jets; separate normal
propagation is still explicit. No single modal-average delay is assigned to an
interfering vector field.

The prepared boundary owner also exposes `field_frequency_dispersion(side,
point_m)` for all six original selections: incident, reflected, their coherent
sum, either transmitted eigenmode, and both transmitted modes. Its
`carrier_field()` returns the original Cartesian E/H and coherent Poynting
exactly; `vector_frequency_jet()` supplies complete local first/second Cartesian
field derivatives. The original half-space, source and phase conventions remain.

Preparation retains original incident bulk derivatives and private DD normal
jets. Public rounded normal outputs are unchanged. Fixed-angle point queries
follow the original signed incident transverse derivative; fixed dimensional
kx has exactly zero transverse derivatives. Reflected normal derivatives change
sign. Repeated point queries perform no material evaluation, root/boundary
solve or allocation.

The native propagation owner factors common complex phase and its derivatives
before coherent interference. DD modal derivative differences are formed before
spatial products; the least-decaying active wave supplies the fixed reference.
All three amplitude orders determine activity, so zero carrier excitation does
not discard a derivative wave. Weighted frequency polynomials precede common
attenuation, reusing native weak-decay and late-tail helpers. Nonfinite
pre-decay polynomials and relative underflow of a lower-order field still needed
by a common-phase product refuse explicitly. The incident curvature retains
the original bulk owner's physical K'' representation.

Two additional source fixtures cover all six selections, both constraints,
signed propagating and evanescent cases, independent original-frequency
evaluations, exact carrier/Poynting parity and physical/scaled derivative units.
They also cover derivative-only tails after carrier extinction, weighted late
tails, zero reflected carriers with nonzero derivative fields, and original
domain/half-space refusal. Changed Rust passed parsing and independent
source/API/math review; these two fixtures were not executed here. This adds
local field derivatives, not a broadband pulse, slab-transfer or
Poynting-derivative model.

Eight numerical controls and four original-source fixtures cover independent
derivatives under both constraints, Maxwell continuity, weak slopes, zero modal
excitation, native carrier identity, prepared/direct parity, source support,
critical/degenerate refusal and grazing limits. The changed Rust was parsed
and independently source/API reviewed. These new numerical and integration
fixtures were not executed here because the available native-build storage
remained below the admitted budget; earlier numerical results are separate.

### Retained vector polarization on the native planar Gauss solve

The default-off `dielectric-hysteresis` feature exposes
`electrostatic::nonlinear::hysteretic::PlanarHystereticElectrostatics` and
`FrozenHystereticDielectric`. Explicit original per-cell fits reuse the native
coupled vector play/spherical saturation law. The existing planar P1 geometry,
electrode reduction, mixed voltage/total-charge controls, matrix-free Newton
tangent, checked CG and backtracking remain the sole field solver. Prepared
geometry and electrode layout are shared across accepted increments.

Every Newton trial freezes the SAME accepted polarization. The voltage
objective is `Phi(E;P0) = E.D - W - Q`, with gradient `D` and native symmetric
positive tangent `dD/dE`. Stable pairwise objective differences retain that
history and avoid subtracting large absolute field energies. Physical stored
energy, per-cell switching heat and endpoint numerical loss remain separate.
The source-work check independently pairs final electrode potentials with
charge changes and final nodal potentials with actual volume-charge changes.
Charge, field work, source work and physical energy gates all precede publication
of the complete new field and memory; cancellation, work budgets and source
domain refusals leave the accepted history replayable. Relative voltages retain
small fields beneath large common electrode biases.

Zero switching branches give ordinary linear cells. The whole-field owner
starts from a compatible virgin zero-field state; poled specimens are prepared
through accepted increments. The lower-level frozen cell adapter also admits
existing native memories but leaves prior global consistency with its caller.
Out-of-plane memory refuses planar binding. Reversible-only force shortcuts and
smooth differential-capacitance extraction refuse hysteretic cells unless their
original admission permits a separate reversible exterior.

This is fixed-geometry, fixed-fit isotropic quasistatics. Finite increments use
endpoint electrical work; the nonnegative numerical remainder is never heat.
Resolve reversals and rotations and refine the path. This lane supplies no time
kinetics, spatial thermal deposition, changing switching coefficients, leakage,
breakdown, switching strain, measured calibration or continuum error estimate.
`tests/hysteretic_electrostatic.rs` compares analytical linear, volume-charge,
floating-conductor and scalar-play accounts, rotated geometry and immutable
refusal/retry behavior through the actual native field owner.

`PlanarHystereticElectrostatics::at_source_queries` refreshes one original query
per physical triangle while retaining exact field, charge, geometry, storage
and polarization history. `HystereticElectrostaticCandidate::with_source_queries`
performs that refresh on an owned unpublished proposal without cloning its
physical field arrays or repeating Gauss. Every fit must remain exactly the
same constant fit through the original source resolver. A short query batch,
source-domain failure or cancellation cannot partly publish changed sources.
The coupled thermal consumer still owns caloric identity, temperature intervals,
receiving volumes and parameter-work admission; these methods do not infer a
temperature-dependent switching law. Two additional spatial regressions retain
poled memory through source refresh and reversal, and exercise refused partial
batches with deterministic retry.

`hysteretic::path` adds bounded continuation along one explicitly authored
affine electrode/volume-charge segment. Every comparison uses the native field
owner: one full increment versus two sequential halves, with only the latter
eligible for accepted history. Dyadic rejection retains a bounded interval
stack and reuses the left coarse trial. One global field-solve/operator budget
counts all attempts, including discarded trials. Separate indicators monitor
branch polarization, the compensated sum of absolute per-cell heat defects
(opposing spatial heating errors cannot cancel), and the larger of the
source-work defect and retained numerical loss. Absolute allowances scale with
interval fraction; accepted indicators and native physical cell accounts sum
over the whole path. Numerical loss remains separate from physical heating.

Relative voltages retain common-bias invariance. Conversion to the original
electrostatic candidate moves cumulative accounts and a compact path summary
to the finite-enthalpy consumer without another solve or history reset. These
indicators are empirical refinement estimates, not certified global errors,
time kinetics or automatic discovery of waveform turns. Supply explicit cuts
at turns and source/control events. The native nonlinear solver may also
accept an Armijo-rejected trial when its actual Gauss residual already meets
the ORIGINAL target; final independent charge and energy checks remain active.
`tests/hysteretic_electrostatic_path.rs` exercises linear work convergence,
charge/bias invariance, exact major-loop heating, changing-direction memory
against an independently subdivided path, and atomic global-budget refusals.

### Dispersive Jones observations on original spectral lines

`AngularSpectrumExposureApertureJonesSpectrum` borrows a validated, strictly
ordered table of `(angular_frequency_rad_s, JonesProjection)` pairs. The four
`measure_{illumination_}aperture_spectral_projection` queries on exposure and
ensemble owners apply the response at each original grouped frequency before
quadratic products. They observe `Re(sum_j a(omega_j)† E_j exp(i omega_j t))`
and the corresponding magnetic signal. Complex conjugation, absolute spectral
phase, real-source carrier sums, proper-factor beats, both physical faces and
the original source clock remain explicit. Frequency-dependent retardance and
common delay therefore change the actual finite-shutter observation.

`AngularSpectrumExposure::angular_frequencies_rad_s` exposes the original sorted
group keys. The table must cover every retained frequency in the mean and all
factors, including equal-frequency groups that combine multiple authored lines.
Coverage is checked before spatial caches or Maxwell reconstruction. Matching is
exact; unused entries are allowed, while missing frequencies refuse. There is
no inferred phase interpolation, extrapolation, material continuation or
frequency-to-wavelength round-trip matching. The caller retains the actual
calibration/material response; finite samples establish no broadband passivity,
causality, absorption or back-action model.

Each response acts on all spatial modes in the original Cartesian basis;
angle-dependent optical scattering remains with the physical field owners.
Admission and coverage are cancellable and allocate nothing. Coverage and
response lookups cost O(log S) per group for S samples, never per spatial mode.
The same double-double accumulator still forms three scalar products per
ordered mode/frequency pair. Existing source/query bounds, reconstruction,
spatial/shutter kernels and factor statistics remain authoritative. Borrowed
table storage belongs to the caller. Constant and real-weight queries retain
their original arithmetic and require no spectral table or lookup.

Focused regressions cover a chromatic quarter-/half-wave analyzer, common-delay
equivalence to shifting the authored shutter, independent spatial/time/phase
quadrature, exact support refusals and cancellation, and bitwise equality for
constant tables. The original-source integration regression exercises all four
query variants through the Maxwell frequency banks and both illumination faces;
numerical-kernel probes alone do not execute that reconstruction.

## Source-bound Lorentz resonance in coherent waves

`plane_wave::DielectricResponse::Lorentz` consumes the original inertial
polarization bank from `fs-material`. Interfaces, finite layers and their
coherent wave consumers retain its actual complex permittivity below, at and
above a damped resonance, including negative real values. Existing outgoing
and evanescent branches, TE/TM conventions and Poynting-flux accounting apply.

`LorentzWithOhmic` adds only the independent `-i sigma/omega` contribution under
`exp(+i omega t)`. Both laws must retain the identical complete original
material snapshot; equal coefficients from different sources do not qualify.
A true undamped constitutive pole refuses instead of acquiring artificial loss.

Three source-level integration regressions specify passive interface/film
scattering through resonance, homogeneous-film subdivision, independent Ohmic
loss, complete-source identity refusal, and lossless negative-permittivity
reflection/tunneling without heat. Their real material-card/quantity fixtures
compile; the full EM tests are authored but unexecuted in the incomplete
constellation dependency environment. Fifteen native production oscillator
tests establish the underlying constitutive dynamics, not full EM integration.

### Finite-window Lorentz heat and reversible storage

The original prepared exposure-layer owner now admits Lorentz and
LorentzWithOhmic media. Each frequency caches the original restoring coordinate
`P/sqrt(delta_epsilon)`, kinetic coordinate `Pdot/(omega0*sqrt(delta_epsilon))`
and friction channel `sqrt(gamma)*Pdot/(omega0*sqrt(delta_epsilon))`. The shared
finite-window products contract complete coherent fields before averaging.
Both restoring and kinetic coordinates contribute their reversible storage
rates; only squared real friction channels become heat. Independent Ohmic heat
enters once. Real and explicitly proper source statistics retain their original
sum/difference-frequency semantics, clock and detector weighting.

`AngularSpectrumExposureHeatChannel::LorentzPole` selects one original pole.
Polarization and Total channels flow through the existing local heat spectrum
and spatially integrated thermal-load owners. Debye and Lorentz selectors never
alias one another. Original material identity, physical volume/support, depth
quadrature, cancellation and spectrum limits remain active. Transfer evaluation
occurs once per source frequency; repeated depths reuse the film/transfer caches.
Pair budgets charge `N*sum(F_r^2)*(3+D+2*L)`, plus existing carrier work, and
workspace admission includes the extra kinetic scratch and retained transfers.
The Debye path does not allocate or contract kinetic scratch.

Three kernel regressions specify original P/current energy balance, proper
slow beats under a large carrier, and lossless/cancellation behavior. Two
original-film integration regressions specify subcycle real-time oracles,
uniform and triangular windows, instantaneous heat spectra, physical volume
loads, zero-loss negative-permittivity storage and the live kinetic work cap.
These five EM regressions are authored but unexecuted without the complete
constellation dependency closure; formatting and bounded source/API review
passed. The separately executed native material tests do not qualify the full
film reconstruction.

This observes settled prescribed harmonics over finite windows. It does not
add switch-on Maxwell transients, temperature-dependent parameter work or
material-frequency derivatives. Steady Poynting absorption is not substituted
for finite-window storage/work.

### Coherent angular spectra in transparent crystals

`uniaxial::oriented::spectrum::PreparedUniaxialAngularSpectrum` binds a sparse
authored incident spectrum to the original `OrientedUniaxialInterface` at one
positive supported frequency. Preparation resolves the original incident law,
both original principal-index laws and permeability once. Each mode supplies
a signed tangential phase index, an incidence-plane azimuth and global complex
Ex/Ey peak coefficients. The fixed global crystal axis is transformed into
each incidence plane before the existing coupled Maxwell solve; the physical
crystal is never rotated independently for different incident directions.
Exact zero coefficients may be omitted without discarding any nonzero mode.

Physical-point queries return all Cartesian E/H components and the full local
cycle-mean Poynting vector. The ordinary and extraordinary waves retain their
own normal wavenumbers, including negative normal phase and evanescent decay.
The original double-double longitudinal phase propagation is reused; transverse
phase products and coherent field sums also retain double-double residuals.
Both eigenmodes and all source directions are summed before computing flux,
preserving cross-polarized reflection, depth-dependent elliptical polarization,
spatial birefringence, interference fringes and incident-side standing waves.
Fields follow Re(E exp(+i omega t)), with z=0 at the interface and crystal z>=0.

Mode/work and numerical-storage caps precede allocations and material solves;
cancellation prevents a partial prepared source or field from being returned.
Preparation is O(N), and each point query uses at most two plane-wave fields
per retained mode with no allocation, material reevaluation or boundary solve.
This sparse monochromatic owner introduces no sampling grid, FFT, paraxial
approximation, spectral interpolation, power normalization or inferred source
coherence. Finite crystals, absorbing anisotropy, films, transient reservoirs
and geometric ray traversal remain outside this half-space capability.

### Finite arbitrarily oriented crystal slabs and coherent beams

`uniaxial::oriented::slab::OrientedUniaxialSlab` binds the original
`OrientedUniaxialInterface`, an original lossless isotropic exit response and a
positive physical thickness. Each frequency query retains both original
principal-index laws, all three original permeabilities and the fixed crystal
axis. The existing arbitrary-axis eigensystem supplies ordinary/extraordinary
waves; exact normal-coordinate reversal supplies the backward family with
polar E and axial H transformation. A tilted extraordinary backward normal is
the other dispersion root, not an assumed negative of its forward normal.

Both faces are solved together through one fixed-size complex 8x8 tangential
E/H system. Forward waves are anchored at the entrance and backward waves at
the exit, so every evanescent propagation factor decays. Coherent internal
reflections retain etalon resonances, polarization conversion and frustrated
total-reflection tunneling. Results expose the four internal vector fields,
the actual exit-plane transmitted field, independently evaluated reflected and
transmitted normal powers, and allocation-free physical-point E/H and Poynting
queries in the incident half-space, crystal interior or exit half-space.
Row-scaled pivoting, coherent boundary-continuity checks and an independent
normal-power balance gate precede publication of the field snapshot.

`PreparedUniaxialSlabSpectrum` reuses the authored sparse incident modes,
resource limits and double-double field accumulator of the half-space owner.
Original source laws are resolved once, each occupied direction receives one
two-face solve, and all four internal waves and all source directions are
summed before forming Poynting flux. The same global crystal axis applies to
every azimuth. Preparation is bounded O(N); point queries allocate nothing and
perform no material or boundary reevaluation. Cancellation returns no partial
source or field.

The shared spatial propagator keeps entrance anchors for incident, reflected
and forward internal waves, and exit anchors for backward internal and
transmitted waves. Backward waves retain their actual signed global normal
wavevector. All physical face distances, relative eigenmode phases and
interdirection phases enter one compensated E/H sum before signed SI Poynting
flux is formed. Scale-aware attenuation preserves representable late
evanescent fields when a standalone exponential would underflow. Physical
region admission occurs once per query; preparation and point queries retain
their original cancellation and storage limits. Focused numerical controls
cover optical face-distance residuals, weak backward retardance, distinct-face
interference, late exit evanescence and sub-ULP interdirection cancellation.

This is one homogeneous transparent uniaxial slab. The exit can propagate,
be critical or be evanescent. Exactly coalesced **internal crystal** critical
modes require a generalized Jordan basis and explicitly refuse, as do
numerically unresolved systems or unrepresentable positive propagating output
power. Existing half-space critical support remains unchanged. Absorbing
anisotropy, crystal multilayer stacks, varying optic axes and causal transient
reservoirs remain outside this capability.

### Analytic thickness sensitivity of coherent crystal slabs

`OrientedUniaxialSlab::thickness_derivatives` and its signed-tangential-index
counterpart differentiate the original two-face Maxwell system at fixed
frequency, material laws, optic axis, incidence and complex incident tangential
E. The result retains the unchanged native carrier and first/second ordinary
derivatives of Cartesian reflected, transmitted and all four internal E/H
fields, together with normal reflected/transmitted power derivatives. The
entrance stays at z=0; transmitted and backward amplitudes follow the moving
exit z=d. These are not fixed-global-observer derivatives.

The original row-scaled pivoted 8x8 factorization serves the carrier and both
derivative right-hand sides, A c'=-A'c and A c''=-A''c-2A'c'. Only the original
bounded forward/backward propagation factors vary. Dimensionless thickness
coordinates keep wave-number scales outside the boundary solves; DD
contractions retain coherent vector products and normal power derivatives
before SI scaling. There are no neighboring material queries, finite-difference
production solves or heap allocations. A reflected carrier zero can retain a
nonzero complex field slope and positive power curvature. Critical/evanescent
isotropic exit normal power remains identically zero.

Fifteen exact native slab numerical tests passed, including four new thickness
tests covering independent five-point carrier controls, homogeneous
moving-exit phase, tilted-axis polarization coupling, signed incidence,
frustrated total reflection, resonance curvature, opaque but representable
tunneling, amplitude scaling and physical units. The additional original-source
spectral-endpoint fixture is authored but unexecuted; full workspace/material
resolver validation remains outstanding. Original internal critical-mode,
boundary-conditioning and carrier-power gates remain authoritative.
Underflowed propagation factors and unrepresentable derivative arithmetic
explicitly refuse instead of fabricating zero sensitivity. This capability
does not differentiate material frequency, optic-axis orientation, multilayer
stack parameters or fixed global field points.

### Coherent stacks of arbitrarily oriented transparent crystals

`uniaxial::oriented::stack::OrientedUniaxialStack` retains the original
lossless isotropic entrance and exit responses and an ordered borrowed slice
of `UniaxialStackLayer`. Each layer binds its original two-axis
`UniaxialOpticalResponse`, same-state scalar permeability, positive physical
thickness and independently authored fixed global optic axis. All source laws
are resolved at the requested positive frequency; their original support and
transparency refusals remain authoritative.

The existing dimensional forward/backward Maxwell basis is shared with the
single-slab owner, including polar-E/axial-H normal-coordinate reversal for
tilted axes. A right-to-left fixed 4x4 Schur reflection sweep eliminates each
remaining load, preserving the noncommuting backward-phase/reflection/
forward-phase order. A forward amplitude sweep then recovers all four
internal eigenmode fields. Every evanescent mode remains anchored at its
decaying face; no inverse exponentially decaying transfer factor is formed.
Preparation and numerical storage scale linearly with the layer count.
The original row-scaled pivoting owner now accepts fixed array sizes and
retains its existing size-eight single-slab arithmetic.

Results retain full reflected and transmitted vector fields, original
physical face positions, four internal waves per layer, independently
evaluated R/T and physical-point E/H/Poynting queries. Every tangential E/H
interface trace and each layer's normal flux must agree before a snapshot is
returned. The propagating incident port and any positive propagating output
power must be representable. Exactly coalesced internal critical modes,
unresolved boundary systems, failed flux/continuity gates and unsupported
source data refuse explicitly. The exit itself may be critical or evanescent.
Global coordinates identify the layer; either adjacent trace is available at
a shared face.

`PreparedUniaxialStackSpectrum` resolves all original source laws once and
solves each nonzero authored `UniaxialIncidentMode` against every original
global layer axis. All four waves in the selected layer and all incident
directions interfere before the common compensated Poynting contraction.
Preparation is O(modes times layers); each nonempty spectrum point query is
O(modes), with no allocation, material reevaluation or boundary solve.
Exactly zero source coefficients may be omitted; nonzero coefficients are
never filtered or renormalized. Layer, mode and logical workspace bounds
precede allocations/solves, and cancellation returns no partial result.

This is a monochromatic stack of parallel, homogeneous, transparent uniaxial
layers. It supplies no absorbing anisotropy, biaxial crystal, rough interface,
continuously varying axis, inferred coherence, transient reservoir or
continuum-error certificate. Large cumulative face coordinates must retain
distinct representable boundaries.

Nine focused numerical regression groups were authored for independent
isotropic TE/TM characteristic matrices with magnetic contrast, noncommuting
resonant waveplates, tilted-layer subdivision and axis-sign invariance,
evanescent tunneling, many-layer Bragg reflection, critical/refusal behavior,
global rotation, coherent fringes, resource bounds and cancellation.
An additional native source/material integration case covers original layer
retention, spectral support, single-slab agreement and every coherent
interface. These new checks have not been executed: the execution environment
became unavailable before the first actual-source probe returned. Earlier
single-slab results do not establish the new stack or its refactored helpers.
## Total-current electrodes in the native moving conductor

`SpatialResistance::field_at_current_with_motion`,
`SpatialConductanceMatrix::field_at_currents_with_motion` and
`SpatialMixedConductance::field_at_controls_with_motion` share their original
electrode admission, conductivity/contact operator and equilibrated terminal
factors. One additional zero-voltage motional solve supplies the affine current
source; prescribed currents close against `G V + I_motion`, including floating
electrodes and a current-controlled reference terminal. Preparation plus motion
costs two native field solves for a pair, or N for N terminals. Caller work
budgets must include that additional solve.

Only native potentials are superposed. The actual field then recovers local
current, bulk/contact Joule heat and Lorentz work together, with independent
terminal-current and signed authored `sum(V I)` checks. Open-circuit EMF,
short-circuit braking and generated electrical power remain distinct from
positive conductor-frame heat; probe/source heat is never deposited. Native
moving-contact restrictions, cancellation and solver limits remain in force.
These are fixed-domain quasistatic snapshots, without induction, sliding-contact
physics, mechanical integration or a regulated-power branch. Focused manufactured
bar regressions cover signed current, mixed floating electrodes, reordered
terminals, contact heat and transactional refusal/reuse.

## Self-consistent inertial polarization on original potential fields

`electrostatic::dispersion::lorentz::{PlanarLorentzDielectric,
AxisymmetricLorentzDielectric}` and
`spatial::electrostatic::dispersion::lorentz::SpatialLorentzDielectric` retain
both polarization P and current J for every original source-bound Lorentz pole.
The same triangle/tetrahedron field kernel solves Gauss/current balance with
`D = epsilon_infinity E + sum P`. Current is an independent inertial state,
not an additional displacement or inferred Ohmic conductivity.

The original material's passive midpoint law condenses each P/J pair locally.
For h=dt/2 its endpoint permittivity contribution is
`delta_epsilon (h omega0)^2 / (1 + h gamma + (h omega0)^2) >= 0`.
This positive time-step coefficient is not a positive-real frequency-domain
permittivity assumption: ringing, field-feedback frequency shifts and finite
lossless resonant growth remain in the retained state. Only original potential
DOFs enter the global solve. Fixed-size steps cache both the material coefficient
bank and original Hodge operator; warm steps do not reprepare pole coefficients.
The constitutive adapter delegates the equations to the material owner and
monomorphizes Debye states without allocating unused Lorentz current vectors.

Physical storage includes instantaneous electric, restoring polarization and
kinetic polarization energy. Independent pole friction supplies nonnegative
heat; no energy residual or kinetic store becomes heat. Original cell volumes,
revolved metrics, electrode traces, free nodal charge, coherent voltage work,
physical residual budgets and cancellation remain owned by the existing field
kernel. Prepared objects require the same immutable original field owner.
Failed candidates publish neither an endpoint nor heat. The common
`DielectricFieldStep` account is the existing dielectric-memory account type.

The new constructors are explicitly insulating, with source-fixed isotropic
coefficients and geometry. They do not add magnetic induction, propagating
Maxwell transients, carrier transport, thermal parameter work or changing
material occupancy. Midpoint is second order and passive but not L-stable;
energy balance does not replace temporal or spatial convergence.

Six focused integration regressions specify held-field second-order dynamics,
isolated-charge longitudinal resonance and voltage reversal, lossless resonant
secular growth, oblique tetrahedral kinetic storage/charge, original revolved
volume accounts, and empty-bank/owner/cancellation boundaries. Two receiving
thermal regressions check physical heat deposition and source/refusal behavior.
These full-crate tests are authored but unexecuted without the complete native
dependency closure. Formatting and bounded source/API/physics review passed;
material numerical tests are reported separately, not as field integration.

### Conducting Lorentz fields and original total-current balance

The planar, axisymmetric and spatial Lorentz owners now expose `new_leaky`
under the existing `resolved-leaky-dielectric` feature. Each cell consumes one
`ResolvedLeakyLorentzMaterial`, preserving the original complete Lorentz source
and independently selected scalar DC conductivity. A private source-extraction
adapter admits both Debye and Lorentz laws into the existing field constructor;
all geometry, Hodge assembly, potential solves, charge recovery and work/loss
accounting remain shared. The insulating constructor retains its original scope.

The field solves the original weak total-current equation. At an unprescribed
node, Delta free charge plus the integrated Ohmic outflow is zero. Lorentz
polarization/current recovery retains the original pole equations. Local
elimination adds `dt*sigma/2` to the positive algorithmic endpoint coefficient
and the corresponding old-field conduction history. This coefficient is never
reported as stored permittivity. Physical storage remains epsilon-infinity
field energy plus each pole's restoring and kinetic energy. Independent
`dt*sigma*|E_mid|^2` Joule heat and original pole friction are reported separately,
as are displacement and conduction contributions to electrode work. An exact
source-backed zero conductivity follows the same numerical path as insulation.

Prepared coefficient banks and native field operators remain reusable at fixed
dt and bound to the original geometry/material owner. Conductivity is visible
through the same original property-name/value getters used by thermal receivers.
The scope is fixed-geometry isotropic electroquasistatics with static DC Ohmic
conduction: no carrier inertia, anisotropic transport, breakdown, temperature
parameter work or propagating Maxwell wave claim is added. Midpoint is passive
and second order; stiff leakage need not be monotone, and resonant phase still
requires time refinement.

Four G1/G3/G4 field regressions are authored: isolated three-state charge/
polarization/current dynamics against independent continuous roots on three
time refinements; equilibrium DC current with unchanged resonant storage;
revolved and oblique tetrahedral physical-volume/current accounting; and exact
zero-conductivity equivalence plus owner/cancellation refusal. Two receiving
thermal tests cover the joint loss and source-validity path in fs-couple. These
full-field tests remain unexecuted while the complete workspace dependency
closure is unavailable; formatting and focused source/API review are the
executed checks for this adapter extension.

### Source-bound anisotropic absorption and local optical heating

`OrientedUniaxialInterface::prepare_absorbing_field` and its signed tangential-
index counterpart admit the original ordered complex principal indices of the
same `UniaxialOpticalResponse`. Both source spectra, complete material snapshot,
original wavelength support, interpolation-gap limit and same-state real scalar
permeability remain authoritative. The source now exposes ordered complex
principal indices and SI permittivities without resampling or discarding loss.
Existing transparent queries retain their original extinction refusals.

The native uniaxial Maxwell eigensystem supplies complex ordinary and
extraordinary waves for the original fixed arbitrary optic axis. Existing
incident-coordinate admission and Maxwell eigenvectors are shared with the
transparent owner. The extraordinary quadratic retains its stable paired roots
until decay or normal-flux admission selects the outgoing branch. Complex
Maxwell residuals, tangential continuity and independent entrance R/T balance
must pass. The existing row-scaled fixed boundary solve and compensated
face-anchored field propagator supply the solved amplitudes and spatial fields.
Exactly zero source extinction delegates to the original transparent solve.

Prepared results retain full complex reflected fields and both transmitted
modes. Point queries return coherent Cartesian E/H, signed cycle-mean Poynting
flux and local dielectric heat `omega/2 E*[-Im(epsilon)]E` in W/m³. The two
eigenmodes interfere before either power contraction; adding their separate
heats or powers is generally incorrect. Positive weighted parallel/perpendicular
squares preserve anisotropic dissipation without subtracting nearly equal
intensities. Pure evanescence in a lossless eigenaxis supplies no heat.
Transmitted interface power means power entering the crystal at z=0, before
bulk absorption. Preparation and each point query have fixed two-mode work and
allocate nothing; point queries neither resolve material sources nor repeat
boundary solves.

`PreparedAbsorbingCrystalField::slab_absorption` integrates this same original
dielectric work over an authored positive face area and finite interval
`0 <= start_depth_m < end_depth_m`. The returned
`AbsorbingCrystalSlabAbsorption` borrows the prepared illumination and its
original interface/material, retaining the area, both depths, physical volume
and cycle-mean absorbed watts. It introduces no exit interface. The original
loss tensor weights coherent ordinary/extraordinary cross terms as well as
modal self terms; a sum of separately attenuated modal powers is insufficient.
Fixed analytic exponential integrals use relative modal phase and compensated
products. Thin intervals contract field Taylor coefficients before squaring
so a coherent node retains its cubic heat, and weak loss never depends on
subtracting nearly equal Poynting powers. Exactly zero original loss and a
pure evanescent lossless eigenaxis deposit zero heat. No allocation, material
query, boundary solve or depth quadrature occurs in a slab query. A receiving
thermal owner must identify this physical material/volume and supply its own
hold interval; the prepared monochromatic source remains frozen.

The absorbing field owner also admits an authored transverse incidence-plane
axis through `prepare_absorbing_field_in_plane` and its signed-tangential-index
counterpart. Jones inputs are global Cartesian `[Ex,Ey]`, with the original
optic axis fixed in the same physical material frame. Preparation normalizes
the plane axis, rotates both material and illumination into the existing local
Maxwell solve, and rotates all electric/magnetic fields back once. Propagation
retains both global transverse phase components; point and analytic volume
heat contract against the original global loss tensor. The stored signed
tangential wavenumber multiplies `incidence_plane_axis()` to give `[kx,ky]`.
Existing x-z queries keep their original arithmetic. Factors with different
azimuths can share the exact same original interface and thermal column.

`AbsorbingCrystalIllumination` borrows the exact original oriented interface
and an explicitly mutually uncorrelated bank of its prepared absorbing fields.
Each factor includes its square-root statistical weight and retains coherent
ordinary/extraordinary interference, its original frequency, signed incidence,
and spectral loss. Factors may differ in frequency and direction; their field
second moments, energy flux and dielectric work add, never their phasors. This
admits partially polarized and unpolarized light without inventing a mean
complex field or converting the crystal to scalar intensity transmission.
The Cartesian electric coherency is `sum(E_i conj(E_j))/2`, with the peak/RMS
factor applied exactly once. Equal transverse p/s powers at oblique incidence
use tangential amplitudes proportional to the incidence cosine and one,
respectively. The bank supplies finite spectral-line mean powers, not a
bandwidth or PSD inferred from unknown phase; cross-frequency coherence and
deterministic beat envelopes remain outside this source model.

Construction checks exact interface identity and an authored factor budget.
The source supplies no implicit field: an empty bank is darkness, with the
same coordinate and positive-volume admission as a nonempty bank. Observations
use compensated sums, fixed work per factor and no allocations, material
queries or repeated boundary solves. The source-bound slab response retains
the complete illumination, area, both depths, physical volume and analytic
sum of the original coherent volume integrals. Cancellation is typed and
polled before, during and after variable work, with no partial observation
published. Illumination and material remain frozen; a thermal consumer owns
the physical hold interval and energy deposition.

`AbsorbingCrystalCoherentIllumination` instead declares mutual coherence
between a borrowed bank of prepared fields from this exact original interface
and exactly the same positive frequency. Directions, incidence azimuths and
polarizations may differ. Every selected native wave enters one shared
compensated relative-phase `FieldSum` before E/H, Poynting flux or anisotropic
dielectric work is finalized. Ordinary/extraordinary interference and
interference between incident directions are both retained. Equal in-phase
beams produce four times the single-beam heat; opposite phases can cancel it.
Crossed beams produce physical spatial absorption fringes, rather than the
sum of independently positive single-beam heat sources.

The coherent bank allocates nothing, reuses the original solved Maxwell modes
and loss tensor, and bounds preparation and each point query by the authored
component count. Exact interface/frequency identity is checked once on
construction; even one frequency ULP of difference requires a temporal
coherence/exposure model and refuses this monochromatic sum. An empty bank is
darkness with no implied carrier and still applies the original half-space and
finite-coordinate admission. Cancellation is typed and polled before, during
and after variable work. Signed Poynting-only queries avoid requiring the
unrelated dielectric-work magnitude to be representable.

Multi-direction heating generally varies in x and y as well as z, and its
energy balance includes lateral flux. A thermal consumer must integrate the
coherent local work over its actual physical volume. The existing area/depth
single-beam slab integrals cannot be added or reused as a coherent volume
integral. This owner does not infer temporal bandwidth or mutual coherence
from an incoherent source bank.

Focused regression cases cover analytic absorbing isotropic TE fringe
contrast and phase translation, destructive cancellation and a weak residual
between large cancelling beams, noncollinear anisotropic three-dimensional
Poynting balance and boundary continuity, and original-source/frequency
admission, darkness and cancellation.

`AbsorbingCrystalCoherentIllumination::box_absorption` integrates this same
source over an authored `AbsorbingCrystalBox` entirely inside the original
half-space. Finite positive x/y/z extents and the original global origin define
the actual volume. Origin and extent remain separate; DD face coordinates
preserve a narrow box whose upper coordinate would round to its lower one.
The source-bound result retains the complete illumination and box, physical
volume, absorbed watts, six signed outward face powers and an independent
power residual. Face order is x_min, x_max, y_min, y_max, z_min, z_max;
negative face power means net entry. The residual is minus the sum of outward
powers minus absorbed watts, and is never redistributed as heat.

Equal complex wavevectors are coalesced in DD before propagation or quadratic
contraction, preserving exact cancellation and weak coherent residuals.
Fields remain anchored to the original interface; only a common irrelevant
phase is removed. General boxes use the existing DD transverse harmonic
integrals and original stable complex depth integrals for all modal products.
Heat contracts the positive original parallel/perpendicular loss channels;
it is not computed by subtracting nearly equal boundary powers. Each of the
six faces integrates the actual coherent Poynting vector independently,
including both lateral directions and signed backflow.

When all varying relative-phase extents and noncommon depth opacity are at
most 0.125, a fixed positive eight-point Gauss rule per varying axis sums
coherent fields before squaring. Exactly shared normal propagation is
integrated analytically, including opaque tails. This preserves narrow
dark-fringe heat that separately rounded diagonal and cross integrals would
erase. The normalized Gauss remainder coefficient and authored modal rates
screen analytic truncation against the retained positive loss integral.
This is not an interval or floating-point accuracy certificate. Unresolved
modal cancellation or unrepresentable geometry, fields or powers refuses;
there is no adaptive quadrature, stochastic sampling or negative-heat clamp.

`AbsorbingCrystalBoxLimits` admits (2 times component count) squared pairs
before allocation or coalescing. Each pair has at most seven fixed analytic
contractions; the short-box route has at most 8 cubed volume plus 6 times
8 squared surface field samples, each linear in the same admitted count.
Peak logical workspace is linear in that count and includes numerical and
public result headers, original-mode buffers, reusable observation buffers
and explicit integration scratch. Borrowed sources, allocator overhead and
ordinary call-stack locals are excluded. Cancellation is checked throughout;
no partial box result is published.

Focused controls include an independent complex-Fresnel six-face fixture with
nonzero lateral power, translated crossed-beam heat fringes, cubic dark-fringe
heat through opaque common depth, anisotropic point-quadrature and x/y/z
partition checks, single-beam slab reduction, large transverse origins,
darkness, zero loss, exact work/storage admission and cancellation. A thermal
consumer must bind the actual rectangular footprint and original material,
choose its hold interval, and retain all external face powers. The box adds
no second optical interface, lateral thermal resolution or optical feedback.
Native execution of the new box regressions remains unverified; the validation
history below predates these box controls and does not report their execution.


`AbsorbingCrystalPartiallyCoherentIllumination` borrows a bank of those original
coherent angular illuminations and declares distinct factors mutually
uncorrelated. Each factor's original peak amplitudes already contain its
square-root statistical weight. All inter-beam and ordinary/extraordinary
products remain inside that factor; the original dielectric work and six
signed face powers add across factors. Equal-frequency factors do not acquire
cross terms. A positive factorization therefore represents partial spatial
and polarization coherence, including contrast between full interference and
independent illumination, without inventing a mean ensemble phasor.

Every factor must retain the exact same original interface object. Frequencies
may differ between independent factors; each retains its own original
spectral loss tensor. A coherent factor still requires one exact frequency.
These are explicitly weighted spectral-line mean powers, not PSD/Hz, an
inferred bandwidth, finite-exposure beat power or a stochastic realization.
A unitary change of the coherent factor basis at one frequency preserves
ensemble heat and all signed face powers.

The bank's `box_absorption` observes the same physical `AbsorbingCrystalBox`
and retains the complete source/factorization. Heat and six face powers are
summed in DD in authored factor order. The reported residual is minus the sum
of returned outward face powers minus returned heat; it is never deposited.
Any child refusal, unrepresentable aggregate or cancellation publishes no
partial result. There is no repeated material query or Maxwell solve.

Construction is allocation-free O(factor_count), checks the authored
maximum_factors even for empty coherent factors, and retains the checked sum
of each factor's (2 times component_count) squared. `AbsorbingCrystalBoxLimits`
admits that TOTAL pair work and the PEAK workspace before any factor query.
The shared coherent admission helper supplies exactly the same requirements
used by each original query. Peak storage is the aggregate result and seven
DD accumulators plus the largest factor's coherent scratch/result headers,
not the sum across factors. Borrowed sources, allocator overhead and ordinary
stack locals are excluded. Sequential work includes factor_count plus the
sum of admitted modal pairs and each factor's fixed short-box rule. Empty
banks are dark, with no implied carrier, field or material heat.

Authored G1/G3 controls cover finite-footprint fringe contrast and phase,
complex unitary factor mixing with nonzero lateral power, mixed original
spectral laws, independent opposite-phase factors, exact aggregate work/peak
storage admission, source identity, darkness and atomic cancellation.
Native execution of these new factor-bank regressions remains unverified.

Each factor uses a monochromatic passive reciprocal uniaxial half-space with a
lossless isotropic incident port. This lane does not advance temperature, infer magnetic loss,
support an absorbing incident port or finite lossy stack, introduce a biaxial
or gain model, or claim a causal transient reservoir. Unresolved/coalesced
complex modes, ambiguous outgoing branches and unrepresentable source loss or
fields refuse explicitly; the original transparent critical support remains.

Focused actual-source rustc validation passed 49 tests: six new absorption
groups and 43 existing principal/arbitrary-axis/propagation cases. The new
groups cover independent complex Fresnel limits with magnetic/metallic
contrast, analytic dichroism/phase, full complex Maxwell curls and continuity,
rotated coupling and coherent cross terms, differential and integrated
Poynting/heating balance, transparent reduction, evanescence and admission.
This probe uses real fs-math, complete production numerical leaves and four
exact native slab helper bodies in a minimal module wrapper; it excludes
material resolution. The native `fs-em` test executable with `plane-wave` subsequently passed all
85 uniaxial numerical/propagation/stack cases. All 14 native original-source
cases also passed, including retained source identities, spectral support/gaps
and depth-dependent absorption heat and flux. The source fixtures create
revision zero and advance through MatDb's original immutable predecessor chain;
they retain their requested revision and source mismatch assertions. Excessive
spectral gaps still refuse at the original optical-source constructor.

### Native three-dimensional harmonic moving-conductor induction

`SpatialInduction::solve_harmonic` solves the fixed-domain Whitney A-phi
magnetoquasistatic equations for peak phasors with the `exp(-i omega t)`
convention. The response borrows the original prepared induction owner: mesh,
relative tree gauge, material sources, physical boundary traces, floating
conductor components and conductivity tensors remain the same authorities as
the transient solve. The harmonic perturbation has zero prescribed tangential
A on physical fixed traces; solving it does not advance or mutate the accepted
transient state.

The native real two-by-two block represents `K - V - i omega M`, with
`q = phi/(-i omega)`. FGMRES uses a scaled complex-diagonal Jacobi inverse.
Each sparse row is traversed once for both phasor components without allocating
full-field work vectors during an operator application. Identically stationary
conductors skip motion assembly. Admission bounds Krylov and recovered-field
storage, cell magnetic Peclet number, and cell diameter relative to skin depth;
the last is a resolution screen, not an a posteriori discretization certificate.

Recovery retains local complex laboratory and conductor-frame electric fields,
current density, scalar potential and magnetic flux. Original quadrature yields
cycle-mean Ohmic heat, Lorentz force and prescribed-motion work, with the peak
phasor factor of one half. Real supply must close against heat plus mechanical
work. Reactive supply must close against reactive mechanical exchange and twice
frequency times mean magnetic energy. Original Ampere rows, including eliminated
tree-gauge rows, original charge rows including floating gauge vertices, and
discrete magnetic-flux closure are checked before publishing a response. A
Krylov residual alone cannot authorize a result.

This lane supports source-backed linear isotropic or orthotropic magnetic and
conductive laws and prescribed Eulerian velocity. Recoil magnets refuse because
a harmonic perturbation alone cannot certify their full-cycle operating
envelope. Saturation, hysteresis, moving mesh/occupancy transport, temperature
evolution and displacement current are not supplied by this harmonic API.
Physical electrodes and winding circuits use the same native field through
their respective terminal owners. Its source-bound powers and
fields are available to receiving models without relabeling mechanical work
as heat.

Focused regressions cover an independent analytic tetrahedral skin-response
problem over four frequencies, signed moving-conductor real and reactive power,
second-order transient-to-phasor agreement, phase and DC-gauge invariance,
resource/skin/material admission and cancellation. A separate dense reference
checks the nonsymmetric motion operator, its actual transpose, and the complex
Jacobi inverse.

### Mean and twice-frequency induction forces on mechanical receivers

`SpatialHarmonicInductionResponse::lorentz_loads` retains the complete
conduction-current force from one AC carrier: mean `Re(J cross conj(B))/2`
and peak complex twice-frequency coefficient `(J cross B)/2`. The original
tetrahedral P1 integrator supplies consistent nodal loads, receiving-frame
wrenches and arbitrary P1 virtual power for all three real coefficient fields.
Its original local-reference torque reduction and compensated scatter remain
authoritative. A zero cycle-mean force does not erase oscillating vibration
forcing.

The accepted response owns the original geometry, current, flux and volume
sources. Body selections are bounded and canonicalized; repeated or foreign
cells refuse. Cancellation publishes no partial transfer. The mean coefficient
paired with the original steady velocity equals the selected cells' mechanical
power; another receiving velocity supplies virtual work, not coupled feedback.

This is fixed-geometry, one-way AC conduction-current forcing. It does not add
DC/AC cross terms, impressed winding-current force, magnetization, interface
traction or magnetostriction, and must not be combined with an overlapping
enclosing-stress load on the same body. Mechanical motion changing the field
requires a coupled evolution.

Focused G1/G3 regressions independently integrate instantaneous `J(t) cross B(t)`
with degree-two tetrahedral quadrature and check nodal force, torque, nonrigid
virtual work, zero-mean oscillation, quadratic amplitude and doubled phase,
receiving frames, canonical selections and cancellation.

`SpatialHarmonicInductionResponse::cross_lorentz_loads` supplies the missing
interaction of two coherent carriers from the same original induction source
and prescribed steady motion. Each carrier keeps its own resolved current,
skin response and magnetic phase. The additional peak force coefficients are
`(J_a cross B_b + J_b cross B_a)/2` at the sum frequency and
`(J_a cross conj(B_b) + conj(J_b) cross B_a)/2` at the nonnegative difference
frequency (a is the higher carrier). Equal frequencies retain their additional
DC force; no incoherent averaging or beat-frequency cutoff is applied. Add
these terms once per unordered pair alongside the existing self-carrier loads.
The four consistent P1 load fields preserve torque and nonrigid virtual work
without another electromagnetic solve or time-domain sampling. Source and
motion mismatches refuse; the accepted response retains its admitted steady
velocity, included in the original field-storage budget. Storage and work stay
linear in mesh vertices and selected-cell data apart from canonical selection
sorting and conductor lookup. Independent real-field superposition followed by
degree-two tetrahedral quadrature checks true harmonics, close inharmonic beats,
equal-frequency coherence, phase signs, force, torque and work (G1/G3), plus
incompatible-source, selection, resource and cancellation refusals (G0/G4).

### Source-bound 3D harmonic winding terminals and original circuits

`spatial::vector_potential::induction::harmonic::winding` binds the original
closed `SpatialWindingBasis` to the exact same linear magnetic problem and
prescribed-edge traces. The shared spatial winding binder integrates resistance
from each original conductor source, winding profile and authored fill
fraction. Stranded support stays stationary and disjoint from resolved eddy
conductors. Separate resolved conductors retain anisotropic conductivity,
skin/proximity resolution admission and prescribed motional induction.

At each positive frequency, native unit field solves construct the complete
small terminal impedance, including resolved mutual flux and explicitly
omitted end storage. Motion may make this matrix nonreciprocal or active;
neither symmetry nor positive terminal resistance is assumed. Current and
voltage drives recover the actual original field and close real and reactive
terminal work against stranded copper, resolved eddy heat, signed motion work,
active magnetic storage and original factored end storage. Source pairings are
checked independently of potentially large end-region power.

Only reduced unit cochains and solver reports are cached, without retaining
redundant unit B/J/heat fields. A combined drive seeds the original real-block
field solve; a sufficiently accurate seed avoids Krylov iterations but never
bypasses the original Ampere, charge, divergence, skin, motion or power gates.
Voltage drive uses the existing scaled small complex LU and checks every
actually reconstructed terminal voltage against its commanded value.

`harmonic::circuit` closes those same terminals on the original RLC, ideal
transformer and source graph. It reuses the existing native descriptor and
terminal closure shared with planar harmonic windings, including original
source indexing, orientation, Kirchhoff constraints and transformer power.
One descriptor factorization serves terminal unit drives and replay. Actual
field voltages, reconstructed from original end factors, drive the replay;
expanded end matrices do not replace those factors. Numerical reference
terminations cancel and contribute no physical heat. Reported field heat and
signed mechanical power remain separate, including electrical generation.

Both paths retain bounded preparation/drive payloads, original source and body
correspondence, exp(-i omega t) peak phasors and cancellation without transient
history changes. They hold magnetic materials, geometry and motion fixed.
They do not infer strand-scale eddies, end resistance, thermal updates,
saturation or evolving mechanical feedback. Focused winding/circuit
regressions cover independent mutual RL references, declared shared end modes,
actual moving-field work, phase covariance, source refusal and retry.

### Physical massive-conductor AC electrodes

`SpatialTerminalInduction::solve_harmonic_voltage` and
`SpatialCurrentInduction::solve_harmonic_current` extend the original physical
electrodes into the native `K - V - i omega M` solve. Voltage control retains
the prescribed charge traces; total-current control retains the grouped
equipotential unknowns. The spatial current distribution and skin/proximity
response are solved rather than imposed as a uniform or stranded profile.
Both reuse the original sparse operators, source-bound conductivity tensors,
gauge, steady Eulerian motion, skin/Peclet limits and cancellation points.

Peak phasors use `exp(-i omega t)`. Each conductor retains its own voltage
reference. Original charge rows, electrode reactions, component current sums
and all free Ampere rows are checked. The receipt separates physical-electrode
and impressed-volume real/reactive supply; their sum closes against magnetic
storage, conductor-frame Joule heat and signed Lorentz work. These same fields
support the existing consistent mean/twice-frequency mechanical loads. Solving
does not change the accepted transient cochain or clock. Materials and geometry
remain fixed; contact impedance, external lead fields, switch-on transients and
temperature or mechanical feedback are not inferred.

G1/G3/G4 regressions in `tests/spatial_induction/terminals/harmonic.rs` compare
independently integrated voltage- and current-driven diffusion modes at four
frequencies, both source power channels, moving-conductor voltage/current
equivalence, complex common-mode invariance, second-order transient refinement,
resource refusal and cancellation/retry. Synthetic source data qualify the
equations and implementation, not an experimentally validated device.

### Finite-interval three-dimensional induction heat

With `harmonic-exposure` enabled,
`SpatialHarmonicInductionResponse::observe_joule_heat` integrates the original
conductor-frame `E(t)^T sigma E(t)` over a positive physical duration and an
authored start phase. It retains the carrier's twice-frequency heating, which
cycle-mean power times duration cannot resolve on a short interval. Each result
borrows the accepted field and the original material sources and reports
immutable per-tetrahedron and total heat.

The shared normalized `fs-fft` temporal Gram and double-double phase arithmetic
retain the positive short-window limit. The same admitted conductivity acts on
the rotated electric field before quadratic evaluation, avoiding artificial
negative heat from independently rounded E/J zero crossings. The original P1
tetrahedral mass is factored into five vector quadratic forms per cell.
Preparation of the temporal kernel is shared; observation is linear in the
number of conducting cells, with no field re-solve, time sampling or per-cell
trigonometry. Output bytes/cell counts are admitted before allocation and
cancellation prevents partial publication.

This observer reports resolved Ohmic heat only. It does not claim complete
finite-window source work, magnetic storage change or mechanical work, and
does not add stranded copper, operating DC/cross terms, switch-on transients
or temperature feedback. A thermal consumer must retain the same source and
volume correspondence, admit the frozen-field approximation and deposit this
heat once. Four focused tests in
`tests/spatial_induction/harmonic/exposure.rs` cover independent analytic and
space/time quadrature references, tiny intervals at a field zero, phase and
interval composition, original-source identity, budgets and cancellation.

## Analytic bulk electromagnetic dispersion and propagation delay

With `plane-wave`, `DielectricResponse::bulk_dispersion` retains the original
material snapshot and computes the complex bulk wave number, its first and
second angular-frequency derivatives, group slowness and GVD. Original Debye
and Lorentz poles and independently sourced Ohmic leakage supply analytic
derivatives; the material owner's permittivity value and native outgoing
Maxwell root remain authoritative. Scaled frequency derivatives and compensated
pole sums avoid dimensional frequency squares and preserve opposing slopes.
Preparation has bounded pole work and no allocation or frequency differencing.

`BulkWaveDispersion::propagation` reuses the native angular-spectrum complex
exponential and double-double phase reduction for exact carrier transmission,
with the same branch's phase delay, group delay, group-delay dispersion and
log-amplitude derivatives. Nondispersive laws recover exactly zero GVD.
Negative or superluminal inverse phase slopes are retained; they are not
signal-front or energy velocities in strongly absorbing/anomalous regions.
The API describes homogeneous bulk propagation at fixed material conditions,
not oblique modal derivatives, interface residence times, a broadband pulse
approximation or a higher-order dispersion bound. Unsupported original laws
refuse rather than treating one monochromatic sample as a constant law.

`tests/plane_wave_dispersion.rs` checks independent frequency differences of
the native interface solve, geometric invariance of bulk Debye dispersion,
single counting of leakage, resonant/anomalous/evanescent limits, native
six-component beam-field propagation and source/domain refusals.

### Fresnel reflection phase and interface frequency dispersion

`PlanarDielectricInterface::frequency_dispersion` reuses the native passive
carrier solution and both original analytic material jets to return the first
and second frequency derivatives of reflected tangential E and transmitted
tangential E / conjugate H. The caller explicitly selects fixed incidence
angle or fixed transverse wavenumber. The latter follows an original spatial
Fourier mode as its angle varies with frequency; these paths have different
phase slopes even for nondispersive total internal reflection. Common frequency
factors cancel before differentiating, retaining exactly zero derivatives for
nondispersive fixed-angle Fresnel amplitudes. Preparation allocates no storage
and adds only fixed-size scalar arithmetic to the original bounded material jets.

Each nonzero amplitude supplies `-Im(H'/H)` group delay, its frequency
derivative, and both log-amplitude derivatives under the original exp(+i omega t)
convention. The fields are anchored at the interface; no propagation distance
or finite-film residence time is inserted. These are signed local narrowband
phase derivatives, not signal-front or energy velocities. Exact Fresnel/Brewster
zeros retain finite complex amplitude derivatives but return no logarithmic
phase/delay. Exact critical transmission, non-lossless incident jets and
unsupported/nonsmooth source frequencies refuse without changing the native
carrier-only response. No broadband pulse approximation is inferred.

Focused G1/G3 cases in `tests/plane_wave_dispersion.rs` compare independent native
frequency and phase sweeps for Debye, Lorentz, Ohmic and dispersive incident laws;
verify closed-form TIR delay/GDD and binary-exact magnetic Brewster cancellation;
and recover the original reflected angular-spectrum beam derivatives at fixed
transverse carrier. Synthetic source fixtures are not measured-material validation.


### Refracted normal propagation and local dispersion

`InterfaceFrequencyDispersion` retains the original transmitted normal
wavenumber's first and second angular-frequency derivatives.
`transmitted_propagation(normal_distance_m)` reuses the native bounded complex
exponential and phase reduction to return the carrier transfer, phase/group
delay, GDD, and attenuation slope/curvature through the same homogeneous
half-space. The fixed-angle or fixed-transverse-wavenumber frequency path is
the one explicitly selected for the interface. There is no extra material
query, finite differencing, allocation, or fitted refractive law.

These are derivatives of exp(-i*k_z*z) at fixed transverse observation
coordinate. At fixed transverse wavenumber in transparent media, the group
slope is n_group/(c0*cos(theta_t)); even nondispersive media retain geometric
normal-transfer GVD. Nondispersive fixed-angle GVD is exactly zero. The exact
normal-incidence limit reuses the original bulk cancellation safeguards.
Other angles use ordinary f64 derivative contractions, without a relative
error guarantee where weak attenuation or curvature terms nearly cancel.

Fresnel factors and their delay remain separate, so no boundary contribution
is counted twice. Evanescent transfer retains decay and its frequency
derivatives, without assigning normal transmitted power or a ray velocity.
Exact critical branches retain the existing derivative refusal. This local
transfer does not supply arbitrary point-to-point travel time, signal-front
speed, finite-film residence time, or a broadband pulse error bound.
Four added G1/G3 cases in `tests/plane_wave_dispersion.rs` specify analytic
Snell/path limits, transparent and absorbing dispersive native frequency
comparisons, TIR attenuation, and original angular-spectrum transfer and
distance composition. Native execution of these additions remains pending
after the implementation environment exhausted storage; no test pass is claimed.

### Coherent coating and etalon frequency dispersion

`PlanarDielectricLayers::second_order_frequency_dispersion` extends the same
original-source frequency jets to finite coherent films. Its reflected E and
transmitted E/H derivatives include material dispersion, propagation and repeated internal
reflection; transmission phase refers to the actual exit face. The original
stack, geometry, frequency, polarization and explicit fixed-angle or fixed-
transverse-wavenumber constraint remain attached. The native carrier response
and positive film-loss calculation are unchanged. Queries perform bounded
O(layers) work without heap allocation or numerical frequency differencing.

The existing bounded homogeneous E/H slab matrices anchor the carrier. Their
local derivatives hold numerical normalization factors fixed and use the
entire cosine/sinc functions of squared normal phase. A fixed Taylor band
retains finite internal critical-layer derivatives without dividing by zero
normal wavenumber or introducing an evanescent-branch cusp. The same normalized
tail recurrence avoids growing transfer exponentials and nearly-unit reflection
subtraction. Native opaque attenuation limits retain their front-interface jet;
transmitted jets are unavailable after propagation or composed-gain underflow,
since a rounded zero field does not establish zero derivatives. The existing
first-order `frequency_dispersion` and driven-resonance APIs retain their cheaper
specialized tangent path, sharing the same near-critical sinc polynomial.
Critical exits and unsupported original material jets refuse; no phase is
assigned to exact zero amplitudes. These local signed delays/GDD do not supply
an incoherent-film phase, broadband pulse model or absorption derivative.

Focused G1/G3 tests recover bulk and empty-stack limits, the Airy-series etalon
dwell time, propagating/evanescent crossing of an exact internal critical film,
full native lossy/dispersive spectral sweeps, subdivision and opaque limits.
The fixtures verify numerical physics for authored laws, not experimental
coating qualification.

### Coherent transverse dispersion and lateral reflection shifts

`PlanarDielectricLayers::transverse_dispersion` returns original reflected E and
transmitted E/H amplitudes with analytic first and second derivatives in signed
transverse angular wavenumber, at fixed prepared frequency. It retains the
original stack, polarization and entrance/exit phase references. It uses only
the admitted monochromatic material values, so no frequency-derivative law or
extra spectral fit is required. The existing compensated angular-spectrum
light-cone check, bounded slab matrices, entire critical-film functions and
normalized E/H jet recurrence are shared with their original owners.

Under exp(+i omega t - i k_parallel x), the logarithmic phase gradient gives
the local +d arg(H)/dk_parallel displacement along the interface; phase
curvature and log-amplitude filtering derivatives are retained separately.
This is the narrow-spectrum Goos–Hänchen/Artmann interpretation, not an exact
finite-width centroid, transverse spin-Hall effect or exterior propagation
distance. Full finite-beam interference remains with `PlanarAngularSpectrum`.
Exact zero amplitudes have no logarithmic shift, exact critical exits refuse,
internal critical films remain regular, and underflowed transmission jets are
unavailable while the reflected jet remains usable. Queries allocate no heap
storage and take bounded O(layers) work without angular differencing.

Five focused tests cover closed-form signed TIR shifts and curvature, native
lossy coating and critical-film angular sweeps, geometrical exit displacement,
mirror/zero/underflow limits, and convergence of an actual finite TE angular-
spectrum beam's entrance-plane intensity centroid to the local phase shift.

### Native AC circuits on physical solid-conductor electrodes

`harmonic::terminals::SpatialHarmonicTerminalFrequency` prepares the complete
complex impedance between original non-reference electrodes and their own
conducting-component references. Each column comes from a native A-phi
total-current trial; currents redistribute spatially and retain resolved
skin/proximity response, source-bound conductivity and prescribed steady
motion. Multi-electrode bodies retain a shared physical reference, and separate
conductors retain their mutual magnetic response. Motion does not authorize
symmetrizing the port matrix or clipping signed mechanical work.

The original harmonic circuit descriptor closes RLC branches, transformers,
sources and external mutual inductances on those ports. Its final replay uses
independently recovered field currents, voltages, Ohmic heat and Lorentz work.
Every occurrence of a physical electrode must map to the same native circuit
node. Numerical reference resistors cancel before physical loss accounting;
field-owned storage cannot also appear in external mutual-inductance links.
The interface is bounded to at most 16 ports, and its additional numerical
payload is admitted before allocation. One field trial per port prepares the
frequency response. Cached cochains seed one original actual-amplitude solve
per subsequent drive; all original residual and power checks still run.

Independent impressed-current excitation is admitted through
`SpatialHarmonicTerminalFrequency::with_impressed_current`. One additional
zero-terminal-current native trial supplies the affine voltage in
`V = V_induced + Z I`; homogeneous port probes prepare Z without subtracting
large driven fields. Both cochain contributions seed actual loaded recovery.
The shared descriptor closes the affine law on the original RLC/transformer
graph, preserving negative terminal power when the field drives a passive load.
Original circuit supply, independently recovered impressed supply and their
total are reported separately; no numerical reference loss becomes heat.
Additional cochain/history storage is admitted before allocation, while an
all-zero impressed source retains homogeneous preparation cost.

Material state, geometry and steady Eulerian motion remain fixed, and accepted
transient histories are unchanged. The focused G1/G3/G4
tests in `tests/spatial_induction/terminals/harmonic_circuit.rs` compare an
independent diffusion/RLC reference, physical node/reference consistency,
orientation and numerical-reference invariance, moving-field power splitting,
and bounded cancellation/retry behavior. Its `harmonic_circuit/impressed.rs`
child compares induced open-circuit voltage, passive loaded secondaries,
simultaneous independent/terminal drive and their separately integrated work.
Synthetic fixtures are implementation evidence, not experimental device validation.

### Original optical phase, group and attenuation dispersion

The bulk dispersion owner accepts original visible and finite-band optical
responses: Cauchy and Sellmeier phase laws, tabulated dielectrics, complex-index
conductors, explicitly lossless IAPWS water, and original molecular-gas optical
responses. Original phase
and absorption samples, fitted bands, interpolation-gap policies, material
snapshots and absolute permeability remain authoritative. Index laws are not
replaced by fits or a monochromatic constant. The native Maxwell carrier root
and propagation exponential remain unchanged.

Analytic wavelength derivatives supply group slowness, GVD and both attenuation
derivatives without frequency sampling or allocation. Tabulated n/k retains
exactly zero local GVD within a linear wavelength segment. Dielectric intensity
absorption is differentiated before conversion to extinction: constant alpha
has exactly zero attenuation slope/curvature even though extinction varies.
Water curvature comes from the original material owner at fixed T and density.
Molecular-mixture derivatives retain the original Cauchy component laws and
nonlinear Lorentz–Lorenz mixing at fixed actual partial number densities. Preparation
caches the wavelength jet; bulk queries at the exact prepared carrier supply
group delay, GVD and GDD without refitting, allocation or frequency differences.
Vacuum and nondispersive composition have exactly zero GVD. Original source-band
endpoints remain valid for monochromatic transport but refuse two-sided
dispersion; no gas absorption or equation-of-state update is inferred.
`MolecularGasOpticalResponse::try_lossless_ciddor_air` also retains the original
NIST Ciddor material law at fixed absolute temperature, pressure and composition.
Its cached phase/group indices and spectral curvature drive the same native
Maxwell interface, layer, angular-mode and bulk dispersion owners. The constructor
explicitly selects a lossless coherent approximation. Magnetic source conditions
must be matching static temperature/pressure or legacy kelvin T; unknown
conditions, including a proxy mass density, refuse. Original atmospheric support,
exact carrier admission and source-state identity remain binding. The material
owner's optical density correction does not replace a fluid equation of state.
Band endpoints and unequal adjacent table slopes refuse two-sided derivatives;
smooth knots remain admitted. These are local original-law derivatives, not
measured derivative accuracy, a causal reservoir or broadband-pulse error bounds.

Focused existing optical/material regressions compare original Maxwell and
permittivity frequency differences, closed Cauchy/linear-index references,
constant and varying power absorption, infrared Sellmeier GVD including its
sign change, water phase/group curvature, source identity and domain refusals.

### Settled AC views of original material-temperature trials

`SpatialInductionMaterialTrial::into_harmonic` exposes the same admitted
source-temperature operators through a dedicated `SpatialHarmonicMaterialTrial`.
Geometry, fixed trace, gauge, conducting support and source selections remain
native, and the original sparse harmonic solver retains its charge, Ampere,
skin/motion and real/reactive-power gates. The view shares the already prepared
operators; it does not assemble a second material matrix or expose a transient
state update. Original magnetic cochains, impressed-current history and clock
remain unchanged.

This view does not supply magnetic parameter-work or thermal accounting. A
consumer must declare its settled-carrier approximation or provide a separate
transient/thermodynamic account. The source-temperature AC heating consumer in
`fs-couple` checks fixed permeability and independently authored caloric/heat
transport before using these responses in its thermal residual.

### Source-temperature AC views of physical electrodes

`SpatialTerminalInduction::prepare_harmonic_materials_trial` and
`SpatialCurrentInduction::prepare_harmonic_materials_trial` refresh the original
source-bound coefficients while retaining their physical voltage or grouped
total-current electrode spaces. Each dedicated harmonic view shares the native
geometry and uses one material preparation. Original electrode references,
charge constraints, source identities and accepted transient histories remain
intact. The view exposes settled AC recovery, including simultaneous impressed
excitation and actual terminal reactions, without publishing a transient
material endpoint or accounting for magnetic parameter work.

The physical-electrode thermal consumers in `fs-couple` retain the same fixed
permeability and slow-thermal restrictions as isolated AC heating. Their focused
regression definitions compare independent current/voltage diffusion and source
work at the accepted temperature, coherent Joule heat, common-mode invariance
and transactional material/thermal refusal.

### Loaded AC views of source-temperature electrode trials

`SpatialHarmonicCurrentMaterialTrial::prepare_harmonic_frequency` exposes the
original affine electrode/circuit response at the admitted material temperature.
Its restricted frequency view retains native port probes, impressed excitation,
cochains and original diagnostics without exposing a refreshed transient owner.
Preparation reports the actual field-trial count; each current/circuit recovery
adds one native solve even when the cached seed already satisfies the residual.
`HarmonicCircuitResponse::map_field` moves the complete circuit receipt alongside
an owned field observation without cloning network vectors or recomputing work.
Existing field/circuit material, topology, residual, storage and cancellation
admissions remain authoritative. A harmonic trial does not supply transient
magnetic parameter work or advance accepted history.

## Physical boxes inside original absorbing crystal films

The existing `AbsorbingSlabCoherentIllumination` and
`AbsorbingSlabIllumination` now integrate their original fields over an
`AbsorbingCrystalBox` through `box_absorption`. The former retains coherent
angular interference; the latter adds the expected powers of its already
declared independent groups, including groups at different original
frequencies. No statistical model, hidden weight or optical solve is added.

Each original beam retains all four forward/backward eigenwaves and both
authored incident ports. The coherent and ensemble results borrow the complete
illumination and original film, and expose physical volume, dielectric watts,
six signed outward face powers, an independent residual and resource counts.
Face order is x_min, x_max, y_min, y_max, z_min, z_max. Negative means entry.
Lateral transport and the residual are never reassigned to material heat.
The exact DD upper depth must remain within the original film thickness,
including for dark illumination; a rounded endpoint cannot silently admit
support outside the film.

The shared native box kernel now retains each wave's original entrance or
exit anchor and attenuation direction. Fields are attenuated from their
nearest box endpoint; a fixed z face uses only its actual physical plane.
Relative phases retain both original anchors. Same-direction products use a
decaying depth integral from their common near endpoint. Opposite directions
use the absolute attenuation-rate difference and the remaining factor
exp(-min(alpha_a,alpha_b)*depth); a backward-selected integral also retains
its exact upper-end relative phase. No exponentially growing propagation or
transfer matrix is introduced, even in an opaque two-sided film.

Short boxes retain the existing positive coherent-field rule and analytical
truncation screen, with backward samples attenuated from the upper box face.
Exactly common normal propagation remains analytically integrated. Original
loss channels, independent six-face Poynting contractions, source-bound
cancellation and representability refusals remain. These screens are not an
interval or floating-point accuracy certificate; unresolved cancellation
still refuses rather than clipping negative heat.

Coherent work admits (4 times component_count) squared pairs before exact
anchored-wave coalescing; an independent bank admits their checked sum.
The existing source constructors retain those counts and bound empty-group
traversal. Shared requirements keep query admission identical to cached
metadata. Storage includes result headers and explicit integration scratch;
a bank retains seven DD accumulators plus the largest child workspace,
evaluating groups sequentially. Borrowed sources, allocator overhead and
ordinary stack locals are excluded. The half-space owner keeps its original
forward-wave integral path and uses the same updated anchored-mode storage
accounting. No partial result escapes cancellation or a failed group.

Authored G1/G3 controls cover opaque two-sided films with all six boundary
powers, the independent TE lateral-work identity, tilted noncollinear
two-port beams against original local E/H and tensor work, x/y/z partition
invariance, positive cubic loss at a short standing-wave node, original
spectral groups, exact work/peak-memory budgets and DD film-support admission.
Native execution of these new film-box regressions remains unverified.
The source remains frozen; thermal time stepping, thermo-optic feedback and
new lateral optical boundaries are outside this electromagnetic observation.

### Cached polarized fields through absorbing crystal stacks

`PreparedAbsorbingStackPolarization::field_at` reuses its two original global
tangential Maxwell basis solves for arbitrary complex Ex/Ey peak amplitudes.
The selected layer's four electric and magnetic eigenwaves combine before
the native compensated propagation, signed Poynting contraction and tensor
heating. Reflection, transmission and incident-side standing fields retain
their original physical face coordinates. Each query takes fixed work and
storage, with no allocation, layer traversal, material query or boundary solve.
The existing layer-integrated heat query retains its electric-only work.

`field_statistics` accepts explicitly independent positive Jones factors at
that same prepared frequency and direction. Each factor includes its own
square-root statistical weight; the shared native accumulator adds Cartesian
electric coherency, magnetic RMS squared, signed flux and local heat without
inventing an ensemble phasor. This represents unpolarized and partially
polarized illumination. Empty/zero factors are dark, but region admission
still applies. The caller bounds factor count; cancellation is checked before,
between and after factors, and no partial observation escapes a refusal.
There is no inferred bandwidth, temporal coherence or thermal state update.

Focused G1/G3/G0 controls cover arbitrary-axis elliptical fields and Maxwell
continuity, an independent matched-isotropic circular/unpolarized limit,
complex unitary factorization invariance, opposite-phase independent factors,
layer heat, darkness, invalid geometry, factor budgets and cancellation.


## Original air-relative catalog glass in Maxwell transport

`SpectralOpticalSource::LosslessAirRelativeSellmeier` borrows the shared material catalog law and an explicitly selected scalar magnetic source. The original `SpectralOpticalResponse` supplies absolute complex index with zero imaginary part, SI permittivity and existing interface/layer consumers. Its frequency jet converts the complete analytic ordinary wavelength curvature through the existing bulk-dispersion machinery, retaining group delay and group-velocity dispersion. Derivatives require an interior source frequency; monochromatic queries retain inclusive mapped support.

The magnetic source must be unconditioned: the catalog's reference-air temperature and pressure do not establish the actual glass state. Refresh through unrelated material-card claims refuses explicitly. This source is an authored lossless approximation, not measured zero absorption or a causal time-domain model. Three shared material controls and five renderer-consumer controls passed. The native EM phase/group/GVD, SI permeability, exact-band and source-ownership regression is authored and signature-reviewed but remains unexecuted in this environment.

### Dispersive coherent-film resonances

`PlanarDielectricLayers::frequency_dispersion` returns analytic first angular-
frequency derivatives of the native reflected and transmitted electric/magnetic
amplitudes, retaining original dielectric source laws and scalar permeability
snapshots. TE/TM and the existing fixed-angle/fixed-transverse-wave-number
constraints use the same original material and incident-port derivatives.

The tangent differentiates the native bounded E/H transfer in a fixed numerical
gauge. Entire functions of the longitudinal phase squared preserve the finite
critical-film limit; no growing propagation exponential or finite difference
kernel is used. Double-double products audit the transfer/tangent rows.
Stack recursion is linear in film count with no per-layer heap workspace;
constitutive evaluations retain their own cost. Original 64-film admission
remains active.

`driven_resonance` reuses the shared peak/half-height search and rebuilds each
original source at the actual sampled frequency. Its evaluation budget and
cancellation remain explicit. Source validity bands, prepared-frequency
restrictions and loss of propagating incidence refuse through native owners.
Transmitted tangents are unavailable after propagation, gain or composed-field
underflow, while the native opaque reflection limit remains usable.

These are selected unit-incident-electric-field amplitude peaks and phase
slopes. A transmitted-amplitude maximum need not maximize transmitted flux;
its half-height quality factor need not equal a pole-decay quality factor.
An exactly critical exit half-space retains the existing Fresnel derivative
refusal, and no finite-band pulse, inferred material law or completeness claim
is added.

### Resonance of the original transmitted energy flux

The shared `PowerFractionDispersion` owner differentiates the original normal
Poynting-flux fraction using both coherent E and H fields, including incident
and exit admittance dispersion. Double-double products retain phase
cancellation. The half-space Fresnel owner and both first/second-order layered
field owners expose the same reflected/transmitted power-fraction slopes.
Original positive carrier fractions remain authoritative; absorption heat and
an optical phase are not inferred from a power ratio.

`LayeredResonanceObservation::TransmittedPowerFraction` resolves a peak and
optional half-height bandwidth of the ACTUAL native transmitted flux. The
shared peak search receives the phase-free observation `sqrt(T)`, with exact
analytic slope `T'/(2 sqrt(T))`; its squared observation is the power fraction.
This retains changing source/port admittances, which can move the flux peak
away from a transmitted-electric-field maximum. Strictly positive resolved
transmission is required at derivative samples. Underflow keeps the original
unavailable-field refusal, and zero-flux samples do not invent a square-root
tangent. Existing source domains, incidence constraints, budgets and
cancellation remain active.


## Physical absorption and heating resonances in coherent films

`PlanarDielectricLayers::frequency_absorption` resolves the complete analytic
first frequency derivative of power dissipated in every original film, together
with the native field/flux response. `LayeredAbsorptionDispersion` exposes each
film's and the total absorbed incident-power fraction. Its `layer_power_density`
and `total_power_density` queries expose actual cycle-mean heat per area and its
frequency slope for a fixed supplied complex peak tangential incident electric
field. `AbsorbedPowerDensityDispersion` reports W/m², `omega*dP/domega`, and
`dP/domega` in W s/m². The incident admittance derivative is included, so a
fixed-field heating peak need not coincide with an absorptance peak.

The existing bounded E/H transfer and exact material-jet sweep is reused. One
additional forward source sweep retains the coherent incoming/opposite fields
of every film. Analytic derivatives act on the original positive even/odd slab
loss channels, their interference, physical attenuation/phase, and dispersive
reference port. Small-depth odd cubic loss is differentiated directly; neither
absorption nor its derivative is obtained by subtracting reflected and
transmitted power. The independent differentiated Poynting balance is measured
and admitted separately. Preparation is O(layer_count), with fixed 64-film
scratch storage and one bounded output vector, without a depth mesh, frequency
differences or a separate stack solve per selected film.

`LayeredResonanceObservation` now selects total or individual-film absorbed power
fractions and physical heat densities, in addition to existing field and
transmitted-power observations. The shared search uses a phase-free square root
of the selected positive power; its peak and half-height crossings therefore
belong to that physical quantity. Each sampled frequency reevaluates original
material laws under the declared fixed-angle or fixed-transverse-wave-number
constraint. Numerical observed linewidth is not a free-mode decay certificate.

Zero-width and lossless films retain zero heat. A lossless carrier with nonzero
constitutive loss slope, unresolved interior propagation/source underflow, a
zero-power derivative observation, material-domain failure or exhausted search
budget refuses explicitly. The ordinary carrier and front-reflection paths
remain available independently of interior derivative admission. Fixed scalar
permeability and immutable material snapshots retain the native domain; this does
not evolve temperature, infer nonlinear heating or couple a thermal reservoir.

Physical source-power contractions retain extreme but representable source,
admittance and loss products before squaring. The native positive-loss amplitude
also retains subnormal loss coefficients before contraction with large fields;
ordinary coefficients preserve the existing numerical path. Six focused native
numerical tests cover independent TE/TM loss-channel slopes, thin weak loss,
near-critical geometry/material directions, and representable power across these
scales. Public stack tests additionally exercise original per-film heat, source
normalization, actual Airy absorption peaks and half-height crossings.

## Independently sourced current absorption in thermo-optic glass

`ThermoOpticOpticalResponse::try_with_current_absorption` combines original
SCHOTT-corrected reference phase with `SpectralDielectricAbsorptionStatePoint`
at the actual temperature. That complete state supplies mu and retains every
selected current requirement. Original PHASE, thermal and current-loss bands
intersect; the current table's own maximum gap is admitted. The Sellmeier fit's
support remains independent of discarded reference alpha support and magnitude.
The resulting epsilon uses n(T)-i alpha(T,lambda) lambda/(4 pi), replacing
reference loss once. Native interfaces, coherent fields, angular spectra,
local heating and frequency jets consume it without another attenuation term.

Analytic complex jets retain current alpha slopes, exact zero attenuation
dispersion for constant alpha, and current-loss or reference-phase derivative
cusps. Discarded reference-loss cusps and permittivity overflow do not limit
independently supported phase/current-loss composition. `at_absorption_state`
rebinds loss/mu together; `at_state` refuses this mode to prevent silently
freezing current absorption. Existing one-coordinate absolute-temperature
admission remains binding, including refusal of pressure/composition axes.
No absorption(T) law or causal time-domain reservoir is inferred beyond sources.
Three focused regressions specify independent support, complex frequency jets,
cusp selection, reference-loss overflow independence and TE/TM Airy/local-work
agreement. Native execution is not claimed; the final source review followed
loss of the executor after the earlier draft's pinned formatting checks.

## Positive short-window film populations and reservoir heat

Uniform-window film observations evaluate admitted short-band electric/magnetic populations and each original Debye/Lorentz friction heat as positive orthogonal temporal squares. Original angular frequencies and source-bank statistics remain attached: proper factors retain their peak-phasor half before products, while real fields retain the full carrier and sum-frequency content. The fixed factor is reused across grid points and reservoirs without observation heap allocation; preparation reserves its bounded projection work. Signed Poynting flux and reversible electromagnetic/restoring/kinetic storage use cross products of the same projected fields and their derivatives. Proper derivatives remove the native common carrier in DD before multiplication; real derivatives retain the full carrier and sum-frequency content. Wider or shaped windows retain the existing pair law.

Focused source-level regressions resolve three-tone T^4 and real T^6 heat nulls through 1 ps, preserve filtered-product residuals and extreme proper power, and retain signed flux/storage and cancellation behavior. This is numerical evidence, not a null-relative certificate. Prepared local and spatially integrated heat spectra are separate owners; this direct-layer change does not alter their observation law.

## Retained factors for reusable local film-heat windows

Prepared local film-heat owners retain the original DD reservoir fields and source/reservoir boundaries already included in their admitted preparation workspace. These buffers live with the owner; they add no preparation heap peak or per-observation allocation. Each admitted source bank uses the short uniform temporal factor and reuses it across all its reservoirs, without another FFT or material solve. Other banks retain their original real/proper DD pair law, sharing each pair kernel across the source's reservoirs within the already admitted ordered-product work. A silent one- or two-tone bank or a wide-carrier real source cannot disable another source's positive contraction. When no bank admits, observation keeps the cached flattened spectrum path. Original conductivity normalization and the proper-complex peak half remain attached to their physical source laws.

The exported flattened spectrum and spatially integrated heat spectra retain their coefficient representation. Instantaneous local evaluation contracts the retained coherent reservoir fields before the real/proper square. This improvement concerns the retained local owner's window observations. Focused regressions exercise actual spectrum preparation, real/proper deep nulls, independent source/reservoir contributions, shifted shutters and native EM phase, conductivity scaling, mixed-bank fallback, finite extreme proper power and cancellation.

## Polarized reflection from correlated thickness ensembles

`PreparedLayerThicknessEnsemble::evaluate_polarized` solves TE and TM for
each same authored joint geometry row, retaining the original cached materials
and coherent Maxwell operator. An explicit work limit includes both solves.
Compensated central moments retain the complex TE/TM response covariance,
rather than inferring a Jones matrix from averaged powers. The result exposes
total, mean-response and fluctuating response dyads, an averaged Mueller
matrix, and constant-work reflected-coherency propagation.

Power-normalized channels use the original fixed tangential-E signs and
positive incident normal admittances. Consumers must explicitly convert to
their own rotating optical frames. Direct coherency propagation retains weak
Brewster channels without reconstructing them from rounded I/Q differences.
Physical finite-film absorption and exit-face transmitted power remain the
original channel-diagonal fluxes; reflected covariance is never heat.
Incident polarization statistics must be independent of the thickness row.
This model supplies neither lateral scattering/diffraction nor an inferred
common-plane transmitted coherency or a physical source-coherence law.

Five G1/G3 regression groups cover rowwise lossy Airy/Jones agreement and
depolarization, single-realization parity, the lossless normal-incidence phase
ensemble, weak channels and small covariance, physical input admission,
deterministic replay, explicit work bounds and cancellation between TE/TM.

The layer query reuses each field projection across its diagonal populations, signed flux and storage, with no positivity clamp on cross moments. Preparation reserves N*(6+3D+5L)*(F+1)*17 bounded vector work. Focused regressions resolve proper three-tone signed flux and real fifth-order storage through 1 ps, retain DD filtered storage residuals, and compare all vector moments and original Debye/Lorentz work identities against independent quadrature. Continuous spatial cell power and volume heat remain separate owners; this temporal refinement is not a null-relative or global conservation certificate.

## Matched-cutoff crystal propagation at grazing incidence

The arbitrary-axis transparent crystal kernel retains the caller's original positive incident normal coordinate at an exactly matched extraordinary cutoff. It evaluates the normal discriminant root from that coordinate, rather than subtracting a rounded tangential square from the cutoff. This preserves propagating fields and finite transmitted power even when the authored cosine is too small for its square to survive subtraction from one. The shared eigenbasis carries this behavior into finite crystal consumers without a cutoff tolerance or an artificial critical mode. Explicit critical and adjacent representable propagating/evanescent inputs retain their distinct meanings.

Native Maxwell and power regressions cover aligned and out-of-plane axes, axis reversal, grazing cosines through 1e-30, and adjacent representable cutoff indices. At a matched principal extraordinary cutoff with incident/extraordinary index 1.5 and ordinary index 1, the transmitted power remains 0.96 rather than collapsing at grazing incidence.

Retained local film heat evaluates each original coherent reservoir field before its real/proper square. Native angular clock differences remain DD; real fields restore the absolute carrier once, while proper factors use their authored half-norm second moment. Conductivity normalization precedes the square, and independent sources/reservoirs add only after contraction. Instantaneous work is linear in admitted retained reservoir lines with fixed scratch, cancellation polls, and no query allocation, FFT or material solve. Exported signed tones and window laws remain unchanged; this does not promise relative accuracy at every arbitrarily deep null.

## Original incidence through finite crystal layers and exit ports

Transparent crystal slabs and multilayer stacks retain the original incident coordinate through every layer and into the isotropic exit basis. Positive authored cosines use the compensated tangential square; equal entrance/exit phase indices retain the original positive normal coordinate directly, with their separate permeabilities still determining E/H and power. Explicit signed tangential-index queries keep their original critical and evanescent meanings. Layer composition no longer reconstructs incidence from a rounded tangential scalar. Transparent and absorbing finite consumers share the same exit-field constructor; absorbing matched-port behavior remains source-bound.

Homogeneous slab/stack controls check unit transmission and finite internal E/H with an observable normal phase at grazing incidence, including multiple authored optic-axis orientations. Original material-source controls match their spectra to the actual resolved incident index rather than assuming nominal decimal indices reconstruct identically. Additional native exit controls retain adjacent propagation/evanescence, exact critical inputs, magnetic contrast and signed-tangent parity. Unresolved boundary systems and exactly grazing incident ports retain their explicit refusals.

## Coherent volume heat over finite source intervals

`PreparedAngularSpectrumExposureLayer::integrated_heat_window` returns mean
irreversible watts from original fields projected before temporal power
contraction. `integrated_heat_interval_from_clock_dd` returns double-double
joules over a positive uniform interval while retaining the original DD source
epoch. The prepared heat-load spectrum exposes this explicit reconstruction as
`observe_window_from_fields`; its cached `observe_window` cost is unchanged.

Each source bank independently selects the shared bounded uniform factor
(3–64 proper tones or 3–32 real tones, at most 17 modes) when its native band is
admitted. Other bands and profiles retain the original spectral law. Proper
sources remove only an irrelevant common carrier; real sources retain carrier
phase and both complex spatial quadratures. Sampled/masked areas, continuous
Parseval/opposite modes and centered rectangle sum/difference integrals retain
the original geometry and carrier. Odd-sized axes pair their signed modes.
Original reservoir filters, proper half-power, volume and duration scales meet
in exponent-scaled DD arithmetic; negative or unrepresentable results refuse.

`AngularSpectrumHeatWindowLimits` admits reconstruction, FFT, temporal/spatial
products, fallback terms and heap before field work. One depth's frequency
planes serve all reservoirs: extra plane storage is O(N*(F+M)), with no retained
depth history. Rectangles retain their exact quadratic spatial contraction.
Fixed bounded projection scratch, borrowed inputs and allocator overhead are
outside the reported heap bound. Unsupported spectral fallbacks retain their
existing coefficient range limits; no relative-error certificate at a null,
depth-convergence claim or material transient is inferred.

Nine actual-source numerical tests pass, including independent time/space
quadrature, fourth-order dark leakage, DD epochs, physical dynamic range,
original depth/area weights and bounded cancellation. Original math, harmonic,
rectangle and volume kernels were compiled; the complete source-owner,
material and FFT adapter graph was inspected but not built or executed.

## Coherent exterior stress and field-momentum modes

Existing deterministic and ensemble stress, momentum and joint observations
select the shared positive temporal field factor independently for each admitted
uniform bank. Six original E/H components are projected before stress and
field-momentum products; momentum storage also projects their derivatives.
Joint observations share those projections and one original field reconstruction
and FFT synthesis per frequency. Proper derivatives remove only an irrelevant
common carrier; actual real fields retain full carrier and sum-frequency terms.
Signed stress, momentum density and its weighted instantaneous rate are never
clipped. Wider bands and shaped shutters retain the original pair fallback.

Continuous cell integration preserves signed opposite modes and unpaired
Nyquist modes. Nonzero real carriers retain the original separable spatial
Gram law, with source/cosine and i-source/i-sine temporal representations;
sampled quadratic products do not replace physical cell integration. Exterior
material admission and final constitutive scaling are unchanged. Mechanical
force still requires the consumer's complete surface/volume balance and any
material momentum or striction law.

Temporal work is O(N*F*M), M<=17, on admitted banks; existing field and carrier
spatial work remains. The helper adds no heap beyond original field planes.
Individual query workspace includes 67,840 bytes of conservative fixed scratch;
the joint composition conservatively counts that allowance twice. Seven
actual-source regressions pass with warnings denied: proper fourth-order and
real cubic dark residuals, carrier shifts, signed time/ensemble tensor oracles,
odd/Nyquist cells, carrier-aware time/space quadrature, cancellation and precise
DD epochs. The real 1ps null retains its documented DD cancellation limit;
no null-relative certificate is asserted. Formatting and a bounded adapter
signature probe passed; full material/reconstruction/FFT execution remains unrun.

## Common-plane transmitted polarization from thickness ensembles

`PreparedLayerThicknessEnsemble::evaluate_polarized_at_plane` evaluates the
original paired TE/TM responses once per illuminated joint geometry and
propagates each transmitted field to an explicit shared downstream z plane.
The plane is measured from the common entrance and must lie beyond every
authored exit, including zero-probability rows. Compensated summed thicknesses
and propagation phase retain sub-ulp exit-position differences at distant
planes; no backwards continuation through an absorbing exit is performed.

Transmitted amplitudes retain tangential-E phase and positive normal-power
normalization. The better-conditioned original E/H chart supplies phase
without forming a singular admittance magnitude ratio. Total, coherent and
fluctuating Jones-response moments use compensated online central covariance;
direct coherency and Mueller queries repeat no film or material solve. Both
polarizations share each authored geometry. Incident statistics must remain
independent of that geometry.

Positive physical absorption in the exit medium between each exit face and
the common plane is reported separately from finite-film heat and polarization
covariance. Reflected power + plane transmission + finite-film absorption +
exit-segment absorption closes the original incident power. Evanescent and
critical exit channels carry zero normal-power coherency, not a claim of zero
near-field electric fields. The independent-column/time model supplies no
lateral diffraction, coherence-length inference, or transient material law.
Work retains the explicit two-polarization layer-realization bound, bounded
geometry admission, no query allocation and transactional cancellation.

Five G1/G3 groups specify common-plane homogeneous propagation (including
distant carriers), rowwise lossy Airy/Jones optics and plane translation,
unequal-port flux/TIR, tiny covariance and subdivision invariance, plus work,
dark-row geometry, input admission, replay and cancellation.

## Current-state absorption for original water and silica refraction

`WaterOpticalResponse::try_with_current_absorption` and
`SilicaOpticalResponse::try_with_current_absorption` compose the original IAPWS
or CHARMS analytic phase with the canonical independently sourced intensity
absorption spectrum. One complete snapshot owns temperature, permeability,
original wavelength/alpha claims and, for water, actual bulk density. No
phase-index table, nonmagnetic default or loss fit is synthesized. The existing
explicit lossless constructors remain separate policies.

Phase and absorption bands intersect in their original frequency coordinates;
zero absorption does not erase source support. Preparation retains a positive
finite wavelength-gap admission. Runtime complex-index and frequency-jet
queries allocate nothing. The shared current-alpha derivative owner preserves
constant-alpha attenuation exactly, original interpolation cusps and two-sided
endpoint refusals while retaining analytic phase/group/GVD. Source matching
includes the loss policy and gap admission, so coherent consumers cannot
silently combine different optical assumptions from the same scalar snapshot.

Existing native interfaces, coherent finite films and angular spectra consume
N=n-i*alpha*lambda/(4*pi) and epsilon=N²/(mu*c0²); absorbed work enters once
through that complete complex permittivity. Absorbing angle-based incident
ports and causal finite-window heat retain their existing refusals. This
composition supplies no scattering, stress, causal reservoir, phase correction
from absorption or empirical qualification of the independently supplied loss.
Three authored regression groups cover oblique TE/TM Airy fields and local
work, analytic complex group/GVD and Beer attenuation, and source/cusp/domain
admission. These tests are defined but were not executed in this environment.

## Coherent film antiresonance and recovery bands

`PlanarDielectricLayers::driven_antiresonance` resolves an isolated minimum of
the original complex reflected electric, transmitted electric or transmitted
magnetic amplitude. Exact reflection cancellation retains finite transmitted
fields and energy flux. Peak and notch queries now share one immutable physical
experiment: original incident/exit and film constitutive sources, original
thicknesses, actual-frequency material evaluation, and either fixed angle or
fixed original transverse wave number. The same native coherent carrier and
analytic material/field tangent preserve phase and source-domain refusals.

Only raw coherent field observations are admitted by the notch adapter.
Phase-free square roots of power/absorption have a cusp at zero and are not
substituted for a differentiable field. Reflected |r|² is actual reflectance;
transmitted E/H observation squares retain their named field-ratio units and
are not generally energy-flux fractions. An optional explicit recovery level
requires both crossings and reports their width without inventing a notch Q
or a phase delay at zero. The shared query/residual/cancellation bounds apply;
recovery queries skip the tangent. Existing absorption and flux-peak behavior
uses the identical shared sampling path.

Three authored public regressions use independent Airy reflection and recovery
formulas for both polarizations, retain unit flux through a half-wave reflection
zero, test a fixed-transverse notch shifted from the prepared carrier, and check
invalid power observations, missing brackets/edges and cancellation. The public
EM tests remain unexecuted while the full dependency build is unavailable; the
extraction and source/constraint/lifetime contracts were independently reviewed.

## Direct sparse tetrahedral conduction electrodes

`SpatialConduction::solve_terminals` identifies complete exterior face electrodes
on the original scalar/tensor conductivity and finite-contact Hodge, then reuses
the geometry-neutral grouped sparse checked-CG solve at the actual drive.
Shared `ConductionTerminalControl` admits voltage, conserved total-current and
mixed/floating controls; `ConductionTerminalLimits` bounds the extra sparse
workspace. All bodies/electrodes need a common bulk/contact/electrode connection.
Original difference-form reactions remove internal equipotential links before
summation while retaining signed anisotropic/contact couplings. Native current,
local bulk/contact heat and Hodge/source-work gates remain independent of the
authored terminal-current/work checks. Absolute field observations retain the
held bias; the terminal report keeps physical differences relative to terminal
zero, including induced drops below the absolute display precision. This is
stationary fixed-geometry DC; it does not add motion, induction or a nonlinear
terminal matrix. Three native regression groups in
`tests/spatial_conduction/direct_terminals.rs` cover analytic contact/floating
networks, rotated tensor sources, gauge invariance, bounds and cancellation.

## Direct tetrahedral terminal controls in nonlinear material response

`SpatialConductionResponseDrive::Terminals` accepts the shared explicit
`ConductionTerminalDrive` and closes original scalar bulk/contact laws at the
actual constrained field with one sparse field solve per constitutive iteration.
Its sparse workspace admission is separate from the existing
`MixedTerminals` dense response budgets, whose meanings are unchanged.
The accepted `SpatialConductionResponseReport::terminals` retains ordered native
reactions, physical voltage differences and independently authored source work.
Original source queries, relaxation, local current/power defects, cancellation
and field-solve budgets remain in the same constitutive loop. This adds neither
a cached nonlinear conductance matrix nor arbitrary nonuniform nonlinear contact
laws. Native absolute field observations and the terminal-relative report keep
their distinct meanings under a large common bias. Regression definitions cover
narrow actual E/J source support, mixed/floating finite-contact loads, both signs,
work, iteration counts and rollback under source/workspace/budget refusal.

## Retained dispersive-film field and polarization energy

`PreparedAngularSpectrumExposureLayer::point_energy` prepares a retained
original-field energy observer at one original grid point and selected depth.
`evaluate_from_clock_dd` returns instantaneous ensemble-mean energy density;
`change_from_clock_dd` returns signed end-minus-start changes over a positive
finite duration. DD results name electromagnetic, polarization restoring,
polarization kinetic and total subaccounts [J/m^3], each included once.

The original causal material law supplies
`U=(epsilon_inf*E^2+mu*H^2)/2+sum(P^2/Delta)/2
+sum_Lorentz(P_dot^2/(omega0^2*Delta))/2`.
Preparation reconstructs each original frequency field once, retaining E/H and
normalized restoring/kinetic coordinates. Native positive-time phasors are
conjugated into the shared negative-time law without rounding filtered products.
Real/proper source identities and the physical proper half-factor remain intact.

Queries allocate no heap, perform no new field/material solve and use linear
coordinate-frequency work. Explicit limits bound original point-line/FFT work,
`3*sum_b F_b*(2+D+2L)` coordinate-frequency terms, original owner workspace,
retained fields/banks and temporary frequency-point storage. Cancellation
publishes no partial observer or result.

Coherent field increments are formed analytically before quadratic changes,
using `delta(q^2)=delta_q*(2*q+delta_q)`. Start epoch and duration remain
separate. Proper factors remove only their irrelevant common carrier; real
fields retain its rotation and increment. The original scaled-DD arithmetic
also retains the phase increment's small quadratic part, allowing finite
changes when bright endpoint energies exceed ordinary DD range. Unrepresentable
source phases or final physical subaccounts are explicitly refused.

Seven numerical regression groups pass with warnings denied: real/proper dark
energy through 1 ps; changes below DD endpoint resolution; precise-clock
partition and carrier invariance; independent mixed ensembles of every store;
polarization heat plus storage change versus independently integrated original
work; physical scale/range; and cancellation/refusal. These execute the new
energy helper, original scaled arithmetic and actual math/harmonic source.
The complete preparation adapter body type-checks with narrow owner metadata
carriers; native material admission, Maxwell reconstruction/FFT and full fs-em/
workspace/DSR integration remain unexecuted.

The observer describes settled prescribed harmonics. Switch-on transients,
thermal energy, parameter work, volume integration, mechanical force and a
complete enclosure balance are not inferred. Original source/material/spatial
uncertainty and arbitrary-null relative accuracy remain unresolved.

## Direct sparse electrodes with prescribed three-dimensional motion

`SpatialConduction::solve_terminals_with_motion` adds the original motional
weak load to the same grouped sparse electrode equations. One actual-amplitude
field solve closes voltage, conserved total-current or mixed/floating controls;
the original conductor-frame current, contact and Lorentz recovery is reused.
Workspace admission also charges the native load, per-cell motional field and
compensated group reduction before preparation.

Retained terminal voltage differences pair with authored current controls and
native voltage-source reactions. Independent terminal work must agree with both
native boundary work and bulk/contact Joule heat plus signed mechanical output.
Open circuits may have nonzero motional voltage and zero heat; generation may
return negative electrical power. Nonzero authored current still requires
representable irreversible heat. Original fixed-occupancy quasistatic limits
and refusal of moving/sliding resistive interfaces remain binding. Three
regression groups in `tests/spatial_conduction/motion/direct_terminals.rs`
cover open/short/generating/motoring states, floating shared electrodes,
stationary contacts with remote motion, huge bias, bounds and cancellation.

## Sampled nonlinear electrical contact laws on planar and revolved traces

`PlanarConductionResponse::with_contact_quadrature` and the axisymmetric facade
accept explicit `ContactResponseQuadrature` resolution and a total point budget.
Positive composite Simpson points include both original trace endpoints. Native
assembly, normal current and contact heat use those same source samples and
physical surface weights, including the exact radial measure for revolution.
The existing exact constant-coefficient stencil and default uniform-trace
admission retain their separate behavior.

Every sample closes the original interface law at its own actual voltage jump
and normal current density. Applied/evaluated full source states, pointwise
R/V/J defects and summed absolute constitutive-work defects remain observable.
Adaptive and fixed relaxation requery every point's owned coordinates; all
actual and damped source evaluations are counted. The point budget is admitted
before sample allocation. Each memoryless solve restarts from its explicitly
refreshed contact seed; that single contact source is a reference state, never
a face-average nonlinear resistance. The complete quadrature policy travels
with continuing frozen-source handoffs.

Five unit tests of the shared production quadrature kernel cover constant-law
planar/radial moments, positive degree-three triangular rules, nonlinear
integration refinement, operator/work agreement and cancellation/refusal.
Native regression definitions cover an analytic nonuniform continuum field,
both signs, multiple resolutions, original source states, relaxation, budgets,
domain refusal and continuation. Resolution remains caller-owned: a fixed
quadrature rule is not an error certificate or empirical material validation.

## Profile-aware EM temporal field products

Normalized Uniform and Triangular detector windows now use their actual shared
positive temporal field factor in interior E/H populations, local reservoir heat,
signed Poynting/storage rates, retained point-heat observations, volume/rectangle
heat observations, and exterior Maxwell stress and field momentum. Proper factors
retain the physical one-half law and original coherence; real fields retain their
carrier, spatial quadratures and DD source epoch. Ineligible bandwidths or tone
counts retain the existing pairwise fallback. General exterior E/H/Poynting and
finite-aperture joint kernels remain separate owners.

A triangular observation is the normalized weighted mean of the physical
quantity. It is not an unweighted interval integral or an endpoint storage
difference divided by duration. Existing physical heat intervals remain Uniform.
The shared profile preparation exposes its fixed scratch payload (6,872 bytes on
x86_64); the exterior stress workspace includes it. Existing volume fixed-stack
exclusions and original source/mode/heap bounds remain explicit.

Validation: 34 focused production numerical tests pass with warnings denied:
18 interior/retained-heat tests, eight exterior stress/momentum tests and eight
integrated-heat tests. They include independent weighted field/ensemble oracles,
deep real/proper nulls, signed storage and momentum rates, dispersive reservoir
transfer, translated rectangles and continuous spatial/carrier products. An
old-source 1 ps triangular null probe reproduced a real result about 3.3e10 times
too large and an incorrectly zero proper result. Validation used actual numerical
leaf sources with narrow outer metadata adapters; the full material/field/FFT
pipeline and native workspace suite were not run in this environment.

## Nonuniform nonlinear contacts on original tetrahedral faces

`SpatialConductionResponse::with_contact_quadrature` uses the same source-owned
sampling policy and constitutive closure on matching original triangular faces.
Each subdivision uses positive vertex/midpoint/centroid weights exact through
degree three. Integer barycentric keys merge shared points deterministically;
`3*n*n+3*n+1` points per face are charged against the global contact budget.

The original triangular contact operator, current and heat share those same
point resistances and physical area weights. Endpoint densities are recovered
from actual endpoint samples; the public sample report retains both applied
and evaluated complete original source states. The original constant-R''
triangle mass matrix remains a separate exact path. All-point relaxation,
source-domain refusal, source-query counts, work acceptance and cancellation
reuse the common nonlinear contact owner. The policy also survives native
frozen handoffs consumed by existing finite-thermal and circuit couplers.

Two focused native regression groups cover spatially varying current-dependent
contacts at multiple resolutions and external temperatures, an independent
pointwise quadratic current law, source identity, native bulk/contact/source
work, exact query counts, budgets, cancellation and continuing source policy.
The rule remains a caller-refined discretization, not a quadrature-error or
measured-material certificate. Moving/sliding contact surfaces and interface
electrical memory remain outside this memoryless fixed-geometry law.

## Sampled-volume reversible film energy and stable storage changes

The original prepared exposure layer now exposes `volume_energy` on explicit
`AngularSpectrumEnergyVolumeSample` quadrature: original grid index, depth in
the selected film and nonnegative physical volume in cubic metres. The retained
owner reports electromagnetic, polarization restoring, polarization kinetic and
total energy in joules, or signed analytic endpoint changes on the original DD
source clock. It shares the point observer's original normalized material
coordinates, real/proper ensemble law and binary-scaled field products.

Preparation groups positive-weight samples deterministically by depth, sharing
each original source-bank frequency reconstruction across all selected points
and reservoirs at that depth. Limits admit authored samples (including zero
weights, hard cap65,536), active depths, original line/FFT work, retained coordinate
terms and explicit workspace before reconstruction. Every authored entry is
validated. The original layer and quadrature remain borrowed; repeated queries
allocate no heap and perform no field/material solve. Volume and constitutive
scales meet original field products before restoration to DD range. Analytic
field increments retain finite changes when bright endpoint energies overflow.

This is sampled physical volume quadrature, not a certificate of spatial
resolution, a continuous Fourier-cell integral, an enclosing force surface,
thermal storage or switch-on dynamics. Original prescribed settled fields and
material ownership remain authoritative.

Validation: five new actual numerical volume tests and seven existing point
energy tests pass. Independent mixed-source ensemble and nonuniform volume
checks cover deep real/proper cancellation, every energy account, extreme
physical scales, signed spatial exchange, precise clocks and partitioned changes.
Four additional narrow adapter checks execute the complete preparation body
against manufactured original frequency fields: batching, point/volume agreement,
all work/memory refusals before reconstruction, zero-weight validation, and late
cancellation/reuse. The adapter preserves the callable owner shape but does not
execute the native Maxwell/material/FFT pipeline. That full pipeline and the
native workspace suite remain unrun here.

## Positive original exterior field observations

Deterministic and ensemble exterior field observations now optionally project
the original E/H vectors before all six field populations and all three signed
Poynting components. Uniform and Triangular normalized windows retain their
actual sensitivity, real carriers and proper half-factors. Continuous cell
watts retain the original opposite-mode, unpaired-Nyquist and nonzero-carrier
laws. Physical area multiplies the completed compensated spatial sum; spatial
cosine/sine fluxes combine before their half-weight, retaining finite aggregate
subnormal power. No material law or finite-aperture model is added.

Preparation reserves optional4*N*(F+1)*17 vector projection terms and50,008 bytes
of fixed logical payload from the existing authored work/workspace headroom.
The ensemble also budgets one retained admission flag for every original bank
when any bank is selected. Banks are selected independently in original order.
Exact headroom admits; insufficient headroom and unsupported spans keep the
original pair law. The original common-owner workspace/pair getters retain their
layer/stress/force meaning; new field_workspace_bytes()/field_projection_terms()
report the additional field-query reservation. Their applicable totals fit the
same original limits. Scalar locals, compiler temporaries and allocator/stack
overhead remain excluded from logical payload accounting.

Admitted banks reuse each original frequency reconstruction and synthesis once,
allocate no extra query heap, and skip the quadratic temporal-kernel table.
The numerical improvement is therefore conditional on field projection admission;
tight legacy budgets retain the original cancellation limits.

Validation: nine committed actual-helper regressions pass with warnings denied,
including independent weighted vector/spatial oracles, real/proper deep nulls,
signed backflow, exact/one-below admission, independent mixed banks, cancellation,
and both spatial subnormal aggregation cases. Three additional scratch checks
compile and execute the verbatim deterministic/ensemble dispatch bodies with
counted manufactured frequency fields, confirming reuse and fallback. Actual
math/cell helpers were executed; native material/FFT constructors and the full
workspace suite were not. The finite-aperture observer remains a separate owner.

## Two-tone field and heat interference

The film positive-field admission and optional exterior field reservation now
include two real or proper source frequencies. Retained heat uses the same film
gate; stress/momentum and volume heat delegate to shared two-tone admission.
Existing maximum capacities, full real-carrier span, source independence,
physical Uniform heat intervals and single-tone fallback remain unchanged.
Film preparation and optical/shared callers admit their corresponding projection
work; tight exterior headroom retains the original pair path.

Eight focused film heat/flux tests and ten exterior field-helper tests pass
against the shared two-tone core with warnings denied. Added cases resolve
cos(t)-cos(2t) fourth-order real leakage and second-order proper beats,
including signed Poynting and a2^49 rad/s common proper carrier. The exterior
two-tone one-point reservation admits exactly204 terms and falls back at203.
This is actual numerical-helper validation; full native material/FFT constructors
and the workspace suite remain unrun.

## Adaptive integration of original nonlinear electrical contacts

`with_contact_quadrature_refinement` closes the original bulk and interface laws
at each successively doubled contact resolution on planar, revolved and original
tetrahedral traces. Every rung starts from the same supplied material/contact
context and the external driver's same accepted history. The finer candidate is
accepted only after consecutive complete solves agree in physical edge voltage
drops, original contact endpoint jumps, individual contact and boundary currents,
bulk/contact and total heat, and current-density components. Separate comparisons
prevent opposing local errors from cancelling. Physical edge drops recovered
from native electric fields remain independent of absolute voltage references.

The common policy bounds subdivision refinements, total contact points, frozen
field trials and original-source queries. Before a rung, its iteration allowance
is narrowed by the remaining global work and its full worst-case queries are
reserved, including relaxed-coordinate reevaluation. This conservative source
reservation may refuse a solve that could converge more cheaply. Successful
reports charge actual work across all discarded and accepted rungs. For the
planar and axisymmetric response owners, the original `maximum_iterations`
also caps total frozen-field trials across refinement, preserving the remaining
budget passed by thermal owners. The tetrahedral `maximum_field_solves` remains
global, including every basis solve for an explicitly requested terminal matrix.
Callback-internal work retains its separate driver bounds.

The accepted native operator retains its actual finer stencil; continuing source
policy retains the declared initial resolution. Frozen and thermal handoffs
therefore repeat the same refinement procedure from their refreshed physical
state. Fixed quadrature explicitly clears automatic refinement. Cancellation,
source failure, unresolved differences or work/point exhaustion returns no
accepted candidate, and refinement differences never enter physical heat.

Three executed shared-helper tests cover componentwise acceptance, conservative
admission with actual accounting, invalid observations and cancellation. Planar
integration regressions cover an exact nonuniform field, explicit sums of
discarded-rung work, continuation and refusals, plus a nonlinear contact whose
current and heat are independently integrated analytically. Tetrahedral
regressions compare the accepted native field and actual cumulative work against
successive fixed-resolution solves, preserve source policy through continuation,
and exercise unresolved, point/work-budget and late-cancellation refusals.
The integration regressions have been authored and formatted, but not compiled
or executed in this environment. Resolution convergence is an estimate, not a
quadrature, mesh-error, uniqueness or empirical material certificate.


### Field-preserving finite-aperture temporal interference

Finite-aperture normalized-window moments now optionally project each original E/H source bank into bounded positive temporal modes before continuous rectangle products. Both Uniform and Triangular shutters retain actual real carrier sums and proper-complex difference statistics, original joint center phases, all spatial modes and Nyquist terms, signed Poynting, Jones-before-product analyses and Maxwell stress. M<=17, source limits 2..32 real or 2..64 proper and half-band<=0.5 remain shared-core admission. Additional logical payload is 31,504+192*N*M bytes; work reserves max(N²F²,N*F*M+N²*M) per selected bank. Tight existing budgets and unsupported bands keep the original pair law. The separate wrench, interferometer and energy-control owners are not covered by this generic measurement change.

Focused verification passed 13 exact numerical tests and two full dispatch-body tests with manufactured reconstruction adapters under -Dwarnings. These cover real/proper dark residuals, both profiles, independent field quadrature, continuous rectangle partitions including unpaired Nyquist modes, exact budget boundaries, cancellation and single reconstruction per original frequency. The centered 1ps real two-tone fixture previously returned 9.183549615799121e-41 instead of the physical leading value 2.8125e-50; the projected route preserves that residual. This is focused numerical and dispatch validation, without a native Maxwell/material/FFT reconstruction or full workspace claim. It does not certify relative error at arbitrary joint spatial nulls.


### Nonlinear conductor-frame response under prescribed motion

The original planar and tetrahedral scalar conductivity/resistivity owners now
close material response at the actual conductor-frame field E' = E + u cross B
and native J. Laboratory E remains a separate observation. Constitutive current
and heat defects use the same E' as the material source; complete source pins,
external coordinates, source-coordinate relaxation, stationary-contact response,
adaptive contact quadrature, cancellation and cumulative field/source work
limits retain their existing owners. No accepted nonlinear operating point is
rescaled, and the private native adapters move existing buffers without cloning.

Planar boundary voltage and distributed-current traces use
conduction::response::motion::MovingConductionResponseCandidate. Its mandatory
field/motion report keeps native Joule heat, signed bulk laboratory electrical
work, stationary-contact heat, signed terminal supply and Lorentz mechanical
output distinct. The 3D response::motion entry points accept original boundary
traces or direct grouped voltage/current/floating/mixed electrodes, returning
the original source owner with SpatialMovingConductionSnapshot. Each direct
electrode constitutive trial uses one actual-amplitude sparse solve.

Positive mechanical work leaves the electrical/thermal system. Native terminal
input equals irreversible Joule heat plus signed mechanical output; generation
may have negative electrical supply, and short-circuit braking may heat with
zero terminal supply. Adaptive contact comparison additionally observes local
and total signed mechanical work. Contacts remain stationary surface laws:
both original adjacent cells must have zero velocity, including when u cross B
vanishes. Motion in remote bulk cells is admitted.

These are prescribed velocity/B snapshots on fixed geometry and material
occupancy. They add no induction evolution, mesh advection, mechanical state
integration, sliding-interface work law or regulated-power moving policy.
Axisymmetric motion and regulated-power moving drives are not exposed by
these entry points. Positive source values do not establish uniqueness or
convergence.

Ten focused native regressions cover original E'/J source coordinates, motoring,
generation, braking, open-circuit screening, zero-motion parity, direct 3D
current/floating/mixed electrodes, planar current traces under a large common
voltage gauge, stationary contact heat and moving-side refusal, source-domain
refusal, budgets, cancellation and retry. Pinned rustfmt parsed/formatted the
changed Rust files. Bounded rustc verification exercised the actual 3D snapshot
getters/observation adapter with unrelated opaque result types substituted and
independent uniform-bar roots; separate analytic signed-work references also
passed. The new native integration regressions have not been compiled or run
in this environment; these checks are not a full crate/workspace validation.


### Interference-preserving continuous aperture traction and torque

The separate finite-aperture wrench observer now shares bounded positive temporal projection before applying its original signed first-spatial-moment kernels. Continuous rectangle lever arms, real/proper source statistics, Uniform/Triangular profiles and torque translation to the requested origin remain intact. Unsupported bands or insufficient extra headroom retain the original pair route, with the same optional temporal work and workspace admission as ordinary aperture moments. The shared linear temporal projection is reused, with each original frequency reconstructed once.

Combined focused verification passed 26 checks: 13 ordinary-aperture numerical tests, ten wrench numerical tests and three complete dispatch-body checks. New cases include analytic two-tone force, an optical three-tone dark fringe with signed torque and translated origin, all six traction/moment components against resolved mixed-polarization products, cancellation and tight-budget fallback. Dispatch checks use nominal reconstruction and medium adapters; no native Maxwell/material solve, complete-body force/torque closure or universal joint-null relative accuracy is claimed.


### Coherent finite-aperture energy controls, gradients and Hessian actions

Prepared aperture energy controls now optionally project the complete coherent field into the shared real temporal basis before applying the original continuous spatial Gram. Energy is contracted within each temporal mode, and each original control bank is coherently projected across all its frequencies before the source pullback. Consequently the value, Cartesian complex-source gradient and matrix-free Hessian action retain the same resolved short-window interference. Fixed illumination remains excluded from Hessian directions, empty controls retain original order, and temporal odd parity retains the real carrier-sum quadrature coupling. No repeated query allocation, Maxwell reconstruction or dense source Hessian is introduced.

For admitted Uniform/Triangular bands the extra work bound is N*M*(F+width+height+1+source_groups+controls), with M<=17. Additional retained payload comprises N*F combined DD fields, N projected fields, F*M DD coefficients and one DD complex gradient per control; preparation scratch is included in admission. Existing pair kernels and workspace remain reserved. Insufficient extra headroom or an unsupported band uses the original pair operator. Sources still carry their original joint aperture/clock center phase, material law and physical units.

Eight actual numerical control tests and one complete constructor/query dispatch check passed, together with the aperture/wrench suite (35 combined checks, -Dwarnings). Independent space/time quadrature covers both the admitted and fallback widths under both profiles; tests also cover Hessian symmetry/finite differences, coherent empty-control ordering, cancellation/replay and exact/one-below optional budgets. The original control operator returns 5.739718509874451e-41 J/m³ for the centered 1ps two-tone fixture whose physical energy is 2.8125e-50; the projected value, gradient and curvature retain that residual. Dispatch validation uses nominal source-reconstruction and medium adapters, without a native Maxwell/material solve or full workspace build. Spatial Gram conditioning and arbitrary deeper nulls are not assigned a universal relative-error certificate.


### Cached temporal interference in prepared aperture interferometers

Prepared aperture interferometers now cache eligible positive real/proper temporal factors and reuse an Mmax*N projected-field buffer with the original Gram workspace. Every read coherently combines complex arm gains and signed delays before frequency projection, then applies the exact separable spatial Gram per temporal mode. All 21 real field moments, signed Poynting, original source clocks and independent-bank statistics remain represented. Repeated reads allocate no heap memory and perform no reconstruction; equal-depth arm preparation continues to share each original frequency reconstruction.

Both deterministic and ensemble owners now expose aperture_interferometer_window_requirements(arms, window), returning enough optional headroom for every eligible bank. The existing arms-only requirements and preparation signature remain available. Per-bank work reserves max(original NF*(Nx+Ny+F), N*M*(F+Nx+Ny+1)+proper*N*F). Extra payload is sum(16*F*(M+2))+192*N*Mmax+12,024 bytes, with existing owner/bank headers counted separately. Unsupported bands and insufficient optional headroom retain the original temporal pair route.

Five new numerical/readout tests and two complete constructor/query checks passed within the 33-check aperture/wrench/interferometer set under -Dwarnings. They exercise real/proper dark fields, all 21 resolved moments, complex gains and signed delays, independent banks, exact window-specific limits and one-below fallbacks, shared-depth reconstruction counts and cancellation/replay. Preparation/dispatch checks use nominal reconstruction owners; native Maxwell/material/FFT reconstruction and a full workspace build were not run. Continuous spatial conditioning retains its existing limitations.


### Shared carrier phase for delayed proper aperture arms

For admitted proper-complex interferometer banks, arm combination now references the first active arm's common carrier delay and composes the retained DD relative-frequency phases. Removing that one global carrier phase leaves every proper second moment unchanged while avoiding independent bright-phase rounding before a dark-envelope projection. Real banks retain their absolute source-clock law. If a delay difference or its carrier product is unrepresentable, the original independent phase path remains available. The original arm reduction order is preserved, no query heap allocation is added, and optional proper work reserves an additional A*(F+1) phase/composition terms.

Three new regression/budget/fallback checks pass within the 36-check aperture/wrench/interferometer set under -Dwarnings. For proper tones 1e15+[0,1,2] rad/s with amplitudes [1,-2,1], centered 1ps Uniform exposure and delay 8e-16s, the old result 6.250162149376229e-51 differs from the independent shifted-fourth-moment value 6.250160000204799e-51 by 3.438586e-7 relative. The repaired unit/coherent-real-gain fixtures agree within 4.45e-16 across both profiles and three delays. Tests also retain nonfirst-active-arm behavior, resolved complex-gain moments, cancellation and previously valid extreme-delay fallback. DD storage still limits arbitrary combined complex-gain, temporal and spatial null conditioning; no universal relative-error claim or native Maxwell/material reconstruction validation is added.


### Direct planar moving electrodes and nonlinear source response

Planar prescribed-motion conduction admits direct simultaneous ideal voltage,
conserved total-current, floating and mixed electrodes through the original
grouped sparse operator. Each actual constitutive trial groups the native
motional weak load with the stiffness and performs one solve. Terminal current
is an integrated reaction; a zero-current electrode may transfer current
between separate traces. Ordered observations preserve small physical voltage
drops beneath a common bias. Independent authored terminal work must equal
conductor-frame Joule/contact heat plus signed Lorentz output, including
motional open circuit and electrical generation.

The existing original-source closure, relaxation and contact refinement retain
their source pins, cancellation, local gates and cumulative native solve cap.
The extra motion/terminal workspace is admitted before allocation. Geometry,
velocity and imposed normal B remain held, and electrical contact-adjacent cells
remain stationary. This adds no induction evolution, moving-interface work,
axisymmetric-motion law or regulated-power policy. Focused native regressions
cover terminal modes, floating transfer, signed work, source coordinates and
transactional refusals. Pinned formatting, bounded exact-production type checks
with dependency stubs and independent reference arithmetic passed; full crate
integration execution remains outstanding.


### Original nonlinear moving-current loads for native mechanics

SpatialMovingConductionResponseCandidate exposes lorentz_loads and
lorentz_loads_for_velocity using its retained original mesh and accepted native
moving snapshot. The same cell P1 integration, compensated scatter, canonical
frame/origin wrench map and local nodal-mean velocity check as the borrowed
fixed-source field remain authoritative. Extraction performs no electrical
solve or material/force re-evaluation and adds no overlapping Maxwell load.
Selection, geometry, cancellation and workspace bounds remain active. The
candidate owns the source/field pairing; an arbitrary external snapshot cannot
be substituted into this handoff. Physical mechanical feedback belongs to the
coupling owner, and heat remains the original conductor-frame dissipation.
