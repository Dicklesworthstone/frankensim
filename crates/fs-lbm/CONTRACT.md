

## Compact borrowed free-surface impulse staging

`FreeSurface3::prepare_impulse_stencil` retains a mutable borrow of the original native owner while sampling a bounded, strictly increasing wet-cell stencil. Native populations, located retained liquid mass, physical Guo-frame velocity and original cell centers feed the shared population impulse kernel. Gas, inactive/static solid, piston, unresolved interface inventory, excessive speed and unrepresentable states refuse. Staging preserves original populations until the final commit; dropping a stage/candidate or final cancellation leaves the original unchanged. Candidate readback exposes each original mass/position and actual endpoint velocity so a finite receiving owner can independently gate reciprocal angular momentum and work before joint publication. No fluid clock, held force, fill/classification, geometry or unlocated carry is changed.

`MaterialFreeSurface3::prepare_native_impulse_stencil` reuses original physical impulse admission and scenario speed, while explicitly returning native-unit staging for a composing owner’s SI transaction. It does not silently reinterpret lattice readback as physical units.

The native stencil's three focused tests execute actual retained populations, force-frame stress, readback and transactional refusal. A further native refinement regression checks retained solid keys outside the physical domain and beyond representable child Morton coordinates: refinement refuses before allocation, preserving source mass, geometry and populations. All six native refinement cases pass after this guard. The prior three SI receipt checks remain unchanged; full streaming/material integration was not executed.

## Temperature-activated molecular wetting

`MolecularKineticLaw::try_new_arrhenius` optionally evolves the equilibrium jump
frequency as `K(T)=K_ref exp[(E_a/k_B)(1/T_ref-1/T)]` inside the original calibrated
temperature domain. Activation energy is explicit joules per molecular event;
the frequency in `MolecularKineticControls` is its value at the declared reference
temperature. This constant-prefactor Arrhenius closure adds no inferred adhesion
measurement, Eyring prefactor, extrapolation or activation-heat source. The
original constructor still selects constant frequency. Reference-temperature and
zero-activation responses preserve that original model exactly.

The same law feeds the existing dense, sparse and enthalpy wall-angle owners,
retaining original interface claims and actual local temperature. Accepted samples
report the local jump frequency separately from static material measurements.
The rate and velocity scales must be finite and positive across both temperature
bounds; thermal domain refusal retains native flow, enthalpy and wetting history.
Per-length molecular loss remains diagnostic, without claiming a resolved line
measure, a global surface-energy balance or heat deposition. G1/G3 unit cases
check independent forward kinetics, low-speed friction and constant-model limits;
the existing thermal integration case also exercises activated rates and rollback.

## Isothermal kinetics from reversible equilibrium adhesion

The dense D2Q9 and sparse D3Q19 material owners can attach an original equilibrium
angle or Young--Dupre work with `with_equilibrium_planar_wetting`, then explicitly
select Cox or molecular kinetics with `with_dynamic_equilibrium_wetting`. The
same frozen liquid tension feeds both adhesion and native pressure. Static
coefficients sampled at fixed external conditions retain their source axes and
receipts; time, loading and rate restrictions cannot be silently frozen. No
advancing/receding measurements or contact-speed query are fabricated from work.

`with_equilibrium_wetting_hysteresis` optionally binds independent static pinning
thresholds from that interface, preserving the equilibrium potential and separate
resistance receipt. The existing full-angle rate domains, geometric speed history
and transactional native publication still apply. Pure equilibrium alone bypasses
the velocity quotient when fill changes on a locally flat wall stencil; adding
kinetics or pinning retains the unresolved-speed refusal. The backward geometric
speed drives the next step's capillary ghosts. Line losses remain diagnostics;
these owners do not infer a line measure, heat deposition or closed energy balance.

Four public integration groups exercise both source origins and rate laws in
both dimensions, independent forward equations and native response, separately
sourced pinning, domain rollback, and the pressure-loaded flat-fill distinction.

## Heated wetting at neutral transverse corners

The height-graph equilibrium-wetting owner now closes the two ends of its
sourced side wall against the explicitly selected neutral Cartesian transverse
walls. The simultaneous conditions `n.n_wall=cos(theta)` and
`n.n_transverse=0` give zero transverse height slope and the existing sourced
normal slope. Graph normals, area, Marangoni projection and public heat patches
share this metric. The original raw/smoothed VOF pressure ghost uses even
reflection at that same neutral wall, including diagonal corner samples.
The internal boundary selection survives thermal angle-field replacement;
ordinary planar wetting retains its original stencil.

This closure requires one sourced stationary wall, no internal solids or piston,
complete monotone liquid/gas columns, one angle through each occupied interface
band, and the original graph-slope envelope. It does not supply an angle for a
second substrate, arbitrary wedges, corner films, corner pinning, capillary heat
or a closed surface-energy balance. Invalid source work, geometry and cancellation
still refuse the complete thermal/momentum candidate. G1 tests cover both
graph orientations on all Cartesian wall pairs, raw affine-plane preservation,
zero artificial Laplace pressure, independent corner traction and an actual
heated Young--Dupre corner response through the public material owner.

## Material-flow enthalpy reconstruction

`MaterialFlowAdvection::with_reconstruction` attaches the original thermal owner's shared `LimitedEnthalpyReconstruction` to solved material face fluxes. Selected contiguous axial cells reconstruct total specific enthalpy, including a declared passive equal-density latent inventory; paired face transfers, source-law inversion, accepted history and the combined conduction/advection timestep restriction remain with `fs-thermal`. Unselected cells and external streams retain donor enthalpy. Empty plans retain the donor-cell arithmetic. The observed timeline retains first-order time integration on stationary control volumes; it does not infer multidimensional junction slopes, phase mobility, density-change flow, pressure work or moving interfaces. The `advection()` accessor exposes the underlying graph and actual transport audit; implicit clients must attach the same `reconstruction()` plan to their own native implicit transport owner.

Pure frozen `MaterialFlowAdvection<_, ()>` owners additionally expose native
`advance_ssprk2_to` and `advance_ssprk2_adaptive_to`, including corresponding
`with_coefficients` entries for repeatable stage sources and link/bath
conductances. These methods retain the solved `rho * face_volume_flux` graph,
original material/density snapshot and selected reconstruction in the returned
accepted wrapper. Original fluid-cell mass, geometry and caloric-law identity
are admitted at stage queries and before publication. The native thermal owner
alone integrates storage and heat/mass accounts. Adaptive comparisons are one
full step versus two fine steps, and publish only the conservative fine state
and its ledgers. Mandatory datum-invariant h error resolves latent storage at
unchanged temperature; optional T error queries the original caloric law.
Native events, CFL/caloric admission, shared work budgets and cancellation
remain authoritative. The smooth-solution local estimate gives no spatial or
certified global accuracy bound. Stateful inner observers and live momentum
feedback are not repeatable coefficient queries and do not expose these methods.
Integration regressions use actual material-driven momentum for analytic
latent-throughflow temporal refinement, reconstruction in both flow directions,
separate supplied/stream heat, map/law admission and transactional refusal.

`MaterialFlowFeedback::with_reconstruction` retains that same shared plan across
live source viscosity, bounded Boussinesq buoyancy, mushy resistance and retained
momentum updates, including both native implicit entries. Each actual source
epoch rebuilds the original mass graph before binding the same reconstruction;
no additional momentum solve or alternate caloric law is introduced. Required
axial interior connections remain present at zero flow, allowing rest and flow
reversal. They count against the original connection budget. Dormant transverse
faces add no artificial junction, while nonzero crossflow remains visible and
refuses an unsupported axial stencil. A missing neighbor at an exact-solid lock
also refuses the candidate; external or cross-solid connectivity is not inferred.

The stationary-skeleton liquid-carrier attachment reconstructs phase-carried
liquid enthalpy, leaving solid latent inventory in original mixture storage.
Its equal-density, pure-substance, source-support and physical Stokes-work gates
still apply. The optional bulk entropy face limiter remains incompatible with
liquid-only carriage. Reconstruction can be selected before the first advance;
accepted feedback preserves it with the original material pins, physical history
and transactional candidate. Empty selection preserves donor arithmetic and CFL.

## Adaptive implicit thermoviscous flow

`MaterialFlowFeedback::advance_implicit_adaptive_to` controls time error in the
existing quasisteady Stokes/material/enthalpy endpoint solve. Each comparison
uses one original backward-Euler interval and two fine intervals; only the fine
thermal states and their actual endpoint flow, source and physical work are
accepted. The original `ViscousHeating` owner integrates both accepted intervals,
retaining wall/pressure/body supply and molecular/Darcy/Forchheimer heat without
turning momentum residuals or temporal error into heat. All coarse, fine,
rejected and finite-difference evaluations share the blueprint's one
`max_solves` budget. Individual one-step reports retain per-interval work counts.

Mandatory enthalpy error scales changes from the comparison's starting state,
so latent inventory remains resolved at fixed temperature and absolute
formation-enthalpy datums do not inflate the allowance. Optional temperature
error queries the original caloric laws. First-order step-doubling uses the
actual representable fine spans and accepts no extrapolated state. The caller's
whole-advance joule allowance bounds cumulative absolute accepted algebraic
residuals and the final compensated, datum-relative storage/heat balance;
native actuator/physical-loss power tolerance remains separate. Reports retain
the final original coupled endpoint and accumulated accepted physical accounts.

This entry retains fixed topology, stateless `()` material feedback, and the
original stationary-heating restrictions. Interior or unprocessed native events
refuse before momentum work; target-time events remain pending. The native
accepted-substeps budget counts both fine intervals. Minimum timestep, trial
budget, source/work/nonlinear refusal and cancellation publish no partial
thermal or flow prefix. The estimate controls local time error on admitted
smooth solutions, not spatial/model error or a certified global bound. Focused
integration cases compare continuous source-viscosity heating, latent mushy
mobility under shifted h datums, and shared speculative-work/transaction gates.

## Hydrostatic rest under material density-dependent gravity

`MaterialFreeSurfaceScaling::free_surface3_hydrostatic` initializes a planar
column using the original material-bound Guo force, collision, streaming and
gas-pressure reconstruction. Gravity points toward one declared Cartesian
supporting wall; zero gravity is also admitted. Complete liquid layers retain
the original stationary solid tiles and at least one layer of gas. The caller
supplies a positive absolute pressure reference and a density-excursion bound.
Capillarity and wetting must be explicitly absent.

For lattice acceleration magnitude `a`, set `h = a/(2*cs^2)`. The nearest wet
layer has `rho = 1/(1-h)`; successive deeper layers multiply by
`(1+h)/(1-h)`. This surface factor follows the native missing-gas-link rule.
A continuum pressure profile evaluated at cell centers alone does not supply
the same discrete rest state. Layer densities are computed once and reused
while the original constructor owns liquid mass, sparse activation and the
physical-zero-velocity Guo frame. Ordinary `step()` retains the same gas datum.

This balances the weakly compressible numerical operator. It is not a liquid
equation of state or an incompressible mass projection: native mass contains
the retained density perturbation. Pressure and submerged traction approach
their incompressible hydrostatic limits as `dt/dx` decreases with the physical
configuration fixed; refining `dx` at fixed `dt/dx` does not remove that
compressibility error. The constructor adds no moving-solid or partially
wetted-body force closure.

G1/G3/G4 regressions in `tests/material_free_surface3/hydrostatic.rs` exercise
all six gravity directions, zero-gravity equality, pressure-datum invariance,
actual reflected-link traction and torque with decreasing compressibility,
mass retention and transactional refusal. The declared density, source-domain,
geometry and cancellation bounds remain with the original owners.

## Original-law shear-free Carreau film base

`film::carreau::CarreauFilm` resolves a positive prescribed flux per width at
fixed temperature under tangential gravity. It owns the original thermal
Carreau-Yasuda law, source-rate ceiling and accepted `CarreauSlitProfile`.
The lower half of a symmetric stationary-wall slit of gap `2h` supplies the
shear-free film exactly: `tau(y)=rho*g_t*(h-y)`. Native constitutive inversion,
positive quadrature, temperature shifting, source support and viscous heat
remain authoritative. No fitted effective viscosity enters the profile.

A bounded thickness solve brackets against the original zero-shear Nusselt
thickness and supported wall stress. Safeguarded Newton uses the exact
`dq/dh=h*gamma_wall` sensitivity and admits both flux residual and relative
thickness correction. Requested and accepted flux remain distinct. Point
samples retain apparent and differential viscosities at the same local shear;
`mu_t*u''=-rho*g_t` supplies curvature. Surface shear is zero, wall speed is
zero, and original heat per wall area equals gravitational work within the
native quadrature tolerance. Density and surface tension are explicit actual
SI inputs, not inferred from a constitutive curve or fabricated material card.

This fully developed, isothermal, planar base does not itself claim stability,
entry-region accuracy, thermal feedback or a generalized turbulent threshold.
All support/work failures and cancellation refuse without a partial film.
Focused regressions compare the exact constant-source Nusselt limit at shifted
temperature, an independently integrated nonlinear profile, gravity/heat balance,
source limits and cancellation.
