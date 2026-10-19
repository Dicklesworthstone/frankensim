# Loaded solids with spatial thermoplastic feedback

`fs_couple::solid_enthalpy::plastic::ThermoplasticSolid` connects the existing
three-dimensional J2 tetrahedral equilibrium solve to the original spatial
`DepositionSimulation`. Material data changes the actual loaded specimen:
temperature-dependent yield changes deformation and support reactions, while
physical plastic dissipation heats the assigned thermal cells. Conductive links
and finite thermal bodies then exchange that heat through the existing graph.

## The shared owners

Assign one `ThermoplasticCell` to each tetrahedron in connectivity order. Its
`SourceQuery` identifies the material card, original source-derived caloric
curve, complete external context, selected hardening convention and strain
bound. The original caloric object must also be registered with the thermal
simulation. Native tetrahedral volumes and sourced density must match its cell
volume and delivered mass. Different elements can use different material cards.

`ThermoplasticSource` shares source admission with the existing uniform
`CaloricPoint`; it owns no alternate thermal inventory or plastic history.
Every temperature query retains all six selected roles: E, Poisson ratio,
yield stress, hardening datum, Cp and density. Native `plastic3::solve_from`
continues the accepted displacement and plastic history through the same
Newton, tangent, load and reaction operators as ordinary load-path solves.

For each physical interval, `advance` accepts endpoint prescribed displacements
and dead nodal forces. Use `geometry().assemble(...)` for existing pressure,
traction and body-force integration, then supply its `nodal_forces_n()`.
The supplied thermal observer composes its existing sources, conductances,
conservative transfers, event handling and time-step limit with plastic heat.

Each coupled proposal restarts from the accepted solid, thermal graph and
observer. The converged candidate provides the actual final-temperature
material fields, native displacement/stress/plastic history, thermal state,
observer and account together. `into_parts()` transfers these accepted owners.
Cancellation, source exit, failed equilibrium, kinematic limits or energy
refusal leave all input owners unchanged.

## Energy and numerical meaning

The energetic hardening convention stores
`psi = 1/2 elastic_strain:C:elastic_strain + 1/2 H alpha^2`.
Only `yield_stress * delta_alpha * reference_volume` is physical plastic heat.
Hardening storage and backward-Euler numerical loss remain distinct.

The coupling compares assembled nodal-force/support-reaction work with the
native element work, then checks the complete interval:

`delta_mechanical_storage + delta_caloric_storage + numerical_loss
 = endpoint_mechanical_work + external_heat + incoming_enthalpy`.

Plastic heat appears in the original thermal-only audit as an input. The
coupled report removes that internal transfer from `external_heat_j`, so it is
not counted twice. The physical heat and the rate supplied to thermal substeps
must agree within an independent joule tolerance.

## Admitted model and remaining boundaries

This is bounded quasistatic, rate-independent, small-strain J2 plasticity with
zero thermal eigenstrain, constant E/Poisson ratio/H/reference density and
single-solid, increasing caloric charts. Yield and Cp may vary with temperature.
The full displacement-gradient norm is bounded by the declared strain limit,
so large rotations cannot bypass a symmetric-strain check. Mechanical and
thermal reference geometry remain fixed. There is no inferred plastic-heat
fraction from an empirical hardening curve.

Each physical interval has one backward-Euler mechanical increment; its heat
is spread uniformly over the actual thermal substeps. Refine both load
intervals and thermal steps. Fixed-point convergence does not establish time
accuracy, mesh convergence or unconditional stability. Variable elastic free
energy, expansion, inertia, viscosity/creep, finite-strain plasticity, moving
thermal contacts and solid/melt support require their own constitutive terms.
This path refuses phase/kinematic exits rather than inventing melt stiffness.

Focused regressions in `solid_plastic_enthalpy.rs` cover an independent
force-loaded uniaxial oracle, material-only Cp changes to temperature and shape,
load/unload history, six-tet localized heating and grip-force relaxation,
conductive temperature redistribution, and atomic late refusal. The numerical
fixtures use synthetic material data; they are not experimental validation of
a named engineering alloy.
