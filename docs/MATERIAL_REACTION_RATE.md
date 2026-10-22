# Material-driven reaction rates on the shared gas state

`GasReaction::with_material_rate` binds an optional source-backed forward rate
scale. The existing `parcel_candidate` evaluates it at the actual input gas
state. Existing reacting piston, nozzle/diffusion graph, material-bound inertia,
and electrical/thermal assembly paths all call that same entry point; no new
scheduler, chemistry integrator, temperature cache or gas inventory is added.
The original fixed-scale `GasReaction::new` path remains unchanged.

## Data and units

The explicit material property is `reaction_progress_rate_scale`, in coherent
SI **mol/(m3 s)**. Its dimensions come from the original thermochemical
`ReactionProgressRateQuantityV1` type. The rate law uses dimensionless partial
pressure activities `p_i/p0`; it does NOT consume an arbitrary concentration
rate constant. For example, a first-order coefficient in 1/s is refused, not
silently treated as an activity-based progress scale.

`MaterialRateScale` retains the original `MaterialCard`, explicit claim selection
and the source's declared activity reference pressure. Binding it to a reaction
with a different reference pressure refuses. The caller remains responsible for
choosing rate data for the declared reaction and activity-order convention.

```rust,ignore
use fs_couple::cavity_phs::piston::reactive::material_rate::MaterialRateScale;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Pressure;

// rate_card contains the explicit property and its valid state domain.
// reaction is the existing certified stoichiometric/caloric binding.
let rate = MaterialRateScale::new(
    rate_card,
    MaterialPropertySelection::SingleClaimOnly,
    Pressure::new(100_000.0), // must match this reaction's retained p0
)?;
let reaction = reaction.with_material_rate(rate)?;
// Pass `&reaction` to the SAME existing parcel/piston/network APIs.
```

The complete parcel `QueryPoint` supplies the source conditions. Temperature,
pressure and any additional axes therefore belong to the same physical state
used by the reaction. Curves and supported multivariate material claims are
resolved by the existing property owner, with its original selection,
interpolation, dimensions and validity checks. Missing/ambiguous sources,
unsupported coordinates, nonpositive scales and invalid units refuse. They do
not fall back to the old declared constant. To retain one specific claim across
steps, use the existing `PinnedByProperty` selection. Rebinding an already bound
reaction also refuses rather than replacing its authority silently.

The material scale is passed to the original `ThermodynamicReverseRateClosureV1`
and `DeclaredStoichiometricMassActionV1`. Forward and reverse rates thus use one
scale and one original equilibrium convention. Stoichiometry, species caloric
models, reaction mass/element checks and the entropy admission are not replaced.
`GasReactionStep::rate_material` retains the actual material snapshot; `rate`
contains the original kinetic evaluation with the effective scale. The
`kinetics()` getter remains the original declaration, not a live coefficient.

## Feedback and energy

Transport is staged before reaction on the existing flow paths. The rate card
therefore sees the actual mixed post-transport gas, once per chemical stage.
The reaction changes composition at conserved total internal energy, including
formation references, and the existing thermal/motion solver then advances.
No material rate or formation energy is credited as external heat. Nonlinear
mechanical/electrical/thermal probes do not react the inventory again.

The NEXT accepted step samples the newly evolved temperature and pressure. This
closes heat/compression -> kinetic scale -> composition/temperature -> pressure,
transport, sound and wall-work feedback through existing owners. It remains a
first-order explicit reaction split. It is not within-step implicit/stiff
kinetics, an Arrhenius fit inferred from a name, a calibrated combustion
mechanism, or a proof of temporal accuracy. Existing consumed-fraction, source,
energy and entropy gates still require timestep refinement. A zero scale is
not supported by the existing positive logarithmic rate-scale type; omitting a
reaction for an interval is an explicit control decision, not a missing-data
fallback.

## Validation status

Native compilation and the existing Rust regressions have not been executed in
the authoring environment: Cargo and rustc are absent. The existing assembly
regressions can be run with:

```sh
cargo test -p fs-couple --test material_exchange -- --test-threads=1
```

An independent Python check exercised 256 manufactured table-rate, reversible
mass-action and conserved-formation-energy cases. Maximum internal-energy
roundoff was 9.10e-13 J. This does not execute Rust or validate any real
reaction's material data.

Concentration-rate conventions and temperature/pressure parameterizations:
https://cantera.org/stable/reference/kinetics/reaction-rates.html

The implementation uses FrankenSim's dimensionless-activity convention, NOT
Cantera's concentration-based constants without conversion.
