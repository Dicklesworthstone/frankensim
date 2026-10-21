# Material-derived conduction in layered self-heating

`fs_couple::strip_enthalpy::transport` connects the existing prepared material
network to the original compatible Maxwell/section/enthalpy solve. Conductivity
is re-evaluated at each proposed thermal endpoint, not copied from the initial
temperature or inferred from a material name. The prepared network keeps its
original claim selection throughout the nonlinear solve and across accepted
steps for which the material/topology bindings remain valid.

## Existing owners, one feedback path

- `MaterialNetwork::prepare` selects the original transport claims.
- `PreparedMaterialNetwork::resolve` evaluates those same claims at each trial.
- `transport::coefficients` checks mechanical/transport source compatibility and
  passes the original conductances and source snapshots to the thermal residual.
- The original strip solve deposits physical Maxwell heat once; the original
  thermal network moves and stores that heat. Temperature-dependent relaxation
  then changes compatible deformation.

The mechanical and conductive properties must share card identity, fixed
coordinates, coordinate conventions and the mechanical temperature domain.
`ResolvedThermoMaxwell::admit_state_point` owns this reusable cross-property check.
It does not reselect the other properties, replace their laws or imply experimental
validation. Numerical damping is still separate from physical heating.

## Usage

Build the existing `MaterialNetwork` with complete `MaterialCell` and
`ConductionLink` bindings, using the same caloric objects as the thermal timeline.
For this thermal/mechanical path its `resistors` list must be empty. Electrical
branches are rejected instead of dropping their work or inventing Joule heat.

```rust,ignore
let transport = network.prepare(&thermal, cancelled)?;
let candidate = fs_couple::strip_enthalpy::transport::step_candidate(
    &thermal, &strip, &layer_bindings, &transport,
    target_s, control, config, &energy_guesses_j, cancelled,
)?;
let (thermal_next, strip_next, report) = candidate.into_parts();
```

For independent extra heaters or cooling sources, use the existing generic
`strip_enthalpy::step_candidate` and call `transport::coefficients` inside its
coefficient callback. Add only those independent watts to `sources_w`; never add
Maxwell loss again. The report's source retains the original complete material
states and actual endpoint conductances. Authored link coefficients are restored
in accepted thermal states; they are not mutated into hidden persistent inputs.

`GrowingStrip::step_transport_candidate` uses the same coefficients while retaining
the existing atomic interval/birth semantics. Its prepared network describes the
pre-event interval. After new cells or contacts are born, explicitly prepare a
complete network for the next interval. The old topology refuses rather than
ignoring a new cell or contact. Use the existing prepared network's
`extend_thermal(&successor_network, born.thermal(), cancelled)` to preserve old
conductivity claims while admitting new cells and newly requested conductivity.
Ordinary `prepare` remains an explicit new selection epoch.
There is no inferred contact or newborn thermal/mechanical aging in Newton trials.

### Source continuity during growth

`PreparedMaterialNetwork::extend_thermal` accepts a complete successor network
after the original timeline has applied births or contact changes. Existing
cells retain their material card, actual caloric object, material ID, fixed
conditions, conductivity key and selected claim. Another claim becoming eligible
at the current temperature cannot change those cells' constitutive model.
An expired original source refuses rather than falling back to a different fit.

Only new cells or previously unrequested conductivity use the successor's
selection policies. Existing pins take precedence by this method's contract;
use `prepare` for an intentional source-selection change. Links and distances
remain explicit successor data. No electrical state is migrated: old or new
resistor branches refuse. No physical history is modified or cached, and late
cancellation leaves the previous prepared network usable.

Five additional integration regressions exercise actual birth and further
self-heating while competing hot-source claims appear, first contact on isolated
cells, expired-source refusal, fixed-context/electrical refusal, and cancellation
with replay. These Rust tests have not been executed in the editing environment.

## Verification and boundary

Six focused Rust regressions cover original-path equivalence, actual endpoint
k(T), conductivity-only changes to both temperature and curvature, same-card
context mismatch, source-domain refusal, cancellation, electrical-branch refusal,
and a later material birth with explicit conductive topology update.

```sh
rch exec -- cargo test -p fs-couple --test strip_material_transport \
  --test strip_enthalpy --test strip_enthalpy_growth
```

Native Rust compilation and tests were not executed in the editing environment:
RCH/Cargo/rustc are absent. Independent NumPy/SciPy reference equations passed 81
conductivity/temperature/duration cases, with maximum energy discrepancy
2.77e-15 J. In one synthetic two-layer case, changing only k from 0.0002 to
0.02 W/(m K) changes curvature from 0.0339724564 to 0.0177450171 per metre.
These are independent equations, not Rust output or measured-material validation.

The original uniform-layer, perfectly bonded, small-strain and fixed-thermal-
geometry assumptions remain. Density, spring energy and free thermal strain must
remain temperature-independent for this caloric coupling. Conductivity and
relaxation may vary within their source support. Scalar series conduction is not
anisotropic heat flow, melting, geometric-contact evolution or a spatial bridge.
