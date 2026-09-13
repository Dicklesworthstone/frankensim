# Fan-driven cooling and speed selection

These extend the existing experimental `cooling-network` command; they are not
native `.fsim` projects or ledger-backed `solve` runs. No new dependencies are
required. The pressure/flow solve uses the existing `FanBank` and `LossGraph`;
convection uses `fs-convection`; the solid/air coupling remains the same FEM path.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/fan-hotspot.json
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/fan-correlated-hotspot.json
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/size-fan-hotspot.json
cargo test -p fs-cli --bin frankensim network_command
```

## Fan boundary drive

Inside `hydraulics`, choose exactly one of the original `boundaries` or a `fan`
object. The complete sample declares the bank's name and source identity, two
terminal node indices, discharge-air temperature, identical-fan count and
series/parallel arrangement, speed ratio and admitted speed domain, minimum
base-curve flow, monotone base-speed pressure/flow points, and explicit scalar
root-search budgets. `pressure_tolerance_rel` is a caller-declared allowance;
it is retained, not propagated into a physical uncertainty interval.

The fan supplies one inlet node from ambient. The outlet node defines zero
pressure. This is not an internal fan edge or a recirculating loop, and a request
cannot simultaneously prescribe pressures and ask the fan to determine them.
Bank arrangement and speed scaling are supplied by the existing fan model.
Disconnected paths, stall-domain violations, or exhausted root-search budgets
refuse. Every leakage or bypass path must still be explicit.

The result's `fan` object records the actual evaluated speed, through-flow,
static pressure, reconstructed fan/network pressure residual and `air_power_w =
Q * delta_p`. Electrical input power is `null`: no efficiency or motor curve was
provided. Neither motor/fan heating nor pressure-loss dissipation is included in
the sensible-heat model. The temperature supplied by the request is the declared
discharge temperature, not an automatically predicted fan temperature rise.

## Flow-derived convection

Each surface requires exactly one of `htc_w_m2_k` or `convection`. A convection
object declares a named card, hydraulic diameter, aggregate free-flow area,
channel length, dynamic viscosity, fluid thermal conductivity and an input
source description. Reynolds uses that surface's owning branch's actual absolute
flow, and Prandtl is computed using the same specific heat as air transport.
Wetted exchange area still comes from the solid mesh, not the free-flow area.

Supported cards are circular and rectangular fully developed laminar CWT,
circular Hausen developing, the retained rectangular developing CWT table slice,
Dittus-Boelter, and Gnielinski. Rectangular cards require `aspect_ratio` in `(0,1]`;
other cards must omit it. Dittus-Boelter requires `thermal_direction` equal to
`heating-fluid` or `cooling-fluid`; other cards must omit it. The solved regional
heat direction is checked against that declaration. Constant-flux, natural and
external-flow cards are not reinterpreted as duct CWT models.

The existing card enforces its validity limits. No automatic regime switch,
extrapolation, or guessed missing geometry is performed. Correlation evaluation
runs again whenever a trial changes flow. Its result includes dimensionless
groups, Nu, h, formula provenance and source-region classification. A declared
engineering bridge remains identified as such. Coefficients and properties stay
frozen within each coupled solve; per-surface lengths do not implement continuous
boundary-layer development across successive surfaces.

Existing effective-h sizing can vary a declared coefficient, but cannot override
a surface whose coefficient is derived from a convection law. Reported inlet
and effective-h gradients remain conditional thermal sensitivities at fixed
hydraulics and properties, not fan-speed or geometry derivatives.

## Fan-speed target search

The optional top-level `fan_speed_design` replaces, rather than combines with,
`design`. It requires `hydraulics.fan` and declares:

```json
{
  "min_speed_ratio": 0.5,
  "max_speed_ratio": 2.0,
  "temperature_limit_k": 302.2,
  "speed_ratio_tolerance": 0.0001,
  "temperature_tolerance_k": 0.00001,
  "max_evaluations": 64
}
```

Every trial solves the fan/graph intersection, recomputes all derived
coefficients, and solves the heterogeneous solid and warmed/mixed air together.
The selected objective can be a mean or the actual discrete hotspot. Bisection
uses evaluated temperatures, not the fixed-flow adjoint as a fictitious fan
sensitivity. `objective.gradient=false` avoids unnecessary adjoints; true retains
the existing conditional thermal gradients on each evaluation.

A feasible declared minimum returns after one evaluation. Otherwise both a
failing lower and passing upper endpoint are required. The returned passing
field must meet the temperature limit, the passing slack tolerance and the
speed-bracket width tolerance. Both endpoints failing means no passing bracket,
not proof that an interior passing design cannot exist. Multiple crossings are
possible for arbitrary models; no global minimum-speed theorem is claimed.

The output includes the selected speed and its matching flow, field,
coefficients, failed lower endpoint, trial history and work count. Domain or
solver failures propagate, rather than masquerading as infeasible designs. The
single invocation wall budget covers all trials; evaluation exhaustion returns
exit 6 with no partial field. There is no motor-speed/RPM conversion without a
supplied reference RPM.

The sample is a two-material solid with one localized watt, a 300 K supply and
an explicit bypass. The curves, channels and materials are illustrative
inputs. Independent NumPy tetrahedral-FEM references put its baseline hotspot
near 302.457593 K at speed ratio 1.0. These are reference calculations, not
retained Rust-run measurements. Rust compilation, formatting and execution of
the added tests have not been performed in the authoring environment.
