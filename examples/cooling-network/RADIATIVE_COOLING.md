# Ambient radiation in coupled cooling

```bash
frankensim --json cooling-network \
  examples/cooling-network/radiative-contact-hotspot.json

frankensim --json cooling-network-uq \
  examples/cooling-network/radiative-contact-hotspot.json \
  examples/cooling-network/uq-radiative-contact.json \
  --checkpoint radiative-study.uqcp
```

The checkpoint destination must be new. The example has a 20 W localized P1
source, two temperature-dependent materials, a matching thermal contact and a
fan-driven split/merge network. Two existing cooling surfaces also exchange
radiation with their declared surroundings. Hardware inputs and uncertainty
laws are illustrative, not measured or calibrated data.

## Declare exposed patches

Add this at the request root, with `objective.gradient=false`:

```json
"radiation": {
  "max_iterations": 128,
  "temperature_tolerance_k": 1e-9,
  "relaxation": 0.5,
  "surfaces": [
    {
      "surface": "first-face",
      "emissivity": 0.85,
      "ambient_temperature_k": 280,
      "source": "Declared exposed gray patch facing a large isothermal reservoir"
    }
  ]
}
```

Each patch identifies exactly one existing `solid.surfaces` entry. It uses that
entry's faces and area; it neither replaces convection nor selects contact faces.
Duplicate patch ownership and unknown names refuse. Emissivity is constant and
strictly in (0,1]; absolute surroundings temperature is positive. Omit a patch
(or the whole radiation option) to omit radiation, rather than using a zero
emissivity to introduce a different boundary mode. The selected emissivity is
admitted through the existing material-card machinery with explicit
caller-declared provenance and unstated uncertainty; an inline card does not
turn a declaration into a measurement.

The first implementation is steady and primal-only. Radiation combined with
transients, requested gradients, mesh studies or nested nominal design searches
refuses at admission. It does not silently drop radiation or return an adjoint
with frozen radiation. Ordinary requests without radiation retain their existing
solver paths. UQ's external finite candidate family can still call the steady
radiating producer, because each candidate is an ordinary complete request.

## Precise physical and spatial model

For the area-mean P1 temperature `Tbar` of a selected patch,

```
Q_rad = epsilon * sigma * A * (Tbar^4 - T_surroundings^4)
h_rad(Tbar) = epsilon * sigma * (Tbar + T_surroundings)
                           * (Tbar^2 + T_surroundings^2)
```

The constant comes from `fs-conduction`'s existing radiation module. The factored
expression avoids subtracting nearly equal fourth powers and has a finite
positive limit at equal temperatures. This is NOT a fixed small-departure
coefficient evaluated once at an arbitrary temperature.

The spatial closure is deliberately explicit: at a fixed patch mean the Robin
flux is `h_rad(Tbar) * (T_h(x) - T_surroundings)`. Integrating this expression
reproduces the stated mean-temperature fourth-power law at convergence. On a
nonuniform patch this is neither a uniform prescribed flux nor the pointwise
integral of `epsilon*sigma*(T_h(x)^4-T_surroundings^4)`. Changing the patch
partition can therefore change the model. The exact isothermal-patch limiting
case is tested against an independent one-dimensional slab solution.

Each patch sees a large isothermal black reservoir with view factor one. There
is no self-shadowing, exchange between modeled patches, gray-enclosure reflection,
participating fluid absorption, solar spectrum or transient radiation model.
Do not assign this reservoir model to mutually facing internal surfaces and
interpret it as a view-factor enclosure calculation.

## Coupling and separate heat ownership

At fixed air references, the solver combines the two Robin laws algebraically:

```
h_total = h_air + h_rad
T_combined = (h_air*T_air_reference + h_rad*T_surroundings) / h_total
```

It solves the actual existing nonlinear-material/contact FEM problem, recomputes
patch means and nonlinear radiation, and iterates with the declared relaxation.
Both the raw patch-temperature update and the nonlinear/applied heat mismatch
must meet their tolerances. Convective coefficients remain the existing declared
or correlation-derived values; radiative h is never sent to a convection card.

The outer mixed-air driver receives ONLY `h_air*A*(Tbar-T_air_reference)`.
Radiative power is not an additional heat source for the air. Both positive
loss to colder surroundings and negative loss (heating from hotter surroundings)
are retained, without clipping. The complete accepted result must pass:

* The existing air/solid interface gates using convective heat alone.
* The independently accumulated FEM Robin split and whole-solid source balance.
* `source_w - convective_out_w - radiative_out_w` within the original `heat_w`,
  with radiative heat recomputed from the accepted surface temperatures.

In a radiating result, `robin_out_w` includes the assembled convection AND
radiation boundary heat. `walls[].outward_heat_w` still means heat exchanged
with the air. The `radiation` object reports these mechanisms separately,
including each patch's temperatures, emissivity, secant coefficient, applied
and recomputed watts, plus the final nonlinear mismatch. A colder external
reservoir can draw more heat than the source supplies, with the air then heating
the solid; negative convective totals are valid in that case.

`radiation.max_iterations` limits inner solid solves per air-reference evaluation
and is capped at 1000. `relaxation` is in (0,1]. Exhaustion returns budget exit 6
and no partial cooling result. The original wall deadline covers all inner and
outer solves. `radiation.solid_solves` counts every completed FEM solve across
these nested iterations, while the existing `coupling_iterations` counts outer
air-interface iterations. Neither counts every inner linear iteration.

## Radiation uncertainty

The new UQ target forms are:

```json
{"kind":"radiation-emissivity","surface":"first-face"}
{"kind":"radiation-ambient-temperature","surface":"last-face"}
```

Their units are 1 and K. Every draw changes the named declared input before the
actual cooling child runs. Neither target changes convective h, inlet temperature,
source footprint, material assignment or contact resistance. Uniform support
must remain in the physical domain. A Gaussian draw outside that domain is a
terminal model refusal, never clipped, skipped or redrawn.

Targets and distributions enter the existing plan/checkpoint identity. Completed
observations use the existing atomic checkpoint, exact-ordinal retry and optional
compliance machinery. The sample example is a small empirical study, not a claim
that 16 solves can resolve a stringent probability target. The probability law
and unknown model discrepancy remain the caller's responsibility.

## Checks and limitations

```bash
cargo test -p fs-cli --test cooling_radiation
cargo test -p fs-cli --test cooling_uq_radiation
cargo test -p fs-cli --bin frankensim uq_command::model::radiation
```

Eleven new Rust tests cover independent hot/cold-surroundings slab equilibria,
nonlinear materials/contact, air-versus-radiation energy ownership, a
vanishing-emissivity control, patch-order/ordinary replay, budget/admission
refusals, actual uncertain solves, exact chunked checkpoint/result replay, and
terminal invalid samples. These tests were not executed in the authoring
environment, which has no Rust toolchain.

Independent NumPy finite-element calculations solve the full reduced equations
with a direct Newton method, eliminating the air references analytically. Twelve
constant/nonlinear/contact cases agreed with a separate nested secant and
least-squares mathematical mirror to within 3.80e-10 K. Three constant-k,
zero-source slabs agreed with their independent two-face algebraic solutions
to within 4.72e-12 K. These are numerical reference comparisons, NOT executions
of the Rust implementation or experimental validation.

For the declared nonlinear contact example, the independent reference peak is
329.419419 K; 3.056382 W leaves radiatively and 16.943618 W enters the air.
Omitting radiation gives 332.435497 K on that same mesh/source. These values
are scoped to the stated mean-patch model, not a pointwise radiation solution,
mesh certificate or physical safety bound.
