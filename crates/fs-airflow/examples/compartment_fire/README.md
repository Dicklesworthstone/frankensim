# Source-owned compartment-fire example

This executable composes the original gas, stratified-opening, Heskestad plume
and thermal-enthalpy modules into an estimated two-zone room simulation. It
initializes a real lower layer and an exactly absent upper layer; entrainment
creates the upper layer without artificial gas seeds. Each gas layer keeps
species mass, NASA calorics and the moving-interface energy account.

## Run

From the repository root:

```bash
cargo run -p fs-airflow --features compartment-fire --example compartment_fire -- \
  data/reference/compartment_fire/steckler_1982_cases.csv 10 10 10 0.5 1800 4
```

Arguments are the retained CSV, case number, inner convection coefficient
[W/(m² K)], outer convection coefficient, maximum timestep [s], end time [s],
and cells through each wall. The last two arguments default to 1800 s and four
cells. Both convection coefficients must be supplied explicitly. Tests 10, 11,
12, 13, 14 and 16 are available. CSV output contains transient temperatures,
interface height and pressure; stderr contains balances and the comparison.

The two h=10 coefficients above are chosen engineering inputs. They are not
measured or fitted parameters. The executable retains the CFAST-authored
62.9 kW heat release, 20% radiant fraction, 0.09 m² burner area and 0.02 m
base height. Its nominal doorway discharge coefficient is 0.7. The retained
raw experiment and authored CFAST inputs remain separately identified in
[the source notes](../../../../data/reference/compartment_fire/steckler_1982_sources.md).

## Original modules and explicit approximations

- Gas state and pressure/interface work: fs-gas two-zone inventory owner, with
  original NASA9 N2/O2 calorics for 200..1000 K. The explicitly chosen dry model
  mixture is 0.79 N2 / 0.21 O2 by mole. Caloric extrapolation is refused.
- Plume: fs-convection Heskestad correlation using the actual lower-layer
  density and mass-specific heat capacity. One donor species/enthalpy vector
  is transferred with identical opposite signs.
- Door: fs-airflow stratified hydrostatic opening, including simultaneous
  directional slabs, solved through the bounded pressure-coupled network.
- Walls: three fixed-mass profiles (vertical walls, ceiling, floor), using the
  original fs-thermal enthalpy kernel, backward-Euler residuals and series
  conductances, including inner and outer half-cell thermal resistance.

Wall k=1 W/(m K), density=200 kg/m³, Cp=1000 J/(kg K), and thickness=0.013 m
are the authored CFAST approximation. They are not measurements of the
experimental multi-layer construction. Interface motion changes contacted
areas without rescaling wall material mass. The prescribed radiant fraction
is allocated uniformly by enclosure area; no gas emissivity is inferred.

The network recomputes duration-dependent wall candidates after every failed
whole-step attempt. Only accepted gas and wall states are published together.
The pressure Jacobian perturbation shrinks as the pressure residual converges,
which avoids unnecessary retries near vent neutral-plane transitions. Gas
source temperatures remain frozen per attempted interval, so this coupling
is first order; it is not a fully temperature-implicit or second-order method.

This executable uses prescribed external heat. It does not consume oxygen or
fuel, generate combustion products/soot, include participating-gas radiation,
apply a weak-plume penetration correction, or model pyrolysis. The separate
InventoryLimitedReaction::apply_split API can be composed through
ZoneNetwork::step_with_stages when an explicit conserved reaction is supplied.
It uses inventories after transport and preserves the physical clock; add its
chemical source exactly once and never add its diagnostic reaction heat again.

## Actual public-data comparison

The following are means over 1500..1800 s at maximum dt=0.5 s, four wall cells
and the two explicit h=10 coefficients above. Each table pair is
**prediction / retained experimental reduction**. No parameters were adjusted
to improve agreement. Reduced-layer experimental standard uncertainty is
unknown; these results do not establish a validation pass.

| Steckler test | Upper temperature (°C) | Lower temperature (°C) | Interface height (m) | Upper temperature-rise error |
| --- | ---: | ---: | ---: | ---: |
| 10 | 186.702 / 180.100 | 49.673 / 65.500 | 0.753202 / 0.480000 | 4.28% |
| 11 | 171.864 / 153.700 | 47.042 / 56.400 | 0.983717 / 0.580000 | 14.45% |
| 12 | 151.400 / 125.900 | 37.914 / 43.800 | 1.109562 / 0.670000 | 24.54% |
| 13 | 144.297 / 118.400 | 37.170 / 40.300 | 1.185806 / 0.710000 | 27.15% |
| 14 | 146.478 / 115.300 | 42.122 / 44.200 | 1.237369 / 0.780000 | 36.13% |
| 16 | 138.159 / 107.400 | 37.890 / 37.800 | 1.278803 / 0.820000 | 37.79% |

The geometric upper-temperature-rise ratio is 1.234814765,
with combined log scatter 0.106556125. Separated model uncertainty
and adjusted FDS bias remain absent because experimental standard uncertainty
is unknown. The interfaces are systematically too high and the hot-layer
depths too small. Lower-layer temperature error is also retained; matching one
upper temperature is not sufficient physical validation.

Across these six runs, maximum per-step mass residual is
1.540434e-15 kg and the largest accumulated gas-plus-wall energy
imbalance is 2.152586e-3 J. Every room receives 113.22 MJ of
prescribed fire heat over the run; gas outflow and outer wall convection remain
explicit external exchanges in the balance.

For Test 10, reducing maximum dt from 0.5 to 0.25 s changes the mean upper
layer temperature by 2.992500e-5 K and the interface by
1.190000e-7 m. This checks the stated late-time quantities,
not a global second-order convergence claim. The adaptive pressure change
reduces Test 14's whole-step halvings from 10,858 to 2; its mean upper-layer
prediction changes by less than 3e-8 K.

[Machine-readable runs](../../../../data/reference/compartment_fire/steckler_1982_two_zone_example.json)
retain all six cases and the timestep refinement. Dalmarnock is not used in
this example or its calibration. The correlation-only rung remains runnable
as fs-convection's compartment_fire_correlations example.
