# Mixed material gas heat laws on one reacting assembly

This patch connects gas-side convection and conduction to the existing prepared
material exchange network. The same assembly retains electrical contact losses,
thermal interface resistance, enclosure radiation, chemical progress, finite
cell enthalpy, gas internal energy and piston momentum. No second integrator,
material catalog, equation of state or heat deposit is introduced.

## Physical entry points

`PreparedMaterialExchangeNetwork::step_with_gas_law_candidate` accepts explicitly
mapped `GasSurface` patches and a typed constitutive callback. The callback gets
the actual gas energy midpoint, proposed thermal history and complete material
sample. It returns total conductance in W/K and its original source payload.
The original owner computes signed heat, sums it once per gas, debits each cell,
and returns all endpoints together. The old `step_candidate` delegates to this
same path using the original series conduction implementation.

`reactive_with_gas_law_candidate` uses the original chemical staging, including
optional bounded chemistry subcycling, then the same exchange solve. Reactions
are not repeated by nonlinear heat/motion probes. Formation energy is not added
again as heat. Original whole-interval chemical storage/entropy and port history
remain intact. The old `reactive_step_candidate` uses the same shared staging.

The concrete data-only adapter is:

```
fs_airflow::graph::thermal::medium::convection::finite_gas::network::exchange::step_candidate
```

Each `ExchangeContact` carries the existing chamber/cell/contact IDs and either
`GasHeatLaw::Conduction(GasConductionGeometry)` or
`GasHeatLaw::ForcedConvection(ForcedGasSurface)`. Empty reaction bindings mean
nonreacting gas. Each patch retains its actual typed source in `GasHeatSample`.
The convection path consumes the reacted trial gas and the actual cell phase;
the conduction path reuses the complete pinned material bundle, not a second
conductivity-only query. The prepared electrical/thermal/radiation sources stay
with their existing owners. Different laws may act on disjoint patches without
changing the solver code.

## Surface, energy and state semantics

All patch areas share one aggregate exposure check across every chamber.
Thermal links and reservoir coverage stay with the original thermal owner.
Nonzero ambient convection on a bound cell refuses, preventing duplicate heat.
Radiation can share an exposed gas patch; two gas heat paths cannot silently
claim overlapping geometry. Area arithmetic does not prove spatial disjointness.
The law's geometry must match the declared area. A zero-area callback is still
invoked, but nonzero accepted conductance on a zero-area patch refuses.

Each signed contact transfer is kept even when summed chamber heat cancels.
Electrical losses are deposited once. Closed-enclosure radiation roundoff is
not external power. Original gas, cell, circuit and mechanical balances remain
in force, with one common interval and one chamber account. Source errors retain
their original type. Refusal, work exhaustion and cancellation return no partial
candidate. Consume every part together; do not apply returned heat or chemistry
a second time.

This extension does not add these laws to `MaterialExchangeState`'s fixed-source
live transport entry points. The new candidate can be advanced repeatedly by
publishing its complete returned states. No gas mass transport is added by this
adapter. Existing transport remains a separate explicit stage/entry point.

## Validation commands and status

```
cargo test -p fs-couple --test material_exchange_network gas_law_tests
cargo test --manifest-path crates/fs-airflow/examples/material_convective_network/Cargo.toml exchange_tests
```

Four shared-owner tests cover exact original conduction/radiation reduction,
typed source retention through one chemical stage, shared area/zero-area
admission, original source errors and cancellation after chemistry. Three
adapter tests cover original-convection reduction, reacting mixed-law history,
data-only kinetic changes, shared coverage, ordering and failed-state atomicity.
They reuse existing manufactured material and circuit fixtures.

These seven Rust tests were authored but NOT executed. Cargo and rustc are
absent in the authoring environment; native compilation is unverified. Baseline
Git hashes, patch application, reconstructed-file bytes, whitespace and lexical
delimiters were checked. Those checks do not establish Rust type correctness,
numerical convergence, native performance or measured-material fidelity.

## Physical limits

The original first-order chemical splitting and thermal-endpoint/gas-midpoint
coupling remain. Cells and gases are uniform, piston walls are independent,
thermal geometry/C/L are fixed, and the radiation model is an opaque gray
nonparticipating-gas enclosure. Forced speed and correlation are explicit;
there is no fan-power, shear-heating, resolved boundary-layer or moving-liquid
model. Melting thermal state does not silently change a rigid wall's mechanics.
Temporal refinement and each source/correlation's validity remain necessary.
