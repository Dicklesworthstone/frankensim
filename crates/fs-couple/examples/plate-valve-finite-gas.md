# Heating the gas in a moving plate-valve performance

[plate-valve-finite-gas.performance](plate-valve-finite-gas.performance) applies
signed heat power to the finite gas already driving the original plate-valve
and DG bore. The source-supported NASA-9 energy law updates temperature and
mean pressure. The same accepted gas then updates acoustic compressibility,
transport, wall memories and aperture density, including the acoustic parameter
work paid by that inventory. The existing moving plate and pressure renderer
continue on their original sample clock.

The example starts with nitrogen at 300 K and 95000 Pa in a 50 mm long,
1 mm radius bore. A constant 5 Pa blowing-pressure control drives the authored
illustrative plate. During the 0.1 s performance, heat power is 10 mW for the
first third, zero for the second, and -5 mW for the final third. Carbon dioxide
is explicitly available in the initial field's source catalogue at zero fraction.
The plate coefficients are illustrative, not a measured instrument specimen.
The initial pressure leaves room for heating within the supplied transport
pack's explicit 75000..101325 Pa application domain.

## Render the example

Build the six nitrogen and carbon-dioxide packs using the commands in
[Gas composition and temperature along a reed bore](reed-axial-gas-field.md#build-the-supplied-source-packs-and-render).
Then run from the repository root:

```sh
cargo run -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/plate-valve-finite-gas.performance \
  target/plate-valve-finite-gas.wav --block 113 --decimate \
  --gas-component target/reed-field-packs/nitrogen.pack \
    target/reed-field-packs/nitrogen-nasa9-mcbride2002.pack \
    target/reed-field-packs/nitrogen-transport-lemmon2004.pack \
  --gas-component target/reed-field-packs/carbon-dioxide.pack \
    target/reed-field-packs/carbon-dioxide-nasa9-mcbride2002.pack \
    target/reed-field-packs/carbon-dioxide-transport-laesecke2017-huber2016.pack \
  --gas-water none --gas-phase frozen
```

The WAV contains the existing **internal inlet-pressure diagnostic**, sampled
at 48 kHz after the original 96 kHz physical trajectory. Use a new output name
for each comparison. To compare with zero external heat, copy the performance
and change all three `heat` powers to zero. Keep geometry, gas sources,
blowing pressure, observation and numerical controls the same.

## Author heat on the physical clock

Place the following optional records after the complete DG declaration and
before `observation`:

```text
finite_gas thermochemical 300 95000 mole 1 N2 1
gas_energy 200 900 1e-8 128
gas_feedback -1e-5 1e-5 1e-15 1e-11 80 2 spatial-viscothermal 0.01
gas_heat 3
heat 0 0.01
heat 3200 0
heat 6400 -0.005
```

`finite_gas` uses the existing thermochemical gas grammar: temperature in K,
pressure in Pa, composition basis, species count, and species/fraction pairs.
It must resolve to the **same complete source state** as the uniform initial
`gas_field`, including material identities. Inventory volume comes from the
actual DG geometry; mass comes from that volume and the original gas density.
Neither volume nor mass is independently authored in this block.
The `ambient` and `finite_gas` records list active species only: omit absent
components. A `gas_field` may retain additional source columns at exactly zero;
its resolved active state still matches the inventory without inventing a trace
concentration.

`gas_energy` supplies minimum and maximum temperature, specific-energy residual
tolerance in J/kg, and the maximum caloric iterations. `gas_feedback` supplies
lower and upper work brackets in J, work and energy tolerances in J, maximum
coupling iterations, maximum temperature change per physical sample in K, and
the existing acoustic loss owner. Select `lossless`, or
`uniform-viscothermal`, `spatial-viscothermal`, or `molecular-diffusion` followed
by its relative error bound. A selected owner must be compatible with the
original load and its source domain; these settings do not add a missing
material law.

Heat powers are signed watts, held until the next assignment. The first sample
must be zero, subsequent samples must strictly increase, and every assignment
must lie inside the performance window. Each accepted physical sample supplies
`power * physical_time_step` joules once. Callback block size and output
decimation do not change this clock. Original source validity, nonlinear
residual limits, work brackets, thermal increments and cancellation remain
active throughout playback.

This route evolves one closed, well-mixed inventory in a uniform circular fixed
bore under its original acoustic boundary conditions. It requires the explicit
frozen-gas phase policy. Computed exterior radiation, independent branch gases,
and an upstream tract require their own mean-pressure and mass exchange and are
refused by this attachment. Spatially varying initial fields and tapered bores
also refuse. The output therefore makes no exterior-microphone or condensation
claim. Acoustic/numerical losses retain their original accounts and are not
automatically deposited as bulk heat.

Library users can inspect `AperturePerformance::finite_gas()` for the accepted
mass, volume and material state, `applied_heat_controls()` for the committed
assignments, and `last_finite_gas_step()` for the original gas/acoustic/mechanical
work report. A refused sample leaves those physical owners at the previous
accepted state and follows the existing renderer's poisoned-callback behavior.
