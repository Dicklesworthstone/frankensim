
### Native finite-gas heat programs in plate-valve performances

Native plate-valve files accept `finite_gas`, `gas_energy`, `gas_feedback`,
`gas_heat` and signed `heat SAMPLE POWER_W` records after the complete DG input.
They retain original source packs and derive inventory volume from the native
geometry. WAV metadata records the initial inventory, original physical clock,
power assignments and solve controls before type erasure. This is closed,
fixed-volume gas heating with the existing internal pressure observer. It does
not supply mean mass/enthalpy transport, phase change, wall heating, an
independent upstream/branch inventory, or a computed exterior mean-pressure
jump. Such inputs refuse instead of changing their physical boundary.

