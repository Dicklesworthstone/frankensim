# One inventory for pressure flow and molecular mixing

`fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::step_candidate`
combines source-backed nozzle connections and Maxwell-Stefan diffusion faces on
one accepted gas snapshot. Both samplers are the existing implementations; neither
reconstructs an endpoint until all rates have been evaluated. No new gas law,
species registry, diffusivity formula, caloric inversion or runtime dependency is
introduced.

The two lists declare **distinct physical paths**. A pressure restriction and a
near-isothermal/isobaric molecular exchange path may connect the same volumes as
separately authored mechanisms. This does not identify a nozzle's bulk velocity
with the diffusion law's molar-average frame or turn them into a resolved
advection-diffusion law for the same aperture. Path geometry and applicability
remain explicit model inputs, not inferred from node IDs or gas names.

## Common conservative update

The original nozzle supplies integrated upstream bulk mass. The existing
composition owner divides it into species masses; the selective transport owner
carries each original species' partial enthalpy. Summing these is the same bulk
enthalpy transport, including flow work and formation references. Molecular
counterflow uses its own signed species rates and the same energy owner.

All species transfers are accumulated before **one reconstruction per affected
volume**. Gross withdrawal is bounded against the ORIGINAL individual species
inventory across both mechanisms. Incoming material or opposite transfers cannot
conceal excessive withdrawal. Individually safe mechanisms can therefore be
refused together. Newly received species cannot leave another path during the
same explicit stage. Reports are ordered by IDs; physical update arithmetic is
not selected by those IDs.

Original source conflicts, local/global species/mass/energy checks, trace-species
admission and cancellation remain. The final combined state must also pass the
existing initially-active nozzle pressure-order check and extensive entropy gate.
These checks are not a time-error estimate or proof of unconditional stability.
There is no clipping, corrective energy deposition or second reaction-heat credit.

## Live electrical, thermal, chemical and moving-wall assembly

The existing `MaterialExchangeState`, under
`fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state`,
now exposes `transport_step_candidate` and `advance_with_transport`.

Supply the prepared original material network, target time, explicit chamber
drives, optional declared reactions and gas-wall patches, both transport lists,
`GasTransportSolve`, the original `NetworkSolve`, and its thermal initial guess.
Per-chamber energy inversion controls come from the existing drive. The shared
thermal clock supplies one duration to every stage:

1. Sample and combine both transport mechanisms on the original inventories.
2. React the reconstructed mixtures once per bound chamber.
3. Run the existing reciprocal wall-heat, circuit/enthalpy and momentum solve.

Nonlinear trials do not repeat steps 1 or 2. The original all-owner accounting
starts BEFORE transport and includes gas internal energy, wall kinetic energy,
electrical storage, finite thermal storage and external work/heat. Partial
enthalpy and formation energy are not supplied again as heat. Publication moves
all owners and accounting histories together, or none. Gross exchange marks the
closed two-port chamber audit unavailable even at zero net mass; the complete
assembly balance remains available.

## Focused checks

```sh
cargo test -p fs-material --lib \
  thermochemical_gas::energy::parcel::transfer::species::network::transport::
cargo test -p fs-couple transport_tests:: -- --test-threads=1
```

Six library regressions cover original-route parity, independent snapshot
sampling, jointly excessive withdrawal, order independence, atomic refusal and
continued use by the original heat/work consumers. Four integration regressions
cover the actual complete assembly, data-only diffusivity changes propagating to
downstream flow/heat/charge/motion/sound, one chemical stage, and late-failure
retry/cancellation. Native compilation and tests were not executed in the
editing environment because Cargo and rustc are absent; authored tests are not
reported as passes.

The method is first-order fixed-volume transport followed by the existing
reaction and coupled heat/motion update. Diffusion keeps its admitted temperature
and pressure differences; heating or motion may invalidate a face on the next
call, which must refuse rather than relax the model. It is not implicit CFD,
pressure/thermal diffusion, stiff chemistry, spatial flame or jet momentum,
shared-crank dynamics, or automatic live acoustic-renderer reconfiguration.
