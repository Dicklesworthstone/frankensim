# fs-nuclear contract

## Prescribed fission, finite inventory and delayed energy

`fs-nuclear` is an L3 source owner. The default kernel uses only `std` and
shares immutable nuclear model data with `Arc`. It accepts a prescribed total
fission rate in expected nuclei per second, finite fuel inventory, explicit
initial delayed-energy reservoirs, SI energy coefficients and caller budgets.
An identified dataset and its domain of applicability remain the caller's
responsibility. There is no built-in isotope or silently assumed irradiation
history. The supplied data identity is a label, not scientific validation.

`EnergyChannels` separates energy deposited in the source, outward radiation
available to another transport owner, and nonrecoverable export. Each model
has a prompt partition and an ordered list of nonnegative exponential decay
groups. A group's coefficients are integrated joules per fission, with a
positive decay constant in inverse seconds. They describe effective energy
reservoirs, not delayed-neutron precursor populations. For impulse-power fit
coefficients `alpha * exp(-lambda * t)`, integrated energy is `alpha/lambda`.

For a constant prescribed rate `F`, each group obeys

```text
dE_i/dt = q_i F - lambda_i E_i.
```

The accepted event count is limited by actual remaining fuel. If it is exhausted
within an interval, fission stops at `N/F` and all groups continue decaying to
the interval endpoint. For active duration `a`, cooling duration `c`, and
`phi(x) = (1-exp(-x))/x`, the endpoint is

```text
E_i,new = E_i,old exp(-lambda_i (a+c))
        + q_i (F a) phi(lambda_i a) exp(-lambda_i c).
```

Integrated release is computed with `expm1` and a short-interval series,
independently of subtracting large stored-energy totals. It therefore retains
small positive decay heat. No explicit decay substeps or endpoint-power
rectangle are needed. Group storage remains nonnegative. Turning the rate to
zero does not erase the accepted reservoirs or stop their heat release.

## Accounting, numerics and work

Every accepted interval reports fissions, active duration, generated energy,
three released channels, surviving delayed-storage change, a fuel residual,
an energy residual and accumulated absolute numerical energy defects. The
energy identity is

```text
generated_energy = delta_delayed_storage
                 + deposited + outward_radiation + nonrecoverable_export.
```

Generated energy sums the prompt and delayed yields of the same event count
that debits fuel. Local group defects cannot cancel one another to satisfy the
energy budget. Compensated accumulators retain small repeated fuel debits and
lifetime release credits. An unresolved reservoir debit consumes the explicit
energy allowance; it cannot silently create heat. Nonfinite inputs, negative
rates/inventories/channels, unrepresentable physical quantities, excess work
and excess accounting defects refuse the candidate. The tolerances bound
numerical account closure, not nuclear-data uncertainty.

Each trial uses O(groups) work and storage, polls the caller checkpoint before
work, at every group, and before publication. `trial` never modifies its input;
`advance` publishes only a fully accepted candidate. The group count is bounded
by the caller. Cloning retains the full accepted inventory and accumulated
accounts. Replay with identical coefficients and interval partition is bitwise
stable on the same ISA/toolchain/libm profile. Different interval partitions
agree to floating-point tolerance; cross-ISA equality is not claimed.

## Physical scope

This is a reduced prescribed-source model. It does not infer neutron flux,
cross sections, criticality, reaction rates, temperature reactivity feedback,
isotopic product yields, activation, poisoning or a nuclide depletion network.
It does not solve radiation transport or assign dose. The radiation channel is
an energy export, not proof of where that radiation subsequently deposits.
Group fits require independent calibration over the intended cooling times
and irradiation histories. No regulatory or reactor-design qualification is
claimed.

Fission cannot be represented by an ordinary chemical reaction preserving
each element. These source accounts therefore do not mutate `fs-thermochem`
species or invent a chemical stoichiometry. Likewise, fuel nuclei consumed by
fission become products; they are not lost bulk material mass. Thermal material
mass, composition evolution and mass-energy changes require their own physical
owners. A fixed material caloric law is a declared reduced-model assumption.

## Evidence and references

The 12 owning source tests cover G0 account closure/admission, G1 independent
constant-source and shutdown solutions, separated channels, multiple decay
timescales, short and stiff intervals, finite fuel exhaustion, G3 interval
partitioning, G4 cancellation/refusal rollback, and G5 identical checkpoint
replay. These passed using the pinned Rust 1.99.0-nightly toolchain through
`rustc --edition=2024 --test crates/fs-nuclear/src/lib.rs`. Coefficients are
explicitly synthetic: equation verification is not isotope validation or a
performance benchmark. Workspace Cargo/DSR verification remains separate.

The exponential source-history formulation and distinction between prompt and
delayed fission energy follow the public [Argonne SAS4A/SASSYS-1 decay-heat
documentation](https://sas-doc.nse.anl.gov/5.7/Part02/Ch04/decay-heat.html),
especially sections 4.4.1.1 and 4.4.1.3. The implementation does not copy or
supply the ANS coefficient tables. [INL BISON's power and burnup
documentation](https://mooseframework.inl.gov/bison/theory/power_burnup.html)
provides context for prescribed fission sources and why energy production,
fuel history and material evolution are distinct quantities.
