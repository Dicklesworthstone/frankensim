# Electrical interfaces on the shared material/enthalpy timeline

`fs_lbm::conjugate::enthalpy::transport::interfaces::
advance_joule_enthalpy_with_material_interfaces_observed` evaluates electrical
interfaces at every actual thermal substep, alongside the existing bulk material
and thermal-link evaluations. It does not introduce another current or heat solver.

One bulk material evaluation per mapped cell supplies the original cell, its
specific enthalpy, temperature, resistivity, conductivity, and concrete material
payload. Both ordered surface samples reach each electrical interface callback.
Its result supplies area-specific resistance and heat partition, retaining the
original constitutive object. Electrical topology remains the explicit ordered
`input.contacts` list; unlisted faces remain perfect electrical contacts.

The existing current solver combines bulk half-cell and interface resistance,
then deposits its combined bulk/contact heat exactly once. Thermal transport uses
its original material laws, reference masses, contact coverage, and actual CFL
step. An additional observer, or a nested `CoupledObservers`, can supply other
heat sources and physical histories on that same time interval. It cannot also
own an already selected thermal link's total conductance.

## Source-backed interfaces, not copied constants

`MaterialContact::into_enthalpy_response(expected_cells)` moves the existing
source-supported surface assignment into the new path. It preserves the original
`ResolvedElectricalInterface`, material/process/finish checks, query coordinates,
and receipts. A changed or reversed spatial assignment refuses. Interface
sources still choose their own explicitly defined temperature, pressure, phase,
finish and history coordinates; the adapter never infers an average temperature
or an electrical contact law from the bulk material name.

Applied electrical interfaces and bulk/thermal epochs are retained separately
from the final endpoint evaluations. Endpoint support produces no additional heat.
A failed endpoint, exhausted work budget, observer refusal or cancellation returns
no accepted thermal/observer candidate. Callbacks must not mutate external physical
history; immutable source objects need not implement `Clone`.

## Run the existing phase-change example with a contact data file

```sh
cargo run -p fs-lbm --features enthalpy-coupling --example joule_phase_change -- \
  crates/fs-lbm/examples/joule_phase_change/transport_phase_dependent.txt \
  10 60 0.05 \
  crates/fs-lbm/examples/joule_phase_change/contact_low.txt
```

Replace only `contact_low.txt` with `contact_high.txt` to change the electrical
interface data. Geometry, bulk caloric/electrical/thermal data, voltage pulse,
clock, and solver remain unchanged. The source turns off halfway through the
requested duration; the stored enthalpy and phase history remain afterward.
The optional output column is endpoint contact resistance in ohm square metres.
Existing arguments, original columns, and no-contact behavior remain available.

The first non-comment contact row declares either `liquid_fraction_b` or
`temperature_b_K`, followed by the fraction of contact heat delivered to the first
surface. The remaining rows are ordered coordinate/resistance pairs. The loader
constructs an immutable interface card; the canonical material resolver evaluates
its bounded linear-inside law. It retains both actual surface temperatures and,
for a phase-axis source, the second surface's liquid fraction. Resistance can then
change across an isothermal melting plateau. A source without support refuses;
it is not replaced by zero resistance, insulation, or an extrapolated value.

These data files are synthetic equation fixtures, **not measured properties of
any named metal, alloy, or surface finish**. The model remains a stationary,
fixed-volume, explicit Ohmic/enthalpy approximation. It does not supply evolving
contact geometry, melt motion, pressure work, magnetic induction, phase kinetics,
or an implicit stability guarantee. Thermal CFL does not bound a nonlinear
source's stiffness: refine `max_dt_s`. Electrical and thermal interface resistance
remain distinct laws; this example adds only electrical interface resistance.
