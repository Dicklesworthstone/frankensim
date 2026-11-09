# Axisymmetric electrothermal feedback

`ConductionThermalNetwork::step_axisymmetric_candidate` advances the existing
finite thermal history and an `AxisymmetricConduction` field to one common
endpoint. It supports the existing `ConductionDrive` voltage, total-current,
regulated-power and distributed-current-density declarations.

Each original electrical triangle maps through `ConductionRegion` to one
thermal receiver. Several triangles can share a receiver, but their complete
swept volume must match that receiver's volume and mass/density. The original
mesh supplies the angular span and radial measure. No planar thickness,
uniform-power redistribution or additional revolution factor is introduced.
Electrical observations use local radial/axial components; watts and joules
refer to the authored sector.

At every thermal trial the adapter refreshes the original complete conductor
sources at the receiver's temperature and explicitly declared liquid fraction.
Source selections, material frames and fixed external coordinates remain
unchanged. The native field determines current redistribution, local E.J
heating and electrical source work. Its cell heat enters the original
backward-Euler enthalpy residual once, together with the existing material
Fourier links and non-electrical thermal sources. Those links retain their
authored transfer areas and half-path distances; this adapter does not infer
radial thermal resistances from the electrical mesh.

The returned candidate contains the accepted thermal timeline, endpoint
electrical material operator and shared energy report. Consume its
`into_parts()` result together for the next interval. Source refusal,
electrical compliance failure, exhausted solve budgets, invalid volume
correspondence and cancellation return neither updated owner.

The original constant-density chart requirement applies through phase change.
An isothermal latent-heat interval still refreshes conductivity through its
declared liquid-fraction coordinate. Geometry motion, expansion, azimuthal
current and a circuit storage model are outside this adapter.

## Finite electrical contacts and two-sided heat

`step_axisymmetric_contact_candidate` additionally accepts the existing
`ConductionContactThermal` declarations, one per original electrical contact
in its canonical order. Each names the original interface card, ordered bulk
cards, thermal/phase coordinates, fixed pressure/history context and authored
fraction of contact heat entering surface A. The remainder enters B. Receiving
thermal cells come from the original adjacent electrical volume routes; neither
a new receiver nor an equal split is inferred.

The same residual refreshes both bulk and contact sources at every trial.
Independent surface temperature and liquid-fraction coordinates follow their
actual caloric receivers, including unequal fractions at one melting
temperature. Source selections, fixed external state and declared phase axes
remain pinned. The native field integrates voltage jumps, current and loss
over the actual swept surface with its exact radial basis moments.

Bulk `cell_power_w` and surface `contact_response` losses are deposited once
through their respective original routes. Their sum equals terminal source
work and the finite thermal energy account. Contact heat fractions are fixed
physical inputs; interface storage, motion, arcing and a thermal resistance
inferred from electrical resistance are outside this model. Authored Fourier
links continue to supply heat transfer independently.

Use the explicit contact entry when surfaces dissipate power. The bulk-only
entry continues to refuse a contact-bearing electrical owner. Both entries
share the existing enthalpy solver and return an atomic electrical/thermal
candidate.

## Multiple terminals and floating electrodes

`step_axisymmetric_terminal_candidate` accepts an
`AxisymmetricTerminalDrive`: complete original electrode vertex sets, native
conductance-response limits and held voltage, total-current or mixed controls.
The original terminal order remains authoritative. Mixed current values follow
the explicitly ordered current-terminal list; voltage values follow its
complement. A zero-total-current terminal is equipotential and can transfer
current between its faces. It is not an insulating boundary.

At every proposed thermal endpoint, the native field rebuilds its N-1 basis
responses from the refreshed bulk and contact material sources. It then recovers
the simultaneous field and local bulk/contact heat, including cross terms
between drives. No conductance matrix from an earlier material state is reused.
The same original enthalpy residual, contact partitions and Fourier links apply.

Terminal source work pairs authored total currents with solved physical voltage
differences; voltage-source reactions supply the remaining work. The pairing
uses retained component-relative potentials and reference differences before
rounding, with one contribution per terminal. Recovered fields report terminal
zero as their voltage gauge.

This entry charges each native basis solve to `maximum_circuit_trials` and
`report().circuit_trials`: N-1 per thermal evaluation, including finite-difference
stencils and the fresh final endpoint. The remaining budget is checked before
those solves start. Native response bounds separately constrain terminal count,
retained storage, dense work and probe/source support. This adds held ideal
terminal control, without circuit storage, finite electrode impedance or a new
contact law.

## Field-dependent bulk and interface transport

`step_axisymmetric_response_candidate` and
`step_axisymmetric_response_contact_candidate` accept the native
`AxisymmetricConductionResponse`. Original scalar conductivity/resistivity
sources can depend on temperature, a declared liquid fraction and actual
electric-field/current-density magnitudes. The same enthalpy residual refreshes
thermal coordinates, then the existing electrical response owner closes its
own coordinates at the actual voltage/current/power drive. Accepted nonlinear
fields are never scaled as if their coefficients were independent of drive.

Explicit nonlinear contact sources likewise close against original surface
voltage jumps/current densities. Their native uniform-trace requirement remains
binding; nonuniform electrically varying interfaces refuse. Thermal/phase
coordinates come from both original caloric receivers, with fixed context and
authored two-sided heat fractions retained. Bulk and surface losses enter the
same energy account once; constitutive defects remain numerical diagnostics.

Every actual native field solve consumes the original thermal trial budget,
including discarded constitutive iterates and final fresh endpoint admission.
Source-derived caloric response/applicability and Fourier conductivity must
remain independent of electrical coordinates. This memoryless scalar model
adds no reversible field storage or thermal-expansion law and promises neither
nonlinear uniqueness nor convergence. Source, budget, compliance or
cancellation failures leave both accepted input histories reusable.

Focused regression definitions cover sector scaling, unequal swept triangle
volumes, local current heating with material Fourier transfer, opposite
voltage/current feedback, regulated-power continuation, latent-phase source
refresh, swept contact work with asymmetric heat partitions, source-bound
surface continuation, latent contact response and late-refusal retry. These
use manufactured material laws and do not qualify measured conductor data or
establish continuum accuracy. Refine
the original spatial and first-order temporal models for the intended case.
