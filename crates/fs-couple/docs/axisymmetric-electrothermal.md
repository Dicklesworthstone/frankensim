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
current and a circuit storage model are outside this adapter. Contact-bearing
electrical owners refuse until an explicit surface heat partition and contact
source refresh are supplied; their surface losses must not be dropped or
counted twice.

Focused regression definitions cover sector scaling, unequal swept triangle
volumes, local current heating with material Fourier transfer, opposite
voltage/current feedback, regulated-power continuation, latent-phase source
refresh, and late-refusal retry. These use manufactured material laws and do
not qualify measured conductor data or establish continuum accuracy. Refine
the original spatial and first-order temporal models for the intended case.
