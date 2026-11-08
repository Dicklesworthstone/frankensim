
## Original density through the shared thermal interface

`PcmComposite`, `CrystallizationCaloricBranch`, `FreezingWaterThermalLaw`,
`IsobaricSuperfluidLaw` and `GasEnthalpyLaw` expose their original density at the
requested specific enthalpy through `ThermalLaw::bulk_density_kg_m3`. These
adapters retain the original domain admission, phase mass fractions, additive
phase/constituent volumes, integrated thermal expansion and the gas mixture's
retained pressure and composition. They do not substitute reference mesh
density, average densities by mass fraction or fit another equation of state.
The existing opt-in material-density conduction layers can therefore derive
their normal distances from conserved mass for these material laws. Density
evaluation itself changes neither inventory nor stored geometry.

Focused G1 checks cover the additive-volume PCM and ice/water plateaus,
crystallization-branch expansion, the analytic superfluid density integral,
the isobaric mixture gas density and bracket, and two explicit PCM conduction
steps with updated resistance and conserved mass and energy. Conductivity,
finite-step accuracy, pressure work and the geometric constraints remain the
responsibilities of their existing consumers.
