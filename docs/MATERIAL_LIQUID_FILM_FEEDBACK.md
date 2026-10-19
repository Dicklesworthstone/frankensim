# Live material temperature in squeeze-film mechanics

`LiquidFilmImpactSystem` connects the existing `LiquidTemperatureLaw` to the
original prepared impact equation and finite material heat stores. A data-only
change to the liquid's viscosity curve changes the solved Reynolds pressure,
force, mechanical trajectory and dissipated work. A data-only change to its
caloric chart changes temperature and, through that same viscosity law, the
subsequent motion.

The public API lives at
`fs_couple::render::plate::impact::material::heating::liquid_film`.
Construct the ordinary `HeatedImpactSystem` with its pressure islands and
explicit `SqueezeFilmHeatRoute`s, then attach `LiquidFilmBinding`s naming an
island, its finite liquid reservoir and the source-pinned temperature law.
Each binding also supplies a maximum relative viscosity change per sample;
an excessive change refuses and requires a smaller mechanical timestep.

## One physical history

At accepted time n, the film coefficient is the source law evaluated at its
reservoir's accepted temperature. The existing fs-phs/fs-tribo solve holds that
coefficient over the mechanical sample while jointly solving all body motion,
pressure, contact and dissipation. Heat is the original film's work at the
actual discrete-gradient velocity, with its existing explicit partition.

The finite thermal owner advances with that heat. The resulting temperature
must resolve the same material claims before either motion or heat publishes.
Only then does the existing pressure graph receive its new viscosity. Its
geometry, work-conjugate ports, limits, coordinates, momentum, contact history
and prepared numerical workspace remain intact. This is safe under load because
the incompressible pressure-eliminated film stores no energy. The update does
not apply to a compressible film's stored gas state.

`step_with_thermal_network` uses the original `DepositionSimulation` and
`DepositionObserver` substeps for conduction and other sources. It shares the
existing mechanical-heating proposal and accounting code. The network's
observer, finite stores, motion, contact sources, film coefficients and liquid
source receipts all publish together. No film work is redeposited, and work from
the completed sample is never recomputed using the hotter endpoint viscosity.

## Scope and evidence

The film has negligible thermal storage and shares the temperature of an
explicit finite, well-mixed liquid reservoir. The reservoir is the thermal
inventory; this is not an additional finite film mass/enthalpy transport model.
The same original material card must identify its caloric and transport data.
The caller supplies their consistent fixed pressure/composition context and
the spatial lumping assumption. Transport Cp is not used for another heat step:
the original caloric chart remains the energy owner.

Only fully liquid, constant-density, fixed-reservoir-volume states are admitted.
Thermal expansion, spatial flash temperature, phase support, cavitation, slip,
rarefaction and compressible storage require their own coupled models. Viscosity
is held for each mechanical sample, so this is first-order constitutive
splitting; the relative-change allowance is not a temporal-error certificate.

Focused G0/G1 tests exercise pressure/force/work/tangent scaling on the unchanged
Poiseuille graph, repeated coefficient refresh without cumulative rescaling,
and unchanged gap/pressure refusals. G1/G3 coupled regressions check independent
midpoint work, data-only heat-capacity feedback into motion, conduction on the
original network's substeps, and total energy/momentum. G4 regressions check
endpoint source-domain refusal and late cancellation preserve accepted owners.
These synthetic materials prove numerical wiring, not measured lubricant accuracy.
