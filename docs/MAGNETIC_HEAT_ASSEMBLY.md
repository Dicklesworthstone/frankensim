# Material heat transfer in a magnetic assembly

`PreparedHeatNetwork::step_nonlinear_winding_candidate` connects the existing
nonlinear winding/core/circuit solver to the complete prepared material heat
assembly. The original bulk-only and complete-heat entries use the same private
field/enthalpy transaction.

At every proposed thermal endpoint, the original owners evaluate bulk thermal
conductivity, area-specific thermal-contact resistance, surface emissivity,
conductor resistivity and the selected magnetic law. The assembly retains the
original source cards, exact electrical property keys, selected claims, material
volumes and heat destinations. Changing these source data changes the actual
current, field and temperature solution; no material-name lookup or independent
temperature effect is added.

Pass the existing `MagneticThermalLaw` bindings and assembly limits through
`magnetic_thermal` to include source-bound permeability, field entropy and
reversible magnetic heat. With `None`, the original fixed magnetic-law path,
including saturation, remains responsible for admitting its temperature domain.

The returned `NonlinearWindingEnthalpyCandidate<HeatNetworkSample>` retains
thermal-interface and radiosity results alongside magnetic sources, current,
field history and the optional physical caloric balance. Commit its thermal and
circuit histories together using `into_parts`.

## Energy ownership

Physical winding, core and ordinary graph-resistor losses enter the original
thermal receivers once. Backward-Euler numerical dissipation stays in its own
account and is never deposited as material heat. The existing magnetic caloric
owner supplies reversible heat from field parameter work and the change in
temperature times entropy.

Radiative exchange in a closed enclosure is internal. Both native free-energy
and physical internal-energy gates exclude its numerical net from external
supply. `caloric_balance()` reports the corrected physical external heat. For
the ordinary report, subtract reversible magnetic heat and
`dt * report.source.network.radiation_net_power_w()` from
`report.non_electrical_heat_j` to recover external non-electrical supply.

## Physical limits and regression

Field geometry, contact areas, enclosure view factors and magnetic-region/cell
correspondence remain fixed. Magnetic cells use the original fully solid,
constant-density caloric charts. Temperature-dependent magnetic storage retains
the existing affine-permeability and positive combined heat-capacity domain.
Electrical series contacts refuse because thermal contact resistance cannot
supply their voltage jump or electrical work. Source and cancellation failures
return no partially advanced history.

The regression module
`crates/fs-couple/tests/nonlinear_winding_enthalpy/heat.rs` compares an
independent two-temperature voltage/current/energy reference, continued physical
steps, source-only contact/emissivity substitutions, exact reduction to the
bulk-only path and refused endpoints. The model still requires time and spatial
refinement; numerical checks do not establish measured-device qualification.
