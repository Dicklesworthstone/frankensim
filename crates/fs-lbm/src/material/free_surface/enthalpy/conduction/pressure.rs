//! Prescribed mechanical pressure on the complete original thermal transaction.

use super::*;

/// External gas pressure held over one accepted thermal/flow interval.
///
/// Pressure is sampled at exposed liquid cell centers in metres. The scaling
/// must already carry an explicit absolute-pressure reference. This is an
/// independent prescribed mechanical load, not the thermodynamic pressure of
/// an evolved atmosphere; the caller owns compatibility with the pinned liquid
/// and thermal gas-source conditions. No pressure work is deposited as heat.
#[derive(Clone, Copy)]
pub struct FreeSurfaceGasPressure<'a> {
    /// Maximum absolute fractional lattice-density excursion, finite in (0,1).
    /// The bound includes the native capillary pressure jump.
    pub maximum_density_excursion: f64,
    /// Positive absolute gas pressure [Pa] at a physical cell center [m].
    pub pressure_pa: &'a dyn Fn([f64; 3]) -> f64,
}

/// Thermal and mechanical boundary inputs for one native liquid interval.
///
/// Empty slices retain insulated boundaries. Existing geometry, exposed-area,
/// source and boundary-overlap admission applies to every supplied slice.
/// A solved interface owns its gas and radiation exchange; do not also assign
/// those same exchanges through `environment` or `radiation`. Surface-stress
/// bindings still require their complete solved interface geometry.
#[derive(Clone, Copy, Default)]
pub struct FreeSurfaceThermalBoundaries<'a> {
    /// Source-bound resistance at the supplied liquid/solid contacts.
    pub interfaces: &'a [FreeSurfaceBodyInterface<'a>],
    /// Prescribed gas-film heat exchange at declared exposed liquid faces.
    pub environment: &'a [FreeSurfaceGasContact<'a>],
    /// Direct liquid-cell radiation from prescribed surroundings.
    pub radiation: &'a [FreeSurfaceRadiationFace],
    /// Solved interface temperatures for gas/radiation and surface stress.
    pub surfaces: &'a [FreeSurfaceInterfaceContact<'a>],
    /// Original explicit or feature-selected implicit surface controls.
    pub surface_solve: FreeSurfaceInterfaceSolve,
    /// Optional normal pressure on the original missing-gas-link boundary.
    pub gas_pressure: Option<FreeSurfaceGasPressure<'a>>,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Advance conducting liquid with combined heat boundaries and gas pressure.
    ///
    /// Uses the original thermal timeline, source refresh, native momentum and
    /// actual mass-route enthalpy transport exactly once. Finite bodies and
    /// their retained radiation boundaries, source-bound contact resistance,
    /// prescribed gas-film and radiation exchange, solved interface temperature,
    /// temperature-dependent capillarity, and bound electrical/Joule/Lorentz
    /// physics retain their existing owners and admission. Selecting viscous
    /// heating adds only the original molecular stress-power source.
    ///
    /// `heat` is supplied uniformly over this interval. Mechanical pressure is
    /// held over its native flow step and uses the same Pascal/coordinate
    /// conversion as [`MaterialFreeSurface3::step_with_gas_pressure`]. A common
    /// datum shift leaves the imposed load unchanged; `None` retains ordinary
    /// thermal arithmetic. This is the existing first-order thermal/flow split,
    /// with unchanged constant-reference-density and fixed source-pressure
    /// caloric scope, not a compressive-heating or total-energy closure.
    ///
    /// Any failed thermal, pressure, source, geometry, mass-transport or
    /// cancellation check leaves accepted liquid, finite-body, electrical and
    /// material history and the liquid clock unchanged together. This entry
    /// does not evolve a finite gas reservoir or bypass surface/electrical
    /// requirements to obtain a pressure-only step.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_thermal_boundaries(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        boundaries: FreeSurfaceThermalBoundaries<'_>,
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        poll(gate)?;
        let pressure =
            boundaries
                .gas_pressure
                .map(|load| {
                    if !load.maximum_density_excursion.is_finite()
                        || load.maximum_density_excursion <= 0.0
                        || load.maximum_density_excursion >= 1.0
                    {
                        return Err(Error::Invalid(
                            "gas-pressure density-excursion bound must be finite and in (0,1)",
                        ));
                    }
                    let reference = self.flow.scaling.reference_pressure_pa.ok_or(
                        Error::Invalid(
                            "absolute gas loading requires an explicit fixed pressure reference",
                        ),
                    )?;
                    Ok((
                        load.maximum_density_excursion,
                        super::super::super::pressure::lattice_pressure_field(
                            reference,
                            self.flow.scaling.cell_size_m,
                            self.flow.scaling.pressure_unit_pa,
                            load.pressure_pa,
                        ),
                    ))
                })
                .transpose()?;
        self.step_with_thermal_pressure(
            faces,
            heat,
            controls,
            bodies,
            contacts,
            boundaries,
            with_viscous_heating,
            pressure
                .as_ref()
                .map(|(maximum, field)| (*maximum, GasPressureField::Cell(field))),
            gate,
        )
    }

    /// Advance all original thermal boundaries with distinct pressure per gas link.
    ///
    /// Positive absolute pressure [Pa] is frozen over the native flow interval
    /// and uses [`super::super::super::MaterialGasPressureLink3`]'s cell centers
    /// and link midpoint. The supplied `boundaries.gas_pressure` must be `None`:
    /// one mechanical pressure field owns each reconstructed gas population.
    ///
    /// Conduction, finite bodies, source viscosity, surface/electrical bindings,
    /// optional molecular viscous heat and mass/enthalpy transport use exactly
    /// the transaction of [`Self::step_with_thermal_boundaries`]. A uniform
    /// link field reproduces its cell-pressure result. All state, source and
    /// clock publication remains atomic; gas topology, gas compression and
    /// pressure-volume heat are not inferred from the numerical free surface.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_thermal_boundaries_and_gas_link_pressure(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        boundaries: FreeSurfaceThermalBoundaries<'_>,
        with_viscous_heating: bool,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(super::super::super::MaterialGasPressureLink3) -> f64,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        poll(gate)?;
        if boundaries.gas_pressure.is_some() {
            return Err(Error::Invalid(
                "cell and link pressure cannot both own the gas boundary",
            ));
        }
        if !maximum_density_excursion.is_finite()
            || maximum_density_excursion <= 0.0
            || maximum_density_excursion >= 1.0
        {
            return Err(Error::Invalid(
                "gas-pressure density-excursion bound must be finite and in (0,1)",
            ));
        }
        let reference = self
            .flow
            .scaling
            .reference_pressure_pa
            .ok_or(Error::Invalid(
                "absolute gas loading requires an explicit fixed pressure reference",
            ))?;
        let field = super::super::super::pressure::lattice_link_pressure_field(
            reference,
            self.flow.scaling.cell_size_m,
            self.flow.scaling.pressure_unit_pa,
            pressure_pa,
        );
        self.step_with_thermal_pressure(
            faces,
            heat,
            controls,
            bodies,
            contacts,
            boundaries,
            with_viscous_heating,
            Some((maximum_density_excursion, GasPressureField::Link(&field))),
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_thermal_pressure(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        boundaries: FreeSurfaceThermalBoundaries<'_>,
        with_viscous_heating: bool,
        gas_pressure: Option<(f64, GasPressureField<'_>)>,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        let applied = if with_viscous_heating {
            self.viscous_sources(gate)?
        } else {
            Vec::new()
        };
        let conduction = self.advance_with_thermal_boundaries_and_pressure(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            boundaries.interfaces,
            boundaries.environment,
            boundaries.radiation,
            boundaries.surfaces,
            boundaries.surface_solve,
            None,
            gas_pressure,
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}
