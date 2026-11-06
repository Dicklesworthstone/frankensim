//! Variable-viscosity turbulent heat transfer in smooth circular tubes.
//!
//! Sieder-Tate uses bulk Re, Pr and conductivity, with only the viscosity
//! denominator evaluated at the fluid-facing wall. Keeping both source owners
//! prevents a film-temperature tuple from silently changing this model.

use super::*;

/// Developing turbulent sections retain the original physical tube entrance.
pub mod developing;
pub use developing::{
    MaterialTurbulentTubeSection, TurbulentTubeSectionDetails,
    evaluate_turbulent_tube_section_from_material,
};

/// Original bulk and wall-temperature fluid sources for one tube evaluation.
/// The surface source is the same FLUID at the wall temperature, not the solid.
#[derive(Debug)]
pub struct MaterialTubeConvection<'a, B: FluidTransport + ?Sized, S: FluidTransport + ?Sized> {
    bulk: MaterialConvection<'a, B>,
    surface: &'a S,
    surface_transport: FluidTransportState,
}

impl<B: FluidTransport + ?Sized, S: FluidTransport + ?Sized> MaterialTubeConvection<'_, B, S> {
    /// Bulk properties, model evidence, coefficient and native Robin lowering.
    pub const fn bulk(&self) -> &MaterialConvection<'_, B> {
        &self.bulk
    }

    /// Original fluid source at the actual wall temperature.
    pub const fn surface_fluid(&self) -> &S {
        self.surface
    }

    /// Snapshot whose viscosity alone enters the wall correction.
    pub const fn surface_transport(&self) -> &FluidTransportState {
        &self.surface_transport
    }

    pub(super) fn into_parts(
        self,
    ) -> (
        FluidTransportState,
        NusseltEvaluation,
        Evidence<HeatTransferCoefficient>,
    ) {
        self.bulk.into_parts()
    }
}

/// Evaluate the original Sieder-Tate turbulent circular-tube card.
///
/// The caller admits smooth, straight, fully developed, single-phase Newtonian
/// flow. `diameter` is the actual circular diameter and `length_over_diameter`
/// is the physical heated length ratio, not a control-volume width. The shared
/// correlation owns its Reynolds, Prandtl, viscosity-ratio and development
/// limits; no entrance, stagnant or transitional-flow correction is inferred.
///
/// Both material owners must describe the same fluid, pressure and composition
/// at bulk and wall temperatures. This generic seam retains their owners but
/// cannot prove identity. Native source-backed cell and link adapters query
/// the SAME resolver at both temperatures at every thermal trial. They retain
/// the actual wall/bulk heat drop separately from property sampling.
pub fn evaluate_tube_from_material<'a, B: FluidTransport + ?Sized, S: FluidTransport + ?Sized>(
    bulk: &'a B,
    surface: &'a S,
    speed_m_s: f64,
    diameter: Length,
    length_over_diameter: f64,
) -> Result<MaterialTubeConvection<'a, B, S>, MaterialConvectionError> {
    let transport = bulk.transport_state()?;
    let surface_transport = surface.transport_state()?;
    let reynolds = transport.reynolds(speed_m_s, diameter.value())?;
    let ratio = transport.dynamic_viscosity_pa_s() / surface_transport.dynamic_viscosity_pa_s();
    let inputs = CorrelationInputs::forced(reynolds, transport.prandtl())
        .with_length_ratio(length_over_diameter)
        .with_viscosity_ratio(ratio);
    let bulk = finish(
        CorrelationId::SiederTateTurbulent,
        bulk,
        transport,
        diameter,
        inputs,
    )?;
    Ok(MaterialTubeConvection {
        bulk,
        surface,
        surface_transport,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_tube_retains_bulk_transport_and_only_wall_viscosity() {
        let bulk = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        // All surface properties deliberately differ; only mu may affect Nu.
        let surface = FluidTransportState::try_new(850.0, 5100.0, 0.8, 0.0005).unwrap();
        let result =
            evaluate_tube_from_material(&bulk, &surface, 2.0, Length::new(0.01), 100.0).unwrap();
        let nu = 0.027 * 20_000.0_f64.powf(0.8) * 8.0_f64.cbrt() * 2.0_f64.powf(0.14);
        assert_eq!(result.bulk().nusselt().groups()["Re"], 20_000.0);
        assert_eq!(result.bulk().nusselt().groups()["Pr"], 8.0);
        assert_eq!(
            result.bulk().nusselt().groups()["mu_bulk_over_surface"],
            2.0
        );
        assert!((result.bulk().coefficient().value.value() - nu * 50.0).abs() < 1e-8);
        assert!(core::ptr::eq(result.bulk().fluid(), &bulk));
        assert!(core::ptr::eq(result.surface_fluid(), &surface));
        assert_eq!(*result.surface_transport(), surface);
        let boundary = result
            .bulk()
            .robin_boundary(Temperature::new(300.0))
            .unwrap();
        assert_eq!(
            boundary.coefficient().value,
            result.bulk().coefficient().value
        );
    }

    #[test]
    fn g0_tube_requires_wall_data_and_admitted_forced_flow() {
        let bulk = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let settings = ForcedConvectionSettings {
            length_over_hydraulic_diameter: Some(100.0),
            ..Default::default()
        };
        assert!(
            evaluate_forced_from_material(
                CorrelationId::SiederTateTurbulent,
                &bulk,
                2.0,
                Length::new(0.01),
                settings,
            )
            .is_err()
        );
        for (speed, length) in [(0.0, 100.0), (0.5, 100.0), (2.0, 5.0)] {
            assert!(
                evaluate_tube_from_material(&bulk, &bulk, speed, Length::new(0.01), length,)
                    .is_err()
            );
        }
        let cool = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.002).unwrap();
        let cooling =
            evaluate_tube_from_material(&bulk, &cool, 2.0, Length::new(0.01), 100.0).unwrap();
        assert_eq!(
            cooling.bulk().nusselt().groups()["mu_bulk_over_surface"],
            0.5
        );
    }
}
