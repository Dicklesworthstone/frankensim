//! Whitaker sphere convection with the source's two property temperatures.
//!
//! Re, Pr and k use the free stream. Only the viscosity correction queries
//! fluid at the sphere surface temperature. A mean-film tuple changes this
//! model; the existing generic single-state forced adapter cannot supply it.

use super::*;

/// Both original material owners and the coefficient computed from them.
/// The surface material is FLUID at the wall temperature, not the solid body.
#[derive(Debug)]
pub struct MaterialSphereConvection<'a, B: FluidTransport + ?Sized, S: FluidTransport + ?Sized> {
    bulk: MaterialConvection<'a, B>,
    surface: &'a S,
    surface_transport: FluidTransportState,
}

impl<B: FluidTransport + ?Sized, S: FluidTransport + ?Sized> MaterialSphereConvection<'_, B, S> {
    /// Bulk properties, source evidence, Nusselt result and Robin lowering.
    pub const fn bulk(&self) -> &MaterialConvection<'_, B> {
        &self.bulk
    }

    /// Original fluid source evaluated at the sphere surface temperature.
    pub const fn surface_fluid(&self) -> &S {
        self.surface
    }

    /// Exact surface snapshot whose viscosity supplied the correction.
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

/// Forced convection over one isolated isothermal sphere, using its diameter.
///
/// Both states must describe the same single-phase Newtonian fluid, pressure
/// and composition, evaluated at bulk and wall temperatures respectively.
/// This generic numerical seam retains both owners but cannot establish their
/// identity; source-pinned gas adapters perform the temperature query through
/// the original material bank. Liquids use their own material law, never gas
/// transport presets. Each supplied transport view is sampled exactly once.
///
/// The shared Whitaker card owns all validity bounds, uncertainty and input
/// derivatives. In particular, `mu_bulk / mu_surface < 1` is refused rather
/// than clamped. No natural-convection term, additional conductive link,
/// evaporation correction, or stagnant-flow extrapolation is added.
pub fn evaluate_sphere_from_material<'a, B: FluidTransport + ?Sized, S: FluidTransport + ?Sized>(
    bulk: &'a B,
    surface: &'a S,
    speed_m_s: f64,
    diameter: Length,
) -> Result<MaterialSphereConvection<'a, B, S>, MaterialConvectionError> {
    let transport = bulk.transport_state()?;
    let surface_transport = surface.transport_state()?;
    let reynolds = transport.reynolds(speed_m_s, diameter.value())?;
    let ratio = transport.dynamic_viscosity_pa_s() / surface_transport.dynamic_viscosity_pa_s();
    let inputs =
        CorrelationInputs::forced(reynolds, transport.prandtl()).with_viscosity_ratio(ratio);
    let bulk = finish(
        CorrelationId::WhitakerSphere,
        bulk,
        transport,
        diameter,
        inputs,
    )?;
    Ok(MaterialSphereConvection {
        bulk,
        surface,
        surface_transport,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_liquid_sphere_uses_bulk_re_pr_k_and_only_surface_viscosity() {
        let bulk = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        // Deliberately different surface rho, Cp and k expose a film substitution.
        let surface = FluidTransportState::try_new(970.0, 4300.0, 0.6, 0.0005).unwrap();
        let result =
            evaluate_sphere_from_material(&bulk, &surface, 0.1, Length::new(0.01)).unwrap();
        let nu = 2.0
            + (0.4 * 1000.0_f64.sqrt() + 0.06 * 1000.0_f64.powf(2.0 / 3.0))
                * 8.0_f64.powf(0.4)
                * 2.0_f64.powf(0.25);
        let evaluation = result.bulk();
        assert!((evaluation.nusselt().groups()["Re"] - 1000.0).abs() < 1e-12);
        assert_eq!(evaluation.nusselt().groups()["Pr"], 8.0);
        assert_eq!(evaluation.nusselt().groups()["mu_bulk_over_surface"], 2.0);
        assert!((evaluation.coefficient().value.value() - nu * 50.0).abs() < 1e-9);
        assert!(core::ptr::eq(evaluation.fluid(), &bulk));
        assert!(core::ptr::eq(result.surface_fluid(), &surface));
        assert_eq!(*result.surface_transport(), surface);
        let robin = evaluation.robin_boundary(Temperature::new(300.0)).unwrap();
        assert_eq!(robin.coefficient().value, evaluation.coefficient().value);
    }

    #[test]
    fn g0_sphere_requires_surface_data_and_refuses_outside_the_source_range() {
        let bulk = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        assert!(
            evaluate_forced_from_material(
                CorrelationId::WhitakerSphere,
                &bulk,
                0.1,
                Length::new(0.01),
                ForcedConvectionSettings::default()
            )
            .is_err()
        );
        for mu in [0.002, 0.0001] {
            let surface = FluidTransportState::try_new(1000.0, 4000.0, 0.5, mu).unwrap();
            assert!(
                evaluate_sphere_from_material(&bulk, &surface, 0.1, Length::new(0.01)).is_err()
            );
        }
        assert!(evaluate_sphere_from_material(&bulk, &bulk, 0.0, Length::new(0.01)).is_err());
    }
}
