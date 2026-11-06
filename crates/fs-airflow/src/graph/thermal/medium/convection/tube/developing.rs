//! Smooth turbulent tube sections using one original bulk transport state.

use super::*;
use fs_convection::developing::{
    TurbulentCircularDuctSection, evaluate_turbulent_circular_duct_section,
};

/// Original bulk material, circular geometry and cumulative entrance transfer.
/// The finite section's coefficient can enter the existing conduction Robin row.
#[derive(Debug)]
pub struct MaterialTurbulentTubeSection<'a, F: FluidTransport + ?Sized> {
    fluid: &'a F,
    transport: FluidTransportState,
    diameter: Length,
    start: Length,
    end: Length,
    section: TurbulentCircularDuctSection,
    coefficient: Evidence<HeatTransferCoefficient>,
}

impl<F: FluidTransport + ?Sized> MaterialTurbulentTubeSection<'_, F> {
    /// The original owner of the bulk transport properties.
    pub const fn fluid(&self) -> &F {
        self.fluid
    }
    /// The single snapshot used for Re, Pr and dimensional conductivity.
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    /// Actual circular inside diameter, not the section length.
    pub const fn diameter(&self) -> Length {
        self.diameter
    }
    /// Upstream distance from the physical entrance.
    pub const fn start(&self) -> Length {
        self.start
    }
    /// Downstream distance from that same entrance.
    pub const fn end(&self) -> Length {
        self.end
    }
    /// Both original cumulative means and their conservative section average.
    pub const fn section(&self) -> &TurbulentCircularDuctSection {
        &self.section
    }
    /// The actual section coefficient, with its own model discrepancy boundary.
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> {
        &self.coefficient
    }

    /// Lower this same section into native conduction. The reference is the
    /// actual bulk temperature; no thermal storage or solver is added here.
    pub fn robin_boundary(
        &self,
        reference_temperature: Temperature,
    ) -> Result<CorrelationRobinBoundary, CorrelationError> {
        self.section.robin_boundary(
            ThermalConductivity::new(self.transport.thermal_conductivity_w_m_k()),
            self.diameter,
            reference_temperature,
        )
    }
}

/// Evaluate a physical `[start,end]` section of a smooth circular turbulent tube.
///
/// The shared Gnielinski developing model uses bulk Re, Pr and conductivity
/// with its explicit constant-property correction K=1. The caller must justify
/// negligible cross-section property variation; this is not Sieder--Tate's
/// wall-viscosity model. The source and its transport tuple are retained, and
/// sampled exactly once. Source updates may re-evaluate a local frozen-property
/// approximation, but additivity requires the SAME bulk properties and thermal
/// driving across every section of the original tube.
///
/// Distances are measured from the actual flow entrance. A new wall cell is
/// not a new entrance, and reversing flow requires rebinding the distances.
/// The original mean-card domain applies at both positive endpoints. Roughness,
/// transition, buoyancy, strong property variation and a local experimental
/// uncertainty bound are outside this model.
pub fn evaluate_turbulent_tube_section_from_material<F: FluidTransport + ?Sized>(
    fluid: &F,
    speed_m_s: f64,
    diameter: Length,
    start: Length,
    end: Length,
) -> Result<MaterialTurbulentTubeSection<'_, F>, MaterialConvectionError> {
    if !start.value().is_finite()
        || start.value() < 0.0
        || !end.value().is_finite()
        || end.value() <= start.value()
    {
        return Err(MaterialError::Parameters {
            what: "turbulent tube section needs finite 0 <= start < end from one physical entrance"
                .into(),
        }
        .into());
    }
    let transport = fluid.transport_state()?;
    let reynolds = transport.reynolds(speed_m_s, diameter.value())?;
    let start_ratio = start.value() / diameter.value();
    if start.value() > 0.0 && start_ratio == 0.0 {
        return Err(MaterialError::Parameters {
            what: "positive tube section start cannot underflow to a new entrance".into(),
        }
        .into());
    }
    let section = evaluate_turbulent_circular_duct_section(
        reynolds,
        transport.prandtl(),
        start_ratio,
        end.value() / diameter.value(),
    )?;
    let coefficient = section.heat_transfer_coefficient(
        ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
        diameter,
    )?;
    Ok(MaterialTurbulentTubeSection {
        fluid,
        transport,
        diameter,
        start,
        end,
        section,
        coefficient,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    struct Counted {
        state: FluidTransportState,
        calls: Cell<usize>,
    }
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.calls.set(self.calls.get() + 1);
            Ok(self.state)
        }
    }

    #[test]
    fn g1_source_and_diameter_reach_native_robin_with_conservative_section_heat() {
        let fluid = Counted {
            state: FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap(),
            calls: Cell::new(0),
        };
        let d = Length::new(0.01);
        let whole = evaluate_turbulent_tube_section_from_material(
            &fluid,
            2.0,
            d,
            Length::new(0.0),
            Length::new(1.0),
        )
        .unwrap();
        assert_eq!(fluid.calls.get(), 1);
        assert!(core::ptr::eq(whole.fluid(), &fluid));
        assert_eq!(*whole.transport(), fluid.state);
        assert_eq!(whole.diameter(), d);
        // Independent Eq. (2), corrected Eq. (6), K=1 at Re=20000, Pr=8.
        let f = (1.8 * 20_000.0_f64.log10() - 1.5).powi(-2);
        let developed = (f / 8.0) * 19_000.0 * 8.0
            / (1.0 + 12.7 * (f / 8.0).sqrt() * (8.0_f64.powf(2.0 / 3.0) - 1.0));
        let expected = developed * (1.0 + 100.0_f64.powf(-2.0 / 3.0)) * 50.0;
        assert!((whole.coefficient().qoi - expected).abs() < 3e-12 * expected);
        let mut sum = 0.0;
        for pair in [0.0, 0.1, 0.2, 0.5, 1.0].windows(2) {
            let section = evaluate_turbulent_tube_section_from_material(
                &fluid,
                2.0,
                d,
                Length::new(pair[0]),
                Length::new(pair[1]),
            )
            .unwrap();
            sum += section.coefficient().qoi * (pair[1] - pair[0]);
            let robin = section.robin_boundary(Temperature::new(300.0)).unwrap();
            assert_eq!(robin.coefficient().qoi, section.coefficient().qoi);
            assert_eq!(section.start().value(), pair[0]);
            assert_eq!(section.end().value(), pair[1]);
        }
        assert!((sum - whole.coefficient().qoi).abs() < 3e-12 * sum);
        // Sampling and Robin lowering never re-query the material source.
        assert_eq!(fluid.calls.get(), 5);
        let reset: f64 = [0.1, 0.1, 0.3, 0.5]
            .into_iter()
            .map(|width| {
                evaluate_turbulent_tube_section_from_material(
                    &fluid,
                    2.0,
                    d,
                    Length::new(0.0),
                    Length::new(width),
                )
                .unwrap()
                .coefficient()
                .qoi * width
            })
            .sum();
        assert!(
            reset > 1.05 * sum,
            "resetting the entrance must overpredict cooling"
        );
    }

    #[test]
    fn g0_tube_sections_keep_actual_material_and_physical_entrance_domain() {
        let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        for (speed, d, start, end) in [
            (0.0, 0.01, 0.0, 1.0),
            (0.1, 0.01, 0.0, 1.0),
            (2.0, 0.0, 0.0, 1.0),
            (2.0, 0.01, 0.0, 0.05),
            (2.0, 0.01, 0.05, 1.0),
            (2.0, 0.01, 1.0, 1.0),
            (2.0, 0.01, -0.1, 1.0),
            (2.0, 0.01, 0.0, f64::NAN),
        ] {
            assert!(
                evaluate_turbulent_tube_section_from_material(
                    &fluid,
                    speed,
                    Length::new(d),
                    Length::new(start),
                    Length::new(end),
                )
                .is_err()
            );
        }
        let high_pr = FluidTransportState::try_new(1000.0, 4000.0, 0.01, 0.001).unwrap();
        assert!(
            evaluate_turbulent_tube_section_from_material(
                &high_pr,
                2.0,
                Length::new(0.01),
                Length::new(0.0),
                Length::new(1.0),
            )
            .is_err()
        );
    }
}
