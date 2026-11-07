//! Thermally developing laminar circular sections from one bulk material state.

use super::*;
use fs_convection::developing::{CircularDuctSection, evaluate_circular_duct_section};

/// Original bulk material, physical heated entrance and conservative Hausen section.
#[derive(Debug)]
pub struct MaterialLaminarTubeSection<'a, F: FluidTransport + ?Sized> {
    fluid: &'a F,
    transport: FluidTransportState,
    diameter: Length,
    start: Length,
    end: Length,
    section: CircularDuctSection,
    coefficient: Evidence<HeatTransferCoefficient>,
}

impl<F: FluidTransport + ?Sized> MaterialLaminarTubeSection<'_, F> {
    pub const fn fluid(&self) -> &F {
        self.fluid
    }
    /// The single snapshot supplying Re, Pr and dimensional conductivity.
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    /// Actual circular inside diameter, not the axial section length.
    pub const fn diameter(&self) -> Length {
        self.diameter
    }
    pub const fn start(&self) -> Length {
        self.start
    }
    pub const fn end(&self) -> Length {
        self.end
    }
    pub const fn section(&self) -> &CircularDuctSection {
        &self.section
    }
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> {
        &self.coefficient
    }

    pub(crate) fn into_parts(
        self,
    ) -> (
        FluidTransportState,
        NusseltEvaluation,
        Evidence<HeatTransferCoefficient>,
        TubeSectionDetails,
    ) {
        let (_, _, upstream, downstream, nusselt) = self.section.into_parts();
        (
            self.transport,
            downstream,
            self.coefficient,
            TubeSectionDetails {
                diameter: self.diameter,
                start: self.start,
                end: self.end,
                upstream,
                nusselt,
            },
        )
    }

    /// Apply this same source-bearing coefficient to the native conduction row.
    /// The reference temperature is the actual bulk fluid, not a separate film.
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

/// Integrate the original Hausen CWT relation over physical `[start,end]`.
///
/// The velocity profile must be hydrodynamically developed at the heated inlet;
/// this is thermal development, not a simultaneously developing inlet model.
/// One bulk source supplies Re, Pr and k. Each positive endpoint must satisfy
/// the original card's domain, and reversing flow requires physical distances
/// from its opposite heated entrance. No wall-viscosity correction is inferred.
///
/// Uniform bulk properties, flow and wall temperature make adjacent sections
/// additive. Refreshing properties locally remains a held-property model, not
/// a solution of nonuniform entrance history. Positive-start sections retain
/// both original means and have no independently validated local discrepancy.
pub fn evaluate_laminar_tube_section_from_material<F: FluidTransport + ?Sized>(
    fluid: &F,
    speed_m_s: f64,
    diameter: Length,
    start: Length,
    end: Length,
) -> Result<MaterialLaminarTubeSection<'_, F>, MaterialConvectionError> {
    if !start.value().is_finite()
        || start.value() < 0.0
        || !end.value().is_finite()
        || end.value() <= start.value()
    {
        return Err(MaterialError::Parameters {
            what: "laminar tube section needs finite 0 <= start < end from one heated entrance"
                .into(),
        }
        .into());
    }
    let transport = fluid.transport_state()?;
    let reynolds = transport.reynolds(speed_m_s, diameter.value())?;
    let start_ratio = start.value() / diameter.value();
    if start.value() > 0.0 && start_ratio == 0.0 {
        return Err(MaterialError::Parameters {
            what: "positive tube section start cannot underflow to a new heated entrance".into(),
        }
        .into());
    }
    let section = evaluate_circular_duct_section(
        reynolds,
        transport.prandtl(),
        start_ratio,
        end.value() / diameter.value(),
    )?;
    let coefficient = section.heat_transfer_coefficient(
        ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
        diameter,
    )?;
    Ok(MaterialLaminarTubeSection {
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
    fn laminar_material_sections_preserve_one_source_and_the_entrance_integral() {
        let fluid = Counted {
            state: FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap(),
            calls: Cell::new(0),
        };
        let d = Length::new(0.01);
        let evaluate = |start, end| {
            evaluate_laminar_tube_section_from_material(
                &fluid,
                0.1,
                d,
                Length::new(start),
                Length::new(end),
            )
            .unwrap()
        };
        let whole = evaluate(0.0, 1.0);
        let mean = |gz: f64| 3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0));
        let expected = mean(80.0) * 50.0; // Re=1000, Pr=8, L/D=100.
        assert!((whole.coefficient().qoi - expected).abs() < 3e-12 * expected);
        assert!(core::ptr::eq(whole.fluid(), &fluid));
        assert_eq!(*whole.transport(), fluid.state);
        assert_eq!(fluid.calls.get(), 1);
        let mut integrated = 0.0;
        for pair in [0.0, 0.05, 0.2, 0.5, 1.0].windows(2) {
            let section = evaluate(pair[0], pair[1]);
            let robin = section.robin_boundary(Temperature::new(300.0)).unwrap();
            assert_eq!(robin.coefficient().qoi, section.coefficient().qoi);
            integrated += section.coefficient().qoi * (pair[1] - pair[0]);
        }
        assert!((integrated - expected).abs() < 3e-12 * expected);
        assert_eq!(fluid.calls.get(), 5); // Robin lowering never re-queries material.
        let downstream = evaluate(0.2, 1.0);
        let expected_downstream = (100.0 * mean(80.0) - 20.0 * mean(400.0)) / 80.0 * 50.0;
        assert!(
            (downstream.coefficient().qoi - expected_downstream).abs()
                < 3e-12 * expected_downstream
        );
        let (_, original_mean, coefficient, details) = downstream.into_parts();
        assert_eq!(original_mean.card().id, CorrelationId::CircularDuctHausen);
        assert_eq!(details.upstream.unwrap().groups()["L_over_Dh"], 20.0);
        assert_eq!((details.start.value(), details.end.value()), (0.2, 1.0));
        assert!((details.nusselt - coefficient.qoi / 50.0).abs() < 3e-12);
        assert!(coefficient.model.discrepancy_rel.is_infinite());
    }

    #[test]
    fn laminar_sections_reject_stagnation_transition_and_restarted_invalid_entrances() {
        let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        for (speed, start, end) in [
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (0.1, -0.1, 1.0),
            (0.1, 0.2, 0.2),
            (0.1, 1e-12, 1.0),
        ] {
            assert!(
                evaluate_laminar_tube_section_from_material(
                    &fluid,
                    speed,
                    Length::new(0.01),
                    Length::new(start),
                    Length::new(end),
                )
                .is_err()
            );
        }
    }
}
