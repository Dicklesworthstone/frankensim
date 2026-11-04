//! Finite rectangular cooling-channel sections from one original fluid state.

use super::*;
use fs_convection::developing::{RectangularDuctSection, evaluate_rectangular_duct_section};

/// Original material, hydraulic geometry, cumulative means and section transfer.
#[derive(Debug)]
pub struct MaterialRectangularDuctSection<'a, F: FluidTransport + ?Sized> {
    fluid: &'a F,
    transport: FluidTransportState,
    hydraulic_diameter: Length,
    start: Length,
    end: Length,
    section: RectangularDuctSection,
    coefficient: Evidence<HeatTransferCoefficient>,
}

impl<F: FluidTransport + ?Sized> MaterialRectangularDuctSection<'_, F> {
    pub const fn fluid(&self) -> &F {
        self.fluid
    }
    /// The single snapshot supplying Re, Pr and dimensional conductivity.
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    pub const fn hydraulic_diameter(&self) -> Length {
        self.hydraulic_diameter
    }
    pub const fn start(&self) -> Length {
        self.start
    }
    pub const fn end(&self) -> Length {
        self.end
    }
    pub const fn section(&self) -> &RectangularDuctSection {
        &self.section
    }
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> {
        &self.coefficient
    }

    /// Use this same coefficient on the original conduction row. The reference
    /// temperature is the actual bulk fluid, not the property evaluation state.
    pub fn robin_boundary(
        &self,
        reference_temperature: Temperature,
    ) -> Result<CorrelationRobinBoundary, CorrelationError> {
        self.section.robin_boundary(
            ThermalConductivity::new(self.transport.thermal_conductivity_w_m_k()),
            self.hydraulic_diameter,
            reference_temperature,
        )
    }
}

/// Integrate the original rectangular CWT developing table over physical
/// `[start,end]` distances from one heated entrance, using actual bulk speed.
///
/// This table admits exactly aspect ratio 0.5 and Pr=0.72. The adapter does not
/// change the source's Pr to fit the table, infer a new entrance at a cell edge,
/// or replace hydraulic diameter with section length. Reversing flow requires
/// distances from its actual opposite entrance. The caller owns that geometry.
/// All source transport properties are sampled once, including the conductivity
/// used to convert Nu to h. At start zero the original whole-length coefficient
/// is preserved. A positive-start coefficient is a conservative section of the
/// mean interpolant, with unbounded local model discrepancy and no local slope
/// validation claim. Variable wall/fluid states remain a local held-property
/// approximation to the uniform-state source relation.
pub fn evaluate_rectangular_duct_section_from_material<F: FluidTransport + ?Sized>(
    fluid: &F,
    speed_m_s: f64,
    hydraulic_diameter: Length,
    aspect_ratio: f64,
    start: Length,
    end: Length,
) -> Result<MaterialRectangularDuctSection<'_, F>, MaterialConvectionError> {
    if !start.value().is_finite()
        || start.value() < 0.0
        || !end.value().is_finite()
        || end.value() <= start.value()
    {
        return Err(MaterialConvectionError::Material(
            MaterialError::Parameters {
                what: "duct section needs finite 0 <= start < end from one heated entrance".into(),
            },
        ));
    }
    let transport = fluid
        .transport_state()
        .map_err(MaterialConvectionError::Material)?;
    let reynolds = transport
        .reynolds(speed_m_s, hydraulic_diameter.value())
        .map_err(MaterialConvectionError::Material)?;
    let start_ratio = start.value() / hydraulic_diameter.value();
    if start.value() > 0.0 && start_ratio == 0.0 {
        return Err(MaterialConvectionError::Material(
            MaterialError::Parameters {
                what: "positive duct section start cannot underflow to a new entrance".into(),
            },
        ));
    }
    let section = evaluate_rectangular_duct_section(
        reynolds,
        transport.prandtl(),
        aspect_ratio,
        start_ratio,
        end.value() / hydraulic_diameter.value(),
    )
    .map_err(MaterialConvectionError::Correlation)?;
    let coefficient = section
        .heat_transfer_coefficient(
            ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
            hydraulic_diameter,
        )
        .map_err(MaterialConvectionError::Correlation)?;
    Ok(MaterialRectangularDuctSection {
        fluid,
        transport,
        hydraulic_diameter,
        start,
        end,
        section,
        coefficient,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn original_fluid_and_hydraulic_geometry_preserve_channel_heat_under_subdivision() {
        let fluid = FluidTransportState::try_new(1.0, 720.0, 1.0, 0.001).unwrap();
        let diameter = Length::new(0.1);
        let whole = evaluate_rectangular_duct_section_from_material(
            &fluid,
            1.0,
            diameter,
            0.5,
            Length::new(0.0),
            Length::new(1.2),
        )
        .unwrap();
        assert!(core::ptr::eq(whole.fluid(), &fluid));
        assert_eq!(whole.transport(), &fluid);
        let original = evaluate_forced_from_material(
            CorrelationId::RectangularDuctLaminarCwtDevelopingPr072,
            &fluid,
            1.0,
            diameter,
            ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(1.2 / 0.1),
                aspect_ratio: Some(0.5),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(whole.coefficient().qoi, original.coefficient().qoi);
        assert_eq!(
            whole.coefficient().provenance,
            original.coefficient().provenance
        );
        let mut sum = 0.0;
        for pair in [0.0, 0.04, 0.1, 0.2, 0.5, 1.2].windows(2) {
            let section = evaluate_rectangular_duct_section_from_material(
                &fluid,
                1.0,
                diameter,
                0.5,
                Length::new(pair[0]),
                Length::new(pair[1]),
            )
            .unwrap();
            sum += section.coefficient().qoi * (pair[1] - pair[0]);
            let robin = section.robin_boundary(Temperature::new(300.0)).unwrap();
            assert_eq!(robin.coefficient().qoi, section.coefficient().qoi);
        }
        assert!((sum - 1.2 * whole.coefficient().qoi).abs() < 2e-12 * sum);
    }

    #[test]
    fn actual_prandtl_and_positive_upstream_distance_cannot_be_fabricated() {
        let fluid = FluidTransportState::try_new(1.0, 720.0, 1.0, 0.001).unwrap();
        let other_pr = FluidTransportState::try_new(1.0, 710.0, 1.0, 0.001).unwrap();
        let evaluate = |fluid, aspect, start| {
            evaluate_rectangular_duct_section_from_material(
                fluid,
                1.0,
                Length::new(0.1),
                aspect,
                Length::new(start),
                Length::new(1.0),
            )
            .map(|receipt| receipt.coefficient().qoi)
        };
        assert!(evaluate(&fluid, 0.5, 0.0).is_ok());
        assert!(evaluate(&fluid, 0.5, 0.01).is_err()); // upstream Gz=720
        assert!(evaluate(&other_pr, 0.5, 0.0).is_err());
        assert!(evaluate(&fluid, 0.4, 0.0).is_err());
    }
}
