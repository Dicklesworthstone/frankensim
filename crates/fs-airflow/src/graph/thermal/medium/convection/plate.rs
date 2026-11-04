//! External plate sections from one original material transport state.
//!
//! Distances retain the physical leading edge. A numerical patch never starts
//! a new boundary layer, and its area remains with the existing thermal owner.

use super::*;
use fs_convection::developing::{
    FlatPlateSection, FlatPlateSectionError, evaluate_flat_plate_section,
};

/// Original material, section geometry or source-correlation refusal.
#[derive(Debug)]
pub enum MaterialFlatPlateSectionError {
    Material(MaterialError),
    Section(FlatPlateSectionError),
    Correlation(CorrelationError),
}

impl fmt::Display for MaterialFlatPlateSectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Material(error) => write!(f, "{error}"),
            Self::Section(error) => write!(f, "{error}"),
            Self::Correlation(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for MaterialFlatPlateSectionError {}

/// A strip's actual source, cumulative endpoint evaluations and coefficient.
#[derive(Debug)]
pub struct MaterialFlatPlateSection<'a, F: FluidTransport + ?Sized> {
    fluid: &'a F,
    transport: FluidTransportState,
    start: Length,
    end: Length,
    section: FlatPlateSection,
    coefficient: Evidence<HeatTransferCoefficient>,
}

impl<F: FluidTransport + ?Sized> MaterialFlatPlateSection<'_, F> {
    /// Original material owner, including its source claims.
    pub const fn fluid(&self) -> &F {
        self.fluid
    }
    /// The single snapshot used for Re, Pr and conductivity.
    pub const fn transport(&self) -> &FluidTransportState {
        &self.transport
    }
    /// Distance from the physical upstream leading edge to the strip start.
    pub const fn start(&self) -> Length {
        self.start
    }
    /// Distance from that same edge to the strip end, not the strip's width.
    pub const fn end(&self) -> Length {
        self.end
    }
    pub const fn section(&self) -> &FlatPlateSection {
        &self.section
    }
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> {
        &self.coefficient
    }

    /// Apply this same strip coefficient to the original solid conduction row.
    /// The reference is the physical bulk temperature, not the property film.
    pub fn robin_boundary(
        &self,
        reference_temperature: Temperature,
    ) -> Result<CorrelationRobinBoundary, CorrelationError> {
        self.section.robin_boundary(
            ThermalConductivity::new(self.transport.thermal_conductivity_w_m_k()),
            self.end,
            reference_temperature,
        )
    }

    pub(super) fn into_coefficient(self) -> Evidence<HeatTransferCoefficient> {
        self.coefficient
    }
}

/// Integrate one selected external isothermal plate correlation over `[start,end]`.
///
/// Supply the actual source fluid at the declared film temperature, with the
/// same pressure and composition as the bulk fluid. `speed_m_s` is the positive
/// relative free-stream speed; both distances run downstream from its physical
/// leading edge. Reversing flow requires distances from the opposite edge.
/// The caller supplies that geometric admission; no velocity projection,
/// transition blend, hydraulic diameter or cell-centre distance is inferred.
///
/// Only the existing laminar and leading-edge-corrected turbulent plate cards
/// are supported. Each positive upstream endpoint must belong to the SAME
/// card as the downstream endpoint. At start zero this is exactly the original
/// whole-length mean. Otherwise the coefficient is a difference of cumulative
/// transfer and has no independently established local model discrepancy.
/// Summed strip heat recovers whole-plate heat only for a common isothermal wall,
/// fluid state and flow. Local temperature/property updates are a held-property
/// approximation; they do not solve a non-isothermal boundary layer.
pub fn evaluate_flat_plate_section_from_material<F: FluidTransport + ?Sized>(
    correlation: CorrelationId,
    fluid: &F,
    speed_m_s: f64,
    start: Length,
    end: Length,
) -> Result<MaterialFlatPlateSection<'_, F>, MaterialFlatPlateSectionError> {
    if !start.value().is_finite()
        || start.value() < 0.0
        || !end.value().is_finite()
        || end.value() <= start.value()
    {
        return Err(MaterialFlatPlateSectionError::Material(
            MaterialError::Parameters {
                what: "plate section needs finite 0 <= start < end measured from one leading edge"
                    .into(),
            },
        ));
    }
    let transport = fluid
        .transport_state()
        .map_err(MaterialFlatPlateSectionError::Material)?;
    let reynolds = transport
        .reynolds(speed_m_s, end.value())
        .map_err(MaterialFlatPlateSectionError::Material)?;
    let fraction = start.value() / end.value();
    if start.value() > 0.0 && fraction == 0.0 {
        return Err(MaterialFlatPlateSectionError::Material(
            MaterialError::Parameters {
                what: "positive plate section start cannot underflow to a new leading edge".into(),
            },
        ));
    }
    let section = evaluate_flat_plate_section(correlation, reynolds, transport.prandtl(), fraction)
        .map_err(MaterialFlatPlateSectionError::Section)?;
    let coefficient = section
        .heat_transfer_coefficient(
            ThermalConductivity::new(transport.thermal_conductivity_w_m_k()),
            end,
        )
        .map_err(MaterialFlatPlateSectionError::Correlation)?;
    Ok(MaterialFlatPlateSection {
        fluid,
        transport,
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
    fn liquid_sections_retain_one_source_and_recover_whole_plate_heat() {
        let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let whole = evaluate_flat_plate_section_from_material(
            CorrelationId::FlatPlateLaminarAverage,
            &fluid,
            0.1,
            Length::new(0.0),
            Length::new(0.2),
        )
        .unwrap();
        let mut integrated = 0.0;
        let mut previous_h = f64::INFINITY;
        for edge in [0.0, 0.005, 0.03, 0.11, 0.2].windows(2) {
            let section = evaluate_flat_plate_section_from_material(
                CorrelationId::FlatPlateLaminarAverage,
                &fluid,
                0.1,
                Length::new(edge[0]),
                Length::new(edge[1]),
            )
            .unwrap();
            let h = section.coefficient().value.value();
            assert!(h < previous_h);
            previous_h = h;
            integrated += h * (edge[1] - edge[0]);
            assert!(core::ptr::eq(section.fluid(), &fluid));
            assert_eq!(*section.transport(), fluid);
            assert_eq!(section.start().value(), edge[0]);
            assert_eq!(section.end().value(), edge[1]);
            let robin = section.robin_boundary(Temperature::new(310.0)).unwrap();
            assert_eq!(robin.coefficient(), section.coefficient());
        }
        assert!((integrated / (0.2 * whole.coefficient().value.value()) - 1.0).abs() < 2e-13);
        let original = evaluate_forced_from_material(
            CorrelationId::FlatPlateLaminarAverage,
            &fluid,
            0.1,
            Length::new(0.2),
            ForcedConvectionSettings::default(),
        )
        .unwrap();
        assert_eq!(whole.coefficient(), original.coefficient());
    }

    #[test]
    fn section_geometry_and_original_regime_never_silently_reset() {
        let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        for (start, end) in [
            (-0.1, 0.2),
            (0.2, 0.2),
            (0.3, 0.2),
            (f64::NAN, 0.2),
            (0.0, f64::INFINITY),
            (1e-12, 0.2),
        ] {
            assert!(
                evaluate_flat_plate_section_from_material(
                    CorrelationId::FlatPlateLaminarAverage,
                    &fluid,
                    0.1,
                    Length::new(start),
                    Length::new(end),
                )
                .is_err()
            );
        }
        // Downstream turbulent support does not authorize an upstream laminar strip.
        assert!(
            evaluate_flat_plate_section_from_material(
                CorrelationId::FlatPlateTurbulentAverage,
                &fluid,
                1.0,
                Length::new(0.1),
                Length::new(0.8),
            )
            .is_err()
        );
        assert!(
            evaluate_flat_plate_section_from_material(
                CorrelationId::FlatPlateLaminarAverage,
                &fluid,
                -0.1,
                Length::new(0.0),
                Length::new(0.2),
            )
            .is_err()
        );
    }
}
