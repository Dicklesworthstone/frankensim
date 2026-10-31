//! Material-driven convection using the existing correlation and conduction owners.
//!
//! A single transport snapshot supplies every dimensionless group AND the
//! conductivity in h = Nu k/L. This is shared by gas, mixtures, liquids and
//! other explicit implementations of FluidTransport. Geometry, correlation
//! selection and expansion remain physical inputs, never material-name switches.
//! Correlation evidence does not include propagated material uncertainty.

use core::fmt;
use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, CorrelationRobinBoundary,
    HeatTransferCoefficient, NusseltEvaluation, ThermalConductivity, ThermalDirection,
};
use fs_evidence::Evidence;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState, ThermalFluidTransport};
use fs_qty::{Length, Temperature};

pub use fs_convection::{HorizontalPlateFacing, NaturalConvectionGeometry};

/// Source-consistent laminar aiding natural/forced vertical-plate convection.
pub mod mixed;
pub use mixed::{AidingVerticalPlateSettings, evaluate_aiding_vertical_plate_from_material};

/// Geometry and heat-flow direction required by the selected forced correlation.
#[derive(Clone, Copy, Debug, Default)]
pub struct ForcedConvectionSettings {
    /// Heated length divided by hydraulic diameter, when required by the card.
    pub length_over_hydraulic_diameter: Option<f64>,
    /// Rectangular minor/major side ratio, when required by the card.
    pub aspect_ratio: Option<f64>,
    /// Dittus-Boelter heat-flow convention.
    pub direction: ThermalDirection,
}

/// Material or existing correlation refusal; no fallback coefficients are used.
#[derive(Debug)]
pub enum MaterialConvectionError {
    /// Inadmissible material transport or physical driving inputs.
    Material(MaterialError),
    /// Existing correlation, domain or dimensional conversion refusal.
    Correlation(CorrelationError),
}
impl fmt::Display for MaterialConvectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Material(error) => write!(f, "{error}"),
            Self::Correlation(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for MaterialConvectionError {}
impl From<MaterialError> for MaterialConvectionError {
    fn from(error: MaterialError) -> Self { Self::Material(error) }
}
impl From<CorrelationError> for MaterialConvectionError {
    fn from(error: CorrelationError) -> Self { Self::Correlation(error) }
}

/// One evaluated material, its exact snapshot, and its existing correlation.
/// Retaining the original source owner avoids promoting a numerical tuple into
/// material evidence. This is a frozen-property evaluation, not a validity claim
/// throughout an evolving temperature, composition or pressure field.
#[derive(Debug)]
pub struct MaterialConvection<'a, F: FluidTransport + ?Sized> {
    fluid: &'a F,
    transport: FluidTransportState,
    length: Length,
    nusselt: NusseltEvaluation,
    coefficient: Evidence<HeatTransferCoefficient>,
}

impl<F: FluidTransport + ?Sized> MaterialConvection<'_, F> {
    /// Original input, retaining any source receipts and model choices.
    #[must_use]
    pub const fn fluid(&self) -> &F { self.fluid }
    /// Exact snapshot used by both the dimensionless and dimensional stages.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// Existing correlation evaluation and its admitted group point.
    #[must_use]
    pub const fn nusselt(&self) -> &NusseltEvaluation { &self.nusselt }
    /// h = Nu k/L from this same material evaluation.
    #[must_use]
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> { &self.coefficient }

    /// Lower to the actual conduction Robin boundary without re-reading the
    /// material or asking the caller to enter a second conductivity. Reference
    /// temperature is explicit and need not equal the property's film temperature.
    ///
    /// # Errors
    /// Propagates the existing boundary's physical-temperature admission.
    pub fn robin_boundary(&self, reference_temperature: Temperature)
        -> Result<CorrelationRobinBoundary, CorrelationError>
    {
        self.nusselt.robin_boundary(
            ThermalConductivity::new(self.transport.thermal_conductivity_w_m_k()),
            self.length, reference_temperature,
        )
    }

    pub(super) fn into_parts(self)
        -> (FluidTransportState, NusseltEvaluation, Evidence<HeatTransferCoefficient>)
    {
        (self.transport, self.nusselt, self.coefficient)
    }
}

/// Evaluate a forced relation from one material state. Speed is a nonnegative
/// magnitude [m/s]; the selected card still owns geometry and regime admission.
///
/// # Errors
/// Refuses invalid material/geometry, missing card inputs and unsupported groups.
pub fn evaluate_forced_from_material<F: FluidTransport + ?Sized>(
    id: CorrelationId,
    fluid: &F,
    speed_m_s: f64,
    characteristic_length: Length,
    settings: ForcedConvectionSettings,
) -> Result<MaterialConvection<'_, F>, MaterialConvectionError> {
    let transport = fluid.transport_state()?;
    let reynolds = transport.reynolds(speed_m_s, characteristic_length.value())?;
    let mut inputs = CorrelationInputs::forced(reynolds, transport.prandtl())
        .with_direction(settings.direction);
    if let Some(ratio) = settings.length_over_hydraulic_diameter {
        inputs = inputs.with_length_ratio(ratio);
    }
    if let Some(ratio) = settings.aspect_ratio {
        inputs = inputs.with_aspect_ratio(ratio);
    }
    finish(id, fluid, transport, characteristic_length, inputs)
}

/// Natural convection for an isothermal VERTICAL plate, using the existing
/// Churchill-Chu correlation. g [m/s²] and volumetric expansion [1/K] are
/// explicit. A liquid's expansion is never guessed as 1/T.
///
/// The vertical-plate coefficient uses the magnitude of the buoyancy driving;
/// reversing hot/cold reverses flow and wall heat, not the sign of h. The actual
/// wall/reference temperatures remain with the caller and Robin boundary.
/// This symmetry is not used to infer horizontal-plate or enclosure behavior.
///
/// # Errors
/// Refuses invalid absolute temperatures, zero/out-of-domain Rayleigh numbers,
/// inadmissible material data and unrepresentable physical driving inputs.
pub fn evaluate_vertical_plate_from_material<F: FluidTransport + ?Sized>(
    fluid: &F,
    gravity_m_s2: f64,
    expansion_per_k: f64,
    wall_temperature: Temperature,
    reference_temperature: Temperature,
    characteristic_length: Length,
) -> Result<MaterialConvection<'_, F>, MaterialConvectionError> {
    evaluate_natural_from_material(
        NaturalConvectionGeometry::VerticalPlate, fluid, gravity_m_s2, expansion_per_k,
        wall_temperature, reference_temperature, characteristic_length,
    )
}

/// Natural convection for an explicitly oriented geometry, from one material
/// snapshot. The signed `beta * (T_wall - T_reference)` selects the horizontal
/// face's stability; the conductivity used for h belongs to this same snapshot.
///
/// Supply plate height, cylinder/sphere diameter, or full horizontal plate
/// area/perimeter as the characteristic length. Contact-masked exposed area cannot reconstruct
/// the planform perimeter. Properties must follow the selected card's declared
/// film convention; use the source-backed cell or film-gas adapter to re-query
/// them inside a thermal solve. The reference remains the actual bulk temperature.
///
/// # Errors
/// Refuses invalid material, geometry, driving inputs and unsupported correlation
/// ranges. Only the sphere admits Ra=0 through its sourced conductive limit Nu=2.
pub fn evaluate_natural_from_material<F: FluidTransport + ?Sized>(
    geometry: NaturalConvectionGeometry,
    fluid: &F,
    gravity_m_s2: f64,
    expansion_per_k: f64,
    wall_temperature: Temperature,
    reference_temperature: Temperature,
    characteristic_length: Length,
) -> Result<MaterialConvection<'_, F>, MaterialConvectionError> {
    let transport = fluid.transport_state()?;
    let rayleigh = natural_rayleigh(
        transport,
        gravity_m_s2,
        expansion_per_k,
        wall_temperature,
        reference_temperature,
        characteristic_length,
    )?;
    let nusselt = geometry.evaluate(rayleigh, transport.prandtl())?;
    finish_evaluation(fluid, transport, characteristic_length, nusselt)
}

/// Natural convection from the original thermal reservoir's complete material
/// state, with no separately authored expansion or ambient temperature.
///
/// The source's temperature is the reservoir temperature; its density, Cp,
/// conductivity, viscosity and volumetric expansion supply the same frozen
/// evaluation. Nonideal gases use their admitted EOS derivative. Liquid
/// expansion remains unknown unless the source owner explicitly supplies it.
/// The original material owner stays attached to the result and all existing
/// Churchill-Chu regime, geometry and model-evidence limits remain unchanged.
/// This does not re-evaluate properties at a film temperature, establish a
/// finite-temperature validity interval or model a finite-capacity reservoir.
///
/// # Errors
/// Refuses unknown expansion, inconsistent thermal/momentum snapshots, invalid
/// temperatures and the same correlation/domain failures as the explicit path.
pub fn evaluate_vertical_plate_from_thermal_material<F: ThermalFluidTransport + ?Sized>(
    fluid: &F,
    gravity_m_s2: f64,
    wall_temperature: Temperature,
    characteristic_length: Length,
) -> Result<MaterialConvection<'_, F>, MaterialConvectionError> {
    let thermal = fluid.thermal_transport_state()?;
    let transport = fluid.transport_state()?;
    if thermal.transport() != transport {
        return Err(MaterialError::Parameters {
            what: "thermal and momentum views must describe the same material state".into(),
        }
        .into());
    }
    let expansion =
        thermal
            .volumetric_expansion_per_k()
            .ok_or_else(|| MaterialError::Parameters {
                what: "natural convection requires the material source's volumetric expansion"
                    .into(),
            })?;
    let inputs = vertical_plate_inputs(
        transport,
        gravity_m_s2,
        expansion,
        wall_temperature,
        Temperature::new(thermal.temperature_k()),
        characteristic_length,
    )?;
    finish(
        CorrelationId::ChurchillChuVerticalPlate,
        fluid,
        transport,
        characteristic_length,
        inputs,
    )
}

fn vertical_plate_inputs(
    transport: FluidTransportState,
    gravity_m_s2: f64,
    expansion_per_k: f64,
    wall_temperature: Temperature,
    reference_temperature: Temperature,
    characteristic_length: Length,
) -> Result<CorrelationInputs, MaterialConvectionError> {
    let rayleigh = natural_rayleigh(transport, gravity_m_s2, expansion_per_k,
        wall_temperature, reference_temperature, characteristic_length)?;
    Ok(CorrelationInputs::natural(rayleigh.abs(), transport.prandtl()))
}

fn natural_rayleigh(
    transport: FluidTransportState,
    gravity_m_s2: f64,
    expansion_per_k: f64,
    wall_temperature: Temperature,
    reference_temperature: Temperature,
    characteristic_length: Length,
) -> Result<f64, MaterialConvectionError> {
    if [wall_temperature.value(), reference_temperature.value()]
        .into_iter().any(|value| !value.is_finite() || value <= 0.0)
    {
        return Err(MaterialError::Parameters {
            what: "natural convection needs positive finite absolute temperatures".into(),
        }.into());
    }
    Ok(transport.rayleigh(
        gravity_m_s2, expansion_per_k,
        wall_temperature.value() - reference_temperature.value(),
        characteristic_length.value(),
    )?)
}

fn finish<F: FluidTransport + ?Sized>(
    id: CorrelationId, fluid: &F, transport: FluidTransportState,
    length: Length, inputs: CorrelationInputs,
) -> Result<MaterialConvection<'_, F>, MaterialConvectionError> {
    let nusselt = fs_convection::evaluate(id, inputs)?;
    finish_evaluation(fluid, transport, length, nusselt)
}

fn finish_evaluation<F: FluidTransport + ?Sized>(
    fluid: &F, transport: FluidTransportState, length: Length, nusselt: NusseltEvaluation,
) -> Result<MaterialConvection<'_, F>, MaterialConvectionError> {
    let coefficient = nusselt.heat_transfer_coefficient(
        ThermalConductivity::new(transport.thermal_conductivity_w_m_k()), length,
    )?;
    Ok(MaterialConvection { fluid, transport, length, nusselt, coefficient })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::gas::{GasSpec, GasState};
    use fs_qty::semantic::{Composition, CompositionBasis};

    #[test]
    fn pressure_recomputes_both_reynolds_and_heat_transfer() {
        let spec = GasSpec::dry_air_ussa1976();
        let a = GasState::try_new(&spec, 300.0, 100000.0).unwrap();
        let b = GasState::try_new(&spec, 300.0, 200000.0).unwrap();
        let first = evaluate_forced_from_material(CorrelationId::FlatPlateLaminarAverage,
            &a, 1.0, Length::new(0.1), ForcedConvectionSettings::default()).unwrap();
        let second = evaluate_forced_from_material(CorrelationId::FlatPlateLaminarAverage,
            &b, 1.0, Length::new(0.1), ForcedConvectionSettings::default()).unwrap();
        assert!((second.nusselt().groups()["Re"] / first.nusselt().groups()["Re"] - 2.0).abs() < 1e-12);
        assert!((second.coefficient().value.value() / first.coefficient().value.value() - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn gas_mixture_reaches_the_existing_conduction_boundary() {
        let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.98, 0.02]).unwrap();
        let gas = GasState::try_new_mixture(
            &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()],
            &composition, 300.0, 100000.0, Some(1),
        ).unwrap();
        let result = evaluate_forced_from_material(CorrelationId::FlatPlateLaminarAverage,
            &gas, 1.0, Length::new(0.1), ForcedConvectionSettings::default()).unwrap();
        let boundary = result.robin_boundary(Temperature::new(290.0)).unwrap();
        assert_eq!(boundary.coefficient().value, result.coefficient().value);
        match boundary.boundary_condition() {
            fs_conduction::ThermalBc::Robin { htc, t_ref } => {
                assert_eq!(htc.at(0), result.coefficient().value.value());
                assert_eq!(t_ref.at(0), 290.0);
            }
            other => panic!("expected Robin boundary, got {other:?}"),
        }
    }

    #[test]
    fn liquid_coefficients_are_not_replaced_by_air_defaults() {
        let liquid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let result = evaluate_forced_from_material(CorrelationId::FlatPlateLaminarAverage,
            &liquid, 1.0, Length::new(0.1), ForcedConvectionSettings::default()).unwrap();
        let direct = fs_convection::evaluate(CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(100000.0, 8.0)).unwrap();
        assert!((result.nusselt().evidence().value - direct.evidence().value).abs() < 1e-10);
        assert!((result.coefficient().value.value() - direct.evidence().value * 5.0).abs() < 1e-9);
    }

    #[test]
    fn vertical_hot_and_cold_walls_keep_positive_h_and_opposite_heat_flow() {
        let fluid = FluidTransportState::try_new(1.0, 1000.0, 0.02, 2e-5).unwrap();
        let hot = evaluate_vertical_plate_from_material(&fluid, 9.81, 1.0 / 300.0,
            Temperature::new(320.0), Temperature::new(300.0), Length::new(0.1)).unwrap();
        let cold = evaluate_vertical_plate_from_material(&fluid, 9.81, 1.0 / 300.0,
            Temperature::new(280.0), Temperature::new(300.0), Length::new(0.1)).unwrap();
        assert_eq!(hot.coefficient().value, cold.coefficient().value);
        let h = cold.coefficient().value.value();
        assert!(h > 0.0);
        assert!(h * (280.0 - 300.0) < 0.0);
        assert!(h * (320.0 - 300.0) > 0.0);
        assert_eq!(hot.nusselt().groups()["Ra"], fluid.rayleigh(9.81, 1.0 / 300.0, 20.0, 0.1).unwrap());
    }

    #[test]
    fn thermal_material_owns_natural_buoyancy_and_ambient_temperature() {
        let source = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 150000.0).unwrap();
        for wall in [280.0, 320.0] {
            let material = evaluate_vertical_plate_from_thermal_material(
                &source,
                9.81,
                Temperature::new(wall),
                Length::new(0.1),
            )
            .unwrap();
            let explicit = evaluate_vertical_plate_from_material(
                &source,
                9.81,
                1.0 / 300.0,
                Temperature::new(wall),
                Temperature::new(300.0),
                Length::new(0.1),
            )
            .unwrap();
            assert_eq!(material.nusselt().groups(), explicit.nusselt().groups());
            assert_eq!(material.coefficient().value, explicit.coefficient().value);
            assert!(core::ptr::eq(material.fluid(), &source));
        }
        let mut stale = source;
        stale.density *= 2.0;
        assert!(
            evaluate_vertical_plate_from_thermal_material(
                &stale,
                9.81,
                Temperature::new(320.0),
                Length::new(0.1)
            )
            .is_err()
        );
    }

    #[test]
    fn thermal_material_does_not_invent_expansion_or_bypass_zero_driving_refusal() {
        use fs_material::fluid_transport::ThermalFluidTransportState;
        let transport = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let unknown = ThermalFluidTransportState::try_new(transport, 300.0, None).unwrap();
        assert!(matches!(
            evaluate_vertical_plate_from_thermal_material(
                &unknown,
                9.81,
                Temperature::new(320.0),
                Length::new(0.1)
            ),
            Err(MaterialConvectionError::Material(_))
        ));
        let zero = ThermalFluidTransportState::try_new(transport, 300.0, Some(0.0)).unwrap();
        assert!(matches!(
            evaluate_vertical_plate_from_thermal_material(
                &zero,
                9.81,
                Temperature::new(320.0),
                Length::new(0.1)
            ),
            Err(MaterialConvectionError::Correlation(_))
        ));
    }

    #[test]
    fn thermal_material_refuses_disagreeing_fluid_views() {
        use fs_material::fluid_transport::ThermalFluidTransportState;
        struct Disagreeing;
        impl FluidTransport for Disagreeing {
            fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
                FluidTransportState::try_new(1.0, 1000.0, 0.02, 2e-5)
            }
        }
        impl ThermalFluidTransport for Disagreeing {
            fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
                ThermalFluidTransportState::try_new(
                    FluidTransportState::try_new(2.0, 1000.0, 0.02, 2e-5)?,
                    300.0,
                    Some(1.0 / 300.0),
                )
            }
        }
        assert!(matches!(
            evaluate_vertical_plate_from_thermal_material(
                &Disagreeing,
                9.81,
                Temperature::new(320.0),
                Length::new(0.1)
            ),
            Err(MaterialConvectionError::Material(_))
        ));
    }

    #[test]
    fn boundary_lowering_uses_one_snapshot_even_for_stateful_adapters() {
        struct Counted(std::cell::Cell<usize>);
        impl FluidTransport for Counted {
            fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
                self.0.set(self.0.get() + 1);
                FluidTransportState::try_new(1.0, 1000.0, 0.02, 2e-5)
            }
        }
        let fluid = Counted(std::cell::Cell::new(0));
        let result = evaluate_forced_from_material(CorrelationId::FlatPlateLaminarAverage,
            &fluid, 1.0, Length::new(0.1), ForcedConvectionSettings::default()).unwrap();
        result.robin_boundary(Temperature::new(300.0)).unwrap();
        assert_eq!(fluid.0.get(), 1);
    }

    #[test]
    fn missing_geometry_and_invalid_driving_still_refuse() {
        let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        assert!(evaluate_forced_from_material(CorrelationId::CircularDuctHausen,
            &fluid, 0.01, Length::new(0.01), ForcedConvectionSettings::default()).is_err());
        assert!(evaluate_forced_from_material(CorrelationId::FlatPlateLaminarAverage,
            &fluid, -1.0, Length::new(0.1), ForcedConvectionSettings::default()).is_err());
        assert!(evaluate_vertical_plate_from_material(&fluid, 9.81, 2e-4,
            Temperature::new(-1.0), Temperature::new(300.0), Length::new(0.1)).is_err());
        assert!(evaluate_vertical_plate_from_material(&fluid, 9.81, 2e-4,
            Temperature::new(300.0), Temperature::new(300.0), Length::new(0.1)).is_err());
    }
}

/// Material-driven boundary powers on the shared finite-cell thermal history.
pub mod deposition;

/// Reciprocal convection with the existing finite gas, enthalpy and wall dynamics.
pub mod finite_gas;
