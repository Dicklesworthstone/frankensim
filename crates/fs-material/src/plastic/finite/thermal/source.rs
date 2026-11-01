//! Exact source cards for bounded affine finite-strain thermoplasticity.
//!
//! `finite-j2-hencky-kirchhoff-thermal-affine`, law/state version 1, selects
//! [`AffineThermalFiniteJ2`] with all seventeen parameters explicit. This is a
//! separate calibrated law from the fixed-temperature finite-J2 model and the
//! small-strain property-based thermoplastic model. No property conversions,
//! aliases, missing-parameter defaults or unused selectors are inferred.
//!
//! The four reference coefficients `shear_modulus_pa`, `bulk_modulus_pa`,
//! `yield_stress_pa`, `hardening_modulus_pa` are in Pa. Their four derivatives
//! `shear_modulus_slope_pa_k`, `bulk_modulus_slope_pa_k`,
//! `yield_stress_slope_pa_k`, `hardening_modulus_slope_pa_k` are in Pa/K.
//! `reference_temperature_k`, `minimum_temperature_k` and
//! `maximum_temperature_k` are absolute kelvin temperatures. Every coefficient
//! is affine about the declared reference temperature over the complete span.
//!
//! `logarithmic_expansion_per_k` and `logarithmic_expansion_slope_per_k2`
//! specify the affine derivative `a(T) = d ln(theta) / dT`, in 1/K and 1/K².
//! Thus `theta(T) = exp(integral(a(T) dT))`. A generic engineering expansion
//! coefficient or mean secant coefficient is not silently this logarithmic law.
//! `minimum_mechanical_jacobian` bounds det(F/theta), not det(F).
//! `maximum_elastic_log_strain`, `maximum_plastic_log_strain` and
//! `maximum_accumulated_plastic_strain` retain their native norm/history meaning.
//! All four domain parameters are dimensionless and remain checked natively.
//!
//! Moduli, hardening storage and dissipation use the ORIGINAL reference volume;
//! yield is in equivalent Kirchhoff stress. There is no additional theta³
//! multiplier on the free energy. At a stress-free expanded F=theta I, the
//! corresponding infinitesimal current-volume shear/bulk moduli are mu/theta³
//! and kappa/theta³. Ordinary measured Cauchy moduli require an explicitly
//! authored conversion or fit before they can calibrate this model card.
//!
//! Exactly one typed absolute-temperature coordinate evolves. Its source domain
//! must cover the entire native temperature span, including the reference.
//! Other declared coordinates remain fixed external conditions. The exact
//! original card, typed reference context and temperature-axis name persist.
//! State schema 1 retains Cp^-1 and accumulated plastic strain; virgin zero
//! plastic strain means Cp^-1=I, alpha=0 at the reference temperature.
//!
//! The native owner supplies entropy -partial(psi)/partial(T) at fixed F/history,
//! internal energy psi+T*s and the mechanical heat-capacity correction. The
//! coupled owner must bind the original caloric chart and mass, admit positive
//! total heat capacity, and declare the caloric pressure/reference convention.
//! Source admission alone neither evolves heat nor converts numerical work
//! defects into heat. It adds no kinematic hardening, phase or creep model.

use super::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters,
    FiniteThermalJ2Material,
};
use crate::MaterialError;
use crate::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

/// Exact affine thermal model, with original-volume and Kirchhoff conventions.
pub const THERMAL_FINITE_J2_LAW_ID: &str = "finite-j2-hencky-kirchhoff-thermal-affine";

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const PA_K: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const K: Dims = Dims([0, 0, 0, 1, 0, 0]);
const INVERSE_K: Dims = Dims([0, 0, 0, -1, 0, 0]);
const INVERSE_K2: Dims = Dims([0, 0, 0, -2, 0, 0]);
const PARAMETERS: [&str; 17] = [
    "reference_temperature_k",
    "minimum_temperature_k",
    "maximum_temperature_k",
    "shear_modulus_pa",
    "bulk_modulus_pa",
    "yield_stress_pa",
    "hardening_modulus_pa",
    "shear_modulus_slope_pa_k",
    "bulk_modulus_slope_pa_k",
    "yield_stress_slope_pa_k",
    "hardening_modulus_slope_pa_k",
    "logarithmic_expansion_per_k",
    "logarithmic_expansion_slope_per_k2",
    "minimum_mechanical_jacobian",
    "maximum_elastic_log_strain",
    "maximum_plastic_log_strain",
    "maximum_accumulated_plastic_strain",
];

/// Original source family, fixed external context and explicit evolving T axis.
///
/// Frozen-temperature trials remain with the native family. Its thermal
/// mapping cannot be bypassed through an isothermal `FiniteJ2Material` adapter.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedAffineThermalFiniteJ2 {
    card: ConstitutiveModelCard,
    context: QueryPoint,
    temperature_axis: String,
    law: AffineThermalFiniteJ2,
}

impl ResolvedAffineThermalFiniteJ2 {
    /// Admit one explicit family over its complete source-supported T interval.
    ///
    /// `context` must name the reference temperature with exact static absolute
    /// temperature semantics. The card must declare that same typed axis and
    /// cover both native endpoints. Its other applicability bounds may constrain
    /// fixed pressure/composition/moisture/concentration; a second temperature,
    /// mechanical coordinate, frequency or rate needs another owner and refuses.
    /// Extra unconstrained context is preserved without granting source coverage.
    pub fn try_from_model_card_at(
        card: &ConstitutiveModelCard,
        context: &QueryPoint,
        temperature_axis: &str,
    ) -> Result<Self, MaterialError> {
        let law = parse_law(card)?;
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if temperature_axis.is_empty()
            || context.axis_quantities().get(temperature_axis) != Some(&absolute)
            || context.axes().get(temperature_axis)
                != Some(&law.parameters().reference_temperature_k)
            || card.validity.axis_quantities().get(temperature_axis) != Some(&absolute)
            || !card.validity.bounds().contains_key(temperature_axis)
        {
            return Err(refusal(
                "finite thermal J2 requires one explicit typed absolute-temperature source axis and its exact reference-temperature context",
            ));
        }
        for temperature in law.parameters().temperature_range_k {
            let endpoint = context
                .clone()
                .with(temperature_axis, temperature)
                .map_err(|error| refusal(error.to_string()))?;
            if !card
                .validity
                .contains_typed(endpoint.axes(), endpoint.axis_quantities())
            {
                return Err(refusal(
                    "finite thermal J2 source context must cover the complete native temperature interval with exact quantity descriptors",
                ));
            }
        }
        for axis in card.validity.bounds().keys() {
            if axis == temperature_axis {
                continue;
            }
            let quantity = context.axis_quantities().get(axis).copied();
            let another_temperature = quantity.is_some_and(|quantity| quantity.dims() == K)
                || axis == "T" && quantity.is_none();
            if another_temperature
                || !crate::state_point::is_fixed_external_coordinate(
                    axis,
                    quantity,
                    context.axes().get(axis).copied(),
                )
            {
                return Err(refusal(format!(
                    "finite thermal J2 source axis {axis} is not a supported fixed external condition; only the named temperature coordinate may evolve"
                )));
            }
        }
        Ok(Self {
            card: card.clone(),
            context: context.clone(),
            temperature_axis: temperature_axis.to_owned(),
            law,
        })
    }

    /// Complete original model, including source hashes and calibration domain.
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.card
    }

    /// Exact typed context at the declared reference temperature.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }

    /// The single explicitly owned evolving absolute-temperature coordinate.
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }

    /// Native thermodynamic family; all original-volume conventions stay intact.
    pub const fn law(&self) -> &AffineThermalFiniteJ2 {
        &self.law
    }
}

impl FiniteThermalJ2Material for ResolvedAffineThermalFiniteJ2 {
    fn thermal_finite_j2(&self) -> &AffineThermalFiniteJ2 {
        &self.law
    }
}

fn refusal(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

fn parameter(card: &ConstitutiveModelCard, name: &str, dims: Dims) -> Result<f64, MaterialError> {
    let parameter = card.parameters.get(name).ok_or_else(|| {
        refusal(format!(
            "{THERMAL_FINITE_J2_LAW_ID}: missing parameter {name}"
        ))
    })?;
    if parameter.dims != dims || !parameter.value.is_finite() {
        return Err(refusal(format!(
            "finite thermal J2 parameter {name}: expected finite {dims:?}, got {parameter:?}"
        )));
    }
    Ok(parameter.value)
}

fn parse_law(card: &ConstitutiveModelCard) -> Result<AffineThermalFiniteJ2, MaterialError> {
    card.validate()
        .map_err(|error| refusal(error.to_string()))?;
    if card.law.0 != THERMAL_FINITE_J2_LAW_ID {
        return Err(refusal(format!(
            "unsupported finite thermal J2 law {}; expected {THERMAL_FINITE_J2_LAW_ID}",
            card.law.0
        )));
    }
    if card.law_version != 1
        || card.state_schema_version != 1
        || card.initial_state != InitialStatePolicy::ZeroInternalState
    {
        return Err(refusal(
            "finite thermal J2 cards require law/state v1 and virgin plastic history at the declared reference temperature",
        ));
    }
    if let Some(name) = card
        .parameters
        .keys()
        .find(|name| !PARAMETERS.contains(&name.as_str()))
    {
        return Err(refusal(format!(
            "unrecognized finite thermal J2 card parameter {name}"
        )));
    }
    let reference = FiniteJ2Parameters {
        shear_modulus_pa: parameter(card, "shear_modulus_pa", PA)?,
        bulk_modulus_pa: parameter(card, "bulk_modulus_pa", PA)?,
        yield_stress_pa: parameter(card, "yield_stress_pa", PA)?,
        hardening_modulus_pa: parameter(card, "hardening_modulus_pa", PA)?,
    };
    let slopes = FiniteJ2TemperatureSlopes {
        shear_modulus_pa_k: parameter(card, "shear_modulus_slope_pa_k", PA_K)?,
        bulk_modulus_pa_k: parameter(card, "bulk_modulus_slope_pa_k", PA_K)?,
        yield_stress_pa_k: parameter(card, "yield_stress_slope_pa_k", PA_K)?,
        hardening_modulus_pa_k: parameter(card, "hardening_modulus_slope_pa_k", PA_K)?,
    };
    let parameters = FiniteJ2ThermalParameters {
        reference_temperature_k: parameter(card, "reference_temperature_k", K)?,
        temperature_range_k: [
            parameter(card, "minimum_temperature_k", K)?,
            parameter(card, "maximum_temperature_k", K)?,
        ],
        reference,
        slopes,
        logarithmic_expansion_per_k: parameter(card, "logarithmic_expansion_per_k", INVERSE_K)?,
        logarithmic_expansion_slope_per_k2: parameter(
            card,
            "logarithmic_expansion_slope_per_k2",
            INVERSE_K2,
        )?,
    };
    let limits = FiniteJ2Limits {
        minimum_jacobian: parameter(card, "minimum_mechanical_jacobian", Dims::NONE)?,
        maximum_elastic_log_strain: parameter(card, "maximum_elastic_log_strain", Dims::NONE)?,
        maximum_plastic_log_strain: parameter(card, "maximum_plastic_log_strain", Dims::NONE)?,
        maximum_accumulated_plastic_strain: parameter(
            card,
            "maximum_accumulated_plastic_strain",
            Dims::NONE,
        )?,
    };
    AffineThermalFiniteJ2::new(parameters, limits)
        .map_err(|error| refusal(format!("inadmissible finite thermal J2 family: {error}")))
}
