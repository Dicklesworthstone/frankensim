//! Exact source cards for bounded affine finite-strain thermoplasticity.
//!
//! `finite-j2-hencky-kirchhoff-thermal-affine`, law/state version 1, selects
//! [`AffineThermalFiniteJ2`] with all seventeen parameters explicit. This is a
//! separate calibrated law from the fixed-temperature finite-J2 model and the
//! small-strain property-based thermoplastic model. No property conversions,
//! aliases, missing-parameter defaults or unused selectors are inferred.
//!
//! `finite-j2-hencky-kirchhoff-thermal-perzyna`, law/state version 1, explicitly
//! adds four kinetic parameters to those seventeen: reference viscosity in
//! Pa s, stress scale in Pa, dimensionless exponent at least one, and an
//! activation temperature Q/R in kelvin. Viscosity is
//! `eta(T)=eta_ref*exp((Q/R)*(1/T-1/T_ref))`, positive and representable over
//! the complete declared interval. The stress scale and exponent stay fixed.
//! These parameters control dissipation, not additional free-energy storage.
//! The selected law requires an explicit physical duration for every return.
//!
//! `finite-norton-hencky-kirchhoff-thermal`, law/state version 1, selects
//! threshold-free secondary creep. Its seventeen explicit parameters are the
//! original thermal set minus the two yield/hardening values and their slopes,
//! plus `creep_reference_viscosity_pa_s`, `creep_stress_scale_pa`,
//! `creep_exponent`, and `creep_activation_temperature_k`, with the same units
//! and Arrhenius convention as above. Yield/hardening and Perzyna kinetic keys
//! are rejected for this selector. No creep/yield combination is implied.
//! Its declared rate is `alpha_dot=(s0/eta(T))*(q/s0)^n` in equivalent KIRCHHOFF
//! stress and ORIGINAL-volume energy. Ordinary Cauchy-stress creep fits require
//! an explicitly authored conversion; they cannot silently select this law.
//! In its native state Cp^-1 denotes the inverse creep metric and alpha the
//! accumulated equivalent creep. The existing `maximum_plastic_log_strain`
//! and `maximum_accumulated_plastic_strain` keys then bound these creep states;
//! no second plastic mechanism or stored hardening is introduced.
//!
//! `finite-norton-j2-hencky-kirchhoff-thermal`, law version 1 and STATE VERSION
//! 2, selects simultaneous Norton creep and rate-independent J2 plasticity.
//! Its twenty-two parameters are the original seventeen plastic/thermal keys,
//! the four `creep_*` kinetics above, and dimensionless positive
//! `maximum_accumulated_creep_strain`. The shared inverse inelastic metric is
//! advanced by both coaxial increments. `alpha` records plastic accumulation,
//! `creep_alpha` records creep accumulation, and only alpha stores hardening
//! energy. The plastic strain bound limits alpha, the separate creep bound
//! limits creep_alpha, and the metric-log bound limits their combined metric.
//! Original state-v1 selections require a zero additional creep counter.
//! Physical yield and creep losses remain separate and both heat the caloric
//! owner once. Accepted combined states satisfy equilibrium yield; the selected
//! viscosity drives Norton creep and does not make yielding rate dependent.
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
//! State schema 1 retains Cp^-1 and accumulated inelastic strain for the selected
//! mechanism; virgin history means Cp^-1=I, alpha=0 at the reference temperature.
//!
//! The native owner supplies entropy -partial(psi)/partial(T) at fixed F/history,
//! internal energy psi+T*s and the mechanical heat-capacity correction. The
//! coupled owner must bind the original caloric chart and mass, admit positive
//! total heat capacity, and declare the caloric pressure/reference convention.
//! Source admission alone neither evolves heat nor converts numerical work
//! defects into heat. It adds no kinematic hardening, phase model, or implicit
//! conversion from a separately calibrated fluid or small-strain viscosity.

use super::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters, FiniteNortonThermalParameters, FiniteThermalJ2Material,
};
use crate::MaterialError;
use crate::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

/// Exact affine thermal model, with original-volume and Kirchhoff conventions.
pub const THERMAL_FINITE_J2_LAW_ID: &str = "finite-j2-hencky-kirchhoff-thermal-affine";
/// Explicit finite Perzyna law with Arrhenius viscosity and affine thermal energy.
pub const THERMAL_PERZYNA_FINITE_J2_LAW_ID: &str = "finite-j2-hencky-kirchhoff-thermal-perzyna";
/// Explicit threshold-free finite Norton creep with Arrhenius thermal kinetics.
pub const THERMAL_FINITE_NORTON_LAW_ID: &str = "finite-norton-hencky-kirchhoff-thermal";
/// Simultaneous finite Norton creep and energetic J2 plasticity, state schema 2.
pub const THERMAL_FINITE_NORTON_J2_LAW_ID: &str = "finite-norton-j2-hencky-kirchhoff-thermal";

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const PA_S: Dims = Dims([-1, 1, -1, 0, 0, 0]);
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
const OVERSTRESS_PARAMETERS: [&str; 4] = [
    "overstress_reference_viscosity_pa_s",
    "overstress_stress_scale_pa",
    "overstress_exponent",
    "overstress_activation_temperature_k",
];
const PLASTIC_PARAMETERS: [&str; 4] = [
    "yield_stress_pa",
    "hardening_modulus_pa",
    "yield_stress_slope_pa_k",
    "hardening_modulus_slope_pa_k",
];
const CREEP_PARAMETERS: [&str; 4] = [
    "creep_reference_viscosity_pa_s",
    "creep_stress_scale_pa",
    "creep_exponent",
    "creep_activation_temperature_k",
];
const CREEP_BOUND: &str = "maximum_accumulated_creep_strain";

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
    let parameter = card
        .parameters
        .get(name)
        .ok_or_else(|| refusal(format!("{}: missing parameter {name}", card.law.0)))?;
    if parameter.dims != dims || !parameter.value.is_finite() {
        return Err(refusal(format!(
            "finite thermal J2 parameter {name}: expected finite {dims:?}, got {parameter:?}"
        )));
    }
    Ok(parameter.value)
}

fn kinetics(
    card: &ConstitutiveModelCard,
    names: &[&str; 4],
) -> Result<FiniteJ2ThermalOverstress, MaterialError> {
    Ok(FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: parameter(card, names[0], PA_S)?,
        stress_scale_pa: parameter(card, names[1], PA)?,
        exponent: parameter(card, names[2], Dims::NONE)?,
        activation_temperature_k: parameter(card, names[3], K)?,
    })
}

fn parse_law(card: &ConstitutiveModelCard) -> Result<AffineThermalFiniteJ2, MaterialError> {
    card.validate()
        .map_err(|error| refusal(error.to_string()))?;
    let (selected_overstress, selected_norton, selected_combined) = match card.law.0.as_str() {
        THERMAL_FINITE_J2_LAW_ID => (false, false, false),
        THERMAL_PERZYNA_FINITE_J2_LAW_ID => (true, false, false),
        THERMAL_FINITE_NORTON_LAW_ID => (false, true, false),
        THERMAL_FINITE_NORTON_J2_LAW_ID => (false, false, true),
        _ => {
            return Err(refusal(format!(
                "unsupported finite thermal law {}; expected {THERMAL_FINITE_J2_LAW_ID}, {THERMAL_PERZYNA_FINITE_J2_LAW_ID}, {THERMAL_FINITE_NORTON_LAW_ID} or {THERMAL_FINITE_NORTON_J2_LAW_ID}",
                card.law.0
            )));
        }
    };
    if card.law_version != 1
        || card.state_schema_version != if selected_combined { 2 } else { 1 }
        || card.initial_state != InitialStatePolicy::ZeroInternalState
    {
        return Err(refusal(
            "finite thermal cards require law v1, state v2 for combined Norton/J2 or v1 otherwise, and virgin inelastic history at the declared reference temperature",
        ));
    }
    if let Some(name) = card.parameters.keys().find(|name| {
        let name = name.as_str();
        if selected_norton {
            !(PARAMETERS.contains(&name) && !PLASTIC_PARAMETERS.contains(&name)
                || CREEP_PARAMETERS.contains(&name))
        } else if selected_combined {
            !(PARAMETERS.contains(&name) || CREEP_PARAMETERS.contains(&name) || name == CREEP_BOUND)
        } else {
            !PARAMETERS.contains(&name)
                && !(selected_overstress && OVERSTRESS_PARAMETERS.contains(&name))
        }
    }) {
        return Err(refusal(format!(
            "unrecognized finite thermal J2 card parameter {name}"
        )));
    }
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
    if selected_norton {
        return AffineThermalFiniteJ2::new_norton(
            FiniteNortonThermalParameters {
                reference_temperature_k: parameter(card, "reference_temperature_k", K)?,
                temperature_range_k: [
                    parameter(card, "minimum_temperature_k", K)?,
                    parameter(card, "maximum_temperature_k", K)?,
                ],
                shear_modulus_pa: parameter(card, "shear_modulus_pa", PA)?,
                bulk_modulus_pa: parameter(card, "bulk_modulus_pa", PA)?,
                shear_modulus_slope_pa_k: parameter(card, "shear_modulus_slope_pa_k", PA_K)?,
                bulk_modulus_slope_pa_k: parameter(card, "bulk_modulus_slope_pa_k", PA_K)?,
                logarithmic_expansion_per_k: parameter(
                    card,
                    "logarithmic_expansion_per_k",
                    INVERSE_K,
                )?,
                logarithmic_expansion_slope_per_k2: parameter(
                    card,
                    "logarithmic_expansion_slope_per_k2",
                    INVERSE_K2,
                )?,
            },
            limits,
            kinetics(card, &CREEP_PARAMETERS)?,
        )
        .map_err(|error| {
            refusal(format!(
                "inadmissible finite thermal Norton family: {error}"
            ))
        });
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
    let family = AffineThermalFiniteJ2::new(parameters, limits)
        .map_err(|error| refusal(format!("inadmissible finite thermal J2 family: {error}")))?;
    if selected_combined {
        family
            .with_norton_creep(
                kinetics(card, &CREEP_PARAMETERS)?,
                parameter(card, CREEP_BOUND, Dims::NONE)?,
            )
            .map_err(|error| {
                refusal(format!(
                    "inadmissible finite thermal Norton/J2 family: {error}"
                ))
            })
    } else if selected_overstress {
        family
            .with_overstress(kinetics(card, &OVERSTRESS_PARAMETERS)?)
            .map_err(|error| {
                refusal(format!(
                    "inadmissible finite thermal Perzyna family: {error}"
                ))
            })
    } else {
        Ok(family)
    }
}
