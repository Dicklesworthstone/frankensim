//! Exact source binding for finite porous thermomechanics.
//!
//! The affine selector requires all eighteen base parameters, with Kirchhoff
//! stress and ORIGINAL-volume energetic coefficients. The Voce selector adds
//! Q_ref [Pa], b [1] and dQ/dT [Pa/K]. Pressure sensitivity and b stay fixed.
//! Both selections retain Cp^-1 and alpha, including plastic volume, under
//! law/state version 1. A virgin state is at the declared thermal reference.
//!
//! Exactly one typed absolute-temperature coordinate evolves, over the entire
//! source-supported interval. The original card/context remain attached; all
//! other constrained coordinates are fixed admitted external conditions. An
//! isothermal porous card or small-strain calibration cannot select this law.
//! Caloric data and reference mass are separate required coupling inputs.

use super::{
    AffineThermalFinitePorous, FinitePorousTemperatureSlopes, FinitePorousThermalParameters,
    FiniteThermalPorousMaterial,
};
use crate::MaterialError;
use crate::plastic::finite::porous::{FinitePorousLimits, FinitePorousParameters};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

pub const THERMAL_FINITE_POROUS_LAW_ID: &str = "finite-porous-hencky-kirchhoff-thermal-affine";
pub const THERMAL_FINITE_POROUS_VOCE_LAW_ID: &str =
    "finite-porous-voce-hencky-kirchhoff-thermal-affine";

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const PA_K: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const K: Dims = Dims([0, 0, 0, 1, 0, 0]);
const INVERSE_K: Dims = Dims([0, 0, 0, -1, 0, 0]);
const INVERSE_K2: Dims = Dims([0, 0, 0, -2, 0, 0]);
const PARAMETERS: [&str; 18] = [
    "reference_temperature_k",
    "minimum_temperature_k",
    "maximum_temperature_k",
    "shear_modulus_pa",
    "bulk_modulus_pa",
    "yield_stress_pa",
    "hardening_modulus_pa",
    "pressure_sensitivity",
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
const VOCE_PARAMETERS: [&str; 3] = [
    "voce_saturation_pa",
    "voce_rate_per_strain",
    "voce_saturation_slope_pa_k",
];

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedAffineThermalFinitePorous {
    card: ConstitutiveModelCard,
    context: QueryPoint,
    temperature_axis: String,
    law: AffineThermalFinitePorous,
}

impl ResolvedAffineThermalFinitePorous {
    /// Admit a complete typed temperature family without querying during Newton.
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
                "finite thermal porous source requires one typed absolute-temperature axis and exact reference-temperature context",
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
                    "finite thermal porous source must cover the complete native temperature interval with exact quantity descriptors",
                ));
            }
        }
        for axis in card.validity.bounds().keys() {
            if axis == temperature_axis {
                continue;
            }
            let quantity = context.axis_quantities().get(axis).copied();
            let another_temperature =
                quantity.is_some_and(|q| q.dims() == K) || axis == "T" && quantity.is_none();
            if another_temperature
                || !crate::state_point::is_fixed_external_coordinate(
                    axis,
                    quantity,
                    context.axes().get(axis).copied(),
                )
            {
                return Err(refusal(format!(
                    "finite thermal porous source axis {axis} is not a supported fixed external condition; only the named temperature may evolve"
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
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.card
    }
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    pub const fn law(&self) -> &AffineThermalFinitePorous {
        &self.law
    }
}
impl FiniteThermalPorousMaterial for ResolvedAffineThermalFinitePorous {
    fn thermal_finite_porous(&self) -> &AffineThermalFinitePorous {
        &self.law
    }
}

fn refusal(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn parameter(card: &ConstitutiveModelCard, name: &str, dims: Dims) -> Result<f64, MaterialError> {
    let value = card
        .parameters
        .get(name)
        .ok_or_else(|| refusal(format!("{}: missing parameter {name}", card.law.0)))?;
    if value.dims != dims || !value.value.is_finite() {
        return Err(refusal(format!(
            "finite thermal porous parameter {name}: expected finite {dims:?}, got {value:?}"
        )));
    }
    Ok(value.value)
}
fn parse_law(card: &ConstitutiveModelCard) -> Result<AffineThermalFinitePorous, MaterialError> {
    card.validate()
        .map_err(|error| refusal(error.to_string()))?;
    let voce = match card.law.0.as_str() {
        THERMAL_FINITE_POROUS_LAW_ID => false,
        THERMAL_FINITE_POROUS_VOCE_LAW_ID => true,
        _ => {
            return Err(refusal(format!(
                "unsupported finite thermal porous law {}",
                card.law.0
            )));
        }
    };
    if card.law_version != 1
        || card.state_schema_version != 1
        || card.initial_state != InitialStatePolicy::ZeroInternalState
    {
        return Err(refusal(
            "finite thermal porous cards require law/state v1 and virgin plastic history",
        ));
    }
    if let Some(name) = card.parameters.keys().find(|name| {
        !PARAMETERS.contains(&name.as_str()) && !(voce && VOCE_PARAMETERS.contains(&name.as_str()))
    }) {
        return Err(refusal(format!(
            "unrecognized finite thermal porous parameter {name}; all selected coefficients and limits must be explicit"
        )));
    }
    let law = AffineThermalFinitePorous::new(
        FinitePorousThermalParameters {
            reference_temperature_k: parameter(card, "reference_temperature_k", K)?,
            temperature_range_k: [
                parameter(card, "minimum_temperature_k", K)?,
                parameter(card, "maximum_temperature_k", K)?,
            ],
            reference: FinitePorousParameters {
                shear_modulus_pa: parameter(card, "shear_modulus_pa", PA)?,
                bulk_modulus_pa: parameter(card, "bulk_modulus_pa", PA)?,
                yield_stress_pa: parameter(card, "yield_stress_pa", PA)?,
                hardening_modulus_pa: parameter(card, "hardening_modulus_pa", PA)?,
                pressure_sensitivity: parameter(card, "pressure_sensitivity", Dims::NONE)?,
            },
            slopes: FinitePorousTemperatureSlopes {
                shear_modulus_pa_k: parameter(card, "shear_modulus_slope_pa_k", PA_K)?,
                bulk_modulus_pa_k: parameter(card, "bulk_modulus_slope_pa_k", PA_K)?,
                yield_stress_pa_k: parameter(card, "yield_stress_slope_pa_k", PA_K)?,
                hardening_modulus_pa_k: parameter(card, "hardening_modulus_slope_pa_k", PA_K)?,
            },
            logarithmic_expansion_per_k: parameter(card, "logarithmic_expansion_per_k", INVERSE_K)?,
            logarithmic_expansion_slope_per_k2: parameter(
                card,
                "logarithmic_expansion_slope_per_k2",
                INVERSE_K2,
            )?,
        },
        FinitePorousLimits {
            minimum_jacobian: parameter(card, "minimum_mechanical_jacobian", Dims::NONE)?,
            maximum_elastic_log_strain: parameter(card, "maximum_elastic_log_strain", Dims::NONE)?,
            maximum_plastic_log_strain: parameter(card, "maximum_plastic_log_strain", Dims::NONE)?,
            maximum_accumulated_plastic_strain: parameter(
                card,
                "maximum_accumulated_plastic_strain",
                Dims::NONE,
            )?,
        },
    )
    .map_err(|error| {
        refusal(format!(
            "inadmissible finite porous thermal family: {error}"
        ))
    })?;
    if voce {
        law.with_voce_hardening(
            parameter(card, "voce_saturation_pa", PA)?,
            parameter(card, "voce_rate_per_strain", Dims::NONE)?,
            parameter(card, "voce_saturation_slope_pa_k", PA_K)?,
        )
        .map_err(|error| {
            refusal(format!(
                "inadmissible finite porous thermal Voce hardening: {error}"
            ))
        })
    } else {
        Ok(law)
    }
}
