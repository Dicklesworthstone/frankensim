//! Original thermal hardening observations on the native finite thermal family.
//!
//! These selectors deliberately omit the four scalar yield/hardening parameters.
//! The selected same-card Kirchhoff stress surface supplies the whole resistance
//! and its temperature derivative. The remaining thirteen thermal parameters are
//! mandatory, with the original Hencky/original-volume conventions. Source
//! preparation happens once; constitutive iterations perform no material queries.

use super::{
    AffineThermalFiniteJ2, CREEP_BOUND, CREEP_PARAMETERS, Dims, FiniteJ2Limits, FiniteJ2Parameters,
    FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters, INVERSE_K, INVERSE_K2,
    InitialStatePolicy, K, MaterialError, OVERSTRESS_PARAMETERS, PA, PA_K, PARAMETERS,
    PLASTIC_PARAMETERS, PRIMARY_PARAMETERS, QueryPoint, ResolvedAffineThermalFiniteJ2, kinetics,
    parameter, refusal, select_primary,
};
use crate::plastic::tabulated::thermal_source::{
    ThermalHardeningCurveRequest, resolve_thermal_tabulated_hardening,
};
use fs_matdb::{ClaimSelection, MaterialCard};
use std::sync::Arc;

/// Temperature-dependent tabulated Kirchhoff yielding, law/state version one.
pub const THERMAL_FINITE_J2_TABULATED_LAW_ID: &str =
    "finite-j2-tabulated-hencky-kirchhoff-thermal-affine";
/// The same energetic table with explicit Arrhenius Perzyna kinetics.
pub const THERMAL_PERZYNA_FINITE_J2_TABULATED_LAW_ID: &str =
    "finite-j2-tabulated-hencky-kirchhoff-thermal-perzyna";
/// Simultaneous tabulated yielding and Norton creep, law one/state two.
pub const THERMAL_FINITE_NORTON_J2_TABULATED_LAW_ID: &str =
    "finite-norton-j2-tabulated-hencky-kirchhoff-thermal";
/// Simultaneous tabulated yielding and strain-hardening Norton creep.
pub const THERMAL_FINITE_PRIMARY_NORTON_J2_TABULATED_LAW_ID: &str =
    "finite-primary-norton-j2-tabulated-hencky-kirchhoff-thermal";

impl ResolvedAffineThermalFiniteJ2 {
    /// Bind an attached model to original temperature/plastic-strain observations.
    ///
    /// All thermal source nodes are retained with their original selection
    /// receipts. Interpolation belongs to this explicitly selected model, with
    /// no relabeling of the original source policy or extrapolated resistance.
    /// `maximum_source_samples` bounds both observations and prepared strain
    /// knots before native preparation. The source owner is shared on cloning.
    pub fn try_from_material_card_at(
        material: &MaterialCard,
        model_index: usize,
        context: &QueryPoint,
        temperature_axis: &str,
        curve_selection: ClaimSelection,
        maximum_source_samples: usize,
    ) -> Result<Self, MaterialError> {
        let card = material.models().get(model_index).ok_or_else(|| {
            refusal("finite thermal tabulated model index is outside the material card")
        })?;
        card.validate()
            .map_err(|error| refusal(error.to_string()))?;
        let (overstress, combined, primary) = match card.law.0.as_str() {
            THERMAL_FINITE_J2_TABULATED_LAW_ID => (false, false, false),
            THERMAL_PERZYNA_FINITE_J2_TABULATED_LAW_ID => (true, false, false),
            THERMAL_FINITE_NORTON_J2_TABULATED_LAW_ID => (false, true, false),
            THERMAL_FINITE_PRIMARY_NORTON_J2_TABULATED_LAW_ID => (false, true, true),
            _ => {
                return Err(refusal(
                    "an exact finite thermal tabulated J2 selector is required",
                ));
            }
        };
        if card.law_version != 1
            || card.state_schema_version != if combined { 2 } else { 1 }
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(refusal(
                "finite thermal tabulated cards require law v1, state v2 for simultaneous creep or v1 otherwise, and virgin history",
            ));
        }
        if let Some(name) = card.parameters.keys().find(|name| {
            let name = name.as_str();
            !(PARAMETERS.contains(&name) && !PLASTIC_PARAMETERS.contains(&name)
                || overstress && OVERSTRESS_PARAMETERS.contains(&name)
                || combined && (CREEP_PARAMETERS.contains(&name) || name == CREEP_BOUND)
                || primary && PRIMARY_PARAMETERS.contains(&name))
        }) {
            return Err(refusal(format!(
                "unrecognized finite thermal tabulated parameter {name}; resistance and its thermal derivative belong to the selected surface"
            )));
        }
        let reference_temperature_k = parameter(card, "reference_temperature_k", K)?;
        let temperature_range_k = [
            parameter(card, "minimum_temperature_k", K)?,
            parameter(card, "maximum_temperature_k", K)?,
        ];
        let source = resolve_thermal_tabulated_hardening(
            material,
            context,
            &ThermalHardeningCurveRequest {
                temperature_axis: temperature_axis.to_owned(),
                reference_temperature_k,
                temperature_range_k,
                maximum_samples: maximum_source_samples,
                maximum_knots: maximum_source_samples,
            },
            curve_selection,
        )
        .map_err(|error| refusal(format!("finite thermal hardening source: {error}")))?;
        let hardening = source.hardening();
        let parameters = FiniteJ2ThermalParameters {
            reference_temperature_k,
            temperature_range_k,
            reference: FiniteJ2Parameters {
                shear_modulus_pa: parameter(card, "shear_modulus_pa", PA)?,
                bulk_modulus_pa: parameter(card, "bulk_modulus_pa", PA)?,
                yield_stress_pa: hardening.reference().initial_yield_stress_pa(),
                hardening_modulus_pa: 0.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: parameter(card, "shear_modulus_slope_pa_k", PA_K)?,
                bulk_modulus_pa_k: parameter(card, "bulk_modulus_slope_pa_k", PA_K)?,
                yield_stress_pa_k: hardening.yield_stress_slopes_pa_k()[0],
                hardening_modulus_pa_k: 0.0,
            },
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
        let family = AffineThermalFiniteJ2::new(parameters, limits)
            .and_then(|family| family.with_tabulated_hardening(hardening.clone()))
            .map_err(|error| refusal(format!("inadmissible finite thermal table: {error}")))?;
        let family = if combined {
            family.with_norton_creep(
                kinetics(card, &CREEP_PARAMETERS)?,
                parameter(card, CREEP_BOUND, Dims::NONE)?,
            )
        } else if overstress {
            family.with_overstress(kinetics(card, &OVERSTRESS_PARAMETERS)?)
        } else {
            Ok(family)
        }
        .map_err(|error| {
            refusal(format!(
                "inadmissible finite thermal tabulated flow: {error}"
            ))
        })?;
        let family = select_primary(card, family, primary)?;
        Self::admit_family(
            card,
            context,
            temperature_axis,
            family,
            Some(Arc::new(source)),
        )
    }
}
