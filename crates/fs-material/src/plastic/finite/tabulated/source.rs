//! Original material-card curves on the native finite J2 owner.
//!
//! The exact `finite-j2-tabulated-hencky-kirchhoff`, law/state v1 selector
//! supplies original-volume Hencky shear/bulk moduli and the four native
//! finite-domain bounds. Its initial yield and complete energetic hardening
//! come from one explicitly selected **Kirchhoff** yield-versus-accumulated-
//! plastic-strain claim on the SAME material card. Cauchy/engineering curves,
//! additional linear hardening and implicit extrapolation are not admitted.
//!
//! The attached model, selected curve, every curve-knot receipt and fixed
//! external query remain available. Source preparation occurs once; all
//! mechanical responses use the original allocation-free finite return.
//! These immutable fixed-condition sources do not authorize a temperature
//! refresh or a transfer of existing history to a different calibration.

use crate::MaterialError;
use crate::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Material, FiniteJ2Parameters};
use crate::plastic::tabulated::source::{
    HardeningCurveRequest, ResolvedTabulatedHardening, YieldStressMeasure,
    resolve_tabulated_hardening,
};
use fs_matdb::{
    ClaimSelection, ConstitutiveModelCard, InitialStatePolicy, MaterialCard, QueryPoint,
};
use fs_qty::{Dims, Pressure};

/// Explicit energetic table model in Kirchhoff stress/original-volume units.
pub const FINITE_TABULATED_J2_LAW_ID: &str = "finite-j2-tabulated-hencky-kirchhoff";

const PARAMETERS: [&str; 6] = [
    "shear_modulus_pa",
    "bulk_modulus_pa",
    "minimum_jacobian",
    "maximum_elastic_log_strain",
    "maximum_plastic_log_strain",
    "maximum_accumulated_plastic_strain",
];

/// Complete native law and the original model/curve selections that define it.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTabulatedFiniteJ2 {
    model: ConstitutiveModelCard,
    model_index: usize,
    context: QueryPoint,
    curve: ResolvedTabulatedHardening,
    law: FiniteJ2,
}

impl ResolvedTabulatedFiniteJ2 {
    /// Select an attached model and one explicit Kirchhoff hardening claim.
    /// The native accumulated-strain domain must fit inside the supplied curve.
    /// No curve fitting, resampling, source query or allocation occurs during
    /// subsequent material-point returns.
    pub fn try_from_material_card_at(
        card: &MaterialCard,
        model_index: usize,
        context: &QueryPoint,
        curve_selection: ClaimSelection,
        maximum_curve_points: usize,
    ) -> Result<Self, MaterialError> {
        let model = card.models().get(model_index).ok_or_else(|| {
            refusal("tabulated finite J2 model index is outside the original material card")
        })?;
        model.validate().map_err(|e| refusal(e.to_string()))?;
        if model.law.0 != FINITE_TABULATED_J2_LAW_ID
            || model.law_version != 1
            || model.state_schema_version != 1
            || model.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(refusal(
                "tabulated finite J2 requires its exact Kirchhoff law/state v1 and virgin native history",
            ));
        }
        if let Some(name) = model
            .parameters
            .keys()
            .find(|name| !PARAMETERS.contains(&name.as_str()))
        {
            return Err(refusal(format!(
                "unrecognized tabulated finite J2 parameter {name}; yield and hardening belong to the selected original curve"
            )));
        }
        if !model
            .validity
            .contains_typed(context.axes(), context.axis_quantities())
        {
            return Err(refusal(
                "tabulated finite J2 context is missing, mistyped or outside the model domain",
            ));
        }
        for axis in model.validity.bounds().keys() {
            if !crate::state_point::is_fixed_external_coordinate(
                axis,
                context.axis_quantities().get(axis).copied(),
                context.axes().get(axis).copied(),
            ) {
                return Err(refusal(format!(
                    "tabulated finite J2 model axis {axis} is not a fixed external condition"
                )));
            }
        }
        let limits = FiniteJ2Limits {
            minimum_jacobian: parameter(model, "minimum_jacobian", Dims::NONE)?,
            maximum_elastic_log_strain: parameter(model, "maximum_elastic_log_strain", Dims::NONE)?,
            maximum_plastic_log_strain: parameter(model, "maximum_plastic_log_strain", Dims::NONE)?,
            maximum_accumulated_plastic_strain: parameter(
                model,
                "maximum_accumulated_plastic_strain",
                Dims::NONE,
            )?,
        };
        let curve = resolve_tabulated_hardening(
            card,
            context,
            &HardeningCurveRequest {
                stress_measure: YieldStressMeasure::Kirchhoff,
                maximum_points: maximum_curve_points,
            },
            curve_selection,
        )
        .map_err(|e| refusal(e.to_string()))?;
        let law = FiniteJ2::new(
            FiniteJ2Parameters {
                shear_modulus_pa: parameter(model, "shear_modulus_pa", Pressure::DIMS)?,
                bulk_modulus_pa: parameter(model, "bulk_modulus_pa", Pressure::DIMS)?,
                yield_stress_pa: curve.hardening().initial_yield_stress_pa(),
                hardening_modulus_pa: 0.0,
            },
            limits,
        )
        .and_then(|law| law.with_tabulated_hardening(curve.hardening().clone()))
        .map_err(|e| refusal(e.to_string()))?;
        Ok(Self {
            model: model.clone(),
            model_index,
            context: context.clone(),
            curve,
            law,
        })
    }

    /// Exact attached model, including source artifacts and original validity.
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.model
    }
    /// Original attached-model index; no same-name model substitution occurs.
    pub const fn model_index(&self) -> usize {
        self.model_index
    }
    /// Original fixed external context with its typed coordinate meanings.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }
    /// Complete selected curve and its original material/claim/knot receipts.
    pub const fn hardening_source(&self) -> &ResolvedTabulatedHardening {
        &self.curve
    }
    /// Original native finite owner, including measured support and strain limits.
    pub const fn law(&self) -> &FiniteJ2 {
        &self.law
    }
}

impl FiniteJ2Material for ResolvedTabulatedFiniteJ2 {
    fn finite_j2(&self) -> &FiniteJ2 {
        &self.law
    }
}

fn refusal(what: impl Into<String>) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

fn parameter(model: &ConstitutiveModelCard, name: &str, dims: Dims) -> Result<f64, MaterialError> {
    let value = model.parameters.get(name).ok_or_else(|| {
        refusal(format!(
            "tabulated finite J2 requires explicit parameter {name}"
        ))
    })?;
    if value.dims != dims || !value.value.is_finite() {
        return Err(refusal(format!(
            "tabulated finite J2 parameter {name} needs finite {dims:?}"
        )));
    }
    Ok(value.value)
}
