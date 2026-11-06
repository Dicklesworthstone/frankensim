//! Source-defined finite porous compaction through the original native law.
//!
//! `finite-porous-hencky-kirchhoff`, law/state version 1, explicitly selects
//! multiplicative pressure-sensitive plasticity with quadratic Hencky storage
//! and linear energetic hardening. Yield and hardening are Kirchhoff stresses;
//! energy and dissipation refer to the original volume. A small-strain foam
//! curve or Cauchy yield value cannot silently select this finite law.
//!
//! All nine coefficients and limits are mandatory. The dimensionless pressure
//! sensitivity is the normalized Deshpande--Fleck ellipse parameter, not
//! porosity. Zero sensitivity explicitly selects its J2 limit; zero hardening
//! explicitly selects perfect plasticity. The retained inverse plastic metric
//! may change volume, and alpha measures accumulated work-conjugate flow.
//! `ZeroInternalState` means the identity metric and zero accumulation.
//!
//! Resolution is cold: native Newton trials borrow the original law and
//! accepted history without querying or cloning source data. External source
//! conditions stay fixed throughout a mechanical path. Re-resolution does not
//! authorize transferring history, changing mass, or discarding parameter work.
//! Evolving porosity, densification, thermal response and measured validation
//! remain outside the original law's claim.

use super::{FinitePorous, FinitePorousLimits, FinitePorousMaterial, FinitePorousParameters};
use crate::MaterialError;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, QueryPoint};
use fs_qty::Dims;

/// Exact model selector, including finite kinematics and stress convention.
pub const FINITE_POROUS_LAW_ID: &str = "finite-porous-hencky-kirchhoff";

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const PARAMETERS: [&str; 9] = [
    "shear_modulus_pa",
    "bulk_modulus_pa",
    "yield_stress_pa",
    "hardening_modulus_pa",
    "pressure_sensitivity",
    "minimum_jacobian",
    "maximum_elastic_log_strain",
    "maximum_plastic_log_strain",
    "maximum_accumulated_plastic_strain",
];

/// Immutable source and fixed context for the same finite compaction owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedFinitePorous {
    card: ConstitutiveModelCard,
    context: QueryPoint,
    law: FinitePorous,
}

impl ResolvedFinitePorous {
    /// Resolve a card with no constrained external coordinates.
    pub fn try_from_model_card(card: &ConstitutiveModelCard) -> Result<Self, MaterialError> {
        Self::try_from_model_card_at(card, &QueryPoint::new())
    }

    /// Resolve the exact typed domain at fixed external conditions.
    ///
    /// Static temperature, external pressure, composition, moisture and
    /// concentration use the shared material-coordinate admission. Evolving
    /// strain, stress, stretch, frequency and rate cannot be frozen as source
    /// conditions. Native parameters separately bound the actual history.
    pub fn try_from_model_card_at(
        card: &ConstitutiveModelCard,
        context: &QueryPoint,
    ) -> Result<Self, MaterialError> {
        let law = parse_law(card)?;
        if !card
            .validity
            .contains_typed(context.axes(), context.axis_quantities())
        {
            return Err(refusal(
                "finite porous source context is missing, mistyped or outside the model domain",
            ));
        }
        for axis in card.validity.bounds().keys() {
            if !crate::state_point::is_fixed_external_coordinate(
                axis,
                context.axis_quantities().get(axis).copied(),
                context.axes().get(axis).copied(),
            ) {
                return Err(refusal(format!(
                    "finite porous source axis {axis} is not a fixed external condition; evolving coordinates require their own owner"
                )));
            }
        }
        Ok(Self {
            card: card.clone(),
            context: context.clone(),
            law,
        })
    }

    /// Complete original card, including calibration identity and validity.
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.card
    }

    /// Exact external context and all its quantity descriptors.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }

    /// Borrow the native compaction, history, tangent and energy owner.
    /// Its `frozen` method supplies the existing finite-element potential API.
    pub const fn law(&self) -> &FinitePorous {
        &self.law
    }
}

impl FinitePorousMaterial for ResolvedFinitePorous {
    fn finite_porous(&self) -> &FinitePorous {
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
            "finite porous parameter {name}: expected finite {dims:?}, got {value:?}"
        )));
    }
    Ok(value.value)
}

fn parse_law(card: &ConstitutiveModelCard) -> Result<FinitePorous, MaterialError> {
    card.validate()
        .map_err(|error| refusal(error.to_string()))?;
    if card.law.0 != FINITE_POROUS_LAW_ID
        || card.law_version != 1
        || card.state_schema_version != 1
        || card.initial_state != InitialStatePolicy::ZeroInternalState
    {
        return Err(refusal(format!(
            "finite porous cards require {FINITE_POROUS_LAW_ID}, law/state v1 and virgin plastic history"
        )));
    }
    if let Some(name) = card
        .parameters
        .keys()
        .find(|name| !PARAMETERS.contains(&name.as_str()))
    {
        return Err(refusal(format!(
            "unrecognized finite porous parameter {name}; coefficients and native limits must be explicit"
        )));
    }
    FinitePorous::new(
        FinitePorousParameters {
            shear_modulus_pa: parameter(card, "shear_modulus_pa", PA)?,
            bulk_modulus_pa: parameter(card, "bulk_modulus_pa", PA)?,
            yield_stress_pa: parameter(card, "yield_stress_pa", PA)?,
            hardening_modulus_pa: parameter(card, "hardening_modulus_pa", PA)?,
            pressure_sensitivity: parameter(card, "pressure_sensitivity", Dims::NONE)?,
        },
        FinitePorousLimits {
            minimum_jacobian: parameter(card, "minimum_jacobian", Dims::NONE)?,
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
            "inadmissible finite porous coefficients/domain: {error}"
        ))
    })
}
