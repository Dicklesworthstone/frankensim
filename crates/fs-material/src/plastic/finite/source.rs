//! Explicit source admission for multiplicative quadratic-Hencky J2 plasticity.
//!
//! Law `finite-j2-hencky-kirchhoff`, version 1, binds the exact native
//! [`FiniteJ2`] law. Its yield stress and isotropic hardening modulus act on
//! equivalent **Kirchhoff** stress, not Cauchy stress. Elastic energy,
//! `H alpha^2 / 2` hardening storage and `yield_stress * delta_alpha`
//! dissipation are per original reference volume. These conventions remain
//! significant when `J != 1`; a small-strain offset-yield property is not a
//! calibration of this law unless an author explicitly creates this model card.
//!
//! All eight parameters are mandatory. `shear_modulus_pa`, `bulk_modulus_pa`,
//! `yield_stress_pa` and `hardening_modulus_pa` have pressure dimensions. The
//! dimensionless limits are `minimum_jacobian`, `maximum_elastic_log_strain`,
//! `maximum_plastic_log_strain` and `maximum_accumulated_plastic_strain`.
//! The elastic moduli and initial yield stress are positive; hardening is
//! nonnegative. `0 < minimum_jacobian < 1`, and all three strain limits are
//! positive. The native constructor also checks numerical representability.
//! The two logarithmic limits bound Frobenius tensor norms; the accumulated
//! limit separately bounds cyclic hardening even when the net plastic metric
//! returns toward identity. The native law checks these limits during response.
//! Names, units, versions and initial-state policy are exact: no aliases,
//! unused selectors, inferred moduli or silent constitutive defaults.
//!
//! State schema 1 means the inverse plastic right Cauchy-Green tensor plus
//! accumulated equivalent plastic strain. `ZeroInternalState` means zero
//! plastic logarithmic strain and zero accumulated strain, represented natively
//! by an identity plastic metric and `alpha = 0`, not an all-zero metric.
//! Declared residual history needs a separate explicit initialization owner.
//!
//! The admitted temperature and other external conditions stay fixed for the
//! entire mechanical path. Re-resolving a card at another temperature does not
//! authorize transferring plastic history or omit coefficient/thermal work.
//! Kinematic/Voce hardening, creep, thermal evolution and changing reference
//! configurations are not variants of this version-1 law.
//!
//! The multiplicative/logarithmic return-map family is described by
//! [Simo (1992)](https://doi.org/10.1016/0045-7825(92)90123-2).
//! The [ANSYS finite-strain plasticity reference](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/ans_mat/finitestrainplas.html)
//! likewise distinguishes Kirchhoff yield and multiplicative elastic/plastic
//! deformation. Source admission preserves supplied calibration; it does not
//! establish experimental validation of this particular parameter set.

use super::{FiniteJ2, FiniteJ2Limits, FiniteJ2Material, FiniteJ2Parameters};
use crate::MaterialError;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, QueryPoint};
use fs_qty::Dims;

/// Exact version-1 model identity, including its finite-strain stress convention.
pub const FINITE_J2_LAW_ID: &str = "finite-j2-hencky-kirchhoff";

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const PARAMETERS: [&str; 8] = [
    "shear_modulus_pa",
    "bulk_modulus_pa",
    "yield_stress_pa",
    "hardening_modulus_pa",
    "minimum_jacobian",
    "maximum_elastic_log_strain",
    "maximum_plastic_log_strain",
    "maximum_accumulated_plastic_strain",
];

/// One finite-strain plastic law retaining its exact source and fixed context.
///
/// Spatial consumers borrow this value through [`FiniteJ2Material`]; no source
/// card, native domain or stress/reference convention needs to be discarded.
/// Construction admits source data once. Constitutive trials remain entirely
/// with the native law and do not clone, query or allocate source objects.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedFiniteJ2 {
    card: ConstitutiveModelCard,
    context: QueryPoint,
    law: FiniteJ2,
}

impl ResolvedFiniteJ2 {
    /// Admit a card with unconstrained external conditions.
    ///
    /// A card with temperature or other external bounds needs an explicit
    /// complete context through [`Self::try_from_model_card_at`].
    pub fn try_from_model_card(card: &ConstitutiveModelCard) -> Result<Self, MaterialError> {
        Self::try_from_model_card_at(card, &QueryPoint::new())
    }

    /// Admit the exact typed source domain at fixed external conditions.
    ///
    /// Supported constrained axes are static temperature, pressure,
    /// composition, moisture and concentration. The existing legacy `T`
    /// coordinate retains its kelvin convention. Even a supplied mechanical
    /// strain, stress, stretch, frequency or rate coordinate cannot be frozen
    /// into this source law. Actual evolving strain/history limits belong to
    /// the explicit native parameters above. Unconstrained extra context is
    /// retained verbatim and grants no additional source coverage.
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
                "finite J2 source context is missing, mistyped or outside the model domain",
            ));
        }
        for axis in card.validity.bounds().keys() {
            if !crate::state_point::is_fixed_external_coordinate(
                axis,
                context.axis_quantities().get(axis).copied(),
                context.axes().get(axis).copied(),
            ) {
                return Err(refusal(format!(
                    "finite J2 source axis {axis} is not a supported fixed external condition; evolving coordinates require an explicit owner"
                )));
            }
        }
        Ok(Self {
            card: card.clone(),
            context: context.clone(),
            law,
        })
    }

    /// Original immutable card, including source identity and calibration domain.
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.card
    }

    /// Original complete fixed external query, with all quantity descriptors.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }

    /// Native constitutive owner, including all eight admitted coefficients/limits.
    pub const fn law(&self) -> &FiniteJ2 {
        &self.law
    }
}

impl FiniteJ2Material for ResolvedFiniteJ2 {
    fn finite_j2(&self) -> &FiniteJ2 {
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
        .ok_or_else(|| refusal(format!("{FINITE_J2_LAW_ID}: missing parameter {name}")))?;
    if parameter.dims != dims || !parameter.value.is_finite() {
        return Err(refusal(format!(
            "finite J2 parameter {name}: expected finite {dims:?}, got {parameter:?}"
        )));
    }
    Ok(parameter.value)
}

fn parse_law(card: &ConstitutiveModelCard) -> Result<FiniteJ2, MaterialError> {
    card.validate()
        .map_err(|error| refusal(error.to_string()))?;
    if card.law.0 != FINITE_J2_LAW_ID {
        return Err(refusal(format!(
            "unsupported finite J2 law {}; expected {FINITE_J2_LAW_ID}",
            card.law.0
        )));
    }
    if card.law_version != 1
        || card.state_schema_version != 1
        || card.initial_state != InitialStatePolicy::ZeroInternalState
    {
        return Err(refusal(
            "finite J2 cards require law/state v1 and virgin zero plastic strain/history",
        ));
    }
    if let Some(name) = card
        .parameters
        .keys()
        .find(|name| !PARAMETERS.contains(&name.as_str()))
    {
        return Err(refusal(format!(
            "unrecognized finite J2 card parameter {name}; this law has fixed isotropic linear hardening and fixed external conditions"
        )));
    }
    let parameters = FiniteJ2Parameters {
        shear_modulus_pa: parameter(card, "shear_modulus_pa", PA)?,
        bulk_modulus_pa: parameter(card, "bulk_modulus_pa", PA)?,
        yield_stress_pa: parameter(card, "yield_stress_pa", PA)?,
        hardening_modulus_pa: parameter(card, "hardening_modulus_pa", PA)?,
    };
    let limits = FiniteJ2Limits {
        minimum_jacobian: parameter(card, "minimum_jacobian", Dims::NONE)?,
        maximum_elastic_log_strain: parameter(card, "maximum_elastic_log_strain", Dims::NONE)?,
        maximum_plastic_log_strain: parameter(card, "maximum_plastic_log_strain", Dims::NONE)?,
        maximum_accumulated_plastic_strain: parameter(
            card,
            "maximum_accumulated_plastic_strain",
            Dims::NONE,
        )?,
    };
    FiniteJ2::new(parameters, limits)
        .map_err(|error| refusal(format!("inadmissible finite J2 parameters/domain: {error}")))
}
