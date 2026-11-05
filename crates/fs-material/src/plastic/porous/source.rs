//! Fixed-condition source coefficients for the native porous plasticity law.
//!
//! All five coefficients are resolved at the same physical point. The
//! hardening datum is the slope against this law's work-conjugate plastic
//! strain, and the sensitivity is its normalized Deshpande--Fleck ellipse
//! parameter. Density, porosity, a friction angle or another law's monotonic
//! tangent cannot supply either coefficient by inference.
//!
//! The complete material bundle retains its original claims, units, query,
//! validity and uncertainty. Conditions remain fixed for a load path;
//! re-resolution does not transfer history or account for thermal/parameter
//! work. The native small-strain, associated-flow and symmetric-yield domain
//! remains binding. No measured foam calibration or densification is implied.

use core::fmt;

use fs_evidence::ModelCard;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure, QuantitySpec};

use super::{PorousPlasticity, PorousResponse, PorousState};
use crate::elastic::IsotropicElastic;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, ResolvedScalarProperty, ScalarAdmissibility,
    ScalarPropertyRequirement, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};

/// Energetic hardening H [Pa] in Y(alpha) = yield_stress + H alpha, where
/// alpha is the native associated flow's accumulated work-conjugate strain.
/// A J2 modulus or total-strain tangent is not selected as a fallback.
pub const POROUS_HARDENING_MODULUS_PROPERTY: &str = "porous_hardening_modulus";
/// Dimensionless a in phi² = (q² + a² p²)/(1 + a²/9).
/// This is neither a friction angle, a porosity nor a hydrostatic yield stress.
pub const POROUS_PRESSURE_SENSITIVITY_PROPERTY: &str = "porous_pressure_sensitivity";

/// Exact five source quantities. Perfect plasticity and the J2 limit require
/// explicit zero H and zero a respectively; absent properties never mean zero.
pub fn porous_requirements() -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    Ok(vec![
        ScalarPropertyRequirement::try_new(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            POISSON_RATIO_PROPERTY,
            Dims([0; 6]),
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        )?,
        ScalarPropertyRequirement::try_new(
            YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            POROUS_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )?,
        ScalarPropertyRequirement::try_new(
            POROUS_PRESSURE_SENSITIVITY_PROPERTY,
            Dims([0; 6]),
            ScalarAdmissibility::NonNegative,
        )?,
    ])
}

/// Native pressure-sensitive response with its complete immutable source bundle.
/// Implements the existing small-strain interface without another return map.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedPorousPlasticity {
    material: ResolvedMaterialStatePoint,
    law: PorousPlasticity,
}

impl ResolvedPorousPlasticity {
    /// Adopt a shared bundle without dropping unrelated material receipts.
    /// Curves/surfaces at supported fixed external conditions are evaluated by
    /// the common resolver. Evolving strain, stress, frequency, rate and history
    /// conditions cannot be frozen into these coefficients.
    ///
    /// # Errors
    /// Missing/wrong quantities, unsupported coefficient conditions and native
    /// parameter refusals. The explicit strain-norm bound must be in (0, 0.05].
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        strain_limit: f64,
    ) -> Result<Self, PorousSourceError> {
        if !(strain_limit.is_finite() && strain_limit > 0.0 && strain_limit <= 0.05) {
            return Err(PorousSourceError::Invalid(
                "small-strain bound must be in (0, 0.05]",
            ));
        }
        let youngs = coefficient(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let poisson = coefficient(&material, POISSON_RATIO_PROPERTY, Dims([0; 6]))?.value_si();
        let yield_stress =
            coefficient(&material, YIELD_STRESS_PROPERTY, Pressure::DIMS)?.value_si();
        let hardening =
            coefficient(&material, POROUS_HARDENING_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let sensitivity = coefficient(
            &material,
            POROUS_PRESSURE_SENSITIVITY_PROPERTY,
            Dims([0; 6]),
        )?
        .value_si();
        let elastic = IsotropicElastic::new(youngs, poisson, strain_limit)?;
        let law = PorousPlasticity::new(elastic, yield_stress, hardening, sensitivity)?;
        Ok(Self { material, law })
    }

    /// Original quantities, source decisions, query, domains and uncertainty.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Original native return, tangent, stored energy and increment work owner.
    pub const fn law(&self) -> &PorousPlasticity {
        &self.law
    }

    /// Typed native trial; the supplied accepted history is never mutated.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &PorousState,
    ) -> Result<PorousResponse, MaterialError> {
        self.law.response(strain, state)
    }

    /// Explicit hardening calibration and its retained usage receipt.
    pub fn hardening_source(&self) -> &ResolvedScalarProperty {
        self.material
            .property(POROUS_HARDENING_MODULUS_PROPERTY)
            .expect("immutable porous bundle retains its admitted hardening source")
    }

    /// Explicit ellipse parameter and its retained usage receipt.
    pub fn pressure_sensitivity_source(&self) -> &ResolvedScalarProperty {
        self.material
            .property(POROUS_PRESSURE_SENSITIVITY_PROPERTY)
            .expect("immutable porous bundle retains its admitted sensitivity source")
    }
}

impl SmallStrainLaw for ResolvedPorousPlasticity {
    type State = PorousState;
    fn initial_state(&self) -> PorousState {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &PorousState) -> Voigt {
        self.law.stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &PorousState) -> Tangent6 {
        self.law.tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &PorousState) -> PorousState {
        self.law.update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

/// Resolve every native coefficient at one explicit material state. No name,
/// missing calibration, density or another law's tangent supplies defaults.
pub fn resolve_porous_plasticity(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    strain_limit: f64,
) -> Result<ResolvedPorousPlasticity, PorousSourceError> {
    let material = resolve_material_state_point(card, point, &porous_requirements()?, selection)?;
    ResolvedPorousPlasticity::try_from_material(material, strain_limit)
}

/// Original source/parameter refusal, or an unsupported coefficient reduction.
#[derive(Clone, Debug, PartialEq)]
pub enum PorousSourceError {
    MaterialState(MaterialStatePointError),
    Parameters(MaterialError),
    Invalid(&'static str),
}
impl fmt::Display for PorousSourceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaterialState(source) => write!(f, "porous plasticity source refusal: {source}"),
            Self::Parameters(source) => write!(f, "porous plasticity parameter refusal: {source}"),
            Self::Invalid(what) => write!(f, "porous plasticity reduction refusal: {what}"),
        }
    }
}
impl std::error::Error for PorousSourceError {}
impl From<MaterialStatePointError> for PorousSourceError {
    fn from(source: MaterialStatePointError) -> Self {
        Self::MaterialState(source)
    }
}
impl From<MaterialError> for PorousSourceError {
    fn from(source: MaterialError) -> Self {
        Self::Parameters(source)
    }
}

fn coefficient<'a>(
    material: &'a ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<&'a ResolvedScalarProperty, PorousSourceError> {
    let property = material.property(name).ok_or(PorousSourceError::Invalid(
        "missing or ambiguous porous coefficient",
    ))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.value_si().is_finite()
        || !property.has_only_fixed_external_conditions(None, false)
        || property.angular_frequency_band_rad_s()?.is_some()
    {
        return Err(PorousSourceError::Invalid(
            "porous plasticity needs absolute SI coefficients at fixed external conditions, without evolving mechanical, history, rate or frequency coordinates",
        ));
    }
    Ok(property)
}
