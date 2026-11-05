//! Original bulk-material data bound to the native brittle crack-band law.
//!
//! These explicitly named fracture coefficients select a constitutive model;
//! ordinary modulus-of-rupture observations, yield stress, hardness, surface
//! adhesion and a material name cannot substitute. The selected source's
//! condition, uncertainty and receipts remain attached. A scalar response is a
//! deterministic estimate, not a weakest-link distribution or a calibrated
//! failure probability. Band width and plane belong to the receiving geometry.
//!
//! Ordinary updates hold external conditions fixed. The explicit
//! [`ResolvedBrittleCrackBand::at_temperature`] and history-transfer methods
//! retain the original source claims, damage and parameter-work accounts.
//! Resolving a second temperature alone never authorizes replacing a law under
//! an accepted history. Rate dependence needs a separate history/work owner.

/// Original-source tensile onset and elastic-to-crack-band handoff.
pub mod initiation;

use super::{
    BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings, BrittleError,
    BrittleSmallStrain,
};
use crate::MaterialError;
use crate::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use crate::elastic::IsotropicElastic;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure, QuantitySpec};

/// Pure opening crack-band peak traction, Pa; not flexural modulus of rupture.
pub const BRITTLE_TENSILE_STRENGTH_PROPERTY: &str = "brittle_tensile_strength";
/// Pure shear crack-band peak traction, Pa; independent of the opening strength.
pub const BRITTLE_SHEAR_STRENGTH_PROPERTY: &str = "brittle_shear_strength";
/// Complete pure-opening fracture work per crack area, J/m²; not toughness K_Ic.
pub const BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY: &str = "brittle_mode_i_fracture_energy";
/// Complete pure-shear fracture work per crack area, J/m².
pub const BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY: &str = "brittle_mode_ii_fracture_energy";
/// Dimensionless Benzeggagh--Kenane exponent for the selected onset-mix law.
pub const BRITTLE_BK_EXPONENT_PROPERTY: &str = "brittle_bk_exponent";

/// Exact source requirements; no absent shear or mixed-mode data are inferred.
pub fn brittle_crack_band_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let positive = ScalarAdmissibility::StrictlyPositive;
    Ok(vec![
        ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, positive)?,
        ScalarPropertyRequirement::try_new(
            POISSON_RATIO_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::OpenInterval {
                lower: -1.0,
                upper: 0.5,
            },
        )?,
        ScalarPropertyRequirement::try_new(
            BRITTLE_TENSILE_STRENGTH_PROPERTY,
            Pressure::DIMS,
            positive,
        )?,
        ScalarPropertyRequirement::try_new(
            BRITTLE_SHEAR_STRENGTH_PROPERTY,
            Pressure::DIMS,
            positive,
        )?,
        ScalarPropertyRequirement::try_new(
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            positive,
        )?,
        ScalarPropertyRequirement::try_new(
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            positive,
        )?,
        ScalarPropertyRequirement::try_new(BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, positive)?,
    ])
}

/// Immutable source bundle and the original prepared native law.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBrittleCrackBand {
    material: ResolvedMaterialStatePoint,
    pub(super) law: BrittleCrackBand,
}

impl ResolvedBrittleCrackBand {
    /// Adopt an existing multiphysics bundle, retaining every additional receipt.
    ///
    /// The complete coefficient set must have scalar SI semantics and supported
    /// fixed external conditions. Mechanical/rate/frequency coordinates refuse,
    /// including restrictions on otherwise constant coefficients. The native
    /// owner admits the geometric width against both fracture-energy envelopes.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        settings: BrittleCrackBandSettings,
    ) -> Result<Self, BrittleSourceError> {
        let youngs = coefficient(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
        let poisson = coefficient(&material, POISSON_RATIO_PROPERTY, Dims::NONE)?;
        let parameters = BrittleCrackBandParameters {
            tensile_strength_pa: coefficient(
                &material,
                BRITTLE_TENSILE_STRENGTH_PROPERTY,
                Pressure::DIMS,
            )?,
            shear_strength_pa: coefficient(
                &material,
                BRITTLE_SHEAR_STRENGTH_PROPERTY,
                Pressure::DIMS,
            )?,
            mode_i_fracture_energy_j_m2: coefficient(
                &material,
                BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
            )?,
            mode_ii_fracture_energy_j_m2: coefficient(
                &material,
                BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
                FRACTURE_ENERGY_PER_AREA_DIMS,
            )?,
            bk_exponent: coefficient(&material, BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE)?,
        };
        let elastic = IsotropicElastic::new(youngs, poisson, settings.maximum_strain)?;
        let law = BrittleCrackBand::new(elastic, parameters, settings)?;
        Ok(Self { material, law })
    }

    /// Complete original material identity, fixed query, claims and receipts.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Checked constitutive owner. A history borrows this particular law.
    #[must_use]
    pub const fn law(&self) -> &BrittleCrackBand {
        &self.law
    }

    /// Borrow the native law for existing uniform/per-element solid solvers.
    /// The source bundle cannot be dropped while its histories are in use.
    #[must_use]
    pub fn small_strain_law(&self) -> BrittleSmallStrain<'_> {
        self.law.small_strain_law()
    }
}

/// Resolve the complete fracture response at one supported material state.
pub fn resolve_brittle_crack_band(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
) -> Result<ResolvedBrittleCrackBand, BrittleSourceError> {
    let material =
        resolve_material_state_point(card, point, &brittle_crack_band_requirements()?, selection)?;
    ResolvedBrittleCrackBand::try_from_material(material, settings)
}

/// Original source, elastic or crack-band admission failure.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleSourceError {
    /// The shared material source owner refused the selected claim/context.
    MaterialState(MaterialStatePointError),
    /// The original isotropic elastic owner refused its input.
    Elastic(MaterialError),
    /// The original crack-band/cohesive owner refused its input.
    Native(BrittleError),
    /// A source is not an admissible fixed scalar for this constitutive role.
    Invalid(&'static str),
}
impl core::fmt::Display for BrittleSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MaterialState(e) => write!(f, "brittle material source: {e}"),
            Self::Elastic(e) => write!(f, "brittle elastic backbone: {e}"),
            Self::Native(e) => write!(f, "brittle crack band: {e}"),
            Self::Invalid(e) => write!(f, "brittle source: {e}"),
        }
    }
}
impl std::error::Error for BrittleSourceError {}
impl From<MaterialStatePointError> for BrittleSourceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::MaterialState(error)
    }
}
impl From<MaterialError> for BrittleSourceError {
    fn from(error: MaterialError) -> Self {
        Self::Elastic(error)
    }
}
impl From<BrittleError> for BrittleSourceError {
    fn from(error: BrittleError) -> Self {
        Self::Native(error)
    }
}

fn coefficient(
    material: &ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<f64, BrittleSourceError> {
    let property = material.property(name).ok_or(BrittleSourceError::Invalid(
        "missing or ambiguous brittle constitutive coefficient",
    ))?;
    let requirement = property.requirement();
    if requirement.quantity() != QuantitySpec::dimensional(dims)
        || requirement.hardness_test().is_some()
        || requirement.elastic_component().is_some()
        || requirement.strain_component().is_some()
        || requirement.stress_component().is_some()
        || !property.value_si().is_finite()
        || !property.is_static_at_fixed_external_conditions()
    {
        return Err(BrittleSourceError::Invalid(
            "brittle coefficients require scalar SI data at fixed external conditions",
        ));
    }
    Ok(property.value_si())
}
