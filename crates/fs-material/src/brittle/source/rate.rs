//! Source-bound viscosity of the surviving fracture-plane bonds.
//!
//! These coefficients are traction per separation rate, Pa s/m. They describe
//! the selected fracture process, not a bulk Kelvin modulus or an arbitrary
//! numerical damage-relaxation time. The geometry supplies jump and volumetric
//! work scaling; it does not silently recalibrate these source coefficients.
//! The integrity evolution still depends on the declared band, so arbitrary-h
//! objectivity of total viscous fracture work is not established. Intrinsic Gc
//! must exclude this additional loss.
//!
//! All coefficients and optional reference density resolve from the SAME card,
//! held external conditions and claim selection. Current strain, separation,
//! damage or rate-conditioned coefficients cannot be frozen into this model.
//! The caller supplies an explicit admissible jump-rate bound. This is a
//! constitutive selection, not evidence of calibration or a failure probability.

use super::{
    BrittleCrackBandSettings, BrittleSourceError, ResolvedBrittleCrackBand,
    brittle_crack_band_requirements, coefficient,
};
use crate::brittle::BrittleError;
use crate::cohesive::CohesiveViscosity;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};

/// Viscous normal traction per fracture-plane opening rate, Pa s/m.
pub const BRITTLE_NORMAL_VISCOSITY_PROPERTY: &str = "brittle_normal_viscosity";
/// Viscous shear traction per fracture-plane sliding rate, Pa s/m.
pub const BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY: &str = "brittle_tangential_viscosity";
/// Interface/fracture-plane viscosity is distinct from bulk dynamic viscosity.
pub use crate::cohesive::source::COHESIVE_VISCOSITY_DIMS as BRITTLE_VISCOSITY_DIMS;

/// Original intrinsic brittle data and both explicit positive viscosities.
pub fn brittle_viscous_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = brittle_crack_band_requirements()?;
    append_viscosities(&mut requirements)?;
    Ok(requirements)
}

/// Original fracture, viscosity and density receipts for reference-mass dynamics.
pub fn brittle_viscous_dynamic_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = super::dynamics::brittle_dynamic_requirements()?;
    append_viscosities(&mut requirements)?;
    Ok(requirements)
}

fn append_viscosities(
    requirements: &mut Vec<ScalarPropertyRequirement>,
) -> Result<(), MaterialStatePointError> {
    for name in [
        BRITTLE_NORMAL_VISCOSITY_PROPERTY,
        BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
    ] {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            BRITTLE_VISCOSITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    Ok(())
}

impl ResolvedBrittleCrackBand {
    /// Adopt the complete original source bundle with physical bond viscosity.
    ///
    /// Both viscosity properties must already be present, scalar and admitted
    /// at the same fixed external conditions. Extra multiphysics receipts are
    /// retained unchanged. No division by the numerical band width, inference
    /// from bulk damping, or fitted total-work substitution for intrinsic Gc is
    /// performed. Mechanical updates of this law require an explicit duration.
    pub fn try_from_viscous_material(
        material: ResolvedMaterialStatePoint,
        settings: BrittleCrackBandSettings,
        maximum_rate_m_s: f64,
    ) -> Result<Self, BrittleSourceError> {
        let normal = coefficient(
            &material,
            BRITTLE_NORMAL_VISCOSITY_PROPERTY,
            BRITTLE_VISCOSITY_DIMS,
        )?;
        let tangential = coefficient(
            &material,
            BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
            BRITTLE_VISCOSITY_DIMS,
        )?;
        let viscosity = CohesiveViscosity::new(normal, tangential, maximum_rate_m_s)
            .map_err(BrittleError::from)?;
        let mut source = Self::try_from_material(material, settings)?;
        source.law = source.law.with_viscosity(viscosity);
        Ok(source)
    }
}

/// Resolve the complete intrinsic fracture and physical rate response once.
pub fn resolve_viscous_brittle_crack_band(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
    maximum_rate_m_s: f64,
) -> Result<ResolvedBrittleCrackBand, BrittleSourceError> {
    let material =
        resolve_material_state_point(card, point, &brittle_viscous_requirements()?, selection)?;
    ResolvedBrittleCrackBand::try_from_viscous_material(material, settings, maximum_rate_m_s)
}

/// Resolve the same response together with original reference density.
/// The mechanical owner retains that mass across damage and later thermal
/// coefficient changes; querying new density never authorizes changing mass.
pub fn resolve_viscous_brittle_dynamic_material(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
    maximum_rate_m_s: f64,
) -> Result<ResolvedBrittleCrackBand, BrittleSourceError> {
    let material = resolve_material_state_point(
        card,
        point,
        &brittle_viscous_dynamic_requirements()?,
        selection,
    )?;
    let resolved =
        ResolvedBrittleCrackBand::try_from_viscous_material(material, settings, maximum_rate_m_s)?;
    resolved.reference_density_kg_m3()?;
    Ok(resolved)
}
