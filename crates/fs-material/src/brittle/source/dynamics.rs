//! Original density and fracture coefficients for reference-mass dynamics.
//!
//! Density is resolved with the same immutable card, complete physical query
//! and claim-selection policy as the brittle response. It uses the existing
//! scalar fixed-external-condition admission: current strain, stress, rate,
//! frequency or field restrictions cannot be frozen into a mass coefficient.
//!
//! A mechanical owner binds this snapshot to its original reference geometry
//! and retains the resulting inventory across motion, fracture and later
//! thermal coefficient transfers. Re-querying density at a changed temperature
//! does not authorize changing that inventory. No current-to-reference volume
//! conversion, porosity evolution, mass source or momentum transfer is inferred.

use super::{
    BrittleCrackBandSettings, BrittleSourceError, ResolvedBrittleCrackBand,
    brittle_crack_band_requirements, coefficient,
};
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::Density;

/// Original brittle requirements plus positive mass density, in kg/m3.
///
/// Quasistatic callers keep the existing seven-property requirement set.
/// Only callers constructing reference mass need this additional property.
pub fn brittle_dynamic_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = brittle_crack_band_requirements()?;
    requirements.push(ScalarPropertyRequirement::try_new(
        DENSITY_PROPERTY,
        Density::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?);
    Ok(requirements)
}

/// Resolve original fracture response and density as one complete source bundle.
///
/// The receiving dynamics owner must bind the density at its chosen reference
/// state and retain that reference inventory. Resolving a second source state
/// does not transfer existing constitutive history, matter or momentum.
pub fn resolve_brittle_dynamic_material(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
) -> Result<ResolvedBrittleCrackBand, BrittleSourceError> {
    let material =
        resolve_material_state_point(card, point, &brittle_dynamic_requirements()?, selection)?;
    let resolved = ResolvedBrittleCrackBand::try_from_material(material, settings)?;
    resolved.reference_density_kg_m3()?;
    Ok(resolved)
}

impl ResolvedBrittleCrackBand {
    /// Read admitted density from this exact original material snapshot.
    ///
    /// No fallback, material-name inference, external query or detached density
    /// source is used. An adopted multiphysics bundle may have weaker scalar
    /// admissibility than the dynamic resolver, so positivity is checked here
    /// independently. Equal-dimensional mass concentration or tensor/property
    /// conventions do not substitute for the declared bulk-density property.
    ///
    /// The caller binds this value once to the original reference geometry.
    /// It must not refresh an existing mass matrix when constitutive temperature
    /// changes without a separate matter/momentum-transfer law.
    pub fn reference_density_kg_m3(&self) -> Result<f64, BrittleSourceError> {
        let density = coefficient(self.material(), DENSITY_PROPERTY, Density::DIMS)?;
        if density <= 0.0 {
            return Err(BrittleSourceError::Invalid(
                "brittle reference mass requires strictly positive original density",
            ));
        }
        Ok(density)
    }
}
