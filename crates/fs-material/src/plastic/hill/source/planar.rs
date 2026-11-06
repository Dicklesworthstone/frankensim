//! Source-retaining planar calibration from explicitly selected Lankford data.
//!
//! The original r-value and transverse-shear claims remain in the material
//! bundle. Derived Hill coefficients are never substituted for those receipts.

use super::{
    HILL_HARDENING_MODULUS_PROPERTY, HILL_REFERENCE_YIELD_STRESS_PROPERTY, HillSourceError,
    ResolvedHill48, coefficient,
};
use crate::elastic::IsotropicElastic;
use crate::plastic::hill::oriented::OrientedHill48;
use crate::plastic::hill::{Hill48, Hill48Coefficients};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, POISSON_RATIO_PROPERTY,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Dims, Pressure};

/// Plastic width/thickness ratios at 0, 45 and 90 degrees to rolling.
pub const HILL_PLANAR_R_PROPERTIES: [&str; 3] = [
    "hill_lankford_r_0",
    "hill_lankford_r_45",
    "hill_lankford_r_90",
];
/// Independent shear yields in planes 23 and 31, divided by `Y_rolling/sqrt(3)`.
pub const HILL_TRANSVERSE_SHEAR_RATIO_PROPERTIES: [&str; 2] =
    ["hill_shear_yield_ratio_23", "hill_shear_yield_ratio_31"];

/// Nine original scalar requirements, without inferred transverse shear data.
/// `hill_reference_yield_stress` is explicitly the rolling-direction yield;
/// the hardening modulus is against the resulting Hill accumulated flow.
pub fn planar_r_requirements() -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = vec![
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
            HILL_REFERENCE_YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            HILL_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            ScalarAdmissibility::NonNegative,
        )?,
    ];
    for property in HILL_PLANAR_R_PROPERTIES
        .into_iter()
        .chain(HILL_TRANSVERSE_SHEAR_RATIO_PROPERTIES)
    {
        requirements.push(ScalarPropertyRequirement::try_new(
            property,
            Dims([0; 6]),
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    Ok(requirements)
}

impl ResolvedHill48 {
    /// Retain original planar r-values and both independent shear-ratio sources.
    /// Columns of `material_to_spatial` are rolling, cross and thickness axes.
    /// The reference yield is rolling-direction yield, fixing `G+H=1`.
    /// Fixed source conditions and the native small-strain limits remain binding.
    pub fn try_from_planar_r_material(
        material: ResolvedMaterialStatePoint,
        strain_limit: f64,
        material_to_spatial: [[f64; 3]; 3],
    ) -> Result<Self, HillSourceError> {
        if !(strain_limit.is_finite() && strain_limit > 0.0 && strain_limit <= 0.05) {
            return Err(HillSourceError::Invalid(
                "small-strain bound must be in (0, 0.05]",
            ));
        }
        let youngs = coefficient(&material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let poisson = coefficient(&material, POISSON_RATIO_PROPERTY, Dims([0; 6]))?.value_si();
        let reference_yield = coefficient(
            &material,
            HILL_REFERENCE_YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
        )?
        .value_si();
        let hardening =
            coefficient(&material, HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS)?.value_si();
        let mut planar_r = [0.0; 3];
        for (value, property) in planar_r.iter_mut().zip(HILL_PLANAR_R_PROPERTIES) {
            *value = coefficient(&material, property, Dims([0; 6]))?.value_si();
        }
        let mut shear_ratios = [0.0; 2];
        for (value, property) in shear_ratios
            .iter_mut()
            .zip(HILL_TRANSVERSE_SHEAR_RATIO_PROPERTIES)
        {
            *value = coefficient(&material, property, Dims([0; 6]))?.value_si();
        }
        let coefficients = Hill48Coefficients::from_planar_r_values(planar_r, shear_ratios)?;
        let local = Hill48::new(
            IsotropicElastic::new(youngs, poisson, strain_limit)?,
            reference_yield,
            hardening,
            coefficients,
        )?;
        Ok(Self {
            material,
            law: OrientedHill48::new(local, material_to_spatial)?,
        })
    }
}

/// Resolve original planar strain-ratio data without supplying missing yields.
pub fn resolve_hill48_from_planar_r_values(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    strain_limit: f64,
    material_to_spatial: [[f64; 3]; 3],
) -> Result<ResolvedHill48, HillSourceError> {
    let material = resolve_material_state_point(card, point, &planar_r_requirements()?, selection)?;
    ResolvedHill48::try_from_planar_r_material(material, strain_limit, material_to_spatial)
}
