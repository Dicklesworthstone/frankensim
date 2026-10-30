//! Source-bound regional Kelvin-Voigt bending, projected on the existing modes.
//! The DKT element owner supplies every strain-energy contribution. No averaged
//! viscosity, damping of conservative prestress, or second modal solve.
//!
//! This is the same diagonal, weak-loss modal approximation used by the uniform
//! plate; its one evaluator is `plate_viscosity::DiagonalPlateViscosity`.
//! Heterogeneous or anisotropic loss is generally nonproportional: off-diagonal
//! modal damping and strongly damped complex modes are not represented by this
//! image (the full tensor image is `plate_viscosity::PlateViscosity`).

use super::*;
use crate::plate_viscosity::{DiagonalPlateViscosity, resolve_region_viscosity};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES;

/// One region's selected Kelvin-Voigt bending loss, in declared region order.
#[derive(Clone, Debug, PartialEq)]
pub enum RegionalBendingViscosity {
    /// Uniaxial eta on an isotropic region: the viscous tensor is (eta/E) times
    /// the region's elastic bending tensor (the proportional legacy law).
    Isotropic(IsotropicPlateBendingViscosity),
    /// Full reduced plane-stress viscosity of an orthotropic region, rotated to
    /// its material axes and integrated over its thickness.
    Orthotropic(OrthotropicPlateBendingViscosity),
}
impl RegionalBendingViscosity {
    /// The proportional isotropic law, or `None` for a full-tensor region.
    #[must_use]
    pub const fn isotropic(&self) -> Option<&IsotropicPlateBendingViscosity> {
        match self {
            Self::Isotropic(law) => Some(law),
            Self::Orthotropic(_) => None,
        }
    }
    /// Applicability of every coefficient of this region [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        match self {
            Self::Isotropic(law) => law.omega_band_rad_s,
            Self::Orthotropic(tensor) => tensor.omega_band_rad_s,
        }
    }
}

/// Rotated DKT bending viscosity of an orthotropic region, resolved from its
/// six explicit `plane_stress_viscosity_*` coefficients by the same owner as the
/// full tensor image. No missing or zero coefficient is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrthotropicPlateBendingViscosity {
    bending_n_m_s: [f64; 9],
    omega_band_rad_s: (f64, f64),
    material_state_identity: ContentHash,
}
impl OrthotropicPlateBendingViscosity {
    /// Bending viscosity tensor [N m s] in mesh coordinates (row-major 3 x 3).
    #[must_use]
    pub const fn bending_matrix_n_m_s(&self) -> &[f64; 9] { &self.bending_n_m_s }
    /// Intersected viscous, elastic and density applicability [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }
    /// Exact material state that produced the tensor; not a validity claim.
    #[must_use]
    pub const fn material_state_identity(&self) -> ContentHash { self.material_state_identity }
}

impl ResolvedPlateChart {
    /// Select Kelvin-Voigt bending in EVERY region, from its retained material
    /// state. An isotropic region needs `kelvin_voigt_bending_viscosity`, with
    /// density, E and nu under the same admission and intersected frequency
    /// policy as a uniform material specimen; its loss tensor is proportional to
    /// its elastic bending tensor. An orthotropic region needs all six
    /// `plane_stress_viscosity_*` coefficients, resolved and rotated exactly as
    /// for the full tensor image. Zero viscosity remains an explicit law with a
    /// checked applicability band. Replacing a region re-resolves this selection;
    /// it cannot silently drop it. Interface/joint damping is not supplied.
    /// Modal realization requires zero authored damping to avoid double counting.
    pub fn with_kelvin_voigt_bending_loss(mut self) -> Result<Self, AcousticRealizeError> {
        let mut viscosities = Vec::with_capacity(self.regions.len());
        for region in &self.regions {
            if region.input.model != PlateMaterialModel::Isotropic {
                let (_, band, bending) = resolve_region_viscosity(region, PLANE_STRESS_VISCOSITY_PROPERTIES)?;
                viscosities.push(RegionalBendingViscosity::Orthotropic(OrthotropicPlateBendingViscosity {
                    bending_n_m_s: bending, omega_band_rad_s: band,
                    material_state_identity: region.input.material.identity(),
                }));
                continue;
            }
            let law = resolve_viscosity(&region.input.material)?;
            let young = plate_property(&region.input.material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
            viscous_time(law, young)?;
            viscosities.push(RegionalBendingViscosity::Isotropic(law));
        }
        self.bending_viscosity = Some(viscosities);
        Ok(self)
    }

    /// Selected loss and source identities in the original region order.
    /// `None` retains the legacy authored-damping path. Geometry/specimen identity
    /// alone is not a complete scenario identity or a claim that loss is selected.
    #[must_use]
    pub fn bending_viscosity(&self) -> Option<&[RegionalBendingViscosity]> {
        self.bending_viscosity.as_deref()
    }

    /// Project selected material loss on an existing mode of this exact chart.
    /// The caller owns the assembly/eigensolve and must supply its matching DOF
    /// map and mode. This keeps valve, body and force ports on the same basis.
    pub(crate) fn modal_bending_damping_ratio(
        &self, model: &PlateModel, phi: &[f64], omega: f64,
    ) -> Result<Option<f64>, AcousticRealizeError> {
        self.bending_viscosity().map(|laws| {
            DiagonalPlateViscosity::new(self, laws, model)?.modal_zeta(model, phi, omega)
        }).transpose()
    }
}

/// Derive body modes, losses and observation ports from one regional specimen.
///
/// With source loss selected, each element contributes `(eta/E) * K_bending`;
/// its contribution is weighted by its ACTUAL modal curvature energy. Neither
/// area nor mass averages decide decay. The denominator includes conservative
/// geometric stiffness, so tension dilutes loss rather than acquiring it.
/// All retained modes must lie inside every region's source-law applicability.
///
/// Single- and multi-port source-bound bodies share the same reduction path.
/// Without selected source loss the original authored chart convention remains.
/// Source bindings remain with the specimen, not the numeric output. This
/// diagonal weak-loss image does not solve full nonproportional damping.
///
/// # Errors
/// Invalid law/controls, nonzero authored damping alongside material loss,
/// source-domain refusal, geometry/modal failure or unrepresentable projection.
pub fn certified_material_chart_radiators(
    specimen: &ResolvedPlateChart, options: &PlateChartRadiation,
) -> Result<Vec<CompactBody>, AcousticRealizeError> {
    ports::reduce(specimen, options, &[], 0).map(|(radiators, _)| radiators)
}

mod ports;

fn refuse(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

// One material-property admission path for uniform and regional plates.
pub(super) fn resolve_viscosity(material: &ResolvedMaterialStatePoint)
    -> Result<IsotropicPlateBendingViscosity, AcousticRealizeError>
{
    let property = material
        .property(KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)
        .ok_or_else(|| refuse("material plate loss needs kelvin_voigt_bending_viscosity"))?;
    let viscosity_pa_s = plate_property(
        material,
        KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
        DynViscosity::DIMS,
    )?;
    let mut band = property
        .angular_frequency_band_rad_s()
        .map_err(|_| refuse("material plate viscosity has invalid frequency-axis semantics or bounds"))?
        .ok_or_else(|| {
            refuse("material plate viscosity needs an explicit omega band in rad/s or a typed frequency band")
        })?;
    for key in [
        DENSITY_PROPERTY,
        YOUNG_MODULUS_PROPERTY,
        POISSON_RATIO_PROPERTY,
        KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
    ] {
        let coefficient = material
            .property(key)
            .ok_or_else(|| refuse("material plate loss needs density, E, nu and viscosity"))?;
        if !coefficient.is_constant_at_fixed_environment() {
            return Err(refuse(
                "Kelvin-Voigt plate coefficients must be constant at the held typed environment, with no frequency-valued or unowned response law",
            ));
        }
        if let Some((lo, hi)) = coefficient.angular_frequency_band_rad_s().map_err(|_| {
            refuse("material plate coefficient has invalid frequency-axis semantics or bounds")
        })? {
            band = (band.0.max(lo), band.1.min(hi));
        }
    }
    Ok(IsotropicPlateBendingViscosity {
        viscosity_pa_s,
        omega_band_rad_s: band,
        material_state_identity: Some(material.identity()),
    })
}

// Same numeric refusal as the original uniform PlateBendingLaw.
pub(crate) fn viscous_time(viscous: IsotropicPlateBendingViscosity, young: f64)
    -> Result<f64, AcousticRealizeError>
{
    let (lo, hi) = viscous.omega_band_rad_s;
    let viscous_time_s = viscous.viscosity_pa_s / young;
    if !viscous.viscosity_pa_s.is_finite()
        || viscous.viscosity_pa_s < 0.0
        || !viscous_time_s.is_finite()
        || (viscous.viscosity_pa_s > 0.0 && viscous_time_s == 0.0)
        || !lo.is_finite()
        || !hi.is_finite()
        || lo < 0.0
        || hi < lo
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "plate viscosity and angular-frequency applicability must be finite and nonnegative with an ordered band",
        });
    }
    Ok(viscous_time_s)
}

#[cfg(test)]
mod tests;
