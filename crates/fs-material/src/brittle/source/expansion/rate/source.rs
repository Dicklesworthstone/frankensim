//! Explicit source selection for thermal strain and physical fracture viscosity.
//!
//! Intrinsic fracture energy, plane viscosity and instantaneous expansion are
//! distinct measured quantities. They resolve from one immutable material card
//! and retain their own selected receipts through the original thermal refresh.
//! No numerical relaxation time, bulk damping coefficient or mean expansion
//! coefficient substitutes for those properties.

use super::super::{BrittleSourceError, ResolvedThermalBrittle, resolve_with_reference};
use crate::brittle::BrittleCrackBandSettings;
use crate::brittle::source::ResolvedBrittleCrackBand;
use crate::brittle::source::rate::{
    brittle_viscous_dynamic_requirements, brittle_viscous_requirements,
};
use crate::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement,
};
use fs_matdb::{MaterialCard, QueryPoint};

/// Original intrinsic fracture data, both physical plane viscosities and
/// instantaneous alpha(T). Zero or negative expansion remains admissible.
pub fn thermal_viscous_brittle_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    with_expansion(brittle_viscous_requirements()?)
}

/// The same complete source bundle with positive original density for an
/// inertial consumer that explicitly binds its reference mass once.
pub fn thermal_viscous_brittle_dynamic_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    with_expansion(brittle_viscous_dynamic_requirements()?)
}

fn with_expansion(
    mut requirements: Vec<ScalarPropertyRequirement>,
) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    requirements.push(ScalarPropertyRequirement::try_new(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        ScalarAdmissibility::Finite,
    )?);
    Ok(requirements)
}

/// Resolve explicit fracture-plane viscosity with the original integrated
/// thermal strain and immutable expansion reference.
///
/// Every selected claim remains attached to the same complete material bundle.
/// Only canonical temperature T may change between queries; all held context,
/// coordinate conventions, alpha-path support and original source checks apply.
/// The physical jump-rate ceiling is an explicit model bound in m/s. Intrinsic
/// fracture energy excludes the additional viscous loss. Mechanical evolution
/// requires the timed adapter; the untimed native law continues to refuse.
pub fn resolve_thermal_viscous_brittle(
    card: &MaterialCard,
    reference_query: &QueryPoint,
    current_query: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
    maximum_rate_m_s: f64,
) -> Result<ResolvedThermalBrittle, BrittleSourceError> {
    resolve_with_reference(
        card,
        reference_query,
        current_query,
        selection,
        &thermal_viscous_brittle_requirements()?,
        |material| {
            ResolvedBrittleCrackBand::try_from_viscous_material(
                material,
                settings,
                maximum_rate_m_s,
            )
        },
    )
}

/// Resolve the same thermal/fracture-rate law with original reference density.
///
/// Density is the CURRENT material snapshot. The receiving motion owner binds
/// an initial snapshot to reference geometry and retains that mass inventory;
/// later density refreshes do not authorize changing mass or momentum.
pub fn resolve_thermal_viscous_brittle_dynamic_material(
    card: &MaterialCard,
    reference_query: &QueryPoint,
    current_query: &QueryPoint,
    selection: MaterialPropertySelection,
    settings: BrittleCrackBandSettings,
    maximum_rate_m_s: f64,
) -> Result<ResolvedThermalBrittle, BrittleSourceError> {
    let resolved = resolve_with_reference(
        card,
        reference_query,
        current_query,
        selection,
        &thermal_viscous_brittle_dynamic_requirements()?,
        |material| {
            ResolvedBrittleCrackBand::try_from_viscous_material(
                material,
                settings,
                maximum_rate_m_s,
            )
        },
    )?;
    resolved.band().reference_density_kg_m3()?;
    Ok(resolved)
}

impl ResolvedThermalBrittle {
    /// Convert an ACTUAL instantaneous temperature rate [K/s] to signed free
    /// linear-strain rate [1/s] using this source's current pinned alpha(T).
    ///
    /// This local derivative neither supplies a temperature trajectory nor
    /// equates endpoint rate with an interval secant. For varying alpha,
    /// alpha(T_end) * (T_end-T_start)/dt generally differs from
    /// (free_end-free_start)/dt. The caller declares which physical path/rate
    /// its force quadrature and accepted endpoint observation represent.
    /// Zero and negative alpha or temperature rates retain their physical sign.
    pub fn free_strain_rate_from_temperature_rate(
        &self,
        temperature_rate_k_s: f64,
    ) -> Result<f64, BrittleSourceError> {
        if !temperature_rate_k_s.is_finite() {
            return Err(BrittleSourceError::Invalid(
                "thermal free-strain rate requires an actual finite temperature rate",
            ));
        }
        let alpha = self.expansion().current().linear_coefficient_per_k();
        let rate = alpha * temperature_rate_k_s;
        if !rate.is_finite() || (alpha != 0.0 && temperature_rate_k_s != 0.0 && rate == 0.0) {
            return Err(BrittleSourceError::Invalid(
                "source-derived instantaneous free-strain rate is unrepresentable",
            ));
        }
        Ok(rate)
    }
}
