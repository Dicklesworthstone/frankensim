//! Local phase response from the same original source-bound caloric owner.
//!
//! Pure-liquid Cp(T) supplies the frozen caloric response; the thermochemical
//! owner supplies the additional response from phase inventory redistribution.
//! Original material phases and their transport claims remain attached. This
//! differential holds pressure and total species amounts fixed and retains
//! the accepted candidate's numerical and local-branch qualifications.

use super::{MulticomponentLiquidCaloricLaw, MulticomponentLiquidError, Result, poll};
use crate::liquid::solution::multicomponent::equilibrium::ResolvedMulticomponentLiquidEquilibrium;
use crate::thermal::transport::isobaric::IsobaricHeatCapacityLaw;
use fs_thermochem::solution::multicomponent::equilibrium::caloric::{
    MulticomponentCaloricResponseControlV1, MulticomponentCaloricResponseV1,
};
use std::sync::Arc;

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Original mixture coordinate explicitly qualified as the shared isobar.
    pub fn pressure_axis(&self) -> &str {
        &self.owner.pressure_axis
    }

    /// Differentiate the already accepted smooth local phase branch using
    /// this owner's original pure-liquid Cp curves and retained source receipts.
    ///
    /// The complete material and thermodynamic owner must match, including
    /// component molar masses, original source receipts and coordinate roles.
    /// Singular phase response, exhausted workspace, or cancellation publishes
    /// no result. No new phase search, heat inversion, global-minimum claim, or
    /// effective transport heat capacity is introduced by this differential.
    pub fn isobaric_caloric_response<'response>(
        &'response self,
        equilibrium: &'response ResolvedMulticomponentLiquidEquilibrium,
        controls: MulticomponentCaloricResponseControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentLiquidCaloricResponse<'response>> {
        poll(&mut cancelled)?;
        if equilibrium.phases().is_empty()
            || equilibrium.phases().iter().any(|phase| {
                !Arc::ptr_eq(&self.owner.liquid, &phase.owner) && self.owner.liquid != phase.owner
            })
        {
            return Err(super::MulticomponentLiquidCaloricError::Invalid(
                "resolved equilibrium belongs to a different complete liquid source owner",
            ));
        }
        let caloric = equilibrium
            .equilibrium_candidate()
            .isobaric_caloric_response(self.model().heat_capacities(), controls, &mut cancelled)?;
        let specific_heat_cp_j_kg_k =
            caloric.isobaric_heat_capacity_j_per_k() / equilibrium.total_mass_kg();
        if !specific_heat_cp_j_kg_k.is_finite() || specific_heat_cp_j_kg_k <= 0.0 {
            return Err(MulticomponentLiquidError::Arithmetic(
                "unrepresentable positive mass-specific equilibrium heat capacity",
            )
            .into());
        }
        poll(&mut cancelled)?;
        Ok(ResolvedMulticomponentLiquidCaloricResponse {
            equilibrium,
            caloric,
            sources: self.sources(),
            pressure_axis: self.pressure_axis(),
            specific_heat_cp_j_kg_k,
        })
    }
}

/// Qualified local Cp retaining the original accepted phases and pure sources.
/// This deliberately supplies no homogeneous or effective transport trait.
pub struct ResolvedMulticomponentLiquidCaloricResponse<'a> {
    equilibrium: &'a ResolvedMulticomponentLiquidEquilibrium,
    caloric: MulticomponentCaloricResponseV1<'a>,
    sources: &'a [IsobaricHeatCapacityLaw],
    pressure_axis: &'a str,
    specific_heat_cp_j_kg_k: f64,
}
impl<'a> ResolvedMulticomponentLiquidCaloricResponse<'a> {
    pub const fn equilibrium(&self) -> &'a ResolvedMulticomponentLiquidEquilibrium {
        self.equilibrium
    }
    /// Total, frozen, redistribution, and species-conserving derivative reports.
    pub const fn caloric(&self) -> &MulticomponentCaloricResponseV1<'a> {
        &self.caloric
    }
    pub const fn sources(&self) -> &'a [IsobaricHeatCapacityLaw] {
        self.sources
    }
    pub const fn pressure_axis(&self) -> &'a str {
        self.pressure_axis
    }
    /// Total equilibrium Cp divided by the original finite component mass.
    pub const fn specific_heat_cp_j_kg_k(&self) -> f64 {
        self.specific_heat_cp_j_kg_k
    }
}
