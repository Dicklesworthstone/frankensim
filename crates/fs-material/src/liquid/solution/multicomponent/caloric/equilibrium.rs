//! Finite local-equilibrium heating with the original material and Cp sources.
//!
//! The existing caloric owner advances conserved phase inventories and their
//! sensible-plus-mixing enthalpy. Only its accepted final equilibrium refreshes
//! material properties. Phase count and support stay fixed: continuation is
//! seeded numerical evidence, not global equilibrium or a phase-onset model.

/// Original-source homogeneous relaxation and unequal-temperature phase mixing.
pub mod transactions;

use super::{
    MulticomponentLiquidCaloricError, MulticomponentLiquidCaloricLaw, MulticomponentLiquidError,
    Result, SourceOwner, poll,
};
use crate::liquid::solution::multicomponent::equilibrium::ResolvedMulticomponentLiquidEquilibrium;
use crate::thermal::transport::isobaric::IsobaricHeatCapacityLaw;
use fs_qty::Temperature;
use fs_thermochem::solution::multicomponent::caloric::equilibrium::{
    MulticomponentEquilibriumCaloricBalanceV1, MulticomponentEquilibriumCaloricSolveV1,
    MulticomponentEquilibriumCaloricStateV1,
};
use fs_thermochem::solution::multicomponent::equilibrium::caloric::MulticomponentCaloricResponseControlV1;
use std::sync::Arc;

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Bind an already resolved equilibrium to this SAME complete source owner.
    /// Original phases are retained without another material query or flash.
    /// Pure Cp support and a resolved physical phase response are required.
    pub fn equilibrium_caloric_state(
        &self,
        equilibrium: ResolvedMulticomponentLiquidEquilibrium,
        response: MulticomponentCaloricResponseControlV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentLiquidEquilibriumCaloricState> {
        poll(&mut cancelled)?;
        if equilibrium.phases().is_empty()
            || equilibrium.phases().iter().any(|phase| {
                !Arc::ptr_eq(&self.owner.liquid, &phase.owner) && self.owner.liquid != phase.owner
            })
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "equilibrium calorics require the complete original liquid source owner",
            ));
        }
        let caloric = self.model().equilibrium_caloric_state(
            equilibrium.equilibrium_candidate().clone(),
            response,
            &mut cancelled,
        )?;
        let state = self.retain_equilibrium_state(caloric, equilibrium)?;
        poll(&mut cancelled)?;
        Ok(state)
    }

    /// Add or remove heat while the original finite phase inventories
    /// re-equilibrate at fixed pressure and conserved total species amounts.
    ///
    /// The shared bounded continuation and caloric inverse own every thermal
    /// trial. Only the final accepted candidate refreshes density, viscosity,
    /// conductivity, mass and phase volume from the original material chart.
    /// A disappearing/unresolved branch, source refusal or cancellation leaves
    /// the original state intact and returns no partial thermal endpoint.
    pub fn add_heat_equilibrium(
        &self,
        initial: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidEquilibriumCaloricBalance> {
        poll(&mut cancelled)?;
        self.admit_equilibrium_caloric(initial)?;
        let balance = self.model().add_heat_equilibrium(
            &initial.caloric,
            heat_input_j,
            controls,
            &mut cancelled,
        )?;
        self.resolve_equilibrium_balance(balance, &mut cancelled)
    }

    fn resolve_equilibrium_balance(
        &self,
        balance: MulticomponentEquilibriumCaloricBalanceV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidEquilibriumCaloricBalance> {
        let equilibrium = self
            .liquid
            .resolve_equilibrium(balance.final_state().candidate().clone(), &mut *cancelled)?;
        let final_state =
            self.retain_equilibrium_state(balance.final_state().clone(), equilibrium)?;
        poll(cancelled)?;
        Ok(MulticomponentLiquidEquilibriumCaloricBalance {
            final_state,
            balance,
        })
    }

    /// Stable original-source enthalpy increment at unchanged total inventory.
    pub fn equilibrium_enthalpy_change_j(
        &self,
        initial: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
        final_state: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
    ) -> Result<f64> {
        self.admit_equilibrium_caloric(initial)?;
        self.admit_equilibrium_caloric(final_state)?;
        Ok(self
            .model()
            .equilibrium_enthalpy_change_j(&initial.caloric, &final_state.caloric)?)
    }

    /// Original pure Cp/T and phase-redistribution entropy increment [J/K].
    /// The caller supplies reservoir entropy for externally delivered heat.
    pub fn equilibrium_entropy_change_j_per_k(
        &self,
        initial: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
        final_state: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
    ) -> Result<f64> {
        self.admit_equilibrium_caloric(initial)?;
        self.admit_equilibrium_caloric(final_state)?;
        Ok(self
            .model()
            .equilibrium_entropy_change_j_per_k(&initial.caloric, &final_state.caloric)?)
    }

    fn admit_equilibrium_caloric(
        &self,
        state: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
    ) -> Result<()> {
        if !Arc::ptr_eq(&self.owner, &state.owner) && self.owner != state.owner {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "equilibrium parcel requires its original caloric and material source binding",
            ));
        }
        Ok(())
    }

    fn retain_equilibrium_state(
        &self,
        caloric: MulticomponentEquilibriumCaloricStateV1,
        equilibrium: ResolvedMulticomponentLiquidEquilibrium,
    ) -> Result<ResolvedMulticomponentLiquidEquilibriumCaloricState> {
        let specific_heat_cp_j_kg_k =
            caloric.isobaric_heat_capacity_j_per_k() / equilibrium.total_mass_kg();
        if !specific_heat_cp_j_kg_k.is_finite() || specific_heat_cp_j_kg_k <= 0.0 {
            return Err(MulticomponentLiquidError::Arithmetic(
                "unrepresentable positive equilibrium specific heat",
            )
            .into());
        }
        Ok(ResolvedMulticomponentLiquidEquilibriumCaloricState {
            owner: Arc::clone(&self.owner),
            caloric,
            equilibrium,
            specific_heat_cp_j_kg_k,
        })
    }
}

/// One accepted finite caloric endpoint with each original material phase.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMulticomponentLiquidEquilibriumCaloricState {
    owner: Arc<SourceOwner>,
    caloric: MulticomponentEquilibriumCaloricStateV1,
    equilibrium: ResolvedMulticomponentLiquidEquilibrium,
    specific_heat_cp_j_kg_k: f64,
}
impl ResolvedMulticomponentLiquidEquilibriumCaloricState {
    pub const fn caloric_state(&self) -> &MulticomponentEquilibriumCaloricStateV1 {
        &self.caloric
    }
    pub const fn equilibrium(&self) -> &ResolvedMulticomponentLiquidEquilibrium {
        &self.equilibrium
    }
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw] {
        &self.owner.sources
    }
    pub fn temperature(&self) -> Temperature {
        self.caloric.temperature()
    }
    pub const fn specific_heat_cp_j_kg_k(&self) -> f64 {
        self.specific_heat_cp_j_kg_k
    }
}

/// Closed fixed-pressure heat balance and its accepted material endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidEquilibriumCaloricBalance {
    final_state: ResolvedMulticomponentLiquidEquilibriumCaloricState,
    balance: MulticomponentEquilibriumCaloricBalanceV1,
}
impl MulticomponentLiquidEquilibriumCaloricBalance {
    pub const fn final_state(&self) -> &ResolvedMulticomponentLiquidEquilibriumCaloricState {
        &self.final_state
    }
    pub const fn balance(&self) -> &MulticomponentEquilibriumCaloricBalanceV1 {
        &self.balance
    }
    /// System entropy change; external heat needs its actual reservoir account.
    pub fn entropy_change_j_per_k(&self) -> f64 {
        self.balance.entropy_change_j_per_k()
    }
    pub fn into_final_state(self) -> ResolvedMulticomponentLiquidEquilibriumCaloricState {
        self.final_state
    }
}
