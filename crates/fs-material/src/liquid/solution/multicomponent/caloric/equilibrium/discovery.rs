//! Original material sources at bounded phase-discovered thermal endpoints.
//!
//! Numerical discovery and caloric accounting stay with the original L1 owner.
//! Only its accepted endpoint resolves actual phase density, conductivity,
//! viscosity, mass and volume. No source query uses an averaged phase state.

use super::super::{MulticomponentLiquidCaloricLaw, Result, poll};
use super::ResolvedMulticomponentLiquidEquilibriumCaloricState;
use super::transactions::MulticomponentLiquidCaloricInventoryRef;
use fs_thermochem::solution::multicomponent::caloric::equilibrium::discovery::{
    MulticomponentDiscoveredCaloricBalanceV1, MulticomponentDiscoveredCaloricSolveV1,
};

/// Original resolved phases and their complete bounded-search thermal balance.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidDiscoveredCaloricBalance {
    final_state: ResolvedMulticomponentLiquidEquilibriumCaloricState,
    balance: MulticomponentDiscoveredCaloricBalanceV1,
}
impl MulticomponentLiquidDiscoveredCaloricBalance {
    pub const fn final_state(&self) -> &ResolvedMulticomponentLiquidEquilibriumCaloricState {
        &self.final_state
    }
    pub const fn balance(&self) -> &MulticomponentDiscoveredCaloricBalanceV1 {
        &self.balance
    }
    /// System entropy change; external heat retains its own reservoir account.
    pub fn entropy_change_j_per_k(&self) -> f64 {
        self.balance.balance().entropy_change_j_per_k()
    }
    pub fn into_final_state(self) -> ResolvedMulticomponentLiquidEquilibriumCaloricState {
        self.final_state
    }
}

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Add or remove heat from an original source-bound parcel while the
    /// declared finite search can select a different final liquid phase count.
    /// The original numerical owner supplies every trial and the complete H/S
    /// balance. Its finite-search qualification is retained without upgrading
    /// it to global equilibrium, a latent plateau or a nucleation model.
    pub fn add_heat_discovered<'a>(
        &self,
        initial: impl Into<MulticomponentLiquidCaloricInventoryRef<'a>>,
        heat_input_j: f64,
        controls: MulticomponentDiscoveredCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidDiscoveredCaloricBalance> {
        poll(&mut cancelled)?;
        let initial = initial.into();
        self.admit_phase_inventory(initial)?;
        let balance = self.model().add_heat_discovered(
            initial.numeric(),
            heat_input_j,
            controls,
            &mut cancelled,
        )?;
        self.resolve_discovered_balance(balance, &mut cancelled)
    }

    /// Mix actual homogeneous or separated parcels at possibly unequal source
    /// temperatures without requiring a fixed-phase-count anchor. Every input
    /// must retain the SAME complete material, caloric and pressure binding.
    pub fn mix_discovered<'a>(
        &self,
        left: impl Into<MulticomponentLiquidCaloricInventoryRef<'a>>,
        right: impl Into<MulticomponentLiquidCaloricInventoryRef<'a>>,
        heat_input_j: f64,
        controls: MulticomponentDiscoveredCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidDiscoveredCaloricBalance> {
        poll(&mut cancelled)?;
        let (left, right) = (left.into(), right.into());
        self.admit_phase_inventory(left)?;
        self.admit_phase_inventory(right)?;
        let balance = self.model().mix_discovered(
            left.numeric(),
            right.numeric(),
            heat_input_j,
            controls,
            &mut cancelled,
        )?;
        self.resolve_discovered_balance(balance, &mut cancelled)
    }

    fn resolve_discovered_balance(
        &self,
        balance: MulticomponentDiscoveredCaloricBalanceV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidDiscoveredCaloricBalance> {
        let equilibrium = self
            .liquid
            .resolve_equilibrium(balance.final_state().candidate().clone(), &mut *cancelled)?;
        let final_state =
            self.retain_equilibrium_state(balance.final_state().clone(), equilibrium)?;
        poll(cancelled)?;
        Ok(MulticomponentLiquidDiscoveredCaloricBalance {
            final_state,
            balance,
        })
    }
}
