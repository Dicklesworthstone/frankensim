//! Actual-source material endpoints for conserved liquid mixing and demixing.
//!
//! Input states retain their selected material and pure-Cp sources. The
//! numerical phase owner closes energy and entropy, then the original material
//! resolver refreshes every actual final phase atomically. No averaged phase
//! transport, independent phase solver or implicit isothermal heat bath appears.

use super::super::{
    MulticomponentLiquidCaloricError, MulticomponentLiquidCaloricLaw,
    ResolvedMulticomponentLiquidCaloricState, Result, poll,
};
use super::{
    MulticomponentLiquidEquilibriumCaloricBalance,
    ResolvedMulticomponentLiquidEquilibriumCaloricState,
};
use fs_thermochem::solution::multicomponent::caloric::equilibrium::MulticomponentEquilibriumCaloricSolveV1;
use fs_thermochem::solution::multicomponent::caloric::equilibrium::transactions::MulticomponentCaloricInventoryRefV1;

/// Borrow an original source-bound parcel without collapsing its actual phases.
#[derive(Clone, Copy, Debug)]
pub enum MulticomponentLiquidCaloricInventoryRef<'a> {
    Homogeneous(&'a ResolvedMulticomponentLiquidCaloricState),
    Equilibrium(&'a ResolvedMulticomponentLiquidEquilibriumCaloricState),
}
impl<'a> From<&'a ResolvedMulticomponentLiquidCaloricState>
    for MulticomponentLiquidCaloricInventoryRef<'a>
{
    fn from(state: &'a ResolvedMulticomponentLiquidCaloricState) -> Self {
        Self::Homogeneous(state)
    }
}
impl<'a> From<&'a ResolvedMulticomponentLiquidEquilibriumCaloricState>
    for MulticomponentLiquidCaloricInventoryRef<'a>
{
    fn from(state: &'a ResolvedMulticomponentLiquidEquilibriumCaloricState) -> Self {
        Self::Equilibrium(state)
    }
}
impl<'a> MulticomponentLiquidCaloricInventoryRef<'a> {
    pub(super) fn numeric(self) -> MulticomponentCaloricInventoryRefV1<'a> {
        match self {
            Self::Homogeneous(state) => state.caloric_state().into(),
            Self::Equilibrium(state) => state.caloric_state().into(),
        }
    }
}

impl MulticomponentLiquidCaloricLaw<'_, '_> {
    /// Relax an actual homogeneous source parcel onto the offered local phase
    /// branch at conserved enthalpy plus signed external heat. In particular
    /// Q=0 preserves the thermal effect of demixing and checks total entropy.
    /// The branch supplies only numerical evidence for the same inventory.
    pub fn relax_homogeneous(
        &self,
        initial: &ResolvedMulticomponentLiquidCaloricState,
        branch: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidEquilibriumCaloricBalance> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        self.admit_equilibrium_caloric(branch)?;
        let balance = self.model().relax_homogeneous(
            initial.caloric_state(),
            branch.caloric_state(),
            heat_input_j,
            controls,
            &mut cancelled,
        )?;
        self.resolve_equilibrium_balance(balance, &mut cancelled)
    }

    /// Mix actual homogeneous or separated parcels at possibly unequal source
    /// temperatures. The admitted branch must have their exact combined species
    /// inventory and complete original material/Cp/pressure binding.
    ///
    /// Every trial stays with the original numerical phase continuation. Only
    /// its accepted endpoint queries the original material sources for each
    /// actual phase. Source refusal or cancellation mutates none of the inputs.
    pub fn mix_equilibrium<'a>(
        &self,
        left: impl Into<MulticomponentLiquidCaloricInventoryRef<'a>>,
        right: impl Into<MulticomponentLiquidCaloricInventoryRef<'a>>,
        branch: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidEquilibriumCaloricBalance> {
        poll(&mut cancelled)?;
        let (left, right) = (left.into(), right.into());
        self.admit_phase_inventory(left)?;
        self.admit_phase_inventory(right)?;
        self.admit_equilibrium_caloric(branch)?;
        let balance = self.model().mix_equilibrium(
            left.numeric(),
            right.numeric(),
            branch.caloric_state(),
            heat_input_j,
            controls,
            &mut cancelled,
        )?;
        self.resolve_equilibrium_balance(balance, &mut cancelled)
    }

    pub(super) fn admit_phase_inventory(
        &self,
        input: MulticomponentLiquidCaloricInventoryRef<'_>,
    ) -> Result<()> {
        match input {
            MulticomponentLiquidCaloricInventoryRef::Homogeneous(state) => self.admit(state),
            MulticomponentLiquidCaloricInventoryRef::Equilibrium(state) => self.admit_equilibrium_caloric(state),
        }.map_err(|error| match error {
            MulticomponentLiquidCaloricError::Invalid(_) => MulticomponentLiquidCaloricError::Invalid(
                "phase transaction requires the complete original material and caloric source binding",
            ),
            other => other,
        })
    }
}
