//! Conserved mixing and demixing into an explicitly admitted local phase branch.
//!
//! The branch anchor is numerical phase evidence for the combined inventory,
//! not an additional physical parcel. Its original predictor/corrector and
//! caloric inverse determine the endpoint. Actual input phases supply all
//! initial energy and entropy, including unequal source temperatures. This
//! permits adiabatic demixing without imposing an isothermal heat sink.

use super::super::{
    MixingInventory, MulticomponentCaloricErrorV1, MulticomponentCaloricStateV1,
    MulticomponentLiquidCaloricModelV1, Sum, Temperature, checked,
};
use super::{
    MulticomponentEquilibriumCaloricBalanceV1, MulticomponentEquilibriumCaloricErrorV1,
    MulticomponentEquilibriumCaloricSolveV1, MulticomponentEquilibriumCaloricStateV1, Result, poll,
};

/// Borrow an actual homogeneous parcel or an already separated phase inventory.
/// Every input retains its original caloric owner and common parcel temperature.
#[derive(Clone, Copy, Debug)]
pub enum MulticomponentCaloricInventoryRefV1<'a> {
    Homogeneous(&'a MulticomponentCaloricStateV1),
    Equilibrium(&'a MulticomponentEquilibriumCaloricStateV1),
}
impl<'a> From<&'a MulticomponentCaloricStateV1> for MulticomponentCaloricInventoryRefV1<'a> {
    fn from(state: &'a MulticomponentCaloricStateV1) -> Self {
        Self::Homogeneous(state)
    }
}
impl<'a> From<&'a MulticomponentEquilibriumCaloricStateV1>
    for MulticomponentCaloricInventoryRefV1<'a>
{
    fn from(state: &'a MulticomponentEquilibriumCaloricStateV1) -> Self {
        Self::Equilibrium(state)
    }
}
impl<'a> MulticomponentCaloricInventoryRefV1<'a> {
    pub fn temperature(self) -> Temperature {
        match self {
            Self::Homogeneous(state) => state.temperature(),
            Self::Equilibrium(state) => state.temperature(),
        }
    }
    pub fn component_moles(self) -> &'a [f64] {
        match self {
            Self::Homogeneous(state) => state.component_moles(),
            Self::Equilibrium(state) => state.component_moles(),
        }
    }
    pub(super) fn phases(self) -> impl Iterator<Item = MixingInventory<'a>> {
        let count = match self {
            Self::Homogeneous(_) => 1,
            Self::Equilibrium(state) => state.candidate().phases().len(),
        };
        (0..count).map(move |phase| match self {
            Self::Homogeneous(state) => MixingInventory::from(state),
            Self::Equilibrium(state) => {
                let phase = &state.candidate().phases()[phase];
                MixingInventory {
                    mixing: phase.thermodynamics(),
                    component_moles: phase.component_moles(),
                    total_moles: phase.total_moles(),
                }
            }
        })
    }
}

impl MulticomponentLiquidCaloricModelV1 {
    /// Relax an actual homogeneous parcel into the explicitly offered local
    /// phase branch with signed external heat. At Q=0, demixing storage changes
    /// temperature while conserving enthalpy; the total entropy may not decrease
    /// beyond the declared allowance. No kinetics or global minimum is inferred.
    ///
    /// `branch` must have exactly the initial conserved inventory and the SAME
    /// complete source owner. It supplies a continuation anchor, not input heat
    /// or matter. Its phase count and support remain fixed during the solve.
    pub fn relax_homogeneous(
        &self,
        initial: &MulticomponentCaloricStateV1,
        branch: &MulticomponentEquilibriumCaloricStateV1,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricBalanceV1> {
        self.equilibrate_inputs(&[initial.into()], branch, heat_input_j, controls, cancelled)
    }

    /// Mix homogeneous or separated parcels at possibly unequal temperatures
    /// into one admitted local phase branch. All initial and final phase
    /// enthalpy and entropy contribute once. The branch's inventory must equal
    /// the representable species-by-species sum of both actual input parcels.
    ///
    /// Every thermal trial uses the existing bounded phase continuation and
    /// exact coupled Cp. A phase loss, unresolved response, source limit, work
    /// exhaustion or cancellation publishes no endpoint and mutates no input.
    pub fn mix_equilibrium<'a>(
        &self,
        left: impl Into<MulticomponentCaloricInventoryRefV1<'a>>,
        right: impl Into<MulticomponentCaloricInventoryRefV1<'a>>,
        branch: &MulticomponentEquilibriumCaloricStateV1,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricBalanceV1> {
        self.equilibrate_inputs(
            &[left.into(), right.into()],
            branch,
            heat_input_j,
            controls,
            cancelled,
        )
    }

    fn equilibrate_inputs(
        &self,
        inputs: &[MulticomponentCaloricInventoryRefV1<'_>],
        branch: &MulticomponentEquilibriumCaloricStateV1,
        heat_input_j: f64,
        controls: MulticomponentEquilibriumCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentEquilibriumCaloricBalanceV1> {
        poll(&mut cancelled)?;
        self.admit_equilibrium_solve(heat_input_j, controls)?;
        self.admit_equilibrium_caloric(branch)?;
        for &input in inputs {
            match input {
                MulticomponentCaloricInventoryRefV1::Homogeneous(state) => self.admit(state)?,
                MulticomponentCaloricInventoryRefV1::Equilibrium(state) => {
                    self.admit_equilibrium_caloric(state)?;
                }
            }
        }
        let amounts = self.combine_component_moles(inputs.iter().map(|x| x.component_moles()))?;
        let n = self.heat_capacities().len();
        if branch.component_moles() != &amounts[..n] {
            return Err(MulticomponentEquilibriumCaloricErrorV1::Invalid(
                "phase anchor must have exactly the combined physical input inventory",
            ));
        }
        let feed = branch.candidate().feed().mole_fractions();
        // Difference through the SAME conserved homogeneous feed. The common
        // feed cancels analytically; no arbitrary pure enthalpy/entropy datum
        // is subtracted, and phase weights remain their actual finite amounts.
        let input_to_feed = self.mixing_changes(
            inputs.iter().flat_map(|input| input.phases()),
            feed,
            &mut cancelled,
        )?;
        let admitted = self.equilibrium_caloric_state(
            branch.candidate().clone(),
            controls.response,
            &mut cancelled,
        )?;
        let (root, completed_phase_solves) =
            self.solve_equilibrium_heat(admitted, controls, &mut cancelled, |state, cancelled| {
                let output_to_feed = self.mixing_changes(
                    MulticomponentCaloricInventoryRefV1::from(state).phases(),
                    feed,
                    cancelled,
                )?;
                let mut change = Sum::default();
                change.add(
                    self.sensible_inventory_change(
                        inputs
                            .iter()
                            .map(|x| (x.temperature(), x.component_moles())),
                        state.temperature(),
                        false,
                        cancelled,
                    )?,
                );
                change.add(input_to_feed[0]);
                change.add(-output_to_feed[0]);
                let change_j = checked(change.value(), "phase transaction enthalpy change")?;
                // Q remains separate from phase storage until the compensated
                // residual is formed. It is never folded into an artificial
                // large branch heat target that can erase a real small Q.
                change.add(-heat_input_j);
                Ok((
                    change_j,
                    checked(change.value(), "phase transaction heat residual")?,
                ))
            })?;
        let final_state = &root.sample.state;
        let output_to_feed = self.mixing_changes(
            MulticomponentCaloricInventoryRefV1::from(final_state).phases(),
            feed,
            &mut cancelled,
        )?;
        let mut entropy = Sum::default();
        entropy.add(
            self.sensible_inventory_change(
                inputs
                    .iter()
                    .map(|x| (x.temperature(), x.component_moles())),
                final_state.temperature(),
                true,
                &mut cancelled,
            )?,
        );
        entropy.add(input_to_feed[1]);
        entropy.add(-output_to_feed[1]);
        let entropy_change_j_per_k = checked(entropy.value(), "phase transaction entropy change")?;
        if heat_input_j == 0.0
            && entropy_change_j_per_k < -controls.thermal.entropy_tolerance_j_per_k
        {
            return Err(MulticomponentCaloricErrorV1::EntropyDecrease {
                change_j_per_k: entropy_change_j_per_k,
            }
            .into());
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentEquilibriumCaloricBalanceV1 {
            final_state: root.sample.state,
            heat_input_j,
            enthalpy_change_j: root.sample.change_j,
            entropy_change_j_per_k,
            residual_j: root.sample.residual_j,
            temperature_bracket: root.bracket.map(Temperature::new),
            iterations: root.iterations,
            completed_phase_solves,
        })
    }
}
