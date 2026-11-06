//! Complete prescribed-exterior durations through the original interval owner.

use super::*;
use super::super::stabilized::{
    StabilizedReactiveEntropyIntervalReport, StabilizedReactiveIntervalControl,
    StabilizedReactiveIntervalReport,
};

/// A completely admitted open reacting duration or represented clock target.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedReactiveExteriorEntropyIntervalReport {
    /// Native interval durations, actual source/boundary balances, and work
    /// consumed by both molecular halves and every discarded coupled attempt.
    pub coupled: StabilizedReactiveIntervalReport,
    /// Direct original-to-final entropy with actual chemical production and
    /// all accepted molecular and advective exterior transfers. One allowance
    /// applies to this whole window; component allowances are not added.
    pub entropy: ReactiveExteriorEntropyReport,
    /// Smallest admitted complete-step generation [J/K], zero for a no-op.
    /// This is diagnostic and does not replace the direct whole-window gate.
    pub minimum_step_entropy_generation_j_k: f64,
    /// Largest relative species continuity residual over accepted whole steps.
    pub maximum_step_relative_species_residual: f64,
}

impl StabilizedReactiveExteriorEntropyIntervalReport {
    fn from_closed(closed: StabilizedReactiveEntropyIntervalReport) -> Self {
        Self {
            coupled: closed.coupled,
            entropy: ReactiveExteriorEntropyReport {
                source_entropy_change_j_k: closed.entropy.source_entropy_change_j_k,
                boundary_entropy_input_j_k: [0.0; 2],
                entropy_generation_j_k: closed.entropy.entropy_change_j_k,
                source_entropy_generation_j_k: closed.entropy.source_entropy_change_j_k,
                reference_entropy_residual_j_k: closed.entropy.reference_entropy_residual_j_k,
                maximum_relative_species_residual: closed.entropy.maximum_relative_species_residual,
                rejected_entropy_steps: closed.entropy.rejected_entropy_steps,
            },
            minimum_step_entropy_generation_j_k: closed.minimum_step_entropy_change_j_k,
            maximum_step_relative_species_residual: closed.maximum_step_relative_species_residual,
        }
    }
}

impl ReactingEuler1d {
    /// Complete a nonnegative reacting duration with source-conjugate RKL2,
    /// original chemistry, and entropy-consistent prescribed exterior transport.
    ///
    /// Every native D-C-A-C-D step uses the same stationary end states for both
    /// molecular halves and the middle inviscid update. They remain fixed
    /// throughout the complete requested window, including rejected attempts.
    /// Each Some supplies the corresponding reflecting or transmissive end;
    /// a None physical end must be a stationary reflecting wall. Two absent
    /// states execute the original closed entropy interval exactly.
    ///
    /// One immutable source-entropy basis is retained from the ORIGINAL state.
    /// Only completely accepted steps contribute their actual chemical species
    /// changes, molecular/advective end transfers and absolute stage throughput.
    /// The complete final source-state entropy is evaluated directly against
    /// the original conserved rollback rows and that same boundary account.
    /// Genuine chemical entropy and species-carried entropy remain included;
    /// only the separately bounded numerical species-source reference defect
    /// is removed. Per-step tolerances and entropy reports are never summed
    /// into the complete-window admission.
    ///
    /// The original interval owner retains its exact duration/target semantics,
    /// one reused conserved backup, compensated physical ledgers, and shared
    /// coupled-attempt and molecular-RHS budgets. Both D halves and abandoned
    /// retries consume those same budgets. Every terminal remainder is solved.
    /// A late source, entropy, work, clock or cancellation refusal restores the
    /// ENTIRE original field and clock, including any accepted internal prefix.
    ///
    /// First-order reconstruction, fixed constant-area geometry, stationary
    /// source-compatible baths and common positive species support (or common
    /// absence) remain required. Bath cells use the original canonical species
    /// order and energy references; these prescribed traces are not depleted
    /// or chemically advanced. Moving baths and finite exterior storage are
    /// outside this entry. Native chemistry and macrostep controls remain
    /// necessary; split symmetry and entropy admission establish neither a
    /// positivity theorem nor a complete temporal-error estimate.
    ///
    /// A zero duration still validates the original molecular operator, bath
    /// states, source and entropy scope. An unrepresentable duration refuses;
    /// use the target-time companion to request an explicit represented clock.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stabilized_entropy_consistent_with_exterior_states_for<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        duration: Time,
        controls: StabilizedReactiveIntervalControl,
        entropy: MixtureEntropyControl,
        exterior_states: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveExteriorEntropyIntervalReport, ReactiveError> {
        if exterior_states.iter().all(Option::is_none) {
            return self.advance_stabilized_entropy_consistent_for(
                operator, duration, controls, entropy, cancelled,
            ).map(StabilizedReactiveExteriorEntropyIntervalReport::from_closed);
        }
        poll(&mut cancelled)?;
        controls.validate()?;
        let mut account = IntervalEntropy::new_with_exterior_states(
            &self.accepted, controls.step.transport, entropy, exterior_states, &mut cancelled,
        )?;
        let coupled = self.advance_stabilized_for_inner(
            operator, duration, controls, exterior_states, Some(&mut account), &mut cancelled,
        )?;
        Ok(account.complete_exterior(coupled))
    }

    /// Reach an explicit represented clock target through the same complete
    /// open reacting interval and original immutable entropy/source account.
    ///
    /// The fixed stationary exterior states and all native component checks
    /// match the duration entry. The original interval loop solves the final
    /// remainder without clock snapping or interpolation. Shared work budgets
    /// include abandoned steps; final admission uses one whole-window entropy
    /// allowance. Failure after any internal accepted prefix restores the
    /// original complete field and clock.
    ///
    /// A target equal to the current time still admits operator, source, bath
    /// and entropy scope. A past or nonfinite target refuses. Two absent
    /// exterior states preserve the original closed target entry exactly.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stabilized_entropy_consistent_with_exterior_states_to<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        target: Time,
        controls: StabilizedReactiveIntervalControl,
        entropy: MixtureEntropyControl,
        exterior_states: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveExteriorEntropyIntervalReport, ReactiveError> {
        if exterior_states.iter().all(Option::is_none) {
            return self.advance_stabilized_entropy_consistent_to(
                operator, target, controls, entropy, cancelled,
            ).map(StabilizedReactiveExteriorEntropyIntervalReport::from_closed);
        }
        poll(&mut cancelled)?;
        controls.validate()?;
        let mut account = IntervalEntropy::new_with_exterior_states(
            &self.accepted, controls.step.transport, entropy, exterior_states, &mut cancelled,
        )?;
        let coupled = self.advance_stabilized_to_inner(
            operator, target, controls, exterior_states, Some(&mut account), &mut cancelled,
        )?;
        Ok(account.complete_exterior(coupled))
    }
}
