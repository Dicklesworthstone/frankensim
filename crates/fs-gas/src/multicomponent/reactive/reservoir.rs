//! Fixed external baths through the original reacting-flow composition.
//!
//! The same borrowed state reaches every advective and molecular stage and
//! every complete-step retry. Chemistry consumes the evolving interior state;
//! the prescribed bath itself is neither reacted nor depleted. Species order
//! and energy gauges follow the original `Cell` reservoir contract.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

pub(super) fn admit(
    flow: &MixtureEuler1d,
    reservoirs: [Option<Cell<'_>>; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ReactiveError> {
    if reservoirs.iter().any(Option::is_some) {
        // Use the same stationary bath/source admission as molecular transport.
        // Check before paying for chemistry or mutating an interval candidate.
        super::super::dissipation::reservoir::admit(
            flow.gas(),
            flow.boundaries(),
            reservoirs,
            cancelled,
        )?;
    }
    Ok(())
}

impl ReactingEuler1d {
    /// Advance finite-rate chemistry and advection against stationary baths.
    ///
    /// Each supplied bath overrides that physical end for the whole accepted
    /// C-A-C step. Its composition, pressure and enthalpy enter the original
    /// characteristic inflow/outflow closure; the source state must use this
    /// gas's canonical species order and per-species energy gauges. Periodic
    /// geometry and moving baths refuse before chemistry starts.
    ///
    /// Acoustic or heat-release retries replay both chemistry halves from the
    /// original field with the same baths. Only a complete step publishes its
    /// state, time and source/boundary balances. Two absent baths preserve
    /// [`Self::step`] exactly.
    pub fn step_with_reservoirs(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        reservoirs: [Option<Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
        if reservoirs.iter().all(Option::is_none) {
            return self.step(maximum_duration, controls, cancelled);
        }
        self.step_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            stabilized::MolecularMethod::Explicit(DissipationControl::default()),
            reservoirs,
            None,
            cancelled,
        )
        .map(|(report, _)| report)
    }

    /// Couple the same baths to advection and both SSP RK2 molecular halves.
    ///
    /// The original D-C-A-C-D owner integrates each molecular half over the
    /// actual accepted duration. Half-cell gradients, bath species enthalpy
    /// and stationary traction work come from the shared reservoir face law.
    /// Advective and molecular boundary transfers enter separate original
    /// accounts and are included once in the complete species/energy balance.
    /// Baths remain prescribed and undepleted throughout the call.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_and_reservoirs<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        reservoirs: [Option<Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveStepReport, ReactiveError> {
        if reservoirs.iter().all(Option::is_none) {
            return self.step_with_dissipation(
                operator,
                maximum_duration,
                controls,
                dissipation_controls,
                cancelled,
            );
        }
        self.step_impl(
            maximum_duration,
            controls,
            Some(operator),
            stabilized::MolecularMethod::Explicit(dissipation_controls),
            reservoirs,
            None,
            cancelled,
        )
        .map(|(report, _)| report)
    }

    /// Use RKL2 for both molecular halves with the same fixed exterior baths.
    ///
    /// The RHS-work counter includes both halves and every discarded coupled
    /// attempt. A shortened acoustic step cannot reset that budget. Chemistry,
    /// boundary/source accounting and final caloric-domain admission retain
    /// their original owners; symmetric splitting supplies no global temporal
    /// error estimate or nonlinear RKL2 stability theorem.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_stabilized_dissipation_and_reservoirs<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        reservoirs: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveStepReport, ReactiveError> {
        if reservoirs.iter().all(Option::is_none) {
            return self.step_with_stabilized_dissipation(
                operator,
                maximum_duration,
                controls,
                dissipation_controls,
                cancelled,
            );
        }
        poll(&mut cancelled)?;
        dissipation_controls.validate()?;
        self.stabilized_step(
            operator,
            maximum_duration,
            controls,
            dissipation_controls,
            reservoirs,
            &mut 0,
            None,
            &mut cancelled,
        )
    }

    /// Complete a duration with fixed baths through the existing RKL2 window.
    ///
    /// The original clock-span admission, complete-window rollback and shared
    /// coupling/RHS budgets apply. Actual advective and molecular bath transfers
    /// accumulate with compensation across every accepted physical step.
    /// Change a prescribed bath between calls to supply a boundary history;
    /// use the finite network/plenum owner when bath depletion is required.
    pub fn advance_stabilized_with_reservoirs_for<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        duration: Time,
        controls: StabilizedReactiveIntervalControl,
        reservoirs: [Option<Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        if reservoirs.iter().all(Option::is_none) {
            return self.advance_stabilized_for(operator, duration, controls, cancelled);
        }
        self.advance_stabilized_for_impl(operator, duration, controls, reservoirs, cancelled)
    }

    /// Reach a represented target with prescribed baths and atomic rollback.
    ///
    /// Every terminal remainder is actually integrated. A zero-length window
    /// still admits the baths and the original molecular operator. Invalid bath
    /// states, late chemistry/source refusals, exhausted work and cancellation
    /// leave the original accepted field and clock intact.
    pub fn advance_stabilized_with_reservoirs_to<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        target: Time,
        controls: StabilizedReactiveIntervalControl,
        reservoirs: [Option<Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        if reservoirs.iter().all(Option::is_none) {
            return self.advance_stabilized_to(operator, target, controls, cancelled);
        }
        self.advance_stabilized_to_impl(operator, target, controls, reservoirs, cancelled)
    }
}
