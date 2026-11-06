//! Open D-C-A-C-D composition through the original molecular and reactor owners.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

/// One complete open reacting step with source-conjugate RKL2 molecular halves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StabilizedReactiveExteriorEntropyStepReport {
    /// Complete source-state/boundary entropy and original conservative reports.
    pub coupled: ReactiveExteriorEntropyStepReport,
    /// Actual accepted molecular half intervals and their original RKL2 work.
    pub dissipation: SplitStabilizedDissipationReport,
    /// All spatial RHS evaluations, including discarded whole-step attempts.
    pub operator_evaluations: usize,
}

impl ReactingEuler1d {
    /// Couple source-conjugate molecular transport, original chemistry and
    /// entropy-consistent advection through stationary prescribed exterior gas.
    ///
    /// The same exterior species, temperature and referenced energy feed both
    /// molecular halves and both inviscid SSP stages. Molecular transport uses
    /// the original Maxwell--Stefan, heat-conduction and viscous energy fluxes,
    /// including the original constitutive cross effects when supplied by the
    /// transport law. No independent reaction or viscous heat is added.
    ///
    /// Every molecular half retains its native state-minus-bath entropy gate.
    /// Actual accepted boundary transfers and absolute stage throughput are
    /// then contracted in the COMPLETE step's immutable pre-D entropy basis.
    /// This keeps genuine chemical entropy and checks each species against
    /// both reaction changes and all molecular/advective boundary exchange.
    ///
    /// The original D(dt/2)-C(dt/2)-A(dt)-C(dt/2)-D(dt/2) restart owner remains
    /// authoritative. If A shortens the duration, both D and C halves restart
    /// from the accepted state. Abandoned transfers and entropy receipts are
    /// cleared; final conservative reporting and cancellation precede publication.
    ///
    /// A supplied exterior must have zero momentum under the existing molecular
    /// bath model. An absent end must reflect; two absent states delegate to the
    /// closed entropy method exactly. First-order space, fixed constant area,
    /// original source compatibility and common species support are required.
    /// Exterior storage is prescribed, not finite. This finite-step admission
    /// does not establish global temporal order, unconditional positivity or
    /// combustion validity.
    #[allow(clippy::too_many_arguments)]
    pub fn step_entropy_consistent_with_dissipation_and_exterior_states<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        entropy: MixtureEntropyControl,
        exterior_states: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveExteriorEntropyStepReport, ReactiveError> {
        if exterior_states.iter().all(Option::is_none) {
            let closed = self.step_entropy_consistent_with_dissipation(
                operator, maximum_duration, controls, dissipation_controls,
                entropy, cancelled,
            )?;
            return Ok(ReactiveExteriorEntropyStepReport::from_closed(closed));
        }
        poll(&mut cancelled)?;
        controls.validate()?;
        dissipation_controls.validate()?;
        let mut account = CoupledEntropy::new_with_exterior_states(
            &self.accepted, controls.transport, entropy, exterior_states, &mut cancelled,
        )?;
        let (coupled, _) = self.step_inner(
            maximum_duration,
            controls,
            Some(operator),
            stabilized::MolecularMethod::EntropyExplicit {
                controls: dissipation_controls, entropy,
            },
            exterior_states,
            None,
            Some(&mut account),
            &mut cancelled,
        )?;
        Ok(account.complete_exterior(coupled))
    }

    /// Use the original source-conjugate RKL2 molecular solver for both halves
    /// of the same open reacting composition and complete entropy admission.
    ///
    /// This has the same stationary exterior, source and geometry requirements
    /// as `step_entropy_consistent_with_dissipation_and_exterior_states`.
    /// One caller-bounded RHS counter spans both molecular halves and every
    /// discarded coupling attempt. Exhaustion is terminal; a shorter retry
    /// cannot renew the work allowance.
    ///
    /// Signed native RKL quadrature supplies the actual conservative exchange;
    /// absolute throughput retains every contributing stage, including negative
    /// anchor weights. Internal Legendre states are not treated as separate
    /// physical entropy-producing substeps. Macrostep, half-interval and complete
    /// composition entropy checks remain distinct and keep their own allowances.
    #[allow(clippy::too_many_arguments)]
    pub fn step_entropy_consistent_with_stabilized_dissipation_and_exterior_states<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        entropy: MixtureEntropyControl,
        exterior_states: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveExteriorEntropyStepReport, ReactiveError> {
        if exterior_states.iter().all(Option::is_none) {
            let closed = self.step_entropy_consistent_with_stabilized_dissipation(
                operator, maximum_duration, controls, dissipation_controls,
                entropy, cancelled,
            )?;
            return Ok(StabilizedReactiveExteriorEntropyStepReport {
                coupled: ReactiveExteriorEntropyStepReport::from_closed(closed.coupled),
                dissipation: closed.dissipation,
                operator_evaluations: closed.operator_evaluations,
            });
        }
        poll(&mut cancelled)?;
        controls.validate()?;
        dissipation_controls.validate()?;
        let mut account = CoupledEntropy::new_with_exterior_states(
            &self.accepted, controls.transport, entropy, exterior_states, &mut cancelled,
        )?;
        let mut evaluations = 0;
        let step = self.stabilized_step_inner(
            operator, maximum_duration, controls, dissipation_controls,
            exterior_states, &mut evaluations, None, Some(&mut account), &mut cancelled,
        )?;
        Ok(StabilizedReactiveExteriorEntropyStepReport {
            coupled: account.complete_exterior(step.coupled),
            dissipation: step.dissipation,
            operator_evaluations: step.operator_evaluations,
        })
    }
}
