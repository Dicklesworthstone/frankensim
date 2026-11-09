//! Held mechanical sources through the original stabilized reacting owner.

use super::*;
use crate::multicomponent::dissipation::StabilizedDissipationControl;

impl ReactingEuler1d {
    /// Couple held body forcing to chemistry, RKL2 diffusion and optional baths.
    ///
    /// Both molecular halves and every discarded complete attempt consume the
    /// original shared RHS budget. The original acoustic/source restarts keep
    /// D/2-C/2-F/2-A-F/2-C/2-D/2 on one actual accepted duration. Acceleration
    /// uses each force stage's real density, and final caloric/Mach admission
    /// follows the last molecular half. Prescribed baths are not accelerated.
    /// RKL2 extends its original parabolic stability interval; it supplies no
    /// unconditional nonlinear positivity or global temporal-error guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_stabilized_dissipation_and_body_force<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: StabilizedDissipationControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        reservoirs: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveStepReport, ReactiveError> {
        poll(&mut cancelled)?;
        dissipation_controls.validate()?;
        self.stabilized_step_inner_forced(
            operator,
            maximum_duration,
            controls,
            dissipation_controls,
            reservoirs,
            &mut 0,
            None,
            None,
            Some((force, force_controls)),
            &mut cancelled,
        )
    }

    /// Integrate a complete held-source duration with atomic whole-window rollback.
    ///
    /// Mechanical impulse, work and representation error accumulate from the
    /// accepted native source maps with compensation. The original-to-final
    /// species, momentum and absolute-energy account includes mechanical work
    /// once alongside distinct chemical and boundary transfers. Every actual
    /// accepted step and discarded attempt shares the original work budgets.
    /// Force values and optional prescribed baths stay fixed throughout this
    /// call; use adjacent calls at actual changes in the authored source history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stabilized_with_body_force_for<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        duration: Time,
        controls: StabilizedReactiveIntervalControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        reservoirs: [Option<Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        self.advance_stabilized_for_inner_forced(
            operator,
            duration,
            controls,
            reservoirs,
            None,
            Some((force, force_controls)),
            cancelled,
        )
    }

    /// Reach an explicit clock target through the same forced reacting window.
    ///
    /// The terminal remainder is actually solved. A zero interval still checks
    /// source shape, force controls, baths and the original molecular operator.
    /// Invalid sources, exhausted global work, late physical refusal and
    /// cancellation restore the entire original field and clock. No force
    /// sample history or per-step trace is retained or allocated.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stabilized_with_body_force_to<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        target: Time,
        controls: StabilizedReactiveIntervalControl,
        force: BodyForce<'_>,
        force_controls: ForceControl,
        reservoirs: [Option<Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedReactiveIntervalReport, ReactiveError> {
        self.advance_stabilized_to_inner_forced(
            operator,
            target,
            controls,
            reservoirs,
            None,
            Some((force, force_controls)),
            cancelled,
        )
    }
}
