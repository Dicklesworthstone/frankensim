//! Steady annular source fields driving actual independent shaft momentum.

use super::super::swirl::wake::WakeModel;
use super::{
    Refusal, RotorState, ShaftEvaluation, ShaftInputs, ShaftMotion, ShaftPropulsion,
    ShaftPropulsionStep,
};
use crate::propcoupling::swirl::wake::{
    ContinuedWakeStep, CoupledWakeStep, WakeControls, WakeTransport, WakeWarmStart,
};
use fs_airscrew::swirl::Controls;
use fs_wing::nonlinear::continuation::{ContinuationControls, WingContinuation};

/// Full downstream angular wash on independently driven, braked or free shafts.
///
/// Each actual axial momentum determines the native blade speed and signed
/// torque. The same annular source fields as `WakePropulsion` then supply the
/// wing's incident air, including material contraction and radial transport
/// when explicitly selected. The retained `ShaftMotion` alone owns carrier
/// recoil, shaft momentum, bearing losses, phase and motor work.
///
/// This is a steady source-field selection. Native wing-branch continuation
/// and numerical field seeds are explicit; neither adds physical wake storage.
/// Native finite-window, disjoint-source and power-budget limits apply.
pub struct WakeShaftPropulsion<'p, 'a, const N: usize> {
    propulsion: &'p ShaftPropulsion<'a, N>,
    model: WakeModel<'a>,
}

impl<'p, 'a, const N: usize> WakeShaftPropulsion<'p, 'a, N> {
    /// Bind the admitted shaft assembly to the existing fixed-cylinder field.
    /// Native rotor and field controls are admitted by the aerodynamic solve.
    #[must_use]
    pub fn new(
        propulsion: &'p ShaftPropulsion<'a, N>,
        native: Controls,
        wake: WakeControls,
    ) -> Self {
        let mechanical_axes = propulsion
            .geometry
            .mounts
            .map(|mount| propulsion.assembly.rotors[mount.rotor_index].axis_body);
        Self {
            propulsion,
            model: WakeModel::new(propulsion.geometry, mechanical_axes, native, wake),
        }
    }

    /// Select fixed-cylinder or material-streamtube transport explicitly.
    #[must_use]
    pub fn with_transport(mut self, transport: WakeTransport) -> Self {
        self.model = self.model.with_transport(transport);
        self
    }

    /// Solve native wing/rotor/field coupling at each actual signed p/J speed.
    /// Spin direction follows the physical shaft axis and blade handedness.
    /// The resulting axial aerodynamic torque belongs only to its external
    /// shaft port; the downstream wing wrench adds no second shaft reaction.
    ///
    /// # Errors
    /// Original stage, wing, rotor, source-field and budget refusals pass
    /// through without advancing mechanical state or modifying caller seeds.
    pub fn evaluate(
        &self,
        state: &RotorState<N>,
        inputs: ShaftInputs<N>,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<&WakeWarmStart>,
    ) -> Result<ShaftEvaluation<N, CoupledWakeStep>, Refusal> {
        self.propulsion
            .evaluate_with(state, inputs, wind_world, |stage| {
                self.model.solve(stage, warm)
            })
    }

    /// Advance the retained native shaft owner from complete initial and
    /// predicted-midpoint annular-field loads. Only an accepted aerodynamic
    /// stage becomes the local numerical seed for the next stage.
    ///
    /// # Errors
    /// Inventory mismatch refuses before callbacks. Any initial, midpoint,
    /// field-budget or late native-mechanics refusal retains both canonical
    /// states and caller seeds. Callback side effects are outside this
    /// transaction. No complete fluid/body energy closure is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn step<C, W>(
        &self,
        motion: &mut ShaftMotion<N>,
        t_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<&WakeWarmStart>,
        inputs: C,
        wind_world: W,
    ) -> Result<ShaftPropulsionStep<N, CoupledWakeStep>, Refusal>
    where
        C: FnMut(f64, &RotorState<N>) -> Result<ShaftInputs<N>, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut stage_warm = warm.cloned();
        self.propulsion.step_with(
            motion,
            t_s,
            angular_origin_world,
            inputs,
            wind_world,
            |state, inputs, wind| {
                let evaluation = self.evaluate(state, inputs, wind, stage_warm.as_ref())?;
                stage_warm = Some(evaluation.coupled.warm_start());
                Ok(evaluation)
            },
        )
    }
    /// Continue the native wing branch under each actual independent p/J speed
    /// and the complete annular field. `previous=None` establishes an initial
    /// accepted branch without advancing mechanical state.
    ///
    /// # Errors
    /// Original field and wing-continuation refusals leave caller inputs intact.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_continued(
        &self,
        state: &RotorState<N>,
        inputs: ShaftInputs<N>,
        wind_world: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
        warm: Option<&WakeWarmStart>,
        previous: Option<&WingContinuation>,
        domain: ContinuationControls,
    ) -> Result<ShaftEvaluation<N, ContinuedWakeStep>, Refusal> {
        self.propulsion
            .evaluate_with(state, inputs, wind_world, |stage| {
                self.model.solve_continued(stage, warm, previous, domain)
            })
    }

    /// Advance native mechanics and the accepted wing anchor transactionally.
    /// Initial and midpoint field solves each hold their wing anchor fixed
    /// through all nested trials. Only a complete accepted initial solution
    /// becomes the local midpoint anchor. Caller state publishes once.
    ///
    /// # Errors
    /// Any field, wing, source-budget or late native-mechanics refusal retains
    /// both mechanical states, caller wing anchor and caller field seeds.
    #[allow(clippy::too_many_arguments)]
    pub fn step_continued<C, W>(
        &self,
        motion: &mut ShaftMotion<N>,
        wing: &mut WingContinuation,
        t_s: f64,
        angular_origin_world: [f64; 3],
        warm: Option<&WakeWarmStart>,
        domain: ContinuationControls,
        inputs: C,
        wind_world: W,
    ) -> Result<ShaftPropulsionStep<N, ContinuedWakeStep>, Refusal>
    where
        C: FnMut(f64, &RotorState<N>) -> Result<ShaftInputs<N>, Refusal>,
        W: Fn(f64, [f64; 3]) -> Result<[f64; 3], Refusal>,
    {
        let mut stage_warm = warm.cloned();
        let mut stage_wing = wing.clone();
        let step = self.propulsion.step_with(
            motion,
            t_s,
            angular_origin_world,
            inputs,
            wind_world,
            |state, inputs, wind| {
                let evaluation = self.evaluate_continued(
                    state,
                    inputs,
                    wind,
                    stage_warm.as_ref(),
                    Some(&stage_wing),
                    domain,
                )?;
                stage_warm = Some(evaluation.coupled.solution.warm_start());
                stage_wing = evaluation.coupled.wing_continuation.clone();
                Ok(evaluation)
            },
        )?;
        // All fallible field/anchor work ran inside the native callbacks.
        *wing = stage_wing;
        Ok(step)
    }
}
