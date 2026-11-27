//! Type erasure: components of different types in one system.
//!
//! The system stores each component behind [`Node`] and each state as
//! `Box<dyn Any>`. The typed wrapper downcasts on every call; the slots are
//! created from the component's own state type, so a failed downcast means the
//! kernel itself is wrong and is reported as [`OwnerError::state_type`], never
//! as a panic.

use core::any::Any;
use core::fmt;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::component::{
    Class, Component, EndOfStep, LoadLine, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};
use fs_qty::Dims;

/// A value that is both printable and recoverable by its concrete type.
trait DebugAny: Any + fmt::Debug {
    fn as_any(&self) -> &dyn Any;
}

impl<T: Any + fmt::Debug> DebugAny for T {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// A component's own refusal, carried unchanged.
///
/// [`OwnerError::downcast_ref`] recovers the owner's error type.
pub struct OwnerError(Box<dyn DebugAny>);

/// Marker error for a state slot that does not hold the component's state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StateTypeMismatch;

/// Marker error for a tangent asked of a component that offers none. The
/// driver checks before it asks, so this means the kernel itself is wrong.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TangentNotOffered;

/// Marker error for a load line asked of a component that was not added on
/// one. The driver checks before it asks, so this means the kernel itself
/// is wrong.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LineNotOffered;

/// Marker error for a stage residual asked of a component that does not offer
/// [`StageResidual`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StageResidualNotOffered;

impl OwnerError {
    fn new<E: fmt::Debug + 'static>(error: E) -> Self {
        Self(Box::new(error))
    }

    fn state_type() -> Self {
        Self(Box::new(StateTypeMismatch))
    }

    fn tangent_not_offered() -> Self {
        Self(Box::new(TangentNotOffered))
    }

    fn line_not_offered() -> Self {
        Self(Box::new(LineNotOffered))
    }

    fn stage_residual_not_offered() -> Self {
        Self(Box::new(StageResidualNotOffered))
    }

    /// The owner's error, if it is of type `E`.
    #[must_use]
    pub fn downcast_ref<E: 'static>(&self) -> Option<&E> {
        // Deref first: calling `as_any` on the box would select the blanket
        // impl for the box itself and the downcast would never match.
        (*self.0).as_any().downcast_ref::<E>()
    }
}

impl fmt::Debug for OwnerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(&self.0, f)
    }
}

/// Object-safe view of a component.
pub(crate) trait Node {
    fn ports(&self) -> &[PortSpec];
    fn state_dimensions(&self) -> Option<Dims>;
    fn class(&self) -> Class;
    fn end_of_step(&self, port: usize) -> EndOfStep;
    fn has_tangent(&self) -> bool;
    fn clone_state(&self, state: &dyn Any) -> Result<Box<dyn Any>, OwnerError>;
    /// The charter of the component's state codec, if it declares one.
    fn state_charter(&self) -> Option<StateIdentityCharterV2>;
    /// `Ok(false)` if the component has no codec; it then wrote nothing.
    fn encode_state(&self, state: &dyn Any, words: &mut Vec<u64>) -> Result<bool, OwnerError>;
    fn decode_state(&self, words: &[u64]) -> Option<Box<dyn Any>>;
    fn identity_words(&self, words: &mut Vec<u64>);
    fn stored_energy_j(&self, state: &dyn Any) -> Result<f64, OwnerError>;
    fn stored_mass_kg(&self, state: &dyn Any) -> Result<Option<f64>, OwnerError>;
    fn stored_entropy_j_per_k(&self, state: &dyn Any) -> Result<Option<f64>, OwnerError>;
    fn stored_carried_quantity(
        &self,
        state: &dyn Any,
        carried: usize,
    ) -> Result<Option<f64>, OwnerError>;
    fn limit_dt_s(&self, state: &dyn Any, proposed_s: f64) -> Result<f64, OwnerError>;
    fn next_event_s(&self, state: &dyn Any) -> Option<f64>;
    fn propose(
        &self,
        old: &dyn Any,
        dt_s: f64,
        inputs: &[f64],
        new: &mut dyn Any,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, OwnerError>;
    fn admit(&self, old: &dyn Any, new: &dyn Any, receipt: &Receipt) -> Result<(), OwnerError>;
    /// Only called when [`Node::has_tangent`] is true.
    fn jvp(
        &self,
        old: &dyn Any,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), OwnerError>;
    /// Adjoint (vector-Jacobian product) of the port map; only called when
    /// [`Node::has_tangent`] is true.
    fn vjp(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _bar_outputs: &[f64],
        _bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        Err(OwnerError::tangent_not_offered())
    }
    /// The ports of its load line, for a component added on one.
    fn line_ports(&self) -> Option<&[usize]> {
        None
    }
    /// The end inputs of its load line's ports, for a component added on
    /// one that names any.
    fn line_ends(&self) -> &[Option<usize>] {
        &[]
    }
    /// Whether its load line reads its ports' flows only through their
    /// means, for a component added on one.
    fn reads_mean_flows_only(&self) -> bool {
        false
    }
    /// Only called when [`Node::line_ports`] is `Some`.
    #[allow(clippy::too_many_arguments)]
    fn propose_on_line(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _compliance: &[f64],
        _new: &mut dyn Any,
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, OwnerError> {
        Err(OwnerError::line_not_offered())
    }
    /// Whether this component offers the [`StageResidual`] monolithic tier.
    fn has_stage_residual(&self) -> bool {
        false
    }
    /// Number of internal stage unknowns `z` (`0` unless [`Self::has_stage_residual`]).
    fn stage_unknowns(&self) -> usize {
        0
    }
    /// Initial guess `z` for monolithic stage solve.
    fn stage_initial_guess(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &mut [f64],
    ) -> Result<(), OwnerError> {
        Err(OwnerError::stage_residual_not_offered())
    }
    /// Evaluate stage residual `r(z)` and port outputs implied by `z`.
    #[allow(clippy::too_many_arguments)]
    fn stage_residual(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        _r: &mut [f64],
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        Err(OwnerError::stage_residual_not_offered())
    }
    /// Directional derivative of stage residual and outputs along `dz` at `z`.
    #[allow(clippy::too_many_arguments)]
    fn stage_jvp_z(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        _dz: &[f64],
        _dr: &mut [f64],
        _d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        Err(OwnerError::stage_residual_not_offered())
    }
    /// Build candidate state and receipt from converged monolithic unknowns `z`.
    fn stage_accept(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _z: &[f64],
        _new: &mut dyn Any,
        _cx: &Cx<'_>,
    ) -> Result<Receipt, OwnerError> {
        Err(OwnerError::stage_residual_not_offered())
    }
}

fn state_of<S: 'static>(state: &dyn Any) -> Result<&S, OwnerError> {
    state.downcast_ref::<S>().ok_or_else(OwnerError::state_type)
}

/// A component offering only the mandatory tier.
pub(crate) struct Plain<C>(pub(crate) C);

/// A component that also offers [`PortTangent`].
pub(crate) struct Tangent<C>(pub(crate) C);

/// A component stepped on its [`LoadLine`].
pub(crate) struct OnLine<C>(pub(crate) C);

/// A component that offers [`StageResidual`].
pub(crate) struct Stage<C>(pub(crate) C);

/// A component that offers both [`PortTangent`] and [`StageResidual`].
pub(crate) struct TangentStage<C>(pub(crate) C);

macro_rules! mandatory_tier {
    () => {
        fn ports(&self) -> &[PortSpec] {
            self.0.ports()
        }

        fn state_dimensions(&self) -> Option<Dims> {
            self.0.state_dimensions()
        }

        fn class(&self) -> Class {
            self.0.class()
        }

        fn end_of_step(&self, port: usize) -> EndOfStep {
            self.0.end_of_step(port)
        }

        fn clone_state(&self, state: &dyn Any) -> Result<Box<dyn Any>, OwnerError> {
            Ok(Box::new(state_of::<C::State>(state)?.clone()))
        }

        fn state_charter(&self) -> Option<StateIdentityCharterV2> {
            C::STATE_CHARTER
        }

        fn encode_state(&self, state: &dyn Any, words: &mut Vec<u64>) -> Result<bool, OwnerError> {
            Ok(self.0.encode_state(state_of::<C::State>(state)?, words))
        }

        fn decode_state(&self, words: &[u64]) -> Option<Box<dyn Any>> {
            self.0
                .decode_state(words)
                .map(|state| Box::new(state) as Box<dyn Any>)
        }

        fn identity_words(&self, words: &mut Vec<u64>) {
            self.0.identity_words(words);
        }

        fn stored_energy_j(&self, state: &dyn Any) -> Result<f64, OwnerError> {
            self.0
                .stored_energy_j(state_of::<C::State>(state)?)
                .map_err(OwnerError::new)
        }

        fn stored_mass_kg(&self, state: &dyn Any) -> Result<Option<f64>, OwnerError> {
            self.0
                .stored_mass_kg(state_of::<C::State>(state)?)
                .map_err(OwnerError::new)
        }

        fn stored_entropy_j_per_k(&self, state: &dyn Any) -> Result<Option<f64>, OwnerError> {
            self.0
                .stored_entropy_j_per_k(state_of::<C::State>(state)?)
                .map_err(OwnerError::new)
        }

        fn stored_carried_quantity(
            &self,
            state: &dyn Any,
            carried: usize,
        ) -> Result<Option<f64>, OwnerError> {
            self.0
                .stored_carried_quantity(state_of::<C::State>(state)?, carried)
                .map_err(OwnerError::new)
        }

        fn limit_dt_s(&self, state: &dyn Any, proposed_s: f64) -> Result<f64, OwnerError> {
            self.0
                .limit_dt_s(state_of::<C::State>(state)?, proposed_s)
                .map_err(OwnerError::new)
        }

        fn next_event_s(&self, state: &dyn Any) -> Option<f64> {
            state
                .downcast_ref::<C::State>()
                .and_then(|state| self.0.next_event_s(state))
        }

        fn propose(
            &self,
            old: &dyn Any,
            dt_s: f64,
            inputs: &[f64],
            new: &mut dyn Any,
            outputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<Proposal, OwnerError> {
            let old = state_of::<C::State>(old)?;
            let new = new
                .downcast_mut::<C::State>()
                .ok_or_else(OwnerError::state_type)?;
            self.0
                .propose(old, dt_s, inputs, new, outputs, cx)
                .map_err(OwnerError::new)
        }

        fn admit(&self, old: &dyn Any, new: &dyn Any, receipt: &Receipt) -> Result<(), OwnerError> {
            self.0
                .admit(
                    state_of::<C::State>(old)?,
                    state_of::<C::State>(new)?,
                    receipt,
                )
                .map_err(OwnerError::new)
        }
    };
}

impl<C> Node for Plain<C>
where
    C: Component + 'static,
    C::State: 'static,
{
    mandatory_tier!();

    fn has_tangent(&self) -> bool {
        false
    }

    fn jvp(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _d_inputs: &[f64],
        _d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        Err(OwnerError::tangent_not_offered())
    }
}

impl<C> Node for OnLine<C>
where
    C: LoadLine + 'static,
    C::State: 'static,
{
    mandatory_tier!();

    fn has_tangent(&self) -> bool {
        false
    }

    fn jvp(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _d_inputs: &[f64],
        _d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        Err(OwnerError::tangent_not_offered())
    }

    fn line_ports(&self) -> Option<&[usize]> {
        Some(self.0.line_ports())
    }

    fn line_ends(&self) -> &[Option<usize>] {
        self.0.line_ends()
    }

    fn reads_mean_flows_only(&self) -> bool {
        self.0.reads_mean_flows_only()
    }

    fn propose_on_line(
        &self,
        old: &dyn Any,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut dyn Any,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, OwnerError> {
        let old = state_of::<C::State>(old)?;
        let new = new
            .downcast_mut::<C::State>()
            .ok_or_else(OwnerError::state_type)?;
        self.0
            .propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)
            .map_err(OwnerError::new)
    }
}

impl<C> Node for Tangent<C>
where
    C: PortTangent + 'static,
    C::State: 'static,
{
    mandatory_tier!();

    fn has_tangent(&self) -> bool {
        true
    }

    fn jvp(
        &self,
        old: &dyn Any,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        self.0
            .jvp(
                state_of::<C::State>(old)?,
                dt_s,
                inputs,
                d_inputs,
                d_outputs,
                cx,
            )
            .map_err(OwnerError::new)
    }

    fn vjp(
        &self,
        old: &dyn Any,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        self.0
            .vjp(
                state_of::<C::State>(old)?,
                dt_s,
                inputs,
                bar_outputs,
                bar_inputs,
                cx,
            )
            .map_err(OwnerError::new)
    }
}

macro_rules! stage_residual_tier {
    () => {
        fn has_stage_residual(&self) -> bool {
            true
        }

        fn stage_unknowns(&self) -> usize {
            self.0.unknowns()
        }

        fn stage_initial_guess(
            &self,
            old: &dyn Any,
            dt_s: f64,
            inputs: &[f64],
            z: &mut [f64],
        ) -> Result<(), OwnerError> {
            self.0
                .initial_guess(state_of::<C::State>(old)?, dt_s, inputs, z)
                .map_err(OwnerError::new)
        }

        fn stage_residual(
            &self,
            old: &dyn Any,
            dt_s: f64,
            inputs: &[f64],
            z: &[f64],
            r: &mut [f64],
            outputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<(), OwnerError> {
            self.0
                .residual(state_of::<C::State>(old)?, dt_s, inputs, z, r, outputs, cx)
                .map_err(OwnerError::new)
        }

        fn stage_jvp_z(
            &self,
            old: &dyn Any,
            dt_s: f64,
            inputs: &[f64],
            z: &[f64],
            dz: &[f64],
            dr: &mut [f64],
            d_outputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<(), OwnerError> {
            self.0
                .jvp_z(
                    state_of::<C::State>(old)?,
                    dt_s,
                    inputs,
                    z,
                    dz,
                    dr,
                    d_outputs,
                    cx,
                )
                .map_err(OwnerError::new)
        }

        fn stage_accept(
            &self,
            old: &dyn Any,
            dt_s: f64,
            inputs: &[f64],
            z: &[f64],
            new: &mut dyn Any,
            cx: &Cx<'_>,
        ) -> Result<Receipt, OwnerError> {
            let old = state_of::<C::State>(old)?;
            let new = new
                .downcast_mut::<C::State>()
                .ok_or_else(OwnerError::state_type)?;
            self.0
                .accept(old, dt_s, inputs, z, new, cx)
                .map_err(OwnerError::new)
        }
    };
}

impl<C> Node for Stage<C>
where
    C: StageResidual + 'static,
    C::State: 'static,
{
    mandatory_tier!();

    fn has_tangent(&self) -> bool {
        false
    }

    fn jvp(
        &self,
        _old: &dyn Any,
        _dt_s: f64,
        _inputs: &[f64],
        _d_inputs: &[f64],
        _d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        Err(OwnerError::tangent_not_offered())
    }

    stage_residual_tier!();
}

impl<C> Node for TangentStage<C>
where
    C: PortTangent + StageResidual + 'static,
    C::State: 'static,
{
    mandatory_tier!();

    fn has_tangent(&self) -> bool {
        true
    }

    fn jvp(
        &self,
        old: &dyn Any,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        self.0
            .jvp(
                state_of::<C::State>(old)?,
                dt_s,
                inputs,
                d_inputs,
                d_outputs,
                cx,
            )
            .map_err(OwnerError::new)
    }

    fn vjp(
        &self,
        old: &dyn Any,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), OwnerError> {
        self.0
            .vjp(
                state_of::<C::State>(old)?,
                dt_s,
                inputs,
                bar_outputs,
                bar_inputs,
                cx,
            )
            .map_err(OwnerError::new)
    }
    stage_residual_tier!();
}
