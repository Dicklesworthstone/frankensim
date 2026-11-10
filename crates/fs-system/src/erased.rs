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
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal, Receipt};

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
    fn class(&self) -> Class;
    fn has_tangent(&self) -> bool;
    fn stored_energy_j(&self, state: &dyn Any) -> Result<f64, OwnerError>;
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
}

fn state_of<S: 'static>(state: &dyn Any) -> Result<&S, OwnerError> {
    state.downcast_ref::<S>().ok_or_else(OwnerError::state_type)
}

/// A component offering only the mandatory tier.
pub(crate) struct Plain<C>(pub(crate) C);

/// A component that also offers [`PortTangent`].
pub(crate) struct Tangent<C>(pub(crate) C);

macro_rules! mandatory_tier {
    () => {
        fn ports(&self) -> &[PortSpec] {
            self.0.ports()
        }

        fn class(&self) -> Class {
            self.0.class()
        }

        fn stored_energy_j(&self, state: &dyn Any) -> Result<f64, OwnerError> {
            self.0
                .stored_energy_j(state_of::<C::State>(state)?)
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
}
