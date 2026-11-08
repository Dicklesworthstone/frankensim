//! # fs-control — sampled plant control
//!
//! The synchronous current controller executes a fixed-period, power-invariant
//! d/q PI law against measured flux, shaft speeds and DC-link voltage. It emits
//! actual averaged three-phase leg duties together with a proposed integrator
//! state. The caller accepts that state only after the corresponding held-input
//! plant interval succeeds, so cancellation or a rejected interval is retryable.
//!
//! The outer speed controller requests signed torque within authored torque
//! and sampled motoring/regeneration power envelopes. A machine's original
//! constitutive law owns torque-to-current allocation and its achieved torque
//! supplies outer-loop anti-windup. Both proposed controller states are accepted
//! only with the successful plant interval. Authored sampled thermal derating
//! adjusts command capability; it does not provide hardware protection.
//! Sample scheduling, observers, gain synthesis and closed-loop stability
//! certification remain separate capabilities.

#![deny(unsafe_code)]

/// Fixed-period synchronous current control with physical inverter modulation.
pub mod synchronous;

/// Relative-shaft speed regulation with signed torque/power limits.
pub mod speed;

/// Authored sampled thermal command derating with an explicit trip latch.
pub mod thermal;

/// A refused controller construction or sample. No proposal is published on an
/// error, and no caller-owned state is modified.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ControlError {
    /// An authored controller parameter is outside its documented domain.
    InvalidParameter {
        /// The physical parameter or constraint that failed.
        what: &'static str,
    },
    /// A sampled measurement, request or restored state is outside its domain.
    InvalidSample {
        /// The measurement or state constraint that failed.
        what: &'static str,
    },
    /// Finite inputs produced an unrepresentable intermediate or result.
    NumericalRange {
        /// The physical computation that cannot be represented.
        what: &'static str,
    },
    /// The cancellation callback requested refusal before publication.
    Cancelled,
}

impl std::fmt::Display for ControlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidParameter { what } => write!(f, "invalid controller parameter: {what}"),
            Self::InvalidSample { what } => write!(f, "invalid controller sample: {what}"),
            Self::NumericalRange { what } => write!(f, "controller numerical range: {what}"),
            Self::Cancelled => f.write_str("controller sample cancelled before publication"),
        }
    }
}

impl std::error::Error for ControlError {}
