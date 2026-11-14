//! Composition kernel: components joined by connections, stepped as one.
//!
//! A coupled device is a set of [`fs_ports::component::Component`]s and the
//! connections between their ports. This crate owns what every pairing used
//! to write for itself: the state store with atomic publication, the
//! equations of the connections, the order of evaluation, the interface
//! solve, step negotiation, and the energy, entropy and declared stream-mass ledger.
//!
//! ```text
//! SystemBuilder::add / add_with_tangent      components and initial states
//!              ::bond / common_effort / common_flow / heat_route / signal
//!              ::transformer / power_map    lossless power conversion
//!              ::power_map_members           one map from several ports
//!              ::common_flow_ratios / common_effort_ratios
//!                                            a junction whose members each
//!                                            sit behind their own ratio
//!              ::boundary                    prescribed inputs
//!              ::seal -> System              admission and plan
//! System::step                               negotiate, solve, audit, admit, publish
//! System::ledger                             window accounts
//! System::checkpoint / restore               accepted-state rewind and replay
//! System::step_adaptive                      state/account error-controlled steps
//! ```
//!
//! [`elements`] holds ideal lumped elements (mass, spring, damper, force,
//! heat capacity, conductance, heater), their power-domain variants, and the
//! dual storage, admittance and flow source for common-effort networks.
//! Mechanical, electrical and thermal devices can be described as data and
//! run through the same kernel.
//!
//! # Version 0
//!
//! The contract in fs-ports is v0 and so is this driver. It solves looped
//! blocks by a Newton iteration on their interface inputs, using each
//! component's port tangent when it offers one and a forward difference
//! otherwise; the inputs that are copies are substituted and the ones that
//! close the loops are eliminated together. It has no
//! non-smooth tier, no multirate
//! stepping, no monolithic solve through `StageResidual` and no snapshot
//! encoding. Step-doubling error control uses an explicitly supplied state
//! norm and convergence order. The remaining work is in beads `.5.3`, `.5.5` and
//! `.5.7` to `.5.10` of `frankensim-rc-1109-root-3gc12`.

#![deny(unsafe_code)]
#![warn(missing_docs)]

mod adaptive;
#[cfg(feature = "constitutive")]
pub mod constitutive;
pub mod elements;
mod erased;
pub mod ledger;
#[cfg(feature = "phs")]
pub mod phs;
#[cfg(feature = "phs")]
pub mod phs_exact;
mod solve;
mod system;
mod wiring;

pub use adaptive::{AdaptiveConfig, AdaptiveError, AdaptiveReport, ErrorNorm, ErrorTolerance};
pub use erased::{LineNotOffered, OwnerError, StateTypeMismatch, TangentNotOffered};
pub use ledger::Ledger;
pub use system::{
    BoundaryFn, Checkpoint, CheckpointError, Phase, SolveConfig, StepReport, System, SystemBuilder,
    SystemError, SystemImage, heat_sink,
};
pub use wiring::{AdmissionError, HeatSink, NodeId, PortRef, RouteTemperature};
