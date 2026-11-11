//! Composition kernel: components joined by connections, stepped as one.
//!
//! A coupled device is a set of [`fs_ports::component::Component`]s and the
//! connections between their ports. This crate owns what every pairing used
//! to write for itself: the state store with atomic publication, the
//! equations of the connections, the order of evaluation, the interface
//! solve, step negotiation, and the energy and entropy ledger.
//!
//! ```text
//! SystemBuilder::add / add_with_tangent      components and initial states
//!              ::bond / common_effort / common_flow / heat_route / signal
//!              ::transformer / power_map    lossless power conversion
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
//! blocks by a dense interface Newton iteration, using each component's port
//! tangent when it offers one and a forward difference otherwise. It has no
//! metric field transfer, no non-smooth tier, no multirate
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
mod solve;
mod system;
mod wiring;

pub use adaptive::{AdaptiveConfig, AdaptiveError, AdaptiveReport, ErrorNorm, ErrorTolerance};
pub use erased::{OwnerError, StateTypeMismatch, TangentNotOffered};
pub use ledger::Ledger;
pub use system::{
    BoundaryFn, Checkpoint, CheckpointError, Phase, SolveConfig, StepReport, System, SystemBuilder,
    SystemError, heat_sink,
};
pub use wiring::{AdmissionError, HeatSink, NodeId, PortRef, RouteTemperature};
