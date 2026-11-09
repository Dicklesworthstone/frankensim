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
//!              ::boundary                    prescribed inputs
//!              ::seal -> System              admission and plan
//! System::step                               negotiate, solve, audit, admit, publish
//! System::ledger                             window accounts
//! ```
//!
//! # Version 0
//!
//! The contract in fs-ports is v0 and so is this driver. It solves looped
//! blocks by a dense interface Newton iteration, using each component's port
//! tangent when it offers one and a forward difference otherwise. It has no
//! transformer or mapped field transfer, no non-smooth tier, no multirate
//! stepping, no monolithic solve through `StageResidual`, no step-size error
//! control and no snapshot encoding. Those are beads `.5.3`, `.5.5` and
//! `.5.7` to `.5.10` of `frankensim-rc-1109-root-3gc12`.

#![deny(unsafe_code)]
#![warn(missing_docs)]

mod erased;
pub mod ledger;
mod solve;
mod system;
mod wiring;

pub use erased::{OwnerError, StateTypeMismatch};
pub use ledger::Ledger;
pub use system::{
    BoundaryFn, Phase, SolveConfig, StepReport, System, SystemBuilder, SystemError, heat_sink,
};
pub use wiring::{AdmissionError, HeatSink, NodeId, PortRef, RouteTemperature};
