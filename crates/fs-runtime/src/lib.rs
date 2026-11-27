//! Generic runtime binding over the FrankenSim composition kernel
//! (`frankensim-rc-1109-root-3gc12.11.9`).

#![deny(unsafe_code)]
#![warn(missing_docs)]

pub mod descriptor;
mod runtime;

pub use descriptor::{
    InputDescriptor, ProbeDescriptor, RuntimeRefusal, SNAPSHOT_DESCRIPTOR_SCHEMA,
    SnapshotDescriptor, dims_to_unit, inputs_to_json,
};
pub use runtime::{AssemblyRuntime, CustomProbeFn, InputBinding, ProbeBinding};
