//! Generic runtime binding over the FrankenSim composition kernel
//! (`frankensim-rc-1109-root-3gc12.11.9`).

#![deny(unsafe_code)]
#![warn(missing_docs)]

pub mod descriptor;
mod draft;
pub mod live_source;
mod probes;
mod runtime;

pub use descriptor::{
    InputDescriptor, ProbeDescriptor, RuntimeRefusal, SNAPSHOT_DESCRIPTOR_SCHEMA,
    SnapshotDescriptor, inputs_to_json,
};
pub use fs_fsim::{MAX_PACKS, UNIT_DIMS, dims_to_unit};
pub use live_source::{LiveInputSpec, LiveSource, RUNTIME_INPUT_FAMILY};
pub use runtime::{AssemblyRuntime, CustomProbeFn, InputBinding, ProbeBinding};
