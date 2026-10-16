//! Target-specific external authority values for the BIPOP full-study replay.
//!
//! This data intentionally lives outside the source snapshot hashed by the
//! replay fixture. Embedding an expected fixture hash in source that contributes
//! to that same fixture hash would demand an accidental cryptographic fixed
//! point. The snapshotted test root owns all comparison semantics; this module
//! supplies only reviewed target data. `None` means that this target has emitted
//! sentinels but has not completed capture plus same-selector reverification.

pub(super) const EXPECTED_FOR_TARGET: Option<[&str; 5]> = None;
