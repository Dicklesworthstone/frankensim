//! Impact observation reuses the fixed-delay owner shared with multirate
//! pressure alignment. Its original physical-travel API remains crate-local.
pub(super) use crate::pcm_wav::delay::FractionalDelay as Propagation;
