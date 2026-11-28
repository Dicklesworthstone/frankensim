//! One-way baffled-pressure Rayleigh microphone observer for modal surfaces (`.11.9.5`).
//! Reference: `fs_couple::pcm_wav::baffled::BaffledPressure::new_multirate`.
//! No radiation backreaction and no shell/room diffraction.

/// Streaming causal Blackman-Harris cascade decimator.
pub mod decimate;

use crate::NodeId;

/// One free-node planar surface quadrature element (`z = 0`).
#[derive(Debug, Clone, PartialEq)]
pub struct RadiatingPatch {
    /// Planar patch coordinates `[x, y]` [m].
    pub xy_m: [f64; 2],
    /// Lumped P1 nodal area [m^2].
    pub area_m2: f64,
    /// Outward-signed mode shapes `outward_sign * phi_k(x, y)` [1/sqrt(kg)].
    pub shapes: Vec<f64>,
}

/// Lowered radiating modal head metadata attached to [`crate::LoweredSystem`].
#[derive(Debug, Clone, PartialEq)]
pub struct RadiatingModalSurface {
    /// Kernel `ModalStructure` node owning modal velocities.
    pub node: NodeId,
    /// Number of retained modes.
    pub mode_count: usize,
    /// Natural angular frequencies `omega_k` [rad/s] of retained modes.
    pub mode_omegas_rad_s: Vec<f64>,
    /// Outward normal sign (`+1.0` or `-1.0`, from touching `PistonFace :sign`).
    pub outward_sign: f64,
    /// Free-node surface patches (fixed rim omitted).
    pub samples: Vec<RadiatingPatch>,
}
