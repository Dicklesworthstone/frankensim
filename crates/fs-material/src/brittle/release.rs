//! Native admission of a completely failed plane before releasing continuity.
//!
//! Complete loss of opening/shear transfer does not imply a traction-free face:
//! a failed band still supports undamaged compression and in-plane elasticity.
//! This observation admits only a face aligned with the declared material plane
//! whose accepted history has fully failed and is not in normal compression.
//! The actual full stress supplies the independent physical face-traction gate.
//!
//! The original bulk history remains the sole owner of its in-plane storage and
//! cumulative fracture loss. No cohesive point, extra fracture energy, impulse,
//! or mechanical work is created. Geometry and motion transfer belong to the
//! structural owner. Partial-damage insertion and contact after a released face
//! closes require their own physical operators; a shifted band coordinate is
//! not a geometric crack-face gap.

use super::{BrittleCrackBand, BrittleError, BrittleState, det, finite, normalize};
use crate::Voigt;

/// Immutable native observation for a proposed physical fracture face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittlePlaneRelease {
    /// Requested physical normal after normalization; its sign is retained.
    pub face_normal: [f64; 3],
    /// Actual accepted full stress in tensor-shear Voigt order, Pa.
    pub stress_pa: Voigt,
    /// Actual `stress * face_normal` in global coordinates, Pa. No component
    /// is clipped to zero to satisfy the caller's numerical traction gate.
    pub traction_pa: [f64; 3],
    /// Original recoverable storage remaining in the in-plane elastic modes.
    /// Retain this in the original bulk state; do not add it as surface energy.
    pub parallel_stored_energy_j_m3: f64,
    /// Original complete fracture history, including earlier thermal transfers.
    /// Retain this loss once in the bulk history, without resetting it to the
    /// current coefficient envelope or duplicating it on a new interface.
    pub fracture_dissipation_j_m3: f64,
}

impl BrittleCrackBand {
    /// Admit a physical face on this ORIGINAL fully failed, unloaded plane.
    ///
    /// `face_normal` may have either orientation and any finite nonzero scale.
    /// `normal_alignment_tolerance` is a dimensionless chord-distance gate in
    /// `[0, 1)`: the Euclidean distance between normalized directions after
    /// choosing the nearer sign. `traction_tolerance_pa` is a finite nonnegative
    /// gate on the Euclidean norm of the actual full stress times that normal.
    ///
    /// Failure and absence of normal compression are exact material-state
    /// requirements. A loose traction tolerance cannot release surviving bonds
    /// or hide active compression. A held native trial verifies owner/domain
    /// and recovers the stress; the supplied state is never changed. This
    /// allocation-free fixed-size operation leaves cancellation to its batch
    /// or element owner.
    pub fn release_failed_plane(
        &self,
        state: &BrittleState<'_>,
        face_normal: [f64; 3],
        normal_alignment_tolerance: f64,
        traction_tolerance_pa: f64,
    ) -> Result<BrittlePlaneRelease, BrittleError> {
        if !normal_alignment_tolerance.is_finite()
            || !(0.0..1.0).contains(&normal_alignment_tolerance)
            || !traction_tolerance_pa.is_finite()
            || traction_tolerance_pa < 0.0
        {
            return Err(BrittleError::Invalid(
                "finite nonnegative traction gate and direction gate in [0, 1) required",
            ));
        }
        let held = self.trial(state.strain(), state)?;
        if !state.failed() {
            return Err(BrittleError::Invalid(
                "plane release requires complete accepted cohesive failure",
            ));
        }
        if state.cohesive_point().jump_m()[0] < 0.0 {
            return Err(BrittleError::Invalid(
                "failed plane still carries normal compression; contact transfer required",
            ));
        }
        let face_normal = normalize(face_normal)?;
        let normal = self.normal();
        let dot = (0..3).map(|i| face_normal[i] * normal[i]).sum::<f64>();
        let sign = if dot < 0.0 { -1.0 } else { 1.0 };
        let direction_error = det::sqrt(
            (0..3)
                .map(|i| (face_normal[i] - sign * normal[i]).powi(2))
                .sum::<f64>(),
        );
        if direction_error > normal_alignment_tolerance {
            return Err(BrittleError::Invalid(
                "physical face normal is not aligned with the failed material plane",
            ));
        }
        let stress_pa = held.stress_pa;
        // Scale before contraction so a finite large in-plane stress does not
        // overflow the intermediate sum of cancelling global tensor entries.
        let scale = stress_pa.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        let mut traction_pa = [0.0; 3];
        if scale > 0.0 {
            let s = stress_pa.map(|x| x / scale);
            let [nx, ny, nz] = face_normal;
            traction_pa = [
                finite(scale * (s[0] * nx + s[3] * ny + s[5] * nz))?,
                finite(scale * (s[3] * nx + s[1] * ny + s[4] * nz))?,
                finite(scale * (s[5] * nx + s[4] * ny + s[2] * nz))?,
            ];
        }
        let scale = traction_pa.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        let magnitude = if scale == 0.0 {
            0.0
        } else {
            finite(scale * det::sqrt(traction_pa.iter().map(|x| (x / scale).powi(2)).sum()))?
        };
        if magnitude > traction_tolerance_pa {
            return Err(BrittleError::Invalid(
                "actual face traction exceeds the release gate",
            ));
        }
        Ok(BrittlePlaneRelease {
            face_normal,
            stress_pa,
            traction_pa,
            parallel_stored_energy_j_m3: finite(self.parallel_energy(state.strain()))?,
            fracture_dissipation_j_m3: finite(state.dissipated_j_m3())?,
        })
    }
}
