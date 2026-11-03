//! Complete damped spectrum of the original linear frame under exp(s t).
//!
//! The same original M, K and full material C used by forced harmonic response
//! enter (s^2 M + s C + K) q = 0. The shared fs-modal quadratic-pencil owner
//! performs binary coordinate/time scaling, mass whitening and the bounded
//! fs-la complex eigensolve. No undamped modal projection or diagonal damping
//! fit is introduced. Real roots remain relaxation or rigid motion, not
//! oscillatory resonances.

use super::linear::FrameLinearError;
use super::{Frame2d, Rayleigh};
use crate::SolidError;
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::{ModalError, quadratic_eigenvalues_with_cancel};

/// Explicit cold-query bounds; full native coordinates include fixed rows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FramePoleControls {
    /// Maximum original coordinate count, BEFORE native dense assembly.
    pub maximum_coordinates: usize,
    /// Conservative simultaneous numerical payload of native assembly and the
    /// full quadratic owner. Excludes source metadata and allocator overhead.
    pub maximum_payload_bytes: usize,
    /// Existing complex QR cap for each active eigenvalue window; positive.
    pub maximum_qr_sweeps_per_window: usize,
    /// Conservative arithmetic-work envelope admitted before assembly:
    /// 64*(2*n)^3*(sweeps+1) + 4096*members, with n the full coordinate count.
    /// This counts a bounded operation envelope, not elapsed time or cycles.
    pub maximum_work: u64,
}

/// Full numerical spectrum with the original support and material-domain maps.
#[derive(Clone, Debug, PartialEq)]
pub struct FramePoleSpectrum {
    /// Original node-major [ux, uy, theta] indices; no massive free coordinate
    /// or rotational coupling was condensed or projected away.
    pub free_dofs: Vec<usize>,
    /// All 2*free_dofs.len() rates under exp(s t), sorted by (real, imaginary).
    /// Real(s) is signed growth/decay; |Im(s)| is angular oscillation rate.
    /// Negative real roots are overdamped relaxation, zero roots are retained,
    /// and small positive real parts are never clipped to claim stability.
    /// These values alone do not establish observability or forced peak height.
    pub rates_per_s: Vec<C64>,
    /// Original intersected real-frequency support of the source viscosities.
    /// None means no source-bound material viscosity is present; authored
    /// Rayleigh damping is still included. This is provenance, not admission
    /// of every complex rate or evidence of a material law outside its band.
    pub material_omega_band_rad_s: Option<(f64, f64)>,
}

/// Source, resource, cancellation or shared numerical refusal.
#[derive(Debug)]
pub enum FramePoleError {
    Frame(SolidError),
    Modal(ModalError),
    Invalid(&'static str),
    Budget {
        quantity: &'static str,
        requested: u64,
        cap: u64,
    },
    Cancelled,
}

impl core::fmt::Display for FramePoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-FRAME-POLES: {self:?}")
    }
}
impl std::error::Error for FramePoleError {}

impl From<FrameLinearError> for FramePoleError {
    fn from(error: FrameLinearError) -> Self {
        match error {
            FrameLinearError::Frame(error) => Self::Frame(error),
            FrameLinearError::Invalid(what) => Self::Invalid(what),
            FrameLinearError::Budget {
                quantity,
                requested,
                cap,
            } => Self::Budget {
                quantity,
                requested: requested as u64,
                cap: cap as u64,
            },
            FrameLinearError::Cancelled => Self::Cancelled,
        }
    }
}
impl From<ModalError> for FramePoleError {
    fn from(error: ModalError) -> Self {
        match error {
            ModalError::Cancelled => Self::Cancelled,
            error => Self::Modal(error),
        }
    }
}

impl Frame2d {
    /// Solve the complete material-damped native frame spectrum.
    ///
    /// This explicitly selects the original constant elastic and Kelvin-Voigt
    /// coefficients, together with the authored Rayleigh model, as a continuous
    /// linear model. The resulting complex rates describe that model. A finite
    /// real-frequency material band cannot validate its continuation to every
    /// complex rate: the original band is returned without rejecting out-of-band
    /// roots or concealing fast relaxation. Pure real roots have no oscillation.
    ///
    /// Fixed linear geometry, elastic sections, stationary supports and positive
    /// reference inertia on EVERY free coordinate are required. Consistent
    /// member mass or physical section spin supplies rotational inertia. Rigid
    /// motion is retained; there is no positive-stiffness filter, mass invention,
    /// material-frequency freezing, modal truncation, or transient Nyquist gate.
    ///
    /// Preparation is shared with native harmonic response. All original
    /// material viscosity and translation/rotation couplings remain intact,
    /// including nonproportional losses across heterogeneous members. An
    /// existing authored Rayleigh reference is honored without mutating it.
    ///
    /// Exactly diagonal pencils reuse the shared analytic scalar path, which
    /// preserves the slow pole of strongly overdamped motion. General pencils
    /// use the shared bounded complex QR owner. This is a numerical spectrum,
    /// not an interval enclosure, stability certificate or mode-shape solve.
    /// Accuracy near defective or extremely separated roots remains the shared
    /// owner's limit; this query adds no clipping or regularization.
    ///
    /// Checkpoints bracket native preparation stages, each member assembly,
    /// scalar pole solves and the existing QR elimination/rotation boundaries.
    /// No partial spectrum or accepted-state mutation is published on refusal.
    pub fn resonance_poles(
        &self,
        rayleigh: Rayleigh,
        controls: FramePoleControls,
        cx: &Cx<'_>,
    ) -> Result<FramePoleSpectrum, FramePoleError> {
        cx.checkpoint().map_err(|_| FramePoleError::Cancelled)?;
        if controls.maximum_qr_sweeps_per_window == 0 {
            return Err(FramePoleError::Invalid("positive QR sweep budget required"));
        }
        let n = self.ndof();
        if n > controls.maximum_coordinates {
            return Err(FramePoleError::Budget {
                quantity: "full frame coordinates",
                requested: n as u64,
                cap: controls.maximum_coordinates as u64,
            });
        }
        let d = (n as u64).checked_mul(2);
        let work = d
            .and_then(|d| d.checked_mul(d).and_then(|d2| d2.checked_mul(d)))
            .and_then(|d3| d3.checked_mul(64))
            .and_then(|v| {
                (controls.maximum_qr_sweeps_per_window as u64)
                    .checked_add(1)
                    .and_then(|sweeps| v.checked_mul(sweeps))
            })
            .and_then(|v| {
                (self.elements.len() as u64)
                    .checked_mul(4096)
                    .and_then(|assembly| v.checked_add(assembly))
            })
            .ok_or(FramePoleError::Invalid("frame pole work envelope overflow"))?;
        if work > controls.maximum_work {
            return Err(FramePoleError::Budget {
                quantity: "frame pole arithmetic work",
                requested: work,
                cap: controls.maximum_work,
            });
        }
        let pencil = self.prepare_linear_pencil(
            rayleigh,
            controls.maximum_coordinates,
            controls.maximum_payload_bytes,
            cx,
        )?;
        let rates_per_s = quadratic_eigenvalues_with_cancel(
            &pencil.mass,
            &pencil.damping,
            &pencil.stiffness,
            pencil.free_dofs.len(),
            controls.maximum_qr_sweeps_per_window,
            || cx.checkpoint().is_err(),
        )?;
        cx.checkpoint().map_err(|_| FramePoleError::Cancelled)?;
        Ok(FramePoleSpectrum {
            free_dofs: pencil.free_dofs,
            rates_per_s,
            material_omega_band_rad_s: pencil.material_omega_band_rad_s,
        })
    }
}
