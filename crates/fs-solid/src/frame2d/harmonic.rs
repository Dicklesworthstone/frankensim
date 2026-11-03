//! Forced resonance of the original linear frame, without a modal truncation.
//!
//! The native reference mass, elastic stiffness and complete material viscosity
//! enter `K - omega^2 M - i omega C` under `exp(-i omega t)`. Authored Rayleigh
//! damping is optional through explicit zero coefficients; it is not inferred
//! from a material name or substituted for the original section viscosity.
//! Stationary supports use the frame's original free-coordinate map.

use super::{Frame2d, Geometry, Rayleigh, SectionModel};
use crate::SolidError;
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponse, HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};

/// Original frame/source refusal or shared forced-response refusal.
#[derive(Debug)]
pub enum FrameHarmonicError {
    Frame(SolidError),
    Harmonic(HarmonicResponseError),
}
impl core::fmt::Display for FrameHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-FRAME-HARMONIC: {self:?}")
    }
}
impl std::error::Error for FrameHarmonicError {}
impl From<SolidError> for FrameHarmonicError {
    fn from(error: SolidError) -> Self {
        match error {
            SolidError::Cancelled => Self::Harmonic(HarmonicResponseError::Cancelled),
            error => Self::Frame(error),
        }
    }
}
impl From<HarmonicResponseError> for FrameHarmonicError {
    fn from(error: HarmonicResponseError) -> Self {
        Self::Harmonic(error)
    }
}
fn invalid(what: &'static str) -> FrameHarmonicError {
    HarmonicResponseError::Invalid(what).into()
}
fn poll(cx: &Cx<'_>) -> Result<(), FrameHarmonicError> {
    cx.checkpoint()
        .map_err(|_| HarmonicResponseError::Cancelled.into())
}

// Native force-based assembly evaluates the two triangles independently.
// Admit only roundoff skew in stiffness-scaled coordinates, then use their
// symmetric mean, as the native modal analysis does. Never drop a coupling.
fn stiffness_entry(
    matrix: &[f64],
    n: usize,
    row: usize,
    col: usize,
) -> Result<f64, FrameHarmonicError> {
    let (a, b) = (matrix[row * n + col], matrix[col * n + row]);
    if !a.is_finite() || !b.is_finite() {
        return Err(invalid("nonfinite frame stiffness"));
    }
    if a == b {
        return Ok(a);
    }
    let si = matrix[row * n + row].abs().sqrt();
    let sj = matrix[col * n + col].abs().sqrt();
    let skew = (a / si / sj - b / si / sj).abs();
    if !skew.is_finite() || skew > 64.0 * n as f64 * f64::EPSILON {
        return Err(invalid(
            "frame stiffness asymmetry exceeds assembly roundoff",
        ));
    }
    Ok(0.5 * a + 0.5 * b)
}

/// One immutable frequency factor, retaining the original frame and supports.
/// Repeated force/receiver queries reuse the same shared native factorization.
pub struct FrameHarmonicOperator<'a> {
    source: &'a Frame2d,
    free_dofs: Vec<usize>,
    response: PreparedHarmonicResponse,
}

impl FrameHarmonicOperator<'_> {
    /// Source geometry, material bindings, supports and accepted state, unchanged.
    pub fn source(&self) -> &Frame2d {
        self.source
    }

    /// Original node-major `[ux, uy, theta]` indices of every retained coordinate.
    /// Translation is in metres and rotation in radians; conjugate forces are
    /// newtons and newton-metres respectively. No massive coordinate is removed.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Solve for peak applied force phasors in `free_dofs()` order. Returned
    /// displacement and velocity use that same map and `exp(-i omega t)`.
    ///
    /// The power account includes the complete original material viscosity AND
    /// the explicitly supplied Rayleigh model once. It is dissipated mechanical
    /// power, not an automatic heat deposit. Fixed coordinates have zero motion;
    /// loads applied there affect reactions, which this response does not report.
    /// This particular periodic solution does not erase an existing transient.
    pub fn solve(
        &self,
        forces: &[C64],
        cx: &Cx<'_>,
    ) -> Result<HarmonicResponse, FrameHarmonicError> {
        Ok(self.response.solve(forces, cx)?)
    }
}

impl Frame2d {
    /// Prepare the full original frame's forced response at one frequency.
    ///
    /// Only fixed linear geometry and elastic sections are admitted. Every free
    /// coordinate requires positive reference inertia: use consistent member
    /// mass or explicit section spin rather than inventing rotational mass or
    /// statically condensing a viscous first-order coordinate. Source viscosity
    /// must admit the actual frequency, including DC when requested. Its full
    /// off-diagonal matrix survives; no fitted modal damping is introduced.
    /// Original rigid modes remain intact at positive forcing frequency. This
    /// query reuses source viscosity assembly without a transient Nyquist test
    /// or an unnecessary natural-mode eigensolve.
    /// Only roundoff-sized native stiffness skew is symmetrized in scratch;
    /// larger asymmetry refuses and the source frame is never rewritten.
    ///
    /// `maximum_coordinates` bounds the FULL frame, including fixed coordinates,
    /// before dense assembly. The payload ceiling covers a conservative envelope
    /// of numerical arrays for native assembly/admission and the shared solve;
    /// source material metadata, caller inputs and allocator overhead are excluded.
    /// Assembly polls per member; existing bounded mass/material preparations are
    /// indivisible calls between checkpoints. No accepted frame state is mutated.
    ///
    /// # Errors
    /// Unsupported geometry/section/mass, source-frequency refusal, resource
    /// limits, nonfinite coefficients, nonpassive operators, singular lossless
    /// poles, cancellation or the original shared solve's numerical refusal.
    pub fn prepare_harmonic_response(
        &self,
        omega_rad_s: f64,
        rayleigh: Rayleigh,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<FrameHarmonicOperator<'_>, FrameHarmonicError> {
        poll(cx)?;
        let n = self.ndof();
        if n == 0
            || !omega_rad_s.is_finite()
            || omega_rad_s < 0.0
            || ![rayleigh.a0, rayleigh.a1]
                .into_iter()
                .all(|x| x.is_finite() && x >= 0.0)
            || !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(invalid(
                "nonempty frame, nonnegative finite frequency/damping and residual budget in (0,1) required",
            ));
        }
        if n > controls.maximum_coordinates {
            return Err(HarmonicResponseError::Budget {
                quantity: "full frame coordinates",
                requested: n,
                cap: controls.maximum_coordinates,
            }
            .into());
        }
        // Includes the shared owner's 128*nf^2 + 128*nf envelope, the native
        // source assembly, full M/K/C, three free matrices and bounded
        // per-element numerical trials (at most seven section stations).
        let bytes = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(512))
            .and_then(|x| n.checked_mul(512).and_then(|y| x.checked_add(y)))
            .and_then(|x| {
                self.elements
                    .len()
                    .checked_mul(2048)
                    .and_then(|y| x.checked_add(y))
            })
            .ok_or_else(|| invalid("frame harmonic numerical payload overflow"))?;
        if bytes > controls.maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "frame harmonic numerical payload bytes",
                requested: bytes,
                cap: controls.maximum_payload_bytes,
            }
            .into());
        }
        self.require_untimed_materials()?;
        if self.geometry != Geometry::Linear
            || self.elements.iter().any(|element| {
                element.sections.iter().any(|section| {
                    !matches!(section, SectionModel::Elastic { ea, ei }
                    if ea.is_finite() && *ea > 0.0 && ei.is_finite() && *ei > 0.0)
                })
            })
        {
            return Err(invalid(
                "harmonic frame requires fixed linear geometry and elastic sections",
            ));
        }
        let mass = self.reference_mass()?;
        let (free_dofs, massless) = self.mass_partition(&mass);
        if free_dofs.is_empty() || !massless.is_empty() {
            return Err(invalid(
                "harmonic frame requires positive inertia on every free coordinate; use consistent member mass",
            ));
        }
        poll(cx)?;
        let (_, stiffness, _) = self.assemble_controlled(&self.u, &mut || {
            cx.checkpoint().map_err(|_| SolidError::Cancelled)
        })?;
        let material = self.assemble_material_damping()?;
        if let Some(operator) = &material {
            let (low, high) = operator.omega_band_rad_s();
            if omega_rad_s < low || omega_rad_s > high {
                return Err(invalid(
                    "harmonic frequency lies outside original frame viscosity support",
                ));
            }
        }
        poll(cx)?;
        let reference = self
            .rayleigh_reference_stiffness
            .as_ref()
            .unwrap_or(&stiffness);
        let nf = free_dofs.len();
        let mut m = vec![0.0; nf * nf];
        let mut k = vec![0.0; nf * nf];
        let mut c = vec![0.0; nf * nf];
        for (i, &row) in free_dofs.iter().enumerate() {
            poll(cx)?;
            for (j, &col) in free_dofs[..=i].iter().enumerate() {
                let full = row * n + col;
                let reduced = i * nf + j;
                m[reduced] = mass.entry(row, col);
                k[reduced] = stiffness_entry(&stiffness, n, row, col)?;
                let reference_entry = if rayleigh.a1 == 0.0 {
                    0.0
                } else {
                    stiffness_entry(reference, n, row, col)?
                };
                c[reduced] = rayleigh.a0 * m[reduced]
                    + rayleigh.a1 * reference_entry
                    + material
                        .as_ref()
                        .map_or(0.0, |operator| operator.entries()[full]);
                m[j * nf + i] = m[reduced];
                k[j * nf + i] = k[reduced];
                c[j * nf + i] = c[reduced];
            }
        }
        let response =
            PreparedHarmonicResponse::try_new(&m, &c, &k, nf, omega_rad_s, controls, cx)?;
        poll(cx)?;
        Ok(FrameHarmonicOperator {
            source: self,
            free_dofs,
            response,
        })
    }
}
