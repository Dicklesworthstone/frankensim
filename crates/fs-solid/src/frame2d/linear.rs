//! Shared native linear frame pencil for forced response and free decay.
//!
//! This is the original assembled coordinate system with stationary supports,
//! not a modal image. Keeping one preparation prevents different resonance
//! queries from silently selecting different inertia, viscosity or couplings.

use super::{Frame2d, Geometry, Rayleigh, SectionModel};
use crate::SolidError;
use fs_exec::Cx;

pub(super) enum FrameLinearError {
    Frame(SolidError),
    Invalid(&'static str),
    Budget {
        quantity: &'static str,
        requested: usize,
        cap: usize,
    },
    Cancelled,
}
impl From<SolidError> for FrameLinearError {
    fn from(error: SolidError) -> Self {
        match error {
            SolidError::Cancelled => Self::Cancelled,
            error => Self::Frame(error),
        }
    }
}
fn invalid(what: &'static str) -> FrameLinearError {
    FrameLinearError::Invalid(what)
}
fn poll(cx: &Cx<'_>) -> Result<(), FrameLinearError> {
    cx.checkpoint().map_err(|_| FrameLinearError::Cancelled)
}

// Native force-based assembly evaluates the two triangles independently.
// Admit only roundoff skew in stiffness-scaled coordinates, then use their
// symmetric mean, as the native modal analysis does. Never drop a coupling.
fn stiffness_entry(
    matrix: &[f64],
    n: usize,
    row: usize,
    col: usize,
) -> Result<f64, FrameLinearError> {
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

pub(super) struct FrameLinearPencil {
    pub free_dofs: Vec<usize>,
    pub mass: Vec<f64>,
    pub stiffness: Vec<f64>,
    pub damping: Vec<f64>,
    pub material_omega_band_rad_s: Option<(f64, f64)>,
}

impl Frame2d {
    /// The payload bound includes native assembly and the numerical arrays of
    /// either the shared harmonic factor or the full quadratic pole solve.
    /// Source metadata, caller inputs and allocator overhead are excluded.
    pub(super) fn prepare_linear_pencil(
        &self,
        rayleigh: Rayleigh,
        maximum_coordinates: usize,
        maximum_payload_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<FrameLinearPencil, FrameLinearError> {
        poll(cx)?;
        let n = self.ndof();
        if n == 0
            || ![rayleigh.a0, rayleigh.a1]
                .into_iter()
                .all(|x| x.is_finite() && x >= 0.0)
        {
            return Err(invalid(
                "nonempty frame and finite nonnegative authored Rayleigh coefficients required",
            ));
        }
        if n > maximum_coordinates {
            return Err(FrameLinearError::Budget {
                quantity: "full frame coordinates",
                requested: n,
                cap: maximum_coordinates,
            });
        }
        // Covers native full M/K/C assembly and three free matrices, plus
        // either shared harmonic preparation or both complex 2*nf companion
        // copies in the quadratic owner. Per-member trials have at most seven
        // section stations. The two consumers do not run simultaneously.
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
            .ok_or_else(|| invalid("frame linear numerical payload overflow"))?;
        if bytes > maximum_payload_bytes {
            return Err(FrameLinearError::Budget {
                quantity: "frame linear numerical payload bytes",
                requested: bytes,
                cap: maximum_payload_bytes,
            });
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
                "linear frame requires fixed linear geometry and elastic sections",
            ));
        }
        let mass = self.reference_mass()?;
        let (free_dofs, massless) = self.mass_partition(&mass);
        if free_dofs.is_empty() || !massless.is_empty() {
            return Err(invalid(
                "linear frame requires positive inertia on every free coordinate; use consistent member mass",
            ));
        }
        poll(cx)?;
        let (_, stiffness, _) = self.assemble_controlled(&self.u, &mut || {
            cx.checkpoint().map_err(|_| SolidError::Cancelled)
        })?;
        let material = self.assemble_material_damping()?;
        let material_omega_band_rad_s = material
            .as_ref()
            .map(|operator| operator.omega_band_rad_s());
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
        Ok(FrameLinearPencil {
            free_dofs,
            mass: m,
            stiffness: k,
            damping: c,
            material_omega_band_rad_s,
        })
    }
}
