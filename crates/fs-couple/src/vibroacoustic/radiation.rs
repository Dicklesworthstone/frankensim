//! Frequency-qualified exterior loading of original physical velocity/effort ports.
//!
//! A radiation impedance is a frequency sample, not a broadband material law.
//! Keep the producer's exact grid and use the shared physical-port Hermitian
//! admission before a sample can enter the coupled equations. No diagonal
//! replacement, passivity repair, interpolation or frozen extrapolation occurs.

use super::{C64, VibroError};
use fs_bem::modal_loading::ModalLoadingBake;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_modal::harmonic::{HarmonicResponseError, admit_passive_port_impedance};

#[derive(Debug, Clone)]
struct RadiationRow {
    omega: f64,
    impedance: Vec<C64>,
}

/// Passive exterior radiation samples in ORIGINAL physical velocity/effort
/// coordinates, with peak `exp(-i omega t)` phasors.
///
/// The structural-only attachment uses mass-normalized modal velocities.
/// [`super::CavityOpeningRadiation`] instead declares a joint basis consisting
/// of those velocities followed by outward opening volume flows. Every matrix
/// row is the work-conjugate opposing effort in the same port order; source
/// traces, SI units and normalization must agree with the chosen attachment.
///
/// Each matrix maps outward port flow to its opposing physical effort. All mutual
/// entries are retained, including reactive coupling and exact loss nullspaces.
/// The full Hermitian part is checked by the existing `fs-modal` physical-port
/// owner, once at preparation. Its bounded dense-work and cancellation rules
/// apply. This is pointwise passivity, not causality or a statement about any
/// frequency between the supplied rows. The producer owns the geometry, medium,
/// mode normalization and spatial convergence; dimensions alone do not establish
/// a match to the structure or openings supplied to [`super::VibroacousticModel`].
#[derive(Debug, Clone)]
pub struct ModalRadiationLoading {
    modes: usize,
    rows: Vec<RadiationRow>,
    source_declaration: Option<ContentHash>,
}

fn poll(cx: &Cx<'_>) -> Result<(), VibroError> {
    cx.checkpoint().map_err(|_| VibroError::Cancelled)
}

fn radiation_error(error: HarmonicResponseError) -> VibroError {
    match error {
        HarmonicResponseError::Cancelled => VibroError::Cancelled,
        error => VibroError::Radiation(error),
    }
}

impl ModalRadiationLoading {
    /// Admit one explicitly evaluated exterior impedance at `omega` [rad/s].
    ///
    /// A subsequent solve at ANY different frequency refuses. For a physical
    /// frequency sweep, use [`Self::from_bem_bake`] to retain every original
    /// frequency row without rebuilding the structural/cavity bases.
    ///
    /// # Errors
    /// Nonpositive/nonfinite frequency, a nonempty square matrix mismatch,
    /// nonfinite/active loading, the shared dense-work cap, or cancellation.
    pub fn try_new(
        omega: f64,
        modes: usize,
        impedance: Vec<C64>,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        poll(cx)?;
        if !omega.is_finite() || omega <= 0.0 {
            return Err(VibroError::BadParameter {
                what: "radiation angular frequency must be finite and positive",
            });
        }
        admit_passive_port_impedance(modes, &impedance, cx).map_err(radiation_error)?;
        poll(cx)?;
        Ok(Self {
            modes,
            rows: vec![RadiationRow { omega, impedance }],
            source_declaration: None,
        })
    }

    /// Retain the complete frequency grid, original matrices and declaration
    /// identity of the native BEM/Rayleigh modal-loading producer.
    ///
    /// Each original row receives the SAME full Hermitian passivity admission
    /// as a directly supplied sample; a producer's approximate reciprocity
    /// allowance cannot hide an active complex direction. No matrix is altered.
    /// Warm FRFs perform only a binary lookup: no repeated material evaluation,
    /// matrix projection or passivity solve. The bake's original aggregate
    /// matrix-entry budget remains the bound on the copied sweep payload.
    ///
    /// # Errors
    /// Any rejected original row or cancellation publishes no partial loading.
    pub fn from_bem_bake(bake: &ModalLoadingBake, cx: &Cx<'_>) -> Result<Self, VibroError> {
        poll(cx)?;
        let mut rows = Vec::with_capacity(bake.rows().len());
        for row in bake.rows() {
            poll(cx)?;
            admit_passive_port_impedance(bake.modes(), row.impedance(), cx)
                .map_err(radiation_error)?;
            rows.push(RadiationRow {
                omega: row.omega(),
                impedance: row.impedance().to_vec(),
            });
        }
        poll(cx)?;
        Ok(Self {
            modes: bake.modes(),
            rows,
            source_declaration: Some(bake.declaration_hash()),
        })
    }

    /// Number of original velocity/effort ports in every impedance matrix.
    #[must_use]
    pub const fn mode_count(&self) -> usize {
        self.modes
    }

    /// The exact producer grid, in its original increasing order [rad/s].
    pub fn angular_frequencies_rad_s(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.rows.iter().map(|row| row.omega)
    }

    /// Original BEM/Rayleigh declaration, when constructed from that producer.
    #[must_use]
    pub const fn source_declaration(&self) -> Option<ContentHash> {
        self.source_declaration
    }

    /// The unmodified full modal matrix at an exact declared frequency.
    ///
    /// # Errors
    /// Refuses unless `omega` is exactly one of the original sample frequencies.
    pub fn at_omega(&self, omega: f64) -> Result<&[C64], VibroError> {
        self.rows
            .binary_search_by(|row| row.omega.total_cmp(&omega))
            .map(|i| self.rows[i].impedance.as_slice())
            .map_err(|_| VibroError::BadParameter {
                what: "radiation load has no sample at the solved angular frequency",
            })
    }
}

/// Contract only the Hermitian work form. Large reactive exchange cancels
/// analytically before multiplication, so it cannot manufacture radiated power.
/// A retained prefix uses the original principal block and original stride.
pub(super) fn outward_power(
    impedance: &[C64],
    modes: usize,
    omega: f64,
    displacement: &[C64],
) -> f64 {
    hermitian_work(impedance, modes, displacement.len(), |r| {
        let b = displacement[r];
        C64::new(omega * b.im, -omega * b.re)
    })
}

/// The same work form for a COMPLETE original physical-port flow vector.
pub(super) fn flow_power(impedance: &[C64], ports: usize, flow: &[C64]) -> f64 {
    hermitian_work(impedance, ports, flow.len(), |r| flow[r])
}

fn hermitian_work(
    impedance: &[C64],
    modes: usize,
    retained: usize,
    flow: impl Fn(usize) -> C64,
) -> f64 {
    let dd = Dd::from_f64;
    let mut work = Dd::ZERO;
    for r in 0..retained {
        let v = flow(r);
        work = work
            + dd(0.5)
                * dd(impedance[r * modes + r].re)
                * (dd(v.re) * dd(v.re) + dd(v.im) * dd(v.im));
        for s in r + 1..retained {
            let w = flow(s);
            let forward = impedance[r * modes + s];
            let reverse = impedance[s * modes + r];
            let re = dd(f64::midpoint(forward.re, reverse.re));
            let im = dd(f64::midpoint(forward.im, -reverse.im));
            work = work
                + re * (dd(v.re) * dd(w.re) + dd(v.im) * dd(w.im))
                + im * (dd(v.im) * dd(w.re) - dd(v.re) * dd(w.im));
        }
    }
    work.to_f64()
}

#[cfg(test)]
mod tests;
