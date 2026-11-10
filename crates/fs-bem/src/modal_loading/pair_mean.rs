//! The frequency sampling of the multi-port radiation law: the loading between two modes is taken
//! at the pair's mean frequency.
//!
//! A [`ModalLoadingBake`] holds frequency samples of the work impedance `Z(w)` between modes
//! (exp(-i w t): resistance `Re Z`, added mass `-Im Z / w`). A modal time step needs constant
//! matrices, so this law takes entry `(r, s)` at `w_rs = (w_r + w_s) / 2`, the mean of the two
//! modes' frequencies. A mode's own entry is taken at its own frequency. The mean lies inside mode
//! `r`'s half-power band only when `|w_r - w_s| <= 2 b_r`: for a bar's mode and its tube's mode
//! 1.9 Hz above it, inside the fast tube mode's band and outside the slow bar mode's. The law takes
//! the reciprocal part `(Z_rs + Z_sr) / 2`, so both matrices are symmetric bit for bit; what it
//! drops is the bake's reciprocity defect, which the bake certified below its bound.
//!
//! Entries between two modes whose radiation one owner already carries are set to zero, so that
//! radiation is not counted twice. The modes of one air column are such a set: the column law loads
//! its open end with its own end correction and loss.
//!
//! The error of the sampling is stated, not hidden, whether or not the pair mean lies in the band.
//! For each mode `r`, [`BandChange`] holds the largest change of `-Im Z_rs / w` and of `Re Z_rs`
//! over the mode's row, between the pair-mean sample and the two edges `w_r - b_r` and `w_r + b_r`
//! of the mode's half-power band, where the mode responds. It samples the two edges only: that is
//! the largest change over the band when each entry is monotone on each half of the band, which
//! the law does not check.
//!
//! No claim: the law does not interpolate between samples (the bake must declare every frequency
//! [`PairMeanLaw::frequencies`] lists); it is not a causal realization of `Z(w)`; a mode's band is
//! supplied by the caller from its loss, not recomputed from the loading.

use super::{ModalLoadingBake, ModalLoadingError, ModalLoadingRow};

/// The modes the law samples for: their frequencies, half-power bands and owners.
#[derive(Debug, Clone, PartialEq)]
pub struct PairMeanLaw {
    /// Angular frequency of each mode [rad/s].
    omegas: Vec<f64>,
    /// Half-width of each mode's half-power band [rad/s]: its amplitude decay rate.
    half_bands: Vec<f64>,
    /// The owner that already carries each mode's radiation, if any.
    owners: Vec<Option<usize>>,
}

/// The constant loading matrices of the law and the statement of their sampling error.
#[derive(Debug, Clone, PartialEq)]
pub struct PairMeanLoading {
    /// `-Im Z_rs / w_rs` from the reciprocal part, row-major `n x n` [kg per kg of modal mass for
    /// mass-normalized traces].
    pub added_mass: Vec<f64>,
    /// `Re Z_rs` from the reciprocal part, row-major `n x n` [the work impedance's unit, 1/s for
    /// mass-normalized traces].
    pub resistance: Vec<f64>,
    /// One entry per mode, in mode order.
    pub band: Vec<BandChange>,
}

/// How much one mode's row of the loading changes across the mode's half-power band, read at the
/// band's two edges.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandChange {
    /// The mode.
    pub mode: usize,
    /// Largest `|(-Im Z_rs / w)(edge) - added_mass[r][s]|` over the row's entries kept and the two
    /// band edges.
    pub added_mass: f64,
    /// Largest `|Re Z_rs(edge) - resistance[r][s]|` over the same.
    pub resistance: f64,
    /// Largest `|added_mass[r][s]|` over the row's entries kept: the scale of the change.
    pub largest_added_mass: f64,
    /// Largest `|resistance[r][s]|` over the row's entries kept.
    pub largest_resistance: f64,
}

impl PairMeanLaw {
    /// The law for modes at `omegas` [rad/s] with half-power half-widths `half_bands` [rad/s] and
    /// `owners`: two modes with the same `Some(owner)` have their mutual entries, and their own
    /// entries, set to zero.
    ///
    /// # Errors
    /// [`ModalLoadingError::Invalid`] for no mode, lists of different lengths, a frequency that is
    /// not positive and finite, or a half-width that is negative, not finite, or reaches zero
    /// frequency.
    #[expect(
        clippy::result_large_err,
        reason = "ModalLoadingError is the module's refusal type (rule 1: no new error type); \
                  boxing it would change every modal_loading signature"
    )]
    pub fn new(
        omegas: &[f64],
        half_bands: &[f64],
        owners: &[Option<usize>],
    ) -> Result<Self, ModalLoadingError> {
        if omegas.is_empty() {
            return Err(ModalLoadingError::Invalid(
                "the law needs at least one mode",
            ));
        }
        if half_bands.len() != omegas.len() || owners.len() != omegas.len() {
            return Err(ModalLoadingError::Invalid(
                "one half-power half-width and one owner entry per mode",
            ));
        }
        for (&w, &b) in omegas.iter().zip(half_bands) {
            if !(w > 0.0 && w.is_finite()) {
                return Err(ModalLoadingError::Invalid(
                    "mode frequencies must be positive and finite",
                ));
            }
            if !(b >= 0.0 && b < w) {
                return Err(ModalLoadingError::Invalid(
                    "a half-power half-width must be finite, nonnegative and below its mode's frequency",
                ));
            }
        }
        Ok(Self {
            omegas: omegas.to_vec(),
            half_bands: half_bands.to_vec(),
            owners: owners.to_vec(),
        })
    }

    /// Number of modes.
    #[must_use]
    pub fn modes(&self) -> usize {
        self.omegas.len()
    }

    /// The pair-mean frequency of modes `r` and `s`, `(w_r + w_s) / 2`, the same bits for `(s, r)`.
    ///
    /// # Panics
    /// If `r` or `s` is not a mode (programmer error).
    #[must_use]
    pub fn pair_mean(&self, r: usize, s: usize) -> f64 {
        0.5 * (self.omegas[r] + self.omegas[s])
    }

    /// The edges `(w_r - b_r, w_r + b_r)` of mode `r`'s half-power band.
    ///
    /// # Panics
    /// If `r` is not a mode (programmer error).
    #[must_use]
    pub fn band_edges(&self, r: usize) -> (f64, f64) {
        (
            self.omegas[r] - self.half_bands[r],
            self.omegas[r] + self.half_bands[r],
        )
    }

    /// Whether one owner carries the radiation of both modes, so their entry is zero.
    ///
    /// # Panics
    /// If `r` or `s` is not a mode (programmer error).
    #[must_use]
    pub fn shared(&self, r: usize, s: usize) -> bool {
        matches!((self.owners[r], self.owners[s]), (Some(a), Some(b)) if a == b)
    }

    /// Every frequency a bake must declare for [`Self::sample`]: each pair's mean and each mode's
    /// band edges, ascending, each once.
    #[must_use]
    pub fn frequencies(&self) -> Vec<f64> {
        let n = self.modes();
        let mut all = Vec::with_capacity(n * (n + 1) / 2 + 2 * n);
        for r in 0..n {
            for s in r..n {
                all.push(self.pair_mean(r, s));
            }
            let (low, high) = self.band_edges(r);
            all.extend([low, high]);
        }
        all.sort_by(f64::total_cmp);
        all.dedup();
        all
    }

    /// The constant matrices and the band statement from a bake that declared every frequency of
    /// [`Self::frequencies`].
    ///
    /// # Errors
    /// [`ModalLoadingError::Invalid`] when the bake's mode count is not the law's or a frequency
    /// the law needs is not a declared sample.
    #[expect(
        clippy::result_large_err,
        reason = "ModalLoadingError is the module's refusal type (rule 1: no new error type); \
                  boxing it would change every modal_loading signature"
    )]
    pub fn sample(&self, bake: &ModalLoadingBake) -> Result<PairMeanLoading, ModalLoadingError> {
        let n = self.modes();
        if bake.modes() != n {
            return Err(ModalLoadingError::Invalid(
                "the bake does not hold the law's modes",
            ));
        }
        let mut added_mass = vec![0.0; n * n];
        let mut resistance = vec![0.0; n * n];
        for r in 0..n {
            for s in r..n {
                if self.shared(r, s) {
                    continue;
                }
                let (mass, real) = reciprocal(bake.at_omega(self.pair_mean(r, s))?, n, r, s);
                added_mass[r * n + s] = mass;
                added_mass[s * n + r] = mass;
                resistance[r * n + s] = real;
                resistance[s * n + r] = real;
            }
        }
        let mut band = Vec::with_capacity(n);
        for r in 0..n {
            let (low, high) = self.band_edges(r);
            let mut change = BandChange {
                mode: r,
                added_mass: 0.0,
                resistance: 0.0,
                largest_added_mass: 0.0,
                largest_resistance: 0.0,
            };
            for edge in [low, high] {
                let row = bake.at_omega(edge)?;
                for s in (0..n).filter(|&s| !self.shared(r, s)) {
                    let (mass, real) = reciprocal(row, n, r, s);
                    let (kept_mass, kept_real) = (added_mass[r * n + s], resistance[r * n + s]);
                    change.added_mass = change.added_mass.max((mass - kept_mass).abs());
                    change.resistance = change.resistance.max((real - kept_real).abs());
                    change.largest_added_mass = change.largest_added_mass.max(kept_mass.abs());
                    change.largest_resistance = change.largest_resistance.max(kept_real.abs());
                }
            }
            band.push(change);
        }
        Ok(PairMeanLoading {
            added_mass,
            resistance,
            band,
        })
    }
}

/// `(-Im / w, Re)` of the reciprocal part `(Z_rs + Z_sr) / 2` of one bake row.
fn reciprocal(row: &ModalLoadingRow, n: usize, r: usize, s: usize) -> (f64, f64) {
    let z = row.impedance();
    let (a, b) = (z[r * n + s], z[s * n + r]);
    (-(0.5 * (a.im + b.im)) / row.omega(), 0.5 * (a.re + b.re))
}
