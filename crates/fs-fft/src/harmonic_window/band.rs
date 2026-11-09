//! Positive uniform-window moments of a short coherent frequency band.
//!
//! On x in [-1,1], exp(-i*a*x) = sum_l (2*l+1)*(-i)^l*j_l(a)*P_l(x).
//! Legendre orthogonality therefore gives <|E|²> = sum_l (2*l+1)*|c_l|²,
//! where c_l = sum_i A_i*exp(-i*delta_omega_i*t_mid)*j_l(a_i). The common
//! (-i)^l cancels from every second moment, including complex cross moments.
//! This avoids subtracting bright pair powers to recover a weak coherent null.
//!
//! Only uniform windows with at most 64 tones and |a_i| <= 1/2 are admitted.
//! Up to 17 modes are retained. Since |j_l(a)| <= |a|^l/(2*l+1)!! here,
//! the exact-arithmetic power in omitted Legendre modes is below
//! 5e-50*(sum_i |A_i|)^2.
//! The small-argument series itself retains 17 terms. Floating arithmetic,
//! phase reduction and cancellation of individual modal coefficients are not
//! covered by that truncation bound; no null-relative certificate is claimed.
//!
//! Actual real-peak fields use their absolute source frequencies, with at most
//! 32 original tones. Conjugate parity evaluates the mirrored +/-frequency
//! expansion without halving amplitudes or dropping its sum-frequency terms.
//! Native Hz and rad/s sources retain their respective phase-clock owners.
//!
//! NIST DLMF 10.60.7 and 10.53.1 supply the expansion and power series:
//! <https://dlmf.nist.gov/10.60.E7>, <https://dlmf.nist.gov/10.53.E1>.

use super::{Dd, HarmonicWindow, HarmonicWindowProfile, WindowError, check_amplitude, checked};

mod basis;
mod workspace;
pub use workspace::{UniformHarmonicBandRef, UniformHarmonicBandWorkspace};

/// Hard frequency-count bound on preparation and fixed observation scratch.
pub const MAX_UNIFORM_HARMONIC_BAND_TONES: usize = 64;
/// Maximum retained temporal degree plus one, including the constant mode.
pub const MAX_UNIFORM_HARMONIC_BAND_MODES: usize = 17;
const HALF: Dd = Dd { hi: 0.5, lo: 0.0 };

/// Cached orthogonal temporal factors using the negative-time convention.
///
/// Source fields retain their original ordering. Equal frequencies should first
/// be coalesced by the caller. Independent source banks must be contracted
/// separately. Analytic constructors return analytic-field power/cross moments;
/// real constructors return actual real-peak products, including sum tones.
/// Proper complex peak factors require the explicit scaled contraction with
/// scale 0.5. Storage is O(m*n), m <= 17 and n <= 64;
/// observations allocate nothing. Each projection has at most 1088 terms, so
/// callers can poll cancellation before/after this bounded leaf.
pub struct PreparedUniformHarmonicBand {
    rotations: Vec<[Dd; 2]>,
    coefficients: Vec<Dd>,
    modes: usize,
    real_field: bool,
}

impl HarmonicWindow {
    /// Allocation-free admission and mode count for the positive band factor.
    ///
    /// `None` selects the original pair kernels: fewer than three or more than
    /// 64 tones, a nonuniform window, or a half-band phase greater than 1/2.
    /// This lets consumers admit projection work BEFORE allocating the factor.
    pub fn uniform_band_mode_count_hz(
        self,
        frequencies_hz: &[f64],
    ) -> Result<Option<usize>, WindowError> {
        basis::mode_count(self, frequencies_hz, basis::Units::Hertz, false)
    }

    /// Prepare the admitted positive factor using original Hz and shutter clock.
    pub fn prepare_uniform_band_hz(
        self,
        frequencies_hz: &[f64],
    ) -> Result<Option<PreparedUniformHarmonicBand>, WindowError> {
        basis::prepare_owned(self, frequencies_hz, basis::Units::Hertz, false)
    }
}

impl PreparedUniformHarmonicBand {
    fn view(&self) -> UniformHarmonicBandRef<'_> {
        UniformHarmonicBandRef {
            rotations: &self.rotations,
            coefficients: &self.coefficients,
            modes: self.modes,
            real_field: self.real_field,
        }
    }

    /// Number of original source frequencies, without a mirrored expansion.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.rotations.len()
    }

    /// Number of retained orthogonal temporal modes.
    #[must_use]
    pub const fn mode_count(&self) -> usize {
        self.modes
    }

    /// Logical complex projection terms plus one squared norm per mode.
    #[must_use]
    pub fn power_work(&self) -> usize {
        (self.rotations.len() + 1) * self.modes
    }

    /// Two field projections plus one complex cross moment per mode.
    #[must_use]
    pub fn cross_work(&self) -> usize {
        (2 * self.rotations.len() + 1) * self.modes
    }

    /// Two projections plus the four real entries of a 2x2 coherency per mode.
    #[must_use]
    pub fn coherency_work(&self) -> usize {
        (2 * self.rotations.len() + 4) * self.modes
    }

    /// Project retained complex fields into the unphased Legendre coefficients.
    /// `output` must have exactly `mode_count()` entries. Only second moments
    /// may omit each mode's common (-i)^l phase. Scratch can change on refusal.
    /// In real preparation, mode zero is also the actual real field mean.
    pub fn project_parts_dd(
        &self,
        fields: &[[Dd; 2]],
        output: &mut [[Dd; 2]],
    ) -> Result<(), WindowError> {
        self.view().project_parts_dd(fields, output)
    }

    /// Weighted complex outer product of one pair of projected coefficients.
    ///
    /// Returns `(2*mode+1)*left*conj(right)` for a retained mode. Inputs are
    /// unweighted coefficients from `project_parts_dd`; do not apply the modal
    /// weight again. There is no real-peak factor of one half. This shares the
    /// scalar/cross contraction's product ordering with Jones coherency users.
    pub fn mode_cross_parts_dd(
        &self,
        mode: usize,
        left: [Dd; 2],
        right: [Dd; 2],
    ) -> Result<[Dd; 2], WindowError> {
        self.view().mode_cross_parts_dd(mode, left, right)
    }

    /// Apply a positive scale before a modal cross product can overflow.
    /// Use scale 0.5 for analytic peak factors with proper complex statistics.
    pub fn mode_cross_parts_dd_scaled(
        &self,
        mode: usize,
        left: [Dd; 2],
        right: [Dd; 2],
        positive_scale: Dd,
    ) -> Result<[Dd; 2], WindowError> {
        self.view()
            .mode_cross_parts_dd_scaled(mode, left, right, positive_scale)
    }

    /// Nonnegative uniform-window average of the prepared field's squared norm.
    pub fn mean_square_parts_dd(&self, fields: &[[Dd; 2]]) -> Result<Dd, WindowError> {
        self.view().mean_square_parts_dd(fields)
    }

    /// Scaled squared norm, applying the positive scale before all products.
    pub fn mean_square_parts_dd_scaled(
        &self,
        fields: &[[Dd; 2]],
        positive_scale: Dd,
    ) -> Result<Dd, WindowError> {
        self.view()
            .mean_square_parts_dd_scaled(fields, positive_scale)
    }

    /// Uniform-window complex cross moment <left(t)*conj(right(t))>.
    pub fn mean_cross_parts_dd(
        &self,
        left: &[[Dd; 2]],
        right: &[[Dd; 2]],
    ) -> Result<[Dd; 2], WindowError> {
        self.view().mean_cross_parts_dd(left, right)
    }

    /// Scaled complex cross moment with the positive scale inside each product.
    pub fn mean_cross_parts_dd_scaled(
        &self,
        left: &[[Dd; 2]],
        right: &[[Dd; 2]],
        positive_scale: Dd,
    ) -> Result<[Dd; 2], WindowError> {
        self.view()
            .mean_cross_parts_dd_scaled(left, right, positive_scale)
    }
}

fn weighted_cross(left: [Dd; 2], right: [Dd; 2], weight: Dd) -> Result<[Dd; 2], WindowError> {
    let real = checked(
        weighted_product(left[0], right[0], weight)? + weighted_product(left[1], right[1], weight)?,
    )?;
    // A norm needs two scalar squares, and its imaginary part is exactly zero.
    // Sharing this path preserves identical scalar and Jones diagonal moments.
    let imaginary = if left == right {
        Dd::ZERO
    } else {
        checked(
            weighted_product(left[1], right[0], weight)?
                - weighted_product(left[0], right[1], weight)?,
        )?
    };
    Ok([real, imaginary])
}

fn weighted_product(a: Dd, b: Dd, weight: Dd) -> Result<Dd, WindowError> {
    // Apply the small positive modal weight before a product can underflow.
    // For a weight below one (notably the proper peak factor), scale the
    // larger operand instead so a subnormal small operand is not erased.
    let scale_left = if weight.lt(Dd::ONE) {
        b.abs().lt(a.abs())
    } else {
        a.abs().lt(b.abs())
    };
    checked(if scale_left {
        (a * weight) * b
    } else {
        a * (b * weight)
    })
}

#[cfg(test)]
mod tests;
