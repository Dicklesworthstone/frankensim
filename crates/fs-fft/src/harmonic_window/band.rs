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
//! NIST DLMF 10.60.7 and 10.53.1 supply the expansion and power series:
//! <https://dlmf.nist.gov/10.60.E7>, <https://dlmf.nist.gov/10.53.E1>.

use super::{Dd, HarmonicWindow, HarmonicWindowProfile, WindowError, check_amplitude, checked};

/// Hard frequency-count bound on preparation and fixed observation scratch.
pub const MAX_UNIFORM_HARMONIC_BAND_TONES: usize = 64;
/// Maximum retained temporal degree plus one, including the constant mode.
pub const MAX_UNIFORM_HARMONIC_BAND_MODES: usize = 17;
const HALF: Dd = Dd { hi: 0.5, lo: 0.0 };

/// Cached orthogonal temporal factors for analytic fields using exp(-i*2*pi*f*t).
///
/// Source fields retain their original ordering. Equal frequencies should first
/// be coalesced by the caller. Independent source banks must be contracted
/// separately. There is no real-peak factor of one half: these methods return
/// analytic-field power/cross moments. Storage is O(m*n), m <= 17 and n <= 64;
/// observations allocate nothing. Each projection has at most 1088 terms, so
/// callers can poll cancellation before/after this bounded leaf.
pub struct PreparedUniformHarmonicBand {
    rotations: Vec<[Dd; 2]>,
    coefficients: Vec<Dd>,
    modes: usize,
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
        if frequencies_hz.len() > MAX_UNIFORM_HARMONIC_BAND_TONES {
            return Ok(None);
        }
        if frequencies_hz.iter().any(|f| !f.is_finite()) {
            return Err(WindowError::InvalidFrequency);
        }
        if self.profile != HarmonicWindowProfile::Uniform
            || !(3..=MAX_UNIFORM_HARMONIC_BAND_TONES).contains(&frequencies_hz.len())
        {
            return Ok(None);
        }
        let (minimum, maximum) = frequency_bounds(frequencies_hz);
        let span = checked(Dd::from_f64(maximum) - Dd::from_f64(minimum))?;
        let width = self.width_hz(span * HALF, 2.0)?;
        if HALF.lt(width) {
            return Ok(None);
        }
        // Keep at least n modes when possible, so short-band cancellation of
        // the first n-1 frequency moments still has a retained temporal mode.
        let minimum_modes = frequencies_hz.len().min(MAX_UNIFORM_HARMONIC_BAND_MODES);
        let radius = width.to_f64().next_up();
        let mut leading = 1.0;
        for modes in 1..=MAX_UNIFORM_HARMONIC_BAND_MODES {
            let odd = (2 * modes + 1) as f64;
            leading *= radius / odd;
            let tail = odd * leading * leading / (1.0 - radius * radius / (odd * (odd + 2.0)));
            // Leave margin for the floating estimate of the tail. At the
            // widest admitted band, 17 modes have tail < 4.15e-50 exactly.
            if modes >= minimum_modes && tail <= 4e-50 {
                return Ok(Some(modes));
            }
        }
        Ok(Some(MAX_UNIFORM_HARMONIC_BAND_MODES))
    }

    /// Prepare the admitted positive factor using original Hz and shutter clock.
    pub fn prepare_uniform_band_hz(
        self,
        frequencies_hz: &[f64],
    ) -> Result<Option<PreparedUniformHarmonicBand>, WindowError> {
        let Some(modes) = self.uniform_band_mode_count_hz(frequencies_hz)? else {
            return Ok(None);
        };
        let (minimum, maximum) = frequency_bounds(frequencies_hz);
        let center = (Dd::from_f64(maximum) - Dd::from_f64(minimum)) * HALF;
        let mut rotations = Vec::with_capacity(frequencies_hz.len());
        let mut coefficients = vec![Dd::ZERO; modes * frequencies_hz.len()];
        for (i, &frequency) in frequencies_hz.iter().enumerate() {
            let relative = (Dd::from_f64(frequency) - Dd::from_f64(minimum)) - center;
            rotations.push(self.rotation_hz(relative)?);
            let a = self.width_hz(relative, 2.0)?;
            let negative_square = -(a * a);
            let mut leading = Dd::ONE;
            for l in 0..modes {
                let mut series = Dd::ONE;
                let mut term = Dd::ONE;
                for k in 1..=16 {
                    term =
                        term * negative_square / Dd::from_f64((2 * k * (2 * l + 2 * k + 1)) as f64);
                    series = series + term;
                }
                coefficients[l * frequencies_hz.len() + i] = leading * series;
                leading = leading * a / Dd::from_f64((2 * l + 3) as f64);
            }
        }
        Ok(Some(PreparedUniformHarmonicBand {
            rotations,
            coefficients,
            modes,
        }))
    }
}

impl PreparedUniformHarmonicBand {
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
    pub fn project_parts_dd(
        &self,
        fields: &[[Dd; 2]],
        output: &mut [[Dd; 2]],
    ) -> Result<(), WindowError> {
        if fields.len() != self.rotations.len() || output.len() != self.modes {
            return Err(WindowError::InvalidBandShape);
        }
        let mut centered = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
        for (i, (&[ar, ai], &[cr, ci])) in fields.iter().zip(&self.rotations).enumerate() {
            check_amplitude(ar)?;
            check_amplitude(ai)?;
            centered[i] = [checked(ar * cr - ai * ci)?, checked(ar * ci + ai * cr)?];
        }
        for (mode, out) in output.iter_mut().enumerate() {
            *out = [Dd::ZERO; 2];
            for (i, field) in centered[..fields.len()].iter().enumerate() {
                let coefficient = self.coefficients[mode * fields.len() + i];
                out[0] = out[0] + coefficient * field[0];
                out[1] = out[1] + coefficient * field[1];
            }
            checked(out[0])?;
            checked(out[1])?;
        }
        Ok(())
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
        if mode >= self.modes {
            return Err(WindowError::InvalidBandShape);
        }
        for amplitude in left.into_iter().chain(right) {
            check_amplitude(amplitude)?;
        }
        weighted_cross(left, right, Dd::from_f64((2 * mode + 1) as f64))
    }

    /// Nonnegative uniform-window average of the squared analytic field.
    pub fn mean_square_parts_dd(&self, fields: &[[Dd; 2]]) -> Result<Dd, WindowError> {
        let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        self.project_parts_dd(fields, &mut modes[..self.modes])?;
        let mut power = Dd::ZERO;
        for (l, &field) in modes[..self.modes].iter().enumerate() {
            let weight = Dd::from_f64((2 * l + 1) as f64);
            let norm = weighted_cross(field, field, weight)?[0];
            power = checked(power + norm)?;
        }
        Ok(power)
    }

    /// Uniform-window complex cross moment <left(t)*conj(right(t))>.
    pub fn mean_cross_parts_dd(
        &self,
        left: &[[Dd; 2]],
        right: &[[Dd; 2]],
    ) -> Result<[Dd; 2], WindowError> {
        let mut a = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        let mut b = a;
        self.project_parts_dd(left, &mut a[..self.modes])?;
        self.project_parts_dd(right, &mut b[..self.modes])?;
        let mut result = [Dd::ZERO; 2];
        for l in 0..self.modes {
            let weight = Dd::from_f64((2 * l + 1) as f64);
            let cross = weighted_cross(a[l], b[l], weight)?;
            result[0] = checked(result[0] + cross[0])?;
            result[1] = checked(result[1] + cross[1])?;
        }
        Ok(result)
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
    // Scaling the smaller operand also avoids an avoidable overflow for
    // differently scaled fields in a complex cross moment.
    checked(if a.abs().lt(b.abs()) {
        (a * weight) * b
    } else {
        a * (b * weight)
    })
}

fn frequency_bounds(frequencies: &[f64]) -> (f64, f64) {
    frequencies
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &f| {
            (lo.min(f), hi.max(f))
        })
}

#[cfg(test)]
mod tests;
