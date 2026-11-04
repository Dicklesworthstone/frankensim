//! Exact-length DFT plans built on the existing Stockham transform.
//!
//! Bluestein's identity factors exp(-2 pi i j k / N) into two chirps and
//! their convolution. Padding that convolution to M >= 2*N-1 does not pad
//! the physical record: the returned N coefficients still have spacing
//! 1/(N*dt). See Rabiner, Schafer and Rader, *The Chirp z-Transform Algorithm*,
//! IEEE Trans. Audio Electroacoustics 17(2), 86--92 (1969),
//! https://web.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/015_czt.pdf.

use crate::{C64, Fft};
use fs_math::det;

/// Explicit record-size and numeric-storage limits for [`GeneralFft`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneralFftBudget {
    /// Largest admitted physical record. Zero admits no transform.
    pub max_samples: usize,
    /// Combined plan and minimum caller-scratch payload, conservatively
    /// including preparation temporaries and existing six-step row buffers.
    /// The transformed data, extra caller scratch, stack and allocator
    /// metadata are excluded. Admission precedes numeric allocation.
    pub max_workspace_bytes: usize,
}

/// Refusal to prepare or evaluate an exact-length DFT.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneralFftError {
    /// N must be positive and exactly representable as an f64 integer.
    InvalidSize { samples: usize },
    /// N exceeds the caller's physical-record limit.
    SampleLimit { samples: usize, maximum: usize },
    /// The convolution size or addressable numeric payload overflowed.
    SizeOverflow,
    /// Required plan/scratch payload exceeds the explicit storage limit.
    WorkspaceBudget {
        required_bytes: usize,
        maximum_bytes: usize,
    },
    /// Data must have exactly N entries and scratch at least scratch_len().
    BufferLength {
        data: usize,
        scratch: usize,
        expected_data: usize,
        minimum_scratch: usize,
    },
    /// An input or output component cannot represent a finite value.
    NonFinite { what: &'static str },
}

impl core::fmt::Display for GeneralFftError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidSize { samples } => {
                write!(f, "FFT size {samples} must be in 1..=2^53")
            }
            Self::SampleLimit { samples, maximum } => {
                write!(f, "FFT record has {samples} samples, limit is {maximum}")
            }
            Self::SizeOverflow => write!(f, "FFT workspace size is unrepresentable"),
            Self::WorkspaceBudget {
                required_bytes,
                maximum_bytes,
            } => write!(
                f,
                "FFT plan and scratch require {required_bytes} bytes, budget is {maximum_bytes}"
            ),
            Self::BufferLength {
                data,
                scratch,
                expected_data,
                minimum_scratch,
            } => write!(
                f,
                "FFT buffers have lengths {data}/{scratch}; need {expected_data}/at least {minimum_scratch}"
            ),
            Self::NonFinite { what } => write!(f, "non-finite FFT {what}"),
        }
    }
}

impl std::error::Error for GeneralFftError {}

#[derive(Debug, Clone)]
enum Route {
    Direct,
    Bluestein {
        chirp: Vec<C64>,
        kernel_spectrum: Vec<C64>,
    },
}

/// Reusable DFT of any admitted positive length, including odd and prime N.
///
/// Power-of-two records dispatch to [`Fft`] with unchanged transform bits.
/// Other lengths use a cached Bluestein convolution of length
/// M = next_power_of_two(2*N-1). Preparation and execution take O(M log M)
/// work and O(M) storage; each execution reuses the chirp, kernel spectrum,
/// underlying FFT plan and caller scratch. No waveform interpolation,
/// truncation or change to the physical sample clock is performed.
///
/// Default transforms allocate nothing after preparation. When the existing
/// `frontier-sixstep` feature selects six-step transforms, its transient row
/// buffers remain in use and are included in workspace_bytes(). As with
/// [`Fft`], a transform is synchronous, without mid-transform cancellation.
/// Explicit sample and storage limits bound admitted work; this is not a
/// deadline, interval-error or arbitrary-null relative-accuracy certificate.
#[derive(Debug, Clone)]
pub struct GeneralFft {
    samples: usize,
    fft: Fft,
    route: Route,
    scratch_len: usize,
    workspace_bytes: usize,
}

impl GeneralFft {
    /// Plan a transform without altering the requested physical record length.
    pub fn new(samples: usize, budget: GeneralFftBudget) -> Result<Self, GeneralFftError> {
        if samples == 0 || samples as u128 > (1_u128 << 53) {
            return Err(GeneralFftError::InvalidSize { samples });
        }
        if samples > budget.max_samples {
            return Err(GeneralFftError::SampleLimit {
                samples,
                maximum: budget.max_samples,
            });
        }
        let layout = layout(samples)?;
        if layout.workspace_bytes > budget.max_workspace_bytes {
            return Err(GeneralFftError::WorkspaceBudget {
                required_bytes: layout.workspace_bytes,
                maximum_bytes: budget.max_workspace_bytes,
            });
        }
        let fft = Fft::new(layout.transform_len);
        let route = if samples.is_power_of_two() {
            Route::Direct
        } else {
            let chirp: Vec<C64> = (0..samples).map(|j| forward_chirp(j, samples)).collect();
            let mut kernel_spectrum = vec![C64::default(); layout.transform_len];
            kernel_spectrum[0] = C64::new(1.0, 0.0);
            for (j, &value) in chirp.iter().enumerate().skip(1) {
                // The negative convolution lags wrap at M, not at N. This
                // distinction is essential for the odd-length DFT.
                kernel_spectrum[j] = value.conj();
                kernel_spectrum[layout.transform_len - j] = value.conj();
            }
            let mut scratch = vec![C64::default(); layout.transform_len];
            fft.forward(&mut kernel_spectrum, &mut scratch);
            Route::Bluestein {
                chirp,
                kernel_spectrum,
            }
        };
        Ok(Self {
            samples,
            fft,
            route,
            scratch_len: layout.scratch_len,
            workspace_bytes: layout.workspace_bytes,
        })
    }

    /// Physical record length, equal to both input and output lengths.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.samples
    }

    /// Plans always contain at least one physical sample.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        false
    }

    /// Minimum reusable scratch length in complex values. Extra entries are
    /// accepted and never read or modified.
    #[must_use]
    pub const fn scratch_len(&self) -> usize {
        self.scratch_len
    }

    /// Numeric plan/scratch payload admitted at preparation.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Determine the required numeric payload without allocating a plan.
    pub fn required_workspace_bytes(samples: usize) -> Result<usize, GeneralFftError> {
        Ok(layout(samples)?.workspace_bytes)
    }

    /// Forward DFT, X[k] = sum_j x[j] exp(-2*pi*i*j*k/N), unnormalized.
    ///
    /// Buffer and input refusals leave both buffers unchanged. A numerical
    /// overflow leaves scratch unspecified; the direct power-of-two path may
    /// also have modified data, matching the underlying in-place transform.
    pub fn forward(&self, data: &mut [C64], scratch: &mut [C64]) -> Result<(), GeneralFftError> {
        self.transform(data, scratch, false)
    }

    /// Inverse DFT with 1/N normalization. Refusal semantics match forward().
    pub fn inverse(&self, data: &mut [C64], scratch: &mut [C64]) -> Result<(), GeneralFftError> {
        self.transform(data, scratch, true)
    }

    fn transform(
        &self,
        data: &mut [C64],
        scratch: &mut [C64],
        inverse: bool,
    ) -> Result<(), GeneralFftError> {
        if data.len() != self.samples || scratch.len() < self.scratch_len {
            return Err(GeneralFftError::BufferLength {
                data: data.len(),
                scratch: scratch.len(),
                expected_data: self.samples,
                minimum_scratch: self.scratch_len,
            });
        }
        if data.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err(GeneralFftError::NonFinite { what: "input" });
        }
        match &self.route {
            Route::Direct => {
                if inverse {
                    self.fft.inverse(data, &mut scratch[..self.samples]);
                } else {
                    self.fft.forward(data, &mut scratch[..self.samples]);
                }
                if data.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
                    return Err(GeneralFftError::NonFinite { what: "output" });
                }
            }
            Route::Bluestein {
                chirp,
                kernel_spectrum,
            } => {
                let (work, fft_scratch) = scratch[..self.scratch_len].split_at_mut(self.fft.len());
                // Normalize before the convolution's larger intermediate
                // spectrum, then restore physical units. This keeps tiny
                // representable records and avoids avoidable intermediate
                // overflow for large finite amplitudes. Relative components
                // beyond f64 dynamic range are not recovered by scaling.
                let scale = data
                    .iter()
                    .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
                if scale == 0.0 {
                    data.fill(C64::default());
                    return Ok(());
                }
                work.fill(C64::default());
                for ((out, &value), &rotation) in work.iter_mut().zip(data.iter()).zip(chirp) {
                    let value = C64::new(value.re / scale, value.im / scale);
                    *out = if inverse { value.conj() } else { value }.mul(rotation);
                }
                self.fft.forward(work, fft_scratch);
                for (value, &kernel) in work.iter_mut().zip(kernel_spectrum) {
                    *value = value.mul(kernel);
                }
                self.fft.inverse(work, fft_scratch);
                for (value, &rotation) in work.iter_mut().zip(chirp) {
                    *value = value.mul(rotation);
                    if inverse {
                        *value = value.conj().scale(1.0 / self.samples as f64);
                    }
                    *value = value.scale(scale);
                    if !value.re.is_finite() || !value.im.is_finite() {
                        return Err(GeneralFftError::NonFinite { what: "output" });
                    }
                }
                data.copy_from_slice(&work[..self.samples]);
            }
        }
        Ok(())
    }
}

struct Layout {
    transform_len: usize,
    scratch_len: usize,
    workspace_bytes: usize,
}

fn layout(samples: usize) -> Result<Layout, GeneralFftError> {
    if samples == 0 || samples as u128 > (1_u128 << 53) {
        return Err(GeneralFftError::InvalidSize { samples });
    }
    let overflow = GeneralFftError::SizeOverflow;
    let direct = samples.is_power_of_two();
    let transform_len = if direct {
        samples
    } else {
        samples
            .checked_mul(2)
            .and_then(|n| n.checked_sub(1))
            .and_then(usize::checked_next_power_of_two)
            .ok_or_else(|| overflow.clone())?
    };
    let scratch_len = transform_len
        .checked_mul(if direct { 1 } else { 2 })
        .ok_or_else(|| overflow.clone())?;
    let retained = if direct {
        0
    } else {
        samples
            .checked_add(transform_len)
            .ok_or_else(|| overflow.clone())?
    };
    let mut entries = retained
        .checked_add(scratch_len)
        .ok_or_else(|| overflow.clone())?;
    let mut n = transform_len;
    loop {
        entries = entries.checked_add(n / 2).ok_or_else(|| overflow.clone())?;
        if !Fft::takes_sixstep(n) {
            break;
        }
        n = 1 << (n.trailing_zeros() / 2);
        // Existing six-step recursion retains one row-scratch plus eight
        // gathered rows per active level. Its sub-plan tables are counted
        // by the next loop iteration.
        entries = n
            .checked_mul(9)
            .and_then(|rows| entries.checked_add(rows))
            .ok_or_else(|| overflow.clone())?;
    }
    let workspace_bytes = entries
        .checked_mul(core::mem::size_of::<C64>())
        .filter(|&bytes| bytes <= isize::MAX as usize)
        .ok_or(overflow)?;
    Ok(Layout {
        transform_len,
        scratch_len,
        workspace_bytes,
    })
}

// exp(-i*pi*j*j/N). Integer reduction happens before any conversion to f64,
// so large j cannot destroy the small phase remaining after whole turns.
// Reduce again to the nearest exact quarter turn: weak real/imaginary parts
// near +/-1 or +/-i are evaluated at their small residual angle directly.
fn forward_chirp(index: usize, samples: usize) -> C64 {
    let n = samples as u128;
    let j = index as u128;
    let residue = (j * j) % (2 * n);
    let quadrant = (4 * residue + n) / (2 * n);
    let remainder = (2 * residue) as i128 - (quadrant * n) as i128;
    let theta = -core::f64::consts::PI * (remainder as f64 / (2.0 * samples as f64));
    let (sin, cos) = det::sin_cos(theta);
    match quadrant % 4 {
        0 => C64::new(cos, sin),
        1 => C64::new(sin, -cos),
        2 => C64::new(-cos, -sin),
        _ => C64::new(-sin, cos),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modular_chirp_retains_small_phase_at_large_indices() {
        for n in [1_000_000_006_usize, 1_000_000_007] {
            let tail = forward_chirp(n - 1, n);
            let theta = core::f64::consts::PI / n as f64;
            let (sin, cos) = det::sin_cos(theta);
            let expected = if n % 2 == 0 {
                C64::new(cos, -sin)
            } else {
                C64::new(-cos, sin)
            };
            assert_eq!(tail, expected);
        }
    }
}
