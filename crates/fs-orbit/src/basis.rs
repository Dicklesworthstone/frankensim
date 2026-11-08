//! Exact integer phase indexing for the uniform AFT quadrature.
//!
//! exp(i * n * 2*pi*j/M) is the (n*j mod M)-th root of unity. Store
//! those M roots once instead of M*N values for N harmonics. Advancing an
//! integer index introduces no accumulated oscillator phase or amplitude drift.

use super::{OrbitError, TAU};
use fs_math::{c64::C64, det};

pub(super) struct PhaseBasis {
    roots: Vec<C64>,
}

impl PhaseBasis {
    pub(super) fn new(samples: usize) -> Result<Self, OrbitError> {
        let invalid = || OrbitError::BadParameter {
            what: "AFT phase grid must have a representable nonempty payload",
        };
        if samples == 0
            || samples as u128 > (1_u128 << 53)
            || samples
                .checked_mul(core::mem::size_of::<C64>())
                .is_none_or(|bytes| bytes > isize::MAX as usize)
        {
            return Err(invalid());
        }
        let mut roots = Vec::new();
        roots.try_reserve_exact(samples).map_err(|_| invalid())?;
        roots.resize(samples, C64::ZERO);
        roots[0] = C64::from_re(1.0);
        for j in 1..=samples / 2 {
            let value = if 2 * j == samples {
                C64::from_re(-1.0)
            } else if samples.is_multiple_of(4) && j == samples / 4 {
                C64::new(0.0, 1.0)
            } else {
                let (sine, cosine) = det::sin_cos(TAU * j as f64 / samples as f64);
                C64::new(cosine, sine)
            };
            roots[j] = value;
            roots[samples - j] = value.conj();
        }
        Ok(Self { roots })
    }

    /// At fixed sample j, cycle(j) visits harmonics n=0,1,... . At fixed
    /// harmonic n, cycle(n) visits samples j=0,1,... . The same roots serve
    /// reconstruction, projection, adaptive checks and exact state tangents.
    pub(super) fn cycle(&self, stride: usize) -> Phases<'_> {
        debug_assert!(stride < self.roots.len());
        Phases {
            roots: &self.roots,
            stride,
            index: 0,
        }
    }
}

pub(super) struct Phases<'a> {
    roots: &'a [C64],
    stride: usize,
    index: usize,
}

impl Iterator for Phases<'_> {
    type Item = C64;

    #[inline]
    fn next(&mut self) -> Option<C64> {
        let value = self.roots[self.index];
        // Modular addition without division, overflowing products or a float
        // recurrence. This also handles DC and non-coprime odd-grid strides.
        let wrap = self.roots.len() - self.stride;
        self.index = if self.index >= wrap {
            self.index - wrap
        } else {
            self.index + self.stride
        };
        Some(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn even_and_odd_grids_match_independent_trigonometric_phases() {
        for samples in [1, 2, 3, 16, 33, 65, 257] {
            let basis = PhaseBasis::new(samples).unwrap();
            for stride in 0..samples {
                for (j, phase) in basis.cycle(stride).take(3 * samples).enumerate() {
                    let angle = TAU * ((j * stride) % samples) as f64 / samples as f64;
                    // The full-angle independent reference rounds its angle
                    // differently from the conjugate half-grid construction.
                    let tolerance = 8.0 * f64::EPSILON;
                    assert!((phase.re - angle.cos()).abs() < tolerance);
                    assert!((phase.im - angle.sin()).abs() < tolerance);
                }
            }
            assert_eq!(basis.roots.len(), samples);
        }
    }

    #[test]
    fn phase_wrap_and_conjugate_symmetry_have_no_accumulated_drift() {
        for samples in [32, 65, 129] {
            let basis = PhaseBasis::new(samples).unwrap();
            for j in 1..samples {
                assert_eq!(basis.roots[j], basis.roots[samples - j].conj());
            }
            for stride in [0, 1, samples - 1] {
                let initial: Vec<_> = basis.cycle(stride).take(samples).collect();
                let later: Vec<_> = basis
                    .cycle(stride)
                    .skip(1000 * samples)
                    .take(samples)
                    .collect();
                assert_eq!(initial, later);
            }
        }
        assert!(PhaseBasis::new(0).is_err());
        assert!(PhaseBasis::new(usize::MAX).is_err());
    }
}
