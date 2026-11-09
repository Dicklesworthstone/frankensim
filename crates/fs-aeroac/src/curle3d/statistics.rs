//! Second moments of a coherently assembled three-dimensional acoustic field.

use crate::AeroacError;
use fs_math::{c64::C64, dd::Dd};

/// Work bounds checked before propagating compact-source coherence factors.
#[derive(Clone, Copy, Debug)]
pub struct CoherenceBudget {
    /// Maximum number of independent factors, including empty factors.
    pub max_modes: usize,
    /// Maximum source-path/factor products; overflow refuses. Free-field
    /// observers count one path per source, ideal half-space observers two.
    /// Orthogonal-plane observers fold their image kernels during preparation,
    /// so their application counts one cached transfer per physical source.
    pub max_source_mode_pairs: usize,
}

/// Physical spectral-line second moments, with the peak-phasor time average 1/2.
///
/// These quantities retain cross-source and cross-mechanism interference within
/// each supplied coherence factor. They are not PSD/Hz or a stochastic confidence
/// interval. No source correlation is inferred from material or distance.
#[derive(Clone, Debug)]
pub struct CompactFieldStatistics {
    /// E[|p|²]/2 [Pa²]. Its square root is pressure RMS, not energy flux.
    pub mean_square_pressure_pa2: f64,
    /// E[u_i conj(u_j)]/2 [m²/s²], Hermitian with real component variances.
    pub velocity_covariance_m2_s2: [[C64; 3]; 3],
    /// Re(E[p conj(u)])/2 [W/m²], with signed Cartesian components.
    pub active_intensity_w_m2: [f64; 3],
    /// Im(E[p conj(u)])/2 [W/m²], in exp(-i omega t) convention.
    pub reactive_intensity_w_m2: [f64; 3],
    /// Prescribed-source and propagation-model boundaries.
    pub scope: &'static str,
}

#[derive(Default)]
pub(super) struct StatisticsSum {
    mean_squares: [Dd; 4],
    // Upper triangle: p ux*, p uy*, p uz*, ux uy*, ux uz*, uy uz*.
    cross: [ComplexSum; 6],
}

impl StatisticsSum {
    pub(super) fn add(&mut self, field: [C64; 4]) -> Result<(), AeroacError> {
        for (sum, value) in self.mean_squares.iter_mut().zip(field) {
            let square = (0.5 * value.re) * value.re + (0.5 * value.im) * value.im;
            if !square.is_finite() {
                return Err(nonfinite("compact field second moment"));
            }
            *sum = *sum + Dd::from_f64(square);
        }
        let mut entry = 0;
        for i in 0..4 {
            for j in i + 1..4 {
                self.cross[entry].add_product(field[i].scale(0.5), field[j].conj())?;
                entry += 1;
            }
        }
        Ok(())
    }

    pub(super) fn finish(self, scope: &'static str) -> Result<CompactFieldStatistics, AeroacError> {
        let squares = self.mean_squares.map(Dd::to_f64);
        let cross = self.cross.map(ComplexSum::value);
        if squares.iter().any(|value| !value.is_finite())
            || cross.iter().any(|value| !finite(*value))
        {
            return Err(nonfinite("accumulated compact field statistics"));
        }
        Ok(CompactFieldStatistics {
            mean_square_pressure_pa2: squares[0],
            velocity_covariance_m2_s2: [
                [C64::from_re(squares[1]), cross[3], cross[4]],
                [cross[3].conj(), C64::from_re(squares[2]), cross[5]],
                [cross[4].conj(), cross[5].conj(), C64::from_re(squares[3])],
            ],
            active_intensity_w_m2: [cross[0].re, cross[1].re, cross[2].re],
            reactive_intensity_w_m2: [cross[0].im, cross[1].im, cross[2].im],
            scope,
        })
    }
}

#[derive(Clone, Copy, Default)]
pub(super) struct ComplexSum {
    re: Dd,
    im: Dd,
}

impl ComplexSum {
    pub(super) fn add(&mut self, value: C64) -> Result<(), AeroacError> {
        if !finite(value) {
            return Err(nonfinite("compact source contribution"));
        }
        self.re = self.re + Dd::from_f64(value.re);
        self.im = self.im + Dd::from_f64(value.im);
        Ok(())
    }

    /// Preserve the product's low parts until all coherent sources have been
    /// combined. Summing already-rounded products loses a weak residual even
    /// when the subsequent addition itself uses double-double arithmetic.
    pub(super) fn add_product(&mut self, left: C64, right: C64) -> Result<(), AeroacError> {
        let product = |a, b| {
            let (hi, lo) = fs_math::eft::two_prod(a, b);
            Dd { hi, lo }
        };
        let real = product(left.re, right.re) - product(left.im, right.im);
        let imaginary = product(left.re, right.im) + product(left.im, right.re);
        if [real.hi, real.lo, imaginary.hi, imaginary.lo]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(nonfinite("compact source product"));
        }
        self.re = self.re + real;
        self.im = self.im + imaginary;
        Ok(())
    }

    pub(super) fn value(self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
    }
}

pub(super) fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

pub(super) fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

pub(super) fn nonfinite(what: &'static str) -> AeroacError {
    AeroacError::NonFinite { what }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// G0/G3: independent factors may transport nearly opposing signed work.
    /// Peak-phasor averaging must keep the small negative flux and its 1/2.
    #[test]
    fn independent_factors_retain_weak_signed_intensity() {
        let epsilon = 2.0_f64.powi(-27);
        let mut moments = StatisticsSum::default();
        moments
            .add([
                C64::from_re(1.0 + epsilon),
                C64::from_re(1.0 - epsilon),
                C64::ZERO,
                C64::ZERO,
            ])
            .unwrap();
        moments
            .add([C64::ONE, -C64::ONE, C64::ZERO, C64::ZERO])
            .unwrap();
        let result = moments.finish("algebraic field factors").unwrap();
        assert_eq!(result.active_intensity_w_m2, [-0.5 * epsilon * epsilon, 0.0, 0.0]);
        assert_eq!(result.reactive_intensity_w_m2, [0.0; 3]);
        assert!(result.mean_square_pressure_pa2 > 0.0);
        assert!(result.velocity_covariance_m2_s2[0][0].re > 0.0);
    }
}
