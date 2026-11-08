//! Positive, fixed-size coherency factors for a calibrated incoherent source.
//!
//! C = sum_j w_j e_j e_j^dagger, with two independent Jones columns. Keeping
//! the columns between optical elements preserves a weak orthogonal population
//! even when the public Stokes I and Q round to the same binary64 number.
//! Columns never interfere with one another. This is an internal numerical
//! representation of a partially polarized source, not inferred source coherence.
//! Cxy = Ex conj(Ey), U = 2 Re(Cxy), V = 2 Im(Cxy), as in the native Mueller owner.

use super::{Mueller, Stokes};
use fs_math::{c64::C64, dd::Dd, det};

/// Principal field responses and their common transverse orientation.
#[derive(Clone, Copy, Debug)]
pub(super) struct JonesResponse {
    cosine: f64,
    sine: f64,
    channels: [C64; 2],
    power_scale: f64,
}

impl JonesResponse {
    pub(super) fn principal(angle: f64, channels: [C64; 2]) -> Self {
        let angle = angle.rem_euclid(core::f64::consts::PI);
        // Cardinal axes are exact: an authored crossed pair must not acquire
        // an artificial sin(pi) leakage floor above its finite extinction.
        let (sine, cosine) = if angle == 0.0 {
            (0.0, 1.0)
        } else if angle == core::f64::consts::FRAC_PI_2 {
            (1.0, 0.0)
        } else {
            det::sin_cos(angle)
        };
        Self { cosine, sine, channels, power_scale: 1.0 }
    }

    pub(super) fn neutral(power_scale: f64) -> Self {
        Self { cosine: 1.0, sine: 0.0, channels: [C64::ONE; 2], power_scale }
    }

    pub(super) fn from_diagonal_mueller(matrix: Mueller) -> Self {
        // A common field phase is irrelevant for an incoherent source. Choose
        // the stronger channel real and reconstruct the other from the retained
        // complex cross product, avoiding subtraction of rounded I/Q entries.
        let difference = matrix.0[0][1];
        let large = (matrix.0[0][0] + difference.abs()).sqrt();
        let cross = C64::new(matrix.0[2][2], matrix.0[3][2]);
        let channels = if large == 0.0 {
            [C64::ZERO; 2]
        } else if difference >= 0.0 {
            [C64::from_re(large), cross.conj().scale(1.0 / large)]
        } else {
            [cross.scale(1.0 / large), C64::from_re(large)]
        };
        Self::principal(0.0, channels)
    }

    fn apply(self, field: [C64; 2]) -> [C64; 2] {
        if self.channels == [C64::ONE; 2] { return field; }
        let [x, y] = field;
        if self.sine == 0.0 {
            return [self.channels[0] * x, self.channels[1] * y];
        }
        if self.cosine == 0.0 {
            return [self.channels[1] * x, self.channels[0] * y];
        }
        let parallel = (x.scale(self.cosine) + y.scale(self.sine)) * self.channels[0];
        let perpendicular = (y.scale(self.cosine) - x.scale(self.sine)) * self.channels[1];
        [parallel.scale(self.cosine) - perpendicular.scale(self.sine),
            parallel.scale(self.sine) + perpendicular.scale(self.cosine)]
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct SourceCoherency {
    columns: [[C64; 2]; 2],
    weights: [f64; 2],
}

impl SourceCoherency {
    pub(super) const UNPOLARIZED: Self = Self {
        columns: [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]],
        weights: [0.5; 2],
    };

    /// Pivoted LDL* factorization of an already admitted unit-intensity state.
    pub(super) fn from_stokes(state: Stokes) -> Self {
        // Only PolarizedSource::try_new calls this after validating I == 1.
        // Thus the pivot is at least 1/2 and never a dark-state 0/0.
        debug_assert_eq!(state.i, 1.0);
        let x = 0.5 * (state.i + state.q);
        let y = 0.5 * (state.i - state.q);
        let cross = C64::new(0.5 * state.u, 0.5 * state.v);
        let pivot = x.max(y);
        let u = Dd::from_f64(cross.re);
        let v = Dd::from_f64(cross.im);
        let residual = (Dd::from_f64(x.min(y)) - (u * u + v * v) / Dd::from_f64(pivot))
            .to_f64().max(0.0);
        let columns = if x >= y {
            [[C64::ONE, cross.conj().scale(1.0 / pivot)], [C64::ZERO, C64::ONE]]
        } else {
            [[cross.scale(1.0 / pivot), C64::ONE], [C64::ONE, C64::ZERO]]
        };
        Self { columns, weights: [pivot, residual] }
    }

    pub(super) fn scale_power(&mut self, factor: f64) {
        for weight in &mut self.weights { *weight *= factor; }
    }

    pub(super) fn through(mut self, response: JonesResponse) -> Self {
        for (column, weight) in self.columns.iter_mut().zip(&mut self.weights) {
            *weight *= response.power_scale;
            if *weight > 0.0 {
                *column = response.apply(*column);
            }
        }
        self
    }

    pub(super) fn stokes(self) -> Stokes {
        let mut axes = [0.0; 2];
        let mut cross = C64::ZERO;
        for ([x, y], weight) in self.columns.into_iter().zip(self.weights) {
            if weight > 0.0 {
                axes[0] += weight * x.norm_sq();
                axes[1] += weight * y.norm_sq();
                cross = cross + (x * y.conj()).scale(weight);
            }
        }
        Stokes { i: axes[0] + axes[1], q: axes[0] - axes[1],
            u: 2.0 * cross.re, v: 2.0 * cross.im }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_2;

    #[test]
    fn admitted_stokes_factorization_retains_mixed_and_circular_states() {
        for (q, u, v) in [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
            (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
            (0.0, 0.0, -1.0), (0.2, -0.4, 0.5)] {
            let state = SourceCoherency::from_stokes(Stokes { i: 1.0, q, u, v }).stokes();
            for (a, b) in [state.i, state.q, state.u, state.v].into_iter().zip([1.0, q, u, v]) {
                assert!((a - b).abs() < 5e-16, "{a} != {b}");
            }
        }
    }

    #[test]
    fn crossed_filters_recover_population_below_stokes_precision() {
        for leakage in [1e-20_f64, 1e-80, 1e-200] {
            let channels = [C64::ONE, C64::from_re(leakage.sqrt())];
            let first = SourceCoherency::UNPOLARIZED
                .through(JonesResponse::principal(0.0, channels));
            assert_eq!(first.stokes().i, first.stokes().q);
            let (s, c) = det::sin_cos(1.7);
            let output = first.through(JonesResponse::neutral(0.3))
                .through(JonesResponse::principal(0.0, [C64::ONE, C64::new(c, s)]))
                .through(JonesResponse::principal(FRAC_PI_2, channels)).stokes();
            assert!((output.i / (0.3 * leakage) - 1.0).abs() < 1e-14);
            assert!((output.q / output.i).abs() < 1e-14);
            assert_eq!(output.u, 0.0);
            assert_eq!(output.v, 0.0);
        }
    }

    #[test]
    fn diagonal_complex_response_retains_either_weak_population() {
        let strong = 0.8_f64;
        let weak = 1e-50_f64;
        let (s, c) = det::sin_cos(1.1);
        let cross = C64::new(c, s).scale(strong.sqrt() * weak.sqrt());
        for sign in [-1.0, 1.0] {
            let a = 0.5 * strong + 0.5 * weak;
            let b = sign * (0.5 * strong - 0.5 * weak);
            let matrix = Mueller([
                [a, b, 0.0, 0.0], [b, a, 0.0, 0.0],
                [0.0, 0.0, cross.re, -cross.im], [0.0, 0.0, cross.im, cross.re],
            ]);
            let angle = if sign > 0.0 { FRAC_PI_2 } else { 0.0 };
            let output = SourceCoherency::UNPOLARIZED
                .through(JonesResponse::from_diagonal_mueller(matrix))
                .through(JonesResponse::principal(angle, [C64::ONE, C64::ZERO])).stokes();
            assert!((output.i / (0.5 * weak) - 1.0).abs() < 1e-14);
            assert_eq!(output.q, -sign * output.i);
        }
    }
}
