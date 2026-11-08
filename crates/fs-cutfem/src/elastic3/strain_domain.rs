//! Allocation-free Frobenius admission of an affine Mandel-strain family.
//!
//! In Mandel coordinates the tensor Frobenius norm is the Euclidean norm.
//! Normalization by the original strain ceiling bounds every arithmetic input.
//! A point or one-dimensional family is checked exactly in real arithmetic;
//! higher-rank families use a rotation-invariant sufficient upper bound.

pub(super) struct FrobeniusStrainBall {
    offset: [f64; 6],
    offset_squared: f64,
    gram: [[f64; 6]; 6],
    coupling_squared: f64,
    directions: usize,
    limit: f64,
}

impl FrobeniusStrainBall {
    pub(super) fn new(offset: [f64; 6], limit: f64) -> Option<Self> {
        let offset = normalized(offset, limit)?;
        let offset_squared = dot(&offset, &offset);
        (offset_squared <= 1.0).then_some(Self {
            offset,
            offset_squared,
            gram: [[0.0; 6]; 6],
            coupling_squared: 0.0,
            directions: 0,
            limit,
        })
    }

    pub(super) fn add(&mut self, direction: [f64; 6]) -> bool {
        let Some(direction) = normalized(direction, self.limit) else {
            return false;
        };
        // max(||e+d||, ||e-d||) >= ||d||; a direction outside the ball
        // cannot be rescued by an offset of the opposite sign.
        if dot(&direction, &direction) > 1.0 {
            return false;
        }
        for i in 0..6 {
            for j in 0..=i {
                self.gram[i][j] += direction[i] * direction[j];
            }
        }
        let coupling = dot(&self.offset, &direction);
        self.coupling_squared += coupling * coupling;
        self.directions += 1;
        true
    }

    pub(super) fn admitted(&self) -> bool {
        if self.directions == 0 {
            return true;
        }
        // For A=B B^T with rank <= r, its largest eigenvalue is bounded
        // by tr(A)/r + sqrt((r-1)/r * (tr(A^2)-tr(A)^2/r)). This is exact
        // for up to two directions (or rank one) and avoids an eigensolve
        // at every quadrature/vertex; rank alone does not set r here.
        // The full affine support obeys
        // max ||e+B a||^2 <= ||e||^2 + lambda_max(A) + 2||B^T e||.
        let rank = self.directions.min(6) as f64;
        let trace: f64 = (0..6).map(|i| self.gram[i][i]).sum();
        let mut trace_squared = 0.0;
        for i in 0..6 {
            trace_squared += self.gram[i][i] * self.gram[i][i];
            for j in 0..i {
                trace_squared += 2.0 * self.gram[i][j] * self.gram[i][j];
            }
        }
        let mean = trace / rank;
        let variance = (trace_squared - trace * mean).max(0.0);
        let spectral_bound = (mean + ((rank - 1.0) / rank * variance).sqrt()).min(trace);
        let upper_squared =
            self.offset_squared + spectral_bound + 2.0 * self.coupling_squared.sqrt();
        upper_squared.is_finite() && upper_squared <= 1.0
    }
}

fn normalized(strain: [f64; 6], limit: f64) -> Option<[f64; 6]> {
    if !limit.is_finite()
        || limit <= 0.0
        || strain.iter().any(|v| !v.is_finite() || v.abs() > limit)
    {
        return None;
    }
    Some(strain.map(|v| v / limit))
}

fn dot(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[cfg(test)]
mod tests {
    use super::FrobeniusStrainBall;

    #[test]
    fn full_tensor_norm_and_extreme_representable_limits() {
        for limit in [1e-200, 0.1, 1e200] {
            assert!(
                FrobeniusStrainBall::new([0.7 * limit, 0.7 * limit, 0.0, 0.0, 0.0, 0.0], limit)
                    .is_some()
            );
            assert!(
                FrobeniusStrainBall::new([0.8 * limit, 0.8 * limit, 0.0, 0.0, 0.0, 0.0], limit)
                    .is_none()
            );
        }
    }

    #[test]
    fn perpendicular_offset_and_direction_use_tensor_geometry() {
        let mut ball = FrobeniusStrainBall::new([0.06, 0.0, 0.0, 0.0, 0.0, 0.0], 0.1).unwrap();
        assert!(ball.add([0.0, 0.07, 0.0, 0.0, 0.0, 0.0]));
        assert!(ball.admitted()); // sqrt(0.06^2+0.07^2), not 0.06+0.07.
        let mut parallel = FrobeniusStrainBall::new([0.06, 0.0, 0.0, 0.0, 0.0, 0.0], 0.1).unwrap();
        assert!(parallel.add([0.07, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert!(!parallel.admitted());
    }

    #[test]
    fn orthogonal_loads_and_correlated_loads_have_distinct_extrema() {
        let mut orthogonal = FrobeniusStrainBall::new([0.0; 6], 0.1).unwrap();
        assert!(orthogonal.add([0.09, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert!(orthogonal.add([0.0, 0.09, 0.0, 0.0, 0.0, 0.0]));
        assert!(orthogonal.admitted());
        let mut correlated = FrobeniusStrainBall::new([0.0; 6], 0.1).unwrap();
        for _ in 0..2 {
            assert!(correlated.add([0.09, 0.0, 0.0, 0.0, 0.0, 0.0]));
        }
        assert!(!correlated.admitted());
    }

    #[test]
    fn full_orthogonal_load_space_and_rotated_offset_stay_admissible() {
        let mut orthogonal = FrobeniusStrainBall::new([0.0; 6], 0.1).unwrap();
        for axis in 0..6 {
            let mut direction = [0.0; 6];
            direction[axis] = 0.09;
            assert!(orthogonal.add(direction));
        }
        assert!(orthogonal.admitted());
        let mut rotated =
            FrobeniusStrainBall::new([0.036, 0.048, 0.0, 0.0, 0.0, 0.0], 0.1).unwrap();
        assert!(rotated.add([-0.056, 0.042, 0.0, 0.0, 0.0, 0.0]));
        assert!(rotated.admitted());
    }
}
