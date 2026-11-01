//! Shared fixed balancing and original-coordinate audits for loaded boundaries.
//! Geometry-specific callers still assemble their own physical equations.
use super::*;

impl PhysicalPencil {
    /// Fixed two-sided scales from one declared reference operator. Frequency
    /// callbacks only APPLY these scales; they never rebalance a contour node.
    pub(super) fn fixed_scales(
        &self,
        n: usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), BemModesError> {
        let mut row_scales = vec![1.0; n];
        let mut column_scales = vec![1.0; n];
        for i in 0..n {
            let row: f64 = self.envelope[i * n..(i + 1) * n].iter().sum();
            if !row.is_finite() || row <= 0.0 {
                return Err(bad("zero or unrepresentable reference row"));
            }
            row_scales[i] = row;
        }
        for j in 0..n {
            let col = (0..n)
                .map(|i| self.envelope[i * n + j] / row_scales[i])
                .fold(0.0, f64::max);
            let scale = 1.0 / col;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(bad("zero or unrepresentable reference coordinate"));
            }
            column_scales[j] = scale;
        }
        for i in 0..n {
            let row: f64 = (0..n)
                .map(|j| self.envelope[i * n + j] * column_scales[j])
                .sum();
            if !row.is_finite() || row <= 0.0 {
                return Err(bad("unrepresentable balanced reference row"));
            }
            row_scales[i] = row;
        }
        poll(cancelled)?;
        Ok((row_scales, column_scales))
    }

    pub(super) fn into_balanced(
        mut self,
        row_scales: &[f64],
        column_scales: &[f64],
    ) -> Result<PoleEvaluation, BemModesError> {
        let n = row_scales.len();
        let mut scale = 0.0_f64;
        for i in 0..n {
            let mut row = 0.0;
            for j in 0..n {
                let k = i * n + j;
                let factor = column_scales[j] / row_scales[i];
                self.value.matrix[k] = self.value.matrix[k].scale(factor);
                self.value.derivative[k] = self.value.derivative[k].scale(factor);
                row += self.envelope[k] * factor;
            }
            scale = scale.max(row);
        }
        if !scale.is_finite()
            || scale <= 0.0
            || self
                .value
                .matrix
                .iter()
                .chain(&self.value.derivative)
                .any(|&z| !finite(z))
        {
            return Err(bad("balanced coupled operator is not representable"));
        }
        self.value.residual_scale = scale;
        Ok(self.value)
    }

    /// Shared compensated audit of original mechanical, boundary and memory
    /// equations. A source RHS acts on the first mechanical block only.
    pub(super) fn physical_residuals(
        &self,
        modes: usize,
        boundary: usize,
        vector: &[C64],
        force: Option<&[C64]>,
        tolerance: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<[f64; 3], BemModesError> {
        let n = vector.len();
        if n.checked_mul(n) != Some(self.value.matrix.len())
            || self.envelope.len() != self.value.matrix.len()
            || modes.checked_add(boundary).is_none_or(|v| v > n)
            || force.is_some_and(|f| f.len() != modes)
        {
            return Err(bad(
                "physical residual needs complete state and force coordinates",
            ));
        }
        let mut residuals = [0.0_f64; 3];
        for i in 0..n {
            poll(cancelled)?;
            let rhs = force.and_then(|f| f.get(i)).copied().unwrap_or(C64::ZERO);
            let mut sum = -rhs;
            let mut correction = C64::ZERO;
            let mut denominator = rhs.abs();
            for j in 0..n {
                let term = self.value.matrix[i * n + j] * vector[j] - correction;
                let next = sum + term;
                correction = (next - sum) - term;
                sum = next;
                denominator += self.envelope[i * n + j] * vector[j].abs();
            }
            if !finite(sum) || !denominator.is_finite() {
                return Err(bad("physical residual contraction overflow"));
            }
            let relative = if denominator == 0.0 {
                if sum.abs() == 0.0 { 0.0 } else { f64::INFINITY }
            } else {
                sum.abs() / denominator
            };
            if !relative.is_finite() {
                return Err(bad("physical residual is not representable"));
            }
            let block = if i < modes {
                0
            } else if i < modes + boundary {
                1
            } else {
                2
            };
            residuals[block] = residuals[block].max(relative);
        }
        if residuals.iter().any(|&r| r > tolerance) {
            return Err(BemModesError::PhysicalResidual {
                mechanical: residuals[0],
                boundary: residuals[1],
                memory: residuals[2],
            });
        }
        poll(cancelled)?;
        Ok(residuals)
    }
}

/// Original physical amplitudes, normalized by maximum structural displacement
/// or (for a zero-motion mode) maximum remaining original coordinate.
pub(super) fn physical_eigenstate(
    pole: &NonlinearPole,
    column_scales: &[f64],
    modes: usize,
) -> Result<Vec<C64>, BemModesError> {
    let mut vector: Vec<_> = pole
        .vector
        .iter()
        .zip(column_scales)
        .map(|(&v, &s)| v.scale(s))
        .collect();
    let mut amplitude = vector[..modes].iter().map(|v| v.abs()).fold(0.0, f64::max);
    if amplitude == 0.0 {
        amplitude = vector.iter().map(|v| v.abs()).fold(0.0, f64::max);
    }
    if !amplitude.is_finite() || amplitude <= 0.0 {
        return Err(bad("original eigenstate scale is not representable"));
    }
    for v in &mut vector {
        *v = divide(*v, amplitude);
    }
    if vector.iter().any(|&v| !finite(v)) {
        return Err(bad("original eigenstate overflow"));
    }
    Ok(vector)
}
