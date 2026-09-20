//! Local numerical resolution of parameter combinations from calibration rows.
//! Assemble the domain's scaled J^T J and reuse fs-la's Jacobi eigensolver.
//! This is not structural/global identifiability, a posterior, or a covariance.
use super::*;

/// Numerical information of the weighted residual Jacobian in declared x units.
#[derive(Clone, Debug, PartialEq)]
pub struct ObservationInformation {
    /// Largest absolute Jacobian entry, factored out before products.
    pub jacobian_scale: f64,
    /// Numerical rank at the requested threshold, None when a mode lies within
    /// the numerical guard of that threshold. Not an interval-certified rank.
    pub numerical_rank: Option<usize>,
    /// Relative singular-value threshold used to classify parameter directions.
    pub relative_threshold: f64,
    /// Singular-value ratios sigma_i/sigma_max, ascending. All zero for J=0.
    /// Computed from a Gram matrix: accuracy is limited by squared conditioning.
    pub relative_singular_values: Vec<f64>,
    /// Unit parameter directions, one vector per ratio, in x coordinates.
    /// Largest-magnitude entry is nonnegative; degenerate subspaces have no
    /// physically distinguished basis or promised cross-ISA vector identity.
    pub directions: Vec<Vec<f64>>,
    /// Condition number only when the threshold resolves every column.
    pub condition_number: Option<f64>,
    /// Recomputed maximum eigenpair infinity residual / Gram infinity norm.
    pub eigen_residual_relative: f64,
    /// Numerical rounding/residual guard divided by the largest eigenvalue.
    /// This is a numerical admission screen, NOT an outward-rounded enclosure.
    pub relative_gram_guard: f64,
}

/// Validate before an expensive producer when this analysis will be requested.
/// The 32-column cap bounds the nonpreemptible existing Jacobi sweep. Gram
/// formation squares conditioning, so sub-1e-5 relative rank requests refuse.
pub fn admit_information(columns: usize, relative_threshold: f64) -> Result<(), DesignError> {
    if !(1..=32).contains(&columns) || !relative_threshold.is_finite()
        || !(1e-5..0.5).contains(&relative_threshold) {
        return Err(bad("local information requires 1..=32 columns and relative threshold in [1e-5,0.5)"));
    }
    Ok(())
}

impl ObservationEvaluation {
    /// Resolve the weighted observation Jacobian's locally observable parameter
    /// combinations. No additional physics or adjoints. Scales and weights are
    /// exactly those of the declared design; changing them changes the question.
    /// No normal-equation inversion or covariance is fabricated for weak modes.
    ///
    /// Cancellation is checked during Gram construction and verification and
    /// around, not inside, the bounded 32-column fs-la eigensolver. Recomputed
    /// eigen residuals and orthogonality gate publication. A mode too close to
    /// the threshold leaves rank unresolved rather than claiming full rank.
    /// Full numerical rank does not establish global uniqueness, uncertainty,
    /// physical validity, or identifiability under active design constraints.
    pub fn information(&self, relative_threshold: f64, gate: &CancelGate)
        -> Result<ObservationInformation, DesignError>
    {
        checkpoint(gate)?;
        let n = self.point.len();
        admit_information(n, relative_threshold)?;
        let scale = self.rows.iter().flat_map(|r| &r.residual_gradient)
            .fold(0.0_f64, |s, x| s.max(x.abs()));
        let divisor = if scale == 0.0 { 1.0 } else { scale };
        let mut gram = vec![0.0; n*n];
        let mut correction = vec![0.0; n*n];
        for row in &self.rows {
            checkpoint(gate)?;
            for j in 0..n {
                for k in 0..=j {
                    let i = j*n+k;
                    let term = (row.residual_gradient[j]/divisor) * (row.residual_gradient[k]/divisor);
                    let adjusted = term - correction[i];
                    let next = gram[i] + adjusted;
                    correction[i] = (next - gram[i]) - adjusted;
                    gram[i] = next;
                }
            }
        }
        for j in 0..n { for k in 0..j { gram[k*n+j] = gram[j*n+k]; } }
        let norm = gram.chunks_exact(n).map(|r| r.iter().map(|v| v.abs()).sum::<f64>()).fold(0.0_f64, f64::max);
        checkpoint(gate)?;
        let (eigenvalues, vectors) = fs_la::eigen::jacobi_eigh(&gram, n);
        checkpoint(gate)?;
        if eigenvalues.iter().chain(&vectors).any(|v| !v.is_finite()) {
            return Err(bad("local information eigensolve returned nonfinite values"));
        }
        let mut residual = 0.0_f64;
        let mut orthogonality = 0.0_f64;
        for k in 0..n {
            checkpoint(gate)?;
            for i in 0..n {
                let applied = (0..n).map(|j| gram[i*n+j]*vectors[j*n+k]).sum::<f64>();
                residual = residual.max((applied - eigenvalues[k]*vectors[i*n+k]).abs());
                let inner = (0..n).map(|j| vectors[j*n+i]*vectors[j*n+k]).sum::<f64>();
                orthogonality = orthogonality.max((inner - if i == k { 1.0 } else { 0.0 }).abs());
            }
        }
        if residual > 1e-10 * norm || orthogonality > 1e-10 {
            return Err(bad("local information eigensolve failed residual/orthogonality admission"));
        }
        let guard = (64.0*n as f64*f64::EPSILON*norm).max(n as f64*residual);
        if eigenvalues[0] < -guard { return Err(bad("local information Gram spectrum is numerically indefinite")); }
        let largest = eigenvalues[n-1].max(0.0);
        let cutoff = relative_threshold*relative_threshold*largest;
        let ambiguous = largest > 0.0 && eigenvalues.iter().any(|v| (*v-cutoff).abs() <= guard);
        let rank = eigenvalues.iter().filter(|v| **v > cutoff).count();
        let ratios: Vec<f64> = eigenvalues.iter().map(|v|
            if largest == 0.0 { 0.0 } else { (v.max(0.0)/largest).sqrt() }).collect();
        let mut directions = Vec::with_capacity(n);
        for k in 0..n {
            let mut vector: Vec<f64> = (0..n).map(|i| vectors[i*n+k]).collect();
            let mut pivot = 0;
            for i in 1..n { if vector[i].abs() > vector[pivot].abs() { pivot = i; } }
            if vector[pivot] < 0.0 { for x in &mut vector { *x = -*x; } }
            directions.push(vector);
        }
        checkpoint(gate)?;
        Ok(ObservationInformation { jacobian_scale: scale,
            numerical_rank: if ambiguous { None } else { Some(rank) }, relative_threshold,
            condition_number: if !ambiguous && rank == n { Some(1.0/ratios[0]) } else { None },
            relative_singular_values: ratios, directions,
            eigen_residual_relative: if norm == 0.0 { 0.0 } else { residual/norm },
            relative_gram_guard: if largest == 0.0 { 0.0 } else { guard/largest } })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn observations(matrix: &[Vec<f64>]) -> ObservationEvaluation {
        let n = matrix[0].len();
        ObservationEvaluation { point: vec![0.0; n], physical_parameters: vec![0.0; n],
            rows: matrix.iter().enumerate().map(|(i, row)| ObservationRow {
                case: 0, target: i, value_m: 0.0, gradient_m: row.clone(), residual: 0.0,
                residual_gradient: row.clone(), adjoint_relative_residual: 0.0,
            }).collect(), equilibria: vec![], unassessed_response_constraints: 0,
            objective: 0.0, gradient: vec![0.0; n] }
    }
    #[test]
    fn duplicate_columns_expose_an_unobservable_combination_despite_zero_loss() {
        let data = observations(&[vec![1.0, 1.0], vec![2.0, 2.0], vec![-1.0, -1.0]]);
        let report = data.information(1e-5, &CancelGate::new()).unwrap();
        assert_eq!(report.numerical_rank, Some(1)); assert_eq!(report.condition_number, None);
        let weak = &report.directions[0]; assert!((weak[0]+weak[1]).abs() < 1e-12);
        assert!(data.gauss_newton_product(weak, &CancelGate::new()).unwrap().iter().all(|x| x.abs() < 1e-10));
        assert_eq!(report, data.information(1e-5, &CancelGate::new()).unwrap());
    }
    #[test]
    fn information_is_scale_stable_and_distinguishes_two_independent_parameters() {
        for scale in [1.0, 1e200, 1e-200] {
            let data = observations(&[vec![scale, 0.0], vec![0.0, 2.0*scale]]);
            let report = data.information(1e-5, &CancelGate::new()).unwrap();
            assert_eq!(report.numerical_rank, Some(2));
            assert!((report.condition_number.unwrap()-2.0).abs() < 1e-12);
            assert_eq!(report.relative_singular_values, vec![0.5, 1.0]);
        }
    }
    #[test]
    fn zero_information_and_threshold_ambiguity_do_not_become_full_rank() {
        let zero = observations(&[vec![0.0, 0.0]]).information(1e-5, &CancelGate::new()).unwrap();
        assert_eq!(zero.numerical_rank, Some(0)); assert_eq!(zero.condition_number, None);
        let borderline = observations(&[vec![1.0, 0.0], vec![0.0, 1e-4]])
            .information(1e-4, &CancelGate::new()).unwrap();
        assert_eq!(borderline.numerical_rank, None); assert_eq!(borderline.condition_number, None);
    }
    #[test]
    fn numerical_envelope_and_cancellation_are_explicit() {
        for tolerance in [0.0, 1e-6, 0.5, f64::NAN] { assert!(admit_information(2, tolerance).is_err()); }
        assert!(admit_information(33, 1e-5).is_err()); assert!(admit_information(0, 1e-5).is_err());
        let gate = CancelGate::new(); gate.request();
        assert!(matches!(observations(&[vec![1.0]]).information(1e-5, &gate), Err(DesignError::Cancelled)));
    }
}
