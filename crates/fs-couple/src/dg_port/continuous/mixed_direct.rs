//! Original mixed-rank reciprocal terminal constraints, with every state kept.
//!
//! x'=A*x+F*q+b*u and H*x=D*q. A fixed complete orthogonal chart separates
//! the positive direct range from its null space. Range flow is algebraic;
//! the original null pressure constraints are differentiated exactly once by
//! the existing zero-direct owner. Only zero-storage flow multipliers disappear.
//! The original D is retained for independent closure and work observations.

use super::super::terminal::TerminalSolve;
use super::zero_direct::{ConstrainedTerminal, eliminate as eliminate_zero};
use super::{AcousticRealizeError, invalid};
use fs_exec::CancelGate;
use fs_la::eigen::{MAX_EIGEN_UNPOLLED_SCALAR_WORK, admit_jacobi_eigh, jacobi_eigh};
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

fn check(
    residual: f64,
    scale: f64,
    tolerance: f64,
    what: &'static str,
) -> Result<(), AcousticRealizeError> {
    if !residual.is_finite() || !scale.is_finite() || residual.abs() > tolerance * scale {
        Err(invalid(what))
    } else {
        Ok(())
    }
}

// A spectral cutoff alone cannot decide physical rank: a small but nonzero
// resistance still carries a fast original state. Inspect original-coordinate
// principal minors through the shared LU owner, without a pivot tolerance.
// One added coordinate detects extra PSD rank; two also detect a hidden tiny
// +/- pair whose individual Schur diagonals vanish. Rounded rank ambiguity is
// a refusal, never permission to replace the supplied D by a lower-rank matrix.
fn admit_original_rank(
    direct: &[f64],
    p: usize,
    rank: usize,
    gate: &CancelGate,
) -> Result<(), AcousticRealizeError> {
    let nullity = p - rank;
    let count = (nullity as u128) * (nullity as u128 + 1) / 2;
    let work = 720u128 * (p as u128).pow(3)
        + 256u128 * (p as u128).pow(2)
        + 16u128 * p as u128 * (rank as u128 + 1).pow(3)
        + 8u128 * count * (rank as u128 + 2).pow(3);
    if work > MAX_EIGEN_UNPOLLED_SCALAR_WORK as u128 {
        return Err(invalid(
            "mixed-direct original-rank checks exceed shared work cap",
        ));
    }
    // Prefer the strongest original diagonal. A dependent candidate is
    // skipped; every selected principal block must pass the existing complete
    // solve's condition gate. These indices are only an admission witness,
    // not a new physical basis or a replacement factorization kernel.
    let mut candidates: Vec<_> = (0..p).collect();
    candidates.sort_by(|&a, &b| {
        direct[b * p + b]
            .total_cmp(&direct[a * p + a])
            .then(a.cmp(&b))
    });
    let mut selected = Vec::with_capacity(rank + 2);
    for candidate in candidates {
        poll(gate)?;
        if selected.len() == rank {
            break;
        }
        selected.push(candidate);
        let block: Vec<_> = selected
            .iter()
            .flat_map(|&i| selected.iter().map(move |&j| direct[i * p + j]))
            .collect();
        if TerminalSolve::new(&block, selected.len()).is_err() {
            selected.pop();
        }
    }
    if selected.len() != rank {
        return Err(invalid(
            "mixed-direct original positive principal block is unresolved",
        ));
    }
    let remaining: Vec<_> = (0..p).filter(|i| !selected.contains(i)).collect();
    let mut single = LuWorkspace::new(rank + 1)
        .map_err(|_| invalid("mixed-direct original-rank workspace refused"))?;
    let mut pair = if nullity > 1 {
        Some(
            LuWorkspace::new(rank + 2)
                .map_err(|_| invalid("mixed-direct original-rank workspace refused"))?,
        )
    } else {
        None
    };
    let mut block = Vec::with_capacity((rank + 2) * (rank + 2));
    for a in 0..nullity {
        for b in a..nullity {
            poll(gate)?;
            selected.truncate(rank);
            selected.push(remaining[a]);
            if b != a {
                selected.push(remaining[b]);
            }
            block.clear();
            for &i in &selected {
                for &j in &selected {
                    block.push(direct[i * p + j]);
                }
            }
            let factor = if b == a {
                &mut single
            } else {
                pair.as_mut()
                    .ok_or_else(|| invalid("mixed-direct pair workspace is missing"))?
            };
            match factor.factor(&block) {
                Err(LuWorkspaceError::Singular { index }) if index == rank => {}
                _ => {
                    return Err(invalid(
                        "mixed-direct original nonzero pivots disagree with spectral null rank",
                    ));
                }
            }
        }
    }
    poll(gate)
}

/// Close an original reciprocal PSD mixed-rank direct term. Callers retain
/// their existing full-rank and exactly-zero paths. Numerical null directions
/// must agree with nonzero pivots of original-coordinate principal minors and
/// independently satisfy original D rows and componentwise range reconstruction; a resolvable small positive/negative direction refuses.
/// There is no diagonal shift, eigenvalue clipping, source rewrite or state
/// projection. Near-rank ambiguity that fails these original checks refuses.
///
/// Shapes/conventions match zero_direct: inputs are p original flow columns,
/// pressure_rows are p rows H, and inlet is one external forcing column. Held
/// input coordinates may already be included in A. Outputs return ORIGINAL
/// port flow rows and inlet factors, not range/null coordinates. Remaining
/// pressure constraints are normalized exactly as the existing zero owner.
///
/// Complete n<=256, p<=n; callers reserve at least256*n² numeric bytes for
/// this cold path, separately from retained/exact-propagation storage. Shared
/// Jacobi workspace/unpolled-work caps apply before allocation. Higher-index,
/// unresolved, nonreciprocal and nonfinite sources refuse without resetting
/// history. The initial state must satisfy the returned original manifold.
pub(crate) fn eliminate(
    source: &[f64],
    inputs: &[Vec<f64>],
    pressure_rows: &[Vec<f64>],
    direct: &[f64],
    inlet: &[f64],
    gate: &CancelGate,
) -> Result<ConstrainedTerminal, AcousticRealizeError> {
    poll(gate)?;
    let n = inlet.len();
    let p = pressure_rows.len();
    if n == 0
        || n > 256
        || p < 2
        || p > n
        || source.len() != n * n
        || direct.len() != p * p
        || inputs.len() != p
        || inputs.iter().chain(pressure_rows).any(|row| row.len() != n)
        || source
            .iter()
            .chain(inlet)
            .chain(direct)
            .chain(inputs.iter().flatten())
            .chain(pressure_rows.iter().flatten())
            .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "mixed-direct terminal needs complete bounded finite original matrices",
        ));
    }
    let work = 720u128 * (p as u128).pow(3) + 256u128 * (p as u128).pow(2);
    if work > MAX_EIGEN_UNPOLLED_SCALAR_WORK as u128 {
        return Err(invalid(
            "mixed-direct spectral chart exceeds shared unpolled work cap",
        ));
    }
    admit_jacobi_eigh(p).map_err(|_| invalid("mixed-direct spectral chart workspace refused"))?;
    let norm = direct.iter().map(|v| v.abs()).fold(0.0, f64::max);
    if norm == 0.0 {
        return Err(invalid(
            "exactly zero direct terms belong to the existing zero-direct owner",
        ));
    }
    for i in 0..p {
        for j in 0..i {
            if direct[i * p + j] != direct[j * p + i] {
                return Err(invalid(
                    "mixed-direct terminal requires an original reciprocal direct matrix",
                ));
            }
        }
    }
    let normalized: Vec<_> = direct.iter().map(|v| v / norm).collect();
    poll(gate)?;
    let (values, vectors) = jacobi_eigh(&normalized, p);
    poll(gate)?;
    let tolerance = 512.0 * f64::EPSILON * p as f64;
    let null_threshold = 32.0 * f64::EPSILON * p as f64;
    let largest = values.iter().copied().fold(0.0, f64::max);
    if !largest.is_finite()
        || largest <= 0.0
        || values.iter().chain(&vectors).any(|v| !v.is_finite())
        || values.iter().any(|v| *v < -null_threshold * largest)
    {
        return Err(invalid(
            "mixed-direct source has no resolved positive semidefinite chart",
        ));
    }
    let null: Vec<_> = (0..p)
        .filter(|&j| values[j].abs() <= null_threshold * largest)
        .collect();
    let range: Vec<_> = (0..p).filter(|j| !null.contains(j)).collect();
    if null.is_empty() || range.is_empty() {
        return Err(invalid(
            "mixed-direct fallback requires a resolved nonempty range and null space",
        ));
    }
    admit_original_rank(direct, p, range.len(), gate)?;
    // Original D*v=0 is tested componentwise. In particular diag(1,tiny)
    // cannot be silently treated as rank one: its isolated tiny row refuses.
    for &column in &null {
        for i in 0..p {
            poll(gate)?;
            let mut residual = 0.0;
            let mut scale = 0.0;
            for j in 0..p {
                let term = direct[i * p + j] * vectors[j * p + column];
                residual += term;
                scale += term.abs();
            }
            check(
                residual,
                scale,
                tolerance,
                "direct null direction is unresolved in the original matrix",
            )?;
        }
    }
    for i in 0..p {
        poll(gate)?;
        for j in 0..p {
            let mut rebuilt = 0.0;
            let mut scale = normalized[i * p + j].abs();
            for &k in &range {
                let term = vectors[i * p + k] * values[k] * vectors[j * p + k];
                rebuilt += term;
                scale += term.abs();
            }
            check(
                rebuilt - normalized[i * p + j],
                scale,
                tolerance,
                "positive direct range fails original componentwise reconstruction",
            )?;
            let orthogonal = (0..p)
                .map(|k| vectors[k * p + i] * vectors[k * p + j])
                .sum::<f64>();
            check(
                orthogonal - if i == j { 1.0 } else { 0.0 },
                1.0,
                tolerance,
                "direct range/null chart is not orthogonal",
            )?;
        }
    }
    let r = range.len();
    let k = null.len();
    // Solve the ACTUAL original range block, rather than replacing it by the
    // computed eigenvalue diagonal. The complete original D stays authoritative.
    let mut range_d = vec![0.0; r * r];
    let mut range_h = vec![vec![0.0; n]; r];
    for (a, &ia) in range.iter().enumerate() {
        poll(gate)?;
        for (b, &ib) in range.iter().enumerate() {
            range_d[a * r + b] = (0..p)
                .map(|i| {
                    vectors[i * p + ia]
                        * (0..p)
                            .map(|j| direct[i * p + j] * vectors[j * p + ib])
                            .sum::<f64>()
                })
                .sum();
        }
        for j in 0..n {
            range_h[a][j] = (0..p)
                .map(|i| vectors[i * p + ia] * pressure_rows[i][j])
                .sum();
        }
    }
    let range_solve = TerminalSolve::new(&range_d, r)
        .map_err(|_| invalid("mixed-direct positive range is unresolved or ill-conditioned"))?;
    let mut range_flow = vec![vec![0.0; n]; p];
    for j in 0..n {
        poll(gate)?;
        let rhs: Vec<_> = range_h.iter().map(|row| row[j]).collect();
        let q = range_solve.apply(&rhs)?;
        for i in 0..p {
            range_flow[i][j] = range
                .iter()
                .enumerate()
                .map(|(a, &ia)| vectors[i * p + ia] * q[a])
                .sum();
        }
    }
    let mut ranged = source.to_vec();
    for i in 0..n {
        poll(gate)?;
        for j in 0..n {
            for a in 0..p {
                ranged[i * n + j] += inputs[a][i] * range_flow[a][j];
            }
        }
    }
    let mut null_h = vec![vec![0.0; n]; k];
    let mut null_f = vec![vec![0.0; n]; k];
    for (a, &ia) in null.iter().enumerate() {
        poll(gate)?;
        for j in 0..n {
            null_h[a][j] = (0..p)
                .map(|i| vectors[i * p + ia] * pressure_rows[i][j])
                .sum();
            null_f[a][j] = (0..p).map(|i| vectors[i * p + ia] * inputs[i][j]).sum();
        }
    }
    let constrained = eliminate_zero(&ranged, &null_f, &null_h, inlet, gate)?;
    let mut flow_columns = range_flow;
    let mut inlet_flow_factors = vec![0.0; p];
    for i in 0..p {
        poll(gate)?;
        for (a, &ia) in null.iter().enumerate() {
            let v = vectors[i * p + ia];
            inlet_flow_factors[i] += v * constrained.inlet_flow_factors[a];
            for j in 0..n {
                flow_columns[i][j] += v * constrained.flow_columns[a][j];
            }
        }
    }
    // H-D*Q must equal precisely the retained null pressure manifold, and
    // D*q_inlet must vanish. These use original port rows after reconstruction.
    for i in 0..p {
        poll(gate)?;
        for j in 0..=n {
            let mut residual = if j < n { pressure_rows[i][j] } else { 0.0 };
            let mut scale = residual.abs();
            for a in 0..p {
                let term = direct[i * p + a]
                    * if j < n {
                        flow_columns[a][j]
                    } else {
                        inlet_flow_factors[a]
                    };
                residual -= term;
                scale += term.abs();
            }
            if j < n {
                for (a, &ia) in null.iter().enumerate() {
                    let term = vectors[i * p + ia] * null_h[a][j];
                    residual -= term;
                    scale += term.abs();
                }
            }
            check(
                residual,
                scale,
                tolerance * 4.0,
                "mixed-direct original pressure closure is unresolved",
            )?;
        }
    }
    // The returned generator/input must use the same original flow maps.
    // This checks the reconstruction sign independently of null elimination.
    for i in 0..n {
        poll(gate)?;
        for j in 0..=n {
            let original = if j < n { source[i * n + j] } else { inlet[i] };
            let actual = if j < n {
                constrained.generator[i * n + j]
            } else {
                constrained.inlet_column[i]
            };
            let mut residual = actual - original;
            let mut scale = actual.abs() + original.abs();
            for a in 0..p {
                let term = inputs[a][i]
                    * if j < n {
                        flow_columns[a][j]
                    } else {
                        inlet_flow_factors[a]
                    };
                residual -= term;
                scale += term.abs();
            }
            check(
                residual,
                scale,
                tolerance * 4.0,
                "mixed-direct original flow reconstruction is unresolved",
            )?;
        }
    }
    if flow_columns
        .iter()
        .flatten()
        .chain(&inlet_flow_factors)
        .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "mixed-direct original flow map is not representable",
        ));
    }
    poll(gate)?;
    Ok(ConstrainedTerminal {
        generator: constrained.generator,
        inlet_column: constrained.inlet_column,
        flow_columns,
        inlet_flow_factors,
        constraints: constrained.constraints,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() < 2e-12 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    }
    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b).map(|(a, b)| a * b).sum()
    }
    #[test]
    fn rank_one_resistance_has_independent_differential_flow_and_correct_common_decay() {
        // Two unit acoustic compliances share P1=P2=R*(U1+U2). The imposed
        // first-node flow gives U1-U2=input; x'=input/2-x/(2R).
        let gate = CancelGate::new_clock_free();
        let r = 2.0;
        let a = eliminate(
            &[0.0; 4],
            &[vec![-1.0, 0.0], vec![0.0, -1.0]],
            &[vec![1.0, 0.0], vec![0.0, 1.0]],
            &[r, r, r, r],
            &[1.0, 0.0],
            &gate,
        )
        .unwrap();
        let x = [3.0, 3.0];
        let input = 4.0;
        let q: Vec<_> = a
            .flow_columns
            .iter()
            .zip(&a.inlet_flow_factors)
            .map(|(row, b)| dot(row, &x) + b * input)
            .collect();
        close(q[0], 3.0 / (2.0 * r) + input / 2.0);
        close(q[1], 3.0 / (2.0 * r) - input / 2.0);
        for i in 0..2 {
            close(
                dot(&a.generator[2 * i..2 * i + 2], &x) + a.inlet_column[i] * input,
                input / 2.0 - 3.0 / (2.0 * r),
            );
        }
        assert_eq!(a.constraints.len(), 1);
        close(dot(&a.constraints[0], &x), 0.0);
        // Multiple held incident coordinates use the SAME elimination once.
        let mut source = vec![0.0; 16];
        source[2] = 1.0;
        source[7] = 1.0;
        let b = eliminate(
            &source,
            &[vec![-1.0, 0.0, 0.0, 0.0], vec![0.0, -1.0, 0.0, 0.0]],
            &[vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 1.0, 0.0, 0.0]],
            &[r, r, r, r],
            &[0.0; 4],
            &gate,
        )
        .unwrap();
        let state = [3.0, 3.0, 4.0, -2.0];
        close(dot(&b.flow_columns[0], &state), 3.0 / (2.0 * r) + 3.0);
        close(dot(&b.flow_columns[1], &state), 3.0 / (2.0 * r) - 3.0);
        for i in 0..2 {
            close(
                dot(&b.generator[4 * i..4 * i + 4], &state),
                1.0 - 3.0 / (2.0 * r),
            );
        }
    }
    #[test]
    fn signed_port_covariance_preserves_original_state_dynamics() {
        let gate = CancelGate::new_clock_free();
        let source = [-0.2, 0.0, 0.0, -0.4];
        let base = eliminate(
            &source,
            &[vec![-1.0, 0.0], vec![0.0, -1.0]],
            &[vec![1.0, 0.0], vec![0.0, 1.0]],
            &[2.0; 4],
            &[1.0, -0.3],
            &gate,
        )
        .unwrap();
        // Signed orthogonal reordering S=[[0,1],[-1,0]], q_old=S*q_new.
        let changed = eliminate(
            &source,
            &[vec![0.0, 1.0], vec![-1.0, 0.0]],
            &[vec![0.0, -1.0], vec![1.0, 0.0]],
            &[2.0, -2.0, -2.0, 2.0],
            &[1.0, -0.3],
            &gate,
        )
        .unwrap();
        for (a, b) in base.generator.iter().zip(&changed.generator) {
            close(*a, *b);
        }
        for (a, b) in base.inlet_column.iter().zip(&changed.inlet_column) {
            close(*a, *b);
        }
        for j in 0..2 {
            close(base.flow_columns[0][j], changed.flow_columns[1][j]);
            close(base.flow_columns[1][j], -changed.flow_columns[0][j]);
        }
        close(base.inlet_flow_factors[0], changed.inlet_flow_factors[1]);
        close(base.inlet_flow_factors[1], -changed.inlet_flow_factors[0]);
    }
    #[test]
    fn unresolved_rank_higher_index_and_nonreciprocal_sources_refuse() {
        let gate = CancelGate::new_clock_free();
        let inputs = [vec![-1.0, 0.0], vec![0.0, -1.0]];
        let h = [vec![1.0, 0.0], vec![0.0, 1.0]];
        for d in [
            [1.0, 0.0, 0.0, 1e-18],
            [1.0, 0.0, 0.0, 1e-320],
            [1.0, 0.0, 0.0, -1e-18],
            [1.0, 0.2, 0.0, 0.0],
            [1.0, 0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0, 1.0],
        ] {
            assert!(eliminate(&[0.0; 4], &inputs, &h, &d, &[0.0; 2], &gate).is_err());
        }
        assert!(
            eliminate(
                &[0.0; 4],
                &[vec![0.0; 2], vec![0.0; 2]],
                &h,
                &[1.0; 4],
                &[0.0; 2],
                &gate
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            eliminate(&[0.0; 4], &inputs, &h, &[1.0; 4], &[0.0; 2], &gate),
            Err(AcousticRealizeError::Cancelled)
        ));
    }
    #[test]
    fn rotated_nonzero_directions_are_never_replaced_by_constraints() {
        let gate = CancelGate::new_clock_free();
        for p in [2, 3] {
            for delta in [1e-14, -1e-14] {
                let mut direct = vec![1.0; p * p];
                for i in 1..p {
                    for j in 1..p {
                        direct[i * p + j] += delta;
                    }
                }
                // With p=3, duplicate rows leave a genuine null in addition to
                // the small positive/negative direction. Full-LU singularity
                // alone must not admit dropping that additional physical mode.
                let h: Vec<Vec<f64>> = (0..p)
                    .map(|i| (0..p).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
                    .collect();
                let f: Vec<Vec<f64>> = h
                    .iter()
                    .map(|row| row.iter().map(|v| -v).collect())
                    .collect();
                let mut inlet = vec![0.0; p];
                inlet[0] = 1.0;
                assert!(eliminate(&vec![0.0; p * p], &f, &h, &direct, &inlet, &gate).is_err());
            }
        }
        let p = 4;
        let h: Vec<Vec<f64>> = (0..p)
            .map(|i| (0..p).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let f: Vec<Vec<f64>> = h
            .iter()
            .map(|row| row.iter().map(|v| -v).collect())
            .collect();
        let mut direct = vec![1.0; p * p];
        direct[p + 2] += 1e-14;
        direct[2 * p + 1] += 1e-14;
        // Every one-coordinate augmentation of principal index0 is singular;
        // the two-coordinate minor exposes the hidden indefinite +/- pair.
        assert!(
            eliminate(
                &vec![0.0; p * p],
                &f,
                &h,
                &direct,
                &[1.0, 0.0, 0.0, 0.0],
                &gate
            )
            .is_err()
        );
    }

    #[test]
    fn independent_positive_channels_and_multiple_exact_nulls_survive_rank_admission() {
        let gate = CancelGate::new_clock_free();
        for p in [3, 4] {
            let h: Vec<Vec<f64>> = (0..p)
                .map(|i| (0..p).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
                .collect();
            let f: Vec<Vec<f64>> = h
                .iter()
                .map(|row| row.iter().map(|v| -v).collect())
                .collect();
            let mut direct = vec![0.0; p * p];
            direct[0] = 1.0;
            direct[1] = 1.0;
            direct[p] = 1.0;
            direct[p + 1] = 1.0;
            direct[2 * p + 2] = 0.5;
            let result =
                eliminate(&vec![0.0; p * p], &f, &h, &direct, &vec![0.0; p], &gate).unwrap();
            assert_eq!(result.constraints.len(), p - 2);
            let mut state = vec![0.0; p];
            state[0] = 3.0;
            state[1] = 3.0;
            state[2] = 2.0;
            close(dot(&result.flow_columns[0], &state), 1.5);
            close(dot(&result.flow_columns[1], &state), 1.5);
            close(dot(&result.flow_columns[2], &state), 4.0);
            for constraint in &result.constraints {
                close(dot(constraint, &state), 0.0);
            }
        }
    }
}
