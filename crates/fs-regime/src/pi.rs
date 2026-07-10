//! Buckingham-Pi machinery: from Qty-typed inputs, compute a basis of
//! the dimensionless products ALGEBRAICALLY — the integer nullspace of
//! the 5×n dimension matrix (SI exponents m, kg, s, K, A), found by
//! exact fraction-free elimination over i128. Groups are dimensionless
//! BY CONSTRUCTION and verified so.

use crate::RegimeError;
use fs_qty::QtyAny;

/// One named, dimensioned problem input.
#[derive(Debug, Clone, PartialEq)]
pub struct Input {
    /// Variable name (for reporting).
    pub name: String,
    /// The dimensioned value.
    pub qty: QtyAny,
}

/// One dimensionless product `Π inputs[i]^exponents[i]`.
#[derive(Debug, Clone, PartialEq)]
pub struct PiGroup {
    /// Integer exponent per input (reduced: gcd 1, first nonzero > 0).
    pub exponents: Vec<i64>,
    /// The group's numerical value (inputs are SI, so this is unit-free).
    pub value: f64,
}

/// The extracted basis.
#[derive(Debug, Clone, PartialEq)]
pub struct PiBasis {
    /// Rank of the dimension matrix.
    pub rank: usize,
    /// Basis groups (count = inputs − rank, the Buckingham count).
    pub groups: Vec<PiGroup>,
}

fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// Fraction-free row echelon over i128; returns the reduced matrix and
/// the pivot (row, col) list.
fn echelon(mut m: Vec<Vec<i128>>, n: usize) -> (Vec<Vec<i128>>, Vec<(usize, usize)>) {
    let mut pivots: Vec<(usize, usize)> = Vec::new();
    let mut row = 0usize;
    for col in 0..n {
        let Some(pivot_row) = (row..5).find(|&r| m[r][col] != 0) else {
            continue;
        };
        m.swap(row, pivot_row);
        let p = m[row][col];
        for r in 0..5 {
            if r != row && m[r][col] != 0 {
                let factor = m[r][col];
                let pivot_row_copy = m[row].clone();
                for (v, &pv) in m[r].iter_mut().zip(&pivot_row_copy) {
                    *v = *v * p - pv * factor;
                }
                let g = m[r].iter().fold(0i128, |acc, &v| gcd(acc, v));
                if g > 1 {
                    for v in &mut m[r] {
                        *v /= g;
                    }
                }
            }
        }
        pivots.push((row, col));
        row += 1;
        if row == 5 {
            break;
        }
    }
    (m, pivots)
}

/// Compute the Pi basis for a set of inputs.
///
/// # Errors
/// [`RegimeError::Degenerate`] on empty input or non-positive values
/// (groups are products of powers — signs/zeros have no regime meaning);
/// [`RegimeError::NotDimensionless`] if the construction self-check fails
/// (impossible unless the elimination is wrong — it is the G0 guard); or
/// [`RegimeError::ExponentOutOfRange`] when an exact exponent cannot be
/// evaluated by the current deterministic i32 power primitive.
pub fn pi_groups(inputs: &[Input]) -> Result<PiBasis, RegimeError> {
    let n = inputs.len();
    if n == 0 {
        return Err(RegimeError::Degenerate {
            what: "no inputs".to_string(),
        });
    }
    for input in inputs {
        if !(input.qty.value.is_finite() && input.qty.value > 0.0) {
            return Err(RegimeError::BadValue {
                what: format!(
                    "input {:?} must be finite and positive (got {})",
                    input.name, input.qty.value
                ),
            });
        }
    }
    // Dimension matrix: 5 rows (SI base dims) × n columns (inputs).
    let m: Vec<Vec<i128>> = (0..5)
        .map(|d| inputs.iter().map(|i| i128::from(i.qty.dims.0[d])).collect())
        .collect();
    let (m, pivots) = echelon(m, n);
    let rank = pivots.len();
    let pivot_cols: Vec<usize> = pivots.iter().map(|&(_, c)| c).collect();
    let free_cols: Vec<usize> = (0..n).filter(|c| !pivot_cols.contains(c)).collect();

    let mut groups = Vec::with_capacity(free_cols.len());
    for &f in &free_cols {
        // Rational solve: x[f] = 1, pivot variables from their rows.
        // x[c_pivot] = -m[row][f] / m[row][c_pivot] (echelon rows only
        // couple pivots to free columns after full elimination).
        let mut numer = vec![0i128; n];
        let mut denom = vec![1i128; n];
        numer[f] = 1;
        for &(r, c) in &pivots {
            numer[c] = -m[r][f];
            denom[c] = m[r][c];
        }
        // Clear denominators: multiply through by lcm.
        let lcm = denom
            .iter()
            .fold(1i128, |acc, &d| acc / gcd(acc, d) * d.abs().max(1));
        let mut exps: Vec<i128> = numer
            .iter()
            .zip(&denom)
            .map(|(&nu, &de)| nu * (lcm / de))
            .collect();
        let g = exps.iter().fold(0i128, |acc, &v| gcd(acc, v));
        if g > 1 {
            for e in &mut exps {
                *e /= g;
            }
        }
        if let Some(first) = exps.iter().find(|&&e| e != 0)
            && *first < 0
        {
            for e in &mut exps {
                *e = -*e;
            }
        }
        // Dimensionless-by-construction self-check (G0).
        let mut residual = [0i128; 5];
        for (input, &e) in inputs.iter().zip(&exps) {
            for (slot, &d) in residual.iter_mut().zip(&input.qty.dims.0) {
                *slot += i128::from(d) * e;
            }
        }
        if residual != [0; 5] {
            return Err(RegimeError::NotDimensionless {
                context: format!("pi group over free column {f}"),
                residual,
            });
        }
        let eval_exps: Vec<i32> = exps
            .iter()
            .enumerate()
            .map(|(index, &e)| {
                i32::try_from(e).map_err(|_| RegimeError::ExponentOutOfRange {
                    context: format!(
                        "pi group over free column {f}, input {:?}",
                        inputs[index].name
                    ),
                    exponent: e,
                })
            })
            .collect::<Result<_, _>>()?;
        let mut value = 1.0f64;
        for (input, &e) in inputs.iter().zip(&eval_exps) {
            value *= fs_math::det::powi(input.qty.value, e);
        }
        groups.push(PiGroup {
            exponents: eval_exps.iter().map(|&e| i64::from(e)).collect(),
            value,
        });
    }
    Ok(PiBasis { rank, groups })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::Dims;

    fn input(name: &str, value: f64, dims: [i8; 5]) -> Input {
        Input {
            name: name.to_string(),
            qty: QtyAny::new(value, Dims(dims)),
        }
    }

    #[test]
    fn oversized_exact_exponents_refuse_instead_of_wrapping() {
        // Valid 5x6 i8 dimension matrix with primitive null exponents
        // [246398765, 1209471509, 174635022, -831044588,
        //  -4691840893, -4000421645]. The final two exceed i32.
        let inputs = [
            input("x0", 1.0, [-67, -88, 107, -71, -83]),
            input("x1", 1.0, [-50, -104, -53, 6, 84]),
            input("x2", 1.0, [-101, -110, 77, 10, -100]),
            input("x3", 1.0, [57, -122, 126, -35, -60]),
            input("x4", 1.0, [-26, -25, 68, -57, -73]),
            input("x5", 1.0, [-5, 13, -112, 72, 114]),
        ];
        assert!(matches!(
            pi_groups(&inputs),
            Err(RegimeError::ExponentOutOfRange { exponent, .. })
                if exponent == -4_691_840_893
        ));
    }

    #[test]
    fn pipe_flow_recovers_reynolds() {
        // (ρ, V, D, μ): rank 3, one group ∝ Re = ρVD/μ.
        let inputs = [
            input("rho", 1000.0, [-3, 1, 0, 0, 0]),
            input("v", 2.0, [1, 0, -1, 0, 0]),
            input("d", 0.05, [1, 0, 0, 0, 0]),
            input("mu", 1e-3, [-1, 1, -1, 0, 0]),
        ];
        let basis = pi_groups(&inputs).expect("pi");
        assert_eq!(basis.rank, 3);
        assert_eq!(basis.groups.len(), 1);
        let g = &basis.groups[0];
        // Exponents ∝ (1, 1, 1, −1) — Re up to normalization.
        let re = 1000.0 * 2.0 * 0.05 / 1e-3;
        let val = g.value;
        assert!(
            (val - re).abs() / re < 1e-12 || (val - 1.0 / re).abs() * re < 1e-12,
            "group must be Re or 1/Re, got {val}"
        );
    }

    #[test]
    fn pendulum_recovers_period_group() {
        // (T, L, g, m): mass is dimensionally isolated ⇒ rank 3, one
        // group ∝ T²g/L (mass cannot appear).
        let inputs = [
            input("t", 2.007_1, [0, 0, 1, 0, 0]),
            input("l", 1.0, [1, 0, 0, 0, 0]),
            input("g", 9.806_65, [1, 0, -2, 0, 0]),
            input("m", 0.3, [0, 1, 0, 0, 0]),
        ];
        let basis = pi_groups(&inputs).expect("pi");
        assert_eq!(basis.rank, 3);
        assert_eq!(basis.groups.len(), 1);
        let g = &basis.groups[0];
        let mass_slot = g.exponents[3];
        assert_eq!(mass_slot, 0, "mass cannot enter a dimensionless group here");
        // T²g/L ≈ (2π)² for a unit pendulum.
        let expect = 2.007_1f64.powi(2) * 9.806_65 / 1.0;
        let val = g.value;
        assert!(
            (val - expect).abs() / expect < 1e-12 || (val - 1.0 / expect) * expect < 1e-12,
            "got {val}, want {expect} (or its inverse)"
        );
    }

    #[test]
    fn drag_problem_has_two_groups() {
        // (F, ρ, V, L, μ) → n − r = 2 (drag coefficient and Reynolds).
        let inputs = [
            input("f", 12.0, [1, 1, -2, 0, 0]),
            input("rho", 1.225, [-3, 1, 0, 0, 0]),
            input("v", 8.0, [1, 0, -1, 0, 0]),
            input("l", 0.12, [1, 0, 0, 0, 0]),
            input("mu", 1.81e-5, [-1, 1, -1, 0, 0]),
        ];
        let basis = pi_groups(&inputs).expect("pi");
        assert_eq!(basis.rank, 3);
        assert_eq!(basis.groups.len(), 2);
    }

    #[test]
    fn degenerate_inputs_refuse() {
        assert!(pi_groups(&[]).is_err());
        assert!(pi_groups(&[input("bad", -1.0, [1, 0, 0, 0, 0])]).is_err());
    }
}
