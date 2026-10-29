//! Exact equal-time Maxwell compression in the original factor coordinates.
//! Orthogonal row reduction preserves A^T A, all work, storage and heat for
//! relaxed/unrelaxed histories driven by the same displacement. It does not
//! reconstruct arbitrary independently authored material-point memories.

use std::collections::BTreeMap;

/// One original relaxation time with at most n complete coupled factor rows.
/// Rows act on physical metre-valued coordinates and have units sqrt(N/m).
#[derive(Clone, Debug, PartialEq)]
pub struct HereditaryBranch {
    relaxation_time_s: f64,
    rows: Vec<Vec<f64>>,
}
impl HereditaryBranch {
    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }
    pub fn rows(&self) -> &[Vec<f64>] {
        &self.rows
    }
}

#[derive(Debug, PartialEq)]
pub(super) enum CompressionError {
    Invalid(&'static str),
    Cancelled,
}

/// Group only bit-identical positive source times. Scaled QR never applies a
/// relative rank threshold: a weak direction remains a physical memory channel.
pub(super) fn compact(
    sources: &[(f64, &[Vec<f64>])],
    n: usize,
    maximum_scalar_terms: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<HereditaryBranch>, CompressionError> {
    use CompressionError::{Cancelled, Invalid};
    if cancelled() {
        return Err(Cancelled);
    }
    if !(1..=256).contains(&n) {
        return Err(Invalid(
            "hereditary compression requires 1..=256 coordinates",
        ));
    }
    let mut families = BTreeMap::<u64, Vec<&[f64]>>::new();
    for &(tau, rows) in sources {
        if cancelled() {
            return Err(Cancelled);
        }
        if !tau.is_finite() || tau <= 0.0 {
            return Err(Invalid(
                "hereditary compression requires positive finite source times",
            ));
        }
        for row in rows {
            if row.len() != n || row.iter().any(|x| !x.is_finite()) {
                return Err(Invalid(
                    "hereditary compression requires complete finite factor rows",
                ));
            }
            if row.iter().any(|x| *x != 0.0) {
                families.entry(tau.to_bits()).or_default().push(row);
            }
        }
    }
    // QR, its original-entry reconstruction and the complete n*n output are
    // bounded before allocating a factor workspace. No source channel is cut
    // to satisfy the caller's budget.
    let mut work = 0usize;
    for rows in families.values() {
        work = rows
            .len()
            .max(n)
            .checked_mul(n)
            .and_then(|x| x.checked_mul(n))
            .and_then(|x| x.checked_mul(8))
            .and_then(|x| work.checked_add(x))
            .ok_or(Invalid("hereditary factor work overflow"))?;
    }
    if work > maximum_scalar_terms {
        return Err(Invalid(
            "hereditary factor work exceeds its explicit budget",
        ));
    }
    let mut result = Vec::with_capacity(families.len());
    for (tau, rows) in families {
        if cancelled() {
            return Err(Cancelled);
        }
        let m = rows.len().max(n);
        let mut scale = vec![0.0_f64; n];
        for row in &rows {
            for j in 0..n {
                scale[j] = scale[j].max(row[j].abs());
            }
        }
        let mut a = vec![0.0; m * n];
        for (i, row) in rows.iter().enumerate() {
            if cancelled() {
                return Err(Cancelled);
            }
            for j in 0..n {
                if scale[j] != 0.0 {
                    a[i * n + j] = row[j] / scale[j];
                    if row[j] != 0.0 && a[i * n + j] == 0.0 {
                        return Err(Invalid(
                            "hereditary factor equilibration lost a nonzero entry",
                        ));
                    }
                }
            }
        }
        let qr = fs_la::factor::qr(&a, m, n);
        let allowance = 256.0 * f64::EPSILON * (m + n) as f64;
        for i in 0..n {
            if cancelled() {
                return Err(Cancelled);
            }
            let left = (0..m)
                .map(|k| a[k * n + i] * a[k * n + i])
                .sum::<f64>()
                .sqrt();
            for j in 0..=i {
                let original = (0..m).map(|k| a[k * n + i] * a[k * n + j]).sum::<f64>();
                let actual = (0..n).map(|k| qr.r(k, i) * qr.r(k, j)).sum::<f64>();
                let right = (0..m)
                    .map(|k| a[k * n + j] * a[k * n + j])
                    .sum::<f64>()
                    .sqrt();
                if !actual.is_finite() || (actual - original).abs() > allowance * left * right {
                    return Err(Invalid(
                        "hereditary factor failed its original Gram reconstruction",
                    ));
                }
            }
        }
        let mut compact_rows = Vec::with_capacity(n);
        for i in 0..n {
            let row = (0..n).map(|j| qr.r(i, j) * scale[j]).collect::<Vec<_>>();
            if row
                .iter()
                .enumerate()
                .any(|(j, x)| !x.is_finite() || (qr.r(i, j) != 0.0 && scale[j] != 0.0 && *x == 0.0))
            {
                return Err(Invalid("hereditary factor scale is unrepresentable"));
            }
            if row.iter().any(|x| *x != 0.0) {
                compact_rows.push(row);
            }
        }
        result.push(HereditaryBranch {
            relaxation_time_s: f64::from_bits(tau),
            rows: compact_rows,
        });
    }
    if cancelled() {
        return Err(Cancelled);
    }
    Ok(result)
}

#[cfg(test)]
#[path = "memory/tests.rs"]
mod tests;
