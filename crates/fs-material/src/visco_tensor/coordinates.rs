//! Exact coordinate sharing for equally timed, displacement-driven histories.
//!
//! For each BIT-IDENTICAL time, A=Q R with Q^T Q=I. The original equations
//! z'=(A q-z)/tau leave z=Q y invariant and become y'=(R q-y)/tau.
//! Thus force, |A q-z|^2/2 and its heat are unchanged. Every original history
//! is reconstructed, including histories that have no direct mechanical force.
//! No rank threshold, pole merging, fitting or constitutive law is introduced.
//!
//! This coordinate image covers relaxed z=A q and unrelaxed z=0 initial states,
//! and histories already supplied IN these coordinates. It does not ingest or
//! project independently authored original histories outside that subspace.
use super::{Result, invalid, poll};
use core::ops::Range;
use std::{collections::BTreeMap, sync::Arc};

/// Fixed work and evolved-coordinate capacities; neither causes truncation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RelaxationCoordinateBudget {
    pub maximum_states: usize,
    pub maximum_scalar_terms: usize,
}

#[derive(Clone, Debug)]
struct MemoryRow {
    coordinates: Range<usize>,
    coefficients: Vec<f64>,
}

#[derive(Clone, Debug)]
struct MemoryTransform {
    coordinates: usize,
    rows: Vec<MemoryRow>,
}

/// Immutable view of the inverse coordinate map for original material points.
/// Clones and point slices share the same complete transform allocation.
#[derive(Clone, Debug)]
pub struct RelaxationMemoryMap {
    transform: Arc<MemoryTransform>,
    original: Range<usize>,
}
impl RelaxationMemoryMap {
    pub fn coordinates(&self) -> usize {
        self.transform.coordinates
    }
    pub fn original_coordinates(&self) -> usize {
        self.original.len()
    }

    /// Recover the complete original history in its original order.
    pub fn reconstruct(&self, memory: &[f64]) -> Result<Vec<f64>> {
        if memory.len() != self.coordinates() || memory.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "shared material history requires its complete finite coordinates",
            ));
        }
        let values: Vec<_> = self.transform.rows[self.original.clone()]
            .iter()
            .map(|row| {
                row.coefficients
                    .iter()
                    .zip(&memory[row.coordinates.clone()])
                    .fold(0.0, |sum, (a, z)| sum + a * z)
            })
            .collect();
        if values.iter().any(|v| !v.is_finite()) {
            return Err(invalid("original material-history reconstruction overflow"));
        }
        Ok(values)
    }
}

/// Orthogonal history coordinates and the original-point reconstruction map.
#[derive(Clone, Debug)]
pub struct RelaxationCoordinates {
    rows: Vec<Vec<f64>>,
    times: Vec<f64>,
    transform: Arc<MemoryTransform>,
}
impl RelaxationCoordinates {
    /// Retain at most `displacements` coordinates for each exact time family.
    /// All QR columns survive, including zero and rank-deficient directions.
    /// Smaller families retain their original rows without QR. Source order is
    /// retained by the inverse map; families use deterministic positive-time order.
    pub fn new(
        rows: &[Vec<f64>],
        times: &[f64],
        displacements: usize,
        budget: RelaxationCoordinateBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let n = displacements;
        if !(1..=256).contains(&n) || rows.len() != times.len() || budget.maximum_scalar_terms == 0
        {
            return Err(invalid(
                "shared material coordinates require complete rows, times and explicit work",
            ));
        }
        let mut groups = BTreeMap::<u64, Vec<usize>>::new();
        for (i, (row, &time)) in rows.iter().zip(times).enumerate() {
            poll(&mut cancelled)?;
            if row.len() != n
                || row.iter().any(|v| !v.is_finite())
                || !time.is_finite()
                || time <= 0.0
                || !time.recip().is_finite()
                || time.recip() <= 0.0
            {
                return Err(invalid(
                    "shared material coordinates require finite factors and positive representable times",
                ));
            }
            groups.entry(time.to_bits()).or_default().push(i);
        }
        let mut count = 0usize;
        let mut work = 0usize;
        for group in groups.values() {
            let m = group.len();
            let k = m.min(n);
            count = count
                .checked_add(k)
                .ok_or_else(|| invalid("shared material extent overflow"))?;
            // QR, thin Q, orthogonality and original-entry reconstruction,
            // plus factor scaling and the retained inverse-map storage.
            let terms = if m > n {
                m.checked_mul(n)
                    .and_then(|v| v.checked_mul(n))
                    .and_then(|v| v.checked_mul(16))
            } else {
                m.checked_mul(n)
                    .and_then(|v| m.checked_mul(m).and_then(|w| v.checked_add(w)))
            };
            work = terms
                .and_then(|v| work.checked_add(v))
                .ok_or_else(|| invalid("shared material work overflow"))?;
        }
        if count > budget.maximum_states || work > budget.maximum_scalar_terms {
            return Err(invalid(
                "complete shared material coordinates exceed explicit state or work capacity",
            ));
        }
        let mut factors = Vec::with_capacity(count);
        let mut rates = Vec::with_capacity(count);
        let mut map = vec![
            MemoryRow {
                coordinates: 0..0,
                coefficients: Vec::new()
            };
            rows.len()
        ];
        for (time, group) in groups {
            poll(&mut cancelled)?;
            let start = factors.len();
            let m = group.len();
            let k = m.min(n);
            let range = start..start + k;
            if m <= n {
                for (local, &source) in group.iter().enumerate() {
                    factors.push(rows[source].clone());
                    let mut coefficients = vec![0.0; k];
                    coefficients[local] = 1.0;
                    map[source] = MemoryRow {
                        coordinates: range.clone(),
                        coefficients,
                    };
                }
            } else {
                let mut scale = vec![0.0_f64; n];
                for &source in &group {
                    for j in 0..n {
                        scale[j] = scale[j].max(rows[source][j].abs());
                    }
                }
                let mut a = vec![0.0; m * n];
                for (i, &source) in group.iter().enumerate() {
                    poll(&mut cancelled)?;
                    for j in 0..n {
                        if scale[j] != 0.0 {
                            a[i * n + j] = rows[source][j] / scale[j];
                            if rows[source][j] != 0.0 && a[i * n + j] == 0.0 {
                                return Err(invalid(
                                    "material coordinate scaling lost an original factor entry",
                                ));
                            }
                        }
                    }
                }
                let qr = fs_la::factor::qr(&a, m, n);
                poll(&mut cancelled)?;
                let mut q = vec![vec![0.0; n]; m];
                for column in 0..n {
                    poll(&mut cancelled)?;
                    let mut unit = vec![0.0; m];
                    unit[column] = 1.0;
                    qr.apply_q(&mut unit);
                    for i in 0..m {
                        q[i][column] = unit[i];
                    }
                }
                let allowance = 256.0 * f64::EPSILON * (m + n) as f64;
                for i in 0..n {
                    for j in 0..=i {
                        let actual = (0..m).fold(0.0, |s, l| s + q[l][i] * q[l][j]);
                        if !actual.is_finite()
                            || (actual - if i == j { 1.0 } else { 0.0 }).abs() > allowance
                        {
                            return Err(invalid(
                                "shared material coordinates failed original orthogonality",
                            ));
                        }
                    }
                }
                let norms: Vec<_> = (0..n)
                    .map(|j| {
                        fs_math::det::sqrt((0..m).fold(0.0, |s, i| s + a[i * n + j] * a[i * n + j]))
                    })
                    .collect();
                for i in 0..m {
                    poll(&mut cancelled)?;
                    for j in 0..n {
                        let actual = (0..n).fold(0.0, |s, l| s + q[i][l] * qr.r(l, j));
                        if !actual.is_finite()
                            || (actual - a[i * n + j]).abs() > allowance * norms[j]
                        {
                            return Err(invalid(
                                "shared material coordinates failed original factor reconstruction",
                            ));
                        }
                    }
                    map[group[i]] = MemoryRow {
                        coordinates: range.clone(),
                        coefficients: q[i].clone(),
                    };
                }
                for i in 0..n {
                    let row: Vec<_> = (0..n).map(|j| qr.r(i, j) * scale[j]).collect();
                    if row.iter().enumerate().any(|(j, v)| {
                        !v.is_finite() || qr.r(i, j) != 0.0 && scale[j] != 0.0 && *v == 0.0
                    }) {
                        return Err(invalid("shared material factor scaling is unrepresentable"));
                    }
                    factors.push(row);
                }
            }
            rates.extend(core::iter::repeat_n(f64::from_bits(time), k));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            rows: factors,
            times: rates,
            transform: Arc::new(MemoryTransform {
                coordinates: count,
                rows: map,
            }),
        })
    }

    pub fn coupling_rows(&self) -> &[Vec<f64>] {
        &self.rows
    }
    pub fn relaxation_times_s(&self) -> &[f64] {
        &self.times
    }
    pub fn original_coordinates(&self) -> usize {
        self.transform.rows.len()
    }
    pub fn coordinates(&self) -> usize {
        self.rows.len()
    }

    /// Select original history rows without losing the shared inverse map.
    pub fn memory_map(&self, original: Range<usize>) -> Result<RelaxationMemoryMap> {
        if original.start > original.end || original.end > self.original_coordinates() {
            return Err(invalid(
                "original material history range is outside its source",
            ));
        }
        Ok(RelaxationMemoryMap {
            transform: Arc::clone(&self.transform),
            original,
        })
    }
}

#[cfg(test)]
mod tests;
