//! Shared bounded mechanical-response admission and supported static reaction.
//!
//! This is the original pivoted factor from the two-axis tangential owner.
//! Factors serve numerical admission/seed selection only. Constitutive owners
//! retain their original physical matrix and independently check compatibility,
//! work and native laws. No regularizing compliance is ever introduced.

use fs_math::det;

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ContactResponseError {
    Invalid(&'static str),
    Cancelled,
}

// Complete-pivot Cholesky admission, normalized by the original largest entry.
// A discarded tail has infinity norm at most the explicit arithmetic tolerance;
// the original matrix, not this factor, is used in every physical residual.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PositiveFactor {
    dimension: usize,
    scale: f64,
    rank: usize,
    permutation: Vec<usize>,
    lower: Vec<f64>,
    gram_lower: Vec<f64>,
}

impl PositiveFactor {
    pub(crate) const fn rank(&self) -> usize {
        self.rank
    }

    pub(crate) fn new(
        matrix: &[f64],
        n: usize,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ContactResponseError> {
        let mut scale = 0.0_f64;
        for row in matrix.chunks(n) {
            poll(cancelled)?;
            for &value in row {
                scale = scale.max(number(value)?.abs());
            }
        }
        let mut result = Self {
            dimension: n,
            scale,
            rank: 0,
            permutation: (0..n).collect(),
            lower: vec![0.0; n * n],
            gram_lower: Vec::new(),
        };
        if scale == 0.0 {
            return Ok(result);
        }
        let mut work = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                let a = matrix[i * n + j] / scale;
                let b = matrix[j * n + i] / scale;
                if (a - b).abs() > tolerance {
                    return Err(ContactResponseError::Invalid(
                        "nonsymmetric full compliance",
                    ));
                }
                work[i * n + j] = 0.5 * a + 0.5 * b;
            }
        }
        for k in 0..n {
            let mut tail_norm = 0.0_f64;
            let mut pivot = k;
            for i in k..n {
                poll(cancelled)?;
                let mut sum = 0.0;
                for j in k..n {
                    sum = number(sum + work[i * n + j].abs())?;
                }
                tail_norm = tail_norm.max(sum);
                if work[i * n + i] > work[pivot * n + pivot] {
                    pivot = i;
                }
            }
            if tail_norm <= tolerance {
                break;
            }
            if work[pivot * n + pivot] <= 0.0 {
                return Err(ContactResponseError::Invalid(
                    "nonpositive-semidefinite full compliance",
                ));
            }
            if pivot != k {
                for j in 0..n {
                    work.swap(k * n + j, pivot * n + j);
                }
                for i in 0..n {
                    work.swap(i * n + k, i * n + pivot);
                }
                for j in 0..k {
                    result.lower.swap(k * n + j, pivot * n + j);
                }
                result.permutation.swap(k, pivot);
            }
            let diagonal = number(det::sqrt(work[k * n + k]))?;
            result.lower[k * n + k] = diagonal;
            for i in k + 1..n {
                result.lower[i * n + k] = number(work[i * n + k] / diagonal)?;
            }
            for i in k + 1..n {
                poll(cancelled)?;
                for j in k + 1..=i {
                    let value = number(
                        work[i * n + j] - result.lower[i * n + k] * result.lower[j * n + k],
                    )?;
                    work[i * n + j] = value;
                    work[j * n + i] = value;
                }
            }
            result.rank += 1;
        }
        // The Schur scratch is no longer needed. Release it before allocating
        // the range Gram factor so preparation stays within three dense arrays
        // including the caller's retained physical compliance.
        drop(work);
        let r = result.rank;
        if r > 0 && r < n {
            result.gram_lower = vec![0.0; r * r];
            for i in 0..r {
                for j in 0..=i {
                    poll(cancelled)?;
                    let mut value = 0.0;
                    for row in 0..n {
                        value =
                            number(value + result.lower[row * n + i] * result.lower[row * n + j])?;
                    }
                    for k in 0..j {
                        value = number(
                            value - result.gram_lower[i * r + k] * result.gram_lower[j * r + k],
                        )?;
                    }
                    result.gram_lower[i * r + j] = if i == j {
                        if value <= 0.0 {
                            return Err(ContactResponseError::Invalid(
                                "unresolved supported compliance factor",
                            ));
                        }
                        det::sqrt(value)
                    } else {
                        number(value / result.gram_lower[j * r + j])?
                    };
                }
            }
        }
        poll(cancelled)?;
        Ok(result)
    }

    pub(crate) fn reaction(
        &self,
        released: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, ContactResponseError> {
        let n = self.dimension;
        let r = self.rank;
        let mut result = vec![0.0; n];
        if r == 0 {
            return Ok(result);
        }
        let mut y = vec![0.0; r];
        for i in 0..r {
            poll(cancelled)?;
            let mut value = released[self.permutation[i]];
            for j in 0..i {
                value -= self.lower[i * n + j] * y[j];
            }
            y[i] = value / self.lower[i * n + i];
        }
        if r == n {
            for i in (0..r).rev() {
                poll(cancelled)?;
                let mut value = y[i];
                for j in i + 1..r {
                    value -= self.lower[j * n + i] * y[j];
                }
                y[i] = value / self.lower[i * n + i];
            }
            for i in 0..n {
                result[self.permutation[i]] = y[i] / self.scale;
            }
        } else {
            // L y = u gives the supported displacement. Minimum-norm force is
            // L (L^T L)^-1 y / scale; no arbitrary independent rows are forced.
            for i in 0..r {
                poll(cancelled)?;
                let mut value = y[i];
                for j in 0..i {
                    value -= self.gram_lower[i * r + j] * y[j];
                }
                y[i] = value / self.gram_lower[i * r + i];
            }
            for i in (0..r).rev() {
                poll(cancelled)?;
                let mut value = y[i];
                for j in i + 1..r {
                    value -= self.gram_lower[j * r + i] * y[j];
                }
                y[i] = value / self.gram_lower[i * r + i];
            }
            for row in 0..n {
                poll(cancelled)?;
                let mut value = 0.0;
                for j in 0..r {
                    value += self.lower[row * n + j] * y[j];
                }
                result[self.permutation[row]] = value / self.scale;
            }
        }
        // An unrepresentable numerical seed need not rule out a bounded kinetic
        // solution. It never changes the physical matrix or accepted history.
        if result.iter().any(|x| !x.is_finite()) {
            result.fill(0.0);
        }
        Ok(result)
    }
}

fn number(value: f64) -> Result<f64, ContactResponseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ContactResponseError::Invalid(
            "nonfinite multiple-contact arithmetic",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ContactResponseError> {
    if cancelled() {
        Err(ContactResponseError::Cancelled)
    } else {
        Ok(())
    }
}
