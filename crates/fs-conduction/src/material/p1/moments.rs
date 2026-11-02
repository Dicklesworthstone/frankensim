//! Scalar moments of an admitted piecewise-linear law over a P1 tetrahedron.
//!
//! For uniform barycentric coordinates, the density of `T = sum(lambda_i T_i)`
//! is `3 B_2(T) / (max(T_i) - min(T_i))`, where B_2 has the four sorted nodal
//! temperatures as knots. Weighting by lambda_i duplicates knot T_i: the
//! weighted density is `B_3(T) / (max(T_i) - min(T_i))`, because multiplying a
//! Dirichlet(1,1,1,1) density by lambda_i gives one quarter of the density with
//! its i-th parameter raised to two. The duplicate knot is that integer
//! parameter. For the simplex interpretation and spline normalization, see
//! de Boor, "Splines as linear combinations of B-splines. A survey" (1976),
//! eq. (4.4), Lemma 9.1 and eq. (9.11), https://ftp.cs.wisc.edu/Approx/survey76.pdf.
//!
//! Thus the average is the integral of k times a quadratic B-spline, and its
//! i-th derivative is the integral of k' times the cubic duplicated-knot
//! B-spline. Two-point Gauss is exact on every interval cut by a nodal or
//! material knot: both integrands have degree at most three. Exact here refers
//! to the polynomial rule, not a floating-point or spatial error certificate.
//! No geometric subdivision, heap allocation or ill-conditioned divided
//! differences are required. Repeated temperatures are ordinary spline knots.

#[derive(Debug, Clone, Copy)]
pub(crate) struct ScalarMoments {
    pub(crate) average: f64,
    pub(crate) derivatives: [f64; 4],
}

#[derive(Debug)]
pub(crate) enum MomentError<E> {
    Numerical(&'static str),
    Interrupted(E),
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }

    fn total(self) -> f64 {
        self.value + self.correction
    }
}

fn checked<E>(average: f64, derivatives: [f64; 4]) -> Result<ScalarMoments, MomentError<E>> {
    if !average.is_finite() || average <= 0.0 || derivatives.iter().any(|d| !d.is_finite()) {
        return Err(MomentError::Numerical(
            "nonfinite or nonpositive integrated conductivity",
        ));
    }
    Ok(ScalarMoments {
        average,
        derivatives,
    })
}

// An interior Gauss point is represented by its interval and fractional
// coordinate. Forming its absolute temperature can round to an endpoint when
// two nodal temperatures are adjacent floats; this representation cannot.
fn bspline(knots: &[f64], lo: f64, hi: f64, alpha: f64) -> f64 {
    let mut basis = [0.0; 4];
    for i in 0..knots.len() - 1 {
        if knots[i] < knots[i + 1] && knots[i] <= lo && hi <= knots[i + 1] {
            basis[i] = 1.0;
        }
    }
    for degree in 1..knots.len() - 1 {
        for i in 0..knots.len() - degree - 1 {
            let left = if basis[i] != 0.0 {
                let distance = (1.0 - alpha).mul_add(lo - knots[i], alpha * (hi - knots[i]));
                distance / (knots[i + degree] - knots[i]) * basis[i]
            } else {
                0.0
            };
            let right = if basis[i + 1] != 0.0 {
                let end = knots[i + degree + 1];
                let distance = (1.0 - alpha).mul_add(end - lo, alpha * (end - hi));
                distance / (end - knots[i + 1]) * basis[i + 1]
            } else {
                0.0
            };
            basis[i] = left + right;
        }
    }
    basis[0]
}

fn segment(knots: &[(f64, f64)], temperature: f64) -> usize {
    knots
        .partition_point(|&(t, _)| t <= temperature)
        .saturating_sub(1)
        .min(knots.len() - 2)
}

/// Inputs are a validated positive table and four supported finite nodal
/// temperatures. Cancellation is checked every 128 crossed intervals. A
/// uniform field at a material kink uses the table's right-segment derivative;
/// this is a semismooth choice, not a unique two-sided classical derivative.
pub(crate) fn integrate_table<E>(
    nodal: [f64; 4],
    knots: &[(f64, f64)],
    with_derivatives: bool,
    mut poll: impl FnMut(usize) -> Result<(), E>,
) -> Result<ScalarMoments, MomentError<E>> {
    if knots.is_empty() || nodal.iter().any(|t| !t.is_finite()) {
        return Err(MomentError::Numerical(
            "invalid P1 conductivity integration input",
        ));
    }
    if knots.len() == 1 {
        return checked(knots[0].1, [0.0; 4]);
    }
    let mut sorted = nodal;
    sorted.sort_by(f64::total_cmp);
    let low = sorted[0];
    let high = sorted[3];
    if low < knots[0].0 || high > knots[knots.len() - 1].0 {
        return Err(MomentError::Numerical(
            "unsupported P1 conductivity integration input",
        ));
    }
    let mut piece = segment(knots, low);
    if high <= knots[piece + 1].0 {
        // Constant/uniform/one-segment fast path retains the production mean
        // interpolation arithmetic. Quartering first is only needed if the
        // sum of four finite, supported temperatures overflows.
        let sum: f64 = nodal.iter().sum();
        let mean = if sum.is_finite() {
            sum / 4.0
        } else {
            nodal.iter().map(|t| t / 4.0).sum()
        };
        let (t0, k0) = knots[piece];
        let (t1, k1) = knots[piece + 1];
        let average = ((mean - t0) / (t1 - t0)).mul_add(k1 - k0, k0);
        let slope = if with_derivatives {
            (k1 - k0) / (t1 - t0) / 4.0
        } else {
            0.0
        };
        return checked(average, [slope; 4]);
    }
    let range = high - low;
    if !range.is_finite() || range <= 0.0 {
        return Err(MomentError::Numerical("invalid P1 temperature range"));
    }
    let mut duplicated = [[0.0; 5]; 4];
    if with_derivatives {
        for (i, row) in duplicated.iter_mut().enumerate() {
            row[..4].copy_from_slice(&sorted);
            row[4] = nodal[i];
            row.sort_by(f64::total_cmp);
        }
    }
    let mut sums = [Sum::default(); 5];
    let mut lo = low;
    let mut node = 0usize;
    let mut interval = 0usize;
    let mut flat = true;
    let first_value = knots[piece].1;
    const ALPHAS: [f64; 2] = [0.211_324_865_405_187_13, 0.788_675_134_594_812_9];
    while lo < high {
        if interval % 128 == 0 {
            poll(interval).map_err(MomentError::Interrupted)?;
        }
        while node < 4 && sorted[node] <= lo {
            node += 1;
        }
        let (t0, k0) = knots[piece];
        let (t1, k1) = knots[piece + 1];
        let hi = high.min(t1).min(sorted[node]);
        if hi <= lo || !((t1 - t0).is_finite() && t1 > t0) {
            return Err(MomentError::Numerical("invalid P1 material interval"));
        }
        flat &= k0 == first_value && k1 == first_value;
        let weight = 0.5 * ((hi - lo) / range);
        let slope = if with_derivatives {
            (k1 - k0) / (t1 - t0)
        } else {
            0.0
        };
        for alpha in ALPHAS {
            let fraction =
                (1.0 - alpha).mul_add((lo - t0) / (t1 - t0), alpha * ((hi - t0) / (t1 - t0)));
            // Interpolate from the nearer endpoint to preserve a small
            // positive k in a strongly decreasing, high-contrast source law.
            let k = if fraction <= 0.5 {
                fraction.mul_add(k1 - k0, k0)
            } else {
                (1.0 - fraction).mul_add(k0 - k1, k1)
            };
            sums[0].add(k * (3.0 * weight * bspline(&sorted, lo, hi, alpha)));
            if with_derivatives {
                for (i, row) in duplicated.iter().enumerate() {
                    sums[i + 1].add(slope * (weight * bspline(row, lo, hi, alpha)));
                }
            }
        }
        lo = hi;
        if lo == t1 && lo < high {
            piece += 1;
        }
        interval += 1;
    }
    if flat {
        return checked(first_value, [0.0; 4]);
    }
    checked(
        sums[0].total(),
        std::array::from_fn(|i| sums[i + 1].total()),
    )
}
