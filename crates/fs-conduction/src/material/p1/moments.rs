//! Scalar moments of an admitted piecewise-linear law over a P1 simplex.
//!
//! The coordinate u may be temperature, specific enthalpy or another admitted
//! scalar. For uniform barycentric coordinates, the density of
//! `u = sum(lambda_i u_i)` is `3 B_2(u) / (max(u_i) - min(u_i))`, where B_2 has
//! the four sorted nodal coordinates as knots. Weighting by lambda_i duplicates
//! knot u_i: the weighted density is `B_3(u) / (max(u_i) - min(u_i))`, because multiplying a
//! Dirichlet(1,1,1,1) density by lambda_i gives one quarter of the density with
//! its i-th parameter raised to two. The duplicate knot is that integer
//! parameter. For the simplex interpretation and spline normalization, see
//! de Boor, "Splines as linear combinations of B-splines. A survey" (1976),
//! eq. (4.4), Lemma 9.1 and eq. (9.11), https://ftp.cs.wisc.edu/Approx/survey76.pdf.
//!
//! A triangle uses three nodal knots, density `2 B_1 / range`, and one
//! duplicated knot for each weighted derivative. The same rule therefore
//! integrates both face and volume laws, without treating a triangle as a
//! tetrahedron with an artificially repeated vertex.
//!
//! Thus the tetrahedral average is the integral of k times a quadratic B-spline, and its
//! i-th derivative is the integral of k' times the cubic duplicated-knot
//! B-spline. Two-point Gauss is exact on every interval cut by a nodal or
//! material knot: both integrands have degree at most three. Exact here refers
//! to the polynomial rule, not a floating-point or spatial error certificate.
//! No geometric subdivision, heap allocation or ill-conditioned divided
//! differences are required. Repeated coordinates are ordinary spline knots.

#[path = "moments/caloric.rs"]
mod caloric;
pub(crate) use caloric::integrate_caloric_table;

#[derive(Debug, Clone, Copy)]
pub(crate) struct ScalarMoments<const N: usize = 4> {
    pub(crate) average: f64,
    pub(crate) derivatives: [f64; N],
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

fn checked<E, const N: usize>(
    average: f64,
    derivatives: [f64; N],
) -> Result<ScalarMoments<N>, MomentError<E>> {
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
// coordinate. Forming its absolute value can round to an endpoint when two
// nodal coordinates are adjacent floats; this representation cannot.
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

fn segment(knots: &[(f64, f64)], coordinate: f64) -> usize {
    knots
        .partition_point(|&(x, _)| x <= coordinate)
        .saturating_sub(1)
        .min(knots.len() - 2)
}

/// Inputs are a validated positive table and three or four supported finite nodal
/// coordinates. Cancellation is checked every 128 crossed intervals. A
/// uniform field at a material kink uses the table's right-segment derivative;
/// this is a semismooth choice, not a unique two-sided classical derivative.
pub(crate) fn integrate_table<E, const N: usize>(
    nodal: [f64; N],
    knots: &[(f64, f64)],
    with_derivatives: bool,
    mut poll: impl FnMut(usize) -> Result<(), E>,
) -> Result<ScalarMoments<N>, MomentError<E>> {
    if !(N == 3 || N == 4) || knots.is_empty() || nodal.iter().any(|t| !t.is_finite()) {
        return Err(MomentError::Numerical(
            "invalid P1 conductivity integration input",
        ));
    }
    if knots.len() == 1 {
        return checked(knots[0].1, [0.0; N]);
    }
    let mut sorted = nodal;
    sorted.sort_by(f64::total_cmp);
    let low = sorted[0];
    let high = sorted[N - 1];
    if low < knots[0].0 || high > knots[knots.len() - 1].0 {
        return Err(MomentError::Numerical(
            "unsupported P1 conductivity integration input",
        ));
    }
    let mut piece = segment(knots, low);
    if high <= knots[piece + 1].0 {
        // Constant/uniform/one-segment fast path retains the production mean
        // interpolation arithmetic. Dividing first is only needed if the
        // sum of the finite, supported coordinates overflows.
        let sum: f64 = nodal.iter().sum();
        let mean = if sum.is_finite() {
            sum / N as f64
        } else {
            nodal.iter().map(|t| t / N as f64).sum()
        };
        let (t0, k0) = knots[piece];
        let (t1, k1) = knots[piece + 1];
        // Preserve exact source endpoints and interpolate from the nearer
        // end, as in the multi-interval rule below. Starting at k0 can erase
        // a small positive k1 when a descending source has strong contrast.
        let interpolate = |t: f64| {
            let fraction = (t - t0) / (t1 - t0);
            if t == t0 {
                k0
            } else if t == t1 {
                k1
            } else if fraction <= 0.5 {
                fraction.mul_add(k1 - k0, k0)
            } else {
                (1.0 - fraction).mul_add(k0 - k1, k1)
            }
        };
        let average = if N == 3 {
            // Dividing a temperature by three can move even a uniform field
            // outside its source endpoint. The mean of affine k at the three
            // ORIGINAL coordinates is the exact face moment and also retains
            // high-contrast laws across adjacent representable temperatures.
            if low == high {
                interpolate(low)
            } else {
                let values = nodal.map(interpolate);
                let scale = values.iter().copied().fold(0.0_f64, f64::max);
                scale * (values.iter().map(|k| k / scale).sum::<f64>() / N as f64)
            }
        } else {
            interpolate(mean)
        };
        let slope = if with_derivatives {
            (k1 - k0) / (t1 - t0) / N as f64
        } else {
            0.0
        };
        return checked(average, [slope; N]);
    }
    let range = high - low;
    if !range.is_finite() || range <= 0.0 {
        return Err(MomentError::Numerical("invalid P1 scalar coordinate range"));
    }
    let mut duplicated = [[0.0; 5]; N];
    if with_derivatives {
        for (i, row) in duplicated.iter_mut().enumerate() {
            row[..N].copy_from_slice(&sorted);
            row[N] = nodal[i];
            row[..=N].sort_by(f64::total_cmp);
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
        while node < N && sorted[node] <= lo {
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
            sums[0].add(k * ((N - 1) as f64 * weight * bspline(&sorted, lo, hi, alpha)));
            if with_derivatives {
                for (i, row) in duplicated.iter().enumerate() {
                    sums[i + 1].add(slope * (weight * bspline(&row[..=N], lo, hi, alpha)));
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
        return checked(first_value, [0.0; N]);
    }
    checked(
        sums[0].total(),
        std::array::from_fn(|i| sums[i + 1].total()),
    )
}

#[cfg(test)]
mod simplex_tests {
    use super::*;

    fn moments<const N: usize>(nodes: [f64; N]) -> ScalarMoments<N> {
        integrate_table(
            nodes,
            &[(300.0, 1.0), (350.0, 2.0), (400.0, 4.0)],
            true,
            |_| Ok::<_, ()>(()),
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-14 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn triangle_kink_moments_match_independent_subtriangle_integrals() {
        // Split the triangle along T=350 and integrate affine k and each
        // original barycentric coordinate on the two resulting polygons.
        let value = moments([300.0, 350.0, 400.0]);
        close(value.average, 13.0 / 6.0);
        for (actual, expected) in
            value
                .derivatives
                .into_iter()
                .zip([1.0 / 120.0, 1.0 / 100.0, 7.0 / 600.0])
        {
            close(actual, expected);
        }
        let repeated = moments([300.0, 300.0, 400.0]);
        close(repeated.average, 7.0 / 4.0);
        for (actual, expected) in
            repeated
                .derivatives
                .into_iter()
                .zip([3.0 / 400.0, 3.0 / 400.0, 1.0 / 100.0])
        {
            close(actual, expected);
        }
    }

    #[test]
    fn triangle_decimal_endpoint_preserves_positive_high_contrast_source() {
        let table = [(350.0, 1e20), (350.1, 1e-20)];
        let uniform = integrate_table([350.1; 3], &table, true, |_| Ok::<_, ()>(())).unwrap();
        assert_eq!(uniform.average, 1e-20);
        let adjacent = 350.1_f64.next_down();
        let endpoint_pair = [(adjacent, 1e20), (350.1, 1e-20)];
        let mixed = integrate_table([adjacent, 350.1, 350.1], &endpoint_pair, true, |_| {
            Ok::<_, ()>(())
        })
        .unwrap();
        close(mixed.average, 1e20 / 3.0);
    }

    #[test]
    fn tetrahedral_moments_keep_the_original_four_node_density() {
        // For s=lambda_3 the uniform tetrahedral density is 3(1-s)^2.
        // k=1+2s+2(s-1/2)+, giving this independent beta-polynomial integral.
        let value = moments([300.0, 300.0, 300.0, 400.0]);
        close(value.average, 49.0 / 32.0);
        for (actual, expected) in value.derivatives.into_iter().zip([
            17.0 / 3200.0,
            17.0 / 3200.0,
            17.0 / 3200.0,
            21.0 / 3200.0,
        ]) {
            close(actual, expected);
        }
    }
}
