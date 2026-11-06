//! Exact original k(T(h)) moments for affine Cp(T) sensible intervals.
//!
//! On each source/nodal cut, h(T) is quadratic. The B2(h) value density,
//! multiplied by affine k(T) and Cp(T), has degree six in T. For a nodal
//! derivative, dk/dh * dh = dk/dT * dT cancels Cp, leaving B3(h(T)) of degree
//! six. Four-point Gauss therefore integrates both original responses exactly
//! up to arithmetic rounding. Authored/latent affine h intervals retain the
//! existing moment owner, including its fast path when no nonlinear interval
//! is crossed. The source and nodal knot cuts remain the only subdivisions.

use super::{MomentError, ScalarMoments, Sum, bspline, checked, integrate_table, segment};

fn temperature_fraction(capacities: [f64; 2], energy_fraction: f64) -> (f64, f64) {
    let scale = capacities[0].max(capacities[1]);
    let a = capacities[0] / scale;
    let b = capacities[1] / scale;
    let capacity = ((1.0 - energy_fraction) * a * a + energy_fraction * b * b).sqrt();
    (energy_fraction * (a + b) / (a + capacity), capacity)
}

fn interpolate(a: f64, b: f64, fraction: f64) -> f64 {
    if fraction == 0.0 {
        a
    } else if fraction == 1.0 {
        b
    } else if fraction <= 0.5 {
        fraction.mul_add(b - a, a)
    } else {
        (1.0 - fraction).mul_add(a - b, b)
    }
}

/// Capacities correspond to the endpoints of each original source interval;
/// None retains its authored affine h law or latent phase path. No source
/// selection, new table approximation, allocation or finite differencing.
pub(crate) fn integrate_caloric_table<E, const N: usize>(
    nodal: [f64; N],
    knots: &[(f64, f64)],
    capacities: &[Option<[f64; 2]>],
    with_derivatives: bool,
    mut poll: impl FnMut(usize) -> Result<(), E>,
) -> Result<ScalarMoments<N>, MomentError<E>> {
    if !(N == 3 || N == 4)
        || knots.len() < 2
        || capacities.len() + 1 != knots.len()
        || nodal.iter().any(|h| !h.is_finite())
    {
        return Err(MomentError::Numerical(
            "invalid P1 caloric conductivity input",
        ));
    }
    let mut sorted = nodal;
    sorted.sort_by(f64::total_cmp);
    let (low, high) = (sorted[0], sorted[N - 1]);
    if low < knots[0].0 || high > knots[knots.len() - 1].0 {
        return Err(MomentError::Numerical(
            "unsupported P1 caloric conductivity input",
        ));
    }
    let first = segment(knots, low);
    let last = segment(knots, high);
    let nonlinear = (first..=last).any(|index| {
        capacities[index].is_some_and(|[a, b]| a != b) && knots[index].1 != knots[index + 1].1
    });
    if !nonlinear {
        return integrate_table(nodal, knots, with_derivatives, poll);
    }
    for cp in capacities[first..=last].iter().flatten() {
        if cp.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(MomentError::Numerical(
                "positive finite original Cp required",
            ));
        }
    }
    if low == high {
        let (a, ka) = knots[first];
        let (b, kb) = knots[first + 1];
        let fraction = (low - a) / (b - a);
        let cp = capacities[first].unwrap_or([1.0, 1.0]);
        let (fraction, local_cp) = temperature_fraction(cp, fraction);
        let scale = cp[0].max(cp[1]);
        let slope = if with_derivatives {
            ((kb - ka) / (b - a)) * ((0.5 * (cp[0] / scale) + 0.5 * (cp[1] / scale)) / local_cp)
                / N as f64
        } else {
            0.0
        };
        return checked(interpolate(ka, kb, fraction), [slope; N]);
    }
    let range = high - low;
    if !range.is_finite() || range <= 0.0 {
        return Err(MomentError::Numerical("unrepresentable P1 caloric range"));
    }
    let mut duplicated = [[0.0; 5]; N];
    if with_derivatives {
        for (i, row) in duplicated.iter_mut().enumerate() {
            row[..N].copy_from_slice(&sorted);
            row[N] = nodal[i];
            row[..=N].sort_by(f64::total_cmp);
        }
    }
    const RULE: [(f64, f64); 4] = [
        (0.069_431_844_202_973_71, 0.173_927_422_568_726_93),
        (0.330_009_478_207_571_87, 0.326_072_577_431_273_07),
        (0.669_990_521_792_428_1, 0.326_072_577_431_273_07),
        (0.930_568_155_797_026_2, 0.173_927_422_568_726_93),
    ];
    let mut sums = [Sum::default(); 5];
    let (mut lo, mut piece, mut node, mut visited) = (low, first, 0, 0);
    while lo < high {
        if visited % 128 == 0 {
            poll(visited).map_err(MomentError::Interrupted)?;
        }
        while node < N && sorted[node] <= lo {
            node += 1;
        }
        let (h0, k0) = knots[piece];
        let (h1, k1) = knots[piece + 1];
        let hi = high.min(h1).min(sorted[node]);
        if hi <= lo || !((h1 - h0).is_finite() && h1 > h0) {
            return Err(MomentError::Numerical("invalid caloric source interval"));
        }
        let cp = capacities[piece].unwrap_or([1.0, 1.0]);
        let (ta, ca) = temperature_fraction(cp, (lo - h0) / (h1 - h0));
        let (tb, cb) = temperature_fraction(cp, (hi - h0) / (h1 - h0));
        let mean_cp = 0.5 * ca + 0.5 * cb;
        if !mean_cp.is_finite() || mean_cp <= 0.0 {
            return Err(MomentError::Numerical(
                "unrepresentable caloric interval capacity",
            ));
        }
        let ka = interpolate(k0, k1, ta);
        let kb = interpolate(k0, k1, tb);
        let interval_weight = (hi - lo) / range;
        // dk/dh dh cancels the capacity. Retain the conductivity change
        // directly from the source fraction, not two rounded absolute values.
        let derivative_weight = if with_derivatives {
            let scale = cp[0].max(cp[1]);
            let original_mean = 0.5 * (cp[0] / scale) + 0.5 * (cp[1] / scale);
            ((k1 - k0) / (h1 - h0)) * (original_mean / mean_cp) * interval_weight
        } else {
            0.0
        };
        for (alpha, gauss_weight) in RULE {
            let heat_fraction = alpha * ((1.0 - 0.5 * alpha) * ca + 0.5 * alpha * cb) / mean_cp;
            let capacity_ratio = ((1.0 - alpha) * ca + alpha * cb) / mean_cp;
            let k = interpolate(ka, kb, alpha);
            let weight = (N - 1) as f64
                * gauss_weight
                * interval_weight
                * capacity_ratio
                * bspline(&sorted, lo, hi, heat_fraction);
            sums[0].add(k * weight);
            if with_derivatives {
                for (i, row) in duplicated.iter().enumerate() {
                    sums[i + 1].add(
                        derivative_weight
                            * gauss_weight
                            * bspline(&row[..=N], lo, hi, heat_fraction),
                    );
                }
            }
        }
        lo = hi;
        if lo == h1 && piece + 2 < knots.len() {
            piece += 1;
        }
        visited += 1;
    }
    poll(visited).map_err(MomentError::Interrupted)?;
    checked(
        sums[0].total(),
        std::array::from_fn(|i| sums[i + 1].total()),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonlinear_source_moments_match_independent_temperature_polynomial_integral() {
        // h=T^2/2, Cp=T, k=T on T in [1,3]. With h=[.5,.5,.5,4.5],
        // the h density is 3*(4.5-h)^2/4^3. Transform dh=T dT.
        let table = [(0.5, 1.0), (4.5, 3.0)];
        let cp = [Some([1.0, 3.0])];
        let nodes = [0.5, 0.5, 0.5, 4.5];
        let result =
            integrate_caloric_table(nodes, &table, &cp, true, |_| Ok::<_, ()>(())).unwrap();
        let primitive =
            |t: f64| 3.0 / 64.0 * (6.75 * t.powi(3) - 0.9 * t.powi(5) + t.powi(7) / 28.0);
        let expected = primitive(3.0) - primitive(1.0);
        assert!(
            (result.average - expected).abs() < 3e-14,
            "{} != {expected}",
            result.average
        );
        // Independent derivative: integrate lambda_i/T(h) in h. For the
        // hot node lambda=(h-.5)/4, and each cold-node conditional mean is
        // (1-lambda)/3. Transformation cancels 1/T exactly.
        let hot_primitive = |t: f64| {
            3.0 / 256.0 * (-10.125 * t + 4.125 * t.powi(3) - 0.475 * t.powi(5) + t.powi(7) / 56.0)
        };
        let hot = hot_primitive(3.0) - hot_primitive(1.0);
        let all_primitive = |t: f64| 3.0 / 64.0 * (20.25 * t - 1.5 * t.powi(3) + t.powi(5) / 20.0);
        let all = all_primitive(3.0) - all_primitive(1.0);
        assert!((result.derivatives[3] - hot).abs() < 3e-14);
        for d in &result.derivatives[..3] {
            assert!((*d - (all - hot) / 3.0).abs() < 3e-14);
        }
    }

    #[test]
    fn nonlinear_moments_cross_source_knots_and_retain_consistent_nodal_tangents() {
        let table = [(0.5, 1.0), (2.0, 2.0), (4.5, 3.0)];
        let cp = [Some([1.0, 2.0]), Some([2.0, 3.0])];
        let nodes = [0.8, 1.7, 3.1, 4.1];
        let run = |h| integrate_caloric_table(h, &table, &cp, true, |_| Ok::<_, ()>(())).unwrap();
        let result = run(nodes);
        for i in 0..4 {
            let mut lo = nodes;
            let mut hi = nodes;
            lo[i] -= 1e-4;
            hi[i] += 1e-4;
            let finite_difference = (run(hi).average - run(lo).average) / 2e-4;
            assert!((result.derivatives[i] - finite_difference).abs() < 2e-10);
        }
        let uniform = run([2.0; 4]);
        assert_eq!(uniform.average, 2.0);
        for d in uniform.derivatives {
            assert!((d - 0.125).abs() < 2e-16);
        }
        let affine_cp = [Some([1000.0, 1000.0]), Some([1000.0, 1000.0])];
        let retained =
            integrate_caloric_table(nodes, &table, &affine_cp, true, |_| Ok::<_, ()>(())).unwrap();
        let native = integrate_table(nodes, &table, true, |_| Ok::<_, ()>(())).unwrap();
        assert_eq!(retained.average, native.average);
        assert_eq!(retained.derivatives, native.derivatives);
    }

    #[test]
    fn adjacent_enthalpies_retain_positive_nodal_conductivity_derivatives() {
        let table = [(0.0, 1.0), (1_500_000.0, 3.0)];
        let cp = [Some([1000.0, 4000.0])];
        let h = 10_000.0_f64;
        let result = integrate_caloric_table([h, h, h, h.next_up()], &table, &cp, true, |_| {
            Ok::<_, ()>(())
        })
        .unwrap();
        let fraction = h / 1_500_000.0;
        let local_cp =
            ((1.0 - fraction) * 1000.0_f64.powi(2) + fraction * 4000.0_f64.powi(2)).sqrt();
        let expected = (2.0 / 1_500_000.0) * (2500.0 / local_cp) / 4.0;
        for derivative in result.derivatives {
            assert!(derivative > 0.0);
            assert!((derivative - expected).abs() < 8.0 * f64::EPSILON * expected);
        }
    }
}
