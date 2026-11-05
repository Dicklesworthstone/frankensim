//! Retarded pressure from a uniform circular aperture in an infinite rigid baffle.
//!
//! Rayleigh's surface integral reduces to a one-dimensional spatial impulse
//! response: `p = rho/A * h * dQ/dt`, where `h = c * arc_fraction`. The fraction
//! is the portion of a circle, centered below the receiver, inside the aperture.
//! This retains the near field and aperture interference without a point-source
//! approximation or a surface quadrature in the audio loop. See Stepanishen,
//! JASA 49 (1971), 1629–1638, equations (15)–(17):
//! <https://doi.org/10.1121/1.1912541>.

use std::f64::consts::PI;

const MAX_SPAN: usize = 4096;
const MAX_DELAY: usize = 1_048_576;
const QUADRATURE_TOLERANCE: f64 = 1e-9;

/// `(first sample, consecutive weights)` multiplying integer-time volume-flow
/// acceleration. The cubic propagation basis is the same as the head receiver's;
/// below one sample of flight it uses the same causal linear basis. Coefficients
/// are aggregated by delay, so the runtime cost depends on aperture transit time,
/// not on a spatial mesh or source-receiver distance.
pub(crate) fn circular_taps(
    area: f64,
    offset: [f64; 3],
    density: f64,
    sound_speed: f64,
    rate: f64,
) -> Result<(usize, Vec<f64>), String> {
    if ![area, density, sound_speed, rate]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
        || !offset.iter().all(|x| x.is_finite())
        || offset[2] == 0.0
    {
        return Err("circular vent receiver needs finite positive area, medium, clock and distance from the baffle".into());
    }
    let a = (area / PI).sqrt();
    let b = offset[0].hypot(offset[1]);
    let z = offset[2].abs();
    let s_min = (b - a).max(0.0);
    let s_max = b + a;
    let r_min = z.hypot(s_min);
    let r_max = z.hypot(s_max);
    let delay = r_min / sound_speed * rate;
    // Rationalize the path-length difference, particularly on-axis in the far field.
    let width = (s_max - s_min) * (s_max + s_min) / (r_max + r_min) / sound_speed * rate;
    let gain = density * sound_speed / (area * rate);
    if ![delay, width, gain]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
        || delay + width > MAX_DELAY as f64
        || width.ceil() > (MAX_SPAN - 4) as f64
        || delay + width == delay
    {
        return Err(
            "circular vent receiver exceeds representable geometry or its delay/tap budget".into(),
        );
    }
    let first_cell = delay.floor() as usize;
    let last_cell = (delay + width).floor() as usize;
    let first = first_cell.saturating_sub(1);
    let mut taps = vec![0.0; last_cell.saturating_sub(1) + 4 - first];
    // Below a point inside the disk, the full-circle interval ends at its near
    // edge. Split there so each integration segment has smooth interior data.
    let full_until = if b < a {
        (a - b).powi(2) / (z.hypot(a - b) + z) / sound_speed * rate
    } else {
        0.0
    };
    for cell in first_cell..=last_cell {
        let lo = (cell as f64 - delay).max(0.0);
        let hi = ((cell + 1) as f64 - delay).min(width);
        if hi <= lo {
            continue;
        }
        let split = full_until.clamp(lo, hi);
        let start = cell.saturating_sub(1);
        for (left, right) in [(lo, split), (split, hi)] {
            if right <= left {
                continue;
            }
            let integrand = |u: f64| {
                // This endpoint map removes the square-root edge singularities.
                let v = (0.5 * PI * u).sin();
                let t = left + (right - left) * v * v;
                let jacobian = (right - left) * 0.5 * PI * (PI * u).sin();
                let dr = t / rate * sound_speed;
                let square_increment = dr * (2.0 * r_min + dr);
                let s = (s_min * s_min + square_increment).sqrt();
                let fraction = if b == 0.0 || (b < a && s <= a - b) {
                    1.0
                } else {
                    // acos(1-epsilon) loses the far-field aperture angle. The
                    // half-angle identity retains it; outside the disk, use the
                    // radial increment from its near edge without subtracting b.
                    let inside = if b >= a {
                        let ds = square_increment / (s + s_min);
                        ds * (2.0 * a - ds)
                    } else {
                        (a - b + s) * (a + b - s)
                    };
                    2.0 / PI * (inside / (4.0 * b * s)).clamp(0.0, 1.0).sqrt().asin()
                };
                let (index, weights) = super::microphone_delay_taps(delay + t);
                let mut value = [0.0; 4];
                for (i, weight) in weights.into_iter().enumerate() {
                    if weight != 0.0 {
                        value[index + i - start] += gain * fraction * jacobian * weight;
                    }
                }
                value
            };
            let integral = adaptive(&integrand, 0.0, 1.0, 12, None)?;
            for (i, value) in integral.into_iter().enumerate() {
                taps[start + i - first] += value;
            }
        }
    }
    if !taps.iter().all(|x| x.is_finite()) || taps.iter().sum::<f64>() <= 0.0 {
        return Err("circular vent pressure coefficients are not representable".into());
    }
    Ok((first, taps))
}

const GAUSS_4: [(f64, f64); 2] = [
    (0.3399810435848563, 0.6521451548625461),
    (0.8611363115940526, 0.3478548451374538),
];
const GAUSS_8: [(f64, f64); 4] = [
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];

fn integrate(f: &impl Fn(f64) -> [f64; 4], lo: f64, hi: f64, rule: &[(f64, f64)]) -> [f64; 4] {
    let mid = 0.5 * (lo + hi);
    let half = 0.5 * (hi - lo);
    let mut sum = [0.0; 4];
    for &(x, weight) in rule {
        let left = f(mid - half * x);
        let right = f(mid + half * x);
        for i in 0..4 {
            sum[i] += half * weight * (left[i] + right[i]);
        }
    }
    sum
}

fn adaptive(
    f: &impl Fn(f64) -> [f64; 4],
    lo: f64,
    hi: f64,
    depth: u32,
    budget: Option<f64>,
) -> Result<[f64; 4], String> {
    let coarse = integrate(f, lo, hi, &GAUSS_4);
    let fine = integrate(f, lo, hi, &GAUSS_8);
    let error: f64 = fine.iter().zip(coarse).map(|(a, b)| (a - b).abs()).sum();
    let scale: f64 = fine.iter().map(|x| x.abs()).sum();
    // Split the parent coefficient-error allowance. Requiring relative
    // accuracy again on a vanishing edge sliver instead chases roundoff.
    let allowance = budget.unwrap_or(QUADRATURE_TOLERANCE * scale);
    if error <= allowance {
        return Ok(fine);
    }
    if depth == 0 {
        return Err(
            "circular vent receiver quadrature did not converge within its setup budget".into(),
        );
    }
    let mid = 0.5 * (lo + hi);
    let a = adaptive(f, lo, mid, depth - 1, Some(0.5 * allowance))?;
    let b = adaptive(f, mid, hi, depth - 1, Some(0.5 * allowance))?;
    Ok(std::array::from_fn(|i| a[i] + b[i]))
}

#[cfg(test)]
mod tests {
    use super::*;

    const RHO: f64 = 1.225;
    const C: f64 = 343.0;

    fn response(first: usize, taps: &[f64], frequency: f64, rate: f64) -> [f64; 2] {
        let mut value = [0.0; 2];
        for (i, tap) in taps.iter().enumerate() {
            let phase = 2.0 * PI * frequency * (first + i) as f64 / rate;
            value[0] += tap * phase.cos();
            value[1] -= tap * phase.sin();
        }
        value
    }

    fn error(a: [f64; 2], b: [f64; 2]) -> f64 {
        (a[0] - b[0]).hypot(a[1] - b[1])
    }

    #[test]
    fn g1_axis_pressure_matches_exact_disk_and_refines_in_time() {
        let (a, z) = (0.08_f64, 0.18_f64);
        let area = PI * a * a;
        let dr = a * a / (z.hypot(a) + z);
        let gain = RHO * dr / area;
        let frequency = 4000.0;
        let half_phase = PI * frequency * dr / C;
        let phase = 2.0 * PI * frequency * (z + 0.5 * dr) / C;
        let amplitude = gain * half_phase.sin() / half_phase;
        let exact = [amplitude * phase.cos(), -amplitude * phase.sin()];
        let mut previous = f64::INFINITY;
        for rate in [24_000.0, 48_000.0, 96_000.0] {
            let (first, taps) = circular_taps(area, [0.0, 0.0, z], RHO, C, rate).unwrap();
            assert!((taps.iter().sum::<f64>() / gain - 1.0).abs() < 3e-12);
            let moment: f64 = taps
                .iter()
                .enumerate()
                .map(|(i, w)| (first + i) as f64 * w)
                .sum();
            assert!((moment / (RHO * rate / (2.0 * PI * C)) - 1.0).abs() < 3e-12);
            let e = error(response(first, &taps, frequency, rate), exact) / gain;
            assert!(
                e < previous / 8.0,
                "clock refinement: rate={rate}, error={e}, previous={previous}"
            );
            previous = e;
        }
        assert!(previous < 2e-4);
        // The near-field finite-disk response differs measurably from a point at its center.
        let phase = 2.0 * PI * frequency * z / C;
        let compact = RHO / (2.0 * PI * z);
        assert!(error(exact, [compact * phase.cos(), -compact * phase.sin()]) / gain > 0.2);
    }

    /// Independent two-dimensional equal-area Rayleigh quadrature. It does not
    /// use the arc formula, time bins, or interpolation of the implementation.
    fn surface_response(a: f64, offset: [f64; 3], frequency: f64, rings: usize) -> [f64; 2] {
        let angles = 4 * rings;
        let mut value = [0.0; 2];
        let directions: Vec<_> = (0..angles)
            .map(|j| {
                let theta = 2.0 * PI * (j as f64 + 0.5) / angles as f64;
                (theta.cos(), theta.sin())
            })
            .collect();
        for i in 0..rings {
            let radius = a * ((i as f64 + 0.5) / rings as f64).sqrt();
            for &(x, y) in &directions {
                let distance = (offset[0] - radius * x)
                    .hypot(offset[1] - radius * y)
                    .hypot(offset[2]);
                let phase = 2.0 * PI * frequency * distance / C;
                value[0] += phase.cos() / distance;
                value[1] -= phase.sin() / distance;
            }
        }
        value.map(|x| x * RHO / (2.0 * PI * (rings * angles) as f64))
    }

    #[test]
    fn g2_off_axis_pressure_matches_surface_integral_inside_at_and_outside_rim() {
        let a = 0.05;
        let rate = 192_000.0;
        for offset in [[0.02_f64, 0.0, 0.08], [a, 0.0, 0.08], [0.09, 0.12, 0.08]] {
            let (first, taps) = circular_taps(PI * a * a, offset, RHO, C, rate).unwrap();
            let moment: f64 = taps
                .iter()
                .enumerate()
                .map(|(i, w)| (first + i) as f64 * w)
                .sum();
            assert!((moment / (RHO * rate / (2.0 * PI * C)) - 1.0).abs() < 2e-10);
            for frequency in [0.0, 6000.0] {
                let reference = surface_response(a, offset, frequency, 256);
                let finer = surface_response(a, offset, frequency, 512);
                let scale = RHO / (2.0 * PI * offset[0].hypot(offset[1]).hypot(offset[2]));
                assert!(
                    error(reference, finer) / scale < 2e-5,
                    "independent surface oracle must resolve"
                );
                assert!(
                    error(response(first, &taps, frequency, rate), finer) / scale < 3e-5,
                    "receiver={offset:?}, frequency={frequency}"
                );
            }
        }
    }

    #[test]
    fn g2_far_field_retains_bessel_directivity_and_compact_limit() {
        let (a, distance, rate) = (0.05_f64, 100.0_f64, 96_000.0);
        let theta = 0.8_f64;
        let offset = [distance * theta.sin(), 0.0, distance * theta.cos()];
        let (first, taps) = circular_taps(PI * a * a, offset, RHO, C, rate).unwrap();
        assert!(taps.len() < 2.0_f64.mul_add(a / C * rate, 6.0).ceil() as usize);
        let gain = RHO / (2.0 * PI * distance);
        for frequency in [10.0, 6000.0] {
            let ka = 2.0 * PI * frequency * a / C;
            let argument = ka * theta.sin();
            let directivity = 2.0 * fs_math::special::bessel_j1(argument) / argument;
            let phase = 2.0 * PI * frequency * distance / C;
            let expected = [
                gain * directivity * phase.cos(),
                -gain * directivity * phase.sin(),
            ];
            assert!(error(response(first, &taps, frequency, rate), expected) / gain < 8e-4);
        }
        assert!((taps.iter().sum::<f64>() / gain - 1.0).abs() < 2e-6);
    }

    #[test]
    fn g0_subsample_geometry_and_invalid_budgets_are_explicit() {
        let (first, taps) = circular_taps(PI * 1e-8, [0.0, 0.0, 0.0001], RHO, C, 48_000.0).unwrap();
        assert_eq!(first, 0);
        assert!(taps.iter().all(|x| x.is_finite()));
        for (area, offset, rate) in [
            (0.0, [0.0, 0.0, 1.0], 48_000.0),
            (0.01, [0.0; 3], 48_000.0),
            (0.01, [0.0, 0.0, 1e8], 48_000.0),
            (0.01, [0.1, 0.0, 1.0], 1e12),
            (0.01, [f64::NAN, 0.0, 1.0], 48_000.0),
        ] {
            assert!(circular_taps(area, offset, RHO, C, rate).is_err());
        }
    }
}
