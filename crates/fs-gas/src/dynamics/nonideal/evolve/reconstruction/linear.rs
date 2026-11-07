//! Source-independent monotone reconstruction and bounded constitutive admission.
//! Geometry uses local cell widths, avoiding cancellation at distant coordinates.

/// A pair of source-admitted endpoint traces and an explicit limiting flag.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Pair<T> {
    pub(crate) traces: [T; 2],
    pub(crate) limited: bool,
}

/// Monotonized-central increments over the CURRENT cell width.
/// Affine data are reproduced on nonuniform meshes. A local extremum gets zero
/// increment. Adjacent center distances are half sums of physical cell widths.
pub(crate) fn increments(
    left: [f64; 3],
    center: [f64; 3],
    right: [f64; 3],
    widths: [f64; 3],
) -> [f64; 3] {
    let [wl, w, wr] = widths;
    let hl = 0.5 * wl + 0.5 * w;
    let hr = 0.5 * w + 0.5 * wr;
    let central_scale = (0.5 * w) / (0.5 * hl + 0.5 * hr);
    core::array::from_fn(|k| {
        minmod(
            2.0 * ((center[k] - left[k]) * (w / hl)),
            (right[k] - left[k]) * central_scale,
            2.0 * ((right[k] - center[k]) * (w / hr)),
        )
    })
}

fn minmod(a: f64, b: f64, c: f64) -> f64 {
    if a > 0.0 && b > 0.0 && c > 0.0 {
        a.min(b).min(c)
    } else if a < 0.0 && b < 0.0 && c < 0.0 {
        a.max(b).max(c)
    } else {
        0.0
    }
}

/// Contract a complete primitive increment toward its original admitted cell.
/// BOTH traces must pass the same original constitutive source before either
/// is used. The one scale retains their coupled direction; individual
/// density/temperature clipping would alter that direction.
/// Twelve contractions bound source work; exhausted admission uses the exact
/// original cell payload. Terminal errors (in particular cancellation) escape.
pub(crate) fn admitted_pair<T: Copy, E>(
    original: T,
    primitive: [f64; 3],
    increment: [f64; 3],
    mut resolve: impl FnMut([f64; 3]) -> Result<T, E>,
    terminal: impl Fn(&E) -> bool,
) -> Result<Pair<T>, E> {
    if increment == [0.0; 3] {
        return Ok(Pair {
            traces: [original; 2],
            limited: false,
        });
    }
    if increment.iter().any(|x| !x.is_finite()) {
        return Ok(Pair {
            traces: [original; 2],
            limited: true,
        });
    }
    let mut scale = 1.0;
    for reduction in 0..=12 {
        let values =
            |side: f64| core::array::from_fn(|k| primitive[k] + side * (scale * increment[k]));
        let left_values = values(-0.5);
        let left_result = if left_values == primitive {
            Ok(original)
        } else {
            resolve(left_values)
        };
        let left = match left_result {
            Ok(value) => value,
            Err(error) => {
                if terminal(&error) {
                    return Err(error);
                }
                scale *= 0.5;
                continue;
            }
        };
        let right_values = values(0.5);
        let right_result = if right_values == primitive {
            Ok(original)
        } else {
            resolve(right_values)
        };
        match right_result {
            Ok(right) => {
                return Ok(Pair {
                    traces: [left, right],
                    limited: reduction != 0,
                });
            }
            Err(error) => {
                if terminal(&error) {
                    return Err(error);
                }
                scale *= 0.5;
            }
        }
    }
    Ok(Pair {
        traces: [original; 2],
        limited: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) {
        assert!((a - b).abs() < 8e-14, "{a} != {b}");
    }

    #[test]
    fn affine_fields_reconstruct_on_nonuniform_cells_without_absolute_coordinates() {
        let widths = [0.2, 0.7, 0.3];
        let center = [13.0, -20.0, 450.0];
        let gradients = [2.0, -3.0, 5.0];
        let left = core::array::from_fn(|k| center[k] - 0.45 * gradients[k]);
        let right = core::array::from_fn(|k| center[k] + 0.5 * gradients[k]);
        let slope = increments(left, center, right, widths);
        for k in 0..3 {
            near(slope[k], 0.7 * gradients[k]);
        }
        let pair = admitted_pair(center, center, slope, |x| Ok::<_, ()>(x), |_| false).unwrap();
        assert!(!pair.limited);
        for k in 0..3 {
            near(pair.traces[0][k], center[k] - 0.35 * gradients[k]);
            near(pair.traces[1][k], center[k] + 0.35 * gradients[k]);
        }
    }

    #[test]
    fn monotone_limiter_adds_no_overshoot_and_flattens_extrema() {
        let s = increments([1.0, 4.0, 1.0], [2.0, 2.0, 3.0], [8.0, 1.0, 2.0], [1.0; 3]);
        assert_eq!(s, [2.0, -1.5, 0.0]);
        for k in 0..2 {
            let a: f64 = [1.0, 4.0][k];
            let b: f64 = [2.0, 2.0][k];
            let c: f64 = [8.0, 1.0][k];
            assert!((a.min(b)..=a.max(b)).contains(&(b - 0.5 * s[k])));
            assert!((b.min(c)..=b.max(c)).contains(&(b + 0.5 * s[k])));
        }
    }

    #[test]
    fn original_source_admission_contracts_both_traces_with_one_scale() {
        let center = [10.0, 20.0, 400.0];
        let increment = [8.0, 80.0, 40.0];
        let source = |x: [f64; 3]| {
            if (9.0..=11.0).contains(&x[0]) && (395.0..=405.0).contains(&x[2]) {
                Ok(x)
            } else {
                Err(())
            }
        };
        let pair = admitted_pair(center, center, increment, source, |_| false).unwrap();
        assert!(pair.limited);
        assert_eq!(pair.traces, [[9.0, 10.0, 395.0], [11.0, 30.0, 405.0]]);
    }

    #[test]
    fn source_exhaustion_reuses_original_payload_instead_of_reencoding() {
        let original = (17_u64, [10.0, 20.0, 400.0]);
        let mut attempts = 0;
        let pair = admitted_pair(
            original,
            original.1,
            [1.0; 3],
            |_| {
                attempts += 1;
                Err::<(u64, [f64; 3]), _>(())
            },
            |_| false,
        )
        .unwrap();
        assert_eq!(attempts, 13);
        assert!(pair.limited);
        assert_eq!(pair.traces, [original; 2]);
    }

    #[test]
    fn cancellation_on_second_trace_is_never_a_limiter_fallback() {
        let mut calls = 0;
        let result = admitted_pair(
            [0.0; 3],
            [1.0; 3],
            [0.1; 3],
            |x| {
                calls += 1;
                if calls == 2 { Err("cancelled") } else { Ok(x) }
            },
            |error| *error == "cancelled",
        );
        assert_eq!(result.unwrap_err(), "cancelled");
        assert_eq!(calls, 2);
    }

    #[test]
    fn constant_state_uses_no_constitutive_calls() {
        let original = [10.0, -30.0, 420.0];
        let pair = admitted_pair(
            original,
            original,
            [0.0; 3],
            |_| -> Result<[f64; 3], ()> { panic!("constant source reevaluation") },
            |_| false,
        )
        .unwrap();
        assert_eq!(pair.traces, [original; 2]);
        assert!(!pair.limited);
    }
    #[test]
    fn unrepresented_trace_increments_keep_original_caloric_payload() {
        let original = (3_u64, [10.0, -30.0, 420.0]);
        let pair = admitted_pair(
            original,
            original.1,
            [1e-100; 3],
            |_| -> Result<(u64, [f64; 3]), ()> { panic!("rounded original was reevaluated") },
            |_| false,
        )
        .unwrap();
        assert_eq!(pair.traces, [original; 2]);
        assert!(!pair.limited);
    }

    #[test]
    fn smooth_periodic_transport_retains_second_order_and_extensive_mean() {
        fn error(n: usize) -> f64 {
            let dx = 1.0 / n as f64;
            let phase_half = core::f64::consts::PI * dx;
            let average = phase_half.sin() / phase_half;
            let mut q: Vec<_> = (0..n)
                .map(|i| {
                    1.0 + 0.2 * average * (core::f64::consts::TAU * (i as f64 + 0.5) * dx).sin()
                })
                .collect();
            let initial_mean = q.iter().sum::<f64>() * dx;
            let forward = |q: &[f64], dt: f64| -> Vec<f64> {
                let traces: Vec<_> = (0..n)
                    .map(|i| {
                        let left = q[(i + n - 1) % n];
                        let right = q[(i + 1) % n];
                        let delta = increments(
                            [1.0, left, 300.0],
                            [1.0, q[i], 300.0],
                            [1.0, right, 300.0],
                            [dx; 3],
                        );
                        q[i] + 0.5 * delta[1]
                    })
                    .collect();
                (0..n)
                    .map(|i| q[i] - (dt / dx) * (traces[i] - traces[(i + n - 1) % n]))
                    .collect()
            };
            let finish = 0.37;
            let mut time: f64 = 0.0;
            while time < finish {
                let dt = (0.3 * dx).min(finish - time);
                let first = forward(&q, dt);
                let second = forward(&first, dt);
                for (q, next) in q.iter_mut().zip(second) {
                    *q = 0.5 * *q + 0.5 * next;
                }
                time += dt;
            }
            assert!((q.iter().sum::<f64>() * dx - initial_mean).abs() < 1e-13);
            q.iter()
                .enumerate()
                .map(|(i, q)| {
                    let exact = 1.0
                        + 0.2
                            * average
                            * (core::f64::consts::TAU * ((i as f64 + 0.5) * dx - finish)).sin();
                    (q - exact).abs() * dx
                })
                .sum()
        }
        let errors = [32, 64, 128].map(error);
        assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.4 * errors[1], "{errors:?}");
    }
}
