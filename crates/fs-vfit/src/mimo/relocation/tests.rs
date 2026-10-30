use super::super::{C64, MimoModel, PoleTerm, det, mimo_vector_fit};
use super::*;

#[test]
fn seventeen_port_fit_retains_a_mode_visible_only_at_the_last_opening() {
    let n = 17;
    let poles = vec![
        PoleTerm::Pair {
            pole: C64::new(-70.0, 500.0),
            residue: C64::ZERO,
        },
        PoleTerm::Pair {
            pole: C64::new(-100.0, 2300.0),
            residue: C64::ZERO,
        },
    ];
    let common: Vec<_> = (0..n * n)
        .map(|i| {
            C64::new(100.0, 20.0)
                .scale((1.0 + (i / n) as f64 / 20.0) * (1.0 + (i % n) as f64 / 20.0))
        })
        .collect();
    let mut last = vec![C64::ZERO; n * n];
    last[n * n - 1] = C64::new(300.0, -30.0);
    let truth = MimoModel {
        ports: n,
        poles,
        residues: vec![common, last],
        d: (0..n * n)
            .map(|i| if i / n == i % n { 0.8 } else { 0.0 })
            .collect(),
        e: vec![0.0; n * n],
        reciprocal: true,
    };
    let omega: Vec<_> = (0..96)
        .map(|i| 80.0 * det::exp(f64::from(i) * det::ln(80.0) / 95.0))
        .collect();
    let data: Vec<_> = omega.iter().map(|&w| truth.eval_iw(w)).collect();
    let options = FitOptions {
        fit_d: true,
        fit_e: false,
        ..FitOptions::new(4)
    };
    let (fit, report) = mimo_vector_fit(&omega, &data, n, &options, true).unwrap();
    assert!(fit.is_stable() && fit.reciprocal);
    assert!(report.max_relative_error < 1e-8, "{report:?}");
    for w in [217.0, 573.0, 2319.0, 3973.0] {
        let expected = truth.eval_iw(w);
        let actual = fit.eval_iw(w);
        for (&a, &b) in actual.iter().zip(&expected) {
            assert!((a - b).abs() < 1e-7 * (1.0 + b.abs()));
        }
        assert_eq!(actual[16], actual[16 * n]);
    }
}

#[test]
fn shared_scratch_bound_scales_with_response_count_and_refuses_overflow() {
    let options = FitOptions {
        fit_e: false,
        ..FitOptions::new(12)
    };
    let bytes = shared_pole_relocation_scratch_bytes(65, 17 * 18 / 2, &options).unwrap();
    assert!(bytes < 16 * 1024 * 1024, "{bytes}");
    let dense_rows = 2 * 65 * 153 + 1;
    let dense_columns = 153 * 13 + 13;
    assert!(8 * dense_rows * dense_columns > 30 * bytes);
    assert!(shared_pole_relocation_scratch_bytes(usize::MAX, usize::MAX, &options).is_none());
}

#[test]
fn shared_scratch_bound_covers_live_fixed_sigma_fallback_vectors() {
    for (samples, responses, order) in [(65usize, 561usize, 3usize), (513, 561, 7)] {
        let options = FitOptions {
            fit_d: true,
            fit_e: false,
            ..FitOptions::new(order)
        };
        // Every reciprocal response of a 33-port array loses its local
        // numerator rows before the common denominator is solved.
        let projected_rows = (2 * samples - order - 1) * responses + 1;
        let columns = order + 1;
        let saved_rows = projected_rows - 1;
        // The original A/RHS and saved sample rows remain live while the
        // fallback owns its A/RHS and QR's copies of both. This independent
        // payload inventory excludes small scale/tau vectors, so it is a
        // lower bound even before allocator bookkeeping.
        let matrices = projected_rows * columns
            + saved_rows * columns
            + 2 * saved_rows * (columns - 1);
        let right_hand_sides = projected_rows + 2 * saved_rows;
        let minimum_bytes = 8 * (matrices + right_hand_sides);
        let bound = shared_pole_relocation_scratch_bytes(samples, responses, &options).unwrap();
        assert!(
            bound >= minimum_bytes,
            "{samples} samples at order {order}: {bound} < {minimum_bytes}"
        );
    }
}
