use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())), "{a} != {b}");
}

#[test]
fn tied_two_designs_have_the_closed_form_value() {
    let result = expected_minimum_improvement(&[4.0, 4.0], &[-2.0, 1.0]).unwrap();
    close(result.expected_improvement(), 3.0 / std::f64::consts::TAU.sqrt(), 2e-14);
    assert_eq!(result.envelope_segments(), 2);
    assert!(result.has_information());
}

#[test]
fn common_mean_motion_has_no_decision_value() {
    for means in [&[1.0, 2.0, 4.0][..], &[1.0, 1.0, 1.0][..]] {
        let result = expected_minimum_improvement(means, &[8.0, 8.0, 8.0]).unwrap();
        assert_eq!(result.expected_improvement(), 0.0);
        assert_eq!(result.log_expected_improvement(), f64::NEG_INFINITY);
        assert!(!result.has_information());
    }
}

#[test]
fn a_worse_mean_third_design_can_carry_the_decision_value() {
    let top_two = expected_minimum_improvement(&[0.0, 0.1], &[0.0, 0.0]).unwrap();
    let full = expected_minimum_improvement(&[0.0, 0.1, 2.0], &[0.0, 0.0, 5.0]).unwrap();
    assert!(!top_two.has_information());
    assert!(full.expected_improvement() > 1.0);
    assert_eq!(full.envelope_segments(), 2);
}

#[test]
fn permutation_common_shift_and_unit_changes_preserve_the_decision_algebra() {
    let a = [0.0, 0.25, 2.0, 1.0];
    let b = [-0.25, 0.0, 3.0, -2.0];
    let base = expected_minimum_improvement(&a, &b).unwrap();
    for order in [[0, 1, 2, 3], [3, 1, 0, 2], [1, 3, 2, 0]] {
        let means = order.map(|i| 8.0 * a[i] + 1024.0);
        let slopes = order.map(|i| 8.0 * b[i] + 32.0);
        let result = expected_minimum_improvement(&means, &slopes).unwrap();
        close(result.expected_improvement(), 8.0 * base.expected_improvement(), 3e-13);
        close(result.log_expected_improvement(), base.log_expected_improvement() + det::ln(8.0), 3e-13);
    }
}

#[test]
fn dominated_parallel_lines_and_coincident_crossings_are_removed() {
    let result = expected_minimum_improvement(&[0.0, 1.0, 0.0, 0.0], &[0.0, 0.0, 1.0, -1.0]).unwrap();
    assert_eq!(result.envelope_segments(), 2);
    close(result.expected_improvement(), (2.0 / std::f64::consts::PI).sqrt(), 2e-14);
}

#[test]
fn far_tail_rankings_survive_linear_underflow_and_large_units() {
    let a = expected_minimum_improvement(&[0.0, 40.0], &[0.0, 1.0]).unwrap();
    let b = expected_minimum_improvement(&[0.0, 41.0], &[0.0, 1.0]).unwrap();
    assert_eq!(a.expected_improvement(), 0.0);
    assert_eq!(b.expected_improvement(), 0.0);
    assert!(a.has_information() && b.has_information());
    assert!(a.log_expected_improvement() > b.log_expected_improvement());
    let scaled = expected_minimum_improvement(&[0.0, 4e301], &[0.0, 1e300]).unwrap();
    assert!(scaled.expected_improvement() > 0.0);
    close(scaled.log_expected_improvement(), a.log_expected_improvement() + det::ln(1e300), 2e-13);
}

#[test]
fn envelope_matches_direct_outcome_integration() {
    // Independent composite Simpson integral of baseline(z) - min_i line_i(z).
    // It does not build the hull or use the hinge formula under test.
    for seed in 0..8 {
        let means: Vec<_> = (0..7).map(|i| ((i * 7 + seed * 3) % 17) as f64 / 5.0).collect();
        let slopes: Vec<_> = (0..7).map(|i| ((i * 11 + seed * 5) % 19) as f64 / 3.0 - 3.0).collect();
        let result = expected_minimum_improvement(&means, &slopes).unwrap();
        let best = (0..means.len()).min_by(|&i, &j| means[i].total_cmp(&means[j])).unwrap();
        let panels = 65_536;
        let step = 24.0 / panels as f64;
        let mut sum = 0.0;
        for j in 0..=panels {
            let z = -12.0 + step * j as f64;
            let minimum = means.iter().zip(&slopes).map(|(a, b)| b.mul_add(z, *a)).fold(f64::INFINITY, f64::min);
            let value = (slopes[best].mul_add(z, means[best]) - minimum) * det::exp(-0.5 * z * z) / std::f64::consts::TAU.sqrt();
            sum += value * if j == 0 || j == panels { 1.0 } else if j % 2 == 1 { 4.0 } else { 2.0 };
        }
        close(result.expected_improvement(), sum * step / 3.0, 2e-7);
    }
}

#[test]
fn cancellation_and_bad_inputs_never_publish_a_partial_value() {
    let mut calls = 0;
    let cancelled = expected_minimum_improvement_controlled(&[0.0, 10.0], &[0.0, 1.0], || {
        calls += 1; calls < 12
    });
    assert_eq!(cancelled, Err(KnowledgeGradientError::Interrupted));
    for (a, b) in [(vec![], vec![]), (vec![1.0], vec![]), (vec![0.0; 257], vec![1.0; 257])] {
        assert_eq!(expected_minimum_improvement(&a, &b), Err(KnowledgeGradientError::InvalidShape));
    }
    assert_eq!(expected_minimum_improvement(&[0.0, f64::NAN], &[0.0, 1.0]), Err(KnowledgeGradientError::NonFiniteInput));
    assert_eq!(expected_minimum_improvement(&[-f64::MAX, f64::MAX], &[1.0, 0.0]), Err(KnowledgeGradientError::NumericalRange));
}
