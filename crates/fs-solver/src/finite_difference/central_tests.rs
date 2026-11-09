use super::*;
use std::cell::Cell;

fn controls() -> CentralDifferenceControls {
    CentralDifferenceControls { initial_step: 0.2, absolute_tolerance: 1e-10,
        relative_tolerance: 1e-8, maximum_refinements: 16, maximum_evaluations: 34 }
}

#[test]
fn g1_refined_central_action_matches_nonlinear_analytic_direction_and_refines() {
    let result = differentiate(&[2.0, -3.0], 2, controls(), |x| {
        Ok::<_, ()>(vec![(0.7+x[0]).sin() + (0.2+x[1]).powi(3),
            (0.7+x[0])*(0.2+x[1])])
    }, || Ok(())).unwrap();
    let exact = [2.0*0.7_f64.cos()-9.0*0.2_f64.powi(2), 0.4-2.1];
    for (&actual, expected) in result.action.iter().zip(exact) {
        assert!((actual-expected).abs() < 2e-8, "{actual} != {expected}");
    }
    assert!(result.refinements > 1);
    assert_eq!(result.evaluations, 2*(result.refinements+1));
}

#[test]
fn g3_direction_scaling_retains_small_and_large_physical_amplitudes() {
    for amplitude in [1e-300, 1e-12, 1.0, 1e150, 1e300] {
        let result = differentiate(&[amplitude, -2.0*amplitude], 2,
            CentralDifferenceControls::default(), |x| {
                Ok::<_, ()>(vec![3.0*x[0]-x[1], x[0]+x[1]])
            }, || Ok(())).unwrap();
        assert!((result.action[0]/amplitude-5.0).abs() < 2e-14);
        assert!((result.action[1]/amplitude+1.0).abs() < 2e-14);
    }
}

#[test]
fn g0_zero_direction_keeps_exact_zero_and_invalid_controls_still_refuse() {
    let result = differentiate(&[0.0, -0.0], 2, controls(), |_| -> Result<Vec<f64>, ()> {
        panic!("zero direction must not query the original source")
    }, || Ok(())).unwrap();
    assert_eq!(result.action, [0.0, 0.0]);
    assert_eq!(result.evaluations, 0);
    assert!(differentiate(&[0.0], 1, CentralDifferenceControls { initial_step: 0.0, ..controls() },
        |_| Ok::<_, ()>(vec![0.0]), || Ok(())).is_err());
}

#[test]
fn g4_original_source_refusals_and_work_ceilings_do_not_publish_partial_actions() {
    let error = differentiate(&[1.0], 1, controls(), |x| {
        if x[0] < 0.0 { Err("left chart") } else { Ok(vec![x[0]]) }
    }, || Ok(())).unwrap_err();
    assert_eq!(error, CentralDifferenceError::Source("left chart"));
    let mut limited = controls(); limited.maximum_evaluations = 4;
    let error = differentiate(&[1.0], 1, limited,
        |x| Ok::<_, ()>(vec![x[0].sin()]), || Ok(())).unwrap_err();
    assert_eq!(error, CentralDifferenceError::EvaluationBudget { used: 4, maximum: 4 });
}

#[test]
fn g4_final_cancellation_is_checked_after_a_resolved_action() {
    let calls = Cell::new(0);
    let result = differentiate(&[1.0], 1, controls(), |x| {
        calls.set(calls.get()+1); Ok::<_, &str>(vec![x[0]])
    }, || if calls.get()==4 { Err("cancel") } else { Ok(()) });
    assert_eq!(result, Err(CentralDifferenceError::Source("cancel")));
}

#[test]
fn g0_bad_source_dimensions_and_erased_offsets_refuse() {
    assert!(matches!(differentiate(&[1.0], 1, controls(),
        |_| Ok::<_, ()>(vec![f64::NAN]), || Ok(())), Err(CentralDifferenceError::Invalid(_))));
    assert!(matches!(differentiate(&[1.0, f64::from_bits(1)], 1, controls(),
        |_| Ok::<_, ()>(vec![0.0]), || Ok(())), Err(CentralDifferenceError::Invalid(_))));
    assert!(matches!(differentiate(&[1.0], 1, controls(),
        |_| Ok::<_, ()>(Vec::new()), || Ok(())), Err(CentralDifferenceError::Invalid(_))));
    assert!(matches!(differentiate(&[1e-300], 1, controls(),
        |x| Ok::<_, ()>(vec![1e-300*x[0]]), || Ok(())), Err(CentralDifferenceError::Invalid(_))));
}
