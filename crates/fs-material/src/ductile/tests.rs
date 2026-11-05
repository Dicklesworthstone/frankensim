use super::*;

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "actual {actual:e}, expected {expected:e}");
}

fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-1.0, 2.0],
        strain_rate_s_inv: [1e-8, 1e8],
        temperature_k: [300.0, 300.0],
    }
}

fn constant_locus() -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(vec![(-1.0, 0.25), (2.0, 0.25)], domain()).unwrap()
}

fn jc_parameters() -> JohnsonCookParameters {
    JohnsonCookParameters {
        d1: 0.05, d2: 0.5, d3: 1.5, d4: 0.03, d5: 0.4,
        reference_strain_rate_s_inv: 1.0,
        transition_temperature_k: 300.0,
        melting_temperature_k: 1500.0,
    }
}

#[test]
fn g0_johnson_cook_formula_sign_rate_and_temperature() {
    let law = JohnsonCookInitiation::new(jc_parameters(), DuctileDomain {
        temperature_k: [200.0, 1600.0], ..domain()
    }).unwrap();
    for (eta, rate, temperature, theta) in [
        (0.0, 1.0, 200.0, 0.0),
        (1.0 / 3.0, 10.0, 300.0, 0.0),
        (0.8, 0.1, 900.0, 0.5),
        (-0.3, 2.0, 1600.0, 1.0),
    ] {
        let expected = (0.05 + 0.5 * (-1.5_f64 * eta).exp())
            * (1.0 + 0.03 * f64::ln(rate)) * (1.0 + 0.4 * theta);
        close(law.fracture_strain(eta, rate, temperature).unwrap(), expected);
    }
    assert!(law.fracture_strain(1.0, 1.0, 300.0).unwrap()
        < law.fracture_strain(0.0, 1.0, 300.0).unwrap());
    assert!(law.fracture_strain(0.0, 10.0, 300.0).unwrap()
        > law.fracture_strain(0.0, 1.0, 300.0).unwrap());
    assert!(law.fracture_strain(0.0, 1.0, 900.0).unwrap()
        > law.fracture_strain(0.0, 1.0, 300.0).unwrap());
}

#[test]
fn g0_signed_coefficients_require_positive_factors_over_the_entire_domain() {
    let p = JohnsonCookParameters {
        d1: -0.2, d2: 0.8, d3: 0.5, d4: -0.05, d5: -0.2,
        ..jc_parameters()
    };
    let support = DuctileDomain {
        triaxiality: [0.0, 2.0], strain_rate_s_inv: [0.1, 10.0],
        temperature_k: [300.0, 1500.0],
    };
    let law = JohnsonCookInitiation::new(p, support).unwrap();
    assert!(law.fracture_strain(2.0, 10.0, 1500.0).unwrap() > 0.0);
    assert!(JohnsonCookInitiation::new(p, DuctileDomain {
        triaxiality: [0.0, 3.0], ..support
    }).is_err());
    assert!(JohnsonCookInitiation::new(p, DuctileDomain {
        strain_rate_s_inv: [0.1, 1e300], ..support
    }).is_err());
    assert!(JohnsonCookInitiation::new(JohnsonCookParameters { d5: -1.0, ..p }, support).is_err());
    assert!(JohnsonCookInitiation::new(JohnsonCookParameters {
        d1: -1.0, d2: 0.0, d5: -2.0, ..p
    }, DuctileDomain { temperature_k: [1500.0, 1500.0], ..support }).is_err());
    assert!(JohnsonCookInitiation::new(JohnsonCookParameters { d3: f64::NAN, ..p }, support).is_err());
    assert!(JohnsonCookInitiation::new(p, DuctileDomain {
        strain_rate_s_inv: [0.0, 1.0], ..support
    }).is_err());
}

#[test]
fn g0_tabulated_locus_keeps_original_knots_and_refuses_extrapolation() {
    let support = DuctileDomain {
        triaxiality: [-1.0, 1.0], strain_rate_s_inv: [0.1, 0.1], ..domain()
    };
    let knots = vec![(-1.0, 0.8), (0.0, 0.5), (1.0, 0.2)];
    let law = TabulatedDuctileLocus::new(knots.clone(), support).unwrap();
    assert_eq!(law.knots(), &knots);
    for &(eta, strain) in &knots {
        assert_eq!(law.fracture_strain(eta, 0.1, 300.0).unwrap(), strain);
    }
    close(law.fracture_strain(0.5, 0.1, 300.0).unwrap(), 0.35);
    close(law.fracture_strain(-0.5, 0.1, 300.0).unwrap(), 0.65);
    for (eta, rate, temperature) in [(1.01, 0.1, 300.0), (0.0, 0.2, 300.0), (0.0, 0.1, 301.0)] {
        assert!(matches!(law.fracture_strain(eta, rate, temperature), Err(DuctileError::Domain { .. })));
    }
    for bad in [vec![], vec![(0.0, 0.2)], vec![(0.0, 0.2), (0.0, 0.3)],
        vec![(-1.0, 0.2), (1.0, 0.0)]] {
        assert!(TabulatedDuctileLocus::new(bad, support).is_err());
    }
    assert!(TabulatedDuctileLocus::new(knots, DuctileDomain {
        triaxiality: [-2.0, 1.0], ..support
    }).is_err());
}

#[test]
fn g1_constant_stress_path_reaches_analytic_initiation_and_retains_overshoot() {
    let law = constant_locus();
    let stress = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let original = DuctileState::virgin();
    let first = original.trial(&law, 0.125, 1.0, 300.0, &stress).unwrap();
    assert_eq!(first.state.utilization(), 0.5);
    assert!(!first.state.initiated());
    assert_eq!(first.initiation_fraction, None);
    let crossing = first.state.trial(&law, 0.375, 2.0, 300.0, &stress).unwrap();
    assert_eq!(crossing.state.utilization(), 1.5);
    assert!(crossing.state.initiated());
    assert_eq!(crossing.initiation_fraction, Some(0.5));
    assert_eq!(crossing.state.accumulated_plastic_strain(), 0.375);
    assert_eq!(crossing.plastic_increment, 0.25);
    assert_eq!(crossing.strain_rate_s_inv, 0.125);
    assert_eq!(crossing.fracture_strain, Some(0.25));
    assert_eq!(original, DuctileState::virgin());
    assert_eq!(first.state.utilization(), 0.5);
}

#[test]
fn g3_triaxiality_changes_capacity_and_unloading_or_reversal_never_heals() {
    let law = TabulatedDuctileLocus::new(
        vec![(-1.0, 0.6), (0.0, 0.4), (1.0, 0.2), (2.0, 0.1)], domain()
    ).unwrap();
    let zero = DuctileState::virgin();
    let uniaxial = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let confined = [200.0, 100.0, 100.0, 0.0, 0.0, 0.0];
    let tensile = zero.trial(&law, 0.0625, 1.0, 300.0, &uniaxial).unwrap();
    let confined = zero.trial(&law, 0.0625, 1.0, 300.0, &confined).unwrap();
    assert!(confined.state.utilization() > tensile.state.utilization());
    let unload = tensile.state.trial(&law, 0.0625, 1.0, 300.0, &[0.0; 6]).unwrap();
    assert_eq!(unload.state, tensile.state);
    assert_eq!(unload.increment_utilization, 0.0);
    assert_eq!(unload.triaxiality, None);
    let reverse = unload.state.trial(&law, 0.125, 1.0, 300.0, &uniaxial.map(|x| -x)).unwrap();
    assert!(reverse.state.utilization() > unload.state.utilization());
    assert!(reverse.triaxiality.unwrap() < 0.0);
}

#[test]
fn g3_triaxiality_is_rotation_and_positive_stress_scale_invariant() {
    close(stress_triaxiality(&[3.0, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap().unwrap(), 1.0/3.0);
    close(stress_triaxiality(&[-3.0, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap().unwrap(), -1.0/3.0);
    assert_eq!(stress_triaxiality(&[0.0, 0.0, 0.0, 7.0, 0.0, 0.0]).unwrap(), Some(0.0));
    assert_eq!(stress_triaxiality(&[7.0, 7.0, 7.0, 0.0, 0.0, 0.0]).unwrap(), None);
    let stress = [3.0, 2.0, 1.0, 0.2, -0.3, 0.4];
    let expected = 2.0 / (3.0_f64 + 3.0 * (0.2*0.2 + 0.3*0.3 + 0.4*0.4)).sqrt();
    close(stress_triaxiality(&stress).unwrap().unwrap(), expected);
    for factor in [1e-300, 1e300] {
        close(stress_triaxiality(&stress.map(|x| x*factor)).unwrap().unwrap(), expected);
    }
    let c = core::f64::consts::FRAC_1_SQRT_2;
    let rotation = [[c,-c,0.0],[c,c,0.0],[0.0,0.0,1.0]];
    let rotated = crate::tensor::rotate(&stress, &rotation);
    close(stress_triaxiality(&rotated).unwrap().unwrap(), expected);
    // Directly squaring this normalized shear would underflow; eta is finite.
    let tiny_deviator = [1.0, 1.0, 1.0, 1e-200, 0.0, 0.0];
    close(stress_triaxiality(&tiny_deviator).unwrap().unwrap(), 1e200 / 3.0_f64.sqrt());
}

#[test]
fn g0_invalid_or_unsupported_plastic_trials_leave_accepted_history_unchanged() {
    let law = constant_locus();
    let stress = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let accepted = DuctileState::virgin().trial(&law, 0.0625, 1.0, 300.0, &stress).unwrap().state;
    let original = accepted;
    assert!(matches!(accepted.trial(&law, 0.0, 1.0, 300.0, &stress),
        Err(DuctileError::NonMonotonePlasticStrain { .. })));
    assert!(matches!(accepted.trial(&law, 0.125, 1.0, 300.0, &[0.0; 6]),
        Err(DuctileError::HydrostaticPlasticFlow)));
    for (alpha, time, temperature) in [
        (0.125, 0.0, 300.0), (0.125, f64::NAN, 300.0),
        (0.125, 1.0, 301.0), (f64::NAN, 1.0, 300.0),
        (0.125, 1e20, 300.0), (0.0625, 1.0, 301.0),
    ] {
        assert!(accepted.trial(&law, alpha, time, temperature, &stress).is_err());
    }
    assert!(matches!(accepted.trial(&law, 0.125, 1.0, 300.0,
        &[1100.0,1000.0,1000.0,0.0,0.0,0.0]), Err(DuctileError::Domain { .. })));
    assert!(accepted.trial(&law, 0.125, 1.0, 300.0, &[f64::NAN; 6]).is_err());
    assert_eq!(accepted, original);
    assert_eq!(accepted.trial(&law, 0.0625, 1.0, 300.0, &[4.0,4.0,4.0,0.0,0.0,0.0])
        .unwrap().state, original);
}

#[test]
fn g3_constant_condition_splitting_and_replay_preserve_native_accumulation() {
    let law = constant_locus();
    let stress = [0.0, 0.0, 0.0, 12.0, 0.0, 0.0];
    let whole = DuctileState::virgin().trial(&law, 0.125, 1.0, 300.0, &stress).unwrap().state;
    let mut split = DuctileState::virgin();
    for step in 1..=4 {
        let alpha = f64::from(step) / 32.0;
        let candidate = split.trial(&law, alpha, 0.25, 300.0, &stress).unwrap();
        assert_eq!(candidate, split.trial(&law, alpha, 0.25, 300.0, &stress).unwrap());
        split = candidate.state;
    }
    assert_eq!(split, whole);
}
