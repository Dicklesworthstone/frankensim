use super::super::TsDriver;
use super::*;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn card() -> TsCard {
    // Deliberately slow, well-scaled analytic test fixture, not a product card.
    TsCard {
        re_ohm: 0.8,
        le_h: 0.6,
        bl: 1.2,
        mms_kg: 0.7,
        cms_m_per_n: 0.9,
        rms_n_s_m: 0.1,
        sd_m2: 0.001,
        x_max_m: 0.1,
    }
}

fn profile(nonlinear: bool) -> VoiceCoilMagnetics {
    VoiceCoilMagnetics::new(
        [-0.05, 0.05],
        0.01,
        if nonlinear {
            [1.2, -0.4, -0.2]
        } else {
            [1.2, 0.0, 0.0]
        },
        if nonlinear {
            [0.6, 0.15, 0.08]
        } else {
            [0.6, 0.0, 0.0]
        },
    )
    .unwrap()
}

#[test]
fn magnetic_force_current_and_back_emf_are_derivatives_of_the_same_energy() {
    let profile = VoiceCoilMagnetics::new(
        [-0.003, 0.003],
        0.002,
        [5.0, 1.0, -0.25],
        [0.002, 0.0005, 0.0002],
    )
    .unwrap();
    let storage = CoilStorage {
        mechanical: QuadraticStorage::new(vec![0.0; 9], 3).unwrap(),
        profile,
        exact_q: None,
    };
    let x = 0.001;
    let current = 0.7;
    let local = profile.sample(x).unwrap();
    let state = [
        local.permanent_linkage_wb + local.inductance_h * current,
        x,
        0.0,
    ];
    let mut gradient = [0.0; 3];
    storage.gradient(&state, &mut gradient);
    let force = current * (5.0 + 0.5 - 0.25 * 0.25)
        + 0.5 * current * current * (0.0005 + 2.0 * 0.0002 * 0.5) / 0.002;
    near(gradient[0], current, 1e-14);
    near(-gradient[1], force, 1e-13);
    for (axis, delta) in [(0, 1e-8), (1, 1e-9)] {
        let mut a = state;
        let mut b = state;
        a[axis] -= delta;
        b[axis] += delta;
        near(
            (storage.hamiltonian(&b) - storage.hamiltonian(&a)) / (2.0 * delta),
            gradient[axis],
            2e-8,
        );
    }
    // Independently differentiate lambda = L(x)i + psi(x), including the
    // often omitted L'(x)i*v voltage. Its power splits into magnetic storage
    // and BOTH permanent-field and reluctance mechanical work.
    let velocity = -0.2;
    let current_rate = 0.4;
    let voltage = local.inductance_h * current_rate
        + (local.force_factor_tm + local.inductance_gradient_h_per_m * current) * velocity;
    let magnetic_rate = local.inductance_h * current * current_rate
        + 0.5 * local.inductance_gradient_h_per_m * current * current * velocity;
    near(voltage * current, magnetic_rate + force * velocity, 2e-15);

    // Reluctance force is even in current, permanent-field force is odd.
    let negative_force = -current * local.force_factor_tm
        + 0.5 * current * current * local.inductance_gradient_h_per_m;
    near(
        force + negative_force,
        current * current * local.inductance_gradient_h_per_m,
        2e-14,
    );
}

#[test]
fn magnetic_profile_checks_interior_inductance_and_refuses_extrapolation() {
    assert!(
        VoiceCoilMagnetics::new([-0.001, 0.001], 0.001, [1.0, 0.0, 0.0], [1e-3, -6e-3, 6e-3])
            .is_err()
    );
    // Both end values and the origin are positive; the minimum at xi=0.5
    // is negative. Testing only knots would admit an active magnetic store.
    let admitted = profile(true);
    assert!(admitted.sample(-0.05000001).is_err());
    assert!(admitted.sample(0.05000001).is_err());
    assert!(admitted.sample(f64::NAN).is_err());
    assert!(VoiceCoilMagnetics::new([-0.1, 0.1], 0.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]).is_err());
    assert!(
        VoiceCoilMagnetics::new([-0.1, 0.1], 0.01, [1.0, f64::NAN, 0.0], [1.0, 0.0, 0.0]).is_err()
    );
    let mismatch =
        VoiceCoilMagnetics::new([-0.05, 0.05], 0.01, [1.3, 0.0, 0.0], [0.6, 0.0, 0.0]).unwrap();
    assert!(MagneticTsDriver::new(card(), None, None, mismatch).is_err());
    let too_far =
        VoiceCoilMagnetics::new([-0.2, 0.2], 0.01, [1.2, 0.0, 0.0], [0.6, 0.0, 0.0]).unwrap();
    assert!(MagneticTsDriver::new(card(), None, None, too_far).is_err());
}

#[test]
fn constant_magnetic_profile_converges_to_original_exact_linear_speaker() {
    let mut exact = TsDriver::new(card(), None, None).unwrap();
    exact.step(0.4, 0.2).unwrap();
    let expected = [exact.coil_current_a(), exact.excursion_m()];
    let mut prior_error: Option<[f64; 2]> = None;
    for count in [10, 20, 40] {
        let mut driver = MagneticTsDriver::new(card(), None, None, profile(false)).unwrap();
        for _ in 0..count {
            driver.step(0.4, 0.2 / f64::from(count)).unwrap();
        }
        let error = [
            (driver.coil_current_a() - expected[0]).abs(),
            (driver.excursion_m() - expected[1]).abs(),
        ];
        if let Some(previous) = prior_error {
            for axis in 0..2 {
                let ratio = previous[axis] / error[axis];
                assert!(
                    ratio > 3.8 && ratio < 4.2,
                    "component {axis} must converge at second order: ratio {ratio}"
                );
            }
        }
        prior_error = Some(error);
    }
    // Accuracy belongs to the refined solution; the deliberately coarse
    // steps establish the order in both electrical and mechanical variables.
    let finest = prior_error.unwrap();
    assert!(finest[0] < 1e-7);
    assert!(finest[1] < 4e-7);
}

#[test]
fn nonlinear_magnetic_speaker_and_cabinet_keep_native_work_and_heat_balance() {
    let panel = PanelMode {
        frequency_hz: 2.0,
        mass_kg: 0.05,
        area_m2: 0.001,
        damping_ratio: 0.02,
    };
    let mut driver = MagneticTsDriver::new(card(), Some(0.02), Some(panel), profile(true)).unwrap();
    let mut constant =
        MagneticTsDriver::new(card(), Some(0.02), Some(panel), profile(false)).unwrap();
    assert!(driver.has_panel());
    let mut supplied = 0.0;
    let mut dissipated = 0.0;
    for tick in 0..100 {
        let voltage = if tick < 60 { 0.4 } else { 0.0 };
        let (record, pressure) = driver.step(voltage, 0.005).unwrap();
        constant.step(voltage, 0.005).unwrap();
        assert!(pressure.is_finite());
        assert!(record.dissipated >= 0.0);
        near(record.delta_h + record.dissipated, record.supplied, 3e-12);
        supplied += record.supplied;
        dissipated += record.dissipated;
        near(driver.stored_energy_j() + dissipated, supplied, 3e-11);
        let local = profile(true).sample(driver.excursion_m()).unwrap();
        let current = driver.coil_current_a();
        near(
            driver.magnetic_force_n(),
            current * local.force_factor_tm
                + 0.5 * current * current * local.inductance_gradient_h_per_m,
            1e-14,
        );
    }
    assert!(
        (driver.excursion_m() - constant.excursion_m()).abs() > 1e-5,
        "nonlinear source coefficients must change actual cone motion"
    );
    assert!(driver.stored_energy_j() > 0.0 && dissipated > 0.0);
}

#[test]
fn magnetic_domain_and_nonfinite_input_failures_leave_the_driver_unchanged() {
    let narrow =
        VoiceCoilMagnetics::new([-1e-6, 1e-6], 0.01, [1.2, -0.4, -0.2], [0.6, 0.15, 0.08]).unwrap();
    let mut driver = MagneticTsDriver::new(card(), None, None, narrow).unwrap();
    let before = driver.x.clone();
    assert!(driver.step(1.0, 0.1).is_err());
    assert_eq!(driver.x, before);
    assert!(driver.step(f64::NAN, 1e-6).is_err());
    assert_eq!(driver.x, before);
    assert!(driver.step(1.0, -1.0).is_err());
    assert_eq!(driver.x, before);
    driver.step(0.01, 1e-5).unwrap();
    assert!(driver.coil_current_a() > 0.0);
}

#[test]
fn magnetic_profile_coordinate_scale_changes_no_physics() {
    let a = profile(true);
    let b =
        VoiceCoilMagnetics::new([-0.05, 0.05], 0.02, [1.2, -0.8, -0.8], [0.6, 0.3, 0.32]).unwrap();
    for x in [-0.01, 0.0, 0.003, 0.02] {
        let first = a.sample(x).unwrap();
        let second = b.sample(x).unwrap();
        near(
            first.permanent_linkage_wb,
            second.permanent_linkage_wb,
            1e-15,
        );
        near(first.force_factor_tm, second.force_factor_tm, 1e-14);
        near(first.inductance_h, second.inductance_h, 1e-14);
        near(
            first.inductance_gradient_h_per_m,
            second.inductance_gradient_h_per_m,
            1e-12,
        );
    }
}
