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
fn distinct_cavity_gas_changes_magnetic_motion_without_changing_exterior_observation() {
    let air = declared_air();
    let fill = GasState::try_new(
        &fs_material::gas_species::GasSpecies::Helium.spec(),
        air.temperature,
        air.pressure,
    )
    .unwrap();
    let panel = Some(PanelMode {
        frequency_hz: 2.0,
        mass_kg: 0.05,
        area_m2: 0.001,
        damping_ratio: 0.02,
    });
    let mut filled =
        MagneticTsDriver::new_with_media(card(), Some(0.02), panel, air, Some(fill), profile(true))
            .unwrap();
    let mut baseline =
        MagneticTsDriver::new_with_gas(card(), Some(0.02), panel, air, profile(true)).unwrap();
    for _ in 0..40 {
        let previous = filled.x[2];
        let (record, pressure) = filled.step(0.2, 0.005).unwrap();
        baseline.step(0.2, 0.005).unwrap();
        near(
            pressure,
            air.density * card().sd_m2 * (record.x[2] - previous)
                / (filled.moving_mass_kg * 0.005 * 4.0 * core::f64::consts::PI),
            1e-15,
        );
        near(record.delta_h + record.dissipated, record.supplied, 3e-12);
    }
    assert!((filled.excursion_m() - baseline.excursion_m()).abs() > 1e-7);
    let refreshed = filled.with_coil_resistance(1.1).unwrap();
    assert_eq!(refreshed.exterior_gas(), &air);
    assert_eq!(refreshed.cavity_gas(), Some(&fill));
    assert_eq!(refreshed.x, filled.x);
    assert_eq!(refreshed.time_s(), filled.time_s());
    assert_eq!(refreshed.storage.mechanical, filled.storage.mechanical);
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

fn chord_storage() -> CoilStorage {
    // Independent cone/panel inertia with a reciprocal cavity cross term.
    let mut q = vec![0.0; 25];
    q[6] = 4.0;
    q[12] = 2.0;
    q[18] = 6.0;
    q[24] = 3.0;
    q[8] = -1.0;
    q[16] = -1.0;
    CoilStorage {
        mechanical: QuadraticStorage::new(q, 5).unwrap(),
        profile: profile(true),
        exact_q: None,
    }
}

#[test]
fn magnetic_chord_keeps_local_work_and_independent_mechanical_efforts() {
    let storage = chord_storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let a = [0.2, -0.006, 0.5, 0.002, -0.2];
    let b = [-0.1, 0.011, -0.3, -0.001, 0.4];
    let mid: [f64; 5] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
    let mut effort = [0.0; 5];
    let mut reverse = [0.0; 5];
    rule.evaluate(&a, &b, &mid, &mut effort);
    rule.evaluate(&b, &a, &mid, &mut reverse);
    let la = storage.profile.sample(a[1]).unwrap();
    let lb = storage.profile.sample(b[1]).unwrap();
    let ia = (a[0] - la.permanent_linkage_wb) / la.inductance_h;
    let ib = (b[0] - lb.permanent_linkage_wb) / lb.inductance_h;
    let current = f64::midpoint(ia, ib);
    let force = current * (lb.permanent_linkage_wb - la.permanent_linkage_wb) / (b[1] - a[1])
        + 0.5 * ia * ib * (lb.inductance_h - la.inductance_h) / (b[1] - a[1]);
    near(effort[0], current, 2e-15);
    near(effort[1], 4.0 * mid[1] - mid[3] - force, 2e-14);
    let magnetic_change = 0.5 * lb.inductance_h * ib * ib - 0.5 * la.inductance_h * ia * ia;
    near(
        current * (b[0] - a[0]) - force * (b[1] - a[1]),
        magnetic_change,
        2e-15,
    );
    let work: f64 = (0..5).map(|i| effort[i] * (b[i] - a[i])).sum();
    near(work, storage.hamiltonian(&b) - storage.hamiltonian(&a), 2e-14);
    for i in 0..5 {
        near(effort[i], reverse[i], 2e-14);
    }
    near(effort[2], 2.0 * mid[2], 2e-15);
    near(effort[3], 6.0 * mid[3] - mid[1], 2e-15);
    near(effort[4], 3.0 * mid[4], 2e-15);

    let mut other_a = a;
    let mut other_b = b;
    other_a[2] = 10.0;
    other_b[2] = -4.0;
    other_a[4] = -6.0;
    other_b[4] = 9.0;
    let other_mid: [f64; 5] = std::array::from_fn(|i| f64::midpoint(other_a[i], other_b[i]));
    let mut other_effort = [0.0; 5];
    rule.evaluate(&other_a, &other_b, &other_mid, &mut other_effort);
    assert_eq!(effort[..2], other_effort[..2]);
    near(other_effort[2], 6.0, 2e-15);
    near(other_effort[4], 4.5, 2e-15);

    // The original mechanical energy still owns its sub-ulp midpoint remainder.
    let mut adjacent = a;
    adjacent[2] = 1.0;
    let mut next = adjacent;
    next[2] += f64::EPSILON;
    let rounded: [f64; 5] = std::array::from_fn(|i| f64::midpoint(adjacent[i], next[i]));
    near(
        storage
            .midpoint_energy_remainder(&adjacent, &next, &rounded)
            .unwrap(),
        f64::EPSILON * f64::EPSILON,
        1e-45,
    );
}

#[test]
fn magnetic_chord_analytic_endpoint_action_matches_the_actual_rule() {
    let storage = chord_storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let a = [0.2, -0.006, 0.5, 0.002, -0.2];
    for b in [
        [-0.1, 0.011, -0.3, -0.001, 0.4],
        [-0.1, a[1], -0.3, -0.001, 0.4],
        a,
    ] {
        for direction in [
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.3, -0.02, 0.4, -0.1, 0.2],
        ] {
            let mut exact = [0.0; 5];
            assert!(rule.end_jacobian_action(&a, &b, &direction, &mut exact));
            let h = 1e-7;
            let minus: [f64; 5] = std::array::from_fn(|i| b[i] - h * direction[i]);
            let plus: [f64; 5] = std::array::from_fn(|i| b[i] + h * direction[i]);
            let minus_mid: [f64; 5] = std::array::from_fn(|i| f64::midpoint(a[i], minus[i]));
            let plus_mid: [f64; 5] = std::array::from_fn(|i| f64::midpoint(a[i], plus[i]));
            let mut gm = [0.0; 5];
            let mut gp = [0.0; 5];
            rule.evaluate(&a, &minus, &minus_mid, &mut gm);
            rule.evaluate(&a, &plus, &plus_mid, &mut gp);
            for i in 0..5 {
                near(
                    (gp[i] - gm[i]) / (2.0 * h),
                    exact[i],
                    5e-7 * (1.0 + exact[i].abs()),
                );
            }
        }
    }
}

#[test]
fn constant_magnetic_chord_retains_the_exact_quadratic_solver_path() {
    let driver = MagneticTsDriver::new(card(), None, None, profile(false)).unwrap();
    assert!(driver.storage.discrete_gradient_rule().is_none());
    assert!(driver.phs.storage().discrete_gradient_rule().is_none());
    assert!(driver.phs.storage().exact_quadratic_matrix().is_some());
    let a = [0.2, -0.006, 0.5];
    let b = [-0.1, 0.011, -0.3];
    let mid: [f64; 3] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
    let mut ga = [0.0; 3];
    let mut gb = [0.0; 3];
    let mut chord = [0.0; 3];
    driver.storage.gradient(&a, &mut ga);
    driver.storage.gradient(&b, &mut gb);
    driver.storage.evaluate(&a, &b, &mid, &mut chord);
    for i in 0..3 {
        near(chord[i], f64::midpoint(ga[i], gb[i]), 2e-14);
    }
    let direction = [0.3, -0.02, 0.4];
    let mut tangent = [0.0; 3];
    let mut analytic = [0.0; 3];
    driver.storage.gradient(&direction, &mut tangent);
    assert!(
        driver
            .storage
            .end_jacobian_action(&a, &b, &direction, &mut analytic)
    );
    for i in 0..3 {
        near(analytic[i], 0.5 * tangent[i], 2e-14);
    }
}

#[test]
fn magnetic_chord_refuses_invalid_support_shape_and_unrepresentable_current() {
    let storage = chord_storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let a = [0.2, -0.006, 0.5, 0.002, -0.2];
    for (axis, value) in [
        (1, 0.051),
        (0, f64::NAN),
        (2, f64::INFINITY),
        (0, f64::MAX),
    ] {
        let mut b = a;
        b[axis] = value;
        let mid: [f64; 5] = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
        let mut out = [0.0; 5];
        rule.evaluate(&a, &b, &mid, &mut out);
        assert!(out.iter().all(|v| v.is_nan()));
        assert!(!rule.end_jacobian_action(&a, &b, &[1.0; 5], &mut out));
        assert!(
            storage
                .midpoint_energy_remainder(&a, &b, &mid)
                .unwrap()
                .is_nan()
        );
    }
    let mut out = [0.0; 5];
    rule.evaluate(&a[..2], &a, &a, &mut out);
    assert!(out.iter().all(|v| v.is_nan()));
    assert!(!rule.end_jacobian_action(&a, &a, &[1.0; 4], &mut out));
    let mut short = [0.0; 4];
    rule.evaluate(&a, &a, &a, &mut short);
    assert!(short.iter().all(|v| v.is_nan()));
}
