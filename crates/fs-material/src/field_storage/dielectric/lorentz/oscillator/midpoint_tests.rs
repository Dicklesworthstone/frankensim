use super::*;

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale.max(1e-305),
        "actual {actual:.17e}, expected {expected:.17e}, scale {scale:.17e}"
    );
}

#[test]
fn g1_condensation_recovery_and_independent_original_midpoint_equations_agree() {
    let before = LorentzState {
        polarization_c_m2: [0.2, -0.1, 0.07],
        current_a_m2: [0.3, -0.8, 1.1],
    };
    let e0 = [0.8, -0.2, 0.05];
    let e1 = [-0.3, 0.6, 0.4];
    for gamma in [0.0, 0.03, 14.0, 40.0] {
        let pole = LorentzPole::try_new(0.7, 7.0, gamma).unwrap();
        for dt in [1e-6, 0.03, 0.4, 10.0, 1e3] {
            let prepared = pole.prepare_midpoint(dt).unwrap();
            let frame = prepared.advance(before, e0, e1).unwrap();
            let history = prepared.history_c_m2(e0, before).unwrap();
            let gain = prepared.permittivity_increment_f_m().unwrap();
            let h = 0.5 * dt;
            let d = 1.0 + h * gamma + (h * 7.0).powi(2);
            near(gain, 0.7 * (h * 7.0).powi(2) / d, gain);
            assert!(gain >= 0.0 && gain <= 0.7);
            let mut heat = 0.0;
            let mut work = 0.0;
            for i in 0..3 {
                let em = 0.5 * (e0[i] + e1[i]);
                let jm = (before.current_a_m2[i]
                    + h * 49.0 * (0.7 * em - before.polarization_c_m2[i]))
                    / d;
                let p1 = before.polarization_c_m2[i] + dt * jm;
                let j1 = 2.0 * jm - before.current_a_m2[i];
                near(frame.state.polarization_c_m2[i], p1, 2.0);
                near(frame.state.current_a_m2[i], j1, 10.0);
                near(
                    frame.state.polarization_c_m2[i],
                    gain * e1[i] + history[i],
                    2.0,
                );
                heat += dt * gamma * jm * jm / (0.7 * 49.0);
                work += em * (frame.state.polarization_c_m2[i] - before.polarization_c_m2[i]);
            }
            near(frame.heat_j_m3, heat, frame.stored_before_j_m3 + heat);
            near(
                frame.supplied_work_j_m3,
                work,
                frame.stored_before_j_m3 + work.abs(),
            );
            near(
                frame.stored_change_j_m3,
                frame.stored_after_j_m3 - frame.stored_before_j_m3,
                frame.stored_before_j_m3 + frame.stored_after_j_m3,
            );
            near(
                frame.energy_residual_j_m3,
                0.0,
                frame.stored_before_j_m3 + frame.stored_after_j_m3 + heat + work.abs(),
            );
            assert_eq!(frame, prepared.advance(before, e0, e1).unwrap());
        }
    }
}

#[test]
fn g1_lossless_ringdown_has_cayley_phase_and_passive_damped_steps() {
    let pole = LorentzPole::try_new(2.0, 7.0, 0.0).unwrap();
    let root = 2.0_f64.sqrt();
    let initial = LorentzState {
        polarization_c_m2: [0.8 * root],
        current_a_m2: [-0.3 * 7.0 * root],
    };
    let mut state = initial;
    let dt = 0.01;
    let prepared = pole.prepare_midpoint(dt).unwrap();
    for _ in 0..4096 {
        let step = prepared.advance(state, [0.0], [0.0]).unwrap();
        assert_eq!(step.heat_j_m3, 0.0);
        assert_eq!(step.supplied_work_j_m3, 0.0);
        state = step.state;
    }
    // Independent exact solution of the discrete Cayley map, not omega0*dt.
    let angle = 4096.0 * 2.0 * (0.5 * 7.0 * dt).atan();
    let (s, c) = angle.sin_cos();
    near(state.polarization_c_m2[0], root * (0.8 * c - 0.3 * s), root);
    near(
        state.current_a_m2[0],
        7.0 * root * (-0.3 * c - 0.8 * s),
        7.0 * root,
    );
    near(
        pole.stored_energy_j_m3(state).unwrap(),
        pole.stored_energy_j_m3(initial).unwrap(),
        1.0,
    );

    for gamma in [0.4, 14.0, 40.0] {
        let pole = LorentzPole::try_new(2.0, 7.0, gamma).unwrap();
        for dt in [1e-4, 0.1, 1.0, 1e3] {
            let prepared = pole.prepare_midpoint(dt).unwrap();
            let mut state = initial;
            for _ in 0..64 {
                let step = prepared.advance(state, [0.0], [0.0]).unwrap();
                assert!(step.heat_j_m3 >= 0.0);
                assert!(step.stored_after_j_m3 <= step.stored_before_j_m3 * (1.0 + 2e-15));
                near(
                    step.stored_change_j_m3 + step.heat_j_m3,
                    0.0,
                    step.stored_before_j_m3,
                );
                state = step.state;
            }
        }
    }
}

#[test]
fn g1_refined_midpoint_converges_to_the_actual_exact_linear_field_owner() {
    let pole = LorentzPole::try_new(0.7, 7.0, 0.4).unwrap();
    let initial = LorentzState {
        polarization_c_m2: [0.07],
        current_a_m2: [-0.2],
    };
    let duration = 0.7;
    let exact = pole
        .prepare_step(duration)
        .unwrap()
        .advance(initial, [0.2], [1.5])
        .unwrap();
    let mut errors = Vec::new();
    for steps in [8, 16, 32, 64] {
        let dt = duration / f64::from(steps);
        let prepared = pole.prepare_midpoint(dt).unwrap();
        let mut state = initial;
        let mut heat = 0.0;
        let mut work = 0.0;
        for i in 0..steps {
            let e0 = 0.2 + 1.3 * f64::from(i) / f64::from(steps);
            let e1 = 0.2 + 1.3 * f64::from(i + 1) / f64::from(steps);
            let frame = prepared.advance(state, [e0], [e1]).unwrap();
            state = frame.state;
            heat += frame.heat_j_m3;
            work += frame.supplied_work_j_m3;
        }
        errors.push(
            (state.polarization_c_m2[0] - exact.state.polarization_c_m2[0])
                .hypot((state.current_a_m2[0] - exact.state.current_a_m2[0]) / 7.0),
        );
        near(
            pole.stored_energy_j_m3(state).unwrap() - pole.stored_energy_j_m3(initial).unwrap()
                + heat,
            work,
            work.abs(),
        );
    }
    for pair in errors.windows(2) {
        assert!(
            (3.6..4.3).contains(&(pair[0] / pair[1])),
            "second-order errors {errors:?}"
        );
    }
}

#[test]
fn g1_extreme_dimensional_products_preserve_condensation_state_and_heat() {
    let reference = LorentzPole::try_new(1.0, 1.0, 1.0)
        .unwrap()
        .prepare_midpoint(1.0)
        .unwrap()
        .advance(LorentzState::default(), [1.0], [1.0])
        .unwrap();
    for (delta, omega, dt, field) in [
        (1e-300, 1e200, 1e-200, 1e150),
        (1e300, 1e-200, 1e200, 1e-150),
    ] {
        let pole = LorentzPole::try_new(delta, omega, omega).unwrap();
        let prepared = pole.prepare_midpoint(dt).unwrap();
        let frame = prepared
            .advance(LorentzState::default(), [field], [field])
            .unwrap();
        near(
            prepared.permittivity_increment_f_m().unwrap(),
            delta / 7.0,
            delta,
        );
        near(
            frame.state.polarization_c_m2[0] / delta.sqrt(),
            reference.state.polarization_c_m2[0],
            1.0,
        );
        near(
            frame.state.current_a_m2[0] / (omega * delta.sqrt()),
            reference.state.current_a_m2[0],
            1.0,
        );
        near(frame.stored_after_j_m3, reference.stored_after_j_m3, 1.0);
        near(frame.heat_j_m3, reference.heat_j_m3, 1.0);
        near(frame.supplied_work_j_m3, reference.supplied_work_j_m3, 1.0);
    }
    // dt*omega and its square overflow, while every retained state is finite.
    // Tiny midpoint current must survive until dt and gamma restore its effect.
    let pole = LorentzPole::try_new(1e-300, 1e200, 1e300).unwrap();
    let before = LorentzState {
        polarization_c_m2: [3e-151],
        current_a_m2: [4e49],
    };
    let frame = pole
        .prepare_midpoint(1e200)
        .unwrap()
        .advance(before, [0.0], [0.0])
        .unwrap();
    near(
        frame.state.polarization_c_m2[0],
        -before.polarization_c_m2[0],
        3e-151,
    );
    near(frame.state.current_a_m2[0], -before.current_a_m2[0], 4e49);
    near(frame.heat_j_m3, 3.6e-301, 3.6e-301);
    assert!(frame.heat_j_m3 > 0.0);

    // gamma*dt overflows but the positive condensed coefficient stays finite.
    let prepared = LorentzPole::try_new(1e300, 1e-200, 1e300)
        .unwrap()
        .prepare_midpoint(1e200)
        .unwrap();
    near(
        prepared.permittivity_increment_f_m().unwrap(),
        5e-201,
        5e-201,
    );
    let frame = prepared
        .advance(LorentzState::default(), [1e150], [1e150])
        .unwrap();
    near(frame.state.polarization_c_m2[0], 1e-50, 1e-50);
    near(frame.state.current_a_m2[0], 2e-250, 2e-250);
    near(frame.heat_j_m3, 1e100, 1e100);
    near(frame.supplied_work_j_m3, 1e100, 1e100);
}

#[test]
fn g1_instantaneous_field_products_keep_scale_and_small_endpoint_changes() {
    near(field_energy(1e-300, [1e150, -1e150]).unwrap(), 1.0, 1.0);
    near(
        field_energy(1e300, [1e-300, -1e-300]).unwrap(),
        1e-300,
        1e-300,
    );
    let before = 1e150_f64;
    let after = f64::from_bits(before.to_bits() + 1);
    let expected = (1e-300 * (after - before)) * f64::midpoint(before, after);
    near(
        field_increment(1e-300, [before], [after]).unwrap(),
        expected,
        expected,
    );
}

#[test]
fn g0_invalid_intervals_vectors_and_unrepresentable_recovery_refuse_transactionally() {
    let pole = LorentzPole::try_new(0.7, 7.0, 0.4).unwrap();
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(pole.prepare_midpoint(dt).is_err());
    }
    let prepared = pole.prepare_midpoint(0.1).unwrap();
    assert!(
        prepared
            .advance(LorentzState::<0>::default(), [], [])
            .is_err()
    );
    assert!(
        prepared
            .advance(LorentzState::<4>::default(), [0.0; 4], [0.0; 4])
            .is_err()
    );
    assert!(
        prepared
            .history_c_m2([f64::NAN], LorentzState::default())
            .is_err()
    );
    let before = LorentzState {
        polarization_c_m2: [0.2, -0.1],
        current_a_m2: [0.3, -0.8],
    };
    let retained = before;
    assert!(prepared.advance(before, [0.0, 0.0], [0.0, 1e300]).is_err());
    assert_eq!(before, retained);
    let invalid = LorentzState {
        polarization_c_m2: [0.0],
        current_a_m2: [f64::NAN],
    };
    assert!(prepared.history_c_m2([0.0], invalid).is_err());
}
