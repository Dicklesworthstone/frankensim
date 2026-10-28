//! Independent numerical checks of stress-controlled Maxwell histories.
//! Synthetic parameters verify the integrator, not material calibration.
#![deny(unsafe_code)]

use fs_material::visco::stress_control::StressControlledMaxwell;
use fs_material::visco::{GeneralizedMaxwell, PronyState};

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn stored(law: &GeneralizedMaxwell, state: &PronyState) -> f64 {
    0.5 * law.e_inf * state.strain * state.strain
        + law
            .terms
            .iter()
            .zip(&state.q)
            .map(|(&(e, _), q)| 0.5 * e * q * q)
            .sum::<f64>()
}
fn same_state(a: &PronyState, b: &PronyState) {
    assert_eq!(
        a.q.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        b.q.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
    );
    for (a, b) in [
        (a.strain, b.strain),
        (a.stress, b.stress),
        (a.work, b.work),
        (a.dissipated, b.dissipated),
    ] {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn g1_standard_linear_solid_creep_and_zero_load_recovery_match_closed_form() {
    let law = GeneralizedMaxwell::new(3.0, vec![(7.0, 2.0)]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut state = law.state();
    let stress = 1.25;
    let tau_creep = 2.0 * 10.0 / 3.0;
    let jump = control.step(&mut state, stress, 0.0).unwrap();
    close(jump.strain, stress / 10.0, 2e-14);
    assert_eq!(jump.incremental_dissipated_j_m3, 0.0);
    assert_eq!(
        jump.strain_rate_bounds_s,
        [f64::NEG_INFINITY, f64::INFINITY]
    );
    let initial = state.clone();
    let time = 4.5;
    let hold = control.step(&mut state, stress, time).unwrap();
    let expected = stress / 3.0 - stress * 7.0 / (3.0 * 10.0) * (-time / tau_creep).exp();
    close(state.strain, expected, 3e-14);
    close(
        hold.incremental_work_j_m3,
        stress * (expected - initial.strain),
        3e-14,
    );
    let heat = 7.0 * initial.q[0] * initial.q[0] * tau_creep / (2.0 * 2.0)
        * -(-2.0 * time / tau_creep).exp_m1();
    close(hold.incremental_dissipated_j_m3, heat, 3e-14);
    assert!(hold.strain_bounds[0] <= initial.strain && hold.strain_bounds[1] >= expected);
    control.step(&mut state, 0.0, 0.0).unwrap();
    let unloaded = state.clone();
    let recovery = control.step(&mut state, 0.0, 9.0).unwrap();
    close(
        state.strain,
        unloaded.strain * (-9.0 / tau_creep).exp(),
        5e-14,
    );
    assert_eq!(recovery.incremental_work_j_m3, 0.0);
    close(
        recovery.incremental_dissipated_j_m3,
        stored(&law, &unloaded) - stored(&law, &state),
        5e-14,
    );
    assert!(state.strain > 0.0 && state.strain < unloaded.strain);
    close(state.work, stored(&law, &state) + state.dissipated, 5e-14);
}

#[test]
fn g1_stress_ramps_match_independent_simpson_work_and_heat_on_both_series_sides() {
    let law = GeneralizedMaxwell::new(3.0, vec![(7.0, 2.0)]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut initial = law.state();
    law.step(&mut initial, 0.08, 0.3);
    let tau_creep = 20.0 / 3.0;
    for ratio in [1e-6, 0.09999, 0.10001, 0.9, 20.0] {
        let dt = ratio * tau_creep;
        let target = -0.3;
        let rate = (target - initial.stress) / dt;
        let q_at = |t: f64| {
            initial.q[0] * (-t / tau_creep).exp() + (rate * 2.0 / 3.0) * -(-t / tau_creep).exp_m1()
        };
        let mut work = 0.0;
        let mut heat = 0.0;
        let intervals = 10_000;
        for i in 0..=intervals {
            let t = dt * i as f64 / intervals as f64;
            let q = q_at(t);
            let eps_rate = (rate + 7.0 * q / 2.0) / 10.0;
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            work += weight * (initial.stress + rate * t) * eps_rate;
            heat += weight * 7.0 * q * q / 2.0;
        }
        work *= dt / (3.0 * intervals as f64);
        heat *= dt / (3.0 * intervals as f64);
        let mut state = initial.clone();
        let response = control.step(&mut state, target, dt).unwrap();
        close(state.q[0], q_at(dt), 1e-10);
        close(response.incremental_work_j_m3, work, 2e-9);
        close(response.incremental_dissipated_j_m3, heat, 2e-9);
        assert!(response.incremental_dissipated_j_m3 > 0.0);
        assert!(response.energy_balance_residual_j_m3.abs() < 1e-13);
    }
}

#[test]
fn g1_multiple_branches_subdivide_and_bound_interior_strain_and_rate() {
    let law = GeneralizedMaxwell::new(3.0, vec![(2.0, 0.05), (7.0, 0.8), (1.5, 8.0)]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut initial = law.state();
    law.step(&mut initial, 0.07, 0.2);
    let target = -0.18;
    let dt = 4.3;
    let mut whole = initial.clone();
    let report = control.step(&mut whole, target, dt).unwrap();
    let mut split = initial.clone();
    for i in 1..=997 {
        let stress = initial.stress + (target - initial.stress) * i as f64 / 997.0;
        control.step(&mut split, stress, dt / 997.0).unwrap();
    }
    close(split.strain, whole.strain, 2e-11);
    close(split.work, whole.work, 2e-11);
    close(split.dissipated, whole.dissipated, 2e-11);
    for (a, b) in split.q.iter().zip(&whole.q) {
        close(*a, *b, 2e-11);
    }
    let e0 = law.e_inf + law.terms.iter().map(|p| p.0).sum::<f64>();
    for i in 0..=37 {
        let fraction = i as f64 / 37.0;
        let mut at = initial.clone();
        control
            .step(
                &mut at,
                initial.stress + (target - initial.stress) * fraction,
                dt * fraction,
            )
            .unwrap();
        let actual_rate = ((target - initial.stress) / dt
            + law
                .terms
                .iter()
                .zip(&at.q)
                .map(|(&(e, tau), q)| e * q / tau)
                .sum::<f64>())
            / e0;
        assert!(
            at.strain >= report.strain_bounds[0] - 1e-13
                && at.strain <= report.strain_bounds[1] + 1e-13
        );
        assert!(
            actual_rate >= report.strain_rate_bounds_s[0] - 1e-13
                && actual_rate <= report.strain_rate_bounds_s[1] + 1e-13
        );
    }
    close(whole.work, stored(&law, &whole) + whole.dissipated, 2e-12);
    let mut replay = initial.clone();
    assert_eq!(report, control.step(&mut replay, target, dt).unwrap());
    same_state(&replay, &whole);
}

#[test]
fn g1_strain_and_stress_control_share_history_through_unloading_and_recovery() {
    let law = GeneralizedMaxwell::new(
        3050.0,
        vec![
            (23900.0, 0.01),
            (19700.0, 0.1),
            (8990.0, 1.0),
            (12400.0, 10.0),
            (8140.0, 100.0),
        ],
    )
    .unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut state = law.state();
    law.step(&mut state, 0.01, 0.2);
    let held_stress = state.stress;
    let previous_strain = state.strain;
    control.step(&mut state, held_stress, 10.0).unwrap();
    assert!(state.strain > previous_strain);
    control.step(&mut state, 0.0, 1.0).unwrap();
    let released_strain = state.strain;
    let before_heat = state.dissipated;
    let recovery = control.step(&mut state, 0.0, 100.0).unwrap();
    assert!(state.strain > 0.0 && state.strain < released_strain);
    assert_eq!(recovery.incremental_work_j_m3, 0.0);
    assert!(state.dissipated > before_heat);
    law.step(&mut state, 0.0, 1.0);
    assert!(state.stress < 0.0); // Clamping the recovered length is not zero-load recovery.
    close(state.work, stored(&law, &state) + state.dissipated, 2e-10);
}

#[test]
fn g1_elastic_and_instantaneous_limits_have_zero_viscous_heat() {
    let law = GeneralizedMaxwell::new(10.0, vec![]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut state = law.state();
    let response = control.step(&mut state, 3.0, 2.0).unwrap();
    close(state.strain, 0.3, 2e-15);
    close(response.incremental_work_j_m3, 0.45, 2e-15);
    assert_eq!(response.incremental_dissipated_j_m3, 0.0);
    assert_eq!(response.strain_rate_bounds_s, [0.15; 2]);
    control.step(&mut state, -3.0, 0.0).unwrap();
    close(state.work, 0.45, 2e-15);
    let law = GeneralizedMaxwell::new(3.0, vec![(7.0, 2.0)]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut state = law.state();
    law.step(&mut state, 0.08, 0.3);
    let before = state.clone();
    let target = -0.2;
    let response = control.step(&mut state, target, 0.0).unwrap();
    let delta = (target - before.stress) / 10.0;
    close(state.strain, before.strain + delta, 2e-14);
    close(state.q[0], before.q[0] + delta, 2e-14);
    assert_eq!(response.incremental_dissipated_j_m3, 0.0);
    close(
        response.incremental_work_j_m3,
        0.5 * (before.stress + target) * delta,
        2e-14,
    );
}

#[test]
fn g1_tiny_steps_keep_positive_heat_and_overflowing_rate_time_keeps_finite_creep_work() {
    let law = GeneralizedMaxwell::new(3.0, vec![(7.0, 2.0)]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut state = law.state();
    control.step(&mut state, 1.0, 0.0).unwrap();
    let response = control.step(&mut state, 1.0, 1e-12).unwrap();
    assert!(response.incremental_dissipated_j_m3 > 0.0);
    close(
        response.incremental_dissipated_j_m3,
        7.0 * 0.1 * 0.1 / 2.0 * 1e-12,
        1e-12,
    );
    let law = GeneralizedMaxwell::new(3.0, vec![(7.0, 1e-200)]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut state = law.state();
    control.step(&mut state, 1.0, 0.0).unwrap();
    let before = state.clone();
    let response = control.step(&mut state, 1.0, 1e200).unwrap();
    close(state.strain, 1.0 / 3.0, 2e-14);
    close(
        response.incremental_work_j_m3,
        1.0 / 3.0 - before.strain,
        3e-14,
    );
    assert!(response.incremental_dissipated_j_m3 > 0.0);
    close(state.work, stored(&law, &state) + state.dissipated, 3e-14);
}

#[test]
fn g0_bad_parameters_and_malformed_states_refuse_without_mutation() {
    for law in [
        GeneralizedMaxwell {
            e_inf: 0.0,
            terms: vec![],
        },
        GeneralizedMaxwell {
            e_inf: 1.0,
            terms: vec![(0.0, 1.0)],
        },
        GeneralizedMaxwell {
            e_inf: 1.0,
            terms: vec![(1.0, 0.0)],
        },
        GeneralizedMaxwell {
            e_inf: 1.0,
            terms: vec![(1.0, 1.0); 65],
        },
    ] {
        assert!(StressControlledMaxwell::new(&law).is_err());
    }
    let law = GeneralizedMaxwell::new(3.0, vec![(7.0, 2.0)]).unwrap();
    let control = StressControlledMaxwell::new(&law).unwrap();
    let mut valid = law.state();
    law.step(&mut valid, 0.1, 0.3);
    for (target, dt) in [
        (f64::NAN, 1.0),
        (1.0, -0.1),
        (1.0, f64::INFINITY),
        (1.0, f64::NAN),
    ] {
        let mut state = valid.clone();
        assert!(control.step(&mut state, target, dt).is_err());
        same_state(&state, &valid);
    }
    let mut invalid = Vec::new();
    let mut state = valid.clone();
    state.q.clear();
    invalid.push(state);
    let mut state = valid.clone();
    state.q[0] = f64::NAN;
    invalid.push(state);
    let mut state = valid.clone();
    state.stress += 1.0;
    invalid.push(state);
    let mut state = valid.clone();
    state.dissipated = -1.0;
    invalid.push(state);
    let mut state = valid.clone();
    state.work = f64::INFINITY;
    invalid.push(state);
    for original in invalid {
        let mut state = original.clone();
        assert!(control.step(&mut state, 0.5, 1.0).is_err());
        same_state(&state, &original);
    }
}
