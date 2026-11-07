//! G1/G3 analytic polarization, physical energy, and exact interval composition.
use super::*;

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

#[test]
fn g1_static_resonance_phase_and_high_frequency_limit() {
    let pole = LorentzPole::try_new(3e-11, 2e15, 4e12).unwrap();
    assert_eq!(pole.susceptibility_f_m(0.0).unwrap(), C64::from_re(3e-11));
    let resonant = pole.susceptibility_f_m(2e15).unwrap();
    assert_eq!(resonant.re, 0.0);
    near(resonant.im, -1.5e-8, 1.5e-8);
    for omega in [1e14, 1.8e15, 2.01e15, 3e15, 1e18] {
        let x = omega / 2e15;
        let denominator = C64::new((1.0 - x) * (1.0 + x), 0.002 * x);
        let reference = C64::from_re(3e-11) / denominator;
        let actual = pole.susceptibility_f_m(omega).unwrap();
        near(actual.re, reference.re, reference.abs());
        near(actual.im, reference.im, reference.abs());
        assert!(actual.im < 0.0);
    }
    assert!(pole.susceptibility_f_m(3e15).unwrap().re < 0.0);
}

#[test]
fn g1_high_q_and_lossless_poles_are_not_regularized() {
    let pole = LorentzPole::try_new(2.0, 1e100, 1e84).unwrap();
    let response = pole.susceptibility_f_m(1e100).unwrap();
    near(response.im, -2e16, 2e16);
    assert_eq!(response.re, 0.0);
    let lossless = LorentzPole::try_new(2.0, 4.0, 0.0).unwrap();
    assert_eq!(
        lossless.susceptibility_f_m(4.0),
        Err(LorentzError::SingularResonance)
    );
    assert_eq!(lossless.susceptibility_f_m(5.0).unwrap().im, 0.0);
    let left = pole.susceptibility_f_m(1e100_f64.next_down()).unwrap();
    let right = pole.susceptibility_f_m(1e100_f64.next_up()).unwrap();
    assert!(left.re > 0.0 && right.re < 0.0);
}

#[test]
fn g1_step_field_matches_damped_closed_form_and_physical_heat() {
    let delta = 2e-6;
    let omega = 100.0;
    let gamma = 14.0;
    let field = 300.0;
    let duration = 0.037;
    let pole = LorentzPole::try_new(delta, omega, gamma).unwrap();
    let frame = pole
        .prepare_step(duration)
        .unwrap()
        .advance(LorentzState::default(), [field], [field])
        .unwrap();
    let a = gamma * 0.5;
    let wd = det::sqrt(omega * omega - a * a);
    let envelope = det::exp(-a * duration);
    let p = delta
        * field
        * (1.0 - envelope * (det::cos(wd * duration) + (a / wd) * det::sin(wd * duration)));
    let j = delta * field * omega * omega / wd * envelope * det::sin(wd * duration);
    near(frame.state.polarization_c_m2[0], p, delta * field);
    near(frame.state.current_a_m2[0], j, delta * field * omega);
    near(frame.supplied_work_j_m3, field * p, (field * p).abs());
    assert!(frame.heat_j_m3 > 0.0);
    near(
        frame.energy_residual_j_m3,
        0.0,
        frame.supplied_work_j_m3.abs(),
    );
}

#[test]
fn g3_linear_field_subdivision_retains_state_work_and_heat() {
    for gamma in [0.0, 0.01, 20.0, 35.0] {
        let pole = LorentzPole::try_new(4e-6, 10.0, gamma).unwrap();
        let initial = LorentzState {
            polarization_c_m2: [1e-3, -2e-4, 7e-4],
            current_a_m2: [-2e-3, 3e-3, 1e-3],
        };
        let start = [400.0, -20.0, 80.0];
        let end = [-20.0, 800.0, -500.0];
        let middle = std::array::from_fn(|i| 0.5 * start[i] + 0.5 * end[i]);
        let whole = pole
            .prepare_step(0.2)
            .unwrap()
            .advance(initial, start, end)
            .unwrap();
        let half = pole.prepare_step(0.1).unwrap();
        let a = half.advance(initial, start, middle).unwrap();
        let b = half.advance(a.state, middle, end).unwrap();
        for i in 0..3 {
            near(
                whole.state.polarization_c_m2[i],
                b.state.polarization_c_m2[i],
                1e-2,
            );
            near(whole.state.current_a_m2[i], b.state.current_a_m2[i], 0.1);
        }
        near(
            whole.supplied_work_j_m3,
            a.supplied_work_j_m3 + b.supplied_work_j_m3,
            10.0,
        );
        near(whole.heat_j_m3, a.heat_j_m3 + b.heat_j_m3, 10.0);
        near(whole.energy_residual_j_m3, 0.0, 10.0);
        assert!(whole.heat_j_m3 >= 0.0);
        if gamma == 0.0 {
            assert_eq!(whole.heat_j_m3, 0.0);
        }
    }
}

#[test]
fn g1_optical_rate_ringdown_has_inertial_storage() {
    let omega = 3e15;
    let pole = LorentzPole::try_new(2e-11, omega, 0.0).unwrap();
    let initial = LorentzState {
        polarization_c_m2: [1e-8],
        current_a_m2: [0.0],
    };
    let frame = pole
        .prepare_step(core::f64::consts::FRAC_PI_2 / omega)
        .unwrap()
        .advance(initial, [0.0], [0.0])
        .unwrap();
    near(frame.state.polarization_c_m2[0], 0.0, 1e-8);
    near(frame.state.current_a_m2[0], -omega * 1e-8, omega * 1e-8);
    near(
        frame.stored_after_j_m3,
        frame.stored_before_j_m3,
        frame.stored_before_j_m3,
    );
    assert_eq!(frame.heat_j_m3, 0.0);
    assert_eq!(frame.supplied_work_j_m3, 0.0);
}

#[test]
fn g0_invalid_sources_states_and_fields_refuse() {
    for (delta, omega, gamma) in [
        (0.0, 1.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, -1.0),
        (f64::NAN, 1.0, 1.0),
        (1.0, f64::INFINITY, 1.0),
    ] {
        assert!(LorentzPole::try_new(delta, omega, gamma).is_err());
    }
    let pole = LorentzPole::try_new(2.0, 10.0, 1.0).unwrap();
    assert!(pole.prepare_step(0.0).is_err());
    assert!(pole.susceptibility_f_m(-1.0).is_err());
    let initial = LorentzState::<2>::default();
    assert!(
        pole.prepare_step(0.1)
            .unwrap()
            .advance(initial, [0.0, 0.0], [0.0, f64::NAN])
            .is_err()
    );
    assert_eq!(initial, LorentzState::default());
    assert!(
        pole.stored_energy_j_m3(LorentzState::<0>::default())
            .is_err()
    );
}
