use super::*;
use crate::broadband_radiation::BroadbandRadiationAuthority;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget};
use fs_vfit::model::{PoleTerm, RationalModel};

// A supplied analytic receiver H(s)=2/(s+3)+1/4. Fitting error is outside
// these trajectory tests; the original continuous realization and exponential
// prepare the held-force part exactly as in the public observer constructor.
fn observer(rate: u32, omega: f64, zeta: f64, ceiling: f64, delay: f64) -> BroadbandModalObserver {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        rate,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 100.0,
            maximum_abs_pressure_pa: ceiling,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap();
    model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.02,
            velocity_m_sqrt_kg_per_s: 0.4,
        }])
        .unwrap();
    let fit = RationalModel {
        terms: vec![PoleTerm::Real {
            pole: -3.0,
            residue: 2.0,
        }],
        d: 0.25,
        e: 0.0,
    };
    let ss = fit.state_space();
    let n = ss.n + 3;
    let dt = model.sample_period_s();
    let q_scale = if omega == 0.0 { 1.0 } else { omega };
    let mut a = vec![0.0; n * n];
    a[1] = q_scale * dt;
    a[n] = -omega * omega / q_scale * dt;
    a[n + 1] = -2.0 * zeta * omega * dt;
    a[n + n - 1] = dt;
    for i in 0..ss.n {
        a[(i + 2) * n + 1] = ss.b[i] * dt;
        for j in 0..ss.n {
            a[(i + 2) * n + j + 2] = ss.a[i * ss.n + j] * dt;
        }
    }
    let filter = ExactFilter {
        n,
        transition: fs_phs::linear_exact::exponential(&a, n).unwrap(),
        state: vec![0.0; ss.n],
        candidate: vec![0.0; ss.n],
        c: ss.c,
        d: ss.d,
        q_scale,
        model: fit,
    };
    let workspace = ModalAcousticWorkspace::new(&model);
    BroadbandModalObserver {
        model,
        workspace,
        filters: vec![filter],
        constant: None,
        delay: Delay::new(delay, 1).unwrap(),
        report: BroadbandModalReport {
            source_id: "analytic one-pole receiver".into(),
            authority: BroadbandRadiationAuthority::EstimateOnly,
            producer_relative_error_estimate: 0.0,
            controls: BroadbandModalControls {
                fit_order: 1,
                fit_iterations: 1,
                delay_kernel_half_len: 1,
                max_held_out_rel_error: 1.0e-10,
            },
            sample_rate_hz: rate,
            held_out_max_rel_error: 0.0,
            delay_samples: delay,
            constant_narrow_band_compatibility: false,
            input_batches: 1,
        },
        linear_force_gains: None,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2.0e-11 * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn free_ramp_receiver_matches_independent_polynomial_convolution() {
    let mut observed = observer(8, 0.0, 0.0, 1.0e6, 0.0);
    observed.prepare_linear_forcing().unwrap();
    let dt = observed.model.sample_period_s();
    let mut q = 0.02;
    let mut v = 0.4;
    let mut z = 0.0;
    for i in 0..8 {
        let f0 = 1.2 - 0.1 * f64::from(i);
        let f1 = -0.4 + 0.2 * f64::from(i);
        let slope = (f1 - f0) / dt;
        let decay = (-3.0 * dt).exp();
        let i0 = (1.0 - decay) / 3.0;
        let i1 = dt / 3.0 - i0 / 3.0;
        let i2 = dt * dt / 3.0 - 2.0 * i1 / 3.0;
        z = decay * z + v * i0 + f0 * i1 + 0.5 * slope * i2;
        q += v * dt + 0.5 * f0 * dt * dt + slope * dt * dt * dt / 6.0;
        v += 0.5 * (f0 + f1) * dt;
        let pressure = observed.step_linear(&[f0], &[f1]).unwrap();
        close(pressure, 2.0 * z + 0.25 * v);
        close(observed.model.states()[0].displacement_m_sqrt_kg, q);
        close(observed.model.states()[0].velocity_m_sqrt_kg_per_s, v);
        assert_eq!(observed.last_modal_frame().viscous_dissipation_j, 0.0);
    }
}

#[test]
fn damped_receiver_work_and_heat_follow_the_same_continuous_ramp() {
    // Independently integrate physical ODEs, input power and viscous power.
    for zeta in [0.0, 0.2, 1.0, 3.0] {
        let omega = 4.0;
        let mut observed = observer(8, omega, zeta, 1.0e6, 0.0);
        observed.prepare_linear_forcing().unwrap();
        let dt = observed.model.sample_period_s();
        let mut x = [0.02, 0.4, 0.0, 0.0, 0.0];
        for (f0, f1) in [(1.2, -0.4), (-0.4, 2.0), (2.0, 0.7)] {
            x[3] = 0.0;
            x[4] = 0.0;
            let derivative = |x: [f64; 5], t: f64| {
                let f = f0 + (f1 - f0) * t / dt;
                let drag = 2.0 * zeta * omega;
                [
                    x[1],
                    f - drag * x[1] - omega * omega * x[0],
                    -3.0 * x[2] + x[1],
                    f * x[1],
                    drag * x[1] * x[1],
                ]
            };
            let add = |x: [f64; 5], a: [f64; 5], h: f64| std::array::from_fn(|i| x[i] + h * a[i]);
            let h = dt / 2048.0;
            for i in 0..2048 {
                let t = f64::from(i) * h;
                let a = derivative(x, t);
                let b = derivative(add(x, a, h / 2.0), t + h / 2.0);
                let c = derivative(add(x, b, h / 2.0), t + h / 2.0);
                let d = derivative(add(x, c, h), t + h);
                x = std::array::from_fn(|j| {
                    x[j] + h * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0
                });
            }
            close(
                observed.step_linear(&[f0], &[f1]).unwrap(),
                2.0 * x[2] + 0.25 * x[1],
            );
            close(observed.model.states()[0].displacement_m_sqrt_kg, x[0]);
            close(observed.model.states()[0].velocity_m_sqrt_kg_per_s, x[1]);
            close(observed.last_modal_frame().input_work_j, x[3]);
            close(observed.last_modal_frame().viscous_dissipation_j, x[4]);
        }
    }
}

#[test]
fn subdividing_the_same_force_ramp_preserves_receiver_and_mechanics() {
    let mut coarse = observer(8, 4.0, 0.6, 1.0e6, 0.0);
    let mut fine = observer(64, 4.0, 0.6, 1.0e6, 0.0);
    coarse.prepare_linear_forcing().unwrap();
    fine.prepare_linear_forcing().unwrap();
    let pressure = coarse.step_linear(&[1.0], &[-2.0]).unwrap();
    let mut refined_pressure = 0.0;
    let mut work = 0.0;
    let mut heat = 0.0;
    for i in 0..8 {
        refined_pressure = fine
            .step_linear(
                &[1.0 - 3.0 * f64::from(i) / 8.0],
                &[1.0 - 3.0 * f64::from(i + 1) / 8.0],
            )
            .unwrap();
        work += fine.last_modal_frame().input_work_j;
        heat += fine.last_modal_frame().viscous_dissipation_j;
    }
    close(pressure, refined_pressure);
    close(
        coarse.model.states()[0].displacement_m_sqrt_kg,
        fine.model.states()[0].displacement_m_sqrt_kg,
    );
    close(
        coarse.model.states()[0].velocity_m_sqrt_kg_per_s,
        fine.model.states()[0].velocity_m_sqrt_kg_per_s,
    );
    close(coarse.last_modal_frame().input_work_j, work);
    close(coarse.last_modal_frame().viscous_dissipation_j, heat);
}

#[test]
fn refusals_preserve_accepted_mechanics_receiver_and_delay() {
    let mut observed = observer(8, 4.0, 0.1, 0.2, 2.0);
    assert!(observed.step_linear(&[0.0], &[1.0]).is_err());
    observed.prepare_linear_forcing().unwrap();
    // Advance delay history before testing failed proposals.
    observed.step_linear(&[0.2], &[-0.1]).unwrap();
    let states = observed.model.states().to_vec();
    let filter = observed.filters[0].state.clone();
    let delay = observed.delay.history.clone();
    let cursor = observed.delay.cursor;
    let frame = observed.last_modal_frame().clone();
    // Delayed pressure has not reached its ceiling, but a huge proposed
    // mechanical displacement refuses after the receiver candidate is staged.
    for (start, end) in [(0.0, f64::NAN), (1.0e12, -1.0e12)] {
        assert!(observed.step_linear(&[start], &[end]).is_err());
        assert_eq!(observed.model.states(), states);
        assert_eq!(observed.filters[0].state, filter);
        assert_eq!(observed.delay.history, delay);
        assert_eq!(observed.delay.cursor, cursor);
        assert_eq!(observed.last_modal_frame(), &frame);
    }
    let mut replay = observer(8, 4.0, 0.1, 0.2, 2.0);
    replay.prepare_linear_forcing().unwrap();
    replay.step_linear(&[0.2], &[-0.1]).unwrap();
    assert_eq!(
        observed.step_linear(&[-0.1], &[0.0]).unwrap(),
        replay.step_linear(&[-0.1], &[0.0]).unwrap()
    );
    assert_eq!(observed.model.states(), replay.model.states());
    assert_eq!(observed.filters[0].state, replay.filters[0].state);

    let mut ceiling = observer(8, 4.0, 0.1, 1.0e-8, 0.0);
    ceiling.prepare_linear_forcing().unwrap();
    let initial = ceiling.model.states().to_vec();
    assert!(ceiling.step_linear(&[0.0], &[1.0]).is_err());
    assert_eq!(ceiling.model.states(), initial);
    assert_eq!(ceiling.filters[0].state, [0.0]);
    assert_eq!(ceiling.delay.history, [0.0]);
}

#[test]
fn equal_endpoints_retain_original_held_force_bits() {
    let mut held = observer(8, 4.0, 0.1, 1.0e6, 2.0);
    let mut ramp = observer(8, 4.0, 0.1, 1.0e6, 2.0);
    ramp.prepare_linear_forcing().unwrap();
    for force in [1.0, -0.5, 0.0, 0.2, 0.2] {
        assert_eq!(
            held.step(&[force]).unwrap().to_bits(),
            ramp.step_linear(&[force], &[force]).unwrap().to_bits()
        );
        assert_eq!(held.model.states(), ramp.model.states());
        assert_eq!(held.last_modal_frame(), ramp.last_modal_frame());
        assert_eq!(held.filters[0].state, ramp.filters[0].state);
    }
}

#[test]
fn real_feedthrough_supports_free_motion_and_complex_quadrature_refuses() {
    for omega in [0.0, 1.0e-200] {
        let mut observed = observer(8, omega, 0.0, 1.0e6, 0.0);
        observed.filters.clear();
        observed.constant = Some(vec![C64::from_re(2.5)]);
        observed.prepare_linear_forcing().unwrap();
        close(
            observed.step_linear(&[1.0], &[1.0]).unwrap(),
            2.5 * (0.4 + 0.125),
        );
        close(
            observed.step_linear(&[1.0], &[-1.0]).unwrap(),
            2.5 * (0.4 + 0.125),
        );
    }
    let mut narrow = observer(8, 4.0, 0.1, 1.0e6, 0.0);
    narrow.filters.clear();
    narrow.constant = Some(vec![C64::new(1.0, 0.1)]);
    let initial = narrow.model.states().to_vec();
    assert!(narrow.prepare_linear_forcing().is_err());
    assert!(narrow.step_linear(&[0.0], &[1.0]).is_err());
    assert_eq!(narrow.model.states(), initial);
}
