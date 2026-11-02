//! G1 affine forcing against an independent ODE; G3 complete source invariance.
use super::super::{ModalColumnConnection, ModalCouplingConfig};
use super::tests::{all_states, fixture};
use super::*;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use broadband::{BroadbandCoupledModalTones, BroadbandModalControls, ModalTransferSamples};
use fs_math::c64::C64;
use fs_ports::PortedModule;

fn budget() -> CoupledToneBudget {
    CoupledToneBudget {
        max_dimension: 256,
        max_dense_bytes: 128 * 1_048_576,
        max_multiply_terms: 100_000_000_000,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-8 * a.abs().max(b.abs()) + 2e-13,
        "{a:e} != {b:e}"
    );
}
fn simple(rate: u32) -> CoupledModalSystem {
    let model = ModalAcousticTimeModel::try_new(
        rate,
        [(500.0, 0.01), (900.0, 0.04)]
            .into_iter()
            .map(|(w, zeta)| ModalAcousticMode {
                angular_frequency_rad_s: w,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect(),
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    CoupledModalSystem::new_with_columns(
        vec![model],
        vec![ModalColumnConnection {
            column: vec![1.0, -0.3],
            stiffness_n_m: 4e4,
            damping_n_s_m: 20.0,
            rest_extension_m: 0.0,
        }],
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 1,
            max_setup_terms: 1_000_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e12,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e8,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn receiver(source: CoupledModalTones, rate: u32) -> BroadbandCoupledModalTones {
    let transfer = |mode: usize, w: f64| {
        C64::from_re([0.25, -0.4][mode])
            + C64::from_re([1.5, -0.6][mode] * 900.0) / C64::new(900.0, -w)
    };
    let omega = vec![100.0, 200.0, 400.0, 650.0, 900.0, 1400.0, 2800.0];
    let samples = ModalTransferSamples {
        source_id: "analytic signed affine-source receiver".into(),
        per_mode: (0..2)
            .map(|m| omega.iter().map(|&w| transfer(m, w)).collect())
            .collect(),
        omega_rad_s: omega,
        held_out: [175.0, 750.0, 1800.0]
            .into_iter()
            .map(|w| (w, (0..2).map(|m| transfer(m, w)).collect()))
            .collect(),
        propagation_delay_s: 0.0,
        producer_relative_error_estimate: 0.0,
    };
    BroadbandCoupledModalTones::build(
        source,
        &samples,
        BroadbandModalControls {
            fit_order: 1,
            fit_iterations: 6,
            delay_kernel_half_len: 1,
            max_held_out_rel_error: 1e-8,
        },
        rate,
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
    .0
}
fn rotate(source: &CoupledModalTones, cosine: &mut [f64], sine: &mut [f64]) {
    let tones = source.tone_count();
    for (i, (c, s)) in cosine.iter_mut().zip(sine).enumerate() {
        (*c, *s) = source.advance_force_phase(i % tones, *c, *s).unwrap();
    }
}

// Independently written physical oscillator/link/filter ODE and its power
// integrals. It calls no prepared flow, heat, force or observer implementation.
fn derivative(
    x: &[f64; 9],
    t: f64,
    h: f64,
    c0: &[f64],
    s0: &[f64],
    c1: &[f64],
    s1: &[f64],
) -> [f64; 9] {
    let fraction = t / h;
    let mut force = [0.0; 2];
    for mode in 0..2 {
        for (tone, omega) in [700.0, 701.3].into_iter().enumerate() {
            let i = 2 * mode + tone;
            force[mode] += (c0[i] * (1.0 - fraction) + c1[i] * fraction) * (omega * t).cos()
                + (s0[i] * (1.0 - fraction) + s1[i] * fraction) * (omega * t).sin();
        }
    }
    let extension = x[0] - 0.3 * x[2];
    let speed = x[1] - 0.3 * x[3];
    let effort = 4e4 * extension + 20.0 * speed;
    [
        x[1],
        -500.0_f64.powi(2) * x[0] - 10.0 * x[1] - effort + force[0],
        x[3],
        -900.0_f64.powi(2) * x[2] - 72.0 * x[3] + 0.3 * effort + force[1],
        -900.0 * x[4] + x[1],
        -900.0 * x[5] + x[3],
        force[0] * x[1] + force[1] * x[3],
        10.0 * x[1] * x[1] + 72.0 * x[3] * x[3],
        20.0 * speed * speed,
    ]
}
fn reference_step(x: &mut [f64; 9], h: f64, c0: &[f64], s0: &[f64], c1: &[f64], s1: &[f64]) {
    let substep = h / 512.0;
    for j in 0..512 {
        let t = j as f64 * substep;
        let a = derivative(x, t, h, c0, s0, c1, s1);
        let b = derivative(
            &std::array::from_fn(|i| x[i] + 0.5 * substep * a[i]),
            t + 0.5 * substep,
            h,
            c0,
            s0,
            c1,
            s1,
        );
        let c = derivative(
            &std::array::from_fn(|i| x[i] + 0.5 * substep * b[i]),
            t + 0.5 * substep,
            h,
            c0,
            s0,
            c1,
            s1,
        );
        let d = derivative(
            &std::array::from_fn(|i| x[i] + substep * c[i]),
            t + substep,
            h,
            c0,
            s0,
            c1,
            s1,
        );
        for i in 0..9 {
            x[i] += substep * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
}

#[test]
fn affine_attack_release_matches_independent_coupled_filter_and_power_ode() {
    let source = CoupledModalTones::new_affine(
        simple(4000),
        &[700.0, 701.3],
        budget(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let mut observed = receiver(source, 4000);
    let (mut cosine, mut sine) = ([0.4, -0.39, -0.2, 0.18], [-0.2, 0.21, 0.3, -0.29]);
    let mut exact = [0.0; 9];
    let gain = |sample: usize| {
        if sample <= 10 {
            sample as f64 / 10.0
        } else {
            (20.0 - sample as f64).max(0.0) / 10.0
        }
    };
    for sample in 0..26 {
        let c0 = cosine.map(|x| gain(sample) * x);
        let s0 = sine.map(|x| gain(sample) * x);
        let c1 = cosine.map(|x| gain(sample + 1) * x);
        let s1 = sine.map(|x| gain(sample + 1) * x);
        let before = exact;
        reference_step(&mut exact, 1.0 / 4000.0, &c0, &s0, &c1, &s1);
        let pressure = observed.step_affine(&c0, &s0, &c1, &s1).unwrap();
        for (actual, expected) in all_states(observed.network()).iter().zip(&exact[..4]) {
            close(*actual, *expected);
        }
        close(
            pressure,
            0.25 * exact[1] - 0.4 * exact[3] + 1350.0 * exact[4] - 540.0 * exact[5],
        );
        let frame = observed.network().last_frame().unwrap();
        close(frame.external_work_j, exact[6] - before[6]);
        close(frame.component_dissipation_j, exact[7] - before[7]);
        close(frame.connection_dissipation_j, exact[8] - before[8]);
        assert_eq!(frame.numerical_dissipation_j, 0.0);
        if sample > 20 {
            assert!(
                pressure.abs() > 1e-9,
                "release must retain physical/filter ringdown"
            );
        }
        rotate(observed.source(), &mut cosine, &mut sine);
    }
}

fn accounts(frame: &super::super::CoupledModalFrame) -> [f64; 5] {
    [
        frame.external_work_j,
        frame.component_dissipation_j,
        frame.connection_dissipation_j,
        frame.material_dissipation_j,
        frame.radiation_dissipation_j,
    ]
}

#[test]
fn affine_coupled_memory_and_observation_share_source_accounts_and_subdivision() {
    let gate = CancelGate::new_clock_free();
    let frequencies = [700.0, 700.3, 0.0];
    let make = |rate| {
        CoupledModalTones::new_affine(fixture(rate, true, 1e-5), &frequencies, budget(), &gate)
            .unwrap()
    };
    let mut physical = make(4000);
    let mut full = receiver(make(4000), 4000);
    let mut half = receiver(make(8000), 8000);
    let mut cosine = [0.4, -0.39, 0.01, -0.2, 0.18, -0.03];
    let mut sine = [-0.2, 0.21, 0.0, 0.3, -0.29, 0.0];
    for sample in 0..8 {
        let scale = if sample % 2 == 0 { 0.7 } else { 1.2 };
        let end_c = cosine.map(|x| scale * x);
        let end_s = sine.map(|x| scale * x);
        let mut mid_c = std::array::from_fn::<_, 6, _>(|i| 0.5 * (cosine[i] + end_c[i]));
        let mut mid_s = std::array::from_fn::<_, 6, _>(|i| 0.5 * (sine[i] + end_s[i]));
        let original = physical
            .step_affine(&cosine, &sine, &end_c, &end_s)
            .unwrap()
            .clone();
        let pressure = full.step_affine(&cosine, &sine, &end_c, &end_s).unwrap();
        let complete = full.network().last_frame().unwrap().clone();
        assert_eq!(all_states(physical.network()), all_states(full.network()));
        assert_eq!(accounts(&original), accounts(&complete));
        assert_eq!(
            physical.network().last_window_transfers().unwrap(),
            full.network().last_window_transfers().unwrap()
        );
        half.step_affine(&cosine, &sine, &mid_c, &mid_s).unwrap();
        let first = half.network().last_frame().unwrap().clone();
        rotate(half.source(), &mut mid_c, &mut mid_s);
        let (mut finish_c, mut finish_s) = (end_c, end_s);
        rotate(half.source(), &mut finish_c, &mut finish_s);
        let divided = half
            .step_affine(&mid_c, &mid_s, &finish_c, &finish_s)
            .unwrap();
        close(pressure, divided);
        for (a, b) in all_states(full.network())
            .iter()
            .zip(all_states(half.network()))
        {
            close(*a, b);
        }
        for ((a, b), c) in accounts(&complete)
            .into_iter()
            .zip(accounts(&first))
            .zip(accounts(half.network().last_frame().unwrap()))
        {
            close(a, b + c);
        }
        assert!(complete.material_dissipation_j > 0.0 && complete.radiation_dissipation_j > 0.0);
        cosine = end_c;
        sine = end_s;
        rotate(full.source(), &mut cosine, &mut sine);
    }
}

#[test]
fn affine_receiver_refusal_preserves_history_and_requires_explicit_preparation() {
    let gate = CancelGate::new_clock_free();
    let make = || {
        let mut network = simple(4000);
        network.config.maximum_abs_pressure_pa = 0.01;
        receiver(
            CoupledModalTones::new_affine(network, &[700.0], budget(), &gate).unwrap(),
            4000,
        )
    };
    let mut actual = make();
    let mut control = make();
    let (c0, s0, c1, s1) = ([1e-6, 0.0], [0.0; 2], [2e-6, 0.0], [0.0; 2]);
    actual.step_affine(&c0, &s0, &c1, &s1).unwrap();
    control.step_affine(&c0, &s0, &c1, &s1).unwrap();
    let state = all_states(actual.network());
    let frame = actual.network().last_frame().unwrap().clone();
    let ports = actual.network().last_window_transfers().unwrap().to_vec();
    let pressure = actual.last_pressure_pa();
    assert!(actual.step_affine(&c0, &s0, &[1e5, 0.0], &s1).is_err());
    assert!(actual.step_affine(&c0, &s0, &[f64::NAN, 0.0], &s1).is_err());
    assert!(actual.step(&c0, &s0).is_err());
    assert_eq!(all_states(actual.network()), state);
    assert_eq!(actual.network().last_frame().unwrap(), &frame);
    assert_eq!(actual.network().last_window_transfers().unwrap(), ports);
    assert_eq!(actual.last_pressure_pa(), pressure);
    assert_eq!(
        actual.step_affine(&c0, &s0, &c1, &s1).unwrap(),
        control.step_affine(&c0, &s0, &c1, &s1).unwrap()
    );
    let mut constant = CoupledModalTones::new(simple(4000), &[700.0], budget(), &gate).unwrap();
    assert!(constant.step_affine(&c0, &s0, &c1, &s1).is_err());
    let cap = CoupledToneBudget {
        max_dimension: 12,
        ..budget()
    };
    assert!(CoupledModalTones::new(simple(4000), &[700.0], cap, &gate).is_ok());
    assert!(CoupledModalTones::new_affine(simple(4000), &[700.0], cap, &gate).is_err());
}
