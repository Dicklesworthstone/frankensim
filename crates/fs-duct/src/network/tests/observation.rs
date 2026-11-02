//! Physical port observations, including real-time and finite-ensemble oracles.
use super::*;
use super::super::observation::*;
use fs_fft::harmonic_window::HarmonicWindow;

fn limits() -> NetworkObservationLimits {
    NetworkObservationLimits {
        maximum_factors: 16,
        maximum_lines: 64,
        maximum_port_line_pairs: 4096,
        maximum_port_frequency_pairs: 1 << 20,
        maximum_workspace_bytes: 2 << 20,
    }
}

fn line(omega: f64, pressure: C64, inflow: C64) -> NetworkCoherentResponse {
    NetworkCoherentResponse {
        omega,
        ports: vec!["observed".into()],
        port_pressures: vec![pressure],
        port_inflows: vec![inflow],
        port_powers_w: vec![0.5 * (pressure * inflow.conj()).re],
        source_power_w: 0.0,
        external_load_power_w: 0.0,
        network_power_w: 0.0,
    }
}

fn prepare(
    mean: &[NetworkCoherentResponse],
    factors: &[&[NetworkCoherentResponse]],
) -> PreparedNetworkObservation {
    let ports = vec!["observed".into()];
    PreparedNetworkObservation::new(
        NetworkObservationScope { ports: &ports, source_clock: "fixture/common-clock" },
        mean, factors, limits(),
    ).unwrap()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2e-10 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}");
}

fn values(point: NetworkPortObservation) -> [f64; 5] {
    [point.mean_pressure_pa, point.mean_inflow_m3_s,
        point.mean_square_pressure_pa2, point.mean_square_inflow_m6_s2,
        point.mean_boundary_power_w]
}

// Independent REAL-signal Simpson integration. The four equally weighted
// phases have exactly E[xi]=E[xi^2]=0 and E[|xi|^2]=1. Every factor chooses its
// own phase, so the finite ensemble has the required joint second moments.
fn real_time_oracle(
    mean: &[NetworkCoherentResponse],
    factors: &[&[NetworkCoherentResponse]],
    window: HarmonicWindow,
    port: usize,
) -> [f64; 5] {
    const STEPS: usize = 20_000;
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let realizations = 4usize.pow(factors.len() as u32);
    let sample = |lines: &[NetworkCoherentResponse], phase: C64, t: f64| {
        let mut output = [0.0; 2];
        for line in lines {
            let (s, c) = (line.omega * t).sin_cos();
            for (sum, a) in output.iter_mut()
                .zip([line.port_pressures[port], line.port_inflows[port]])
            {
                let a = a * phase;
                *sum += a.re * c + a.im * s;
            }
        }
        output
    };
    let mut total = [0.0; 5];
    for realization in 0..realizations {
        for step in 0..=STEPS {
            let t = window.start_s() + window.duration_s() * step as f64 / STEPS as f64;
            let mut signal = sample(mean, C64::ONE, t);
            let mut code = realization;
            for factor in factors {
                let part = sample(factor, phases[code % 4], t);
                code /= 4;
                signal[0] += part[0];
                signal[1] += part[1];
            }
            let [p, u] = signal;
            let weight = if step == 0 || step == STEPS { 1.0 }
                else if step % 2 == 0 { 2.0 } else { 4.0 };
            for (sum, value) in total.iter_mut().zip([p, u, p*p, u*u, p*u]) {
                *sum += weight * value;
            }
        }
    }
    total.map(|value| value / (3 * STEPS * realizations) as f64)
}

#[test]
fn solved_matched_duct_lines_match_analytic_propagation_and_real_time_observation() {
    let gas = gas();
    let radius = 0.01;
    let length = 0.23;
    let resistance = gas.characteristic_impedance / (core::f64::consts::PI * radius * radius);
    let net = DuctNetwork {
        nodes: vec![port("drive"), port("load")],
        edges: vec![edge("tube", 0, 1, radius, length)],
    };
    let drive = [C64::new(1e-6, -0.3e-6), C64::new(-0.7e-6, 0.2e-6)];
    let mut solved = Vec::new();
    let mut analytic = Vec::new();
    for (omega, flow) in [1700.0, 1960.0].into_iter().zip(drive) {
        let response = network_drive_coherent(
            &net, &gas, omega, LossModel::Lossless,
            &[PortFlowSource { label: "drive".into(), volume_flow_m3_s: flow }],
            &[PortTermination { label: "load".into(),
                load: NetworkLoad::Impedance(C64::from_re(resistance)) }],
            None,
        ).unwrap();
        let phase = omega * length / gas.sound_speed;
        let transmitted = flow * C64::new(phase.cos(), phase.sin());
        let expected_p = [flow.scale(resistance), transmitted.scale(resistance)];
        let expected_u = [flow, -transmitted];
        for i in 0..2 {
            assert!((response.port_pressures[i] - expected_p[i]).abs() < 2e-10);
            assert!((response.port_inflows[i] - expected_u[i]).abs() < 2e-16);
        }
        let mut expected = response.clone();
        expected.port_pressures = expected_p.to_vec();
        expected.port_inflows = expected_u.to_vec();
        analytic.push(expected);
        solved.push(response);
    }
    let prepared = PreparedNetworkObservation::new(
        NetworkObservationScope { ports: &solved[0].ports, source_clock: "matched/source-clock" },
        &solved, &[], limits(),
    ).unwrap();
    let window = HarmonicWindow::new(0.0037, 0.009).unwrap();
    let actual = prepared.observe_window(window).unwrap();
    assert_eq!(actual.ports(), &["drive".to_owned(), "load".to_owned()]);
    assert_eq!(actual.source_clock(), "matched/source-clock");
    for port in 0..2 {
        let expected = real_time_oracle(&analytic, &[], window, port);
        for ((value, reference), scale) in values(actual.port_values()[port]).into_iter()
            .zip(expected).zip([5.0, 3e-6, 10.0, 5e-12, 1e-5])
        {
            near(value, reference, scale);
        }
    }
    // Work at the two ends includes the actual delayed envelope and can differ
    // over a finite window even though the matched line is lossless.
    near(actual.boundary_power_w(),
        actual.port_values().iter().map(|value| value.mean_boundary_power_w).sum(), 1e-5);
}

#[test]
fn deterministic_mean_plus_proper_factors_matches_four_phase_real_time_ensemble() {
    let mean = [line(2.0, C64::new(0.2, 0.3), C64::new(-0.1, 0.2))];
    let first = [
        line(3.0, C64::new(1.0, 0.2), C64::new(0.4, -0.1)),
        line(5.0, C64::new(-0.7, 0.1), C64::new(0.2, 0.5)),
    ];
    let second = [line(3.0, C64::new(0.3, -0.8), C64::new(0.1, 0.2))];
    let factors = [&first[..], &second[..]];
    let prepared = prepare(&mean, &factors);
    let window = HarmonicWindow::new(0.17, 0.83).unwrap();
    let actual = prepared.observe_window(window).unwrap();
    let expected = real_time_oracle(&mean, &factors, window, 0);
    for (a, b) in values(actual.port_values()[0]).into_iter().zip(expected) {
        near(a, b, 4.0);
    }
}

#[test]
fn tiny_destructive_beats_retain_the_midpoint_and_width_contributions() {
    let delta = 2.0_f64.powi(-20);
    let duration = delta;
    let lines = [
        line(1.0, C64::ONE, C64::ONE),
        line(1.0 + delta, -C64::ONE, -C64::ONE),
    ];
    let prepared = prepare(&[], &[&lines]);
    let result = prepared.observe_window(HarmonicWindow::new(0.0, duration).unwrap()).unwrap();
    let expected = (delta * duration).powi(2) / 6.0;
    let point = result.port_values()[0];
    for actual in [point.mean_square_pressure_pa2, point.mean_square_inflow_m6_s2,
        point.mean_boundary_power_w]
    {
        assert!((actual / expected - 1.0).abs() < 2e-12,
            "weak beat lost: {actual:e}, expected {expected:e}");
    }
    assert_eq!(point.mean_pressure_pa, 0.0);
    assert_eq!(point.mean_inflow_m3_s, 0.0);
}

#[test]
fn grouped_low_components_survive_later_frequency_pair_cancellation() {
    let large = 2.0_f64.powi(54);
    let lines = [
        line(2.0, C64::from_re(large), C64::ONE),
        line(2.0, C64::ONE, C64::ZERO),
        line(3.0, C64::from_re(-large), C64::ONE),
    ];
    let prepared = prepare(&[], &[&lines]);
    assert_eq!(prepared.port_frequency_pairs(), 4);
    let result = prepared.observe_window(HarmonicWindow::new(0.0, 1.0).unwrap()).unwrap();
    // The grouped pressure large+1 cannot be rounded to C64 first. The large
    // diagonal and cross-frequency pressure/flow terms cancel, leaving this.
    let expected = 0.5 * (1.0 + 1.0_f64.sin());
    near(result.port_values()[0].mean_boundary_power_w, expected, 1.0);
}

#[test]
fn factors_are_not_added_as_fields_and_empty_input_is_dark() {
    let a = [line(7.0, C64::ONE, C64::from_re(-2.0))];
    let b = [line(7.0, -C64::ONE, C64::from_re(2.0))];
    let combined = [a[0].clone(), b[0].clone()];
    let coherent = prepare(&[], &[&combined]);
    let independent = prepare(&[], &[&a, &b]);
    let dark = prepare(&[], &[]);
    let window = HarmonicWindow::new(0.13, 0.0001).unwrap();
    assert_eq!(values(coherent.observe_window(window).unwrap().port_values()[0]), [0.0; 5]);
    let value = independent.observe_window(window).unwrap().port_values()[0];
    near(value.mean_square_pressure_pa2, 1.0, 1.0);
    near(value.mean_square_inflow_m6_s2, 4.0, 4.0);
    near(value.mean_boundary_power_w, -2.0, 2.0);
    assert_eq!(values(dark.observe_window(window).unwrap().port_values()[0]), [0.0; 5]);
}

#[test]
fn finite_reactive_port_work_keeps_its_sign_and_is_not_heat() {
    let lines = [line(1.0, C64::ONE, C64::new(0.0, 1.0))];
    let prepared = prepare(&lines, &[]);
    let window = HarmonicWindow::new(core::f64::consts::FRAC_PI_2,
        core::f64::consts::FRAC_PI_2).unwrap();
    let result = prepared.observe_window(window).unwrap();
    near(result.port_values()[0].mean_boundary_power_w, -1.0 / core::f64::consts::PI, 1.0);
    assert!(result.boundary_power_w() < 0.0);
    assert_eq!(lines[0].port_powers_w[0], 0.0);
}

#[test]
fn authored_clock_shift_preserves_all_real_observables() {
    let mean = [line(2.3, C64::new(0.2, -0.1), C64::new(0.4, 0.2))];
    let factor = [
        line(3.7, C64::new(0.8, 0.3), C64::new(0.2, -0.1)),
        line(5.2, C64::new(-0.4, 0.5), C64::new(0.1, 0.6)),
    ];
    let shift = 0.27;
    let translate = |lines: &[NetworkCoherentResponse]| {
        lines.iter().map(|line| {
            let mut result = line.clone();
            let angle = line.omega * shift;
            let rotation = C64::new(angle.cos(), -angle.sin());
            for z in result.port_pressures.iter_mut().chain(&mut result.port_inflows) {
                *z = *z * rotation;
            }
            result
        }).collect::<Vec<_>>()
    };
    let moved_mean = translate(&mean);
    let moved_factor = translate(&factor);
    let original = prepare(&mean, &[&factor]);
    let moved = PreparedNetworkObservation::new(
        NetworkObservationScope {
            ports: &moved_mean[0].ports, source_clock: "fixture/shifted-clock",
        },
        &moved_mean, &[&moved_factor], limits(),
    ).unwrap();
    let a = original.observe_window(HarmonicWindow::new(0.41, 0.73).unwrap()).unwrap();
    let b = moved.observe_window(HarmonicWindow::new(0.41 - shift, 0.73).unwrap()).unwrap();
    for (x, y) in values(a.port_values()[0]).into_iter().zip(values(b.port_values()[0])) {
        near(x, y, 5.0);
    }
}

#[test]
fn scope_budgets_invalid_fields_and_cancellation_publish_no_partial_result() {
    let lines = [line(2.0, C64::ONE, C64::ONE), line(3.0, C64::ONE, -C64::ONE)];
    let ports = lines[0].ports.clone();
    let scope = NetworkObservationScope { ports: &ports, source_clock: "authored/clock" };
    let prepared = PreparedNetworkObservation::new(scope, &lines, &[], limits()).unwrap();
    let exact = NetworkObservationLimits {
        maximum_lines: 2, maximum_port_line_pairs: 2, maximum_port_frequency_pairs: 4,
        maximum_workspace_bytes: prepared.workspace_bytes(), ..limits()
    };
    assert!(PreparedNetworkObservation::new(scope, &lines, &[], exact).is_ok());
    for bad_limits in [
        NetworkObservationLimits { maximum_lines: 1, ..exact },
        NetworkObservationLimits { maximum_port_line_pairs: 1, ..exact },
        NetworkObservationLimits { maximum_port_frequency_pairs: 3, ..exact },
        NetworkObservationLimits { maximum_workspace_bytes: exact.maximum_workspace_bytes - 1, ..exact },
    ] {
        assert!(PreparedNetworkObservation::new(scope, &lines, &[], bad_limits).is_err());
    }
    assert!(PreparedNetworkObservation::new(scope, &[], &[&lines], limits()).is_ok());
    assert!(PreparedNetworkObservation::new(scope, &[], &[&lines],
        NetworkObservationLimits { maximum_factors: 0, ..limits() }).is_err());
    assert!(PreparedNetworkObservation::new(
        NetworkObservationScope { ports: &ports, source_clock: " " },
        &lines, &[], limits()).is_err());
    let mut wrong = lines.clone();
    wrong[0].ports[0] = "different".into();
    assert!(PreparedNetworkObservation::new(scope, &wrong, &[], limits()).is_err());
    wrong = lines.clone();
    wrong[0].port_pressures[0].re = f64::NAN;
    assert!(PreparedNetworkObservation::new(scope, &wrong, &[], limits()).is_err());
    wrong = lines.clone();
    wrong[0].omega = 0.0;
    assert!(PreparedNetworkObservation::new(scope, &wrong, &[], limits()).is_err());
    assert!(matches!(PreparedNetworkObservation::new_cancellable(
        scope, &lines, &[], limits(), || true), Err(DuctError::Cancelled { .. })));
    let window = HarmonicWindow::new(0.1, 0.7).unwrap();
    let before = prepared.observe_window(window).unwrap();
    assert!(prepared.observe_window(HarmonicWindow::new(1e308, 1.0).unwrap()).is_err());
    let mut calls = 0;
    assert!(matches!(prepared.observe_window_cancellable(window, || {
        calls += 1; calls > 4
    }), Err(DuctError::Cancelled { .. })));
    let after = prepared.observe_window(window).unwrap();
    assert_eq!(before.port_values(), after.port_values());
}
