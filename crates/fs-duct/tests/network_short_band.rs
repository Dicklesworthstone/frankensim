//! Real/proper finite-window port interference on the authored rad/s clock.
use fs_duct::DuctError;
use fs_duct::network::{
    NetworkCoherentResponse,
    observation::{NetworkObservationLimits, NetworkObservationScope, PreparedNetworkObservation},
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::c64::C64;

fn line(omega: f64, pressure: C64, inflow: C64) -> NetworkCoherentResponse {
    NetworkCoherentResponse {
        omega,
        ports: vec!["receiver".into()],
        port_pressures: vec![pressure],
        port_inflows: vec![inflow],
        port_powers_w: vec![0.0],
        source_power_w: 0.0,
        external_load_power_w: 0.0,
        network_power_w: 0.0,
    }
}

fn limits() -> NetworkObservationLimits {
    NetworkObservationLimits {
        maximum_factors: 4,
        maximum_lines: 16,
        maximum_port_line_pairs: 64,
        maximum_port_frequency_pairs: 256,
        maximum_workspace_bytes: 100_000,
    }
}

fn prepare(
    mean: &[NetworkCoherentResponse],
    factors: &[&[NetworkCoherentResponse]],
) -> PreparedNetworkObservation {
    PreparedNetworkObservation::new(
        NetworkObservationScope {
            ports: &["receiver".into()],
            source_clock: "original/source-clock",
        },
        mean,
        factors,
        limits(),
    )
    .unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite()
            && actual.signum() == expected.signum()
            && (actual / expected - 1.0).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_real_short_window_keeps_weak_variance_mean_and_negative_boundary_work() {
    let lines = [1.0, -2.0, 1.0]
        .into_iter()
        .enumerate()
        .map(|(i, amplitude)| {
            line(
                (i + 1) as f64,
                C64::new(amplitude, 0.0),
                C64::new(-2.0 * amplitude, 0.0),
            )
        })
        .collect::<Vec<_>>();
    let owner = prepare(&lines, &[]);
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
        let sample = owner.observe_window(window).unwrap();
        let value = sample.port_values()[0];
        // p(t)=cos(t)-2cos(2t)+cos(3t)=-4sin(t/2)^2 cos(2t).
        // <p>=-T^2/12+O(T^4), <p^2>=T^4/80+O(T^6); U=-2p.
        let square = duration.powi(4) / 80.0;
        relative(value.mean_pressure_pa, -duration.powi(2) / 12.0, 3e-6);
        relative(value.mean_inflow_m3_s, duration.powi(2) / 6.0, 3e-6);
        relative(value.mean_square_pressure_pa2, square, 3e-6);
        relative(value.mean_square_inflow_m6_s2, 4.0 * square, 3e-6);
        relative(value.mean_boundary_power_w, -2.0 * square, 3e-6);
        assert!(value.mean_pressure_pa.powi(2) <= value.mean_square_pressure_pa2);
        assert_eq!(sample.boundary_power_w(), value.mean_boundary_power_w);
    }
}

#[test]
fn g1_proper_native_angular_carrier_keeps_half_power_and_source_independence() {
    let carrier = 2.0_f64.powi(49);
    let lines = [1.0, -2.0, 1.0]
        .into_iter()
        .enumerate()
        .map(|(i, amplitude)| {
            line(
                carrier + i as f64,
                C64::new(amplitude, 0.0),
                C64::new(-2.0 * amplitude, 0.0),
            )
        })
        .collect::<Vec<_>>();
    let shared = prepare(&[], &[&lines]);
    let independent = prepare(&[], &[&lines[0..1], &lines[1..2], &lines[2..3]]);
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
        let sample = shared.observe_window(window).unwrap();
        let value = sample.port_values()[0];
        let square = duration.powi(4) / 160.0;
        relative(value.mean_square_pressure_pa2, square, 3e-6);
        relative(value.mean_square_inflow_m6_s2, 4.0 * square, 3e-6);
        relative(value.mean_boundary_power_w, -2.0 * square, 3e-6);
        assert_eq!(value.mean_pressure_pa, 0.0);
        let incoherent = independent.observe_window(window).unwrap();
        relative(
            incoherent.port_values()[0].mean_square_pressure_pa2,
            3.0,
            2e-15,
        );
    }
}

#[test]
fn g1_real_quadrature_fields_keep_carrier_sum_terms_at_the_null() {
    let lines = [1.0, -2.0, 1.0]
        .into_iter()
        .enumerate()
        .map(|(i, amplitude)| {
            line(
                (i + 1) as f64,
                C64::new(0.0, amplitude),
                C64::new(0.0, -amplitude),
            )
        })
        .collect::<Vec<_>>();
    let owner = prepare(&lines, &[]);
    for duration in [1e-4_f64, 1e-7, 1e-10] {
        let sample = owner
            .observe_window(HarmonicWindow::new(-duration / 2.0, duration).unwrap())
            .unwrap();
        let value = sample.port_values()[0];
        // sin(t)-2sin(2t)+sin(3t)=-4sin(t/2)^2 sin(2t)=-2t^3+O(t^5).
        // The real square therefore starts at T^6/112, unlike analytic power.
        relative(
            value.mean_square_pressure_pa2,
            duration.powi(6) / 112.0,
            3e-6,
        );
        relative(value.mean_boundary_power_w, -duration.powi(6) / 112.0, 3e-6);
        assert_eq!(value.mean_pressure_pa, 0.0);
    }
}

#[test]
fn g0_g4_proper_half_range_and_cancelled_reuse_are_preserved() {
    let lines = [
        line(1.0, C64::new(1.5e154, 0.0), C64::ZERO),
        line(2.0, C64::ZERO, C64::ZERO),
        line(3.0, C64::ZERO, C64::ZERO),
    ];
    let owner = prepare(&[], &[&lines]);
    let window = HarmonicWindow::new(0.0, 0.01).unwrap();
    let baseline = owner.observe_window(window).unwrap();
    relative(
        baseline.port_values()[0].mean_square_pressure_pa2,
        1.125e308,
        2e-15,
    );
    let mut polls = 0;
    assert!(matches!(
        owner.observe_window_cancellable(window, || {
            polls += 1;
            polls == 6
        }),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(
        owner.observe_window(window).unwrap().port_values(),
        baseline.port_values()
    );
    assert_eq!(owner.port_frequency_pairs(), 9);
    assert!(
        PreparedNetworkObservation::new(
            NetworkObservationScope {
                ports: &["receiver".into()],
                source_clock: "original/source-clock",
            },
            &[],
            &[&lines],
            NetworkObservationLimits {
                maximum_port_frequency_pairs: 8,
                ..limits()
            }
        )
        .is_err()
    );
}
