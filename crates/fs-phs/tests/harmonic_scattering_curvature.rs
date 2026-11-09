//! Original loaded equations and power curvature, without dummy parameter solves.
use fs_math::c64::C64;
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringModel};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, HarmonicModel, HarmonicParameterDirection};
use fs_phs::{PortHamiltonian, QuadraticStorage};

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 16,
        max_dense_bytes: 4_000_000,
        max_dense_work: 8_000_000,
        relative_residual_tolerance: 2e-11,
    }
}
fn wave_budget() -> ScatteringBudget {
    ScatteringBudget {
        maximum_ports: 4,
        maximum_payload_bytes: 1_000_000,
        maximum_work: 2_000_000,
    }
}
fn oscillator(c: f64, wave: ScatteringBudget) -> ScatteringModel {
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., c],
        vec![0., 0.9],
        Box::new(QuadraticStorage::new(vec![4.2, 0., 0., 1. / 1.3], 2).unwrap()),
    )
    .unwrap();
    ScatteringModel::new(
        HarmonicModel::from_ode(&system, budget(), || false).unwrap(),
        &[(0, 2.4)],
        wave,
        || false,
    )
    .unwrap()
}
fn source(w: f64) -> [C64; 3] {
    let phase = C64::new((0.37 * w).cos(), (0.37 * w).sin());
    [
        phase.scale(1. + 0.1 * w),
        phase * C64::new(0.1, 0.37 * (1. + 0.1 * w)),
        phase * C64::new(-0.37 * 0.37 * (1. + 0.1 * w), 2. * 0.37 * 0.1),
    ]
}
fn close(a: C64, b: C64) {
    assert!((a - b).abs() < 3e-10 * b.abs().max(1e-5), "{a:?} != {b:?}");
}
fn real(a: f64, b: f64) {
    // The independent conjugate-work formula cancels O(1) terms at a
    // lossless bare resonance; its near-zero remainder needs an absolute bound.
    assert!((a - b).abs() < 3e-10 * b.abs().max(0.1), "{a:e} != {b:e}");
}
fn norm_curvature(v: [C64; 3]) -> f64 {
    2. * (v[1].norm_sq() + (v[0].conj() * v[2]).re)
}

#[test]
fn three_solve_curvature_matches_original_loaded_oscillator_and_wave_work() {
    for c in [0., 0.17] {
        let model = oscillator(c, wave_budget());
        for w in [0.03, 1., (4.2_f64 / 1.3).sqrt(), 3.1] {
            let a = source(w);
            let d = C64::new(4.2 - 1.3 * w * w, -w * (c + 0.9 * 0.9 / 2.4));
            let dw = C64::new(-2. * 1.3 * w, -(c + 0.9 * 0.9 / 2.4));
            let q = a[0].scale(2. * 0.9 / 2.4) / d;
            let qw = (a[1].scale(2. * 0.9 / 2.4) - dw * q) / d;
            let qww = (a[2].scale(2. * 0.9 / 2.4) - dw * qw.scale(2.) + q.scale(2. * 1.3)) / d;
            let s = C64::new(0., -w);
            let v = [
                s * q,
                s * qw + C64::new(0., -1.) * q,
                s * qww + C64::new(0., -2.) * qw,
            ];
            let y = v.map(|v| v.scale(0.9));
            let u = [0, 1, 2].map(|i| (a[i].scale(2.) - y[i]).scale(1. / 2.4));
            let b = [0, 1, 2].map(|i| y[i] - a[i]);
            let r = model
                .factor(w, || false)
                .unwrap()
                .response_with_frequency_curvature(&[a[0]], Some(&[a[1]]), Some(&[a[2]]), || false)
                .unwrap();
            close(r.physical.state_second_derivative[0], qww);
            close(r.physical.state_second_derivative[1], v[2].scale(1.3));
            close(r.physical.output_second_derivative[0], y[2]);
            close(r.physical_inputs[0].second_frequency, u[2]);
            close(r.waves[0].outgoing.second_frequency, b[2]);
            real(
                r.physical.dissipated_power_second_derivative,
                0.5 * c * norm_curvature(v),
            );
            real(
                r.waves[0].incident_power.second_frequency,
                norm_curvature(a) / (2. * 2.4),
            );
            real(
                r.waves[0].outgoing_power.second_frequency,
                norm_curvature(b) / (2. * 2.4),
            );
            let work = 0.5
                * ((u[2].conj() * y[0]).re
                    + 2. * (u[1].conj() * y[1]).re
                    + (u[0].conj() * y[2]).re);
            real(r.physical.port_power_second_derivative[0], work);
            assert!(r.waves[0].power_residual.second_frequency.abs() < 3e-9);
            assert!(r.wave_power_residual.second_frequency.abs() < 3e-9);
            assert!(r.physical.relative_second_derivative_residual < 2e-11);
        }
    }
}

#[test]
fn curvature_preserves_small_coherent_source_remainder_and_matches_shared_mixed_path() {
    let system = PortHamiltonian::new(
        1,
        4,
        vec![0.],
        vec![1.],
        vec![1., 1., 1., 0.],
        Box::new(QuadraticStorage::new(vec![1.], 1).unwrap()),
    )
    .unwrap();
    let model = ScatteringModel::new(
        HarmonicModel::from_ode(&system, budget(), || false).unwrap(),
        &[(3, 1.)],
        wave_budget(),
        || false,
    )
    .unwrap();
    let factor = model.factor(1., || false).unwrap();
    let input = [C64::ZERO, C64::ONE, C64::ZERO, C64::ZERO];
    let large = C64::from_re(2.0_f64.powi(54));
    let second = [large, C64::ONE, -large, C64::ZERO];
    let r = factor
        .response_with_frequency_curvature(&input, None, Some(&second), || false)
        .unwrap();
    // x'' = 1/(1-i) - 2/(1-i)^3 = 1 exactly in the original scalar equation.
    close(r.physical.state_second_derivative[0], C64::ONE);
    let full = factor
        .response_with_frequency_parameter_derivative(
            &input,
            HarmonicParameterDirection::default(),
            None,
            Some(&second),
            None,
            None,
            None,
            || false,
        )
        .unwrap();
    assert_eq!(r.physical, full.physical.frequency_response);
    real(
        r.wave_power_residual.second_frequency,
        full.wave_power_residual.second_frequency,
    );
}

#[test]
fn curvature_has_its_own_lower_budget_and_cancelled_factor_remains_reusable() {
    let mut limits = wave_budget();
    // n=2,m=k=1: curvature requires2032 additional bytes, mixed requires4064.
    limits.maximum_payload_bytes = 4000;
    let model = oscillator(0.17, limits);
    let factor = model.factor(1.7, || false).unwrap();
    let input = [C64::ONE];
    let r = factor
        .response_with_frequency_curvature(&input, None, None, || false)
        .unwrap();
    assert!(matches!(
        factor.response_with_frequency_parameter_derivative(
            &input,
            HarmonicParameterDirection::default(),
            None,
            None,
            None,
            None,
            None,
            || false
        ),
        Err(HarmonicError::Invalid { .. })
    ));
    assert!(matches!(
        factor.response_with_frequency_curvature(&input, None, None, || true),
        Err(HarmonicError::Cancelled)
    ));
    assert!(
        factor
            .response_with_frequency_curvature(&input, None, Some(&[]), || false)
            .is_err()
    );
    assert_eq!(
        r,
        factor
            .response_with_frequency_curvature(&input, None, None, || false)
            .unwrap()
    );
}
