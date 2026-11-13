//! G1/G3: nonlinear saturated supply, original relaxation and generated harmonics.
use super::*;
use fs_math::c64::C64;
use fs_orbit::{HbAnchor, HbOrbit};
use fs_phs::{ContinuousStorageTangent, Storage};

const TAU: f64 = core::f64::consts::TAU;
const M: f64 = 1.2;
const K: f64 = 1.6;
const ALPHA: f64 = 0.1;
const A: f64 = 0.45;
const RELAX: f64 = 0.7;
const C: f64 = 0.15;
const BETA: f64 = 0.4;
const SIGMA: f64 = 0.18;
const ETA: f64 = 0.2;
const GAIN: f64 = 0.9;

#[derive(Debug)]
struct Memory;
impl Storage for Memory {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * K * x[0].powi(2)
            + 0.25 * ALPHA * x[0].powi(4)
            + 0.5 * x[1].powi(2) / M
            + 0.5 * (A * x[0] - x[2]).powi(2)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = K * x[0] + ALPHA * x[0].powi(3) + A * (A * x[0] - x[2]);
        out[1] = x[1] / M;
        out[2] = x[2] - A * x[0];
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
}
impl ContinuousStorageTangent for Memory {
    fn apply(&self, x: &[f64], d: &[f64], out: &mut [f64]) -> bool {
        out[0] = (K + 3.0 * ALPHA * x[0].powi(2) + A * A) * d[0] - A * d[2];
        out[1] = d[1] / M;
        out[2] = d[2] - A * d[0];
        true
    }
}
struct Drive;
impl AutonomousPortFeedback for Drive {
    fn port_count(&self) -> usize {
        1
    }
    fn input(&self, _: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        out[0] = GAIN * flows[0];
        true
    }
}
fn feedback_tangent(_: &[f64], _: &[f64], out: &mut [f64]) -> bool {
    out.copy_from_slice(&[0., GAIN / M, 0.]);
    true
}
fn system() -> PortHamiltonian {
    PortHamiltonian::new(
        3,
        1,
        vec![0., 1., 0., -1., 0., 0., 0., 0., 0.],
        vec![0., 0., 0., 0., C, 0., 0., 0., 1.0 / RELAX],
        vec![0., 1., 0.],
        Box::new(Memory),
    )
    .unwrap()
}
fn loss(x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
    let s = e[1] + ETA * e[2];
    let d = BETA * (1.0 + SIGMA * x[0].powi(2)) * s.powi(3);
    out.copy_from_slice(&[0., d, ETA * d]);
    true
}
fn loss_action(x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]) -> bool {
    let s = e[1] + ETA * e[2];
    let d = 2.0 * BETA * SIGMA * x[0] * dx[0] * s.powi(3)
        + 3.0 * BETA * (1.0 + SIGMA * x[0].powi(2)) * s * s * (de[1] + ETA * de[2]);
    out.copy_from_slice(&[0., d, ETA * d]);
    true
}
fn law() -> NonlinearDissipation<'static> {
    NonlinearDissipation {
        flow: &loss,
        tangent: &loss_action,
        max_tangent_actions: 3,
    }
}
fn budget(order: usize) -> HbBudget {
    HbBudget {
        harmonics: order,
        max_newton: 30,
        tolerance: 2e-11,
        deflate_equilibrium: true,
    }
}
fn scales() -> HbScales {
    HbScales {
        states: vec![1.; 3],
        equations: vec![1.; 3],
        angular_frequency: 1.2,
    }
}
fn controls() -> VariationalShootControls {
    VariationalShootControls {
        state_scales: vec![1.; 3],
        period_scale: 5.2,
        minimum_autonomous_excursion: 0.01,
        maximum_dimension: 3,
        maximum_payload_bytes: 4_000_000,
        maximum_work: 5_000_000_000,
        maximum_period_maps: 40,
    }
}
fn close(a: f64, b: f64, eps: f64) {
    assert!((a - b).abs() < eps, "{a:e} != {b:e}, eps {eps:e}");
}

// Independent original mechanical equations, including the memory force and
// its nonlinear resisting flow. No production field/storage/tangent is called.
fn physical(x: [f64; 3]) -> [f64; 3] {
    let (q, v, stress) = (x[0], x[1] / M, A * x[0] - x[2]);
    let speed = v - ETA * stress;
    let drag = BETA * (1.0 + SIGMA * q * q) * speed.powi(3);
    [
        v,
        GAIN * v - K * q - ALPHA * q.powi(3) - A * stress - C * v - drag,
        stress / RELAX - ETA * drag,
    ]
}
fn rk4(x: [f64; 3], h: f64) -> [f64; 3] {
    let a = physical(x);
    let b = physical(std::array::from_fn(|j| x[j] + h * a[j] / 2.));
    let c = physical(std::array::from_fn(|j| x[j] + h * b[j] / 2.));
    let d = physical(std::array::from_fn(|j| x[j] + h * c[j]));
    std::array::from_fn(|j| x[j] + h * (a[j] + 2. * b[j] + 2. * c[j] + d[j]) / 6.)
}
fn independent(steps: usize, order: usize) -> (HbOrbit, [f64; 3]) {
    let h = 5.2 / steps as f64;
    let mut x = [0.1, 0., 0.];
    let mut last = (0_usize, 0.);
    let mut count = 0;
    let mut period = 0.;
    for i in 0..140 * steps {
        let next = rk4(x, h);
        if x[0] <= 0. && next[0] > 0. {
            let (mut lo, mut hi) = (0., h);
            for _ in 0..40 {
                let middle = 0.5 * (lo + hi);
                if rk4(x, middle)[0] < 0. {
                    lo = middle;
                } else {
                    hi = middle;
                }
            }
            let offset = 0.5 * (lo + hi);
            period = (i - last.0) as f64 * h + offset - last.1;
            last = (i, offset);
            count += 1;
            if count == 120 {
                x = rk4(x, offset);
                break;
            }
        }
        x = next;
    }
    assert_eq!(count, 120);
    let mut coeffs = vec![vec![C64::ZERO; order + 1]; 3];
    let mut powers = [0.; 3];
    for i in 0..steps {
        for (j, row) in coeffs.iter_mut().enumerate() {
            for (n, value) in row.iter_mut().enumerate() {
                let phase = TAU * (n * i) as f64 / steps as f64;
                *value = *value + C64::new(phase.cos(), -phase.sin()).scale(x[j] / steps as f64);
            }
        }
        let v = x[1] / M;
        let stress = A * x[0] - x[2];
        powers[0] += GAIN * v * v / steps as f64;
        powers[1] += (C * v * v + stress * stress / RELAX) / steps as f64;
        powers[2] += BETA * (1. + SIGMA * x[0] * x[0]) * (v - ETA * stress).powi(4) / steps as f64;
        x = rk4(x, period / steps as f64);
    }
    for row in &mut coeffs {
        row[0].im = 0.;
    }
    (
        HbOrbit {
            omega: TAU / period,
            coeffs,
            residual_trace: vec![],
            residual: 0.,
            aft_samples: steps,
            sampling_defect: 0.,
        },
        powers,
    )
}

#[test]
fn original_feedback_and_internal_loss_keep_complete_state_effort_tangent_and_work() {
    let system = system();
    let source =
        DissipativeAutonomousHamiltonian::new(&system, &Drive, &feedback_tangent, law()).unwrap();
    for x in [[0.7, -1.3, 0.2], [-0.4, 1.2, -0.5]] {
        let mut field = [0.; 3];
        source.vector_field_into(&x, &mut field).unwrap();
        for (a, b) in field.into_iter().zip(physical(x)) {
            close(a, b, 2e-15);
        }
        let mut tangent = [0.; 9];
        source.state_jacobian_into(&x, &mut tangent).unwrap();
        for column in 0..3 {
            let (mut plus, mut minus) = (x, x);
            plus[column] += 2e-6;
            minus[column] -= 2e-6;
            let (a, b) = (physical(plus), physical(minus));
            for row in 0..3 {
                close(tangent[3 * row + column], (a[row] - b[row]) / 4e-6, 3e-10);
            }
        }
        let power = source.power(&x).unwrap();
        let v = x[1] / M;
        let stress = A * x[0] - x[2];
        close(power.port_power, GAIN * v * v, 1e-15);
        close(
            power.resistance_power,
            C * v * v + stress * stress / RELAX,
            1e-15,
        );
        close(
            power.nonlinear_power,
            BETA * (1. + SIGMA * x[0] * x[0]) * (v - ETA * stress).powi(4),
            1e-15,
        );
        close(
            power.energy_rate,
            power.port_power - power.resistance_power - power.nonlinear_power,
            2e-15,
        );
    }
}

#[test]
fn self_sustained_memory_cycle_resolves_frequency_harmonics_power_and_variational_period() {
    let system = system();
    let source =
        DissipativeAutonomousHamiltonian::new(&system, &Drive, &feedback_tangent, law()).unwrap();
    let guess = independent(256, 3).0;
    let seed = guess.packed_seed(&HbAnchor::Autonomous {
        omega_guess: guess.omega,
    });
    let adaptive = source
        .solve_harmonic_balance_adaptive_scaled_seeded(
            seed,
            &budget(3),
            &HbRefinementBudget {
                max_harmonics: 31,
                residual_tolerances: vec![2e-7; 3],
                state_tolerances: vec![2e-7; 3],
                frequency_tolerance: 2e-7,
            },
            &scales(),
        )
        .unwrap();
    let orbit = &adaptive.orbit;
    assert!(adaptive.refinement_trace.len() >= 3);
    assert!(orbit.first_harmonic_amplitude(0) > 1.);
    assert!((orbit.omega - (K / M).sqrt()).abs() > 0.03);
    assert!(orbit.coeffs[0][3].abs() > 0.005 && orbit.coeffs[2][3].abs() > 0.001);
    assert!(orbit.coeffs[1][5].abs() > 1e-4);
    let order = orbit.coeffs[0].len() - 1;
    let refs: Vec<_> = [256, 512, 1024]
        .into_iter()
        .map(|n| independent(n, order))
        .collect();
    let errors: Vec<_> = refs
        .iter()
        .map(|(other, _)| (orbit.omega - other.omega).abs())
        .collect();
    assert!(
        errors[2] < 3e-9 && errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    for (row, reference) in orbit.coeffs.iter().zip(&refs[2].0.coeffs) {
        for (a, b) in row.iter().zip(reference) {
            close(a.abs(), b.abs(), 5e-8);
        }
    }
    let mut powers = [0.; 4];
    for i in 0..2048 {
        let x = std::array::from_fn::<_, 3, _>(|j| orbit.sample(j, TAU * i as f64 / 2048.));
        let p = source.power(&x).unwrap();
        for (mean, v) in powers.iter_mut().zip([
            p.port_power,
            p.resistance_power,
            p.nonlinear_power,
            p.energy_rate,
        ]) {
            *mean += v / 2048.;
        }
    }
    for (a, b) in powers.iter().zip(refs[2].1) {
        close(*a, b, 5e-8);
    }
    close(powers[3], 0., 2e-10);
    close(powers[0], powers[1] + powers[2], 2e-10);
    let x: Vec<_> = (0..3).map(|j| orbit.sample(j, 0.)).collect();
    let period = TAU / orbit.omega;
    let mut shooting_errors = Vec::new();
    for steps in [128, 256, 512] {
        let shot = source
            .solve_shooting_variational(
                &x,
                period,
                &ShootBudget {
                    steps_per_period: steps,
                    max_newton: 16,
                    tolerance: 2e-12,
                },
                &controls(),
                || false,
            )
            .unwrap();
        shooting_errors.push((shot.orbit.period - period).abs());
        assert!(shot.maximum_scaled_excursion > 1. && shot.period_derivative.is_some());
        assert!(
            shot.orbit
                .multipliers
                .iter()
                .any(|&(re, im)| C64::new(re - 1., im).abs() < 2e-6)
        );
    }
    assert!(
        shooting_errors[2] < 2e-8
            && shooting_errors[1] < 0.3 * shooting_errors[0]
            && shooting_errors[2] < 0.3 * shooting_errors[1],
        "{shooting_errors:?}"
    );
}

#[test]
fn incomplete_source_tangents_active_loss_and_cancellation_refuse_without_poisoning_reuse() {
    let system = system();
    let mut insufficient = law();
    insufficient.max_tangent_actions = 2;
    assert!(matches!(
        DissipativeAutonomousHamiltonian::new(&system, &Drive, &feedback_tangent, insufficient),
        Err(DissipativeError::Invalid(_))
    ));
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.;
        true
    };
    let source =
        DissipativeAutonomousHamiltonian::new(&system, &Drive, &incomplete, law()).unwrap();
    assert!(matches!(
        source.state_jacobian_into(&[0.2, 0.3, -0.1], &mut [0.; 9]),
        Err(DissipativeError::SourceRefused)
    ));
    let active = |_: &[f64], e: &[f64], out: &mut [f64]| {
        for (a, b) in out.iter_mut().zip(e) {
            *a = -*b;
        }
        true
    };
    let source = DissipativeAutonomousHamiltonian::new(
        &system,
        &Drive,
        &feedback_tangent,
        NonlinearDissipation {
            flow: &active,
            ..law()
        },
    )
    .unwrap();
    assert!(matches!(
        source.power(&[0.2, 0.3, -0.1]),
        Err(DissipativeError::Dissipation(_))
    ));
    let source =
        DissipativeAutonomousHamiltonian::new(&system, &Drive, &feedback_tangent, law()).unwrap();
    assert!(matches!(
        source.solve_shooting_variational(
            &[0., 1.8, -0.38],
            5.2,
            &ShootBudget {
                steps_per_period: 128,
                max_newton: 16,
                tolerance: 1e-10
            },
            &controls(),
            || true
        ),
        Err(DissipativeError::Shooting(VariationalShootError::Cancelled))
    ));
    let mut x = [0.; 3];
    source.vector_field_into(&[0.2, 0.3, -0.1], &mut x).unwrap();
    for (a, b) in x.into_iter().zip(physical([0.2, 0.3, -0.1])) {
        close(a, b, 2e-15);
    }
    assert!(
        source
            .solve_shooting_variational(
                &[0.; 3],
                5.2,
                &ShootBudget {
                    steps_per_period: 128,
                    max_newton: 16,
                    tolerance: 1e-10
                },
                &controls(),
                || false
            )
            .is_err()
    );
}
