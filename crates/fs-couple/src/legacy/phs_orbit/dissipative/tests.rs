use super::*;
use fs_math::c64::C64;
use fs_orbit::HbOrbit;
use fs_phs::{ContinuousStorageTangent, QuadraticStorage, Storage};

const TAU: f64 = core::f64::consts::TAU;
const M: f64 = 1.2;
const K: f64 = 1.6;
const ALPHA: f64 = 0.15;
const A: f64 = 0.7;
const RELAX: f64 = 0.8;
const C: f64 = 0.25;
const BETA: f64 = 0.4;
const SIGMA: f64 = 0.3;
const ETA: f64 = 0.2;
const OMEGA: f64 = 0.9;

#[derive(Debug)]
struct Memory;
impl Storage for Memory {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * K * x[0] * x[0]
            + 0.25 * ALPHA * x[0].powi(4)
            + 0.5 * x[1] * x[1] / M
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
        out[0] = (K + 3.0 * ALPHA * x[0] * x[0] + A * A) * d[0] - A * d[2];
        out[1] = d[1] / M;
        out[2] = d[2] - A * d[0];
        true
    }
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
fn force(t: f64) -> f64 {
    0.14 + 0.5 * (OMEGA * t).cos() + 0.1 * (2.0 * OMEGA * t + 0.37).cos()
}
fn input() -> PeriodicPortInput {
    PeriodicPortInput::new(
        OMEGA,
        vec![vec![
            C64::new(0.14, 0.),
            C64::new(0.25, 0.),
            C64::new(0.05 * 0.37_f64.cos(), 0.05 * 0.37_f64.sin()),
        ]],
    )
    .unwrap()
}
fn loss(x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
    let s = e[1] + ETA * e[2];
    let d = BETA * (1.0 + SIGMA * x[0] * x[0]) * s.powi(3);
    out.copy_from_slice(&[0., d, ETA * d]);
    true
}
fn loss_action(x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]) -> bool {
    let s = e[1] + ETA * e[2];
    let d = 2.0 * BETA * SIGMA * x[0] * dx[0] * s.powi(3)
        + 3.0 * BETA * (1.0 + SIGMA * x[0] * x[0]) * s * s * (de[1] + ETA * de[2]);
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
fn scales() -> HbScales {
    HbScales {
        states: vec![1.; 3],
        equations: vec![1.; 3],
        angular_frequency: OMEGA,
    }
}
fn budget(order: usize) -> HbBudget {
    HbBudget {
        harmonics: order,
        max_newton: 24,
        tolerance: 2e-12,
        deflate_equilibrium: false,
    }
}
fn shooting_controls() -> VariationalShootControls {
    VariationalShootControls {
        state_scales: vec![1.; 3],
        period_scale: TAU / OMEGA,
        minimum_autonomous_excursion: 1e-5,
        maximum_dimension: 3,
        maximum_payload_bytes: 4_000_000,
        maximum_work: 5_000_000_000,
        maximum_period_maps: 32,
    }
}

// Independently written physical equations: no production field, loss,
// storage, tangent, time integrator, HB, or observation callback is used.
fn physical(t: f64, x: &[f64; 3]) -> [f64; 3] {
    let (q, v, stress) = (x[0], x[1] / M, A * x[0] - x[2]);
    let s = v - ETA * stress;
    let drag = BETA * (1.0 + SIGMA * q * q) * s * s * s;
    [
        v,
        force(t) - K * q - ALPHA * q * q * q - A * stress - C * v - drag,
        stress / RELAX - ETA * drag,
    ]
}
fn rk4(t: f64, x: [f64; 3], h: f64) -> [f64; 3] {
    let k1 = physical(t, &x);
    let k2 = physical(t + h / 2., &std::array::from_fn(|j| x[j] + h * k1[j] / 2.));
    let k3 = physical(t + h / 2., &std::array::from_fn(|j| x[j] + h * k2[j] / 2.));
    let k4 = physical(t + h, &std::array::from_fn(|j| x[j] + h * k3[j]));
    std::array::from_fn(|j| x[j] + h * (k1[j] + 2. * k2[j] + 2. * k3[j] + k4[j]) / 6.)
}
fn independent(steps: usize, order: usize) -> (Vec<Vec<C64>>, [f64; 3]) {
    let h = TAU / (OMEGA * steps as f64);
    let mut x = [0.; 3];
    for i in 0..60 * steps {
        x = rk4(i as f64 * h, x, h);
    }
    let mut coefficients = vec![vec![C64::ZERO; order + 1]; 3];
    let mut powers = [0.; 3];
    for i in 0..steps {
        let t = i as f64 * h;
        for (j, row) in coefficients.iter_mut().enumerate() {
            for (n, value) in row.iter_mut().enumerate() {
                let phase = n as f64 * OMEGA * t;
                *value = *value + C64::new(phase.cos(), -phase.sin()).scale(x[j] / steps as f64);
            }
        }
        let v = x[1] / M;
        let stress = A * x[0] - x[2];
        let s = v - ETA * stress;
        powers[0] += force(t) * v / steps as f64;
        powers[1] += (C * v * v + stress * stress / RELAX) / steps as f64;
        powers[2] += BETA * (1.0 + SIGMA * x[0] * x[0]) * s.powi(4) / steps as f64;
        x = rk4(t, x, h);
    }
    (coefficients, powers)
}
fn error(orbit: &HbOrbit, other: &[Vec<C64>]) -> f64 {
    orbit
        .coeffs
        .iter()
        .zip(other)
        .flat_map(|(a, b)| a.iter().zip(b))
        .map(|(a, b)| (*a - *b).abs())
        .fold(0., f64::max)
}

#[test]
fn full_nonlinear_loss_and_memory_generate_resolved_harmonics_and_separate_heat() {
    let (system, input) = (system(), input());
    let source = DissipativeHamiltonian::new(&system, &input, law()).unwrap();
    let adaptive = source
        .solve_harmonic_balance_adaptive_scaled_seeded(
            vec![0.; 3 * 7],
            &budget(3),
            &HbRefinementBudget {
                max_harmonics: 31,
                residual_tolerances: vec![2e-7; 3],
                state_tolerances: vec![2e-7; 3],
                frequency_tolerance: 1e-10,
            },
            &scales(),
        )
        .unwrap();
    assert!(adaptive.refinement_trace.len() >= 3);
    let orbit = &adaptive.orbit;
    let order = orbit.coeffs[0].len() - 1;
    assert!(orbit.coeffs[1][3].abs() > 1e-4);
    assert!(orbit.coeffs[2][3].abs() > 1e-6);
    assert!(orbit.coeffs[1][5].abs() > 1e-6);
    let coarse = independent(512, order);
    let middle = independent(1024, order);
    let fine = independent(2048, order);
    let errors = [
        error(orbit, &coarse.0),
        error(orbit, &middle.0),
        error(orbit, &fine.0),
    ];
    assert!(
        errors[2] < 2e-9 && errors[1] < 0.2 * errors[0] && errors[2] < 0.2 * errors[1],
        "{errors:?}"
    );
    let mut means = [0.; 4];
    for i in 0..2048 {
        let t = TAU * i as f64 / (OMEGA * 2048.);
        let x = std::array::from_fn::<_, 3, _>(|j| orbit.sample(j, OMEGA * t));
        let p = source.power(t, &x).unwrap();
        assert!(p.nonlinear_power >= 0. && p.resistance_power >= 0.);
        let balance = p.port_power - p.resistance_power - p.nonlinear_power - p.energy_rate;
        assert!(balance.abs() < 1e-14);
        for (mean, value) in means.iter_mut().zip([
            p.port_power,
            p.resistance_power,
            p.nonlinear_power,
            p.energy_rate,
        ]) {
            *mean += value / 2048.;
        }
    }
    for j in 0..3 {
        assert!(
            (means[j] - fine.1[j]).abs() < 2e-9,
            "{means:?} {:?}",
            fine.1
        );
    }
    assert!(means[2] > 1e-4 && means[1] > means[2]);
    assert!(means[3].abs() < 1e-10);
    assert!((means[0] - means[1] - means[2]).abs() < 1e-10);
}

#[test]
fn full_state_effort_chain_matches_independent_tangent_and_variational_shooting() {
    let (system, input) = (system(), input());
    let source = DissipativeHamiltonian::new(&system, &input, law()).unwrap();
    let x = [0.6, -0.7, 0.2];
    let mut jac = [0.; 9];
    source.state_jacobian_into(0.31, &x, &mut jac).unwrap();
    let v = x[1] / M;
    let s = v + ETA * (x[2] - A * x[0]);
    let coefficient = BETA * (1. + SIGMA * x[0] * x[0]);
    let drag_q = 2. * BETA * SIGMA * x[0] * s.powi(3) - 3. * coefficient * s * s * ETA * A;
    let drag_p = 3. * coefficient * s * s / M;
    let drag_z = 3. * coefficient * s * s * ETA;
    let expected = [
        0.,
        1. / M,
        0.,
        -(K + 3. * ALPHA * x[0] * x[0] + A * A) - drag_q,
        -C / M - drag_p,
        A - drag_z,
        A / RELAX - ETA * drag_q,
        -ETA * drag_p,
        -1. / RELAX - ETA * drag_z,
    ];
    for (a, b) in jac.iter().zip(expected) {
        assert!((a - b).abs() < 1e-14, "{jac:?} {expected:?}");
    }
    let orbit = source
        .solve_harmonic_balance_scaled_seeded(vec![0.; 3 * 31], &budget(15), &scales())
        .unwrap()
        .orbit;
    let seed: Vec<_> = (0..3).map(|j| orbit.sample(j, 0.)).collect();
    let mut shooting = Vec::new();
    for steps in [256, 512, 1024] {
        let shot = source
            .solve_shooting_variational(
                &seed,
                &ShootBudget {
                    steps_per_period: steps,
                    max_newton: 12,
                    tolerance: 2e-12,
                },
                &shooting_controls(),
                || false,
            )
            .unwrap();
        assert!(shot.jacobian_evaluations > 0);
        let discrepancy = shot
            .orbit
            .x0
            .iter()
            .zip(&seed)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        shooting.push(discrepancy);
    }
    assert!(
        shooting[2] < 3e-9 && shooting[1] < 0.2 * shooting[0] && shooting[2] < 0.2 * shooting[1],
        "{shooting:?}"
    );
}

#[test]
fn passive_negative_differential_slope_is_kept_instead_of_a_secant_damper() {
    let system = PortHamiltonian::new(
        1,
        1,
        vec![0.],
        vec![0.],
        vec![1.],
        Box::new(QuadraticStorage::new(vec![1.], 1).unwrap()),
    )
    .unwrap();
    let input = PeriodicPortInput::new(1., vec![vec![C64::new(0.4, 0.)]]).unwrap();
    let flow = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = e[0] / (1. + e[0] * e[0]);
        true
    };
    let tangent = |_: &[f64], e: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = (1. - e[0] * e[0]) / (1. + e[0] * e[0]).powi(2) * de[0];
        true
    };
    let source = DissipativeHamiltonian::new(
        &system,
        &input,
        NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 1,
        },
    )
    .unwrap();
    let mut field = [0.];
    let mut jac = [0.];
    source.vector_field_into(0., &[2.], &mut field).unwrap();
    source.state_jacobian_into(0., &[2.], &mut jac).unwrap();
    assert_eq!(field, [0.]);
    assert!((jac[0] - 0.12).abs() < 1e-15);
    let p = source.power(0., &[2.]).unwrap();
    assert!((p.nonlinear_power - 0.8).abs() < 1e-15);
    assert_eq!(p.energy_rate, 0.);
}

#[test]
fn source_admission_refusals_and_cancellation_do_not_poison_reuse() {
    let (system, input) = (system(), input());
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.;
        true
    };
    let partial = DissipativeHamiltonian::new(
        &system,
        &input,
        NonlinearDissipation {
            flow: &incomplete,
            tangent: &loss_action,
            max_tangent_actions: 3,
        },
    )
    .unwrap();
    assert!(matches!(
        partial.vector_field_into(0., &[1., 1., 1.], &mut [0.; 3]),
        Err(DissipativeError::Dissipation(
            LinearizationError::DissipationRefused { column: None }
        ))
    ));
    let active = |_: &[f64], e: &[f64], out: &mut [f64]| {
        for (o, e) in out.iter_mut().zip(e) {
            *o = -*e;
        }
        true
    };
    let source = DissipativeHamiltonian::new(
        &system,
        &input,
        NonlinearDissipation {
            flow: &active,
            tangent: &loss_action,
            max_tangent_actions: 3,
        },
    )
    .unwrap();
    assert!(matches!(
        source.vector_field_into(0., &[1., 1., 1.], &mut [0.; 3]),
        Err(DissipativeError::Dissipation(
            LinearizationError::ActiveDissipation { .. }
        ))
    ));
    assert!(
        DissipativeHamiltonian::new(
            &system,
            &input,
            NonlinearDissipation {
                flow: &loss,
                tangent: &loss_action,
                max_tangent_actions: 2
            }
        )
        .is_err()
    );
    let source = DissipativeHamiltonian::new(&system, &input, law()).unwrap();
    assert!(matches!(
        source.solve_harmonic_balance_scaled_seeded(vec![0.; 9], &budget(1), &scales()),
        Err(DissipativeError::Invalid(_))
    ));
    let b = ShootBudget {
        steps_per_period: 64,
        max_newton: 10,
        tolerance: 1e-10,
    };
    assert!(matches!(
        source.solve_shooting_variational(&[0.; 3], &b, &shooting_controls(), || true),
        Err(DissipativeError::Shooting(VariationalShootError::Cancelled))
    ));
    let mut f = [0.; 3];
    source
        .vector_field_into(0.3, &[0.1, 0.2, 0.3], &mut f)
        .unwrap();
    let mut retry = [0.; 3];
    source
        .vector_field_into(0.3, &[0.1, 0.2, 0.3], &mut retry)
        .unwrap();
    assert_eq!(f, retry);
}

#[test]
fn shared_physical_contractions_preserve_weak_flow_and_internal_work() {
    let big = 2.0_f64.powi(54);
    let system = PortHamiltonian::new(
        1,
        2,
        vec![0.],
        vec![0.],
        vec![big, 1.],
        Box::new(QuadraticStorage::new(vec![1.], 1).unwrap()),
    )
    .unwrap();
    let mut flow = [0.];
    system
        .vector_field_with_dissipation_from_effort_into(&[1.], &[1., 1.], &[big], &mut flow)
        .unwrap();
    assert_eq!(flow, [1.]);
    let malformed = PortHamiltonian::from_raw_parts(
        1,
        0,
        vec![],
        vec![0.],
        vec![],
        Box::new(QuadraticStorage::new(vec![1.], 1).unwrap()),
    );
    assert!(
        malformed
            .vector_field_with_dissipation_from_effort_into(&[1.], &[], &[1.], &mut flow,)
            .is_err()
    );
    let d = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[big, 1., -big]);
        true
    };
    let unused = |_: &[f64], _: &[f64], _: &[f64], _: &[f64], _: &mut [f64]| false;
    let law = NonlinearDissipation {
        flow: &d,
        tangent: &unused,
        max_tangent_actions: 3,
    };
    assert_eq!(
        law.evaluate_into(&[0.; 3], &[1.; 3], &mut [0.; 3], || false)
            .unwrap(),
        1.
    );
}
