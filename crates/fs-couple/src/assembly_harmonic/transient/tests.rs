use super::*;
use fs_math::c64::C64;
use fs_phs::linear_exact::LinearSplitBudget;
use fs_phs::{PortHamiltonian, QuadraticStorage};
use fs_vfit::mimo::{MimoModel, load::MimoLoadLimits};
use fs_vfit::model::PoleTerm;

fn budget() -> AssemblyMimoTransientBudget {
    AssemblyMimoTransientBudget {
        structure: AssemblyResponseBudget {
            maximum_projection_terms: 1_000_000,
            maximum_projection_bytes: 1 << 20,
            linear: LinearSplitBudget {
                max_modes: 8,
                max_omega_dt: 1.,
            },
        },
        exact: DenseExponentialBudget {
            max_dimension: 32,
            max_dense_bytes: 16 << 20,
            max_multiply_terms: 100_000_000,
        },
        maximum_ports: 2,
        maximum_real_states: 12,
        maximum_connection_bytes: 8 << 20,
        maximum_connection_terms: 1_000_000,
        maximum_relative_balance_error: 1e-10,
    }
}

fn system(k: &[f64], inverse_mass: &[f64], c: &[f64], n: usize) -> PortHamiltonian {
    let s = 2 * n;
    let (mut q, mut j, mut r, mut g) = (
        vec![0.; s * s],
        vec![0.; s * s],
        vec![0.; s * s],
        vec![0.; s * n],
    );
    for i in 0..n {
        j[i * s + n + i] = 1.;
        j[(n + i) * s + i] = -1.;
        g[(n + i) * n + i] = 1.;
        for a in 0..n {
            q[i * s + a] = k[i * n + a];
            q[(n + i) * s + n + a] = inverse_mass[i * n + a];
            r[(n + i) * s + n + a] = c[i * n + a];
        }
    }
    PortHamiltonian::new(
        s,
        n,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, s).unwrap()),
    )
    .unwrap()
}

fn source(model: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        model,
        dt,
        MimoLoadLimits {
            band_rad_s: [0.01, 100_000.],
            sampled_passivity_points: 33,
            max_states: 8,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}
fn maxwell(dt: f64) -> PreparedMimoLoad {
    source(
        &MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -3.,
                residue: 0.,
            }],
            residues: vec![vec![C64::from_re(5.)]],
            d: vec![0.],
            e: vec![0.],
            reciprocal: true,
        },
        dt,
    )
}
fn near(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a} != {b}; scale={scale}"
    );
}

// m=1, k=6, Z(s)=5/(s+3), B=1 has characteristic
// (s+2)(s^2+s+9). This is the physical loaded resonance, not the dry sqrt(6).
fn maxwell_oracle(t: f64) -> [f64; 3] {
    let b = 35_f64.sqrt() / 2.;
    let (a, c, d) = (3. / 11., 8. / 11., 10. / (11. * b));
    let (sin, cos) = (b * t).sin_cos();
    let decay = (-0.5 * t).exp();
    let fast = a * (-2. * t).exp();
    let oscillation = c * cos + d * sin;
    let derivative = b * (-c * sin + d * cos);
    let q = fast + decay * oscillation;
    let v = -2. * fast + decay * (-0.5 * oscillation + derivative);
    let acceleration = 4. * fast + decay * ((0.25 - b * b) * oscillation - derivative);
    [q, v, -(acceleration + 6. * q) / 5.]
}

#[test]
fn exact_maxwell_ringdown_retains_the_loaded_frequency_and_decay_with_zero_direct_term() {
    let sys = system(&[6.], &[1.], &[0.], 1);
    for duration in [0.02, 0.8, 2.4, 7.] {
        let load = maxwell(duration);
        let joint = PreparedJointFlow::new(
            &sys,
            &[6.],
            &[1.],
            &[0.],
            &[1.],
            &load,
            budget(),
            &mut || false,
        )
        .unwrap();
        let step = joint.flow.evaluate(&[1., 0., 0.], &[0.]).unwrap();
        let expected = maxwell_oracle(duration);
        for (a, b) in step.state.iter().zip(expected) {
            near(*a, b, 1., 3e-12);
        }
        near(step.integrals[0], 0., 1., 1e-14);
        let delta_h = joint::energy_change(&[6.], &[1.], &[1., 0.], &step.state[..2]).unwrap();
        near(delta_h + step.integrals[1], 0., 3., 3e-12);
        near(step.integrals[2], step.state[0] - 1., 1., 3e-12);
        let staged = load
            .continuous_operator()
            .candidate_from_coordinates(
                &step.state[2..],
                &joint.endpoint_flows(&step.state).unwrap(),
            )
            .unwrap();
        assert_eq!(staged.source_model(), load.source_model());
        assert_eq!(staged.limits(), load.limits());
        assert_eq!(staged.time_step_s().to_bits(), duration.to_bits());
        assert_eq!(staged.reports().0.class, load.reports().0.class);
        assert_eq!(staged.reports().1.class, load.reports().1.class);
        near(staged.pole_states()[0].re, expected[2], 1., 3e-12);
    }
}

#[test]
fn held_force_and_existing_load_history_exchange_signed_work_without_resetting_memory() {
    let sys = system(&[6.], &[1.], &[0.], 1);
    let load = maxwell(0.4);
    let joint = PreparedJointFlow::new(
        &sys,
        &[6.],
        &[1.],
        &[0.],
        &[1.],
        &load,
        budget(),
        &mut || false,
    )
    .unwrap();
    let forced = joint.flow.evaluate(&[1., 0., 0.], &[0.7]).unwrap();
    let dry = maxwell_oracle(0.4);
    let equilibrium = 0.7 / 6.;
    for (actual, expected) in forced.state.iter().zip([
        equilibrium + (1. - equilibrium) * dry[0],
        (1. - equilibrium) * dry[1],
        (1. - equilibrium) * dry[2],
    ]) {
        near(*actual, expected, 1., 3e-12);
    }
    let delta_h = joint::energy_change(&[6.], &[1.], &[1., 0.], &forced.state[..2]).unwrap();
    near(
        delta_h + forced.integrals[1],
        0.7 * (forced.state[0] - 1.),
        3.,
        3e-12,
    );

    let release = joint.flow.evaluate(&[0., 0., 1.], &[0.]).unwrap();
    assert!(release.integrals[1] < 0.); // original memory supplies structural energy
    let structural_energy = sys.hamiltonian(&release.state[..2]);
    near(
        structural_energy + release.integrals[1],
        0.,
        structural_energy,
        3e-12,
    );
    assert!(release.state[2] != 0.);
}

#[test]
fn full_mass_mutual_load_and_nonproportional_glue_match_physical_ode_and_semigroup() {
    let (m, k, c, b) = (
        [2., 0.2, 0.2, 1.],
        [8., -2., -2., 5.],
        [0.26, -0.06, -0.06, 0.2],
        [1., 0.3, -0.2, 0.7],
    );
    let determinant = m[0] * m[3] - m[1] * m[2];
    let inverse = [
        m[3] / determinant,
        -m[1] / determinant,
        -m[2] / determinant,
        m[0] / determinant,
    ];
    let residue = [5., 1., 1., 2.];
    let direct = [0.3, 0.08, -0.02, 0.2]; // retains a passive nonreciprocal part
    let model = MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Real {
            pole: -3.,
            residue: 0.,
        }],
        residues: vec![residue.map(C64::from_re).to_vec()],
        d: direct.to_vec(),
        e: vec![0.; 4],
        reciprocal: false,
    };
    let sys = system(&k, &inverse, &c, 2);
    let initial = [0.2, -0.1, 0.61, 0.11, 0.1, -0.2]; // canonical p=M*[0.3,0.05]
    let forcing = [0.7, -0.4];
    let mut reference = vec![0.2, -0.1, 0.3, 0.05, 0.1, -0.2, 0., 0., 0.];
    // Independent physical q,v,z equations and direct trajectory work. RK4
    // exists only in this test; production uses the original exact owner.
    let derivative = |x: &[f64]| {
        let q = [x[0], x[1]];
        let v = [x[2], x[3]];
        let z = [x[4], x[5]];
        let u = [b[0] * v[0] + b[2] * v[1], b[1] * v[0] + b[3] * v[1]];
        let e = [
            residue[0] * z[0] + residue[1] * z[1] + direct[0] * u[0] + direct[1] * u[1],
            residue[2] * z[0] + residue[3] * z[1] + direct[2] * u[0] + direct[3] * u[1],
        ];
        let force = [
            forcing[0]
                - k[0] * q[0]
                - k[1] * q[1]
                - c[0] * v[0]
                - c[1] * v[1]
                - b[0] * e[0]
                - b[1] * e[1],
            forcing[1]
                - k[2] * q[0]
                - k[3] * q[1]
                - c[2] * v[0]
                - c[3] * v[1]
                - b[2] * e[0]
                - b[3] * e[1],
        ];
        vec![
            v[0],
            v[1],
            (m[3] * force[0] - m[1] * force[1]) / determinant,
            (-m[2] * force[0] + m[0] * force[1]) / determinant,
            -3. * z[0] + u[0],
            -3. * z[1] + u[1],
            v[0] * (c[0] * v[0] + c[1] * v[1]) + v[1] * (c[2] * v[0] + c[3] * v[1]),
            e[0] * u[0] + e[1] * u[1],
            forcing[0] * v[0] + forcing[1] * v[1],
        ]
    };
    let duration = 0.8;
    let dt = duration / 4096.;
    for _ in 0..4096 {
        let a = derivative(&reference);
        let shifted = |slope: &[f64], h: f64| {
            reference
                .iter()
                .zip(slope)
                .map(|(x, dx)| x + h * dx)
                .collect::<Vec<_>>()
        };
        let b2 = derivative(&shifted(&a, 0.5 * dt));
        let c2 = derivative(&shifted(&b2, 0.5 * dt));
        let d = derivative(&shifted(&c2, dt));
        for i in 0..reference.len() {
            reference[i] += dt / 6. * (a[i] + 2. * b2[i] + 2. * c2[i] + d[i]);
        }
    }
    let mut expected = reference[..6].to_vec();
    expected[2] = m[0] * reference[2] + m[1] * reference[3];
    expected[3] = m[2] * reference[2] + m[3] * reference[3];
    let mut one_step = Vec::new();
    let mut one_integrals = Vec::new();
    for subdivisions in [1, 2, 4] {
        let load = source(&model, duration / subdivisions as f64);
        let joint =
            PreparedJointFlow::new(&sys, &k, &inverse, &c, &b, &load, budget(), &mut || false)
                .unwrap();
        let mut state = initial.to_vec();
        let mut integrals = vec![0.; 6];
        for _ in 0..subdivisions {
            let next = joint.flow.evaluate(&state, &forcing).unwrap();
            state = next.state;
            for (sum, value) in integrals.iter_mut().zip(next.integrals) {
                *sum += value;
            }
        }
        for (a, b) in state.iter().zip(&expected) {
            near(*a, *b, 1., 3e-11);
        }
        near(integrals[0], reference[6], 1., 3e-11);
        near(integrals[1], reference[7], 1., 3e-11);
        near(
            forcing[0] * integrals[2] + forcing[1] * integrals[3],
            reference[8],
            1.,
            3e-11,
        );
        let change = joint::energy_change(&k, &inverse, &initial[..4], &state[..4]).unwrap();
        near(
            change + integrals[0] + integrals[1] - reference[8],
            0.,
            1.,
            3e-11,
        );
        if subdivisions == 1 {
            one_step = state;
            one_integrals = integrals;
        } else {
            for (a, b) in state.iter().zip(&one_step) {
                near(*a, *b, 1., 3e-12);
            }
            for (a, b) in integrals.iter().zip(&one_integrals) {
                near(*a, *b, 1., 3e-12);
            }
        }
    }
}

#[test]
fn original_pair_chart_is_preserved_and_complete_capacity_and_cancellation_refuse() {
    let sys = system(&[6.], &[1.], &[0.2], 1);
    let mut load = source(
        &MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Pair {
                pole: C64::new(-0.3, 4.),
                residue: C64::ZERO,
            }],
            residues: vec![vec![C64::from_re(1.7)]],
            d: vec![0.],
            e: vec![0.],
            reciprocal: true,
        },
        0.3,
    );
    load.step_mean_into(&[0.8], &mut [0.]).unwrap();
    let original_history = load.pole_states().to_vec();
    let original = load.continuous_operator().real_coordinates();
    assert_eq!(original[0], original_history[0].re);
    assert_eq!(original[1], -original_history[0].im);
    let joint = PreparedJointFlow::new(
        &sys,
        &[6.],
        &[1.],
        &[0.2],
        &[1.],
        &load,
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut initial = vec![0.1, -0.2];
    initial.extend_from_slice(&original);
    let endpoint = joint.flow.evaluate(&initial, &[0.]).unwrap();
    let candidate = load
        .continuous_operator()
        .candidate_from_coordinates(
            &endpoint.state[2..],
            &joint.endpoint_flows(&endpoint.state).unwrap(),
        )
        .unwrap();
    assert_eq!(candidate.pole_states()[0].re, endpoint.state[2]);
    assert_eq!(candidate.pole_states()[0].im, -endpoint.state[3]);
    for fail_at in [0, 1, 4, 8] {
        let mut calls = 0;
        let result = PreparedJointFlow::new(
            &sys,
            &[6.],
            &[1.],
            &[0.2],
            &[1.],
            &load,
            budget(),
            &mut || {
                calls += 1;
                calls > fail_at
            },
        );
        assert!(result.is_err());
    }
    for kind in 0..6 {
        let mut cap = budget();
        match kind {
            0 => cap.maximum_real_states = 3,
            1 => cap.maximum_ports = 0,
            2 => cap.maximum_connection_bytes = 1,
            3 => cap.maximum_connection_terms = 0,
            4 => cap.exact.max_dense_bytes = 1,
            _ => cap.exact.max_multiply_terms = 0,
        }
        assert!(
            PreparedJointFlow::new(&sys, &[6.], &[1.], &[0.2], &[1.], &load, cap, &mut || false,)
                .is_err()
        );
    }
    assert_eq!(load.pole_states(), original_history);
}
