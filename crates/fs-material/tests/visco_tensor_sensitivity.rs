//! G0/G1/G3 derivatives of original anisotropic damping memory and heat.
//! Authored tensor fixtures; no measured-material or coupled-solver claim.

use fs_material::visco_tensor::{
    TensorProny, TensorPronyBranch, TensorPronyError, TensorPronyStep, TensorPronyStepDerivative,
    TensorPronyStepDirection, TensorPronyStepLinearization,
};
use fs_material::{Tangent6, Voigt};
use fs_math::det;

const FACTORS: [Tangent6; 2] = [
    [
        [2.0, -0.75, 0.5, 1.5, 0.0, -0.25],
        [0.0, 1.0, -0.5, 0.0, 0.75, 0.4],
        [0.6, 0.0, 0.0, 0.3, -0.9, 0.2],
        [0.0; 6],
        [0.0; 6],
        [0.0; 6],
    ],
    [
        [0.25, 0.0, 1.5, -0.5, 0.25, 1.0],
        [0.5, -0.5, 0.0, 0.25, 0.0, 0.75],
        [0.0; 6],
        [0.0; 6],
        [0.0; 6],
        [0.0; 6],
    ],
];
const STRAINS: [Voigt; 5] = [
    [0.01, -0.02, 0.015, 0.025, -0.01, 0.012],
    [-0.015, 0.012, 0.02, -0.008, 0.015, 0.006],
    [-0.015, 0.012, 0.02, -0.008, 0.015, 0.006],
    [0.018, -0.004, -0.007, 0.016, 0.008, -0.015],
    [0.0; 6],
];
const TIMES: [f64; 5] = [0.03, 0.08, 0.2, 0.017, 0.4];

fn dot(a: &Voigt, b: &Voigt) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn apply(c: &Tangent6, v: &Voigt) -> Voigt {
    c.map(|row| dot(&row, v))
}
fn gram(b: &Tangent6) -> Tangent6 {
    core::array::from_fn(|i| core::array::from_fn(|j| (0..6).map(|k| b[k][i] * b[k][j]).sum()))
}
fn mandel(mut v: Voigt) -> Voigt {
    for x in &mut v[3..] {
        *x *= core::f64::consts::SQRT_2;
    }
    v
}
fn material(times: [f64; 2], kelvin: bool) -> TensorProny {
    let mut eq = core::array::from_fn(|i| {
        core::array::from_fn(|j| if i == j { 12.0 + i as f64 } else { 0.0 })
    });
    eq[0][3] = 0.7;
    eq[3][0] = 0.7;
    eq[1][4] = -0.4;
    eq[4][1] = -0.4;
    let branches = FACTORS
        .iter()
        .zip(times)
        .map(|(b, tau)| TensorPronyBranch::from_gram_rows(*b, tau).unwrap())
        .collect();
    let law = TensorProny::new_mandel(eq, branches, 0.2).unwrap();
    if kelvin {
        let mut factor = [[0.0; 6]; 6];
        factor[0] = [0.2, -0.1, 0.0, 0.3, 0.05, -0.08];
        factor[1] = [0.0, 0.12, 0.1, 0.0, -0.2, 0.04];
        law.with_kelvin_viscosity_mandel(gram(&factor)).unwrap()
    } else {
        law
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * (a.abs().max(b.abs()).max(1e-12)),
        "{a:e} != {b:e} (relative tolerance {tolerance:e})"
    );
}
fn pair(a: TensorPronyStepDerivative, b: TensorPronyStepDerivative) -> f64 {
    dot(&a.stress, &b.stress)
        + dot(&a.mean_stress, &b.mean_stress)
        + a.stored_energy * b.stored_energy
        + a.work * b.work
        + a.dissipated * b.dissipated
}
fn observations(step: &TensorPronyStep) -> TensorPronyStepDerivative {
    TensorPronyStepDerivative {
        stress: step.stress,
        mean_stress: step.mean_stress,
        stored_energy: step.stored_energy,
        work: step.work,
        dissipated: step.dissipated,
    }
}
fn seed(step: usize) -> TensorPronyStepDerivative {
    TensorPronyStepDerivative {
        stress: [0.3 + step as f64 * 0.07, -0.2, 0.4, -0.5, 0.15, 0.08],
        mean_stress: [-0.1, 0.25, -0.3, 0.6, -0.12, 0.3],
        stored_energy: 0.4,
        work: -0.25,
        dissipated: 1.7,
    }
}
fn objective(law: &TensorProny, strains: &[Voigt; 5]) -> f64 {
    let mut history = law.history();
    strains
        .iter()
        .zip(TIMES)
        .enumerate()
        .map(|(k, (&strain, dt))| {
            pair(
                seed(k),
                observations(&history.advance(strain, dt, || false).unwrap()),
            )
        })
        .sum()
}

#[test]
fn g0_full_tensor_tangent_pullback_duality_keeps_original_primal_bits() {
    let law = material([0.05, 0.7], true);
    let mut history = law.history();
    history.advance(STRAINS[0], TIMES[0], || false).unwrap();
    let mut original = history.clone();
    let expected = original.advance(STRAINS[1], TIMES[1], || false).unwrap();
    let tape = history
        .advance_linearized(STRAINS[1], TIMES[1], || false)
        .unwrap();
    assert_eq!(tape.response(), &expected);
    assert_eq!(history, original);
    let previous_strain = [0.03, -0.01, 0.02, 0.04, -0.02, 0.01];
    let strain = [-0.02, 0.05, -0.03, 0.01, 0.03, -0.04];
    let memory_direction = [
        [0.02, -0.03, 0.01, 0.05, -0.01, 0.04],
        [-0.01, 0.04, -0.02, 0.03, -0.05, 0.01],
    ];
    let log_times = [0.7, -0.4];
    let mut next_memory = [[0.0; 6]; 2];
    let pushed = tape
        .push_forward(
            TensorPronyStepDirection {
                previous_strain,
                previous_memory: &memory_direction,
                strain,
                log_relaxation_times: &log_times,
            },
            &mut next_memory,
            || false,
        )
        .unwrap();
    let outgoing_seed = [
        [0.3, -0.2, 0.1, 0.4, -0.5, 0.6],
        [0.2, 0.1, -0.3, 0.5, -0.4, 0.3],
    ];
    let mut previous_memory = [[0.0; 6]; 2];
    let mut time_gradient = [0.0; 2];
    let strains = tape
        .pull_back(
            seed(1),
            &outgoing_seed,
            &mut previous_memory,
            &mut time_gradient,
            || false,
        )
        .unwrap();
    let lhs = pair(seed(1), pushed)
        + (0..2)
            .map(|j| dot(&outgoing_seed[j], &next_memory[j]))
            .sum::<f64>();
    let rhs = dot(&strains.previous_strain, &previous_strain)
        + dot(&strains.strain, &strain)
        + (0..2)
            .map(|j| {
                dot(&previous_memory[j], &memory_direction[j]) + time_gradient[j] * log_times[j]
            })
            .sum::<f64>();
    close(lhs, rhs, 5e-13);
    // The endpoint block also agrees with the already exposed native tangent,
    // including normal/shear cross terms and physical tensor-shear weights.
    let endpoint = tape
        .push_forward(
            TensorPronyStepDirection {
                previous_strain: [0.0; 6],
                previous_memory: &[[0.0; 6]; 2],
                strain,
                log_relaxation_times: &[0.0; 2],
            },
            &mut next_memory,
            || false,
        )
        .unwrap();
    let expected = apply(&law.tangent(TIMES[1]).unwrap(), &strain);
    for i in 0..6 {
        close(endpoint.stress[i], expected[i], 3e-13);
    }
}

#[test]
fn g1_multistep_loading_and_log_time_adjoints_match_actual_history_differences() {
    let times = [0.05, 0.7];
    let law = material(times, true);
    let mut history = law.history();
    let tapes: Vec<TensorPronyStepLinearization> = STRAINS
        .iter()
        .zip(TIMES)
        .map(|(&strain, dt)| history.advance_linearized(strain, dt, || false).unwrap())
        .collect();
    let loading: [Voigt; 5] = core::array::from_fn(|k| {
        core::array::from_fn(|i| (k as f64 - 1.5) * (i as f64 - 2.2) * 0.003)
    });
    let time_direction = [0.6, -0.35];
    let mut previous_strain = [0.0; 6];
    let mut memory = [[0.0; 6]; 2];
    let mut next_memory = memory;
    let (mut forward, mut work, mut heat, mut final_storage) = (0.0, 0.0, 0.0, 0.0);
    for (k, tape) in tapes.iter().enumerate() {
        let pushed = tape
            .push_forward(
                TensorPronyStepDirection {
                    previous_strain,
                    previous_memory: &memory,
                    strain: loading[k],
                    log_relaxation_times: &time_direction,
                },
                &mut next_memory,
                || false,
            )
            .unwrap();
        memory = next_memory;
        previous_strain = loading[k];
        forward += pair(seed(k), pushed);
        work += pushed.work;
        heat += pushed.dissipated;
        final_storage = pushed.stored_energy;
    }
    close(work, final_storage + heat, 3e-12);
    let mut memory_seed = [[0.0; 6]; 2];
    let mut incoming_seed = memory_seed;
    let mut total_times = [0.0; 2];
    let mut loading_gradient = [[0.0; 6]; 5];
    for k in (0..tapes.len()).rev() {
        let mut times = [0.0; 2];
        let pulled = tapes[k]
            .pull_back(
                seed(k),
                &memory_seed,
                &mut incoming_seed,
                &mut times,
                || false,
            )
            .unwrap();
        memory_seed = incoming_seed;
        for j in 0..2 {
            total_times[j] += times[j];
        }
        for i in 0..6 {
            loading_gradient[k][i] += pulled.strain[i];
            if k > 0 {
                loading_gradient[k - 1][i] += pulled.previous_strain[i];
            }
        }
    }
    let reverse = (0..2)
        .map(|j| total_times[j] * time_direction[j])
        .sum::<f64>()
        + (0..5)
            .map(|k| dot(&loading_gradient[k], &loading[k]))
            .sum::<f64>();
    close(forward, reverse, 2e-12);
    // Perturbed constitutive histories are a test oracle only. Production
    // propagates the exact retained-memory derivatives without rerunning them.
    let h = 1e-5;
    for j in 0..2 {
        let mut plus = times;
        plus[j] *= det::exp(h);
        let mut minus = times;
        minus[j] *= det::exp(-h);
        let difference = (objective(&material(plus, true), &STRAINS)
            - objective(&material(minus, true), &STRAINS))
            / (2.0 * h);
        close(total_times[j], difference, 3e-8);
    }
    for (k, i) in [(0, 3), (1, 0), (2, 5), (4, 4)] {
        let mut plus = STRAINS;
        plus[k][i] += h;
        let mut minus = STRAINS;
        minus[k][i] -= h;
        let difference = (objective(&law, &plus) - objective(&law, &minus)) / (2.0 * h);
        close(loading_gradient[k][i], difference, 2e-9);
    }
}

#[test]
fn g1_anisotropic_heat_and_log_times_match_independent_hold_solution() {
    let strain = STRAINS[0];
    let epsilon = mandel(strain);
    for times in [[1e-12, 0.2], [1e18, 0.7]] {
        let law = material(times, false);
        let mut history = law.history();
        history.advance(strain, 0.0, || false).unwrap();
        let mut elapsed = 0.0;
        let tapes: Vec<_> = [0.01, 0.07, 0.23, 0.4]
            .into_iter()
            .map(|dt| {
                elapsed += dt;
                let tape = history.advance_linearized(strain, dt, || false).unwrap();
                assert_eq!(tape.response().work, 0.0);
                tape
            })
            .collect();
        let mut memory_seed = [[0.0; 6]; 2];
        let mut previous_seed = memory_seed;
        let mut gradient = [0.0; 2];
        for tape in tapes.iter().rev() {
            let mut times = [0.0; 2];
            tape.pull_back(
                TensorPronyStepDerivative {
                    dissipated: 1.0,
                    ..Default::default()
                },
                &memory_seed,
                &mut previous_seed,
                &mut times,
                || false,
            )
            .unwrap();
            memory_seed = previous_seed;
            for j in 0..2 {
                gradient[j] += times[j];
            }
        }
        let mut expected_heat = 0.0;
        for j in 0..2 {
            let energy_twice = dot(
                &epsilon,
                &apply(law.branches()[j].stiffness_mandel(), &epsilon),
            );
            let r = elapsed / times[j];
            expected_heat -= 0.5 * energy_twice * det::expm1(-2.0 * r);
            let expected = -energy_twice * r * det::exp(-2.0 * r);
            if expected == 0.0 {
                assert_eq!(gradient[j], 0.0);
            } else {
                assert!(gradient[j] < 0.0);
                assert!(
                    (gradient[j] / expected - 1.0).abs() < 4e-13,
                    "weak physical heat gradient: {} != {expected}",
                    gradient[j]
                );
            }
        }
        close(history.dissipated(), expected_heat, 3e-13);
    }
}

#[test]
fn g1_tensor_ramp_heat_gradient_matches_independent_continuous_integral() {
    let times = [0.05, 0.7];
    let law = material(times, false);
    let mut history = law.history();
    history.advance(STRAINS[0], 0.0, || false).unwrap();
    let dt = 0.08;
    let tape = history
        .advance_linearized(STRAINS[1], dt, || false)
        .unwrap();
    let mut previous_memory = [[0.0; 6]; 2];
    let mut gradient = [0.0; 2];
    tape.pull_back(
        TensorPronyStepDerivative {
            dissipated: 1.0,
            ..Default::default()
        },
        &[[0.0; 6]; 2],
        &mut previous_memory,
        &mut gradient,
        || false,
    )
    .unwrap();
    let initial = mandel(STRAINS[0]);
    let next = mandel(STRAINS[1]);
    let delta: Voigt = core::array::from_fn(|i| next[i] - initial[i]);
    let n = 8192;
    let mut heat = 0.0;
    for j in 0..2 {
        let c = law.branches()[j].stiffness_mandel();
        let r = dt / times[j];
        let mut derivative = 0.0;
        let mut branch_heat = 0.0;
        for i in 0..=n {
            let u = f64::from(i) / f64::from(n);
            let a = det::exp(-r * u);
            let b = -det::expm1(-r * u) / r;
            // Physical branch strain, independently contracted with the
            // original full stiffness tensor, not its work-factor channels.
            let q = core::array::from_fn(|k| a * initial[k] + b * delta[k]);
            let dq = core::array::from_fn(|k| r * u * a * initial[k] + (b - u * a) * delta[k]);
            let c_q = apply(c, &q);
            let power = dot(&q, &c_q);
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            branch_heat += weight * r * power;
            derivative += weight * r * (2.0 * dot(&dq, &c_q) - power);
        }
        heat += branch_heat / (3.0 * f64::from(n));
        close(gradient[j], derivative / (3.0 * f64::from(n)), 4e-11);
    }
    close(tape.response().dissipated, heat, 4e-12);
}

#[test]
fn g0_tensor_derivative_refusals_and_cancellation_publish_nothing() {
    let law = material([0.05, 0.7], true);
    let mut initial = law.history();
    initial.advance(STRAINS[0], TIMES[0], || false).unwrap();
    let mut unchanged = initial.clone();
    assert!(
        unchanged
            .advance_linearized(STRAINS[0], 0.0, || false)
            .is_err()
    );
    assert_eq!(unchanged, initial);
    assert!(unchanged.advance(STRAINS[0], 0.0, || false).is_ok());
    for stop in 1..=5 {
        let mut history = initial.clone();
        let mut polls = 0;
        assert!(matches!(
            history.advance_linearized(STRAINS[1], TIMES[1], || {
                polls += 1;
                polls == stop
            }),
            Err(TensorPronyError::Cancelled)
        ));
        assert_eq!(history, initial);
    }
    let tape = initial
        .advance_linearized(STRAINS[1], TIMES[1], || false)
        .unwrap();
    let zero_memory = [[0.0; 6]; 2];
    let direction = TensorPronyStepDirection {
        previous_strain: [0.0; 6],
        previous_memory: &zero_memory,
        strain: [0.0; 6],
        log_relaxation_times: &[0.0; 2],
    };
    for stop in 1..=4 {
        let mut output = [[123.0; 6]; 2];
        let mut polls = 0;
        assert_eq!(
            tape.push_forward(direction, &mut output, || {
                polls += 1;
                polls == stop
            }),
            Err(TensorPronyError::Cancelled)
        );
        assert_eq!(output, [[123.0; 6]; 2]);
        let mut times = [456.0; 2];
        let mut polls = 0;
        assert_eq!(
            tape.pull_back(seed(0), &zero_memory, &mut output, &mut times, || {
                polls += 1;
                polls == stop
            }),
            Err(TensorPronyError::Cancelled)
        );
        assert_eq!(output, [[123.0; 6]; 2]);
        assert_eq!(times, [456.0; 2]);
    }
    let mut output = [[123.0; 6]; 2];
    for direction in [
        TensorPronyStepDirection {
            log_relaxation_times: &[0.0],
            ..direction
        },
        TensorPronyStepDirection {
            strain: [f64::NAN; 6],
            ..direction
        },
        TensorPronyStepDirection {
            strain: [f64::MAX; 6],
            ..direction
        },
    ] {
        assert!(tape.push_forward(direction, &mut output, || false).is_err());
        assert_eq!(output, [[123.0; 6]; 2]);
    }
    let mut times = [456.0; 2];
    assert!(
        tape.pull_back(
            TensorPronyStepDerivative {
                stress: [f64::MAX; 6],
                ..Default::default()
            },
            &zero_memory,
            &mut output,
            &mut times,
            || false
        )
        .is_err()
    );
    assert_eq!(output, [[123.0; 6]; 2]);
    assert_eq!(times, [456.0; 2]);
}

#[test]
fn g0_zero_duration_and_pure_elastic_derivatives_retain_exact_limits() {
    let law = material([0.05, 0.7], false);
    let mut history = law.history();
    let tape = history
        .advance_linearized(STRAINS[0], 0.0, || false)
        .unwrap();
    let mut memory = [[0.0; 6]; 2];
    let time_only = tape
        .push_forward(
            TensorPronyStepDirection {
                previous_strain: [0.0; 6],
                previous_memory: &[[0.0; 6]; 2],
                strain: [0.0; 6],
                log_relaxation_times: &[1.0, -1.0],
            },
            &mut memory,
            || false,
        )
        .unwrap();
    assert_eq!(time_only, TensorPronyStepDerivative::default());
    assert_eq!(memory, [[0.0; 6]; 2]);
    let elastic = TensorProny::new_mandel(*law.equilibrium_mandel(), vec![], 0.2).unwrap();
    let tape = elastic
        .history()
        .advance_linearized(STRAINS[0], 0.0, || false)
        .unwrap();
    let pushed = tape
        .push_forward(
            TensorPronyStepDirection {
                previous_strain: [0.0; 6],
                previous_memory: &[],
                strain: STRAINS[1],
                log_relaxation_times: &[],
            },
            &mut [],
            || false,
        )
        .unwrap();
    assert_eq!(pushed.dissipated, 0.0);
    close(pushed.work, pushed.stored_energy, 2e-14);
    let pulled = tape
        .pull_back(seed(0), &[], &mut [], &mut [], || false)
        .unwrap();
    close(
        pair(seed(0), pushed),
        dot(&pulled.strain, &STRAINS[1]),
        2e-14,
    );
}
