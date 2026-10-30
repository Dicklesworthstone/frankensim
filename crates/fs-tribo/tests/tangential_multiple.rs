//! G1/G3: simultaneous retained contacts against independent global equations.
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::tangential::multiple::{
    MultipleTangentialConfig, MultipleTangentialError, PreparedTangentialContacts,
    TangentialContactInput,
};
use fs_tribo::{
    FrictionLaw, FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
};

fn controls() -> MultipleTangentialConfig {
    MultipleTangentialConfig {
        maximum_contacts: 8,
        maximum_matrix_terms: 256,
        maximum_sweeps: 256,
        maximum_law_evaluations: 65536,
        local: MatrixTangentialConfig {
            maximum_evaluations: 512,
            displacement_tolerance_m: 2e-13,
            energy_tolerance_j: 2e-12,
            matrix_relative_tolerance: 1e-12,
        },
        total_energy_tolerance_j: 4e-12,
    }
}
fn input(free: [f64; 2]) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 1.0,
        slip_increment_m: free,
        dt_s: 0.02,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 1000.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-12,
        },
    }
}
fn matrix() -> Vec<f64> {
    vec![
        0.003, 0.0002, 0.001, 0.0004, 0.0002, 0.005, -0.0003, 0.001, 0.001, -0.0003, 0.004, 0.0001,
        0.0004, 0.001, 0.0001, 0.002,
    ]
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * actual.abs().max(expected.abs()).max(1.0),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}
fn close_slice(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        close(a, b);
    }
}
fn product(matrix: &[f64], vector: &[f64]) -> Vec<f64> {
    matrix
        .chunks(vector.len())
        .map(|row| row.iter().zip(vector).map(|(a, b)| a * b).sum())
        .collect()
}
// Independent dense Gaussian elimination; production uses local constitutive
// solves and block sweeps, never this complete contact equilibrium solve.
fn solve(mut matrix: Vec<f64>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| matrix[i * n + k].abs().total_cmp(&matrix[j * n + k].abs()))
            .unwrap();
        for j in 0..n {
            matrix.swap(k * n + j, pivot * n + j);
        }
        rhs.swap(k, pivot);
        let d = matrix[k * n + k];
        assert!(d.abs() > 1e-18);
        for i in k + 1..n {
            let a = matrix[i * n + k] / d;
            for j in k..n {
                matrix[i * n + j] -= a * matrix[k * n + j];
            }
            rhs[i] -= a * rhs[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            rhs[i] -= matrix[i * n + j] * rhs[j];
        }
        rhs[i] /= matrix[i * n + i];
    }
    rhs
}
fn prepare(matrix: Vec<f64>) -> PreparedTangentialContacts {
    PreparedTangentialContacts::new(2, matrix, controls(), || false).unwrap()
}
fn sites<'a, L: TangentialLaw>(
    laws: &'a [L],
    states: &'a [L::State],
    free: &[f64],
) -> Vec<TangentialContactInput<'a, L>> {
    laws.iter()
        .zip(states)
        .enumerate()
        .map(|(i, (law, state))| TangentialContactInput {
            law,
            state,
            input: input([free[2 * i], free[2 * i + 1]]),
        })
        .collect()
}
fn layer(k: f64, y: f64, resistance: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(k, y, resistance, 300.0, 10.0).unwrap()
}

#[test]
fn simultaneous_elastic_preloads_match_an_independent_complete_equilibrium_solve() {
    let laws = [layer(300.0, 100.0, 0.8), layer(700.0, 100.0, 1.2)];
    let elastic = [0.02, -0.01, -0.015, 0.025];
    let states: Vec<_> = laws
        .iter()
        .enumerate()
        .map(|(i, law)| {
            law.advance(
                &law.initial_state(),
                &input([elastic[2 * i], elastic[2 * i + 1]]),
            )
            .unwrap()
            .state
        })
        .collect();
    let saved = states.clone();
    let c = matrix();
    let prepared = prepare(c.clone());
    let free = [0.001, -0.002, 0.003, 0.0005];
    let result = prepared
        .advance(&sites(&laws, &states, &free), || false)
        .unwrap();
    let mut a = c.clone();
    let mut b = free.to_vec();
    for i in 0..4 {
        a[4 * i + i] += 2.0 / laws[i / 2].stiffness_n_m().unwrap();
        b[i] += 2.0 * elastic[i];
    }
    let expected = solve(a.clone(), b.clone());
    let force: Vec<_> = result
        .contacts
        .iter()
        .flat_map(|step| step.traction_n)
        .collect();
    close_slice(&force, &expected);
    let recoil = product(&c, &force);
    let mut expected_change = 0.0;
    for i in 0..2 {
        let e = result.contacts[i].state.elastic_slip_m();
        for axis in 0..2 {
            close(
                e[axis],
                elastic[2 * i + axis] + free[2 * i + axis] - recoil[2 * i + axis],
            );
        }
        expected_change += 0.5
            * laws[i].stiffness_n_m().unwrap()
            * (e[0] * e[0] + e[1] * e[1] - elastic[2 * i].powi(2) - elastic[2 * i + 1].powi(2));
    }
    close(result.stored_change_j, expected_change);
    assert_eq!(result.dissipated_j, 0.0);
    assert!(result.sweeps > 1);
    assert!(result.body_work_residual_j.abs() <= controls().total_energy_tolerance_j);
    // Removing cross-contact recoil is a measurably different physical answer.
    for i in 0..4 {
        for j in 0..4 {
            if i / 2 != j / 2 {
                a[4 * i + j] = 0.0;
            }
        }
    }
    let independent = solve(a, b);
    assert!(
        force
            .iter()
            .zip(independent)
            .any(|(&a, b)| (a - b).abs() > 0.1)
    );
    assert_eq!(states, saved);
    assert_eq!(prepared.clone().into_compliance_m_per_n(), c);
}

#[test]
fn global_static_priority_and_mixed_stick_slide_use_all_cross_blocks() {
    let laws: [_; 2] = std::array::from_fn(|_| FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: 0.5,
    });
    let states: Vec<_> = laws.iter().map(TangentialLaw::initial_state).collect();
    let c = vec![
        0.5, 0.0, 0.4, 0.0, 0.0, 0.2, 0.0, 0.05, 0.4, 0.0, 0.5, 0.0, 0.0, 0.05, 0.0, 0.2,
    ];
    // Each isolated contact would slide: .72/.5 > mu_s. Together both stick.
    let prepared = prepare(c);
    let free = [0.72, 0.0, 0.72, 0.0];
    let result = prepared
        .advance(&sites(&laws, &states, &free), || false)
        .unwrap();
    for contact in &result.contacts {
        close_slice(&contact.traction_n, &[0.8, 0.0]);
        assert_eq!(contact.regime, FrictionRegime::Sticking);
        assert_eq!(contact.slip_increment_m, [0.0; 2]);
    }
    assert_eq!(result.sweeps, 1);
    assert_eq!(result.law_evaluations, 2);
    assert_eq!(result.dissipated_j, 0.0);
    let mut limited = sites(&laws, &states, &free);
    limited[0].input.solve.maximum_force_n = 0.7;
    assert!(matches!(
        prepared.advance(&limited, || false),
        Err(MultipleTangentialError::Contact { index: 0, .. })
    ));
    let c = vec![
        0.003, 0.0, 0.001, 0.0, 0.0, 0.002, 0.0, 0.0005, 0.001, 0.0, 0.004, 0.0, 0.0, 0.0005, 0.0,
        0.003,
    ];
    let result = prepare(c)
        .advance(&sites(&laws, &states, &[0.0005, 0.0, 0.02, 0.0]), || false)
        .unwrap();
    assert_eq!(result.contacts[0].regime, FrictionRegime::Sticking);
    assert_eq!(result.contacts[1].regime, FrictionRegime::Sliding);
    close_slice(&result.contacts[0].traction_n, &[0.0, 0.0]);
    close_slice(&result.contacts[1].traction_n, &[0.5, 0.0]);
    close(result.dissipated_j, 0.5 * 0.018);
}

// Independent global Newton solve of the convex force-compliance equation:
// (C + diag(2/k)) F + flow*(1 - Y/|F|)_+ F = du_free + 2 e_old.
fn nonlinear_reference(
    c: &[f64],
    free: &[f64],
    elastic: &[f64],
    laws: &[ThirdBodyLayerContact],
) -> Vec<f64> {
    let n = free.len();
    let mut a = c.to_vec();
    let b: Vec<_> = free.iter().zip(elastic).map(|(u, e)| u + 2.0 * e).collect();
    for i in 0..n {
        a[n * i + i] += 2.0 / laws[i / 2].stiffness_n_m().unwrap();
    }
    let mut force = vec![0.0; n];
    for _ in 0..64 {
        let mut residual = product(&a, &force);
        let mut jacobian = a.clone();
        for i in 0..n {
            residual[i] -= b[i];
        }
        for (site, law) in laws.iter().enumerate() {
            let row = 2 * site;
            let length = force[row].hypot(force[row + 1]);
            if length > law.yield_force_n() {
                let flow = 0.02 / law.plastic_resistance_n_s_m();
                let ratio = law.yield_force_n() / length;
                let direction = [force[row] / length, force[row + 1] / length];
                for i in 0..2 {
                    residual[row + i] += flow * (1.0 - ratio) * force[row + i];
                    for j in 0..2 {
                        jacobian[(row + i) * n + row + j] += flow
                            * (if i == j { 1.0 - ratio } else { 0.0 }
                                + ratio * direction[i] * direction[j]);
                    }
                }
            }
        }
        if residual.iter().all(|r| r.abs() < 1e-14) {
            return force;
        }
        let correction = solve(jacobian, residual);
        for i in 0..n {
            force[i] -= correction[i];
        }
    }
    panic!("independent global layer reference did not converge");
}

#[test]
fn multisite_nonproportional_layer_history_matches_independent_global_newton() {
    let laws = [layer(300.0, 0.7, 0.8), layer(500.0, 1.2, 1.4)];
    let c = matrix();
    let prepared = prepare(c.clone());
    let mut states: Vec<_> = laws.iter().map(TangentialLaw::initial_state).collect();
    let mut elastic = vec![0.0; 4];
    let mut plastic = vec![0.0; 4];
    let mut total_heat = 0.0;
    for step in 0..48 {
        let phase = 0.31 * f64::from(step);
        let free = [
            0.025 * phase.cos(),
            0.018 * phase.sin(),
            -0.017 * (1.7 * phase).sin(),
            0.023 * (0.8 * phase).cos(),
        ];
        let expected = nonlinear_reference(&c, &free, &elastic, &laws);
        let result = prepared
            .advance(&sites(&laws, &states, &free), || false)
            .unwrap();
        let force: Vec<_> = result.contacts.iter().flat_map(|s| s.traction_n).collect();
        close_slice(&force, &expected);
        let recoil = product(&c, &expected);
        let mut heat = 0.0;
        for i in 0..2 {
            let length = expected[2 * i].hypot(expected[2 * i + 1]);
            let flow = 0.02 / laws[i].plastic_resistance_n_s_m();
            let factor = if length > laws[i].yield_force_n() {
                flow * (1.0 - laws[i].yield_force_n() / length)
            } else {
                0.0
            };
            let distance = factor * length;
            heat += laws[i].yield_force_n() * distance + distance * distance / flow;
            for axis in 0..2 {
                let j = 2 * i + axis;
                let dp = factor * expected[j];
                elastic[j] += free[j] - recoil[j] - dp;
                plastic[j] += dp;
            }
            close_slice(
                &result.contacts[i].state.elastic_slip_m(),
                &elastic[2 * i..2 * i + 2],
            );
            close_slice(
                &result.contacts[i].state.plastic_slip_m(),
                &plastic[2 * i..2 * i + 2],
            );
            assert!(
                result.contact_work_residuals_j[i].abs() <= controls().local.energy_tolerance_j
            );
        }
        close(result.dissipated_j, heat);
        total_heat += result.dissipated_j;
        states = result.contacts.into_iter().map(|s| s.state).collect();
    }
    assert!(total_heat > 0.1);
}

#[test]
fn independent_frame_rotations_and_contact_permutation_preserve_physical_history() {
    let laws = [layer(300.0, 0.7, 0.8), layer(500.0, 1.2, 1.4)];
    let reordered = [laws[1].clone(), laws[0].clone()];
    let c = matrix();
    let mut transform = vec![0.0; 16];
    for (row, (old, angle)) in [(1, 0.73_f64), (0, -0.41_f64)].into_iter().enumerate() {
        let (s, c) = angle.sin_cos();
        transform[(2 * row) * 4 + 2 * old] = c;
        transform[(2 * row) * 4 + 2 * old + 1] = -s;
        transform[(2 * row + 1) * 4 + 2 * old] = s;
        transform[(2 * row + 1) * 4 + 2 * old + 1] = c;
    }
    let mut transported = vec![0.0; 16];
    for i in 0..4 {
        for j in 0..4 {
            for a in 0..4 {
                for b in 0..4 {
                    transported[4 * i + j] +=
                        transform[4 * i + a] * c[4 * a + b] * transform[4 * j + b];
                }
            }
        }
    }
    let left = prepare(c);
    let right = prepare(transported);
    let mut original: Vec<_> = laws.iter().map(TangentialLaw::initial_state).collect();
    let mut rotated: Vec<_> = reordered.iter().map(TangentialLaw::initial_state).collect();
    for free in [
        [0.02, -0.01, 0.005, 0.03],
        [-0.03, 0.01, 0.02, -0.02],
        [0.01, 0.04, -0.03, 0.0],
    ] {
        let a = left
            .advance(&sites(&laws, &original, &free), || false)
            .unwrap();
        let b = right
            .advance(
                &sites(&reordered, &rotated, &product(&transform, &free)),
                || false,
            )
            .unwrap();
        let forces: Vec<_> = a.contacts.iter().flat_map(|s| s.traction_n).collect();
        let actual: Vec<_> = b.contacts.iter().flat_map(|s| s.traction_n).collect();
        close_slice(&actual, &product(&transform, &forces));
        let elastic: Vec<_> = a
            .contacts
            .iter()
            .flat_map(|s| s.state.elastic_slip_m())
            .collect();
        let actual: Vec<_> = b
            .contacts
            .iter()
            .flat_map(|s| s.state.elastic_slip_m())
            .collect();
        close_slice(&actual, &product(&transform, &elastic));
        close(a.dissipated_j, b.dissipated_j);
        close(a.stored_change_j, b.stored_change_j);
        original = a.contacts.into_iter().map(|s| s.state).collect();
        rotated = b.contacts.into_iter().map(|s| s.state).collect();
    }
}

#[test]
fn redundant_contact_rows_use_minimum_norm_static_seed_and_allow_incompatible_drives() {
    let laws: [_; 2] = std::array::from_fn(|_| FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: 0.5,
    });
    let states: Vec<_> = laws.iter().map(TangentialLaw::initial_state).collect();
    let c = vec![
        0.003, 0.0, 0.003, 0.0, 0.0, 0.003, 0.0, 0.003, 0.003, 0.0, 0.003, 0.0, 0.0, 0.003, 0.0,
        0.003,
    ];
    let prepared = prepare(c);
    let still = prepared
        .advance(
            &sites(&laws, &states, &[0.0012, -0.0006, 0.0012, -0.0006]),
            || false,
        )
        .unwrap();
    for contact in &still.contacts {
        close_slice(&contact.traction_n, &[0.2, -0.1]);
        assert_eq!(contact.regime, FrictionRegime::Sticking);
    }
    assert_eq!(still.dissipated_j, 0.0);
    let moving = prepared
        .advance(&sites(&laws, &states, &[0.01, 0.0, -0.01, 0.0]), || false)
        .unwrap();
    close_slice(&moving.contacts[0].traction_n, &[0.5, 0.0]);
    close_slice(&moving.contacts[1].traction_n, &[-0.5, 0.0]);
    close_slice(&moving.contacts[0].slip_increment_m, &[0.01, 0.0]);
    close_slice(&moving.contacts[1].slip_increment_m, &[-0.01, 0.0]);
    close(moving.dissipated_j, 0.01);
}

#[test]
fn full_matrix_admission_and_collective_refusals_leave_every_original_state_retryable() {
    let indefinite = vec![
        1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 1.0,
    ];
    assert!(PreparedTangentialContacts::new(2, indefinite, controls(), || false).is_err());
    assert!(PreparedTangentialContacts::new(0, vec![], controls(), || false).is_err());
    assert!(PreparedTangentialContacts::new(2, vec![0.0; 15], controls(), || false).is_err());
    let mut bad = matrix();
    bad[2] *= 2.0;
    assert!(PreparedTangentialContacts::new(2, bad, controls(), || false).is_err());
    let mut bad = matrix();
    bad[0] = f64::NAN;
    assert!(PreparedTangentialContacts::new(2, bad, controls(), || false).is_err());
    let laws = [layer(300.0, 0.7, 0.8), layer(500.0, 1.2, 1.4)];
    let states: Vec<_> = laws
        .iter()
        .map(|law| {
            law.advance(&law.initial_state(), &input([0.01, -0.02]))
                .unwrap()
                .state
        })
        .collect();
    let before = states.clone();
    let prepared = prepare(matrix());
    let free = [0.02, 0.01, -0.01, 0.015];
    let short = PreparedTangentialContacts::new(
        2,
        matrix(),
        MultipleTangentialConfig {
            maximum_sweeps: 1,
            ..controls()
        },
        || false,
    )
    .unwrap();
    assert!(matches!(
        short.advance(&sites(&laws, &states, &free), || false),
        Err(MultipleTangentialError::Unresolved { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        prepared.advance(&sites(&laws, &states, &free), || {
            polls += 1;
            polls == 12
        }),
        Err(MultipleTangentialError::Cancelled)
    ));
    let mut wrong = sites(&laws, &states, &free);
    wrong[1].input.contact_temperature_k = Some(301.0);
    assert!(matches!(
        prepared.advance(&wrong, || false),
        Err(MultipleTangentialError::Contact { index: 1, .. })
    ));
    wrong[1].input = input([0.0; 2]);
    wrong[1].input.dt_s = 0.01;
    assert!(matches!(
        prepared.advance(&wrong, || false),
        Err(MultipleTangentialError::Invalid(_))
    ));
    assert_eq!(states, before);
    assert_eq!(
        prepared
            .advance(&sites(&laws, &states, &free), || false)
            .unwrap(),
        prepared
            .advance(&sites(&laws, &before, &free), || false)
            .unwrap()
    );
    let dry: [_; 2] = std::array::from_fn(|_| FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: 0.5,
    });
    let dry_states: Vec<_> = dry.iter().map(TangentialLaw::initial_state).collect();
    for budget in [1, 2, 3] {
        let c = PreparedTangentialContacts::new(
            2,
            matrix(),
            MultipleTangentialConfig {
                maximum_law_evaluations: budget,
                ..controls()
            },
            || false,
        )
        .unwrap();
        let result = c.advance(
            &sites(&dry, &dry_states, &[0.0001, 0.0, -0.0001, 0.0]),
            || false,
        );
        if budget < 2 {
            assert!(result.is_err());
        } else {
            assert_eq!(result.unwrap().law_evaluations, 2);
        }
    }
}
