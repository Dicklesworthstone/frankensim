use super::*;
use fs_material::Tangent6;
use fs_material::visco_tensor::TensorPronyBranch;
use fs_math::c64::C64;

const W: [f64; 6] = [
    1.0,
    1.0,
    1.0,
    core::f64::consts::SQRT_2,
    core::f64::consts::SQRT_2,
    core::f64::consts::SQRT_2,
];
const VOLUME: f64 = 0.17;
const TIMES: [f64; 2] = [0.08, 0.63];

fn gram(rows: &Tangent6) -> Tangent6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..6).map(|k| rows[k][i] * rows[k][j]).sum())
    })
}

fn factors() -> [Tangent6; 2] {
    let mut a = [[0.0; 6]; 6];
    let mut b = [[0.0; 6]; 6];
    a[0] = [3.0, -1.0, 4.0, 0.7, 2.0, -0.8];
    a[2] = [0.2, 1.4, -0.8, 1.0, -0.6, 1.2];
    b[1] = [-1.0, 2.0, 2.5, 1.0, -1.4, 3.0];
    b[4] = [0.5, -0.2, 1.0, -0.4, 1.5, 0.3];
    [a, b]
}

fn law() -> TensorProny {
    let mut eq = [[0.0; 6]; 6];
    for (i, row) in eq.iter_mut().enumerate() {
        row[i] = 7.0 + i as f64;
    }
    eq[0][2] = 1.4;
    eq[1][4] = -0.7;
    eq[2][4] = 1.1;
    eq[2][5] = -0.4;
    TensorProny::new_mandel(
        gram(&eq),
        factors()
            .into_iter()
            .zip(TIMES)
            .map(|(b, tau)| TensorPronyBranch::from_gram_rows(b, tau).unwrap())
            .collect(),
        0.1,
    )
    .unwrap()
}

fn budget() -> TensorProjectionBudget {
    TensorProjectionBudget {
        max_coordinates: 16,
        max_scalar_terms: 4_000_000,
    }
}

fn body_budget() -> TensorMaterialBudget {
    TensorMaterialBudget {
        max_state_coordinates: 128,
        max_material_points: 4,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

fn columns() -> [[f64; 5]; 5] {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
}

fn retained_strain(columns: &[[f64; 5]], u: &[f64]) -> [f64; 6] {
    let eng: [f64; 5] =
        core::array::from_fn(|i| columns.iter().zip(u).map(|(column, u)| column[i] * u).sum());
    [
        eng[0],
        eng[1],
        0.0,
        eng[2] / W[3],
        eng[3] / W[4],
        eng[4] / W[5],
    ]
}

// Independent original 3D constitutive law in Mandel units. This neither reads
// the point's projected/condensed matrices nor uses its relaxation image.
fn complete(law: &TensorProny, strain: &[f64; 6], a: &[f64]) -> ([f64; 6], f64, Vec<f64>) {
    let mut stress: [f64; 6] = core::array::from_fn(|i| {
        (0..6)
            .map(|j| law.equilibrium_mandel()[i][j] * strain[j])
            .sum()
    });
    let mut energy = 0.5 * VOLUME * stress.iter().zip(strain).map(|(s, e)| s * e).sum::<f64>();
    let mut residuals = Vec::new();
    for (branch, rows) in factors().iter().enumerate() {
        for (row, values) in rows.iter().enumerate() {
            let residual = VOLUME.sqrt()
                * values.iter().zip(strain).map(|(b, e)| b * e).sum::<f64>()
                - a[6 * branch + row];
            energy += 0.5 * residual * residual;
            for i in 0..6 {
                stress[i] += values[i] * residual / VOLUME.sqrt();
            }
            residuals.push(residual);
        }
    }
    (stress, energy, residuals)
}

fn constrained_strain(law: &TensorProny, mut strain: [f64; 6], a: &[f64]) -> [f64; 6] {
    let mut instant = *law.equilibrium_mandel();
    for rows in factors() {
        let c = gram(&rows);
        for i in 0..6 {
            for j in 0..6 {
                instant[i][j] += c[i][j];
            }
        }
    }
    let history: f64 = factors()
        .iter()
        .enumerate()
        .map(|(j, rows)| {
            (0..6)
                .map(|r| rows[r][2] * a[6 * j + r] / VOLUME.sqrt())
                .sum::<f64>()
        })
        .sum();
    strain[2] = (history - (0..6).map(|k| instant[2][k] * strain[k]).sum::<f64>()) / instant[2][2];
    strain
}

#[test]
fn shell_keeps_physical_shear_and_complete_source_energy() {
    let gate = CancelGate::new_clock_free();
    let law = law();
    let columns = columns();
    let point = PlaneStressPoint::new_shell(
        &law,
        &columns,
        VOLUME,
        TensorProjectionBudget {
            max_coordinates: 6,
            ..budget()
        },
        &gate,
    )
    .unwrap();
    assert_eq!(
        point.source.coordinates(),
        6,
        "only zz consumes a constraint coordinate"
    );
    assert_eq!(point.normal_map.len(), 1);
    assert_eq!(
        point.memory_coordinates(),
        12,
        "every six-row source history survives"
    );
    let u = [0.001, -0.002, 0.003, -0.0017, 0.0023];
    let a: Vec<_> = (0..12).map(|i| 0.0003 * (i as f64 - 4.0)).collect();
    let strain = constrained_strain(&law, retained_strain(&columns, &u), &a);
    let (stress, energy, _) = complete(&law, &strain, &a);
    for (actual, expected) in point
        .strain(&u, &a)
        .unwrap()
        .iter()
        .zip((0..6).map(|i| strain[i] / W[i]))
    {
        close(*actual, expected, 2e-13);
    }
    let tractions = point.normal_stress_pa(&u, &a).unwrap();
    for (actual, i) in tractions.into_iter().zip([2, 4, 5]) {
        close(actual, stress[i] / W[i], 2e-12);
    }
    close(
        point.thickness_normal_stress_pa(&u, &a).unwrap(),
        0.0,
        2e-12,
    );
    assert!(
        tractions[1].abs() > 0.03 && tractions[2].abs() > 0.03,
        "shell shear cannot be silently constrained to zero: {tractions:?}"
    );
    let x: Vec<_> = u.into_iter().chain(a.iter().copied()).collect();
    let m = x.len();
    let actual_energy: f64 = (0..m)
        .map(|i| {
            (0..m)
                .map(|j| 0.5 * x[i] * point.energy_matrix()[i * m + j] * x[j])
                .sum::<f64>()
        })
        .sum();
    close(actual_energy, energy, 2e-14);
    for i in 0..5 {
        let force: f64 = (0..m)
            .map(|j| point.energy_matrix()[i * m + j] * x[j])
            .sum();
        let axis = [0, 1, 3, 4, 5][i];
        close(force, VOLUME * stress[axis] / W[axis], 2e-13);
    }
    assert!(
        PlaneStressPoint::new_shell(
            &law,
            &columns,
            VOLUME,
            TensorProjectionBudget {
                max_coordinates: 5,
                ..budget()
            },
            &gate
        )
        .is_err()
    );
    assert!(
        PlaneStressPoint::new_shell(
            &law,
            &columns,
            VOLUME,
            TensorProjectionBudget {
                max_scalar_terms: 1,
                ..budget()
            },
            &gate
        )
        .is_err()
    );
    assert!(
        point.strain(&[0.2; 5], &a).is_err(),
        "original 3D strain limit still applies"
    );
}

#[test]
fn shell_relaxation_rates_and_cross_couplings_match_full_single_traction_schur() {
    let gate = CancelGate::new_clock_free();
    let law = law();
    let point = PlaneStressPoint::new_shell(&law, &columns(), VOLUME, budget(), &gate).unwrap();
    let image = point.clone().relaxation_image(4_000_000, &gate).unwrap();
    assert_eq!(image.coupling_rows().len(), 12);
    assert!(
        image
            .relaxation_rates_per_s()
            .iter()
            .any(|rate| TIMES.iter().all(|tau| (rate - tau.recip()).abs() > 1e-4)),
        "the thickness constraint must couple full histories and change their rates"
    );
    for omega in [0.0, 0.01, 0.3, 1.0, 10.0, 1000.0] {
        let (storage, loss) = law.modulus_mandel(omega).unwrap();
        let c: [[C64; 6]; 6] =
            core::array::from_fn(|i| core::array::from_fn(|j| C64::new(storage[i][j], loss[i][j])));
        for (i, a) in [0, 1, 3, 4, 5].into_iter().enumerate() {
            for (j, b) in [0, 1, 3, 4, 5].into_iter().enumerate() {
                // Eliminate sigma_zz of the COMPLETE complex modulus once.
                let expected = (c[a][b] - c[a][2] * c[2][b] / c[2][2]).scale(VOLUME / W[a] / W[b]);
                let actual = image
                    .coupling_rows()
                    .iter()
                    .zip(image.relaxation_rates_per_s())
                    .fold(
                        C64::from_re(image.equilibrium_stiffness()[5 * i + j]),
                        |s, (row, rate)| {
                            s + (C64::new(0.0, omega) / C64::new(*rate, omega))
                                .scale(row[i] * row[j])
                        },
                    );
                close(actual.re, expected.re, 3e-10);
                close(actual.im, expected.im, 3e-10);
            }
        }
    }
    let u = [0.001, -0.002, 0.003, -0.001, 0.002];
    let z: Vec<_> = (0..12).map(|i| 0.0002 * (i as f64 - 3.0)).collect();
    let a = image.original_memory(&z).unwrap();
    let strain = constrained_strain(&law, retained_strain(&columns(), &u), &a);
    let (_, original_energy, _) = complete(&law, &strain, &a);
    let elastic: f64 = (0..5)
        .map(|i| {
            (0..5)
                .map(|j| 0.5 * u[i] * image.equilibrium_stiffness()[5 * i + j] * u[j])
                .sum::<f64>()
        })
        .sum();
    let spring: f64 = image
        .coupling_rows()
        .iter()
        .zip(&z)
        .map(|(row, z)| 0.5 * (row.iter().zip(u).map(|(c, u)| c * u).sum::<f64>() - z).powi(2))
        .sum();
    close(elastic + spring, original_energy, 2e-13);
    close(
        image.thickness_normal_stress_pa(&u, &z).unwrap(),
        0.0,
        2e-12,
    );
    assert!(image.normal_stress_pa(&u, &z).unwrap()[1].abs() > 1e-3);
}

#[test]
fn shell_complete_history_matches_independent_three_dimensional_dynamics_and_heat() {
    let gate = CancelGate::new_clock_free();
    let law = law();
    let columns = [[1.0, -0.2, 0.3, 0.7, -0.4], [-0.1, 0.8, -0.5, 0.2, 0.9]];
    let point = PlaneStressPoint::new_shell(&law, &columns, VOLUME, budget(), &gate).unwrap();
    let mass = [1.4, 0.31, 0.31, 0.9];
    let determinant = mass[0] * mass[3] - mass[1] * mass[2];
    let inverse = [
        mass[3] / determinant,
        -mass[1] / determinant,
        -mass[2] / determinant,
        mass[0] / determinant,
    ];
    let body = TensorMaterialBody::from_plane_stress_mass_matrix(
        2,
        mass.to_vec(),
        vec![point],
        body_budget(),
        &gate,
    )
    .unwrap();
    for initial in [InitialMemory::Relaxed, InitialMemory::Instantaneous] {
        let u = [0.001, -0.002];
        let velocity = [0.003, -0.002];
        let mut history = body.history(&u, &velocity, initial, &gate).unwrap();
        let mut reference: Vec<_> = u.into_iter().chain(velocity).chain([0.0; 14]).collect();
        if initial == InitialMemory::Relaxed {
            let mut strain = retained_strain(&columns, &u);
            strain[2] = -(0..6)
                .map(|i| law.equilibrium_mandel()[2][i] * strain[i])
                .sum::<f64>()
                / law.equilibrium_mandel()[2][2];
            for (j, rows) in factors().iter().enumerate() {
                for (r, row) in rows.iter().enumerate() {
                    reference[4 + 6 * j + r] =
                        VOLUME.sqrt() * row.iter().zip(strain).map(|(b, e)| b * e).sum::<f64>();
                }
            }
        }
        for (actual, expected) in history.state()[4..].iter().zip(&reference[4..16]) {
            close(*actual, *expected, 2e-13);
        }
        for tick in 0..12 {
            let force = if tick < 4 { [0.03, -0.02] } else { [0.0; 2] };
            let derivative = |x: &[f64]| {
                let mut out = vec![0.0; 18];
                out[..2].copy_from_slice(&x[2..4]);
                let strain =
                    constrained_strain(&law, retained_strain(&columns, &x[..2]), &x[4..16]);
                let (stress, _, residuals) = complete(&law, &strain, &x[4..16]);
                let restoring: [f64; 2] = core::array::from_fn(|i| {
                    let column = retained_strain(
                        &columns,
                        &[
                            if i == 0 { 1.0 } else { 0.0 },
                            if i == 1 { 1.0 } else { 0.0 },
                        ],
                    );
                    VOLUME * column.iter().zip(stress).map(|(e, s)| e * s).sum::<f64>()
                });
                for i in 0..2 {
                    out[2 + i] = inverse[2 * i] * (force[0] - restoring[0])
                        + inverse[2 * i + 1] * (force[1] - restoring[1]);
                }
                for (i, residual) in residuals.iter().enumerate() {
                    out[4 + i] = residual / TIMES[i / 6];
                    out[17] += residual * residual / TIMES[i / 6];
                }
                out[16] = force[0] * x[2] + force[1] * x[3];
                out
            };
            let dt = 0.008;
            let h = dt / 64.0;
            for _ in 0..64 {
                let k1 = derivative(&reference);
                let stage = |k: &[f64], t: f64| {
                    reference
                        .iter()
                        .zip(k)
                        .map(|(x, k)| x + t * k)
                        .collect::<Vec<_>>()
                };
                let k2 = derivative(&stage(&k1, 0.5 * h));
                let k3 = derivative(&stage(&k2, 0.5 * h));
                let k4 = derivative(&stage(&k3, h));
                for i in 0..reference.len() {
                    reference[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
                }
            }
            history.advance(&force, dt, &gate).unwrap();
            for (actual, expected) in history
                .displacement_m()
                .iter()
                .chain(&history.velocity_m_s())
                .chain(&history.state()[4..])
                .zip(&reference[..16])
            {
                close(*actual, *expected, 2e-11);
            }
            close(history.supplied_work_j(), reference[16], 2e-13);
            close(history.dissipated_heat_j(), reference[17], 2e-13);
            close(history.balance_residual_j(), 0.0, 2e-13);
            close(
                history.plane_stress_thickness_normal_stress_pa(0).unwrap(),
                0.0,
                2e-12,
            );
        }
        assert!(history.dissipated_heat_j() > 0.0);
    }
}

#[test]
fn shell_keeps_one_kelvin_normal_rate_and_actual_transverse_shear_stress() {
    let gate = CancelGate::new_clock_free();
    let source = law();
    let mut rows = [[0.0; 6]; 6];
    rows[0] = [0.4, -0.2, 1.1, 0.1, 0.7, -0.5];
    rows[1] = [0.0, 0.0, 0.0, 0.0, 0.8, 0.3];
    // Resolve both retained loss directions: the unchanged exact heat owner
    // may refuse roundoff-unresolved singular resistance at propagation.
    rows[2] = [0.6, 0.8, 0.0, -0.2, 0.1, -0.3];
    let kv = gram(&rows);
    let law = TensorProny::new(
        *source.c_eq(),
        Some(kv),
        source.branches().to_vec(),
        (0.0, 1e6),
    )
    .unwrap();
    let columns = [[1.0, -0.2, 0.3, 0.7, -0.4], [-0.1, 0.8, -0.5, 0.2, 0.9]];
    let point = PlaneStressPoint::new_shell(&law, &columns, VOLUME, budget(), &gate).unwrap();
    assert_eq!(point.normal_strain_coordinates(), 1);
    assert_eq!(point.memory_coordinates(), 13);
    assert!(point.clone().relaxation_image(4_000_000, &gate).is_err());
    let body =
        TensorMaterialBody::from_plane_stress(vec![1.0, 1.4], vec![point], body_budget(), &gate)
            .unwrap();
    let point = &body.plane_stress_points()[0];
    assert!(
        body.history(
            &[0.001, -0.002],
            &[0.0; 2],
            InitialMemory::Instantaneous,
            &gate
        )
        .is_err()
    );
    let mut history = body
        .history(
            &[0.001, -0.002],
            &[0.003, -0.002],
            InitialMemory::Relaxed,
            &gate,
        )
        .unwrap();
    for _ in 0..8 {
        history.advance(&[0.03, -0.02], 0.003, &gate).unwrap();
        let a = &history.state()[4..];
        let strain = history.plane_stress_strain(0).unwrap();
        let mandel = core::array::from_fn(|i| W[i] * strain[i]);
        let (mut stress, material_energy, _) = complete(&law, &mandel, &a[..12]);
        let velocity = history.velocity_m_s();
        let mut rate = retained_strain(&columns, &velocity);
        rate[2] = (-stress[2] - (0..6).map(|j| kv[2][j] * rate[j]).sum::<f64>()) / kv[2][2];
        for i in 0..6 {
            stress[i] += (0..6).map(|j| kv[i][j] * rate[j]).sum::<f64>();
        }
        let actual = history.plane_stress_normal_stress_pa(0).unwrap();
        for (s, i) in actual.into_iter().zip([2, 4, 5]) {
            close(s, stress[i] / W[i], 2e-11);
        }
        assert!(actual[1].abs() > 1e-4);
        assert!(point.normal_stress_pa(history.displacement_m(), a).is_err());
        close(
            history.stored_energy_j(),
            material_energy + 0.5 * (velocity[0] * velocity[0] + 1.4 * velocity[1] * velocity[1]),
            2e-13,
        );
        close(history.balance_residual_j(), 0.0, 2e-13);
    }
    assert!(history.dissipated_heat_j() > 0.0);
    // Tangential Kelvin viscosity needs a velocity for actual shear stress,
    // while sigma_zz alone remains an algebraic observation.
    let mut tangential = [[0.0; 6]; 6];
    tangential[4][4] = 0.3;
    let law = TensorProny::new(*source.c_eq(), Some(tangential), Vec::new(), (0.0, 1e6)).unwrap();
    let point = PlaneStressPoint::new_shell(&law, &columns, VOLUME, budget(), &gate).unwrap();
    assert_eq!(point.normal_strain_coordinates(), 0);
    assert!(point.normal_stress_pa(&[0.001, -0.002], &[]).is_err());
    close(
        point
            .thickness_normal_stress_pa(&[0.001, -0.002], &[])
            .unwrap(),
        0.0,
        2e-12,
    );
}

#[test]
fn shell_image_projects_volume_and_strains_without_recomputing_constitutive_rates() {
    let gate = CancelGate::new_clock_free();
    let law = law();
    let unit = PlaneStressPoint::new_shell(&law, &columns(), 1.0, budget(), &gate)
        .unwrap()
        .relaxation_image(4_000_000, &gate)
        .unwrap();
    let columns = [
        [1.0, -0.2, 0.3, 0.7, -0.4],
        [-0.1, 0.8, -0.5, 0.2, 0.9],
        [0.0; 5],
    ];
    let z: Vec<_> = (0..12).map(|i| 0.0002 * (i as f64 - 3.0)).collect();
    for volume in [0.003, VOLUME, 1.7] {
        let image = unit
            .project_shell(&columns, volume, budget(), &gate)
            .unwrap();
        assert_eq!(
            image.relaxation_rates_per_s(),
            unit.relaxation_rates_per_s(),
            "quadrature geometry must not round originally identical constitutive rates differently"
        );
        assert_eq!(
            image.original_memory(&z).unwrap(),
            unit.original_memory(&z).unwrap(),
            "original and image histories both scale by sqrt(volume)"
        );
        assert!(image.coupling_rows().iter().all(|row| row[2] == 0.0));
        assert!(
            image.equilibrium_stiffness()[6..9]
                .iter()
                .all(|v| *v == 0.0)
        );
        for omega in [0.0, 0.3, 10.0, 1000.0] {
            let (storage, loss) = law.modulus_mandel(omega).unwrap();
            let c: [[C64; 6]; 6] = core::array::from_fn(|i| {
                core::array::from_fn(|j| C64::new(storage[i][j], loss[i][j]))
            });
            for i in 0..3 {
                for j in 0..3 {
                    let mut expected = C64::from_re(0.0);
                    for (a, ia) in [0, 1, 3, 4, 5].into_iter().enumerate() {
                        for (b, ib) in [0, 1, 3, 4, 5].into_iter().enumerate() {
                            expected = expected
                                + (c[ia][ib] - c[ia][2] * c[2][ib] / c[2][2])
                                    .scale(volume * columns[i][a] * columns[j][b] / W[ia] / W[ib]);
                        }
                    }
                    let actual = image
                        .coupling_rows()
                        .iter()
                        .zip(image.relaxation_rates_per_s())
                        .fold(
                            C64::from_re(image.equilibrium_stiffness()[3 * i + j]),
                            |s, (row, rate)| {
                                s + (C64::new(0.0, omega) / C64::new(*rate, omega))
                                    .scale(row[i] * row[j])
                            },
                        );
                    close(actual.re, expected.re, 3e-10);
                    close(actual.im, expected.im, 3e-10);
                }
            }
        }
        let u = [0.001, -0.002, 0.003];
        let a = image.original_memory(&z).unwrap();
        let x: Vec<_> = u.into_iter().chain(a.iter().copied()).collect();
        let m = x.len();
        let original: f64 = (0..m)
            .map(|i| {
                (0..m)
                    .map(|j| 0.5 * x[i] * image.point().energy_matrix()[i * m + j] * x[j])
                    .sum::<f64>()
            })
            .sum();
        let projected: f64 = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| 0.5 * u[i] * image.equilibrium_stiffness()[3 * i + j] * u[j])
                    .sum::<f64>()
            })
            .sum::<f64>()
            + image
                .coupling_rows()
                .iter()
                .zip(&z)
                .map(|(row, z)| {
                    0.5 * (row.iter().zip(u).map(|(c, u)| c * u).sum::<f64>() - z).powi(2)
                })
                .sum::<f64>();
        close(original, projected, 2e-13);
        close(
            image.thickness_normal_stress_pa(&u, &z).unwrap(),
            0.0,
            2e-12,
        );
        assert_eq!(
            image.strain(&u, &z).unwrap(),
            image.point().strain(&u, &a).unwrap()
        );
        assert!(
            image
                .project_shell(&columns, volume, budget(), &gate)
                .is_err()
        );
    }
    assert!(
        unit.project_shell(
            &columns,
            VOLUME,
            TensorProjectionBudget {
                max_scalar_terms: 1,
                ..budget()
            },
            &gate
        )
        .is_err()
    );
    assert!(
        unit.project_shell(
            &columns,
            VOLUME,
            TensorProjectionBudget {
                max_coordinates: 3,
                ..budget()
            },
            &gate
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        unit.project_shell(&columns, VOLUME, budget(), &cancelled)
            .is_err()
    );
}
