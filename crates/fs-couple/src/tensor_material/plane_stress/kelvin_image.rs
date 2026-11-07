use super::*;
use fs_material::{Tangent6, visco_tensor::TensorPronyBranch};
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;

const W: [f64; 6] = [
    1.0,
    1.0,
    1.0,
    core::f64::consts::SQRT_2,
    core::f64::consts::SQRT_2,
    core::f64::consts::SQRT_2,
];
const TIMES: [f64; 2] = [0.13, 0.8];

fn gram(rows: &Tangent6) -> Tangent6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..6).map(|k| rows[k][i] * rows[k][j]).sum())
    })
}
fn factors() -> [Tangent6; 2] {
    let mut a = [[0.0; 6]; 6];
    let mut b = [[0.0; 6]; 6];
    a[0] = [2.0, -0.7, 3.0, 0.4, 1.4, -0.8];
    a[3] = [0.3, 1.2, -0.5, 1.0, -0.4, 0.7];
    b[1] = [-0.6, 2.0, 1.5, 0.5, -1.1, 2.0];
    [a, b]
}
fn law(kv: Option<Tangent6>, branches: usize) -> TensorProny {
    let mut eq = [[0.0; 6]; 6];
    for (i, value) in [80.0, 95.0, 120.0, 55.0, 45.0, 65.0]
        .into_iter()
        .enumerate()
    {
        eq[i][i] = value;
    }
    for (i, j, value) in [
        (0, 2, 3.2),
        (1, 2, -2.1),
        (2, 4, 1.7),
        (0, 4, 2.4),
        (1, 5, -0.9),
    ] {
        eq[i][j] = value;
        eq[j][i] = value;
    }
    let branches = factors()
        .into_iter()
        .zip(TIMES)
        .take(branches)
        .map(|(row, tau)| TensorPronyBranch::from_gram_rows(row, tau).unwrap())
        .collect();
    TensorProny::new(eq, kv, branches, (0.0, 1e9)).unwrap()
}
fn viscosity() -> Tangent6 {
    let mut rows = [[0.0; 6]; 6];
    for (i, value) in [0.7, 0.8, 1.1, 0.5, 0.6, 0.9].into_iter().enumerate() {
        rows[i][i] = value;
    }
    rows[0][2] = 0.25;
    rows[0][4] = -0.13;
    rows[2][4] = -0.2;
    rows[2][5] = 0.17;
    rows[1][5] = 0.23;
    gram(&rows)
}
fn tangential_viscosity() -> Tangent6 {
    let mut c = viscosity();
    for i in 0..6 {
        c[2][i] = 0.0;
        c[i][2] = 0.0;
    }
    c
}
fn budget() -> TensorProjectionBudget {
    TensorProjectionBudget {
        max_coordinates: 16,
        max_scalar_terms: 4_000_000,
    }
}
fn unit() -> [[f64; 5]; 5] {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn retained(columns: &[[f64; 5]], q: &[f64]) -> [f64; 6] {
    let e: [f64; 5] = core::array::from_fn(|i| columns.iter().zip(q).map(|(c, q)| c[i] * q).sum());
    [e[0], e[1], 0.0, e[2] / W[3], e[3] / W[4], e[4] / W[5]]
}
fn source_energy_stress(
    law: &TensorProny,
    e: &[f64; 6],
    a: &[f64],
    volume: f64,
) -> (f64, [f64; 6], Vec<f64>) {
    let mut stress: [f64; 6] =
        core::array::from_fn(|i| (0..6).map(|j| law.c_eq()[i][j] * e[j]).sum());
    let mut energy = 0.5 * volume * stress.iter().zip(e).map(|(s, e)| s * e).sum::<f64>();
    let mut residuals = Vec::new();
    for (j, rows) in factors().iter().take(law.branches().len()).enumerate() {
        for (r, row) in rows.iter().enumerate() {
            let value =
                volume.sqrt() * row.iter().zip(e).map(|(b, e)| b * e).sum::<f64>() - a[6 * j + r];
            energy += 0.5 * value * value;
            for i in 0..6 {
                stress[i] += row[i] * value / volume.sqrt();
            }
            residuals.push(value);
        }
    }
    (energy, stress, residuals)
}
fn raw_modulus(law: &TensorProny, omega: f64) -> [[C64; 6]; 6] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let mut c = C64::new(
                law.c_eq()[i][j],
                omega * law.kv_pa_s().map_or(0.0, |d| d[i][j]),
            );
            for b in law.branches() {
                c = c
                    + (C64::new(0.0, omega * b.tau_s) / C64::new(1.0, omega * b.tau_s))
                        .scale(b.c_pa[i][j]);
            }
            c
        })
    })
}
fn image_modulus(image: &PlaneStressRelaxation, i: usize, j: usize, omega: f64) -> C64 {
    let n = image.point().coordinates();
    image
        .coupling_rows()
        .iter()
        .zip(image.relaxation_rates_per_s())
        .fold(
            C64::new(
                image.equilibrium_stiffness()[n * i + j],
                omega * image.feedthrough_viscosity()[n * i + j],
            ),
            |c, (row, rate)| {
                c + (C64::new(0.0, omega) / C64::new(*rate, omega)).scale(row[i] * row[j])
            },
        )
}
fn image_energy(image: &PlaneStressRelaxation, q: &[f64], z: &[f64]) -> f64 {
    let n = q.len();
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| 0.5 * q[i] * image.equilibrium_stiffness()[n * i + j] * q[j])
                .sum::<f64>()
        })
        .sum::<f64>()
        + image
            .coupling_rows()
            .iter()
            .zip(z)
            .map(|(row, z)| 0.5 * (row.iter().zip(q).map(|(c, q)| c * q).sum::<f64>() - z).powi(2))
            .sum::<f64>()
}

#[test]
fn combined_image_matches_raw_tensor_schur_storage_and_original_kelvin_maxwell_power() {
    let gate = CancelGate::new_clock_free();
    let volume = 0.17;
    for (kv, branches) in [
        (viscosity(), 2),
        (tangential_viscosity(), 2),
        (viscosity(), 0),
        (tangential_viscosity(), 0),
    ] {
        let law = law(Some(kv), branches);
        let point = PlaneStressPoint::new_shell(&law, &unit(), volume, budget(), &gate).unwrap();
        assert!(point.clone().relaxation_image(4_000_000, &gate).is_err());
        let normal = point.normal_strain_coordinates();
        let image = point.viscoelastic_image(4_000_000, &gate).unwrap();
        assert_eq!(image.coupling_rows().len(), 6 * branches + normal);
        for omega in [0.0, 1e-4, 0.1, 1.0, 30.0, 1e5] {
            let c = raw_modulus(&law, omega);
            for (i, a) in [0, 1, 3, 4, 5].into_iter().enumerate() {
                for (j, b) in [0, 1, 3, 4, 5].into_iter().enumerate() {
                    let expected =
                        (c[a][b] - c[a][2] * c[2][b] / c[2][2]).scale(volume / W[a] / W[b]);
                    let actual = image_modulus(&image, i, j, omega);
                    near(actual.re, expected.re, 3e-10);
                    near(actual.im, expected.im, 3e-10);
                }
            }
        }
        let q = [0.001, -0.002, 0.0007, 0.0014, -0.0011];
        let v = [0.003, -0.001, 0.002, -0.004, 0.0017];
        let z: Vec<_> = (0..image.coupling_rows().len())
            .map(|i| 0.0002 * (i as f64 - 3.0))
            .collect();
        let a = image.original_memory_at(&q, &z).unwrap();
        if normal > 0 {
            assert!(image.original_memory(&z).is_err());
        }
        let e = image.strain(&q, &z).unwrap();
        let e = core::array::from_fn(|i| e[i] * W[i]);
        let (energy, mut stress, residuals) = source_energy_stress(&law, &e, &a, volume);
        let mut rate = retained(&unit(), &v);
        if normal > 0 {
            rate[2] = (-stress[2] - (0..6).map(|j| kv[2][j] * rate[j]).sum::<f64>()) / kv[2][2];
        }
        let kv_power = volume
            * (0..6)
                .map(|i| (0..6).map(|j| rate[i] * kv[i][j] * rate[j]).sum::<f64>())
                .sum::<f64>();
        let branch_power: f64 = residuals
            .iter()
            .enumerate()
            .map(|(i, z)| z * z / TIMES[i / 6])
            .sum();
        let expected_heat = kv_power + branch_power;
        let actual_heat = (0..5)
            .map(|i| {
                (0..5)
                    .map(|j| v[i] * image.feedthrough_viscosity()[5 * i + j] * v[j])
                    .sum::<f64>()
            })
            .sum::<f64>()
            + image
                .coupling_rows()
                .iter()
                .zip(&z)
                .zip(image.relaxation_rates_per_s())
                .map(|((row, z), rate)| {
                    rate * (row.iter().zip(q).map(|(c, q)| c * q).sum::<f64>() - z).powi(2)
                })
                .sum::<f64>();
        near(image_energy(&image, &q, &z), energy, 2e-13);
        near(actual_heat, expected_heat, 2e-12);
        assert!(expected_heat > 0.0 && actual_heat > 0.0);
        for i in 0..6 {
            stress[i] += (0..6).map(|j| kv[i][j] * rate[j]).sum::<f64>();
        }
        let actual = image.normal_stress_pa_with_velocity(&q, &v, &z).unwrap();
        for (s, i) in actual.into_iter().zip([2, 4, 5]) {
            near(s, stress[i] / W[i], 2e-11);
        }
        near(
            image
                .thickness_normal_stress_pa_with_velocity(&q, &v, &z)
                .unwrap(),
            0.0,
            2e-11,
        );
    }
}

#[test]
fn combined_projection_retains_exact_rates_physical_normal_units_and_feedthrough() {
    let gate = CancelGate::new_clock_free();
    let law = law(Some(viscosity()), 2);
    let unit_image = PlaneStressPoint::new_shell(&law, &unit(), 1.0, budget(), &gate)
        .unwrap()
        .viscoelastic_image(4_000_000, &gate)
        .unwrap();
    let columns = [
        [1.0, -0.2, 0.3, 0.7, -0.4],
        [-0.1, 0.8, -0.5, 0.2, 0.9],
        [0.0; 5],
    ];
    let q = [0.001, -0.002, 0.0];
    let v = [0.003, -0.001, 0.0];
    let q_unit: Vec<_> = (0..5)
        .map(|i| columns.iter().zip(q).map(|(s, q)| s[i] * q).sum())
        .collect();
    let z_unit: Vec<_> = (0..13).map(|i| 0.0002 * (i as f64 - 4.0)).collect();
    let a_unit = unit_image.original_memory_at(&q_unit, &z_unit).unwrap();
    for volume in [0.003_f64, 0.17, 1.7] {
        let image = unit_image
            .project_shell(&columns, volume, budget(), &gate)
            .unwrap();
        assert_eq!(
            image.relaxation_rates_per_s(),
            unit_image.relaxation_rates_per_s()
        );
        let z: Vec<_> = z_unit.iter().map(|z| volume.sqrt() * z).collect();
        let a = image.original_memory_at(&q, &z).unwrap();
        for i in 0..12 {
            near(a[i], volume.sqrt() * a_unit[i], 2e-13);
        }
        near(a[12], a_unit[12], 2e-13);
        assert!(a[12].abs() > 1e-5, "the normal offset must be exercised");
        near(
            image_energy(&image, &q, &z),
            volume * image_energy(&unit_image, &q_unit, &z_unit),
            2e-13,
        );
        for (got, want) in image
            .strain(&q, &z)
            .unwrap()
            .iter()
            .zip(unit_image.strain(&q_unit, &z_unit).unwrap())
        {
            near(*got, want, 2e-13);
        }
        for i in 0..3 {
            for j in 0..3 {
                let expected = volume
                    * (0..5)
                        .map(|a| {
                            (0..5)
                                .map(|b| {
                                    columns[i][a]
                                        * unit_image.feedthrough_viscosity()[5 * a + b]
                                        * columns[j][b]
                                })
                                .sum::<f64>()
                        })
                        .sum::<f64>();
                near(image.feedthrough_viscosity()[3 * i + j], expected, 2e-13);
            }
        }
        assert!(image.coupling_rows().iter().all(|row| row[2] == 0.0));
        assert!(
            image.feedthrough_viscosity()[6..9]
                .iter()
                .all(|x| *x == 0.0)
        );
        near(
            image
                .thickness_normal_stress_pa_with_velocity(&q, &v, &z)
                .unwrap(),
            0.0,
            2e-11,
        );
        let direct = PlaneStressPoint::new_shell(&law, &columns, volume, budget(), &gate)
            .unwrap()
            .viscoelastic_image(4_000_000, &gate)
            .unwrap();
        for omega in [0.01, 3.0, 300.0] {
            for i in 0..3 {
                for j in 0..3 {
                    let got = image_modulus(&image, i, j, omega);
                    let want = image_modulus(&direct, i, j, omega);
                    near(got.re, want.re, 3e-10);
                    near(got.im, want.im, 3e-10);
                }
            }
        }
    }
}

#[test]
fn combined_image_keeps_zero_kelvin_bits_and_explicit_rank_and_capacity_boundaries() {
    let gate = CancelGate::new_clock_free();
    for kv in [None, Some([[0.0; 6]; 6])] {
        let source = law(kv, 2);
        let point = PlaneStressPoint::new_shell(&source, &unit(), 0.17, budget(), &gate).unwrap();
        let pure = point.clone().relaxation_image(4_000_000, &gate).unwrap();
        let full = point.viscoelastic_image(4_000_000, &gate).unwrap();
        let bits = |values: &[f64]| values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(
            bits(pure.equilibrium_stiffness()),
            bits(full.equilibrium_stiffness())
        );
        assert_eq!(
            bits(pure.relaxation_rates_per_s()),
            bits(full.relaxation_rates_per_s())
        );
        for (a, b) in pure.coupling_rows().iter().zip(full.coupling_rows()) {
            assert_eq!(bits(a), bits(b));
        }
        for i in 0..12 {
            let mut z = [0.0; 12];
            z[i] = 1.0;
            assert_eq!(
                bits(&pure.original_memory(&z).unwrap()),
                bits(&full.original_memory(&z).unwrap())
            );
        }
        assert!(full.feedthrough_viscosity().iter().all(|v| *v == 0.0));
    }
    // Rank-one Kelvin loss may leave exactly zero direct damping while its
    // thickness relaxation remains finite, passive and mechanically coupled.
    let mut row = [[0.0; 6]; 6];
    row[0] = [1.0, 0.0, 2.0, 0.0, 0.0, 0.0];
    let source = law(Some(gram(&row)), 0);
    let point = PlaneStressPoint::new_shell(&source, &unit(), 1.0, budget(), &gate).unwrap();
    let image = point.clone().viscoelastic_image(4_000_000, &gate).unwrap();
    assert_eq!(image.coupling_rows().len(), 1);
    assert!(image.feedthrough_viscosity().iter().all(|v| *v == 0.0));
    assert!(image.relaxation_rates_per_s()[0] > 0.0);
    assert!(image.coupling_rows()[0][0].abs() > 1.0);
    assert!(image.original_memory(&[0.0]).is_err());
    assert!(image.original_memory_at(&[0.0; 4], &[0.0]).is_err());
    assert!(point.clone().viscoelastic_image(1, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        point.viscoelastic_image(4_000_000, &cancelled),
        Err(TensorMaterialError::Cancelled)
    ));
    let mut weak = [[0.0; 6]; 6];
    weak[2][2] = 1e-14;
    let image = PlaneStressPoint::new_shell(&law(Some(weak), 0), &unit(), 1.0, budget(), &gate)
        .unwrap()
        .viscoelastic_image(4_000_000, &gate)
        .unwrap();
    assert_eq!(image.coupling_rows().len(), 1);
    assert!(
        image.relaxation_rates_per_s()[0] > 1e15,
        "weak positive viscosity is not replaced by zero"
    );
    let old = PlaneStressPoint::new(
        &law(Some(viscosity()), 0),
        &[[1.0, 0.0, 0.0]],
        1.0,
        budget(),
        &gate,
    )
    .unwrap();
    assert_eq!(old.normal_strain_coordinates(), 3);
    assert!(
        old.viscoelastic_image(4_000_000, &gate).is_err(),
        "multiple coupled normal mobilities remain explicit"
    );
}

fn image_system(image: &PlaneStressRelaxation, inverse: &[f64; 4]) -> PortHamiltonian {
    let mut q = vec![0.0; 16];
    let mut j = vec![0.0; 16];
    let mut r = vec![0.0; 16];
    let mut g = vec![0.0; 8];
    for i in 0..2 {
        for k in 0..2 {
            q[4 * i + k] = image.equilibrium_stiffness()[2 * i + k];
            q[4 * (2 + i) + 2 + k] = inverse[2 * i + k];
            r[4 * (2 + i) + 2 + k] = image.feedthrough_viscosity()[2 * i + k];
        }
        j[4 * i + 2 + i] = 1.0;
        j[4 * (2 + i) + i] = -1.0;
        g[(2 + i) * 2 + i] = 1.0;
    }
    PortHamiltonian::new(
        4,
        2,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(
        image
            .coupling_rows()
            .iter()
            .zip(image.relaxation_rates_per_s())
            .map(|(row, rate)| RelaxationBranch {
                projection: vec![row[0], row[1], 0.0, 0.0],
                stiffness: 1.0,
                relaxation_time_s: rate.recip(),
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn combined_image_and_complete_body_match_independent_physical_dynamics_and_heat() {
    let gate = CancelGate::new_clock_free();
    let kv = viscosity();
    let source = law(Some(kv), 2);
    let volume = 0.17;
    let columns = [[1.0, -0.2, 0.3, 0.7, -0.4], [-0.1, 0.8, -0.5, 0.2, 0.9]];
    let point = PlaneStressPoint::new_shell(&source, &columns, volume, budget(), &gate).unwrap();
    let image = point.clone().viscoelastic_image(4_000_000, &gate).unwrap();
    let mass = [1.3, 0.2, 0.2, 0.85];
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
        TensorMaterialBudget {
            max_state_coordinates: 128,
            max_material_points: 4,
        },
        &gate,
    )
    .unwrap();
    let system = image_system(&image, &inverse);
    let split = split_quadratic(
        &system,
        LinearSplitBudget {
            max_modes: 17,
            max_omega_dt: core::f64::consts::PI,
        },
    )
    .unwrap();
    for (q, initial) in [
        ([0.001, -0.002], InitialMemory::Relaxed),
        ([0.0; 2], InitialMemory::Instantaneous),
    ] {
        let v = [0.003, -0.002];
        let mut history = body.history(&q, &v, initial, &gate).unwrap();
        let z: Vec<_> = image
            .coupling_rows()
            .iter()
            .map(|row| row.iter().zip(q).map(|(c, q)| c * q).sum())
            .collect();
        let mut x = vec![
            q[0],
            q[1],
            mass[0] * v[0] + mass[1] * v[1],
            mass[2] * v[0] + mass[3] * v[1],
        ];
        x.extend_from_slice(&z);
        let mut reference: Vec<_> = q.into_iter().chain(v).chain([0.0; 15]).collect();
        let mut e = retained(&columns, &q);
        e[2] = -(0..6).map(|j| source.c_eq()[2][j] * e[j]).sum::<f64>() / source.c_eq()[2][2];
        reference[16] = e[2];
        for (j, rows) in factors().iter().enumerate() {
            for (r, row) in rows.iter().enumerate() {
                reference[4 + 6 * j + r] =
                    volume.sqrt() * row.iter().zip(e).map(|(b, e)| b * e).sum::<f64>();
            }
        }
        for (got, want) in image
            .original_memory_at(&q, &z)
            .unwrap()
            .iter()
            .zip(&reference[4..17])
        {
            near(*got, *want, 2e-12);
        }
        let mut heat = 0.0;
        let mut work = 0.0;
        for tick in 0..10 {
            let force = if tick < 4 { [0.03, -0.02] } else { [0.0; 2] };
            let rhs = |state: &[f64]| {
                let mut out = vec![0.0; 19];
                out[..2].copy_from_slice(&state[2..4]);
                let mut e = retained(&columns, &state[..2]);
                e[2] = state[16];
                let (_, mut stress, residuals) =
                    source_energy_stress(&source, &e, &state[4..16], volume);
                let mut rate = retained(&columns, &state[2..4]);
                rate[2] = (-stress[2] - (0..6).map(|j| kv[2][j] * rate[j]).sum::<f64>()) / kv[2][2];
                out[16] = rate[2];
                out[18] = volume
                    * (0..6)
                        .map(|i| (0..6).map(|j| rate[i] * kv[i][j] * rate[j]).sum::<f64>())
                        .sum::<f64>();
                for i in 0..6 {
                    stress[i] += (0..6).map(|j| kv[i][j] * rate[j]).sum::<f64>();
                }
                let restoring: [f64; 2] = core::array::from_fn(|i| {
                    let column = retained(
                        &columns,
                        &[
                            if i == 0 { 1.0 } else { 0.0 },
                            if i == 1 { 1.0 } else { 0.0 },
                        ],
                    );
                    volume * column.iter().zip(stress).map(|(c, s)| c * s).sum::<f64>()
                });
                for i in 0..2 {
                    out[2 + i] = inverse[2 * i] * (force[0] - restoring[0])
                        + inverse[2 * i + 1] * (force[1] - restoring[1]);
                }
                for (i, z) in residuals.iter().enumerate() {
                    out[4 + i] = z / TIMES[i / 6];
                    out[18] += z * z / TIMES[i / 6];
                }
                out[17] = force[0] * state[2] + force[1] * state[3];
                out
            };
            let dt = 0.003;
            let h = dt / 64.0;
            for _ in 0..64 {
                let k1 = rhs(&reference);
                let stage = |k: &[f64], t: f64| {
                    reference
                        .iter()
                        .zip(k)
                        .map(|(x, k)| x + t * k)
                        .collect::<Vec<_>>()
                };
                let k2 = rhs(&stage(&k1, 0.5 * h));
                let k3 = rhs(&stage(&k2, 0.5 * h));
                let k4 = rhs(&stage(&k3, h));
                for i in 0..19 {
                    reference[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
                }
            }
            let record = step_render(&system, &split, &x, &force, dt).unwrap();
            x = record.x;
            heat += record.dissipated;
            work += record.supplied;
            history.advance(&force, dt, &gate).unwrap();
            let velocity = [
                inverse[0] * x[2] + inverse[1] * x[3],
                inverse[2] * x[2] + inverse[3] * x[3],
            ];
            let original = image.original_memory_at(&x[..2], &x[4..]).unwrap();
            for (got, want) in x[..2]
                .iter()
                .chain(&velocity)
                .chain(&original)
                .zip(&reference[..17])
            {
                near(*got, *want, 3e-11);
            }
            for (got, want) in history
                .displacement_m()
                .iter()
                .chain(&history.velocity_m_s())
                .chain(&history.state()[4..])
                .zip(&reference[..17])
            {
                near(*got, *want, 3e-11);
            }
            near(heat, reference[18], 3e-13);
            near(work, reference[17], 3e-13);
            near(history.dissipated_heat_j(), reference[18], 3e-13);
            near(history.supplied_work_j(), reference[17], 3e-13);
            near(system.hamiltonian(&x), history.stored_energy_j(), 3e-13);
            near(
                image
                    .thickness_normal_stress_pa_with_velocity(&x[..2], &velocity, &x[4..])
                    .unwrap(),
                0.0,
                2e-11,
            );
        }
        assert!(heat > 0.0 && history.dissipated_heat_j() > 0.0);
    }
}
