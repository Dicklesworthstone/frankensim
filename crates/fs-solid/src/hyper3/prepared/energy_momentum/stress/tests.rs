use super::*;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};

fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 3.0,
            c01: 2.0,
            kappa: 40.0,
        },
        4.0,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}; tolerance {tol:e}");
}
fn states() -> ([f64; 9], [f64; 9]) {
    (
        [1.18, 0.12, -0.08, 0.06, 0.92, 0.09, 0.03, -0.04, 1.06],
        [1.04, -0.18, 0.11, 0.23, 1.03, -0.07, -0.08, 0.15, 0.96],
    )
}
// Independent closed-form Mooney-Rivlin potential, using columns and the
// cofactor norm rather than the production AD invariant evaluation.
fn energy(f: &[f64; 9]) -> f64 {
    let c: [[f64; 3]; 3] = core::array::from_fn(|i| core::array::from_fn(|j| f[3 * j + i]));
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let cof = [cross(c[1], c[2]), cross(c[2], c[0]), cross(c[0], c[1])];
    let j = (0..3).map(|i| c[0][i] * cof[0][i]).sum::<f64>();
    let i1 = f.iter().map(|x| x * x).sum::<f64>();
    let i2 = cof.iter().flatten().map(|x| x * x).sum::<f64>();
    3.0 * (i1 * j.powf(-2.0 / 3.0) - 3.0)
        + 2.0 * (i2 * j.powf(-4.0 / 3.0) - 3.0)
        + 20.0 * (j - 1.0).powi(2)
}

#[test]
fn finite_strain_work_and_nodal_force_torque_close() {
    let (a, b) = states();
    let response = evaluate(&law(), &a, &b, false).unwrap();
    let work = (0..9)
        .map(|i| response.piola[i] * (b[i] - a[i]))
        .sum::<f64>();
    near(work, energy(&b) - energy(&a), 2e-13);
    let mid: [f64; 9] = core::array::from_fn(|i| 0.5 * a[i] + 0.5 * b[i]);
    let grads = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let x = [
        [0.0; 3],
        [mid[0], mid[3], mid[6]],
        [mid[1], mid[4], mid[7]],
        [mid[2], mid[5], mid[8]],
    ];
    let mut force = [0.0; 3];
    let mut torque = [0.0; 3];
    for n in 0..4 {
        let f: [f64; 3] = core::array::from_fn(|i| {
            (0..3)
                .map(|j| response.piola[3 * i + j] * grads[n][j] / 6.0)
                .sum()
        });
        for i in 0..3 {
            force[i] += f[i];
            torque[i] += x[n][(i + 1) % 3] * f[(i + 2) % 3] - x[n][(i + 2) % 3] * f[(i + 1) % 3];
        }
    }
    for i in 0..3 {
        near(force[i], 0.0, 2e-15);
        near(torque[i], 0.0, 2e-14);
    }
}

#[test]
fn analytic_endpoint_tangent_matches_directional_difference() {
    let (a, b) = states();
    let material = law();
    let response = evaluate(&material, &a, &b, true).unwrap();
    let tangent = response.tangent_after.unwrap();
    for column in 0..9 {
        let mut high = b;
        let mut low = b;
        let h = 2e-6;
        high[column] += h;
        low[column] -= h;
        let hp = evaluate(&material, &a, &high, false).unwrap().piola;
        let lp = evaluate(&material, &a, &low, false).unwrap().piola;
        for row in 0..9 {
            near(tangent[row][column], (hp[row] - lp[row]) / (2.0 * h), 2e-7);
        }
    }
}

#[test]
fn superposed_rotation_and_endpoint_reversal_preserve_work_and_stress() {
    let (a, b) = states();
    let material = law();
    let forward = evaluate(&material, &a, &b, false).unwrap().piola;
    let reverse = evaluate(&material, &b, &a, false).unwrap().piola;
    for i in 0..9 {
        near(forward[i], reverse[i], 3e-14);
    }
    let q = [0.36, -0.8, 0.48, 0.8, 0.0, -0.6, 0.48, 0.6, 0.64];
    let ra = product(&q, &a);
    let rb = product(&q, &b);
    let rotated = evaluate(&material, &ra, &rb, false).unwrap().piola;
    let expected = product(&q, &forward);
    for i in 0..9 {
        near(rotated[i], expected[i], 2e-13);
    }
}

#[test]
fn continuous_limit_and_original_material_domain_are_retained() {
    let (a, _) = states();
    let material = law();
    let response = evaluate(&material, &a, &a, true).unwrap();
    let physical = material.piola(&a).unwrap();
    let tangent = material.tangent(&a).unwrap();
    for i in 0..9 {
        near(response.piola[i], physical[i], 3e-13);
        for j in 0..9 {
            near(
                response.tangent_after.unwrap()[i][j],
                0.5 * tangent[i][j],
                2e-12,
            );
        }
    }
    let mut inverted = a;
    inverted[0] = -4.0;
    assert!(evaluate(&material, &a, &inverted, false).is_err());
    let huge = [10.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert!(evaluate(&material, &a, &huge, false).is_err());
}

#[test]
fn original_material_increment_survives_an_arbitrary_storage_datum() {
    struct Shifted(Hyperelastic);
    impl HyperelasticMaterial for Shifted {
        fn hyperelastic(&self) -> &Hyperelastic {
            &self.0
        }
        fn energy_density(&self, f: &[f64; 9]) -> f64 {
            1e30 + self.0.energy(f)
        }
        fn energy_density_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
            self.0.energy_increment(before, after)
        }
    }
    let material = Shifted(law());
    let (before, after) = states();
    assert_eq!(
        material.energy_density(&before),
        material.energy_density(&after)
    );
    let response = evaluate(&material, &before, &after, true).unwrap();
    let work = (0..9)
        .map(|i| response.piola[i] * (after[i] - before[i]))
        .sum::<f64>();
    // The independent original Mooney-Rivlin expression has no arbitrary datum.
    near(work, energy(&after) - energy(&before), 2e-13);
    let native = evaluate(&material.0, &before, &after, true).unwrap();
    assert_eq!(response.piola, native.piola);
    assert_eq!(response.tangent_after, native.tangent_after);
    let pressure = 7.5;
    let loaded =
        evaluate_with_uniform_pressure(&material, &before, &after, pressure, true).unwrap();
    let native_loaded =
        evaluate_with_uniform_pressure(&material.0, &before, &after, pressure, true).unwrap();
    assert_eq!(loaded.piola, native_loaded.piola);
    assert_eq!(loaded.tangent_after, native_loaded.tangent_after);
}

#[test]
fn weak_isochoric_work_and_tangent_survive_large_bulk_storage() {
    let material = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 900.0,
            c01: 600.0,
            kappa: 1e12,
        },
        4.0,
    )
    .unwrap();
    let before = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
    for shear in [4e-7, 8e-7] {
        let mut after = before;
        after[1] = shear;
        // J=8 at both endpoints. The independent reduced invariants give
        // Delta I1_bar=Delta I2_bar=shear²/4; bulk does no work at all.
        // Absolute storage is 2.45e13 Pa, while shear work is only 6e-11 Pa
        // at the first rung. These increments are in the resolved DG branch.
        let expected = 375.0 * shear * shear;
        assert_eq!(material.energy(&before), material.energy(&after));
        let response = evaluate(&material, &before, &after, true).unwrap();
        near(response.piola[1] * shear, expected, 2e-14 * expected);
        near(response.tangent_after.unwrap()[1][1], 375.0, 2e-11);
        let reversed = evaluate(&material, &after, &before, false).unwrap();
        assert_eq!(response.piola, reversed.piola);

        // Symmetric algorithmic S retains zero midpoint torque. Check the
        // returned f64 stress through compensated products, with an allowance
        // for rounding its opposing ~1e6 Pa shear/normal contributions.
        let midpoint: [Dd; 9] = core::array::from_fn(|i| {
            Dd::from_f64(0.5) * (Dd::from_f64(before[i]) + Dd::from_f64(after[i]))
        });
        let spatial = product_dd(
            &response.piola.map(Dd::from_f64),
            &core::array::from_fn(|ij| midpoint[3 * (ij % 3) + ij / 3]),
        );
        for (i, j) in [(1, 3), (2, 6), (5, 7)] {
            near((spatial[i] - spatial[j]).to_f64(), 0.0, 3e-9);
        }
    }
}

fn pressure_states() -> ([f64; 9], [f64; 9]) {
    let before = [1.30, 0.18, -0.08, 0.06, 0.87, 0.12, 0.04, -0.09, 1.08];
    let increment = [0.16, -0.23, 0.09, 0.15, -0.08, 0.21, -0.17, 0.11, 0.07];
    (before, core::array::from_fn(|i| before[i] + increment[i]))
}
fn nodes(f: &[f64; 9]) -> [[f64; 3]; 4] {
    [
        [0.0; 3],
        [f[0], f[3], f[6]],
        [f[1], f[4], f[7]],
        [f[2], f[5], f[8]],
    ]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn volume(f: &[f64; 9]) -> f64 {
    let x = nodes(f);
    let c = cross(x[2], x[3]);
    (0..3).map(|i| x[1][i] * c[i]).sum::<f64>() / 6.0
}
// Independently integrate native P1 face tractions on outward oriented faces.
fn face_force(f: &[f64; 9], p: f64) -> [[f64; 3]; 4] {
    let x = nodes(f);
    let mut force = [[0.0; 3]; 4];
    for [a, b, c] in [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]] {
        let area = cross(
            core::array::from_fn(|i| x[b][i] - x[a][i]),
            core::array::from_fn(|i| x[c][i] - x[a][i]),
        );
        for n in [a, b, c] {
            for i in 0..3 {
                force[n][i] -= p * area[i] / 6.0;
            }
        }
    }
    force
}
fn nodal_torque(x: [[f64; 3]; 4], f: [[f64; 3]; 4]) -> [f64; 3] {
    core::array::from_fn(|i| (0..4).map(|n| cross(x[n], f[n])[i]).sum())
}

#[test]
fn uniform_pressure_closes_volume_work_and_midpoint_torque_together() {
    let (a, b) = pressure_states();
    let material = law();
    let p = 3.0;
    let mid = core::array::from_fn(|i| 0.5 * a[i] + 0.5 * b[i]);
    let original = evaluate(&material, &a, &b, false).unwrap().piola;
    let combined = evaluate_with_uniform_pressure(&material, &a, &b, p, false)
        .unwrap()
        .piola;
    let grads = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let force: [[f64; 3]; 4] = core::array::from_fn(|n| {
        core::array::from_fn(|i| {
            -(0..3)
                .map(|j| (combined[3 * i + j] - original[3 * i + j]) * grads[n][j] / 6.0)
                .sum::<f64>()
        })
    });
    let x0 = nodes(&a);
    let x1 = nodes(&b);
    let work = (0..4)
        .flat_map(|n| (0..3).map(move |i| force[n][i] * (x1[n][i] - x0[n][i])))
        .sum::<f64>();
    near(work, -p * (volume(&b) - volume(&a)), 5e-14);
    for i in 0..3 {
        near(force.iter().map(|f| f[i]).sum(), 0.0, 3e-15);
        near(nodal_torque(nodes(&mid), force)[i], 0.0, 3e-14);
    }
    let energy_work = (0..9)
        .map(|i| combined[i] * (b[i] - a[i]) / 6.0)
        .sum::<f64>();
    near(
        energy_work,
        (energy(&b) - energy(&a)) / 6.0 + p * (volume(&b) - volume(&a)),
        5e-14,
    );
    // Exact straight-segment work does not imply a midpoint angular balance.
    let (f0, fm, f1) = (face_force(&a, p), face_force(&mid, p), face_force(&b, p));
    let simpson = core::array::from_fn(|n| {
        core::array::from_fn(|i| (f0[n][i] + 4.0 * fm[n][i] + f1[n][i]) / 6.0)
    });
    let simpson_work = (0..4)
        .flat_map(|n| (0..3).map(move |i| simpson[n][i] * (x1[n][i] - x0[n][i])))
        .sum::<f64>();
    near(simpson_work, work, 5e-14);
    assert!(
        nodal_torque(nodes(&mid), simpson)
            .iter()
            .any(|v| v.abs() > 1e-3)
    );
    let limit = evaluate_with_uniform_pressure(&material, &a, &a, p, false)
        .unwrap()
        .piola;
    let physical = material.piola(&a).unwrap();
    let endpoint_force: [[f64; 3]; 4] = core::array::from_fn(|n| {
        core::array::from_fn(|i| {
            -(0..3)
                .map(|j| (limit[3 * i + j] - physical[3 * i + j]) * grads[n][j] / 6.0)
                .sum::<f64>()
        })
    });
    for n in 0..4 {
        for i in 0..3 {
            near(endpoint_force[n][i], f0[n][i], 3e-14);
        }
    }
}

#[test]
fn signed_pressure_tangent_objectivity_and_exact_zero_delegation() {
    let (a, b) = pressure_states();
    let material = law();
    let original = evaluate(&material, &a, &b, true).unwrap();
    let zero = evaluate_with_uniform_pressure(&material, &a, &b, 0.0, true).unwrap();
    assert_eq!(zero.piola, original.piola);
    assert_eq!(zero.tangent_after, original.tangent_after);
    for p in [-2.0, 3.0] {
        let response = evaluate_with_uniform_pressure(&material, &a, &b, p, true).unwrap();
        for column in 0..9 {
            let (mut high, mut low) = (b, b);
            let h = 2e-6;
            high[column] += h;
            low[column] -= h;
            let hp = evaluate_with_uniform_pressure(&material, &a, &high, p, false)
                .unwrap()
                .piola;
            let lp = evaluate_with_uniform_pressure(&material, &a, &low, p, false)
                .unwrap()
                .piola;
            for row in 0..9 {
                near(
                    response.tangent_after.unwrap()[row][column],
                    (hp[row] - lp[row]) / (2.0 * h),
                    3e-7,
                );
            }
        }
        let q = [0.36, -0.8, 0.48, 0.8, 0.0, -0.6, 0.48, 0.6, 0.64];
        let rotated =
            evaluate_with_uniform_pressure(&material, &product(&q, &a), &product(&q, &b), p, false)
                .unwrap()
                .piola;
        let expected = product(&q, &response.piola);
        let reverse = evaluate_with_uniform_pressure(&material, &b, &a, p, false)
            .unwrap()
            .piola;
        for i in 0..9 {
            near(rotated[i], expected[i], 2e-13);
            near(reverse[i], response.piola[i], 5e-14);
        }
    }
    assert!(evaluate_with_uniform_pressure(&material, &a, &b, f64::NAN, false).is_err());
}

#[test]
fn uniform_pressure_preserves_resolved_equal_volume_shear_work() {
    let material = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 500.0,
            c01: 0.0,
            kappa: 1e6,
        },
        4.0,
    )
    .unwrap();
    let before: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    for shear in [2e-7_f64, 5e-7, 1e-6, 1e-5] {
        let mut after = before;
        after[1] = shear;
        // J0=J1=1 exactly. Pressure does no work, and the independent
        // Mooney-Rivlin shear potential is W=c10*gamma². These increments
        // enter the resolved correction branch, not its midpoint limit.
        let delta_c_squared = 2.0 * shear * shear + shear.powi(4);
        let mean_c_squared = 3.0 + 1.5 * shear * shear + 0.25 * shear.powi(4);
        assert!(delta_c_squared > 64.0 * f64::EPSILON * mean_c_squared);
        let expected = 500.0 * shear * shear;
        for pressure in [-1e5_f64, 1e5] {
            let response =
                evaluate_with_uniform_pressure(&material, &before, &after, pressure, false)
                    .unwrap();
            near(response.piola[1] * shear, expected, expected * 5e-10);
        }
    }
}

#[test]
fn uniform_pressure_volume_work_coefficient_and_endpoint_tangent_agree() {
    let material = law();
    let before: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let after: [f64; 9] = [
        17.0 / 16.0,
        1.0 / 8.0,
        -1.0 / 16.0,
        0.0,
        31.0 / 32.0,
        1.0 / 32.0,
        0.0,
        0.0,
        65.0 / 64.0,
    ];
    let delta = core::array::from_fn::<_, 9, _>(|i| after[i] - before[i]);
    // The triangular dyadic endpoints give this independent exact volume
    // increment; neither the production determinant helper nor its cofactors
    // are needed by the oracle.
    let volume_change = 1487.0 / 32768.0;
    let material_change = energy(&after) - energy(&before);
    let coefficient = volume_stress(&material, &before, &after).unwrap();
    near(dot(&coefficient, &delta), volume_change, 2e-14);
    let original = evaluate(&material, &before, &after, false).unwrap();
    for pressure in [-1e5_f64, 1e5] {
        let response =
            evaluate_with_uniform_pressure(&material, &before, &after, pressure, true).unwrap();
        near(
            dot(&response.piola, &delta),
            pressure.mul_add(volume_change, material_change),
            5e-10,
        );
        for i in 0..9 {
            near(
                response.piola[i],
                pressure.mul_add(coefficient[i], original.piola[i]),
                5e-9,
            );
        }
        let tangent = response.tangent_after.unwrap();
        let scale = tangent
            .iter()
            .flatten()
            .map(|v| v.abs())
            .fold(1.0, f64::max);
        for column in 0..9 {
            let (mut high, mut low) = (after, after);
            let step = 2e-6;
            high[column] += step;
            low[column] -= step;
            let high_force =
                evaluate_with_uniform_pressure(&material, &before, &high, pressure, false)
                    .unwrap()
                    .piola;
            let low_force =
                evaluate_with_uniform_pressure(&material, &before, &low, pressure, false)
                    .unwrap()
                    .piola;
            for row in 0..9 {
                near(
                    tangent[row][column],
                    (high_force[row] - low_force[row]) / (2.0 * step),
                    5e-8 * scale,
                );
            }
        }
    }
}
