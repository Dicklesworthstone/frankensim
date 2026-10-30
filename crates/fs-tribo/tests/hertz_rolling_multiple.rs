//! G1/G3: joint surface and rolling work from the original full 5N equations.

use core::f64::consts::PI;

use fs_tribo::sliding_patch::contact::multiple::rolling::{
    HertzRollingContactInput, MultipleHertzRollingConfig, MultipleHertzRollingError,
    PreparedHertzRollingContacts,
};
use fs_tribo::sliding_patch::contact::multiple::{
    HertzContactInput, MultipleHertzConfig, PreparedHertzContacts,
};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchMode, HertzWrenchParameters,
};

fn law(radius: f64, normal: f64, static_mu: f64, kinetic_mu: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: normal,
        contact_radius_m: radius,
        static_mu,
        kinetic_mu,
    })
    .unwrap()
}
fn config() -> MultipleHertzRollingConfig {
    MultipleHertzRollingConfig {
        hertz: MultipleHertzConfig {
            maximum_contacts: 8,
            maximum_matrix_terms: 576,
            maximum_sweeps: 256,
            maximum_evaluations: 65536,
            local: HertzWrenchConfig {
                maximum_static_iterations: 96,
                maximum_iterations: 256,
                maximum_evaluations: 1024,
                force_tolerance_n: 2e-8,
                moment_tolerance_nm: 2e-9,
                displacement_tolerance_m: 2e-12,
                energy_tolerance_j: 2e-9,
                relative_matrix_tolerance: 1e-13,
                maximum_force_n: 30.0,
                maximum_spin_moment_nm: 5.0,
            },
            total_energy_tolerance_j: 4e-9,
        },
        maximum_matrix_terms: 1600,
        rolling_force_tolerance_n: 2e-9,
        total_energy_tolerance_j: 4e-9,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}
fn diagonal(d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        result[i * n + i] = d[i];
    }
    result
}
fn inputs<'a>(
    laws: &'a [HertzWrenchContact],
    free: &[[f64; 5]],
) -> Vec<HertzRollingContactInput<'a>> {
    laws.iter()
        .zip(free)
        .map(|(law, &free_coordinates_m)| HertzRollingContactInput {
            law,
            free_coordinates_m,
        })
        .collect()
}
fn matrix() -> Vec<f64> {
    // Direct strict diagonal dominance establishes SPD independently of the
    // production factor. Every surface/rolling and inter-patch coupling exists.
    let mut c = vec![0.0; 100];
    for i in 0..10 {
        for j in 0..=i {
            let value = if i == j {
                0.06 + 0.003 * i as f64
            } else {
                (if (i + j) % 3 == 0 { -1.0 } else { 1.0 })
                    * (0.001 + 0.0001 * ((3 * i + j) % 5) as f64)
            };
            c[i * 10 + j] = value;
            c[j * 10 + i] = value;
        }
    }
    c
}
fn free_for(c: &[f64], motion: &[[f64; 5]], wrench: &[[f64; 5]]) -> Vec<[f64; 5]> {
    let n = motion.len() * 5;
    let mut result = motion.to_vec();
    for i in 0..n {
        for j in 0..n {
            result[i / 5][i % 5] -= c[i * n + j] * wrench[j / 5][j % 5];
        }
    }
    result
}

// Independent pivoted Gaussian elimination of the ORIGINAL equations. This
// deliberately uses neither Cholesky nor the production symmetric Schur form.
fn dense_solve(mut a: Vec<f64>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))
            .unwrap();
        assert!(a[pivot * n + k].abs() > 1e-20);
        for j in 0..n {
            a.swap(k * n + j, pivot * n + j);
        }
        b.swap(k, pivot);
        for i in k + 1..n {
            let q = a[i * n + k] / a[k * n + k];
            for j in k + 1..n {
                a[i * n + j] -= q * a[k * n + j];
            }
            b[i] -= q * b[k];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            b[i] -= a[i * n + j] * b[j];
        }
        b[i] /= a[i * n + i];
    }
    b
}

#[test]
fn shared_frictionless_surfaces_retain_all_rolling_recoil_and_bulk_work() {
    let laws = [law(0.08, 7.0, 0.0, 0.0), law(0.13, 11.0, 0.0, 0.0)];
    let c = matrix();
    let damping = [5.0, 13.0];
    let free = [[0.3, -0.2, 0.1, 0.4, -0.3], [-0.1, 0.2, -0.2, -0.25, 0.6]];
    let mut full = vec![0.0; 100];
    for i in 0..10 {
        for j in 0..10 {
            full[i * 10 + j] = if i == j { 1.0 } else { 0.0 }
                + if j % 5 >= 3 {
                    c[i * 10 + j] * damping[j / 5]
                } else {
                    0.0
                };
        }
    }
    let reference = dense_solve(full, free.iter().flatten().copied().collect());
    let prepared =
        PreparedHertzRollingContacts::new(2, c.clone(), damping.to_vec(), config(), || false)
            .unwrap();
    assert_eq!(prepared.contact_count(), 2);
    assert_eq!(prepared.scaled_compliance_m_per_n(), c);
    assert_eq!(prepared.rolling_damping_n_per_m(), damping);
    let step = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    let mut expected_heat = 0.0;
    for i in 0..2 {
        assert_eq!(step.hertz.contacts[i].mode, HertzWrenchMode::Frictionless);
        for j in 0..5 {
            close(step.scaled_motion_m[i][j], reference[5 * i + j], 2e-15);
        }
        for j in 3..5 {
            let x = reference[5 * i + j];
            close(step.scaled_acting_wrench_n[i][j], -damping[i] * x, 3e-15);
            expected_heat += damping[i] * x * x;
        }
    }
    assert_eq!(step.surface_heat_j, 0.0);
    close(step.bulk_heat_j, expected_heat, 4e-15);
    close(step.contact_work_j, expected_heat, 4e-15);
    assert!(step.maximum_compatibility_residual_m < 3e-16);
    assert_eq!(
        prepared.advance(&inputs(&laws, &free), || false).unwrap(),
        step
    );
}

#[test]
fn simultaneous_static_traction_and_nonzero_rolling_match_manufactured_full_wrench() {
    let laws = [law(0.1, 10.0, 0.6, 0.4), law(0.2, 8.0, 0.7, 0.5)];
    let damping = [4.0, 9.0];
    let c = matrix();
    let motion = [[0.0, 0.0, 0.0, 0.08, -0.03], [0.0, 0.0, 0.0, -0.04, 0.06]];
    let wrench = [
        [0.4, -0.2, 0.1, -0.32, 0.12],
        [-0.3, 0.1, -0.15, 0.36, -0.54],
    ];
    let free = free_for(&c, &motion, &wrench);
    let step = PreparedHertzRollingContacts::new(2, c, damping.to_vec(), config(), || false)
        .unwrap()
        .advance(&inputs(&laws, &free), || false)
        .unwrap();
    for i in 0..2 {
        assert_eq!(step.hertz.contacts[i].mode, HertzWrenchMode::Sticking);
        for j in 0..5 {
            close(step.scaled_motion_m[i][j], motion[i][j], 2e-16);
            close(step.scaled_acting_wrench_n[i][j], wrench[i][j], 2e-15);
        }
    }
    assert_eq!(step.surface_heat_j, 0.0);
    close(
        step.bulk_heat_j,
        4.0 * (0.08 * 0.08 + 0.03 * 0.03) + 9.0 * (0.04 * 0.04 + 0.06 * 0.06),
        3e-16,
    );
    close(step.work_residual_j, 0.0, 3e-16);
}

#[test]
fn mixed_shear_spin_and_rolling_match_independent_active_full_system_and_permutation() {
    let laws = [law(0.08, 10.0, 0.6, 0.4), law(0.21, 7.0, 0.5, 0.3)];
    let damping = [7.0, 2.5];
    let c = matrix();
    let motion = [[3.0, 0.0, 0.0, 0.04, -0.06], [0.0, 0.0, 0.8, -0.03, 0.02]];
    let spin = 0.3 * 7.0 * 3.0 * PI / 16.0;
    let wrench = [
        [-4.0, 0.0, 0.0, -0.28, 0.42],
        [0.0, 0.0, -spin, 0.075, -0.05],
    ];
    let free = free_for(&c, &motion, &wrench);
    // The manufactured coordinate check is tighter than the default native
    // effort budget propagated through this response. Request corresponding
    // native effort accuracy instead of hiding it in a looser oracle.
    let mut controls = config();
    controls.hertz.local.force_tolerance_n = 1e-10;
    controls.hertz.local.moment_tolerance_nm = 5e-12;
    let prepared =
        PreparedHertzRollingContacts::new(2, c.clone(), damping.to_vec(), controls, || false)
            .unwrap();
    let step = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    for i in 0..2 {
        assert_eq!(step.hertz.contacts[i].mode, HertzWrenchMode::Sliding);
        for j in 0..5 {
            close(step.scaled_motion_m[i][j], motion[i][j], 2e-11);
            close(step.scaled_acting_wrench_n[i][j], wrench[i][j], 3e-9);
        }
    }
    close(step.surface_heat_j, 4.0 * 3.0 + spin * 0.8, 3e-9);
    close(
        step.bulk_heat_j,
        7.0 * (0.04 * 0.04 + 0.06 * 0.06) + 2.5 * (0.03 * 0.03 + 0.02 * 0.02),
        3e-10,
    );
    // Independent original linear system on this analytically known kinetic
    // branch. A sequential rolling-after-surface solve misses its cross terms.
    let mut a = vec![0.0; 100];
    let mut rhs: Vec<_> = free.iter().flatten().copied().collect();
    for i in 0..10 {
        for j in 0..10 {
            a[i * 10 + j] = if i == j { 1.0 } else { 0.0 };
            if j % 5 >= 3 {
                a[i * 10 + j] += c[i * 10 + j] * damping[j / 5];
            } else {
                rhs[i] += c[i * 10 + j] * wrench[j / 5][j % 5];
            }
        }
    }
    let reference = dense_solve(a, rhs);
    for i in 0..10 {
        close(step.scaled_motion_m[i / 5][i % 5], reference[i], 2e-11);
    }
    let order = [5, 6, 7, 8, 9, 0, 1, 2, 3, 4];
    let mut swapped = vec![0.0; 100];
    for i in 0..10 {
        for j in 0..10 {
            swapped[i * 10 + j] = c[order[i] * 10 + order[j]];
        }
    }
    let permuted = PreparedHertzRollingContacts::new(
        2,
        swapped,
        vec![damping[1], damping[0]],
        controls,
        || false,
    )
    .unwrap()
    .advance(&inputs(&[laws[1], laws[0]], &[free[1], free[0]]), || false)
    .unwrap();
    for i in 0..2 {
        for j in 0..5 {
            close(
                step.scaled_motion_m[i][j],
                permuted.scaled_motion_m[1 - i][j],
                2e-11,
            );
        }
    }
    close(step.surface_heat_j, permuted.surface_heat_j, 3e-9);
    close(step.bulk_heat_j, permuted.bulk_heat_j, 3e-10);
}

#[test]
fn surface_wrench_can_induce_rolling_from_zero_released_rolling_speed() {
    let laws = [law(0.1, 10.0, 0.6, 0.4)];
    let mut c = diagonal(&[0.2, 0.3, 0.4, 0.3, 0.4]);
    c[3] = 0.1;
    c[15] = 0.1;
    let free = [[2.0, 0.0, 0.0, 0.0, 0.0]];
    let d = 5.0;
    let x_roll = -0.1 * 4.0 / (1.0 + 0.3 * d);
    let f_roll = -d * x_roll;
    let x_surface = 2.0 - 0.2 * 4.0 + 0.1 * f_roll;
    let step = PreparedHertzRollingContacts::new(1, c, vec![d], config(), || false)
        .unwrap()
        .advance(&inputs(&laws, &free), || false)
        .unwrap();
    close(step.scaled_motion_m[0][3], x_roll, 2e-15);
    close(step.scaled_acting_wrench_n[0][3], f_roll, 2e-15);
    close(step.scaled_motion_m[0][0], x_surface, 2e-15);
    close(step.bulk_heat_j, d * x_roll * x_roll, 3e-16);
    close(step.surface_heat_j, 4.0 * x_surface, 3e-15);
    assert!(step.bulk_heat_j > 0.0);
}

#[test]
fn zero_damping_exact_prescribed_motion_and_zero_load_keep_original_limits() {
    let laws = [law(0.08, 10.0, 0.6, 0.4), law(0.15, 8.0, 0.7, 0.5)];
    let c = matrix();
    let free = [[2.0, -0.1, 0.2, 0.4, -0.3], [-1.0, 0.2, -0.3, -0.1, 0.6]];
    let prepared =
        PreparedHertzRollingContacts::new(2, c.clone(), vec![0.0; 2], config(), || false).unwrap();
    let step = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    let indices = [0, 1, 2, 5, 6, 7];
    let mut hc = vec![0.0; 36];
    for i in 0..6 {
        for j in 0..6 {
            hc[i * 6 + j] = c[indices[i] * 10 + indices[j]];
        }
    }
    let reference = PreparedHertzContacts::new(2, hc, config().hertz, || false)
        .unwrap()
        .advance(
            &[
                HertzContactInput {
                    law: &laws[0],
                    free_displacement_m: [free[0][0], free[0][1]],
                    free_spin_angle_rad: free[0][2] / 0.08,
                },
                HertzContactInput {
                    law: &laws[1],
                    free_displacement_m: [free[1][0], free[1][1]],
                    free_spin_angle_rad: free[1][2] / 0.15,
                },
            ],
            || false,
        )
        .unwrap();
    // The joint owner also spends a native support evaluation at the original
    // full-C endpoint. Physical Hertz values remain exactly the zero-damping
    // original values; numerical work counters intentionally include that audit.
    for i in 0..2 {
        let joint = step.hertz.contacts[i];
        let original = reference.contacts[i];
        assert_eq!(joint.force_n, original.force_n);
        assert_eq!(joint.spin_moment_nm, original.spin_moment_nm);
        assert_eq!(joint.slip_displacement_m, original.slip_displacement_m);
        assert_eq!(joint.spin_angle_rad, original.spin_angle_rad);
        assert_eq!(joint.dissipated_energy_j, original.dissipated_energy_j);
        assert_eq!(joint.mode, original.mode);
    }
    assert_eq!(step.bulk_heat_j, 0.0);
    for rolling in &step.rolling {
        assert_eq!(rolling.scaled_force_n, [0.0; 2]);
        assert_eq!(rolling.dissipated_energy_j, 0.0);
    }
    let prescribed = [[1e-20, 0.0, 0.0, 0.2, -0.3]];
    let p = PreparedHertzRollingContacts::new(1, vec![0.0; 25], vec![2.0], config(), || false)
        .unwrap()
        .advance(&inputs(&laws[..1], &prescribed), || false)
        .unwrap();
    assert_eq!(p.hertz.contacts[0].mode, HertzWrenchMode::Sliding);
    assert_eq!(p.scaled_motion_m[0], prescribed[0]);
    close(p.surface_heat_j, 4e-20, 1e-34);
    close(p.bulk_heat_j, 2.0 * (0.2 * 0.2 + 0.3 * 0.3), 2e-16);
    let released = [law(0.1, 0.0, 0.6, 0.4)];
    let r =
        PreparedHertzRollingContacts::new(1, diagonal(&[0.1; 5]), vec![0.0], config(), || false)
            .unwrap()
            .advance(&inputs(&released, &[[0.2, -0.1, 0.5, 0.4, -0.2]]), || false)
            .unwrap();
    assert_eq!(r.scaled_acting_wrench_n, [[0.0; 5]]);
    assert_eq!(r.surface_heat_j, 0.0);
    assert_eq!(r.bulk_heat_j, 0.0);
}

#[test]
fn dependent_surface_rows_redistribute_static_load_with_shared_rolling_recoil() {
    let laws = [law(0.1, 5.0, 0.6, 0.4), law(0.15, 20.0, 0.6, 0.4)];
    let mut c = diagonal(&[0.1, 0.1, 0.1, 0.3, 0.4, 0.1, 0.1, 0.1, 0.2, 0.3]);
    for i in 0..3 {
        c[i * 10 + i + 5] = 0.1;
        c[(i + 5) * 10 + i] = 0.1;
    }
    c[3 * 10 + 8] = -0.05;
    c[8 * 10 + 3] = -0.05;
    let free = [[1.2, 0.0, 0.0, 0.3, -0.2], [1.2, 0.0, 0.0, -0.1, 0.4]];
    let step = PreparedHertzRollingContacts::new(2, c, vec![2.0, 5.0], config(), || false)
        .unwrap()
        .advance(&inputs(&laws, &free), || false)
        .unwrap();
    assert_eq!(step.hertz.contacts[0].mode, HertzWrenchMode::Sticking);
    assert_eq!(step.hertz.contacts[1].mode, HertzWrenchMode::Sticking);
    close(
        step.hertz.contacts[0].force_n[0] + step.hertz.contacts[1].force_n[0],
        -12.0,
        3e-11,
    );
    assert!(step.hertz.contacts[0].force_n[0].abs() <= 3.0 + 1e-8);
    assert!(step.hertz.contacts[1].force_n[0].abs() <= 12.0 + 1e-8);
    assert_eq!(step.surface_heat_j, 0.0);
    assert!(step.bulk_heat_j > 0.0);
    close(step.contact_work_j, step.bulk_heat_j, 3e-10);
}

#[test]
fn original_matrix_budgets_unresolved_stiffness_cancellation_and_retry() {
    let c = matrix();
    let laws = [law(0.08, 10.0, 0.6, 0.4), law(0.15, 8.0, 0.7, 0.5)];
    let free = [[0.3, -0.1, 0.2, 0.4, -0.3], [-0.2, 0.2, -0.3, -0.1, 0.6]];
    let prep = PreparedHertzRollingContacts::new(2, c.clone(), vec![2.0, 4.0], config(), || false)
        .unwrap();
    let expected = prep.advance(&inputs(&laws, &free), || false).unwrap();
    assert_eq!(
        prep.advance(&inputs(&laws, &free), || true),
        Err(MultipleHertzRollingError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        prep.advance(&inputs(&laws, &free), || {
            polls += 1;
            polls == 25
        }),
        Err(MultipleHertzRollingError::Cancelled)
    );
    assert_eq!(
        prep.advance(&inputs(&laws, &free), || false).unwrap(),
        expected
    );
    assert_eq!(
        PreparedHertzRollingContacts::new(2, c.clone(), vec![2.0, 4.0], config(), || true),
        Err(MultipleHertzRollingError::Cancelled)
    );
    let mut limited = config();
    limited.maximum_matrix_terms = 99;
    assert!(
        PreparedHertzRollingContacts::new(2, c.clone(), vec![2.0, 4.0], limited, || false).is_err()
    );
    assert!(
        PreparedHertzRollingContacts::new(2, c.clone(), vec![2.0], config(), || false).is_err()
    );
    assert!(
        PreparedHertzRollingContacts::new(2, c.clone(), vec![2.0, -1.0], config(), || false)
            .is_err()
    );
    let mut indefinite = c.clone();
    indefinite[1] = 1.0;
    indefinite[10] = 1.0;
    assert!(
        PreparedHertzRollingContacts::new(2, indefinite, vec![2.0, 4.0], config(), || false)
            .is_err()
    );
    let mut nonfinite = free;
    nonfinite[1][4] = f64::NAN;
    assert!(prep.advance(&inputs(&laws, &nonfinite), || false).is_err());
    assert!(
        prep.advance(&inputs(&laws[..1], &free[..1]), || false)
            .is_err()
    );
    // Perfect h/r correlation with enormous finite damping leaves a physically
    // positive Schur diagonal below resolvable subtraction precision. It must
    // refuse, rather than turn rounded zero into a prescribed-motion model.
    let mut correlated = diagonal(&[1.0; 5]);
    correlated[3] = 1.0;
    correlated[15] = 1.0;
    assert!(matches!(
        PreparedHertzRollingContacts::new(1, correlated, vec![1e30], config(), || false),
        Err(MultipleHertzRollingError::Invalid(
            "unresolved eliminated Hertz diagonal"
        ))
    ));
    // A tiny admitted antisymmetry can reverse a very small prescribed slip.
    // Work and displacement tolerances alone cannot detect the O(mu*N) wrong
    // kinetic force. The final original-full-C native support must reject it.
    let mut almost_symmetric = diagonal(&[0.0, 0.0, 0.0, 1.0, 1.0]);
    almost_symmetric[3] = 1e-15;
    almost_symmetric[15] = -1e-15;
    let reversal =
        PreparedHertzRollingContacts::new(1, almost_symmetric, vec![1.0], config(), || false)
            .unwrap()
            .advance(&inputs(&laws[..1], &[[1e-16, 0.0, 0.0, 1.0, 0.0]]), || {
                false
            });
    assert!(matches!(
        reversal,
        Err(MultipleHertzRollingError::Audit { maximum_surface_force_residual_n, .. })
            if maximum_surface_force_residual_n > 7.9
    ));
    assert_eq!(
        prep.advance(&inputs(&laws, &free), || false).unwrap(),
        expected
    );
}
