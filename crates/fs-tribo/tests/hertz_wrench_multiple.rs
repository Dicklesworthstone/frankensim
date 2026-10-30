//! G1/G3: simultaneous original Hertz surfaces, full recoil and physical work.

use core::f64::consts::PI;

use fs_tribo::sliding_patch::contact::multiple::{
    HertzContactInput, MultipleHertzConfig, MultipleHertzError, PreparedHertzContacts,
};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
    HertzWrenchParameters,
};

fn law(radius: f64, normal: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: normal,
        contact_radius_m: radius,
        static_mu: 0.6,
        kinetic_mu: 0.4,
    })
    .unwrap()
}
fn config() -> MultipleHertzConfig {
    MultipleHertzConfig {
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
            maximum_force_n: 20.0,
            maximum_spin_moment_nm: 5.0,
        },
        total_energy_tolerance_j: 4e-9,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}
fn diagonal(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        c[i * n + i] = values[i];
    }
    c
}
fn product(c: &[f64], x: &[f64]) -> Vec<f64> {
    c.chunks(x.len())
        .map(|row| row.iter().zip(x).map(|(&a, &b)| a * b).sum())
        .collect()
}
fn inputs<'a>(laws: &'a [HertzWrenchContact], free: &[f64]) -> Vec<HertzContactInput<'a>> {
    laws.iter()
        .enumerate()
        .map(|(i, law)| HertzContactInput {
            law,
            free_displacement_m: [free[3 * i], free[3 * i + 1]],
            free_spin_angle_rad: free[3 * i + 2] / law.parameters().contact_radius_m,
        })
        .collect()
}
fn full_matrix() -> Vec<f64> {
    // Symmetric strict diagonal dominance gives an independent SPD fixture,
    // with genuine force/spin and cross-patch response in every coordinate.
    vec![
        0.040, 0.002, -0.001, -0.006, 0.003, 0.002, 0.002, 0.050, 0.003, 0.001, -0.008, -0.002,
        -0.001, 0.003, 0.035, 0.004, -0.002, -0.005, -0.006, 0.001, 0.004, 0.045, 0.003, -0.002,
        0.003, -0.008, -0.002, 0.003, 0.055, 0.001, 0.002, -0.002, -0.005, -0.002, 0.001, 0.030,
    ]
}

#[test]
fn full_static_seed_retains_cross_patch_force_spin_and_zero_heat() {
    let laws = [law(0.25, 10.0), law(0.12, 8.0)];
    let c = full_matrix();
    let requested = [1.1, -0.4, 0.2, -0.7, 0.3, -0.15];
    let free = product(&c, &requested);
    let prepared = PreparedHertzContacts::new(2, c.clone(), config(), || false).unwrap();
    assert_eq!(prepared.contact_count(), 2);
    assert_eq!(prepared.scaled_compliance_m_per_n(), c);
    let result = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    assert_eq!(result.sweeps, 0);
    for (i, step) in result.contacts.iter().enumerate() {
        assert_eq!(step.mode, HertzWrenchMode::Sticking);
        close(step.force_n[0], -requested[3 * i], 2e-14);
        close(step.force_n[1], -requested[3 * i + 1], 2e-14);
        close(
            step.spin_moment_nm,
            -requested[3 * i + 2] * laws[i].parameters().contact_radius_m,
            2e-14,
        );
        assert_eq!(step.slip_displacement_m, [0.0; 2]);
        assert_eq!(step.spin_angle_rad, 0.0);
        assert_eq!(step.dissipated_energy_j, 0.0);
    }
    close(result.work_residual_j, 0.0, 1e-15);
    assert_eq!(result.dissipated_energy_j, 0.0);
    assert_eq!(
        prepared.advance(&inputs(&laws, &free), || false).unwrap(),
        result
    );
}

#[test]
fn simultaneous_sliding_and_static_contact_match_dense_scalar_reference() {
    let laws = [law(0.25, 10.0), law(0.12, 8.0)];
    let mut c = diagonal(&[0.2, 0.3, 0.4, 0.3, 0.2, 0.4]);
    c[3] = -0.1;
    c[18] = -0.1;
    let free = [3.0, 0.0, 0.0, 0.05, 0.0, 0.0];
    let first_force = 4.0;
    let second_force = (free[3] + 0.1 * first_force) / 0.3;
    let first_motion = free[0] - 0.2 * first_force + 0.1 * second_force;
    let result = PreparedHertzContacts::new(2, c.clone(), config(), || false)
        .unwrap()
        .advance(&inputs(&laws, &free), || false)
        .unwrap();
    assert_eq!(result.contacts[0].mode, HertzWrenchMode::Sliding);
    assert_eq!(result.contacts[1].mode, HertzWrenchMode::Sticking);
    close(result.contacts[0].force_n[0], -first_force, 2e-12);
    close(result.contacts[1].force_n[0], -second_force, 2e-12);
    close(
        result.contacts[0].slip_displacement_m[0],
        first_motion,
        2e-12,
    );
    close(result.contacts[1].dissipated_energy_j, 0.0, 0.0);
    close(
        result.dissipated_energy_j,
        first_force * first_motion,
        2e-11,
    );
    close(result.contact_work_j, first_force * first_motion, 2e-11);
    // Permuting both physical blocks gives the same solution in the new order.
    let order = [3, 4, 5, 0, 1, 2];
    let swapped: Vec<_> = order
        .iter()
        .flat_map(|&i| order.iter().map(move |&j| (i, j)))
        .map(|(i, j)| c[i * 6 + j])
        .collect();
    let swapped_free: Vec<_> = order.iter().map(|&i| free[i]).collect();
    let permuted = PreparedHertzContacts::new(2, swapped, config(), || false)
        .unwrap()
        .advance(&inputs(&[laws[1], laws[0]], &swapped_free), || false)
        .unwrap();
    close(
        permuted.contacts[0].force_n[0],
        result.contacts[1].force_n[0],
        2e-11,
    );
    close(
        permuted.contacts[1].force_n[0],
        result.contacts[0].force_n[0],
        2e-11,
    );
    close(
        permuted.dissipated_energy_j,
        result.dissipated_energy_j,
        2e-11,
    );
}

#[test]
fn pure_spin_coupling_uses_each_original_pressure_radius_and_full_work() {
    let laws = [law(0.3, 10.0), law(0.08, 7.0)];
    let mut c = diagonal(&[0.1, 0.1, 0.2, 0.1, 0.1, 0.3]);
    c[2 * 6 + 5] = -0.04;
    c[5 * 6 + 2] = -0.04;
    let free = [0.0, 0.0, 1.5, 0.0, 0.0, 1.2];
    let spin = [4.0 * 3.0 * PI / 16.0, 2.8 * 3.0 * PI / 16.0];
    let actual = [
        free[2] - 0.2 * spin[0] + 0.04 * spin[1],
        free[5] + 0.04 * spin[0] - 0.3 * spin[1],
    ];
    let result = PreparedHertzContacts::new(2, c, config(), || false)
        .unwrap()
        .advance(&inputs(&laws, &free), || false)
        .unwrap();
    for i in 0..2 {
        assert_eq!(result.contacts[i].mode, HertzWrenchMode::Sliding);
        assert_eq!(result.contacts[i].force_n, [0.0; 2]);
        close(
            result.contacts[i].spin_moment_nm,
            -spin[i] * laws[i].parameters().contact_radius_m,
            2e-13,
        );
        close(
            result.contacts[i].spin_angle_rad,
            actual[i] / laws[i].parameters().contact_radius_m,
            2e-12,
        );
        close(
            result.contacts[i].dissipated_energy_j,
            spin[i] * actual[i],
            2e-12,
        );
    }
    close(
        result.dissipated_energy_j,
        spin[0] * actual[0] + spin[1] * actual[1],
        2e-12,
    );
    close(result.work_residual_j, 0.0, 2e-13);
}

// Independent polar midpoint integration of local Coulomb traction and work.
fn spatial(law: HertzWrenchContact, motion: [f64; 3], radial: usize) -> [f64; 4] {
    let p = law.parameters();
    let angular = 4 * radial;
    let dt = (PI / 2.0) / radial as f64;
    let mut result = [0.0; 4];
    for i in 0..radial {
        let (r, z) = ((i as f64 + 0.5) * dt).sin_cos();
        let weight = 3.0 * p.normal_force_n * p.kinetic_mu * r * z * z * dt / angular as f64;
        for j in 0..angular {
            let (s, c) = (2.0 * PI * (j as f64 + 0.5) / angular as f64).sin_cos();
            let x = r * c;
            let y = r * s;
            let v = [motion[0] - motion[2] * y, motion[1] + motion[2] * x];
            let speed = v[0].hypot(v[1]);
            if speed > 0.0 {
                let traction = v.map(|v| weight * v / speed);
                result[0] += traction[0];
                result[1] += traction[1];
                result[2] += x * traction[1] - y * traction[0];
                result[3] += weight * speed;
            }
        }
    }
    result
}

#[test]
fn mixed_six_axis_recoil_satisfies_independent_local_traction_and_heat_integrals() {
    let laws = [law(0.25, 10.0), law(0.12, 8.0)];
    let c = full_matrix();
    let motions = [[0.28, -0.08, 0.18], [-0.22, 0.1, -0.25]];
    let physical = [
        spatial(laws[0], motions[0], 256),
        spatial(laws[1], motions[1], 256),
    ];
    let forces = [
        physical[0][0],
        physical[0][1],
        physical[0][2],
        physical[1][0],
        physical[1][1],
        physical[1][2],
    ];
    let free: Vec<_> = product(&c, &forces)
        .into_iter()
        .enumerate()
        .map(|(i, value)| value + motions[i / 3][i % 3])
        .collect();
    let result = PreparedHertzContacts::new(2, c.clone(), config(), || false)
        .unwrap()
        .advance(&inputs(&laws, &free), || false)
        .unwrap();
    let accepted: Vec<_> = result
        .contacts
        .iter()
        .zip(laws)
        .flat_map(|(step, law)| {
            [
                -step.force_n[0],
                -step.force_n[1],
                -step.spin_moment_nm / law.parameters().contact_radius_m,
            ]
        })
        .collect();
    let recoil = product(&c, &accepted);
    let mut expected_work = 0.0;
    for i in 0..2 {
        let actual = core::array::from_fn(|axis| free[3 * i + axis] - recoil[3 * i + axis]);
        let coarse = spatial(laws[i], actual, 128);
        let fine = spatial(laws[i], actual, 256);
        let step = result.contacts[i];
        assert_eq!(step.mode, HertzWrenchMode::Sliding);
        for axis in 0..3 {
            let error = (fine[axis] - coarse[axis]).abs();
            close(accepted[3 * i + axis], fine[axis], 4.0 * error + 2e-5);
            close(actual[axis], motions[i][axis], 2e-5);
            expected_work += accepted[3 * i + axis] * actual[axis];
        }
        close(
            step.dissipated_energy_j,
            fine[3],
            4.0 * (fine[3] - coarse[3]).abs() + 2e-7,
        );
        assert!(step.force_residual_n <= config().local.force_tolerance_n);
        assert!(step.spin_moment_residual_nm <= config().local.moment_tolerance_nm);
    }
    close(result.contact_work_j, expected_work, 2e-14);
    close(
        result.dissipated_energy_j,
        expected_work,
        config().total_energy_tolerance_j,
    );
}

#[test]
fn exact_zero_response_preserves_tiny_prescribed_slip_and_released_modes() {
    let laws = [law(0.25, 10.0), law(0.12, 0.0)];
    let prepared = PreparedHertzContacts::new(2, vec![0.0; 36], config(), || false).unwrap();
    for free in [[1e-20, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0; 6]] {
        let result = prepared.advance(&inputs(&laws, &free), || false).unwrap();
        for i in 0..2 {
            let single = laws[i]
                .advance(
                    &HertzWrenchInput {
                        free_displacement_m: [free[3 * i], free[3 * i + 1]],
                        free_spin_angle_rad: free[3 * i + 2]
                            / laws[i].parameters().contact_radius_m,
                        scaled_compliance_m_per_n: [[0.0; 3]; 3],
                    },
                    config().local,
                    || false,
                )
                .unwrap();
            assert_eq!(result.contacts[i].mode, single.mode);
            assert_eq!(result.contacts[i].force_n, single.force_n);
            assert_eq!(
                result.contacts[i].dissipated_energy_j,
                single.dissipated_energy_j
            );
        }
    }
    let changed = [law(0.25, 20.0), laws[1]];
    let original = prepared
        .advance(&inputs(&laws, &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]), || false)
        .unwrap();
    let updated = prepared
        .advance(&inputs(&changed, &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]), || false)
        .unwrap();
    close(
        updated.contacts[0].force_n[0],
        2.0 * original.contacts[0].force_n[0],
        0.0,
    );
    // A tiny global static seed must not impose a reaction on an exactly
    // frictionless patch, even when it lies below the numerical force budget.
    let mobile = PreparedHertzContacts::new(2, diagonal(&[1.0; 6]), config(), || false).unwrap();
    let tiny = mobile
        .advance(&inputs(&laws, &[1e-20, 0.0, 0.0, 1e-22, 0.0, 0.0]), || {
            false
        })
        .unwrap();
    assert_eq!(tiny.contacts[1].mode, HertzWrenchMode::Frictionless);
    assert_eq!(tiny.contacts[1].force_n, [0.0; 2]);
    assert_eq!(tiny.contacts[1].spin_moment_nm, 0.0);
    assert_eq!(tiny.contacts[1].slip_displacement_m, [1e-22, 0.0]);
}

#[test]
fn matrix_effort_iteration_cancellation_refusals_preserve_reusable_owner() {
    let laws = [law(0.25, 10.0), law(0.12, 8.0)];
    let c = full_matrix();
    for broken in [
        {
            let mut a = c.clone();
            a[0] = f64::NAN;
            a
        },
        {
            let mut a = c.clone();
            a[1] += 0.01;
            a
        },
        {
            let mut a = c.clone();
            a[0] = -0.1;
            a
        },
        diagonal(&[1.0, 1.0, 1.0, 1.0, 1.0, 0.0]),
    ] {
        assert!(matches!(
            PreparedHertzContacts::new(2, broken, config(), || false),
            Err(MultipleHertzError::Invalid(_))
        ));
    }
    let prepared = PreparedHertzContacts::new(2, c.clone(), config(), || false).unwrap();
    let saved = prepared.clone();
    let free = [0.8, -0.2, 0.25, -0.6, 0.15, -0.3];
    let clean = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.advance(&inputs(&laws, &free), || {
            polls += 1;
            polls == 12
        }),
        Err(MultipleHertzError::Cancelled)
    ));
    assert_eq!(prepared, saved);
    assert_eq!(
        prepared.advance(&inputs(&laws, &free), || false).unwrap(),
        clean
    );
    let limited = MultipleHertzConfig {
        maximum_evaluations: 1,
        ..config()
    };
    let error = PreparedHertzContacts::new(2, c.clone(), limited, || false)
        .unwrap()
        .advance(&inputs(&laws, &free), || false)
        .unwrap_err();
    assert!(matches!(
        error,
        MultipleHertzError::Unresolved { .. } | MultipleHertzError::Contact { .. }
    ));
    let one_sweep = MultipleHertzConfig {
        maximum_sweeps: 1,
        ..config()
    };
    assert!(matches!(
        PreparedHertzContacts::new(2, c.clone(), one_sweep, || false)
            .unwrap()
            .advance(&inputs(&laws, &free), || false),
        Err(MultipleHertzError::Unresolved { .. })
    ));
    let effort = MultipleHertzConfig {
        local: HertzWrenchConfig {
            maximum_force_n: 0.01,
            ..config().local
        },
        ..config()
    };
    assert!(matches!(
        PreparedHertzContacts::new(2, c, effort, || false)
            .unwrap()
            .advance(&inputs(&laws, &free), || false),
        Err(MultipleHertzError::Contact {
            error: HertzWrenchError::EffortBudget { .. },
            ..
        })
    ));
    assert_eq!(
        prepared.advance(&inputs(&laws, &free), || false).unwrap(),
        clean
    );
}
