//! G1: independent uniform-disk limits and manufactured full-matrix response.
//! G3: physical reorder/tangent covariance. G4: cancellation/refusal rollback.

use core::f64::consts::PI;

use fs_tribo::adhesive::spinning::contact::multiple::{
    AdhesiveContactInput, MultipleAdhesiveConfig, MultipleAdhesiveError, PreparedAdhesiveContacts,
};
use fs_tribo::adhesive::spinning::contact::{
    AdhesiveWrenchConfig, AdhesiveWrenchContact, AdhesiveWrenchMode, AdhesiveWrenchParameters,
    AdhesiveWrenchPhase,
};

// Manufacture explicit SI shear strengths from desired tau*A capacities. No
// normal load or pressure distribution enters these adhesive fixtures.
fn law(radius: f64, static_capacity: f64, kinetic_capacity: f64) -> AdhesiveWrenchContact {
    let area = PI * radius * radius;
    AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
        contact_radius_m: radius,
        static_shear_strength_pa: static_capacity / area,
        kinetic_shear_strength_pa: kinetic_capacity / area,
    })
    .unwrap()
}

fn config() -> MultipleAdhesiveConfig {
    MultipleAdhesiveConfig {
        maximum_contacts: 8,
        maximum_matrix_terms: 576,
        maximum_sweeps: 256,
        maximum_evaluations: 65536,
        local: AdhesiveWrenchConfig {
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

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}

fn product(matrix: &[f64], vector: &[f64]) -> Vec<f64> {
    matrix
        .chunks(vector.len())
        .map(|row| row.iter().zip(vector).map(|(a, b)| a * b).sum())
        .collect()
}

fn inputs<'a>(laws: &'a [AdhesiveWrenchContact], free: &[f64]) -> Vec<AdhesiveContactInput<'a>> {
    laws.iter()
        .enumerate()
        .map(|(i, law)| AdhesiveContactInput {
            law,
            free_displacement_m: [free[3 * i], free[3 * i + 1]],
            free_spin_angle_rad: free[3 * i + 2] / law.parameters().contact_radius_m,
        })
        .collect()
}

fn diagonal(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        result[i * n + i] = values[i];
    }
    result
}

fn duplicate_response(response: f64) -> Vec<f64> {
    let mut result = vec![0.0; 36];
    for i in 0..2 {
        for j in 0..2 {
            for axis in 0..3 {
                result[(3 * i + axis) * 6 + 3 * j + axis] = response;
            }
        }
    }
    result
}

fn full_response() -> Vec<f64> {
    // Positive diagonal plus two PSD outer products: every three-axis block
    // has actual shear/spin and cross-junction coupling, without a solver-made
    // symmetry, regularizer or pressure-based approximation.
    let u = [1.0, -0.7, 0.8, 0.3, -0.9, 0.6, -0.4, 0.5, 1.2];
    let v = [0.4, 0.2, -0.3, 0.8, 0.6, -0.5, 0.7, -0.9, 0.1];
    let mut c = vec![0.0; 81];
    for i in 0..9 {
        for j in 0..9 {
            c[i * 9 + j] = 0.002 * u[i] * u[j]
                + 0.001 * v[i] * v[j]
                + if i == j { 0.03 + 0.002 * i as f64 } else { 0.0 };
        }
    }
    c
}

fn manufactured() -> (
    [AdhesiveWrenchContact; 3],
    Vec<f64>,
    Vec<f64>,
    [f64; 9],
    [f64; 9],
) {
    let laws = [
        law(0.30, 3.0, 2.0),
        law(0.08, 2.0, 1.2),
        law(0.11, 1.5, 0.9),
    ];
    let motion = [0.3, 0.4, 0.5, -0.12, 0.16, -0.2, 0.0, 0.0, 0.3];
    // Independent exact disk integrals at |u|=|a*spin| (first two), then
    // pure spin. These constants differ from the Hertz pressure-weighted law.
    let force_ratio = 8.0 / (3.0 * PI);
    let spin_ratio = 8.0 / (9.0 * PI);
    let wrench = [
        2.0 * 0.6 * force_ratio,
        2.0 * 0.8 * force_ratio,
        2.0 * spin_ratio,
        -1.2 * 0.6 * force_ratio,
        1.2 * 0.8 * force_ratio,
        -1.2 * spin_ratio,
        0.0,
        0.0,
        (2.0 / 3.0) * 0.9,
    ];
    let c = full_response();
    let free = product(&c, &wrench)
        .into_iter()
        .zip(motion)
        .map(|(response, slip)| response + slip)
        .collect();
    (laws, c, free, motion, wrench)
}

#[test]
fn three_uniform_disks_match_independent_nine_coordinate_solution_and_gauge_covariance() {
    let (laws, c, free, expected_motion, expected_wrench) = manufactured();
    let prepared = PreparedAdhesiveContacts::new(3, c.clone(), config(), || false).unwrap();
    assert_eq!(prepared.contact_count(), 3);
    assert_eq!(prepared.scaled_compliance_m_per_n(), c);
    assert_eq!(prepared.config(), config());
    let result = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    let mut actual_wrench = [0.0; 9];
    for (i, step) in result.contacts.iter().enumerate() {
        assert_eq!(step.mode, AdhesiveWrenchMode::Sliding);
        actual_wrench[3 * i..3 * i + 3].copy_from_slice(&[
            -step.force_n[0],
            -step.force_n[1],
            -step.spin_moment_nm / laws[i].parameters().contact_radius_m,
        ]);
        let actual_motion = [
            step.slip_displacement_m[0],
            step.slip_displacement_m[1],
            step.spin_angle_rad * laws[i].parameters().contact_radius_m,
        ];
        for axis in 0..3 {
            close(
                actual_wrench[3 * i + axis],
                expected_wrench[3 * i + axis],
                6e-8,
            );
            close(actual_motion[axis], expected_motion[3 * i + axis], 6e-9);
        }
    }
    // Independent work uses the unchanged COMPLETE matrix, including all
    // cross-contact terms. No diagonal-only or per-sweep heat may replace it.
    let recoil = product(&c, &actual_wrench);
    let work: f64 = actual_wrench
        .iter()
        .enumerate()
        .map(|(i, w)| w * (free[i] - recoil[i]))
        .sum();
    let analytic_work: f64 = expected_wrench
        .iter()
        .zip(expected_motion)
        .map(|(w, x)| w * x)
        .sum();
    close(result.contact_work_j, work, 2e-14);
    close(
        result.dissipated_energy_j,
        work,
        config().total_energy_tolerance_j,
    );
    close(result.dissipated_energy_j, analytic_work, 6e-8);

    // Reorder physical junctions and rotate the first/third tangent gauges by
    // +/-90 degrees. Transform C, released motion and dual wrench together.
    let order = [7, 6, 8, 0, 1, 2, 4, 3, 5];
    let signs = [1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0];
    let mut transformed = vec![0.0; 81];
    for i in 0..9 {
        for j in 0..9 {
            transformed[i * 9 + j] = signs[i] * signs[j] * c[order[i] * 9 + order[j]];
        }
    }
    let transformed_free: Vec<_> = (0..9).map(|i| signs[i] * free[order[i]]).collect();
    let reordered = [laws[2], laws[0], laws[1]];
    let transformed_result = PreparedAdhesiveContacts::new(3, transformed, config(), || false)
        .unwrap()
        .advance(&inputs(&reordered, &transformed_free), || false)
        .unwrap();
    for i in 0..3 {
        let step = transformed_result.contacts[i];
        let actual = [
            -step.force_n[0],
            -step.force_n[1],
            -step.spin_moment_nm / reordered[i].parameters().contact_radius_m,
        ];
        for axis in 0..3 {
            close(
                actual[axis],
                signs[3 * i + axis] * actual_wrench[order[3 * i + axis]],
                9e-8,
            );
        }
    }
    close(
        transformed_result.dissipated_energy_j,
        result.dissipated_energy_j,
        8e-8,
    );
}

#[test]
fn dependent_torsion_redistributes_static_junction_capacity_without_false_sliding() {
    let laws = [law(0.25, 3.0, 2.0), law(0.25, 12.0, 8.0)];
    let prepared =
        PreparedAdhesiveContacts::new(2, duplicate_response(0.2), config(), || false).unwrap();
    // Exact rank three. Total scaled torque 8 N is feasible; a minimum-norm
    // 4/4 split exceeds the weak junction's 2 N static spin capacity. Combined
    // kinetic capacity is only 20/3 N, so false kinetic selection would move.
    let free = [0.0, 0.0, 1.6, 0.0, 0.0, 1.6];
    let result = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    assert!(result.sweeps > 0);
    for (i, step) in result.contacts.iter().enumerate() {
        assert_eq!(step.mode, AdhesiveWrenchMode::Sticking);
        assert_eq!(step.force_n, [0.0; 2]);
        assert_eq!(step.spin_angle_rad, 0.0);
        close(step.spin_moment_nm, -[2.0, 6.0][i] * 0.25, 2e-13);
    }
    assert_eq!(result.dissipated_energy_j, 0.0);
    close(result.contact_work_j, 0.0, 2e-13);
    // A dependent static allocation may vary with order; its net physical
    // torque, zero motion and lack of heat must be unchanged.
    let reversed = prepared
        .advance(&inputs(&[laws[1], laws[0]], &free), || false)
        .unwrap();
    close(
        reversed.contacts.iter().map(|s| s.spin_moment_nm).sum(),
        -2.0,
        2e-13,
    );
    assert!(
        reversed
            .contacts
            .iter()
            .all(|s| s.mode == AdhesiveWrenchMode::Sticking)
    );
    assert_eq!(reversed.dissipated_energy_j, 0.0);
}

#[test]
fn mixed_accepted_phases_keep_stationary_strength_and_prevent_moving_capture() {
    let laws = [law(0.25, 3.0, 1.0), law(0.12, 5.0, 0.2)];
    let original = laws;
    let mut c = diagonal(&[0.2, 0.3, 0.4, 0.3, 0.2, 0.4]);
    c[3] = -0.1;
    c[18] = -0.1;
    let prepared = PreparedAdhesiveContacts::new(2, c, config(), || false).unwrap();
    let sites = inputs(&laws, &[0.3, 0.0, 0.0, 0.2, 0.0, 0.0]);
    let rest = prepared.advance(&sites, || false).unwrap();
    assert!(
        rest.contacts
            .iter()
            .all(|s| s.mode == AdhesiveWrenchMode::Sticking)
    );
    close(rest.contacts[0].force_n[0], -2.2, 2e-13);
    close(rest.contacts[1].force_n[0], -1.4, 2e-13);
    let mixed = prepared
        .advance_phases(
            &sites,
            &[
                AdhesiveWrenchPhase::Kinetic,
                AdhesiveWrenchPhase::StaticAdmissible,
            ],
            || false,
        )
        .unwrap();
    assert_eq!(mixed.contacts[0].mode, AdhesiveWrenchMode::Sliding);
    assert_eq!(mixed.contacts[1].mode, AdhesiveWrenchMode::Sticking);
    close(mixed.contacts[0].force_n[0], -1.0, 2e-13);
    close(mixed.contacts[0].slip_displacement_m[0], 0.2, 2e-13);
    close(mixed.contacts[1].force_n[0], -1.0, 2e-13); // greater than its kinetic .2 N
    close(mixed.dissipated_energy_j, 0.2, 2e-12);
    let moving = prepared
        .advance_phases(&sites, &[AdhesiveWrenchPhase::Kinetic; 2], || false)
        .unwrap();
    assert!(
        moving
            .contacts
            .iter()
            .all(|s| s.mode == AdhesiveWrenchMode::Sliding)
    );
    close(moving.contacts[0].slip_displacement_m[0], 0.12, 2e-13);
    close(moving.contacts[1].slip_displacement_m[0], 0.24, 2e-13);
    close(moving.dissipated_energy_j, 0.12 + 0.2 * 0.24, 2e-12);
    assert_eq!(laws, original);
}

#[test]
fn unsupported_opposing_spin_and_zero_strength_rows_cannot_create_stiction() {
    let prepared =
        PreparedAdhesiveContacts::new(2, duplicate_response(0.2), config(), || false).unwrap();
    let laws = [law(0.25, 3.0, 3.0); 2];
    let free = [0.0, 0.0, 0.3, 0.0, 0.0, -0.3];
    let result = prepared.advance(&inputs(&laws, &free), || false).unwrap();
    for (i, step) in result.contacts.iter().enumerate() {
        assert_eq!(step.mode, AdhesiveWrenchMode::Sliding);
        close(step.spin_moment_nm, [-0.5, 0.5][i], 2e-13);
        close(step.spin_angle_rad, [1.2, -1.2][i], 2e-13);
    }
    close(result.dissipated_energy_j, 1.2, 2e-12);
    let zero = law(0.25, 0.0, 0.0);
    let mixed = prepared
        .advance(
            &inputs(&[laws[0], zero], &[0.0, 0.0, 0.2, 0.0, 0.0, 0.6]),
            || false,
        )
        .unwrap();
    assert_eq!(mixed.contacts[0].mode, AdhesiveWrenchMode::Sticking);
    close(mixed.contacts[0].spin_moment_nm, -0.25, 2e-13);
    assert_eq!(mixed.contacts[1].mode, AdhesiveWrenchMode::Frictionless);
    assert_eq!(mixed.contacts[1].spin_moment_nm, 0.0);
    close(mixed.contacts[1].spin_angle_rad, 1.6, 2e-13);
    assert_eq!(mixed.dissipated_energy_j, 0.0);
}

#[test]
fn cancellation_budgets_and_unresolved_static_decisions_leave_the_whole_owner_retryable() {
    let (laws, c, free, _, _) = manufactured();
    let prepared = PreparedAdhesiveContacts::new(3, c.clone(), config(), || false).unwrap();
    let saved = prepared.clone();
    let sites = inputs(&laws, &free);
    let mut polls = 0;
    let reference = prepared
        .advance(&sites, || {
            polls += 1;
            false
        })
        .unwrap();
    for allowance in [0, 3, polls - 1] {
        let mut calls = 0;
        assert_eq!(
            prepared.advance(&sites, || {
                calls += 1;
                calls > allowance
            }),
            Err(MultipleAdhesiveError::Cancelled)
        );
        assert_eq!(prepared, saved);
    }
    assert_eq!(prepared.advance(&sites, || false).unwrap(), reference);
    assert!(PreparedAdhesiveContacts::new(3, c.clone(), config(), || true).is_err());
    assert!(
        prepared
            .advance_phases(&sites, &[AdhesiveWrenchPhase::Kinetic], || false)
            .is_err()
    );
    let mut insufficient = config();
    insufficient.maximum_evaluations = 1;
    assert!(
        PreparedAdhesiveContacts::new(3, c, insufficient, || false)
            .unwrap()
            .advance(&sites, || false)
            .is_err()
    );
    let mut indefinite = diagonal(&[1.0; 6]);
    indefinite[3] = 2.0;
    indefinite[18] = 2.0;
    assert!(PreparedAdhesiveContacts::new(2, indefinite, config(), || false).is_err());

    let laws = [law(0.25, 3.0, 2.0), law(0.25, 12.0, 8.0)];
    let mut controls = config();
    controls.maximum_sweeps = 8;
    let dependent =
        PreparedAdhesiveContacts::new(2, duplicate_response(0.2), controls, || false).unwrap();
    // The physical drive slightly exceeds the 10 N total static spin bound,
    // but the separating work margin is below the requested decision budget.
    // Neither a lower kinetic capacity nor numerical stiction is justified.
    for spin in [[2.0 + 1e-10, 2.0 + 1e-10], [1e-20, -1e-20]] {
        assert!(matches!(
            dependent.advance(
                &inputs(&laws, &[0.0, 0.0, spin[0], 0.0, 0.0, spin[1]]),
                || false
            ),
            Err(MultipleAdhesiveError::StaticFeasibilityUnresolved { .. })
        ));
    }
}
