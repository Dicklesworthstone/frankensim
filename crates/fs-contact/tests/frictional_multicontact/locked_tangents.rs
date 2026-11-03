//! Exact locked tangents retain physical Coulomb traction and work authority.

use super::*;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContact, solve_elliptic_frictional_contacts,
};

fn normal_mobility(normal: &[f64], contacts: usize) -> Vec<f64> {
    let dimension = 3 * contacts;
    let mut matrix = vec![0.0; dimension * dimension];
    for row in 0..contacts {
        for column in 0..contacts {
            matrix[(3 * row) * dimension + 3 * column] = normal[row * contacts + column];
        }
    }
    matrix
}

#[test]
fn locked_sliding_pairs_keep_both_tractions_and_original_relative_work() {
    let matrix = normal_mobility(&[0.5], 1);
    let mut original = contact(7, [-2.0, 3.0, 4.0], 0.25);
    original.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 };
    for angle in [0.0_f64, 0.63] {
        let (sine, cosine) = angle.sin_cos();
        let mut c = original;
        c.relative_velocity[1] = 3.0 * cosine + 4.0 * sine;
        c.relative_velocity[2] = -3.0 * sine + 4.0 * cosine;
        let result = solve(&[c], &matrix);
        close(result.impulses[0][0], 6.0);
        close(result.impulses[0][1], -0.3 * c.relative_velocity[1]);
        close(result.impulses[0][2], -0.3 * c.relative_velocity[2]);
        close(result.relative_velocities_after[0][0], 1.0);
        assert_eq!(
            result.relative_velocities_after[0][1..],
            c.relative_velocity[1..]
        );
        close(result.relative_work, -10.5);
        close(result.dissipated_energy, 10.5);
        let FrictionalContactDeterminacy::UniqueVelocity(witness) = &result.determinacy else {
            panic!("a singular mobility must retain its range authority")
        };
        assert_eq!(witness.dimension, 3);
        assert_eq!(witness.rank, 1);
        close(witness.friction_feedback_bound, 0.0);
        let work = result.impulse_work_witness.as_ref().unwrap();
        assert_eq!(work.locked_tangent_substitutions.len(), 1);
        let substitution = &work.locked_tangent_substitutions[0];
        assert_eq!(substitution.active_contact, 0);
        close(substitution.normal_velocity_shift, -1.25);
        for axis in 0..2 {
            close(
                substitution.weighted_slip_direction[axis],
                0.05 * c.relative_velocity[axis + 1],
            );
        }
        assert!(work.range_velocity_residual <= work.allowed_range_velocity_residual);
        let bound = work.locked_tangent_work_admission.as_ref().unwrap();
        assert!(bound.impulse_norm_bound.is_some());
        assert!(bound.maximum_work_difference <= bound.allowed_work_difference);
    }
}

#[test]
fn every_nonzero_locked_slip_gets_its_finite_coulomb_traction() {
    let matrix = normal_mobility(&[0.5], 1);
    for speed in [1e-12, 1e-200] {
        let result = solve(&[contact(1, [-1.0, speed, -2.0 * speed], 0.5)], &matrix);
        close(result.impulses[0][0], 2.0);
        close(result.impulses[0][1], -1.0 / 5.0_f64.sqrt());
        close(result.impulses[0][2], 2.0 / 5.0_f64.sqrt());
        close(result.impulses[0][1].hypot(result.impulses[0][2]), 1.0);
        assert_eq!(result.relative_velocities_after[0][1], speed);
        assert_eq!(result.relative_velocities_after[0][2], -2.0 * speed);
    }
}

#[test]
fn zero_slip_retains_a_static_load_gauge_and_frictionless_axes_stay_zero() {
    let matrix = normal_mobility(&[0.5], 1);
    let result = solve(&[contact(1, [-2.0, 0.0, 0.0], 0.25)], &matrix);
    close(result.impulses[0][0], 4.0);
    close(result.impulses[0][1], 0.0);
    close(result.impulses[0][2], 0.0);
    close(result.relative_work, -4.0);
    assert!(matches!(
        result.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    let work = result.impulse_work_witness.as_ref().unwrap();
    assert!(work.locked_tangent_substitutions.is_empty());
    assert_eq!(work.impulse_coordinates, [0, 1, 2]);
    // The independent static traction (0.3,0.4) is also inside the unit disk.
    // Its zero-mobility image and zero-slip work are identical, so neither the
    // solver's zero representative nor the local support reaction is unique.
    assert!(0.3_f64.hypot(0.4) < 0.25 * result.impulses[0][0]);
    close(0.5 * result.impulses[0][0], 2.0);

    let frictionless = solve(&[contact(2, [-2.0, 3.0, -5.0], 0.0)], &matrix);
    close(frictionless.impulses[0][0], 4.0);
    assert_eq!(frictionless.impulses[0][1..], [0.0, 0.0]);
    assert_eq!(
        frictionless
            .impulse_work_witness
            .unwrap()
            .impulse_coordinates,
        [0]
    );
    let separating = solve(&[contact(3, [1.0, 0.0, 0.0], 0.7)], &matrix);
    assert_eq!(separating.impulses, [[0.0; 3]]);
    close(separating.relative_work, 0.0);

    let mixed = solve(
        &[
            contact(7, [-2.0, 3.0, 0.0], 0.25),
            contact(8, [-2.0, 10.0, -5.0], 0.0),
        ],
        &normal_mobility(&[0.5, 0.0, 0.0, 0.5], 2),
    );
    close(mixed.dissipated_energy, 11.0);
    assert_eq!(
        mixed
            .impulse_work_witness
            .unwrap()
            .locked_tangent_work_admission
            .unwrap()
            .zero_friction_coordinates,
        [4, 5]
    );
}

fn ambiguous_elastic_pressure_gauge(mobility: f64, mu: f64, slip: f64) -> FrictionalContactError {
    let matrix = normal_mobility(&[mobility; 4], 2);
    let contacts = [
        FrictionalContact {
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 },
            ..contact(1, [-1.0, slip, 0.0], mu)
        },
        FrictionalContact {
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 },
            ..contact(2, [-1.0, 0.0, 0.0], mu)
        },
    ];
    let total = 2.0 / mobility;
    // Independent exact family: pn=(a,total-a), pt1=-mu*a, pt2=0.
    // Both endpoint normal velocities are +1, tangents stay fixed, both
    // Coulomb laws hold, and every member is passive. Its loss ranges from
    // zero to this finite value; a selected representative cannot identify it.
    close(-1.0 + mobility * total, 1.0);
    let maximum_loss = mu * slip * total;
    assert!(maximum_loss > 1e-4);
    let mut config = FrictionalContactConfig::default();
    config.newton.relative_tolerance = 1e-18;
    solve_frictional_contacts(&contacts, &matrix, config, || false)
        .expect_err("a finite same-motion loss gauge must not publish arbitrary work")
}

#[test]
fn tiny_raw_slip_does_not_hide_finite_loss_with_large_friction_and_pressure() {
    // Raw slip is smaller than the original absolute velocity allowance.
    // Nevertheless the same-velocity family has losses ranging from 0 to .01 J.
    let error = ambiguous_elastic_pressure_gauge(1e-6, 10.0, 5e-10);
    assert!(
        matches!(error, FrictionalContactError::UnresolvedWork { .. }),
        "{error:?}"
    );
}

#[test]
fn small_work_covector_projection_is_not_a_small_energy_gauge() {
    // The adjusted null velocity is below 1e-9, but loss varies by .001 J.
    let error = ambiguous_elastic_pressure_gauge(1e-6, 0.1, 5e-9);
    assert!(
        matches!(error, FrictionalContactError::UnresolvedWork { .. }),
        "{error:?}"
    );
}

#[test]
fn a_sub_ulp_friction_coefficient_cannot_disappear_beside_normal_speed() {
    let mu = 0.1;
    let slip = 1e-16;
    assert_eq!(-1.0 - mu * slip, -1.0);
    match ambiguous_elastic_pressure_gauge(2.0_f64.powi(-50), mu, slip) {
        FrictionalContactError::UnresolvedWork { .. } => {}
        FrictionalContactError::ResidualNotSatisfied(residual) => {
            // At this enormous pressure the numerical natural map can stop
            // at a rounded point before reaching work admission. Preserve
            // the original physical refusal; an exact-root private test
            // separately exercises the real DD work helper at this scale.
            assert!(residual.normal_complementarity.is_finite());
            assert!(
                residual.normal_complementarity
                    > FrictionalContactConfig::default().complementarity_tolerance
            );
        }
        other => panic!("unexpected refusal of sub-ULP loss gauge: {other:?}"),
    }
}

#[test]
fn coupled_normal_mobility_is_retained_with_distinct_fixed_slips() {
    let matrix = normal_mobility(&[1.0, -0.5, -0.5, 1.0], 2);
    let contacts = [
        contact(20, [-1.0, 3.0, 0.0], 0.2),
        contact(10, [-0.5, -4.0, 0.0], 0.3),
    ];
    for reverse in [false, true] {
        let mut input = contacts;
        if reverse {
            input.reverse();
        }
        let result = solve(&input, &matrix);
        for (index, c) in input.iter().enumerate() {
            let expected = if c.id == 20 {
                [5.0 / 3.0, -1.0 / 3.0, 0.0]
            } else {
                [4.0 / 3.0, 0.4, 0.0]
            };
            for axis in 0..3 {
                close(result.impulses[index][axis], expected[axis]);
            }
            close(result.relative_velocities_after[index][0], 0.0);
            assert_eq!(
                result.relative_velocities_after[index][1..],
                c.relative_velocity[1..]
            );
        }
        close(result.dissipated_energy, 113.0 / 30.0);
    }
}

#[test]
fn redundant_normals_require_compatible_total_sliding_work() {
    let matrix = normal_mobility(&[1.0; 4], 2);
    let contacts = [
        contact(8, [-1.0, 3.0, 0.0], 0.2),
        contact(2, [-1.0, 1.5, 0.0], 0.4),
    ];
    for input in [contacts, [contacts[1], contacts[0]]] {
        let result = solve(&input, &matrix);
        close(result.impulses.iter().map(|p| p[0]).sum(), 1.0);
        close(result.dissipated_energy, 1.1);
        for (index, c) in input.iter().enumerate() {
            close(
                result.impulses[index][1],
                -c.friction_coefficient * result.impulses[index][0],
            );
        }
        let work = result.impulse_work_witness.unwrap();
        assert_eq!(work.locked_tangent_substitutions.len(), 2);
        for substitution in &work.locked_tangent_substitutions {
            close(substitution.normal_velocity_shift, -0.6);
        }
    }
    let incompatible = [
        contacts[0],
        FrictionalContact {
            relative_velocity: [-1.0, 2.0, 0.0],
            ..contacts[1]
        },
    ];
    assert!(matches!(
        solve_frictional_contacts(
            &incompatible,
            &matrix,
            FrictionalContactConfig::default(),
            || false,
        ),
        Err(FrictionalContactError::UnresolvedWork { .. })
    ));

    // Different normal restitution losses can be offset by their known
    // locked-slip losses. Both rows demand the same mobility image, and the
    // COMPLETE work covector, rather than normal work alone, is determined.
    let compensated = [
        FrictionalContact {
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 },
            ..contact(1, [-1.0, 3.0, 0.0], 0.4)
        },
        contact(2, [-2.0, 1.0, 0.0], 0.2),
    ];
    let result = solve(&compensated, &matrix);
    close(result.impulses.iter().map(|p| p[0]).sum(), 2.0);
    close(result.dissipated_energy, 2.4);
}

#[test]
fn normal_and_nonstructural_locks_are_not_regularized_or_finite_mobility_dropped() {
    let c = contact(1, [-1.0, 0.0, 0.0], 0.5);
    let matrices = [
        [0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5],
        [0.5, 1e-14, 0.0, 1e-14, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    for matrix in matrices {
        let result =
            solve_frictional_contacts(&[c], &matrix, FrictionalContactConfig::default(), || false);
        assert!(
            matches!(
                result,
                Err(FrictionalContactError::MobilityNotAdmitted(_)
                    | FrictionalContactError::NonFinite)
            ),
            "{result:?}"
        );
    }
    // Preserve the original range owner's admitted partial-axis static case.
    // The new full-pair slip identity is not needed and does not replace it.
    let partial = [0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5];
    let result = solve(&[c], &partial);
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = &result.determinacy else {
        panic!("lost original partial-axis range admission")
    };
    assert_eq!(witness.rank, 2);
    close(result.dissipated_energy, 1.0);
    assert!(
        result
            .impulse_work_witness
            .unwrap()
            .locked_tangent_substitutions
            .is_empty()
    );
    // Genuine finite tangent mobility, however small, keeps the original SPD
    // path and its unique-impulse authority; it is never declared a lock.
    let weak = [0.5, 0.0, 0.0, 0.0, 1e-24, 0.0, 0.0, 0.0, 1e-24];
    let result = solve(&[c], &weak);
    assert_eq!(
        result.determinacy,
        FrictionalContactDeterminacy::UniqueImpulses
    );
    assert!(result.minimum_normalized_pivot > 0.0);
    assert!(result.impulse_work_witness.is_none());
}

#[test]
fn locked_directional_slip_keeps_the_original_physical_ellipse() {
    let matrix = normal_mobility(&[0.5], 1);
    let c = EllipticFrictionalContact {
        friction: EllipticFriction {
            coefficients_mu: [0.2, 0.5],
            principal_axis: [1.0, 0.0],
        },
        ..contact(1, [-2.0, 3.0, 4.0], 0.0).into()
    };
    let result = solve_elliptic_frictional_contacts(
        &[c],
        &matrix,
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap();
    let speed = 0.6_f64.hypot(2.0);
    close(result.impulses[0][0], 4.0);
    close(result.impulses[0][1], -0.48 / speed);
    close(result.impulses[0][2], -4.0 / speed);
    close(result.dissipated_energy, 4.0 + 4.0 * speed);
    assert!(matches!(
        result.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    assert_eq!(
        result
            .transformed
            .impulse_work_witness
            .unwrap()
            .locked_tangent_substitutions
            .len(),
        1
    );
}

#[test]
fn locked_pairs_keep_cancellation_solver_caps_and_passivity_gates() {
    let matrix = normal_mobility(&[0.5], 1);
    let c = contact(1, [-1.0, 3.0, 4.0], 0.5);
    let mut polls = 0;
    solve_frictional_contacts(&[c], &matrix, FrictionalContactConfig::default(), || {
        polls += 1;
        false
    })
    .unwrap();
    for checkpoint in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            solve_frictional_contacts(&[c], &matrix, FrictionalContactConfig::default(), || {
                seen += 1;
                seen == checkpoint
            }),
            Err(FrictionalContactError::Cancelled)
        ));
    }
    let config = FrictionalContactConfig {
        maximum_evaluations: 1,
        ..FrictionalContactConfig::default()
    };
    assert!(matches!(
        solve_frictional_contacts(&[c], &matrix, config, || false),
        Err(FrictionalContactError::EvaluationBudget)
    ));
    let mut config = FrictionalContactConfig::default();
    config.matrix.max_variables = 2;
    assert!(matches!(
        solve_frictional_contacts(&[c], &matrix, config, || false),
        Err(FrictionalContactError::MobilityNotAdmitted(_))
    ));
    let energetic = FrictionalContact {
        normal_law: NormalVelocityLaw::Target { speed: 2.0 },
        ..contact(2, [-1.0, 0.0, 0.0], 0.5)
    };
    assert!(matches!(
        solve_frictional_contacts(
            &[energetic],
            &matrix,
            FrictionalContactConfig::default(),
            || false,
        ),
        Err(FrictionalContactError::EnergyIncrease { .. })
    ));
}
