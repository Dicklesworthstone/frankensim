use super::{close, contact, normal_chain, shared_body_case, solve};
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, FrictionalContactError,
    NormalVelocityLaw,
    elliptic::{
        EllipticFriction, EllipticFrictionalContact, EllipticFrictionalContactError,
        EllipticFrictionalContactResult, solve_elliptic_frictional_contacts,
    },
};

fn directional(id: u64, velocity: [f64; 3], mu: [f64; 2], angle: f64) -> EllipticFrictionalContact {
    EllipticFrictionalContact {
        friction: EllipticFriction {
            coefficients_mu: mu,
            principal_axis: [angle.cos(), angle.sin()],
        },
        ..contact(id, velocity, 0.0).into()
    }
}

fn solve_ellipse(c: &[EllipticFrictionalContact], w: &[f64]) -> EllipticFrictionalContactResult {
    solve_elliptic_frictional_contacts(c, w, Default::default(), || false).unwrap()
}

// Independent support function in physical principal coordinates, without
// using the production transform or its scalar result.
fn sliding_impulse(f: EllipticFriction, pn: f64, velocity: [f64; 3]) -> [f64; 3] {
    let [c, s] = f.principal_axis;
    let mu = f.coefficients_mu;
    let v = [
        c * velocity[1] + s * velocity[2],
        -s * velocity[1] + c * velocity[2],
    ];
    let support = (mu[0] * v[0]).hypot(mu[1] * v[1]);
    let p = [
        -pn * mu[0] * mu[0] * v[0] / support,
        -pn * mu[1] * mu[1] * v[1] / support,
    ];
    [pn, c * p[0] - s * p[1], s * p[0] + c * p[1]]
}

fn coupled_ellipses() -> (Vec<EllipticFrictionalContact>, Vec<f64>, Vec<[f64; 3]>) {
    let (_, w, _) = shared_body_case();
    let endpoint = [[0.0, 2.0, 0.7], [0.0, -0.4, 1.0]];
    let mut c = vec![
        directional(31, endpoint[0], [0.25, 0.15], 0.31),
        directional(7, endpoint[1], [0.18, 0.3], -0.43),
    ];
    let expected = vec![
        sliding_impulse(c[0].friction, 1.0, endpoint[0]),
        sliding_impulse(c[1].friction, 0.8, endpoint[1]),
    ];
    for i in 0..6 {
        for j in 0..6 {
            c[i / 3].relative_velocity[i % 3] -= w[6 * i + j] * expected[j / 3][j % 3];
        }
    }
    (c, w, expected)
}

#[test]
fn directional_sliding_retains_normal_tangent_feedback_and_physical_work() {
    let w = [1.0, 0.3, 0.0, 0.3, 1.0, 0.0, 0.0, 0.0, 1.4];
    let mut c = directional(1, [-1.0, 3.0, 0.0], [0.4, 0.15], 0.0);
    c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.2 };
    let out = solve_ellipse(&[c], &w);
    let pn = 1.2 / 0.88;
    close(out.impulses[0][0], pn);
    close(out.impulses[0][1], -0.4 * pn);
    close(out.impulses[0][2], 0.0);
    close(out.relative_velocities_after[0][0], 0.2);
    close(out.relative_velocities_after[0][1], 3.0 - 0.1 * pn);
    assert!((out.impulses[0][0] - 1.2).abs() > 0.1);
    let independent_work = 0.5 * pn * (-1.0 + 0.2) + 0.5 * (-0.4 * pn) * (6.0 - 0.1 * pn);
    close(out.relative_work, independent_work);
    close(out.dissipated_energy, -independent_work);
    assert!(matches!(
        out.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueImpulses
    ));
    assert!(out.transformed.contraction_bound < 1.0);
    assert_eq!(out.coordinates[0].scalar_friction_coefficient, 0.4);
    close(out.coordinates[0].tangent_impulse_map[1][1], 0.375);
    close(out.transformed_config.velocity_tolerance, 0.375e-9);
    assert!(out.residuals.coordinate_work_residual_j.abs() < 1e-12);
    assert!(out.residuals.ellipse_impulse_violation < 1e-10);
    assert!(out.residuals.sliding_impulse_residual < 1e-10);
}

#[test]
fn independent_coupled_ellipse_root_is_permutation_invariant() {
    let (contacts, w, expected) = coupled_ellipses();
    let out = solve_ellipse(&contacts, &w);
    assert_eq!(out.active_indices, vec![1, 0]);
    let mut work = 0.0;
    for i in 0..2 {
        for a in 0..3 {
            close(out.impulses[i][a], expected[i][a]);
            let v = contacts[i].relative_velocity[a]
                + (0..6)
                    .map(|j| w[(3 * i + a) * 6 + j] * expected[j / 3][j % 3])
                    .sum::<f64>();
            close(out.relative_velocities_after[i][a], v);
            work += expected[i][a] * (contacts[i].relative_velocity[a] + v) * 0.5;
        }
    }
    close(out.relative_work, work);
    close(out.dissipated_energy, -work);
    let mut permuted = vec![0.0; 36];
    for i in 0..6 {
        for j in 0..6 {
            permuted[6 * i + j] = w[6 * ((i + 3) % 6) + (j + 3) % 6];
        }
    }
    let swapped = solve_ellipse(&[contacts[1], contacts[0]], &permuted);
    assert_eq!(out.impulses[0], swapped.impulses[1]);
    assert_eq!(out.impulses[1], swapped.impulses[0]);
    assert_eq!(out.relative_work.to_bits(), swapped.relative_work.to_bits());
    assert_eq!(out.transformed.nonlinear, swapped.transformed.nonlinear);
}

#[test]
fn equal_coefficients_and_zero_friction_recover_original_owner_exactly() {
    let (contacts, w, _) = shared_body_case();
    let original = solve(&contacts, &w);
    let mut directional: Vec<_> = contacts
        .iter()
        .copied()
        .map(EllipticFrictionalContact::from)
        .collect();
    for c in &mut directional {
        c.friction.principal_axis = [0.3, -0.8];
    }
    let out = solve_ellipse(&directional, &w);
    assert_eq!(out.transformed, original);
    assert_eq!(out.impulses, original.impulses);
    assert_eq!(
        out.relative_velocities_after,
        original.relative_velocities_after
    );
    assert_eq!(out.relative_work, original.relative_work);
    for coordinates in out.coordinates {
        assert_eq!(coordinates.tangent_impulse_map, [[1.0, 0.0], [0.0, 1.0]]);
    }

    let mut c = [
        contact(1, [-100.0, 7.0, -2.0], 0.0),
        contact(2, [-1.0, 0.5, 0.1], 0.0),
    ];
    c[0].gap = 0.1;
    let w = normal_chain();
    let original = solve(&c, &w);
    let out = solve_ellipse(&c.map(EllipticFrictionalContact::from), &w);
    assert_eq!(out.transformed, original);
    assert_eq!(out.impulses, original.impulses);
    assert_eq!(
        out.relative_velocities_after,
        original.relative_velocities_after
    );
    close(out.relative_velocities_after[0][0], -100.6);
    assert_eq!(out.impulses[0], [0.0; 3]);
    assert!(solve_ellipse(&[], &[]).impulses.is_empty());
}

#[test]
fn independent_tangent_frame_rotations_preserve_the_physical_solution() {
    let (contacts, w, _) = coupled_ellipses();
    let original = solve_ellipse(&contacts, &w);
    let angles: [f64; 2] = [0.61, -0.47];
    let q = angles.map(|a| {
        [
            [1.0, 0.0, 0.0],
            [0.0, a.cos(), -a.sin()],
            [0.0, a.sin(), a.cos()],
        ]
    });
    let mut rotated = contacts.clone();
    for i in 0..2 {
        for a in 0..3 {
            rotated[i].relative_velocity[a] = (0..3)
                .map(|b| q[i][a][b] * contacts[i].relative_velocity[b])
                .sum();
        }
        let axis = contacts[i].friction.principal_axis;
        rotated[i].friction.principal_axis = [
            q[i][1][1] * axis[0] + q[i][1][2] * axis[1],
            q[i][2][1] * axis[0] + q[i][2][2] * axis[1],
        ];
    }
    let mut rotated_w = vec![0.0; 36];
    for i in 0..6 {
        for j in i..6 {
            for a in 0..3 {
                for b in 0..3 {
                    rotated_w[6 * i + j] += q[i / 3][i % 3][a]
                        * w[(3 * (i / 3) + a) * 6 + 3 * (j / 3) + b]
                        * q[j / 3][j % 3][b];
                }
            }
            rotated_w[6 * j + i] = rotated_w[6 * i + j];
        }
    }
    let out = solve_ellipse(&rotated, &rotated_w);
    for i in 0..2 {
        for a in 0..3 {
            close(
                out.impulses[i][a],
                (0..3).map(|b| q[i][a][b] * original.impulses[i][b]).sum(),
            );
            close(
                out.relative_velocities_after[i][a],
                (0..3)
                    .map(|b| q[i][a][b] * original.relative_velocities_after[i][b])
                    .sum(),
            );
        }
    }
    close(out.relative_work, original.relative_work);
    close(out.dissipated_energy, original.dissipated_energy);
}

#[test]
fn singular_axes_budgets_and_cancellation_refuse_without_changing_inputs() {
    let (contacts, w, _) = coupled_ellipses();
    let before = contacts.clone();
    let matrix_before = w.clone();
    let mut invalid = contacts.clone();
    invalid[0].friction.coefficients_mu[1] = 0.0;
    assert!(matches!(
        solve_elliptic_frictional_contacts(&invalid, &w, Default::default(), || false),
        Err(EllipticFrictionalContactError::UnsupportedDegenerateEllipse { index: 0 })
    ));
    invalid[0].friction.coefficients_mu[1] = -0.1;
    assert!(matches!(
        solve_elliptic_frictional_contacts(&invalid, &w, Default::default(), || false),
        Err(EllipticFrictionalContactError::Invalid(_))
    ));
    invalid = contacts.clone();
    invalid[0].friction.principal_axis = [0.0; 2];
    assert!(matches!(
        solve_elliptic_frictional_contacts(&invalid, &w, Default::default(), || false),
        Err(EllipticFrictionalContactError::Invalid(_))
    ));
    invalid = contacts.clone();
    invalid[0].id = invalid[1].id;
    assert!(matches!(
        solve_elliptic_frictional_contacts(&invalid, &w, Default::default(), || false),
        Err(EllipticFrictionalContactError::Invalid(_))
    ));
    assert!(matches!(
        solve_elliptic_frictional_contacts(&contacts, &w[..35], Default::default(), || false),
        Err(EllipticFrictionalContactError::Invalid(_))
    ));
    let small = FrictionalContactConfig {
        maximum_contacts: 1,
        ..Default::default()
    };
    assert!(matches!(
        solve_elliptic_frictional_contacts(&contacts, &w, small, || false),
        Err(EllipticFrictionalContactError::Invalid(_))
    ));
    let small = FrictionalContactConfig {
        maximum_evaluations: 1,
        ..Default::default()
    };
    assert!(matches!(
        solve_elliptic_frictional_contacts(&contacts, &w, small, || false),
        Err(EllipticFrictionalContactError::Transformed {
            error: FrictionalContactError::EvaluationBudget,
            ..
        })
    ));
    let mut polls = 0;
    let good = solve_elliptic_frictional_contacts(&contacts, &w, Default::default(), || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 10);
    for stop in [1, polls / 2, polls - 1] {
        let mut calls = 0;
        assert!(matches!(
            solve_elliptic_frictional_contacts(&contacts, &w, Default::default(), || {
                calls += 1;
                calls >= stop
            }),
            Err(EllipticFrictionalContactError::Cancelled)
        ));
    }
    assert_eq!(contacts, before);
    assert_eq!(w, matrix_before);
    assert_eq!(solve_ellipse(&contacts, &w), good);
}
