use super::{close, contact, solve};
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, FrictionalContactError,
    FrictionalContactResult, solve_frictional_contacts,
};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Actual point Jacobians of one mass-one, identity-inertia rigid body. The
// point frame is [world z normal, world x shear, world y shear]. Three points
// give nine contact rows but only six native mechanical degrees of freedom.
fn gram(arms: &[[f64; 3]]) -> Vec<f64> {
    let axes = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let rows: Vec<_> = arms
        .iter()
        .flat_map(|&r| {
            axes.map(|axis| {
                let angular = cross(r, axis);
                [
                    axis[0], axis[1], axis[2], angular[0], angular[1], angular[2],
                ]
            })
        })
        .collect();
    rows.iter()
        .flat_map(|a| {
            rows.iter()
                .map(|b| a.iter().zip(b).map(|(x, y)| x * y).sum())
        })
        .collect()
}

fn triangle() -> [[f64; 3]; 3] {
    [[-1.0, -1.0, -0.25], [1.0, -1.0, -0.25], [0.0, 1.0, -0.25]]
}

fn native_increment(arms: &[[f64; 3]], impulses: &[[f64; 3]]) -> [f64; 6] {
    let mut output = [0.0; 6];
    for (&arm, &p) in arms.iter().zip(impulses) {
        let force = [p[1], p[2], p[0]];
        let torque = cross(arm, force);
        for a in 0..3 {
            output[a] += force[a];
            output[3 + a] += torque[a];
        }
    }
    output
}

fn velocity_witness(out: &FrictionalContactResult, dimension: usize) {
    let FrictionalContactDeterminacy::UniqueVelocity(w) = &out.determinacy else {
        panic!("a rank-deficient manifold cannot claim unique contact loads");
    };
    assert_eq!(w.dimension, dimension);
    assert_eq!(w.rank, 6);
    assert!(w.friction_feedback_bound < 1.0);
    assert!(w.spectral_residual < 1e-13);
    assert!(w.orthogonality_residual < 1e-13);
    assert!(w.weighted_normal_nullspace_residual < 1e-13);
    assert!(w.minimum_retained_eigenvalue > 1e-3);
    assert!(w.maximum_null_eigenvalue_magnitude < 1e-13);
    assert_eq!(out.contraction_bound, 1.0);
    assert!(out.contraction_weights.is_empty());
    assert_eq!(out.minimum_normalized_pivot, 0.0);
}

#[test]
fn rigid_tripod_sticks_with_unique_motion_and_an_explicit_shear_impulse_gauge() {
    let arms = triangle();
    let w = gram(&arms);
    let contacts = [
        contact(51, [-1.0, 0.1, 0.0], 0.2),
        contact(7, [-1.0, 0.1, 0.0], 0.2),
        contact(23, [-1.0, 0.1, 0.0], 0.2),
    ];
    let out = solve(&contacts, &w);
    velocity_witness(&out, 9);
    for v in &out.relative_velocities_after {
        for &a in v {
            close(a, 0.0);
        }
    }
    // Independent force/torque balance: sum pn=1, sum y*pn=0,
    // sum x*pn=0.25*0.1; no assumed shear allocation is needed.
    for (p, expected) in out.impulses.iter().zip([0.2375, 0.2625, 0.5]) {
        close(p[0], expected);
        assert!(p[1].hypot(p[2]) < 0.2 * p[0]);
    }
    let expected_increment = [-0.1, 0.0, 1.0, 0.0, 0.0, 0.0];
    for (a, b) in native_increment(&arms, &out.impulses)
        .iter()
        .zip(expected_increment)
    {
        close(*a, b);
    }
    close(out.dissipated_energy, 0.505);

    // A second, independently constructed valid contact-load allocation:
    // opposing collinear shear on the first edge changes neither force nor
    // torque. Both endpoints remain strictly inside their Coulomb disks.
    let mut alternative = out.impulses.clone();
    alternative[0][1] += 0.001;
    alternative[1][1] -= 0.001;
    assert_ne!(alternative, out.impulses);
    for p in &alternative {
        assert!(p[1].hypot(p[2]) < 0.2 * p[0]);
    }
    for i in 0..9 {
        let response: f64 = (0..9)
            .map(|j| w[9 * i + j] * (alternative[j / 3][j % 3] - out.impulses[j / 3][j % 3]))
            .sum();
        close(response, 0.0);
    }
    for (a, b) in native_increment(&arms, &alternative)
        .iter()
        .zip(expected_increment)
    {
        close(*a, b);
    }
}

#[test]
fn rigid_tripod_slides_with_coupled_load_transfer_and_canonical_order() {
    let arms = triangle();
    let w = gram(&arms);
    let contacts = [
        contact(51, [-1.0, 3.0, 0.0], 0.2),
        contact(7, [-1.0, 3.0, 0.0], 0.2),
        contact(23, [-1.0, 3.0, 0.0], 0.2),
    ];
    let out = solve(&contacts, &w);
    velocity_witness(&out, 9);
    // Friction at finite height shifts the normal resultant by h*mu.
    // A normal-first solve would keep [0.25,0.25,0.5] and induce rotation.
    for (i, pn) in [0.225, 0.275, 0.5].into_iter().enumerate() {
        close(out.impulses[i][0], pn);
        close(out.impulses[i][1], -0.2 * pn);
        close(out.impulses[i][2], 0.0);
        close(out.relative_velocities_after[i][0], 0.0);
        close(out.relative_velocities_after[i][1], 2.8);
        close(out.relative_velocities_after[i][2], 0.0);
    }
    close(out.dissipated_energy, 1.08);
    let permutation = [2, 0, 1];
    let mut reordered_w = vec![0.0; 81];
    for i in 0..9 {
        for j in 0..9 {
            reordered_w[9 * i + j] =
                w[9 * (3 * permutation[i / 3] + i % 3) + 3 * permutation[j / 3] + j % 3];
        }
    }
    let other = solve(&permutation.map(|i| contacts[i]), &reordered_w);
    assert_eq!(out.determinacy, other.determinacy);
    assert_eq!(out.relative_work.to_bits(), other.relative_work.to_bits());
    for (i, &original) in permutation.iter().enumerate() {
        assert_eq!(out.impulses[original], other.impulses[i]);
        assert_eq!(
            out.relative_velocities_after[original],
            other.relative_velocities_after[i]
        );
    }
}

#[test]
fn redundant_frictionless_normal_support_has_unique_motion_without_unique_loads() {
    let arms = [
        [-1.0, -1.0, -0.25],
        [1.0, -1.0, -0.25],
        [1.0, 1.0, -0.25],
        [-1.0, 1.0, -0.25],
    ];
    let w = gram(&arms);
    let contacts = [1, 2, 3, 4].map(|id| contact(id, [-1.0, 3.0, -0.2], 0.0));
    let out = solve(&contacts, &w);
    velocity_witness(&out, 12);
    for v in &out.relative_velocities_after {
        close(v[0], 0.0);
        close(v[1], 3.0);
        close(v[2], -0.2);
    }
    for (a, b) in native_increment(&arms, &out.impulses)
        .iter()
        .zip([0.0, 0.0, 1.0, 0.0, 0.0, 0.0])
    {
        close(*a, b);
    }
    close(out.dissipated_energy, 0.5);
    let mut frictionful = contacts;
    for c in &mut frictionful {
        c.friction_coefficient = 0.01;
    }
    // The normal-observable proof cannot admit nonzero friction here. The
    // candidate-direction energy proof can, without choosing physical pressure.
    let frictionful = solve(&frictionful, &w);
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = frictionful.determinacy else {
        panic!("normal loads remain indeterminate");
    };
    assert!(matches!(
        witness.criterion,
        fs_contact::frictional_multicontact::PsdVelocityCriterion::CandidateFrictionDirection { .. }
    ));
}

#[test]
fn mobility_admission_distinguishes_bad_rank_from_a_large_friction_coefficient() {
    let mut contacts = [1, 2, 3].map(|id| contact(id, [-1.0, 0.1, 0.0], 0.2));
    let w = gram(&triangle());
    for eigenvalue in [-1e-5, 1e-11] {
        let mut hostile = w.clone();
        // [t1 at point 0, -t1 at point 1]/sqrt(2) is an exact force/torque
        // null vector. Give it an inadmissible negative or unresolved eigenvalue.
        hostile[9 + 1] += 0.5 * eigenvalue;
        hostile[9 + 4] -= 0.5 * eigenvalue;
        hostile[36 + 1] -= 0.5 * eigenvalue;
        hostile[36 + 4] += 0.5 * eigenvalue;
        assert!(matches!(
            solve_frictional_contacts(&contacts, &hostile, Default::default(), || false),
            Err(FrictionalContactError::MobilityNotAdmitted(_))
        ));
    }
    for c in &mut contacts {
        c.friction_coefficient = 1.0;
    }
    // A large available friction coefficient can defeat the worst-case norm
    // bound while the actual sticking traction has a positive global margin.
    let result = solve(&contacts, &w);
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = result.determinacy else {
        panic!("redundant shear loads remain indeterminate");
    };
    assert!(matches!(
        witness.criterion,
        fs_contact::frictional_multicontact::PsdVelocityCriterion::CandidateFrictionDirection { .. }
    ));
    for velocity in result.relative_velocities_after {
        for value in velocity {
            close(value, 0.0);
        }
    }
}

#[test]
fn redundant_contact_cancellation_and_work_exhaustion_still_return_no_authority() {
    let contacts = [1, 2, 3].map(|id| contact(id, [-1.0, 0.1, 0.0], 0.2));
    let w = gram(&triangle());
    let mut polls = 0;
    assert!(matches!(
        solve_frictional_contacts(&contacts, &w, Default::default(), || {
            polls += 1;
            polls == 55
        }),
        Err(FrictionalContactError::Cancelled)
    ));
    assert!(matches!(
        solve_frictional_contacts(
            &contacts,
            &w,
            FrictionalContactConfig {
                maximum_evaluations: 1,
                ..Default::default()
            },
            || false
        ),
        Err(FrictionalContactError::EvaluationBudget)
    ));
    let first = solve(&contacts, &w);
    let retry = solve(&contacts, &w);
    assert_eq!(first, retry);
}
