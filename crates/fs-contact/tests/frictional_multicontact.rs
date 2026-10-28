use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, NormalVelocityLaw,
    solve_frictional_contacts,
};

fn contact(id: u64, velocity: [f64; 3], mu: f64) -> FrictionalContact {
    FrictionalContact {
        id,
        gap: 0.0,
        activation_tolerance: 0.0,
        relative_velocity: velocity,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8, "{a:.15e} != {b:.15e}");
}
fn solve(
    c: &[FrictionalContact],
    w: &[f64],
) -> fs_contact::frictional_multicontact::FrictionalContactResult {
    solve_frictional_contacts(c, w, FrictionalContactConfig::default(), || false).unwrap()
}

#[test]
fn oblique_impact_solves_normal_and_shear_together() {
    let w = [1.0, 0.3, 0.0, 0.3, 1.0, 0.0, 0.0, 0.0, 1.4];
    let mut c = contact(1, [-1.0, 3.0, 0.0], 0.4);
    c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.2 };
    let out = solve(&[c], &w);
    // Independent sliding branch: pn + 0.3*pt = 1.2; pt = -0.4*pn.
    let pn = 1.2 / 0.88;
    close(out.impulses[0][0], pn);
    close(out.impulses[0][1], -0.4 * pn);
    close(out.relative_velocities_after[0][0], 0.2);
    close(out.relative_velocities_after[0][1], 3.0 - 0.1 * pn);
    assert!(
        (out.impulses[0][0] - 1.2).abs() > 0.1,
        "sequential normal solve is wrong"
    );
    let independent_work = 0.5 * pn * (-1.0 + 0.2) + 0.5 * (-0.4 * pn) * (6.0 - 0.1 * pn);
    close(out.relative_work, independent_work);
    close(out.dissipated_energy, -independent_work);
    assert!(out.contraction_bound < 1.0);
    assert!(out.residuals.normal_complementarity < 1e-9);
    assert!(out.residuals.tangential_duality_gap < 1e-9);
    // Rotate the tangential basis: both shear components now slide. A box or
    // component-wise clip would change the physical impulse and work.
    let angle: f64 = 0.72;
    let q = [
        [1.0, 0.0, 0.0],
        [0.0, angle.cos(), angle.sin()],
        [0.0, -angle.sin(), angle.cos()],
    ];
    let mut rotated_w = [0.0; 9];
    let mut rotated = c;
    rotated.relative_velocity = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            rotated.relative_velocity[i] += q[i][j] * c.relative_velocity[j];
        }
        for j in 0..3 {
            for a in 0..3 {
                for b in 0..3 {
                    rotated_w[3 * i + j] += q[i][a] * w[3 * a + b] * q[j][b];
                }
            }
        }
    }
    let rotated_out = solve(&[rotated], &rotated_w);
    for i in 0..3 {
        let expected: f64 = (0..3).map(|j| q[i][j] * out.impulses[0][j]).sum();
        close(rotated_out.impulses[0][i], expected);
    }
    close(rotated_out.relative_work, out.relative_work);
}

#[test]
fn coupled_sticking_and_explicit_target_have_distinct_normal_laws() {
    let w = [1.0, 0.2, 0.1, 0.2, 1.5, -0.1, 0.1, -0.1, 2.0];
    let expected = [1.0, -0.15, 0.12];
    let mut before = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            before[i] -= w[3 * i + j] * expected[j];
        }
    }
    let out = solve(&[contact(9, before, 0.5)], &w);
    for i in 0..3 {
        close(out.impulses[0][i], expected[i]);
        close(out.relative_velocities_after[0][i], 0.0);
    }
    assert!(out.impulses[0][1].hypot(out.impulses[0][2]) < 0.5 * out.impulses[0][0]);

    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let mut target = contact(10, [-2.0, 0.0, 0.0], 2.0);
    target.normal_law = NormalVelocityLaw::Target { speed: 0.25 };
    let target_out = solve(&[target], &identity);
    close(target_out.impulses[0][0], 2.25);
    close(target_out.normal_targets[0], 0.25);
    assert!(
        target_out.contraction_bound < 1.0,
        "large mu alone does not imply ambiguity"
    );
}

fn shared_body_case() -> (Vec<FrictionalContact>, Vec<f64>, Vec<[f64; 3]>) {
    // Three translation masses 1,10,2: each contact shares the middle body.
    // Rotating the second contact frame creates cross-contact normal/shear terms.
    let angle: f64 = 0.6;
    let q = [
        [angle.cos(), angle.sin(), 0.0],
        [-angle.sin(), angle.cos(), 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut b = [[0.0; 9]; 6];
    for a in 0..3 {
        b[a][a] = 1.0;
        b[a][3 + a] = -1.0 / 10.0_f64.sqrt();
    }
    for a in 0..3 {
        for j in 0..3 {
            b[3 + a][3 + j] = q[a][j] / 10.0_f64.sqrt();
            b[3 + a][6 + j] = -q[a][j] / 2.0_f64.sqrt();
        }
    }
    let mut w = vec![0.0; 36];
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..9 {
                w[6 * i + j] += b[i][k] * b[j][k];
            }
        }
    }
    let expected = vec![[1.0, -0.25, 0.0], [0.8, 0.0, -0.2]];
    let endpoint = [[0.0, 2.0, 0.0], [0.0, 0.0, 1.0]];
    let mut initial = endpoint;
    for i in 0..6 {
        for j in 0..6 {
            initial[i / 3][i % 3] -= w[6 * i + j] * expected[j / 3][j % 3];
        }
    }
    (
        vec![contact(31, initial[0], 0.25), contact(7, initial[1], 0.25)],
        w,
        expected,
    )
}

#[test]
fn simultaneous_shared_body_sliding_is_permutation_invariant() {
    let (contacts, w, expected) = shared_body_case();
    let out = solve(&contacts, &w);
    assert_eq!(out.active_indices, vec![1, 0]);
    for i in 0..2 {
        for a in 0..3 {
            close(out.impulses[i][a], expected[i][a]);
        }
    }
    close(out.relative_velocities_after[0][0], 0.0);
    close(out.relative_velocities_after[1][0], 0.0);
    close(out.relative_velocities_after[0][1], 2.0);
    close(out.relative_velocities_after[1][2], 1.0);
    let mut permuted = vec![0.0; 36];
    for i in 0..6 {
        for j in 0..6 {
            permuted[6 * i + j] = w[6 * ((i + 3) % 6) + (j + 3) % 6];
        }
    }
    let swapped = solve(&[contacts[1], contacts[0]], &permuted);
    assert_eq!(out.impulses[0], swapped.impulses[1]);
    assert_eq!(out.impulses[1], swapped.impulses[0]);
    assert_eq!(out.relative_work.to_bits(), swapped.relative_work.to_bits());
    assert_eq!(out.nonlinear, swapped.nonlinear);
}

fn normal_chain() -> Vec<f64> {
    let mut w = vec![0.0; 36];
    for i in 0..6 {
        w[6 * i + i] = 2.5;
    }
    w[3] = -1.5;
    w[18] = -1.5;
    w
}

#[test]
fn zero_friction_matches_coupled_normal_solution_and_inactive_cross_response() {
    let w = normal_chain();
    let mut contacts = [
        contact(1, [-1.2, 0.0, 0.0], 0.0),
        contact(2, [-0.8, 0.0, 0.0], 0.0),
    ];
    for c in &mut contacts {
        c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 };
    }
    let out = solve(&contacts, &w);
    close(out.impulses[0][0], 1.575);
    close(out.impulses[1][0], 1.425);
    close(out.dissipated_energy, 0.7575);
    contacts[0].gap = 0.1;
    contacts[0].relative_velocity[0] = -100.0;
    contacts[1].relative_velocity[0] = -1.0;
    contacts[1].normal_law = NormalVelocityLaw::Target { speed: 0.0 };
    let inactive = solve(&contacts, &w);
    assert_eq!(inactive.impulses[0], [0.0; 3]);
    close(inactive.impulses[1][0], 0.4);
    close(inactive.relative_velocities_after[0][0], -100.6);
}

#[test]
fn restitution_cannot_inject_global_energy_through_an_initially_separating_contact() {
    let mut c = [
        contact(1, [-1.0, 0.0, 0.0], 0.0),
        contact(2, [0.5, 0.0, 0.0], 0.0),
    ];
    for contact in &mut c {
        contact.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 };
    }
    match solve_frictional_contacts(&c, &normal_chain(), Default::default(), || false).unwrap_err()
    {
        FrictionalContactError::EnergyIncrease { work, .. } => close(work, 0.109375),
        e => panic!("wrong refusal: {e:?}"),
    }
    for contact in &mut c {
        contact.normal_law = NormalVelocityLaw::Target { speed: 0.0 };
    }
    let plastic = solve(&c, &normal_chain());
    close(plastic.impulses[0][0], 0.4375);
    close(plastic.impulses[1][0], 0.0625);
    close(plastic.dissipated_energy, 0.203125);
}

#[test]
fn two_passive_coulomb_roots_are_reported_instead_of_arbitrary_selection() {
    let w = [1.0, 0.8, 0.0, 0.8, 1.0, 0.0, 0.0, 0.0, 1.0];
    let c = contact(1, [0.1, 1.0, 0.0], 2.0);
    match solve_frictional_contacts(&[c], &w, Default::default(), || false).unwrap_err() {
        FrictionalContactError::Nonunique { first, second } => {
            assert_eq!(first[0], [0.0; 3]);
            // There are THREE independently derived passive roots: zero,
            // sliding [1/6,-1/3,0], and sticking [35/18,-23/9,0]. Either
            // nonzero witness suffices; no particular branch is authoritative.
            if second[0][0] > 1.0 {
                close(second[0][0], 35.0 / 18.0);
                close(second[0][1], -23.0 / 9.0);
                close(1.0 + 0.8 * second[0][0] + second[0][1], 0.0);
                assert!(second[0][1].abs() < 2.0 * second[0][0]);
            } else {
                close(second[0][0], 1.0 / 6.0);
                close(second[0][1], -1.0 / 3.0);
                close(1.0 + 0.8 * second[0][0] + second[0][1], 0.8);
            }
            close(0.1 + second[0][0] + 0.8 * second[0][1], 0.0);
            let vt_after = 1.0 + 0.8 * second[0][0] + second[0][1];
            assert!(0.05 * second[0][0] + 0.5 * (1.0 + vt_after) * second[0][1] < 0.0);
        }
        e => panic!("missed competing physical branches: {e:?}"),
    }
    // This case has a root but lies outside the sufficient contraction bound.
    let unresolved = solve_frictional_contacts(
        &[contact(2, [-1.0, 2.0, 0.0], 0.8)],
        &w,
        Default::default(),
        || false,
    );
    assert!(
        matches!(
            unresolved,
            Err(FrictionalContactError::UnresolvedDeterminacy { .. })
        ),
        "unexpected determinacy classification: {unresolved:?}"
    );
}

#[test]
fn degeneracy_penetration_and_work_exhaustion_never_become_physical_authority() {
    let c = contact(1, [-1.0, 2.0, 0.0], 0.3);
    let singular = [1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert!(matches!(
        solve_frictional_contacts(&[c], &singular, Default::default(), || false),
        Err(FrictionalContactError::MobilityNotAdmitted(_))
    ));
    let mut penetrating = c;
    penetrating.gap = -1e-15;
    penetrating.activation_tolerance = 1e-3;
    assert!(matches!(
        solve_frictional_contacts(&[penetrating], &singular, Default::default(), || false),
        Err(FrictionalContactError::PenetratingContact { .. })
    ));
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let limited = FrictionalContactConfig {
        maximum_evaluations: 1,
        ..Default::default()
    };
    assert!(matches!(
        solve_frictional_contacts(&[c], &identity, limited, || false),
        Err(FrictionalContactError::EvaluationBudget)
    ));
    let (coupled, w, _) = shared_body_case();
    let limited_iterations = FrictionalContactConfig {
        maximum_newton_iterations: 1,
        ..Default::default()
    };
    assert!(matches!(
        solve_frictional_contacts(&coupled, &w, limited_iterations, || false),
        Err(FrictionalContactError::NotConverged(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        solve_frictional_contacts(&[c], &identity, Default::default(), || {
            polls += 1;
            polls == 20
        }),
        Err(FrictionalContactError::Cancelled)
    ));
    let empty = solve(&[], &[]);
    assert_eq!(empty.relative_work, 0.0);
    assert!(empty.impulses.is_empty());
    assert!(matches!(
        solve_frictional_contacts(&[], &[], Default::default(), || true),
        Err(FrictionalContactError::Cancelled)
    ));
}
