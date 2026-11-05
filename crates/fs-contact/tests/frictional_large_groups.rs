use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContact, EllipticFrictionalContactError,
    solve_elliptic_frictional_contacts,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactError, MAX_FRICTIONAL_CONTACTS, NormalVelocityLaw, solve_frictional_contacts,
};

fn controls(contacts: usize) -> FrictionalContactConfig {
    let mut config = FrictionalContactConfig::default();
    config.maximum_contacts = contacts;
    config.matrix.max_variables = 3 * contacts;
    config.newton.linear_restart = 64;
    config.newton.max_linear_cycles = 16;
    config
}

fn contact(id: usize, velocity: [f64; 3], mu: f64) -> FrictionalContact {
    FrictionalContact {
        id: id as u64,
        gap: 0.0,
        activation_tolerance: 0.0,
        relative_velocity: velocity,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8 * (1.0 + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}

// Actual Delassus matrix of n+1 equal unit-mass solid spheres in a straight
// touching chain: I=2*r^2/5. A contact's normal mobility is 2; its shear mobility
// is 2+2*(5/2)=7. Adjacent contacts share one sphere, giving -1 normal coupling
// and -1+5/2=3/2 shear coupling. The complete 3n-row matrix is SPD and connected.
fn sphere_chain(n: usize) -> Vec<f64> {
    let d = 3 * n;
    let mut w = vec![0.0; d * d];
    for i in 0..n {
        for (axis, diagonal, adjacent) in [(0, 2.0, -1.0), (1, 7.0, 1.5), (2, 7.0, 1.5)] {
            let row = 3 * i + axis;
            w[row * d + row] = diagonal;
            if i + 1 < n {
                w[row * d + row + 3] = adjacent;
                w[(row + 3) * d + row] = adjacent;
            }
        }
    }
    w
}

// Independent local momentum balance, without multiplying the supplied matrix.
fn chain_response(p: &[[f64; 3]]) -> Vec<[f64; 3]> {
    (0..p.len())
        .map(|i| {
            std::array::from_fn(|axis| {
                let (diagonal, adjacent) = if axis == 0 { (2.0, -1.0) } else { (7.0, 1.5) };
                diagonal * p[i][axis]
                    + adjacent * if i > 0 { p[i - 1][axis] } else { 0.0 }
                    + adjacent * if i + 1 < p.len() { p[i + 1][axis] } else { 0.0 }
            })
        })
        .collect()
}

fn normal_impulse(i: usize, n: usize) -> f64 {
    // Strictly concave loads make every initial normal velocity approaching.
    4.0 * ((i + 1) * (n - i)) as f64 / (n * n) as f64
}

fn circular_chain(
    n: usize,
    restitution: f64,
) -> (Vec<FrictionalContact>, Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut impulses = Vec::with_capacity(n);
    let mut after = Vec::with_capacity(n);
    for i in 0..n {
        let pn = normal_impulse(i, n);
        let (tangent, velocity) = match i % 3 {
            0 => ([0.1 * pn, -0.075 * pn], [0.0, 0.0]),
            1 => ([-0.24 * pn, -0.18 * pn], [2.0, 1.5]),
            _ => ([0.18 * pn, -0.24 * pn], [-1.5, 2.0]),
        };
        impulses.push([pn, tangent[0], tangent[1]]);
        after.push([0.0, velocity[0], velocity[1]]);
    }
    let response = chain_response(&impulses);
    let contacts = (0..n)
        .map(|i| {
            let normal_before = -response[i][0] / (1.0 + restitution);
            after[i][0] = -restitution * normal_before;
            let mut c = contact(
                i,
                [
                    normal_before,
                    after[i][1] - response[i][1],
                    after[i][2] - response[i][2],
                ],
                0.3,
            );
            c.normal_law = NormalVelocityLaw::NewtonRestitution {
                coefficient: restitution,
            };
            c
        })
        .collect();
    (contacts, impulses, after)
}

#[test]
fn sixty_four_connected_contacts_preserve_coulomb_restitution_and_input_order() {
    let n = 64;
    let w = sphere_chain(n);
    let (contacts, expected, after) = circular_chain(n, 0.35);
    let result = solve_frictional_contacts(&contacts, &w, controls(n), || false).unwrap();
    assert_eq!(
        result.determinacy,
        FrictionalContactDeterminacy::UniqueImpulses
    );
    assert!(result.minimum_normalized_pivot > 0.0);
    assert!(result.contraction_bound < 1.0 - 1e-8);
    assert_eq!(result.contraction_weights.len(), 2 * n);
    assert!(result.impulse_energy_witness.is_none());
    assert!(result.impulse_work_witness.is_none());
    let mut work = 0.0;
    for i in 0..n {
        for axis in 0..3 {
            close(result.impulses[i][axis], expected[i][axis]);
            close(result.relative_velocities_after[i][axis], after[i][axis]);
            work +=
                0.5 * expected[i][axis] * (contacts[i].relative_velocity[axis] + after[i][axis]);
        }
        close(result.normal_targets[i], after[i][0]);
    }
    close(result.relative_work, work);
    close(result.dissipated_energy, -work);
    assert!(work < 0.0);

    let order: Vec<_> = (0..n).map(|i| (17 * i + 11) % n).collect();
    let shuffled: Vec<_> = order.iter().map(|&i| contacts[i]).collect();
    let d = 3 * n;
    let mut shuffled_w = vec![0.0; d * d];
    for row in 0..d {
        for column in 0..d {
            shuffled_w[row * d + column] =
                w[(3 * order[row / 3] + row % 3) * d + 3 * order[column / 3] + column % 3];
        }
    }
    let shuffled_result =
        solve_frictional_contacts(&shuffled, &shuffled_w, controls(n), || false).unwrap();
    for (i, &original) in order.iter().enumerate() {
        assert_eq!(shuffled_result.impulses[i], result.impulses[original]);
        assert_eq!(
            shuffled_result.relative_velocities_after[i],
            result.relative_velocities_after[original]
        );
    }
    assert_eq!(shuffled_result.relative_work, result.relative_work);
}

#[test]
fn connected_directional_contacts_keep_the_physical_ellipse_and_exact_free_axis() {
    let n = 32;
    let axis = [0.6, 0.8];
    let rotate = |v: [f64; 2]| {
        [
            axis[0] * v[0] - axis[1] * v[1],
            axis[1] * v[0] + axis[0] * v[1],
        ]
    };
    let mut expected = Vec::with_capacity(n);
    let mut after = Vec::with_capacity(n);
    let mut coefficients = Vec::with_capacity(n);
    for i in 0..n {
        let mu: [f64; 2] = if i % 4 == 0 {
            [0.0, 0.15]
        } else {
            [0.15, 0.12]
        };
        let slip = if i % 2 == 0 { [1.5, -2.0] } else { [-2.0, 1.5] };
        let pn = normal_impulse(i, n);
        let support = (mu[0] * slip[0]).hypot(mu[1] * slip[1]);
        let tangent = rotate([
            -pn * mu[0] * mu[0] * slip[0] / support,
            -pn * mu[1] * mu[1] * slip[1] / support,
        ]);
        let velocity = rotate(slip);
        expected.push([pn, tangent[0], tangent[1]]);
        after.push([0.0, velocity[0], velocity[1]]);
        coefficients.push(mu);
    }
    let response = chain_response(&expected);
    let contacts: Vec<_> = (0..n)
        .map(|i| EllipticFrictionalContact {
            friction: EllipticFriction {
                coefficients_mu: coefficients[i],
                principal_axis: axis,
            },
            ..contact(
                i,
                std::array::from_fn(|a| after[i][a] - response[i][a]),
                0.0,
            )
            .into()
        })
        .collect();
    let result =
        solve_elliptic_frictional_contacts(&contacts, &sphere_chain(n), controls(n), || false)
            .unwrap();
    assert!(result.transformed.contraction_bound < 1.0 - 1e-8);
    let mut work = 0.0;
    for i in 0..n {
        for a in 0..3 {
            close(result.impulses[i][a], expected[i][a]);
            close(result.relative_velocities_after[i][a], after[i][a]);
            work += 0.5 * expected[i][a] * (contacts[i].relative_velocity[a] + after[i][a]);
        }
        if coefficients[i][0] == 0.0 {
            assert_eq!(result.transformed.impulses[i][1], 0.0);
            assert!(
                (axis[0] * result.impulses[i][1] + axis[1] * result.impulses[i][2]).abs() < 1e-14
            );
        }
    }
    close(result.relative_work, work);
    close(result.dissipated_energy, -work);
}

fn identity(n: usize) -> Vec<f64> {
    let d = 3 * n;
    let mut w = vec![0.0; d * d];
    for i in 0..d {
        w[i * d + i] = 1.0;
    }
    w
}

#[test]
fn large_singular_and_indefinite_mobility_refuse_without_softening() {
    let n = 24;
    let d = 3 * n;
    let contacts: Vec<_> = (0..n).map(|i| contact(i, [-1.0, 0.0, 0.0], 0.0)).collect();
    for coupling in [1.0, 2.0] {
        let mut w = identity(n);
        w[3] = coupling;
        w[3 * d] = coupling;
        assert!(matches!(
            solve_frictional_contacts(&contacts, &w, controls(n), || false),
            Err(FrictionalContactError::MobilityNotAdmitted(_))
        ));
    }

    // The exact free axis's zero coordinate map cannot hide negative ORIGINAL
    // mobility: transformed auxiliary mobility alone would look positive.
    let mut w = identity(n);
    w[d + 1] = -1.0;
    let directional: Vec<_> = contacts
        .iter()
        .map(|&c| EllipticFrictionalContact {
            friction: EllipticFriction {
                coefficients_mu: [0.0, 0.3],
                principal_axis: [1.0, 0.0],
            },
            ..c.into()
        })
        .collect();
    assert!(matches!(
        solve_elliptic_frictional_contacts(&directional, &w, controls(n), || false),
        Err(EllipticFrictionalContactError::OriginalMobility(
            FrictionalContactError::MobilityNotAdmitted(_)
        ))
    ));
}

#[test]
fn large_input_with_small_active_psd_group_retains_velocity_determinacy() {
    let n = 32;
    let d = 3 * n;
    let mut w = identity(n);
    for a in 0..3 {
        w[a * d + 3 + a] = 1.0;
        w[(3 + a) * d + a] = 1.0;
    }
    let contacts: Vec<_> = (0..n)
        .map(|i| {
            let mut c = contact(i, if i < 2 { [-1.0, 0.0, 0.0] } else { [0.0; 3] }, 0.0);
            if i >= 2 {
                c.gap = 1.0;
            }
            c
        })
        .collect();
    let result = solve_frictional_contacts(&contacts, &w, controls(n), || false).unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(witness) = &result.determinacy else {
        panic!("redundant active impulses cannot claim unique local loads");
    };
    assert_eq!(witness.dimension, 6);
    assert_eq!(witness.rank, 3);
    close(result.impulses[0][0] + result.impulses[1][0], 1.0);
    close(result.dissipated_energy, 0.5);
    assert_eq!(result.active_indices, [0, 1]);
    for p in &result.impulses[2..] {
        assert_eq!(*p, [0.0; 3]);
    }
}

#[test]
fn large_spd_candidate_without_global_certificate_remains_unresolved() {
    let n = 24;
    let d = 3 * n;
    let mut w = identity(n);
    for i in 0..n {
        let a = 3 * i;
        w[a * d + a + 1] = 0.8;
        w[(a + 1) * d + a] = 0.8;
        if i + 1 < n {
            w[a * d + a + 3] = 0.005;
            w[(a + 3) * d + a] = 0.005;
        }
    }
    // All rows are connected. The sliding candidate is exact and passive, but
    // strong normal/shear feedback defeats the global block contraction bound.
    let expected = [25.0 / 9.0, -20.0 / 9.0, 0.0];
    let contacts: Vec<_> = (0..n)
        .map(|i| {
            let neighbors = usize::from(i > 0) + usize::from(i + 1 < n);
            contact(
                i,
                [-1.0 - 0.005 * neighbors as f64 * expected[0], 2.0, 0.0],
                0.8,
            )
        })
        .collect();
    let error = solve_frictional_contacts(&contacts, &w, controls(n), || false).unwrap_err();
    let FrictionalContactError::UnresolvedDeterminacy { candidate_impulses } = error else {
        panic!("an uncertified large candidate must remain unresolved: {error:?}");
    };
    for p in candidate_impulses {
        for a in 0..3 {
            close(p[a], expected[a]);
        }
    }
}

#[test]
fn large_group_keeps_original_energy_cap_and_explicit_allocation_limits() {
    let n = 24;
    let contacts: Vec<_> = (0..n)
        .map(|i| {
            let mut c = contact(i, [-0.1, 0.0, 0.0], 0.0);
            c.normal_law = NormalVelocityLaw::Target { speed: 0.3 };
            c
        })
        .collect();
    let w = sphere_chain(n);
    assert!(matches!(
        solve_frictional_contacts(&contacts, &w, controls(n), || false),
        Err(FrictionalContactError::EnergyIncrease { .. })
    ));
    assert!(matches!(
        solve_frictional_contacts(&contacts, &w, FrictionalContactConfig::default(), || false),
        Err(FrictionalContactError::Invalid(_))
    ));
    let mut config = controls(n);
    config.matrix.max_variables = 64;
    assert!(matches!(
        solve_frictional_contacts(&contacts, &w, config, || false),
        Err(FrictionalContactError::MobilityNotAdmitted(_))
    ));
    config.maximum_contacts = MAX_FRICTIONAL_CONTACTS + 1;
    assert!(matches!(
        solve_frictional_contacts(&[], &[], config, || false),
        Err(FrictionalContactError::Invalid(_))
    ));
}

#[test]
fn cancellation_and_newton_budget_never_return_a_large_partial_response() {
    let n = 24;
    let (contacts, _, _) = circular_chain(n, 0.25);
    let w = sphere_chain(n);
    let mut polls = 0;
    solve_frictional_contacts(&contacts, &w, controls(n), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, 4 * n, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            solve_frictional_contacts(&contacts, &w, controls(n), || {
                seen += 1;
                seen == stop
            }),
            Err(FrictionalContactError::Cancelled)
        ));
    }
    let mut config = controls(n);
    config.maximum_evaluations = 1;
    assert!(matches!(
        solve_frictional_contacts(&contacts, &w, config, || false),
        Err(FrictionalContactError::EvaluationBudget)
    ));
}
