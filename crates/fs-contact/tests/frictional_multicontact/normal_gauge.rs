use super::{close, contact, solve};
use fs_contact::frictional_multicontact::{
    FrictionalContactDeterminacy, FrictionalContactError, FrictionalContactResult,
    PsdVelocityCriterion, solve_frictional_contacts,
};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn arms(height: f64) -> [[f64; 3]; 4] {
    [
        [-1.0, -1.0, -height],
        [1.0, -1.0, -height],
        [1.0, 1.0, -height],
        [-1.0, 1.0, -height],
    ]
}

fn gram(arms: &[[f64; 3]], inertia: [f64; 3]) -> Vec<f64> {
    let axes = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let rows: Vec<_> = arms
        .iter()
        .flat_map(|&r| {
            axes.map(|a| {
                let angular = cross(r, a);
                [a[0], a[1], a[2], angular[0], angular[1], angular[2]]
            })
        })
        .collect();
    rows.iter()
        .flat_map(|a| {
            rows.iter().map(|b| {
                (0..3)
                    .map(|i| a[i] * b[i] / 3.0 + a[i + 3] * b[i + 3] / inertia[i])
                    .sum()
            })
        })
        .collect()
}

fn direction_witness(result: &FrictionalContactResult) {
    let FrictionalContactDeterminacy::UniqueVelocity(w) = &result.determinacy else {
        panic!("four-foot support has a physical pressure gauge");
    };
    assert_eq!((w.dimension, w.rank), (12, 6));
    assert!(w.weighted_normal_nullspace_residual > 1e-3);
    assert!(w.friction_feedback_bound < 0.2);
    let PsdVelocityCriterion::CandidateFrictionDirection {
        nullspace_coupling_residual,
        energy_coercivity_lower_bound,
        energy_spectral_residual,
        energy_orthogonality_residual,
        traction_residual,
        subgradient_velocity_residual,
        ..
    } = &w.criterion
    else {
        panic!("normal observability does not hold here")
    };
    assert!(*nullspace_coupling_residual < 1e-13);
    assert!(*energy_coercivity_lower_bound > 0.8);
    assert!(*energy_spectral_residual < 1e-13);
    assert!(*energy_orthogonality_residual < 1e-13);
    assert!(*traction_residual < 1e-10);
    assert!(*subgradient_velocity_residual < 1e-10);
}

#[test]
fn four_point_cube_stick_threshold_and_circular_slip_have_unique_motion() {
    // Actual mass-three uniform cube of side two: Ixx=Iyy=Izz=2, with
    // four coplanar lower corners one metre below its centre of mass.
    let points = arms(1.0);
    let w = gram(&points, [2.0; 3]);
    for direction in [[1.0, 0.0], [0.6, 0.8]] {
        for speed in [0.02_f64, 0.1, 3.0] {
            let contacts = [51, 7, 23, 17]
                .map(|id| contact(id, [-1.0, speed * direction[0], speed * direction[1]], 0.1));
            let result = solve(&contacts, &w);
            direction_witness(&result);
            let traction = speed.min(0.1);
            let after = speed - traction;
            let mut force = [0.0; 3];
            let mut torque = [0.0; 3];
            let mut normal_moment = [0.0; 2];
            for (i, p) in result.impulses.iter().enumerate() {
                let f = [p[1], p[2], p[0]];
                let t = cross(points[i], f);
                for a in 0..3 {
                    force[a] += f[a];
                    torque[a] += t[a];
                }
                for a in 0..2 {
                    normal_moment[a] += points[i][a] * p[0];
                }
                close(result.relative_velocities_after[i][0], 0.0);
                close(result.relative_velocities_after[i][1], after * direction[0]);
                close(result.relative_velocities_after[i][2], after * direction[1]);
            }
            close(force[2], 3.0);
            for a in 0..2 {
                close(force[a], -3.0 * traction * direction[a]);
                close(normal_moment[a], 3.0 * traction * direction[a]);
            }
            for value in torque {
                close(value, 0.0);
            }
            close(
                result.dissipated_energy,
                1.5 * (1.0 + speed * speed - after * after),
            );

            // Two independent pressure distributions, both obeying the SAME
            // Coulomb law and resulting in the SAME body velocity. The diagonal
            // load transfer has zero force and torque, including its shear.
            let mut previous = None;
            for gauge in [-0.1, 0.1] {
                let candidate: Vec<_> = points
                    .iter()
                    .enumerate()
                    .map(|(i, r)| {
                        let pn = 0.75
                            * (1.0 + traction * (r[0] * direction[0] + r[1] * direction[1]))
                            + if i % 2 == 0 { gauge } else { -gauge };
                        assert!(pn > 0.0);
                        [
                            pn,
                            -traction * direction[0] * pn,
                            -traction * direction[1] * pn,
                        ]
                    })
                    .collect();
                if let Some(previous) = previous {
                    assert_ne!(candidate, previous);
                }
                for i in 0..12 {
                    let expected = if i % 3 == 0 {
                        1.0
                    } else {
                        -traction * direction[i % 3 - 1]
                    };
                    close(
                        (0..12)
                            .map(|j| w[12 * i + j] * candidate[j / 3][j % 3])
                            .sum(),
                        expected,
                    );
                }
                previous = Some(candidate);
            }
        }
    }
}

#[test]
fn four_point_admission_is_covariant_under_each_tangent_frame_and_input_order() {
    let w = gram(&arms(1.0), [2.0; 3]);
    let contacts = [51, 7, 23, 17].map(|id| contact(id, [-1.0, 0.012, 0.016], 0.1));
    let reference = solve(&contacts, &w);
    let rotations = [0.3_f64, -0.7, 1.1, 2.2].map(|a| {
        [
            [1.0, 0.0, 0.0],
            [0.0, a.cos(), a.sin()],
            [0.0, -a.sin(), a.cos()],
        ]
    });
    let mut rotated = contacts;
    let mut rotated_w = vec![0.0; 144];
    for i in 0..12 {
        rotated[i / 3].relative_velocity[i % 3] = (0..3)
            .map(|a| rotations[i / 3][i % 3][a] * contacts[i / 3].relative_velocity[a])
            .sum();
        for j in 0..12 {
            for a in 0..3 {
                for b in 0..3 {
                    rotated_w[12 * i + j] += rotations[i / 3][i % 3][a]
                        * w[12 * (3 * (i / 3) + a) + 3 * (j / 3) + b]
                        * rotations[j / 3][j % 3][b];
                }
            }
        }
    }
    let result = solve(&rotated, &rotated_w);
    direction_witness(&result);
    close(result.relative_work, reference.relative_work);
    for i in 0..4 {
        for a in 0..3 {
            close(result.relative_velocities_after[i][a], 0.0);
        }
    }
    let permutation = [2, 0, 3, 1];
    let mut permuted_w = vec![0.0; 144];
    for i in 0..12 {
        for j in 0..12 {
            permuted_w[12 * i + j] =
                rotated_w[12 * (3 * permutation[i / 3] + i % 3) + 3 * permutation[j / 3] + j % 3];
        }
    }
    let other = solve(&permutation.map(|i| rotated[i]), &permuted_w);
    assert_eq!(result.determinacy, other.determinacy);
    assert_eq!(
        result.relative_work.to_bits(),
        other.relative_work.to_bits()
    );
    for (i, &j) in permutation.iter().enumerate() {
        assert_eq!(other.impulses[i], result.impulses[j]);
        assert_eq!(
            other.relative_velocities_after[i],
            result.relative_velocities_after[j]
        );
    }
}

fn independent_sliding_root(w: &[f64], mu: &[f64; 4], p: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let mut velocity = vec![[-1.0, 3.0, 0.0]; 4];
    let mut work = 0.0;
    for i in 0..12 {
        velocity[i / 3][i % 3] += (0..12)
            .map(|j| w[12 * i + j] * p[j / 3][j % 3])
            .sum::<f64>();
    }
    for i in 0..4 {
        let v = velocity[i];
        close(v[0], 0.0);
        assert!(p[i][0] >= 0.0 && v[1] > 0.0);
        close(p[i][1], -mu[i] * p[i][0]);
        close(p[i][2], 0.0);
        work += 0.5 * (-p[i][0] + p[i][1] * (3.0 + v[1]) + p[i][2] * v[2]);
    }
    assert!(work < 0.0);
    velocity
}

#[test]
fn heterogeneous_friction_exposes_two_passive_roots_with_different_motion() {
    // A mass-three thin square plate: Ixx=Iyy=1, Izz=2. Two diagonal
    // normal-load allocations create different total sliding friction while
    // each has zero torque. No numerical uniqueness test can pick a load law.
    let w = gram(&arms(0.0), [1.0, 1.0, 2.0]);
    let mu = [0.1, 0.3, 0.1, 0.3];
    for (diagonal, expected) in [(0, 2.9), (1, 2.7)] {
        let p: Vec<_> = (0..4)
            .map(|i| {
                let pn = if i % 2 == diagonal { 1.5 } else { 0.0 };
                [pn, -mu[i] * pn, 0.0]
            })
            .collect();
        for v in independent_sliding_root(&w, &mu, &p) {
            close(v[1], expected);
        }
    }
    let contacts = [0, 1, 2, 3].map(|i| contact(i as u64, [-1.0, 3.0, 0.0], mu[i]));
    let error = solve_frictional_contacts(&contacts, &w, Default::default(), || false).unwrap_err();
    let FrictionalContactError::Nonunique { first, second } = error else {
        panic!("the loaded probe should expose the independently known ambiguity: {error:?}");
    };
    let first_velocity = independent_sliding_root(&w, &mu, &first);
    let second_velocity = independent_sliding_root(&w, &mu, &second);
    assert!((first_velocity[0][1] - second_velocity[0][1]).abs() > 1e-3);
}

#[test]
fn tiny_heterogeneous_friction_is_not_mistaken_for_a_null_feedback_identity() {
    let w = gram(&arms(0.0), [1.0, 1.0, 2.0]);
    let mu = [1e-12, 3e-12, 1e-12, 3e-12];
    let contacts = [0, 1, 2, 3].map(|i| contact(i as u64, [-1.0, 3.0, 0.0], mu[i]));
    // Distinct motions are below the configured velocity resolution, so the
    // bounded probe may remain unresolved. It still must not claim uniqueness
    // by using an absolute zero test on a physically nonzero friction coupling.
    assert!(matches!(
        solve_frictional_contacts(&contacts, &w, Default::default(), || false),
        Err(FrictionalContactError::UnresolvedDeterminacy { .. })
            | Err(FrictionalContactError::Nonunique { .. })
    ));
}
