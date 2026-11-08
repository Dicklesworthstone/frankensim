//! G1/G3: original-force Laurent residues of the full Maxwell pencil.

use fs_math::c64::C64;
use fs_modal::relaxation_modes::{
    MaxwellBranchDirection, MaxwellDampingBranch, MaxwellModeError, MaxwellPoleDirection,
    PoleOptions, PreparedMaxwellModes,
};

fn options(rate: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: rate,
        step_tolerance: 2e-13,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 8,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 200_000_000,
    }
}
fn branch(stiffness: f64, relaxation_time_s: f64, strain: &[f64]) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        stiffness,
        relaxation_time_s,
        strain: strain.to_vec(),
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn all_scalar_material_residues_reconstruct_original_force_response() {
    // T(s)^-1=(s+5)/(4*(s+1)*(s+2)*(s+4)). These signed partial
    // fractions follow independently from the original scalar equations.
    let branches = [branch(9.6, 0.2, &[1.0])];
    let prepared =
        PreparedMaxwellModes::new(&[4.0], &[8.0], &[6.4], 1, &branches, options(5.0), || false)
            .unwrap();
    let force = C64::new(2.0, -0.7);
    let mut residues = Vec::new();
    for (pole, coefficient) in [(-1.0, 1.0 / 3.0), (-2.0, -3.0 / 8.0), (-4.0, 1.0 / 24.0)] {
        let residue = prepared
            .refine_residue(
                C64::from_re(pole - 0.02),
                &[C64::new(0.5, 0.3)],
                &[C64::new(-2.0, 0.4)],
                1e-8,
                || false,
            )
            .unwrap();
        close(
            residue.right_mode().displacement[0],
            C64::from_re(0.5),
            1e-14,
        );
        close(
            residue.transpose_displacement()[0],
            C64::from_re(0.5),
            1e-14,
        );
        close(
            residue.derivative_normalization(),
            C64::from_re(0.25 / coefficient),
            3e-12,
        );
        assert!(residue.derivative_alignment() > 0.99);
        assert!(residue.transpose_source_residuals().0 < 2e-12);
        assert!(residue.transpose_source_residuals().1 < 2e-12);
        let applied = residue.force_residue(&[force], || false).unwrap();
        close(applied.displacement[0], force.scale(coefficient), 3e-12);
        close(applied.velocity[0], force.scale(pole * coefficient), 3e-12);
        close(
            applied.branch_stress[0],
            force.scale(9.6 * pole / (pole + 5.0) * coefficient),
            4e-12,
        );
        let zero = residue.force_residue(&[C64::ZERO], || false).unwrap();
        assert_eq!(zero.displacement, [C64::ZERO]);
        assert_eq!(zero.velocity, [C64::ZERO]);
        assert_eq!(zero.branch_stress, [C64::ZERO]);
        residues.push(applied);
    }
    for s in [C64::new(0.0, 0.3), C64::new(0.0, 2.0), C64::new(0.2, 12.0)] {
        let pencil = s * s.scale(4.0)
            + s.scale(8.0)
            + C64::from_re(6.4)
            + (s / (s + C64::from_re(5.0))).scale(9.6);
        let displacement = force / pencil;
        let stress = (s / (s + C64::from_re(5.0))).scale(9.6) * displacement;
        let mut reconstructed = [C64::ZERO; 3];
        for r in &residues {
            for (i, value) in [r.displacement[0], r.velocity[0], r.branch_stress[0]]
                .into_iter()
                .enumerate()
            {
                reconstructed[i] = reconstructed[i] + value / (s - r.laplace_frequency_per_s);
            }
        }
        close(reconstructed[0], displacement, 3e-11);
        close(reconstructed[1], s * displacement, 3e-11);
        close(reconstructed[2], stress, 3e-11);
    }
}

fn coupled_source(gyro: f64) -> ([f64; 4], [f64; 4], [f64; 4], Vec<MaxwellDampingBranch>) {
    // Plant s=-1+2i and q=[1,i] in a full real mechanical pencil.
    // The independent branch compensation fixes T at that point only;
    // its complete frequency derivative is still needed for the residue.
    let mass = [1.0, 0.0, 0.0, 1.0];
    let mut damping = [1.0, 0.5, 0.5, 3.0];
    let mut stiffness = [5.0, 2.5, 2.5, 5.0];
    let branches = vec![branch(0.8, 0.4, &[1.0, 0.5]), branch(0.3, 0.7, &[0.3, 1.0])];
    let pole = C64::new(-1.0, 2.0);
    for branch in &branches {
        let st = pole.scale(branch.relaxation_time_s);
        let response = (st / (C64::ONE + st)).scale(branch.stiffness);
        let dc = -response.im / pole.im;
        let dk = -response.re - pole.re * dc;
        for i in 0..2 {
            for j in 0..2 {
                damping[2 * i + j] += dc * branch.strain[i] * branch.strain[j];
                stiffness[2 * i + j] += dk * branch.strain[i] * branch.strain[j];
            }
        }
    }
    damping[1] -= gyro;
    damping[2] += gyro;
    for i in [0, 3] {
        damping[i] -= gyro / 2.0;
        stiffness[i] -= 2.5 * gyro;
    }
    (mass, damping, stiffness, branches)
}
fn pencil(
    mass: &[f64; 4],
    damping: &[f64; 4],
    stiffness: &[f64; 4],
    branches: &[MaxwellDampingBranch],
    s: C64,
) -> ([C64; 4], [C64; 4]) {
    let mut t = core::array::from_fn(|i| {
        s * s.scale(mass[i]) + s.scale(damping[i]) + C64::from_re(stiffness[i])
    });
    let mut derivative =
        core::array::from_fn(|i| s.scale(2.0 * mass[i]) + C64::from_re(damping[i]));
    for branch in branches {
        let st = s.scale(branch.relaxation_time_s);
        let denominator = C64::ONE + st;
        let weight = st / denominator;
        let slope = C64::from_re(branch.relaxation_time_s) / (denominator * denominator);
        for i in 0..2 {
            for j in 0..2 {
                let coefficient = branch.stiffness * branch.strain[i] * branch.strain[j];
                t[2 * i + j] = t[2 * i + j] + weight.scale(coefficient);
                derivative[2 * i + j] = derivative[2 * i + j] + slope.scale(coefficient);
            }
        }
    }
    (t, derivative)
}
fn adjugate_action(t: &[C64; 4], f: &[C64; 2]) -> [C64; 2] {
    [t[3] * f[0] - t[1] * f[1], t[0] * f[1] - t[2] * f[0]]
}

#[test]
fn nonsymmetric_complex_residue_matches_adjugate_and_nearby_forced_response() {
    let pole = C64::new(-1.0, 2.0);
    let force = [C64::new(0.8, 0.3), C64::new(-0.4, 0.2)];
    for gyro in [-0.3, 0.3] {
        let (mass, damping, stiffness, branches) = coupled_source(gyro);
        let prepared = PreparedMaxwellModes::new(
            &mass,
            &damping,
            &stiffness,
            2,
            &branches,
            options(3.0),
            || false,
        )
        .unwrap();
        let (t, dt) = pencil(&mass, &damping, &stiffness, &branches, pole);
        let left_seed = [-t[2], t[0]];
        let residue = prepared
            .refine_residue(
                C64::new(-0.96, 2.04),
                &[C64::ONE, C64::new(0.0, 0.9)],
                &left_seed,
                1e-8,
                || false,
            )
            .unwrap();
        let applied = residue.force_residue(&force, || false).unwrap();
        let determinant_slope = dt[0] * t[3] + t[0] * dt[3] - dt[1] * t[2] - t[1] * dt[2];
        let expected = adjugate_action(&t, &force).map(|x| x / determinant_slope);
        for (i, &value) in expected.iter().enumerate() {
            close(applied.displacement[i], value, 2e-11);
            close(applied.velocity[i], pole * value, 2e-11);
        }
        // Reciprocal/real-mode substitutions are observably wrong here.
        let right = &residue.right_mode().displacement;
        let left = residue.transpose_displacement();
        assert!((left[1] / left[0] - right[1] / right[0]).abs() > 0.1);
        for (branch, &actual) in branches.iter().zip(&applied.branch_stress) {
            let strain = expected[0].scale(branch.strain[0]) + expected[1].scale(branch.strain[1]);
            let st = pole.scale(branch.relaxation_time_s);
            close(
                actual,
                (st / (C64::ONE + st)).scale(branch.stiffness) * strain,
                3e-11,
            );
        }
        let mut previous = f64::INFINITY;
        for delta in [1e-2, 1e-3, 1e-4] {
            let s = pole + C64::from_re(delta);
            let (near, _) = pencil(&mass, &damping, &stiffness, &branches, s);
            let determinant = near[0] * near[3] - near[1] * near[2];
            let forced = adjugate_action(&near, &force).map(|x| (x / determinant).scale(delta));
            let error = (forced[0] - expected[0]).abs() + (forced[1] - expected[1]).abs();
            assert!(
                error < previous * 0.13,
                "near-pole convergence: {error} after {previous}"
            );
            previous = error;
        }
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One physical congruence/rate map checks force, stress and pole derivatives.
fn original_force_and_stress_residues_survive_coordinate_and_rate_scaling() {
    let (mass, damping, stiffness, branches) = coupled_source(-0.2);
    let force = [C64::new(0.8, 0.3), C64::new(-0.4, 0.2)];
    let pole = C64::new(-1.0, 2.0);
    let (t, dt) = pencil(&mass, &damping, &stiffness, &branches, pole);
    let derivative = dt[0] * t[3] + t[0] * dt[3] - dt[1] * t[2] - t[1] * dt[2];
    let expected = adjugate_action(&t, &force).map(|x| x / derivative);
    let left_seed = [-t[2], t[0]];
    for (coordinates, rate) in [
        ([2.0_f64.powi(300), 2.0_f64.powi(-300)], 1.0),
        ([1.0, 1.0], 2.0_f64.powi(400)),
        ([1.0, 1.0], 2.0_f64.powi(-400)),
    ] {
        let mut matrices = [[0.0; 4]; 3];
        for (order, a) in [mass, damping, stiffness].iter().enumerate() {
            let factor = match order {
                0 => 1.0,
                1 => rate,
                _ => rate * rate,
            };
            for i in 0..2 {
                for j in 0..2 {
                    matrices[order][2 * i + j] =
                        a[2 * i + j] * coordinates[i] * coordinates[j] * factor;
                }
            }
        }
        let transformed: Vec<_> = branches
            .iter()
            .map(|b| {
                branch(
                    (b.stiffness * rate) * rate,
                    b.relaxation_time_s / rate,
                    &[b.strain[0] * coordinates[0], b.strain[1] * coordinates[1]],
                )
            })
            .collect();
        let prepared = PreparedMaxwellModes::new(
            &matrices[0],
            &matrices[1],
            &matrices[2],
            2,
            &transformed,
            options(3.0 * rate),
            || false,
        )
        .unwrap();
        let right = [
            C64::from_re(1.0 / coordinates[0]),
            C64::new(0.0, 1.0 / coordinates[1]),
        ];
        let left = core::array::from_fn::<_, 2, _>(|i| left_seed[i].scale(1.0 / coordinates[i]));
        let f = core::array::from_fn::<_, 2, _>(|i| force[i].scale(coordinates[i]));
        let residue = prepared
            .refine_residue(pole.scale(rate), &right, &left, 1e-8, || false)
            .unwrap();
        let response = residue.force_residue(&f, || false).unwrap();
        // Uniform physical rate scaling: C -> a*C, K,k -> a^2*(K,k),
        // tau -> tau/a. Thus d(s)/dln(a)=s, even in extreme source units.
        let stiffness_direction = matrices[2].map(|x| 2.0 * x);
        let branch_directions: Vec<_> = transformed
            .iter()
            .map(|branch| MaxwellBranchDirection {
                stiffness: 2.0 * branch.stiffness,
                log_relaxation_time: -1.0,
                strain: &[0.0, 0.0],
            })
            .collect();
        close(
            residue
                .pole_derivative(
                    &MaxwellPoleDirection {
                        mass: &[0.0; 4],
                        damping: &matrices[1],
                        relaxed_stiffness: &stiffness_direction,
                        branches: &branch_directions,
                    },
                    || false,
                )
                .unwrap(),
            pole.scale(rate),
            4e-11,
        );
        for i in 0..2 {
            close(
                response.displacement[i].scale(coordinates[i]).scale(rate),
                expected[i],
                3e-11,
            );
            close(
                response.velocity[i].scale(coordinates[i]),
                pole * expected[i],
                3e-11,
            );
        }
        for (branch, &actual) in branches.iter().zip(&response.branch_stress) {
            let strain = expected[0].scale(branch.strain[0]) + expected[1].scale(branch.strain[1]);
            let st = pole.scale(branch.relaxation_time_s);
            close(
                actual.scale(1.0 / rate),
                (st / (C64::ONE + st)).scale(branch.stiffness) * strain,
                4e-11,
            );
        }
    }
}

#[test]
fn weak_relaxation_retains_large_signed_residues() {
    let branches = [branch(1e-20, 1.0, &[1.0])];
    let prepared =
        PreparedMaxwellModes::new(&[1.0], &[0.0], &[0.0], 1, &branches, options(1e-20), || {
            false
        })
        .unwrap();
    let residue = prepared
        .refine_residue(
            C64::from_re(-1.03e-20),
            &[C64::ONE],
            &[C64::ONE],
            1e-8,
            || false,
        )
        .unwrap();
    let response = residue.force_residue(&[C64::ONE], || false).unwrap();
    close(
        response.laplace_frequency_per_s,
        C64::from_re(-1e-20),
        6e-13,
    );
    close(response.displacement[0], C64::from_re(-1e20), 7e-13);
    close(response.velocity[0], C64::ONE, 7e-13);
    close(response.branch_stress[0], C64::from_re(1e-20), 8e-13);
}

#[test]
fn nearly_defective_derivative_cancellation_and_alignment_refuse() {
    let epsilon = 1e-3;
    let stiffness = [1.0 - epsilon * epsilon];
    let prepared =
        PreparedMaxwellModes::new(&[1.0], &[2.0], &stiffness, 1, &[], options(1.0), || false)
            .unwrap();
    assert!(matches!(prepared.refine_residue(
        C64::from_re(-1.0+epsilon),&[C64::ONE],&[C64::ONE],1e-2,|| false,
    ),Err(MaxwellModeError::ResidueNormalization { cancellation_ratio,.. }) if cancellation_ratio<1e-2));

    // Both derivative terms have the same sign; the source cancellation
    // ratio alone cannot detect nearly orthogonal right/transpose modes.
    let stiffness = [4.0, -1.0, -epsilon * epsilon, 4.0];
    let prepared = PreparedMaxwellModes::new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0; 4],
        &stiffness,
        2,
        &[],
        options(2.0),
        || false,
    )
    .unwrap();
    let right = [C64::ONE, C64::from_re(epsilon)];
    let left = [C64::from_re(epsilon), C64::ONE];
    assert!(matches!(prepared.refine_residue(
        C64::new(0.0,(4.0-epsilon).sqrt()),&right,&left,1e-2,|| false,
    ),Err(MaxwellModeError::ResidueNormalization { alignment,cancellation_ratio,.. }) if alignment<1e-2 && cancellation_ratio>0.9));
}

#[test]
fn paired_budgets_invalid_inputs_and_cancellation_publish_no_partial_residue() {
    let branches = [branch(9.6, 0.2, &[1.0])];
    for (bytes, work, quantity) in [
        (2000, 200_000_000, "paired numerical payload bytes"),
        (16 * 1024 * 1024, 50_000, "paired arithmetic work"),
    ] {
        let mut controls = options(5.0);
        controls.max_dense_bytes = bytes;
        controls.max_work = work;
        let prepared =
            PreparedMaxwellModes::new(&[4.0], &[8.0], &[6.4], 1, &branches, controls, || false)
                .unwrap();
        assert!(
            prepared
                .refine(C64::from_re(-1.03), &[C64::ONE], || false)
                .is_ok()
        );
        let mut polls = 0;
        assert!(
            matches!(prepared.refine_residue(C64::from_re(-1.03),&[C64::ONE],&[C64::ONE],1e-8,|| {polls+=1;false}),Err(MaxwellModeError::Budget { quantity:q,.. }) if q==quantity)
        );
        assert_eq!(polls, 1);
    }
    let prepared =
        PreparedMaxwellModes::new(&[4.0], &[8.0], &[6.4], 1, &branches, options(5.0), || false)
            .unwrap();
    for floor in [0.0, 1.0, -1.0, f64::NAN] {
        assert!(matches!(
            prepared.refine_residue(C64::from_re(-1.03), &[C64::ONE], &[C64::ONE], floor, || {
                false
            }),
            Err(MaxwellModeError::Invalid(_))
        ));
    }
    assert!(matches!(
        prepared.refine_residue(C64::from_re(-5.0), &[C64::ONE], &[C64::ONE], 1e-8, || false),
        Err(MaxwellModeError::MaterialPole { .. })
    ));
    let mut count = 0;
    let baseline = prepared
        .refine_residue(C64::from_re(-1.03), &[C64::ONE], &[C64::ONE], 1e-8, || {
            count += 1;
            false
        })
        .unwrap();
    for stop in [1, count / 2, count] {
        let mut calls = 0;
        assert!(matches!(
            prepared.refine_residue(C64::from_re(-1.03), &[C64::ONE], &[C64::ONE], 1e-8, || {
                calls += 1;
                calls >= stop
            }),
            Err(MaxwellModeError::Cancelled)
        ));
    }
    let force = [C64::new(0.7, -0.2)];
    let expected = baseline.force_residue(&force, || false).unwrap();
    assert!(matches!(
        baseline.force_residue(&[], || false),
        Err(MaxwellModeError::Invalid(_))
    ));
    assert!(matches!(
        baseline.force_residue(&[C64::new(f64::NAN, 0.0)], || false),
        Err(MaxwellModeError::Invalid(_))
    ));
    let mut calls = 0;
    assert!(matches!(
        baseline.force_residue(&force, || {
            calls += 1;
            calls >= 6
        }),
        Err(MaxwellModeError::Cancelled)
    ));
    let replay = prepared
        .refine_residue(C64::from_re(-1.03), &[C64::ONE], &[C64::ONE], 1e-8, || {
            false
        })
        .unwrap()
        .force_residue(&force, || false)
        .unwrap();
    assert_eq!(replay.displacement, expected.displacement);
    assert_eq!(replay.velocity, expected.velocity);
    assert_eq!(replay.branch_stress, expected.branch_stress);
}

#[test]
fn scalar_pole_derivatives_retain_every_original_source_parameter() {
    let branches = [branch(9.6, 0.2, &[1.0])];
    let prepared =
        PreparedMaxwellModes::new(&[4.0], &[8.0], &[6.4], 1, &branches, options(5.0), || false)
            .unwrap();
    let residue = prepared
        .refine_residue(C64::from_re(-2.03), &[C64::ONE], &[C64::ONE], 1e-8, || {
            false
        })
        .unwrap();
    // At s=-2: T'=-8/3. Differentiate m*s²+c*s+K+k*s*tau*b²/(1+s*tau)
    // directly in m,c,K,k,ln(tau),b, without any eigenvector calculation.
    for (index, expected) in [1.5, -0.75, 0.375, -0.25, -4.0, -4.8]
        .into_iter()
        .enumerate()
    {
        let mut d = [0.0; 6];
        d[index] = 1.0;
        let dstrain = [d[5]];
        let dbranch = [MaxwellBranchDirection {
            stiffness: d[3],
            log_relaxation_time: d[4],
            strain: &dstrain,
        }];
        let dm = [d[0]];
        let dc = [d[1]];
        let dk = [d[2]];
        let direction = MaxwellPoleDirection {
            mass: &dm,
            damping: &dc,
            relaxed_stiffness: &dk,
            branches: &dbranch,
        };
        close(
            residue.pole_derivative(&direction, || false).unwrap(),
            C64::from_re(expected),
            5e-12,
        );
        assert!(matches!(
            residue.pole_derivative(&direction, || true),
            Err(MaxwellModeError::Cancelled)
        ));
    }
    let missing = MaxwellPoleDirection {
        mass: &[],
        damping: &[0.0],
        relaxed_stiffness: &[0.0],
        branches: &[],
    };
    assert!(matches!(
        residue.pole_derivative(&missing, || false),
        Err(MaxwellModeError::Invalid(_))
    ));
    let bad_branch = [MaxwellBranchDirection {
        stiffness: 0.0,
        log_relaxation_time: f64::NAN,
        strain: &[0.0],
    }];
    let invalid = MaxwellPoleDirection {
        mass: &[0.0],
        damping: &[0.0],
        relaxed_stiffness: &[0.0],
        branches: &bad_branch,
    };
    assert!(matches!(
        residue.pole_derivative(&invalid, || false),
        Err(MaxwellModeError::Invalid(_))
    ));
}

#[test]
fn coupled_directional_pole_derivative_matches_original_source_perturbations() {
    let (mass, damping, stiffness, branches) = coupled_source(0.2);
    let pole = C64::new(-1.0, 2.0);
    let (t, _) = pencil(&mass, &damping, &stiffness, &branches, pole);
    let right = [C64::ONE, C64::new(0.0, 1.0)];
    let left = [-t[2], t[0]];
    let prepared = PreparedMaxwellModes::new(
        &mass,
        &damping,
        &stiffness,
        2,
        &branches,
        options(3.0),
        || false,
    )
    .unwrap();
    let residue = prepared
        .refine_residue(pole, &right, &left, 1e-8, || false)
        .unwrap();
    let dm = [0.2, 0.04, 0.04, -0.1];
    let dc = [0.07, -0.04, 0.08, -0.05];
    let dk = [-0.2, 0.03, -0.06, 0.1];
    let db = [
        MaxwellBranchDirection {
            stiffness: 0.15,
            log_relaxation_time: 0.08,
            strain: &[0.12, -0.09],
        },
        MaxwellBranchDirection {
            stiffness: -0.05,
            log_relaxation_time: -0.12,
            strain: &[-0.02, 0.07],
        },
    ];
    let direction = MaxwellPoleDirection {
        mass: &dm,
        damping: &dc,
        relaxed_stiffness: &dk,
        branches: &db,
    };
    let derivative = residue.pole_derivative(&direction, || false).unwrap();
    let perturbed = |h: f64| {
        let m: [f64; 4] = core::array::from_fn(|i| mass[i] + h * dm[i]);
        let c: [f64; 4] = core::array::from_fn(|i| damping[i] + h * dc[i]);
        let k: [f64; 4] = core::array::from_fn(|i| stiffness[i] + h * dk[i]);
        let memory: Vec<_> = branches
            .iter()
            .zip(db)
            .map(|(b, d)| {
                branch(
                    b.stiffness + h * d.stiffness,
                    b.relaxation_time_s * (h * d.log_relaxation_time).exp(),
                    &[b.strain[0] + h * d.strain[0], b.strain[1] + h * d.strain[1]],
                )
            })
            .collect();
        PreparedMaxwellModes::new(&m, &c, &k, 2, &memory, options(3.0), || false)
            .unwrap()
            .refine(pole + derivative.scale(h), &right, || false)
            .unwrap()
            .laplace_frequency_per_s
    };
    let mut previous = f64::INFINITY;
    for h in [0.002, 0.001, 0.0005] {
        let difference = (perturbed(h) - perturbed(-h)).scale(0.5 / h);
        let error = (difference - derivative).abs();
        assert!(
            error < previous * 0.3,
            "central derivative convergence: {error} after {previous}"
        );
        previous = error;
    }
    assert!(
        previous < 5e-7 * derivative.abs(),
        "derivative {derivative:?}, final error {previous}, relative {}",
        previous / derivative.abs()
    );
    let nonsymmetric_mass = [0.0, 1.0, 0.0, 0.0];
    let bad = MaxwellPoleDirection {
        mass: &nonsymmetric_mass,
        ..direction
    };
    assert!(matches!(
        residue.pole_derivative(&bad, || false),
        Err(MaxwellModeError::Invalid(_))
    ));
}
