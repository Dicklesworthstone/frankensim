use super::*;
use crate::rotordynamics::RotorBearingModel;
use crate::rotordynamics::periodic::periodic_floquet;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn spectrum_close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let mut used = vec![false; actual.len()];
    for target in expected {
        let (index, error) = actual
            .iter()
            .enumerate()
            .filter(|(i, _)| !used[*i])
            .map(|(i, value)| (i, (*value - *target).abs()))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        assert!(
            error <= tolerance * (1.0 + target.abs()),
            "error={error:e}, target={target:?}"
        );
        used[index] = true;
    }
}

fn arm(strain: &[f64], stiffness: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: strain.to_vec(),
        stiffness,
        relaxation_time_s: tau,
    }
}

fn oscillator(stiffness: f64, damping: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[1.0], &[damping], &[stiffness], &[0.0], || false).unwrap()
}

#[test]
fn three_exact_physical_eigenmotions_survive_every_segment_boundary() {
    // det(sI-A) = s^3+6s^2+11s+6 = (s+1)(s+2)(s+3).
    // M=1, C=2, K_relaxed=1.5, k=1.5, tau=1/4, b=1.
    let model = oscillator(1.5, 2.0);
    let branches = [arm(&[1.0], 1.5, 0.25)];
    let cycle = [
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.2,
        },
        RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.5,
        },
    ];
    let prepared = prepare_maxwell_cycle(
        7.0,
        &cycle,
        &branches,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let report = prepared.floquet_report();
    spectrum_close(
        &report.multipliers,
        &[-1.0_f64, -2.0, -3.0].map(|rate| C64::from_re((0.7 * rate).exp())),
        2e-12,
    );
    assert_eq!(report.dofs, 1);
    assert_eq!(report.memory_states, 1);
    assert_eq!(report.monodromy.len(), 9);
    close(report.growth_rate_per_s, -1.0, 2e-12);
    assert!(report.maximum_mass_residual < 1e-12);
    assert!(report.maximum_memory_residual < 1e-12);
    for rate in [-1.0, -2.0, -3.0] {
        let stress = -(rate * rate + 2.0 * rate + 1.5);
        let initial = RotorMaxwellState {
            displacement: vec![1.0],
            velocity: vec![rate],
            branch_stress: vec![stress],
        };
        let trajectory = prepared.advance(&initial, 1e-11, || false).unwrap();
        assert!(trajectory.relative_map_residual < 1e-11);
        for (&time, state) in trajectory
            .boundary_times_seconds
            .iter()
            .zip(&trajectory.boundary_states)
        {
            let amplitude = (rate * time).exp();
            close(state.displacement[0], amplitude, 3e-12);
            close(state.velocity[0], rate * amplitude, 3e-12);
            close(state.branch_stress[0], stress * amplitude, 3e-12);
        }
    }
    let one = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.7,
    }];
    let unsplit = periodic_maxwell_floquet(
        7.0,
        &one,
        &branches,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    for (&actual, &expected) in report.monodromy.iter().zip(&unsplit.monodromy) {
        close(actual, expected, 3e-12);
    }
    assert_eq!(
        *report,
        periodic_maxwell_floquet(
            7.0,
            &cycle,
            &branches,
            1e-12,
            RotorFloquetBudget::default(),
            || false
        )
        .unwrap()
    );
}

#[test]
fn equal_rate_dark_combination_and_uncoupled_stress_keep_their_own_decay() {
    let model = oscillator(4.0, 0.2);
    let branches = [
        arm(&[1.0], 2.0, 0.4),
        arm(&[1.0], 5.0, 0.4),
        arm(&[0.0], 9.0, 0.7),
    ];
    let cycle = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.3,
    }];
    let prepared = prepare_maxwell_cycle(
        0.0,
        &cycle,
        &branches,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let report = prepared.floquet_report();
    assert_eq!(report.multipliers.len(), 5);
    assert_eq!(report.memory_states, 3);
    for tau in [0.4, 0.7] {
        let expected = (-0.3_f64 / tau).exp();
        assert!(
            report
                .multipliers
                .iter()
                .any(|value| (*value - C64::from_re(expected)).abs() < 1e-12)
        );
    }
    // Opposing physical stresses cancel mechanically despite different k.
    let initial = RotorMaxwellState {
        displacement: vec![0.0],
        velocity: vec![0.0],
        branch_stress: vec![3.0, -3.0, 2.0],
    };
    let trajectory = prepared.advance(&initial, 1e-11, || false).unwrap();
    let final_state = trajectory.boundary_states.last().unwrap();
    close(final_state.displacement[0], 0.0, 2e-13);
    close(final_state.velocity[0], 0.0, 2e-13);
    close(final_state.branch_stress[0], 3.0 * (-0.75_f64).exp(), 2e-12);
    close(
        final_state.branch_stress[1],
        -3.0 * (-0.75_f64).exp(),
        2e-12,
    );
    close(
        final_state.branch_stress[2],
        2.0 * (-0.3_f64 / 0.7).exp(),
        2e-12,
    );
}

// Independent time-domain oracle: integrate the three ORIGINAL scalar ODEs
// with classical RK4, not the production exponential or any reduced damping.
fn rk4(
    mut state: [f64; 3],
    stiffness: f64,
    damping: f64,
    branch: &MaxwellDampingBranch,
    duration: f64,
    steps: usize,
) -> [f64; 3] {
    let derivative = |x: [f64; 3]| {
        [
            x[1],
            -stiffness * x[0] - damping * x[1] - branch.strain[0] * x[2],
            branch.stiffness * branch.strain[0] * x[1] - x[2] / branch.relaxation_time_s,
        ]
    };
    let h = duration / steps as f64;
    for _ in 0..steps {
        let a = derivative(state);
        let b = derivative(std::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = derivative(std::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = derivative(std::array::from_fn(|i| state[i] + h * c[i]));
        for i in 0..3 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn original_memory_survives_modulation_and_reveals_parametric_growth() {
    let soft = oscillator(1.0, 0.08);
    let stiff = oscillator(9.0, 0.08);
    let branches = [arm(&[1.0], 0.15, 0.6)];
    for model in [&soft, &stiff] {
        let frozen = model
            .at_speed_with_maxwell_memory(0.0, &branches, || false)
            .unwrap()
            .poles(|| false)
            .unwrap();
        assert!(frozen.poles_per_s.iter().all(|value| value.re < 0.0));
    }
    let h1 = std::f64::consts::FRAC_PI_2;
    let h2 = std::f64::consts::PI / 6.0;
    let cycle = [
        RotorPeriodicSegment {
            model: &soft,
            duration_seconds: h1,
        },
        RotorPeriodicSegment {
            model: &stiff,
            duration_seconds: h2,
        },
    ];
    let prepared = prepare_maxwell_cycle(
        0.0,
        &cycle,
        &branches,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let report = prepared.floquet_report();
    assert!(report.spectral_radius > 2.0);
    assert!(report.growth_rate_per_s > 0.0);
    for column in 0..3 {
        let basis = std::array::from_fn(|i| if i == column { 1.0 } else { 0.0 });
        let reference = |steps| {
            let first = rk4(basis, 1.0, 0.08, &branches[0], h1, steps);
            let last = rk4(first, 9.0, 0.08, &branches[0], h2, steps);
            (first, last)
        };
        let (first, last) = reference(1024);
        let (_, coarse) = reference(512);
        for i in 0..3 {
            close(coarse[i], last[i], 3e-10);
            close(report.monodromy[i * 3 + column], last[i], 3e-10);
        }
        let trajectory = prepared
            .advance(
                &RotorMaxwellState {
                    displacement: vec![basis[0]],
                    velocity: vec![basis[1]],
                    branch_stress: vec![basis[2]],
                },
                1e-10,
                || false,
            )
            .unwrap();
        for (state, expected) in trajectory.boundary_states[1..].iter().zip([first, last]) {
            close(state.displacement[0], expected[0], 3e-10);
            close(state.velocity[0], expected[1], 3e-10);
            close(state.branch_stress[0], expected[2], 3e-10);
        }
    }
    // Reversing chronological order changes the physical map even though its
    // nonzero spectrum is phase-origin invariant.
    let reverse = [cycle[1], cycle[0]];
    let reversed = periodic_maxwell_floquet(
        0.0,
        &reverse,
        &branches,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    assert!((report.monodromy[0] - reversed.monodromy[0]).abs() > 1.0);
    spectrum_close(&report.multipliers, &reversed.multipliers, 3e-11);
}

#[test]
fn full_mass_signed_gyro_and_original_strain_are_covariant_under_units() {
    let mass = [2.0, 0.4, 0.4, 1.5];
    let damping = [0.12, -0.02, 0.07, 0.18];
    let k1 = [4.0, 0.3, -0.15, 7.0];
    let k2 = [9.0, -0.5, 0.2, 3.0];
    let gyro = [0.0, 0.8, -0.8, 0.0];
    let branches = [arm(&[1.2, -0.4], 3.5, 0.7), arm(&[-0.3, 0.8], 2.0, 0.2)];
    let units = [1e5, 1e-4];
    let congruence = |matrix: &[f64; 4]| -> [f64; 4] {
        std::array::from_fn(|i| matrix[i] * units[i / 2] * units[i % 2])
    };
    let build = |k: &[f64]| RotorBearingModel::new(2, &mass, &damping, k, &gyro, || false).unwrap();
    let transformed = |k: &[f64; 4]| {
        RotorBearingModel::new(
            2,
            &congruence(&mass),
            &congruence(&damping),
            &congruence(k),
            &congruence(&gyro),
            || false,
        )
        .unwrap()
    };
    let first = build(&k1);
    let second = build(&k2);
    let first_scaled = transformed(&k1);
    let second_scaled = transformed(&k2);
    let scaled_branches: Vec<_> = branches
        .iter()
        .map(|branch| {
            arm(
                &[units[0] * branch.strain[0], units[1] * branch.strain[1]],
                branch.stiffness,
                branch.relaxation_time_s,
            )
        })
        .collect();
    let cycle = [
        RotorPeriodicSegment {
            model: &first,
            duration_seconds: 0.17,
        },
        RotorPeriodicSegment {
            model: &second,
            duration_seconds: 0.29,
        },
    ];
    let scaled_cycle = [
        RotorPeriodicSegment {
            model: &first_scaled,
            duration_seconds: 0.17,
        },
        RotorPeriodicSegment {
            model: &second_scaled,
            duration_seconds: 0.29,
        },
    ];
    let transform = [units[0], units[1], units[0], units[1], 1.0, 1.0];
    for spin in [-2.7, 2.7] {
        let original = periodic_maxwell_floquet(
            spin,
            &cycle,
            &branches,
            1e-11,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        let scaled = periodic_maxwell_floquet(
            spin,
            &scaled_cycle,
            &scaled_branches,
            1e-11,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        for i in 0..6 {
            for j in 0..6 {
                close(
                    scaled.monodromy[i * 6 + j],
                    original.monodromy[i * 6 + j] * transform[j] / transform[i],
                    2e-11,
                );
            }
        }
        spectrum_close(&original.multipliers, &scaled.multipliers, 3e-11);
        assert!(scaled.maximum_mass_residual < 1e-11);
        assert!(scaled.maximum_memory_residual < 1e-11);
    }
}

#[test]
fn empty_memory_agrees_with_original_owner_and_dense_budget_counts_all_states() {
    let model = oscillator(4.0, 0.2);
    let cycle = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.1,
    }];
    let original =
        periodic_floquet(-3.0, &cycle, 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    let empty = periodic_maxwell_floquet(
        -3.0,
        &cycle,
        &[],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    assert_eq!(original.monodromy, empty.monodromy);
    assert_eq!(original.multipliers, empty.multipliers);
    assert_eq!(empty.maximum_memory_residual, 0.0);
    // 66 original states exceed the mechanical-only 64-state ceiling.
    // Every uncoupled arm nevertheless remains a separate physical channel.
    let branches: Vec<_> = (0..MAX_ROTOR_MAXWELL_BRANCHES)
        .map(|i| arm(&[0.0], 1.0 + i as f64, 0.5 + 0.02 * i as f64))
        .collect();
    let full = periodic_maxwell_floquet(
        0.0,
        &cycle,
        &branches,
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    assert_eq!(full.multipliers.len(), 66);
    assert_eq!(full.monodromy.len(), 66 * 66);
    for (j, branch) in branches.iter().enumerate() {
        close(
            full.monodromy[(2 + j) * 66 + 2 + j],
            (-0.1 / branch.relaxation_time_s).exp(),
            2e-13,
        );
    }
    let tiny = RotorFloquetBudget {
        max_multiply_terms: 1_000_000,
        ..RotorFloquetBudget::default()
    };
    assert!(periodic_maxwell_floquet(0.0, &cycle, &[], 1e-12, tiny, || false).is_ok());
    assert!(matches!(
        periodic_maxwell_floquet(0.0, &cycle, &branches, 1e-12, tiny, || false),
        Err(RotorError::PeriodicBudget(_))
    ));
    let bytes = RotorFloquetBudget {
        max_dense_bytes: 1,
        ..RotorFloquetBudget::default()
    };
    assert!(matches!(
        prepare_maxwell_cycle(0.0, &cycle, &branches, 1e-12, bytes, || false),
        Err(RotorError::PeriodicBudget(_))
    ));
}

#[test]
fn cancellation_and_invalid_original_states_publish_no_partial_result() {
    let model = oscillator(4.0, 0.2);
    let cycle = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.2,
    }];
    let branches = [arm(&[1.0], 2.0, 0.5)];
    let budget = RotorFloquetBudget::default();
    assert!(matches!(
        prepare_maxwell_cycle(0.0, &cycle, &branches, 1e-12, budget, || true),
        Err(RotorError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        prepare_maxwell_cycle(0.0, &cycle, &branches, 1e-12, budget, || {
            polls += 1;
            polls == 20
        }),
        Err(RotorError::Cancelled)
    ));
    for bad in [
        arm(&[], 2.0, 0.5),
        arm(&[1.0], 2.0, 0.0),
        arm(&[1.0], -2.0, 0.5),
        arm(&[f64::NAN], 2.0, 0.5),
    ] {
        assert!(prepare_maxwell_cycle(0.0, &cycle, &[bad], 1e-12, budget, || false).is_err());
    }
    let overfull = vec![arm(&[0.0], 1.0, 1.0); MAX_ROTOR_MAXWELL_BRANCHES + 1];
    assert!(prepare_maxwell_cycle(0.0, &cycle, &overfull, 1e-12, budget, || false).is_err());
    let other_mass = RotorBearingModel::new(1, &[2.0], &[0.2], &[4.0], &[0.0], || false).unwrap();
    let changed = [
        cycle[0],
        RotorPeriodicSegment {
            model: &other_mass,
            duration_seconds: 0.2,
        },
    ];
    assert!(prepare_maxwell_cycle(0.0, &changed, &branches, 1e-12, budget, || false).is_err());
    let prepared = prepare_maxwell_cycle(0.0, &cycle, &branches, 1e-12, budget, || false).unwrap();
    for invalid in [
        RotorMaxwellState {
            displacement: vec![0.0],
            velocity: vec![0.0],
            branch_stress: vec![],
        },
        RotorMaxwellState {
            displacement: vec![0.0],
            velocity: vec![0.0],
            branch_stress: vec![f64::NAN],
        },
    ] {
        assert!(prepared.advance(&invalid, 1e-10, || false).is_err());
    }
    let valid = RotorMaxwellState {
        displacement: vec![0.0],
        velocity: vec![0.0],
        branch_stress: vec![1.0],
    };
    assert!(matches!(
        prepared.advance(&valid, 1e-10, || true),
        Err(RotorError::Cancelled)
    ));
    assert!(prepared.advance(&valid, 0.0, || false).is_err());
    // A vanished material direction cannot become an infinite-decay claim.
    let unresolved = [arm(&[0.0], 1.0, 0.0001)];
    assert!(periodic_maxwell_floquet(0.0, &cycle, &unresolved, 1e-12, budget, || false).is_err());
}
