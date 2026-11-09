//! G1/G3 original Maxwell clocks, continuous memory and independently observed heat.
use fs_mbd::rotordynamics::periodic::power::RotorCyclePowerError;
use fs_mbd::rotordynamics::periodic::{
    RotorFloquetBudget, RotorFourierForce, RotorHarmonicForce, RotorMaxwellKinetics,
    RotorMaxwellState, RotorPeriodicForcing, RotorPeriodicSegment, prepare_maxwell_cycle,
    prepare_maxwell_forced_power_cycle, prepare_maxwell_fourier_power_cycle,
    prepare_maxwell_harmonic_power_cycle,
};
use fs_mbd::rotordynamics::relaxation::MaxwellDampingBranch;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::f64::consts::TAU;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn model(mass: f64, damping: f64, stiffness: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[mass], &[damping], &[stiffness], &[0.0], || false).unwrap()
}

fn arm(strain: f64, stiffness: f64, tau: f64) -> MaxwellDampingBranch {
    MaxwellDampingBranch {
        strain: vec![strain],
        stiffness,
        relaxation_time_s: tau,
    }
}

fn state(values: &[f64; 4]) -> RotorMaxwellState {
    RotorMaxwellState {
        displacement: vec![values[0]],
        velocity: vec![values[1]],
        branch_stress: vec![values[2], values[3]],
    }
}

fn budget() -> RotorFloquetBudget {
    RotorFloquetBudget {
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 500_000_000,
        ..RotorFloquetBudget::default()
    }
}

#[test]
fn each_original_clock_controls_dark_memory_and_heat_without_compounding_or_switch_work() {
    let body = model(1.0, 0.3, 2.0);
    let branches = [arm(0.0, 2.0, 0.25), arm(0.0, 3.0, 0.5)];
    let times = [0.15, 0.25, 0.3];
    let segments = times.map(|duration_seconds| RotorPeriodicSegment {
        model: &body,
        duration_seconds,
    });
    let scales = [vec![2.0, 0.25], vec![0.5, 4.0], vec![3.0, 0.8]];
    let clocks = RotorMaxwellKinetics {
        relaxation_time_scales: &scales,
    };
    let prepared = clocks
        .prepare_power_cycle(
            7.0,
            &segments,
            &branches,
            RotorPeriodicForcing::Homogeneous,
            1e-12,
            budget(),
            || false,
        )
        .unwrap();
    let initial = state(&[0.0, 0.0, 4.0, -3.0]);
    let (trajectory, power) = prepared.advance(&initial, 1e-10, 1e-10, || false).unwrap();
    // These times are authored from the original (.25,.5) s references on
    // EVERY row, not from the already shifted times of the preceding row.
    let actual_times = [[0.5, 0.125], [0.125, 2.0], [0.75, 0.4]];
    let mut stress = [4.0_f64, -3.0_f64];
    let mut accumulated_rate = [0.0; 2];
    let mut total_heat = [0.0; 2];
    for segment in 0..3 {
        for branch in 0..2 {
            let rate_time: f64 = times[segment] / actual_times[segment][branch];
            let stored = stress[branch].powi(2) / (2.0 * branches[branch].stiffness);
            let heat = stored * -(-2.0 * rate_time).exp_m1();
            stress[branch] *= (-rate_time).exp();
            accumulated_rate[branch] += rate_time;
            total_heat[branch] += heat;
            near(
                trajectory.boundary_states[segment + 1].branch_stress[branch],
                stress[branch],
                5e-12,
            );
            near(power.segments[segment].branch_heat_j[branch], heat, 5e-12);
            near(
                power.segments[segment].branch_elastic_energy_change_j[branch],
                -heat,
                5e-12,
            );
        }
        assert_eq!(power.segments[segment].stiffness_switch_work_j, 0.0);
    }
    for branch in 0..2 {
        near(power.branch_heat_j[branch], total_heat[branch], 5e-12);
        near(
            power.branch_elastic_energy_change_j[branch],
            -total_heat[branch],
            5e-12,
        );
    }
    assert_eq!(power.applied_work_j, 0.0);
    assert_eq!(power.velocity_force_work_j, 0.0);
    assert_eq!(power.kinetic_energy_change_j, 0.0);
    assert_eq!(power.relaxed_elastic_energy_change_j, 0.0);
    assert_eq!(power.stiffness_switch_work_j, 0.0);
    near(power.balance_residual_j, 0.0, 1e-10);
    let floquet = clocks
        .floquet(7.0, &segments, &branches, 1e-12, budget(), || false)
        .unwrap();
    assert_eq!(floquet.memory_states, 2);
    for branch in 0..2 {
        near(
            floquet.monodromy[(2 + branch) * 4 + 2 + branch],
            (-accumulated_rate[branch]).exp(),
            5e-12,
        );
    }
    assert_eq!(branches[0].relaxation_time_s, 0.25);
    assert_eq!(branches[1].relaxation_time_s, 0.5);
}

// Independent original physical equations, with no production flow, Gramian,
// eigensolver or fixed-point solver. The last four slots integrate actual work
// and the two distinct branch heats alongside [q,v,sigma_0,sigma_1].
fn rhs(t: f64, x: &[f64; 8], tau: [f64; 2]) -> [f64; 8] {
    let force = 0.3 + 0.7 * (TAU * t).cos() + 0.2 * (TAU * t).sin() + 0.1 * (3.0 * TAU * t).cos()
        - 0.05 * (3.0 * TAU * t).sin();
    [
        x[1],
        (force - 0.08 * x[1] - 41.0 * x[0] - x[2] + 0.6 * x[3]) / 1.2,
        8.0 * x[1] - x[2] / tau[0],
        -2.1 * x[1] - x[3] / tau[1],
        force * x[1],
        0.08 * x[1] * x[1],
        x[2] * x[2] / (8.0 * tau[0]),
        x[3] * x[3] / (3.5 * tau[1]),
    ]
}

fn reference_segment(
    mut x: [f64; 8],
    start: f64,
    duration: f64,
    tau: [f64; 2],
    steps: usize,
) -> [f64; 8] {
    let h = duration / steps as f64;
    for step in 0..steps {
        let t = start + step as f64 * h;
        let a = rhs(t, &x, tau);
        let b = rhs(
            t + 0.5 * h,
            &std::array::from_fn(|i| x[i] + 0.5 * h * a[i]),
            tau,
        );
        let c = rhs(
            t + 0.5 * h,
            &std::array::from_fn(|i| x[i] + 0.5 * h * b[i]),
            tau,
        );
        let d = rhs(t + h, &std::array::from_fn(|i| x[i] + h * c[i]), tau);
        for i in 0..8 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn switched_clocks_change_periodic_motion_and_heat_with_the_original_fourier_phase() {
    let body = model(1.2, 0.08, 41.0);
    let branches = [arm(1.0, 8.0, 0.1), arm(-0.6, 3.5, 0.25)];
    let segments = [0.37, 0.63].map(|duration_seconds| RotorPeriodicSegment {
        model: &body,
        duration_seconds,
    });
    let scales = [vec![0.25, 3.0], vec![4.0, 0.4]];
    let clocks = RotorMaxwellKinetics {
        relaxation_time_scales: &scales,
    };
    let tones = [
        RotorHarmonicForce {
            cosine: &[0.7],
            sine: &[0.2],
        },
        RotorHarmonicForce {
            cosine: &[0.1],
            sine: &[-0.05],
        },
    ];
    let sources = [RotorFourierForce {
        steady: &[0.3],
        harmonics: &tones,
    }; 2];
    let prepared = clocks
        .prepare_power_cycle(
            0.0,
            &segments,
            &branches,
            RotorPeriodicForcing::Fourier {
                rows: &sources,
                orders: &[1, 3],
            },
            1e-12,
            budget(),
            || false,
        )
        .unwrap();
    let (orbit, power) = prepared.periodic_orbit(1e-9, 1e8, 1e-9, || false).unwrap();
    // Independent positive-loss dynamics settle from zero under the actual
    // clocks (.025,.75) then (.4,.1) s. Both pieces use the same global phase.
    let mut reference = [0.0; 8];
    let mut settled = false;
    for _ in 0..256 {
        let before = reference;
        reference[4..].fill(0.0);
        reference = reference_segment(reference, 0.0, 0.37, [0.025, 0.75], 1024);
        reference = reference_segment(reference, 0.37, 0.63, [0.4, 0.1], 1024);
        if (0..4).all(|i| (reference[i] - before[i]).abs() < 2e-12) {
            settled = true;
            break;
        }
    }
    assert!(
        settled,
        "independent original dynamics did not reach a periodic state"
    );
    reference[4..].fill(0.0);
    let first = reference_segment(reference, 0.0, 0.37, [0.025, 0.75], 4096);
    let last = reference_segment(first, 0.37, 0.63, [0.4, 0.1], 4096);
    for (actual, expected) in orbit
        .trajectory
        .boundary_states
        .iter()
        .zip([reference, first, last])
    {
        near(actual.displacement[0], expected[0], 3e-8);
        near(actual.velocity[0], expected[1], 3e-8);
        near(actual.branch_stress[0], expected[2], 3e-8);
        near(actual.branch_stress[1], expected[3], 3e-8);
    }
    near(power.applied_work_j, last[4], 3e-8);
    near(power.velocity_force_work_j, last[5], 3e-8);
    for branch in 0..2 {
        near(power.branch_heat_j[branch], last[6 + branch], 3e-8);
        near(
            power.segments[0].branch_heat_j[branch],
            first[6 + branch],
            3e-8,
        );
        near(
            power.segments[1].branch_heat_j[branch],
            last[6 + branch] - first[6 + branch],
            3e-8,
        );
        near(power.branch_elastic_energy_change_j[branch], 0.0, 3e-9);
    }
    assert_eq!(power.stiffness_switch_work_j, 0.0);
    near(power.balance_residual_j, 0.0, 1e-9);
    let unshifted = prepare_maxwell_fourier_power_cycle(
        0.0,
        &segments,
        &branches,
        &sources,
        &[1, 3],
        1e-12,
        budget(),
        || false,
    )
    .unwrap()
    .periodic_orbit(1e-9, 1e8, 1e-9, || false)
    .unwrap();
    let current = &orbit.trajectory.boundary_states[0];
    let original = &unshifted.0.trajectory.boundary_states[0];
    assert!((current.displacement[0] - original.displacement[0]).abs() > 1e-4);
    assert!((current.velocity[0] - original.velocity[0]).abs() > 1e-3);
    assert!(
        (power.branch_heat_j.iter().sum::<f64>() - unshifted.1.branch_heat_j.iter().sum::<f64>())
            .abs()
            > 1e-4
    );
}

#[test]
fn unit_clocks_reuse_original_held_harmonic_fourier_and_homogeneous_motion() {
    let body = model(1.2, 0.2, 2.4);
    let branches = [arm(1.0, 0.7, 0.3), arm(-0.6, 0.35, 0.8)];
    let segments = [0.2, 0.3].map(|duration_seconds| RotorPeriodicSegment {
        model: &body,
        duration_seconds,
    });
    let scales = [vec![1.0; 2], vec![1.0; 2]];
    let clocks = RotorMaxwellKinetics {
        relaxation_time_scales: &scales,
    };
    let initial = state(&[0.1, -0.4, 0.2, -0.3]);
    let homogeneous = clocks
        .prepare_cycle(
            -2.0,
            &segments,
            &branches,
            RotorPeriodicForcing::Homogeneous,
            1e-12,
            budget(),
            || false,
        )
        .unwrap();
    let original =
        prepare_maxwell_cycle(-2.0, &segments, &branches, 1e-12, budget(), || false).unwrap();
    assert_eq!(homogeneous.floquet_report(), original.floquet_report());
    assert_eq!(
        homogeneous.advance(&initial, 1e-9, || false).unwrap(),
        original.advance(&initial, 1e-9, || false).unwrap()
    );
    let held: [&[f64]; 2] = [&[0.7], &[-0.2]];
    let harmonics = [
        RotorHarmonicForce {
            cosine: &[0.4],
            sine: &[0.2],
        },
        RotorHarmonicForce {
            cosine: &[-0.3],
            sine: &[0.5],
        },
    ];
    let fourier = [
        RotorFourierForce {
            steady: held[0],
            harmonics: &harmonics,
        },
        RotorFourierForce {
            steady: held[1],
            harmonics: &harmonics,
        },
    ];
    for source in [
        RotorPeriodicForcing::Held(&held),
        RotorPeriodicForcing::Harmonic {
            rows: &harmonics,
            order: 1,
        },
        RotorPeriodicForcing::Fourier {
            rows: &fourier,
            orders: &[1, 3],
        },
    ] {
        let expected = match source {
            RotorPeriodicForcing::Held(rows) => prepare_maxwell_forced_power_cycle(
                -2.0,
                &segments,
                &branches,
                rows,
                1e-12,
                budget(),
                || false,
            ),
            RotorPeriodicForcing::Harmonic { rows, order } => prepare_maxwell_harmonic_power_cycle(
                -2.0,
                &segments,
                &branches,
                rows,
                order,
                1e-12,
                budget(),
                || false,
            ),
            RotorPeriodicForcing::Fourier { rows, orders } => prepare_maxwell_fourier_power_cycle(
                -2.0,
                &segments,
                &branches,
                rows,
                orders,
                1e-12,
                budget(),
                || false,
            ),
            RotorPeriodicForcing::Homogeneous => unreachable!(),
        }
        .unwrap();
        let actual = clocks
            .prepare_power_cycle(-2.0, &segments, &branches, source, 1e-12, budget(), || {
                false
            })
            .unwrap();
        assert_eq!(
            actual.motion_cycle().floquet_report(),
            expected.motion_cycle().floquet_report()
        );
        assert_eq!(
            actual.advance(&initial, 1e-9, 1e-9, || false).unwrap(),
            expected.advance(&initial, 1e-9, 1e-9, || false).unwrap()
        );
    }
}

#[test]
fn incomplete_or_unrepresentable_clocks_and_cancellation_publish_no_partial_cycle() {
    let body = model(1.0, 0.2, 2.0);
    let branches = [arm(0.0, 2.0, 2.0), arm(0.0, 3.0, 0.25)];
    let segments = [0.2, 0.3].map(|duration_seconds| RotorPeriodicSegment {
        model: &body,
        duration_seconds,
    });
    for scales in [
        vec![],
        vec![vec![1.0, 1.0]],
        vec![vec![1.0], vec![1.0, 1.0]],
        vec![vec![0.0, 1.0], vec![1.0, 1.0]],
        vec![vec![-1.0, 1.0], vec![1.0, 1.0]],
        vec![vec![f64::NAN, 1.0], vec![1.0, 1.0]],
        vec![vec![1.0, f64::INFINITY], vec![1.0, 1.0]],
        vec![vec![f64::MAX, 1.0], vec![1.0, 1.0]],
        vec![vec![1.0, 1e-309], vec![1.0, 1.0]], // Positive time, infinite inverse.
        vec![vec![1.0, f64::from_bits(1)], vec![1.0, 1.0]],
    ] {
        let clocks = RotorMaxwellKinetics {
            relaxation_time_scales: &scales,
        };
        assert!(matches!(
            clocks.validate(2, &branches, || false),
            Err(RotorError::Invalid(_))
        ));
        assert!(
            clocks
                .prepare_cycle(
                    0.0,
                    &segments,
                    &branches,
                    RotorPeriodicForcing::Homogeneous,
                    1e-12,
                    budget(),
                    || false
                )
                .is_err()
        );
        assert!(
            clocks
                .prepare_power_cycle(
                    0.0,
                    &segments,
                    &branches,
                    RotorPeriodicForcing::Homogeneous,
                    1e-12,
                    budget(),
                    || false
                )
                .is_err()
        );
    }
    let scales = [vec![0.5, 2.0], vec![3.0, 0.25]];
    let clocks = RotorMaxwellKinetics {
        relaxation_time_scales: &scales,
    };
    assert!(matches!(
        clocks.validate(2, &branches, || true),
        Err(RotorError::Cancelled)
    ));
    let mut polls = 0;
    let prepared = clocks
        .prepare_power_cycle(
            0.0,
            &segments,
            &branches,
            RotorPeriodicForcing::Homogeneous,
            1e-12,
            budget(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            clocks.prepare_power_cycle(
                0.0,
                &segments,
                &branches,
                RotorPeriodicForcing::Homogeneous,
                1e-12,
                budget(),
                || {
                    count += 1;
                    count == stop
                }
            ),
            Err(RotorError::Cancelled)
        ));
    }
    let initial = state(&[0.0, 0.0, 0.2, -0.3]);
    let expected = prepared.advance(&initial, 1e-9, 1e-10, || false).unwrap();
    assert!(matches!(
        prepared.advance(&initial, 1e-9, 1e-10, || true),
        Err(RotorCyclePowerError::Dynamics(RotorError::Cancelled))
    ));
    assert_eq!(
        expected,
        prepared.advance(&initial, 1e-9, 1e-10, || false).unwrap()
    );
    assert!(
        clocks
            .prepare_power_cycle(
                0.0,
                &segments,
                &branches,
                RotorPeriodicForcing::Homogeneous,
                1e-12,
                RotorFloquetBudget {
                    max_dense_bytes: 1,
                    ..budget()
                },
                || false
            )
            .is_err()
    );
    assert_eq!(branches[0].relaxation_time_s, 2.0);
    assert_eq!(branches[1].relaxation_time_s, 0.25);
}
