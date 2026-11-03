//! Native affine rotor cycles and periodic forced response, with analytic physics oracles.
use fs_mbd::rotordynamics::periodic::{
    RotorFloquetBudget, RotorPeriodicSegment, periodic_floquet, prepare_forced_cycle,
};
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::cell::Cell;
use std::f64::consts::PI;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn oscillator(mass: f64, damping: f64, stiffness: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[mass], &[damping], &[stiffness], &[0.0], || false).unwrap()
}

fn scalar_state(q: f64, v: f64) -> RotorTransientState {
    RotorTransientState {
        displacement: vec![q],
        velocity: vec![v],
    }
}

// Independent analytic underdamped transition, in original displacement/velocity.
fn damped_transition(mass: f64, damping: f64, stiffness: f64, time: f64) -> [f64; 4] {
    let alpha = damping / (2.0 * mass);
    let frequency = (stiffness / mass - alpha * alpha).sqrt();
    let (sin, cos) = (frequency * time).sin_cos();
    let decay = (-alpha * time).exp();
    [
        decay * (cos + alpha * sin / frequency),
        decay * sin / frequency,
        -stiffness / mass * decay * sin / frequency,
        decay * (cos - alpha * sin / frequency),
    ]
}

#[test]
fn damped_square_wave_resonance_matches_analytic_orbit_across_force_scales() {
    let model = oscillator(1.0, 0.02, 1.0);
    let cycle = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: PI,
    }; 2];
    let e = damped_transition(1.0, 0.02, 1.0, PI);
    // Half-period positive forcing has offset g=(I-E)[F/K,0].
    // Antisymmetry gives x_half=-x_0, hence (I+E)x_0=-g.
    let (gq, gv) = (1.0 - e[0], -e[2]);
    let determinant = (1.0 + e[0]) * (1.0 + e[3]) - e[1] * e[2];
    let expected = [
        (-(1.0 + e[3]) * gq + e[1] * gv) / determinant,
        (e[2] * gq - (1.0 + e[0]) * gv) / determinant,
    ];
    assert!(expected[0].abs() > 50.0);
    for amplitude in [1e-80, 1.0, 1e80] {
        let positive = [amplitude];
        let negative = [-amplitude];
        let prepared = prepare_forced_cycle(
            0.0,
            &cycle,
            &[&positive, &negative],
            1e-12,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        let orbit = prepared.periodic_orbit(1e-10, 1e8, || false).unwrap();
        assert_eq!(orbit.trajectory.boundary_states.len(), 3);
        for (state, sign) in orbit
            .trajectory
            .boundary_states
            .iter()
            .zip([1.0, -1.0, 1.0])
        {
            // Divide out the force: an absolute floor must not hide small loads.
            near(state.displacement[0] / amplitude, sign * expected[0], 2e-11);
            near(state.velocity[0] / amplitude, sign * expected[1], 2e-11);
        }
        assert!(orbit.relative_fixed_point_residual < 1e-10);
        assert!(orbit.trajectory.relative_affine_residual < 1e-10);
        assert!(orbit.relative_closure_residual < 1e-10);
        assert!(prepared.floquet_report().spectral_radius < 1.0);
    }
    // Splitting a true held load leaves the physical waveform unchanged.
    let refined = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: PI / 2.0,
    }; 4];
    let prepared = prepare_forced_cycle(
        0.0,
        &refined,
        &[&[1.0], &[1.0], &[-1.0], &[-1.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let orbit = prepared.periodic_orbit(1e-10, 1e8, || false).unwrap();
    near(
        orbit.trajectory.boundary_states[0].displacement[0],
        expected[0],
        2e-11,
    );
    near(
        orbit.trajectory.boundary_states[0].velocity[0],
        expected[1],
        2e-11,
    );
}

#[test]
fn nonbinary_static_equilibrium_does_not_refuse_a_zero_physical_velocity() {
    let model = oscillator(3.1, 0.7, 7.3);
    let prepared = prepare_forced_cycle(
        -2.0,
        &[RotorPeriodicSegment {
            model: &model,
            duration_seconds: 0.7,
        }],
        &[&[2.7]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let orbit = prepared.periodic_orbit(1e-11, 1e8, || false).unwrap();
    for state in &orbit.trajectory.boundary_states {
        near(state.displacement[0], 2.7 / 7.3, 2e-12);
        near(state.velocity[0], 0.0, 2e-12);
    }
    assert!(orbit.relative_closure_residual < 1e-11);
    let direct = prepared
        .advance(&scalar_state(2.7 / 7.3, 0.0), 1e-11, || false)
        .unwrap();
    near(direct.boundary_states[1].displacement[0], 2.7 / 7.3, 2e-12);
    near(direct.boundary_states[1].velocity[0], 0.0, 2e-12);
}

#[test]
fn unstable_stiffness_cycle_has_a_periodic_particular_orbit_but_growing_perturbations() {
    let soft = oscillator(1.0, 0.0, 1.0);
    let stiff = oscillator(1.0, 0.0, 4.0);
    let prepared = prepare_forced_cycle(
        7.0,
        &[
            RotorPeriodicSegment {
                model: &soft,
                duration_seconds: PI / 2.0,
            },
            RotorPeriodicSegment {
                model: &stiff,
                duration_seconds: PI / 4.0,
            },
        ],
        &[&[1.0], &[1.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    // Quarter-period maps give P=diag(-1/2,-2), b=[3/4,-3/2].
    near(prepared.forcing_offset()[0], 0.75, 2e-12);
    near(prepared.forcing_offset()[1], -1.5, 2e-12);
    near(prepared.floquet_report().spectral_radius, 2.0, 2e-12);
    let orbit = prepared.periodic_orbit(1e-11, 1e8, || false).unwrap();
    for (state, velocity) in orbit
        .trajectory
        .boundary_states
        .iter()
        .zip([-0.5, 0.5, -0.5])
    {
        near(state.displacement[0], 0.5, 2e-12);
        near(state.velocity[0], velocity, 2e-12);
    }
    let perturbation = 0.125;
    let first = prepared
        .advance(&scalar_state(0.5, -0.5 + perturbation), 1e-11, || false)
        .unwrap();
    let first_end = first.boundary_states.last().unwrap();
    near(first_end.velocity[0], -0.5 - 2.0 * perturbation, 2e-12);
    let second = prepared.advance(first_end, 1e-11, || false).unwrap();
    near(
        second.boundary_states.last().unwrap().velocity[0],
        -0.5 + 4.0 * perturbation,
        2e-12,
    );
}

#[test]
fn held_forces_use_full_mutual_inertia_and_preserve_physical_boundary_continuity() {
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.5, 0.5, 1.0],
        &[0.0; 4],
        &[0.0; 4],
        &[0.0; 4],
        || false,
    )
    .unwrap();
    let prepared = prepare_forced_cycle(
        9.0,
        &[
            RotorPeriodicSegment {
                model: &model,
                duration_seconds: 0.25,
            },
            RotorPeriodicSegment {
                model: &model,
                duration_seconds: 0.5,
            },
        ],
        &[&[1.0, -0.5], &[-0.25, 0.75]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let initial = RotorTransientState {
        displacement: vec![0.2, -0.3],
        velocity: vec![-0.4, 0.5],
    };
    let trace = prepared.advance(&initial, 1e-11, || false).unwrap();
    assert_eq!(trace.boundary_times_seconds, vec![0.0, 0.25, 0.75]);
    assert_eq!(trace.boundary_states[0], initial);
    let mut expected = initial;
    for (index, (duration, acceleration)) in [
        (0.25, [5.0 / 7.0, -6.0 / 7.0]),
        (0.5, [-5.0 / 14.0, 13.0 / 14.0]),
    ]
    .into_iter()
    .enumerate()
    {
        for (i, a) in acceleration.into_iter().enumerate() {
            expected.displacement[i] +=
                duration * expected.velocity[i] + 0.5 * duration * duration * a;
            expected.velocity[i] += duration * a;
            near(
                trace.boundary_states[index + 1].displacement[i],
                expected.displacement[i],
                2e-12,
            );
            near(
                trace.boundary_states[index + 1].velocity[i],
                expected.velocity[i],
                2e-12,
            );
        }
    }
    assert!(matches!(
        prepared.periodic_orbit(1e-11, 1e8, || false),
        Err(RotorError::PeriodicLinear(_))
    ));
}

fn congruence(matrix: &[f64], h: &[f64; 4]) -> Vec<f64> {
    let mut result = vec![0.0; 4];
    for i in 0..2 {
        for j in 0..2 {
            for p in 0..2 {
                for q in 0..2 {
                    result[2 * i + j] += h[2 * p + i] * matrix[2 * p + q] * h[2 * q + j];
                }
            }
        }
    }
    result
}

fn physical_pair(h: &[f64; 4], values: &[f64]) -> [f64; 2] {
    [
        h[0] * values[0] + h[1] * values[1],
        h[2] * values[0] + h[3] * values[1],
    ]
}

#[test]
fn forced_orbit_respects_dense_coordinate_congruence_and_cycle_origin() {
    let mass = [1.0, 0.0, 0.0, 1.0];
    let gyro = [0.0, 0.5, -0.5, 0.0];
    let damping = [[0.125, 0.25, -0.125, -0.0625], [0.25, 0.125, 0.0625, 0.125]];
    let stiffness = [[4.0, 0.5, -0.25, 2.0], [3.0, -0.25, 0.75, 5.0]];
    let forces = [[0.5, -0.25], [0.125, 0.75]];
    let h = [2.0, 0.25, -0.5, 1.5];
    let original: Vec<_> = (0..2)
        .map(|i| {
            RotorBearingModel::new(2, &mass, &damping[i], &stiffness[i], &gyro, || false).unwrap()
        })
        .collect();
    let transformed: Vec<_> = (0..2)
        .map(|i| {
            RotorBearingModel::new(
                2,
                &congruence(&mass, &h),
                &congruence(&damping[i], &h),
                &congruence(&stiffness[i], &h),
                &congruence(&gyro, &h),
                || false,
            )
            .unwrap()
        })
        .collect();
    let transformed_forces: Vec<_> = forces
        .iter()
        .map(|force| {
            [
                h[0] * force[0] + h[2] * force[1],
                h[1] * force[0] + h[3] * force[1],
            ]
        })
        .collect();
    let analyze = |models: &[RotorBearingModel], loads: &[[f64; 2]]| {
        let prepared = prepare_forced_cycle(
            -1.25,
            &[
                RotorPeriodicSegment {
                    model: &models[0],
                    duration_seconds: 0.3,
                },
                RotorPeriodicSegment {
                    model: &models[1],
                    duration_seconds: 0.2,
                },
            ],
            &[&loads[0], &loads[1]],
            1e-12,
            RotorFloquetBudget::default(),
            || false,
        )
        .unwrap();
        prepared.periodic_orbit(1e-10, 1e8, || false).unwrap()
    };
    let a = analyze(&original, &forces);
    let b = analyze(&transformed, &transformed_forces);
    for (a, b) in a
        .trajectory
        .boundary_states
        .iter()
        .zip(&b.trajectory.boundary_states)
    {
        for (actual, expected) in a
            .displacement
            .iter()
            .zip(physical_pair(&h, &b.displacement))
        {
            near(*actual, expected, 2e-10);
        }
        for (actual, expected) in a.velocity.iter().zip(physical_pair(&h, &b.velocity)) {
            near(*actual, expected, 2e-10);
        }
    }
    let shifted = prepare_forced_cycle(
        -1.25,
        &[
            RotorPeriodicSegment {
                model: &original[1],
                duration_seconds: 0.2,
            },
            RotorPeriodicSegment {
                model: &original[0],
                duration_seconds: 0.3,
            },
        ],
        &[&forces[1], &forces[0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap()
    .periodic_orbit(1e-10, 1e8, || false)
    .unwrap();
    for (actual, expected) in shifted.trajectory.boundary_states[0]
        .displacement
        .iter()
        .zip(&a.trajectory.boundary_states[1].displacement)
    {
        near(*actual, *expected, 2e-10);
    }
    for (actual, expected) in shifted.trajectory.boundary_states[0]
        .velocity
        .iter()
        .zip(&a.trajectory.boundary_states[1].velocity)
    {
        near(*actual, *expected, 2e-10);
    }
}

#[test]
fn singular_periodic_families_and_unresolved_lossless_resonance_refuse_unique_orbits() {
    let free = oscillator(1.0, 0.0, 0.0);
    let family = prepare_forced_cycle(
        0.0,
        &[RotorPeriodicSegment {
            model: &free,
            duration_seconds: 1.0,
        }],
        &[&[0.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        family.periodic_orbit(1e-11, 1e8, || false),
        Err(RotorError::PeriodicLinear(_))
    ));
    // Perfectly compatible zero forcing still leaves an arbitrary displacement.
    let trace = family
        .advance(&scalar_state(0.25, 0.5), 1e-11, || false)
        .unwrap();
    near(trace.boundary_states[1].displacement[0], 0.75, 1e-14);
    near(trace.boundary_states[1].velocity[0], 0.5, 1e-14);

    let lossless = oscillator(1.0, 0.0, 1.0);
    let resonant = prepare_forced_cycle(
        0.0,
        &[RotorPeriodicSegment {
            model: &lossless,
            duration_seconds: PI,
        }; 2],
        &[&[1.0], &[-1.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        resonant.periodic_orbit(1e-11, 1e8, || false),
        Err(RotorError::PeriodicLinear(_)) | Err(RotorError::PeriodicOrbitUnresolved { .. })
    ));
}

#[test]
fn a_small_force_is_not_hidden_by_a_large_unchanged_velocity() {
    let free = oscillator(1.0, 0.0, 0.0);
    let prepared = prepare_forced_cycle(
        0.0,
        &[RotorPeriodicSegment {
            model: &free,
            duration_seconds: 1.0,
        }],
        &[&[1e-100]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    let resolved = prepared
        .advance(&scalar_state(0.0, 0.0), 1e-11, || false)
        .unwrap();
    near(
        resolved.boundary_states[1].displacement[0] / 1e-100,
        0.5,
        1e-14,
    );
    near(resolved.boundary_states[1].velocity[0] / 1e-100, 1.0, 1e-14);
    let initial = scalar_state(0.0, 1e100);
    let original = initial.clone();
    assert!(matches!(
        prepared.advance(&initial, 1e-11, || false),
        Err(RotorError::Residual { .. })
    ));
    assert_eq!(initial, original);
}

#[test]
fn zero_forcing_reuses_the_homogeneous_owner_and_zero_state_closes() {
    let model = oscillator(3.1, 0.7, 7.3);
    let cycle = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.35,
    }; 2];
    let homogeneous =
        periodic_floquet(2.0, &cycle, 1e-12, RotorFloquetBudget::default(), || false).unwrap();
    let prepared = prepare_forced_cycle(
        2.0,
        &cycle,
        &[&[0.0], &[0.0]],
        1e-12,
        RotorFloquetBudget::default(),
        || false,
    )
    .unwrap();
    assert_eq!(prepared.floquet_report(), &homogeneous);
    assert_eq!(prepared.forcing_offset(), &[0.0, 0.0]);
    let orbit = prepared.periodic_orbit(1e-11, 1e8, || false).unwrap();
    assert!(
        orbit
            .trajectory
            .boundary_states
            .iter()
            .all(|state| *state == scalar_state(0.0, 0.0))
    );
    assert_eq!(orbit.relative_fixed_point_residual, 0.0);
    assert_eq!(orbit.relative_closure_residual, 0.0);
}

#[test]
fn force_shapes_budgets_and_cancellation_refuse_without_partial_publication() {
    let model = oscillator(1.0, 0.2, 1.0);
    let cycle = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.3,
    }; 2];
    let defaults = RotorFloquetBudget::default();
    for loads in [
        vec![],
        vec![&[1.0][..]],
        vec![&[1.0, 2.0][..], &[0.0][..]],
        vec![&[f64::NAN][..], &[0.0][..]],
    ] {
        assert!(prepare_forced_cycle(0.0, &cycle, &loads, 1e-12, defaults, || false).is_err());
    }
    let small = RotorFloquetBudget {
        max_dense_bytes: 1200,
        ..defaults
    };
    assert!(periodic_floquet(0.0, &cycle, 1e-12, small, || false).is_ok());
    assert!(matches!(
        prepare_forced_cycle(0.0, &cycle, &[&[1.0], &[-0.5]], 1e-12, small, || false),
        Err(RotorError::PeriodicBudget(_))
    ));
    let polls = Cell::new(0);
    let prepared = prepare_forced_cycle(0.0, &cycle, &[&[1.0], &[-0.5]], 1e-12, defaults, || {
        polls.set(polls.get() + 1);
        false
    })
    .unwrap();
    let count = polls.get();
    for stop in [1, count / 2, count] {
        let calls = Cell::new(0);
        assert!(matches!(
            prepare_forced_cycle(0.0, &cycle, &[&[1.0], &[-0.5]], 1e-12, defaults, || {
                calls.set(calls.get() + 1);
                calls.get() >= stop
            }),
            Err(RotorError::Cancelled)
        ));
    }
    let initial = scalar_state(0.2, -0.3);
    let original = initial.clone();
    assert!(matches!(
        prepared.advance(&initial, 1e-11, || true),
        Err(RotorError::Cancelled)
    ));
    assert!(matches!(
        prepared.periodic_orbit(1e-11, 1e8, || true),
        Err(RotorError::Cancelled)
    ));
    assert_eq!(initial, original);
    for limit in [0.0, f64::NAN, f64::INFINITY] {
        assert!(prepared.periodic_orbit(1e-11, limit, || false).is_err());
    }
}
