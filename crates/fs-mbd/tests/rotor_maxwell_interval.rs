//! G1/G3 open physical-frequency motion, original memory and finite-window work.
use fs_mbd::rotordynamics::periodic::{
    RotorFloquetBudget, RotorFourierForce, RotorHarmonicForce, RotorIntervalForcing,
    RotorMaxwellKinetics, RotorMaxwellState, RotorPeriodicSegment, prepare_maxwell_cycle,
    prepare_maxwell_interval, prepare_maxwell_power_interval,
};
use fs_mbd::rotordynamics::relaxation::MaxwellDampingBranch;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn model(c: f64, k: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[1.2], &[c], &[k], &[0.0], || false).unwrap()
}

fn branches() -> [MaxwellDampingBranch; 2] {
    [
        MaxwellDampingBranch {
            strain: vec![1.0],
            stiffness: 8.0,
            relaxation_time_s: 0.1,
        },
        MaxwellDampingBranch {
            strain: vec![-0.6],
            stiffness: 3.5,
            relaxation_time_s: 0.25,
        },
    ]
}

fn state() -> RotorMaxwellState {
    RotorMaxwellState {
        displacement: vec![0.2],
        velocity: vec![-0.1],
        branch_stress: vec![0.35, -0.2],
    }
}

fn budget() -> RotorFloquetBudget {
    RotorFloquetBudget {
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 500_000_000,
        ..RotorFloquetBudget::default()
    }
}

fn same_state(a: &RotorMaxwellState, b: &RotorMaxwellState, tolerance: f64) {
    for (&actual, &expected) in a
        .displacement
        .iter()
        .chain(&a.velocity)
        .chain(&a.branch_stress)
        .zip(
            b.displacement
                .iter()
                .chain(&b.velocity)
                .chain(&b.branch_stress),
        )
    {
        near(actual, expected, tolerance);
    }
}

// Independent physical ODE, without the production exponentials, modal poles,
// maps or quadratic-flow integrals. Slots retain both original stresses and
// the four separate physical work/heat integrals.
fn rhs(t: f64, y: [f64; 8], c: f64, k: f64, tau: [f64; 2], a: [f64; 5]) -> [f64; 8] {
    let force = a[0]
        + a[1] * (2.7 * t).cos()
        + a[2] * (2.7 * t).sin()
        + a[3] * (6.25 * t).cos()
        + a[4] * (6.25 * t).sin();
    [
        y[1],
        (force - c * y[1] - k * y[0] - y[2] + 0.6 * y[3]) / 1.2,
        8.0 * y[1] - y[2] / tau[0],
        -2.1 * y[1] - y[3] / tau[1],
        force * y[1],
        c * y[1] * y[1],
        y[2] * y[2] / (8.0 * tau[0]),
        y[3] * y[3] / (3.5 * tau[1]),
    ]
}

fn integrate(
    start: f64,
    dt: f64,
    mut y: [f64; 8],
    c: f64,
    k: f64,
    tau: [f64; 2],
    a: [f64; 5],
    steps: usize,
) -> [f64; 8] {
    let h = dt / steps as f64;
    for i in 0..steps {
        let t = start + h * i as f64;
        let k1 = rhs(t, y, c, k, tau, a);
        let k2 = rhs(
            t + h / 2.0,
            std::array::from_fn(|j| y[j] + h * k1[j] / 2.0),
            c,
            k,
            tau,
            a,
        );
        let k3 = rhs(
            t + h / 2.0,
            std::array::from_fn(|j| y[j] + h * k2[j] / 2.0),
            c,
            k,
            tau,
            a,
        );
        let k4 = rhs(
            t + h,
            std::array::from_fn(|j| y[j] + h * k3[j]),
            c,
            k,
            tau,
            a,
        );
        y = std::array::from_fn(|j| y[j] + h * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0);
    }
    y
}

#[test]
fn switched_physical_tones_match_original_ode_work_heat_and_open_storage() {
    let c = [0.17, 0.08, 0.31];
    let k = [14.0, 19.0, 11.0];
    let durations = [0.17, 0.23, 0.31];
    let models: Vec<_> = (0..3).map(|i| model(c[i], k[i])).collect();
    let segments: Vec<_> = (0..3)
        .map(|i| RotorPeriodicSegment {
            model: &models[i],
            duration_seconds: durations[i],
        })
        .collect();
    let branches = branches();
    let scales = vec![vec![0.25, 3.0], vec![1.6, 0.4], vec![0.75, 1.1]];
    let amplitudes = [
        [0.3, 0.7, 0.2, 0.1, -0.05],
        [-0.1, 0.4, -0.3, 0.2, 0.15],
        [0.2, -0.5, 0.1, 0.3, 0.2],
    ];
    let tones: Vec<_> = amplitudes
        .iter()
        .map(|a| {
            [
                RotorHarmonicForce {
                    cosine: &a[1..2],
                    sine: &a[2..3],
                },
                RotorHarmonicForce {
                    cosine: &a[3..4],
                    sine: &a[4..5],
                },
            ]
        })
        .collect();
    let rows: Vec<_> = (0..3)
        .map(|i| RotorFourierForce {
            steady: &amplitudes[i][0..1],
            harmonics: &tones[i],
        })
        .collect();
    let initial = state();
    let prepared = prepare_maxwell_power_interval(
        0.0,
        0.43,
        &segments,
        &branches,
        Some(RotorMaxwellKinetics {
            relaxation_time_scales: &scales,
        }),
        RotorIntervalForcing::Fourier {
            rows: &rows,
            angular_frequencies_rad_s: &[2.7, 6.25],
        },
        1e-11,
        budget(),
        || false,
    )
    .unwrap();
    let (trajectory, work) = prepared.advance(&initial, 1e-10, 1e-9, || false).unwrap();
    let mut y = [0.2, -0.1, 0.35, -0.2, 0.0, 0.0, 0.0, 0.0];
    let mut time = 0.43;
    let mut switch_work = 0.0;
    for i in 0..3 {
        let before = y;
        let tau = [0.1 * scales[i][0], 0.25 * scales[i][1]];
        y = integrate(
            time,
            durations[i],
            before,
            c[i],
            k[i],
            tau,
            amplitudes[i],
            4096,
        );
        let coarse = integrate(
            time,
            durations[i],
            before,
            c[i],
            k[i],
            tau,
            amplitudes[i],
            2048,
        );
        for j in 0..8 {
            near(y[j], coarse[j], 1e-11);
        }
        time += durations[i];
        near(trajectory.boundary_times_seconds[i + 1], time, 2e-15);
        same_state(
            &trajectory.boundary_states[i + 1],
            &RotorMaxwellState {
                displacement: vec![y[0]],
                velocity: vec![y[1]],
                branch_stress: vec![y[2], y[3]],
            },
            2e-9,
        );
        near(work.segments[i].applied_work_j, y[4] - before[4], 2e-9);
        near(
            work.segments[i].velocity_force_work_j,
            y[5] - before[5],
            2e-9,
        );
        for j in 0..2 {
            near(
                work.segments[i].branch_heat_j[j],
                y[6 + j] - before[6 + j],
                2e-9,
            );
            near(
                work.segments[i].branch_elastic_energy_change_j[j],
                (y[2 + j].powi(2) - before[2 + j].powi(2)) / (2.0 * branches[j].stiffness),
                2e-9,
            );
        }
        if i < 2 {
            let jump = (k[i + 1] - k[i]) * y[0].powi(2) / 2.0;
            switch_work += jump;
            near(work.segments[i].stiffness_switch_work_j, jump, 2e-9);
        } else {
            assert_eq!(work.segments[i].stiffness_switch_work_j, 0.0);
        }
    }
    near(work.applied_work_j, y[4], 2e-9);
    near(work.velocity_force_work_j, y[5], 2e-9);
    near(work.stiffness_switch_work_j, switch_work, 2e-9);
    near(
        work.relaxed_elastic_energy_change_j,
        (k[2] * y[0].powi(2) - k[0] * 0.2_f64.powi(2)) / 2.0,
        2e-9,
    );
    near(
        work.kinetic_energy_change_j,
        0.6 * (y[1].powi(2) - 0.01),
        2e-9,
    );
    for j in 0..2 {
        near(work.branch_heat_j[j], y[6 + j], 2e-9);
    }
    near(work.balance_residual_j, 0.0, 1e-9);
    assert_eq!(trajectory.boundary_times_seconds[0], 0.43);
    assert_eq!(initial, state());
}

#[test]
fn unequal_windows_keep_global_phase_and_agree_with_harmonic_fourier_ingress() {
    let model = model(0.17, 14.0);
    let branches = branches();
    let initial = state();
    let durations = [0.13, 0.29, 0.17];
    let total = durations.iter().sum();
    let start = -0.37;
    let single = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: total,
    }];
    let tones = [RotorHarmonicForce {
        cosine: &[0.7],
        sine: &[-0.3],
    }];
    let direct = prepare_maxwell_power_interval(
        0.0,
        start,
        &single,
        &branches,
        None,
        RotorIntervalForcing::Harmonic {
            rows: &tones,
            angular_frequency_rad_s: 2.7,
        },
        1e-11,
        budget(),
        || false,
    )
    .unwrap();
    let (whole, whole_power) = direct.advance(&initial, 1e-10, 1e-9, || false).unwrap();
    let mut current = initial.clone();
    let mut time = start;
    let mut heat = [0.0; 2];
    let mut applied = 0.0;
    for duration in durations {
        let segment = [RotorPeriodicSegment {
            model: &model,
            duration_seconds: duration,
        }];
        let rows = [RotorFourierForce {
            steady: &[0.0],
            harmonics: &tones,
        }];
        let next = prepare_maxwell_power_interval(
            0.0,
            time,
            &segment,
            &branches,
            None,
            RotorIntervalForcing::Fourier {
                rows: &rows,
                angular_frequencies_rad_s: &[2.7],
            },
            1e-11,
            budget(),
            || false,
        )
        .unwrap();
        let (trajectory, power) = next.advance(&current, 1e-10, 1e-9, || false).unwrap();
        assert_eq!(trajectory.boundary_times_seconds[0], time);
        time += duration;
        near(
            *trajectory.boundary_times_seconds.last().unwrap(),
            time,
            2e-15,
        );
        current = trajectory.boundary_states.last().unwrap().clone();
        applied += power.applied_work_j;
        for j in 0..2 {
            heat[j] += power.branch_heat_j[j];
        }
    }
    same_state(&current, whole.boundary_states.last().unwrap(), 2e-10);
    near(applied, whole_power.applied_work_j, 2e-10);
    for j in 0..2 {
        near(heat[j], whole_power.branch_heat_j[j], 2e-10);
    }
    // A shifted start is a physically different drive, not a new phase origin.
    let shifted = prepare_maxwell_interval(
        0.0,
        start + 0.2,
        &single,
        &branches,
        None,
        RotorIntervalForcing::Harmonic {
            rows: &tones,
            angular_frequency_rad_s: 2.7,
        },
        1e-11,
        budget(),
        || false,
    )
    .unwrap()
    .advance(&initial, 1e-10, || false)
    .unwrap();
    assert!((shifted.boundary_states[1].displacement[0] - current.displacement[0]).abs() > 1e-4);

    // The original binary64 omega=1.3 multiplied by 1e9 has a nonzero low
    // product. Rounding it before phase reduction loses 4.44e-8 radians.
    // These free-coordinate endpoints were computed at 90 decimal digits
    // from the exact binary64 omega and integer start, then rounded once:
    // v=(sin(w*(t+h))-sin(w*t))/w,
    // q=(cos(w*t)-cos(w*(t+h)))/w²-h*sin(w*t)/w.
    let free = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let free_segment = [RotorPeriodicSegment {
        model: &free,
        duration_seconds: 0.5,
    }];
    let free_rows = [RotorHarmonicForce {
        cosine: &[1.0],
        sine: &[0.0],
    }];
    let free_motion = prepare_maxwell_interval(
        0.0,
        1e9,
        &free_segment,
        &[],
        None,
        RotorIntervalForcing::Harmonic {
            rows: &free_rows,
            angular_frequency_rad_s: 1.3,
        },
        1e-11,
        budget(),
        || false,
    )
    .unwrap()
    .advance(
        &RotorMaxwellState {
            displacement: vec![0.0],
            velocity: vec![0.0],
            branch_stress: vec![],
        },
        1e-10,
        || false,
    )
    .unwrap();
    near(
        free_motion.boundary_states[1].displacement[0],
        0.116_524_462_757_239_5,
        5e-13,
    );
    near(
        free_motion.boundary_states[1].velocity[0],
        0.442_915_743_990_642_45,
        5e-13,
    );
}

#[test]
fn open_window_needs_no_floquet_spectrum_even_after_complete_dark_decay() {
    let model = model(0.4, 2.0);
    let branches = [MaxwellDampingBranch {
        strain: vec![0.0],
        stiffness: 2.0,
        relaxation_time_s: 0.001,
    }];
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 1.0,
    }];
    let initial = RotorMaxwellState {
        displacement: vec![0.0],
        velocity: vec![0.0],
        branch_stress: vec![4.0],
    };
    let mut no_qr = budget();
    no_qr.max_qr_sweeps_per_window = 0;
    let plan = prepare_maxwell_power_interval(
        0.0,
        3.0,
        &segments,
        &branches,
        None,
        RotorIntervalForcing::Homogeneous,
        1e-11,
        no_qr,
        || false,
    )
    .unwrap();
    let (trajectory, power) = plan.advance(&initial, 1e-10, 1e-9, || false).unwrap();
    assert_eq!(trajectory.boundary_times_seconds, vec![3.0, 4.0]);
    near(trajectory.boundary_states[1].branch_stress[0], 0.0, 1e-14);
    near(power.branch_heat_j[0], 4.0, 1e-10);
    near(power.branch_elastic_energy_change_j[0], -4.0, 1e-10);
    assert_eq!(power.applied_work_j, 0.0);
    assert_eq!(power.stiffness_switch_work_j, 0.0);
    assert!(prepare_maxwell_cycle(0.0, &segments, &branches, 1e-11, no_qr, || false).is_err());
}

#[test]
fn physical_clock_ranges_budgets_cancellation_and_replay_are_enforced() {
    let model = model(0.17, 14.0);
    let branches = branches();
    let segments = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.2,
    }];
    let harmonic = RotorHarmonicForce {
        cosine: &[0.7],
        sine: &[0.3],
    };
    let rows = [harmonic];
    for start in [f64::NAN, f64::INFINITY, f64::MAX, 1e16] {
        assert!(
            prepare_maxwell_interval(
                0.0,
                start,
                &segments,
                &branches,
                None,
                RotorIntervalForcing::Homogeneous,
                1e-11,
                budget(),
                || false
            )
            .is_err()
        );
    }
    // The endpoint is distinct, but its rounded increment is .0009765625 s,
    // not the declared .001 s. Refuse an underresolved physical time boundary.
    let short = [RotorPeriodicSegment {
        model: &model,
        duration_seconds: 0.001,
    }];
    assert!(
        prepare_maxwell_interval(
            0.0,
            1e12,
            &short,
            &branches,
            None,
            RotorIntervalForcing::Homogeneous,
            1e-11,
            budget(),
            || false,
        )
        .is_err()
    );
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e9] {
        assert!(
            prepare_maxwell_interval(
                0.0,
                0.0,
                &segments,
                &branches,
                None,
                RotorIntervalForcing::Harmonic {
                    rows: &rows,
                    angular_frequency_rad_s: omega
                },
                1e-11,
                budget(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        prepare_maxwell_interval(
            0.0,
            1e13,
            &segments,
            &branches,
            None,
            RotorIntervalForcing::Harmonic {
                rows: &rows,
                angular_frequency_rad_s: 2.7
            },
            1e-11,
            budget(),
            || false
        )
        .is_err()
    );
    for frequencies in [
        vec![2.7, 2.7],
        vec![6.25, 2.7],
        (1..=33).map(f64::from).collect(),
    ] {
        let harmonics = vec![harmonic; frequencies.len()];
        let fourier = [RotorFourierForce {
            steady: &[0.0],
            harmonics: &harmonics,
        }];
        assert!(
            prepare_maxwell_interval(
                0.0,
                0.0,
                &segments,
                &branches,
                None,
                RotorIntervalForcing::Fourier {
                    rows: &fourier,
                    angular_frequencies_rad_s: &frequencies
                },
                1e-11,
                budget(),
                || false
            )
            .is_err()
        );
    }
    let forcing = RotorIntervalForcing::Harmonic {
        rows: &rows,
        angular_frequency_rad_s: 2.7,
    };
    for limited in [
        RotorFloquetBudget {
            max_segments: 0,
            ..budget()
        },
        RotorFloquetBudget {
            max_dense_bytes: 0,
            ..budget()
        },
        RotorFloquetBudget {
            max_multiply_terms: 0,
            ..budget()
        },
    ] {
        assert!(
            prepare_maxwell_interval(
                0.0,
                0.2,
                &segments,
                &branches,
                None,
                forcing,
                1e-11,
                limited,
                || false
            )
            .is_err()
        );
    }
    assert!(matches!(
        prepare_maxwell_interval(
            0.0,
            0.2,
            &segments,
            &branches,
            None,
            forcing,
            1e-11,
            budget(),
            || true
        ),
        Err(RotorError::Cancelled)
    ));
    let plan = prepare_maxwell_interval(
        0.0,
        0.2,
        &segments,
        &branches,
        None,
        forcing,
        1e-11,
        budget(),
        || false,
    )
    .unwrap();
    let initial = state();
    let before = initial.clone();
    let first = plan.advance(&initial, 1e-10, || false).unwrap();
    assert!(matches!(
        plan.advance(&initial, 1e-10, || true),
        Err(RotorError::Cancelled)
    ));
    assert_eq!(first, plan.advance(&initial, 1e-10, || false).unwrap());
    assert_eq!(initial, before);
}
