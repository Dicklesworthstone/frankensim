use fs_euler_disc_e2e::coupled_runner::{
    ADAPTER_INTEGRATION_PLAN, CoupledChannelFactors, CoupledControls, CoupledError, CoupledFactors,
    CoupledInitialState, initial_contact_point_velocity, initial_qois, run_closed_profile_reduced,
    run_closed_reduced,
};
use fs_euler_disc_e2e::specimen::DiscProfileSpec;
use fs_exec::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_frep::SquatDiscEdgeTreatment;

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4555_4c45_525f_434f,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn factors(radius_m: f64, density: f64) -> CoupledFactors {
    let mass = std::f64::consts::PI * radius_m * radius_m * 0.006 * density;
    CoupledFactors {
        mass_kg: mass,
        radius_m,
        thickness_m: 0.006,
        transverse_inertia_kg_m2: mass * (3.0 * radius_m * radius_m + 0.006f64.powi(2)) / 12.0,
        axial_inertia_kg_m2: 0.5 * mass * radius_m * radius_m,
        gravity_m_per_s2: 9.806_65,
        sliding_friction_coefficient: 0.42,
        rolling_resistance_m: 4e-5,
        contact_stiffness_n_per_m: 8e4,
        contact_damping_n_s_per_m: 3.0,
        base_effective_mass_kg: 0.25,
        base_stiffness_n_per_m: 4e4,
        base_damping_n_s_per_m: 4.0,
        gas_rotational_damping_n_m_s: 2e-7,
        gas_translation_damping_n_s_per_m: 4e-4,
    }
}
fn controls() -> CoupledControls {
    CoupledControls {
        timestep_s: 2e-5,
        maximum_steps: 200,
        terminal_inclination_rad: 0.002,
        reimpact_limit: 32,
    }
}
fn initial() -> CoupledInitialState {
    CoupledInitialState {
        inclination_rad: 0.08,
        precession_rad_per_s: 16.0,
        spin_rad_per_s: 120.0,
    }
}

#[test]
fn evolves_and_replays_from_a_deterministic_checkpoint() {
    let first =
        run_closed_reduced(factors(0.038, 2680.0), controls(), initial(), None).expect("run");
    let second =
        run_closed_reduced(factors(0.038, 2680.0), controls(), initial(), None).expect("replay");
    assert_eq!(first, second);
    assert!(!first.samples.is_empty());
    assert!(
        first
            .samples
            .windows(2)
            .any(|pair| pair[0].spin_rad_per_s != pair[1].spin_rad_per_s
                || pair[0].inclination_rad != pair[1].inclination_rad)
    );
    assert!(first.samples[0].contact_active);
    assert_eq!(first.samples[0].reimpact_count, 0);
    assert!(
        first
            .samples
            .iter()
            .all(|s| s.energy_defect_j.is_finite() && s.mechanical_energy_j.is_finite())
    );
    assert!(
        first
            .samples
            .iter()
            .all(|sample| sample.channels.rolling.work_j <= 1.0e-14)
    );
    let split = run_closed_reduced(
        factors(0.038, 2680.0),
        CoupledControls {
            maximum_steps: 50,
            ..controls()
        },
        initial(),
        None,
    )
    .expect("prefix");
    let resumed = run_closed_reduced(
        factors(0.038, 2680.0),
        CoupledControls {
            maximum_steps: 150,
            ..controls()
        },
        initial(),
        Some(split.checkpoint),
    )
    .expect("resume");
    assert_eq!(first.checkpoint.state, resumed.checkpoint.state);
    assert_eq!(
        first.checkpoint.accumulated_channel_work_j,
        resumed.checkpoint.accumulated_channel_work_j
    );
}

#[test]
fn restart_rejects_changed_identity_and_documents_adapter_integration() {
    let prefix = run_closed_reduced(
        factors(0.038, 2680.0),
        CoupledControls {
            maximum_steps: 10,
            ..controls()
        },
        initial(),
        None,
    )
    .expect("prefix");
    let assert_mismatch = |factors, controls, checkpoint| {
        assert_eq!(
            run_closed_reduced(factors, controls, initial(), Some(checkpoint)),
            Err(CoupledError::CheckpointMismatch)
        );
    };
    assert_mismatch(
        factors(0.038, 2680.0),
        CoupledControls {
            timestep_s: 1e-5,
            maximum_steps: 10,
            ..controls()
        },
        prefix.checkpoint.clone(),
    );
    assert_mismatch(
        factors(0.038, 2680.0),
        CoupledControls {
            terminal_inclination_rad: 0.003,
            maximum_steps: 10,
            ..controls()
        },
        prefix.checkpoint.clone(),
    );
    assert_mismatch(
        factors(0.038, 2680.0),
        CoupledControls {
            reimpact_limit: 31,
            maximum_steps: 10,
            ..controls()
        },
        prefix.checkpoint.clone(),
    );
    let mut changed_energy = prefix.checkpoint.clone();
    changed_energy.initial_total_energy_j += 1.0;
    assert_mismatch(
        factors(0.038, 2680.0),
        CoupledControls {
            maximum_steps: 10,
            ..controls()
        },
        changed_energy,
    );
    let mut changed_factors = factors(0.038, 2680.0);
    changed_factors.gas_rotational_damping_n_m_s *= 2.0;
    assert_mismatch(
        changed_factors,
        CoupledControls {
            maximum_steps: 10,
            ..controls()
        },
        prefix.checkpoint.clone(),
    );
    let spin_only = CoupledInitialState {
        precession_rad_per_s: 0.0,
        ..initial()
    };
    let spin_prefix = run_closed_reduced(
        factors(0.038, 2680.0),
        CoupledControls {
            maximum_steps: 10,
            ..controls()
        },
        spin_only,
        None,
    )
    .expect("spin prefix");
    let opposite_spin = CoupledInitialState {
        spin_rad_per_s: -spin_only.spin_rad_per_s,
        ..spin_only
    };
    let opposite_spin_fresh = run_closed_reduced(
        factors(0.038, 2680.0),
        CoupledControls {
            maximum_steps: 10,
            ..controls()
        },
        opposite_spin,
        None,
    )
    .expect("opposite-spin run");
    assert_eq!(
        spin_prefix.checkpoint.initial_total_energy_j.to_bits(),
        opposite_spin_fresh
            .checkpoint
            .initial_total_energy_j
            .to_bits()
    );
    assert_eq!(
        run_closed_reduced(
            factors(0.038, 2680.0),
            CoupledControls {
                maximum_steps: 10,
                ..controls()
            },
            opposite_spin,
            Some(spin_prefix.checkpoint),
        ),
        Err(CoupledError::CheckpointMismatch)
    );
    assert!(ADAPTER_INTEGRATION_PLAN.contains("mechanics=>contact/base"));
    assert!(ADAPTER_INTEGRATION_PLAN.contains("rolling_contact=>rolling"));
    assert!(ADAPTER_INTEGRATION_PLAN.contains("air=>gas"));
}

#[test]
fn radius_factor_changes_the_computed_trajectory_without_a_rank_assertion() {
    let small =
        run_closed_reduced(factors(0.030, 2680.0), controls(), initial(), None).expect("small");
    let large =
        run_closed_reduced(factors(0.052, 2680.0), controls(), initial(), None).expect("large");
    assert_ne!(
        small.samples.last().unwrap().spin_rad_per_s,
        large.samples.last().unwrap().spin_rad_per_s
    );
}

#[test]
fn initializer_qoi_round_trip_and_gravity_only_energy_closure() {
    let declared = initial();
    let initial_qoi = initial_qois(factors(0.038, 2680.0), declared).expect("initial qoi");
    assert!((initial_qoi.0 - declared.inclination_rad).abs() < 1.0e-12);
    assert!((initial_qoi.1 - declared.precession_rad_per_s).abs() < 1.0e-12);
    assert!((initial_qoi.2 - declared.spin_rad_per_s).abs() < 1.0e-12);
    assert!(
        initial_contact_point_velocity(factors(0.038, 2680.0), declared)
            .expect("initial contact velocity")
            .norm_squared()
            .sqrt()
            < 1.0e-12
    );

    let mut conservative = factors(0.038, 2680.0);
    conservative.sliding_friction_coefficient = 0.0;
    conservative.rolling_resistance_m = 0.0;
    conservative.contact_stiffness_n_per_m = 0.0;
    conservative.contact_damping_n_s_per_m = 0.0;
    conservative.base_stiffness_n_per_m = 0.0;
    conservative.base_damping_n_s_per_m = 0.0;
    conservative.gas_rotational_damping_n_m_s = 0.0;
    conservative.gas_translation_damping_n_s_per_m = 0.0;
    let run = run_closed_reduced(
        conservative,
        CoupledControls {
            maximum_steps: 1,
            ..controls()
        },
        declared,
        None,
    )
    .expect("gravity-only run");
    assert!(run.checkpoint.accumulated_energy_defect_j.abs() < 1.0e-10);
}

#[test]
fn true_profile_run_uses_profile_support_and_binds_restart_to_chart_identity() {
    with_cx(|cx| {
        let cylinder = DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        }
        .resolve(2_680.0, cx)
        .expect("cylinder profile");
        let tapered = DiscProfileSpec::SymmetricTapered {
            outer_radius_m: 0.038,
            face_radius_m: 0.015,
            thickness_m: 0.006,
        }
        .resolve(2_680.0, cx)
        .expect("tapered profile");
        let channels: CoupledChannelFactors = factors(0.038, 2_680.0).channel_factors();
        let short_controls = CoupledControls {
            maximum_steps: 20,
            ..controls()
        };
        let cylinder_run =
            run_closed_profile_reduced(&cylinder, channels, short_controls, initial(), None, cx)
                .expect("profile cylinder run");
        let tapered_run =
            run_closed_profile_reduced(&tapered, channels, short_controls, initial(), None, cx)
                .expect("profile tapered run");
        assert!(
            cylinder_run
                .samples
                .iter()
                .all(|sample| sample.support_source_feature.is_some())
        );
        assert_ne!(
            cylinder_run
                .samples
                .last()
                .unwrap()
                .spin_rad_per_s
                .to_bits(),
            tapered_run.samples.last().unwrap().spin_rad_per_s.to_bits()
        );
        assert_eq!(
            run_closed_profile_reduced(
                &tapered,
                channels,
                short_controls,
                initial(),
                Some(cylinder_run.checkpoint),
                cx,
            ),
            Err(CoupledError::CheckpointMismatch)
        );
    });
}
