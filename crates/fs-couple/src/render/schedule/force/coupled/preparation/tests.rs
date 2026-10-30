use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::plate::impact::radiation::{Model, Pole};
use crate::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;
use crate::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort,
};
use fs_math::c64::C64;

fn budget(bytes: usize) -> PreparedModalBudget {
    PreparedModalBudget {
        max_modes: 8,
        max_dense_bytes: bytes,
        max_multiply_terms: 100_000_000_000,
    }
}

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 8,
        max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9,
    }
}

fn models() -> Vec<ModalAcousticTimeModel> {
    vec![
        ModalAcousticTimeModel::try_new(
            48_000,
            vec![
                ModalAcousticMode {
                    angular_frequency_rad_s: 500.0,
                    damping_ratio: 0.02,
                    pressure_per_modal_velocity: C64::new(0.3, -0.1),
                },
                ModalAcousticMode {
                    angular_frequency_rad_s: 900.0,
                    damping_ratio: 0.03,
                    pressure_per_modal_velocity: C64::new(-0.2, 0.04),
                },
            ],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap(),
    ]
}

fn network(
    choice: ModalPreparation,
    gate: &CancelGate,
) -> Result<CoupledModalSystem, ModalCouplingError> {
    CoupledModalSystem::new_with_columns_preparation(
        models(),
        vec![
            ModalColumnConnection {
                column: vec![1.0, -0.4],
                stiffness_n_m: 50_000.0,
                damping_n_s_m: 12.0,
                rest_extension_m: 1e-5,
            },
            ModalColumnConnection {
                column: vec![0.3, 0.6],
                stiffness_n_m: 10_000.0,
                damping_n_s_m: 3.0,
                rest_extension_m: 0.0,
            },
        ],
        &[ConstantLossConnection {
            connection: 0,
            reference_omega_rad_s: 700.0,
        }],
        config(),
        choice,
        gate,
    )?
    .with_material_relaxation(
        vec![
            fs_phs::RelaxationBranch {
                projection: vec![1.0, 0.0, 0.2, 0.0],
                stiffness: 30_000.0,
                relaxation_time_s: 0.002,
            },
            fs_phs::RelaxationBranch {
                projection: vec![-0.3, 0.0, 1.0, 0.0],
                stiffness: 7_000.0,
                relaxation_time_s: 0.03,
            },
        ],
        Some(vec![0.004, -0.007]),
        gate,
    )?
    .with_radiation_load(
        &Model {
            ports: 2,
            poles: vec![Pole {
                omega: 700.0,
                zeta: 0.1,
                coupling: vec![100.0, -30.0],
            }],
        },
        &[0, 1],
        1,
        gate,
    )
}

#[test]
fn harmonic_only_retains_prepared_physics_without_a_time_propagator() {
    // G1/G3: the complete continuous response, including coupled constant-loss
    // damping directions, signed Maxwell arms and air, matches playback's
    // original image. The small analysis budget cannot admit its exponential.
    let gate = CancelGate::new_clock_free();
    let small = budget(32 * 1024);
    assert!(network(ModalPreparation::Prepared(small), &gate).is_err());
    let mut analysis = network(ModalPreparation::HarmonicOnly(small), &gate).unwrap();
    let mut playback =
        network(ModalPreparation::Prepared(budget(128 * 1024 * 1024)), &gate).unwrap();
    assert!(analysis.exact.is_none());
    assert!(playback.exact.is_some());
    assert_eq!(
        analysis.connection_damping_column(0),
        playback.connection_damping_column(0)
    );
    assert_ne!(
        analysis.connection_damping_column(0).unwrap(),
        analysis.columns[0]
    );
    let force = [C64::new(0.2, 0.03), C64::new(-0.1, 0.05)];
    for omega in [70.0, 550.0, 700.0, 1500.0] {
        assert_eq!(
            analysis
                .harmonic_response(omega, &force, 10_000, &gate)
                .unwrap(),
            playback
                .harmonic_response(omega, &force, 10_000, &gate)
                .unwrap()
        );
    }
    let update = [coefficients::ConnectionCoefficients {
        connection: 1,
        stiffness_n_m: 16_000.0,
        damping_n_s_m: 4.0,
        rest_extension_m: -1e-5,
    }];
    assert_eq!(
        analysis
            .set_connection_coefficients(&update, &gate)
            .unwrap(),
        playback
            .set_connection_coefficients(&update, &gate)
            .unwrap()
    );
    assert!(analysis.exact.is_none());
    assert!(analysis.is_harmonic_only());
    assert_eq!(
        analysis.connection_damping_column(0),
        playback.connection_damping_column(0)
    );
    assert_eq!(
        analysis
            .harmonic_response(650.0, &force, 10_000, &gate)
            .unwrap(),
        playback
            .harmonic_response(650.0, &force, 10_000, &gate)
            .unwrap()
    );
}

#[test]
fn harmonic_only_refuses_time_and_preserves_histories_on_every_rejection() {
    // G0/G4: no absent time propagator may fall back to uncoupled oscillators.
    let gate = CancelGate::new_clock_free();
    let mut analysis = network(ModalPreparation::HarmonicOnly(budget(32 * 1024)), &gate).unwrap();
    let state = analysis.models[0].states().to_vec();
    let memory = analysis.material_memory().to_vec();
    let air = analysis.radiation_state().to_vec();
    assert!(analysis.step(&[0.1, 0.0]).is_err());
    assert!(analysis.step_under_gate(&[0.1, 0.0], &gate).is_err());
    assert!(
        analysis
            .prepare_prediction(&[0.1, 0.0], Some(&gate))
            .is_err()
    );
    assert!(analysis.linear_response(&[1.0, 0.0], &gate).is_err());
    assert_eq!(analysis.models[0].states(), state);
    assert_eq!(analysis.material_memory(), memory);
    assert_eq!(analysis.radiation_state(), air);
    assert_eq!(analysis.samples_rendered(), 0);
    assert!(analysis.last_frame().is_none());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        network(
            ModalPreparation::HarmonicOnly(budget(32 * 1024)),
            &cancelled
        ),
        Err(ModalCouplingError::Cancelled)
    ));
    assert!(network(ModalPreparation::HarmonicOnly(budget(1)), &gate).is_err());
    assert!(
        network(
            ModalPreparation::HarmonicOnly(PreparedModalBudget {
                max_modes: 1,
                ..budget(32 * 1024)
            }),
            &gate
        )
        .is_err()
    );
}

#[test]
fn harmonic_only_friction_ports_keep_pressure_without_preparing_compliance() {
    // G1/G4: released pressure still observes the original connected solution;
    // single and multiple friction owners cannot consume analysis placeholders.
    let gate = CancelGate::new_clock_free();
    let analysis = network(ModalPreparation::HarmonicOnly(budget(32 * 1024)), &gate).unwrap();
    let port = ModalFrictionPort {
        left: ModalAttachment {
            component: 0,
            shapes: vec![0.8, -0.2],
        },
        right: None,
        law: fs_tribo::FrictionLaw::Coulomb {
            static_mu: 0.4,
            kinetic_mu: 0.3,
        },
    };
    let config = ModalFrictionConfig {
        max_iterations: 80,
        maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-10,
        sticking_tolerance_m: 1e-12,
    };
    let mut friction = FrictionModalSystem::new(analysis, port.clone(), config, &gate).unwrap();
    let force = [C64::from_re(0.2), C64::ZERO];
    let observed = friction
        .released_harmonic_response(650.0, &force, 10_000, &gate)
        .unwrap();
    assert!(observed.observer_pressure_pa.abs() > 0.0);
    assert!(friction.is_harmonic_only());
    let drive = FrictionDrive {
        speed_m_s: 0.1,
        normal_force_n: 1.0,
    };
    assert!(friction.step(&[0.1, 0.0], drive).is_err());
    assert_eq!(friction.samples_rendered(), 0);
    let mut multiple = friction
        .with_additional_ports(
            vec![(port, config)],
            MultiFrictionConfig {
                max_ports: 2,
                max_sweeps: 16,
                max_setup_terms: 10_000,
            },
            &gate,
        )
        .unwrap();
    assert_eq!(
        multiple
            .released_harmonic_response(650.0, &force, 10_000, &gate)
            .unwrap(),
        observed
    );
    assert!(multiple.step(&[0.1, 0.0], &[drive, drive]).is_err());
    assert_eq!(multiple.samples_rendered(), 0);
}
