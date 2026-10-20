//! One actual bow, shared stretching, both transverse planes, twist and fluid.
use super::support::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{
    SharedMaterialPlateBody, SharedPlateBridge,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    ReciprocalStringBank, StringBankOptions, StringBowLocation,
};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, ModalFrictionConfig,
};
use fs_couple::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig, PreparedModalBudget};
use fs_couple::thin_plate::{
    PlateMaterialModel, PlateThicknessConstraint, with_uniform_plate_material_state,
};
use fs_exec::{CancelGate, Cx};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::{
    ThinPlate,
    gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue},
};
use fs_tribo::FrictionLaw;

mod multiple;

fn nonlinear_config() -> NonlinearModalConfig {
    NonlinearModalConfig {
        max_iterations: 64,
        maximum_generalized_force: 1e5,
        force_absolute_tolerance: 1e-11,
        force_relative_tolerance: 1e-10,
    }
}
fn resolution() -> CylinderFlowResolution {
    CylinderFlowResolution {
        minimum_wavenumber_radius: 0.001,
        maximum_wavenumber_radius: 30.0,
        intervals: 4,
    }
}
fn physical_bank(
    cx: &Cx<'_>,
    gate: &CancelGate,
    young: f64,
    medium: GasState,
    strings: usize,
) -> ReciprocalStringBank {
    let sources = (0..strings)
        .map(|_| {
            source(7800.0, young, 1)
                .with_torsion("G", "eta_s", 1)
                .unwrap()
                .with_transverse_polarizations(
                    TransverseDirections::try_new([0.6, 0.8], [1.0, 0.0]).unwrap(),
                )
                .unwrap()
        })
        .collect();
    let body = SharedMaterialPlateBody {
        specimen: with_uniform_plate_material_state(
            ThinPlate {
                length_m: 0.4,
                width_m: 0.3,
                thickness_m: 1.0,
                density_kg_m3: 1.0,
                e1_pa: 1.0,
                e2_pa: 1.0,
                nu12: 0.0,
                g12_pa: 1.0,
                material_angle_rad: 0.0,
                damping_ratio: 0.0,
                thermoelastic: None,
                kelvin_voigt_bending: None,
                n_modes: 2,
                geometric_nonlinearity: false,
                pretension_n_m: 0.0,
                clamped: false,
            },
            &material(7800.0, 200e9),
            PlateMaterialModel::Isotropic,
            PlateThicknessConstraint::FixedThickness(0.001),
        )
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap(),
        bridges: vec![
            SharedPlateBridge {
                station_fraction: 0.8,
                stiffness_n_m: 3000.0,
                damping_n_s_m: 0.1,
                rest_extension_m: 0.0
            };
            strings
        ],
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    };
    let options = StringBankOptions {
        sample_rate_hz: 48_000,
        bow: StringBowLocation {
            string: 0,
            station_fraction: 0.2,
        },
        string_budget: ModalAcousticTimeBudget::audible_reference(),
        network: ModalCouplingConfig {
            max_modes: 64,
            max_connections: 64,
            max_setup_terms: 1_048_576,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1000.0,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
        },
        friction: ModalFrictionConfig {
            max_iterations: 96,
            maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-9,
            sticking_tolerance_m: 1e-12,
        },
    };
    ReciprocalStringBank::from_material_plate(
        cx,
        sources,
        medium,
        body,
        FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.4,
        },
        options,
        gate,
    )
    .unwrap()
    .with_prepared_linear_propagation(
        PreparedModalBudget {
            max_modes: 64,
            max_dense_bytes: 67_108_864,
            max_multiply_terms: 200_000_000_000,
        },
        gate,
    )
    .unwrap()
    .with_cylinder_gas_loading(resolution(), gate)
    .unwrap()
}
fn active() -> FrictionDrive {
    FrictionDrive {
        speed_m_s: 0.04,
        normal_force_n: 0.15,
    }
}
fn released() -> FrictionDrive {
    FrictionDrive {
        speed_m_s: 0.0,
        normal_force_n: 0.0,
    }
}

#[test]
fn joint_bow_stretching_and_fluid_memory_retain_all_coordinates_and_close_work() {
    with_cx(|cx, gate| {
        let original = physical_bank(cx, gate, 200e9, gas(), 2);
        let modes: Vec<_> = original
            .system()
            .components()
            .iter()
            .map(|m| m.modes().to_vec())
            .collect();
        let shapes = original.system().port().left.shapes.clone();
        let weights = original
            .system()
            .acceleration_pressure_weights()
            .unwrap()
            .to_vec();
        let mut bank = original
            .with_geometric_stretching(0.2, nonlinear_config(), gate)
            .unwrap();
        assert_eq!(
            bank.system()
                .components()
                .iter()
                .map(|m| m.modes().to_vec())
                .collect::<Vec<_>>(),
            modes
        );
        assert_eq!(bank.system().port().left.shapes, shapes);
        assert_eq!(
            bank.system().acceleration_pressure_weights(),
            Some(weights.as_slice())
        );
        assert_eq!(
            bank.source()
                .geometric_stretching()
                .unwrap()
                .potentials()
                .len(),
            2
        );
        assert_eq!(
            bank.source()
                .cylinder_gas_loading()
                .unwrap()
                .strings()
                .len(),
            2
        );
        let zero = vec![0.0; bank.system().mode_count()];
        let mut work = 0.0;
        let mut loss = 0.0;
        let mut stretch = 0.0_f64;
        for k in 0..320 {
            let frame = *bank
                .step_under_gate(&zero, if k < 192 { active() } else { released() }, gate)
                .unwrap();
            work += frame.drive_work_j + frame.external_work_j;
            loss += frame.network_dissipation_j + frame.friction_dissipation_j;
            let nonlinear = bank.system().nonlinear_frame().unwrap();
            assert!(nonlinear.reaction_residual <= nonlinear.reaction_tolerance);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert_eq!(nonlinear.friction, frame);
            stretch = stretch.max(nonlinear.nonlinear_energy_j);
            assert!((bank.system().total_energy_j().unwrap() + loss - work).abs() < 1e-10);
            for index in 0..2 {
                let card = &bank.cards()[index];
                let q = bank.system().components()[index].states();
                let mass = card.linear_density_kg_m * card.length_m * 0.5;
                let strain = 0.25
                    * (core::f64::consts::PI / card.length_m).powi(2)
                    * (q[0].displacement_m_sqrt_kg.powi(2) + q[1].displacement_m_sqrt_kg.powi(2))
                    / mass;
                let observed = bank.stretching_observation(index).unwrap();
                assert!((observed.additional_strain - strain).abs() < 1e-17);
                let ea = 200e9 * core::f64::consts::PI * 0.0003_f64.powi(2);
                assert!((observed.tension_n - (card.tension_n + ea * strain)).abs() < 1e-12);
            }
        }
        assert!(stretch > 0.0);
        assert!(loss > 0.0);
        assert!(
            bank.system().components()[1]
                .states()
                .iter()
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0),
            "unplayed string lost shared-body reaction"
        );
        assert!(
            bank.system().components()[0].states()[2].velocity_m_sqrt_kg_per_s != 0.0,
            "nonlinear string lost material twist"
        );
        assert!(
            bank.system().components()[3..]
                .iter()
                .flat_map(|c| c.states())
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0),
            "gas memory became detached metadata"
        );
    });
}
#[test]
fn source_stiffness_and_gas_viscosity_change_the_actual_combined_bow_pressure() {
    with_cx(|cx, gate| {
        let mut specification = GasSpec::dry_air_ussa1976();
        specification.sutherland_beta *= 4.0;
        let thick = GasState::try_new(&specification, 293.15, 101_325.0).unwrap();
        assert_eq!(thick.density, gas().density);
        let mut runs: Vec<_> = [(200e9, gas()), (50e9, gas()), (200e9, thick)]
            .into_iter()
            .map(|(young, medium)| {
                physical_bank(cx, gate, young, medium, 1)
                    .with_geometric_stretching(0.2, nonlinear_config(), gate)
                    .unwrap()
            })
            .collect();
        let zero = vec![0.0; runs[0].system().mode_count()];
        let mut pressure = vec![Vec::new(); runs.len()];
        for k in 0..192 {
            for (run, out) in runs.iter_mut().zip(&mut pressure) {
                out.push(
                    run.step(&zero, if k < 128 { active() } else { released() })
                        .unwrap()
                        .observer_pressure_pa,
                );
            }
        }
        assert_eq!(
            runs[0].system().components()[1].modes(),
            runs[1].system().components()[1].modes()
        );
        assert_eq!(
            runs[0].cards()[0].fundamental_hz(),
            runs[2].cards()[0].fundamental_hz()
        );
        for index in 1..3 {
            assert!(
                pressure[0]
                    .iter()
                    .zip(&pressure[index])
                    .any(|(a, b)| (a - b).abs() > 1e-15)
            );
            assert_ne!(
                runs[0].system().components()[1].states(),
                runs[index].system().components()[1].states()
            );
            assert!(runs[index].system().nonlinear_energy_j().unwrap() > 0.0);
        }
    });
}
#[test]
fn original_gesture_scheduler_executes_the_same_nonlinear_bank_for_every_block_partition() {
    with_cx(|cx, gate| {
        let schedule = GestureSchedule::try_new(
            700,
            vec![GestureTrack {
                id: "bow".into(),
                target: GestureTarget::BowStroke { string: 0 },
                initial: GestureValue::Bow {
                    velocity_m_per_s: active().speed_m_s,
                    normal_force_n: active().normal_force_n,
                    station: 0.2,
                },
                events: vec![],
            }],
        )
        .unwrap();
        let build = || {
            physical_bank(cx, gate, 200e9, gas(), 1)
                .with_geometric_stretching(0.2, nonlinear_config(), gate)
                .unwrap()
        };
        let mut direct = build();
        let zero = vec![0.0; direct.system().mode_count()];
        let expected: Vec<_> = (0..129)
            .map(|_| direct.step(&zero, active()).unwrap().observer_pressure_pa)
            .collect();
        assert!(expected.iter().any(|p| p.abs() > 1e-15));
        for block in [1, 37, 128] {
            let mut performance = build()
                .into_scheduled(
                    &schedule,
                    "bow",
                    FrictionGestureConfig {
                        sample_rate_hz: 48_000,
                        samples: 129,
                        max_block: 128,
                        max_work: 1024,
                        max_events: 16,
                        station_fraction: 0.2,
                    },
                )
                .unwrap();
            assert!(performance.source().geometric_stretching().is_some());
            assert!(performance.source().cylinder_gas_loading().is_some());
            let mut output = vec![0.0; 129];
            for chunk in output.chunks_mut(block) {
                performance.renderer_mut().block(chunk).unwrap();
            }
            assert_eq!(output, expected);
        }
    });
}
#[test]
fn unconverged_and_out_of_envelope_trials_keep_complete_bank_state_and_clock() {
    with_cx(|cx, gate| {
        for (slope, config) in [
            (
                0.2,
                NonlinearModalConfig {
                    max_iterations: 1,
                    force_absolute_tolerance: 1e-30,
                    force_relative_tolerance: 1e-15,
                    ..nonlinear_config()
                },
            ),
            (1e-15, nonlinear_config()),
        ] {
            let mut bank = physical_bank(cx, gate, 200e9, gas(), 1)
                .with_geometric_stretching(slope, config, gate)
                .unwrap();
            let before: Vec<_> = bank
                .system()
                .components()
                .iter()
                .map(|c| c.states().to_vec())
                .collect();
            let zero = vec![0.0; bank.system().mode_count()];
            assert!(bank.step_under_gate(&zero, active(), gate).is_err());
            assert_eq!(
                bank.system()
                    .components()
                    .iter()
                    .map(|c| c.states().to_vec())
                    .collect::<Vec<_>>(),
                before
            );
            assert_eq!(bank.system().samples_rendered(), 0);
            assert!(bank.system().last_frame().is_none());
            assert!(bank.system().nonlinear_frame().is_none());
            assert!(bank.system().has_nonlinear_potentials());
        }
    });
}
