//! G1/G3/G4: sourced junction events retain complete original elastic networks.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;
use common::{Z, close, sources};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::SourceDmtSphereContact;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::{
    LinearNetworkEventConfig, LinearNetworkEventKind, LinearNetworkPhase,
};
use fs_couple::render::schedule::force::coupled::contact::adhesive_network_events::*;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
    PreparedModalBudget,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_ports::{ModuleWindowAudit, PortedModule};

const MASS: [f64; 2] = [2e-9, 3e-9];
const DT: f64 = 1e-4;
fn controls() -> AdhesiveNetworkEventModalConfig {
    AdhesiveNetworkEventModalConfig {
        events: LinearNetworkEventConfig {
            maximum_modes: 8,
            maximum_connections: 4,
            maximum_workspace_bytes: 64 * 1024,
            maximum_cells: 4096,
            maximum_events: 64,
            maximum_operator_evaluations: 200_000,
            maximum_root_evaluations: 200_000,
            maximum_root_iterations: 80,
            maximum_step_s: DT,
            event_time_tolerance_s: 1e-16,
            sticking_velocity_tolerance_m_s: 1e-14,
            maximum_abs_modal_coordinate: 1.0,
            maximum_abs_modal_velocity: 1.0,
            energy_absolute_tolerance_j: 2e-24,
            energy_relative_tolerance: 2e-9,
            maximum_numerical_work_uncertainty_j: 1e-21,
        },
        maximum_contact_force_n: 1.0,
        maximum_tangent_speed_m_s: 1000.0,
        maximum_accumulated_slip_m: 0.1,
    }
}
fn column() -> Vec<f64> {
    vec![
        1.0 / det::sqrt(MASS[0]),
        0.25 / det::sqrt(MASS[0]),
        -1.0 / det::sqrt(MASS[1]),
        -0.4 / det::sqrt(MASS[1]),
    ]
}
fn network(
    moving: bool,
    displacement_limit: f64,
    absolute: f64,
    relative: f64,
) -> CoupledModalSystem {
    network_with_prepared(moving, displacement_limit, absolute, relative, false)
}
fn network_with_prepared(
    moving: bool,
    displacement_limit: f64,
    absolute: f64,
    relative: f64,
    prepared: bool,
) -> CoupledModalSystem {
    let models = MASS
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let mut budget = ModalAcousticTimeBudget::audible_reference();
            budget.maximum_abs_displacement_m_sqrt_kg = displacement_limit;
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                10_000,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    },
                    ModalAcousticMode {
                        angular_frequency_rad_s: [1800.0, 2700.0][i],
                        damping_ratio: [0.025, 0.04][i],
                        pressure_per_modal_velocity: C64 { re: 0.5, im: 0.1 },
                    },
                ],
                budget,
            )
            .unwrap();
            model
                .restore_states(&[
                    ModalAcousticState {
                        displacement_m_sqrt_kg: det::sqrt(mass) * [2e-8, -1e-8][i],
                        velocity_m_sqrt_kg_per_s: if moving {
                            det::sqrt(mass) * [0.0015, 0.0007][i]
                        } else {
                            0.0
                        },
                    },
                    ModalAcousticState {
                        displacement_m_sqrt_kg: det::sqrt(mass) * [3e-8, -2e-8][i],
                        velocity_m_sqrt_kg_per_s: if moving {
                            det::sqrt(mass) * [0.0003, -0.0005][i]
                        } else {
                            0.0
                        },
                    },
                ])
                .unwrap();
            model
        })
        .collect();
    let links = [
        ([1.0, 0.2], [1.0, -0.1], 0.04, 2e-7, 1e-8),
        ([1.0, 1.0], [1.0, -0.7], 0.015, 1e-7, -0.5e-8),
    ]
    .into_iter()
    .map(|(left, right, k, c, rest)| ModalConnection {
        left: ModalAttachment {
            component: 0,
            shapes: left.map(|x| x / det::sqrt(MASS[0])).to_vec(),
        },
        right: ModalAttachment {
            component: 1,
            shapes: right.map(|x| x / det::sqrt(MASS[1])).to_vec(),
        },
        stiffness_n_m: k,
        damping_n_s_m: c,
        rest_extension_m: rest,
    })
    .collect();
    let config = ModalCouplingConfig {
        max_modes: 4,
        max_connections: 2,
        max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1e-5,
        maximum_abs_pressure_pa: 1.0,
        maximum_abs_connection_force_n: 1.0,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: absolute,
        energy_relative_tolerance: relative,
    };
    if prepared {
        CoupledModalSystem::new_prepared(
            models,
            links,
            config,
            PreparedModalBudget {
                max_modes: 4,
                max_dense_bytes: 1024 * 1024,
                max_multiply_terms: 100_000_000,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    } else {
        CoupledModalSystem::new(models, links, config, &CancelGate::new_clock_free()).unwrap()
    }
}
fn owner<'a>(
    normal: &'a SourceDmtSphereContact,
    shear: &'a SourceAdhesiveShear,
    moving: bool,
    displacement_limit: f64,
    absolute: f64,
    relative: f64,
) -> AdhesiveNetworkEventModalSystem<'a> {
    let materials = normal
        .normal_state()
        .elastic()
        .interface()
        .surface_materials();
    AdhesiveNetworkEventModalSystem::new(
        network(moving, displacement_limit, absolute, relative),
        normal,
        shear,
        [&materials[0], &materials[1]],
        0.5 * Z,
        column(),
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn states(models: &[ModalAcousticTimeModel]) -> Vec<ModalAcousticState> {
    models
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn momentum(models: &[ModalAcousticTimeModel]) -> f64 {
    MASS.into_iter()
        .zip(models)
        .map(|(mass, model)| det::sqrt(mass) * model.states()[0].velocity_m_sqrt_kg_per_s)
        .sum()
}
fn compare(a: f64, b: f64) {
    close(a, b, 2e-24, 4e-9);
}

#[test]
fn every_original_elastic_mode_and_nonparallel_link_matches_the_existing_exact_owner() {
    let (normal, shear) = sources(1e-5, 0.0, 0.0);
    let mut system = owner(&normal, &shear, true, 1.0, 2e-24, 3e-9);
    let mut reference = network_with_prepared(true, 1.0, 2e-24, 3e-9, true);
    let external = [2e-5, -3e-6, -1e-5, 7e-6];
    for _ in 0..12 {
        let frame = *system.step(&external).unwrap();
        let expected = reference.step(&external).unwrap().clone();
        for (actual, expected) in states(system.components())
            .iter()
            .zip(states(reference.components()))
        {
            compare(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
            );
            compare(
                actual.velocity_m_sqrt_kg_per_s,
                expected.velocity_m_sqrt_kg_per_s,
            );
        }
        compare(frame.network_energy_j, expected.total_energy_j());
        compare(
            frame.modal_damping_dissipation_j,
            expected.component_dissipation_j,
        );
        compare(
            frame.connection_damping_dissipation_j,
            expected.connection_dissipation_j,
        );
        compare(frame.external_work_j, expected.external_work_j);
        compare(frame.observer_pressure_pa, expected.observer_pressure_pa);
        assert_eq!(frame.friction_dissipation_j, 0.0);
        assert!(frame.modal_damping_dissipation_j > 0.0);
        assert!(frame.connection_damping_dissipation_j > 0.0);
        assert_ne!(frame.observer_pressure_pa, 0.0);
        for (actual, expected) in system
            .network()
            .last_frame()
            .unwrap()
            .connection_forces_n
            .iter()
            .zip(&expected.connection_forces_n)
        {
            compare(*actual, *expected);
        }
    }
}

#[test]
fn sourced_stopping_and_loaded_restart_retain_all_modes_momentum_and_original_ports() {
    let (normal, shear) = sources(1e-5, 2e7, 1e7);
    let mut system = owner(&normal, &shear, true, 1.0, 2e-24, 3e-9);
    let initial = system.total_energy_j().unwrap();
    let p0 = momentum(system.components());
    let frame = *system.step(&[0.0; 4]).unwrap();
    assert!(
        frame
            .path
            .events
            .iter()
            .flatten()
            .any(|event| event.kind == LinearNetworkEventKind::Stopped)
    );
    assert_eq!(system.phase(), LinearNetworkPhase::Sticking);
    assert!(frame.friction_dissipation_j > 0.0);
    compare(momentum(system.components()), p0);
    compare(
        system.total_energy_j().unwrap()
            + frame.friction_dissipation_j
            + frame.modal_damping_dissipation_j
            + frame.connection_damping_dissipation_j,
        initial,
    );
    let ModuleWindowAudit::Available(audit) = system.network().audit_last_window().unwrap() else {
        panic!("original ports unavailable")
    };
    assert!(audit.violations().is_empty());
    let external: Vec<_> = column()
        .into_iter()
        .map(|b| 3.0 * frame.static_capacity_n * b)
        .collect();
    let before = system.total_energy_j().unwrap();
    let restart = *system.step(&external).unwrap();
    assert!(restart.path.events.iter().flatten().any(|event| matches!(
        event.kind,
        LinearNetworkEventKind::SlipStarted { direction: 1 }
    )));
    assert!(restart.friction_dissipation_j > 0.0);
    compare(momentum(system.components()), p0);
    compare(
        system.total_energy_j().unwrap()
            + restart.friction_dissipation_j
            + restart.modal_damping_dissipation_j
            + restart.connection_damping_dissipation_j
            - restart.external_work_j,
        before,
    );
    assert_eq!(system.samples_rendered(), 2);
}

#[test]
fn a_stuck_junction_constrains_only_its_coordinate_while_internal_elastic_motion_continues() {
    let (normal, shear) = sources(1e-5, 2e10, 1e10);
    let mut system = owner(&normal, &shear, false, 1.0, 2e-24, 3e-9);
    let before = states(system.components());
    let x0: f64 = column()
        .iter()
        .zip(&before)
        .map(|(b, state)| b * state.displacement_m_sqrt_kg)
        .sum();
    let frame = *system.step(&[0.0; 4]).unwrap();
    assert_eq!(frame.path.phase, LinearNetworkPhase::Sticking);
    assert_eq!(frame.friction_dissipation_j, 0.0);
    assert_eq!(frame.path.accumulated_slip_m, 0.0);
    close(frame.tangent_position_m, x0, 1e-20, 1e-12);
    assert!(frame.tangent_velocity_m_s.abs() < 1e-14);
    assert_ne!(states(system.components()), before);
    assert!(
        system.components()[0].states()[1]
            .velocity_m_sqrt_kg_per_s
            .abs()
            > 0.0
    );
    assert!(frame.modal_damping_dissipation_j > 0.0);
    compare(momentum(system.components()), 0.0);
}

#[test]
fn late_original_component_refusal_and_cancellation_leave_every_state_and_phase_retryable() {
    let (normal, shear) = sources(1e-5, 2e10, 1e10);
    let (mut system, mut reference) = (
        owner(&normal, &shear, false, 2e-12, 2e-24, 3e-9),
        owner(&normal, &shear, false, 2e-12, 2e-24, 3e-9),
    );
    system.step(&[0.0; 4]).unwrap();
    reference.step(&[0.0; 4]).unwrap();
    let old = (
        states(system.components()),
        system.phase(),
        *system.last_frame().unwrap(),
        system.network().last_frame().cloned(),
    );
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(system.step_under_gate(&[0.0; 4], &cancel).is_err());
    // Equal physical translation acceleration leaves the tangent unforced,
    // but exceeds an ORIGINAL component's position budget after propagation.
    let common = [
        1000.0 * det::sqrt(MASS[0]),
        0.0,
        1000.0 * det::sqrt(MASS[1]),
        0.0,
    ];
    assert!(matches!(
        system.step(&common),
        Err(AdhesiveNetworkEventModalError::Network(
            ModalCouplingError::Component { .. }
        ))
    ));
    assert_eq!(
        (
            states(system.components()),
            system.phase(),
            *system.last_frame().unwrap(),
            system.network().last_frame().cloned()
        ),
        old
    );
    assert_eq!(system.samples_rendered(), 1);
    assert_eq!(
        *system.step(&[0.0; 4]).unwrap(),
        *reference.step(&[0.0; 4]).unwrap()
    );
    assert_eq!(states(system.components()), states(reference.components()));
}

#[test]
fn numerical_trajectory_uncertainty_must_fit_the_original_native_energy_budget() {
    let (normal, shear) = sources(1e-5, 0.0, 0.0);
    let mut ordinary = owner(&normal, &shear, true, 1.0, 2e-24, 3e-9);
    let frame = *ordinary.step(&[0.0; 4]).unwrap();
    assert!(frame.numerical_work_uncertainty_j > 0.0);
    let mut strict = owner(
        &normal,
        &shear,
        true,
        1.0,
        0.5 * frame.numerical_work_uncertainty_j,
        1e-30,
    );
    let before = states(strict.components());
    assert!(matches!(
        strict.step(&[0.0; 4]),
        Err(AdhesiveNetworkEventModalError::EnergyAccuracy { .. })
    ));
    assert_eq!(states(strict.components()), before);
    assert_eq!(strict.phase(), LinearNetworkPhase::Fresh);
    assert_eq!(strict.samples_rendered(), 0);
    assert!(strict.last_frame().is_none());
    assert!(strict.network().last_frame().is_none());
}

#[test]
fn source_identity_and_complete_contact_column_are_required_at_admission() {
    let (normal, shear) = sources(1e-5, 2e7, 1e7);
    let (_, alien) = sources(2e-5, 2e7, 1e7);
    let materials = normal
        .normal_state()
        .elastic()
        .interface()
        .surface_materials();
    for (wrong_source, swapped, bad_column) in [
        (true, false, false),
        (false, true, false),
        (false, false, true),
    ] {
        let pair = if swapped {
            [&materials[1], &materials[0]]
        } else {
            [&materials[0], &materials[1]]
        };
        assert!(
            AdhesiveNetworkEventModalSystem::new(
                network(false, 1.0, 2e-24, 3e-9),
                &normal,
                if wrong_source { &alien } else { &shear },
                pair,
                0.5 * Z,
                if bad_column { vec![0.0; 4] } else { column() },
                controls(),
                &CancelGate::new_clock_free()
            )
            .is_err()
        );
    }
}
