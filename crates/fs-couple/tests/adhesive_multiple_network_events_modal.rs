//! G1/G3/G4: original bodies and independently sourced junctions share events.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;
use common::{Z, close, sources};
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::dmt::SourceDmtSphereContact;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::{
    LinearNetworkEventConfig, LinearNetworkEventKind, LinearNetworkPhase};
use fs_couple::render::schedule::force::coupled::contact::adhesive_events::linear_network::multiple::MultipleLinearNetworkEventConfig;
use fs_couple::render::schedule::force::coupled::contact::adhesive_network_events::{
    AdhesiveNetworkEventModalConfig, AdhesiveNetworkEventModalError,
    AdhesiveNetworkEventModalSystem, multiple::*};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment,
    ModalConnection, ModalCouplingConfig, ModalCouplingError, PreparedModalBudget};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_ports::{ModuleWindowAudit, PortedModule};

const MASS: [f64; 3] = [2e-9, 3e-9, 4e-9];
const DT: f64 = 1e-4;
fn controls() -> MultipleAdhesiveNetworkEventModalConfig {
    MultipleAdhesiveNetworkEventModalConfig {
        events: MultipleLinearNetworkEventConfig {
            network: LinearNetworkEventConfig {
                maximum_modes: 8,
                maximum_connections: 4,
                maximum_workspace_bytes: 256 * 1024,
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
            maximum_contacts: 2,
            minimum_column_pivot: 1e-10,
            maximum_phase_trials: 128,
        },
        maximum_contact_force_n: 1.0,
        maximum_tangent_speed_m_s: 1000.0,
        maximum_accumulated_slip_m: 0.1,
    }
}
fn columns(elastic: bool) -> Vec<Vec<f64>> {
    let width = if elastic { 2 } else { 1 };
    let mut result = vec![vec![0.0; 3 * width]; 2];
    for p in 0..2 {
        result[p][p * width] = 1.0 / det::sqrt(MASS[p]);
        result[p][(p + 1) * width] = -1.0 / det::sqrt(MASS[p + 1]);
        if elastic {
            result[p][p * width + 1] = [0.25, -0.3][p] / det::sqrt(MASS[p]);
            result[p][(p + 1) * width + 1] = [-0.4, 0.2][p] / det::sqrt(MASS[p + 1]);
        }
    }
    result
}
fn network(
    elastic: bool,
    moving: bool,
    displacement_limit: f64,
    absolute: f64,
    relative: f64,
) -> CoupledModalSystem {
    let models = MASS
        .into_iter()
        .enumerate()
        .map(|(i, mass)| {
            let mut budget = ModalAcousticTimeBudget::audible_reference();
            budget.maximum_abs_displacement_m_sqrt_kg = displacement_limit;
            let mut modes = vec![ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }];
            if elastic {
                modes.push(ModalAcousticMode {
                    angular_frequency_rad_s: [1800.0, 2300.0, 2700.0][i],
                    damping_ratio: 0.03,
                    pressure_per_modal_velocity: C64 { re: 0.5, im: 0.1 },
                });
            }
            let mut model =
                ModalAcousticTimeModel::try_new_with_free_coordinates(10_000, modes, budget)
                    .unwrap();
            let mut states = vec![ModalAcousticState {
                displacement_m_sqrt_kg: if elastic {
                    det::sqrt(mass) * [2e-8, -1e-8, 1e-8][i]
                } else {
                    0.0
                },
                velocity_m_sqrt_kg_per_s: if moving {
                    det::sqrt(mass) * [0.002, 0.0, -0.001][i]
                } else {
                    0.0
                },
            }];
            if elastic {
                states.push(ModalAcousticState {
                    displacement_m_sqrt_kg: det::sqrt(mass) * [3e-8, -2e-8, 1e-8][i],
                    velocity_m_sqrt_kg_per_s: if moving {
                        det::sqrt(mass) * [0.0003, -0.0005, 0.0002][i]
                    } else {
                        0.0
                    },
                });
            }
            model.restore_states(&states).unwrap();
            model
        })
        .collect();
    let links = if elastic {
        (0..2)
            .map(|i| ModalConnection {
                left: ModalAttachment {
                    component: i,
                    shapes: vec![1.0 / det::sqrt(MASS[i]), 0.2 / det::sqrt(MASS[i])],
                },
                right: ModalAttachment {
                    component: i + 1,
                    shapes: vec![1.0 / det::sqrt(MASS[i + 1]), -0.1 / det::sqrt(MASS[i + 1])],
                },
                stiffness_n_m: [0.04, 0.025][i],
                damping_n_s_m: [2e-7, 1e-7][i],
                rest_extension_m: [1e-8, -0.5e-8][i],
            })
            .collect()
    } else {
        vec![]
    };
    CoupledModalSystem::new_prepared(
        models,
        links,
        ModalCouplingConfig {
            max_modes: 6,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-5,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: absolute,
            energy_relative_tolerance: relative,
        },
        PreparedModalBudget {
            max_modes: 6,
            max_dense_bytes: 1024 * 1024,
            max_multiply_terms: 100_000_000,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn junctions<'a>(
    normal: [&'a SourceDmtSphereContact; 2],
    shear: [&'a SourceAdhesiveShear; 2],
    elastic: bool,
) -> Vec<AdhesiveNetworkEventJunction<'a>> {
    columns(elastic)
        .into_iter()
        .enumerate()
        .map(|(i, column)| {
            let materials = normal[i]
                .normal_state()
                .elastic()
                .interface()
                .surface_materials();
            AdhesiveNetworkEventJunction {
                normal: normal[i],
                shear: shear[i],
                surface_materials: [&materials[0], &materials[1]],
                indentation_m: 0.5 * Z,
                column,
            }
        })
        .collect()
}
fn owner<'a>(
    normal: [&'a SourceDmtSphereContact; 2],
    shear: [&'a SourceAdhesiveShear; 2],
    elastic: bool,
    moving: bool,
) -> MultipleAdhesiveNetworkEventModalSystem<'a> {
    MultipleAdhesiveNetworkEventModalSystem::new(
        network(elastic, moving, 1.0, 2e-24, 3e-9),
        junctions(normal, shear, elastic),
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
fn simultaneous_stops_and_coupled_restart_keep_original_momentum_energy_and_ports() {
    let (n0, s0) = sources(1e-5, 2e8, 1e8);
    let (n1, s1) = sources(2e-5, 2e8, 1e8);
    let mut system = owner([&n0, &n1], [&s0, &s1], false, true);
    let before = system.total_energy_j().unwrap();
    let p0 = momentum(system.components());
    let frame = *system.step(&[0.0; 3]).unwrap();
    assert_eq!(system.phases(), &[LinearNetworkPhase::Sticking; 2]);
    for i in 0..2 {
        assert!(
            frame
                .path
                .events
                .iter()
                .flatten()
                .any(|e| e.contact_index == i && e.kind == LinearNetworkEventKind::Stopped)
        );
        assert!(frame.path.contact_friction_dissipation_j[i] > 0.0);
        assert!(frame.tangent_velocities_m_s[i].abs() < 1e-14);
    }
    compare(momentum(system.components()), p0);
    compare(
        system.total_energy_j().unwrap() + frame.friction_dissipation_j,
        before,
    );
    compare(
        frame.path.contact_friction_dissipation_j[..2].iter().sum(),
        frame.friction_dissipation_j,
    );
    let ModuleWindowAudit::Available(audit) = system.network().audit_last_window().unwrap() else {
        panic!("native ports missing")
    };
    assert!(audit.violations().is_empty());
    let external: Vec<_> = columns(false)[0]
        .iter()
        .map(|b| b * 3.0 * frame.static_capacities_n[0])
        .collect();
    let before = system.total_energy_j().unwrap();
    let restart = *system.step(&external).unwrap();
    assert!(
        restart
            .path
            .events
            .iter()
            .flatten()
            .any(|e| matches!(e.kind, LinearNetworkEventKind::SlipStarted { .. }))
    );
    assert!(restart.friction_dissipation_j > 0.0);
    compare(momentum(system.components()), p0);
    compare(
        system.total_energy_j().unwrap() + restart.friction_dissipation_j - restart.external_work_j,
        before,
    );
}

#[test]
fn both_held_junctions_preserve_internal_elastic_motion_and_each_original_loss_channel() {
    let (n0, s0) = sources(1e-5, 2e8, 1e8);
    let (n1, s1) = sources(2e-5, 3e8, 1.5e8);
    let mut system = owner([&n0, &n1], [&s0, &s1], true, false);
    let initial = states(system.components());
    let b = columns(true);
    let x0: Vec<f64> = b
        .iter()
        .map(|column| {
            column
                .iter()
                .zip(&initial)
                .map(|(b, s)| b * s.displacement_m_sqrt_kg)
                .sum()
        })
        .collect();
    for _ in 0..8 {
        let frame = *system.step(&[0.0; 6]).unwrap();
        assert_eq!(system.phases(), &[LinearNetworkPhase::Sticking; 2]);
        assert_eq!(frame.friction_dissipation_j, 0.0);
        assert!(frame.modal_damping_dissipation_j > 0.0);
        assert!(frame.connection_damping_dissipation_j > 0.0);
        for i in 0..2 {
            close(frame.tangent_positions_m[i], x0[i], 2e-20, 1e-11);
        }
        compare(momentum(system.components()), 0.0);
    }
    assert_ne!(states(system.components()), initial);
}

#[test]
fn zero_junction_strength_recovers_the_original_complete_prepared_network() {
    let (n0, s0) = sources(1e-5, 0.0, 0.0);
    let (n1, s1) = sources(2e-5, 0.0, 0.0);
    let mut system = owner([&n0, &n1], [&s0, &s1], true, true);
    let mut reference = network(true, true, 1.0, 2e-24, 3e-9);
    let external = [2e-5, -3e-6, -1e-5, 7e-6, 4e-6, -2e-6];
    for _ in 0..12 {
        let frame = *system.step(&external).unwrap();
        let expected = reference.step(&external).unwrap().clone();
        for (a, b) in states(system.components())
            .iter()
            .zip(states(reference.components()))
        {
            compare(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
            compare(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
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
        compare(frame.observer_pressure_pa, expected.observer_pressure_pa);
        compare(frame.external_work_j, expected.external_work_j);
        assert_eq!(frame.friction_dissipation_j, 0.0);
    }
}

#[test]
fn contact_order_changes_indices_without_changing_the_joint_physical_motion() {
    let (n0, s0) = sources(1e-5, 2e7, 1e7);
    let (n1, s1) = sources(2e-5, 3e7, 1.5e7);
    let mut a = owner([&n0, &n1], [&s0, &s1], true, true);
    let mut sites = junctions([&n0, &n1], [&s0, &s1], true);
    sites.swap(0, 1);
    let mut b = MultipleAdhesiveNetworkEventModalSystem::new(
        network(true, true, 1.0, 2e-24, 3e-9),
        sites,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for _ in 0..6 {
        let fa = *a.step(&[0.0; 6]).unwrap();
        let fb = *b.step(&[0.0; 6]).unwrap();
        for (a, b) in states(a.components()).iter().zip(states(b.components())) {
            compare(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
            compare(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
        }
        compare(fa.friction_dissipation_j, fb.friction_dissipation_j);
        for i in 0..2 {
            compare(
                fa.path.contact_friction_dissipation_j[i],
                fb.path.contact_friction_dissipation_j[1 - i],
            );
            compare(fa.tangent_positions_m[i], fb.tangent_positions_m[1 - i]);
        }
    }
}

#[test]
fn adding_an_independent_zero_capacity_junction_preserves_the_original_single_junction_owner() {
    let (n0, s0) = sources(1e-5, 2e7, 1e7);
    let (n1, s1) = sources(2e-5, 0.0, 0.0);
    let mut multiple = owner([&n0, &n1], [&s0, &s1], true, true);
    let materials = n0.normal_state().elastic().interface().surface_materials();
    let c = controls();
    let mut single = AdhesiveNetworkEventModalSystem::new(
        network(true, true, 1.0, 2e-24, 3e-9),
        &n0,
        &s0,
        [&materials[0], &materials[1]],
        0.5 * Z,
        columns(true)[0].clone(),
        AdhesiveNetworkEventModalConfig {
            events: c.events.network,
            maximum_contact_force_n: c.maximum_contact_force_n,
            maximum_tangent_speed_m_s: c.maximum_tangent_speed_m_s,
            maximum_accumulated_slip_m: c.maximum_accumulated_slip_m,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    for _ in 0..6 {
        let a = *multiple.step(&[0.0; 6]).unwrap();
        let b = *single.step(&[0.0; 6]).unwrap();
        for (a, b) in states(multiple.components())
            .iter()
            .zip(states(single.components()))
        {
            compare(a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg);
            compare(a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s);
        }
        compare(a.friction_dissipation_j, b.friction_dissipation_j);
        compare(a.path.contact_impulse_n_s[0], b.path.contact_impulse_n_s);
        assert_eq!(a.path.contact_friction_dissipation_j[1], 0.0);
    }
}

#[test]
fn late_original_body_failure_and_cancellation_preserve_all_phases_and_identical_retry() {
    let (n0, s0) = sources(1e-5, 2e8, 1e8);
    let (n1, s1) = sources(2e-5, 2e8, 1e8);
    let make = || {
        MultipleAdhesiveNetworkEventModalSystem::new(
            network(true, false, 2e-12, 2e-24, 3e-9),
            junctions([&n0, &n1], [&s0, &s1], true),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    };
    let (mut system, mut reference) = (make(), make());
    system.step(&[0.0; 6]).unwrap();
    reference.step(&[0.0; 6]).unwrap();
    let before = (
        states(system.components()),
        system.phases().to_vec(),
        *system.last_frame().unwrap(),
        system.network().last_frame().cloned(),
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(system.step_under_gate(&[0.0; 6], &gate).is_err());
    let common: Vec<_> = MASS
        .into_iter()
        // Exceed the original displacement limit while remaining within the
        // independent event-propagation work-uncertainty budget.
        .flat_map(|m| [100.0 * det::sqrt(m), 0.0])
        .collect();
    let refusal = system.step(&common).unwrap_err();
    assert!(
        matches!(
            refusal,
            AdhesiveNetworkEventModalError::Network(ModalCouplingError::Component { .. })
        ),
        "unexpected refusal: {refusal:?}"
    );
    assert_eq!(
        (
            states(system.components()),
            system.phases().to_vec(),
            *system.last_frame().unwrap(),
            system.network().last_frame().cloned()
        ),
        before
    );
    assert_eq!(system.samples_rendered(), 1);
    assert_eq!(
        *system.step(&[0.0; 6]).unwrap(),
        *reference.step(&[0.0; 6]).unwrap()
    );
    assert_eq!(states(system.components()), states(reference.components()));
}

#[test]
fn original_energy_accuracy_rejects_every_junction_history_together() {
    let (n0, s0) = sources(1e-5, 0.0, 0.0);
    let (n1, s1) = sources(2e-5, 0.0, 0.0);
    let mut ordinary = owner([&n0, &n1], [&s0, &s1], true, true);
    let frame = *ordinary.step(&[0.0; 6]).unwrap();
    assert!(frame.numerical_work_uncertainty_j > 0.0);
    let mut strict = MultipleAdhesiveNetworkEventModalSystem::new(
        network(
            true,
            true,
            1.0,
            0.5 * frame.numerical_work_uncertainty_j,
            1e-30,
        ),
        junctions([&n0, &n1], [&s0, &s1], true),
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let before = states(strict.components());
    assert!(matches!(
        strict.step(&[0.0; 6]),
        Err(AdhesiveNetworkEventModalError::EnergyAccuracy { .. })
    ));
    assert_eq!(states(strict.components()), before);
    assert_eq!(strict.phases(), &[LinearNetworkPhase::Fresh; 2]);
    assert_eq!(strict.samples_rendered(), 0);
    assert!(strict.last_frame().is_none());
    assert!(strict.network().last_frame().is_none());
}

#[test]
fn each_original_source_and_independent_complete_work_column_is_required() {
    let (n0, s0) = sources(1e-5, 2e7, 1e7);
    let (n1, s1) = sources(2e-5, 2e7, 1e7);
    for defect in 0..4 {
        let mut sites = junctions([&n0, &n1], [&s0, &s1], true);
        match defect {
            0 => sites[1].shear = &s0,
            1 => sites[1].surface_materials.swap(0, 1),
            2 => sites[1].column = sites[0].column.clone(),
            _ => {
                sites[1].column.pop();
            }
        }
        assert!(
            MultipleAdhesiveNetworkEventModalSystem::new(
                network(true, false, 1.0, 2e-24, 3e-9),
                sites,
                controls(),
                &CancelGate::new_clock_free()
            )
            .is_err()
        );
    }
}
