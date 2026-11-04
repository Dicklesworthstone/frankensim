//! Original native spring geometry, work ports and inertial JKR branch events.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use fs_contact::normal_patch::adhesive::motion::{JkrNormalMotionBudget, MAX_JKR_NORMAL_EVENTS};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrEnvelope, JkrJumpKind, JkrState, SourceJkrSphereContact,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::adhesive_normal_events::*;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalPreparation,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use fs_ports::{ModuleWindowAudit, PortedModule};

const MA: f64 = 2e-5;
const MB: f64 = 3e-5;
const MU: f64 = MA * MB / (MA + MB);
const DT: f64 = 1e-4;

fn source() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}
fn scales(law: &SourceJkrSphereContact) -> (f64, f64, f64) {
    let d = law.indentation_scale_m();
    let f = core::f64::consts::PI * 0.05 * 0.01 / 6.0;
    (d, f, det::sqrt(f * d / MU))
}
fn config(law: &SourceJkrSphereContact, divisions: f64) -> JkrNormalEventModalConfig {
    let (d, f, v) = scales(law);
    JkrNormalEventModalConfig {
        motion: JkrNormalMotionBudget {
            maximum_step_s: DT,
            maximum_substep_s: DT / divisions,
            maximum_substeps: 2048,
            maximum_root_iterations: 64,
            maximum_events: MAX_JKR_NORMAL_EVENTS,
            minimum_indentation_m: -500.0 * d,
            maximum_accumulated_travel_m: 200.0 * d,
            maximum_normal_speed_m_s: 100.0 * v,
            maximum_force_n: 1.0,
            coordinate_absolute_tolerance_m: 1e-10 * d,
            coordinate_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-10 * f * d,
            energy_relative_tolerance: 1e-10,
            impulse_absolute_tolerance_n_s: 1e-10 * MU * v,
            impulse_relative_tolerance: 1e-10,
        },
        coordinate_absolute_tolerance_m: 2e-10 * d,
        coordinate_relative_tolerance: 2e-10,
        velocity_absolute_tolerance_m_s: 2e-10 * v,
        velocity_relative_tolerance: 2e-10,
    }
}

#[derive(Clone, Copy)]
struct Setup {
    k: f64,
    x: f64,
    extension: f64,
    velocity: f64,
    common_velocity: f64,
    damping: f64,
    links: usize,
    force_limit: f64,
    displacement_limit: f64,
}
impl Setup {
    fn new(k: f64, x: f64, extension: f64) -> Self {
        Self {
            k,
            x,
            extension,
            velocity: 0.0,
            common_velocity: 0.003,
            damping: 0.0,
            links: 1,
            force_limit: 1.0,
            displacement_limit: 1.0,
        }
    }
}
fn network(s: Setup) -> CoupledModalSystem {
    let models = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(side, mass)| {
            let factor = if side == 0 {
                MB / (MA + MB)
            } else {
                -MA / (MA + MB)
            };
            let mut budget = ModalAcousticTimeBudget::audible_reference();
            budget.maximum_abs_displacement_m_sqrt_kg = s.displacement_limit;
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                10_000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                budget,
            )
            .unwrap();
            model
                .restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: det::sqrt(mass) * factor * s.x,
                    velocity_m_sqrt_kg_per_s: det::sqrt(mass)
                        * (s.common_velocity + factor * s.velocity),
                }])
                .unwrap();
            model
        })
        .collect();
    let links = (0..s.links)
        .map(|_| ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0 / det::sqrt(MA)],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![1.0 / det::sqrt(MB)],
            },
            stiffness_n_m: s.k,
            damping_n_s_m: s.damping,
            rest_extension_m: s.x - s.extension,
        })
        .collect();
    CoupledModalSystem::new_with_preparation(
        models,
        links,
        &[],
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-4,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: s.force_limit,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-20,
            energy_relative_tolerance: 5e-10,
        },
        ModalPreparation::Ordinary,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn make<'a>(
    law: &'a SourceJkrSphereContact,
    contact: JkrState,
    setup: Setup,
    divisions: f64,
) -> JkrNormalEventModalSystem<'a> {
    let ids = law.normal_state().elastic().interface().surface_materials();
    JkrNormalEventModalSystem::new_with_elastic_connection(
        network(setup),
        law,
        [&ids[0], &ids[1]],
        contact,
        config(law, divisions),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}, scaled error {:e}",
        (actual - expected).abs() / scale
    );
}
fn states(owner: &JkrNormalEventModalSystem<'_>) -> Vec<ModalAcousticState> {
    owner
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn momentum(owner: &JkrNormalEventModalSystem<'_>) -> f64 {
    [MA, MB]
        .into_iter()
        .zip(owner.components())
        .map(|(m, body)| det::sqrt(m) * body.states()[0].velocity_m_sqrt_kg_per_s)
        .sum()
}
fn ports(owner: &JkrNormalEventModalSystem<'_>, spring_change: f64, scale: f64) {
    let f = owner.last_frame().unwrap();
    let native = owner.network().last_frame().unwrap();
    close(native.connection_energy_j, f.spring_energy_j, scale, 2e-12);
    assert_eq!(native.connection_forces_n, vec![f.connection_force_n]);
    close(
        native.external_work_j,
        f.external_work_j + f.normal_work_on_network_j,
        scale,
        2e-12,
    );
    let transfers = owner.network().last_window_transfers().unwrap();
    assert_eq!(transfers.len(), 1);
    close(
        transfers[0].energy_into().value(),
        -spring_change,
        scale,
        2e-12,
    );
    let ModuleWindowAudit::Available(audit) = owner.network().audit_last_window().unwrap() else {
        panic!("original native connection audit unavailable");
    };
    assert!(audit.violations().is_empty());
}

#[test]
fn detached_native_spring_converges_to_original_oscillator_and_retains_its_ports() {
    let law = source();
    let (d, f, v) = scales(&law);
    let k = 0.4 * f / d;
    let setup = Setup::new(k, 2.0 * d, 5.0 * d);
    let initial = law.detached_state(-20.0 * d).unwrap();
    let mut errors = Vec::new();
    for divisions in [2.0, 4.0] {
        let mut owner = make(&law, initial, setup, divisions);
        let energy = owner.total_energy_j().unwrap();
        let mut spring = 0.5 * k * setup.extension.powi(2);
        for _ in 0..8 {
            let frame = owner.step(&[0.0; 2]).unwrap().clone();
            assert_eq!(frame.path.event_count(), 0);
            assert!(frame.normal_work_on_network_j.abs() < 1e-10 * f * d);
            close(owner.total_energy_j().unwrap(), energy, energy, 2e-9);
            close(
                momentum(&owner),
                (MA + MB) * setup.common_velocity,
                (MA + MB) * v,
                2e-11,
            );
            ports(&owner, frame.spring_energy_j - spring, f * d);
            spring = frame.spring_energy_j;
        }
        let omega = det::sqrt(k / MU);
        let time = 8.0 * DT;
        let exact = -25.0 * d + 5.0 * d * det::cos(omega * time);
        let exact_velocity = -5.0 * d * omega * det::sin(omega * time);
        let error = ((owner.state().contact_state().indentation_m() - exact) / d).abs()
            + ((owner.state().velocity_m_s() - exact_velocity) / v).abs();
        errors.push(error);
        assert_eq!(owner.state().contact_state().branch(), JkrBranch::Detached);
    }
    assert!(
        errors[0] > 3.5 * errors[1],
        "second-order oscillator refinement: {errors:?}"
    );
    assert!(errors[1] < 0.02, "refined oscillator error: {errors:?}");
}

#[test]
fn original_preload_balances_tensile_adhesion_and_anchors_actual_native_geometry() {
    let law = source();
    let (d, f, v) = scales(&law);
    let contact = law.evaluate_attached(0.0).unwrap();
    let k = 2.0 * f / d;
    let extension = -contact.normal_force_n / k;
    let setup = Setup::new(k, 7.0 * d, extension);
    let mut owner = make(&law, contact.state(), setup, 20.0);
    close(
        owner.support_zero_force_indentation_m().unwrap(),
        -extension,
        d,
        2e-14,
    );
    assert_eq!(owner.support_stiffness_n_per_m(), k);
    close(owner.effective_mass_kg(), MU, MU, 2e-15);
    let before = owner.total_energy_j().unwrap();
    let frame = owner.step(&[0.0; 2]).unwrap().clone();
    close(owner.state().contact_state().indentation_m(), 0.0, d, 3e-9);
    close(owner.state().velocity_m_s(), 0.0, v, 3e-9);
    close(
        frame.path.normal_impulse_n_s,
        contact.normal_force_n * DT,
        f * DT,
        3e-9,
    );
    close(
        frame.path.support_impulse_n_s,
        -contact.normal_force_n * DT,
        f * DT,
        3e-9,
    );
    close(owner.total_energy_j().unwrap(), before, before.abs(), 3e-10);
    ports(&owner, 0.0, f * d);
    let mut shifted = make(
        &law,
        contact.state(),
        Setup {
            extension: extension + d,
            ..setup
        },
        20.0,
    );
    shifted.step(&[0.0; 2]).unwrap();
    assert!(shifted.state().velocity_m_s() < 0.0);
    close(
        shifted.support_zero_force_indentation_m().unwrap(),
        -extension - d,
        d,
        2e-14,
    );
}

#[test]
fn original_spring_load_changes_capture_and_pull_off_with_native_energy_and_momentum() {
    let law = source();
    let (d, f, v) = scales(&law);
    for (attached, extension, expected) in [
        (false, -4.0 * d, JkrJumpKind::Attach),
        (true, 30.0 * d, JkrJumpKind::Detach),
    ] {
        let state = if attached {
            law.evaluate_attached(0.0).unwrap().state()
        } else {
            law.detached_state(-d).unwrap()
        };
        let setup = Setup::new(0.5 * f / d, 2.0 * d, extension);
        let mut owner = make(&law, state, setup, 20.0);
        let energy = owner.total_energy_j().unwrap();
        let mut released = 0.0;
        let mut found = false;
        let mut spring = 0.5 * setup.k * extension * extension;
        for _ in 0..80 {
            let frame = owner.step(&[0.0; 2]).unwrap().clone();
            released += frame.released_energy_j;
            close(
                owner.total_energy_j().unwrap() + released,
                energy,
                energy,
                3e-8,
            );
            close(
                momentum(&owner),
                (MA + MB) * setup.common_velocity,
                (MA + MB) * v,
                5e-11,
            );
            close(
                frame.normal_work_on_network_j,
                -frame.path.normal_work_j,
                f * d,
                2e-8,
            );
            ports(&owner, frame.spring_energy_j - spring, f * d);
            spring = frame.spring_energy_j;
            if frame.path.events().any(|event| event.jump.kind == expected) {
                found = true;
                break;
            }
        }
        assert!(found, "original spring did not produce {expected:?}");
        // Under the same zero external force, the free detached initial state
        // stays in its gap: its new contact event is supplied by native spring work.
        if !attached {
            let ids = law.normal_state().elastic().interface().surface_materials();
            let mut free = JkrNormalEventModalSystem::new(
                network(Setup { links: 0, ..setup }),
                &law,
                [&ids[0], &ids[1]],
                state,
                vec![1.0 / det::sqrt(MA), -1.0 / det::sqrt(MB)],
                config(&law, 20.0),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
            for _ in 0..owner.samples_rendered() {
                assert_eq!(free.step(&[0.0; 2]).unwrap().path.event_count(), 0);
            }
            assert_eq!(free.state().contact_state().branch(), JkrBranch::Detached);
        }
    }
}

#[test]
fn unsupported_links_source_and_native_force_caps_refuse_without_publishing() {
    let law = source();
    let (d, f, _) = scales(&law);
    let setup = Setup::new(0.5 * f / d, 2.0 * d, -4.0 * d);
    let state = law.detached_state(-d).unwrap();
    let ids = law.normal_state().elastic().interface().surface_materials();
    for bad in [
        Setup { links: 0, ..setup },
        Setup { links: 2, ..setup },
        Setup {
            damping: 1e-4,
            ..setup
        },
        Setup { k: 0.0, ..setup },
    ] {
        assert!(
            JkrNormalEventModalSystem::new_with_elastic_connection(
                network(bad),
                &law,
                [&ids[0], &ids[1]],
                state,
                config(&law, 20.0),
                &CancelGate::new_clock_free()
            )
            .is_err()
        );
    }
    assert!(
        JkrNormalEventModalSystem::new_with_elastic_connection(
            network(setup),
            &law,
            [&ids[1], &ids[0]],
            state,
            config(&law, 20.0),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut owner = make(
        &law,
        law.detached_state(-20.0 * d).unwrap(),
        Setup {
            extension: 0.0,
            velocity: -0.02,
            force_limit: 0.1 * f,
            ..setup
        },
        20.0,
    );
    let old = states(&owner);
    let token = owner.state();
    assert!(owner.step(&[0.0; 2]).is_err());
    assert_eq!(states(&owner), old);
    assert_eq!(owner.state(), token);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
    let mut owner = make(&law, state, setup, 20.0);
    let old = states(&owner);
    let token = owner.state();
    let rejected: Result<(), JkrNormalEventModalError> =
        owner.step_with_admission(&[0.0; 2], &CancelGate::new_clock_free(), |_, _| {
            Err(JkrNormalEventModalError::Compatibility {
                what: "receiving owner refused",
                residual: 1.0,
                tolerance: 0.0,
            })
        });
    assert!(rejected.is_err());
    assert_eq!(states(&owner), old);
    assert_eq!(owner.state(), token);
    assert_eq!(owner.samples_rendered(), 0);
    let frame = owner.step(&[0.0; 2]).unwrap().clone();
    let mut replay = make(&law, state, setup, 20.0);
    assert_eq!(frame, replay.step(&[0.0; 2]).unwrap().clone());
}

#[path = "adhesive_jkr_elastic_modal/viscous.rs"]
mod viscous;
