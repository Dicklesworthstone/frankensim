//! Native contact/spring motion, event-timed caloric cooling and joint rollback.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

#[path = "adhesive_jkr_thermal_modal/damping.rs"]
mod damping;

use fs_contact::normal_patch::adhesive::{motion::*, *};
use fs_couple::friction_heat::adhesive::snap::{timed::*, viscous::*, *};
use fs_couple::modal_acoustic_time::*;
use fs_couple::render::schedule::force::coupled::contact::adhesive_normal_events::{thermal::*, *};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalPreparation,
};
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_math::{c64::C64, det};
use fs_thermal::enthalpy::*;
use fs_tribo::HeatPartition;

const MA: f64 = 2e-5;
const MB: f64 = 3e-5;
const MU: f64 = MA * MB / (MA + MB);
const DT: f64 = 1e-4;
const AMBIENT: f64 = 293.15;
const ALPHA: f64 = 1500.0;

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
    (d, f * d, det::sqrt(f * d / MU))
}
fn release(law: &SourceJkrSphereContact) -> f64 {
    -law.evaluate_attached(0.0).unwrap().reversible_energy_j
}
fn mechanical<'a>(
    law: &'a SourceJkrSphereContact,
    state: JkrState,
    incoming: f64,
    support: Option<(f64, f64)>,
) -> JkrNormalEventModalSystem<'a> {
    mechanical_with_damping(
        law,
        state,
        incoming,
        support.map(|(k, extension)| (k, extension, 0.0)),
    )
}
fn mechanical_with_damping<'a>(
    law: &'a SourceJkrSphereContact,
    state: JkrState,
    incoming: f64,
    support: Option<(f64, f64, f64)>,
) -> JkrNormalEventModalSystem<'a> {
    let (d, e, v) = scales(law);
    let models = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(side, mass)| {
            let factor = if side == 0 {
                MB / (MA + MB)
            } else {
                -MA / (MA + MB)
            };
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                10_000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model
                .restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: det::sqrt(mass) * (0.003 + factor * incoming),
                }])
                .unwrap();
            model
        })
        .collect();
    let links = support
        .into_iter()
        .map(|(k, extension, damping)| ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0 / det::sqrt(MA)],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![1.0 / det::sqrt(MB)],
            },
            stiffness_n_m: k,
            damping_n_s_m: damping,
            rest_extension_m: -extension,
        })
        .collect();
    let native = CoupledModalSystem::new_with_preparation(
        models,
        links,
        &[],
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 1,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-4,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-20,
            energy_relative_tolerance: 5e-10,
        },
        ModalPreparation::Ordinary,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let config = JkrNormalEventModalConfig {
        motion: JkrNormalMotionBudget {
            maximum_step_s: DT,
            maximum_substep_s: DT / 20.0,
            maximum_substeps: 2048,
            maximum_root_iterations: 64,
            maximum_events: MAX_JKR_NORMAL_EVENTS,
            minimum_indentation_m: -500.0 * d,
            maximum_accumulated_travel_m: 200.0 * d,
            maximum_normal_speed_m_s: 100.0 * v,
            maximum_force_n: 1.0,
            coordinate_absolute_tolerance_m: 1e-10 * d,
            coordinate_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-10 * e,
            energy_relative_tolerance: 1e-10,
            impulse_absolute_tolerance_n_s: 1e-10 * MU * v,
            impulse_relative_tolerance: 1e-10,
        },
        coordinate_absolute_tolerance_m: 2e-10 * d,
        coordinate_relative_tolerance: 2e-10,
        velocity_absolute_tolerance_m_s: 2e-10 * v,
        velocity_relative_tolerance: 2e-10,
    };
    let ids = law.normal_state().elastic().interface().surface_materials();
    if support.is_some_and(|(_, _, c)| c > 0.0) {
        JkrNormalEventModalSystem::new_with_viscoelastic_connection(
            native,
            law,
            [&ids[0], &ids[1]],
            state,
            config,
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    } else if support.is_some() {
        JkrNormalEventModalSystem::new_with_elastic_connection(
            native,
            law,
            [&ids[0], &ids[1]],
            state,
            config,
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    } else {
        JkrNormalEventModalSystem::new(
            native,
            law,
            [&ids[0], &ids[1]],
            state,
            vec![1.0 / det::sqrt(MA), -1.0 / det::sqrt(MB)],
            config,
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    }
}
fn curves(law: &SourceJkrSphereContact, widths: [f64; 2]) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let e = law.normal_state().elastic();
    let ids = [
        e.surface_a_state().resolved().card_identity(),
        e.surface_b_state().resolved().card_identity(),
    ];
    // Manufactured constant heat capacities C_i = release * widths[i], J/K.
    std::array::from_fn(|side| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ids[side],
            SolidLiquidPhase::Solid,
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: AMBIENT,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: [1000.0, 2000.0][side],
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 100.0 * release(law) * widths[side] / [MA, MB][side],
                    temperature_k: AMBIENT + 100.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: [1000.0, 2000.0][side],
                },
            ],
        )
        .unwrap()
    })
}
fn controls() -> JkrTimedSnapHeatControls {
    JkrTimedSnapHeatControls {
        heat: JkrSnapHeatControls {
            maximum_segments: MAX_JKR_NORMAL_EVENTS + 1,
            cell_energy_tolerance_j: 1e-21,
            thermal_energy_tolerance_j: 2e-21,
            transfer_energy_tolerance_j: 4e-21,
        },
        transport: StepConfig {
            max_dt_s: DT / 128.0,
            energy_tolerance_j: 1e-20,
        },
        maximum_transport_substeps: 4096,
        total_energy_tolerance_j: 5e-20,
    }
}
fn timeline<'a>(
    law: &SourceJkrSphereContact,
    caloric: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    widths: [f64; 2],
    origin: f64,
) -> DepositionSimulation<'a> {
    let births = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(side, mass)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(side as u64),
                segment_id: side as u64,
                body_id: side as u64,
                material_id: MaterialId(side as u64),
                mass_kg: mass,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [side as f64, 0.0, 0.0],
                    end_m: [side as f64 + 0.001, 0.0, 0.0],
                    volume_m3: mass / [1000.0, 2000.0][side],
                    surface_area_m2: 1e-4,
                },
                boundary: Boundary {
                    ambient_k: AMBIENT,
                    htc_w_m2_k: ALPHA * release(law) * widths[side] / 1e-4,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        births,
        &[(MaterialId(0), &caloric[0]), (MaterialId(1), &caloric[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 4096,
        },
    )
    .unwrap()
    .advance_to(origin, controls().transport, || Ok(()))
    .unwrap()
    .0
}
fn joint<'s, 'm>(
    mechanical: JkrNormalEventModalSystem<'s>,
    thermal: DepositionSimulation<'m>,
    caloric: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    origin: f64,
) -> JkrThermalModalSystem<'s, 'm> {
    JkrThermalModalSystem::new(
        mechanical,
        thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
        Some("outgoing elastic waves".into()),
        origin,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn joint_with_support_heat<'s, 'm>(
    mechanical: JkrNormalEventModalSystem<'s>,
    thermal: DepositionSimulation<'m>,
    caloric: &'s [EquilibriumEnthalpyPhaseCurve; 2],
    origin: f64,
) -> JkrThermalModalSystem<'s, 'm> {
    JkrThermalModalSystem::new_with_support_heat(
        mechanical,
        thermal,
        [CellId(0), CellId(1)],
        [&caloric[0], &caloric[1]],
        HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
        Some("outgoing elastic waves".into()),
        origin,
        controls(),
        JkrViscousHeatRoute::new(
            HeatPartition::new(0.6, 0.1, 0.3).unwrap(),
            Some("unresolved support reservoir".into()),
        )
        .unwrap(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn states(owner: &JkrNormalEventModalSystem<'_>) -> Vec<ModalAcousticState> {
    owner
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} != {b:e}; scaled error {:e}",
        (a - b).abs() / scale
    );
}
fn independent_temperature(time: f64, side: usize, pulses: &[(f64, f64)], capacity: f64) -> f64 {
    AMBIENT
        + pulses
            .iter()
            .map(|&(at, released)| {
                [0.25, 0.5][side] * released / capacity * (-ALPHA * (time - at)).exp()
            })
            .sum::<f64>()
}

#[test]
fn native_rebound_snap_temperatures_follow_event_exponentials_and_full_energy_closes() {
    let law = source();
    let (d, e, v) = scales(&law);
    let state = law.detached_state(-0.25 * d).unwrap();
    let caloric = curves(&law, [1.0; 2]);
    let mut reference = mechanical(&law, state, 10.0 * v, None);
    let mut owner = joint(
        mechanical(&law, state, 10.0 * v, None),
        timeline(&law, &caloric, [1.0; 2], 0.0),
        &caloric,
        0.0,
    );
    let initial = owner.stored_energy_j().unwrap();
    let (mut outside, mut transport, mut incoming, mut work) = (0.0, 0.0, 0.0, 0.0);
    let mut pulses = Vec::new();
    let mut kinds = Vec::new();
    let mut after_rebound = 0;
    let mut first_cooled = false;
    for _ in 0..60 {
        let frame = owner.step(&[0.0; 2]).unwrap().clone();
        let bare = reference.step(&[0.0; 2]).unwrap().clone();
        assert_eq!(frame.mechanical, bare);
        assert_eq!(states(owner.mechanical()), states(&reference));
        for event in frame.mechanical.path.events() {
            let at = frame.mechanical.path.initial_state().time_s() + event.offset_s;
            pulses.push((at, event.jump.released_energy_j));
            kinds.push(event.jump.kind);
        }
        let time = owner.thermal().time_s();
        assert_eq!(time, owner.mechanical().state().time_s());
        assert_eq!(
            owner.heat_history().state(),
            owner.mechanical().state().contact_state()
        );
        for side in 0..2 {
            close(
                owner.thermal().temperature_k(CellId(side as u64)).unwrap(),
                independent_temperature(time, side, &pulses, release(&law)),
                1.0,
                1e-3,
            );
        }
        if owner.mechanical().samples_rendered() == 1 {
            assert_eq!(pulses.len(), 1);
            first_cooled = AMBIENT + 0.5 - owner.thermal().temperature_k(CellId(1)).unwrap() > 0.01;
        }
        outside += frame.heat.snap.work.other_work_j();
        transport += frame.heat.transport.external_heat_j;
        incoming += frame.heat.transport.incoming_energy_j;
        work += frame.mechanical.external_work_j;
        close(
            owner.stored_energy_j().unwrap() + outside - transport - incoming - work,
            initial,
            e,
            3e-7,
        );
        assert!(frame.combined_energy_residual_j.abs() <= frame.combined_energy_tolerance_j);
        if kinds.len() == 2 {
            after_rebound += 1;
        }
        if after_rebound == 5 {
            break;
        }
    }
    assert_eq!(kinds, vec![JkrJumpKind::Attach, JkrJumpKind::Detach]);
    assert!(
        first_cooled,
        "first snap must cool before the first interval ends"
    );
    assert!(transport < 0.0);
    assert_eq!(
        owner.mechanical().state().contact_state().branch(),
        JkrBranch::Detached
    );
}

#[test]
fn original_spring_capture_and_pull_off_feed_timed_heat_on_an_explicit_clock() {
    let law = source();
    let (d, e, _) = scales(&law);
    let caloric = curves(&law, [1.0; 2]);
    let origin = 3.0 * DT;
    for (attached, extension, expected) in [
        (false, -4.0 * d, JkrJumpKind::Attach),
        (true, 30.0 * d, JkrJumpKind::Detach),
    ] {
        let initial = if attached {
            law.evaluate_attached(0.0).unwrap().state()
        } else {
            law.detached_state(-d).unwrap()
        };
        let mut owner = joint(
            mechanical(&law, initial, 0.0, Some((0.5 * e / (d * d), extension))),
            timeline(&law, &caloric, [1.0; 2], origin),
            &caloric,
            origin,
        );
        let energy = owner.stored_energy_j().unwrap();
        let (mut outside, mut transport) = (0.0, 0.0);
        let mut event_seen = false;
        let mut spring_work = 0.0;
        let mut pulses = Vec::new();
        for _ in 0..80 {
            let frame = owner.step(&[0.0; 2]).unwrap().clone();
            outside += frame.heat.snap.work.other_work_j();
            transport += frame.heat.transport.external_heat_j;
            spring_work += frame.mechanical.path.support_work_j.abs();
            assert_eq!(frame.mechanical.external_work_j, 0.0);
            assert!(frame.mechanical.spring_energy_j > 0.0);
            assert_eq!(
                owner.thermal().time_s(),
                origin + owner.mechanical().state().time_s()
            );
            close(
                owner.stored_energy_j().unwrap() + outside - transport,
                energy,
                e,
                5e-7,
            );
            for event in frame.mechanical.path.events() {
                let at = origin + (frame.mechanical.path.initial_state().time_s() + event.offset_s);
                pulses.push((at, event.jump.released_energy_j));
                event_seen |= event.jump.kind == expected;
            }
            for side in 0..2 {
                close(
                    owner.thermal().temperature_k(CellId(side as u64)).unwrap(),
                    independent_temperature(owner.thermal().time_s(), side, &pulses, release(&law)),
                    1.0,
                    1e-3,
                );
            }
            if event_seen {
                break;
            }
        }
        assert!(event_seen);
        assert!(spring_work > 0.0);
        assert!(owner.heat_history().ledger().dissipated_work_j() > 0.0);
        assert!(owner.thermal().energy_j() > 0.0);
        assert_eq!(
            owner.heat_history().state(),
            owner.mechanical().state().contact_state()
        );
    }
}

#[test]
fn caloric_failure_cancellation_retry_and_clock_refusal_preserve_all_native_owners() {
    let law = source();
    let (d, _, v) = scales(&law);
    let state = law.detached_state(-0.25 * d).unwrap();
    let narrow = curves(&law, [1.0, 1e-4]);
    let mut failed = joint(
        mechanical(&law, state, 10.0 * v, None),
        timeline(&law, &narrow, [1.0, 1e-4], 0.0),
        &narrow,
        0.0,
    );
    let before = states(failed.mechanical());
    let motion = failed.mechanical().state();
    let stamp = failed.heat_history().stamp();
    let cells = failed.thermal().cells().clone();
    assert!(matches!(
        failed.step(&[0.0; 2]),
        Err(JkrThermalModalError::Thermal(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    ));
    assert_eq!(states(failed.mechanical()), before);
    assert_eq!(failed.mechanical().state(), motion);
    assert_eq!(failed.heat_history().stamp(), stamp);
    assert_eq!(failed.thermal().cells(), &cells);
    assert_eq!(failed.thermal().time_s(), 0.0);
    assert_eq!(failed.mechanical().samples_rendered(), 0);
    assert!(failed.last_frame().is_none());
    assert!(failed.mechanical().last_frame().is_none());
    assert!(failed.mechanical().network().last_frame().is_none());
    let caloric = curves(&law, [1.0; 2]);
    let mut owner = joint(
        mechanical(&law, state, 10.0 * v, None),
        timeline(&law, &caloric, [1.0; 2], 0.0),
        &caloric,
        0.0,
    );
    let stamp = owner.heat_history().stamp();
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(owner.step_under_gate(&[0.0; 2], &gate).is_err());
    assert_eq!(owner.heat_history().stamp(), stamp);
    assert_eq!(owner.mechanical().samples_rendered(), 0);
    assert_eq!(owner.thermal().time_s(), 0.0);
    assert!(owner.last_frame().is_none());
    let accepted = owner.step(&[0.0; 2]).unwrap().clone();
    let mut fresh = joint(
        mechanical(&law, state, 10.0 * v, None),
        timeline(&law, &caloric, [1.0; 2], 0.0),
        &caloric,
        0.0,
    );
    assert_eq!(accepted, *fresh.step(&[0.0; 2]).unwrap());
    assert_eq!(owner.thermal().cells(), fresh.thermal().cells());
    assert_eq!(states(owner.mechanical()), states(fresh.mechanical()));
    assert!(
        JkrThermalModalSystem::new(
            mechanical(&law, state, 10.0 * v, None),
            timeline(&law, &caloric, [1.0; 2], 3.0 * DT),
            [CellId(0), CellId(1)],
            [&caloric[0], &caloric[1]],
            HeatPartition::new(0.25, 0.5, 0.25).unwrap(),
            Some("outgoing elastic waves".into()),
            0.0,
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
}
