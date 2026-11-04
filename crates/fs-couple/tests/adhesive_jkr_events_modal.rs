//! Original JKR branch events, native modal bodies and finite thermal receivers.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::motion::{
    JkrNormalMotion, JkrNormalMotionBudget, MAX_JKR_NORMAL_EVENTS,
};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrEnvelope, JkrJumpKind, SourceJkrSphereContact,
};
use fs_couple::friction_heat::adhesive::snap::{
    JkrSnapHeat, JkrSnapHeatControls, JkrSnapHeatError,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::adhesive_normal_events::*;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalCouplingConfig, ModalCouplingError, ModalPreparation,
    PreparedModalBudget,
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

fn source(work: f64) -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, work).normal_state().clone(),
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
fn scales(source: &SourceJkrSphereContact) -> (f64, f64, f64) {
    let h = source.indentation_scale_m();
    let e = PI * 0.05 * 0.01 * h / 6.0;
    (h, e, det::sqrt(e / MU))
}
fn column() -> Vec<f64> {
    vec![1.0 / det::sqrt(MA), 0.0, -1.0 / det::sqrt(MB), 0.0]
}
fn config(source: &SourceJkrSphereContact) -> JkrNormalEventModalConfig {
    let (h, e, v) = scales(source);
    JkrNormalEventModalConfig {
        motion: JkrNormalMotionBudget {
            maximum_step_s: DT,
            maximum_substep_s: DT / 20.0,
            maximum_substeps: 2048,
            maximum_root_iterations: 64,
            maximum_events: MAX_JKR_NORMAL_EVENTS,
            minimum_indentation_m: -500.0 * h,
            maximum_accumulated_travel_m: 200.0 * h,
            maximum_normal_speed_m_s: 100.0 * v,
            maximum_force_n: 1.0,
            coordinate_absolute_tolerance_m: 1e-10 * h,
            coordinate_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-10 * e,
            energy_relative_tolerance: 1e-10,
            impulse_absolute_tolerance_n_s: 1e-10 * MU * v,
            impulse_relative_tolerance: 1e-10,
        },
        coordinate_absolute_tolerance_m: 2e-10 * h,
        coordinate_relative_tolerance: 2e-10,
        velocity_absolute_tolerance_m_s: 2e-10 * v,
        velocity_relative_tolerance: 2e-10,
    }
}
fn network(
    vrel: f64,
    common_v: f64,
    displacement_limit: f64,
    harmonic_only: bool,
) -> CoupledModalSystem {
    let models = [MA, MB]
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
                        pressure_per_modal_velocity: C64::ZERO
                    };
                    2
                ],
                budget,
            )
            .unwrap();
            let physical_v = common_v
                + if i == 0 {
                    MB / (MA + MB) * vrel
                } else {
                    -MA / (MA + MB) * vrel
                };
            model
                .restore_states(&[
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: det::sqrt(mass) * physical_v,
                    },
                    ModalAcousticState {
                        displacement_m_sqrt_kg: 0.0,
                        velocity_m_sqrt_kg_per_s: det::sqrt(mass) * (i as f64 + 1.0) * 0.002,
                    },
                ])
                .unwrap();
            model
        })
        .collect();
    let preparation = if harmonic_only {
        ModalPreparation::HarmonicOnly(PreparedModalBudget {
            max_modes: 4,
            max_dense_bytes: 1 << 20,
            max_multiply_terms: 1 << 30,
        })
    } else {
        ModalPreparation::Ordinary
    };
    CoupledModalSystem::new_with_preparation(
        models,
        vec![],
        &[],
        ModalCouplingConfig {
            max_modes: 4,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-4,
            maximum_abs_pressure_pa: 1.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-20,
            energy_relative_tolerance: 5e-10,
        },
        preparation,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn make<'a>(
    source: &'a SourceJkrSphereContact,
    factor: f64,
    common_v: f64,
) -> JkrNormalEventModalSystem<'a> {
    let (h, _, v) = scales(source);
    let materials = source
        .normal_state()
        .elastic()
        .interface()
        .surface_materials();
    JkrNormalEventModalSystem::new(
        network(factor * v, common_v, 1.0, false),
        source,
        [&materials[0], &materials[1]],
        source.detached_state(-0.25 * h).unwrap(),
        column(),
        config(source),
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
fn momentum(owner: &JkrNormalEventModalSystem<'_>) -> f64 {
    let s = states(owner);
    det::sqrt(MA) * s[0].velocity_m_sqrt_kg_per_s + det::sqrt(MB) * s[2].velocity_m_sqrt_kg_per_s
}
fn center(owner: &JkrNormalEventModalSystem<'_>) -> f64 {
    let s = states(owner);
    (det::sqrt(MA) * s[0].displacement_m_sqrt_kg + det::sqrt(MB) * s[2].displacement_m_sqrt_kg)
        / (MA + MB)
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}; scaled error {:e}",
        (actual - expected).abs() / scale
    );
}

#[test]
fn incoming_energy_selects_capture_or_rebound_on_the_original_native_bodies() {
    let source = source(0.05);
    let (h, energy_scale, vscale) = scales(&source);
    let attach_release = -source.evaluate_attached(0.0).unwrap().reversible_energy_j;
    let detach_release = source
        .evaluate_attached(source.minimum_attached_indentation_m())
        .unwrap()
        .reversible_energy_j;
    for factor in [5.0, 10.0] {
        let vin = factor * vscale;
        let escapes = 0.5 * MU * vin * vin > attach_release + detach_release;
        assert_eq!(escapes, factor == 10.0);
        let common_v = 0.003;
        let mut owner = make(&source, factor, common_v);
        let initial_energy = owner.total_energy_j().unwrap();
        let initial_p = momentum(&owner);
        let mut released = 0.0;
        let mut jumps = Vec::new();
        let mut reversed = false;
        for _ in 0..60 {
            let frame = owner.step(&[0.0; 4]).unwrap().clone();
            released += frame.released_energy_j;
            jumps.extend(frame.path.events().copied());
            reversed |= owner.state().velocity_m_s() < 0.0;
            close(
                owner.total_energy_j().unwrap() + released,
                initial_energy,
                energy_scale,
                3e-7,
            );
            close(momentum(&owner), initial_p, (MA + MB) * vscale, 2e-10);
            let time = owner.samples_rendered() as f64 * DT;
            close(
                center(&owner),
                common_v * time,
                h.max(common_v * time),
                2e-10,
            );
            for (side, mass) in [MA, MB].into_iter().enumerate() {
                let untouched = owner.components()[side].states()[1];
                let speed = det::sqrt(mass) * (side as f64 + 1.0) * 0.002;
                close(untouched.velocity_m_sqrt_kg_per_s, speed, speed, 1e-14);
                close(
                    untouched.displacement_m_sqrt_kg,
                    speed * time,
                    speed * time,
                    1e-13,
                );
            }
            assert_eq!(frame.path.state().steps(), owner.samples_rendered());
            assert_eq!(
                owner
                    .network()
                    .last_frame()
                    .unwrap()
                    .numerical_dissipation_j,
                0.0
            );
            if jumps.len() == 2 {
                break;
            }
        }
        assert!(reversed);
        assert_eq!(jumps[0].jump.kind, JkrJumpKind::Attach);
        close(jumps[0].velocity_m_s, vin, vscale, 2e-12);
        close(
            jumps[0].jump.released_energy_j,
            attach_release,
            energy_scale,
            1e-14,
        );
        if escapes {
            assert_eq!(jumps.len(), 2);
            assert_eq!(jumps[1].jump.kind, JkrJumpKind::Detach);
            assert_eq!(jumps[1].before().branch(), JkrBranch::Attached);
            assert_eq!(jumps[1].after().branch(), JkrBranch::Detached);
            assert_eq!(
                jumps[1].before().indentation_m(),
                jumps[1].after().indentation_m()
            );
            let vout = -det::sqrt(vin * vin - 2.0 * (attach_release + detach_release) / MU);
            close(owner.state().velocity_m_s(), vout, vscale, 2e-7);
            close(jumps[1].velocity_m_s, vout, vscale, 2e-7);
            assert_eq!(owner.state().contact_state().branch(), JkrBranch::Detached);
        } else {
            assert_eq!(jumps.len(), 1);
            assert_eq!(owner.state().contact_state().branch(), JkrBranch::Attached);
        }
    }
}

#[test]
fn common_acceleration_and_uncontacted_loads_keep_their_original_work_and_motion() {
    let source = source(0.05);
    let (h, e, _) = scales(&source);
    let common_v = 0.003;
    let mut owner = make(&source, 10.0, common_v);
    let initial_energy = owner.total_energy_j().unwrap();
    let initial_p = momentum(&owner);
    let acceleration = [0.02, 0.03, 0.02, -0.01];
    let external = [
        det::sqrt(MA) * acceleration[0],
        det::sqrt(MA) * acceleration[1],
        det::sqrt(MB) * acceleration[2],
        det::sqrt(MB) * acceleration[3],
    ];
    let mut release = 0.0;
    let mut work = 0.0;
    for _ in 0..8 {
        let frame = owner.step(&external).unwrap().clone();
        release += frame.released_energy_j;
        work += frame.external_work_j;
        let t = owner.samples_rendered() as f64 * DT;
        close(
            momentum(&owner),
            initial_p + (MA + MB) * acceleration[0] * t,
            initial_p.abs(),
            2e-10,
        );
        close(
            center(&owner),
            common_v * t + 0.5 * acceleration[0] * t * t,
            h,
            2e-9,
        );
        let state = states(&owner);
        for (side, mass) in [MA, MB].into_iter().enumerate() {
            let index = 2 * side + 1;
            let v0 = (side as f64 + 1.0) * 0.002;
            close(
                state[index].displacement_m_sqrt_kg,
                det::sqrt(mass) * (v0 * t + 0.5 * acceleration[index] * t * t),
                det::sqrt(mass) * h,
                1e-10,
            );
            close(
                state[index].velocity_m_sqrt_kg_per_s,
                det::sqrt(mass) * (v0 + acceleration[index] * t),
                det::sqrt(mass) * v0,
                1e-12,
            );
        }
        let independent_work: f64 = external
            .iter()
            .zip(&state)
            .map(|(force, state)| force * state.displacement_m_sqrt_kg)
            .sum();
        close(work, independent_work, e, 2e-10);
        close(
            owner.total_energy_j().unwrap() + release - work,
            initial_energy,
            e,
            3e-7,
        );
    }
}

#[derive(Debug)]
enum AdmissionError {
    Native(JkrNormalEventModalError),
    Heat(JkrSnapHeatError),
    Receiver,
}
impl From<JkrNormalEventModalError> for AdmissionError {
    fn from(value: JkrNormalEventModalError) -> Self {
        Self::Native(value)
    }
}
impl From<JkrSnapHeatError> for AdmissionError {
    fn from(value: JkrSnapHeatError) -> Self {
        Self::Heat(value)
    }
}

#[test]
fn receiving_owner_refusal_and_late_cancellation_preserve_bodies_branch_and_clock() {
    let source = source(0.05);
    let mut owner = make(&source, 10.0, 0.003);
    let initial = states(&owner);
    let initial_contact = owner.state();
    let refused: Result<(), AdmissionError> = owner.step_with_admission(
        &[0.0; 4],
        &CancelGate::new_clock_free(),
        |frame, candidate| {
            assert_eq!(frame.path.event_count(), 1);
            assert_ne!(candidate[0].states()[0], initial[0]);
            Err(AdmissionError::Receiver)
        },
    );
    assert!(matches!(refused, Err(AdmissionError::Receiver)));
    assert_eq!(states(&owner), initial);
    assert_eq!(owner.state(), initial_contact);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
    assert!(owner.network().last_frame().is_none());
    let gate = CancelGate::new_clock_free();
    let cancelled: Result<(), AdmissionError> =
        owner.step_with_admission(&[0.0; 4], &gate, |_, _| {
            gate.request();
            Ok(())
        });
    assert!(matches!(
        cancelled,
        Err(AdmissionError::Native(JkrNormalEventModalError::Network(
            ModalCouplingError::Cancelled
        )))
    ));
    assert_eq!(states(&owner), initial);
    assert_eq!(owner.state(), initial_contact);
    assert_eq!(owner.samples_rendered(), 0);
    let accepted = owner.step(&[0.0; 4]).unwrap().clone();
    let mut reference = make(&source, 10.0, 0.003);
    assert_eq!(accepted, reference.step(&[0.0; 4]).unwrap().clone());
    assert_eq!(states(&owner), states(&reference));
}

#[test]
fn original_receiver_budgets_source_identity_and_harmonic_scope_remain_binding() {
    let source = source(0.05);
    let (h, _, v) = scales(&source);
    let materials = source
        .normal_state()
        .elastic()
        .interface()
        .surface_materials();
    let initial = source.detached_state(-0.25 * h).unwrap();
    let mut limited = JkrNormalEventModalSystem::new(
        network(10.0 * v, 0.003, 1e-15, false),
        &source,
        [&materials[0], &materials[1]],
        initial,
        column(),
        config(&source),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let before = states(&limited);
    let before_contact = limited.state();
    assert!(matches!(
        limited.step(&[0.0; 4]),
        Err(JkrNormalEventModalError::Network(
            ModalCouplingError::Component { .. }
        ))
    ));
    assert_eq!(states(&limited), before);
    assert_eq!(limited.state(), before_contact);
    assert_eq!(limited.samples_rendered(), 0);
    let foreign = self::source(0.06);
    assert!(
        JkrNormalEventModalSystem::new(
            network(10.0 * v, 0.003, 1.0, false),
            &source,
            [&materials[0], &materials[1]],
            foreign.detached_state(-0.25 * h).unwrap(),
            column(),
            config(&source),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    assert!(matches!(
        JkrNormalEventModalSystem::new(
            network(10.0 * v, 0.003, 1.0, true),
            &source,
            [&materials[0], &materials[1]],
            initial,
            column(),
            config(&source),
            &CancelGate::new_clock_free()
        ),
        Err(JkrNormalEventModalError::Network(
            ModalCouplingError::Invalid { .. }
        ))
    ));
}

#[test]
fn detached_touch_is_retained_until_the_original_zero_time_attachment_event() {
    let source = source(0.05);
    let (_, _, vscale) = scales(&source);
    let materials = source
        .normal_state()
        .elastic()
        .interface()
        .surface_materials();
    let incoming = 10.0 * vscale;
    let mut owner = JkrNormalEventModalSystem::new(
        network(incoming, 0.0, 1.0, false),
        &source,
        [&materials[0], &materials[1]],
        source.detached_state(0.0).unwrap(),
        column(),
        config(&source),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(owner.state().contact_state().branch(), JkrBranch::Detached);
    let frame = owner.step(&[0.0; 4]).unwrap().clone();
    let event = *frame.path.events().next().unwrap();
    assert_eq!(event.offset_s, 0.0);
    assert_eq!(event.before().branch(), JkrBranch::Detached);
    assert_eq!(event.after().branch(), JkrBranch::Attached);
    assert_eq!(
        event.before().indentation_m(),
        event.after().indentation_m()
    );
    close(event.velocity_m_s, incoming, vscale, 1e-14);
}

fn caloric(source: &SourceJkrSphereContact) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let elastic = source.normal_state().elastic();
    let ids = [
        elastic.surface_a_state().resolved().card_identity(),
        elastic.surface_b_state().resolved().card_identity(),
    ];
    let energy = -source.evaluate_attached(0.0).unwrap().reversible_energy_j;
    std::array::from_fn(|side| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ids[side],
            SolidLiquidPhase::Solid,
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 293.15,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: [1000.0, 2000.0][side],
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 100.0 * energy / [MA, MB][side],
                    temperature_k: 303.15,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: [1000.0, 2000.0][side],
                },
            ],
        )
        .unwrap()
    })
}
fn thermal(curves: &[EquilibriumEnthalpyPhaseCurve; 2]) -> DepositionSimulation<'_> {
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
                    ambient_k: 293.15,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        births,
        &[(MaterialId(0), &curves[0]), (MaterialId(1), &curves[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-22,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn native_admission_stages_original_snap_heat_and_combined_energy_once() {
    let source = source(0.05);
    let curves = caloric(&source);
    let mut thermal = thermal(&curves);
    let mut owner = make(&source, 10.0, 0.003);
    let mut heat = JkrSnapHeat::new(
        &source,
        owner.state().contact_state(),
        &thermal,
        [CellId(0), CellId(1)],
        [&curves[0], &curves[1]],
        HeatPartition::new(0.4, 0.4, 0.2).unwrap(),
        Some("resolved elastic waves".into()),
    )
    .unwrap();
    let initial_energy = owner.total_energy_j().unwrap();
    // A real caloric-domain refusal occurs after native mechanics was staged.
    // The finite receiver is already at its original chart's upper endpoint.
    let capacity = -100.0 * source.evaluate_attached(0.0).unwrap().reversible_energy_j;
    let (full, _) = thermal
        .apply_heat_j(
            &std::collections::BTreeMap::from([(CellId(0), capacity), (CellId(1), capacity)]),
            1e-20,
            2e-20,
            || Ok(()),
        )
        .unwrap();
    let initial_states = states(&owner);
    let initial_contact = owner.state();
    let initial_heat = heat.stamp();
    let refused =
        owner.step_with_admission(&[0.0; 4], &CancelGate::new_clock_free(), |frame, _| {
            heat.propose_normal_motion(
                &full,
                heat.stamp(),
                &frame.path,
                JkrSnapHeatControls {
                    maximum_segments: MAX_JKR_NORMAL_EVENTS + 1,
                    cell_energy_tolerance_j: 1e-22,
                    thermal_energy_tolerance_j: 2e-22,
                    transfer_energy_tolerance_j: 4e-22,
                },
                || Ok(()),
            )
            .map_err(AdmissionError::from)
        });
    assert!(matches!(
        refused,
        Err(AdmissionError::Heat(JkrSnapHeatError::Thermal(
            DepositionError::MaterialDomain
        )))
    ));
    assert_eq!(states(&owner), initial_states);
    assert_eq!(owner.state(), initial_contact);
    assert_eq!(owner.samples_rendered(), 0);
    assert_eq!(heat.stamp(), initial_heat);
    assert_eq!(thermal.energy_j(), 0.0);
    let mut event_count = 0;
    for _ in 0..60 {
        let candidate = owner
            .step_with_admission(&[0.0; 4], &CancelGate::new_clock_free(), |frame, _| {
                let endpoint_thermal = thermal
                    .advance_to(
                        frame.path.state().time_s(),
                        StepConfig {
                            max_dt_s: DT,
                            energy_tolerance_j: 1e-22,
                        },
                        || Ok(()),
                    )
                    .map_err(JkrSnapHeatError::from)?
                    .0;
                // Public observations cannot invent heat. Original sealed
                // trajectory history still supplies every actual snap.
                let mut observations = frame.path;
                observations.released_energy_j = 1e100;
                observations.normal_work_j = 1e100;
                observations.contact.reversible_energy_j = -1e100;
                heat.propose_normal_motion(
                    &endpoint_thermal,
                    heat.stamp(),
                    &observations,
                    JkrSnapHeatControls {
                        maximum_segments: MAX_JKR_NORMAL_EVENTS + 1,
                        cell_energy_tolerance_j: 1e-22,
                        thermal_energy_tolerance_j: 2e-22,
                        transfer_energy_tolerance_j: 4e-22,
                    },
                    || Ok(()),
                )
                .map_err(AdmissionError::from)
            })
            .unwrap();
        let (next_heat, next_thermal, report) = candidate.into_parts();
        event_count += report.jumps.len();
        heat = next_heat;
        thermal = next_thermal;
        assert_eq!(heat.state(), owner.state().contact_state());
        assert_eq!(thermal.time_s(), owner.state().time_s());
        close(
            report.work.dissipated_work_j(),
            owner.last_frame().unwrap().released_energy_j,
            scales(&source).1,
            1e-14,
        );
        close(
            owner.total_energy_j().unwrap() + thermal.energy_j() + heat.ledger().other_work_j(),
            initial_energy,
            scales(&source).1,
            3e-7,
        );
        if event_count == 2 {
            break;
        }
    }
    assert_eq!(event_count, 2);
    assert!(thermal.temperature_k(CellId(0)).unwrap() > 293.15);
    assert!(thermal.temperature_k(CellId(1)).unwrap() > 293.15);
}

#[test]
fn sealed_motion_chronology_refuses_replay_even_when_contact_state_is_unchanged() {
    let source = source(0.05);
    let curves = caloric(&source);
    let thermal = thermal(&curves);
    let equilibrium = source.evaluate_attached(0.0).unwrap();
    let motion = JkrNormalMotion::new(&source, MU, config(&source).motion).unwrap();
    let initial = motion.initialize(equilibrium.state(), 0.0).unwrap();
    let gate = CancelGate::new_clock_free();
    let step = motion
        .advance(initial, equilibrium.normal_force_n, DT, &gate)
        .unwrap();
    assert_eq!(
        step.initial_state().contact_state(),
        step.state().contact_state()
    );
    assert_ne!(step.initial_state(), step.state());
    assert_eq!(step.event_count(), 0);
    let heat = JkrSnapHeat::new(
        &source,
        equilibrium.state(),
        &thermal,
        [CellId(0), CellId(1)],
        [&curves[0], &curves[1]],
        HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
        None,
    )
    .unwrap();
    let controls = JkrSnapHeatControls {
        maximum_segments: MAX_JKR_NORMAL_EVENTS + 1,
        cell_energy_tolerance_j: 1e-22,
        thermal_energy_tolerance_j: 2e-22,
        transfer_energy_tolerance_j: 4e-22,
    };
    let (accepted, thermal, _) = heat
        .propose_normal_motion(&thermal, heat.stamp(), &step, controls, || Ok(()))
        .unwrap()
        .into_parts();
    assert!(matches!(
        accepted.propose_normal_motion(&thermal, accepted.stamp(), &step, controls, || Ok(())),
        Err(JkrSnapHeatError::HistoryMismatch)
    ));
    // The raw endpoint API cannot erase the stronger chronology once bound.
    assert!(matches!(
        accepted.propose(
            &thermal,
            accepted.stamp(),
            &[0.0],
            equilibrium.state(),
            controls,
            || Ok(())
        ),
        Err(JkrSnapHeatError::Invalid(_))
    ));
    let next = motion
        .advance(step.state(), equilibrium.normal_force_n, DT, &gate)
        .unwrap();
    assert!(
        accepted
            .propose_normal_motion(&thermal, accepted.stamp(), &next, controls, || Ok(()))
            .is_ok()
    );
    // Same geometric history and elapsed time from a different relative mass
    // still belongs to another physical owner and cannot continue this ledger.
    let other = JkrNormalMotion::new(&source, 2.0 * MU, config(&source).motion).unwrap();
    let other_initial = other.initialize(equilibrium.state(), 0.0).unwrap();
    let other_first = other
        .advance(other_initial, equilibrium.normal_force_n, DT, &gate)
        .unwrap();
    let other_next = other
        .advance(other_first.state(), equilibrium.normal_force_n, DT, &gate)
        .unwrap();
    assert!(matches!(
        accepted.propose_normal_motion(
            &thermal,
            accepted.stamp(),
            &other_next,
            controls,
            || Ok(())
        ),
        Err(JkrSnapHeatError::HistoryMismatch)
    ));
}
