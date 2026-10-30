//! Native source-backed gases, original impulse, complete DG and rational memory.
//! Synthetic caloric/transport laws test the implementation, not measured gases.

use super::*;
use crate::cavity_phs::piston::GasPiston;
use crate::cavity_phs::piston::tests::{controls as single_controls, drive, gas};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_duct::{Duct, Segment};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_math::c64::C64;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::model::PoleTerm;

pub(super) fn mechanics() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: single_controls(),
        opposing_gas: single_controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}
pub(super) fn motion() -> AcousticMotionSolve {
    AcousticMotionSolve {
        maximum_faces: 2,
        pressure_tolerance_pa: 1e-7,
        velocity_tolerance_m_s: 1e-10,
        energy_tolerance_j: 1e-7,
    }
}
pub(super) fn forcing() -> DoublePistonInput {
    DoublePistonInput {
        duration_s: drive(0.0).duration_s,
        heat_j: [0.2, 0.1],
        external_force_n: 0.4,
    }
}
pub(super) fn pair(velocity: f64) -> DoubleActingGasPiston {
    DoubleActingGasPiston::new(
        GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, velocity).unwrap(),
        gas(0.25, 0.001),
        0.01,
    )
    .unwrap()
}
pub(super) fn load(fraction: f64, dt: f64) -> DgPortLoad {
    let parcel = gas(fraction, 0.001);
    let field = AxialGasField::new(
        vec![0.0, 0.05],
        parcel.state().state().pressure,
        vec![
            PrimitiveGasState {
                temperature_k: parcel.state().state().temperature,
                composition: parcel.model().thermodynamics().mole_composition().clone(),
            };
            2
        ],
        parcel.model().clone(),
        QueryPoint::new().with("specimen-lot", 42.0).unwrap(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    assert_eq!(field.uniform_state(), Some(parcel.state()));
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.05,
            }],
        },
        &field,
        DgMesh {
            elements: 2,
            order: 1,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap();
    let terminal = PreparedMimoLoad::new(
        &MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -400.0,
                residue: 0.0,
            }],
            residues: vec![vec![C64::from_re(4e10)]],
            d: vec![1e7],
            e: vec![0.0],
            reciprocal: true,
        },
        dt,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / dt],
            sampled_passivity_points: 129,
            max_states: 4,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    DgPortLoad::new_with_terminal(duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal).unwrap()
}
pub(super) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}
pub(super) fn assert_piston_same(actual: &DoubleActingGasPiston, expected: &DoubleActingGasPiston) {
    assert_eq!(actual, expected);
    assert_eq!(actual.ports, expected.ports);
    assert_eq!(actual.primary().ports, expected.primary().ports);
}
pub(super) fn assert_load_same(actual: &DgPortLoad, expected: &DgPortLoad) {
    assert_eq!(actual.state(), expected.state());
    assert_eq!(actual.duct().gas_field(), expected.duct().gas_field());
    assert_eq!(
        actual.port_response().unwrap(),
        expected.port_response().unwrap()
    );
    assert_eq!(
        actual.terminal().unwrap().pole_states(),
        expected.terminal().unwrap().pole_states()
    );
    assert_eq!(
        actual.last_terminal_frames(),
        expected.last_terminal_frames()
    );
    assert_eq!(actual.stored_energy_j(), expected.stored_energy_j());
    assert_eq!(
        actual.cumulative_dissipated_j(),
        expected.cumulative_dissipated_j()
    );
}
// Integrate the fixture's two Cp/R polynomials independently of native reports.
pub(super) fn caloric_energy(parcel: &ClosedGasParcel) -> f64 {
    let x = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let t = parcel.state().state().temperature;
    parcel.mass_kg() * R / (x[0] * 0.028 + x[1] * 0.004)
        * (x[0] * (2.5 * t + 0.0005 * t * t) + x[1] * 1.5 * t)
}
fn faces(loads: &[DgPortLoad], areas: [f64; 2]) -> [AcousticMotionFace<'_>; 2] {
    [
        AcousticMotionFace {
            load: &loads[0],
            signed_area_m2: areas[0],
        },
        AcousticMotionFace {
            load: &loads[1],
            signed_area_m2: areas[1],
        },
    ]
}

#[test]
fn g1_independent_media_share_actual_wall_work_and_retain_every_original_history() {
    let dt = forcing().duration_s;
    let areas = [0.001, -0.0007];
    let mut live = pair(0.03);
    let mut loads = vec![load(0.15, dt), load(0.85, dt)];
    assert_ne!(
        loads[0].duct().inlet_gas().density,
        loads[1].duct().inlet_gas().density
    );
    assert_ne!(
        loads[0].duct().inlet_gas().sound_speed,
        loads[1].duct().inlet_gas().sound_speed
    );
    let fields = [
        loads[0].duct().gas_field().cloned(),
        loads[1].duct().gas_field().cloned(),
    ];
    let initial = live.clone();
    let mut external = 0.0;
    let mut exported = 0.0;
    for sample in 0..12 {
        let input = DoublePistonInput {
            heat_j: if sample < 3 {
                forcing().heat_j
            } else {
                [0.0; 2]
            },
            external_force_n: if sample < 4 {
                forcing().external_force_n
            } else {
                0.0
            },
            ..forcing()
        };
        let previous = live.clone();
        let previous_loads = loads.clone();
        let relations: Vec<_> = loads
            .iter()
            .map(|load| load.port_response().unwrap())
            .collect();
        let candidate = live
            .acoustic_candidate(input, mechanics(), &faces(&loads, areas), motion(), || {
                false
            })
            .unwrap();
        assert_piston_same(&live, &previous);
        for (actual, expected) in loads.iter().zip(&previous_loads) {
            assert_load_same(actual, expected);
        }
        let mechanical = candidate.piston_step();
        let exchange = candidate.motion();
        let dx = mechanical.displacement_m;
        let velocity = dx / dt;
        near(
            velocity,
            0.5 * (mechanical.velocity_m_s[0] + mechanical.velocity_m_s[1]),
            1e-10,
        );
        near(exchange.midpoint_velocity_m_s, velocity, 1e-10);
        let mut actual_force = 0.0;
        let mut actual_work = 0.0;
        for i in 0..2 {
            let q = areas[i] * velocity;
            let frame = &exchange.faces[i];
            let port = frame
                .duct
                .port_outputs
                .iter()
                .find(|port| port.id == loads[i].inlet())
                .unwrap();
            near(frame.inlet_flow_m3_s, q, 1e-17);
            near(frame.mean_pressure_pa, port.conjugate_rate, 1e-7);
            near(
                frame.mean_pressure_pa,
                relations[i].history_pressure_pa + relations[i].impedance_pa_s_m3 * q,
                1e-7,
            );
            near(frame.inlet_work_j, q * port.delta_x, 1e-13);
            near(frame.inlet_work_j, dt * q * port.conjugate_rate, 1e-13);
            near(frame.port_work_residual_j, 0.0, 1e-12);
            actual_force -= areas[i] * port.conjugate_rate;
            actual_work += q * port.delta_x;
            // Replay the ORIGINAL complete owner at the actual geometric flow.
            // This checks acoustic and rational states independently of the
            // mechanical candidate, including each native port-work ledger.
            let mut reference = previous_loads[i].clone();
            let output = reference.preview_step(q).unwrap().accept();
            assert_eq!(frame.duct, output);
            assert_load_same(&candidate.loads()[i], &reference);
            assert_eq!(candidate.loads()[i].duct().gas_field(), fields[i].as_ref());
            assert_eq!(candidate.loads()[i].state().step_count(), sample + 1);
        }
        near(exchange.force_n, actual_force, 1e-9);
        near(exchange.mechanical_work_j, actual_force * dx, 1e-12);
        near(actual_work + actual_force * dx, 0.0, 1e-12);
        near(exchange.transfer_residual_j, 0.0, 1e-12);
        let gas_force = 0.01
            * (mechanical.midpoint[0].state().pressure - mechanical.midpoint[1].state().pressure);
        near(
            live.primary().mass_kg()
                * (candidate.piston().primary().velocity_m_s() - live.primary().velocity_m_s()),
            dt * (gas_force + input.external_force_n + actual_force),
            2e-10,
        );
        if sample == 0 {
            let mut unloaded = live.clone();
            unloaded.step(input, mechanics(), || false).unwrap();
            assert!(
                candidate.piston().primary().velocity_m_s() < unloaded.primary().velocity_m_s()
            );
            assert!(actual_work > 1e-8);
        }
        near(
            exchange.acoustic_energy_residual_j,
            0.0,
            motion().energy_tolerance_j,
        );
        near(
            candidate.total_energy_residual_j(),
            0.0,
            motion().energy_tolerance_j,
        );
        external += input.heat_j.iter().sum::<f64>() + input.external_force_n * dx;
        exported += exchange.terminal_work_j;
        let (next, next_loads, _, _, _) = candidate.into_parts();
        live = next;
        loads = next_loads;
    }
    let gas_change = caloric_energy(live.primary().gas()) - caloric_energy(initial.primary().gas())
        + caloric_energy(live.opposing())
        - caloric_energy(initial.opposing());
    let kinetic_change = 0.5
        * live.primary().mass_kg()
        * (live.primary().velocity_m_s().powi(2) - initial.primary().velocity_m_s().powi(2));
    let acoustics: f64 = loads
        .iter()
        .map(|load| load.stored_energy_j() + load.cumulative_dissipated_j())
        .sum();
    near(
        gas_change + kinetic_change + acoustics + exported,
        external,
        12.0 * motion().energy_tolerance_j,
    );
    assert!(loads.iter().all(|load| {
        load.terminal()
            .unwrap()
            .pole_states()
            .iter()
            .any(|z| *z != C64::ZERO)
    }));
}

#[test]
fn g3_zero_area_reduces_exactly_to_native_motion_and_continues_live_acoustic_memory() {
    let dt = forcing().duration_s;
    let original = pair(0.03);
    let cold = vec![load(0.15, dt), load(0.85, dt)];
    let warm = original
        .acoustic_candidate(
            forcing(),
            mechanics(),
            &faces(&cold, [0.001, -0.0007]),
            motion(),
            || false,
        )
        .unwrap();
    let (live, loads, _, _, _) = warm.into_parts();
    let before = live.clone();
    let mut native = live.clone();
    let native_step = native.step(forcing(), mechanics(), || false).unwrap();
    let candidate = live
        .acoustic_candidate(
            forcing(),
            mechanics(),
            &faces(&loads, [0.0; 2]),
            motion(),
            || false,
        )
        .unwrap();
    assert_piston_same(candidate.piston(), &native);
    assert_eq!(candidate.piston_step(), &native_step);
    assert_piston_same(&live, &before);
    assert_eq!(candidate.motion().force_n, 0.0);
    assert_eq!(candidate.motion().mechanical_work_j, 0.0);
    for i in 0..2 {
        let mut free = loads[i].clone();
        let output = free.preview_step(0.0).unwrap().accept();
        assert_load_same(&candidate.loads()[i], &free);
        assert_eq!(candidate.motion().faces[i].duct, output);
        assert_ne!(candidate.loads()[i].state(), loads[i].state());
        assert_ne!(
            candidate.loads()[i].terminal().unwrap().pole_states(),
            loads[i].terminal().unwrap().pole_states()
        );
        assert_eq!(candidate.motion().faces[i].inlet_work_j, 0.0);
    }
    near(
        candidate.motion().acoustic_energy_residual_j,
        0.0,
        motion().energy_tolerance_j,
    );
}

#[test]
fn g4_refusal_and_final_cancellation_publish_nothing_and_retry_preserves_the_suffix() {
    let dt = forcing().duration_s;
    let live = pair(0.03);
    let loads = vec![load(0.15, dt), load(0.85, dt)];
    let original = live.clone();
    let original_loads = loads.clone();
    let areas = [0.001, -0.0007];
    let wrong_clock = load(0.15, dt * 2.0);
    assert!(
        live.acoustic_candidate(
            forcing(),
            mechanics(),
            &[
                AcousticMotionFace {
                    load: &loads[0],
                    signed_area_m2: areas[0]
                },
                AcousticMotionFace {
                    load: &wrong_clock,
                    signed_area_m2: areas[1]
                },
            ],
            motion(),
            || false
        )
        .is_err()
    );
    assert!(
        live.acoustic_candidate(
            forcing(),
            mechanics(),
            &[
                AcousticMotionFace {
                    load: &loads[0],
                    signed_area_m2: areas[0]
                },
                AcousticMotionFace {
                    load: &loads[0],
                    signed_area_m2: areas[1]
                },
            ],
            motion(),
            || false
        )
        .is_err()
    );
    for bad_areas in [[f64::NAN, areas[1]], [areas[0], f64::INFINITY]] {
        assert!(
            live.acoustic_candidate(
                forcing(),
                mechanics(),
                &faces(&loads, bad_areas),
                motion(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        live.acoustic_candidate(
            forcing(),
            mechanics(),
            &faces(&loads, areas),
            AcousticMotionSolve {
                maximum_faces: 1,
                ..motion()
            },
            || false
        )
        .is_err()
    );
    let mut narrow = mechanics();
    narrow.piston.volume_bracket_m3 = [0.001 - 1e-12, 0.001 + 1e-12];
    assert!(
        live.acoustic_candidate(forcing(), narrow, &faces(&loads, areas), motion(), || false)
            .is_err()
    );
    let mut polls = 0usize;
    let expected = live
        .acoustic_candidate(
            forcing(),
            mechanics(),
            &faces(&loads, areas),
            motion(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 2);
    let final_poll = polls;
    polls = 0;
    let refused = live.acoustic_candidate(
        forcing(),
        mechanics(),
        &faces(&loads, areas),
        motion(),
        || {
            polls += 1;
            polls == final_poll
        },
    );
    assert!(matches!(refused, Err(AcousticMotionError::Cancelled)));
    assert_eq!(polls, final_poll);
    assert_piston_same(&live, &original);
    for (actual, before) in loads.iter().zip(&original_loads) {
        assert_load_same(actual, before);
    }
    let retried = live
        .acoustic_candidate(
            forcing(),
            mechanics(),
            &faces(&loads, areas),
            motion(),
            || false,
        )
        .unwrap();
    assert_piston_same(retried.piston(), expected.piston());
    assert_eq!(retried.piston_step(), expected.piston_step());
    assert_eq!(retried.motion().force_n, expected.motion().force_n);
    assert_eq!(
        retried.motion().transfer_residual_j,
        expected.motion().transfer_residual_j
    );
    for (actual, reference) in retried.loads().iter().zip(expected.loads()) {
        assert_load_same(actual, reference);
    }
    let (a, loads_a, _, _, _) = retried.into_parts();
    let (b, loads_b, _, _, _) = expected.into_parts();
    let input = DoublePistonInput {
        heat_j: [0.0; 2],
        external_force_n: -0.3,
        ..forcing()
    };
    let next_a = a
        .acoustic_candidate(
            input,
            mechanics(),
            &faces(&loads_a, areas),
            motion(),
            || false,
        )
        .unwrap();
    let next_b = b
        .acoustic_candidate(
            input,
            mechanics(),
            &faces(&loads_b, areas),
            motion(),
            || false,
        )
        .unwrap();
    assert_piston_same(next_a.piston(), next_b.piston());
    assert_eq!(next_a.piston_step(), next_b.piston_step());
    for (actual, reference) in next_a.loads().iter().zip(next_b.loads()) {
        assert_load_same(actual, reference);
    }
}
