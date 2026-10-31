//! Original PR storage and DG histories share one reciprocal motion account.

use super::{assert_load_same, faces, load, motion, near};
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput,
    acoustic::{AcousticMotionError, AcousticMotionSolve},
    real_gas::{real_controls, real_pair},
};
use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
    UniformEnthalpyStepInput,
};
use fs_material::thermochemical_gas::nonideal::{
    ResolvedThermochemicalRealGasState, energy::parcel::IsochoricRealGasParcel,
};
use std::convert::Infallible;

const AREAS: [f64; 2] = [5e-5, -3e-5];

fn input() -> DoublePistonInput {
    DoublePistonInput {
        duration_s: 2e-5,
        heat_j: [0.3, 0.1],
        external_force_n: 30.0,
    }
}
fn acoustic() -> AcousticMotionSolve {
    AcousticMotionSolve {
        energy_tolerance_j: 2e-6,
        ..motion()
    }
}
fn exchange() -> DoubleHeatExchangeSolve {
    DoubleHeatExchangeSolve {
        faces: [HeatExchangeSolve {
            heat_bracket_j: [-0.25, 0.25],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 2e-6,
            maximum_iterations: 128,
        }; 2],
        maximum_joint_trials: 128,
    }
}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        fs_blake3::hash_bytes(b"synthetic real-gas acoustic thermal store"),
        SolidLiquidPhase::Solid,
        [200.0, 800.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 500.0 * (temperature_k - 200.0),
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 4000.0,
            })
            .to_vec(),
    )
    .unwrap()
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve.state_at_specific_enthalpy(62_500.0).unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.005,
        duration_s: input().duration_s,
    }
}
fn interface(
    gases: &[ResolvedThermochemicalRealGasState; 2],
    _: EquilibriumPhaseState,
) -> Result<[ConductanceSample<ResolvedThermochemicalRealGasState>; 2], Infallible> {
    Ok(std::array::from_fn(|i| ConductanceSample {
        conductance_w_k: [50.0, 30.0][i],
        source: gases[i].clone(),
    }))
}
fn assert_real_same(
    actual: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    expected: &DoubleActingGasPiston<IsochoricRealGasParcel>,
) {
    assert_eq!(actual, expected);
    assert_eq!(actual.ports, expected.ports);
    assert_eq!(actual.primary().ports, expected.primary().ports);
}
fn stored_change(
    initial: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    next: &DoubleActingGasPiston<IsochoricRealGasParcel>,
) -> f64 {
    let gas: f64 = [
        (initial.primary().gas(), next.primary().gas()),
        (initial.opposing(), next.opposing()),
    ]
    .into_iter()
    .map(|(a, b)| {
        assert_eq!(a.mass_kg(), b.mass_kg());
        assert_eq!(a.model(), b.model());
        a.mass_kg()
            * a.model()
                .specific_internal_energy_change_between_densities(
                    a.density_kg_m3(),
                    a.state().temperature_k(),
                    b.density_kg_m3(),
                    b.state().temperature_k(),
                )
                .unwrap()
    })
    .sum();
    let (v0, v1) = (
        initial.primary().velocity_m_s(),
        next.primary().velocity_m_s(),
    );
    gas + 0.5 * initial.primary().mass_kg() * (v1 - v0) * (v1 + v0)
}

#[test]
fn dense_chambers_drive_each_original_acoustic_history_once_and_close_energy() {
    let initial = real_pair();
    let mut piston = initial.clone();
    let mut loads = vec![
        load(0.15, input().duration_s),
        load(0.85, input().duration_s),
    ];
    let initial_loads = loads.clone();
    let mut supplied = 0.0;
    let mut exported = 0.0;
    for step in 0..3 {
        let previous = piston.clone();
        let previous_loads = loads.clone();
        let candidate = piston
            .acoustic_candidate(
                input(),
                real_controls(),
                &faces(&loads, AREAS),
                acoustic(),
                || false,
            )
            .unwrap();
        assert_real_same(&piston, &previous);
        let report = candidate.piston_step();
        let observed_force: f64 = candidate
            .motion()
            .faces
            .iter()
            .map(|face| -face.signed_area_m2 * face.mean_pressure_pa)
            .sum();
        let force = piston.primary().area_m2() * report.midpoint[0].pressure_pa()
            - piston.opposing_area_m2() * report.midpoint[1].pressure_pa()
            + input().external_force_n
            + observed_force;
        near(
            piston.primary().mass_kg() * (report.velocity_m_s[1] - report.velocity_m_s[0]),
            input().duration_s * force,
            real_controls().piston.momentum_tolerance_n_s,
        );
        for i in 0..2 {
            assert_load_same(&loads[i], &previous_loads[i]);
            let flow = AREAS[i] * (report.displacement_m / input().duration_s);
            let mut replay = previous_loads[i].clone();
            let output = replay.preview_step(flow).unwrap().accept();
            assert_eq!(candidate.motion().faces[i].duct, output);
            assert_load_same(&candidate.loads()[i], &replay);
            assert_eq!(candidate.loads()[i].state().step_count(), step + 1);
            assert_eq!(
                candidate.loads()[i].duct().gas_field(),
                initial_loads[i].duct().gas_field()
            );
        }
        near(
            candidate.motion().mechanical_work_j,
            observed_force * report.displacement_m,
            acoustic().energy_tolerance_j,
        );
        near(
            candidate.motion().transfer_residual_j,
            0.0,
            acoustic().energy_tolerance_j,
        );
        near(
            candidate.total_energy_residual_j(),
            0.0,
            acoustic().energy_tolerance_j,
        );
        if step == 0 {
            let mut unloaded = piston.clone();
            unloaded.step(input(), real_controls(), || false).unwrap();
            assert!(
                candidate.piston().primary().velocity_m_s() < unloaded.primary().velocity_m_s()
            );
            assert!(
                candidate
                    .motion()
                    .faces
                    .iter()
                    .map(|f| f.inlet_work_j)
                    .sum::<f64>()
                    > 0.0
            );
        }
        supplied += input().heat_j.iter().sum::<f64>()
            + input().external_force_n
                * ((candidate.piston().primary().gas().volume_m3()
                    - piston.primary().gas().volume_m3())
                    / piston.primary().area_m2());
        exported += candidate.motion().terminal_work_j;
        let (next, next_loads, _, _, _) = candidate.into_parts();
        piston = next;
        loads = next_loads;
    }
    let acoustic_energy: f64 = loads
        .iter()
        .zip(&initial_loads)
        .map(|(a, b)| {
            (a.stored_energy_j() - b.stored_energy_j())
                + (a.cumulative_dissipated_j() - b.cumulative_dissipated_j())
        })
        .sum();
    near(
        stored_change(&initial, &piston) + acoustic_energy + exported,
        supplied,
        3.0 * acoustic().energy_tolerance_j,
    );
}

#[test]
fn real_heat_callbacks_and_late_refusal_preserve_all_thermal_acoustic_owners() {
    let piston = real_pair();
    let saved = piston.clone();
    let loads = vec![
        load(0.15, input().duration_s),
        load(0.85, input().duration_s),
    ];
    let saved_loads = loads.clone();
    let curve = curve();
    let thermal = body(&curve);
    let mut polls = 0;
    let mut callbacks = 0;
    let candidate = piston
        .acoustic_thermal_candidate(
            &faces(&loads, AREAS),
            input(),
            thermal,
            real_controls(),
            exchange(),
            acoustic(),
            |gases, wall| {
                callbacks += 1;
                interface(gases, wall)
            },
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let native = candidate.native();
    let report = native.step();
    assert!(callbacks > 1);
    assert_eq!(callbacks, report.joint_trials);
    assert_eq!(native.thermal().report, report.piston.midpoint);
    assert!(report.exchanged_heat_j[0] > 0.0 && report.exchanged_heat_j[1] < 0.0);
    for i in 0..2 {
        near(
            report.exchanged_heat_j[i],
            input().duration_s
                * [50.0, 30.0][i]
                * (report.thermal_midpoint.temperature_k()
                    - report.piston.midpoint[i].temperature_k()),
            exchange().faces[i].source_tolerance_j,
        );
        assert_eq!(candidate.loads()[i].state().step_count(), 1);
        let mut replay = loads[i].clone();
        replay
            .preview_step(candidate.motion().faces[i].inlet_flow_m3_s)
            .unwrap()
            .accept();
        assert_load_same(&candidate.loads()[i], &replay);
    }
    let stored_body = thermal.mass_kg
        * (native.thermal().state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
    let acoustic_energy: f64 = candidate
        .loads()
        .iter()
        .zip(&loads)
        .map(|(a, b)| {
            (a.stored_energy_j() - b.stored_energy_j())
                + (a.cumulative_dissipated_j() - b.cumulative_dissipated_j())
        })
        .sum();
    near(
        stored_change(&piston, native.piston())
            + stored_body
            + acoustic_energy
            + candidate.motion().terminal_work_j
            - input().heat_j.iter().sum::<f64>()
            - thermal.internal_heat_j
            - input().external_force_n * report.piston.displacement_m,
        0.0,
        acoustic().energy_tolerance_j,
    );
    near(
        candidate.total_energy_residual_j(),
        0.0,
        acoustic().energy_tolerance_j,
    );

    let source_failure = piston.acoustic_thermal_candidate(
        &faces(&loads, AREAS),
        input(),
        thermal,
        real_controls(),
        exchange(),
        acoustic(),
        |_, _| Err::<[ConductanceSample<()>; 2], _>("opposing heat source unavailable"),
        || false,
    );
    assert!(matches!(
        source_failure,
        Err(AcousticMotionError::Mechanical(
            ThermalPistonError::Interface("opposing heat source unavailable")
        ))
    ));
    let mut count = 0;
    let refused = piston.acoustic_thermal_candidate(
        &faces(&loads, AREAS),
        input(),
        thermal,
        real_controls(),
        exchange(),
        acoustic(),
        interface,
        || {
            count += 1;
            count == polls
        },
    );
    assert!(matches!(refused, Err(AcousticMotionError::Cancelled)));
    assert_eq!(count, polls);
    assert_real_same(&piston, &saved);
    for (actual, before) in loads.iter().zip(&saved_loads) {
        assert_load_same(actual, before);
    }
    assert_eq!(thermal.initial, body(&curve).initial);
}

#[test]
fn zero_acoustic_area_preserves_native_real_motion_and_advances_existing_memory() {
    let initial = real_pair();
    let cold = vec![
        load(0.15, input().duration_s),
        load(0.85, input().duration_s),
    ];
    let warm = initial
        .acoustic_candidate(
            input(),
            real_controls(),
            &faces(&cold, AREAS),
            acoustic(),
            || false,
        )
        .unwrap();
    let (piston, loads, _, _, _) = warm.into_parts();
    let mut expected = piston.clone();
    let native = expected.step(input(), real_controls(), || false).unwrap();
    let candidate = piston
        .acoustic_candidate(
            input(),
            real_controls(),
            &faces(&loads, [0.0; 2]),
            acoustic(),
            || false,
        )
        .unwrap();
    assert_real_same(candidate.piston(), &expected);
    assert_eq!(candidate.piston_step(), &native);
    assert_eq!(candidate.motion().force_n, 0.0);
    assert_eq!(candidate.motion().mechanical_work_j, 0.0);
    for i in 0..2 {
        let mut free = loads[i].clone();
        let output = free.preview_step(0.0).unwrap().accept();
        assert_eq!(candidate.motion().faces[i].duct, output);
        assert_load_same(&candidate.loads()[i], &free);
        assert_ne!(candidate.loads()[i].state(), loads[i].state());
        assert_eq!(candidate.loads()[i].state().step_count(), 2);
    }
}

#[test]
fn material_density_and_sources_survive_loaded_heat_continuation_and_last_wrapper_poll() {
    use crate::cavity_phs::piston::material::{
        MaterialDoubleActingGasPiston, tests::body as specimen,
    };

    let gas_pair = real_pair();
    let curve = curve();
    // Resolve the combined gas spring even for the lighter source specimen.
    let drive = DoublePistonInput {
        duration_s: 2e-6,
        ..input()
    };
    let thermal = UniformEnthalpyStepInput {
        duration_s: drive.duration_s,
        ..body(&curve)
    };
    let heat = DoubleHeatExchangeSolve {
        faces: exchange().faces.map(|face| HeatExchangeSolve {
            heat_bracket_j: [-0.025, 0.025],
            ..face
        }),
        ..exchange()
    };
    let mut velocities = Vec::new();
    for density in [1000.0, 8000.0] {
        let source_body = specimen(density);
        let wall = MaterialDoubleActingGasPiston::new(
            gas_pair.primary().gas().clone(),
            gas_pair.opposing().clone(),
            source_body.clone(),
            [gas_pair.primary().area_m2(), gas_pair.opposing_area_m2()],
            0.0,
        )
        .unwrap();
        let loads = vec![load(0.15, drive.duration_s), load(0.85, drive.duration_s)];
        let native = wall
            .piston()
            .acoustic_candidate(
                drive,
                real_controls(),
                &faces(&loads, AREAS),
                acoustic(),
                || false,
            )
            .unwrap();
        let bound = wall
            .acoustic_candidate(
                drive,
                real_controls(),
                &faces(&loads, AREAS),
                acoustic(),
                || false,
            )
            .unwrap();
        assert_eq!(bound.piston().body(), &source_body);
        assert_eq!(
            bound.piston().piston().primary().mass_kg(),
            source_body.inertia().mass()
        );
        assert_real_same(bound.piston().piston(), native.piston());
        assert_eq!(bound.piston_step(), native.piston_step());
        assert_eq!(bound.motion(), native.motion());
        assert_eq!(
            bound.total_energy_residual_j(),
            native.total_energy_residual_j()
        );
        for (a, b) in bound.loads().iter().zip(native.loads()) {
            assert_load_same(a, b);
        }
        let (wall, loads, _, _, _) = bound.into_parts();
        let saved = wall.clone();
        let saved_loads = loads.clone();

        let native = wall
            .piston()
            .acoustic_thermal_candidate(
                &faces(&loads, AREAS),
                drive,
                thermal,
                real_controls(),
                heat,
                acoustic(),
                interface,
                || false,
            )
            .unwrap();
        let mut polls = 0;
        let bound = wall
            .acoustic_thermal_candidate(
                &faces(&loads, AREAS),
                drive,
                thermal,
                real_controls(),
                heat,
                acoustic(),
                interface,
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        assert_eq!(bound.native().piston().body(), &source_body);
        assert_real_same(bound.native().piston().piston(), native.native().piston());
        assert_eq!(bound.native().thermal(), native.native().thermal());
        assert_eq!(bound.native().step(), native.native().step());
        assert_eq!(bound.motion(), native.motion());
        assert_eq!(
            bound.total_energy_residual_j(),
            native.total_energy_residual_j()
        );
        for (a, b) in bound.loads().iter().zip(native.loads()) {
            assert_load_same(a, b);
            assert_eq!(a.state().step_count(), 2);
        }
        for (gas, original) in bound
            .native()
            .piston()
            .gases()
            .into_iter()
            .zip([gas_pair.primary().gas(), gas_pair.opposing()])
        {
            assert_eq!(gas.mass_kg(), original.mass_kg());
            assert_eq!(gas.model(), original.model());
        }
        velocities.push(bound.native().piston().piston().primary().velocity_m_s());

        let mut count = 0;
        let cancelled = wall.acoustic_thermal_candidate(
            &faces(&loads, AREAS),
            drive,
            thermal,
            real_controls(),
            heat,
            acoustic(),
            interface,
            || {
                count += 1;
                count == polls
            },
        );
        assert!(matches!(cancelled, Err(AcousticMotionError::Cancelled)));
        assert_eq!(count, polls);
        assert_eq!(wall.body(), saved.body());
        assert_real_same(wall.piston(), saved.piston());
        for (a, b) in loads.iter().zip(&saved_loads) {
            assert_load_same(a, b);
        }
        assert_eq!(thermal.initial, body(&curve).initial);
    }
    assert!(velocities[0] > velocities[1] && velocities[1] > 0.0);
}
