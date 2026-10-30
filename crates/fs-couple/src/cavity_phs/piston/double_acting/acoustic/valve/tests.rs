//! Native mass, finite-body heat and acoustic reaction share one moving wall.
use super::super::tests::{assert_load_same, load, motion};
use super::*;
use crate::cavity_phs::piston::double_acting::valve::tests::{
    assert_unchanged, controls, exchange as flow_solve, geometry, input, internal_energy, near,
    pair, species_masses,
};
use crate::cavity_phs::piston::double_acting::valve::thermal::tests::{
    body, conductance, curve, exchange as heat_solve, links,
};
use crate::cavity_phs::piston::thermal::ThermalPistonError;
use crate::dg_port::DgPortLoad;
use fs_material::phase::UniformEnthalpyStepInput;

const AREAS: [f64; 2] = [0.001, -0.0007];

fn faces(loads: &[DgPortLoad], areas: [f64; 2]) -> [AcousticMotionFace<'_>; 2] {
    std::array::from_fn(|i| AcousticMotionFace {
        load: &loads[i],
        signed_area_m2: areas[i],
    })
}

fn gas_and_kinetic_change(original: &DoubleActingGasPiston, next: &DoubleActingGasPiston) -> f64 {
    let before = [original.primary().gas(), original.opposing()];
    let after = [next.primary().gas(), next.opposing()];
    let gas: f64 = (0..2)
        .map(|i| internal_energy(after[i], 0.001) - internal_energy(before[i], 0.001))
        .sum();
    let kinetic = 0.5
        * original.primary().mass_kg()
        * (next.primary().velocity_m_s().powi(2) - original.primary().velocity_m_s().powi(2));
    gas + kinetic
}

// Original storage and loss getters plus actual external terminal work, without
// reconstructing the result from the new coupling's enclosing residual.
fn acoustic_change(before: &[DgPortLoad], after: &[DgPortLoad]) -> f64 {
    before
        .iter()
        .zip(after)
        .map(|(old, new)| {
            let terminal: f64 = new
                .last_terminal_frames()
                .iter()
                .map(|frame| frame.load_work_j)
                .sum();
            (new.stored_energy_j() - old.stored_energy_j())
                + (new.cumulative_dissipated_j() - old.cumulative_dissipated_j())
                + terminal
        })
        .sum()
}

#[test]
fn g1_bidirectional_metering_drives_original_acoustics_and_closes_species_and_power() {
    for heat in [[100.0, 0.0], [0.0, 100.0]] {
        let original = pair(heat, [0.001; 2], 0.03);
        let saved = original.clone();
        let forcing = DoublePistonInput {
            heat_j: [0.002, 0.001],
            external_force_n: 3.0,
            ..input()
        };
        let loads = [
            load(0.75, forcing.duration_s),
            load(0.25, forcing.duration_s),
        ];
        let old_loads = loads.clone();
        let native = original
            .valve_candidate(
                geometry(),
                forcing,
                flow_solve(),
                controls(),
                motion().energy_tolerance_j,
                || false,
            )
            .unwrap();
        let candidate = original
            .acoustic_valve_candidate(
                &faces(&loads, AREAS),
                geometry(),
                forcing,
                flow_solve(),
                controls(),
                motion(),
                || false,
            )
            .unwrap();
        let next = candidate.native().piston();
        let report = candidate.native().step();
        // Acoustic reaction changes the actual trajectory, but this partition
        // must not sample or debit the original nozzle more than once.
        assert_eq!(report.flow, native.step().flow);
        assert_eq!(report.transfer, native.step().transfer);
        assert!(report.transfer.transferred_mass_kg > 0.0);
        assert_ne!(
            next.primary().velocity_m_s(),
            native.piston().primary().velocity_m_s()
        );
        let before = [original.primary().gas(), original.opposing()];
        let after = [next.primary().gas(), next.opposing()];
        for species in 0..2 {
            near(
                species_masses(after[0])[species] + species_masses(after[1])[species],
                species_masses(before[0])[species] + species_masses(before[1])[species],
                1e-14,
            );
        }
        let displacement = (next.primary().gas().volume_m3()
            - original.primary().gas().volume_m3())
            / original.primary().area_m2();
        let velocity = displacement / forcing.duration_s;
        let mut inlet_work = 0.0;
        for i in 0..2 {
            let flow = AREAS[i] * velocity;
            let mut replay = loads[i].clone();
            let trial = replay.preview_step(flow).unwrap();
            inlet_work += trial.inlet_work_j();
            let expected = trial.accept();
            assert_eq!(candidate.motion().faces[i].duct, expected);
            assert_load_same(&candidate.loads()[i], &replay);
            near(candidate.motion().faces[i].inlet_flow_m3_s, flow, 1e-18);
            assert!(candidate.loads()[i].stored_energy_j() > 0.0);
            assert_eq!(candidate.loads()[i].state().step_count(), 1);
        }
        near(
            report.piston.external_force_work_j - forcing.external_force_n * displacement
                + inlet_work,
            0.0,
            motion().energy_tolerance_j,
        );
        let residual = gas_and_kinetic_change(&original, next)
            + acoustic_change(&loads, candidate.loads())
            - forcing.heat_j.iter().sum::<f64>()
            - forcing.external_force_n * displacement;
        near(residual, 0.0, 2e-7);
        near(candidate.total_energy_residual_j(), residual, 1e-8);
        assert_unchanged(&original, &saved);
        for i in 0..2 {
            assert_load_same(&loads[i], &old_loads[i]);
        }
    }
}

#[test]
fn g3_sensible_and_latent_heat_reach_loaded_motion_and_zero_area_is_exactly_native() {
    for latent in [false, true] {
        let original = pair([100.0, 0.0], [0.001; 2], 0.03);
        let chart = curve(latent);
        let thermal = UniformEnthalpyStepInput {
            internal_heat_j: 0.004,
            ..body(&chart, latent)
        };
        let forcing = DoublePistonInput {
            heat_j: [0.002, 0.001],
            external_force_n: 3.0,
            ..input()
        };
        let loads = [
            load(0.75, forcing.duration_s),
            load(0.25, forcing.duration_s),
        ];
        let dry = original
            .acoustic_valve_candidate(
                &faces(&loads, AREAS),
                geometry(),
                forcing,
                flow_solve(),
                controls(),
                motion(),
                || false,
            )
            .unwrap();
        let candidate = original
            .acoustic_thermal_valve_candidate(
                &faces(&loads, AREAS),
                geometry(),
                forcing,
                thermal,
                controls(),
                flow_solve(),
                heat_solve(),
                motion(),
                links([1.0, 0.3]),
                || false,
            )
            .unwrap();
        let native = candidate.native();
        let next = native.piston();
        let report = native.step();
        assert_eq!(report.valve.flow, dry.native().step().flow);
        assert_eq!(report.valve.transfer, dry.native().step().transfer);
        assert_ne!(
            next.primary().velocity_m_s(),
            dry.native().piston().primary().velocity_m_s()
        );
        for i in 0..2 {
            let gas = &report.valve.piston.midpoint[i];
            let g = conductance([1.0, 0.3][i], gas, report.thermal_midpoint);
            let q = forcing.duration_s
                * g
                * (report.thermal_midpoint.temperature_k() - gas.state().temperature);
            near(report.conductance_w_k[i], g, 1e-14);
            near(
                report.exchanged_heat_j[i],
                q,
                heat_solve().faces[i].source_tolerance_j,
            );
            assert_eq!(
                native.thermal().report[i],
                (gas.clone(), report.thermal_midpoint)
            );
        }
        let stored = thermal.mass_kg
            * (native.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        near(
            stored + report.exchanged_heat_j.iter().sum::<f64>(),
            thermal.internal_heat_j,
            1e-9,
        );
        let displacement = (next.primary().gas().volume_m3()
            - original.primary().gas().volume_m3())
            / original.primary().area_m2();
        let residual = gas_and_kinetic_change(&original, next)
            + stored
            + acoustic_change(&loads, candidate.loads())
            - forcing.heat_j.iter().sum::<f64>()
            - thermal.internal_heat_j
            - forcing.external_force_n * displacement;
        near(residual, 0.0, 2e-7);
        near(candidate.total_energy_residual_j(), residual, 1e-8);
        if latent {
            assert_eq!(native.thermal().state.temperature_k(), 350.0);
            assert!(
                native.thermal().state.liquid_mass_fraction()
                    > thermal.initial.liquid_mass_fraction()
            );
        } else {
            assert!(native.thermal().state.temperature_k() < thermal.initial.temperature_k());
        }
        // Removing radiator area preserves the original heat/transport/motion
        // exactly while the acoustic clocks still advance under zero flow.
        let reference = original
            .thermal_valve_candidate(
                geometry(),
                forcing,
                thermal,
                controls(),
                flow_solve(),
                heat_solve(),
                motion().energy_tolerance_j,
                links([1.0, 0.3]),
                || false,
            )
            .unwrap();
        let disconnected = original
            .acoustic_thermal_valve_candidate(
                &faces(&loads, [0.0; 2]),
                geometry(),
                forcing,
                thermal,
                controls(),
                flow_solve(),
                heat_solve(),
                motion(),
                links([1.0, 0.3]),
                || false,
            )
            .unwrap();
        assert_eq!(disconnected.native().step(), reference.step());
        assert_eq!(disconnected.native().thermal(), reference.thermal());
        assert_unchanged(disconnected.native().piston(), reference.piston());
        for i in 0..2 {
            let mut replay = loads[i].clone();
            replay.preview_step(0.0).unwrap().accept();
            assert_load_same(&disconnected.loads()[i], &replay);
        }

        // A closed metering window reduces to the direct finite-body/acoustic
        // path, without a distinct mechanical, thermal or propagation model.
        let closed = RectangularMeteringValve {
            discharge_coefficient: 0.0,
            ..geometry()
        };
        let direct = original
            .acoustic_thermal_candidate(
                &faces(&loads, AREAS),
                forcing,
                thermal,
                controls(),
                heat_solve(),
                motion(),
                links([1.0, 0.3]),
                || false,
            )
            .unwrap();
        let closed_valve = original
            .acoustic_thermal_valve_candidate(
                &faces(&loads, AREAS),
                closed,
                forcing,
                thermal,
                controls(),
                flow_solve(),
                heat_solve(),
                motion(),
                links([1.0, 0.3]),
                || false,
            )
            .unwrap();
        assert_eq!(
            closed_valve.native().step().valve.piston,
            direct.native().step().piston
        );
        assert_eq!(closed_valve.native().thermal(), direct.native().thermal());
        assert_unchanged(closed_valve.native().piston(), direct.native().piston());
        for i in 0..2 {
            assert_load_same(&closed_valve.loads()[i], &direct.loads()[i]);
        }
    }
}

#[test]
fn g4_late_thermal_refusal_and_final_cancel_preserve_all_owners_and_exact_suffix() {
    let original = pair([100.0, 0.0], [0.001; 2], 0.03);
    let saved = original.clone();
    let chart = curve(false);
    let thermal = body(&chart, false);
    let forcing = input();
    let loads = [
        load(0.75, forcing.duration_s),
        load(0.25, forcing.duration_s),
    ];
    let old_loads = loads.clone();
    let mut calls = 0;
    let mut source = links([1.0, 0.3]);
    let refused = original.acoustic_thermal_valve_candidate(
        &faces(&loads, AREAS),
        geometry(),
        forcing,
        thermal,
        controls(),
        flow_solve(),
        heat_solve(),
        motion(),
        |gas, body| {
            calls += 1;
            if calls == 2 {
                Err("late thermal source refusal")
            } else {
                source(gas, body)
            }
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(AcousticMotionError::Mechanical(
            ThermalPistonError::Interface("late thermal source refusal")
        ))
    ));
    assert_eq!(calls, 2);
    let mut limited = heat_solve();
    limited.maximum_joint_trials = 1;
    assert!(matches!(
        original.acoustic_thermal_valve_candidate(
            &faces(&loads, AREAS),
            geometry(),
            forcing,
            thermal,
            controls(),
            flow_solve(),
            limited,
            motion(),
            links([1.0, 0.3]),
            || false,
        ),
        Err(AcousticMotionError::Mechanical(
            ThermalPistonError::Invalid("shared thermal joint trial budget exhausted")
        ))
    ));
    let mut polls = 0;
    let expected = original
        .acoustic_thermal_valve_candidate(
            &faces(&loads, AREAS),
            geometry(),
            forcing,
            thermal,
            controls(),
            flow_solve(),
            heat_solve(),
            motion(),
            links([1.0, 0.3]),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    let cancelled = original.acoustic_thermal_valve_candidate(
        &faces(&loads, AREAS),
        geometry(),
        forcing,
        thermal,
        controls(),
        flow_solve(),
        heat_solve(),
        motion(),
        links([1.0, 0.3]),
        || {
            seen += 1;
            seen == polls
        },
    );
    assert!(matches!(cancelled, Err(AcousticMotionError::Cancelled)));
    assert_eq!(seen, polls);
    assert_unchanged(&original, &saved);
    assert_eq!(thermal.initial, body(&chart, false).initial);
    for i in 0..2 {
        assert_load_same(&loads[i], &old_loads[i]);
    }
    let retry = original
        .acoustic_thermal_valve_candidate(
            &faces(&loads, AREAS),
            geometry(),
            forcing,
            thermal,
            controls(),
            flow_solve(),
            heat_solve(),
            motion(),
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    assert_eq!(retry.native().step(), expected.native().step());
    assert_eq!(retry.native().thermal(), expected.native().thermal());
    assert_unchanged(retry.native().piston(), expected.native().piston());
    for i in 0..2 {
        assert_load_same(&retry.loads()[i], &expected.loads()[i]);
    }

    let (retry, retry_loads, _, _) = retry.into_parts();
    let (expected, expected_loads, _, _) = expected.into_parts();
    let (retry, accepted, _) = retry.into_parts();
    let (expected, reference, _) = expected.into_parts();
    assert_eq!(accepted, reference);
    let next_body = UniformEnthalpyStepInput {
        initial: accepted.state,
        ..thermal
    };
    let next = retry
        .acoustic_thermal_valve_candidate(
            &faces(&retry_loads, AREAS),
            geometry(),
            forcing,
            next_body,
            controls(),
            flow_solve(),
            heat_solve(),
            motion(),
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    let clean = expected
        .acoustic_thermal_valve_candidate(
            &faces(&expected_loads, AREAS),
            geometry(),
            forcing,
            next_body,
            controls(),
            flow_solve(),
            heat_solve(),
            motion(),
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    assert_eq!(next.native().step(), clean.native().step());
    assert_eq!(next.native().thermal(), clean.native().thermal());
    assert_unchanged(next.native().piston(), clean.native().piston());
    for i in 0..2 {
        assert_load_same(&next.loads()[i], &clean.loads()[i]);
    }
}
