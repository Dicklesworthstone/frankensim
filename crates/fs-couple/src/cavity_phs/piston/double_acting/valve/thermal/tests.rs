//! Original valve transport, finite phase storage and both native gas faces.
use super::super::tests::{
    assert_unchanged, controls, exchange as flow_solve, geometry, input, internal_energy, near,
    pair, species_masses,
};
use super::*;
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
};

pub(in crate::cavity_phs::piston::double_acting) type Source = (ResolvedThermochemicalGasState, EquilibriumPhaseState);

pub(in crate::cavity_phs::piston::double_acting) fn exchange() -> DoubleHeatExchangeSolve {
    DoubleHeatExchangeSolve {
        faces: [HeatExchangeSolve {
            heat_bracket_j: [-0.5, 0.5],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        }; 2],
        maximum_joint_trials: 4096,
    }
}

pub(in crate::cavity_phs::piston::double_acting) fn curve(latent: bool) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, t, x| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: x,
        bulk_density_kg_m3: 4000.0,
    };
    if latent {
        EquilibriumEnthalpyPhaseCurve::try_new(
            fs_blake3::hash_bytes(b"synthetic valve latent store"),
            vec![
                knot(0.0, 200.0, 0.0),
                knot(100000.0, 350.0, 0.0),
                knot(160000.0, 350.0, 1.0),
                knot(400000.0, 800.0, 1.0),
            ],
        )
        .unwrap()
    } else {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            fs_blake3::hash_bytes(b"synthetic valve body cp=500"),
            SolidLiquidPhase::Solid,
            vec![knot(0.0, 200.0, 0.0), knot(300000.0, 800.0, 0.0)],
        )
        .unwrap()
    }
}

pub(in crate::cavity_phs::piston::double_acting) fn body(curve: &EquilibriumEnthalpyPhaseCurve, latent: bool) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: curve
            .state_at_specific_enthalpy(if latent { 130000.0 } else { 200000.0 })
            .unwrap(),
        mass_kg: 0.001,
        volume_m3: 0.001 / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.0,
        duration_s: input().duration_s,
    }
}

// Synthetic state-dependent conductances keep the physical midpoint feedback live.
pub(in crate::cavity_phs::piston::double_acting) fn conductance(
    base: f64,
    gas: &ResolvedThermochemicalGasState,
    body: EquilibriumPhaseState,
) -> f64 {
    base * gas.state().temperature / 400.0 * (1.0 + 0.1 * body.liquid_mass_fraction())
}

pub(in crate::cavity_phs::piston::double_acting) fn links(
    base: [f64; 2],
) -> impl FnMut(
    &[ResolvedThermochemicalGasState; 2],
    EquilibriumPhaseState,
) -> Result<[ConductanceSample<Source>; 2], &'static str> {
    move |gases, body| {
        Ok(std::array::from_fn(|i| ConductanceSample {
            conductance_w_k: conductance(base[i], &gases[i], body),
            source: (gases[i].clone(), body),
        }))
    }
}

#[test]
fn g1_transport_and_both_faces_share_one_finite_sensible_or_latent_body() {
    for latent in [false, true] {
        let original = pair([100.0, 0.0], [0.001; 2], 0.03);
        let saved = original.clone();
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
        let native = original
            .valve_candidate(geometry(), forcing, flow_solve(), controls(), 1e-7, || {
                false
            })
            .unwrap();
        let candidate = original
            .thermal_valve_candidate(
                geometry(),
                forcing,
                thermal,
                controls(),
                flow_solve(),
                exchange(),
                1e-7,
                links([1.0, 0.3]),
                || false,
            )
            .unwrap();
        let next = candidate.piston();
        let report = candidate.step();
        // Heat feedback must not resample transport from a trial midpoint.
        assert_eq!(report.valve.flow, native.step().flow);
        assert_eq!(report.valve.transfer, native.step().transfer);
        assert!(report.valve.transfer.transferred_mass_kg > 0.0);
        assert_ne!(
            next.primary().velocity_m_s(),
            native.piston().primary().velocity_m_s()
        );
        assert_unchanged(&original, &saved);
        let old_gas = [original.primary().gas(), original.opposing()];
        let new_gas = [next.primary().gas(), next.opposing()];
        let mut gas_change = 0.0;
        for i in 0..2 {
            let delta = internal_energy(new_gas[i], 0.001) - internal_energy(old_gas[i], 0.001);
            near(report.valve.gas_energy_change_j[i], delta, 1e-8);
            gas_change += delta;
            let gas_midpoint = &report.valve.piston.midpoint[i];
            let g = conductance([1.0, 0.3][i], gas_midpoint, report.thermal_midpoint);
            near(report.conductance_w_k[i], g, 1e-14);
            let q = forcing.duration_s
                * g
                * (report.thermal_midpoint.temperature_k() - gas_midpoint.state().temperature);
            near(
                report.exchanged_heat_j[i],
                q,
                exchange().faces[i].source_tolerance_j,
            );
            assert_eq!(
                candidate.thermal().report[i],
                (gas_midpoint.clone(), report.thermal_midpoint)
            );
            assert_eq!(
                report.valve.piston.gas[i].supplied_heat_j,
                forcing.heat_j[i] + report.exchanged_heat_j[i]
            );
        }
        for species in 0..2 {
            near(
                species_masses(new_gas[0])[species] + species_masses(new_gas[1])[species],
                species_masses(old_gas[0])[species] + species_masses(old_gas[1])[species],
                1e-14,
            );
        }
        assert_ne!(
            new_gas[1].model().thermodynamics().mole_composition(),
            old_gas[1].model().thermodynamics().mole_composition()
        );
        let stored = thermal.mass_kg
            * (candidate.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let shared_heat: f64 = report.exchanged_heat_j.iter().sum();
        near(stored + shared_heat, thermal.internal_heat_j, 1e-9);
        assert_eq!(candidate.thermal().external_heat_j, -shared_heat);
        near(
            report.thermal_midpoint.specific_enthalpy_j_kg(),
            thermal.initial.specific_enthalpy_j_kg() + 0.5 * stored / thermal.mass_kg,
            1e-8,
        );
        let kinetic = 0.5
            * next.primary().mass_kg()
            * (next.primary().velocity_m_s().powi(2) - original.primary().velocity_m_s().powi(2));
        let displacement = (next.primary().gas().volume_m3()
            - original.primary().gas().volume_m3())
            / original.primary().area_m2();
        let residual = gas_change + kinetic + stored
            - forcing.heat_j.iter().sum::<f64>()
            - thermal.internal_heat_j
            - forcing.external_force_n * displacement;
        near(residual, 0.0, 2e-7);
        near(report.total_energy_residual_j, residual, 1e-8);
        assert!(report.joint_trials > 1 && report.joint_trials <= exchange().maximum_joint_trials);
        if latent {
            assert!(report.exchanged_heat_j.iter().all(|q| *q < 0.0));
            assert_eq!(candidate.thermal().state.temperature_k(), 350.0);
            assert_eq!(report.thermal_midpoint.temperature_k(), 350.0);
            assert!(
                candidate.thermal().state.liquid_mass_fraction()
                    > thermal.initial.liquid_mass_fraction()
            );
        } else {
            assert!(report.exchanged_heat_j.iter().all(|q| *q > 0.0));
            assert!(candidate.thermal().state.temperature_k() < thermal.initial.temperature_k());
        }
    }
}

#[test]
fn g3_insulation_and_closed_metering_reduce_to_the_original_native_owners() {
    let original = pair([100.0, 0.0], [0.001; 2], 0.03);
    let chart = curve(false);
    let thermal = body(&chart, false);
    let forcing = DoublePistonInput {
        heat_j: [0.02, 0.01],
        ..input()
    };
    let native = original
        .valve_candidate(geometry(), forcing, flow_solve(), controls(), 1e-7, || {
            false
        })
        .unwrap();
    let insulated = original
        .thermal_valve_candidate(
            geometry(),
            forcing,
            thermal,
            controls(),
            flow_solve(),
            exchange(),
            1e-7,
            links([0.0; 2]),
            || false,
        )
        .unwrap();
    assert_unchanged(insulated.piston(), native.piston());
    assert_eq!(insulated.step().valve, *native.step());
    assert_eq!(insulated.thermal().state, thermal.initial);
    assert_eq!(insulated.step().exchanged_heat_j, [0.0; 2]);
    assert_eq!(insulated.step().joint_trials, 1);

    let closed = RectangularMeteringValve {
        edge_height_at_reference_m: -0.01,
        ..geometry()
    };
    let native = original
        .thermal_candidate(
            forcing,
            thermal,
            controls(),
            exchange(),
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    let coupled = original
        .thermal_valve_candidate(
            closed,
            forcing,
            thermal,
            controls(),
            flow_solve(),
            exchange(),
            1e-7,
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    assert_unchanged(coupled.piston(), native.piston());
    assert_eq!(coupled.thermal(), native.thermal());
    assert_eq!(coupled.step().valve.piston, native.step().piston);
    assert_eq!(coupled.step().valve.direction, None);
    assert_eq!(coupled.step().valve.transfer.transferred_mass_kg, 0.0);
    assert_eq!(
        coupled.step().thermal_midpoint,
        native.step().thermal_midpoint
    );
    assert_eq!(
        coupled.step().exchanged_heat_j,
        native.step().exchanged_heat_j
    );
    assert_eq!(
        coupled.step().source_residual_j,
        native.step().source_residual_j
    );
    assert_eq!(
        coupled.step().material_trials,
        native.step().material_trials
    );
}

#[test]
fn g4_late_source_budget_and_final_cancel_return_no_partial_system_and_retry_exactly() {
    let original = pair([100.0, 0.0], [0.001; 2], 0.03);
    let saved = original.clone();
    let chart = curve(false);
    let thermal = body(&chart, false);
    let mut called = 0;
    let mut source = links([1.0, 0.3]);
    let refused = original.thermal_valve_candidate(
        geometry(),
        input(),
        thermal,
        controls(),
        flow_solve(),
        exchange(),
        1e-7,
        |gas, body| {
            called += 1;
            if called == 2 {
                Err("second face source unavailable")
            } else {
                source(gas, body)
            }
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(ThermalPistonError::Interface(
            "second face source unavailable"
        ))
    ));
    assert_eq!(called, 2);
    let refused = original.thermal_valve_candidate(
        geometry(),
        input(),
        thermal,
        controls(),
        flow_solve(),
        DoubleHeatExchangeSolve {
            maximum_joint_trials: 1,
            ..exchange()
        },
        1e-7,
        links([1.0, 0.3]),
        || false,
    );
    assert!(matches!(
        refused,
        Err(ThermalPistonError::Invalid(
            "shared thermal joint trial budget exhausted"
        ))
    ));
    assert_unchanged(&original, &saved);
    assert_eq!(thermal.initial, body(&chart, false).initial);

    let mut polls = 0;
    let expected = original
        .thermal_valve_candidate(
            geometry(),
            input(),
            thermal,
            controls(),
            flow_solve(),
            exchange(),
            1e-7,
            links([1.0, 0.3]),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    let cancelled = original.thermal_valve_candidate(
        geometry(),
        input(),
        thermal,
        controls(),
        flow_solve(),
        exchange(),
        1e-7,
        links([1.0, 0.3]),
        || {
            seen += 1;
            seen == polls
        },
    );
    assert!(matches!(cancelled, Err(ThermalPistonError::Cancelled)));
    assert_eq!(seen, polls);
    assert_unchanged(&original, &saved);
    let retried = original
        .thermal_valve_candidate(
            geometry(),
            input(),
            thermal,
            controls(),
            flow_solve(),
            exchange(),
            1e-7,
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    assert_eq!(retried.step(), expected.step());
    assert_eq!(retried.thermal(), expected.thermal());
    assert_unchanged(retried.piston(), expected.piston());
    let (retried, accepted, _) = retried.into_parts();
    let (expected, reference, _) = expected.into_parts();
    assert_eq!(accepted, reference);
    let next_body = UniformEnthalpyStepInput {
        initial: accepted.state,
        ..thermal
    };
    let next = retried
        .thermal_valve_candidate(
            geometry(),
            input(),
            next_body,
            controls(),
            flow_solve(),
            exchange(),
            1e-7,
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    let clean = expected
        .thermal_valve_candidate(
            geometry(),
            input(),
            next_body,
            controls(),
            flow_solve(),
            exchange(),
            1e-7,
            links([1.0, 0.3]),
            || false,
        )
        .unwrap();
    assert_eq!(next.step(), clean.step());
    assert_eq!(next.thermal(), clean.thermal());
    assert_unchanged(next.piston(), clean.piston());
}
