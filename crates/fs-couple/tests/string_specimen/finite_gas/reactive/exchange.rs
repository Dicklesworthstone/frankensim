//! G1/G3/G4: one reacting endpoint closes chemistry AND reciprocal heat.

use super::super::exchange::quiet_body;
use super::*;
use fs_couple::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, layers::MaterialConductionLayers,
};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

const G: f64 = 0.1;
fn exchange_solve() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-8e-5, 8e-5],
        source_tolerance_j: 1e-13,
        energy_tolerance_j: 2e-10,
        maximum_iterations: 64,
    }
}
fn conductance(
    _: &ResolvedThermochemicalGasState,
    _: EquilibriumPhaseState,
) -> Result<ConductanceSample<()>, &'static str> {
    Ok(ConductanceSample {
        conductance_w_k: G,
        source: (),
    })
}

// Independent elimination of the TWO constant-capacity thermal balances,
// followed by scalar bisection of the exact isomer mass-action equation.
fn reference(solid: f64, gas: &IsochoricGasParcel, cs: f64, heat: FiniteGasStringHeat) -> [f64; 4] {
    let n = gas.mass_kg() / 0.028;
    let cg = n * 2.5 * R;
    let endpoint = |xi: f64| {
        let source_s = solid + heat.string_heat_j / cs;
        let source_g = gas.state().state().temperature + (heat.gas_heat_j + RELEASE * xi) / cg;
        let q = DT * G * (source_s - source_g) / (1.0 + DT * G * (1.0 / cs + 1.0 / cg));
        let ts = source_s - q / cs;
        let tg = source_g + q / cg;
        let k = (RELEASE / (R * tg)).exp();
        let activity_factor = R * tg / (gas.volume_m3() * 100000.0);
        let residual = xi
            - DT * gas.volume_m3() * RATE * activity_factor * (0.9 * n - xi - (0.1 * n + xi) / k);
        [residual, q, ts, tg]
    };
    let (mut lo, mut hi) = (0.0, 0.9 * n);
    assert!(endpoint(lo)[0] < 0.0 && endpoint(hi)[0] > 0.0);
    for _ in 0..80 {
        let xi = 0.5 * (lo + hi);
        if endpoint(xi)[0] > 0.0 {
            hi = xi;
        } else {
            lo = xi;
        }
    }
    let xi = 0.5 * (lo + hi);
    let [_, q, ts, tg] = endpoint(xi);
    [xi, q, ts, tg]
}

#[test]
fn g1_reacting_exchange_solves_both_endpoint_equations_in_both_heat_directions() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        for solid in [300.0, 330.0] {
            let mut gas = gas(&reaction);
            let (mut live, _, _, _) = quiet_body(cx, gas.state().state(), solid);
            let initial = live.clone();
            let old_gas = gas.clone();
            let cs = live.mechanical().specimen().mass_kg() * 0.001;
            let heat = FiniteGasStringHeat {
                string_heat_j: 2.0 * cs,
                gas_heat_j: 0.125,
            };
            let [xi, q, ts, tg] = reference(solid, &gas, cs, heat);
            let result = live
                .step_with_reacting_isochoric_gas_exchange(
                    cx,
                    0,
                    &[0.0],
                    &mut gas,
                    heat,
                    &reaction,
                    controls(),
                    exchange_solve(),
                    conductance,
                )
                .unwrap();
            let frame = &result.coupled;
            near(result.chemistry.reaction.extent_mol, xi, 3e-12);
            near(frame.solid_to_gas_heat_j, q, 2e-12);
            near(live.thermal().temperature_k(), ts, 1e-6);
            near(gas.state().state().temperature, tg, 2e-8);
            assert!(frame.solid_to_gas_heat_j * (ts - tg) > 0.0);
            assert_eq!(frame.coupled.drag_heat_j, 0.0);
            assert_eq!(frame.coupled.gas, result.chemistry.heat);
            assert_eq!(
                frame.coupled.gas.supplied_heat_j,
                heat.gas_heat_j + frame.solid_to_gas_heat_j
            );
            assert_eq!(
                frame.coupled.string.external_heat_j,
                heat.string_heat_j - frame.solid_to_gas_heat_j
            );
            assert_eq!(frame.external_string_heat_j, heat.string_heat_j);
            let solid_storage = live.mechanical().specimen().mass_kg()
                * (live.thermal().specific_enthalpy_j_kg()
                    - initial.thermal().specific_enthalpy_j_kg());
            let gas_storage = gas.mass_kg()
                * (gas.state().specific_internal_energy_j_kg()
                    - old_gas.state().specific_internal_energy_j_kg());
            assert!(
                (solid_storage + gas_storage - heat.string_heat_j - heat.gas_heat_j).abs()
                    <= exchange_solve().energy_tolerance_j
            );
            assert!(frame.source_residual_j.abs() <= exchange_solve().source_tolerance_j);
            assert!(
                frame.coupled.energy_balance_residual_j.abs() <= frame.coupled.energy_tolerance_j
            );
            assert_eq!(live.mechanical().ambient(), gas.state().state());
            assert_eq!(live.mechanical().epoch(), 1);
            assert_eq!(live.mechanical().accepted_samples(), 1);
            assert_eq!(
                result
                    .chemistry
                    .reaction
                    .implicit_mechanism
                    .as_ref()
                    .unwrap()
                    .duration_s,
                DT
            );
        }
    });
}

#[test]
fn g3_conduction_uses_final_composition_and_actual_solid_source_conditions() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let mut gas = gas(&reaction);
        let old_gas = gas.clone();
        let (mut live, card, point, key) = quiet_body(cx, gas.state().state(), 300.0);
        let radius = live.mechanical().specimen().radius_m();
        let area = 2.0
            * core::f64::consts::PI
            * radius
            * (live.mechanical().specimen().string().length_m + radius);
        let layers = MaterialConductionLayers {
            material: &card,
            context: &point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
            conductivity: &key,
            selection: MaterialPropertySelection::SingleClaimOnly,
            area_m2: area,
            gas_distance_m: 0.001,
            body_distance_m: radius / 2.0,
        };
        let result = live
            .step_with_reacting_isochoric_gas_conduction(
                cx,
                0,
                &[0.0],
                &mut gas,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
                exchange_solve(),
                &layers,
            )
            .unwrap();
        let frame = &result.coupled;
        let ks = 20.0 + 0.5 * (live.thermal().temperature_k() - 300.0);
        let kg = gas.state().state().thermal_conductivity;
        let expected_g = area / (layers.gas_distance_m / kg + layers.body_distance_m / ks);
        near(frame.conductance_w_k, expected_g, 2e-13);
        near(
            frame.source.property_by_key(&key).unwrap().value_si(),
            ks,
            1e-11,
        );
        assert_eq!(frame.source.card_identity(), card.content_hash());
        assert!(
            frame
                .source
                .query_point()
                .contains(&("T".into(), live.thermal().temperature_k()))
        );
        assert!(frame.source.query_point().contains(&("omega".into(), 1.0)));
        near(
            frame.solid_to_gas_heat_j,
            DT * expected_g * (live.thermal().temperature_k() - gas.state().state().temperature),
            exchange_solve().source_tolerance_j,
        );
        // Isolate composition from temperature: frozen reactants heated to the
        // SAME endpoint cannot supply the mixture's actual conductivity.
        let mut frozen = old_gas.clone();
        let cg = old_gas.mass_kg() / 0.028 * 2.5 * R;
        frozen
            .add_heat(
                cg * (gas.state().state().temperature - old_gas.state().state().temperature),
                solve().gas,
                || false,
            )
            .unwrap();
        assert!(kg > frozen.state().state().thermal_conductivity * 1.01);
        assert!(live.thermal().temperature_k() > 300.0);
        assert_eq!(live.mechanical().ambient(), gas.state().state());
        assert!(frame.coupled.energy_balance_residual_j.abs() <= frame.coupled.energy_tolerance_j);
    });
}

#[test]
fn g3_insulated_exchange_is_the_same_native_reacting_sound_trajectory() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let (mut direct, mut direct_gas) = prepared(cx, &reaction);
        let mut exchange = direct.clone();
        let mut exchange_gas = direct_gas.clone();
        for epoch in 8..12 {
            let heat = FiniteGasStringHeat {
                string_heat_j: 5e-8,
                gas_heat_j: 0.01,
            };
            let native = direct
                .step_with_reacting_isochoric_gas(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut direct_gas,
                    heat,
                    &reaction,
                    controls(),
                )
                .unwrap();
            let isolated = exchange
                .step_with_reacting_isochoric_gas_exchange(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut exchange_gas,
                    heat,
                    &reaction,
                    controls(),
                    exchange_solve(),
                    |_, _| {
                        Ok::<_, &'static str>(ConductanceSample {
                            conductance_w_k: 0.0,
                            source: (),
                        })
                    },
                )
                .unwrap();
            assert_eq!(isolated.coupled.solid_to_gas_heat_j, 0.0);
            assert_eq!(isolated.coupled.trials, 1);
            assert_eq!(isolated.chemistry, native.chemistry);
            assert_eq!(exchange_gas, direct_gas);
            assert_eq!(exchange.thermal(), direct.thermal());
            assert_eq!(exchange.mechanical().states(), direct.mechanical().states());
            assert_eq!(
                isolated.coupled.coupled.string.vibration,
                native.coupled.string.vibration
            );
            assert_eq!(exchange.mechanical().epoch(), epoch + 1);
            assert_eq!(exchange.mechanical().accepted_samples(), epoch + 1);
        }
    });
}

#[test]
fn g4_late_reacting_exchange_refusal_and_cancellation_preserve_retry_state() {
    with_string_cx(|cx, gate| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let old_gas = gas(&reaction);
        let (before, _, _, _) = quiet_body(cx, old_gas.state().state(), 300.0);
        let mut live = before.clone();
        let mut gas = old_gas.clone();
        let mut calls = 0;
        let tiny = ReactiveFiniteGasStringSolve {
            maximum_gas_temperature_increment_k: 1e-8,
            ..controls()
        };
        assert!(
            live.step_with_reacting_isochoric_gas_exchange(
                cx,
                0,
                &[0.0],
                &mut gas,
                FiniteGasStringHeat::default(),
                &reaction,
                tiny,
                exchange_solve(),
                |gas, solid| {
                    calls += 1;
                    conductance(gas, solid)
                },
            )
            .is_err()
        );
        assert!(calls > 1, "failure must follow actual coupled proposals");
        assert_eq!(gas, old_gas);
        assert_eq!(live.thermal(), before.thermal());
        assert_eq!(live.mechanical().states(), before.mechanical().states());
        assert_eq!(live.mechanical().ambient(), before.mechanical().ambient());
        assert_eq!(live.mechanical().epoch(), 0);
        assert_eq!(live.mechanical().accepted_samples(), 0);
        let mut reference = before.clone();
        let mut reference_gas = old_gas.clone();
        let retry = live
            .step_with_reacting_isochoric_gas_exchange(
                cx,
                0,
                &[0.0],
                &mut gas,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
                exchange_solve(),
                conductance,
            )
            .unwrap();
        let expected = reference
            .step_with_reacting_isochoric_gas_exchange(
                cx,
                0,
                &[0.0],
                &mut reference_gas,
                FiniteGasStringHeat::default(),
                &reaction,
                controls(),
                exchange_solve(),
                conductance,
            )
            .unwrap();
        assert_eq!(retry.chemistry, expected.chemistry);
        assert_eq!(
            retry.coupled.coupled.string,
            expected.coupled.coupled.string
        );
        assert_eq!(gas, reference_gas);
        // Cancellation arrives from the source callback AFTER a native
        // chemical proposal, before the shared owner can publish anything.
        let mut cancelled = before.clone();
        let mut cancelled_gas = old_gas.clone();
        let mut proposals = 0;
        assert!(
            cancelled
                .step_with_reacting_isochoric_gas_exchange(
                    cx,
                    0,
                    &[0.0],
                    &mut cancelled_gas,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    controls(),
                    exchange_solve(),
                    |gas, solid| {
                        proposals += 1;
                        gate.request();
                        conductance(gas, solid)
                    },
                )
                .is_err()
        );
        assert_eq!(proposals, 1);
        assert_eq!(cancelled_gas, old_gas);
        assert_eq!(cancelled.thermal(), before.thermal());
        assert_eq!(
            cancelled.mechanical().states(),
            before.mechanical().states()
        );
        assert_eq!(cancelled.mechanical().epoch(), 0);
        assert_eq!(cancelled.mechanical().accepted_samples(), 0);
    });
}
