//! G1 finite-gas Stefan/heat event balances and G3 common-clock owner parity.
use super::*;
use crate::evaporation::gas::{
    equilibrium::{
        heat::advance_equilibrium_surface_with_heat,
        tests::{CP, FUSION, MELTING, TRANSFER, VAPORIZATION, bulk, catalog, curve, near},
    },
    tests::{Phase, binding, controls, gas},
};
use crate::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::evaporation::EvaporationProperties;

const EVENT: FreezingHeatSolve = FreezingHeatSolve {
    maximum_iterations: 128,
};

fn heat_controls(chart: &EquilibriumEnthalpyPhaseCurve) -> EquilibriumGasHeatSolve {
    EquilibriumGasHeatSolve {
        bulk_enthalpy_range_j_kg: [
            chart.knots()[0].specific_enthalpy_j_kg,
            chart.knots().last().unwrap().specific_enthalpy_j_kg,
        ],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}

fn phase_controls() -> GasEvaporationSolve {
    let mut result = controls();
    result.maximum_inventory_fraction = 0.5;
    result
}

#[test]
fn evaporative_freezing_replays_both_owners_at_one_common_duration() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, 280.0, 0.0, 1.0);
    let maximum = 20_000.0;
    let result = advance_to_freezing_with_heat(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &original_gas,
        &source,
        TRANSFER,
        1.0,
        maximum,
        phase_controls(),
        heat_controls(&chart),
        EVENT,
        || false,
    )
    .unwrap();
    let time = result.boundary().elapsed_duration_s;
    assert!(result.boundary().freezing_completed && time < maximum);
    assert!(result.trials() > 1 && result.trials() <= EVENT.maximum_iterations);
    assert_eq!(result.step().duration_s, time);
    assert_eq!(
        result.transfer().unwrap().exchange().step().duration_s,
        time
    );
    assert_eq!(
        result.transfer().unwrap().boundary().remaining_duration_s,
        0.0
    );
    near(
        time + result.boundary().remaining_duration_s,
        maximum,
        1e-11,
    );

    // Independent constant-Cv gas, plateau heat and Stefan-film event oracle.
    // At complete freezing: m_liquid,0*Lf + Q = dm*(Lf+Lv).
    let capacity = 100_000.0 / 280.0 * 2.5;
    let heat =
        |t: f64| initial.area_m2 * t * (280.0 - MELTING) / (1.0 + initial.area_m2 * t / capacity);
    let mass = |t: f64| {
        let gas_t = 280.0 - heat(t) / capacity;
        let pressure = 100_000.0 * gas_t / 280.0;
        initial.area_m2 * TRANSFER * pressure / ((CP / 4.0) * gas_t)
            * (pressure / (pressure - 1500.0)).ln()
            * t
    };
    let mut lo = 0.0;
    let mut hi = maximum;
    for _ in 0..100 {
        let t = lo + 0.5 * (hi - lo);
        if 0.005 * FUSION + heat(t) - mass(t) * (FUSION + VAPORIZATION) > 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    near(time, lo + 0.5 * (hi - lo), 2e-4);
    near(result.heat().heat_into_liquid_j, heat(time), 2e-6);
    near(result.step().transferred_mass_kg, mass(time), 2e-11);
    assert!((result.heat().heat_into_liquid_j - heat(maximum)).abs() > 100.0);
    assert_eq!(
        result.bulk().specific_enthalpy_j_kg,
        inventory.solid_at_melting_h_j_kg
    );
    assert_eq!(result.phase().final_state.liquid_mass_kg, 0.0);
    near(
        result.phase().fusion_energy_j,
        result.step().transferred_mass_kg * VAPORIZATION - result.heat().heat_into_liquid_j,
        1e-8,
    );
    near(
        result.step().liquid_enthalpy_change_j + result.step().gas_internal_energy_change_j,
        0.0,
        1e-8,
    );
    assert!(result.sum_absolute_energy_residual_j() <= phase_controls().gas.energy_tolerance_j);
    assert!(std::ptr::eq(result.catalog(), &source));

    // The exact published heat receipt can be replayed at the SAME duration
    // through the original sensible owner and its admitted terminal bound.
    let mut native = heat_controls(&chart)
        .original_controls(inventory, initial.specific_enthalpy_j_kg)
        .unwrap();
    native.pure_liquid_enthalpy_range_j_kg[0] = inventory.solid_at_melting_h_j_kg;
    let (heated, heated_gas, replayed_heat) = sensible_heat(
        initial,
        &original_gas,
        MELTING,
        &|h| {
            chart
                .temperature_at_specific_enthalpy(h)
                .map_err(|_| EvaporationError::Invalid("test source"))
        },
        1.0,
        time,
        phase_controls(),
        native,
        &mut || false,
    )
    .unwrap();
    let replayed_phase = advance_to_freezing(
        inventory,
        &Phase,
        &binding(0.0),
        heated,
        &heated_gas,
        &source,
        TRANSFER,
        time,
        phase_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(&replayed_heat, result.heat());
    assert_eq!(&replayed_phase, result.transfer().unwrap());
    let FreezingHeatTrial::Accepted(single) = single_trial_to_freezing_with_heat(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &original_gas,
        &source,
        TRANSFER,
        1.0,
        time,
        phase_controls(),
        heat_controls(&chart),
        || false,
    )
    .unwrap() else {
        panic!("the common-duration event must remain an accepted single trial");
    };
    assert_eq!(single.bulk(), result.bulk());
    assert_eq!(single.gas(), result.gas());
    assert_eq!(single.heat(), result.heat());
    assert_eq!(single.step(), result.step());
    assert_eq!(single.boundary().remaining_duration_s, 0.0);
    assert_eq!(single.trials(), 1);
}

#[test]
fn the_published_heat_is_continuous_across_the_ulp_where_the_bisection_path_flipped() {
    // At the event of the test above the heat solve of two ADJACENT durations,
    // 5554.513257088766 s (0x40b5b28364d10ab1) and the next double, took 31 and
    // 34 bisection steps and published heats that differ by about 2.8e-8 J: the
    // first midpoint inside the 5e-9 K band was a different point of the band.
    // That jump moves the phase owner's mass by 9.7e-14 kg, nine times its
    // 1e-14 kg completion window, which is 2.9e-9 J of heat. The published heat
    // is the root to the secant's accuracy, so it follows the true rate
    // (about 0.07 W, 6e-14 J per ulp of the duration) and the two agree to
    // far less than that window.
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, 280.0, 0.0, 1.0);
    let mut native = heat_controls(&chart)
        .original_controls(inventory, initial.specific_enthalpy_j_kg)
        .unwrap();
    native.pure_liquid_enthalpy_range_j_kg[0] = inventory.solid_at_melting_h_j_kg;
    let heat_at = |duration: f64| {
        sensible_heat(
            initial,
            &original_gas,
            MELTING,
            &|h| {
                chart
                    .temperature_at_specific_enthalpy(h)
                    .map_err(|_| EvaporationError::Invalid("test source"))
            },
            1.0,
            duration,
            phase_controls(),
            native,
            &mut || false,
        )
        .unwrap()
        .2
        .heat_into_liquid_j
    };
    let before = f64::from_bits(0x40b5_b283_64d1_0ab1);
    let after = before.next_up();
    assert_eq!(before, 5554.513257088766);
    let jump = (heat_at(after) - heat_at(before)).abs();
    assert!(jump < 1e-9, "{jump:e} J across one ulp of the duration");
}

#[test]
fn reciprocal_heat_alone_completes_isothermal_freezing_conservatively() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, 260.0, 0.0, 1.0);
    let result = advance_to_freezing_with_heat(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &original_gas,
        &source,
        0.0,
        60.0,
        500.0,
        phase_controls(),
        heat_controls(&chart),
        EVENT,
        || false,
    )
    .unwrap();
    let expected_heat = -0.005 * FUSION;
    let capacity = 100_000.0 / 260.0 * 2.5;
    let expected_gas_t = 260.0 - expected_heat / capacity;
    let expected_time = expected_heat / (60.0 * initial.area_m2 * (expected_gas_t - MELTING));
    assert!(result.boundary().freezing_completed);
    assert!(result.transfer().is_none());
    assert_eq!(result.bulk().mass_kg, initial.mass_kg);
    assert_eq!(result.gas().mass_kg(), original_gas.mass_kg());
    assert_eq!(
        result.bulk().specific_enthalpy_j_kg,
        inventory.solid_at_melting_h_j_kg
    );
    assert_eq!(result.phase().final_state.liquid_mass_kg, 0.0);
    assert_eq!(result.step().surface_temperature_k, [MELTING, MELTING]);
    assert_eq!(result.step().transferred_mass_kg, 0.0);
    assert_eq!(result.step().latent_energy_j, 0.0);
    near(result.heat().heat_into_liquid_j, expected_heat, 1e-10);
    near(result.boundary().elapsed_duration_s, expected_time, 1e-6);
    near(result.step().gas_temperature_k[1], expected_gas_t, 1e-9);
    near(result.phase().fusion_energy_j, -expected_heat, 1e-10);
    near(
        result.step().liquid_enthalpy_change_j + result.step().gas_internal_energy_change_j,
        0.0,
        1e-8,
    );
    assert!(result.sum_absolute_energy_residual_j() <= phase_controls().gas.energy_tolerance_j);
    assert!(inventory.state(result.bulk()).is_err());
}

#[test]
fn exact_saturation_allows_thermal_completion_with_nonzero_vapor_coefficient() {
    struct Saturated {
        pressure: f64,
    }
    impl EvaporationLaw for Saturated {
        fn temperature_range_k(&self) -> [f64; 2] {
            Phase.temperature_range_k()
        }
        fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
            let mut phase = Phase.properties(t)?;
            phase.saturation_pressure_pa =
                self.pressure * (VAPORIZATION / (CP / 4.0) * (1.0 / MELTING - 1.0 / t)).exp();
            phase.saturation_pressure_derivative_pa_k =
                phase.saturation_pressure_pa * VAPORIZATION / ((CP / 4.0) * t * t);
            Ok(phase)
        }
    }
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, 260.0, 1000.0, 1.0);
    let baseline = advance_to_freezing_with_heat(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &original_gas,
        &source,
        0.0,
        60.0,
        500.0,
        phase_controls(),
        heat_controls(&chart),
        EVENT,
        || false,
    )
    .unwrap();
    let phase = Saturated {
        pressure: baseline
            .gas()
            .species_partial_pressure_pa(&binding(0.0).species)
            .unwrap(),
    };
    let result = advance_to_freezing_with_heat(
        inventory,
        &phase,
        &binding(0.0),
        initial,
        &original_gas,
        &source,
        TRANSFER,
        60.0,
        500.0,
        phase_controls(),
        heat_controls(&chart),
        EVENT,
        || false,
    )
    .unwrap();
    assert!(result.boundary().freezing_completed);
    assert!(result.transfer().is_none());
    assert_eq!(result.step().transferred_mass_kg, 0.0);
    assert_eq!(result.bulk(), baseline.bulk());
    assert_eq!(result.gas(), baseline.gas());
    near(
        result.boundary().elapsed_duration_s,
        baseline.boundary().elapsed_duration_s,
        1e-6,
    );
}

#[test]
fn zero_heat_keeps_exact_phase_event_parity_and_ignores_unused_controls() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, MELTING, 0.0, 1.0);
    let unused = EquilibriumGasHeatSolve {
        bulk_enthalpy_range_j_kg: [f64::NAN; 2],
        temperature_tolerance_k: f64::NAN,
        maximum_iterations: 0,
    };
    for duration in [0.0, 1.0, 20_000.0] {
        let phase = advance_to_freezing(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &original_gas,
            &source,
            TRANSFER,
            duration,
            phase_controls(),
            || false,
        )
        .unwrap();
        let result = advance_to_freezing_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &original_gas,
            &source,
            TRANSFER,
            0.0,
            duration,
            phase_controls(),
            unused,
            FreezingHeatSolve {
                maximum_iterations: 0,
            },
            || false,
        )
        .unwrap();
        assert_eq!(result.transfer(), Some(&phase));
        assert_eq!(result.step(), phase.exchange().step());
        assert_eq!(result.phase(), phase.phase());
        assert_eq!(result.boundary(), phase.boundary());
        assert_eq!(result.heat(), &GasFilmHeatStep::default());
        assert_eq!(result.trials(), 0);
        let single = single_trial_to_freezing_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &original_gas,
            &source,
            TRANSFER,
            0.0,
            duration,
            phase_controls(),
            unused,
            || false,
        )
        .unwrap();
        if phase.boundary().elapsed_duration_s < duration {
            assert!(
                matches!(single, FreezingHeatTrial::Bracket { residual_energy_j }
                if residual_energy_j < 0.0)
            );
        } else {
            let FreezingHeatTrial::Accepted(single) = single else {
                panic!("a full-duration phase endpoint must remain available");
            };
            assert_eq!(single.transfer(), Some(&phase));
        }
    }
}

#[test]
fn ordinary_partial_steps_retain_the_original_heat_and_phase_owners() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    for (fraction, t, pressure, duration) in [
        (0.5, 330.0, 0.0, 1.0),
        (0.5, 280.0, 5000.0, 0.5),
        (1.0, 330.0, 5000.0, 0.5),
        (0.5, 280.0, 0.0, 0.0),
    ] {
        let initial = bulk(&chart, fraction);
        let original_gas = gas(&source, t, pressure, 0.01);
        let old = advance_equilibrium_surface_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &original_gas,
            &source,
            TRANSFER,
            30.0,
            duration,
            phase_controls(),
            heat_controls(&chart),
            || false,
        )
        .unwrap();
        let new = advance_to_freezing_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &original_gas,
            &source,
            TRANSFER,
            30.0,
            duration,
            phase_controls(),
            heat_controls(&chart),
            EVENT,
            || false,
        )
        .unwrap();
        assert!(!new.boundary().freezing_completed);
        assert_eq!(new.boundary().elapsed_duration_s, duration);
        assert_eq!(new.boundary().remaining_duration_s, 0.0);
        assert_eq!(new.bulk(), old.bulk());
        assert_eq!(new.gas(), old.gas());
        assert_eq!(new.heat(), old.heat());
        assert_eq!(new.step(), old.step());
        assert_eq!(new.phase(), old.phase());
    }
}

#[test]
fn original_donor_temperature_cancellation_and_work_bounds_remain_atomic() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, 280.0, 0.0, 1.0);
    let saved = original_gas.clone();
    let run = |phase, event, cancelled| {
        advance_to_freezing_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &original_gas,
            &source,
            TRANSFER,
            1.0,
            20_000.0,
            phase,
            heat_controls(&chart),
            event,
            cancelled,
        )
    };
    let never: fn() -> bool = || false;
    assert!(matches!(
        run(
            phase_controls(),
            FreezingHeatSolve {
                maximum_iterations: 1
            },
            never
        ),
        Err(GasEvaporationError::StepSizeLimit {
            reason: "common-duration heat/phase freezing event did not converge",
            ..
        })
    ));
    assert!(matches!(
        run(controls(), EVENT, never),
        Err(GasEvaporationError::StepSizeLimit {
            reason: "initial liquid or vapor inventory withdrawal limit",
            ..
        })
    ));
    let mut tight = phase_controls();
    tight.maximum_temperature_increment_k = 1e-7;
    assert!(matches!(
        run(tight, EVENT, never),
        Err(GasEvaporationError::StepSizeLimit { .. })
    ));
    assert!(matches!(
        run(
            phase_controls(),
            FreezingHeatSolve {
                maximum_iterations: 0
            },
            never
        ),
        Err(GasEvaporationError::Invalid(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        advance_to_freezing_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &original_gas,
            &source,
            TRANSFER,
            1.0,
            20_000.0,
            phase_controls(),
            heat_controls(&chart),
            EVENT,
            || {
                polls += 1;
                polls >= 40
            },
        ),
        Err(GasEvaporationError::Cancelled)
    ));
    assert!(polls >= 40);
    assert_eq!(original_gas, saved);
    assert_eq!(initial, bulk(&chart, 0.5));
}
