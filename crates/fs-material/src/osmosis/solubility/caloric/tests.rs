use super::*;
use fs_qty::AmountConcentration;
use fs_thermochem::{
    MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    solution::{
        caloric::{MolarHeatCapacityCurveV1, PureLiquidHeatCapacityV1},
        solubility::{IdealSolubilityV1, caloric::SolubilityCaloricComponentV1},
    },
};

fn law(solvent: &str, solute: &str, saturation: f64, delta_h: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new(solvent).unwrap(),
        SpeciesId::new(solute).unwrap(),
        Temperature::new(290.0),
        AmountConcentration::new(saturation),
        delta_h,
        [Temperature::new(280.0), Temperature::new(340.0)],
    )
    .unwrap()
}

fn thermal_model(molar_volume: f64, laws: &[IdealSolubilityV1]) -> SolubilityCaloricModelV1 {
    let knots = [280.0, 340.0].map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(50.0)));
    let solvent = PureLiquidHeatCapacityV1::new(
        laws[0].solvent().clone(),
        Temperature::new(300.0),
        knots.to_vec(),
    )
    .unwrap();
    let components = laws
        .iter()
        .map(|law| {
            if law.phase_heat_capacities().is_some() {
                return SolubilityCaloricComponentV1::from_phase_heat_capacities(law.clone())
                    .unwrap();
            }
            SolubilityCaloricComponentV1::new(
                law.clone(),
                MolarHeatCapacityCurveV1::new(
                    law.solute().clone(),
                    Temperature::new(300.0),
                    knots.to_vec(),
                )
                .unwrap(),
            )
            .unwrap()
        })
        .collect();
    SolubilityCaloricModelV1::new(molar_volume, solvent, components).unwrap()
}

fn controls(low: f64, high: f64) -> SolubilityCaloricSolveV1 {
    SolubilityCaloricSolveV1 {
        temperature_bracket: [Temperature::new(low), Temperature::new(high)],
        temperature_tolerance_k: 1e-10,
        energy_tolerance_j: 1e-8,
        maximum_iterations: 128,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}

fn consistent(initial: &FiniteSolubilitySolution, step: &FiniteSolubilityHeatStep) {
    let material = step.solution();
    let thermal = step.balance().final_state();
    assert_eq!(
        material.solution().temperature_k(),
        thermal.temperature().value()
    );
    assert_eq!(
        material.solution().solvent_volume_m3(),
        initial.solution().solvent_volume_m3(),
    );
    assert_eq!(
        thermal.solvent_volume().value(),
        initial.solution().solvent_volume_m3(),
    );
    assert_eq!(material.phase_partitions(), thermal.equilibrium().unwrap(),);
    for ((before, after), total) in initial
        .components()
        .iter()
        .zip(material.components())
        .zip(thermal.total_solute())
    {
        assert_eq!(before.law(), after.law());
        assert_eq!(before.total(), after.total());
        assert_eq!(after.total(), *total);
    }
    for ((amount, phase), aqueous_moles) in thermal
        .dissolved_solute()
        .iter()
        .zip(material.phase_partitions())
        .zip(material.solution().solute_moles())
    {
        assert_eq!(*amount, phase.dissolved);
        assert_eq!(amount.value(), *aqueous_moles);
        near(
            phase.dissolved.value() + phase.solid.value(),
            phase.total_solute.value(),
            8.0 * f64::EPSILON * phase.total_solute.value(),
        );
    }
    near(step.balance().energy_residual_j(), 0.0, 1e-8);
}

#[test]
fn unequal_phase_heat_capacities_follow_live_material_stocks_and_osmotic_activity() {
    let original_law = law("solvent", "A", 250.0, 1000.0);
    let phase_cp = |capacity| {
        MolarHeatCapacityCurveV1::new(
            original_law.solute().clone(),
            Temperature::new(300.0),
            [280.0, 340.0]
                .into_iter()
                .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(capacity)))
                .collect(),
        )
        .unwrap()
    };
    let source = original_law
        .clone()
        .with_phase_heat_capacities(phase_cp(30.0), phase_cp(90.0))
        .unwrap();
    let model = thermal_model(2e-5, std::slice::from_ref(&source));
    let original = FiniteSolubilitySolution::new(
        Temperature::new(290.0),
        2e-5,
        Volume::new(0.0002),
        vec![(source.clone(), Amount::new(0.2))],
    )
    .unwrap();
    let mut transported = original.independent_clone();
    transported.replace_inventories(0.00011, &[0.06]).unwrap();
    let before = transported.clone();
    // Independent Kirchhoff/van't Hoff closed form; both endpoints retain solid.
    let target = 320.0;
    let saturation = 250.0
        * ((1000.0 * (1.0 / 290.0 - 1.0 / target)
            + 60.0 * ((target / 290.0_f64).ln() + 290.0 / target - 1.0))
            / R)
            .exp();
    let final_dissolved = saturation * 0.00011;
    let initial_dissolved = 250.0 * 0.00011;
    let heat =
        (0.00011 / 2e-5 * 50.0 + (0.06 - initial_dissolved) * 30.0 + initial_dissolved * 90.0)
            * 30.0
            + (final_dissolved - initial_dissolved) * (1000.0 + 60.0 * 30.0);
    let step = transported
        .add_heat(&model, heat, controls(310.0, 330.0), || false)
        .unwrap();
    consistent(&transported, &step);
    near(step.solution().solution().temperature_k(), target, 1e-9);
    near(
        step.solution().phase_partitions()[0].dissolved.value(),
        final_dissolved,
        2e-12,
    );
    near(
        step.solution().solution().osmotic_pressure_pa().unwrap(),
        R * target * saturation,
        1e-7,
    );
    near(
        step.solution().solution().water_activity().unwrap(),
        (-2e-5 * saturation).exp(),
        1e-14,
    );
    assert_eq!(step.solution().components()[0].law(), &source);
    assert_eq!(transported, before);

    // A numerically identical reference saturation without the original phase
    // Cp sources cannot stand in for this material's actual thermal owner.
    let wrong_model = thermal_model(2e-5, &[original_law]);
    assert_eq!(
        transported.add_heat(&wrong_model, heat, controls(310.0, 330.0), || false),
        Err(FiniteSolubilityHeatError::SourceMismatch)
    );
}

#[test]
fn heating_exhausts_crystals_and_reverse_cooling_restores_phases_with_one_heat_balance() {
    let law = law("solvent", "A", 250.0, 30_000.0);
    let model = thermal_model(2e-5, std::slice::from_ref(&law));
    let initial = FiniteSolubilitySolution::new(
        Temperature::new(290.0),
        2e-5,
        Volume::new(0.0002),
        vec![(law, Amount::new(0.1))],
    )
    .unwrap();
    let before = initial.clone();
    near(initial.phase_partitions()[0].dissolved.value(), 0.05, 1e-15);
    // Independent constant-Cp balance: (10+0.1)*50*30 J sensible heat,
    // plus 0.05*30000 J to dissolve the remaining solid. At 320 K the
    // original van't Hoff capacity exceeds the complete 0.1 mol stock.
    let heat = 505.0 * 30.0 + 0.05 * 30_000.0;
    let heated = initial
        .add_heat(&model, heat, controls(315.0, 325.0), || false)
        .unwrap();
    consistent(&initial, &heated);
    near(heated.solution().solution().temperature_k(), 320.0, 1e-9);
    assert_eq!(
        heated.solution().phase_partitions()[0].dissolved.value(),
        0.1
    );
    assert_eq!(heated.solution().phase_partitions()[0].solid.value(), 0.0);
    near(heated.balance().enthalpy_change_j(), heat, 1e-8);
    near(
        heated.solution().solution().osmotic_pressure_pa().unwrap(),
        R * 320.0 * 500.0,
        1e-7,
    );
    near(
        heated.solution().solution().water_activity().unwrap(),
        (-0.01_f64).exp(),
        2e-15,
    );
    let (hot_material, balance) = heated.into_parts();
    assert_eq!(balance.final_state().dissolved_solute()[0].value(), 0.1);
    let cooled = hot_material
        .add_heat(&model, -heat, controls(285.0, 295.0), || false)
        .unwrap();
    consistent(&hot_material, &cooled);
    near(cooled.solution().solution().temperature_k(), 290.0, 1e-9);
    near(
        cooled.solution().phase_partitions()[0].dissolved.value(),
        0.05,
        1e-12,
    );
    near(
        cooled.solution().phase_partitions()[0].solid.value(),
        0.05,
        1e-12,
    );
    near(cooled.balance().enthalpy_change_j(), -heat, 1e-8);
    assert_eq!(initial, before);
}

#[test]
fn thermal_step_uses_transported_totals_and_preserves_the_exact_original_volume() {
    let law = law("solvent", "A", 250.0, 30_000.0);
    let molar_volume = 1.8e-5;
    let volume = 0.00011;
    assert_ne!((volume / molar_volume) * molar_volume, volume);
    let model = thermal_model(molar_volume, std::slice::from_ref(&law));
    let original = FiniteSolubilitySolution::new(
        Temperature::new(290.0),
        molar_volume,
        Volume::new(0.0002),
        vec![(law, Amount::new(0.2))],
    )
    .unwrap();
    let mut transported = original.independent_clone();
    transported.replace_inventories(volume, &[0.06]).unwrap();
    let before = transported.clone();
    // The CURRENT bath has V/Vbar solvent moles and 0.06 mol total solute.
    // At 290 K its dissolved capacity is 0.00011*250=0.0275 mol; heating
    // to 320 K exhausts the 0.0325 mol solid. No construction-time stock
    // or rounded quotient/product volume may enter this balance.
    let heat = (volume / molar_volume + 0.06) * 50.0 * 30.0 + (0.06 - 0.0275) * 30_000.0;
    let heated = transported
        .add_heat(&model, heat, controls(315.0, 325.0), || false)
        .unwrap();
    consistent(&transported, &heated);
    near(heated.solution().solution().temperature_k(), 320.0, 1e-9);
    assert_eq!(heated.solution().components()[0].total().value(), 0.06);
    assert_eq!(
        heated.solution().phase_partitions()[0].dissolved.value(),
        0.06
    );
    assert_eq!(heated.solution().phase_partitions()[0].solid.value(), 0.0);
    assert_eq!(heated.solution().solution().solvent_volume_m3(), volume);
    assert_eq!(
        heated.balance().final_state().solvent_volume().value(),
        volume
    );
    assert_eq!(
        heated.balance().final_state().solvent_moles().value(),
        volume / molar_volume,
    );
    near(heated.balance().enthalpy_change_j(), heat, 1e-8);
    assert_eq!(transported, before);
    assert_eq!(original.components()[0].total().value(), 0.2);
    assert_eq!(original.solution().solvent_volume_m3(), 0.0002);
}

#[test]
fn exact_ordered_source_binding_and_unsupported_heat_refuse_without_changing_material() {
    let a = law("solvent", "A", 250.0, 30_000.0);
    let b = law("solvent", "B", 100.0, 0.0);
    let initial = FiniteSolubilitySolution::new(
        Temperature::new(290.0),
        2e-5,
        Volume::new(0.0002),
        vec![
            (a.clone(), Amount::new(0.1)),
            (b.clone(), Amount::new(0.05)),
        ],
    )
    .unwrap();
    let before = initial.clone();
    let matching = thermal_model(2e-5, &[a.clone(), b.clone()]);
    let incompatible = [
        thermal_model(2e-5, &[b.clone(), a.clone()]),
        thermal_model(2e-5, &[law("solvent", "A", 251.0, 30_000.0), b.clone()]),
        thermal_model(1.8e-5, &[a, b]),
        thermal_model(
            2e-5,
            &[
                law("other-solvent", "A", 250.0, 30_000.0),
                law("other-solvent", "B", 100.0, 0.0),
            ],
        ),
    ];
    for model in incompatible {
        assert_eq!(
            initial.add_heat(&model, 0.0, controls(285.0, 325.0), || false),
            Err(FiniteSolubilityHeatError::SourceMismatch),
        );
        assert_eq!(initial, before);
    }
    assert!(matches!(
        initial.add_heat(&matching, 1e12, controls(285.0, 325.0), || false),
        Err(FiniteSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::HeatOutsideBracket { .. }
        )),
    ));
    assert!(matches!(
        initial.add_heat(&matching, 0.0, controls(279.0, 325.0), || false),
        Err(FiniteSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::Invalid(_)
        )),
    ));
    assert_eq!(initial, before);
}

#[test]
fn early_and_final_publication_cancellation_preserve_all_material_views() {
    let law = law("solvent", "A", 1000.0, 0.0);
    let model = thermal_model(2e-5, std::slice::from_ref(&law));
    let initial = FiniteSolubilitySolution::new(
        Temperature::new(300.0),
        2e-5,
        Volume::new(0.0002),
        vec![(law, Amount::new(0.1))],
    )
    .unwrap();
    let before = initial.clone();
    let heat = 505.0 * 1.3;
    let budget = controls(290.0, 310.0);
    let mut calls = 0;
    let accepted = initial
        .add_heat(&model, heat, budget, || {
            calls += 1;
            false
        })
        .unwrap();
    consistent(&initial, &accepted);
    near(accepted.solution().solution().temperature_k(), 301.3, 1e-9);
    assert!(calls > 2);
    let mut early_calls = 0;
    assert_eq!(
        initial.add_heat(&model, heat, budget, || {
            early_calls += 1;
            true
        }),
        Err(FiniteSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::Cancelled
        )),
    );
    assert_eq!(early_calls, 1);
    let mut late_calls = 0;
    assert_eq!(
        initial.add_heat(&model, heat, budget, || {
            late_calls += 1;
            late_calls == calls
        }),
        Err(FiniteSolubilityHeatError::Caloric(
            SolubilityCaloricErrorV1::Cancelled
        )),
    );
    assert_eq!(late_calls, calls);
    assert_eq!(initial, before);
    assert_eq!(
        initial.add_heat(&model, heat, budget, || false).unwrap(),
        accepted,
    );
}
