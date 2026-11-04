use super::*;
use crate::{MolarThermalQuantityV1, SpeciesId};
use fs_qty::AmountConcentration;

fn species(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn common_cp(name: &str, capacity: f64) -> MolarHeatCapacityCurveV1 {
    MolarHeatCapacityCurveV1::new(
        species(name),
        Temperature::new(300.0),
        [280.0, 340.0]
            .into_iter()
            .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(capacity)))
            .collect(),
    )
    .unwrap()
}
fn solvent_cp(capacity: f64) -> PureLiquidHeatCapacityV1 {
    PureLiquidHeatCapacityV1::new(
        species("solvent"),
        Temperature::new(300.0),
        [280.0, 340.0]
            .into_iter()
            .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(capacity)))
            .collect(),
    )
    .unwrap()
}
fn component(
    name: &str,
    reference_temperature: f64,
    concentration: f64,
    delta_h: f64,
    capacity: f64,
) -> SolubilityCaloricComponentV1 {
    SolubilityCaloricComponentV1::new(
        IdealSolubilityV1::new(
            species("solvent"),
            species(name),
            Temperature::new(reference_temperature),
            AmountConcentration::new(concentration),
            delta_h,
            [Temperature::new(280.0), Temperature::new(340.0)],
        )
        .unwrap(),
        common_cp(name, capacity),
    )
    .unwrap()
}
fn model(components: Vec<SolubilityCaloricComponentV1>) -> SolubilityCaloricModelV1 {
    SolubilityCaloricModelV1::new(2e-5, solvent_cp(50.0), components).unwrap()
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
fn conserved(initial: &SolubilityCaloricStateV1, final_state: &SolubilityCaloricStateV1) {
    assert_eq!(initial.solvent_moles(), final_state.solvent_moles());
    assert_eq!(initial.solvent_volume(), final_state.solvent_volume());
    assert_eq!(initial.total_solute(), final_state.total_solute());
    for partition in final_state.equilibrium().unwrap() {
        near(
            partition.dissolved.value() + partition.solid.value(),
            partition.total_solute.value(),
            4.0 * f64::EPSILON * partition.total_solute.value(),
        );
        assert!(partition.dissolved.value() >= 0.0);
        assert!(partition.solid.value() >= 0.0);
        if partition.solid.value() > 0.0 {
            assert_eq!(partition.dissolved, partition.capacity);
        } else {
            assert!(partition.dissolved.value() <= partition.capacity.value());
        }
    }
}

#[test]
fn g1_adiabatic_dissolution_and_precipitation_close_the_same_latent_energy() {
    // Ten solvent moles and two solute moles each carry 50 J/(mol K):
    // C_sensible=600 J/K, independently of the phase allocation.
    // At the target temperature V*c_sat=0.0002*500=0.1 mol.
    // The three exact balances are respectively -600+600, +600-600,
    // and +600-600 J. They cover endothermic/exothermic dissolution
    // and heat release on precipitation from an initial supersaturation.
    for (delta_h, target, initially_dissolved) in [
        (6000.0, 299.0, 0.0),
        (-6000.0, 301.0, 0.0),
        (6000.0, 301.0, 0.2),
    ] {
        let model = model(vec![component("A", target, 500.0, delta_h, 50.0)]);
        let initial = model
            .prepare(
                Temperature::new(300.0),
                Amount::new(10.0),
                &[Amount::new(2.0)],
                &[Amount::new(initially_dissolved)],
                || false,
            )
            .unwrap();
        let original = initial.clone();
        assert!(initial.equilibrium().is_none());
        let result = model
            .add_heat(&initial, 0.0, controls(target - 1.0, target + 1.0), || {
                false
            })
            .unwrap();
        let state = result.final_state();
        near(state.temperature().value(), target, 1e-9);
        near(state.dissolved_solute()[0].value(), 0.1, 1e-12);
        near(state.solid_solute().next().unwrap().value(), 1.9, 1e-12);
        near(result.enthalpy_change_j(), 0.0, 1e-8);
        near(result.energy_residual_j(), 0.0, 1e-8);
        conserved(&initial, state);
        assert_eq!(initial, original);
    }
}

#[test]
fn g1_heat_and_cooling_cross_complete_solid_exhaustion_without_losing_stock() {
    let model = model(vec![component("A", 290.0, 250.0, 30000.0, 50.0)]);
    let initial = model
        .equilibrate(
            Temperature::new(290.0),
            Amount::new(10.0),
            &[Amount::new(0.1)],
            || false,
        )
        .unwrap();
    near(initial.dissolved_solute()[0].value(), 0.05, 1e-15);
    // C_sensible=505 J/K, and 0.05 mol remaining solid dissolves.
    // At 320 K the van't Hoff capacity exceeds the full 0.1 mol stock.
    let heat = 505.0 * 30.0 + 30000.0 * 0.05;
    let heated = model
        .add_heat(&initial, heat, controls(315.0, 325.0), || false)
        .unwrap();
    near(heated.final_state().temperature().value(), 320.0, 1e-9);
    assert_eq!(heated.final_state().dissolved_solute()[0].value(), 0.1);
    assert_eq!(
        heated.final_state().solid_solute().next().unwrap().value(),
        0.0
    );
    near(heated.enthalpy_change_j(), heat, 1e-8);
    conserved(&initial, heated.final_state());

    let cooled = model
        .add_heat(heated.final_state(), -heat, controls(285.0, 295.0), || {
            false
        })
        .unwrap();
    near(cooled.final_state().temperature().value(), 290.0, 1e-9);
    near(
        cooled.final_state().dissolved_solute()[0].value(),
        0.05,
        1e-12,
    );
    conserved(&initial, cooled.final_state());
    assert_eq!(
        model
            .enthalpy_change_j(&initial, heated.final_state())
            .unwrap(),
        -model
            .enthalpy_change_j(heated.final_state(), &initial)
            .unwrap(),
    );
}

#[test]
fn g1_opposed_phase_transfers_share_one_bath_and_count_solvent_once() {
    let a = component("A", 299.0, 500.0, 6000.0, 50.0);
    let b = component("B", 299.0, 500.0, -5500.0, 72.5);
    // C=10*50 + 0.05*50 + 1*72.5 =575 J/K.
    // A exhausts its solid: 0.05*6000=300 J.
    // B precipitates: (0.1-0.15)*(-5500)=275 J.
    // The shared adiabatic bath therefore cools from300 to299 K.
    for (components, totals, dissolved, expected) in [
        (
            vec![a.clone(), b.clone()],
            [0.05, 1.0],
            [0.0, 0.15],
            [0.05, 0.1],
        ),
        (vec![b, a], [1.0, 0.05], [0.15, 0.0], [0.1, 0.05]),
    ] {
        let model = model(components);
        let initial = model
            .prepare(
                Temperature::new(300.0),
                Amount::new(10.0),
                &totals.map(Amount::new),
                &dissolved.map(Amount::new),
                || false,
            )
            .unwrap();
        let result = model
            .add_heat(&initial, 0.0, controls(298.0, 300.0), || false)
            .unwrap();
        near(result.final_state().temperature().value(), 299.0, 1e-9);
        for (actual, expected) in result.final_state().dissolved_solute().iter().zip(expected) {
            near(actual.value(), expected, 1e-12);
        }
        conserved(&initial, result.final_state());
    }
}

#[test]
fn g1_piecewise_cp_uses_original_integral_for_liquid_and_common_solute_phases() {
    let knots = [(280.0, 20.0), (300.0, 60.0), (320.0, 100.0)]
        .into_iter()
        .map(|(t, cp)| (Temperature::new(t), MolarThermalQuantityV1::new(cp)))
        .collect::<Vec<_>>();
    let liquid =
        PureLiquidHeatCapacityV1::new(species("solvent"), Temperature::new(300.0), knots.clone())
            .unwrap();
    let numeric =
        MolarHeatCapacityCurveV1::new(species("solvent"), Temperature::new(300.0), knots).unwrap();
    let from = Temperature::new(295.0);
    let to = Temperature::new(315.0);
    // Cp(T)=2T-540; mean Cp=70 over [295,315], crossing the300 K knot.
    near(
        numeric.enthalpy_change_j_per_mol(from, to).unwrap(),
        1400.0,
        1e-12,
    );
    assert_eq!(
        liquid.enthalpy_change_j_per_mol(from, to).unwrap(),
        numeric.enthalpy_change_j_per_mol(from, to).unwrap(),
    );
    assert_eq!(
        liquid.entropy_change_j_per_mol_k(from, to).unwrap(),
        numeric.entropy_change_j_per_mol_k(from, to).unwrap(),
    );
    let model =
        SolubilityCaloricModelV1::new(2e-5, liquid, vec![component("A", 300.0, 1000.0, 0.0, 80.0)])
            .unwrap();
    let initial = model
        .equilibrate(from, Amount::new(10.0), &[Amount::new(0.1)], || false)
        .unwrap();
    // Solvent: 10*1400 J; solute: 0.1*80*20 J, all dissolved throughout.
    let result = model
        .add_heat(&initial, 14160.0, controls(310.0, 320.0), || false)
        .unwrap();
    near(result.final_state().temperature().value(), 315.0, 1e-9);
    near(result.energy_residual_j(), 0.0, 1e-8);
    conserved(&initial, result.final_state());
}

#[test]
fn g1_tiny_heat_survives_large_reference_enthalpy() {
    let knots = [1.0, 1e17].map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(1e-10)));
    let solvent =
        PureLiquidHeatCapacityV1::new(species("solvent"), Temperature::new(1e17), knots.to_vec())
            .unwrap();
    let law = IdealSolubilityV1::new(
        species("solvent"),
        species("A"),
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        0.0,
        [Temperature::new(1.0), Temperature::new(1e17)],
    )
    .unwrap();
    let component = SolubilityCaloricComponentV1::new(
        law,
        MolarHeatCapacityCurveV1::new(species("A"), Temperature::new(1e17), knots.to_vec())
            .unwrap(),
    )
    .unwrap();
    let model = SolubilityCaloricModelV1::new(2e-5, solvent, vec![component]).unwrap();
    let t0 = 300.0_f64;
    let t1 = t0.next_up();
    let initial = model
        .equilibrate(
            Temperature::new(t0),
            Amount::new(10.0),
            &[Amount::new(0.1)],
            || false,
        )
        .unwrap();
    let heat = (10.0 * 1e-10 + 0.1 * 1e-10) * (t1 - t0);
    let mut budget = controls(t0, t1);
    budget.energy_tolerance_j = 1e-35;
    let result = model.add_heat(&initial, heat, budget, || false).unwrap();
    assert_eq!(result.final_state().temperature().value(), t1);
    assert_eq!(result.final_state().enthalpy_j(), initial.enthalpy_j());
    assert!(result.enthalpy_change_j() > 0.0);
    near(result.enthalpy_change_j(), heat, 1e-35);
    near(result.energy_residual_j(), 0.0, 1e-35);
}

#[test]
fn g1_empty_solute_stock_recovers_pure_solvent_and_source_domains_remain_binding() {
    let c = component("A", 300.0, 500.0, 6000.0, 50.0);
    let model = model(vec![c.clone()]);
    let initial = model
        .equilibrate(
            Temperature::new(300.0),
            Amount::new(10.0),
            &[Amount::new(0.0)],
            || false,
        )
        .unwrap();
    let result = model
        .add_heat(&initial, 1000.0, controls(301.0, 303.0), || false)
        .unwrap();
    near(result.final_state().temperature().value(), 302.0, 1e-9);
    assert_eq!(result.final_state().dissolved_solute()[0].value(), 0.0);
    conserved(&initial, result.final_state());
    assert!(SolubilityCaloricComponentV1::new(c.law().clone(), common_cp("B", 50.0)).is_err());
    assert!(SolubilityCaloricModelV1::new(2e-5, solvent_cp(50.0), vec![c.clone(), c]).is_err());
    assert!(SolubilityCaloricModelV1::new(0.0, solvent_cp(50.0), vec![]).is_err());
    assert!(
        model
            .prepare(
                Temperature::new(300.0),
                Amount::new(10.0),
                &[Amount::new(0.1)],
                &[Amount::new(0.2)],
                || false,
            )
            .is_err()
    );
    assert!(
        model
            .equilibrate(
                Temperature::new(279.0),
                Amount::new(10.0),
                &[Amount::new(0.0)],
                || false,
            )
            .is_err()
    );
    assert!(
        model
            .equilibrate(
                Temperature::new(300.0),
                Amount::new(0.0),
                &[Amount::new(0.0)],
                || false,
            )
            .is_err()
    );
    assert!(
        model
            .equilibrate(Temperature::new(300.0), Amount::new(10.0), &[], || false,)
            .is_err()
    );
    assert!(
        model
            .equilibrate(
                Temperature::new(300.0),
                Amount::new(10.0),
                &[Amount::new(f64::NAN)],
                || false,
            )
            .is_err()
    );
    assert!(matches!(
        model.add_heat(&initial, 0.0, controls(279.0, 310.0), || false),
        Err(SolubilityCaloricErrorV1::Invalid(_))
    ));
    let other = SolubilityCaloricModelV1::new(
        2e-5,
        solvent_cp(51.0),
        vec![component("A", 300.0, 500.0, 6000.0, 50.0)],
    )
    .unwrap();
    assert!(other.enthalpy_change_j(&initial, &initial).is_err());
    assert!(
        other
            .add_heat(&initial, 0.0, controls(290.0, 310.0), || false)
            .is_err()
    );
}

#[test]
fn g1_iteration_heat_refusal_and_late_cancellation_leave_original_state_unchanged() {
    let model = model(vec![component("A", 300.0, 1000.0, 0.0, 50.0)]);
    let initial = model
        .equilibrate(
            Temperature::new(300.0),
            Amount::new(10.0),
            &[Amount::new(0.1)],
            || false,
        )
        .unwrap();
    let original = initial.clone();
    let budget = controls(290.0, 310.0);
    let heat = 505.0 * 1.3;
    let short = SolubilityCaloricSolveV1 {
        maximum_iterations: 1,
        ..budget
    };
    assert!(matches!(
        model.add_heat(&initial, heat, short, || false),
        Err(SolubilityCaloricErrorV1::IterationLimit)
    ));
    assert!(matches!(
        model.add_heat(&initial, 1e9, budget, || false),
        Err(SolubilityCaloricErrorV1::HeatOutsideBracket { .. })
    ));
    let mut calls = 0;
    let accepted = model
        .add_heat(&initial, heat, budget, || {
            calls += 1;
            false
        })
        .unwrap();
    near(accepted.final_state().temperature().value(), 301.3, 1e-9);
    assert!(calls > 2);
    let mut cancelled_calls = 0;
    assert_eq!(
        model.add_heat(&initial, heat, budget, || {
            cancelled_calls += 1;
            cancelled_calls == calls
        }),
        Err(SolubilityCaloricErrorV1::Cancelled),
    );
    assert_eq!(cancelled_calls, calls);
    assert_eq!(initial, original);
    assert_eq!(
        model.add_heat(&initial, heat, budget, || false).unwrap(),
        accepted,
    );
}

#[test]
fn g1_component_work_bound_is_explicit() {
    let components = (0..=MAX_SOLUBILITY_CALORIC_COMPONENTS_V1)
        .map(|index| component(&format!("S{index}"), 300.0, 1000.0, 0.0, 50.0))
        .collect::<Vec<_>>();
    assert!(SolubilityCaloricModelV1::new(2e-5, solvent_cp(50.0), components.clone()).is_err());
    let model = model(components[..MAX_SOLUBILITY_CALORIC_COMPONENTS_V1].to_vec());
    let amounts = vec![Amount::new(0.001); MAX_SOLUBILITY_CALORIC_COMPONENTS_V1];
    let state = model
        .equilibrate(Temperature::new(300.0), Amount::new(10.0), &amounts, || {
            false
        })
        .unwrap();
    assert_eq!(
        state.equilibrium().unwrap().len(),
        MAX_SOLUBILITY_CALORIC_COMPONENTS_V1
    );
    assert_eq!(state.dissolved_solute(), amounts);
}

#[test]
fn g1_material_volume_is_preserved_without_quotient_product_drift() {
    let volume = 0.00011;
    let molar_volume = 1.8e-5;
    assert_ne!((volume / molar_volume) * molar_volume, volume);
    let model = SolubilityCaloricModelV1::new(
        molar_volume,
        solvent_cp(50.0),
        vec![component("A", 300.0, 1000.0, 0.0, 50.0)],
    )
    .unwrap();
    let total = [Amount::new(0.1)];
    let initial = model
        .prepare_at_solvent_volume(
            Temperature::new(300.0),
            Volume::new(volume),
            &total,
            &total,
            || false,
        )
        .unwrap();
    assert_eq!(initial.solvent_volume().value(), volume);
    assert_eq!(initial.solvent_moles().value(), volume / molar_volume);
    let heat = (volume / molar_volume + 0.1) * 50.0;
    let result = model
        .add_heat(&initial, heat, controls(300.0, 302.0), || false)
        .unwrap();
    assert_eq!(result.final_state().solvent_volume().value(), volume);
    near(result.final_state().temperature().value(), 301.0, 1e-9);
    conserved(&initial, result.final_state());
    let by_moles = model
        .prepare(
            initial.temperature(),
            initial.solvent_moles(),
            &total,
            &total,
            || false,
        )
        .unwrap();
    assert!(model.enthalpy_change_j(&initial, &by_moles).is_err());
}
