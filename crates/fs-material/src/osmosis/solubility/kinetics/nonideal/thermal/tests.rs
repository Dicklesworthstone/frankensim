//! Independent potential derivatives, trace increments and thermal dynamics.

use super::*;
use fs_qty::{Amount, Time, VolumetricFlowRate};
use fs_thermochem::{
    MolarThermalQuantityV1, SpeciesId,
    solution::{
        caloric::{MolarHeatCapacityCurveV1, PureLiquidHeatCapacityV1},
        solubility::{
            IdealSolubilityV1,
            caloric::{
                SolubilityCaloricSolveV1,
                kinetics::{SolubilityThermalKineticControlsV1, SolubilityThermalKineticsV1},
            },
            kinetics::PrecipitationPolicyV1,
        },
    },
};

const T: f64 = 300.0;
const V: f64 = 1e-3;
const VBAR: f64 = 18e-6;
const CP_WATER: f64 = 75.0;

fn temperature(t: f64) -> Temperature {
    Temperature::new(t)
}

fn knots(cp: f64) -> Vec<(Temperature, MolarThermalQuantityV1)> {
    vec![
        (temperature(280.0), MolarThermalQuantityV1::new(cp)),
        (temperature(340.0), MolarThermalQuantityV1::new(cp)),
    ]
}

fn sources(count: usize, unequal_cp: bool) -> SolubilityCaloricModelV1 {
    let water = SpeciesId::new("water").unwrap();
    let components = (0..count)
        .map(|i| {
            let solute = SpeciesId::new(format!("thermal-component-{i}")).unwrap();
            let law = IdealSolubilityV1::new(
                water.clone(),
                solute.clone(),
                temperature(T),
                AmountConcentration::new([100.0, 80.0][i]),
                [6000.0, -4000.0][i],
                [temperature(280.0), temperature(340.0)],
            )
            .unwrap();
            if unequal_cp {
                let solid = MolarHeatCapacityCurveV1::new(
                    solute.clone(),
                    temperature(T),
                    knots(40.0 + i as f64 * 10.0),
                )
                .unwrap();
                let dissolved = MolarHeatCapacityCurveV1::new(
                    solute,
                    temperature(T),
                    knots(70.0 + i as f64 * 10.0),
                )
                .unwrap();
                SolubilityCaloricComponentV1::from_phase_heat_capacities(
                    law.with_phase_heat_capacities(solid, dissolved).unwrap(),
                )
                .unwrap()
            } else {
                SolubilityCaloricComponentV1::new(
                    law,
                    MolarHeatCapacityCurveV1::new(
                        solute,
                        temperature(T),
                        knots(50.0 + i as f64 * 10.0),
                    )
                    .unwrap(),
                )
                .unwrap()
            }
        })
        .collect();
    SolubilityCaloricModelV1::new(
        VBAR,
        PureLiquidHeatCapacityV1::new(water, temperature(T), knots(CP_WATER)).unwrap(),
        components,
    )
    .unwrap()
}

fn family(b: Vec<f64>, c: Vec<f64>, unequal_cp: bool) -> VirialSolubilityThermalModel {
    let count = if b.len() == 1 { 1 } else { 2 };
    VirialSolubilityThermalModel::new(
        sources(count, unequal_cp),
        Arc::new(VirialSolutionModel::new(T, VBAR, 1.0, vec![500.0; count], b).unwrap()),
        c,
        [temperature(290.0), temperature(320.0)],
    )
    .unwrap()
}

fn phases(amounts: &[f64], total: f64) -> Vec<SolubilityKineticStateV1> {
    amounts.iter().map(|&n| {
        SolubilityKineticStateV1::new(Amount::new(total), Amount::new(n)).unwrap()
    }).collect()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e}"
    );
}

// Independent absolute expressions are used ONLY as finite-difference oracles.
// Production thermal work never subtracts these large absolute energies.
fn free_energy(source: &VirialSolubilityThermalModel, t: f64, volume: f64, amounts: &[f64]) -> f64 {
    let view = source.at_temperature(temperature(t)).unwrap();
    let b = view.interaction_matrix_m3_mol();
    let n = amounts.len();
    let mut reduced = 0.0;
    for i in 0..n {
        for j in 0..n {
            reduced += 0.5 * amounts[i] * b[i * n + j] * amounts[j] / volume;
        }
        let sat = source.caloric().components()[i].law()
            .saturation_concentration(temperature(t)).unwrap().value();
        reduced -= amounts[i] * b[i * n + i] * sat;
    }
    R * t * reduced
}

fn absolute_heat(source: &VirialSolubilityThermalModel, t: f64, volume: f64, amounts: &[f64]) -> f64 {
    let view = source.at_temperature(temperature(t)).unwrap();
    let b = view.interaction_matrix_m3_mol();
    let c = source.inverse_temperature_m3_k_mol();
    let n = amounts.len();
    let mut heat = 0.0;
    for i in 0..n {
        for j in 0..n {
            heat += 0.5 * R * amounts[i] * c[i * n + j] * amounts[j] / volume;
        }
        let law = source.caloric().components()[i].law();
        let sat = law.saturation_concentration(temperature(t)).unwrap().value();
        let h = law.dissolution_enthalpy_at(temperature(t)).unwrap();
        heat += amounts[i] * sat * (b[i * n + i] * h - R * c[i * n + i]);
    }
    heat
}

fn capacity(source: &VirialSolubilityThermalModel, t: f64, amounts: &[f64]) -> f64 {
    let view = source.at_temperature(temperature(t)).unwrap();
    let b = view.interaction_matrix_m3_mol();
    let c = source.inverse_temperature_m3_k_mol();
    let count = amounts.len();
    amounts.iter().enumerate().map(|(i, &n)| {
        let component = &source.caloric().components()[i];
        let law = component.law();
        let sat = law.saturation_concentration(temperature(t)).unwrap().value();
        let h = law.dissolution_enthalpy_at(temperature(t)).unwrap();
        let dcp = component.dissolved_heat_capacity().heat_capacity_at(temperature(t)).unwrap().value()
            - component.solid_heat_capacity().heat_capacity_at(temperature(t)).unwrap().value();
        n * sat * (b[i * count + i] * (h * h / (R * t * t) + dcp)
            - 2.0 * c[i * count + i] * h / (t * t))
    }).sum()
}

fn derivative(mut f: impl FnMut(f64) -> f64, t: f64, h: f64) -> f64 {
    (f(t - 2.0 * h) - 8.0 * f(t - h) + 8.0 * f(t + h) - f(t + 2.0 * h)) / (12.0 * h)
}

#[test]
fn original_views_and_temperature_dependent_conditional_activity_share_one_family() {
    let source = family(vec![0.001, 0.0007, 0.0007, 0.0015], vec![0.4, -0.1, -0.1, 0.3], true);
    assert!(Arc::ptr_eq(&source.at_temperature(temperature(T)).unwrap(), source.reference_model()));
    let stocks = phases(&[0.02, 0.12], 0.3);
    for t in [290.0, 307.0, 320.0] {
        let view = source.at_temperature(temperature(t)).unwrap();
        for (i, &actual) in view.interaction_matrix_m3_mol().iter().enumerate() {
            let expected = source.reference_model().interaction_matrix_m3_mol()[i]
                + source.inverse_temperature_m3_k_mol()[i] * (1.0 / t - 1.0 / T);
            close(actual, expected, 3e-15);
        }
        let conditional = source.conditional_saturation(
            temperature(t), Volume::new(V), &stocks, 0,
        ).unwrap().value();
        let pure = source.caloric().components()[0].law()
            .saturation_concentration(temperature(t)).unwrap().value();
        let b = view.interaction_matrix_m3_mol();
        let force = (conditional / pure).ln() + b[0] * (conditional - pure) + b[1] * 120.0;
        assert!(force.abs() < 3e-14, "{force:e}");
        let pure_stocks = phases(&[0.02, 0.0], 0.3);
        assert_eq!(source.conditional_saturation(
            temperature(t), Volume::new(V), &pure_stocks, 0,
        ).unwrap().value(), pure);
    }
    assert!(source.at_temperature(temperature(289.0)).is_err());
}

#[test]
fn excess_heat_and_capacity_follow_gibbs_helmholtz_including_constant_b_reference() {
    for c in [vec![0.0; 4], vec![0.4, -0.1, -0.1, 0.3]] {
        let source = family(vec![0.001, 0.0007, 0.0007, 0.0015], c, true);
        let amounts = [0.04, 0.08];
        let all_solid = phases(&[0.0, 0.0], 0.3);
        let dissolved = phases(&amounts, 0.3);
        let t = 307.0;
        let expected = absolute_heat(&source, t, V, &amounts);
        let actual = source.enthalpy_increment_j(
            Volume::new(V), temperature(t), &all_solid, temperature(t), &dissolved,
        ).unwrap();
        close(actual, expected, 3e-14);
        assert!(actual.abs() > 1.0);
        let gibbs_helmholtz = free_energy(&source, t, V, &amounts)
            - t * derivative(|x| free_energy(&source, x, V, &amounts), t, 0.02);
        close(actual, gibbs_helmholtz, 3e-9);
        close(
            derivative(|x| absolute_heat(&source, x, V, &amounts), t, 0.02),
            capacity(&source, t, &amounts),
            3e-9,
        );
    }
}

#[test]
fn adjacent_temperature_and_trace_phase_heat_do_not_subtract_large_datums() {
    let source = family(vec![0.002], vec![0.4], false);
    let initial = phases(&[0.08], 0.2);
    let next = f64::from_bits(T.to_bits() + 1);
    let heat = source.enthalpy_increment_j(
        Volume::new(V), temperature(T), &initial, temperature(next), &initial,
    ).unwrap();
    let expected = capacity(&source, T, &[0.08]) * (next - T);
    assert_ne!(heat, 0.0);
    close(heat, expected, 2e-13);

    let before = [SolubilityKineticStateV1::from_phases(Amount::new(1.0), Amount::new(1e-18)).unwrap()];
    let after = [SolubilityKineticStateV1::from_phases(Amount::new(1.0), Amount::new(0.5e-18)).unwrap()];
    assert_eq!(before[0].dissolved(), after[0].dissolved());
    let delta = before[0].solid().value() - after[0].solid().value();
    let actual = source.enthalpy_increment_j(
        Volume::new(0.01), temperature(T), &before, temperature(T), &after,
    ).unwrap();
    let expected = delta * (R * 0.4 * 100.0 + 100.0 * (0.002 * 6000.0 - R * 0.4));
    assert_ne!(actual, 0.0);
    close(actual, expected, 5e-15);
}

#[test]
fn complete_temperature_and_capacity_support_is_required() {
    let original = sources(1, false);
    let reference = Arc::new(VirialSolutionModel::new(T, VBAR, 1.0, vec![500.0], vec![0.001]).unwrap());
    assert!(VirialSolubilityThermalModel::new(
        original.clone(), Arc::clone(&reference), vec![0.0], [temperature(279.0), temperature(320.0)],
    ).is_err());
    assert!(VirialSolubilityThermalModel::new(
        original.clone(), Arc::clone(&reference), vec![100.0], [temperature(290.0), temperature(320.0)],
    ).is_err());
    let unsupported_cap = Arc::new(VirialSolutionModel::new(T, VBAR, 1.0, vec![101.0], vec![0.001]).unwrap());
    assert!(VirialSolubilityThermalModel::new(
        original.clone(), unsupported_cap, vec![0.0], [temperature(290.0), temperature(320.0)],
    ).is_err());

    let water = SpeciesId::new("water").unwrap();
    let solute = SpeciesId::new("interior-valley").unwrap();
    let solid = MolarHeatCapacityCurveV1::new(solute.clone(), temperature(T), knots(100.0)).unwrap();
    let dissolved = MolarHeatCapacityCurveV1::new(solute.clone(), temperature(T), vec![
        (temperature(280.0), MolarThermalQuantityV1::new(100.0)),
        (temperature(290.0), MolarThermalQuantityV1::new(100.0)),
        (temperature(300.0), MolarThermalQuantityV1::new(1.0)),
        (temperature(320.0), MolarThermalQuantityV1::new(100.0)),
        (temperature(340.0), MolarThermalQuantityV1::new(100.0)),
    ]).unwrap();
    let law = IdealSolubilityV1::new(
        water.clone(), solute, temperature(T), AmountConcentration::new(100.0), 0.0,
        [temperature(280.0), temperature(340.0)],
    ).unwrap().with_phase_heat_capacities(solid, dissolved).unwrap();
    let caloric = SolubilityCaloricModelV1::new(
        VBAR,
        PureLiquidHeatCapacityV1::new(water, temperature(T), knots(CP_WATER)).unwrap(),
        vec![SolubilityCaloricComponentV1::from_phase_heat_capacities(law).unwrap()],
    ).unwrap();
    // Both endpoint dissolved Cp values are 100, but the interior effective
    // value is 1 + 100*0.02*(1-100) = -197 J/(mol K).
    let rejected = VirialSolubilityThermalModel::new(
        caloric,
        Arc::new(VirialSolutionModel::new(T, VBAR, 1.0, vec![500.0], vec![0.02]).unwrap()),
        vec![0.0],
        [temperature(290.0), temperature(320.0)],
    );
    assert!(matches!(rejected, Err(SolubilityCaloricErrorV1::Invalid(
        "positive effective dissolved heat capacity cannot be certified over thermal source"
    ))));
}

fn controls(step: f64) -> SolubilityThermalKineticControlsV1 {
    SolubilityThermalKineticControlsV1 {
        root: SolubilityCaloricSolveV1 {
            temperature_bracket: [temperature(290.0), temperature(320.0)],
            temperature_tolerance_k: 1e-11,
            energy_tolerance_j: 1e-8,
            maximum_iterations: 128,
        },
        initial_step: Time::new(step),
        minimum_step: Time::new(1e-6),
        maximum_step: Time::new(step),
        relative_tolerance: 1e-2,
        amount_absolute_tolerance_mol: 1e-9,
        temperature_absolute_tolerance_k: 1e-7,
        maximum_attempts: 8192,
    }
}

// Independent continuous particle/temperature ODE with the full calibrated
// potential. One component's pure saturation is analytic; no production phase
// step, conditional root or thermal projection is used by this reference.
fn continuous_particle(duration: f64) -> [f64; 2] {
    let rhs = |state: [f64; 2]| {
        let [amount, t] = state;
        let solid = 0.2 - amount;
        let saturation = 100.0 * (6000.0 / R * (1.0 / T - 1.0 / t)).exp();
        let b = 0.002 + 0.4 * (1.0 / t - 1.0 / T);
        let q = saturation * (b * 6000.0 - R * 0.4);
        let dq = saturation * (b * 6000.0 * 6000.0 / (R * t * t) - 2.0 * 0.4 * 6000.0 / (t * t));
        let amount_rate = 2e-4 * (solid / 0.18).powf(2.0 / 3.0) * (saturation - amount / V);
        let heat_capacity = V / VBAR * CP_WATER + 0.2 * 50.0 + amount * dq;
        let transfer_enthalpy = 6000.0 + q + R * 0.4 * amount / V;
        [amount_rate, -transfer_enthalpy * amount_rate / heat_capacity]
    };
    let h = duration / 16384.0;
    let mut state = [0.02, T];
    let shifted = |s: [f64; 2], d: [f64; 2], h: f64| [s[0] + h*d[0], s[1] + h*d[1]];
    for _ in 0..16384 {
        let a = rhs(state);
        let b = rhs(shifted(state, a, 0.5*h));
        let c = rhs(shifted(state, b, 0.5*h));
        let d = rhs(shifted(state, c, h));
        for i in 0..2 {
            state[i] += h * (a[i] + 2.0*b[i] + 2.0*c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn shared_particle_heat_controller_refines_to_independent_continuous_dynamics() {
    let source = Arc::new(family(vec![0.002], vec![0.4], false));
    let kinetic = SolubilityKineticsV1::new(
        source.caloric().components()[0].law().clone(),
        VolumetricFlowRate::new(2e-4),
        PrecipitationPolicyV1::ExistingSolidOnly,
    ).unwrap().with_particle_surface(Amount::new(0.18)).unwrap();
    let thermal = SolubilityThermalKineticsV1::new(source.caloric().clone(), vec![kinetic])
        .unwrap().with_interactions(source.clone()).unwrap();
    let initial = thermal.prepare_at_solvent_volume(
        temperature(T), Volume::new(V), &phases(&[0.02], 0.2), || false,
    ).unwrap();
    let expected = continuous_particle(2.0);
    let mut errors = Vec::new();
    for step in [0.2, 0.05] {
        let outcome = thermal.advance(&initial, Time::new(2.0), 0.0, controls(step), || false).unwrap();
        let actual = outcome.final_state();
        assert!(actual.temperature().value() < T);
        assert_eq!(actual.phases()[0].total_solute(), Amount::new(0.2));
        assert!(outcome.energy_residual_j().abs() <= 1e-8);
        errors.push((actual.phases()[0].dissolved().value() - expected[0]).abs()
            + (actual.temperature().value() - expected[1]).abs() * 0.01);
    }
    assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
    assert!(errors[1] < 2e-6, "{errors:?}");
}
