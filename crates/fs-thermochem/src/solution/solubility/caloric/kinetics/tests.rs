use super::*;
use crate::solution::caloric::{MolarHeatCapacityCurveV1, PureLiquidHeatCapacityV1};
use crate::solution::solubility::IdealSolubilityV1;
use crate::solution::solubility::caloric::SolubilityCaloricComponentV1;
use crate::solution::solubility::kinetics::PrecipitationPolicyV1;
use crate::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};
use fs_qty::{AmountConcentration, VolumetricFlowRate};

fn species(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn phase_cp(name: &str, value: f64) -> MolarHeatCapacityCurveV1 {
    MolarHeatCapacityCurveV1::new(
        species(name),
        Temperature::new(300.0),
        [280.0, 340.0]
            .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(value)))
            .to_vec(),
    )
    .unwrap()
}
pub(super) fn component(
    name: &str,
    concentration: f64,
    enthalpy: f64,
    solid_cp: f64,
    dissolved_cp: f64,
) -> SolubilityCaloricComponentV1 {
    let law = IdealSolubilityV1::new(
        species("solvent"),
        species(name),
        Temperature::new(300.0),
        AmountConcentration::new(concentration),
        enthalpy,
        [Temperature::new(280.0), Temperature::new(340.0)],
    )
    .unwrap();
    if solid_cp == dissolved_cp {
        SolubilityCaloricComponentV1::new(law, phase_cp(name, solid_cp)).unwrap()
    } else {
        SolubilityCaloricComponentV1::from_phase_heat_capacities(
            law.with_phase_heat_capacities(phase_cp(name, solid_cp), phase_cp(name, dissolved_cp))
                .unwrap(),
        )
        .unwrap()
    }
}
pub(super) fn model(
    components: Vec<SolubilityCaloricComponentV1>,
    conductances: &[f64],
    policy: PrecipitationPolicyV1,
) -> SolubilityThermalKineticsV1 {
    let kinetics = components
        .iter()
        .zip(conductances)
        .map(|(c, &g)| {
            SolubilityKineticsV1::new(c.law().clone(), VolumetricFlowRate::new(g), policy).unwrap()
        })
        .collect();
    let solvent = PureLiquidHeatCapacityV1::new(
        species("solvent"),
        Temperature::new(300.0),
        [280.0, 340.0]
            .map(|t| (Temperature::new(t), MolarThermalQuantityV1::new(50.0)))
            .to_vec(),
    )
    .unwrap();
    SolubilityThermalKineticsV1::new(
        SolubilityCaloricModelV1::new(2e-5, solvent, components).unwrap(),
        kinetics,
    )
    .unwrap()
}
fn one(enthalpy: f64, conductance: f64) -> SolubilityThermalKineticsV1 {
    model(
        vec![component("A", 600.0, enthalpy, 30.0, 30.0)],
        &[conductance],
        PrecipitationPolicyV1::SeededSurface,
    )
}
pub(super) fn state(
    model: &SolubilityThermalKineticsV1,
    phases: &[(f64, f64)],
) -> SolubilityThermalKineticStateV1 {
    let phases: Vec<_> = phases
        .iter()
        .map(|&(total, dissolved)| {
            SolubilityKineticStateV1::new(Amount::new(total), Amount::new(dissolved)).unwrap()
        })
        .collect();
    model
        .prepare_at_solvent_volume(Temperature::new(300.0), Volume::new(1e-3), &phases, || {
            false
        })
        .unwrap()
}
pub(super) fn controls(step: f64) -> SolubilityThermalKineticControlsV1 {
    SolubilityThermalKineticControlsV1 {
        root: SolubilityCaloricSolveV1 {
            temperature_bracket: [Temperature::new(280.0), Temperature::new(340.0)],
            temperature_tolerance_k: 1e-10,
            energy_tolerance_j: 1e-7,
            maximum_iterations: 128,
        },
        initial_step: Time::new(step),
        minimum_step: Time::new(1e-6),
        maximum_step: Time::new(step),
        relative_tolerance: 2e-5,
        amount_absolute_tolerance_mol: 1e-10,
        temperature_absolute_tolerance_k: 1e-7,
        maximum_attempts: 8192,
    }
}
pub(super) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}
fn conserved(
    initial: &SolubilityThermalKineticStateV1,
    final_state: &SolubilityThermalKineticStateV1,
) {
    assert_eq!(initial.solvent_moles(), final_state.solvent_moles());
    assert_eq!(initial.solvent_volume(), final_state.solvent_volume());
    for (a, b) in initial.phases().iter().zip(final_state.phases()) {
        assert_eq!(a.total_solute(), b.total_solute());
        assert!(b.solid().value() >= 0.0 && b.dissolved().value() >= 0.0);
        assert!(
            b.balance_residual().value().abs() <= 8.0 * f64::EPSILON * b.total_solute().value()
        );
    }
}

#[test]
fn particle_geometry_survives_thermal_stages_without_resetting_its_reference() {
    let fixed = model(
        vec![component("A", 1000.0, 0.0, 30.0, 30.0)],
        &[1e-3],
        PrecipitationPolicyV1::ExistingSolidOnly,
    );
    let particles = SolubilityThermalKineticsV1::new(
        fixed.caloric().clone(),
        vec![
            fixed.kinetics()[0]
                .clone()
                .with_particle_surface(Amount::new(0.2))
                .unwrap(),
        ],
    )
    .unwrap();
    let initial = state(&particles, &[(1.0, 0.8)]);
    let mut ctl = controls(10.0);
    ctl.relative_tolerance = 1e-8;
    let result = particles
        .advance(&initial, Time::new(50.0), 0.0, ctl, || false)
        .unwrap();
    let expected = 0.2 / (1.0 + (2.0 / 3.0) * 50.0_f64).powf(1.5);
    near(
        result.final_state().phases()[0].solid().value(),
        expected,
        1e-12 * expected,
    );
    assert_eq!(result.final_state().temperature(), initial.temperature());
    assert_eq!(result.enthalpy_change_j(), 0.0);
    assert_eq!(result.root_iterations(), 0);
    assert_eq!(
        particles.kinetics()[0].particle_reference_solid(),
        Some(Amount::new(0.2))
    );
    conserved(&initial, result.final_state());
}

#[test]
fn zero_enthalpy_recovers_exact_isothermal_rate_and_positive_donor_tail() {
    let model = model(
        vec![component("A", 1000.0, 0.0, 30.0, 30.0)],
        &[1e-3],
        PrecipitationPolicyV1::ExistingSolidOnly,
    );
    let initial = state(&model, &[(1.0, 0.8)]);
    let mut ctl = controls(50.0);
    ctl.relative_tolerance = 1e-8;
    ctl.amount_absolute_tolerance_mol = 1e-35;
    let result = model
        .advance(&initial, Time::new(50.0), 0.0, ctl, || false)
        .unwrap();
    let expected = initial.phases()[0].solid().value() * (-50.0_f64).exp();
    assert_eq!(result.final_state().temperature(), initial.temperature());
    assert_eq!(result.final_state().phases()[0].dissolved().value(), 1.0);
    near(
        result.final_state().phases()[0].solid().value(),
        expected,
        1e-13 * expected,
    );
    assert!(result.final_state().phases()[0].solid().value() > 0.0);
    assert_eq!(result.enthalpy_change_j(), 0.0);
    assert_eq!(result.root_iterations(), 0);
    conserved(&initial, result.final_state());
    let later = model
        .advance(result.final_state(), Time::new(10.0), 0.0, ctl, || false)
        .unwrap();
    near(
        later.final_state().phases()[0].solid().value(),
        expected * (-10.0_f64).exp(),
        1e-13 * expected * (-10.0_f64).exp(),
    );
    assert!(later.final_state().phases()[0].solid().value() > 0.0);
}

#[test]
fn zero_conductance_heats_actual_unequal_phases_without_equilibrating() {
    let model = model(
        vec![component("A", 600.0, 1000.0, 30.0, 90.0)],
        &[0.0],
        PrecipitationPolicyV1::SeededSurface,
    );
    let initial = state(&model, &[(1.0, 0.25)]);
    let capacity = 50.0 * 50.0 + 0.75 * 30.0 + 0.25 * 90.0;
    let result = model
        .advance(
            &initial,
            Time::new(5.0),
            capacity * 10.0,
            controls(2.0),
            || false,
        )
        .unwrap();
    near(result.final_state().temperature().value(), 310.0, 1e-9);
    assert_eq!(result.final_state().phases(), initial.phases());
    near(result.enthalpy_change_j(), capacity * 10.0, 1e-7);
    near(
        model
            .enthalpy_change_j(&initial, result.final_state())
            .unwrap(),
        capacity * 10.0,
        1e-7,
    );
    assert_eq!(result.heat_input_j(), capacity * 10.0);
    conserved(&initial, result.final_state());
}

#[test]
fn adiabatic_dissolution_cools_or_warms_and_feeds_back_on_the_rate() {
    for enthalpy in [20_000.0, -20_000.0] {
        let model = one(enthalpy, 1e-4);
        let initial = state(&model, &[(2.0, 0.05)]);
        let frozen = model.kinetics()[0]
            .advance(
                initial.temperature(),
                initial.solvent_volume(),
                &initial.phases()[0],
                Time::new(10.0),
            )
            .unwrap();
        let result = model
            .advance(&initial, Time::new(10.0), 0.0, controls(3.0), || false)
            .unwrap();
        let final_state = result.final_state();
        let dissolved = final_state.phases()[0].dissolved().value();
        assert!(dissolved > 0.05 && dissolved < frozen.state().dissolved().value());
        assert!((final_state.temperature().value() - 300.0) * enthalpy < 0.0);
        let independent = (50.0 * 50.0 + 2.0 * 30.0) * (final_state.temperature().value() - 300.0)
            + (dissolved - 0.05) * enthalpy;
        near(independent, 0.0, 1e-7);
        near(result.energy_residual_j(), independent, 1e-10);
        assert!(result.maximum_accepted_error_norm() <= 1.0);
        conserved(&initial, final_state);
    }
}

#[test]
fn adiabatic_precipitation_respects_growth_surfaces_and_thermal_feedback() {
    for enthalpy in [12_000.0, -12_000.0] {
        for policy in [
            PrecipitationPolicyV1::ExistingSolidOnly,
            PrecipitationPolicyV1::SeededSurface,
        ] {
            let model = model(
                vec![component("A", 200.0, enthalpy, 30.0, 30.0)],
                &[1e-3],
                policy,
            );
            let initial = state(&model, &[(0.9, 0.9)]);
            let result = model
                .advance(&initial, Time::new(2.0), 0.0, controls(1.0), || false)
                .unwrap();
            if policy == PrecipitationPolicyV1::ExistingSolidOnly {
                assert_eq!(result.final_state(), &initial);
            } else {
                let frozen = 0.2 + 0.7 * (-2.0_f64).exp();
                let final_state = result.final_state();
                assert!(final_state.phases()[0].dissolved().value() > frozen);
                assert!(final_state.phases()[0].solid().value() > 0.0);
                assert!((final_state.temperature().value() - 300.0) * enthalpy > 0.0);
                let independent = (2500.0 + 0.9 * 30.0)
                    * (final_state.temperature().value() - 300.0)
                    + (final_state.phases()[0].dissolved().value() - 0.9) * enthalpy;
                near(independent, 0.0, 1e-7);
            }
            conserved(&initial, result.final_state());
        }
    }
}

#[test]
fn finite_time_lags_but_long_time_reaches_the_same_adiabatic_equilibrium() {
    let model = one(4000.0, 1e-3);
    let initial = state(&model, &[(2.0, 0.0)]);
    let declared = model
        .caloric()
        .prepare_at_solvent_volume(
            initial.temperature(),
            initial.solvent_volume(),
            &[Amount::new(2.0)],
            &[Amount::new(0.0)],
            || false,
        )
        .unwrap();
    let equilibrium = model
        .caloric()
        .add_heat(&declared, 0.0, controls(1.0).root, || false)
        .unwrap();
    let short = model
        .advance(&initial, Time::new(0.1), 0.0, controls(0.1), || false)
        .unwrap();
    let long = model
        .advance(&initial, Time::new(40.0), 0.0, controls(10.0), || false)
        .unwrap();
    assert!(
        short.final_state().phases()[0].dissolved().value()
            < 0.25 * equilibrium.final_state().dissolved_solute()[0].value()
    );
    near(
        long.final_state().phases()[0].dissolved().value(),
        equilibrium.final_state().dissolved_solute()[0].value(),
        2e-8,
    );
    near(
        long.final_state().temperature().value(),
        equilibrium.final_state().temperature().value(),
        1e-7,
    );
    conserved(&initial, long.final_state());
}

#[test]
fn fixed_step_refinement_converges_at_first_order_to_an_independent_ode() {
    let model = one(40_000.0, 1e-3);
    let initial = state(&model, &[(2.0, 0.02)]);
    // Independent ODE reduction: constant phase Cp and the first law eliminate
    // T analytically; RK4 supplies a much finer reference for dn/dt alone.
    let derivative = |n: f64| {
        let t = 300.0 - 40_000.0 * (n - 0.02) / 2560.0;
        0.6 * (40_000.0 / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (1.0 / 300.0 - 1.0 / t)).exp() - n
    };
    let mut reference = 0.02;
    let dt = 2.0 / 8192.0;
    for _ in 0..8192 {
        let k1 = derivative(reference);
        let k2 = derivative(reference + 0.5 * dt * k1);
        let k3 = derivative(reference + 0.5 * dt * k2);
        let k4 = derivative(reference + dt * k3);
        reference += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let errors = [1.0, 0.5, 0.25].map(|h| {
        let mut ctl = controls(h);
        ctl.minimum_step = Time::new(h);
        ctl.relative_tolerance = 1.0;
        ctl.amount_absolute_tolerance_mol = 1.0;
        ctl.temperature_absolute_tolerance_k = 10.0;
        let result = model
            .advance(&initial, Time::new(2.0), 0.0, ctl, || false)
            .unwrap();
        assert_eq!(result.rejected_steps(), 0);
        (result.final_state().phases()[0].dissolved().value() - reference).abs()
    });
    assert!(errors[0] > 1e-6, "{errors:?}");
    assert!(
        errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn two_components_share_one_bath_with_the_original_distinct_phase_capacities() {
    let model = model(
        vec![
            component("A", 400.0, 1000.0, 30.0, 90.0),
            component("B", 100.0, -8000.0, 20.0, 45.0),
        ],
        &[3e-4, 2e-4],
        PrecipitationPolicyV1::SeededSurface,
    );
    let initial = state(&model, &[(1.0, 0.1), (0.7, 0.65)]);
    let result = model
        .advance(&initial, Time::new(3.0), 1500.0, controls(1.0), || false)
        .unwrap();
    let final_state = result.final_state();
    let dt = final_state.temperature().value() - 300.0;
    let da = final_state.phases()[0].dissolved().value() - 0.1;
    let db = final_state.phases()[1].dissolved().value() - 0.65;
    assert!(da > 0.0 && db < 0.0);
    let independent = (2500.0 + 0.9 * 30.0 + 0.1 * 90.0 + 0.05 * 20.0 + 0.65 * 45.0) * dt
        + da * (1000.0 + 60.0 * dt)
        + db * (-8000.0 + 25.0 * dt);
    near(independent, 1500.0, 1e-7);
    near(result.enthalpy_change_j(), independent, 1e-10);
    assert_eq!(result.duration(), Time::new(3.0));
    conserved(&initial, final_state);
}

#[test]
fn source_refusals_cancellation_and_work_limits_leave_the_input_unchanged() {
    let model = model(
        vec![component("A", 400.0, 1000.0, 30.0, 90.0)],
        &[3e-4],
        PrecipitationPolicyV1::SeededSurface,
    );
    let wrong_law = component("A", 400.0, 1000.0, 31.0, 90.0).law().clone();
    let wrong_kinetics = SolubilityKineticsV1::new(
        wrong_law,
        VolumetricFlowRate::new(3e-4),
        PrecipitationPolicyV1::SeededSurface,
    )
    .unwrap();
    assert!(matches!(
        SolubilityThermalKineticsV1::new(model.caloric().clone(), vec![wrong_kinetics]),
        Err(SolubilityCaloricErrorV1::Invalid(_))
    ));
    let initial = state(&model, &[(1.0, 0.1)]);
    let before = initial.clone();
    let mut polls = 0;
    let cancelled = model.advance(&initial, Time::new(4.0), 0.0, controls(2.0), || {
        polls += 1;
        polls >= 12
    });
    assert_eq!(cancelled, Err(SolubilityCaloricErrorV1::Cancelled));
    assert!(polls >= 12);
    let mut strict = controls(2.0);
    strict.relative_tolerance = 1e-12;
    strict.amount_absolute_tolerance_mol = 1e-14;
    strict.temperature_absolute_tolerance_k = 1e-12;
    strict.maximum_attempts = 1;
    assert_eq!(
        model.advance(&initial, Time::new(4.0), 0.0, strict, || false),
        Err(SolubilityCaloricErrorV1::IntegrationLimit)
    );
    strict.minimum_step = strict.initial_step;
    assert!(matches!(
        model.advance(&initial, Time::new(4.0), 0.0, strict, || false),
        Err(SolubilityCaloricErrorV1::TemporalTolerance { .. })
    ));
    let mut outside = controls(1.0);
    outside.minimum_step = outside.initial_step;
    outside.root.temperature_bracket = [Temperature::new(300.0), Temperature::new(301.0)];
    assert!(matches!(
        model.advance(&initial, Time::new(1.0), 0.0, outside, || false),
        Err(SolubilityCaloricErrorV1::HeatOutsideBracket { .. })
    ));
    assert!(matches!(
        model.advance(&initial, Time::new(0.0), 1.0, controls(1.0), || false),
        Err(SolubilityCaloricErrorV1::Invalid(_))
    ));
    let zero = model
        .advance(&initial, Time::new(0.0), 0.0, controls(1.0), || false)
        .unwrap();
    assert_eq!(zero.final_state(), &initial);
    assert_eq!(zero.accepted_steps(), 0);
    assert_eq!(initial, before);
}
