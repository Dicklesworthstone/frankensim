//! G0/G1/G3: retained gas inventories, analytic equilibria and native dynamics.
use super::*;
use crate::cavitation::{
    BubbleGas, BubbleLimits, BubbleNumerics, KellerMiksisBubble, KellerMiksisLimits,
    PopulationEvolutionNumerics,
};
use fs_material::cavitation::CavitationState;

fn liquid(sigma: f64) -> CavitationState {
    CavitationState::try_new(300.0, 1000.0, 0.001, sigma, 2000.0, 1500.0).unwrap()
}

fn class(material: CavitationState, radius: f64, number: f64, kappa: f64) -> BubblePopulationClass {
    BubblePopulationClass::try_new(
        KellerMiksisBubble::try_new(
            material,
            BubbleGas::equilibrium(material, radius, 100_000.0, kappa).unwrap(),
            KellerMiksisLimits {
                bubble: BubbleLimits {
                    minimum_radius_m: 0.25 * radius,
                    maximum_radius_m: 2.0 * radius,
                    maximum_wall_mach_number: 0.01,
                    maximum_vapor_density_ratio: 0.01,
                },
                maximum_pressure_ratio: 0.01,
            },
        )
        .unwrap(),
        number,
    )
    .unwrap()
}

fn population(sigma: f64) -> BubblePopulation {
    let material = liquid(sigma);
    BubblePopulation::try_new(
        material,
        vec![
            class(material, 40e-6, 4e6, 1.0),
            class(material, 120e-6, 1e6, 1.4),
        ],
        2,
    )
    .unwrap()
}

fn controls() -> EquilibriumControls {
    EquilibriumControls {
        relative_radius_tolerance: 1e-13,
        absolute_radius_tolerance_m: 1e-19,
        absolute_pressure_tolerance_pa: 1e-8,
        maximum_iterations: 100,
    }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn pressure_changes_radius_volume_and_stiffness_without_changing_gas_or_number() {
    let original = population(0.0);
    for pressure in [25_000.0, 250_000.0] {
        let prepared = original
            .at_equilibrium_pressure(pressure, &controls(), &mut || false)
            .unwrap();
        assert_eq!(prepared.population.liquid(), original.liquid());
        assert_eq!(
            prepared.population.number_density_per_m3(),
            original.number_density_per_m3()
        );
        assert_eq!(
            prepared.population.equilibrium_absolute_pressure_pa(),
            Some(pressure)
        );
        let mut expected_volume = 0.0;
        for ((old, new), root) in original
            .classes()
            .iter()
            .zip(prepared.population.classes())
            .zip(&prepared.report.equilibria)
        {
            assert_eq!(old.bubble(), new.bubble());
            assert_eq!(old.number_density_per_m3(), new.number_density_per_m3());
            let gas = old.bubble().gas();
            let gas_pressure = pressure - original.liquid().saturation_pressure_pa();
            let radius = gas.reference_radius_m()
                * (gas.reference_pressure_pa() / gas_pressure)
                    .powf(1.0 / (3.0 * gas.polytropic_exponent()));
            let linear = new.linearization().unwrap();
            near(linear.equilibrium_radius_m, radius, 3e-13);
            near(
                linear.stiffness_coefficient_per_s2,
                3.0 * gas.polytropic_exponent() * gas_pressure / (1000.0 * radius * radius),
                1e-12,
            );
            let fraction =
                old.number_density_per_m3() * 4.0 * core::f64::consts::PI * radius.powi(3) / 3.0;
            near(new.void_fraction(), fraction, 1e-12);
            expected_volume += fraction;
            assert!(
                root.unwrap().pressure_residual_pa.abs()
                    <= controls().absolute_pressure_tolerance_pa
            );
        }
        near(prepared.population.void_fraction(), expected_volume, 1e-12);
        assert!(prepared.report.total_iterations <= 2 * controls().maximum_iterations);
    }
}

#[test]
fn pressure_cycles_have_no_reference_drift_and_split_bins_keep_the_same_physics() {
    let original = population(0.072);
    let lowered = original
        .at_equilibrium_pressure(40_000.0, &controls(), &mut || false)
        .unwrap();
    let cycled = lowered
        .population
        .at_equilibrium_pressure(180_000.0, &controls(), &mut || false)
        .unwrap();
    let direct = original
        .at_equilibrium_pressure(180_000.0, &controls(), &mut || false)
        .unwrap();
    assert_eq!(cycled, direct);

    let source = original.classes()[0];
    let part =
        BubblePopulationClass::try_new(*source.bubble(), source.number_density_per_m3() / 8.0)
            .unwrap();
    let merged = BubblePopulation::try_new(original.liquid(), vec![source], 1).unwrap();
    let split = BubblePopulation::try_new(original.liquid(), vec![part; 8], 8).unwrap();
    let a = merged
        .at_equilibrium_pressure(40_000.0, &controls(), &mut || false)
        .unwrap();
    let b = split
        .at_equilibrium_pressure(40_000.0, &controls(), &mut || false)
        .unwrap();
    assert_eq!(a.population.void_fraction(), b.population.void_fraction());
    assert_eq!(
        a.population.number_density_per_m3(),
        b.population.number_density_per_m3()
    );
    for member in b.population.classes() {
        assert_eq!(
            member.linearization(),
            a.population.classes()[0].linearization()
        );
    }
}

#[test]
fn continued_population_starts_the_original_nonlinear_models_at_the_actual_equilibria() {
    let original = population(0.072);
    let pressure = 40_000.0;
    let prepared = original
        .at_equilibrium_pressure(pressure, &controls(), &mut || false)
        .unwrap();
    let radii: Vec<_> = prepared
        .population
        .classes()
        .iter()
        .map(|class| class.linearization().unwrap().equilibrium_radius_m)
        .collect();
    let expected_fraction = prepared.population.void_fraction();
    let mut evolving = prepared
        .into_evolving(
            0.0,
            PopulationEvolutionLimits {
                maximum_void_fraction: 0.001,
            },
        )
        .unwrap();
    for ((state, radius), old) in evolving.states().iter().zip(&radii).zip(original.classes()) {
        assert_eq!(state.radius_m(), *radius);
        assert_ne!(state.radius_m(), old.bubble().gas().reference_radius_m());
    }
    let numerics = PopulationEvolutionNumerics {
        per_class: BubbleNumerics {
            relative_tolerance: 1e-9,
            radius_absolute_tolerance_m: 1e-16,
            velocity_absolute_tolerance_m_s: 1e-10,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-13,
            maximum_step_s: 1e-6,
            maximum_attempts: 10_000,
        },
        maximum_total_attempts: 20_000,
    };
    let step = evolving
        .advance_controlled(
            50e-6,
            &numerics,
            &|_| BubblePressureSample {
                absolute_pressure_pa: pressure,
                time_derivative_pa_s: 0.0,
            },
            &mut || false,
        )
        .unwrap();
    near(step.observation.void_fraction, expected_fraction, 2e-11);
    for ((state, radius), (current, old)) in evolving.states().iter().zip(&radii).zip(
        evolving
            .population()
            .classes()
            .iter()
            .zip(original.classes()),
    ) {
        near(state.radius_m(), *radius, 1e-11);
        assert!(state.wall_velocity_m_s().abs() < 1e-8);
        assert_eq!(current.bubble(), old.bubble());
    }
}

#[test]
fn failure_of_any_class_or_final_cancellation_publishes_no_partial_distribution() {
    let original = population(0.072);
    let saved = original.clone();
    let mut polls = 0;
    let good = original
        .at_equilibrium_pressure(40_000.0, &controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(
        original
            .at_equilibrium_pressure(40_000.0, &controls(), &mut || {
                calls += 1;
                calls == polls
            })
            .is_err()
    );
    assert_eq!(original, saved);
    assert_eq!(
        original
            .at_equilibrium_pressure(40_000.0, &controls(), &mut || false)
            .unwrap(),
        good
    );
    assert!(matches!(
        original.at_equilibrium_pressure(-1e6, &controls(), &mut || false),
        Err(BubblePopulationError::Equilibrium { class_index: 0, .. })
    ));
    let tight = EquilibriumControls {
        maximum_iterations: 1,
        ..controls()
    };
    assert!(
        original
            .at_equilibrium_pressure(40_000.0, &tight, &mut || false)
            .is_err()
    );
    assert_eq!(original, saved);

    let material = liquid(0.0);
    let radius: f64 = 1e-4;
    let number = 0.3 / (4.0 * core::f64::consts::PI * radius.powi(3) / 3.0);
    let crowded =
        BubblePopulation::try_new(material, vec![class(material, radius, number, 1.0)], 1).unwrap();
    let next_pressure = material.saturation_pressure_pa() + 98_000.0 / 6.0;
    assert!(matches!(
        crowded.at_equilibrium_pressure(next_pressure, &controls(), &mut || false),
        Err(BubblePopulationError::InvalidInput(_))
    ));
}

#[test]
fn absent_nuclei_do_not_manufacture_equilibria_or_work() {
    let material = liquid(0.072);
    let reference = class(material, 1e-4, 1.0, 1.0);
    let empty_model = KellerMiksisBubble::try_new(
        material,
        BubbleGas::try_new(1e-4, 0.0, 1.0).unwrap(),
        reference.bubble().compressibility_limits(),
    )
    .unwrap();
    let empty = BubblePopulation::try_new(
        material,
        vec![BubblePopulationClass::try_new(empty_model, 0.0).unwrap()],
        1,
    )
    .unwrap();
    let prepared = empty
        .at_equilibrium_pressure(-100_000.0, &controls(), &mut || false)
        .unwrap();
    assert_eq!(prepared.population, empty);
    assert_eq!(prepared.report.equilibria, vec![None]);
    assert_eq!(prepared.report.total_iterations, 0);
    assert_eq!(prepared.report.maximum_absolute_pressure_residual_pa, 0.0);
    assert!(
        empty
            .at_equilibrium_pressure(f64::NAN, &controls(), &mut || false)
            .is_err()
    );
}

#[test]
fn static_pressure_excursions_cannot_bypass_frozen_liquid_compressibility_limits() {
    let material = CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 50.0).unwrap();
    let original = class(material, 1e-4, 1e6, 1.4);
    // rho*c^2 = 2.5 MPa; the authored 0.01 ratio allows at most 25 kPa
    // away from the original 100 kPa operating point, not per query.
    let (first, _) = original
        .at_equilibrium_pressure(120_000.0, &controls(), &mut || false)
        .unwrap();
    let error = first
        .at_equilibrium_pressure(140_000.0, &controls(), &mut || false)
        .unwrap_err();
    assert!(matches!(
        error,
        crate::cavitation::inception::InceptionError::Input(
            CavitationError::PressureRatioExceeded { .. }
        )
    ));
    assert_eq!(first.bubble(), original.bubble());
    assert!(
        original
            .at_equilibrium_pressure(140_000.0, &controls(), &mut || false)
            .is_err()
    );
}
