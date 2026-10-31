//! Independent manufactured-mixture caloric identities and conserved closures.
//! These are algebraic validation fixtures, not empirical gas qualifications.
use super::*;
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::nonideal::PengRobinsonSpeciesBinding;
use crate::thermochemical_gas::tests::{model, point};
use fs_thermochem::SpeciesId;

fn bound(source: ThermochemicalGasModel) -> ThermochemicalPengRobinsonModel {
    ThermochemicalPengRobinsonModel::new(
        source,
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap()
}

fn controls() -> GasEntropySolve {
    GasEntropySolve {
        minimum_temperature_k: 280.0,
        maximum_temperature_k: 450.0,
        entropy_tolerance_j_kg_k: 1.0e-9,
        maximum_iterations: 128,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_isentropic_compression_matches_independent_helmholtz_entropy_oracle() {
    let model = bound(model(vec![0.75, 0.25], 0.001));
    let rho0 = 47.889_003_167_633_395;
    let rho1 = rho0 / 0.8;
    let context = point(300.0, 5.0e6).with("material-lot", 12.0).unwrap();
    let initial = model.evaluate_density_temperature(&context, rho0).unwrap();
    // Independent 75-digit Decimal solution, with M=.022 and ideal Cv/R =
    // 2.25+.00075*T. Entropy difference is the Cv/T integral plus
    // R ln((v1-b)/(v0-b)) + a1' I(v1)-a0' I(v0), all divided by M.
    let solution = model
        .reconstruct_density_entropy_change(&context, rho0, rho1, 0.0, controls(), || false)
        .unwrap();
    let state = solution.state();
    near(state.temperature_k(), 331.606_449_235_535_6, 2.0e-9);
    near(state.pressure_pa(), 6_981_183.016_655_789, 1.0e-4);
    assert_eq!(state.density_kg_m3().to_bits(), rho1.to_bits());
    assert_eq!(
        state.source().components()[0]
            .query()
            .axes()
            .get("material-lot"),
        Some(&12.0)
    );
    assert!(solution.residual_j_kg_k().abs() <= controls().entropy_tolerance_j_kg_k);
    near(
        state.specific_entropy_j_kg_k().unwrap(),
        initial.specific_entropy_j_kg_k().unwrap(),
        2.0e-9,
    );
    near(
        model
            .specific_internal_energy_change_between_densities(
                rho0,
                300.0,
                rho1,
                state.temperature_k(),
            )
            .unwrap(),
        24_622.834_659_551_524,
        2.0e-6,
    );
    // A constant gamma from the same initial material misses the PR adiabat.
    assert!(
        (300.0 * (rho1 / rho0).powf(initial.gamma() - 1.0) - state.temperature_k()).abs() > 1.0
    );
}

#[test]
fn g1_entropy_caloric_and_acoustic_derivatives_share_the_density_eos() {
    let model = bound(model(vec![0.75, 0.25], 0.001));
    let (rho, temperature) = (50.0, 300.0);
    let context = point(temperature, 5.0e6);
    let initial = model.evaluate_density_temperature(&context, rho).unwrap();
    let dt = 1.0e-3;
    let entropy_change = model
        .specific_entropy_change_between_densities(rho, temperature - dt, rho, temperature + dt)
        .unwrap();
    near(
        entropy_change / (2.0 * dt),
        initial.specific_heat_cv_j_kg_k() / temperature,
        2.0e-8,
    );
    let epsilon = 1.0e-4;
    let solve = |density| {
        model
            .reconstruct_density_entropy_change(&context, rho, density, 0.0, controls(), || false)
            .unwrap()
    };
    let lower = solve(rho * (1.0 - epsilon));
    let upper = solve(rho * (1.0 + epsilon));
    let dp_drho_s =
        (upper.state().pressure_pa() - lower.state().pressure_pa()) / (2.0 * rho * epsilon);
    assert!((dp_drho_s / initial.sound_speed_m_s().powi(2) - 1.0).abs() < 5.0e-8);
}

#[test]
fn g1_general_entropy_increment_is_reversible_additive_and_reconstructs() {
    let model = bound(model(vec![0.75, 0.25], 0.001));
    let a = (40.0, 300.0);
    let b = (50.0, 330.0);
    let c = (45.0, 370.0);
    let increment = |a: (f64, f64), b: (f64, f64)| {
        model
            .specific_entropy_change_between_densities(a.0, a.1, b.0, b.1)
            .unwrap()
    };
    let target = increment(a, c);
    near(increment(a, b) + increment(b, c), target, 2.0e-12);
    assert_eq!(increment(a, a), 0.0);
    assert_eq!(increment(c, a), -target);
    let solved = model
        .reconstruct_density_entropy_change(
            &point(a.1, 5.0e6),
            a.0,
            c.0,
            target,
            controls(),
            || false,
        )
        .unwrap();
    near(solved.state().temperature_k(), c.1, 2.0e-9);
    // Fixed-density heating uses complete Cv/T, including the PR departure.
    let fixed = model
        .reconstruct_density_entropy_change(
            &point(a.1, 5.0e6),
            a.0,
            a.0,
            1.0e-7,
            controls(),
            || false,
        )
        .unwrap();
    assert!(fixed.state().temperature_k() > a.1);
    assert!(fixed.residual_j_kg_k().abs() <= controls().entropy_tolerance_j_kg_k);
}

fn with_entropy_reference(offset: f64) -> ThermochemicalGasModel {
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};
    let original = model(vec![0.75, 0.25], 0.001);
    let components = original
        .thermodynamics()
        .components()
        .iter()
        .map(|source| {
            let regions = source
                .regions()
                .iter()
                .map(|region| {
                    let mut card = region.card().clone();
                    card.parameters.get_mut("a8").unwrap().value += offset;
                    Nasa9RegionV1::from_card(card).unwrap()
                })
                .collect();
            Nasa9StandardStateModelV1::new(
                source.species().clone(),
                source.molar_mass(),
                source.convention().clone(),
                regions,
            )
            .unwrap()
        })
        .collect();
    let mixture = FrozenIdealGasMixtureModelV1::new(
        components,
        original.thermodynamics().declared_composition().clone(),
    )
    .unwrap();
    ThermochemicalGasModel::new(mixture, original.transport_materials().to_vec(), None).unwrap()
}

#[test]
fn g1_large_species_entropy_references_cancel_before_adiabatic_inversion() {
    let plain = bound(with_entropy_reference(0.0));
    let shifted = bound(with_entropy_reference(1.0e24));
    let context = point(300.0, 5.0e6);
    let solve = |model: &ThermochemicalPengRobinsonModel| {
        model
            .reconstruct_density_entropy_change(&context, 40.0, 50.0, 0.0, controls(), || false)
            .unwrap()
    };
    let a = solve(&plain);
    let b = solve(&shifted);
    assert_eq!(a.state().root(), b.state().root());
    assert_eq!(a.state().source().state(), b.state().source().state());
    assert_eq!(a.residual_j_kg_k(), b.residual_j_kg_k());
    assert_eq!(a.iterations(), b.iterations());
    assert!(b.state().temperature_k() > 300.0);
}

#[test]
fn g0_bad_states_sources_work_bounds_and_cancellation_return_no_candidate() {
    let model = bound(model(vec![0.75, 0.25], 0.001));
    let context = point(300.0, 5.0e6);
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .reconstruct_density_entropy_change(
                    &context,
                    density,
                    50.0,
                    0.0,
                    controls(),
                    || false
                )
                .is_err()
        );
        assert!(
            model
                .reconstruct_density_entropy_change(
                    &context,
                    40.0,
                    density,
                    0.0,
                    controls(),
                    || false
                )
                .is_err()
        );
    }
    for change in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            model
                .reconstruct_density_entropy_change(
                    &context,
                    40.0,
                    50.0,
                    change,
                    controls(),
                    || false
                )
                .is_err()
        );
    }
    assert!(matches!(
        model.reconstruct_density_entropy_change(&context, 40.0, 50.0, 1.0e9, controls(), || false),
        Err(GasEntropyError::OutsideBracket { .. })
    ));
    assert!(matches!(
        model.reconstruct_density_entropy_change(
            &context,
            40.0,
            50.0,
            0.0,
            GasEntropySolve {
                maximum_iterations: 1,
                ..controls()
            },
            || false
        ),
        Err(GasEntropyError::NotConverged { .. })
    ));
    assert!(matches!(
        model.reconstruct_density_entropy_change(&context, 40.0, 50.0, 0.0, controls(), || true),
        Err(GasEntropyError::Cancelled)
    ));
    let mut total_polls = 0;
    model
        .reconstruct_density_entropy_change(&context, 40.0, 50.0, 0.0, controls(), || {
            total_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        model.reconstruct_density_entropy_change(&context, 40.0, 50.0, 0.0, controls(), || {
            polls += 1;
            polls == total_polls
        }),
        Err(GasEntropyError::Cancelled)
    ));
    // The pure constitutive endpoint is valid, but its pressure exceeds the
    // retained transport cards' 10 MPa domain. Closure must refuse publication.
    let change = model
        .specific_entropy_change_between_densities(40.0, 300.0, 100.0, 440.0)
        .unwrap();
    assert!(matches!(
        model.reconstruct_density_entropy_change(&context, 40.0, 100.0, change, controls(), || {
            false
        }),
        Err(GasEntropyError::Material(_))
    ));
    let zero = model
        .reconstruct_density_entropy_change(&context, 40.0, 40.0, 0.0, controls(), || false)
        .unwrap();
    assert_eq!(zero.state().temperature_k(), 300.0);
    assert_eq!(zero.iterations(), 0);
}
