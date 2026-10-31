//! Public production-crate regressions; no copied EOS or private test modules.
#[path = "common/nonideal_equilibrium.rs"]
mod common;
use common::{MOLAR_MASS, controls, model};
use fs_material::thermochemical_gas::nonideal::equilibrium::{
    PureFluidEquilibriumError as E, PureFluidEquilibriumPhase as Phase,
};
use fs_qty::{Density, Pressure, Temperature};

#[test]
fn g1_conserved_density_energy_recovers_independent_phase_states() {
    let source = model(0.0, [100.0, 250.0]);
    // NumPy cubic roots + SciPy brentq on Maxwell's area equation;
    // u_ideal=2.5 R T/M and original analytical PR departure. Independent
    // density/T endpoints include both sides of the convex phase envelope.
    for (volume, temperature, energy, pressure, quality) in [
        (
            0.0005,
            140.0,
            -101_904.284_402_097_46,
            645_034.343_438_627_5,
            0.301_267_706_509_236_6,
        ),
        (
            0.0005,
            150.0,
            -6705.825_981_239_993,
            1_047_062.532_316_559_9,
            0.493_351_640_620_488_4,
        ),
        (
            0.0005,
            170.0,
            178_185.720_263_928_7,
            2_043_249.266_435_749_8,
            1.0,
        ),
        (
            0.00005,
            170.0,
            -95_436.407_775_541_81,
            2_532_958.544_984_899_5,
            0.0,
        ),
    ] {
        let density = MOLAR_MASS / volume;
        let result = source
            .reconstruct_pure_density_energy(density, energy, controls(), || false)
            .unwrap();
        let state = result.state();
        assert_eq!(state.density_kg_m3().to_bits(), density.to_bits());
        assert!((state.temperature_k() - temperature).abs() < 3e-8);
        assert!((state.pressure_pa() - pressure).abs() < 0.01);
        assert!((state.vapor_mass_fraction() - quality).abs() < 3e-9);
        assert!(result.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
        assert_eq!(result.accepted_increment_j_kg(), None);
        assert!(
            (state.specific_enthalpy_j_kg()
                - state.specific_internal_energy_j_kg()
                - state.pressure_pa() / density)
                .abs()
                < 1e-8
        );
        assert!(result.temperature_evaluations() <= controls().maximum_temperature_evaluations);
        assert!(result.eos_evaluations() <= controls().maximum_eos_evaluations);
    }
}

#[test]
fn g0_isochoric_heating_and_cooling_cross_both_dome_boundaries() {
    let source = model(0.0, [100.0, 250.0]);
    for (volume, heat, vapor) in [
        (0.0005, 280_090.004_666_026_16, true),
        (0.00005, 119_577.435_965_805_97, false),
    ] {
        let density = MOLAR_MASS / volume;
        let initial = source
            .evaluate_pure_equilibrium_density_temperature(
                density,
                Temperature::new(140.0),
                controls(),
                || false,
            )
            .unwrap();
        assert!(matches!(initial.phase(), Phase::Coexistence(_)));
        let heated = source
            .reconstruct_pure_density_internal_energy_increment(
                density,
                Temperature::new(140.0),
                heat,
                controls(),
                || false,
            )
            .unwrap();
        assert!((heated.state().temperature_k() - 170.0).abs() < 3e-8);
        assert!(matches!(heated.state().phase(), Phase::Vapor(_)) == vapor);
        assert!(matches!(heated.state().phase(), Phase::Liquid(_)) != vapor);
        assert!(
            (heated.accepted_increment_j_kg().unwrap() - heat).abs()
                <= controls().energy_tolerance_j_kg
        );
        let cooled = source
            .reconstruct_pure_density_internal_energy_increment(
                density,
                Temperature::new(heated.state().temperature_k()),
                -heat,
                controls(),
                || false,
            )
            .unwrap();
        assert!((cooled.state().temperature_k() - initial.temperature_k()).abs() < 3e-8);
        assert!(matches!(cooled.state().phase(), Phase::Coexistence(_)));
        assert!(
            (cooled.state().vapor_mass_fraction() - initial.vapor_mass_fraction()).abs() < 3e-9
        );
    }
}

#[test]
fn g3_equilibrium_replaces_metastable_homogeneous_states_and_cancels_large_datums() {
    let base = model(0.0, [100.0, 250.0]);
    let density = MOLAR_MASS / 0.0005;
    let initial = base
        .evaluate_pure_equilibrium_density_temperature(
            density,
            Temperature::new(140.0),
            controls(),
            || false,
        )
        .unwrap();
    let homogeneous = base
        .mechanics()
        .evaluate_density_temperature(Temperature::new(140.0), Density::new(density))
        .unwrap();
    assert!(homogeneous.root().pressure_derivatives().1 < 0.0);
    assert!(matches!(initial.phase(), Phase::Coexistence(_)));
    assert!((homogeneous.pressure().value() - initial.pressure_pa()).abs() > 1e5);
    let unshifted = base
        .reconstruct_pure_density_internal_energy_increment(
            density,
            Temperature::new(140.0),
            10.0,
            controls(),
            || false,
        )
        .unwrap();
    for datum in [-1e20, 1e20, 1e28] {
        let shifted = model(datum, [100.0, 250.0]);
        let before = shifted
            .evaluate_pure_equilibrium_density_temperature(
                density,
                Temperature::new(140.0),
                controls(),
                || false,
            )
            .unwrap();
        assert!(matches!(
            shifted.reconstruct_pure_density_energy(
                density,
                before.specific_internal_energy_j_kg(),
                controls(),
                || false
            ),
            Err(E::Invalid(_))
        ));
        if datum == 1e28 {
            // Even loose numerical tolerances cannot choose an arbitrary
            // temperature from an entirely collapsed represented energy span.
            let lower = shifted
                .evaluate_pure_equilibrium_density_temperature(
                    density,
                    controls().temperature_bracket[0],
                    controls(),
                    || false,
                )
                .unwrap();
            let upper = shifted
                .evaluate_pure_equilibrium_density_temperature(
                    density,
                    controls().temperature_bracket[1],
                    controls(),
                    || false,
                )
                .unwrap();
            assert_eq!(
                lower.specific_internal_energy_j_kg(),
                upper.specific_internal_energy_j_kg()
            );
            let mut loose = controls();
            loose.energy_tolerance_j_kg = 1e20;
            loose.temperature_tolerance_k = 1e20;
            assert!(matches!(
                shifted.reconstruct_pure_density_energy(
                    density,
                    before.specific_internal_energy_j_kg(),
                    loose,
                    || false
                ),
                Err(E::Invalid(_))
            ));
        }
        let after = shifted
            .reconstruct_pure_density_internal_energy_increment(
                density,
                Temperature::new(140.0),
                10.0,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(
            after.state().temperature_k(),
            unshifted.state().temperature_k()
        );
        assert_eq!(after.state().phase(), unshifted.state().phase());
        assert_eq!(
            after.accepted_increment_j_kg(),
            unshifted.accepted_increment_j_kg()
        );
        assert!(after.state().temperature_k() > before.temperature_k());
        assert!((after.accepted_increment_j_kg().unwrap() - 10.0).abs() < 1e-6);
    }
}

#[test]
fn g4_domain_unresolved_phase_cancellation_and_whole_call_caps_refuse() {
    let source = model(0.0, [100.0, 250.0]);
    let density = MOLAR_MASS / 0.0005;
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            source
                .reconstruct_pure_density_energy(bad, 0.0, controls(), || false)
                .is_err()
        );
    }
    assert!(matches!(
        source.reconstruct_pure_density_energy(density, 1e9, controls(), || false),
        Err(E::OutsideBracket { .. })
    ));
    let mut work = controls();
    work.maximum_eos_evaluations = 4;
    assert!(matches!(
        source.reconstruct_pure_density_energy(density, 0.0, work, || false),
        Err(E::WorkLimit {
            eos_evaluations: 4,
            ..
        })
    ));
    let mut tight = controls();
    tight.maximum_temperature_evaluations = 2;
    assert!(matches!(
        source.reconstruct_pure_density_energy(density, 0.0, tight, || false),
        Err(E::NotConverged {
            temperature_evaluations: 2,
            ..
        })
    ));
    let mut wrong_phase = controls();
    wrong_phase.pressure_bracket = [Pressure::new(2e6), Pressure::new(4.4e6)];
    assert!(matches!(
        source.reconstruct_pure_density_energy(density, 0.0, wrong_phase, || false),
        Err(E::Phase(_))
    ));
    assert!(
        model(0.0, [140.0, 160.0])
            .reconstruct_pure_density_energy(density, 0.0, controls(), || false)
            .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        source.reconstruct_pure_density_internal_energy_increment(
            density,
            Temperature::new(140.0),
            10.0,
            controls(),
            || {
                polls += 1;
                polls >= 7
            }
        ),
        Err(E::Cancelled)
    );
    assert_eq!(polls, 7);
    let mut strict = controls();
    strict.energy_tolerance_j_kg = 1e-8;
    strict.temperature_tolerance_k = 100.0;
    strict.maximum_temperature_evaluations = 2;
    assert!(matches!(
        source.reconstruct_pure_density_energy(density, 0.123456789, strict, || false),
        Err(E::NotConverged { .. })
    ));
}
