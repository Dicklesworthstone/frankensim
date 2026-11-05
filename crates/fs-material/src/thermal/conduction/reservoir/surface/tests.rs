use super::*;
use std::cell::Cell;

fn config() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g0_constant_film_matches_original_series_owner_in_both_directions_and_equilibrium() {
    let curve = ConductivityCurve::new(vec![[200.0, 4.0], [1000.0, 4.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    let film_g = 5.0;
    let series_g = 1.0 / (0.04 / (0.02 * 4.0) + 1.0 / film_g);
    for [cell, bulk] in [[700.0, 300.0], [300.0, 700.0], [450.0, 450.0]] {
        let response = geometry
            .evaluate_convective_material_to_bath(
                &curve,
                [cell, bulk],
                config(),
                |surface| {
                    Ok::<_, &'static str>((film_g * (surface - bulk), film_g, surface.to_bits()))
                },
                || Ok(()),
            )
            .unwrap();
        let original = geometry
            .evaluate_material_to_bath(&curve, [cell, bulk], 1.0 / film_g, config(), || Ok(()))
            .unwrap();
        let expected = series_g * (cell - bulk);
        close(response.power_to_bath_w, expected, 1e-8);
        close(response.solid_power_to_bath_w, expected, 1e-8);
        close(
            response.surface_temperature_k,
            bulk + expected / film_g,
            1e-8,
        );
        close(response.conductance_w_k, series_g, 1e-10);
        close(
            response.power_to_bath_w,
            original.power_to_reservoir_w,
            1e-8,
        );
        assert_eq!(response.film, response.surface_temperature_k.to_bits());
        assert!(response.heat_residual_w.abs() <= 1e-9);
    }
    let neutral = geometry
        .evaluate_convective_material_to_bath(
            &curve,
            [450.0; 2],
            config(),
            |_| Ok::<_, ()>((0.0, 0.0, "neutral source")),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(neutral.conductance_w_k, 0.0);
    assert_eq!(neutral.power_to_bath_w, 0.0);
    assert_eq!(neutral.iterations, 0);
    assert_eq!(neutral.film, "neutral source");
}

#[test]
fn g1_nonlinear_film_and_crossed_conductivity_knots_recover_manufactured_surfaces() {
    // k(T)=0.01*T, C=0.5 m. Its exact integral gives +/-500 W between
    // 400 and 600 K. G(Ts)=1+0.01*Ts sets the independently chosen baths.
    let curve = ConductivityCurve::new(vec![[200.0, 2.0], [450.0, 4.5], [800.0, 8.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    for (cell, bulk, face, heat) in [
        (600.0, 300.0, 400.0, 500.0),
        (400.0, 600.0 + 500.0 / 7.0, 600.0, -500.0),
    ] {
        let calls = Cell::new(0);
        let response = geometry
            .evaluate_convective_material_to_bath(
                &curve,
                [cell, bulk],
                config(),
                |surface| {
                    calls.set(calls.get() + 1);
                    let g = 1.0 + 0.01 * surface;
                    Ok::<_, ()>((g * (surface - bulk), g, [surface, bulk, g]))
                },
                || Ok(()),
            )
            .unwrap();
        close(response.surface_temperature_k, face, 1e-8);
        close(response.power_to_bath_w, heat, 1e-8);
        close(response.solid_power_to_bath_w, heat, 1e-8);
        close(response.conductance_w_k, heat / (cell - bulk), 1e-10);
        assert_eq!(
            response.film[0].to_bits(),
            response.surface_temperature_k.to_bits()
        );
        assert_eq!(
            response.power_to_bath_w,
            response.film[2] * (response.film[0] - bulk)
        );
        assert!(response.heat_residual_w.abs() <= 1e-9);
        assert!(calls.get() <= config().max_iterations + 2);
    }
}

#[test]
fn g0_g1_only_actual_material_profile_needs_support_and_boundary_roots_stay_exact() {
    let curve = ConductivityCurve::new(vec![[400.0, 4.0], [500.0, 5.0], [600.0, 6.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    for (cell, bulk, face, heat, film_g) in [
        (550.0, 250.0, 450.0, 250.0, 1.25),
        (450.0, 750.0, 550.0, -250.0, 1.25),
        (550.0, 250.0, 400.0, 356.25, 2.375),
        (450.0, 750.0, 600.0, -393.75, 2.625),
    ] {
        let response = geometry
            .evaluate_convective_material_to_bath(
                &curve,
                [cell, bulk],
                config(),
                |surface| Ok::<_, ()>((film_g * (surface - bulk), film_g, surface)),
                || Ok(()),
            )
            .unwrap();
        close(response.surface_temperature_k, face, 1e-8);
        close(response.power_to_bath_w, heat, 1e-8);
        if face == 400.0 || face == 600.0 {
            assert_eq!(response.surface_temperature_k, face);
            assert_eq!(response.iterations, 0);
        }
        let unsupported = geometry.evaluate_convective_material_to_bath(
            &curve,
            [cell, bulk],
            config(),
            |surface| Ok::<_, ()>((100.0 * (surface - bulk), 100.0, ())),
            || Ok(()),
        );
        assert!(matches!(
            unsupported,
            Err(ConvectiveLayerError::Thermal(
                DepositionError::MaterialDomain
            ))
        ));
    }
}

#[test]
fn g1_thin_metal_retains_flux_accuracy_without_reconstructing_from_rounded_surface() {
    let curve = ConductivityCurve::new(vec![[200.0, 200.0], [1000.0, 200.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.1, 0.0001).unwrap();
    let controls = LayerConductionConfig {
        absolute_heat_tolerance_w: 1e-12,
        ..config()
    };
    for cell in [300.003, 299.997] {
        let bulk = 300.0;
        let response = geometry
            .evaluate_convective_material_to_bath(
                &curve,
                [cell, bulk],
                controls,
                |surface| Ok::<_, ()>((0.05 * (surface - bulk), 0.05, ())),
                || Ok(()),
            )
            .unwrap();
        let expected = (cell - bulk) / (20.0 + 0.0001 / (0.1 * 200.0));
        close(response.power_to_bath_w, expected, 1e-12);
        close(response.solid_power_to_bath_w, expected, 1e-12);
        assert!(response.heat_residual_w.abs() <= 1e-12);
        assert_ne!(response.surface_temperature_k, cell);
        let rounded_temperature_power =
            0.1 / 0.0001 * 200.0 * (cell - response.surface_temperature_k);
        assert!((rounded_temperature_power - expected).abs() > 1e-11);
    }
    // A still thinner layer cannot represent its nonzero face drop at all;
    // preserve the original inverse's refusal instead of flooring that drop.
    let unrepresentable = ReservoirLayerGeometry::planar(1.0, 1e-12)
        .unwrap()
        .evaluate_convective_material_to_bath(
            &curve,
            [300.003, 300.0],
            controls,
            |surface| Ok::<_, ()>((0.05 * (surface - 300.0), 0.05, ())),
            || Ok(()),
        );
    assert!(matches!(
        unrepresentable,
        Err(ConvectiveLayerError::Thermal(DepositionError::Arithmetic))
    ));
}

#[test]
fn g0_g4_source_errors_cancellation_and_exhausted_surface_budget_remain_typed() {
    #[derive(Debug, PartialEq)]
    enum SourceFailure {
        MissingProperty(u8),
    }
    let curve = ConductivityCurve::new(vec![[200.0, 2.0], [800.0, 8.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    let refused = geometry.evaluate_convective_material_to_bath(
        &curve,
        [600.0, 300.0],
        config(),
        |_| Err::<(f64, f64, ()), _>(SourceFailure::MissingProperty(17)),
        || Ok(()),
    );
    assert!(matches!(
        refused,
        Err(ConvectiveLayerError::Film(SourceFailure::MissingProperty(
            17
        )))
    ));
    let calls = Cell::new(0);
    let cancelled = geometry.evaluate_convective_material_to_bath(
        &curve,
        [600.0, 300.0],
        config(),
        |surface| {
            calls.set(calls.get() + 1);
            Ok::<_, ()>((surface - 300.0, 1.0, ()))
        },
        || {
            if calls.get() > 0 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(ConvectiveLayerError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(calls.get(), 1);
    let budget = LayerConductionConfig {
        max_iterations: 1,
        absolute_heat_tolerance_w: 0.0,
        relative_heat_tolerance: 0.0,
    };
    let incomplete = geometry.evaluate_convective_material_to_bath(
        &curve,
        [600.0, 300.0],
        budget,
        |surface| {
            let g = 1.0 + 0.01 * surface;
            Ok::<_, ()>((g * (surface - 300.0), g, ()))
        },
        || Ok(()),
    );
    assert!(
        matches!(incomplete, Err(ConvectiveLayerError::SurfaceBalance {
        residual_w, iterations: 1,
    }) if residual_w.is_finite() && residual_w.abs() > 1.0)
    );
}

#[test]
fn g0_driven_zero_or_active_film_cannot_become_a_passive_conducting_layer() {
    let curve = ConductivityCurve::new(vec![[200.0, 4.0], [800.0, 4.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    for (power, g) in [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (-1.0, 1.0)] {
        let result = geometry.evaluate_convective_material_to_bath(
            &curve,
            [600.0, 300.0],
            config(),
            |_| Ok::<_, ()>((power, g, ())),
            || Ok(()),
        );
        assert!(matches!(
            result,
            Err(ConvectiveLayerError::Thermal(
                DepositionError::InvalidInput(_)
            ))
        ));
    }
    let equilibrium_power = geometry.evaluate_convective_material_to_bath(
        &curve,
        [450.0; 2],
        config(),
        |_| Ok::<_, ()>((1.0, 1.0, ())),
        || Ok(()),
    );
    assert!(matches!(
        equilibrium_power,
        Err(ConvectiveLayerError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
}
