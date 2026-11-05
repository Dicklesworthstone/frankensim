use super::*;
use std::cell::Cell;

fn controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g1_three_reservoir_oracle_includes_opposing_drives_neutral_body_and_insulation() {
    let curve = ConductivityCurve::new(vec![[200.0, 4.0], [1000.0, 4.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    let gs = 2.0;
    for (cell, baths, g) in [
        (700.0, [300.0, 900.0], [5.0, 3.0]),
        (300.0, [700.0, 200.0], [5.0, 3.0]),
        (450.0, [450.0, 900.0], [5.0, 3.0]),
        (450.0, [450.0, 200.0], [5.0, 3.0]),
        (450.0, [450.0, 450.0], [5.0, 3.0]),
        (500.0, [300.0, 700.0], [5.0, 5.0]),
        (700.0, [300.0, 900.0], [5.0, 0.0]),
        (300.0, [700.0, 200.0], [5.0, 0.0]),
    ] {
        let response = geometry
            .evaluate_material_to_two_baths(
                &curve,
                cell,
                baths,
                [200.0, 1000.0],
                controls(),
                |ts| {
                    Ok::<_, ()>((
                        [g[0] * (ts - baths[0]), g[1] * (ts - baths[1])],
                        g,
                        ts.to_bits(),
                    ))
                },
                || Ok(()),
            )
            .unwrap();
        let total = gs + g[0] + g[1];
        let ts = (gs * cell + g[0] * baths[0] + g[1] * baths[1]) / total;
        near(response.surface_temperature_k, ts, 1e-9);
        near(response.solid_power_to_surface_w, gs * (cell - ts), 1e-8);
        for i in 0..2 {
            near(response.powers_to_baths_w[i], g[i] * (ts - baths[i]), 1e-8);
            near(
                response.material_bath_conductances_w_k[i],
                gs * g[i] / total,
                1e-12,
            );
        }
        near(
            response.bath_pair_conductance_w_k,
            g[0] * g[1] / total,
            1e-12,
        );
        assert!(response.heat_residual_w.abs() <= 1e-8);
        assert_eq!(response.film, response.surface_temperature_k.to_bits());
        if g == [5.0, 5.0] {
            assert_eq!(response.surface_temperature_k, cell);
            assert_eq!(response.solid_power_to_surface_w, 0.0);
            assert_eq!(response.powers_to_baths_w, [1000.0, -1000.0]);
            assert_eq!(response.heat_residual_w, 0.0);
            assert!(response.bath_pair_conductance_w_k > 0.0);
        }
        if cell == baths[0] && cell != baths[1] {
            assert_ne!(response.powers_to_baths_w[0], 0.0);
            assert_ne!(response.solid_power_to_surface_w, 0.0);
        }
        if g[1] == 0.0 {
            let original = geometry
                .evaluate_convective_material_to_bath(
                    &curve,
                    [cell, baths[0]],
                    controls(),
                    |ts| Ok::<_, ()>((g[0] * (ts - baths[0]), g[0], ())),
                    || Ok(()),
                )
                .unwrap();
            near(
                response.powers_to_baths_w[0],
                original.power_to_bath_w,
                1e-8,
            );
            near(
                response.material_bath_conductances_w_k[0],
                original.conductance_w_k,
                1e-12,
            );
            assert_eq!(response.bath_pair_conductance_w_k, 0.0);
        }
    }
}

#[test]
fn g1_affine_material_and_nonlinear_convection_radiation_recover_manufactured_surfaces() {
    // C=0.5 m and k=.01*T give Qs=.0025*(Tc²-Ts²). Select Ts first,
    // then derive a fourth-power reservoir analytically, without a root solve.
    let curve = ConductivityCurve::new(vec![[200.0, 2.0], [450.0, 4.5], [1000.0, 10.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    let radiation = 5e-8;
    for (cell, surface, gas) in [(600.0_f64, 400.0_f64, 450.0), (400.0, 600.0, 500.0)] {
        let solid = 0.0025 * (cell * cell - surface * surface);
        let gas_heat = (1.0 + 0.01 * surface) * (surface - gas);
        let radiant = (surface.powi(4) - (solid - gas_heat) / radiation)
            .sqrt()
            .sqrt();
        let calls = Cell::new(0);
        let response = geometry
            .evaluate_material_to_two_baths(
                &curve,
                cell,
                [gas, radiant],
                [200.0, 1000.0],
                controls(),
                |ts| {
                    calls.set(calls.get() + 1);
                    let gg = 1.0 + 0.01 * ts;
                    let gr = radiation * (ts + radiant) * (ts * ts + radiant * radiant);
                    Ok::<_, ()>(([gg * (ts - gas), gr * (ts - radiant)], [gg, gr], ts))
                },
                || Ok(()),
            )
            .unwrap();
        near(response.surface_temperature_k, surface, 1e-8);
        near(response.solid_power_to_surface_w, solid, 1e-8);
        near(response.powers_to_baths_w[0], gas_heat, 1e-8);
        near(response.powers_to_baths_w[1], solid - gas_heat, 1e-8);
        assert!(response.powers_to_baths_w[0] * response.powers_to_baths_w[1] < 0.0);
        assert_eq!(response.film, response.surface_temperature_k);
        assert!(calls.get() <= controls().max_iterations + 2);
    }
}

#[test]
fn g0_g1_supported_surface_need_not_include_either_bath_or_the_material_cell() {
    let curve = ConductivityCurve::new(vec![[400.0, 4.0], [600.0, 6.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    for span in [[425.0, 475.0], [450.0, 450.0], [455.0, 475.0]] {
        let result = geometry.evaluate_material_to_two_baths(
            &curve,
            550.0,
            [250.0, 850.0],
            span,
            controls(),
            |ts| {
                assert!((span[0]..=span[1]).contains(&ts));
                Ok::<_, ()>(([2.0 * (ts - 250.0), 0.375 * (ts - 850.0)], [2.0, 0.375], ()))
            },
            || Ok(()),
        );
        if span[0] > 450.0 {
            assert!(matches!(
                result,
                Err(ConvectiveLayerError::Thermal(
                    DepositionError::MaterialDomain
                ))
            ));
        } else {
            let response = result.unwrap();
            near(response.surface_temperature_k, 450.0, 1e-8);
            near(response.solid_power_to_surface_w, 250.0, 1e-8);
            if span[0] == span[1] {
                assert_eq!(response.surface_temperature_k, 450.0);
                assert_eq!(response.iterations, 0);
            }
        }
    }
}

#[test]
fn g1_thin_metal_retains_each_bath_power_and_material_integral_accuracy() {
    let curve = ConductivityCurve::new(vec![[200.0, 200.0], [1000.0, 200.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.1, 0.0001).unwrap();
    let config = LayerConductionConfig {
        absolute_heat_tolerance_w: 1e-12,
        ..controls()
    };
    for cell in [300.003, 299.997] {
        let response = geometry
            .evaluate_material_to_two_baths(
                &curve,
                cell,
                [300.0; 2],
                [200.0, 1000.0],
                config,
                |ts| Ok::<_, ()>(([0.05 * (ts - 300.0), 0.03 * (ts - 300.0)], [0.05, 0.03], ())),
                || Ok(()),
            )
            .unwrap();
        let expected = (cell - 300.0) / (1.0 / 0.08 + 0.0001 / (0.1 * 200.0));
        near(response.solid_power_to_surface_w, expected, 1e-12);
        near(response.powers_to_baths_w[0], expected * 0.05 / 0.08, 1e-12);
        near(response.powers_to_baths_w[1], expected * 0.03 / 0.08, 1e-12);
        assert!(response.heat_residual_w.abs() <= 1e-12);
        let reconstructed = 200000.0 * (cell - response.surface_temperature_k);
        assert!((reconstructed - expected).abs() > 1e-11);
    }
}

#[test]
fn g0_g4_bath_passivity_original_errors_cancellation_and_work_budget_remain_authoritative() {
    let curve = ConductivityCurve::new(vec![[200.0, 2.0], [1000.0, 10.0]]).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.02, 0.04).unwrap();
    for (power, g) in [(1.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (1.0, 1.0)] {
        let result = geometry.evaluate_material_to_two_baths(
            &curve,
            600.0,
            [300.0; 2],
            [200.0, 1000.0],
            controls(),
            |_| Ok::<_, ()>(([power, 0.0], [g, 0.0], ())),
            || Ok(()),
        );
        assert!(matches!(
            result,
            Err(ConvectiveLayerError::Thermal(
                DepositionError::InvalidInput(_)
            ))
        ));
    }
    let refused = geometry.evaluate_material_to_two_baths(
        &curve,
        600.0,
        [300.0; 2],
        [200.0, 1000.0],
        controls(),
        |_| Err::<([f64; 2], [f64; 2], ()), _>("original missing source"),
        || Ok(()),
    );
    assert!(matches!(
        refused,
        Err(ConvectiveLayerError::Film("original missing source"))
    ));
    let calls = Cell::new(0);
    let cancelled = geometry.evaluate_material_to_two_baths(
        &curve,
        600.0,
        [300.0; 2],
        [200.0, 1000.0],
        controls(),
        |ts| {
            calls.set(calls.get() + 1);
            Ok::<_, ()>(([ts - 300.0, 0.0], [1.0, 0.0], ()))
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
    let incomplete = geometry.evaluate_material_to_two_baths(
        &curve,
        600.0,
        [300.0, 800.0],
        [200.0, 1000.0],
        LayerConductionConfig {
            max_iterations: 1,
            absolute_heat_tolerance_w: 0.0,
            relative_heat_tolerance: 0.0,
        },
        |ts| {
            let gg = 1.0 + 0.01 * ts;
            let gr = 1e-8 * (ts + 800.0) * (ts * ts + 640000.0);
            Ok::<_, ()>(([gg * (ts - 300.0), gr * (ts - 800.0)], [gg, gr], ()))
        },
        || Ok(()),
    );
    assert!(
        matches!(incomplete, Err(ConvectiveLayerError::SurfaceBalance { iterations: 1, residual_w }) if residual_w.is_finite() && residual_w != 0.0)
    );
}
