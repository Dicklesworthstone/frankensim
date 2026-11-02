//! Exact radial Fourier balances, actual interface areas and source-supported profiles.
use super::*;
use fs_conduction::layered::{NonlinearGenerationOptions, RadialLayeredWall};

fn controls() -> NonlinearGenerationOptions {
    NonlinearGenerationOptions {
        temperature_interval_k: [280.0, 700.0],
        solver: options(),
    }
}

fn shell(cylindrical: bool, wall: LayeredWall, inner: f64) -> RadialLayeredWall {
    if cylindrical {
        RadialLayeredWall::cylindrical(wall, inner, 0.5).unwrap()
    } else {
        RadialLayeredWall::spherical(wall, inner).unwrap()
    }
}

fn affine_wall(thickness: f64) -> LayeredWall {
    LayeredWall::new(
        vec![curve(
            "radial k=T/300",
            thickness,
            vec![(280.0, 280.0 / 300.0), (700.0, 700.0 / 300.0)],
        )],
        vec![],
        [0.0; 2],
    )
    .unwrap()
}

// Independent dimensional integrals for moderate radii, with no source-weight
// or inverse-volume implementation shared with the production geometry leaf.
fn integrals(cylinder: bool, a: f64, b: f64, reference: f64) -> [f64; 3] {
    if cylinder {
        let log = (b / a).ln();
        [
            reference * log,
            (b * b - a * a) / (2.0 * reference),
            (b * b - a * a) / 4.0 - a * a * log / 2.0,
        ]
    } else {
        [
            reference * reference * (1.0 / a - 1.0 / b),
            (b.powi(3) - a.powi(3)) / (3.0 * reference * reference),
            (b - a).powi(2) * (b + 2.0 * a) / (6.0 * b),
        ]
    }
}

#[test]
fn radial_zero_source_recovers_logarithmic_and_inverse_radius_resistance() {
    for cylinder in [true, false] {
        let model = shell(cylinder, affine_wall(1.5), 0.5);
        let answer = with_cx(|cx| {
            model.solve_temperature_dependent_with_generation(
                cx,
                [600.0, 300.0],
                &[0.0],
                controls(),
            )
        })
        .unwrap();
        let area = model.reference_area_m2();
        let [distance, _, _] = integrals(cylinder, 0.5, 2.0, 0.5);
        // integral_300^600 T/300 dT = 450 W/m.
        let q = 450.0 / distance;
        let row = &answer.layer_temperatures[0];
        close(row.radial_heat_rates_w[0], q * area);
        close(row.radial_heat_rates_w[1], q * area);
        close(answer.reservoir_heat_rates_w[0], -q * area);
        close(answer.reservoir_heat_rates_w[1], q * area);
        close(answer.generated_heat_rate_w, 0.0);
        for x in [0.0, 0.2, 0.75, 1.5] {
            let [d, _, _] = integrals(cylinder, 0.5, 0.5 + x, 0.5);
            close(
                with_cx(|cx| answer.temperature_at(cx, 0, x)).unwrap(),
                (600.0_f64.powi(2) - 600.0 * q * d).sqrt(),
            );
            let ratio = if cylinder {
                (0.5 + x) / 0.5
            } else {
                ((0.5 + x) / 0.5_f64).powi(2)
            };
            close(
                with_cx(|cx| answer.heat_flux_at(cx, 0, x)).unwrap(),
                q / ratio,
            );
        }
        assert_eq!(answer.wall(), model.wall());
    }
}

#[test]
fn radial_generation_recovers_nonlinear_hot_and_cold_spots_and_full_energy() {
    for cylinder in [true, false] {
        let model = shell(cylinder, affine_wall(2.0), 1.0);
        let [d, volume, moment] = integrals(cylinder, 1.0, 3.0, 1.0);
        for source in [60.0, -60.0] {
            let reservoir: f64 = if source > 0.0 { 300.0 } else { 600.0 };
            let inlet = -source * moment / d;
            let outlet = inlet + source * volume;
            let radius = if cylinder {
                (1.0 - 2.0 * inlet / source).sqrt()
            } else {
                (1.0 - 3.0 * inlet / source).cbrt()
            };
            let [di, _, bi] = integrals(cylinder, 1.0, radius, 1.0);
            let extremum = (reservoir.powi(2) - 600.0 * (inlet * di + source * bi)).sqrt();
            let answer = with_cx(|cx| {
                model.solve_temperature_dependent_with_generation(
                    cx,
                    [reservoir; 2],
                    &[source],
                    controls(),
                )
            })
            .unwrap();
            let area = model.reference_area_m2();
            let row = &answer.layer_temperatures[0];
            close(row.radial_heat_rates_w[0], inlet * area);
            close(row.radial_heat_rates_w[1], outlet * area);
            close(answer.generated_heat_rate_w, source * volume * area);
            close(answer.energy_residual_w, 0.0);
            let actual = if source > 0.0 {
                [row.maximum_radius_m, row.maximum_temperature_k]
            } else {
                [row.minimum_radius_m, row.minimum_temperature_k]
            };
            close(actual[0], radius);
            close(actual[1], extremum);
            assert!(
                with_cx(|cx| answer.heat_flux_at(cx, 0, radius - 1.0))
                    .unwrap()
                    .abs()
                    < 1e-8
            );
            for x in [0.25, 0.7, 1.25, 1.8] {
                let [dx, _, bx] = integrals(cylinder, 1.0, 1.0 + x, 1.0);
                close(
                    with_cx(|cx| answer.temperature_at(cx, 0, x)).unwrap(),
                    (reservoir.powi(2) - 600.0 * (inlet * dx + source * bx)).sqrt(),
                );
            }
            if source > 0.0 {
                // Supported faces at 300 K cannot admit an unsupported interior maximum.
                let mut limited = controls();
                limited.temperature_interval_k[1] = extremum - 1.0;
                assert!(
                    with_cx(|cx| model.solve_temperature_dependent_with_generation(
                        cx,
                        [reservoir; 2],
                        &[source],
                        limited,
                    ))
                    .is_err()
                );
            }
        }
    }
}

#[test]
fn multilayer_generation_uses_each_contact_and_film_actual_area() {
    for cylinder in [true, false] {
        let wall = LayeredWall::new(
            vec![
                curve(
                    "heater",
                    0.5,
                    vec![(280.0, 280.0 / 300.0), (700.0, 700.0 / 300.0)],
                ),
                PlanarThermalLayer::new("sink", 0.75, ConductivityTable::declared(2.0).unwrap())
                    .unwrap(),
            ],
            vec![0.4],
            [0.2, 0.6],
        )
        .unwrap();
        let model = shell(cylinder, wall.clone(), 1.0);
        let [da, va, ba] = integrals(cylinder, 1.0, 1.5, 1.0);
        let [db, vb, bb] = integrals(cylinder, 1.5, 2.25, 1.0);
        let ratio = |r: f64| if cylinder { r } else { r * r };
        let inlet = -2.0;
        let junction = inlet + 20.0 * va;
        let outlet = junction - 5.0 * vb;
        let t0 = 300.0 - inlet * 0.2;
        let t1 = (t0 * t0 - 600.0 * (inlet * da + 20.0 * ba)).sqrt();
        let t2 = t1 - junction / ratio(1.5) * 0.4;
        let t3 = t2 - (junction * db - 5.0 * bb) / 2.0;
        let reservoir = t3 - outlet / ratio(2.25) * 0.6;
        let answer = with_cx(|cx| {
            model.solve_temperature_dependent_with_generation(
                cx,
                [300.0, reservoir],
                &[20.0, -5.0],
                controls(),
            )
        })
        .unwrap();
        let area = model.reference_area_m2();
        let a = &answer.layer_temperatures[0];
        let b = &answer.layer_temperatures[1];
        for (actual, expected) in a
            .face_temperatures_k
            .into_iter()
            .chain(b.face_temperatures_k)
            .zip([t0, t1, t2, t3])
        {
            close(actual, expected);
        }
        close(a.radial_heat_rates_w[1], junction * area);
        close(b.radial_heat_rates_w[0], junction * area);
        close(a.radial_heat_fluxes_w_per_m2[1], junction / ratio(1.5));
        close(b.radial_heat_fluxes_w_per_m2[1], outlet / ratio(2.25));
        close(a.face_areas_m2[1], area * ratio(1.5));
        close(b.face_areas_m2[1], area * ratio(2.25));
        close(answer.reservoir_heat_rates_w[0], -inlet * area);
        close(answer.reservoir_heat_rates_w[1], outlet * area);
        close(answer.generated_heat_rate_w, (20.0 * va - 5.0 * vb) * area);
        close(answer.energy_residual_w, 0.0);
        assert_eq!(answer.wall(), &wall);
        assert_eq!(answer.volumetric_generation_w_m3(), &[20.0, -5.0]);
    }
}

#[test]
fn radial_subdivision_length_scaling_and_thin_shell_limit_preserve_physics() {
    let single = affine_wall(1.0);
    let split = LayeredWall::new(
        vec![
            curve(
                "a",
                0.25,
                single.layers()[0].conductivity().knots().to_vec(),
            ),
            curve(
                "b",
                0.75,
                single.layers()[0].conductivity().knots().to_vec(),
            ),
        ],
        vec![0.0],
        [0.0; 2],
    )
    .unwrap();
    for cylinder in [true, false] {
        let full = shell(cylinder, single.clone(), 2.0);
        let divided = shell(cylinder, split.clone(), 2.0);
        let a = with_cx(|cx| {
            full.solve_temperature_dependent_with_generation(
                cx,
                [300.0, 310.0],
                &[30.0],
                controls(),
            )
        })
        .unwrap();
        let b = with_cx(|cx| {
            divided.solve_temperature_dependent_with_generation(
                cx,
                [300.0, 310.0],
                &[30.0; 2],
                controls(),
            )
        })
        .unwrap();
        for side in 0..2 {
            close(
                a.reservoir_heat_rates_w[side],
                b.reservoir_heat_rates_w[side],
            );
        }
        close(a.generated_heat_rate_w, b.generated_heat_rate_w);
        close(
            with_cx(|cx| a.temperature_at(cx, 0, 0.6)).unwrap(),
            with_cx(|cx| b.temperature_at(cx, 1, 0.35)).unwrap(),
        );
        if cylinder {
            let long = RadialLayeredWall::cylindrical(single.clone(), 2.0, 1.5).unwrap();
            let c = with_cx(|cx| {
                long.solve_temperature_dependent_with_generation(
                    cx,
                    [300.0, 310.0],
                    &[30.0],
                    controls(),
                )
            })
            .unwrap();
            assert_eq!(
                a.layer_temperatures[0].face_temperatures_k,
                c.layer_temperatures[0].face_temperatures_k
            );
            close(c.generated_heat_rate_w, 3.0 * a.generated_heat_rate_w);
            close(
                c.reservoir_heat_rates_w[0],
                3.0 * a.reservoir_heat_rates_w[0],
            );
        }
        let coating = LayeredWall::new(
            vec![
                PlanarThermalLayer::new("coating", 1.0, ConductivityTable::declared(2.0).unwrap())
                    .unwrap(),
            ],
            vec![],
            [0.0; 2],
        )
        .unwrap();
        let thin = shell(cylinder, coating, 1e6);
        let answer = with_cx(|cx| {
            thin.solve_temperature_dependent_with_generation(cx, [300.0; 2], &[10.0], controls())
        })
        .unwrap();
        let peak = &answer.layer_temperatures[0];
        assert!((peak.maximum_temperature_k - 300.625).abs() < 1e-7);
        assert!((peak.maximum_radius_m - 1e6 - 0.5).abs() < 2e-7);
    }
}

#[test]
fn radial_invalid_geometry_sources_budget_and_cancellation_refuse() {
    let wall = affine_wall(1.0);
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e20] {
        assert!(RadialLayeredWall::spherical(wall.clone(), radius).is_err());
    }
    assert!(RadialLayeredWall::cylindrical(wall.clone(), 1.0, 0.0).is_err());
    let model = shell(true, wall, 1.0);
    let before = model.clone();
    for source in [vec![], vec![f64::NAN], vec![1.0, 2.0]] {
        assert!(
            with_cx(|cx| model.solve_temperature_dependent_with_generation(
                cx,
                [300.0, 500.0],
                &source,
                controls(),
            ))
            .is_err()
        );
    }
    let mut limited = controls();
    limited.solver.max_iterations = 1;
    assert!(
        with_cx(|cx| model.solve_temperature_dependent_with_generation(
            cx,
            [300.0, 500.0],
            &[12.0],
            limited,
        ))
        .is_err()
    );
    assert!(matches!(
        with_cancelled_cx(|cx| model.solve_temperature_dependent_with_generation(
            cx,
            [300.0, 500.0],
            &[12.0],
            controls(),
        )),
        Err(ConductionError::Cancelled { .. })
    ));
    assert_eq!(model, before);
}
