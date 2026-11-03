//! Independent solid-centre balances, nonlinear profiles and radial interfaces.
use super::*;
use fs_conduction::layered::{NonlinearGenerationOptions, SolidRadialLayeredBody};

const PI: f64 = core::f64::consts::PI;

fn controls() -> NonlinearGenerationOptions {
    NonlinearGenerationOptions {
        temperature_interval_k: [200.0, 800.0],
        solver: options(),
    }
}
fn affine_layer(name: &str, thickness: f64) -> PlanarThermalLayer {
    curve(
        name,
        thickness,
        vec![(200.0, 200.0 / 300.0), (800.0, 800.0 / 300.0)],
    )
}
fn solid(cylinder: bool, wall: LayeredWall, length: f64) -> SolidRadialLayeredBody {
    if cylinder {
        SolidRadialLayeredBody::cylindrical(wall, length).unwrap()
    } else {
        SolidRadialLayeredBody::spherical(wall).unwrap()
    }
}
fn physical_area(cylinder: bool, r: f64, length: f64) -> f64 {
    if cylinder {
        2.0 * PI * r * length
    } else {
        4.0 * PI * r * r
    }
}
fn physical_volume(cylinder: bool, r: f64, length: f64) -> f64 {
    if cylinder {
        PI * r * r * length
    } else {
        4.0 * PI * r.powi(3) / 3.0
    }
}
// Independent integrals for the shell 1..r, referenced to its inner area:
// D=int dr/a and B=int int(a du)/a dr, with a=r for cylinders and r² for spheres.
fn shell_moments(cylinder: bool, r: f64) -> (f64, f64) {
    if cylinder {
        (r.ln(), (r * r - 1.0) / 4.0 - r.ln() / 2.0)
    } else {
        (1.0 - 1.0 / r, (r * r - 3.0 + 2.0 / r) / 6.0)
    }
}

#[test]
fn generated_solid_core_and_shell_keep_separate_source_support() {
    for cylinder in [true, false] {
        let dimension = if cylinder { 2.0 } else { 3.0 };
        let (distance, _) = shell_moments(cylinder, 2.0);
        let cold_slope = 100.0 * distance / 35_000.0;
        let make = |core_maximum: f64| {
            LayeredWall::new(
                vec![
                    curve(
                        "hot core",
                        1.0,
                        vec![(550.0, 0.55), (core_maximum, core_maximum / 1000.0)],
                    ),
                    curve(
                        "cold shell",
                        1.0,
                        vec![(250.0, 250.0 * cold_slope), (450.0, 450.0 * cold_slope)],
                    ),
                ],
                vec![2.0],
                [0.0; 2],
            )
            .unwrap()
        };
        let model = solid(cylinder, make(750.0), 3.0);
        let answer = with_cx(|cx| {
            model.solve_temperature_dependent_with_generation(
                cx,
                300.0,
                &[dimension * 100.0, 0.0],
                controls(),
            )
        })
        .unwrap();
        // q at the unit-radius core surface is 100 W/m2. The source
        // integral to the centre is q*r/2=50 W/m, while the contact
        // contributes 200 K and the shell's Kirchhoff drop gives 100 K.
        close(
            answer.layer_temperatures[0].face_temperatures_k[0],
            460_000.0_f64.sqrt(),
        );
        close(answer.layer_temperatures[0].face_temperatures_k[1], 600.0);
        close(answer.layer_temperatures[1].face_temperatures_k[0], 400.0);
        close(answer.layer_temperatures[1].face_temperatures_k[1], 300.0);
        close(
            answer.reservoir_heat_rate_w,
            100.0 * model.reference_area_m2(),
        );
        close(answer.energy_residual_w, 0.0);
        // The actual 678 K centre must still be refused by a 670 K source
        // ceiling, even though the interface is supported at 600 K.
        let narrow = solid(cylinder, make(670.0), 3.0);
        assert!(
            with_cx(|cx| narrow.solve_temperature_dependent_with_generation(
                cx,
                300.0,
                &[dimension * 100.0, 0.0],
                controls(),
            ))
            .is_err()
        );
    }
}

#[test]
fn solid_core_affine_conductivity_has_regular_signed_exact_profiles() {
    for cylinder in [true, false] {
        let m = if cylinder { 2.0 } else { 3.0 };
        let wall = LayeredWall::new(vec![affine_layer("core", 2.0)], vec![], [0.0; 2]).unwrap();
        let model = solid(cylinder, wall, 3.0);
        for source in [-60.0, 60.0] {
            for reservoir in [300.0_f64, 600.0] {
                let answer = with_cx(|cx| {
                    model.solve_temperature_dependent_with_generation(
                        cx,
                        reservoir,
                        &[source],
                        controls(),
                    )
                })
                .unwrap();
                let heat = source * physical_volume(cylinder, 2.0, 3.0);
                close(answer.reservoir_heat_rate_w, heat);
                close(answer.generated_heat_rate_w, heat);
                close(answer.energy_residual_w, 0.0);
                let row = &answer.layer_temperatures[0];
                assert_eq!(row.face_radii_m, [0.0, 2.0]);
                assert_eq!(row.face_areas_m2[0], 0.0);
                assert_eq!(row.radial_heat_rates_w[0], 0.0);
                assert_eq!(row.radial_heat_fluxes_w_per_m2[0], 0.0);
                close(row.face_areas_m2[1], physical_area(cylinder, 2.0, 3.0));
                close(row.radial_heat_rates_w[1], heat);
                close(row.radial_heat_fluxes_w_per_m2[1], source * 2.0 / m);
                for r in [0.0_f64, 0.25, 0.9, 1.75, 2.0] {
                    let expected_squared =
                        reservoir.powi(2) + 600.0 * source * (4.0 - r * r) / (2.0 * m);
                    close(
                        with_cx(|cx| answer.temperature_at(cx, 0, r)).unwrap(),
                        expected_squared.sqrt(),
                    );
                    close(
                        with_cx(|cx| answer.heat_flux_at(cx, 0, r)).unwrap(),
                        source * r / m,
                    );
                }
                if source > 0.0 {
                    close(row.maximum_radius_m, 0.0);
                    close(row.minimum_radius_m, 2.0);
                } else {
                    close(row.minimum_radius_m, 0.0);
                    close(row.maximum_radius_m, 2.0);
                }
                assert!(answer.maximum_heat_rate_residual_w < 1e-7);
            }
        }
    }
}

#[test]
fn solid_generated_core_and_sink_shell_keep_physical_contacts_films_and_sources() {
    for cylinder in [true, false] {
        let m = if cylinder { 2.0 } else { 3.0 };
        let wall = LayeredWall::new(
            vec![
                affine_layer("core", 1.0),
                curve("shell", 1.0, vec![(200.0, 2.0), (800.0, 2.0)]),
            ],
            vec![0.4],
            [0.0, 0.2],
        )
        .unwrap();
        let model = solid(cylinder, wall, 3.0);
        let answer = with_cx(|cx| {
            model.solve_temperature_dependent_with_generation(cx, 300.0, &[20.0, -2.0], controls())
        })
        .unwrap();
        let reference_area = physical_area(cylinder, 1.0, 3.0);
        let outer_area = physical_area(cylinder, 2.0, 3.0);
        let core_volume = physical_volume(cylinder, 1.0, 3.0);
        let shell_volume = physical_volume(cylinder, 2.0, 3.0) - core_volume;
        let core_heat = 20.0 * core_volume;
        let shell_heat = -2.0 * shell_volume;
        let total_heat = core_heat + shell_heat;
        let inlet = 20.0 / m;
        let (d, b) = if cylinder {
            (2.0_f64.ln(), 0.75 - 0.5 * 2.0_f64.ln())
        } else {
            (0.5, 1.0 / 3.0)
        };
        let outer_face = 300.0 + 0.2 * total_heat / outer_area;
        let shell_inner = outer_face + (inlet * d - 2.0 * b) / 2.0;
        let core_outer = shell_inner + 0.4 * core_heat / reference_area;
        let core_center = (core_outer * core_outer + 600.0 * 20.0 / (2.0 * m)).sqrt();
        let core = &answer.layer_temperatures[0];
        let shell = &answer.layer_temperatures[1];
        close(core.face_temperatures_k[0], core_center);
        close(core.face_temperatures_k[1], core_outer);
        close(shell.face_temperatures_k[0], shell_inner);
        close(shell.face_temperatures_k[1], outer_face);
        close(core.generated_heat_rate_w, core_heat);
        close(shell.generated_heat_rate_w, shell_heat);
        close(answer.generated_heat_rate_w, total_heat);
        close(answer.reservoir_heat_rate_w, total_heat);
        close(core.radial_heat_rates_w[1], shell.radial_heat_rates_w[0]);
        close(shell.radial_heat_rates_w[1], total_heat);
        close(core.face_areas_m2[1], reference_area);
        close(shell.face_areas_m2[0], reference_area);
        close(shell.face_areas_m2[1], outer_area);
        close(
            core_outer - shell_inner,
            0.4 * core.radial_heat_fluxes_w_per_m2[1],
        );
        close(
            outer_face - 300.0,
            0.2 * shell.radial_heat_fluxes_w_per_m2[1],
        );
        let r = 1.37;
        let (local_d, local_b) = shell_moments(cylinder, r);
        close(
            with_cx(|cx| answer.temperature_at(cx, 1, r - 1.0)).unwrap(),
            shell_inner - (inlet * local_d - 2.0 * local_b) / 2.0,
        );
        let local_heat = core_heat - 2.0 * (physical_volume(cylinder, r, 3.0) - core_volume);
        close(
            with_cx(|cx| answer.heat_flux_at(cx, 1, r - 1.0)).unwrap(),
            local_heat / physical_area(cylinder, r, 3.0),
        );
        close(answer.energy_residual_w, 0.0);
    }
}

#[test]
fn solid_shell_interior_cold_spot_is_found_and_original_source_support_is_enforced() {
    for cylinder in [true, false] {
        let m = if cylinder { 2.0 } else { 3.0 };
        let build = |lower| {
            LayeredWall::new(
                vec![
                    curve("core", 1.0, vec![(lower, 2.0), (350.0, 2.0)]),
                    curve("shell", 1.0, vec![(lower, 2.0), (350.0, 2.0)]),
                ],
                vec![0.0],
                [0.0; 2],
            )
            .unwrap()
        };
        let wide = NonlinearGenerationOptions {
            temperature_interval_k: [250.0, 350.0],
            ..controls()
        };
        let model = solid(cylinder, build(250.0), 3.0);
        let answer = with_cx(|cx| {
            model.solve_temperature_dependent_with_generation(cx, 300.0, &[20.0, -20.0], wide)
        })
        .unwrap();
        let inlet = 20.0 / m;
        let (d, b) = shell_moments(cylinder, 2.0);
        let shell_inner = 300.0 + (inlet * d - 20.0 * b) / 2.0;
        let minimum_radius = 2.0_f64.powf(1.0 / m);
        let (d_min, b_min) = shell_moments(cylinder, minimum_radius);
        let minimum = shell_inner - (inlet * d_min - 20.0 * b_min) / 2.0;
        let row = &answer.layer_temperatures[1];
        close(row.minimum_radius_m, minimum_radius);
        close(row.minimum_temperature_k, minimum);
        close(
            with_cx(|cx| answer.temperature_at(cx, 1, minimum_radius - 1.0)).unwrap(),
            minimum,
        );
        close(
            with_cx(|cx| answer.heat_flux_at(cx, 1, minimum_radius - 1.0)).unwrap(),
            0.0,
        );
        // Both face temperatures and the centre remain inside the narrowed
        // original source chart; only the true interior cold spot leaves it.
        let lower = 0.5 * (minimum + shell_inner.min(300.0));
        assert!(minimum < lower);
        for layer in &answer.layer_temperatures {
            assert!(layer.face_temperatures_k.iter().all(|&t| t > lower));
        }
        let narrowed = solid(cylinder, build(lower), 3.0);
        let controls = NonlinearGenerationOptions {
            temperature_interval_k: [lower, 350.0],
            ..controls()
        };
        assert!(
            with_cx(|cx| narrowed.solve_temperature_dependent_with_generation(
                cx,
                300.0,
                &[20.0, -20.0],
                controls,
            ))
            .is_err()
        );
    }
}

#[test]
fn solid_homogeneous_partition_and_cylinder_length_preserve_profiles_and_heat() {
    for cylinder in [true, false] {
        let whole_wall =
            LayeredWall::new(vec![affine_layer("whole", 2.0)], vec![], [0.0, 0.2]).unwrap();
        let split_wall = LayeredWall::new(
            vec![affine_layer("core", 0.75), affine_layer("shell", 1.25)],
            vec![0.0],
            [0.0, 0.2],
        )
        .unwrap();
        let whole = solid(cylinder, whole_wall.clone(), 3.0);
        let split = solid(cylinder, split_wall, 3.0);
        for source in [0.0, 30.0] {
            let a = with_cx(|cx| {
                whole.solve_temperature_dependent_with_generation(cx, 400.0, &[source], controls())
            })
            .unwrap();
            let b = with_cx(|cx| {
                split.solve_temperature_dependent_with_generation(
                    cx,
                    400.0,
                    &[source, source],
                    controls(),
                )
            })
            .unwrap();
            close(a.reservoir_heat_rate_w, b.reservoir_heat_rate_w);
            close(a.generated_heat_rate_w, b.generated_heat_rate_w);
            for r in [0.0, 0.2, 0.74, 0.75, 1.25, 2.0] {
                let (layer, offset) = if r <= 0.75 { (0, r) } else { (1, r - 0.75) };
                let temperature = with_cx(|cx| a.temperature_at(cx, 0, r)).unwrap();
                close(
                    temperature,
                    with_cx(|cx| b.temperature_at(cx, layer, offset)).unwrap(),
                );
                close(
                    with_cx(|cx| a.heat_flux_at(cx, 0, r)).unwrap(),
                    with_cx(|cx| b.heat_flux_at(cx, layer, offset)).unwrap(),
                );
                if source == 0.0 {
                    close(temperature, 400.0);
                }
            }
            if source == 0.0 {
                close(a.reservoir_heat_rate_w, 0.0);
                close(b.reservoir_heat_rate_w, 0.0);
            }
        }
        if cylinder {
            let longer = solid(true, whole_wall, 12.0);
            let a = with_cx(|cx| {
                whole.solve_temperature_dependent_with_generation(cx, 400.0, &[30.0], controls())
            })
            .unwrap();
            let b = with_cx(|cx| {
                longer.solve_temperature_dependent_with_generation(cx, 400.0, &[30.0], controls())
            })
            .unwrap();
            close(b.reservoir_heat_rate_w, 4.0 * a.reservoir_heat_rate_w);
            for r in [0.0, 0.75, 2.0] {
                close(
                    with_cx(|cx| a.temperature_at(cx, 0, r)).unwrap(),
                    with_cx(|cx| b.temperature_at(cx, 0, r)).unwrap(),
                );
                close(
                    with_cx(|cx| a.heat_flux_at(cx, 0, r)).unwrap(),
                    with_cx(|cx| b.heat_flux_at(cx, 0, r)).unwrap(),
                );
            }
        }
    }
}

#[test]
fn solid_radial_invalid_symmetry_sources_cancellation_and_queries_refuse() {
    let wall = LayeredWall::new(vec![affine_layer("core", 1.0)], vec![], [0.0; 2]).unwrap();
    let bad_film = LayeredWall::new(wall.layers().to_vec(), vec![], [0.1, 0.2]).unwrap();
    assert!(SolidRadialLayeredBody::spherical(bad_film.clone()).is_err());
    assert!(SolidRadialLayeredBody::cylindrical(bad_film, 3.0).is_err());
    assert!(SolidRadialLayeredBody::cylindrical(wall.clone(), 0.0).is_err());
    let body = SolidRadialLayeredBody::spherical(wall).unwrap();
    let before = body.clone();
    for sources in [vec![], vec![20.0, 0.0], vec![f64::NAN]] {
        assert!(
            with_cx(|cx| body.solve_temperature_dependent_with_generation(
                cx,
                300.0,
                &sources,
                controls(),
            ))
            .is_err()
        );
    }
    assert!(matches!(
        with_cancelled_cx(|cx| body.solve_temperature_dependent_with_generation(
            cx,
            300.0,
            &[20.0],
            controls(),
        )),
        Err(ConductionError::Cancelled { .. })
    ));
    let answer = with_cx(|cx| {
        body.solve_temperature_dependent_with_generation(cx, 300.0, &[20.0], controls())
    })
    .unwrap();
    for (layer, offset) in [
        (1, 0.0),
        (0, -f64::EPSILON),
        (0, 1.0_f64.next_up()),
        (0, f64::NAN),
    ] {
        assert!(with_cx(|cx| answer.temperature_at(cx, layer, offset)).is_err());
        assert!(with_cx(|cx| answer.heat_flux_at(cx, layer, offset)).is_err());
    }
    assert!(matches!(
        with_cancelled_cx(|cx| answer.temperature_at(cx, 0, 0.5)),
        Err(ConductionError::Cancelled { .. })
    ));
    assert_eq!(body, before);
}
