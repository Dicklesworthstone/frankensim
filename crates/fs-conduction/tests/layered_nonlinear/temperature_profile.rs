//! Original nonlinear source derivatives for refractive thermal snapshots.
use super::*;
use fs_conduction::layered::{NonlinearGenerationOptions, NonlinearTemperatureProfileError};

fn generated_options(interval: [f64; 2]) -> NonlinearGenerationOptions {
    NonlinearGenerationOptions {
        temperature_interval_k: interval,
        solver: options(),
    }
}

fn affine_wall() -> LayeredWall {
    LayeredWall::new(
        vec![curve(
            "k=T/300",
            2.0,
            vec![(280.0, 280.0 / 300.0), (800.0, 800.0 / 300.0)],
        )],
        vec![],
        [0.0; 2],
    )
    .unwrap()
}

#[test]
fn nonlinear_profile_derivatives_follow_the_integrated_conductivity_law() {
    with_cx(|cx| {
        for (boundary, source) in [(300.0, 600.0), (600.0, -600.0)] {
            let solution = affine_wall()
                .solve_temperature_dependent_with_generation(
                    cx,
                    [boundary; 2],
                    2.0,
                    &[source],
                    generated_options([280.0, 800.0]),
                )
                .unwrap();
            let profile = solution.temperature_profile(cx, 0).unwrap();
            assert!(core::ptr::eq(profile.solution(), &solution));
            assert_eq!(profile.layer_index(), 0);
            assert_eq!(profile.thickness_m(), 2.0);
            for x in [0.0, 0.25, 0.5, 1.0, 1.75, 2.0] {
                // Independent Kirchhoff solution for k=T/300 and L=2:
                // T²=Tb²+300s(2x-x²), T'=300s(1-x)/T.
                let temperature = (boundary * boundary + 300.0 * source * (2.0 * x - x * x)).sqrt();
                let numerator = 300.0 * source * (1.0 - x);
                let gradient = numerator / temperature;
                let curvature =
                    -300.0 * source / temperature - numerator * numerator / temperature.powi(3);
                let jet = profile.temperature_jet_at(cx, x).unwrap();
                close(jet.sample.temperature_k, temperature);
                close(jet.sample.heat_flux_w_per_m2, source * (x - 1.0));
                close(jet.sample.gradient_k_per_m, gradient);
                close(jet.curvature_k_per_m2, curvature);
                assert_eq!(profile.temperature_gradient_at(cx, x).unwrap(), jet.sample);
                assert_eq!(solution.temperature_jet_at(cx, 0, x).unwrap(), jet);
            }
            for (x, face) in [(0.0, 0), (2.0, 1)] {
                let sample = profile.temperature_gradient_at(cx, x).unwrap();
                assert_eq!(
                    sample.temperature_k.to_bits(),
                    solution.layer_temperatures[0].face_temperatures_k[face].to_bits()
                );
                assert_eq!(
                    sample.heat_flux_w_per_m2.to_bits(),
                    solution.layer_temperatures[0].face_fluxes_w_per_m2[face].to_bits()
                );
            }
        }
        let solution = affine_wall()
            .solve_temperature_dependent_with_generation(
                cx,
                [600.0, 300.0],
                1.0,
                &[0.0],
                generated_options([280.0, 800.0]),
            )
            .unwrap();
        for x in [0.0, 0.5, 1.0, 2.0] {
            // With no generation, q=225 and T²=360000-135000x.
            // Curvature remains nonzero solely because conductivity depends on T.
            let temperature = (360_000.0_f64 - 135_000.0 * x).sqrt();
            let jet = solution.temperature_jet_at(cx, 0, x).unwrap();
            close(jet.sample.heat_flux_w_per_m2, 225.0);
            close(jet.sample.temperature_k, temperature);
            close(jet.sample.gradient_k_per_m, -67_500.0 / temperature);
            close(
                jet.curvature_k_per_m2,
                -67_500.0_f64.powi(2) / temperature.powi(3),
            );
            assert!(jet.curvature_k_per_m2 < 0.0);
        }
    });
}

#[test]
fn nonlinear_profile_distinguishes_source_curvature_jumps_from_zero_flux_knots() {
    let wall = LayeredWall::new(
        vec![curve(
            "slope changes",
            2.0,
            vec![(290.0, 1.0), (302.0, 1.0), (304.0, 3.0), (310.0, 3.0)],
        )],
        vec![],
        [0.0; 2],
    )
    .unwrap();
    with_cx(|cx| {
        for (reservoirs, knot_x, face_x) in [
            ([304.0, 300.0], 4.0 / 3.0, 0.0),
            ([300.0, 304.0], 2.0 / 3.0, 2.0),
        ] {
            let solution = wall
                .solve_temperature_dependent_with_generation(
                    cx,
                    reservoirs,
                    1.0,
                    &[0.0],
                    generated_options([290.0, 304.0]),
                )
                .unwrap();
            let first = solution.temperature_gradient_at(cx, 0, knot_x).unwrap();
            assert_eq!(first.temperature_k, 302.0);
            close(first.heat_flux_w_per_m2.abs(), 3.0);
            close(first.gradient_k_per_m.abs(), 3.0);
            assert_eq!(
                solution.temperature_jet_at(cx, 0, knot_x),
                Err(NonlinearTemperatureProfileError::UndefinedCurvature {
                    temperature_k: 302.0
                })
            );
            // The304 K physical face must use the rising source segment from
            // inside the layer, although the table's right-hand slope is zero.
            let face = solution.temperature_jet_at(cx, 0, face_x).unwrap();
            close(face.sample.temperature_k, 304.0);
            close(face.curvature_k_per_m2, -1.0 / 3.0);
        }
        let zero_flux = wall
            .solve_temperature_dependent_with_generation(
                cx,
                [300.0; 2],
                1.0,
                &[4.0],
                generated_options([290.0, 310.0]),
            )
            .unwrap();
        let peak = zero_flux.temperature_jet_at(cx, 0, 1.0).unwrap();
        assert_eq!(peak.sample.temperature_k, 302.0);
        assert_eq!(peak.sample.heat_flux_w_per_m2, 0.0);
        assert_eq!(peak.sample.gradient_k_per_m, 0.0);
        assert_eq!(peak.curvature_k_per_m2, -4.0);
    });
}

#[test]
fn nonlinear_profile_retains_physical_support_and_source_cancellation() {
    let solution = with_cx(|cx| {
        affine_wall().solve_temperature_dependent_with_generation(
            cx,
            [300.0; 2],
            1.0,
            &[600.0],
            generated_options([280.0, 800.0]),
        )
    })
    .unwrap();
    with_cx(|cx| {
        assert!(solution.temperature_profile(cx, 1).is_err());
        let profile = solution.temperature_profile(cx, 0).unwrap();
        for x in [-0.01, 2.01, f64::NAN, f64::INFINITY] {
            assert!(profile.temperature_gradient_at(cx, x).is_err());
            assert!(profile.temperature_jet_at(cx, x).is_err());
        }
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                profile.temperature_gradient_at(cancelled, 1.0),
                Err(ConductionError::Cancelled { .. })
            ));
            assert!(matches!(
                profile.temperature_jet_at(cancelled, 1.0),
                Err(NonlinearTemperatureProfileError::Conduction(
                    ConductionError::Cancelled { .. }
                ))
            ));
        });
    });
}
