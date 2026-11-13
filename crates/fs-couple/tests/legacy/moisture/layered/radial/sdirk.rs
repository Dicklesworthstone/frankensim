//! Time accuracy and complete physical accounts on actual curved layer paths.
use super::*;

const GAMMA: f64 = 1.0 - core::f64::consts::FRAC_1_SQRT_2;

#[test]
fn coated_sdirk_converges_with_unequal_shell_storage_and_finite_curved_contact() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.2 },
        ] {
            let rho = [300.0, 900.0];
            let slope = [0.2, 0.4];
            let diffusion = [0.001, 0.0002];
            let initial = [0.8, 0.1];
            let radius = [[0.0, 0.03], [0.03, 0.08]];
            let mass: [f64; 2] =
                std::array::from_fn(|i| rho[i] * measures(symmetry, radius[i][0], radius[i][1]).0);
            let capacity: [f64; 2] = std::array::from_fn(|i| mass[i] * slope[i]);
            let h = 0.03;
            let half = [
                rho[0] * diffusion[0] * slope[0] / measures(symmetry, 0.015, 0.03).1,
                rho[1] * diffusion[1] * slope[1] / measures(symmetry, 0.03, 0.055).1,
            ];
            let g = 1.0 / (1.0 / half[0] + 1.0 / (h * area(symmetry, 0.03)) + 1.0 / half[1]);
            let rate = g * (1.0 / capacity[0] + 1.0 / capacity[1]);
            let duration = 1.0 / rate;
            let mean =
                (capacity[0] * initial[0] + capacity[1] * initial[1]) / (capacity[0] + capacity[1]);
            let exact_difference = (initial[0] - initial[1]) * (-rate * duration).exp();
            let exact = [
                mean + capacity[1] / (capacity[0] + capacity[1]) * exact_difference,
                mean - capacity[0] / (capacity[0] + capacity[1]) * exact_difference,
            ];
            let original = MoistureLayeredRadial::new(
                (0..2)
                    .map(|i| {
                        radial_layer(
                            linear_card(rho[i], slope[i], diffusion[i]),
                            symmetry,
                            radius[i].to_vec(),
                            vec![initial[i]],
                            cx,
                        )
                    })
                    .collect(),
                2,
                cx,
            )
            .unwrap()
            .with_transfer(
                [MoistureActivityTransfer::IdealContact; 2],
                vec![MoistureActivityTransfer::ConductanceKgM2S(h)],
                cx,
            )
            .unwrap();
            let mut errors = Vec::new();
            for count in [4, 8, 16] {
                let dt = duration / f64::from(count);
                let mut state = original.clone();
                for _ in 0..count {
                    let old_water: [f64; 2] =
                        std::array::from_fn(|i| state.layers()[i].radial().water_mass_kg());
                    let (next, report) = state
                        .advance_sdirk2(dt, &sealed(), config(), nonlinear(), cx)
                        .unwrap();
                    let first: [f64; 2] = std::array::from_fn(|i| {
                        report.first_stage.transport_end[i].query().axes()["moisture_content"]
                            / slope[i]
                    });
                    let first_rate = g * (first[0] - first[1]);
                    let last_rate = g * (next.activity()[0] - next.activity()[1]);
                    let transferred = dt * ((1.0 - GAMMA) * first_rate + GAMMA * last_rate);
                    near(
                        report.endpoint.interfaces[0].water_left_to_right_kg,
                        transferred,
                        2e-12,
                    );
                    near(
                        report.first_stage.interfaces[0].water_left_to_right_kg,
                        GAMMA * dt * first_rate,
                        2e-12,
                    );
                    near(
                        next.layers()[0].radial().water_mass_kg() - old_water[0],
                        -transferred,
                        2e-12,
                    );
                    near(
                        next.layers()[1].radial().water_mass_kg() - old_water[1],
                        transferred,
                        2e-12,
                    );
                    near(
                        report.endpoint.interfaces[0].relative_activity,
                        next.activity()[0] - last_rate / half[0],
                        3e-12,
                    );
                    near(
                        report.endpoint.interfaces[0].right_relative_activity,
                        next.activity()[1] + last_rate / half[1],
                        3e-12,
                    );
                    near(next.water_mass_kg(), original.water_mass_kg(), 3e-12);
                    for i in 0..2 {
                        near(next.layers()[i].radial().dry_mass_kg(), mass[i], 2e-14);
                        assert_eq!(
                            next.layers()[i].radial().dry_masses_kg(),
                            original.layers()[i].radial().dry_masses_kg()
                        );
                    }
                    assert_eq!(report.endpoint.boundary_water_kg, [0.0; 2]);
                    state = next;
                }
                errors.push(
                    (state.activity()[0] - exact[0])
                        .abs()
                        .max((state.activity()[1] - exact[1]).abs()),
                );
            }
            for pair in errors.windows(2) {
                assert!(
                    pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.3,
                    "{symmetry:?}: {errors:?}"
                );
            }
            assert_eq!(original.time_s(), 0.0);
            assert!(
                original
                    .advance_sdirk2(
                        0.1,
                        &[
                            MoistureLayeredFace::PrescribedRelativeHumidity(0.5),
                            MoistureLayeredFace::Sealed
                        ],
                        config(),
                        nonlinear(),
                        cx
                    )
                    .is_err()
            );
            assert_eq!(original.time_s(), 0.0);
        }
    });
}

#[test]
fn coated_sdirk_preserves_steady_throughflow_and_reports_a_full_interval_once() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
        ] {
            let radii = [[0.02, 0.04], [0.04, 0.09]];
            let center = [0.03, 0.065];
            let rho = [450.0, 800.0];
            let slope = [0.2, 0.5];
            let diffusion = [0.001, 0.0004];
            let conductivity: [f64; 2] = std::array::from_fn(|i| rho[i] * slope[i] * diffusion[i]);
            let exterior_h = [0.04, 0.02];
            let interface_h = 0.03;
            let film = [
                1.0 / (area(symmetry, 0.02) * exterior_h[0]),
                1.0 / (area(symmetry, 0.09) * exterior_h[1]),
            ];
            let join = 1.0 / (area(symmetry, 0.04) * interface_h);
            let bulk: [f64; 2] = std::array::from_fn(|i| {
                measures(symmetry, radii[i][0], radii[i][1]).1 / conductivity[i]
            });
            let flux = (0.85 - 0.15) / (film[0] + bulk[0] + join + bulk[1] + film[1]);
            let left_surface = 0.85 - flux * film[0];
            let contact = [
                left_surface - flux * bulk[0],
                left_surface - flux * (bulk[0] + join),
            ];
            let initial = [
                left_surface - flux * measures(symmetry, 0.02, center[0]).1 / conductivity[0],
                contact[1] - flux * measures(symmetry, 0.04, center[1]).1 / conductivity[1],
            ];
            let original = MoistureLayeredRadial::new(
                (0..2)
                    .map(|i| {
                        radial_layer(
                            linear_card(rho[i], slope[i], diffusion[i]),
                            symmetry,
                            radii[i].to_vec(),
                            vec![initial[i]],
                            cx,
                        )
                    })
                    .collect(),
                2,
                cx,
            )
            .unwrap()
            .with_transfer(
                exterior_h.map(MoistureActivityTransfer::ConductanceKgM2S),
                vec![MoistureActivityTransfer::ConductanceKgM2S(interface_h)],
                cx,
            )
            .unwrap();
            let dt = 37.5;
            let (next, report) = original
                .advance_sdirk2(
                    dt,
                    &[
                        MoistureLayeredFace::PrescribedRelativeHumidity(0.85),
                        MoistureLayeredFace::PrescribedRelativeHumidity(0.15),
                    ],
                    config(),
                    nonlinear(),
                    cx,
                )
                .unwrap();
            for i in 0..2 {
                near(next.activity()[i], initial[i], 4e-12);
                near(
                    next.layers()[i].radial().water_mass_kg(),
                    original.layers()[i].radial().water_mass_kg(),
                    3e-12,
                );
            }
            near(report.endpoint.boundary_water_kg[0], dt * flux, 3e-12);
            near(report.endpoint.boundary_water_kg[1], -dt * flux, 3e-12);
            near(
                report.endpoint.interfaces[0].water_left_to_right_kg,
                dt * flux,
                3e-12,
            );
            near(
                report.first_stage.interfaces[0].water_left_to_right_kg,
                GAMMA * dt * flux,
                3e-12,
            );
            near(
                report.endpoint.interfaces[0].relative_activity,
                contact[0],
                3e-12,
            );
            near(
                report.endpoint.interfaces[0].right_relative_activity,
                contact[1],
                3e-12,
            );
            near(
                report.endpoint.interfaces[0].moisture[0],
                slope[0] * contact[0],
                3e-12,
            );
            near(
                report.endpoint.interfaces[0].moisture[1],
                slope[1] * contact[1],
                3e-12,
            );
            near(report.endpoint.water_balance_residual_kg, 0.0, 3e-12);
            assert_eq!(report.endpoint.interval_s, [0.0, dt]);
            assert_eq!(report.first_stage.interval_s, [0.0, GAMMA * dt]);
            assert_eq!(original.time_s(), 0.0);
        }
    });
}

#[test]
fn layered_sdirk_reduces_to_the_native_nonlinear_radial_owner_under_subdivision() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.2 },
        ] {
            let source = material(400.0, vec![(0.0, 0.0), (1.0, 0.25)], [0.0002, 0.0014]);
            let full = radial_layer(
                source.clone(),
                symmetry,
                vec![0.0, 0.005, 0.02, 0.05],
                vec![0.5, 0.1, 0.3],
                cx,
            );
            let native = full.radial().clone();
            let one = MoistureLayeredRadial::new(vec![full], 3, cx).unwrap();
            let split = MoistureLayeredRadial::new(
                vec![
                    radial_layer(
                        source.clone(),
                        symmetry,
                        vec![0.0, 0.005, 0.02],
                        vec![0.5, 0.1],
                        cx,
                    ),
                    radial_layer(source, symmetry, vec![0.02, 0.05], vec![0.3], cx),
                ],
                3,
                cx,
            )
            .unwrap();
            let dt = 0.001;
            let boundary = [
                MoistureLayeredFace::Sealed,
                MoistureLayeredFace::PrescribedRelativeHumidity(0.8),
            ];
            let (expected, native_report) = native
                .advance_sdirk2(
                    dt,
                    &[MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)],
                    config(),
                    nonlinear(),
                    cx,
                )
                .unwrap();
            let (one_next, one_report) = one
                .advance_sdirk2(dt, &boundary, config(), nonlinear(), cx)
                .unwrap();
            let (split_next, split_report) = split
                .advance_sdirk2(dt, &boundary, config(), nonlinear(), cx)
                .unwrap();
            for ((a, b), &expected) in contents(&one_next)
                .into_iter()
                .zip(contents(&split_next))
                .zip(expected.moisture())
            {
                near(a, expected, 5e-12);
                near(b, expected, 5e-12);
            }
            near(
                one_report.endpoint.boundary_water_kg[1],
                native_report.endpoint.boundary_water_kg[1],
                3e-12,
            );
            near(
                split_report.endpoint.boundary_water_kg[1],
                native_report.endpoint.boundary_water_kg[1],
                3e-12,
            );
            near(
                split_report.first_stage.boundary_water_kg[1],
                native_report.first_stage.boundary_water_kg[1],
                3e-12,
            );
            let transfer = split_report.endpoint.interfaces[0].water_left_to_right_kg;
            near(
                split_next.layers()[0].radial().water_mass_kg()
                    - split.layers()[0].radial().water_mass_kg(),
                -transfer,
                3e-12,
            );
            near(
                split_next.layers()[1].radial().water_mass_kg()
                    - split.layers()[1].radial().water_mass_kg(),
                transfer + split_report.endpoint.boundary_water_kg[1],
                3e-12,
            );
            assert_eq!(
                split_next.layers()[0].radial().dry_masses_kg(),
                split.layers()[0].radial().dry_masses_kg()
            );
            assert_eq!(split_next.time_s(), dt);
            assert_eq!(split.time_s(), 0.0);
        }
    });
}
