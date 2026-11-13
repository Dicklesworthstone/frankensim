//! G1/G3 curved composites: original sources, partitions and physical storage.
use super::*;
use fs_couple::moisture::radial::{
    MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry,
    layered::{MoistureLayeredRadial, MoistureRadialLayer},
};

#[path = "radial/adaptive.rs"]
mod adaptive;
#[path = "radial/atmosphere.rs"]
mod atmosphere;
#[path = "radial/sdirk.rs"]
mod sdirk;

fn radial_layer(
    card: MaterialCard,
    symmetry: RadialMoistureSymmetry,
    radii: Vec<f64>,
    activity: Vec<f64>,
    cx: &Cx<'_>,
) -> MoistureRadialLayer {
    let sorption = isotherm(&card);
    let initial = activity
        .iter()
        .map(|&a| sorption.moisture(a).unwrap())
        .collect();
    let radial = MoistureRadial::new(
        card,
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        MoistureRadialGeometry {
            symmetry,
            faces_m: radii,
            maximum_moisture: 0.6,
            maximum_cells: 128,
        },
        initial,
        cx,
    )
    .unwrap();
    MoistureRadialLayer::new(radial, sorption, cx).unwrap()
}

fn linear_card(rho: f64, slope: f64, diffusion: f64) -> MaterialCard {
    material(rho, vec![(0.0, 0.0), (1.0, slope)], [diffusion; 2])
}

fn measures(symmetry: RadialMoistureSymmetry, a: f64, b: f64) -> (f64, f64) {
    match symmetry {
        RadialMoistureSymmetry::Sphere => (
            4.0 * core::f64::consts::PI * (b.powi(3) - a.powi(3)) / 3.0,
            (1.0 / a - 1.0 / b) / (4.0 * core::f64::consts::PI),
        ),
        RadialMoistureSymmetry::Cylinder { length_m } => (
            core::f64::consts::PI * length_m * (b * b - a * a),
            (b / a).ln() / (2.0 * core::f64::consts::PI * length_m),
        ),
    }
}

fn area(symmetry: RadialMoistureSymmetry, radius: f64) -> f64 {
    match symmetry {
        RadialMoistureSymmetry::Sphere => 4.0 * core::f64::consts::PI * radius * radius,
        RadialMoistureSymmetry::Cylinder { length_m } => {
            2.0 * core::f64::consts::PI * length_m * radius
        }
    }
}

fn contents(state: &MoistureLayeredRadial) -> Vec<f64> {
    state
        .layers()
        .iter()
        .flat_map(|l| l.radial().moisture().iter().copied())
        .collect()
}

#[test]
fn layered_radial_reduces_to_native_nonlinear_diffusion_and_material_subdivision() {
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
            let original = full.radial().clone();
            let homogeneous = MoistureLayeredRadial::new(vec![full], 3, cx).unwrap();
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
            let dt = 0.01;
            let boundary = [
                MoistureLayeredFace::Sealed,
                MoistureLayeredFace::PrescribedRelativeHumidity(0.8),
            ];
            let (native, native_report) = original
                .advance(
                    dt,
                    &[MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)],
                    config(),
                    nonlinear(),
                    cx,
                )
                .unwrap();
            let (one, one_report) = homogeneous
                .advance(dt, &boundary, config(), nonlinear(), cx)
                .unwrap();
            let (two, two_report) = split
                .advance(dt, &boundary, config(), nonlinear(), cx)
                .unwrap();
            for ((a, b), &expected) in contents(&one)
                .into_iter()
                .zip(contents(&two))
                .zip(native.moisture())
            {
                near(a, expected, 5e-12);
                near(b, expected, 5e-12);
            }
            near(
                one_report.boundary_water_kg[1],
                native_report.boundary_water_kg[1],
                3e-12,
            );
            near(
                two_report.boundary_water_kg[1],
                native_report.boundary_water_kg[1],
                3e-12,
            );
            near(two.water_mass_kg(), one.water_mass_kg(), 3e-12);
            let transfer = two_report.interfaces[0].water_left_to_right_kg;
            near(
                two.layers()[0].radial().water_mass_kg()
                    - split.layers()[0].radial().water_mass_kg(),
                -transfer,
                3e-12,
            );
            near(
                two.layers()[1].radial().water_mass_kg()
                    - split.layers()[1].radial().water_mass_kg(),
                transfer + two_report.boundary_water_kg[1],
                3e-12,
            );
            assert_eq!(
                two.layers()[0].radial().identity(),
                split.layers()[0].radial().identity()
            );
            assert_eq!(two.time_s(), dt);
            assert_eq!(split.time_s(), 0.0);
            assert_eq!(
                one.layers()[0].radial().dry_masses_kg(),
                original.dry_masses_kg()
            );
        }
    });
}

#[test]
fn layered_radial_curved_steady_flow_keeps_partition_and_each_actual_transfer_area() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
        ] {
            for finite_join in [false, true] {
                let radii = [vec![0.02, 0.027, 0.04], vec![0.04, 0.052, 0.075, 0.09]];
                let rho = [450.0, 800.0];
                let slope = [0.2, 0.5];
                let diffusion = [0.001, 0.0004];
                let conductivity: [f64; 2] =
                    std::array::from_fn(|i| rho[i] * slope[i] * diffusion[i]);
                let exterior_h = [0.04, 0.02];
                let interface_h = 0.03;
                let film_r = [
                    1.0 / (area(symmetry, 0.02) * exterior_h[0]),
                    1.0 / (area(symmetry, 0.09) * exterior_h[1]),
                ];
                let join_r = if finite_join {
                    1.0 / (area(symmetry, 0.04) * interface_h)
                } else {
                    0.0
                };
                let bulk_r = [
                    measures(symmetry, 0.02, 0.04).1 / conductivity[0],
                    measures(symmetry, 0.04, 0.09).1 / conductivity[1],
                ];
                let rate = (0.85 - 0.15) / (film_r[0] + bulk_r[0] + join_r + bulk_r[1] + film_r[1]);
                let inner_surface = 0.85 - rate * film_r[0];
                let interface = [
                    inner_surface - rate * bulk_r[0],
                    inner_surface - rate * (bulk_r[0] + join_r),
                ];
                let initial: [Vec<f64>; 2] = std::array::from_fn(|k| {
                    radii[k]
                        .windows(2)
                        .map(|r| {
                            let center = 0.5 * (r[0] + r[1]);
                            let start = if k == 0 { inner_surface } else { interface[1] };
                            start
                                - rate * measures(symmetry, radii[k][0], center).1 / conductivity[k]
                        })
                        .collect()
                });
                let stack = MoistureLayeredRadial::new(
                    (0..2)
                        .map(|k| {
                            radial_layer(
                                linear_card(rho[k], slope[k], diffusion[k]),
                                symmetry,
                                radii[k].clone(),
                                initial[k].clone(),
                                cx,
                            )
                        })
                        .collect(),
                    5,
                    cx,
                )
                .unwrap()
                .with_transfer(
                    exterior_h.map(MoistureActivityTransfer::ConductanceKgM2S),
                    vec![if finite_join {
                        MoistureActivityTransfer::ConductanceKgM2S(interface_h)
                    } else {
                        MoistureActivityTransfer::IdealContact
                    }],
                    cx,
                )
                .unwrap();
                let dt = 0.05;
                let (next, report) = stack
                    .advance(
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
                for (&a, &old) in next.activity().iter().zip(stack.activity()) {
                    near(a, old, 3e-12);
                }
                near(report.boundary_water_kg[0], dt * rate, 2e-12);
                near(report.boundary_water_kg[1], -dt * rate, 2e-12);
                near(
                    report.interfaces[0].water_left_to_right_kg,
                    dt * rate,
                    2e-12,
                );
                near(report.surface_activities[0].unwrap(), inner_surface, 3e-12);
                near(
                    report.surface_activities[1].unwrap(),
                    0.15 + rate * film_r[1],
                    3e-12,
                );
                near(report.interfaces[0].relative_activity, interface[0], 3e-12);
                near(
                    report.interfaces[0].right_relative_activity,
                    interface[1],
                    3e-12,
                );
                near(
                    report.interfaces[0].moisture[0],
                    slope[0] * interface[0],
                    3e-12,
                );
                near(
                    report.interfaces[0].moisture[1],
                    slope[1] * interface[1],
                    3e-12,
                );
                assert_ne!(
                    report.interfaces[0].moisture[0],
                    report.interfaces[0].moisture[1]
                );
                assert_eq!(
                    report.interfaces[0].relative_activity
                        > report.interfaces[0].right_relative_activity,
                    finite_join
                );
                near(next.water_mass_kg(), stack.water_mass_kg(), 3e-12);
                assert_eq!(report.diffusivity_sources.len(), 2);
                assert_eq!(report.sorption_sources.len(), 2);
            }
        }
    });
}

#[test]
fn layered_radial_unequal_density_partition_has_exact_closed_two_cell_transient() {
    with_cx(|cx| {
        let rho = [300.0, 900.0];
        let slopes = [0.2, 0.4];
        let diffusion = [1e-4, 2e-5];
        let initial = [0.8, 0.1];
        let h = 0.002;
        let dt = 0.3;
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.2 },
        ] {
            let radii = [[0.0, 0.03], [0.03, 0.08]];
            let mass: [f64; 2] =
                std::array::from_fn(|i| rho[i] * measures(symmetry, radii[i][0], radii[i][1]).0);
            let capacity: [f64; 2] = std::array::from_fn(|i| mass[i] * slopes[i]);
            let gl = rho[0] * diffusion[0] * slopes[0] / measures(symmetry, 0.015, 0.03).1;
            let gr = rho[1] * diffusion[1] * slopes[1] / measures(symmetry, 0.03, 0.055).1;
            let gh = h * area(symmetry, 0.03);
            let g = 1.0 / (1.0 / gl + 1.0 / gh + 1.0 / gr);
            let mean =
                (capacity[0] * initial[0] + capacity[1] * initial[1]) / (capacity[0] + capacity[1]);
            let difference = (initial[0] - initial[1])
                / (1.0 + dt * g * (1.0 / capacity[0] + 1.0 / capacity[1]));
            let expected = [
                mean + capacity[1] / (capacity[0] + capacity[1]) * difference,
                mean - capacity[0] / (capacity[0] + capacity[1]) * difference,
            ];
            let rate = g * difference;
            let stack = MoistureLayeredRadial::new(
                (0..2)
                    .map(|i| {
                        radial_layer(
                            linear_card(rho[i], slopes[i], diffusion[i]),
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
                [MoistureActivityTransfer::IdealContact; 2],
                vec![MoistureActivityTransfer::ConductanceKgM2S(h)],
                cx,
            )
            .unwrap();
            let (next, report) = stack
                .advance(dt, &sealed(), config(), nonlinear(), cx)
                .unwrap();
            for i in 0..2 {
                near(next.activity()[i], expected[i], 4e-12);
                near(
                    next.layers()[i].radial().moisture()[0],
                    slopes[i] * expected[i],
                    4e-12,
                );
                near(next.layers()[i].radial().dry_mass_kg(), mass[i], 2e-14);
            }
            near(
                report.interfaces[0].water_left_to_right_kg,
                dt * rate,
                2e-12,
            );
            near(
                report.interfaces[0].relative_activity,
                expected[0] - rate / gl,
                3e-12,
            );
            near(
                report.interfaces[0].right_relative_activity,
                expected[1] + rate / gr,
                3e-12,
            );
            near(next.water_mass_kg(), stack.water_mass_kg(), 3e-12);
            assert_eq!(report.boundary_water_kg, [0.0; 2]);
            near(report.water_balance_residual_kg, 0.0, 2e-12);
            assert_eq!(stack.time_s(), 0.0);
        }
    });
}

#[test]
fn layered_radial_nonlinear_potential_retains_curved_steady_partition_and_wets_dry_coating() {
    with_cx(|cx| {
        let symmetry = RadialMoistureSymmetry::Sphere;
        let radii = [vec![0.02, 0.03, 0.045], vec![0.045, 0.065, 0.09]];
        let rho = [500.0, 750.0];
        let slope: [f64; 2] = [0.25, 0.5];
        let alpha = [0.01, 0.025];
        let coefficient: [f64; 2] =
            std::array::from_fn(|i| 0.5 * rho[i] * alpha[i] * slope[i].powi(2));
        // D_i(M)=alpha_i*M and M_i=s_i*a give K_i rho_i=c_i*a².
        // The exact steady squared activity is affine in spherical resistance.
        let resistance = [
            measures(symmetry, 0.02, 0.045).1 / coefficient[0],
            measures(symmetry, 0.045, 0.09).1 / coefficient[1],
        ];
        let rate = 0.8_f64.powi(2) / (resistance[0] + resistance[1]);
        let interface_squared = 0.8_f64.powi(2) - rate * resistance[0];
        let source = |i| {
            material(
                rho[i],
                vec![(0.0, 0.0), (1.0, slope[i])],
                [0.0, 0.6 * alpha[i]],
            )
        };
        let initial: [Vec<f64>; 2] = std::array::from_fn(|i| {
            radii[i]
                .windows(2)
                .map(|r| {
                    let center = 0.5 * (r[0] + r[1]);
                    let start = if i == 0 {
                        0.8_f64.powi(2)
                    } else {
                        interface_squared
                    };
                    (start - rate * measures(symmetry, radii[i][0], center).1 / coefficient[i])
                        .sqrt()
                })
                .collect()
        });
        let stack = MoistureLayeredRadial::new(
            (0..2)
                .map(|i| {
                    radial_layer(
                        source(i),
                        symmetry,
                        radii[i].clone(),
                        initial[i].clone(),
                        cx,
                    )
                })
                .collect(),
            4,
            cx,
        )
        .unwrap();
        let dt = 0.1;
        let (next, report) = stack
            .advance(
                dt,
                &[
                    MoistureLayeredFace::PrescribedRelativeHumidity(0.8),
                    MoistureLayeredFace::PrescribedRelativeHumidity(0.0),
                ],
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        for (&a, &old) in next.activity().iter().zip(stack.activity()) {
            near(a, old, 3e-12);
        }
        near(report.boundary_water_kg[0], dt * rate, 3e-12);
        near(report.boundary_water_kg[1], -dt * rate, 3e-12);
        near(
            report.interfaces[0].relative_activity,
            interface_squared.sqrt(),
            3e-12,
        );
        near(
            report.interfaces[0].water_left_to_right_kg,
            dt * rate,
            3e-12,
        );

        let dry = MoistureLayeredRadial::new(
            vec![
                radial_layer(source(0), symmetry, vec![0.0, 0.03], vec![0.0], cx),
                radial_layer(source(1), symmetry, vec![0.03, 0.05], vec![0.0], cx),
            ],
            2,
            cx,
        )
        .unwrap();
        let (wet, uptake) = dry
            .advance(
                0.01,
                &[
                    MoistureLayeredFace::Sealed,
                    MoistureLayeredFace::PrescribedRelativeHumidity(0.8),
                ],
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        assert!(wet.activity()[0] > 0.0 && wet.activity()[1] > wet.activity()[0]);
        assert!(uptake.interfaces[0].water_left_to_right_kg < 0.0);
        near(wet.water_mass_kg(), uptake.boundary_water_kg[1], 2e-12);
        near(
            wet.layers()[0].radial().water_mass_kg(),
            -uptake.interfaces[0].water_left_to_right_kg,
            2e-12,
        );
        near(
            wet.layers()[1].radial().water_mass_kg(),
            uptake.boundary_water_kg[1] + uptake.interfaces[0].water_left_to_right_kg,
            2e-12,
        );
        assert_eq!(contents(&dry), vec![0.0, 0.0]);
    });
}

#[test]
fn layered_radial_preserves_partition_equilibrium_and_refuses_bad_joins_and_centers() {
    with_cx(|cx| {
        let symmetry = RadialMoistureSymmetry::Sphere;
        let core = radial_layer(
            linear_card(200.0, 0.25, 0.001),
            symmetry,
            vec![0.0, 0.03],
            vec![0.5],
            cx,
        );
        let coating = radial_layer(
            linear_card(800.0, 0.5, 0.004),
            symmetry,
            vec![0.03, 0.05],
            vec![0.5],
            cx,
        );
        let stack = MoistureLayeredRadial::new(vec![core.clone(), coating.clone()], 2, cx).unwrap();
        let (same, report) = stack
            .advance(1e8, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        assert_eq!(contents(&same), vec![0.125, 0.25]);
        assert_eq!(same.activity(), &[0.5; 2]);
        assert_eq!(report.interfaces[0].water_left_to_right_kg, 0.0);
        assert_eq!(report.interfaces[0].moisture, [0.125, 0.25]);
        assert!(MoistureLayeredRadial::new(vec![core.clone(), coating.clone()], 1, cx).is_err());
        let gap = radial_layer(
            linear_card(800.0, 0.5, 0.004),
            symmetry,
            vec![0.031, 0.05],
            vec![0.5],
            cx,
        );
        assert!(matches!(
            MoistureLayeredRadial::new(vec![core.clone(), gap], 2, cx),
            Err(MoistureError::Invalid(_))
        ));
        let different_context = isotherm(&linear_card(300.0, 0.25, 0.001));
        assert!(matches!(
            MoistureRadialLayer::new(core.radial().clone(), different_context, cx),
            Err(MoistureError::Invalid(_))
        ));
        assert!(matches!(
            stack.advance(
                0.1,
                &[
                    MoistureLayeredFace::PrescribedRelativeHumidity(0.8),
                    MoistureLayeredFace::Sealed
                ],
                config(),
                nonlinear(),
                cx
            ),
            Err(MoistureError::Invalid(_))
        ));
        assert!(
            stack
                .clone()
                .with_transfer(
                    [
                        MoistureActivityTransfer::ConductanceKgM2S(0.01),
                        MoistureActivityTransfer::IdealContact
                    ],
                    vec![MoistureActivityTransfer::IdealContact],
                    cx
                )
                .is_err()
        );
        let mut limited = config();
        limited.maximum_content_increment = 1e-12;
        assert!(matches!(
            stack.advance(
                0.1,
                &[
                    MoistureLayeredFace::Sealed,
                    MoistureLayeredFace::PrescribedRelativeHumidity(0.8)
                ],
                limited,
                nonlinear(),
                cx
            ),
            Err(MoistureError::ContentIncrement { .. })
        ));
        assert!(matches!(
            stack.advance(
                0.1,
                &[
                    MoistureLayeredFace::Sealed,
                    MoistureLayeredFace::PrescribedRelativeHumidity(1.1)
                ],
                config(),
                nonlinear(),
                cx
            ),
            Err(MoistureError::Invalid(_))
        ));
        assert_eq!(contents(&stack), vec![0.125, 0.25]);
        assert_eq!(stack.time_s(), 0.0);

        let isolated = MoistureLayeredRadial::new(
            vec![
                radial_layer(
                    linear_card(200.0, 0.25, 0.001),
                    symmetry,
                    vec![0.0, 0.03],
                    vec![0.8],
                    cx,
                ),
                radial_layer(
                    linear_card(800.0, 0.5, 0.004),
                    symmetry,
                    vec![0.03, 0.05],
                    vec![0.1],
                    cx,
                ),
            ],
            2,
            cx,
        )
        .unwrap()
        .with_transfer(
            [MoistureActivityTransfer::IdealContact; 2],
            vec![MoistureActivityTransfer::Sealed],
            cx,
        )
        .unwrap();
        let (closed, account) = isolated
            .advance(1.0, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        for (&a, &old) in closed.activity().iter().zip(isolated.activity()) {
            near(a, old, 1e-15);
        }
        assert_eq!(account.interfaces[0].water_left_to_right_kg, 0.0);
        near(closed.water_mass_kg(), isolated.water_mass_kg(), 1e-14);
    });
}
