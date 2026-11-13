//! Authored radial diffusion laws, checked against independent continuum limits.
use super::*;
use fs_couple::moisture::radial::*;

#[path = "radial/adaptive.rs"]
mod adaptive;
#[path = "radial/atmosphere.rs"]
mod finite_atmosphere;
#[path = "radial/solution.rs"]
mod finite_solution;
#[path = "radial/sdirk.rs"]
mod sdirk;

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn specimen(
    symmetry: RadialMoistureSymmetry,
    radii: Vec<f64>,
    initial: Vec<f64>,
    diffusion: [f64; 2],
    cx: &Cx<'_>,
) -> MoistureRadial {
    MoistureRadial::new(
        card(diffusion, 0.3),
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        MoistureRadialGeometry {
            symmetry,
            faces_m: radii,
            maximum_moisture: 0.3,
            maximum_cells: 256,
        },
        initial,
        cx,
    )
    .unwrap()
}

#[test]
fn radial_nonuniform_shells_retain_exact_nonlinear_steady_flux() {
    with_cx(|cx| {
        let a: f64 = 0.02;
        let b: f64 = 0.1;
        let left: f64 = 0.05;
        let right: f64 = 0.2;
        let radii: Vec<_> = (0..=17)
            .map(|i| a + (b - a) * (i as f64 / 17.0).powi(2))
            .collect();
        for symmetry in [
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
            RadialMoistureSymmetry::Sphere,
        ] {
            let (factor, resistance) = match symmetry {
                RadialMoistureSymmetry::Cylinder { length_m } => {
                    (2.0 * core::f64::consts::PI * length_m, (b / a).ln())
                }
                RadialMoistureSymmetry::Sphere => (4.0 * core::f64::consts::PI, 1.0 / a - 1.0 / b),
            };
            // D(M)=0.1*M, hence K(M)=0.05*M^2. Its steady potential is
            // affine in log(r) for cylinders and in -1/r for spheres.
            let initial: Vec<_> = radii
                .windows(2)
                .map(|r| {
                    let center = 0.5 * (r[0] + r[1]);
                    let coordinate = match symmetry {
                        RadialMoistureSymmetry::Cylinder { .. } => (center / a).ln(),
                        RadialMoistureSymmetry::Sphere => 1.0 / a - 1.0 / center,
                    } / resistance;
                    (left * left + coordinate * (right * right - left * left)).sqrt()
                })
                .collect();
            let old = specimen(symmetry, radii.clone(), initial.clone(), [0.0, 0.03], cx);
            let boundary = [
                MoistureFace::PrescribedDryBasis(left),
                MoistureFace::PrescribedDryBasis(right),
            ];
            let dt = 0.2;
            let (next, report) = old
                .advance(dt, &boundary, config(), nonlinear(), cx)
                .unwrap();
            for (&actual, &expected) in next.moisture().iter().zip(&initial) {
                near(actual, expected, 2e-13);
            }
            let incoming = dt * 500.0 * factor * 0.05 * (right * right - left * left) / resistance;
            near(report.boundary_water_kg[0], -incoming, 1e-11);
            near(report.boundary_water_kg[1], incoming, 1e-11);
            near(next.water_mass_kg(), old.water_mass_kg(), 1e-12);
            assert!(report.kirchhoff_source.is_some());
        }
    });
}

#[test]
fn radial_sphere_dry_front_matches_nonlinear_one_cell_balance() {
    with_cx(|cx| {
        let radius: f64 = 0.1;
        let dt = 0.01;
        let equilibrium: f64 = 0.2;
        let old = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![0.0],
            [0.0, 0.03],
            cx,
        );
        let boundary = [
            MoistureFace::Sealed,
            MoistureFace::PrescribedDryBasis(equilibrium),
        ];
        let (next, report) = old
            .advance(dt, &boundary, config(), nonlinear(), cx)
            .unwrap();
        // For center R/2, resistance=(1/(R/2)-1/R)/(4*pi).
        // Storage gives M = a*(Meq^2-M^2), a=3*dt*0.05/R^2.
        let a = 0.15 * dt / radius.powi(2);
        let expected = 2.0 * a * equilibrium.powi(2)
            / (1.0 + (1.0 + 4.0 * a.powi(2) * equilibrium.powi(2)).sqrt());
        near(next.moisture()[0], expected, 1e-13);
        assert!(next.moisture()[0] > 0.0);
        near(report.boundary_water_kg[0], 0.0, 0.0);
        near(report.boundary_water_kg[1], next.water_mass_kg(), 1e-12);
        assert_eq!(old.moisture(), &[0.0]);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn radial_finite_surface_exchange_uses_actual_sphere_area() {
    with_cx(|cx| {
        let radius: f64 = 0.2;
        let h = 0.004;
        let d = 0.01;
        let dt = 2.0;
        let old = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![0.1],
            [d; 2],
            cx,
        );
        for rh in [0.2, 0.8] {
            let boundary = [
                MoistureFace::Sealed,
                MoistureFace::SourceSorption {
                    law: surface::law([h; 2]),
                    relative_humidity: rh,
                },
            ];
            let (next, report) = old
                .advance(dt, &boundary, config(), nonlinear(), cx)
                .unwrap();
            let equilibrium = 0.25 * rh;
            let bulk_speed = d / radius;
            let speed = 1.0 / (1.0 / bulk_speed + 1.0 / h);
            let alpha = dt * 3.0 * speed / radius;
            let expected = (0.1 + alpha * equilibrium) / (1.0 + alpha);
            let surface = (bulk_speed * expected + h * equilibrium) / (bulk_speed + h);
            near(next.moisture()[0], expected, 2e-14);
            near(report.surface_moisture[1].unwrap(), surface, 2e-14);
            let exchange = dt
                * 500.0
                * 4.0
                * core::f64::consts::PI
                * radius.powi(2)
                * h
                * (equilibrium - surface);
            near(report.boundary_water_kg[1], exchange, 2e-13);
            near(next.water_mass_kg() - old.water_mass_kg(), exchange, 2e-13);
            assert!(report.surface_sorption[1].is_some());
            assert!(report.surface_transfer[1].is_some());
        }
        for h in [0.0, 1e300] {
            let boundary = [
                MoistureFace::Sealed,
                MoistureFace::SourceSorption {
                    law: surface::law([h; 2]),
                    relative_humidity: 0.8,
                },
            ];
            let (next, report) = old
                .advance(dt, &boundary, config(), nonlinear(), cx)
                .unwrap();
            if h == 0.0 {
                assert_eq!(next.moisture(), old.moisture());
                assert_eq!(report.boundary_water_kg, [0.0; 2]);
            } else {
                let prescribed = [MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)];
                let expected = old
                    .advance(dt, &prescribed, config(), nonlinear(), cx)
                    .unwrap()
                    .0;
                near(next.moisture()[0], expected.moisture()[0], 2e-14);
                near(report.surface_moisture[1].unwrap(), 0.2, 2e-14);
            }
        }
    });
}

#[test]
fn radial_sealed_nonuniform_shells_conserve_volume_weighted_water() {
    with_cx(|cx| {
        let radii = vec![0.0, 0.01, 0.03, 0.06, 0.1];
        let values = vec![0.26, 0.1, 0.06, 0.18];
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
        ] {
            let old = specimen(symmetry, radii.clone(), values.clone(), [0.001, 0.01], cx);
            let (next, report) = old
                .advance(
                    0.2,
                    &[MoistureFace::Sealed, MoistureFace::Sealed],
                    config(),
                    nonlinear(),
                    cx,
                )
                .unwrap();
            near(next.water_mass_kg(), old.water_mass_kg(), 1e-11);
            assert_eq!(report.boundary_water_kg, [0.0; 2]);
            assert!(next.moisture().iter().all(|m| *m >= 0.06 && *m <= 0.26));
            assert!(next.moisture()[0] < values[0]);
            assert!(next.moisture()[2] > values[2]);
            let volume = match symmetry {
                RadialMoistureSymmetry::Sphere => {
                    4.0 * core::f64::consts::PI * 0.1_f64.powi(3) / 3.0
                }
                RadialMoistureSymmetry::Cylinder { length_m } => {
                    core::f64::consts::PI * 0.1_f64.powi(2) * length_m
                }
            };
            near(next.dry_mass_kg(), 500.0 * volume, 1e-13);
        }
    });
}

#[test]
fn radial_spherical_eigenmode_converges_under_space_and_time_refinement() {
    with_cx(|cx| {
        let radius: f64 = 0.1;
        let d = 0.01;
        let duration = 0.05;
        let exact = |r: f64, t: f64| {
            let x = core::f64::consts::PI * r / radius;
            0.05 + 0.1 * x.sin() / x
                * (-d * core::f64::consts::PI.powi(2) * t / radius.powi(2)).exp()
        };
        let mut errors = Vec::new();
        for n in [12, 24] {
            let radii: Vec<_> = (0..=n).map(|i| radius * i as f64 / n as f64).collect();
            let initial: Vec<_> = radii
                .windows(2)
                .map(|r| exact(0.5 * (r[0] + r[1]), 0.0))
                .collect();
            let mut state = specimen(RadialMoistureSymmetry::Sphere, radii, initial, [d; 2], cx);
            let steps = n * n / 6;
            for _ in 0..steps {
                state = state
                    .advance(
                        duration / steps as f64,
                        &[MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.05)],
                        config(),
                        nonlinear(),
                        cx,
                    )
                    .unwrap()
                    .0;
            }
            let error = state
                .moisture()
                .iter()
                .zip(state.centers_m())
                .zip(state.dry_masses_kg())
                .map(|((&m, &r), &mass)| mass * (m - exact(r, duration)).powi(2))
                .sum::<f64>();
            errors.push((error / state.dry_mass_kg()).sqrt());
        }
        assert!(
            errors[1] < errors[0] / 3.3,
            "radial refinement errors: {errors:?}"
        );
        assert!(errors[1] < 1e-4, "radial continuum error: {errors:?}");
    });
}

#[test]
fn radial_center_and_nonlinear_refusals_leave_the_checkpoint_unchanged() {
    with_cx(|cx| {
        let old = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.0],
            [0.0, 0.03],
            cx,
        );
        assert!(matches!(
            old.advance(0.01, &faces(0.2), config(), nonlinear(), cx),
            Err(MoistureError::Invalid(_))
        ));
        let mut limited = nonlinear();
        limited.maximum_iterations = 1;
        assert!(matches!(
            old.advance(
                0.01,
                &[MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.2)],
                config(),
                limited,
                cx
            ),
            Err(MoistureError::NonlinearNotConverged(_))
        ));
        assert_eq!(old.moisture(), &[0.0]);
        assert_eq!(old.water_mass_kg(), 0.0);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn radial_exact_equilibrium_keeps_zero_flux_at_extreme_time_scale() {
    with_cx(|cx| {
        let old = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.3_f64.sqrt()],
            vec![0.1],
            [1.0; 2],
            cx,
        );
        let (next, report) = old
            .advance(
                1e307,
                &[MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.1)],
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        assert_eq!(next.moisture(), old.moisture());
        assert_eq!(report.boundary_water_kg, [0.0; 2]);
        assert_eq!(report.water_balance_residual_kg, 0.0);
        assert_eq!(next.time_s(), 1e307);
    });
}

#[test]
fn radial_unrepresentable_total_wet_mass_refuses_transactionally() {
    with_cx(|cx| {
        let radius = (1.6e308 / 500.0 * (3.0 / (4.0 * core::f64::consts::PI))).cbrt();
        let old = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![0.0],
            [1.0; 2],
            cx,
        );
        // Half-cell spherical resistance gives alpha=3*D*dt/R^2=1.
        // The endpoint M=0.15 has finite water, but dry + water overflows.
        assert!(old.dry_mass_kg().is_finite());
        assert!((0.15 * old.dry_mass_kg()).is_finite());
        assert!(!(1.15 * old.dry_mass_kg()).is_finite());
        let mut step = config();
        step.water_balance_tolerance_kg = 1e296;
        assert!(matches!(
            old.advance(
                radius.powi(2) / 3.0,
                &[MoistureFace::Sealed, MoistureFace::PrescribedDryBasis(0.3)],
                step,
                nonlinear(),
                cx,
            ),
            Err(MoistureError::Invalid(_))
        ));
        assert_eq!(old.moisture(), &[0.0]);
        assert_eq!(old.wet_mass_kg(), old.dry_mass_kg());
        assert_eq!(old.time_s(), 0.0);
    });
}
