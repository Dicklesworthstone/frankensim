//! Principal-axis source potentials, dry fronts and real nonlinear surfaces.
//! All coefficients are authored verification data, not measured wood claims.
use super::*;

const PROPERTIES: [&str; 3] = ["nonlinear_x", "nonlinear_y", "nonlinear_z"];

fn material(slopes: [f64; 3]) -> MaterialCard {
    let base = card([0.0, 0.03], 0.3);
    let template = base.claims_for("moisture_diffusivity")[0].1.clone();
    let mut claims = base.claims().clone();
    for (key, slope) in PROPERTIES.into_iter().zip(slopes) {
        let mut claim = template.clone();
        claim.key = PropertyKey::new(key, MOISTURE_DIFFUSIVITY_DIMS);
        claim.value = PropertyValue::Curve {
            abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.0, 0.0), (0.3, 0.3 * slope)],
            dims: MOISTURE_DIFFUSIVITY_DIMS,
        };
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

#[test]
fn g1_six_anisotropic_faces_wet_a_dry_cell_with_the_correct_potential_flux() {
    with_cx(|cx| {
        let slopes = [0.1, 0.2, 0.05];
        let lengths = [1.0, 0.5, 0.25];
        let dt = 4.0;
        let equilibrium = [0.14, 0.2, 0.1, 0.18, 0.06, 0.16];
        let old = MoistureCartesian::new_kirchhoff(
            material(slopes),
            point(0.1),
            PROPERTIES,
            MaterialPropertySelection::SingleClaimOnly,
            grid_geometry(lengths, [1; 3]),
            vec![0.0],
            cx,
        )
        .unwrap();
        let boundaries = equilibrium.map(MoistureFace::PrescribedDryBasis);
        let (next, report) = old
            .advance(dt, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        // Phi_k(M)=slope_k*M²/2. Each physical half-cell contributes
        // dt*slope_k/length_k² * (M_reservoir²-M²) to this ONE inventory.
        let weights: [f64; 6] = std::array::from_fn(|side| {
            let axis = side / 2;
            dt * slopes[axis] / lengths[axis].powi(2)
        });
        let a = weights.iter().sum::<f64>();
        let b = weights
            .iter()
            .zip(equilibrium)
            .map(|(w, m)| w * m * m)
            .sum::<f64>();
        let expected = 2.0 * b / (1.0 + (1.0 + 4.0 * a * b).sqrt());
        near(next.moisture()[0], expected, 2e-12);
        for side in 0..6 {
            near(
                report.boundary_water_kg[side],
                old.dry_mass_kg() * weights[side] * (equilibrium[side].powi(2) - expected.powi(2)),
                2e-10,
            );
            assert_eq!(report.surface_moisture[side], vec![(0, equilibrium[side])]);
        }
        for (axis, key) in PROPERTIES.into_iter().enumerate() {
            assert_eq!(
                report.transport_start[0].property(key).unwrap().value_si(),
                0.0
            );
            near(
                report.transport_end[0].property(key).unwrap().value_si(),
                slopes[axis] * expected,
                3e-13,
            );
        }
        let sources = report.kirchhoff_sources.unwrap();
        assert_ne!(sources[0], sources[1]);
        assert_ne!(sources[1], sources[2]);
        assert!(report.linear.converged_euclidean());
        assert!(report.nonlinear_relative_residual < nonlinear().relative_tolerance);
        near(
            report.water_balance_residual_kg,
            0.0,
            config().water_balance_tolerance_kg,
        );
        assert_eq!(old.moisture(), &[0.0]);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(
            next.moisture(),
            old.advance(dt, &boundaries, config(), nonlinear(), cx)
                .unwrap()
                .0
                .moisture()
        );
    });
}

#[test]
fn g1_anisotropic_three_dimensional_interior_transport_wets_a_dry_corner() {
    with_cx(|cx| {
        let slopes = [0.1, 0.2, 0.05];
        let lengths = [1.0, 0.5, 0.25];
        let squared_base: f64 = 0.0004;
        let squared_increments = [0.0016, 0.0032, 0.0064];
        // M² is affine in all three coordinates. With two sealed cells per
        // axis, its Kirchhoff divergence is the signed sum of three constants.
        let rate: [f64; 3] = std::array::from_fn(|k| {
            0.5 * slopes[k] * squared_increments[k] / (lengths[k] / 2.0_f64).powi(2)
        });
        let dt = squared_base.sqrt() / rate.iter().sum::<f64>();
        let expected = (0..8)
            .map(|i| {
                (squared_base
                    + (0..3)
                        .map(|k| squared_increments[k] * ((i >> k) & 1) as f64)
                        .sum::<f64>())
                .sqrt()
            })
            .collect::<Vec<_>>();
        let initial = (0..8)
            .map(|i| {
                if i == 0 {
                    return 0.0;
                }
                expected[i]
                    - dt * (0..3)
                        .map(|k| rate[k] * (1.0 - 2.0 * ((i >> k) & 1) as f64))
                        .sum::<f64>()
            })
            .collect::<Vec<_>>();
        let old = MoistureCartesian::new_kirchhoff(
            material(slopes),
            point(0.1),
            PROPERTIES,
            MaterialPropertySelection::SingleClaimOnly,
            grid_geometry(lengths, [2; 3]),
            initial,
            cx,
        )
        .unwrap();
        let (next, report) = old
            .advance(dt, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        for (&actual, &expected) in next.moisture().iter().zip(&expected) {
            near(actual, expected, 3e-12);
        }
        assert!(next.moisture()[0] > 0.0);
        assert_eq!(report.boundary_water_kg, [0.0; 6]);
        near(next.water_mass_kg(), old.water_mass_kg(), 2e-12);
        assert!(report.linear.converged_euclidean());
    });
}

#[test]
fn g1_zero_transverse_laws_leave_independent_dry_fronts_equal_to_native_slabs() {
    with_cx(|cx| {
        let shape = grid_geometry([1.0, 0.2, 0.1], [8, 3, 2]);
        let [nx, ny, nz] = shape.cells;
        let initial = (0..nx * ny * nz)
            .map(|i| {
                if i % nx < 3 {
                    0.12 + 0.01 * (i / nx) as f64
                } else {
                    0.0
                }
            })
            .collect::<Vec<_>>();
        let old = MoistureCartesian::new_kirchhoff(
            material([0.1, 0.0, 0.0]),
            point(0.1),
            PROPERTIES,
            MaterialPropertySelection::SingleClaimOnly,
            shape,
            initial.clone(),
            cx,
        )
        .unwrap();
        let mut boundaries = sealed();
        boundaries[0] = MoistureFace::PrescribedDryBasis(0.2);
        // Strong transverse reservoirs must not couple rows through a law
        // that vanishes identically, including finite-emission half-cells.
        boundaries[2] = MoistureFace::PrescribedDryBasis(0.3);
        boundaries[3] = MoistureFace::PrescribedDryBasis(0.0);
        boundaries[4] = MoistureFace::SourceSorption {
            law: super::super::surface::law([0.004; 2]),
            relative_humidity: 0.8,
        };
        boundaries[5] = MoistureFace::SourceSorption {
            law: super::super::surface::law([0.0; 2]),
            relative_humidity: 0.2,
        };
        let (next, report) = old
            .advance(0.2, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        let mut received = 0.0;
        for row in 0..ny * nz {
            let slab = MoistureSlab::new_kirchhoff(
                card([0.0, 0.03], 0.3),
                point(0.1),
                "moisture_diffusivity",
                MaterialPropertySelection::SingleClaimOnly,
                MoistureSlabGeometry {
                    length_m: 1.0,
                    cross_section_m2: shape.lengths_m[1] * shape.lengths_m[2] / (ny * nz) as f64,
                    maximum_moisture: 0.3,
                    maximum_cells: 256,
                },
                initial[row * nx..(row + 1) * nx].to_vec(),
                nonlinear(),
                cx,
            )
            .unwrap();
            let (expected, slab_report) = slab
                .advance(
                    0.2,
                    &[boundaries[0].clone(), MoistureFace::Sealed],
                    config(),
                    cx,
                )
                .unwrap();
            for i in 0..nx {
                near(next.moisture()[row * nx + i], expected.moisture()[i], 3e-12);
            }
            assert!(
                next.moisture()[row * nx + 3] > 0.0,
                "interior dry front must wet"
            );
            received += slab_report.boundary_water_kg[0];
        }
        near(report.boundary_water_kg[0], received, 3e-12);
        assert_eq!(&report.boundary_water_kg[1..], &[0.0; 5]);
        near(next.water_mass_kg() - old.water_mass_kg(), received, 3e-12);
        assert!(
            next.moisture()
                .iter()
                .all(|m| m.is_finite() && *m >= 0.0 && *m <= 0.2)
        );
    });
}

#[test]
fn g1_each_principal_half_cell_balances_nonlinear_diffusion_and_surface_emission() {
    with_cx(|cx| {
        let slopes = [0.1, 0.2, 0.05];
        let lengths = [1.0, 0.5, 0.25];
        let endpoint: f64 = 0.12;
        let surfaces: [f64; 6] = [0.10, 0.14, 0.11, 0.13, 0.105, 0.135];
        let h = 0.2;
        let dt = 0.5;
        let flux: [f64; 6] = std::array::from_fn(|side| {
            let axis = side / 2;
            slopes[axis] / lengths[axis] * (surfaces[side].powi(2) - endpoint.powi(2))
        });
        let equilibrium: [f64; 6] = std::array::from_fn(|side| surfaces[side] + flux[side] / h);
        let increment = (0..6)
            .map(|side| dt / lengths[side / 2] * flux[side])
            .sum::<f64>();
        // Manufacture one common endpoint and six distinct physical interfaces.
        // Old moisture comes directly from their signed conservation equation.
        let old = MoistureCartesian::new(
            material(slopes),
            point(0.1),
            PROPERTIES,
            MaterialPropertySelection::SingleClaimOnly,
            grid_geometry(lengths, [1; 3]),
            vec![endpoint - increment],
            cx,
        )
        .unwrap()
        .with_kirchhoff_diffusivity();
        let boundaries = equilibrium.map(|m| MoistureFace::SourceSorption {
            law: super::super::surface::law([h; 2]),
            relative_humidity: m / 0.25,
        });
        let (next, report) = old
            .advance(dt, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        near(next.moisture()[0], endpoint, 2e-12);
        for side in 0..6 {
            let actual = report.surface_moisture[side][0].1;
            near(actual, surfaces[side], 2e-12);
            let expected = old.dry_mass_kg() * dt / lengths[side / 2] * flux[side];
            near(report.boundary_water_kg[side], expected, 2e-11);
            near(
                report.boundary_water_kg[side],
                old.dry_mass_kg() * dt / lengths[side / 2] * h * (equilibrium[side] - actual),
                2e-11,
            );
            assert!(report.surface_sorption[side].is_some());
            assert!(report.surface_transfer[side].is_some());
        }
        near(
            report.water_balance_residual_kg,
            0.0,
            config().water_balance_tolerance_kg,
        );
    });
}

#[test]
fn g0_bounded_cartesian_nonlinear_refusal_keeps_the_original_dry_inventory() {
    with_cx(|cx| {
        let old = MoistureCartesian::new_kirchhoff(
            material([0.1, 0.2, 0.05]),
            point(0.1),
            PROPERTIES,
            MaterialPropertySelection::SingleClaimOnly,
            grid_geometry([1.0; 3], [2; 3]),
            vec![0.0; 8],
            cx,
        )
        .unwrap();
        let boundaries = std::array::from_fn(|_| MoistureFace::PrescribedDryBasis(0.2));
        let controls = MoistureNonlinearConfig {
            maximum_iterations: 1,
            ..nonlinear()
        };
        assert!(matches!(
            old.advance(1.0, &boundaries, config(), controls, cx),
            Err(MoistureError::SpatialNotConverged { iterations: 1, .. })
        ));
        assert_eq!(old.moisture(), &[0.0; 8]);
        assert_eq!(old.water_mass_kg(), 0.0);
        assert_eq!(old.time_s(), 0.0);
    });
}
