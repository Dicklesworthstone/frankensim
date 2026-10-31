//! Actual 3D eigenmodes, six physical faces and native differential swelling.
//! All material coefficients are authored verification data, not measured wood.
use super::*;
use fs_blake3::ContentHash;
use fs_couple::moisture::cartesian::*;
use fs_couple::moisture::solid::*;
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig, source::SolidStiffnessModel};

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}
fn grid_geometry(lengths_m: [f64; 3], cells: [usize; 3]) -> MoistureCartesianGeometry {
    MoistureCartesianGeometry {
        lengths_m,
        cells,
        frame: ContentHash([4; 32]),
        maximum_moisture: 0.3,
        maximum_cells: 256,
    }
}
fn directional_card(values: [f64; 3]) -> MaterialCard {
    let base = card([0.01; 2], 0.3);
    let template = base.claims_for("moisture_diffusivity")[0].1.clone();
    let mut claims = base.claims().clone();
    for (key, value) in ["diffusion_x", "diffusion_y", "diffusion_z"]
        .into_iter()
        .zip(values)
    {
        let mut claim = template.clone();
        claim.key = PropertyKey::new(key, MOISTURE_DIFFUSIVITY_DIMS);
        claim.value = PropertyValue::Curve {
            abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.0, value), (0.3, value)],
            dims: MOISTURE_DIFFUSIVITY_DIMS,
        };
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}
fn grid(
    material: MaterialCard,
    geometry: MoistureCartesianGeometry,
    initial: Vec<f64>,
    cx: &Cx<'_>,
) -> MoistureCartesian {
    MoistureCartesian::new(
        material,
        point(0.1),
        ["moisture_diffusivity"; 3],
        MaterialPropertySelection::SingleClaimOnly,
        geometry,
        initial,
        cx,
    )
    .unwrap()
}
fn sealed() -> [MoistureFace; 6] {
    std::array::from_fn(|_| MoistureFace::Sealed)
}

#[test]
fn g1_three_dimensional_principal_diffusion_has_exact_discrete_cosine_decay() {
    with_cx(|cx| {
        let geometry = grid_geometry([0.12, 0.09, 0.07], [4, 3, 5]);
        let d = [5e-6, 2e-5, 5e-5];
        let [nx, ny, nz] = geometry.cells;
        let wave = |i: usize| {
            [i % nx, (i / nx) % ny, i / (nx * ny)]
                .into_iter()
                .zip(geometry.cells)
                .map(|(j, n)| (std::f64::consts::PI * (j as f64 + 0.5) / n as f64).cos())
                .product::<f64>()
        };
        let old = MoistureCartesian::new(
            directional_card(d),
            point(0.1),
            ["diffusion_x", "diffusion_y", "diffusion_z"],
            MaterialPropertySelection::SingleClaimOnly,
            geometry,
            (0..nx * ny * nz).map(|i| 0.15 + 0.03 * wave(i)).collect(),
            cx,
        )
        .unwrap();
        let dt = 0.7;
        let eigenvalue = (0..3)
            .map(|k| {
                let dx = geometry.lengths_m[k] / geometry.cells[k] as f64;
                4.0 * d[k] / (dx * dx)
                    * (std::f64::consts::PI / (2.0 * geometry.cells[k] as f64))
                        .sin()
                        .powi(2)
            })
            .sum::<f64>();
        let (next, report) = old
            .advance(dt, &sealed(), config(), nonlinear(), cx)
            .unwrap();
        for (i, &m) in next.moisture().iter().enumerate() {
            near(m, 0.15 + 0.03 * wave(i) / (1.0 + dt * eigenvalue), 3e-13);
        }
        near(next.water_mass_kg(), old.water_mass_kg(), 1e-13);
        assert_eq!(report.boundary_water_kg, [0.0; 6]);
        assert_eq!(report.nonlinear_iterations, 1);
        assert!(report.linear.converged_euclidean());
        assert_eq!(next.cell_index([3, 2, 4]), Some(59));
        assert_eq!(next.cell_index([4, 0, 0]), None);
    });
}

#[test]
fn g1_single_cell_resolves_all_six_half_cell_fluxes_without_double_counting_storage() {
    with_cx(|cx| {
        let geometry = grid_geometry([1.0, 2.0, 0.5], [1; 3]);
        let d = [0.01, 0.02, 0.04];
        let old = MoistureCartesian::new(
            directional_card(d),
            point(0.1),
            ["diffusion_x", "diffusion_y", "diffusion_z"],
            MaterialPropertySelection::SingleClaimOnly,
            geometry,
            vec![0.1],
            cx,
        )
        .unwrap();
        let reservoir = [0.08, 0.21, 0.12, 0.17, 0.09, 0.24];
        let dt = 2.0;
        let coupling: [f64; 6] = std::array::from_fn(|side| {
            2.0 * dt * d[side / 2] / geometry.lengths_m[side / 2].powi(2)
        });
        let expected = (0.1
            + coupling
                .iter()
                .zip(reservoir)
                .map(|(a, m)| a * m)
                .sum::<f64>())
            / (1.0 + coupling.iter().sum::<f64>());
        let (next, report) = old
            .advance(
                dt,
                &reservoir.map(MoistureFace::PrescribedDryBasis),
                config(),
                nonlinear(),
                cx,
            )
            .unwrap();
        near(next.moisture()[0], expected, 1e-14);
        near(next.dry_mass_kg(), 500.0, 1e-12);
        for side in 0..6 {
            near(
                report.boundary_water_kg[side],
                old.dry_mass_kg() * coupling[side] * (reservoir[side] - expected),
                2e-12,
            );
            assert_eq!(report.surface_moisture[side], vec![(0, reservoir[side])]);
        }
        near(
            next.water_mass_kg() - old.water_mass_kg(),
            report.boundary_water_kg.iter().sum(),
            2e-12,
        );
    });
}

#[test]
fn g1_three_dimensional_nonlinear_endpoint_and_refusal_use_original_sources() {
    with_cx(|cx| {
        let old = grid(
            card([0.01, 0.04], 0.3),
            grid_geometry([1.0; 3], [1; 3]),
            vec![0.1],
            cx,
        );
        let faces = std::array::from_fn(|_| MoistureFace::PrescribedDryBasis(0.2));
        // Six half-cell faces give M-.1=12*dt*(.01+.1*M)*(.2-M).
        // With dt=10/3, the positive root of 4M²+.6M-.18 is .15.
        let (next, report) = old
            .advance(10.0 / 3.0, &faces, config(), nonlinear(), cx)
            .unwrap();
        near(next.moisture()[0], 0.15, 1e-12);
        near(
            report.transport_end[0]
                .property("moisture_diffusivity")
                .unwrap()
                .value_si(),
            0.025,
            1e-13,
        );
        assert!(report.nonlinear_iterations > 1);
        assert!(report.nonlinear_relative_residual < nonlinear().relative_tolerance);
        assert!(matches!(
            old.advance(
                10.0 / 3.0,
                &faces,
                config(),
                MoistureNonlinearConfig {
                    maximum_iterations: 1,
                    ..nonlinear()
                },
                cx
            ),
            Err(MoistureError::SpatialNotConverged { iterations: 1, .. })
        ));
        assert_eq!(old.moisture(), &[0.1]);
        assert_eq!(old.time_s(), 0.0);
        let replay = old
            .advance(10.0 / 3.0, &faces, config(), nonlinear(), cx)
            .unwrap();
        assert_eq!(replay.0.moisture(), next.moisture());
    });
}

#[test]
fn g1_each_face_retains_source_sorption_and_finite_surface_resistance() {
    with_cx(|cx| {
        let geometry = grid_geometry([0.2, 0.1, 0.05], [1; 3]);
        let diffusion = 0.01;
        let old = grid(card([diffusion; 2], 0.3), geometry, vec![0.1], cx);
        let h = 0.004;
        let dt = 2.0;
        let faces = std::array::from_fn(|side| MoistureFace::SourceSorption {
            law: super::surface::law([if side == 0 { 0.0 } else { h }; 2]),
            relative_humidity: 0.8,
        });
        let coupling: [f64; 6] = std::array::from_fn(|side| {
            if side == 0 {
                return 0.0;
            }
            let dx = geometry.lengths_m[side / 2];
            dt / dx / (1.0 / h + dx / (2.0 * diffusion))
        });
        let alpha = coupling.iter().sum::<f64>();
        let expected = (0.1 + alpha * 0.2) / (1.0 + alpha);
        let (next, report) = old.advance(dt, &faces, config(), nonlinear(), cx).unwrap();
        near(next.moisture()[0], expected, 1e-14);
        for side in 0..6 {
            let dx = geometry.lengths_m[side / 2];
            let speed = 2.0 * diffusion / dx;
            let actual_h = if side == 0 { 0.0 } else { h };
            let surface = (speed * expected + actual_h * 0.2) / (speed + actual_h);
            near(report.surface_moisture[side][0].1, surface, 1e-14);
            near(
                report.boundary_water_kg[side],
                old.dry_mass_kg() * dt / dx * actual_h * (0.2 - surface),
                1e-13,
            );
            assert!(report.surface_sorption[side].is_some());
            assert!(report.surface_transfer[side].is_some());
        }
        assert_eq!(report.boundary_water_kg[0], 0.0);
        near(report.water_balance_residual_kg, 0.0, 1e-13);
    });
}

#[test]
fn g1_transversely_uniform_cartesian_transport_reduces_to_the_original_slab() {
    with_cx(|cx| {
        let material = card([0.01, 0.04], 0.3);
        let geometry = grid_geometry([0.1, 0.05, 0.02], [4, 3, 2]);
        let contents = vec![0.1, 0.12, 0.15, 0.18];
        let old = grid(
            material.clone(),
            geometry,
            (0..24).map(|i| contents[i % 4]).collect(),
            cx,
        );
        let slab = slab(material, 0.1, contents, cx)
            .with_implicit_diffusivity(nonlinear())
            .unwrap();
        let mut boundaries = sealed();
        boundaries[0] = MoistureFace::PrescribedDryBasis(0.2);
        boundaries[1] = MoistureFace::PrescribedDryBasis(0.05);
        let (next, report) = old
            .advance(0.001, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        let (slab_next, slab_report) = slab
            .advance(
                0.001,
                &[
                    MoistureFace::PrescribedDryBasis(0.2),
                    MoistureFace::PrescribedDryBasis(0.05),
                ],
                config(),
                cx,
            )
            .unwrap();
        for (i, &m) in next.moisture().iter().enumerate() {
            near(m, slab_next.moisture()[i % 4], 2e-12);
        }
        near(next.water_mass_kg(), slab_next.water_mass_kg(), 2e-13);
        for side in 0..2 {
            near(
                report.boundary_water_kg[side],
                slab_report.boundary_water_kg[side],
                2e-13,
            );
        }
    });
}

#[test]
fn g3_through_thickness_uptake_reaches_each_original_solid_cell_and_compatible_shape() {
    with_cx(|cx| {
        let geometry = grid_geometry([0.02; 3], [2; 3]);
        let old = grid(card([0.01; 2], 0.3), geometry, vec![0.1; 8], cx);
        let axes = HygroscopicAxes {
            principal_to_mesh: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            principal_frame: ContentHash([3; 32]),
            mesh_frame: geometry.frame,
        };
        // Clamp the actual x=0 face in the native x-fastest 3x3x3 node lattice.
        let fixed = (0..27)
            .filter(|i| i % 3 == 0)
            .flat_map(|i| [3 * i, 3 * i + 1, 3 * i + 2])
            .collect();
        let controls = MoistureSolidControls {
            assembly: TetAssemblyBudget {
                maximum_nodes: 64,
                maximum_tetrahedra: 96,
                maximum_free_dofs: 192,
                minimum_scaled_jacobian: 1e-8,
            },
            equilibrium: TetStaticSolveConfig {
                relative_residual_tolerance: 1e-10,
                maximum_iterations: 4096,
            },
            maximum_absolute_mandel_strain: 0.02,
            maximum_nodal_displacement_m: 0.01,
            mass_tolerance_kg: 1e-12,
        };
        let solid = MoistureSolid::from_cartesian(
            &old,
            SolidStiffnessModel::Isotropic,
            MaterialPropertySelection::SingleClaimOnly,
            MaterialPropertySelection::SingleClaimOnly,
            axes,
            fixed,
            controls,
            cx,
        )
        .unwrap();
        let mut boundaries = sealed();
        boundaries[5] = MoistureFace::PrescribedDryBasis(0.2);
        let step = solid
            .advance_cartesian(&old, 0.005, &boundaries, config(), nonlinear(), cx)
            .unwrap();
        assert_eq!(solid.reference_nodes_m().len(), 27);
        assert_eq!(solid.tetrahedra().len(), 48);
        // Independent two-layer solve: [1.5,-.5;-.5,2.5] M = [.1,.3].
        for i in 0..8 {
            let expected = if i < 4 { 0.4 / 3.5 } else { 0.5 / 3.5 };
            near(step.moisture.moisture()[i], expected, 1e-13);
            near(
                step.solid.expansion_paths()[i].free_principal_strains()[2],
                0.02 * (expected - 0.1),
                1e-13,
            );
            for tet in 6 * i..6 * (i + 1) {
                assert!(step.solid.current_density_kg_m3()[tet].is_finite());
            }
        }
        assert!(
            step.solid
                .solution()
                .displacement_m()
                .iter()
                .flatten()
                .any(|u| u.abs() > 1e-7)
        );
        assert!(
            step.solid
                .stress()
                .elements()
                .iter()
                .any(|e| e.stress_mandel_pa().iter().any(|s| s.abs() > 1.0))
        );
        near(step.solid.wet_mass_kg(), step.moisture.wet_mass_kg(), 1e-13);
        assert_eq!(step.solid.time_s(), step.moisture.time_s());
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.moisture(), &[0.1; 8]);
    });
}
