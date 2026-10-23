//! G0/G1 analytic and manufactured-field checks; synthetic source claims only.
mod common;
use common::*;
use fs_em::{
    planar::{MagneticError, PlanarMagnetostatics},
    transfer::{CellSample, PlanarFrame},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{PropertyValue, QueryPoint};
use fs_material::field_storage::{MAGNETIC_PERMEABILITY_DIMS, ResolvedInductor};
use fs_rep_mesh::{Metric2, TriComplex2};

#[test]
fn affine_field_and_energy_share_the_source_used_by_circuit_storage() {
    let mesh = TriComplex2::from_indexed_triangles(
        "affine-field",
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        vec![[0, 1, 2]],
        Metric2::planar(0.2).unwrap(),
    )
    .unwrap();
    let source = material(2e-6);
    let circuit = ResolvedInductor::try_from_material(source.clone(), 1.0, 1.0, 1).unwrap();
    let problem =
        PlanarMagnetostatics::new(mesh, vec![source.clone()], limits(), &mut || false).unwrap();
    near(problem.permeability_h_m()[0], circuit.inductance_h());
    assert_eq!(&problem.materials()[0], &source);
    let field = problem
        .solve(
            &[0.0],
            &[(0, 7.0), (1, 9.0), (2, 10.0)],
            config(),
            &mut || false,
        )
        .unwrap();
    near(field.magnetic_flux_density_t()[0][0], 3.0);
    near(field.magnetic_flux_density_t()[0][1], -2.0);
    near(field.energy_j(), 0.5 * 0.1 * 13.0 / 2e-6);
    assert!(field.report().is_none());
    let frame = PlanarFrame::new([1.0, 2.0, 3.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let sampled = field
        .sample_cells(
            &frame,
            &[CellSample {
                triangle: 0,
                position_m: [1.1, 2.2, 3.2],
            }],
            1e-12,
            1,
            &mut || false,
        )
        .unwrap();
    for (got, expected) in sampled[0].into_iter().zip([0.0, 3.0, -2.0]) {
        near(got, expected);
    }
}

#[test]
fn changing_only_one_material_changes_the_solved_flux_and_preserves_tangential_h() {
    let solve = |left| {
        PlanarMagnetostatics::new(
            layered(1.0, 1.0, 0.5),
            vec![
                material(left),
                material(left),
                material(1e-6),
                material(1e-6),
            ],
            limits(),
            &mut || false,
        )
        .unwrap()
    };
    let a = solve(1e-6);
    let b = solve(2e-6);
    let boundary = ends(6.0);
    let fa = a
        .solve(&[0.0; 4], &boundary, config(), &mut || false)
        .unwrap();
    let fb = b
        .solve(&[0.0; 4], &boundary, config(), &mut || false)
        .unwrap();
    for cell in 0..4 {
        near(fa.magnetic_flux_density_t()[cell][0], 0.0);
        near(fa.magnetic_flux_density_t()[cell][1], -3.0);
        near(
            fb.magnetic_flux_density_t()[cell][1],
            if cell < 2 { -4.0 } else { -2.0 },
        );
        near(fb.magnetic_field_a_m()[cell][1], -2e6);
    }
    near(fa.energy_j(), 4.5e6);
    near(fb.energy_j(), 3e6);
    assert!(fb.report().unwrap().converged_euclidean());
}

#[test]
fn impressed_current_uses_consistent_cell_loads_and_permeability() {
    let mesh = || {
        TriComplex2::from_indexed_triangles(
            "impressed-current",
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
            vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
            Metric2::planar(0.25).unwrap(),
        )
        .unwrap()
    };
    for mu in [1e-6, 2e-6] {
        let problem =
            PlanarMagnetostatics::new(mesh(), vec![material(mu); 4], limits(), &mut || false)
                .unwrap();
        let bc = [(0, 0.0), (1, 0.0), (2, 0.0), (3, 0.0)];
        let field = problem
            .solve(&[12e6; 4], &bc, config(), &mut || false)
            .unwrap();
        // Four triangles: K_44=4*depth/mu and b_4=J*depth/3.
        near(field.potential_wb_m()[4], mu * 1e6);
        assert!(matches!(
            problem.solve(
                &[12e6; 4],
                &bc,
                fs_solver::CheckedCgConfig {
                    max_iterations: 0,
                    ..config()
                },
                &mut || false
            ),
            Err(MagneticError::Nonconvergence(_))
        ));
    }
}

fn manufactured_error(n: usize) -> f64 {
    let mut vertices = Vec::new();
    for y in 0..=n {
        for x in 0..=n {
            vertices.push([x as f64 / n as f64, y as f64 / n as f64]);
        }
    }
    let mut faces = Vec::new();
    for y in 0..n {
        for x in 0..n {
            let a = (y * (n + 1) + x) as u32;
            let b = a + 1;
            let d = a + (n + 1) as u32;
            let c = d + 1;
            faces.extend([[a, b, c], [a, c, d]]);
        }
    }
    let currents: Vec<_> = faces
        .iter()
        .map(|face| {
            let [x, y] =
                [0, 1].map(|k| face.iter().map(|&v| vertices[v as usize][k]).sum::<f64>() / 3.0);
            2e6 * (y * (1.0 - y) + x * (1.0 - x))
        })
        .collect();
    let boundary: Vec<_> = vertices
        .iter()
        .enumerate()
        .filter_map(|(i, &[x, y])| {
            (x == 0.0 || x == 1.0 || y == 0.0 || y == 1.0).then_some((i, 0.0))
        })
        .collect();
    let mesh = TriComplex2::from_indexed_triangles(
        "manufactured-magnetic",
        vertices,
        faces,
        Metric2::planar(0.3).unwrap(),
    )
    .unwrap();
    let materials = vec![material(1e-6); mesh.faces().len()];
    let problem = PlanarMagnetostatics::new(mesh, materials, limits(), &mut || false).unwrap();
    let field = problem
        .solve(&currents, &boundary, config(), &mut || false)
        .unwrap();
    field
        .potential_wb_m()
        .iter()
        .zip(problem.mesh().vertices())
        .map(|(&a, &[x, y])| (a - x * (1.0 - x) * y * (1.0 - y)).abs())
        .fold(0.0_f64, f64::max)
}
#[test]
fn manufactured_poisson_field_converges_with_mesh_refinement() {
    let coarse = manufactured_error(4);
    let fine = manufactured_error(8);
    assert!(fine < coarse / 3.0, "coarse={coarse}, fine={fine}");
}

#[test]
fn gauge_source_domain_metric_and_transfer_refusals_are_explicit() {
    let mesh = layered(1.0, 1.0, 0.5);
    let p = PlanarMagnetostatics::new(mesh.clone(), vec![material(1e-6); 4], limits(), &mut || {
        false
    })
    .unwrap();
    assert!(p.solve(&[0.0; 4], &[], config(), &mut || false).is_err());
    assert!(matches!(
        p.solve(&[0.0; 4], &ends(1.0), config(), &mut || true),
        Err(MagneticError::Cancelled)
    ));
    let field = p
        .solve(&[0.0; 4], &ends(1.0), config(), &mut || false)
        .unwrap();
    let frame = PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    for position in [[1.5, 0.5, 0.2], [0.5, 0.5, 0.6], [f64::NAN, 0.5, 0.2]] {
        assert!(
            field
                .sample_cells(
                    &frame,
                    &[CellSample {
                        triangle: 0,
                        position_m: position
                    }],
                    1e-12,
                    1,
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(PlanarFrame::new([0.0; 3], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]).is_err());
    let nonlinear = source(
        PropertyValue::Curve {
            abscissa: "B".into(),
            abscissa_dims: fs_qty::Dims([0, 1, -2, 0, -1, 0]),
            knots: vec![(0.0, 1e-6), (1.0, 2e-6)],
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        ValidityDomain::unconstrained().with("B", 0.0, 1.0),
    );
    let state = resolve(&nonlinear, &QueryPoint::new().with("B", 0.5).unwrap());
    assert!(matches!(
        PlanarMagnetostatics::new(mesh, vec![state; 4], limits(), &mut || false),
        Err(MagneticError::Material(_))
    ));
    let axisymmetric = TriComplex2::from_indexed_triangles(
        "axisymmetric-refusal",
        vec![[1.0, 0.0], [2.0, 0.0], [1.0, 1.0]],
        vec![[0, 1, 2]],
        Metric2::axisymmetric(1.0).unwrap(),
    )
    .unwrap();
    assert!(
        PlanarMagnetostatics::new(axisymmetric, vec![material(1e-6)], limits(), &mut || false)
            .is_err()
    );
}

#[test]
fn temperature_source_refresh_changes_actual_fields_without_changing_the_law() {
    let card = source(
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, 1e-6), (400.0, 2e-6)],
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        ValidityDomain::unconstrained().with("T", 300.0, 400.0),
    );
    let prepare = |temperature| {
        let state = resolve(&card, &QueryPoint::new().with("T", temperature).unwrap());
        PlanarMagnetostatics::new(
            layered(1.0, 1.0, 0.5),
            vec![state.clone(), state, material(1e-6), material(1e-6)],
            limits(),
            &mut || false,
        )
        .unwrap()
    };
    let cold = prepare(300.0);
    let hot = prepare(400.0);
    let bc = ends(6.0);
    let a = cold.solve(&[0.0; 4], &bc, config(), &mut || false).unwrap();
    let b = hot.solve(&[0.0; 4], &bc, config(), &mut || false).unwrap();
    assert_eq!(
        cold.materials()[0].card_identity(),
        hot.materials()[0].card_identity()
    );
    assert_ne!(cold.materials()[0], hot.materials()[0]);
    near(a.magnetic_flux_density_t()[0][1], -3.0);
    near(b.magnetic_flux_density_t()[0][1], -4.0);
}

#[test]
fn disconnected_components_cannot_hide_an_unfixed_gauge() {
    let mesh = TriComplex2::from_indexed_triangles(
        "disconnected-gauges",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [2.0, 1.0],
        ],
        vec![[0, 1, 2], [3, 4, 5]],
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap();
    let problem =
        PlanarMagnetostatics::new(mesh, vec![material(1e-6); 2], limits(), &mut || false).unwrap();
    assert!(
        problem
            .solve(&[0.0; 2], &[(0, 0.0)], config(), &mut || false)
            .is_err()
    );
    assert!(
        problem
            .solve(&[0.0; 2], &[(0, 0.0), (1, 0.0)], config(), &mut || false)
            .is_err()
    );
    let field = problem
        .solve(
            &[0.0; 2],
            &[(0, 0.0), (1, 0.0), (3, 0.0), (4, 0.0)],
            config(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(field.energy_j(), 0.0);
}
