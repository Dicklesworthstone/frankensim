//! A bath on one solid cannot set another insulated solid's temperature.
mod support;

use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::{TransientConfig, TransientProblem, VolumetricHeatCapacity, march};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    DofMap, InitialGuess, Nonlinearity, ScalarField, SolveConfig, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder, solve,
};
use fs_rep_mesh::TetComplex;
use support::with_cx;

fn separated_cubes() -> ConductionMesh {
    let (first, mut positions) = box_grid([1, 1, 1], [1.0; 3]);
    let offset = positions.len() as u32;
    let extra: Vec<_> = positions.iter().map(|p| [p[0] + 3.0, p[1], p[2]]).collect();
    positions.extend(extra);
    let mut tets = first.tets.clone();
    tets.extend(first.tets.into_iter().map(|tet| tet.map(|v| v + offset)));
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}

fn left_bath(mesh: &ConductionMesh, robin: bool) -> ThermalBoundary {
    let bath = if robin {
        ThermalBc::robin(2.0, 300.0).unwrap()
    } else {
        ThermalBc::dirichlet(300.0).unwrap()
    };
    ThermalBoundaryBuilder::new(mesh)
        .region("left bath", |f| on_box_face(f.centroid[0], 0.0), bath)
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

#[test]
fn one_bath_cannot_hide_a_floating_body_at_zero_residual() {
    let mesh = separated_cubes();
    let source = ScalarField::Uniform(0.0);
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 2.0), (450.0, 4.0)]).unwrap(),
    );
    for robin in [false, true] {
        let boundary = left_bath(&mesh, robin);
        let dofs = DofMap::new(&boundary, mesh.vertex_count()).unwrap();
        for newton in [false, true] {
            for floating_temperature in [300.0, 340.0] {
                let mut config = SolveConfig::default();
                if newton {
                    config.nonlinearity = Nonlinearity::Newton {
                        line_search: Default::default(),
                    };
                }
                config.initial = InitialGuess::Free(
                    dofs.free()
                        .iter()
                        .map(|&v| {
                            if mesh.positions()[v][0] < 2.0 {
                                300.0
                            } else {
                                floating_temperature
                            }
                        })
                        .collect(),
                );
                let error = with_cx(|cx| {
                    solve(
                        cx,
                        ConductionProblem {
                            mesh: &mesh,
                            boundary: &boundary,
                            material: &material,
                            element_materials: None,
                            source: &source,
                        },
                        config,
                    )
                    .unwrap_err()
                });
                assert!(matches!(error, ConductionError::SingularPureNeumann));
            }
        }
    }
}

#[test]
fn conductivity_design_cannot_invert_a_floating_component() {
    let mesh = separated_cubes();
    let boundary = left_bath(&mesh, false);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let design = ConductivityDesign::new(
        ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        },
        SolveConfig::default().linear,
    )
    .unwrap();
    let error = with_cx(|cx| {
        design
            .solve(cx, &vec![1.0; mesh.element_count()])
            .unwrap_err()
    });
    assert!(matches!(error, ConductionError::SingularPureNeumann));
}

#[test]
fn separate_dirichlet_and_spatial_robin_baths_set_independent_temperatures() {
    let mesh = separated_cubes();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    for robin in [false, true] {
        let first = if robin {
            ThermalBc::robin(2.0, 300.0).unwrap()
        } else {
            ThermalBc::dirichlet(300.0).unwrap()
        };
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("left bath", |f| on_box_face(f.centroid[0], 0.0), first)
            .unwrap()
            .region(
                "right bath",
                |f| on_box_face(f.centroid[0], 4.0),
                ThermalBc::Robin {
                    htc: ScalarField::Nodal(
                        mesh.positions()
                            .iter()
                            .map(|p| 1.0 + p[1] + 2.0 * p[2])
                            .collect(),
                    ),
                    t_ref: ScalarField::Uniform(340.0),
                },
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let solution = with_cx(|cx| {
            solve(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                },
                SolveConfig {
                    initial: InitialGuess::Uniform(320.0),
                    ..Default::default()
                },
            )
            .unwrap()
        });
        for (p, &actual) in mesh.positions().iter().zip(&solution.temperature) {
            let expected = if p[0] < 2.0 { 300.0 } else { 340.0 };
            assert!(
                (actual - expected).abs() < 1e-7,
                "{actual} versus {expected}"
            );
        }
        assert!(solution.report.energy.closure_w.abs() < 1e-7);
    }
}

#[test]
fn a_robin_tag_without_represented_conductance_is_not_an_anchor() {
    let mesh = separated_cubes();
    // Positive as supplied, but every integrated face entry rounds to zero.
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "left bath",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(300.0).unwrap(),
        )
        .unwrap()
        .region(
            "vanished bath",
            |f| on_box_face(f.centroid[0], 4.0),
            ThermalBc::robin(f64::from_bits(1), 300.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let error = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            },
            SolveConfig {
                initial: InitialGuess::Uniform(300.0),
                ..Default::default()
            },
        )
        .unwrap_err()
    });
    assert!(matches!(error, ConductionError::SingularPureNeumann));
}

#[test]
fn isolated_transient_bodies_retain_their_own_caloric_initial_data() {
    let mesh = separated_cubes();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(1.0);
    let initial: Vec<_> = mesh
        .positions()
        .iter()
        .map(|p| if p[0] < 2.0 { 300.0 } else { 340.0 })
        .collect();
    let config = TransientConfig::backward_euler(0.1, SolveConfig::default().linear).unwrap();
    let result = with_cx(|cx| {
        march(
            cx,
            TransientProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
                capacity: VolumetricHeatCapacity::declared(2.0).unwrap(),
            },
            &config,
            &initial,
            4,
        )
        .unwrap()
    });
    // Each 1 m^3 body receives 0.4 J and has 2 J/K heat capacity.
    for (&before, &after) in initial.iter().zip(&result.temperature) {
        assert!((after - before - 0.2).abs() < 1e-9);
    }
}
