//! Heterogeneous per-element conductivity (bead s93ej.2).
//!
//! One-material assignment must reproduce the uniform path. A two-layer
//! slab with k=1 then k=2 has the series-resistance interface
//! temperature 1/3 when T(0)=1 and T(2)=0.

mod support;

use fs_adjoint::verify_gradient;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::assemble::{DofMap, residual};
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::{
    ConductivityModel, ConductivityTable, ElementMaterials, MaterialId, MaterialTable,
};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    CONDUCTIVITY_DIMS, ConductionError, ConductionProblem, InitialGuess, InterfaceFacePair,
    InterfaceResistance, InterfaceSurface, LineSearch, LinearConfig, Nonlinearity, ProvenanceClass,
    SolveConfig, StopRule, TEMPERATURE_DIMS, ThermalInterfaces,
    assemble_jacobian_with_element_materials, assemble_operator_with_element_materials,
    element_heat_flux_assigned, solve, solve_with_element_materials,
    solve_with_element_materials_and_interfaces, solve_with_interfaces,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_mesh::{
    RegionId, RegionKind, RegionSpec, UnverifiedPlc, VolumetricPolicy, box_triangles, box_vertices,
    volumetricize,
};
use fs_rep_mesh::TetComplex;
use std::collections::BTreeMap;
use support::with_cx;

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 4,
        },
        stop: StopRule {
            residual_rtol: 1e-12,
            residual_atol: 1e-14,
            step_atol: 1e-12,
            max_iterations: 8,
        },
        linear: LinearConfig {
            tolerance: 1e-14,
            max_iterations: 200,
            restart: 40,
        },
        initial: InitialGuess::DirichletMean,
    }
}

#[test]
fn empty_table_and_unknown_id_refuse() {
    assert!(matches!(
        MaterialTable::new(Vec::new()),
        Err(ConductionError::MaterialAssignment { .. })
    ));
    let k = ConductivityModel::isotropic_declared(1.0).expect("k");
    assert!(matches!(
        MaterialTable::new([(MaterialId(1), k.clone()), (MaterialId(1), k.clone()),]),
        Err(ConductionError::MaterialAssignment { .. })
    ));
    let table = MaterialTable::new([(MaterialId(1), k)]).expect("table");
    let err = ElementMaterials::new(table, vec![MaterialId(2)]).expect_err("unknown");
    assert!(matches!(err, ConductionError::MaterialAssignment { .. }));
}

#[test]
fn one_material_assignment_matches_the_uniform_path() {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(20.0).expect("k");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(340.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
    };
    let uniform = with_cx(|cx| solve(cx, problem, config()).expect("uniform"));
    let table = MaterialTable::new([(MaterialId(1), material.clone())]).expect("table");
    let assigned =
        ElementMaterials::new(table, vec![MaterialId(1); mesh.element_count()]).expect("assign");
    let hetero = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("hetero")
    });
    let worst = uniform
        .temperature
        .iter()
        .zip(&hetero.temperature)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-12,
        "uniform vs one-material assignment differed by {worst:e}"
    );
}

#[test]
fn two_layer_slab_hits_the_series_interface_temperature() {
    let length = 2.0;
    let (complex, positions) = box_grid([8, 2, 2], [length, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let right = ConductivityModel::isotropic_declared(2.0).expect("k2");
    let fallback = left.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], length),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let ids: Vec<MaterialId> = (0..mesh.element_count())
        .map(|e| {
            let tet = mesh.complex().tets[e];
            let cx = tet
                .iter()
                .map(|&v| mesh.positions()[v as usize][0])
                .sum::<f64>()
                / 4.0;
            if cx < 1.0 {
                MaterialId(1)
            } else {
                MaterialId(2)
            }
        })
        .collect();
    let assigned = ElementMaterials::new(table, ids).expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    // Equal thicknesses, k1=1, k2=2, T(0)=1, T(2)=0. Heat flows +x:
    // q = 1 / (1/1 + 1/2) = 2/3, so T(1) = 1 − q/k1 = 1/3.
    let mut worst = 0.0f64;
    let mut counted = 0usize;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if (p[0] - 1.0).abs() < 1e-9 {
            worst = worst.max((solution.temperature[v] - 1.0 / 3.0).abs());
            counted += 1;
        }
    }
    assert!(counted > 0, "no interface nodes at x=1");
    assert!(
        worst < 2e-3,
        "interface temperature off by {worst:e} from 1/3"
    );
    let flux = element_heat_flux_assigned(&mesh, &assigned, &solution.temperature).expect("flux");
    let want = 2.0 / 3.0;
    let worst_qx = flux
        .iter()
        .map(|q| (q[0] - want).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst_qx < 0.05,
        "recovered q_x off by {worst_qx} from series {want}"
    );
}

fn two_layer_ids(mesh: &ConductionMesh) -> Vec<u32> {
    (0..mesh.element_count())
        .map(|e| {
            let tet = mesh.complex().tets[e];
            let cx = tet
                .iter()
                .map(|&v| mesh.positions()[v as usize][0])
                .sum::<f64>()
                / 4.0;
            if cx < 1.0 { 1 } else { 2 }
        })
        .collect()
}

#[test]
fn region_map_refuses_an_unmapped_label_and_accepts_a_complete_map() {
    let table = MaterialTable::new([
        (
            MaterialId(10),
            ConductivityModel::isotropic_declared(1.0).expect("k1"),
        ),
        (
            MaterialId(20),
            ConductivityModel::isotropic_declared(2.0).expect("k2"),
        ),
    ])
    .expect("table");
    let mut map = BTreeMap::new();
    map.insert(1, MaterialId(10));
    let err =
        ElementMaterials::from_region_ids(table.clone(), &[1, 2], &map).expect_err("unmapped");
    assert!(matches!(err, ConductionError::MaterialAssignment { .. }));
    map.insert(2, MaterialId(20));
    let assigned = ElementMaterials::from_region_ids(table, &[1, 2, 1], &map).expect("map");
    assert_eq!(
        assigned.of_element(),
        &[MaterialId(10), MaterialId(20), MaterialId(10)]
    );
}

#[test]
fn assignment_length_mismatch_refuses_at_bind() {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(1.0).expect("k");
    let table = MaterialTable::new([(MaterialId(1), material.clone())]).expect("table");
    let assigned = ElementMaterials::new(table, vec![MaterialId(1)]).expect("one");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
    };
    let err = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect_err("len")
    });
    assert!(
        matches!(err, ConductionError::FieldLength { .. }),
        "length mismatch as {err:?}"
    );
}

#[test]
fn temperature_dependent_layer_refuses_outside_its_span() {
    let (complex, positions) = box_grid([4, 1, 1], [2.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let right = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(0.0, 2.0), (1.0, 2.0)]).expect("curve"),
    );
    let fallback = left.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(400.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let regions = two_layer_ids(&mesh);
    let mut map = BTreeMap::new();
    map.insert(1, MaterialId(1));
    map.insert(2, MaterialId(2));
    let assigned = ElementMaterials::from_region_ids(table, &regions, &map).expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let err = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect_err("span")
    });
    assert!(
        matches!(err, ConductionError::OutsideTemperatureSpan { .. }),
        "out-of-span as {err:?}"
    );
}

#[test]
fn anisotropic_layers_keep_the_x_series_interface() {
    let length = 2.0;
    let (complex, positions) = box_grid([8, 2, 2], [length, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left =
        ConductivityModel::constant_tensor([[1.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
            .expect("k1");
    let right =
        ConductivityModel::constant_tensor([[2.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
            .expect("k2");
    let fallback = left.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], length),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let regions = two_layer_ids(&mesh);
    let mut map = BTreeMap::new();
    map.insert(1, MaterialId(1));
    map.insert(2, MaterialId(2));
    let assigned = ElementMaterials::from_region_ids(table, &regions, &map).expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    let mut worst = 0.0f64;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if (p[0] - 1.0).abs() < 1e-9 {
            worst = worst.max((solution.temperature[v] - 1.0 / 3.0).abs());
        }
    }
    assert!(
        worst < 2e-3,
        "anisotropic x-series interface off by {worst:e}"
    );
}

fn weld(parts: &[Vec<[f64; 3]>]) -> (Vec<[f64; 3]>, Vec<Vec<u32>>) {
    let mut verts = Vec::new();
    let mut index = BTreeMap::new();
    let mut remaps = Vec::new();
    for part in parts {
        let mut remap = Vec::with_capacity(part.len());
        for p in part {
            let key = [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
            let id = if let Some(&existing) = index.get(&key) {
                existing
            } else {
                let id = u32::try_from(verts.len()).expect("vertex count");
                verts.push(*p);
                index.insert(key, id);
                id
            };
            remap.push(id);
        }
        remaps.push(remap);
    }
    (verts, remaps)
}

fn remap_tris(tris: &[[u32; 3]], remap: &[u32]) -> Vec<[u32; 3]> {
    tris.iter()
        .map(|t| {
            [
                remap[t[0] as usize],
                remap[t[1] as usize],
                remap[t[2] as usize],
            ]
        })
        .collect()
}

#[test]
fn labeled_adjacent_volumes_solve_as_two_materials() {
    let left = box_vertices(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    let right = box_vertices(1.0, 2.0, 0.0, 1.0, 0.0, 1.0);
    let (verts, remaps) = weld(&[left, right]);
    let regions = vec![
        RegionSpec {
            id: RegionId(1),
            kind: RegionKind::Solid,
            seed: [0.5, 0.5, 0.5],
            triangles: remap_tris(&box_triangles(0), &remaps[0]),
        },
        RegionSpec {
            id: RegionId(2),
            kind: RegionKind::Solid,
            seed: [1.5, 0.5, 0.5],
            triangles: remap_tris(&box_triangles(0), &remaps[1]),
        },
    ];
    let audited = with_cx(|cx| {
        volumetricize(
            UnverifiedPlc::new(verts, regions),
            VolumetricPolicy::fixture_default("m"),
            cx,
        )
        .expect("volume")
    });
    let labeled = audited.labeled();
    let complex = TetComplex::from_tets(labeled.positions().len(), labeled.tets().to_vec());
    let mesh = ConductionMesh::new(complex, labeled.positions().to_vec()).expect("mesh");
    let k1 = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let k2 = ConductivityModel::isotropic_declared(2.0).expect("k2");
    let fallback = k1.clone();
    let table = MaterialTable::new([(MaterialId(1), k1), (MaterialId(2), k2)]).expect("table");
    let region_ids: Vec<u32> = labeled.region_of_tet().iter().map(|r| r.0).collect();
    let mut map = BTreeMap::new();
    map.insert(1, MaterialId(1));
    map.insert(2, MaterialId(2));
    let assigned = ElementMaterials::from_region_ids(table, &region_ids, &map).expect("assign");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    let mut worst = 0.0f64;
    let mut counted = 0usize;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if (p[0] - 1.0).abs() < 1e-9 {
            worst = worst.max((solution.temperature[v] - 1.0 / 3.0).abs());
            counted += 1;
        }
    }
    assert!(counted > 0, "no interface nodes from the labeled mesh");
    assert!(
        worst < 5e-2,
        "labeled-mesh two-material interface off by {worst:e}"
    );
}

#[test]
fn assigned_linear_adjoint_refuses_k_of_t_and_matches_finite_differences() {
    let (complex, positions) = box_grid([4, 2, 2], [2.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let right = ConductivityModel::isotropic_declared(2.0).expect("k2");
    let k_of_t = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.0), (450.0, 2.0)]).expect("curve"),
    );
    let fallback = left.clone();
    let source = ScalarField::Uniform(4.0e3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(290.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table_bad =
        MaterialTable::new([(MaterialId(1), left.clone()), (MaterialId(2), k_of_t)]).expect("bad");
    let regions = two_layer_ids(&mesh);
    let mut map = BTreeMap::new();
    map.insert(1, MaterialId(1));
    map.insert(2, MaterialId(2));
    let assigned_bad =
        ElementMaterials::from_region_ids(table_bad, &regions, &map).expect("assign-bad");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let refused = ConductivityDesign::new_with_element_materials(
        problem,
        &assigned_bad,
        LinearConfig {
            tolerance: 1e-14,
            max_iterations: 20_000,
            restart: 40,
        },
    );
    match refused {
        Err(ConductionError::Conductivity { .. }) => {}
        Err(other) => panic!("expected Conductivity refusal for k(T), got {other}"),
        Ok(_) => panic!("k(T) assigned model must refuse the linear IFT hook"),
    }

    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let assigned = ElementMaterials::from_region_ids(table, &regions, &map).expect("assign");
    let design = ConductivityDesign::new_with_element_materials(
        problem,
        &assigned,
        LinearConfig {
            tolerance: 1e-14,
            max_iterations: 20_000,
            restart: 40,
        },
    )
    .expect("design");
    let np = design.parameter_count();
    let nf = design.dofs().n();
    let rho: Vec<f64> = (0..np)
        .map(|e| 0.75 + 0.5 * ((e % 7) as f64) / 7.0)
        .collect();
    let weights = vec![1.0 / nf as f64; nf];
    let (gradient, report) = with_cx(|cx| design.gradient(cx, &rho, &weights).expect("grad"));
    assert!(report.converged);
    let objective =
        |p: &[f64]| -> f64 { with_cx(|cx| design.objective(cx, p, &weights).expect("J")) };
    let dirs = vec![
        {
            let mut d = vec![0.0; np];
            d[0] = 1.0;
            d
        },
        {
            let mut d = vec![0.0; np];
            d[np / 2] = 1.0;
            d
        },
        (0..np).map(|i| (i as f64 + 1.0) / np as f64).collect(),
    ];
    let verdict = verify_gradient(&objective, &rho, &gradient, &dirs, 1e-6, 8e-6);
    assert!(
        verdict.pass,
        "assigned IFT vs FD failed: max_rel_err={:e} informative={} pairs={:?}",
        verdict.max_rel_err, verdict.informative_directions, verdict.pairs
    );
    assert_eq!(verdict.informative_directions, dirs.len());
}

fn three_layer_ids(mesh: &ConductionMesh) -> Vec<MaterialId> {
    (0..mesh.element_count())
        .map(|e| {
            let tet = mesh.complex().tets[e];
            let cx = tet
                .iter()
                .map(|&v| mesh.positions()[v as usize][0])
                .sum::<f64>()
                / 4.0;
            if cx < 1.0 {
                MaterialId(1)
            } else if cx < 2.0 {
                MaterialId(2)
            } else {
                MaterialId(3)
            }
        })
        .collect()
}

fn newton_config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::Newton {
            line_search: LineSearch::default(),
        },
        stop: StopRule {
            residual_rtol: 1e-12,
            residual_atol: 1e-14,
            step_atol: 1e-12,
            max_iterations: 16,
        },
        linear: LinearConfig {
            tolerance: 1e-14,
            max_iterations: 400,
            restart: 40,
        },
        initial: InitialGuess::DirichletMean,
    }
}

#[test]
fn three_layer_slab_hits_the_series_interface_temperatures() {
    let length = 3.0;
    let (complex, positions) = box_grid([9, 2, 2], [length, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let k1 = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let k2 = ConductivityModel::isotropic_declared(2.0).expect("k2");
    let k3 = ConductivityModel::isotropic_declared(4.0).expect("k3");
    let fallback = k1.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], length),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([
        (MaterialId(1), k1),
        (MaterialId(2), k2),
        (MaterialId(3), k3),
    ])
    .expect("table");
    let assigned = ElementMaterials::new(table, three_layer_ids(&mesh)).expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    // Equal thicknesses, k = 1, 2, 4, T(0)=1, T(3)=0.
    // R = 1 + 1/2 + 1/4 = 7/4, q = 4/7, T(1)=3/7, T(2)=1/7.
    let mut worst_1 = 0.0f64;
    let mut worst_2 = 0.0f64;
    let mut counted_1 = 0usize;
    let mut counted_2 = 0usize;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if (p[0] - 1.0).abs() < 1e-9 {
            worst_1 = worst_1.max((solution.temperature[v] - 3.0 / 7.0).abs());
            counted_1 += 1;
        }
        if (p[0] - 2.0).abs() < 1e-9 {
            worst_2 = worst_2.max((solution.temperature[v] - 1.0 / 7.0).abs());
            counted_2 += 1;
        }
    }
    assert!(counted_1 > 0 && counted_2 > 0, "missing series interfaces");
    assert!(
        worst_1 < 3e-3 && worst_2 < 3e-3,
        "three-layer interfaces off by {worst_1:e} and {worst_2:e}"
    );
    let flux = element_heat_flux_assigned(&mesh, &assigned, &solution.temperature).expect("flux");
    let want = 4.0 / 7.0;
    let worst_qx = flux
        .iter()
        .map(|q| (q[0] - want).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst_qx < 0.06,
        "recovered q_x off by {worst_qx} from series {want}"
    );
    // Source-free Dirichlet–Dirichlet: net Dirichlet inflow cancels, so
    // `relative_closure` is vacuous. The assembled residual must still
    // close in Watts.
    assert!(
        solution.report.energy.closure_w.abs() < 1e-12,
        "three-layer energy did not close: {:?}",
        solution.report.energy
    );
}

#[test]
fn assigned_temperature_dependent_newton_solves_inside_span() {
    let (complex, positions) = box_grid([6, 2, 2], [2.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.0), (450.0, 1.0)]).expect("left"),
    );
    let right = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.5), (450.0, 2.5)]).expect("right"),
    );
    let fallback = left.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(340.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let regions = two_layer_ids(&mesh);
    let mut map = BTreeMap::new();
    map.insert(1, MaterialId(1));
    map.insert(2, MaterialId(2));
    let assigned = ElementMaterials::from_region_ids(table, &regions, &map).expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, newton_config()).expect("newton")
    });
    assert!(
        solution.report.energy.closure_w.abs() < 1e-8,
        "assigned k(T) energy did not close: {:?}",
        solution.report.energy
    );
    for &t in &solution.temperature {
        assert!(
            (300.0..=340.0).contains(&t),
            "assigned k(T) temperature {t} left the Dirichlet interval"
        );
    }
    let flux = element_heat_flux_assigned(&mesh, &assigned, &solution.temperature).expect("flux");
    assert!(
        flux.iter().all(|q| q[0] > 0.0),
        "heat must flow +x from the hot face"
    );
}

#[test]
fn assigned_newton_jacobian_matches_central_differences() {
    let (complex, positions) = box_grid([3, 2, 2], [2.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let right = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.5), (450.0, 4.5)]).expect("curve"),
    );
    let fallback = left.clone();
    let source = ScalarField::Uniform(2.0e3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(320.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let regions = two_layer_ids(&mesh);
    let mut map = BTreeMap::new();
    map.insert(1, MaterialId(1));
    map.insert(2, MaterialId(2));
    let assigned = ElementMaterials::from_region_ids(table, &regions, &map).expect("assign");
    let dofs = DofMap::new(&boundary, mesh.vertex_count()).expect("dofs");
    let free: Vec<f64> = (0..dofs.n())
        .map(|i| 305.0 + 7.0 * ((i % 5) as f64) - 3.0 * ((i % 3) as f64))
        .collect();
    let residual_at = |free: &[f64]| -> Vec<f64> {
        let full = dofs.scatter(free);
        with_cx(|cx| {
            let system = assemble_operator_with_element_materials(
                cx, &mesh, &boundary, &fallback, &source, &full, &assigned,
            )
            .expect("asm");
            residual(&system, &dofs, &full)
        })
    };
    let base = residual_at(&free);
    assert!(
        base.iter().any(|v| v.abs() > 1e-6),
        "the probe iterate must produce a non-trivial residual"
    );
    let jacobian = with_cx(|cx| {
        let full = dofs.scatter(&free);
        assemble_jacobian_with_element_materials(cx, &mesh, &boundary, &fallback, &full, &assigned)
            .expect("jacobian")
    });
    let (j_ff, _) = fs_conduction::assemble::reduce_matrix_and_lift(&jacobian, &dofs);
    let eps = 1e-4;
    let mut worst = 0.0f64;
    let mut scale = 0.0f64;
    for c in 0..dofs.n() {
        let mut plus = free.clone();
        let mut minus = free.clone();
        plus[c] += eps;
        minus[c] -= eps;
        let rp = residual_at(&plus);
        let rm = residual_at(&minus);
        for row in 0..dofs.n() {
            let fd = (rp[row] - rm[row]) / (2.0 * eps);
            let analytic = j_ff.get(row, c);
            scale = scale.max(fd.abs()).max(analytic.abs());
            worst = worst.max((fd - analytic).abs());
        }
    }
    assert!(
        worst <= 1e-5 * scale,
        "assigned Jacobian vs FD off by {worst:e} against scale {scale:e}"
    );
}

#[test]
fn assigned_materials_do_not_erase_an_undeclared_interface() {
    // Two tets share a geometrically coincident triangle on distinct
    // vertices. Region materials cannot stand in for the contact card.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(positions.len(), vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        positions,
    )
    .expect("two-tet contact mesh");
    let left = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let right = ConductivityModel::isotropic_declared(2.0).expect("k2");
    let fallback = left.clone();
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let assigned =
        ElementMaterials::new(table, vec![MaterialId(1), MaterialId(2)]).expect("assign");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| f.centroid[2] < -0.2 && f.centroid[0] > 0.2 && f.centroid[1] > 0.2,
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| f.centroid[2] > 0.2 && f.centroid[0] > 0.2 && f.centroid[1] > 0.2,
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let err = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect_err("interface")
    });
    assert!(
        matches!(err, ConductionError::Interface { .. }),
        "undeclared coincident faces must refuse even with a material map: {err}"
    );
}

fn two_slab_contact_mesh(n: usize) -> (ConductionMesh, usize) {
    let mut tets = Vec::new();
    let mut positions = Vec::new();
    let mut vertices_per_slab = 0usize;
    for slab in 0..2 {
        let (complex, slab_positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
        if slab == 0 {
            vertices_per_slab = slab_positions.len();
        }
        let offset = u32::try_from(positions.len()).expect("vertex count");
        tets.extend(
            complex
                .tets
                .into_iter()
                .map(|tet| tet.map(|vertex| vertex + offset)),
        );
        positions.extend(
            slab_positions
                .into_iter()
                .map(|[x, y, z]| [x + slab as f64, y, z]),
        );
    }
    (
        ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions)
            .expect("two-slab mesh"),
        vertices_per_slab,
    )
}

fn contact_card(r_m2k_per_w: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: r_m2k_per_w,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "s93ej.2 assigned-contact fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("contact claim");
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: MaterialStateId {
                chemistry: "solid-a".to_string(),
                phase: "solid".to_string(),
                process: "as-fixtured".to_string(),
                revision: 0,
            },
            texture_frame: "interface-normal-plus-x".to_string(),
        },
        SurfaceSpec {
            material: MaterialStateId {
                chemistry: "solid-b".to_string(),
                phase: "solid".to_string(),
                process: "as-fixtured".to_string(),
                revision: 0,
            },
            texture_frame: "interface-normal-minus-x".to_string(),
        },
        SystemContext {
            medium: "dry".to_string(),
            third_body: Some("declared-contact-layer".to_string()),
            environment: "vacuum".to_string(),
            history: "unaged".to_string(),
        },
        claims,
        Vec::new(),
    )
    .expect("interface card")
}

fn oriented_pairs(mesh: &ConductionMesh) -> Vec<InterfaceFacePair> {
    ThermalInterfaces::coincident_face_pairs(mesh)
        .expect("coincident pairs")
        .into_iter()
        .map(|pair| {
            if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 {
                pair
            } else {
                InterfaceFacePair {
                    side_a: pair.side_b,
                    side_b: pair.side_a,
                }
            }
        })
        .collect()
}

fn bind_bondline(
    mesh: &ConductionMesh,
    boundary: &fs_conduction::ThermalBoundary,
) -> ThermalInterfaces {
    let resistance = InterfaceResistance::from_card(
        "bondline",
        &contact_card(0.1),
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("resistance");
    let surface =
        InterfaceSurface::new("bondline", oriented_pairs(mesh), resistance).expect("surface");
    ThermalInterfaces::new(mesh, boundary, vec![surface]).expect("interfaces")
}

fn slab_ids(mesh: &ConductionMesh) -> Vec<MaterialId> {
    (0..mesh.element_count())
        .map(|e| {
            let tet = mesh.complex().tets[e];
            let cx = tet
                .iter()
                .map(|&v| mesh.positions()[v as usize][0])
                .sum::<f64>()
                / 4.0;
            if cx < 1.0 {
                MaterialId(1)
            } else {
                MaterialId(2)
            }
        })
        .collect()
}

#[test]
fn assigned_one_material_contact_matches_the_uniform_contact_path() {
    let (mesh, _) = two_slab_contact_mesh(3);
    let k = ConductivityModel::isotropic_declared(10.0).expect("k");
    let fallback = ConductivityModel::isotropic_declared(1.0).expect("wrong");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(330.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let interfaces = bind_bondline(&mesh, &boundary);
    let table = MaterialTable::new([(MaterialId(1), k.clone())]).expect("table");
    let assigned =
        ElementMaterials::new(table, vec![MaterialId(1); mesh.element_count()]).expect("assign");
    let assigned_problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let uniform_problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &k,
        source: &source,
    };
    let assigned_sol = with_cx(|cx| {
        solve_with_element_materials_and_interfaces(
            cx,
            assigned_problem,
            &assigned,
            &interfaces,
            config(),
        )
        .expect("assigned")
    });
    let uniform_sol = with_cx(|cx| {
        solve_with_interfaces(cx, uniform_problem, &interfaces, config()).expect("uniform")
    });
    let worst = assigned_sol
        .temperature
        .iter()
        .zip(&uniform_sol.temperature)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-12,
        "assigned one-material contact differed from uniform by {worst:e}"
    );
    assert_eq!(assigned_sol.report.interface_fluxes.len(), 1);
    assert_eq!(
        assigned_sol.report.interface_fluxes[0].heat_rate_a_to_b_w,
        uniform_sol.report.interface_fluxes[0].heat_rate_a_to_b_w
    );
}

#[test]
fn assigned_two_material_contact_hits_the_series_jump() {
    let (mesh, left_n) = two_slab_contact_mesh(3);
    let left = ConductivityModel::isotropic_declared(10.0).expect("k1");
    let right = ConductivityModel::isotropic_declared(20.0).expect("k2");
    let fallback = ConductivityModel::isotropic_declared(1.0).expect("wrong");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(330.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let interfaces = bind_bondline(&mesh, &boundary);
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let assigned = ElementMaterials::new(table, slab_ids(&mesh)).expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials_and_interfaces(cx, problem, &assigned, &interfaces, config())
            .expect("solve")
    });
    // R = 1/10 + 0.1 + 1/20 = 0.25 K/W, Q = 30/0.25 = 120 W.
    // Left: T = 330 − 12 x. Right: T = 312 − 6 x.
    let mut worst = 0.0f64;
    for (v, &p) in mesh.positions().iter().enumerate() {
        let exact = if v < left_n {
            330.0 - 12.0 * p[0]
        } else {
            312.0 - 6.0 * p[0]
        };
        worst = worst.max((solution.temperature[v] - exact).abs());
    }
    assert!(
        worst < 1e-6,
        "two-material contact profile off by {worst:e}"
    );
    let flux = solution
        .report
        .interface_fluxes
        .first()
        .expect("interface flux");
    assert!((flux.heat_rate_a_to_b_w - 120.0).abs() < 1e-6);
    assert!((flux.mean_jump_k - 12.0).abs() < 1e-6);
}

fn conductivity_claims(knots: Vec<(f64, f64)>) -> ClaimSet {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".to_string(),
                abscissa_dims: TEMPERATURE_DIMS,
                knots,
                dims: CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(),
            provenance: Provenance {
                source: "s93ej.2 assigned-receipt fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("claim");
    claims
}

#[test]
fn assigned_matdb_receipts_travel_with_the_solve() {
    let claims_a = conductivity_claims(vec![(250.0, 10.0), (500.0, 10.0)]);
    let claims_b = conductivity_claims(vec![(250.0, 20.0), (350.0, 20.0), (500.0, 20.0)]);
    let table_a = ConductivityTable::from_claims(
        &claims_a,
        "thermal_conductivity",
        &[280.0, 400.0],
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("table-a");
    let table_b = ConductivityTable::from_claims(
        &claims_b,
        "thermal_conductivity",
        &[280.0, 340.0, 400.0],
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("table-b");
    assert_eq!(table_a.receipts().len(), 2);
    assert_eq!(table_b.receipts().len(), 3);
    let left = ConductivityModel::isotropic(table_a);
    let right = ConductivityModel::isotropic(table_b);
    let fallback = ConductivityModel::isotropic_declared(1.0).expect("declared fallback");
    let (complex, positions) = box_grid([4, 2, 2], [2.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(340.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let materials =
        MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    assert_eq!(materials.receipts().len(), 5);
    let assigned = ElementMaterials::from_region_ids(
        materials,
        &two_layer_ids(&mesh),
        &BTreeMap::from([(1, MaterialId(1)), (2, MaterialId(2))]),
    )
    .expect("assign");
    assert_eq!(assigned.receipts().len(), 5);
    assert_eq!(assigned.provenance(), ProvenanceClass::MatdbReceipts);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    assert_eq!(solution.report.material_receipts, 5);
    assert_eq!(
        solution.report.material_provenance,
        ProvenanceClass::MatdbReceipts
    );
    // The unused fallback is declared and must not overwrite the assignment.
    assert_eq!(fallback.provenance(), ProvenanceClass::Declared);
}

#[test]
fn unused_table_entries_and_mixed_provenance_stay_honest() {
    let claims = conductivity_claims(vec![(250.0, 10.0), (500.0, 10.0)]);
    let sampled = ConductivityTable::from_claims(
        &claims,
        "thermal_conductivity",
        &[280.0, 400.0],
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("sampled");
    let matdb = ConductivityModel::isotropic(sampled);
    let declared = ConductivityModel::isotropic_declared(20.0).expect("declared");
    let table = MaterialTable::new([
        (MaterialId(1), matdb.clone()),
        (MaterialId(2), declared.clone()),
    ])
    .expect("table");
    assert_eq!(table.receipts().len(), 2);
    assert_eq!(table.provenance(), ProvenanceClass::Declared);

    let only_matdb = ElementMaterials::new(table.clone(), vec![MaterialId(1); 4]).expect("only");
    assert_eq!(only_matdb.receipts().len(), 2);
    assert_eq!(only_matdb.provenance(), ProvenanceClass::MatdbReceipts);

    let mixed = ElementMaterials::new(table, vec![MaterialId(1), MaterialId(2)]).expect("mixed");
    assert_eq!(mixed.receipts().len(), 2);
    assert_eq!(mixed.provenance(), ProvenanceClass::Declared);

    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(340.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    // Keep the unused declared model in the table; it must not change the report.
    let assigned = ElementMaterials::new(
        MaterialTable::new([(MaterialId(1), matdb), (MaterialId(2), declared)]).expect("unused"),
        vec![MaterialId(1); mesh.element_count()],
    )
    .expect("assign-unused");
    assert_eq!(assigned.receipts().len(), 2);
    assert_eq!(assigned.provenance(), ProvenanceClass::MatdbReceipts);
    let fallback = ConductivityModel::isotropic_declared(1.0).expect("fallback");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    assert_eq!(solution.report.material_receipts, 2);
    assert_eq!(
        solution.report.material_provenance,
        ProvenanceClass::MatdbReceipts
    );
}

fn rotated_kx(kx: f64, k_perp: f64) -> ConductivityModel {
    // Principal axis 0 = +y, 1 = −x, 2 = +z. Conductivity along x is table 1.
    ConductivityModel::orthotropic(
        [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [
            ConductivityTable::declared(k_perp).expect("k0"),
            ConductivityTable::declared(kx).expect("k1"),
            ConductivityTable::declared(k_perp).expect("k2"),
        ],
    )
    .expect("orthotropic")
}

#[test]
fn rotated_orthotropic_layers_keep_the_x_series_interface() {
    let length = 2.0;
    let (complex, positions) = box_grid([8, 2, 2], [length, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = rotated_kx(1.0, 10.0);
    let right = rotated_kx(2.0, 10.0);
    let fallback = left.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], length),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let assigned = ElementMaterials::from_region_ids(
        table,
        &two_layer_ids(&mesh),
        &BTreeMap::from([(1, MaterialId(1)), (2, MaterialId(2))]),
    )
    .expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    let mut worst = 0.0f64;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if (p[0] - 1.0).abs() < 1e-9 {
            worst = worst.max((solution.temperature[v] - 1.0 / 3.0).abs());
        }
    }
    assert!(
        worst < 2e-3,
        "rotated orthotropic x-series interface off by {worst:e}"
    );
}

fn two_layer_assigned(mesh: &ConductionMesh, k1: f64, k2: f64) -> ElementMaterials {
    let table = MaterialTable::new([
        (
            MaterialId(1),
            ConductivityModel::isotropic_declared(k1).expect("k1"),
        ),
        (
            MaterialId(2),
            ConductivityModel::isotropic_declared(k2).expect("k2"),
        ),
    ])
    .expect("table");
    ElementMaterials::from_region_ids(
        table,
        &two_layer_ids(mesh),
        &BTreeMap::from([(1, MaterialId(1)), (2, MaterialId(2))]),
    )
    .expect("assign")
}

fn mean_interface_t(mesh: &ConductionMesh, temperature: &[f64], x: f64) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if (p[0] - x).abs() < 1e-9 {
            sum += temperature[v];
            n += 1;
        }
    }
    assert!(n > 0, "no nodes at x={x}");
    sum / n as f64
}

fn mean_qx(mesh: &ConductionMesh, assigned: &ElementMaterials, temperature: &[f64]) -> f64 {
    let flux = element_heat_flux_assigned(mesh, assigned, temperature).expect("flux");
    flux.iter().map(|q| q[0]).sum::<f64>() / flux.len() as f64
}

#[test]
fn raising_the_right_layer_conductivity_raises_flux_and_drops_interface_t() {
    let fallback = ConductivityModel::isotropic_declared(1.0).expect("fallback");
    let mut last_q = f64::NEG_INFINITY;
    let mut last_t = f64::INFINITY;
    for k2 in [1.0, 2.0, 4.0] {
        let (complex, positions) = box_grid([6, 2, 2], [2.0, 1.0, 1.0]);
        let mesh = ConductionMesh::new(complex, positions).expect("mesh");
        let assigned = two_layer_assigned(&mesh, 1.0, k2);
        let source = ScalarField::Uniform(0.0);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "hot",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(1.0).expect("bc"),
            )
            .expect("hot")
            .region(
                "cold",
                |f| on_box_face(f.centroid[0], 2.0),
                ThermalBc::dirichlet(0.0).expect("bc"),
            )
            .expect("cold")
            .adiabatic_remainder()
            .finish()
            .expect("boundary");
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &fallback,
            source: &source,
        };
        let solution = with_cx(|cx| {
            solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
        });
        let q = mean_qx(&mesh, &assigned, &solution.temperature);
        let t = mean_interface_t(&mesh, &solution.temperature, 1.0);
        assert!(
            q > last_q + 1e-4,
            "q should rise with k2={k2}: {q} vs previous {last_q}"
        );
        assert!(
            t < last_t - 1e-4,
            "interface T should fall with k2={k2}: {t} vs previous {last_t}"
        );
        last_q = q;
        last_t = t;
    }
}

#[test]
fn two_ids_with_the_same_k_match_one_material() {
    let (complex, positions) = box_grid([4, 2, 2], [2.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let k = ConductivityModel::isotropic_declared(20.0).expect("k");
    let fallback = ConductivityModel::isotropic_declared(1.0).expect("wrong");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let one = ElementMaterials::new(
        MaterialTable::new([(MaterialId(1), k.clone())]).expect("one"),
        vec![MaterialId(1); mesh.element_count()],
    )
    .expect("one");
    let two = ElementMaterials::new(
        MaterialTable::new([(MaterialId(7), k.clone()), (MaterialId(9), k)]).expect("two"),
        two_layer_ids(&mesh)
            .into_iter()
            .map(|r| if r == 1 { MaterialId(7) } else { MaterialId(9) })
            .collect(),
    )
    .expect("two");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let a = with_cx(|cx| solve_with_element_materials(cx, problem, &one, config()).expect("one"));
    let b = with_cx(|cx| solve_with_element_materials(cx, problem, &two, config()).expect("two"));
    let worst = a
        .temperature
        .iter()
        .zip(&b.temperature)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-12,
        "identical-k two-id assignment differed by {worst:e}"
    );
}

#[test]
fn swapped_region_map_is_the_same_physics() {
    let (complex, positions) = box_grid([6, 2, 2], [2.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let k1 = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let k2 = ConductivityModel::isotropic_declared(2.0).expect("k2");
    let fallback = k1.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let regions = two_layer_ids(&mesh);
    let forward = ElementMaterials::from_region_ids(
        MaterialTable::new([(MaterialId(1), k1.clone()), (MaterialId(2), k2.clone())])
            .expect("fwd"),
        &regions,
        &BTreeMap::from([(1, MaterialId(1)), (2, MaterialId(2))]),
    )
    .expect("fwd");
    let swapped = ElementMaterials::from_region_ids(
        MaterialTable::new([(MaterialId(1), k2), (MaterialId(2), k1)]).expect("swp"),
        &regions,
        &BTreeMap::from([(1, MaterialId(2)), (2, MaterialId(1))]),
    )
    .expect("swp");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let a =
        with_cx(|cx| solve_with_element_materials(cx, problem, &forward, config()).expect("fwd"));
    let b =
        with_cx(|cx| solve_with_element_materials(cx, problem, &swapped, config()).expect("swp"));
    let worst = a
        .temperature
        .iter()
        .zip(&b.temperature)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-12,
        "swapped region map changed the field by {worst:e}"
    );
}

#[test]
fn two_layer_interface_tightens_under_refinement() {
    let fallback = ConductivityModel::isotropic_declared(1.0).expect("fallback");
    let mut last_err = f64::INFINITY;
    for n in [4usize, 8, 12] {
        let (complex, positions) = box_grid([n, 2, 2], [2.0, 1.0, 1.0]);
        let mesh = ConductionMesh::new(complex, positions).expect("mesh");
        let assigned = two_layer_assigned(&mesh, 1.0, 2.0);
        let source = ScalarField::Uniform(0.0);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "hot",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(1.0).expect("bc"),
            )
            .expect("hot")
            .region(
                "cold",
                |f| on_box_face(f.centroid[0], 2.0),
                ThermalBc::dirichlet(0.0).expect("bc"),
            )
            .expect("cold")
            .adiabatic_remainder()
            .finish()
            .expect("boundary");
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &fallback,
            source: &source,
        };
        let solution = with_cx(|cx| {
            solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
        });
        let err = (mean_interface_t(&mesh, &solution.temperature, 1.0) - 1.0 / 3.0).abs();
        assert!(
            err <= last_err + 1e-12,
            "interface error grew under refinement n={n}: {err:e} vs {last_err:e}"
        );
        last_err = err;
    }
    assert!(
        last_err < 2e-3,
        "finest two-layer interface still off by {last_err:e}"
    );
}
