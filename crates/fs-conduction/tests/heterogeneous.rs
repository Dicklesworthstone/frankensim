//! Heterogeneous per-element conductivity (bead s93ej.2).
//!
//! One-material assignment must reproduce the uniform path. A two-layer
//! slab with k=1 then k=2 has the series-resistance interface
//! temperature 1/3 when T(0)=1 and T(2)=0.

mod support;

use fs_adjoint::verify_gradient;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::{
    ConductivityModel, ConductivityTable, ElementMaterials, MaterialId, MaterialTable,
};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::{
    ConductionError, ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig,
    StopRule, element_heat_flux_assigned, solve, solve_with_element_materials,
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
    let table = MaterialTable::new([(
        MaterialId(1),
        ConductivityModel::isotropic_declared(1.0).expect("k"),
    )])
    .expect("table");
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
