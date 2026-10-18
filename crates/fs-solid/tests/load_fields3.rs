//! Physical field-mass integration and native equilibrium, not centroid lumping.
#[path = "tet10/support.rs"]
mod support;
use support::{block, with_cx};
use fs_blake3::ContentHash;
use fs_solid::linear3::{TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState};
use fs_solid::loads3::{LoadError3, SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::fields::{SurfaceField3, VolumeField3};

fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }
fn component(values: &[(usize, f64)], dof: usize) -> f64 {
    values.binary_search_by_key(&dof, |v| v.0).map_or(0.0, |i| values[i].1)
}

#[test]
fn g1_pressure_gradient_keeps_its_force_and_bending_moment() {
    for n in [1, 2, 4] {
        let mesh = block([n; 3], [1.0; 3]);
        with_cx(|cx| {
            let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 2000, cx).unwrap();
            let right: Vec<_> = g.faces().iter().enumerate()
                .filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
            let values: Vec<_> = right.iter().map(|&i| g.faces()[i].nodes.map(|node| {
                let p = mesh.nodes[node]; 10.0 + 6.0*p[1] + 12.0*p[2]
            })).collect();
            let case = g.assemble_p1_fields(&[SurfaceField3::Pressure { faces: &right, nodal_pascals: &values }],
                VolumeField3::None, &[], cx).unwrap();
            let (force, moment) = case.resultant([0.0; 3], cx).unwrap();
            for (a, b) in force.into_iter().zip([-19.0, 0.0, 0.0]) { close(a, b, 1e-12); }
            for (a, b) in moment.into_iter().zip([0.0, -10.5, 10.0]) { close(a, b, 1e-12); }
        });
    }
}

#[test]
fn g0_surface_fields_use_consistent_not_equal_nodal_shares() {
    let mesh = block([1; 3], [1.0; 3]);
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 100, cx).unwrap();
        let face = g.faces()[0];
        let values = [[[1.0, 2.0, -1.0], [7.0, -4.0, 3.0], [-2.0, 8.0, 5.0]]];
        let case = g.assemble_p1_fields(&[SurfaceField3::Traction { faces: &[0], nodal_pascals: &values }],
            VolumeField3::None, &[], cx).unwrap();
        for i in 0..3 { for c in 0..3 {
            let sum: f64 = values[0].iter().map(|v| v[c]).sum();
            close(component(case.nodal_forces_n(), 3*face.nodes[i]+c),
                face.area_m2*(sum+values[0][i][c])/12.0, 1e-14);
        } }
        assert!((component(case.nodal_forces_n(), 3*face.nodes[0])
            - component(case.nodal_forces_n(), 3*face.nodes[1])).abs() > 0.1);
        let pressures = [[1.0, 3.0, 9.0]];
        let equivalent = [pressures[0].map(|p| face.outward_normal.map(|n| -p*n))];
        let pressure = g.assemble_p1_fields(&[SurfaceField3::Pressure { faces: &[0], nodal_pascals: &pressures }],
            VolumeField3::None, &[], cx).unwrap();
        let traction = g.assemble_p1_fields(&[SurfaceField3::Traction { faces: &[0], nodal_pascals: &equivalent }],
            VolumeField3::None, &[], cx).unwrap();
        assert_eq!(pressure.nodal_forces_n(), traction.nodal_forces_n());
    });
}

#[test]
fn g1_body_fields_preserve_first_moments_and_element_discontinuities() {
    for n in [1, 2] {
        let mesh = block([n; 3], [1.0; 3]);
        let values: Vec<_> = mesh.nodes.iter().map(|p| [p[0]+2.0*p[1], 2.0*p[2], -3.0*p[0]]).collect();
        let element_values: Vec<_> = mesh.tets.iter().map(|tet| tet.map(|node| values[node])).collect();
        with_cx(|cx| {
            let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 1000, cx).unwrap();
            let case = g.assemble_p1_fields(&[], VolumeField3::PerNode(&values), &[], cx).unwrap();
            let per_element = g.assemble_p1_fields(&[], VolumeField3::PerElement(&element_values), &[], cx).unwrap();
            assert_eq!(case.nodal_forces_n(), per_element.nodal_forces_n());
            let (force, moment) = case.resultant([0.0; 3], cx).unwrap();
            for (a, b) in force.into_iter().zip([1.5, 1.0, -1.5]) { close(a, b, 1e-12); }
            for (a, b) in moment.into_iter().zip([-17.0/12.0, 7.0/4.0, -5.0/12.0]) { close(a, b, 1e-12); }
            let discontinuous: Vec<_> = mesh.tets.iter().enumerate().map(|(i, _)| [[(i+1) as f64, 0.0, 0.0]; 4]).collect();
            let discontinuous = g.assemble_p1_fields(&[], VolumeField3::PerElement(&discontinuous), &[], cx).unwrap();
            let expected: f64 = g.element_volumes_m3().iter().enumerate().map(|(i, v)| (i+1) as f64*v).sum();
            close(discontinuous.resultant([0.0; 3], cx).unwrap().0[0], expected, 1e-12);
        });
    }
}

#[test]
fn g3_constant_field_limit_reaches_the_same_native_equilibrium() {
    let mesh = block([2; 3], [1.0; 3]);
    let mat = TetElasticMaterial::try_new(2.0, 1e6, 0.25, [7; 32]).unwrap();
    let thermal = TetThermalStrainState::try_new(300.0, 290.0, [0.0; 6], [7; 32], ContentHash([8; 32]), 0.1).unwrap();
    let fixed: Vec<_> = mesh.nodes.iter().enumerate().filter(|(_, p)| p[0] == 0.0)
        .flat_map(|(i, _)| [3*i, 3*i+1, 3*i+2]).collect();
    let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&mat), fixed_dofs: &fixed, budget: Default::default() };
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 1000, cx).unwrap();
        let right: Vec<_> = g.faces().iter().enumerate()
            .filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
        let values = vec![[[10.0, 2.0, 3.0]; 3]; right.len()];
        let body = vec![[0.0, 0.0, -4.0]; mesh.nodes.len()];
        let old = g.assemble(&[SurfaceLoad3::Traction { faces: &right, pascals: [10.0, 2.0, 3.0] }],
            VolumeLoad3::Uniform([0.0, 0.0, -4.0]), &[(0, 5.0)], cx).unwrap();
        let fields = g.assemble_p1_fields(&[SurfaceField3::Traction { faces: &right, nodal_pascals: &values }],
            VolumeField3::PerNode(&body), &[(0, 5.0)], cx).unwrap();
        for dof in 0..3*mesh.nodes.len() {
            close(component(old.nodal_forces_n(), dof), component(fields.nodal_forces_n(), dof), 1e-13);
        }
        let a = old.solve(&problem, TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([9; 32]), 0.01, cx).unwrap();
        let b = fields.solve(&problem, TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([9; 32]), 0.01, cx).unwrap();
        for (u, v) in a.displacement.displacement_m().iter().flatten().zip(b.displacement.displacement_m().iter().flatten()) { close(*u, *v, 1e-11); }
        assert!(b.free_force_residual_n < 1e-6);
        let force = fields.resultant([0.0; 3], cx).unwrap().0;
        for c in 0..3 { close(b.support_reactions_n.iter().map(|r| r[c]).sum::<f64>() + force[c], 0.0, 1e-6); }
    });
}

#[test]
fn g0_bad_field_shapes_values_and_partition_refuse_without_altering_geometry() {
    let mesh = block([1; 3], [1.0; 3]);
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 100, cx).unwrap();
        let faces = g.faces().to_vec();
        for (ids, values) in [(vec![], vec![]), (vec![0], vec![]), (vec![0, 0], vec![[1.0; 3]; 2]),
            (vec![g.faces().len()], vec![[1.0; 3]]), (vec![0], vec![[f64::NAN; 3]])] {
            assert!(g.assemble_p1_fields(&[SurfaceField3::Pressure { faces: &ids, nodal_pascals: &values }], VolumeField3::None, &[], cx).is_err());
        }
        assert!(g.assemble_p1_fields(&[], VolumeField3::PerNode(&[]), &[], cx).is_err());
        assert!(g.assemble_p1_fields(&[], VolumeField3::PerElement(&[]), &[], cx).is_err());
        let bad = vec![[f64::INFINITY; 3]; mesh.nodes.len()];
        assert!(g.assemble_p1_fields(&[], VolumeField3::PerNode(&bad), &[], cx).is_err());
        assert!(g.assemble_p1_fields(&[], VolumeField3::None, &[(0, 1.0), (0, 2.0)], cx).is_err());
        assert_eq!(g.faces(), faces);
    });
}

#[test]
fn g4_cancelled_fields_publish_nothing_and_preserve_reuse() {
    let mesh = block([1; 3], [1.0; 3]);
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 100, cx).unwrap();
        let values = vec![[1.0; 3]; mesh.nodes.len()];
        let gate = fs_exec::CancelGate::new_clock_free(); gate.request();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let stopped = fs_exec::Cx::new(&gate, arena, fs_exec::StreamKey {
                seed: 58, kernel_id: 58, tile: 0, iteration: 0,
            }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
            assert!(matches!(g.assemble_p1_fields(&[], VolumeField3::PerNode(&values), &[], &stopped), Err(LoadError3::Cancelled)));
        });
        let case = g.assemble_p1_fields(&[], VolumeField3::PerNode(&values), &[(0, 2.0)], cx).unwrap();
        let force = case.resultant([0.0; 3], cx).unwrap().0;
        for (a, b) in force.into_iter().zip([3.0, 1.0, 1.0]) { close(a, b, 1e-13); }
    });
}
