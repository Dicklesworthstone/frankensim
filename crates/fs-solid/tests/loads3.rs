//! Physical G0/G1/G3 load, equilibrium and reaction checks on P1 tetrahedra.
#[path = "tet10/support.rs"]
mod support;
use support::{block, with_cx};
use fs_blake3::ContentHash;
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem,
    TetMaterialField, TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState};
use fs_solid::loads3::{LoadError3, SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};

fn material() -> TetElasticMaterial { TetElasticMaterial::try_new(2.0, 1e6, 0.25, [7; 32]).unwrap() }
fn thermal(strain: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(300.0, 290.0, [strain, strain, strain, 0.0, 0.0, 0.0],
        [7; 32], ContentHash([8; 32]), 0.1).unwrap()
}
fn rollers(nodes: &[[f64; 3]]) -> Vec<usize> {
    nodes.iter().enumerate().flat_map(|(i, p)| (0..3).filter_map(move |c| (p[c] == 0.0).then_some(3*i+c))).collect()
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }

#[test]
fn g0_exterior_faces_pressure_and_volume_loads_preserve_force_and_moment() {
    let mesh = block([2; 3], [2.0, 1.0, 1.0]);
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 1000, cx).unwrap();
        assert_eq!(g.faces().len(), 48);
        close(g.element_volumes_m3().iter().sum(), 2.0, 1e-13);
        let all: Vec<_> = (0..g.faces().len()).collect();
        let closed = g.assemble(&[SurfaceLoad3::Pressure { faces: &all, pascals: 12.0 }],
            VolumeLoad3::None, &[], cx).unwrap();
        let (force, moment) = closed.resultant([0.25, -4.0, 8.0], cx).unwrap();
        for value in force.into_iter().chain(moment) { close(value, 0.0, 1e-12); }
        let right: Vec<_> = g.faces().iter().enumerate()
            .filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
        let case = g.assemble(&[SurfaceLoad3::Traction { faces: &right, pascals: [0.0, 0.0, 3.0] }],
            VolumeLoad3::Uniform([0.0, -4.0, 0.0]), &[(0, 2.0)], cx).unwrap();
        let (force, moment) = case.resultant([0.0; 3], cx).unwrap();
        for (a, b) in force.into_iter().zip([2.0, -8.0, 3.0]) { close(a, b, 1e-12); }
        // Surface centroid (2,.5,.5), body centroid (1,.5,.5), node 0 at origin.
        for (a, b) in moment.into_iter().zip([5.5, -6.0, -8.0]) { close(a, b, 1e-12); }
        let densities = vec![[0.0, -4.0, 0.0]; mesh.tets.len()];
        let per_element = g.assemble(&[], VolumeLoad3::PerElement(&densities), &[], cx).unwrap();
        let uniform = g.assemble(&[], VolumeLoad3::Uniform([0.0, -4.0, 0.0]), &[], cx).unwrap();
        assert_eq!(per_element.nodal_forces_n(), uniform.nodal_forces_n());
    });
}

#[test]
fn g1_surface_tension_and_thermal_expansion_reach_native_stress_and_displacement() {
    for n in [1, 2] { for expansion in [0.0, 0.0003] {
        let mesh = block([n; 3], [2.0, 1.0, 1.0]);
        let material = material(); let thermal = thermal(expansion); let fixed = rollers(&mesh.nodes);
        let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
            materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget: Default::default() };
        with_cx(|cx| {
            let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 1000, cx).unwrap();
            let right: Vec<_> = g.faces().iter().enumerate()
                .filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
            let case = g.assemble(&[SurfaceLoad3::Traction { faces: &right, pascals: [100.0, 0.0, 0.0] }],
                VolumeLoad3::None, &[], cx).unwrap();
            let answer = case.solve(&problem, TetThermalStrainField::Uniform(&thermal),
                TetStaticSolveConfig::default(), ContentHash([9; 32]), 0.01, cx).unwrap();
            for (p, u) in mesh.nodes.iter().zip(answer.displacement.displacement_m()) {
                close(u[0], (expansion + 1e-4)*p[0], 1e-10);
                close(u[1], (expansion - 2.5e-5)*p[1], 1e-10);
                close(u[2], (expansion - 2.5e-5)*p[2], 1e-10);
            }
            for stress in answer.stress.elements() {
                for (a, b) in stress.stress_mandel_pa().iter().zip([100.0, 0.0, 0.0, 0.0, 0.0, 0.0]) { close(*a, b, 1e-5); }
            }
            close(answer.support_reactions_n.iter().map(|r| r[0]).sum(), -100.0, 1e-6);
            assert!(answer.free_force_residual_n < 1e-5);
            let native = problem.solve_thermal_displacement_with_nodal_forces(
                TetThermalStrainField::Uniform(&thermal), case.nodal_forces_n(), Default::default(), cx).unwrap();
            assert_eq!(native, answer.displacement);
        });
    } }
}

#[test]
fn g1_closed_pressure_produces_hydrostatic_compression_and_fixed_loads_remain_reactions() {
    let mesh = block([1; 3], [1.0; 3]);
    let material = material(); let thermal = thermal(0.0); let fixed = rollers(&mesh.nodes);
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 100, cx).unwrap();
        let faces: Vec<_> = (0..g.faces().len()).collect();
        let case = g.assemble(&[SurfaceLoad3::Pressure { faces: &faces, pascals: 80.0 }], VolumeLoad3::None, &[], cx).unwrap();
        let mut problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
            materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget: Default::default() };
        let answer = case.solve(&problem, TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([9; 32]), 0.01, cx).unwrap();
        for (p, u) in mesh.nodes.iter().zip(answer.displacement.displacement_m()) { for c in 0..3 { close(u[c], -4e-5*p[c], 1e-11); } }
        for stress in answer.stress.elements() { for (a, b) in stress.stress_mandel_pa().iter().zip([-80.0, -80.0, -80.0, 0.0, 0.0, 0.0]) { close(*a, b, 1e-6); } }
        for x in answer.support_reactions_n.iter().flatten() { close(*x, 0.0, 1e-6); }
        let all_fixed: Vec<_> = (0..mesh.nodes.len()*3).collect(); problem.fixed_dofs = &all_fixed;
        let clamped = case.solve(&problem, TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([9; 32]), 0.01, cx).unwrap();
        assert_eq!(clamped.displacement.iterations(), 0);
        assert!(clamped.displacement.displacement_m().iter().flatten().all(|v| *v == 0.0));
        for &(dof, value) in case.nodal_forces_n() { close(clamped.support_reactions_n[dof/3][dof%3], -value, 1e-12); }
    });
}

#[test]
fn g3_face_winding_and_tet_order_do_not_change_pressure_resultants() {
    let mesh = block([2; 3], [1.0; 3]);
    let mut reversed = mesh.tets.clone(); reversed.reverse(); for tet in &mut reversed { tet.swap(0, 1); }
    with_cx(|cx| {
        let a = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 1000, cx).unwrap();
        let b = TetLoadGeometry3::new(&mesh.nodes, &reversed, Default::default(), 1000, cx).unwrap();
        assert_eq!(a.faces().len(), b.faces().len());
        for (x, y) in a.faces().iter().zip(b.faces()) {
            assert_eq!(x.nodes, y.nodes);
            for c in 0..3 { close(x.outward_normal[c], y.outward_normal[c], 1e-14); }
        }
        let face = [3];
        let a = a.assemble(&[SurfaceLoad3::Pressure { faces: &face, pascals: 17.0 }], VolumeLoad3::None, &[], cx).unwrap();
        let b = b.assemble(&[SurfaceLoad3::Pressure { faces: &face, pascals: 17.0 }], VolumeLoad3::None, &[], cx).unwrap();
        for ((i, x), (j, y)) in a.nodal_forces_n().iter().zip(b.nodal_forces_n()) { assert_eq!(i, j); close(*x, *y, 1e-14); }
    });
}

#[test]
fn g0_bad_geometry_duplicate_loads_and_foreign_meshes_refuse() {
    let mesh = block([1; 3], [1.0; 3]);
    with_cx(|cx| {
        let budget = TetAssemblyBudget::default();
        assert!(matches!(TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, budget, 23, cx), Err(LoadError3::Budget)));
        let mut bad = mesh.tets.clone(); bad.push(mesh.tets[0]);
        assert!(matches!(TetLoadGeometry3::new(&mesh.nodes, &bad, budget, 100, cx), Err(LoadError3::Face { .. })));
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, budget, 100, cx).unwrap();
        for ids in [vec![], vec![0, 0], vec![g.faces().len()]] {
            assert!(g.assemble(&[SurfaceLoad3::Pressure { faces: &ids, pascals: 1.0 }], VolumeLoad3::None, &[], cx).is_err());
        }
        assert!(g.assemble(&[], VolumeLoad3::PerElement(&[]), &[], cx).is_err());
        assert!(g.assemble(&[], VolumeLoad3::Uniform([f64::NAN, 0.0, 0.0]), &[], cx).is_err());
        assert!(g.assemble(&[], VolumeLoad3::None, &[(0, 1.0), (0, 2.0)], cx).is_err());
        let material = material(); let thermal = thermal(0.0); let fixed = rollers(&mesh.nodes);
        let mut moved = mesh.nodes.clone(); moved[0][0] = -0.1;
        let p = TetLinearElasticProblem { nodes_m: &moved, tetrahedra: &mesh.tets,
            materials: TetMaterialField::Uniform(&material), fixed_dofs: &fixed, budget };
        let case = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        assert!(case.solve(&p, TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([9; 32]), 0.1, cx).is_err());
    });
}

#[test]
fn g4_cancelled_assembly_publishes_no_partial_load_and_can_resume() {
    let mesh = block([1; 3], [1.0; 3]);
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 100, cx).unwrap();
        let gate = fs_exec::CancelGate::new_clock_free(); gate.request();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let stopped = fs_exec::Cx::new(&gate, arena, fs_exec::StreamKey {
                seed: 58, kernel_id: 58, tile: 0, iteration: 0,
            }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
            assert!(matches!(TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 100, &stopped), Err(LoadError3::Cancelled)));
            assert!(matches!(g.assemble(&[], VolumeLoad3::Uniform([1.0; 3]), &[], &stopped), Err(LoadError3::Cancelled)));
        });
        let resumed = g.assemble(&[], VolumeLoad3::Uniform([1.0; 3]), &[], cx).unwrap();
        for v in resumed.resultant([0.0; 3], cx).unwrap().0 { close(v, 1.0, 1e-14); }
    });
}
