//! Same-time force and support histories, through the original P1 mechanics.
#[path = "tet10/support.rs"]
mod support;
use support::{block, material, with_cx};
use fs_blake3::ContentHash;
use fs_ga::Vec3;
use fs_qty::{Dims, QtyAny};
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, TimeSignal, WORLD};
use fs_scenario::payload::{OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload};
use fs_solid::linear3::*;
use fs_solid::loads3::{LoadError3, SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::history::{LoadHistoryError3, SurfaceHistory3, TetLoadHistory3,
    VolumeHistory3, TRACTION_DIMS, FORCE_DENSITY_DIMS};
use fs_solid::loads3::history::supports::{DISPLACEMENT_DIMS, SupportMotion3};

fn source(dims: Dims, frame: FrameId, a: [f64; 3], b: [f64; 3], outside: OutsideDomainPolicy)
    -> VectorPayload {
    VectorPayload::new(PayloadMeta::new(QuantityContract::Dimensions(dims),
        PayloadId::new(if frame == WORLD { "basis/world-cartesian" } else { "basis/frame-cartesian" }).unwrap(),
        frame, OrientationParity::Even, ReferenceSemantics::Continuous).unwrap(),
        SampleSource::table(vec![QtyAny::new(0.0, Dims([0,0,1,0,0,0])), QtyAny::new(1.0, Dims([0,0,1,0,0,0]))],
            vec![a.map(|v| QtyAny::new(v, dims)).to_vec(), b.map(|v| QtyAny::new(v, dims)).to_vec()],
            TableInterpolation::Linear, outside).unwrap()).unwrap()
}
fn motion(a: [f64; 3], b: [f64; 3]) -> VectorPayload {
    source(DISPLACEMENT_DIMS, WORLD, a, b, OutsideDomainPolicy::Refuse)
}
fn state(mat: &TetElasticMaterial, alpha: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(310.0, 300.0, [alpha, alpha, alpha, 0.0, 0.0, 0.0],
        mat.material_state_identity, ContentHash([17; 32]), 0.1).unwrap()
}
fn face_nodes(mesh: &support::Mesh, axis: usize, coordinate: f64) -> Vec<usize> {
    mesh.nodes.iter().enumerate().filter_map(|(i, p)| (p[axis] == coordinate).then_some(i)).collect()
}
fn problem<'a>(mesh: &'a support::Mesh, mat: &'a TetElasticMaterial, fixed: &'a [usize])
    -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(mat), fixed_dofs: fixed, budget: Default::default() }
}
fn empty_history(frames: &FrameTree) -> TetLoadHistory3<'_> {
    TetLoadHistory3 { surfaces: &[], volume: VolumeHistory3::Fixed(VolumeLoad3::None),
        nodal_forces_n: &[], frames }
}
fn rotating() -> FrameTree {
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), parent: WORLD, name: "actuator-vector".into(),
        motion: FrameMotion::Rotating { axis: [0.0, 0.0, 1.0], center: Vec3::new(1e6, -2e6, 0.0),
            rate: QtyAny::new(std::f64::consts::FRAC_PI_2, Dims([0,0,-1,0,0,0])) } });
    frames
}

#[test]
fn g3_load_and_motion_histories_match_manual_combined_endpoints() {
    with_cx(|cx| {
        let mesh = block([2; 3], [1.0; 3]); let mat = material(1e6, 0.25);
        let left = face_nodes(&mesh, 0, 0.0);
        let fixed: Vec<_> = left.iter().flat_map(|&n| [3*n, 3*n+1, 3*n+2]).collect();
        let p = problem(&mesh, &mat, &fixed); let frames = FrameTree::new();
        let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 4*mesh.tets.len(), cx).unwrap();
        let right: Vec<_> = geometry.faces().iter().enumerate().filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
        let top: Vec<_> = geometry.faces().iter().enumerate().filter_map(|(i, f)| (f.outward_normal[2] > 0.99).then_some(i)).collect();
        let traction = source(TRACTION_DIMS, WORLD, [2.0, 1.0, 0.0], [4.0, 3.0, 2.0], OutsideDomainPolicy::Refuse);
        let body = source(FORCE_DENSITY_DIMS, WORLD, [0.0, 0.0, -0.5], [0.0, 0.0, -1.0], OutsideDomainPolicy::Refuse);
        let pressure = TimeSignal::Ramp { t_start: 0.0, t_end: 1.0,
            from: QtyAny::new(1.0, TRACTION_DIMS), to: QtyAny::new(2.0, TRACTION_DIMS) };
        let moving = motion([0.0; 3], [1.0/1024.0, -1.0/2048.0, 1.0/4096.0]);
        let rows = [SupportMotion3 { nodes: &left, components: [true; 3], displacement: &moving }];
        let surfaces = [SurfaceHistory3::Traction { faces: &right, value: &traction },
            SurfaceHistory3::Pressure { faces: &top, value: &pressure }];
        let histories = TetLoadHistory3 { surfaces: &surfaces, volume: VolumeHistory3::ForceDensity(&body),
            nodal_forces_n: &[(0, 5.0)], frames: &frames };
        let thermal = state(&mat, 0.0001);
        for t in [0.0, 0.25, 0.5, 1.0] {
            let sampled = histories.at_with_supports(&geometry, &p, &rows, t, 0.1, cx).unwrap();
            let manual = geometry.assemble(&[
                SurfaceLoad3::Traction { faces: &right, pascals: [2.0+2.0*t, 1.0+2.0*t, 2.0*t] },
                SurfaceLoad3::Pressure { faces: &top, pascals: 1.0+t },
            ], VolumeLoad3::Uniform([0.0, 0.0, -0.5-0.5*t]), &[(0, 5.0)], cx).unwrap();
            let displacement = if t == 0.0 { [0.0;3] } else { [t/1024.0, -t/2048.0, t/4096.0] };
            let g: Vec<_> = left.iter().flat_map(|&n| displacement
                .into_iter().enumerate().map(move |(c, v)| (3*n+c, v))).collect();
            assert_eq!(sampled.time_s(), t);
            assert_eq!(sampled.prescribed_dofs_m(), g);
            assert_eq!(sampled.nodal_forces_n(), manual.nodal_forces_n());
            let a = sampled.solve(TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([18;32]), 0.1, cx).unwrap();
            let b = manual.solve_with_prescribed_dofs(&p, TetThermalStrainField::Uniform(&thermal),
                &g, 0.1, Default::default(), ContentHash([18;32]), 0.1, cx).unwrap();
            assert_eq!(a.displacement, b.displacement); assert_eq!(a.stress, b.stress);
            assert_eq!(a.support_reactions_n, b.support_reactions_n);
            assert!(a.free_force_residual_n < 1e-6);
        }
    });
}

#[test]
fn g1_driven_support_extension_without_external_force_changes_actual_stress() {
    with_cx(|cx| {
        let mesh = block([2;3], [1.0;3]); let mat = material(2000.0, 0.25);
        let sets = [face_nodes(&mesh, 0, 0.0), face_nodes(&mesh, 0, 1.0),
            face_nodes(&mesh, 1, 0.0), face_nodes(&mesh, 2, 0.0)];
        let mut fixed = Vec::new();
        for (s, c) in [(0,0), (1,0), (2,1), (3,2)] { fixed.extend(sets[s].iter().map(|n| 3*n+c)); }
        fixed.sort_unstable();
        let p = problem(&mesh, &mat, &fixed); let frames = FrameTree::new();
        let history = empty_history(&frames); let thermal = state(&mat, 0.0);
        let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 4*mesh.tets.len(), cx).unwrap();
        let zero = motion([0.0;3], [0.0;3]); let stretch = motion([0.0;3], [1.0/512.0, 0.0, 0.0]);
        let rows = [
            SupportMotion3 { nodes: &sets[0], components: [true,false,false], displacement: &zero },
            SupportMotion3 { nodes: &sets[1], components: [true,false,false], displacement: &stretch },
            SupportMotion3 { nodes: &sets[2], components: [false,true,false], displacement: &zero },
            SupportMotion3 { nodes: &sets[3], components: [false,false,true], displacement: &zero },
        ];
        for t in [0.0, 0.25, 0.5, 1.0] {
            let endpoint = history.at_with_supports(&geometry, &p, &rows, t, 0.1, cx).unwrap();
            assert!(endpoint.nodal_forces_n().is_empty());
            let result = endpoint.solve(TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([18;32]), 0.1, cx).unwrap();
            for (point, u) in mesh.nodes.iter().zip(result.displacement.displacement_m()) {
                let expected = [point[0]*t/512.0, -0.25*point[1]*t/512.0, -0.25*point[2]*t/512.0];
                for c in 0..3 { assert!((u[c]-expected[c]).abs() < 1e-10); }
            }
            for s in result.stress.elements() { assert!((s.stress_mandel_pa()[0]-2000.0*t/512.0).abs() < 1e-6); }
        }
    });
}

#[test]
fn g3_clamped_and_periodic_motion_values_keep_the_frames_physical_time() {
    with_cx(|cx| {
        let mesh = block([2;3], [1.0;3]); let mat = material(1000.0, 0.25); let frames = rotating();
        let left = face_nodes(&mesh, 0, 0.0);
        let fixed: Vec<_> = left.iter().flat_map(|&n| [3*n,3*n+1,3*n+2]).collect();
        let p = problem(&mesh, &mat, &fixed); let history = empty_history(&frames);
        let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 4*mesh.tets.len(), cx).unwrap();
        let thermal = state(&mat, 0.0);
        for (policy, t, phase) in [(OutsideDomainPolicy::Clamp, 2.0, 1.0),
            (OutsideDomainPolicy::Periodic, 2.5, 0.5), (OutsideDomainPolicy::Periodic, -0.5, 0.5)] {
            let value = source(DISPLACEMENT_DIMS, FrameId(1), [1.0/1024.0,0.0,0.0],
                [2.0/1024.0,0.0,0.0], policy);
            let rows = [SupportMotion3 { nodes: &left, components: [true;3], displacement: &value }];
            let sample = history.at_with_supports(&geometry, &p, &rows, t, 0.1, cx).unwrap();
            let angle: f64 = std::f64::consts::FRAC_PI_2*t;
            let expected = [(1.0+phase)/1024.0*angle.cos(), (1.0+phase)/1024.0*angle.sin(), 0.0];
            for &(dof, v) in sample.prescribed_dofs_m() { assert!((v-expected[dof%3]).abs() < 1e-14); }
            let result = sample.solve(TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([18;32]), 0.1, cx).unwrap();
            for u in result.displacement.displacement_m() {
                for c in 0..3 { assert!((u[c]-expected[c]).abs() < 1e-10); }
            }
            assert!(result.stress.elastic_energy_j() < 1e-14);
        }
    });
}

#[test]
fn g0_component_selection_does_not_create_rotating_constraints_or_move_other_fixed_dofs() {
    with_cx(|cx| {
        let mesh = block([1;3], [1.0;3]); let mat = material(1000.0, 0.0);
        let fixed: Vec<_> = (0..3*mesh.nodes.len()).collect(); let p = problem(&mesh, &mat, &fixed);
        let frames = FrameTree::new(); let history = empty_history(&frames);
        let geometry = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 24, cx).unwrap();
        let v = motion([0.001,1.0,0.0], [0.001,1.0,0.0]);
        let rows = [SupportMotion3 { nodes: &[0], components: [true,false,false], displacement: &v }];
        let sample = history.at_with_supports(&geometry, &p, &rows, 0.0, 0.01, cx).unwrap();
        assert_eq!(sample.prescribed_dofs_m(), &[(0,0.001)]);
        let result = sample.solve(TetThermalStrainField::Uniform(&state(&mat,0.0)), Default::default(), ContentHash([18;32]), 0.1, cx).unwrap();
        assert_eq!(result.displacement.displacement_m()[0], [0.001,0.0,0.0]);
        assert_eq!(result.displacement.iterations(), 0);
    });
}

#[test]
fn g0_conflicts_wrong_units_unfixed_components_and_excessive_work_refuse() {
    with_cx(|cx| {
        let mesh = block([1;3], [1.0;3]); let mat = material(1000.0,0.25);
        let p = problem(&mesh,&mat,&[0,1,2]); let frames = FrameTree::new(); let h = empty_history(&frames);
        let g = TetLoadGeometry3::new(&mesh.nodes,&mesh.tets,p.budget,24,cx).unwrap();
        let value = motion([0.001;3],[0.002;3]);
        let wrong = source(TRACTION_DIMS,WORLD,[0.0;3],[0.0;3],OutsideDomainPolicy::Refuse);
        let valid = SupportMotion3 { nodes: &[0], components: [true,false,false], displacement: &value };
        let before = h.at_with_supports(&g,&p,&[valid],0.5,0.1,cx).unwrap().prescribed_dofs_m().to_vec();
        for bad in [SupportMotion3 { nodes: &[], ..valid }, SupportMotion3 { nodes: &[1], ..valid },
            SupportMotion3 { nodes: &[8], ..valid }, SupportMotion3 { nodes: &[0,0], ..valid },
            SupportMotion3 { components: [false;3], ..valid }, SupportMotion3 { displacement: &wrong, ..valid }] {
            assert!(h.at_with_supports(&g,&p,&[bad],0.5,0.1,cx).is_err());
        }
        assert!(h.at_with_supports(&g,&p,&[valid,valid],0.5,0.1,cx).is_err());
        assert!(h.at_with_supports(&g,&p,&vec![valid;1025],0.5,0.1,cx).is_err());
        for (time,limit) in [(f64::NAN,0.1),(2.0,0.1),(0.5,0.0001),(0.5,0.0),(0.5,f64::INFINITY)] {
            assert!(h.at_with_supports(&g,&p,&[valid],time,limit,cx).is_err());
        }
        assert_eq!(h.at_with_supports(&g,&p,&[valid],0.5,0.1,cx).unwrap().prescribed_dofs_m(),before);
    });
}

#[test]
fn g4_cancelled_binding_or_solve_does_not_consume_the_history() {
    with_cx(|cx| {
        let mesh = block([1;3],[1.0;3]); let mat = material(1000.0,0.25);
        let fixed: Vec<_> = (0..3*mesh.nodes.len()).collect(); let p = problem(&mesh,&mat,&fixed);
        let frames = FrameTree::new(); let h = empty_history(&frames); let value = motion([0.0;3],[0.001;3]);
        let rows = [SupportMotion3 { nodes: &[0], components: [true;3], displacement: &value }];
        let g = TetLoadGeometry3::new(&mesh.nodes,&mesh.tets,p.budget,24,cx).unwrap();
        let sample = h.at_with_supports(&g,&p,&rows,0.5,0.1,cx).unwrap();
        let gate = fs_exec::CancelGate::new_clock_free(); gate.request();
        let pool = fs_alloc::ArenaPool::new(Default::default());
        pool.scope(|arena| {
            let stopped = fs_exec::Cx::new(&gate,arena,fs_exec::StreamKey { seed:1,kernel_id:1,tile:0,iteration:0 },
                fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
            assert!(matches!(h.at_with_supports(&g,&p,&rows,0.5,0.1,&stopped),Err(LoadHistoryError3::Cancelled)));
            assert!(matches!(sample.solve(TetThermalStrainField::Uniform(&state(&mat,0.0)),Default::default(),
                ContentHash([18;32]),0.1,&stopped),Err(LoadError3::Cancelled)));
        });
        let again = h.at_with_supports(&g,&p,&rows,0.5,0.1,cx).unwrap();
        assert_eq!(sample.prescribed_dofs_m(),again.prescribed_dofs_m());
        assert!(again.solve(TetThermalStrainField::Uniform(&state(&mat,0.0)),Default::default(),
            ContentHash([18;32]),0.1,cx).is_ok());
    });
}
