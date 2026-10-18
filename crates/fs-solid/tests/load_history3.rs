//! Typed sources reach the same native tetrahedral equilibrium and reactions.
#[path = "tet10/support.rs"]
mod support;
use support::{block, with_cx};
use fs_blake3::ContentHash;
use fs_ga::Vec3;
use fs_qty::{Dims, QtyAny};
use fs_scenario::{BcKind, BcValue, BoundaryCondition, Frame, FrameId, FrameMotion, FrameTree,
    Physics, TimeSignal, WORLD};
use fs_scenario::payload::{OrientationParity, OutsideDomainPolicy, Payload, PayloadId,
    PayloadMeta, QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload};
use fs_solid::linear3::{TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetThermalStrainField, TetThermalStrainState};
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::history::{ACCELERATION_DIMS, FORCE_DENSITY_DIMS, TRACTION_DIMS,
    LoadHistoryError3, SurfaceHistory3, TetLoadHistory3, VolumeHistory3};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
fn vector_source(dims: Dims, frame: FrameId, a: [f64; 3], b: [f64; 3], outside: OutsideDomainPolicy) -> VectorPayload {
    let meta = PayloadMeta::new(QuantityContract::Dimensions(dims),
        PayloadId::new(if frame == WORLD { "basis/world-cartesian" } else { "basis/frame-cartesian" }).unwrap(),
        frame, OrientationParity::Even, ReferenceSemantics::Continuous).unwrap();
    VectorPayload::new(meta, SampleSource::table(vec![QtyAny::new(0.0, TIME), QtyAny::new(1.0, TIME)],
        vec![a.map(|v| QtyAny::new(v, dims)).to_vec(), b.map(|v| QtyAny::new(v, dims)).to_vec()],
        TableInterpolation::Linear, outside).unwrap()).unwrap()
}
fn material(rho: f64) -> TetElasticMaterial { TetElasticMaterial::try_new(rho, 1e6, 0.25, [7; 32]).unwrap() }
fn thermal() -> TetThermalStrainState {
    TetThermalStrainState::try_new(300.0, 300.0, [0.0; 6], [7; 32], ContentHash([8; 32]), 0.1).unwrap()
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a:e} != {b:e}"); }
fn rotating() -> FrameTree {
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), name: "loading-direction".into(), parent: WORLD,
        motion: FrameMotion::Rotating { axis: [0.0, 0.0, 1.0], center: Vec3::new(1e6, -2e6, 0.0),
            rate: QtyAny::new(std::f64::consts::FRAC_PI_2, Dims([0, 0, -1, 0, 0, 0])) } });
    frames
}

#[test]
fn g3_frame_keeps_physical_time_after_the_vector_table_clamps() {
    let mesh = block([1; 3], [1.0; 3]); let mat = material(2.0); let frames = rotating();
    let value = vector_source(TRACTION_DIMS, FrameId(1), [10.0, 0.0, 0.0], [20.0, 0.0, 0.0], OutsideDomainPolicy::Clamp);
    let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&mat), fixed_dofs: &[], budget: Default::default() };
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 100, cx).unwrap();
        let right: Vec<_> = g.faces().iter().enumerate().filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
        let rows = [SurfaceHistory3::Traction { faces: &right, value: &value }];
        let history = TetLoadHistory3 { surfaces: &rows, volume: VolumeHistory3::Fixed(VolumeLoad3::None), nodal_forces_n: &[], frames: &frames };
        for (t, expected) in [(0.0, [10.0, 0.0, 0.0]), (1.0, [0.0, 20.0, 0.0]), (2.0, [-20.0, 0.0, 0.0])] {
            let sample = history.at(&g, &p, t, cx).unwrap(); assert_eq!(sample.time_s(), t);
            for (a, b) in sample.resultant([0.0; 3], cx).unwrap().0.into_iter().zip(expected) { close(a, b, 1e-10); }
        }
    });
}

#[test]
fn g1_scheduled_traction_pressure_and_body_force_match_manual_native_endpoints() {
    let mesh = block([1; 3], [1.0; 3]); let mat = material(2.0); let frames = FrameTree::new();
    let fixed: Vec<_> = mesh.nodes.iter().enumerate().filter(|(_, p)| p[0] == 0.0).flat_map(|(i, _)| [3*i, 3*i+1, 3*i+2]).collect();
    let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets, materials: TetMaterialField::Uniform(&mat), fixed_dofs: &fixed, budget: Default::default() };
    let traction = vector_source(TRACTION_DIMS, WORLD, [10.0, 2.0, 0.0], [30.0, 6.0, 8.0], OutsideDomainPolicy::Refuse);
    let body = vector_source(FORCE_DENSITY_DIMS, WORLD, [0.0, 0.0, -2.0], [0.0, 0.0, -4.0], OutsideDomainPolicy::Refuse);
    let pressure = TimeSignal::Ramp { t_start: 0.0, t_end: 1.0, from: QtyAny::new(4.0, TRACTION_DIMS), to: QtyAny::new(12.0, TRACTION_DIMS) };
    let thermal = thermal();
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 100, cx).unwrap();
        let right: Vec<_> = g.faces().iter().enumerate().filter_map(|(i, f)| (f.outward_normal[0] > 0.99).then_some(i)).collect();
        let top: Vec<_> = g.faces().iter().enumerate().filter_map(|(i, f)| (f.outward_normal[2] > 0.99).then_some(i)).collect();
        let rows = [SurfaceHistory3::Traction { faces: &right, value: &traction }, SurfaceHistory3::Pressure { faces: &top, value: &pressure }];
        let history = TetLoadHistory3 { surfaces: &rows, volume: VolumeHistory3::ForceDensity(&body), nodal_forces_n: &[(0, 5.0)], frames: &frames };
        for t in [0.0, 0.25, 0.5, 1.0] {
            let sampled = history.at(&g, &p, t, cx).unwrap();
            let manual = g.assemble(&[
                SurfaceLoad3::Traction { faces: &right, pascals: [10.0+20.0*t, 2.0+4.0*t, 8.0*t] },
                SurfaceLoad3::Pressure { faces: &top, pascals: 4.0+8.0*t },
            ], VolumeLoad3::Uniform([0.0, 0.0, -2.0-2.0*t]), &[(0, 5.0)], cx).unwrap();
            assert_eq!(sampled.nodal_forces_n(), manual.nodal_forces_n());
            let a = sampled.solve(TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([9; 32]), 0.1, cx).unwrap();
            let b = manual.solve(&p, TetThermalStrainField::Uniform(&thermal), Default::default(), ContentHash([9; 32]), 0.1, cx).unwrap();
            assert_eq!(a.displacement, b.displacement); assert_eq!(a.stress, b.stress);
            assert_eq!(a.support_reactions_n, b.support_reactions_n);
            assert!(a.free_force_residual_n < 1e-7);
        }
    });
}

#[test]
fn g0_acceleration_uses_the_actual_uniform_and_per_element_material_mass() {
    let mesh = block([1; 3], [1.0; 3]); let frames = FrameTree::new();
    let acceleration = vector_source(ACCELERATION_DIMS, WORLD, [0.0, 0.0, -2.0], [0.0, 0.0, -2.0], OutsideDomainPolicy::Clamp);
    let history = TetLoadHistory3 { surfaces: &[], volume: VolumeHistory3::Acceleration(&acceleration), nodal_forces_n: &[], frames: &frames };
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, Default::default(), 100, cx).unwrap();
        for rho in [2.0, 5.0] {
            let mat = material(rho); let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets, materials: TetMaterialField::Uniform(&mat), fixed_dofs: &[], budget: Default::default() };
            let force = history.at(&g, &p, 0.5, cx).unwrap().resultant([0.0; 3], cx).unwrap().0;
            close(force[2], -2.0*rho, 1e-12);
        }
        let mats: Vec<_> = (1..=mesh.tets.len()).map(|i| material(i as f64)).collect();
        let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets, materials: TetMaterialField::PerElement(&mats), fixed_dofs: &[], budget: Default::default() };
        close(history.at(&g, &p, 0.0, cx).unwrap().resultant([0.0; 3], cx).unwrap().0[2], -7.0, 1e-12);
    });
}

#[test]
fn g0_typed_neumann_row_remains_distinct_from_scalar_traction() {
    let value = vector_source(TRACTION_DIMS, WORLD, [1.0, 2.0, 3.0], [2.0, 3.0, 4.0], OutsideDomainPolicy::Refuse);
    let mut row = BoundaryCondition { region: "loaded-face".into(), physics: Physics::Elasticity, kind: BcKind::Neumann,
        value: Some(BcValue::Typed(Payload::Vector(value))), compatibility: None, frame: 0 };
    assert!(SurfaceHistory3::scenario_traction(&[0], &row).is_ok());
    row.kind = BcKind::Traction; assert!(SurfaceHistory3::scenario_traction(&[0], &row).is_err());
    row.kind = BcKind::Neumann; row.frame = 7; assert!(SurfaceHistory3::scenario_traction(&[0], &row).is_err());
}

#[test]
fn g0_invalid_time_units_and_duplicate_faces_never_publish_partial_histories() {
    let mesh = block([1; 3], [1.0; 3]); let mat = material(2.0); let frames = FrameTree::new();
    let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets, materials: TetMaterialField::Uniform(&mat), fixed_dofs: &[], budget: Default::default() };
    let wrong = vector_source(ACCELERATION_DIMS, WORLD, [1.0; 3], [2.0; 3], OutsideDomainPolicy::Refuse);
    let value = vector_source(TRACTION_DIMS, WORLD, [1.0; 3], [2.0; 3], OutsideDomainPolicy::Refuse);
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 100, cx).unwrap();
        let rows = [SurfaceHistory3::Traction { faces: &[0], value: &wrong }];
        let mut h = TetLoadHistory3 { surfaces: &rows, volume: VolumeHistory3::Fixed(VolumeLoad3::None), nodal_forces_n: &[], frames: &frames };
        assert!(matches!(h.at(&g, &p, 0.5, cx), Err(LoadHistoryError3::Quantity { .. })));
        let good = [SurfaceHistory3::Traction { faces: &[0], value: &value }]; h.surfaces = &good;
        assert!(h.at(&g, &p, f64::NAN, cx).is_err()); assert!(h.at(&g, &p, 2.0, cx).is_err());
        let before = h.at(&g, &p, 0.5, cx).unwrap().nodal_forces_n().to_vec();
        let dup = [good[0], good[0]]; h.surfaces = &dup; assert!(h.at(&g, &p, 0.5, cx).is_err());
        h.surfaces = &good; assert_eq!(h.at(&g, &p, 0.5, cx).unwrap().nodal_forces_n(), before);
    });
}

#[test]
fn g4_cancelled_history_leaves_sources_reusable() {
    let mesh = block([1; 3], [1.0; 3]); let mat = material(2.0); let frames = FrameTree::new();
    let value = vector_source(FORCE_DENSITY_DIMS, WORLD, [0.0, 0.0, -1.0], [0.0, 0.0, -2.0], OutsideDomainPolicy::Clamp);
    let history = TetLoadHistory3 { surfaces: &[], volume: VolumeHistory3::ForceDensity(&value), nodal_forces_n: &[], frames: &frames };
    let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets, materials: TetMaterialField::Uniform(&mat), fixed_dofs: &[], budget: Default::default() };
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 100, cx).unwrap();
        let gate = fs_exec::CancelGate::new_clock_free(); gate.request();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let stopped = fs_exec::Cx::new(&gate, arena, fs_exec::StreamKey {
                seed: 58, kernel_id: 58, tile: 0, iteration: 0,
            }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
            assert!(matches!(history.at(&g, &p, 0.5, &stopped), Err(LoadHistoryError3::Cancelled)));
        });
        close(history.at(&g, &p, 0.5, cx).unwrap().resultant([0.0; 3], cx).unwrap().0[2], -1.5, 1e-12);
    });
}
