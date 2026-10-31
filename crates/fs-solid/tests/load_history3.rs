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
use fs_solid::loads3::fields::AccelerationField3;
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

fn force_component(values: &[(usize, f64)], dof: usize) -> f64 {
    values.binary_search_by_key(&dof, |v| v.0).map_or(0.0, |i| values[i].1)
}

#[test]
fn g0_spatial_acceleration_applies_full_consistent_material_mass_before_constraints() {
    let mesh = block([1; 3], [1.0; 3]);
    let frames = FrameTree::new();
    let materials: Vec<_> = (0..mesh.tets.len()).map(|i| material(2.0+i as f64)).collect();
    let field: Vec<_> = mesh.nodes.iter().map(|p| [-p[0], -2.0*p[1], -3.0*p[2]]).collect();
    let per_element: Vec<_> = mesh.tets.iter().map(|tet| tet.map(|node| field[node])).collect();
    let problem = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::PerElement(&materials), fixed_dofs: &[], budget: Default::default() };
    let history = TetLoadHistory3 { surfaces: &[], nodal_forces_n: &[], frames: &frames,
        volume: VolumeHistory3::FixedAcceleration(AccelerationField3::PerNode(&field)) };
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, problem.budget, 100, cx).unwrap();
        let sampled = history.at(&g, &problem, 0.0, cx).unwrap();
        let assembly = problem.assemble(cx).unwrap();
        for dof in 0..mesh.nodes.len()*3 {
            let (columns, values) = assembly.mass.row(dof);
            let expected: f64 = columns.iter().zip(values)
                .map(|(&column, value)| value*field[column/3][column%3]).sum();
            close(force_component(sampled.nodal_forces_n(), dof), expected, 1e-13);
        }
        let element_history = TetLoadHistory3 {
            volume: VolumeHistory3::FixedAcceleration(AccelerationField3::PerElement(&per_element)), ..history
        };
        assert_eq!(element_history.at(&g, &problem, 0.0, cx).unwrap().nodal_forces_n(), sampled.nodal_forces_n());
        // Node zero's acceleration vanishes, yet the consistent off-diagonal
        // terms carry a finite load there. Supports must retain that load.
        assert!(force_component(sampled.nodal_forces_n(), 0).abs() > 0.01);
        let constrained = TetLinearElasticProblem { fixed_dofs: &[0, 1, 2], ..problem };
        assert_eq!(history.at(&g, &constrained, 0.0, cx).unwrap().nodal_forces_n(), sampled.nodal_forces_n());

        // The very same history reaches the actual sparse inertia owner.
        let zero = vec![[0.0; 3]; mesh.nodes.len()];
        let run = history.begin_transient(&g, &problem, None, 0.0, &zero, &zero, Default::default(), cx).unwrap();
        for (actual, expected) in run.nodal_motion(cx).unwrap()[2].iter().flatten().zip(field.iter().flatten()) {
            close(*actual, *expected, 1e-9);
        }
    });
}

#[test]
fn g1_spatial_gravity_preserves_analytic_resultant_moment_and_rotated_support_reactions() {
    let mesh = block([1; 3], [1.0; 3]); let mat = material(2.0); let thermal = thermal();
    let frames = FrameTree::new();
    let fixed: Vec<_> = mesh.nodes.iter().enumerate().filter(|(_, p)| p[0] == 0.0)
        .flat_map(|(i, _)| [3*i, 3*i+1, 3*i+2]).collect();
    for turned in [false, true] {
        let rotate = |v: [f64; 3]| if turned { [-v[1], v[0], v[2]] } else { v };
        let nodes: Vec<_> = mesh.nodes.iter().copied().map(rotate).collect();
        // g = -grad(x²/2 + y² + 3z²/2), sampled in the reference frame.
        let field: Vec<_> = mesh.nodes.iter().map(|p| rotate([-p[0], -2.0*p[1], -3.0*p[2]])).collect();
        let p = TetLinearElasticProblem { nodes_m: &nodes, tetrahedra: &mesh.tets,
            materials: TetMaterialField::Uniform(&mat), fixed_dofs: &fixed, budget: Default::default() };
        let history = TetLoadHistory3 { surfaces: &[], nodal_forces_n: &[], frames: &frames,
            volume: VolumeHistory3::FixedAcceleration(AccelerationField3::PerNode(&field)) };
        with_cx(|cx| {
            let g = TetLoadGeometry3::new(&nodes, &mesh.tets, p.budget, 100, cx).unwrap();
            let sample = history.at(&g, &p, 0.0, cx).unwrap();
            let (force, moment) = sample.resultant([0.0; 3], cx).unwrap();
            for (actual, expected) in force.into_iter().zip(rotate([-1.0, -2.0, -3.0])) { close(actual, expected, 1e-13); }
            for (actual, expected) in moment.into_iter().zip(rotate([-0.5, 1.0, -0.5])) { close(actual, expected, 1e-13); }
            let solution = sample.solve(TetThermalStrainField::Uniform(&thermal), Default::default(),
                ContentHash([9; 32]), 0.01, cx).unwrap();
            let mut support_moment = [0.0; 3];
            for (x, r) in nodes.iter().zip(&solution.support_reactions_n) {
                support_moment[0] += x[1]*r[2] - x[2]*r[1];
                support_moment[1] += x[2]*r[0] - x[0]*r[2];
                support_moment[2] += x[0]*r[1] - x[1]*r[0];
            }
            for c in 0..3 {
                close(solution.support_reactions_n.iter().map(|r| r[c]).sum::<f64>() + force[c], 0.0, 1e-6);
                close(support_moment[c] + moment[c], 0.0, 1e-6);
            }
        });
    }
}

#[test]
fn g0_fixed_gravity_shares_uniform_history_and_refuses_malformed_spatial_fields() {
    let mesh = block([1; 3], [1.0; 3]); let mat = material(2.0); let frames = FrameTree::new();
    let p = TetLinearElasticProblem { nodes_m: &mesh.nodes, tetrahedra: &mesh.tets,
        materials: TetMaterialField::Uniform(&mat), fixed_dofs: &[], budget: Default::default() };
    let acceleration = [1.0, -2.0, -9.80665];
    let source = vector_source(ACCELERATION_DIMS, WORLD, acceleration, acceleration, OutsideDomainPolicy::Clamp);
    let history = TetLoadHistory3 { surfaces: &[], nodal_forces_n: &[], frames: &frames,
        volume: VolumeHistory3::FixedAcceleration(AccelerationField3::Uniform(acceleration)) };
    with_cx(|cx| {
        let g = TetLoadGeometry3::new(&mesh.nodes, &mesh.tets, p.budget, 100, cx).unwrap();
        let expected = history.at(&g, &p, 0.5, cx).unwrap();
        let typed = TetLoadHistory3 { volume: VolumeHistory3::Acceleration(&source), ..history };
        assert_eq!(typed.at(&g, &p, 0.5, cx).unwrap().nodal_forces_n(), expected.nodal_forces_n());
        let zero = TetLoadHistory3 { volume: VolumeHistory3::FixedAcceleration(AccelerationField3::Uniform([0.0; 3])), ..history };
        assert!(zero.at(&g, &p, 0.5, cx).unwrap().nodal_forces_n().is_empty());
        let bad_values = vec![[f64::INFINITY; 3]; mesh.nodes.len()];
        for bad in [AccelerationField3::PerNode(&[]), AccelerationField3::PerElement(&[]),
            AccelerationField3::Uniform([f64::NAN; 3]), AccelerationField3::PerNode(&bad_values)] {
            let malformed = TetLoadHistory3 { volume: VolumeHistory3::FixedAcceleration(bad), ..history };
            assert!(malformed.at(&g, &p, 0.0, cx).is_err());
        }
        let invalid_material_count = TetLinearElasticProblem { materials: TetMaterialField::PerElement(&[]), ..p };
        assert!(history.at(&g, &invalid_material_count, 0.0, cx).is_err());
        let mut invalid_material = mat.clone(); invalid_material.density_kg_m3 = -1.0;
        let invalid_density = TetLinearElasticProblem { materials: TetMaterialField::Uniform(&invalid_material), ..p };
        assert!(history.at(&g, &invalid_density, 0.0, cx).is_err());
        assert_eq!(history.at(&g, &p, 0.5, cx).unwrap().nodal_forces_n(), expected.nodal_forces_n());
    });
}
