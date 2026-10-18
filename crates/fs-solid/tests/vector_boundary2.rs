//! G1/G3: typed histories drive the actual P1/Q1 elasticity owner.
use fs_ga::{Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{BcKind, BcValue, BoundaryCondition, Frame, FrameId, FrameMotion, FrameTree, Physics, WORLD};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_solid::{Formulation, LinearProblem, Mesh2, Patch, PlaneKind};
use fs_solid::vector_boundary2::{VectorBoundary2 as Bc, VectorBoundaryError2 as Error};
use std::ops::ControlFlow;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const STRESS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
fn payload(dims: Dims, frame: FrameId, source: SampleSource<Vec<QtyAny>>) -> VectorPayload {
    let meta = PayloadMeta::new(QuantityContract::Dimensions(dims), PayloadId::new("basis/frame-cartesian").unwrap(),
        frame, OrientationParity::Even, ReferenceSemantics::Continuous).unwrap();
    VectorPayload::new(meta, source).unwrap()
}
fn values(v: [f64; 3], dims: Dims) -> Vec<QtyAny> { v.map(|x| QtyAny::new(x, dims)).to_vec() }
fn fixed(v: [f64; 3], dims: Dims) -> VectorPayload {
    payload(dims, WORLD, SampleSource::fixed(values(v, dims)))
}
fn ramp(a: [f64; 3], b: [f64; 3], dims: Dims) -> VectorPayload {
    payload(dims, WORLD, SampleSource::table(
        vec![QtyAny::new(0.0, TIME), QtyAny::new(2.0, TIME)], vec![values(a, dims), values(b, dims)],
        TableInterpolation::Linear, OutsideDomainPolicy::Refuse).unwrap())
}
fn problem(mesh: &Mesh2) -> LinearProblem<'_> {
    LinearProblem { mesh, youngs: 1000.0, poisson: 0.0, plane: PlaneKind::Stress,
        formulation: Formulation::Standard, body_force: None, dirichlet: vec![],
        traction: vec![], symmetry: vec![] }
}
fn run(problem: &LinearProblem<'_>, frames: &FrameTree, rows: &[Bc<'_>], time: f64) -> Vec<[f64; 2]> {
    problem.solve_with_vector_boundaries(frames, rows, time, &mut || ControlFlow::Continue(())).unwrap()
}

#[test]
fn rotated_scheduled_traction_recovers_analytical_tension_on_both_element_families() {
    let mut frames = FrameTree::new();
    frames.add(Frame { id: FrameId(1), name: "loading frame".into(), parent: WORLD,
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2),
            translation: Vec3::new(100.0, -200.0, 300.0),
        } });
    let load = payload(STRESS, FrameId(1), SampleSource::table(
        vec![QtyAny::new(0.0, TIME), QtyAny::new(2.0, TIME)],
        vec![values([0.0, -2.0, 0.0], STRESS), values([0.0, -4.0, 0.0], STRESS)],
        TableInterpolation::Linear, OutsideDomainPolicy::Refuse).unwrap());
    let zero = fixed([0.0; 3], LENGTH);
    for mesh in [Mesh2::quads(2.0, 1.0, 4, 2), Mesh2::triangles(2.0, 1.0, 4, 2)] {
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let mut p = problem(&mesh);
            p.formulation = formulation;
            let rows = [Bc::Displacement { patch: Patch::Left, value: &zero },
                Bc::Traction { patch: Patch::Right, value: &load }];
            let u = run(&p, &frames, &rows, 1.0);
            for (point, displacement) in mesh.nodes.iter().zip(&u) {
                assert!((displacement[0] - 3.0 * point[0] / 1000.0).abs() < 1e-9);
                assert!(displacement[1].abs() < 1e-9);
            }
            assert_eq!(run(&p, &frames, &rows, 1.0), u);
            let doubled = run(&p, &frames, &rows, 2.0);
            for (a, b) in u.iter().zip(doubled) { assert!((4.0 / 3.0 * a[0] - b[0]).abs() < 1e-9); }
        }
    }
}

#[test]
fn nonzero_scheduled_supports_use_the_actual_dirichlet_rhs_lift() {
    let mesh = Mesh2::quads(2.0, 1.0, 4, 2);
    let p = problem(&mesh);
    let left = fixed([0.001, 0.0, 0.0], LENGTH);
    let right = ramp([0.001, 0.0, 0.0], [0.009, 0.0, 0.0], LENGTH);
    let rows = [Bc::Displacement { patch: Patch::Left, value: &left },
        Bc::Displacement { patch: Patch::Right, value: &right }];
    let u = run(&p, &FrameTree::new(), &rows, 1.0);
    for (point, displacement) in mesh.nodes.iter().zip(u) {
        assert!((displacement[0] - (0.001 + 0.002 * point[0])).abs() < 1e-9);
        assert!(displacement[1].abs() < 1e-9);
    }
    let translated = run(&p, &FrameTree::new(), &rows, 0.0);
    for value in translated { assert!((value[0] - 0.001).abs() < 1e-9); }
}

#[test]
fn existing_body_static_loads_and_symmetry_survive_the_vector_binding() {
    let mesh = Mesh2::quads(2.0, 1.0, 4, 2);
    let zero = |_, _| [0.0, 0.0];
    let body = |_, _| [0.0, -0.1];
    let top = |_, _| [0.2, 0.0];
    let right = |_, _| [3.0, 0.0];
    let mut p = problem(&mesh);
    p.body_force = Some(&body);
    p.dirichlet.push((Patch::Left, &zero));
    p.traction.push((Patch::Top, &top));
    p.symmetry.push((Patch::Bottom, 1));
    let load = fixed([3.0, 0.0, 0.0], STRESS);
    let actual = run(&p, &FrameTree::new(), &[Bc::Traction { patch: Patch::Right, value: &load }], 0.0);
    assert_eq!(p.traction.len(), 1, "the template is not changed");
    p.traction.push((Patch::Right, &right));
    assert_eq!(actual, p.solve().unwrap());
}

#[test]
fn canonical_typed_neumann_row_reaches_equilibrium_without_scalar_guessing() {
    let load = fixed([2.0, 0.0, 0.0], STRESS);
    let row = BoundaryCondition { region: "loaded-end".into(), physics: Physics::Elasticity,
        kind: BcKind::Neumann, value: Some(BcValue::Typed(Payload::Vector(load))),
        compatibility: None, frame: 0 };
    let mesh = Mesh2::quads(2.0, 1.0, 4, 2);
    let zero = fixed([0.0; 3], LENGTH);
    let rows = [Bc::Displacement { patch: Patch::Left, value: &zero }, Bc::scenario_traction(Patch::Right, &row).unwrap()];
    let u = run(&problem(&mesh), &FrameTree::new(), &rows, 0.0);
    for node in mesh.patch_nodes(Patch::Right) { assert!((u[node][0] - 0.004).abs() < 1e-9); }
    let mut bad = row.clone();
    bad.frame = 1;
    assert!(Bc::scenario_traction(Patch::Right, &bad).is_err());
    bad = row;
    bad.kind = BcKind::Traction;
    assert!(Bc::scenario_traction(Patch::Right, &bad).is_err());
}

#[test]
fn conflicting_shared_supports_and_duplicate_roles_refuse_without_mutation() {
    let mesh = Mesh2::quads(2.0, 1.0, 4, 2);
    let p = problem(&mesh);
    let zero = fixed([0.0; 3], LENGTH);
    let shifted = fixed([0.001, 0.0, 0.0], LENGTH);
    let rows = [Bc::Displacement { patch: Patch::Left, value: &zero },
        Bc::Displacement { patch: Patch::Top, value: &shifted }];
    assert!(matches!(p.solve_with_vector_boundaries(&FrameTree::new(), &rows, 0.0,
        &mut || ControlFlow::Continue(())), Err(Error::ConflictingDisplacement { .. })));
    let duplicate = [rows[0], rows[0]];
    assert!(matches!(p.solve_with_vector_boundaries(&FrameTree::new(), &duplicate, 0.0,
        &mut || ControlFlow::Continue(())), Err(Error::Invalid(_))));
    assert!(p.dirichlet.is_empty() && p.traction.is_empty());
}

#[test]
fn dimensions_out_of_plane_time_and_cancellation_are_not_silently_ignored() {
    let mesh = Mesh2::quads(2.0, 1.0, 4, 2);
    let p = problem(&mesh);
    let wrong_units = fixed([1.0, 0.0, 0.0], LENGTH);
    let out_of_plane = fixed([1.0, 0.0, 1e-20], STRESS);
    let frames = FrameTree::new();
    assert!(matches!(p.solve_with_vector_boundaries(&frames,
        &[Bc::Traction { patch: Patch::Right, value: &wrong_units }], 0.0,
        &mut || ControlFlow::Continue(())), Err(Error::Quantity { .. })));
    assert!(matches!(p.solve_with_vector_boundaries(&frames,
        &[Bc::Traction { patch: Patch::Right, value: &out_of_plane }], 0.0,
        &mut || ControlFlow::Continue(())), Err(Error::OutOfPlane { .. })));
    assert_eq!(p.solve_with_vector_boundaries(&frames, &[], 0.0,
        &mut || ControlFlow::Break(())), Err(Error::Cancelled));
    assert!(p.solve_with_vector_boundaries(&frames, &[], f64::NAN,
        &mut || ControlFlow::Continue(())).is_err());
}
