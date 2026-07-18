//! G0/G3 conformance for fs-ascent's migration to the fs-opt manifold runtime.
//!
//! Fixed fixtures distinguish packed point storage from optimizer-parameter
//! storage, exercise an SO(3) Study finite-difference step, and pin the public
//! Riemannian state dimensions for SO(3) and Stiefel. Constraint adapters are
//! required to fail closed outside their retained Euclidean point-coordinate
//! contract.
//!
//! This target makes no constrained-manifold, Stiefel transport-isometry,
//! superlinear-convergence, arbitrary-conditioning, cross-ISA, cancellation,
//! or performance claim.

#![deny(unsafe_code)]

use fs_ascent::{
    Packing, RiemannianLbfgs, StopReason, StopRule, Study, StudyError, retract as root_retract,
    tangent_project as root_tangent_project,
};
use fs_opt::{Manifold, NodeId, Problem, ProblemBuilder, Sense};
use fs_qty::Dims;
use std::panic::{AssertUnwindSafe, catch_unwind};

const D0: Dims = Dims([0, 0, 0, 0, 0, 0]);
const STIEFEL: Manifold = Manifold::Stiefel { n: 4, p: 2 };
const STIEFEL_POINT: [f64; 8] = [
    0.5, 0.5, 0.5, 0.5, // first column
    0.5, -0.5, 0.5, -0.5, // second column
];

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn sum_nodes(builder: &mut ProblemBuilder, nodes: &[NodeId]) -> NodeId {
    let mut total = nodes[0];
    for node in &nodes[1..] {
        total = builder.add(total, *node).expect("dimensionless scalar sum");
    }
    total
}

fn mixed_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let prefix = builder
        .var("prefix", Manifold::Rn { dim: 1 }, D0)
        .expect("leading Euclidean variable");
    let rotation = builder
        .var("rotation", Manifold::So3, D0)
        .expect("SO(3) variable");
    let frame = builder.var("frame", STIEFEL, D0).expect("Stiefel variable");
    let euclidean = builder
        .var("offset", Manifold::Rn { dim: 2 }, D0)
        .expect("Euclidean variable");
    let prefix_ref = builder.var_ref(prefix).expect("prefix reference");
    let rotation_ref = builder.var_ref(rotation).expect("rotation reference");
    let frame_ref = builder.var_ref(frame).expect("frame reference");
    let euclidean_ref = builder.var_ref(euclidean).expect("offset reference");
    let roots = [
        builder.component(prefix_ref, 0).expect("prefix component"),
        builder.component(rotation_ref, 0).expect("quaternion w"),
        builder.component(frame_ref, 0).expect("frame component"),
        builder
            .component(euclidean_ref, 0)
            .expect("offset component"),
    ];
    let objective = sum_nodes(&mut builder, &roots);
    builder
        .objective(objective, Sense::Minimize, 1.0)
        .expect("mixed objective");
    builder.finish()
}

fn mixed_point() -> Vec<f64> {
    let mut point = Vec::from([1.5, 1.0, 0.0, 0.0, 0.0]);
    point.extend(STIEFEL_POINT);
    point.extend([2.0, -1.0]);
    point
}

fn so3_antipodal_problem() -> Problem {
    let mut builder = ProblemBuilder::new();
    let rotation = builder
        .var("rotation", Manifold::So3, D0)
        .expect("SO(3) variable");
    let rotation_ref = builder.var_ref(rotation).expect("rotation reference");
    let quaternion_w = builder
        .component(rotation_ref, 0)
        .expect("quaternion w component");
    let quaternion_x = builder
        .component(rotation_ref, 1)
        .expect("quaternion x component");
    let objective = builder
        .mul(quaternion_w, quaternion_x)
        .expect("antipodally invariant quaternion product");
    builder
        .objective(objective, Sense::Minimize, 1.0)
        .expect("antipodally invariant quaternion objective");
    builder.finish()
}

fn euclidean_problem(dim: u32) -> Problem {
    let mut builder = ProblemBuilder::new();
    let variable = builder
        .var("euclidean", Manifold::Rn { dim }, D0)
        .expect("Euclidean variable");
    let variable_ref = builder.var_ref(variable).expect("Euclidean reference");
    let objective = builder
        .component(variable_ref, 0)
        .expect("Euclidean component");
    builder
        .objective(objective, Sense::Minimize, 1.0)
        .expect("Euclidean objective");
    builder.finish()
}

#[test]
fn g0_mixed_packing_keeps_point_and_parameter_offsets_distinct() {
    let problem = mixed_problem();
    let packing = Packing::new(&problem);
    assert_eq!(packing.dim, 15);
    assert_eq!(packing.point_dim(), 15);
    assert_eq!(packing.param_dim, 14);

    let point = mixed_point();
    let bindings = packing.unpack(&point);
    assert_eq!(bindings.len(), 4);
    assert_eq!(bits(&bindings[0]), bits(&point[..1]));
    assert_eq!(bits(&bindings[1]), bits(&point[1..5]));
    assert_eq!(bits(&bindings[2]), bits(&point[5..13]));
    assert_eq!(bits(&bindings[3]), bits(&point[13..]));

    let ambient = [
        -2.0, // leading Rn gradient
        0.0, 2.0, 0.0, 0.0, // SO(3) ambient quaternion gradient
        0.75, -0.5, 0.25, 1.0, -0.25, 0.5, 1.25, -0.75, // Stiefel ambient gradient
        3.0, -4.0, // Rn ambient gradient
    ];
    let parameter = packing.project(&point, &ambient);
    assert_eq!(parameter.len(), packing.param_dim);
    assert_eq!(parameter[0].to_bits(), (-2.0_f64).to_bits());
    assert_eq!(bits(&parameter[1..4]), bits(&[1.0, 0.0, 0.0]));
    assert_eq!(
        bits(&root_tangent_project(
            &Manifold::So3,
            &point[1..5],
            &ambient[1..5],
        )),
        bits(&parameter[1..4]),
    );
    assert_eq!(
        bits(&parameter[4..12]),
        bits(
            &STIEFEL
                .parameter_gradient(&point[5..13], &ambient[5..13])
                .expect("direct Stiefel gradient pullback")
        )
    );
    STIEFEL
        .validate_parameter_tangent(&point[5..13], &parameter[4..12])
        .expect("packed Stiefel gradient is tangent");
    assert_eq!(bits(&parameter[12..]), bits(&[3.0, -4.0]));

    let mut step = Vec::from([0.25, 0.2, -0.1, 0.3]);
    step.extend([0.01, -0.02, 0.0, 0.0, 0.015, 0.0, -0.01, 0.02]);
    step.extend([0.5, -0.25]);
    let landed = packing.retract(&point, &step);
    assert_eq!(landed.len(), packing.dim);
    assert_eq!(landed[0].to_bits(), 1.75_f64.to_bits());
    assert_eq!(
        bits(&landed[1..5]),
        bits(
            &Manifold::So3
                .retract(&point[1..5], &step[1..4])
                .expect("direct SO(3) authority landing")
        )
    );
    assert_eq!(
        bits(&root_retract(&Manifold::So3, &point[1..5], &step[1..4],)),
        bits(&landed[1..5]),
    );
    assert_eq!(
        bits(&landed[5..13]),
        bits(
            &STIEFEL
                .retract(&point[5..13], &step[4..12])
                .expect("direct Stiefel authority landing")
        )
    );
    STIEFEL
        .validate_parameter_tangent(&landed[5..13], &[0.0; 8])
        .expect("packed Stiefel landing remains a valid point");
    assert_eq!(bits(&landed[13..]), bits(&[2.5, -1.25]));

    let mut antipodal = point;
    for coordinate in &mut antipodal[1..5] {
        *coordinate = -*coordinate;
    }
    let zero_step = vec![0.0; packing.param_dim];
    let canonical = packing.retract(&antipodal, &zero_step);
    assert_eq!(bits(&canonical[1..5]), bits(&[1.0, 0.0, 0.0, 0.0]));
}

#[test]
fn g3_so3_study_uses_three_parameter_probes_for_four_point_lanes() {
    let problem = so3_antipodal_problem();
    let origin = [1.0, 0.0, 0.0, 0.0];
    let fd_h = 1.0e-4;
    let learning_rate = 0.25;
    let mut study = Study::new(&problem, &origin, fd_h, learning_rate);
    let report = study.run(&problem, &StopRule::GradNorm(0.0), 1);

    assert_eq!(report.reason, StopReason::IterationCap);
    assert_eq!(study.steps, 1);
    assert_eq!(study.x.len(), 4);
    assert_eq!(study.evals, 8, "one objective + 2*3 probes + one landing");
    assert_eq!(report.evals, study.evals);
    assert!(report.f < 0.0);
    assert!(study.x[1] < 0.0);
    Manifold::So3
        .validate_parameter_tangent(&study.x, &[0.0; 3])
        .expect("Study retains an authoritative quaternion point");

    let gradient: [f64; 3] = core::array::from_fn(|index| {
        let mut plus_step = [0.0; 3];
        plus_step[index] = fd_h;
        let mut minus_step = [0.0; 3];
        minus_step[index] = -fd_h;
        let plus = Manifold::So3
            .retract(&origin, &plus_step)
            .expect("positive direct probe");
        let minus = Manifold::So3
            .retract(&origin, &minus_step)
            .expect("negative direct probe");
        (plus[0] * plus[1] - minus[0] * minus[1]) / (2.0 * fd_h)
    });
    let expected_step = gradient.map(|value| -learning_rate * value);
    let expected = Manifold::So3
        .retract(&origin, &expected_step)
        .expect("direct expected Study landing");
    assert_eq!(bits(&study.x), bits(&expected));

    let mut repeat = Study::new(&problem, &origin, fd_h, learning_rate);
    let repeat_report = repeat.run(&problem, &StopRule::GradNorm(0.0), 1);
    assert_eq!(bits(&repeat.x), bits(&study.x));
    assert_eq!(repeat.evals, study.evals);
    assert_eq!(repeat_report.f.to_bits(), report.f.to_bits());
}

#[test]
fn g0_study_canonicalizes_so3_starts_and_types_invalid_point_refusals() {
    let problem = so3_antipodal_problem();
    let study = Study::new(&problem, &[-1.0, 0.0, 0.0, 0.0], 1.0e-4, 0.25);
    assert_eq!(bits(&study.x), bits(&[1.0, 0.0, 0.0, 0.0]));

    let refusal = Study::try_new(&problem, &[2.0, 0.0, 0.0, 0.0], 1.0e-4, 0.25)
        .expect_err("off-manifold start must fail before objective evaluation");
    assert!(matches!(
        refusal,
        StudyError::ManifoldPointInvalid { variable: 0, .. }
    ));
}

#[test]
fn g0_riemannian_state_uses_parameter_sized_gradients_for_so3_and_stiefel() {
    let mut so3_objective =
        |point: &[f64]| (point[0] * point[1], vec![point[1], point[0], 0.0, 0.0]);
    let so3 = RiemannianLbfgs::new(Manifold::So3, &[1.0, 0.0, 0.0, 0.0], 4, &mut so3_objective);
    assert_eq!(so3.x.len(), 4);
    assert_eq!(so3.g.len(), 3);
    Manifold::So3
        .validate_parameter_tangent(&so3.x, &so3.g)
        .expect("SO(3) state gradient uses body parameters");

    let ambient = [0.75, -0.5, 0.25, 1.0, -0.25, 0.5, 1.25, -0.75];
    let mut stiefel_objective = |point: &[f64]| {
        let value = point.iter().zip(ambient).map(|(x, g)| x * g).sum::<f64>();
        (value, ambient.to_vec())
    };
    let stiefel = RiemannianLbfgs::new(STIEFEL, &STIEFEL_POINT, 4, &mut stiefel_objective);
    assert_eq!(stiefel.x.len(), 8);
    assert_eq!(stiefel.g.len(), 8);
    STIEFEL
        .validate_parameter_tangent(&stiefel.x, &stiefel.g)
        .expect("Stiefel state gradient uses an embedded tangent parameter");
}

#[test]
fn g0_riemannian_runs_authoritative_so3_and_stiefel_curves() {
    let mut so3_objective =
        |point: &[f64]| (point[0] * point[1], vec![point[1], point[0], 0.0, 0.0]);
    let mut so3 = RiemannianLbfgs::new(Manifold::So3, &[1.0, 0.0, 0.0, 0.0], 4, &mut so3_objective);
    let so3_report = so3.run(&mut so3_objective, &StopRule::GradNorm(0.0), 1);
    assert_eq!(so3_report.reason, StopReason::IterationCap);
    assert_eq!(so3.iters, 1);
    assert_eq!(so3.x.len(), 4);
    assert_eq!(so3.g.len(), 3);
    assert!(so3.f < 0.0);

    let ambient = [0.75, -0.5, 0.25, 1.0, -0.25, 0.5, 1.25, -0.75];
    let tangent = STIEFEL
        .parameter_gradient(&STIEFEL_POINT, &ambient)
        .expect("Stiefel tangent direction");
    let target: Vec<f64> = STIEFEL_POINT
        .iter()
        .zip(&tangent)
        .map(|(point, direction)| point + 0.05 * direction)
        .collect();
    let mut stiefel_objective = |point: &[f64]| {
        let gradient: Vec<f64> = point
            .iter()
            .zip(&target)
            .map(|(x, goal)| x - goal)
            .collect();
        let value = 0.5 * gradient.iter().map(|value| value * value).sum::<f64>();
        (value, gradient)
    };
    let mut stiefel = RiemannianLbfgs::new(STIEFEL, &STIEFEL_POINT, 4, &mut stiefel_objective);
    let initial = stiefel.f;
    let stiefel_report = stiefel.run(&mut stiefel_objective, &StopRule::GradNorm(0.0), 1);
    assert_eq!(stiefel_report.reason, StopReason::IterationCap);
    assert_eq!(stiefel.iters, 1);
    assert!(stiefel.f < initial);
    STIEFEL
        .validate_parameter_tangent(&stiefel.x, &stiefel.g)
        .expect("accepted Stiefel state remains authoritative");
}

#[test]
fn g0_constraint_adapters_fail_closed_for_non_euclidean_packing() {
    let problem = mixed_problem();
    let packing = Packing::new(&problem);
    let refusal = catch_unwind(AssertUnwindSafe(|| {
        let _ = Study::constraint_adapters(&problem, &packing, 1.0e-6);
    }));
    assert!(refusal.is_err());

    let euclidean = euclidean_problem(15);
    let unrelated_euclidean_packing = Packing::new(&euclidean);
    let mismatch_refusal = catch_unwind(AssertUnwindSafe(|| {
        let _ = Study::constraint_adapters(&problem, &unrelated_euclidean_packing, 1.0e-6);
    }));
    assert!(mismatch_refusal.is_err());
}
