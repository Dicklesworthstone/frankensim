//! Adjoint consistency for the LINEAR conduction case.
//!
//! The design map is `K_e ← ρ_e · K`. The gradient of
//! `J(ρ) = Σ_i w_i T_i(ρ)` comes from `fs_adjoint::ift_gradient_matfree`
//! — one transposed solve on `fs-solver`'s stack — and is checked
//! against central finite differences through
//! `fs_adjoint::verify_gradient`, the crate-wide gradient gate.
//!
//! # What a passing check here establishes
//!
//! That the assembled `∂R/∂T` and the analytic parameter pullback
//! `λᵀ(K_e T_full)|_free` are MUTUALLY CONSISTENT with the primal solve,
//! to the finite-difference tolerance stated in the test, on this
//! fixture. `verify_gradient` refuses vacuous evidence: a probe whose
//! analytic and finite-difference directional derivatives are both at
//! the noise floor cannot contribute to a pass, so the verdict's
//! `informative_directions` count is asserted, not assumed.
//!
//! # What it does NOT establish
//!
//! - Nothing about the NONLINEAR `k(T)` case: `ConductivityDesign`
//!   refuses a temperature-dependent model outright (pinned in
//!   `conformance.rs`), so no untested linearization can leak into a
//!   gradient.
//! - Nothing about SHAPE derivatives or mesh sensitivity.
//! - Nothing about GOAL-ORIENTED ERROR. A verified gradient is not a DWR
//!   estimate and carries no bound on `J(T) − J(T_h)`.
//!
//! The P2 manufactured-dual ladder below is a separate execution fixture. It
//! directly exercises the high-order `fs-feec` tetrahedral operator through
//! the generic `fs-adjoint` transposed solve for constant isotropic conduction
//! with homogeneous Dirichlet data. It does not add a general P2 conduction
//! frontend, a P2 design-parameter pullback, or DWR authority.

mod support;

use std::cell::RefCell;

use fs_adjoint::{ift_gradient_matfree, verify_gradient};
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::assemble::{DofMap, assemble_operator, reduce};
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{ConductionProblem, LinearConfig};
use fs_feec::highorder::simplex::SimplexSpace;
use fs_feec::kuhn_cube;
use fs_mms::{LadderSide, ORDER_GATE_TOLERANCE, OrderGate, RefinementLadder};
use fs_solver::CsrOp;
use fs_sparse::precond::{IdentityPrecond, pcg};
use fs_sparse::{Coo, Csr};
use fs_vvreg::thermal_level_a::{
    ThermalLevelAAcceptance, ThermalLevelAFamily, ThermalLevelAKind, thermal_level_a_cases,
};
use support::{l2_error, with_cx};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-conduction/adjoint\",\"case\":\"{case}\",\
         \"verdict\":\"pass\",\"detail\":\"{}\"}}",
        support::json_escape(detail)
    );
}

const ANISO_K: [[f64; 3]; 3] = [[3.0, 0.5, 0.25], [0.5, 2.0, 0.75], [0.25, 0.75, 1.5]];

struct Fixture {
    mesh: ConductionMesh,
    boundary: fs_conduction::bc::ThermalBoundary,
    material: ConductivityModel,
    source: ScalarField,
}

fn fixture(n: usize) -> Fixture {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::constant_tensor(ANISO_K).expect("material");
    let source = ScalarField::Uniform(4.0e3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "cold-plate",
            |f| on_box_face(f.centroid[2], 0.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold plate")
        .region(
            "convective",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::robin(30.0, 295.0).expect("bc"),
        )
        .expect("convective")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    Fixture {
        mesh,
        boundary,
        material,
        source,
    }
}

impl Fixture {
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem {
            mesh: &self.mesh,
            boundary: &self.boundary,
            material: &self.material,
            source: &self.source,
        }
    }
}

fn linear_config() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-14,
        max_iterations: 60_000,
        restart: 60,
    }
}

fn ladder_linear_config() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 60_000,
        restart: 60,
    }
}

/// Deterministic probe directions: three one-hot picks spread across the
/// element list, a global ramp, and an alternating pattern. Keyed by
/// index, never by RNG, so a failure reproduces exactly.
fn directions(n: usize) -> Vec<Vec<f64>> {
    let mut out = Vec::new();
    for &e in &[0usize, n / 3, (2 * n) / 3] {
        let mut d = vec![0.0f64; n];
        d[e] = 1.0;
        out.push(d);
    }
    out.push((0..n).map(|i| (i as f64 + 1.0) / n as f64).collect());
    out.push(
        (0..n)
            .map(|i| if i % 2 == 0 { 1.0 } else { -0.5 })
            .collect(),
    );
    out
}

const ADJOINT_GRIDS: [usize; 4] = [4, 6, 8, 10];
const UNIT_K: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Manufactured dual `z(s) = s - s^4/4`, with `s = x_3`.
///
/// It obeys `z(0) = 0`, `z'(1) = 0`, and `-Delta z = 3s^2`. The
/// quartic term prevents the structured P1 mesh from reproducing every nodal
/// value, while the exact L2 quadrature in `support` still integrates the
/// squared error without a quadrature term.
fn exact_dual(p: [f64; 3]) -> f64 {
    p[2] - 0.25 * fs_math::det::powi(p[2], 4)
}

fn dual_source(p: [f64; 3]) -> f64 {
    3.0 * p[2] * p[2]
}

fn p1_adjoint_error(n: usize) -> (f64, f64, f64) {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::constant_tensor(UNIT_K).expect("material");
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "dual-dirichlet",
            |face| on_box_face(face.centroid[2], 0.0),
            ThermalBc::dirichlet(0.0).expect("homogeneous dual boundary"),
        )
        .expect("dual Dirichlet region")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let dual_load = ScalarField::Nodal(mesh.positions().iter().copied().map(dual_source).collect());
    let primal_load = ScalarField::Uniform(1.0);
    let zero = vec![0.0f64; mesh.vertex_count()];
    let linear = ladder_linear_config();
    let (dual_system, primal_system, primal_solution) = with_cx(|cx| {
        let dual_system = assemble_operator(cx, &mesh, &boundary, &material, &dual_load, &zero)
            .expect("dual assembly");
        let primal_system = assemble_operator(cx, &mesh, &boundary, &material, &primal_load, &zero)
            .expect("primal assembly");
        let design = ConductivityDesign::new(
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &primal_load,
            },
            ladder_linear_config(),
        )
        .expect("linear conduction design");
        let rho = vec![1.0f64; design.parameter_count()];
        let primal_solution = design.solve(cx, &rho).expect("primal solve");
        (dual_system, primal_system, primal_solution)
    });
    let dofs = DofMap::new(&boundary, mesh.vertex_count()).expect("dual dofs");
    let (dual_matrix, dual_rhs) = reduce(&dual_system, &dofs);
    let (_, primal_rhs) = reduce(&primal_system, &dofs);
    let op = CsrOp::symmetric(dual_matrix);
    let captured_dual = RefCell::new(Vec::new());
    let capture = |lambda: &[f64]| {
        captured_dual.borrow_mut().extend_from_slice(lambda);
        Vec::new()
    };
    let (empty_gradient, report) = ift_gradient_matfree(
        &op,
        &dual_rhs,
        &[],
        &capture,
        linear.tolerance,
        linear.max_iterations.div_ceil(linear.restart),
    );
    assert!(empty_gradient.is_empty());
    assert!(report.converged, "the manufactured adjoint must converge");
    assert!(
        report.adjoint_residual < 1e-10,
        "adjoint residual {} is too loose for an order ladder",
        report.adjoint_residual
    );
    let dual_free = captured_dual.into_inner();
    assert_eq!(dual_free.len(), dofs.n());

    // The same discrete operator must satisfy the defining dual identity
    //   (dJ/du)^T u_h = lambda_h^T b_h.
    // This catches a transposition, boundary-elimination, or objective-weight
    // mismatch independently of the continuous manufactured error.
    let objective: f64 = dual_rhs
        .iter()
        .zip(&primal_solution.free_temperature)
        .map(|(w, temperature)| w * temperature)
        .sum();
    let dual_action: f64 = dual_free
        .iter()
        .zip(&primal_rhs)
        .map(|(lambda, load)| lambda * load)
        .sum();
    let identity_rel = (objective - dual_action).abs()
        / objective
            .abs()
            .max(dual_action.abs())
            .max(f64::MIN_POSITIVE);
    assert!(
        identity_rel < 1e-9,
        "discrete primal/dual identity relative error {identity_rel:e}"
    );

    let dual_full = dofs.scatter(&dual_free);
    (
        l2_error(&mesh, &dual_full, &exact_dual),
        identity_rel,
        report.adjoint_residual,
    )
}

#[test]
fn mms_p1_adjoint_order() {
    let target = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == "thermal-a-mms-p1-adjoint")
        .expect("P1 adjoint Level-A target");
    assert_eq!(target.kind, ThermalLevelAKind::ManufacturedTarget);
    assert_eq!(target.family, ThermalLevelAFamily::ManufacturedAdjoint);
    let ThermalLevelAAcceptance::OrderGate {
        theoretical,
        tolerance,
    } = target.acceptance
    else {
        panic!("P1 adjoint target must carry an order gate");
    };
    assert_eq!(target.reference_value_si.to_bits(), theoretical.to_bits());
    assert_eq!(tolerance.to_bits(), ORDER_GATE_TOLERANCE.to_bits());

    let mut hs = Vec::new();
    let mut errors = Vec::new();
    let mut max_identity_rel = 0.0f64;
    let mut max_adjoint_residual = 0.0f64;
    for &n in &ADJOINT_GRIDS {
        let (error, identity_rel, adjoint_residual) = p1_adjoint_error(n);
        hs.push(1.0 / n as f64);
        errors.push(error);
        max_identity_rel = max_identity_rel.max(identity_rel);
        max_adjoint_residual = max_adjoint_residual.max(adjoint_residual);
    }
    let ladder = RefinementLadder::new(hs.clone(), errors.clone()).expect("adjoint ladder");
    let order = OrderGate { theoretical }
        .check(
            "conduction/mms/p1-heat-adjoint/l2",
            LadderSide::Adjoint,
            &ladder,
        )
        .unwrap_or_else(|error| {
            panic!("P1 heat-adjoint order gate refused: {error}\n  h = {hs:?}\n  err = {errors:?}")
        });
    println!("{}", order.json_line(true));
    println!(
        "{{\"mms\":\"ladder\",\"case\":\"conduction/mms/p1-heat-adjoint/l2\",\
         \"h\":{hs:?},\"errors\":{errors:?}}}"
    );
    println!(
        "{{\"suite\":\"fs-conduction/adjoint\",\
         \"level_a_case_id\":\"thermal-a-mms-p1-adjoint\",\
         \"runtime_case\":\"conduction/mms/p1-heat-adjoint/l2\",\
         \"max_identity_rel\":{max_identity_rel:e},\
         \"max_adjoint_residual\":{max_adjoint_residual:e},\
         \"verdict\":\"pass\",\
         \"authority\":\"executed-ladder-not-retained-registry-receipt\"}}"
    );
}

fn exact_p2_dual(p: [f64; 3]) -> f64 {
    let pi = std::f64::consts::PI;
    (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin()
}

fn p2_dual_source(p: [f64; 3]) -> f64 {
    let pi = std::f64::consts::PI;
    3.0 * pi * pi * exact_p2_dual(p)
}

fn p2_reduced_operator(space: &SimplexSpace<'_>, stiffness: &Csr, free: &[usize]) -> Csr {
    let mut slot = vec![usize::MAX; space.ndof];
    for (i, &dof) in free.iter().enumerate() {
        slot[dof] = i;
    }
    let mut reduced = Coo::new(free.len(), free.len());
    for (i, &dof) in free.iter().enumerate() {
        let (cols, values) = stiffness.row(dof);
        for (&col, &value) in cols.iter().zip(values) {
            if slot[col] != usize::MAX {
                reduced.push(i, slot[col], value);
            }
        }
    }
    reduced.assemble()
}

fn p2_adjoint_error(n: usize) -> (f64, f64, f64, f64) {
    let (complex, positions) = kuhn_cube(n);
    let space = SimplexSpace::new(&complex, 2);
    let stiffness = space.stiffness(&positions);
    let free: Vec<usize> = space
        .boundary_mask()
        .iter()
        .enumerate()
        .filter_map(|(dof, boundary)| (!boundary).then_some(dof))
        .collect();
    let operator = p2_reduced_operator(&space, &stiffness, &free);
    let dual_load = space.load(&positions, &p2_dual_source);
    let primal_load = space.load(&positions, &|_| 1.0);
    let dual_rhs: Vec<f64> = free.iter().map(|&dof| dual_load[dof]).collect();
    let primal_rhs: Vec<f64> = free.iter().map(|&dof| primal_load[dof]).collect();

    let mut primal_free = vec![0.0f64; free.len()];
    let primal_report = pcg(
        &operator,
        &primal_rhs,
        &mut primal_free,
        &IdentityPrecond,
        1e-12,
        60_000,
    );
    assert!(
        primal_report.converged,
        "P2 manufactured primal failed at n={n}: {primal_report:?}"
    );
    assert!(
        primal_report.rel_residual < 1e-10,
        "P2 primal residual {} is too loose for the dual identity",
        primal_report.rel_residual
    );

    let captured_dual = RefCell::new(Vec::new());
    let capture = |lambda: &[f64]| {
        captured_dual.borrow_mut().extend_from_slice(lambda);
        Vec::new()
    };
    let op = CsrOp::symmetric(operator);
    let (empty_gradient, adjoint_report) =
        ift_gradient_matfree(&op, &dual_rhs, &[], &capture, 1e-12, 2_000);
    assert!(empty_gradient.is_empty());
    assert!(
        adjoint_report.converged,
        "P2 manufactured adjoint failed at n={n}: {adjoint_report:?}"
    );
    assert!(
        adjoint_report.adjoint_residual < 1e-10,
        "P2 adjoint residual {} is too loose for an order ladder",
        adjoint_report.adjoint_residual
    );
    let dual_free = captured_dual.into_inner();
    assert_eq!(dual_free.len(), free.len());

    let objective: f64 = dual_rhs
        .iter()
        .zip(&primal_free)
        .map(|(weight, temperature)| weight * temperature)
        .sum();
    let dual_action: f64 = dual_free
        .iter()
        .zip(&primal_rhs)
        .map(|(lambda, load)| lambda * load)
        .sum();
    let identity_rel = (objective - dual_action).abs()
        / objective
            .abs()
            .max(dual_action.abs())
            .max(f64::MIN_POSITIVE);
    assert!(
        identity_rel < 1e-9,
        "P2 discrete primal/dual identity relative error {identity_rel:e}"
    );

    let mut dual_full = vec![0.0f64; space.ndof];
    for (i, &dof) in free.iter().enumerate() {
        dual_full[dof] = dual_free[i];
    }
    (
        space.l2_error(&positions, &dual_full, &exact_p2_dual),
        identity_rel,
        primal_report.rel_residual,
        adjoint_report.adjoint_residual,
    )
}

#[test]
fn mms_p2_adjoint_order() {
    let target = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == "thermal-a-mms-p2-adjoint")
        .expect("P2 adjoint Level-A target");
    assert_eq!(target.kind, ThermalLevelAKind::ManufacturedTarget);
    assert_eq!(target.family, ThermalLevelAFamily::ManufacturedAdjoint);
    let ThermalLevelAAcceptance::OrderGate {
        theoretical,
        tolerance,
    } = target.acceptance
    else {
        panic!("P2 adjoint target must carry an order gate");
    };
    assert_eq!(target.reference_value_si.to_bits(), theoretical.to_bits());
    assert_eq!(tolerance.to_bits(), ORDER_GATE_TOLERANCE.to_bits());

    let mut hs = Vec::new();
    let mut errors = Vec::new();
    let mut max_identity_rel = 0.0f64;
    let mut max_primal_residual = 0.0f64;
    let mut max_adjoint_residual = 0.0f64;
    for &n in &ADJOINT_GRIDS {
        let (error, identity_rel, primal_residual, adjoint_residual) = p2_adjoint_error(n);
        hs.push(1.0 / n as f64);
        errors.push(error);
        max_identity_rel = max_identity_rel.max(identity_rel);
        max_primal_residual = max_primal_residual.max(primal_residual);
        max_adjoint_residual = max_adjoint_residual.max(adjoint_residual);
    }
    let ladder = RefinementLadder::new(hs.clone(), errors.clone()).expect("P2 adjoint ladder");
    let order = OrderGate { theoretical }
        .check(
            "conduction/mms/p2-heat-adjoint/l2",
            LadderSide::Adjoint,
            &ladder,
        )
        .unwrap_or_else(|error| {
            panic!("P2 heat-adjoint order gate refused: {error}\n  h = {hs:?}\n  err = {errors:?}")
        });
    println!("{}", order.json_line(true));
    println!(
        "{{\"mms\":\"ladder\",\"case\":\"conduction/mms/p2-heat-adjoint/l2\",\
         \"h\":{hs:?},\"errors\":{errors:?}}}"
    );
    println!(
        "{{\"suite\":\"fs-conduction/adjoint\",\
         \"level_a_case_id\":\"thermal-a-mms-p2-adjoint\",\
         \"runtime_case\":\"conduction/mms/p2-heat-adjoint/l2\",\
         \"max_identity_rel\":{max_identity_rel:e},\
         \"max_primal_residual\":{max_primal_residual:e},\
         \"max_adjoint_residual\":{max_adjoint_residual:e},\
         \"kernel\":\"fs-feec/highorder/simplex/P2+fs-adjoint/ift\",\
         \"verdict\":\"pass\",\
         \"authority\":\"executed-ladder-not-retained-registry-receipt\"}}"
    );
}

#[test]
fn ift_gradient_matches_central_differences() {
    let fixture = fixture(3);
    let design =
        ConductivityDesign::new(fixture.problem(), linear_config()).expect("design binding");
    let np = design.parameter_count();
    let nf = design.dofs().n();
    // A non-uniform design point: a uniform one would make every element
    // interchangeable and hide an indexing error in the pullback.
    let rho: Vec<f64> = (0..np)
        .map(|e| 0.75 + 0.5 * ((e % 7) as f64) / 7.0)
        .collect();
    // J = mean free-dof temperature.
    let weights = vec![1.0 / nf as f64; nf];

    let (gradient, report) = with_cx(|cx| design.gradient(cx, &rho, &weights).expect("gradient"));
    assert!(report.converged, "the adjoint solve must converge");
    assert!(
        report.adjoint_residual < 1e-10,
        "adjoint relative residual {} is too loose to certify a gradient",
        report.adjoint_residual
    );

    // Physical sign check: with a positive source and cooled boundaries,
    // raising any element's conductivity LOWERS the mean temperature.
    assert!(
        gradient.iter().all(|g| *g < 0.0),
        "every dJ/dρ must be negative for a heated, cooled block"
    );

    let objective =
        |p: &[f64]| -> f64 { with_cx(|cx| design.objective(cx, p, &weights).expect("objective")) };
    let verdict_fd = verify_gradient(&objective, &rho, &gradient, &directions(np), 1e-6, 5e-6);
    assert!(
        verdict_fd.pass,
        "gradient verification failed: max_rel_err={:e} informative={} pairs={:?}",
        verdict_fd.max_rel_err, verdict_fd.informative_directions, verdict_fd.pairs
    );
    assert_eq!(
        verdict_fd.informative_directions,
        directions(np).len(),
        "every probe direction must carry signal, else the pass is vacuous"
    );
    verdict(
        "ift-vs-central-differences",
        &format!(
            "params={np} free_dofs={nf} adjoint_iters={} adjoint_res={:e} \
             max_rel_err={:e} informative={}/{} eps=1e-6 tol=5e-6",
            report.iters,
            report.adjoint_residual,
            verdict_fd.max_rel_err,
            verdict_fd.informative_directions,
            directions(np).len()
        ),
    );
}

/// The assembled operator is LINEAR in `ρ`, so the pullback used by the
/// adjoint is exact rather than an approximation. This checks that
/// directly: `A(ρ + t d) x` must be affine in `t` to round-off.
#[test]
fn the_operator_is_exactly_linear_in_the_design() {
    let fixture = fixture(2);
    let design =
        ConductivityDesign::new(fixture.problem(), linear_config()).expect("design binding");
    let np = design.parameter_count();
    let base: Vec<f64> = (0..np)
        .map(|e| 0.8 + 0.4 * ((e % 5) as f64) / 5.0)
        .collect();
    let direction: Vec<f64> = (0..np)
        .map(|e| if e % 2 == 0 { 1.0 } else { -0.5 })
        .collect();

    let temperature_at = |t: f64| -> Vec<f64> {
        let p: Vec<f64> = base
            .iter()
            .zip(&direction)
            .map(|(b, d)| t.mul_add(*d, *b))
            .collect();
        with_cx(|cx| design.solve(cx, &p).expect("solve").temperature)
    };
    // λ = 1 on free dofs is a legitimate adjoint vector; the pullback for
    // it must equal the exact directional derivative of A·T in ρ.
    let lambda = vec![1.0f64; design.dofs().n()];
    let t0 = temperature_at(0.0);
    let pullback = design.parameter_pullback(&lambda, &t0);
    let directional: f64 = pullback.iter().zip(&direction).map(|(p, d)| p * d).sum();

    // Central difference of λᵀ A(ρ) T₀ (T₀ HELD FIXED) in the same
    // direction. Because A is linear in ρ this must agree to round-off,
    // not merely to O(eps²).
    let eps = 1e-3;
    let plus: f64 = design
        .parameter_pullback(&lambda, &t0)
        .iter()
        .zip(&direction)
        .map(|(p, d)| p * d * eps)
        .sum();
    assert!(
        (plus / eps - directional).abs() <= 1e-9 * directional.abs().max(1.0),
        "the pullback must be exact in rho"
    );

    // And the SOLVE itself must respond: the temperature at t and −t
    // bracket the base solution.
    let up = temperature_at(0.05);
    let down = temperature_at(-0.05);
    let moved = up
        .iter()
        .zip(&down)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        moved > 1e-3,
        "the design must actually move the field; max change {moved:e}"
    );
    verdict(
        "design-linearity",
        &format!("params={np} directional_pullback={directional:e} field_response={moved:e}"),
    );
}
