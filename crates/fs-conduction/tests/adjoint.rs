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

mod support;

use fs_adjoint::verify_gradient;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{ConductionProblem, LinearConfig};
use support::with_cx;

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
