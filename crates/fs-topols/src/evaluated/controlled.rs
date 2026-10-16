//! Controlled final evaluations using the canonical residual-admitted solver.
use super::*;
use fs_cutfem::ControlledElasticitySolution;
use std::ops::ControlFlow;

/// Cooperative boundaries shared by final compliance and stress evaluations.
///
/// CG counts include all residual-correction passes. Assembly, area quadrature,
/// and one cell's stress probes remain indivisible; no hard deadline is implied.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DesignEvaluationStage {
    /// Before input validation or discretization allocation.
    Prepare,
    /// Before building the grid and assembling the canonical operator.
    Assemble,
    /// Cumulative CG iterations within this independent solve.
    Solve(usize),
    /// Before evaluating the numerical material area.
    Area,
    /// Before one cell's stress quadrature; includes cells later found empty.
    StressCell(usize),
    /// All requested quantities are complete, but no result has been returned.
    Publish,
}

// Retain the admitted field only for the duration of this final evaluation.
// It is never serialized or returned as a partially solved optimizer state.
pub(crate) struct EvaluationFields {
    pub grid: Quadtree,
    pub solution: ControlledElasticitySolution,
    pub lambda: f64,
    pub mu: f64,
    pub state: EvaluatedFinalState,
}

pub(crate) fn solve_fields_controlled<B>(
    phi: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
    poll_iters: usize,
    control: &mut impl FnMut(DesignEvaluationStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, EvaluationFields>, CutFemError> {
    if poll_iters == 0 {
        return Err(invalid_input("final-evaluation CG poll interval must be positive"));
    }
    if let ControlFlow::Break(reason) = control(DesignEvaluationStage::Prepare) {
        return Ok(ControlFlow::Break(reason));
    }
    let (support, material) = validate_inputs(phi, fixture, settings)?;
    let (lambda, mu) = material.lame();
    if let ControlFlow::Break(reason) = control(DesignEvaluationStage::Assemble) {
        return Ok(ControlFlow::Break(reason));
    }
    let grid = Quadtree::uniform(settings.level);
    let clamp = |x: f64, _y: f64| x < 1e-9;
    let traction = |_: f64, _: f64| [0.0, -fixture.load];
    let solver = CutElasticity {
        grid: &grid,
        sdf: phi,
        material: &material,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        stabilization_scaling: fs_cutfem::CutStabilizationScaling::LongitudinalModulus,
        quad_depth: 2,
        clamp: Some(&clamp),
        boundary_traction: None,
        traction_free_interface: true,
        solver_tol: SOLVER_TOL,
        solver_max_iters: SOLVER_MAX_ITERS,
    };
    let operator = solver.assemble_with_boundary_traction(
        &|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0],
        BoundaryTraction::EdgeBand { support, value: &traction },
    )?;
    let solution = match operator.solve_controlled(
        SOLVER_TOL, SOLVER_MAX_ITERS, poll_iters,
        |iters| control(DesignEvaluationStage::Solve(iters)),
    )? {
        ControlFlow::Continue(solution) => solution,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    if let ControlFlow::Break(reason) = control(DesignEvaluationStage::Area) {
        return Ok(ControlFlow::Break(reason));
    }
    let compliance = solution.compliance();
    let volume = material_volume(&grid, phi);
    if !(compliance.is_finite() && compliance >= 0.0 && volume.is_finite() && volume > 0.0) {
        return Err(invalid_input("final design evaluation produced invalid compliance or area"));
    }
    let state = EvaluatedFinalState {
        compliance, volume, snapshot: snapshot(phi),
        trajectory_compliance: None, trajectory_volume: None,
        compliance_delta: 0.0, volume_delta: 0.0,
    };
    Ok(ControlFlow::Continue(EvaluationFields { grid, solution, lambda, mu, state }))
}

/// Independently evaluate an exact field with interruption inside the CG solve.
///
/// Uses the same operator, material law, true-residual gate, area functional and
/// arithmetic ordering as [`evaluate_compliance_design`]. Cancellation returns
/// the caller's reason, never an approximate result or a numerical refusal.
/// The input is borrowed and no unfinished state is published.
///
/// # Errors
/// Refuses a zero polling interval, malformed input, or canonical PDE failure.
pub fn evaluate_compliance_design_controlled<B>(
    phi: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
    poll_iters: usize,
    mut control: impl FnMut(DesignEvaluationStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, EvaluatedFinalState>, CutFemError> {
    let fields = match solve_fields_controlled(phi, fixture, settings, poll_iters, &mut control)? {
        ControlFlow::Continue(fields) => fields,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    if let ControlFlow::Break(reason) = control(DesignEvaluationStage::Publish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue(fields.state))
}
