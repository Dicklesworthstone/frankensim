//! Transactional, independently re-evaluated publication for topology optimization.
//!
//! The legacy descent trajectory is useful evolution evidence, but its iteration
//! rows historically couple a pre-evolution elasticity solve to a post-evolution
//! level-set snapshot. This module provides a publication boundary that always
//! re-solves the returned geometry with the canonical CutFEM elasticity operator.
//! The caller's geometry is replaced only after both evolution and that final
//! evaluation succeed.

use crate::GridSdf;
use crate::optimize::{Cantilever, OptimizeReport, OptimizeSettings, material_volume, optimize_compliance};
use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutFemError, DesignBoxEdge, EdgeBand,
    MAX_PLANE_STRAIN_STIFFNESS_RATIO, Quadtree,
};
use fs_material::IsotropicElastic;

const MATERIAL_STRAIN_LIMIT: f64 = 1.0;
const SOLVER_TOL: f64 = 1e-12;
const SOLVER_MAX_ITERS: usize = 60_000;

/// One authoritative evaluation of the geometry actually returned to the caller.
#[derive(Debug, Clone, Copy)]
pub struct EvaluatedFinalState {
    /// Canonical discrete external work `b^T u` on the returned geometry.
    pub compliance: f64,
    /// Cut-quadrature material area of the returned geometry.
    pub volume: f64,
    /// FNV-64 fingerprint of the returned nodal level-set bits.
    pub snapshot: u64,
    /// Compliance stored by the final legacy trajectory row, when one exists.
    pub trajectory_compliance: Option<f64>,
    /// Volume stored by the final legacy trajectory row, when one exists.
    pub trajectory_volume: Option<f64>,
    /// `compliance - trajectory_compliance`, or zero when there is no trajectory row.
    pub compliance_delta: f64,
    /// `volume - trajectory_volume`, or zero when there is no trajectory row.
    pub volume_delta: f64,
}

/// Evolution evidence plus an independently solved final publication state.
#[derive(Debug, Clone)]
pub struct EvaluatedOptimizeReport {
    /// Existing per-iteration evolution evidence. Treat the independently
    /// evaluated `final_state` as authoritative for the returned geometry.
    pub trajectory: OptimizeReport,
    /// Independently re-solved state for the exact geometry returned to the caller.
    pub final_state: EvaluatedFinalState,
}

fn invalid_input(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn snapshot(phi: &GridSdf) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for value in phi.nodes() {
        for byte in value.to_bits().to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

fn validate_inputs(
    phi: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
) -> Result<(EdgeBand, IsotropicElastic), CutFemError> {
    let expected = 1usize.checked_shl(settings.level).ok_or_else(|| {
        invalid_input("optimizer grid level exceeds the platform lattice address space")
    })?;
    if phi.n() != expected {
        return Err(invalid_input(format!(
            "SDF lattice has {} cells per side but level {} requires {expected}",
            phi.n(), settings.level
        )));
    }
    if phi.nodes().iter().any(|value| !value.is_finite()) {
        return Err(invalid_input("optimizer level set must contain only finite nodal values"));
    }
    if !(fixture.load.is_finite() && fixture.load > 0.0) {
        return Err(invalid_input("cantilever load must be finite and strictly positive"));
    }
    let support = EdgeBand::new(
        DesignBoxEdge::Right,
        0.5 - fixture.band,
        0.5 + fixture.band,
    )
    .map_err(|error| invalid_input(format!("invalid cantilever load band: {error}")))?;
    if !(settings.youngs.is_finite() && settings.youngs > 0.0) {
        return Err(invalid_input("Young's modulus must be finite and strictly positive"));
    }
    if !(settings.poisson.is_finite() && settings.poisson > -1.0 && settings.poisson < 0.5) {
        return Err(invalid_input("Poisson ratio must lie in (-1, 0.5)"));
    }
    let material = IsotropicElastic::new(settings.youngs, settings.poisson, MATERIAL_STRAIN_LIMIT)
        .map_err(|error| invalid_input(format!("material card refused final evaluation: {error}")))?;
    let (lambda, mu) = material.lame();
    let ratio = (lambda + 2.0 * mu) / mu;
    if !(lambda.is_finite()
        && mu.is_finite()
        && mu > 0.0
        && ratio.is_finite()
        && ratio <= MAX_PLANE_STRAIN_STIFFNESS_RATIO)
    {
        return Err(invalid_input(format!(
            "plane-strain stiffness ratio {ratio} exceeds the admitted final-evaluation regime"
        )));
    }
    Ok((support, material))
}

/// Independently solve compliance and material area for one exact level set.
///
/// This does not evolve the geometry and makes no optimality, convergence, or
/// physical-validation claim.
///
/// # Errors
/// Returns a typed elasticity refusal for invalid geometry/material/load inputs or
/// whenever the canonical CutFEM solve refuses the design.
pub fn evaluate_compliance_design(
    phi: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
) -> Result<EvaluatedFinalState, CutFemError> {
    let (support, material) = validate_inputs(phi, fixture, settings)?;
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
    let solution = solver.solve_with_boundary_traction(
        &|_, _| [0.0, 0.0],
        &|_, _| [0.0, 0.0],
        BoundaryTraction::EdgeBand { support, value: &traction },
    )?;
    let compliance = solution.compliance();
    let volume = material_volume(&grid, phi);
    if !(compliance.is_finite() && compliance >= 0.0 && volume.is_finite() && volume > 0.0) {
        return Err(invalid_input("final design evaluation produced invalid compliance or area"));
    }
    Ok(EvaluatedFinalState {
        compliance,
        volume,
        snapshot: snapshot(phi),
        trajectory_compliance: None,
        trajectory_volume: None,
        compliance_delta: 0.0,
        volume_delta: 0.0,
    })
}

/// Run the existing level-set optimizer transactionally and independently
/// re-evaluate the exact final geometry before publishing it.
///
/// The input geometry is untouched if admission, evolution, or final evaluation
/// fails. This is the recommended publication boundary for consumers that need
/// objective/area values bound to the returned geometry.
///
/// # Errors
/// Propagates typed admission, evolution, and final CutFEM evaluation refusals.
pub fn optimize_compliance_evaluated(
    phi: &mut GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
) -> Result<EvaluatedOptimizeReport, CutFemError> {
    validate_inputs(phi, fixture, settings)?;
    let mut candidate = phi.clone();
    let trajectory = optimize_compliance(&mut candidate, fixture, settings)?;
    let mut final_state = evaluate_compliance_design(&candidate, fixture, settings)?;
    final_state.trajectory_compliance = trajectory.compliance.last().copied();
    final_state.trajectory_volume = trajectory.volume.last().copied();
    final_state.compliance_delta = final_state
        .trajectory_compliance
        .map_or(0.0, |reported| final_state.compliance - reported);
    final_state.volume_delta = final_state
        .trajectory_volume
        .map_or(0.0, |reported| final_state.volume - reported);
    *phi = candidate;
    Ok(EvaluatedOptimizeReport { trajectory, final_state })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn beam(level: u32) -> GridSdf {
        GridSdf::from_fn(1usize << level, &|_, y| (y - 0.5).abs() - 0.35)
    }

    #[test]
    fn dimension_mismatch_refuses_without_mutation() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let settings = OptimizeSettings { level: 4, iterations: 1, ..OptimizeSettings::default() };
        assert!(optimize_compliance_evaluated(
            &mut phi,
            Cantilever { load: 1.0, band: 0.125 },
            settings,
        ).is_err());
        assert_eq!(phi.nodes(), before.as_slice());
    }

    #[test]
    fn zero_iteration_publication_still_evaluates_exact_returned_geometry() {
        let mut phi = beam(3);
        let settings = OptimizeSettings { level: 3, iterations: 0, ..OptimizeSettings::default() };
        let report = optimize_compliance_evaluated(
            &mut phi,
            Cantilever { load: 1.0, band: 0.125 },
            settings,
        ).expect("final geometry should solve");
        assert!(report.trajectory.rows.is_empty());
        assert_eq!(report.final_state.snapshot, snapshot(&phi));
        assert!(report.final_state.compliance.is_finite());
        assert!(report.final_state.volume > 0.0);
        assert_eq!(report.final_state.trajectory_compliance, None);
        assert_eq!(report.final_state.trajectory_volume, None);
    }

    #[test]
    fn bad_material_refuses_without_mutating_caller_geometry() {
        let mut phi = beam(3);
        let before = phi.nodes().to_vec();
        let settings = OptimizeSettings {
            level: 3,
            iterations: 1,
            youngs: f64::NAN,
            ..OptimizeSettings::default()
        };
        assert!(optimize_compliance_evaluated(
            &mut phi,
            Cantilever { load: 1.0, band: 0.125 },
            settings,
        ).is_err());
        assert_eq!(phi.nodes(), before.as_slice());
    }
}
