//! Final-design goal-error assessment using the canonical vector DWR owner.
//! The level set is unchanged on both grids. This is an estimate, never a
//! certified continuum-error bound; absolute indicator mass is for marking.

use super::*;
use std::convert::Infallible;
use std::ops::ControlFlow;

/// Two solves on the exact declared geometry, with its identity and area.
#[derive(Debug, Clone)]
pub struct ComplianceDwrAssessment {
    /// Fingerprint of the original nodal level-set bits.
    pub snapshot: u64,
    /// Material area on the original grid.
    pub volume: f64,
    /// Original grid level; the enriched solve uses one extra level.
    pub level: u32,
    /// Actual coarse/enriched solves and the signed residual decomposition.
    pub estimate: fs_dwr::ElasticityDwrEstimate,
}

/// Checkpoints bracketing an indivisible two-solve DWR assessment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComplianceDwrStage {
    BeforeEstimate,
    BeforePublish,
}

/// Assess the final material design without advancing or projecting it.
///
/// The canonical load, clamp, material, quadrature and stabilization settings
/// match [`evaluate_compliance_design`]. Each solve has the same 60,000-iteration
/// cap and 1e-12 true-Euclidean residual gate. The enriched grid has one extra
/// level but samples the original bilinear field, not an evolved geometry.
///
/// # Errors
/// Returns typed admission/physics refusals; assessment admits levels 1..=5
/// before allocating either grid, bounding enrichment at level 6.
pub fn assess_compliance_dwr(
    phi: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
) -> Result<ComplianceDwrAssessment, CutFemError> {
    match assess_compliance_dwr_controlled(phi, fixture, settings, |_| {
        ControlFlow::<Infallible>::Continue(())
    })? {
        ControlFlow::Continue(assessment) => Ok(assessment),
        ControlFlow::Break(never) => match never {},
    }
}

/// Phase-boundary cancellation for final-design assessment.
///
/// A break returns no partial assessment and never mutates the input. DWR
/// assembly, both linear solves and residual integration remain indivisible:
/// this API does not promise intra-kernel cancellation or a wall-time bound.
///
/// # Errors
/// Propagates the same admission and physics refusals as [`assess_compliance_dwr`].
pub fn assess_compliance_dwr_controlled<B>(
    phi: &GridSdf,
    fixture: Cantilever,
    settings: OptimizeSettings,
    mut control: impl FnMut(ComplianceDwrStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, ComplianceDwrAssessment>, CutFemError> {
    if !(1..=5).contains(&settings.level) {
        return Err(invalid_input(
            "final elasticity DWR assessment admits levels 1..=5, enriched through level 6",
        ));
    }
    let (support, material) = validate_inputs(phi, fixture, settings)?;
    if let ControlFlow::Break(reason) = control(ComplianceDwrStage::BeforeEstimate) {
        return Ok(ControlFlow::Break(reason));
    }
    let grid = Quadtree::uniform(settings.level);
    let clamp = |x: f64, _y: f64| x < 1e-9;
    let traction = |_: f64, _: f64| [0.0, -fixture.load];
    let problem = CutElasticity {
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
    let estimate = fs_dwr::estimate_elasticity_compliance_with_boundary_traction(
        &problem,
        &|_, _| [0.0, 0.0],
        &|_, _| [0.0, 0.0],
        BoundaryTraction::EdgeBand { support, value: &traction },
    )?;
    let volume = material_volume(&grid, phi);
    if !(volume.is_finite() && volume > 0.0) {
        return Err(invalid_input("DWR assessment produced invalid material area"));
    }
    if let ControlFlow::Break(reason) = control(ComplianceDwrStage::BeforePublish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue(ComplianceDwrAssessment {
        snapshot: snapshot(phi),
        volume,
        level: settings.level,
        estimate,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings() -> OptimizeSettings {
        OptimizeSettings { level: 3, iterations: 1, ..OptimizeSettings::default() }
    }

    fn beam(width: f64) -> GridSdf {
        GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - width)
    }

    fn fixture(load: f64) -> Cantilever { Cantilever { load, band: 0.125 } }

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 1e-8 * a.abs().max(b.abs()).max(1e-8),
            "{a:e} != {b:e}",
        );
    }

    #[test]
    fn dwr_binds_exact_design_and_pde_with_quadratic_load_scaling() {
        let phi = beam(0.35);
        let before = phi.nodes().to_vec();
        let assessed = assess_compliance_dwr(&phi, fixture(1.0), settings()).unwrap();
        let canonical = evaluate_compliance_design(&phi, fixture(1.0), settings()).unwrap();
        assert_eq!(assessed.snapshot, canonical.snapshot);
        assert_eq!(assessed.volume.to_bits(), canonical.volume.to_bits());
        assert_eq!(assessed.estimate.j_primal.to_bits(), canonical.compliance.to_bits());
        assert_eq!(phi.nodes(), before);
        let e = &assessed.estimate;
        close(e.eta_signed, e.terms.total());
        assert!(e.eta_abs + 1e-14 >= e.eta_signed.abs());
        assert!(e.enriched_dofs > e.dofs);
        assert!(e.primal_relative_residual < SOLVER_TOL);
        assert!(e.enriched_relative_residual < SOLVER_TOL);
        assert!(e.primal_iterations <= SOLVER_MAX_ITERS);
        assert!(e.enriched_iterations <= SOLVER_MAX_ITERS);
        let doubled = assess_compliance_dwr(&phi, fixture(2.0), settings()).unwrap();
        for (a, b) in [
            (e.j_primal, doubled.estimate.j_primal),
            (e.j_enriched, doubled.estimate.j_enriched),
            (e.eta_signed, doubled.estimate.eta_signed),
            (e.eta_abs, doubled.estimate.eta_abs),
        ] {
            close(4.0 * a, b);
        }
        let thinner = assess_compliance_dwr(&beam(0.30), fixture(1.0), settings()).unwrap();
        assert_ne!(assessed.snapshot, thinner.snapshot);
        assert!(thinner.volume < assessed.volume);
        assert!(thinner.estimate.j_primal > e.j_primal);
    }

    #[test]
    fn dwr_phase_stops_do_not_publish_or_mutate_the_design() {
        let phi = beam(0.35);
        let before = phi.nodes().to_vec();
        for stop in [ComplianceDwrStage::BeforeEstimate, ComplianceDwrStage::BeforePublish] {
            let result = assess_compliance_dwr_controlled(
                &phi, fixture(1.0), settings(), |stage| {
                    if stage == stop { ControlFlow::Break(stage) }
                    else { ControlFlow::Continue(()) }
                },
            ).unwrap();
            assert!(matches!(result, ControlFlow::Break(stage) if stage == stop));
            assert_eq!(phi.nodes(), before);
        }
    }
}
