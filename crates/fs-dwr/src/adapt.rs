//! The octree h-refinement loop (mechanism 1 of 4): solve → estimate →
//! Dörfler-mark → constraint-aware split admission → rebalance. Scalar
//! ghost stabilization advances the whole interface only when a marked
//! candidate raises the declared one-cell-halo target. A mismatch that leaves
//! that target unchanged is deferred in favor of the next deterministic
//! indicator. The accuracy-per-DOF trajectory is the ledgered evidence —
//! goal-oriented refinement must beat uniform on localized QoIs or the
//! estimator is decoration.

use crate::estimate::{DwrEstimate, GoalContext, estimate};
use crate::mark::{dorfler, indicator_order};
use fs_cutfem::{CellKey, CutFemError, CutSdf, FemParams, Quadtree, Space};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// One adaptive iteration's evidence.
#[derive(Debug, Clone)]
pub struct AdaptStep {
    /// Primal free DOFs at this step.
    pub dofs: usize,
    /// J(u_h).
    pub j: f64,
    /// Signed estimate.
    pub eta_signed: f64,
    /// Marking mass Σ|η_K|.
    pub eta_abs: f64,
    /// Indicator cells actually refined by the constrained marking plan
    /// (0 on the final, estimate-only step).
    pub marked: usize,
}

impl AdaptStep {
    /// Ledger-style JSON row.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"dofs\":{},\"j\":{:.10e},\"eta_signed\":{:.4e},\
             \"eta_abs\":{:.4e},\"marked\":{}}}",
            self.dofs, self.j, self.eta_signed, self.eta_abs, self.marked
        );
        s
    }
}

fn refinement_ceiling(max_level: u32) -> Result<u32, CutFemError> {
    max_level
        .checked_sub(1)
        .ok_or_else(|| CutFemError::InvalidFemInput {
            what: "scalar DWR adaptivity requires refinement headroom when another marked iteration remains"
                .to_string(),
        })
}

fn validate_theta(theta: f64) -> Result<(), CutFemError> {
    if theta.is_finite() && (0.0..=1.0).contains(&theta) {
        Ok(())
    } else {
        Err(CutFemError::InvalidFemInput {
            what: format!("DWR marking fraction theta must be finite and in [0, 1], got {theta}"),
        })
    }
}

fn interface_target_level(grid: &Quadtree, sdf: &dyn CutSdf) -> u32 {
    let mut target = 0u32;
    for cell in grid.leaves() {
        let (lo, hi) = grid.rect(cell);
        let h = hi[0] - lo[0];
        let inflated_lo = [(lo[0] - h).max(0.0), (lo[1] - h).max(0.0)];
        let inflated_hi = [(hi[0] + h).min(1.0), (hi[1] + h).min(1.0)];
        if sdf.enclose(inflated_lo, inflated_hi).contains_zero() {
            target = target.max(cell.0);
        }
    }
    target
}

fn constrained_refinement(
    grid: &Quadtree,
    sdf: &dyn CutSdf,
    params: FemParams,
    indicators: &BTreeMap<CellKey, f64>,
    theta: f64,
    ceiling: u32,
) -> Result<(Quadtree, usize), CutFemError> {
    for &cell in indicators.keys() {
        if !grid.is_leaf(cell) {
            return Err(CutFemError::InvalidFemInput {
                what: format!(
                    "DWR indicator cell {cell:?} is not a leaf of the admitted analysis grid"
                ),
            });
        }
    }
    let total: f64 = indicators.values().map(|value| value.abs()).sum();
    let target = theta * total;
    if total <= 0.0 || target <= 0.0 {
        return Ok((grid.clone(), 0));
    }

    let mut planned = grid.clone();
    let mut refined = BTreeSet::new();
    let mut refined_mass = 0.0;
    let mut protected = 0usize;
    let mut capped = 0usize;
    let mut admitted_band = interface_target_level(&planned, sdf);
    for (cell, _) in indicator_order(indicators) {
        if refined_mass >= target {
            break;
        }
        if !planned.is_leaf(cell) {
            continue;
        }
        if cell.0 >= ceiling {
            capped += 1;
            continue;
        }

        let mut trial = planned.clone();
        trial.split(cell);
        trial.balance();
        let requested_band = interface_target_level(&trial, sdf);
        let admitted = match Space::build(&trial, sdf, params).map(|_| ()) {
            Ok(()) => true,
            Err(CutFemError::CutBandNotUniform { .. }) => {
                if requested_band > admitted_band {
                    // A target-raising mismatch deliberately advances the
                    // interface (the original, quality-gated 3 -> 4 policy).
                    // An unchanged target is merely an incidental mismatch.
                    trial.refine_toward_interface(sdf, requested_band);
                    Space::build(&trial, sdf, params).map(|_| ())?;
                    true
                } else {
                    protected += 1;
                    false
                }
            }
            Err(error) => return Err(error),
        };
        if !admitted {
            continue;
        }

        planned = trial;
        admitted_band = interface_target_level(&planned, sdf);
        for (&indicator_cell, &indicator) in indicators {
            if !planned.is_leaf(indicator_cell) && refined.insert(indicator_cell) {
                refined_mass += indicator.abs();
            }
        }
    }

    if refined_mass < target {
        return Err(CutFemError::InvalidFemInput {
            what: format!(
                "constrained DWR marking found only {refined_mass:.6e} admissible mass toward \
                 {target:.6e}; {protected} candidate(s) would violate the scalar ghost cut band \
                 and {capped} candidate(s) lacked refinement headroom"
            ),
        });
    }
    Ok((planned, refined.len()))
}

/// Run `iters` adaptive cycles (the last records without refining).
/// The grid must carry enough `with_room` headroom for the splits.
///
/// # Errors
/// Propagates fs-cutfem build/solve errors and the DWR estimator's structured
/// refusals for non-nested active coverage or missing/non-finite field
/// evidence. A non-final iteration with marked cells also refuses when the
/// grid has no reserved refinement headroom or the requested Dörfler mass
/// cannot be refined without violating scalar ghost-band admission. A
/// non-finite marking fraction or one outside `[0, 1]` refuses before solving.
#[allow(clippy::too_many_arguments)] // the PDE problem statement is the argument list
pub fn adapt_loop(
    grid: &mut Quadtree,
    sdf: &dyn CutSdf,
    params: FemParams,
    f: &dyn Fn(f64, f64) -> f64,
    g: &dyn Fn(f64, f64) -> f64,
    goal: &GoalContext<'_>,
    theta: f64,
    iters: usize,
) -> Result<(Vec<AdaptStep>, DwrEstimate), CutFemError> {
    validate_theta(theta)?;
    let mut steps = Vec::new();
    loop {
        let est = estimate(grid, sdf, params, f, g, goal)?;
        let last = steps.len() + 1 >= iters;
        if last {
            steps.push(AdaptStep {
                dofs: est.dofs,
                j: est.j_primal,
                eta_signed: est.eta_signed,
                eta_abs: est.eta_abs,
                marked: 0,
            });
            return Ok((steps, est));
        }
        let requested = dorfler(&est.indicators, theta);
        let admitted = if requested.is_empty() {
            0
        } else {
            let ceiling = refinement_ceiling(grid.max_level())?;
            let (planned, admitted) =
                constrained_refinement(grid, sdf, params, &est.indicators, theta, ceiling)?;
            *grid = planned;
            admitted
        };
        steps.push(AdaptStep {
            dofs: est.dofs,
            j: est.j_primal,
            eta_signed: est.eta_signed,
            eta_abs: est.eta_abs,
            marked: admitted,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::{
        adapt_loop, constrained_refinement, interface_target_level, refinement_ceiling,
        validate_theta,
    };
    use crate::estimate::{GoalContext, estimate};
    use crate::mark::dorfler;
    use fs_cutfem::{Circle, CutFemError, CutSdf, FemParams, Quadtree, Space};
    use std::collections::BTreeMap;

    #[test]
    fn constrained_refinement_defers_inside_halo_without_global_band_promotion() {
        let sdf = Circle {
            center: [0.5, 0.5],
            radius: 0.42,
        };
        let grid = Quadtree::with_room(4, 6);
        let protected_halo = (4, 13, 7);
        let safe_interior = (4, 7, 7);
        for cell in [protected_halo, safe_interior] {
            let (lo, hi) = grid.rect(cell);
            assert!(
                sdf.enclose(lo, hi).hi() < 0.0,
                "fixture cell {cell:?} must be certified inside rather than cut"
            );
        }
        let indicators = BTreeMap::from([(protected_halo, 0.6), (safe_interior, 0.4)]);
        assert_eq!(dorfler(&indicators, 0.4), vec![protected_halo]);
        let ceiling = refinement_ceiling(grid.max_level()).expect("fixture has headroom");

        let mut raw_trial = grid.clone();
        raw_trial.split(protected_halo);
        raw_trial.balance();
        assert!(matches!(
            Space::build(&raw_trial, &sdf, FemParams::default()),
            Err(CutFemError::CutBandNotUniform { cell, neighbor })
                if cell == (4, 14, 7) && neighbor == (5, 27, 15)
        ));

        let (planned, marked) =
            constrained_refinement(&grid, &sdf, FemParams::default(), &indicators, 0.4, ceiling)
                .expect("a lower-ranked interior mark preserves scalar admission");
        assert_eq!(marked, 1);
        assert!(grid.is_leaf(protected_halo) && grid.is_leaf(safe_interior));
        assert!(planned.is_leaf(protected_halo));
        assert!(!planned.is_leaf(safe_interior));
        assert_eq!(planned.leaves().count(), grid.leaves().count() + 3);
        Space::build(&planned, &sdf, FemParams::default())
            .expect("constrained plan must remain independently admissible");

        let (replayed, replay_marked) =
            constrained_refinement(&grid, &sdf, FemParams::default(), &indicators, 0.4, ceiling)
                .expect("deterministic replay");
        assert_eq!(replay_marked, marked);
        assert_eq!(
            replayed.leaves().collect::<Vec<_>>(),
            planned.leaves().collect::<Vec<_>>()
        );

        let protected_only = BTreeMap::from([(protected_halo, 1.0)]);
        let error = constrained_refinement(
            &grid,
            &sdf,
            FemParams::default(),
            &protected_only,
            0.5,
            ceiling,
        )
        .expect_err("an all-protected marked mass must refuse without a partial plan");
        assert!(matches!(
            error,
            CutFemError::InvalidFemInput { what }
                if what.contains("admissible mass") && what.contains("scalar ghost cut band")
        ));

        let ghost_free = FemParams {
            ghost_gamma: 0.0,
            ..FemParams::default()
        };
        let (graded, ghost_free_marked) =
            constrained_refinement(&grid, &sdf, ghost_free, &protected_only, 0.5, ceiling)
                .expect("ghost-free scalar admission accepts the same local grade");
        assert_eq!(ghost_free_marked, 1);
        assert!(!graded.is_leaf(protected_halo));
        Space::build(&graded, &sdf, ghost_free).expect("ghost-free graded plan builds");
    }

    #[test]
    fn constrained_refinement_retains_deliberate_band_target_advance() {
        let sdf = Circle {
            center: [0.5, 0.5],
            radius: 0.42,
        };
        let grid = Quadtree::with_room(3, 6);
        let target_raising_halo = (3, 6, 3);
        let indicators = BTreeMap::from([(target_raising_halo, 1.0)]);
        let ceiling = refinement_ceiling(grid.max_level()).expect("fixture has headroom");
        assert_eq!(interface_target_level(&grid, &sdf), 3);

        let (planned, marked) =
            constrained_refinement(&grid, &sdf, FemParams::default(), &indicators, 0.5, ceiling)
                .expect("a marked halo cell may deliberately advance the band target");
        assert_eq!(marked, 1);
        assert_eq!(interface_target_level(&planned, &sdf), 4);
        assert!(!planned.is_leaf(target_raising_halo));
        assert!(
            planned.leaves().count() > grid.leaves().count() + 3,
            "authorized band advancement must be distinct from one local split"
        );
        Space::build(&planned, &sdf, FemParams::default())
            .expect("deliberately advanced band must remain admissible");

        let ghost_free = FemParams {
            ghost_gamma: 0.0,
            ..FemParams::default()
        };
        let (local, local_marked) =
            constrained_refinement(&grid, &sdf, ghost_free, &indicators, 0.5, ceiling)
                .expect("raw-admissible target-raising split must remain local");
        assert_eq!(local_marked, 1);
        assert_eq!(local.leaves().count(), grid.leaves().count() + 3);
        Space::build(&local, &sdf, ghost_free).expect("ghost-free local target raiser builds");
    }

    #[test]
    fn adapt_loop_refuses_zero_level_marked_headroom_without_bypassing_safe_paths() {
        for theta in [f64::NAN, -f64::EPSILON, 1.0 + f64::EPSILON] {
            let error = validate_theta(theta).expect_err("invalid theta must refuse");
            assert!(matches!(
                error,
                CutFemError::InvalidFemInput { what } if what.contains("theta")
            ));
        }
        validate_theta(0.0).expect("zero theta is the documented empty marking");
        validate_theta(1.0).expect("full-mass theta is valid");

        let error = refinement_ceiling(0).expect_err("zero-level grid has no refinement headroom");
        assert!(matches!(
            error,
            CutFemError::InvalidFemInput { what }
                if what.contains("refinement headroom")
        ));
        assert_eq!(
            refinement_ceiling(1).expect("one level has a zero ceiling"),
            0
        );
        assert_eq!(refinement_ceiling(4).expect("positive headroom"), 3);

        let sdf = Circle {
            center: [0.5, 0.5],
            radius: 0.35,
        };
        let params = FemParams {
            ghost_gamma: 0.0,
            ..FemParams::default()
        };
        let source = |x: f64, y: f64| 1.0 + x + 2.0 * y;
        let boundary = |_: f64, _: f64| 0.0;
        let weight = |x: f64, y: f64| 1.0 + x * x + y;
        let goal = GoalContext { weight: &weight };
        let probe_grid = Quadtree::uniform(0);
        let probe = estimate(&probe_grid, &sdf, params, &source, &boundary, &goal)
            .expect("level-zero fixture has a valid DWR estimate");
        assert!(
            !dorfler(&probe.indicators, 0.5).is_empty(),
            "nonzero fixture must reach the marked headroom gate"
        );

        let mut final_only = Quadtree::uniform(0);
        let (steps, _) = adapt_loop(
            &mut final_only,
            &sdf,
            params,
            &source,
            &boundary,
            &goal,
            0.5,
            1,
        )
        .expect("the final estimate-only iteration requires no headroom");
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].marked, 0);

        let mut no_headroom = Quadtree::uniform(0);
        let error = adapt_loop(
            &mut no_headroom,
            &sdf,
            params,
            &source,
            &boundary,
            &goal,
            0.5,
            2,
        )
        .expect_err("a marked continuation must refuse zero-level headroom");
        assert!(matches!(
            error,
            CutFemError::InvalidFemInput { what }
                if what.contains("refinement headroom")
        ));

        let zero = |_: f64, _: f64| 0.0;
        let empty_goal = GoalContext { weight: &zero };
        let mut empty_marking = Quadtree::uniform(0);
        let (steps, _) = adapt_loop(
            &mut empty_marking,
            &sdf,
            params,
            &zero,
            &zero,
            &empty_goal,
            0.5,
            2,
        )
        .expect("an empty marked set requires no refinement headroom");
        assert_eq!(steps.len(), 2);
        assert!(steps.iter().all(|step| step.marked == 0));
    }
}
