//! Uniform spatial refinement of complete native thermal trajectories.
//!
//! Every rung starts from the declared physical initial state and uses the
//! same two event-aligned time grids. The forward-step allowance belongs to
//! the whole spatial study. The final temperature's spatial and temporal
//! estimates remain separate until both support the stated additive estimate.

use std::cell::{Cell, RefCell};

use fs_project::ProjectSpec;

use super::super::{
    LadderDiscretization, LadderEstimate, SolveRefusal, canonical_f64, json_string,
};
use super::{TimeGrids, bad, policy};

pub(in crate::solve) struct SpatialWork {
    steps_per_rung: usize,
    limit: usize,
    reserved: Cell<usize>,
    temporal: RefCell<Vec<Option<f64>>>,
}

impl SpatialWork {
    pub(in crate::solve) fn new(spec: &ProjectSpec) -> Result<Option<Self>, SolveRefusal> {
        let Some(policy) = policy(spec) else {
            return Ok(None);
        };
        if !spec
            .solver
            .as_ref()
            .is_some_and(|solver| solver.fidelity == "ladder")
        {
            return Ok(None);
        }
        let grids = TimeGrids::new(policy)?;
        let steps_per_rung = grids
            .coarse
            .len()
            .checked_add(grids.fine.len())
            .ok_or_else(|| bad("spatial/temporal work count overflows"))?;
        let limit = (policy.max_steps as usize).min(10_000);
        if steps_per_rung == 0 || steps_per_rung > limit {
            return Err(bad(
                "one complete coarse/fine trajectory exceeds the shared spatial-ladder time-step allowance",
            ));
        }
        Ok(Some(Self {
            steps_per_rung,
            limit,
            reserved: Cell::new(0),
            temporal: RefCell::new(Vec::new()),
        }))
    }

    pub(in crate::solve) fn can_start_rung(&self) -> bool {
        self.reserved
            .get()
            .checked_add(self.steps_per_rung)
            .is_some_and(|next| next <= self.limit)
    }

    /// Reserve the complete coarse/fine pair before any trajectory work.
    /// Failed work cannot be followed by a fresh rung within this invocation.
    pub(in crate::solve) fn begin_rung(&self) -> Result<(), SolveRefusal> {
        if self
            .temporal
            .borrow()
            .len()
            .checked_mul(self.steps_per_rung)
            != Some(self.reserved.get())
        {
            return Err(bad(
                "the previous spatial rung did not complete its trajectory",
            ));
        }
        let next = self
            .reserved
            .get()
            .checked_add(self.steps_per_rung)
            .filter(|next| *next <= self.limit)
            .ok_or_else(|| {
                bad("the spatial ladder exhausted its shared forward time-step allowance")
            })?;
        self.reserved.set(next);
        Ok(())
    }

    pub(in crate::solve) fn finish_rung(
        &self,
        temporal_half_width_k: Option<f64>,
    ) -> Result<(), SolveRefusal> {
        let mut rows = self.temporal.borrow_mut();
        if rows
            .len()
            .checked_add(1)
            .and_then(|count| count.checked_mul(self.steps_per_rung))
            != Some(self.reserved.get())
            || temporal_half_width_k.is_some_and(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(bad(
                "the completed spatial rung has inconsistent time work or temporal evidence",
            ));
        }
        rows.push(temporal_half_width_k);
        Ok(())
    }

    pub(in crate::solve) fn report(
        &self,
        spatial: &LadderEstimate,
    ) -> Result<(String, Option<LadderDiscretization>), SolveRefusal> {
        let temporal = self.temporal.borrow();
        if temporal.len() != spatial.rungs
            || temporal.len().checked_mul(self.steps_per_rung) != Some(self.reserved.get())
        {
            return Err(bad(
                "spatial and temporal convergence rows do not describe the same completed rungs",
            ));
        }
        let finite = |value| {
            canonical_f64(value).ok_or_else(|| bad("space-time convergence evidence is nonfinite"))
        };
        let optional = |value: Option<f64>| value.map_or(Ok("null".to_string()), finite);
        let final_temporal = temporal.last().copied().flatten();
        // A data-range fallback measures variation rather than asymptotic
        // error. Bitwise agreement does not establish zero spatial error.
        // Neither may become a full space-time error estimate by addition.
        let spatial_error = if spatial.status == "observed-order" {
            spatial
                .half_width_k
                .filter(|value| value.is_finite() && *value >= 0.0)
        } else {
            None
        };
        let combined = spatial_error
            .zip(final_temporal)
            .map(|(space, time)| space + time);
        if combined.is_some_and(|value| !value.is_finite()) {
            return Err(bad("combined spatial and temporal estimate overflows"));
        }
        let status = if combined.is_some() {
            "estimated-space-time"
        } else if spatial.rungs < 3 {
            "incomplete-spatial-ladder"
        } else if spatial_error.is_none() {
            "unresolved-spatial-order"
        } else {
            "unresolved-temporal-error"
        };
        let rows = temporal
            .iter()
            .enumerate()
            .map(|(rung, &estimate)| {
                Ok(format!(
                    "{{\"rung\":{rung},\"forward_steps\":{},\"temporal_half_width_k\":{}}}",
                    self.steps_per_rung,
                    optional(estimate)?,
                ))
            })
            .collect::<Result<Vec<_>, SolveRefusal>>()?
            .join(",");
        let json = format!(
            "{{\"status\":{},\"forward_steps_per_rung\":{},\"forward_steps_completed\":{},\"max_forward_steps\":{},\"rungs\":[{}],\"spatial_status\":{},\"spatial_half_width_k\":{},\"spatial_order\":{},\"final_temporal_half_width_k\":{},\"temporal_assumed_order\":1,\"temporal_safety_factor\":1.25,\"estimated_half_width_k\":{},\"authority\":\"Estimated\",\"scope\":{}}}",
            json_string(status),
            self.steps_per_rung,
            self.reserved.get(),
            self.limit,
            rows,
            json_string(spatial.status),
            optional(spatial_error)?,
            optional(spatial.order)?,
            optional(final_temporal)?,
            optional(combined)?,
            json_string(
                "Each spatial rung independently advances the same declared initial state and event-aligned coarse/fine time grids. The shared forward-step allowance counts every completed trajectory on every rung. Spatial Richardson/GCI compares final maxima at the common fine time grid; temporal error is the final mesh's coarse/fine comparison at assumed backward-Euler order one. Their sum is an additive asymptotic estimate, not a continuum, continuous-time maximum, phase-fraction, or enthalpy bound. Missing temporal evidence, insufficient spatial rungs, bitwise spatial agreement, or an unobserved spatial order leaves the full Discretization term unmeasured. Adjoint reverse work retains its existing separately bounded trajectory policy."
            ),
        );
        let term = combined.map(|half_width_k| LadderDiscretization {
            half_width_k,
            status: "transient-space-time",
            order: spatial.order,
            rungs: spatial.rungs,
        });
        Ok((json, term))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn work(steps_per_rung: usize, limit: usize) -> SpatialWork {
        SpatialWork {
            steps_per_rung,
            limit,
            reserved: Cell::new(0),
            temporal: RefCell::new(Vec::new()),
        }
    }

    fn spatial(status: &'static str, rungs: usize) -> LadderEstimate {
        LadderEstimate {
            status,
            stop: "complete",
            rungs,
            order: Some(2.0),
            extrapolated_k: Some(300.0),
            half_width_k: Some(0.25),
            detail: "test-only spatial observation".into(),
        }
    }

    #[test]
    fn g0_every_complete_trajectory_consumes_the_original_shared_allowance() {
        let work = work(12, 25);
        for _ in 0..2 {
            assert!(work.can_start_rung());
            work.begin_rung().unwrap();
            work.finish_rung(Some(0.1)).unwrap();
        }
        assert!(!work.can_start_rung());
        assert!(work.begin_rung().is_err());
        let (json, term) = work.report(&spatial("two-rungs", 2)).unwrap();
        assert!(json.contains("\"forward_steps_completed\":24"));
        assert!(json.contains("\"max_forward_steps\":25"));
        assert!(term.is_none());
    }

    #[test]
    fn g0_incomplete_or_invalid_rungs_cannot_publish_or_reset_work() {
        let work = work(12, 36);
        work.begin_rung().unwrap();
        assert!(work.begin_rung().is_err());
        assert!(work.report(&spatial("single-rung", 1)).is_err());
        for bad in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
            assert!(work.finish_rung(Some(bad)).is_err());
        }
        work.finish_rung(Some(0.1)).unwrap();
        assert!(work.finish_rung(Some(0.1)).is_err());
    }

    #[test]
    fn g0_full_discretization_combines_final_mesh_time_error_with_observed_space_error() {
        let work = work(12, 36);
        for estimate in [Some(0.8), Some(0.4), Some(0.1)] {
            work.begin_rung().unwrap();
            work.finish_rung(estimate).unwrap();
        }
        let (json, term) = work.report(&spatial("observed-order", 3)).unwrap();
        let term = term.unwrap();
        assert_eq!(term.half_width_k, 0.35);
        assert_eq!(term.status, "transient-space-time");
        assert!(json.contains("\"status\":\"estimated-space-time\""));
        assert!(json.contains("\"forward_steps_completed\":36"));
    }

    #[test]
    fn g0_missing_or_non_asymptotic_evidence_never_becomes_complete_discretization() {
        for status in ["observed-order", "data-range", "converged-exactly"] {
            let work = work(12, 36);
            for _ in 0..3 {
                work.begin_rung().unwrap();
                work.finish_rung(if status == "observed-order" {
                    None
                } else {
                    Some(0.1)
                })
                .unwrap();
            }
            assert!(work.report(&spatial(status, 3)).unwrap().1.is_none());
        }
    }
}
