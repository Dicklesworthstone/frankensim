//! Bounded normal-load control on the original accepted plastic history.
//!
//! Every pressure trial uses the same predecessor plastic depth and work. Only
//! its prescribed approach changes; trial plastic states are discarded. A
//! recession-at-load solve supplies the new recession directly to those trials,
//! so no fictitious held-approach plastic intermediate is committed.
//!
//! Total force is monotone for the convex pressure box problem. A safeguarded
//! secant/bisection search retains tested approach bounds and reports actual
//! force residuals. A positive perfectly plastic force plateau does not select
//! a unique approach and is refused. Zero load uses the explicitly supplied
//! lower approach bound, unless the accepted open state can be retained.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;

use super::{
    FiniteGapPlasticAdvance, FiniteGapPlasticContinuation, FiniteGapPlasticError,
    FiniteGapPlasticModel, FiniteGapPlasticState, checkpoint,
};

/// Explicit bounded search controls for one prescribed compressive resultant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticLoadControls {
    /// Lower admitted approach, also the unloading position for zero load [m].
    pub minimum_approach_m: f64,
    /// Upper admitted approach [m]. Both tested endpoints must be admissible.
    pub maximum_approach_m: f64,
    /// Absolute accepted normal-force residual [N].
    pub absolute_force_tolerance_n: f64,
    /// Additional accepted residual divided by requested force, in [0, 1).
    pub relative_force_tolerance: f64,
    /// Smallest resolvable approach bracket or saturation interval [m].
    ///
    /// Reaching this width never substitutes for satisfying force tolerance.
    /// This is not a rigorous position enclosure beyond the native contact
    /// discretization and complementarity tolerance.
    pub approach_tolerance_m: f64,
    /// Maximum native equilibrium trials, including the two bracket endpoints.
    pub maximum_equilibrium_evaluations: usize,
}

/// How the final force-controlled position was selected.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum FiniteGapPlasticLoadDisposition {
    /// A converged native response satisfies the requested positive resultant.
    ForceBalanced = 0,
    /// An unchanged accepted equilibrium already satisfies the load tolerance.
    /// No assertion of unique position on a plastic plateau is made.
    RetainedExistingEquilibrium = 1,
    /// Zero force was realized at the caller's explicit lower approach bound.
    UnloadedAtLowerApproach = 2,
}

/// Candidate physical contact response and the bounded load-solver receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLoadAdvance {
    /// Original normal/plastic response from the accepted predecessor.
    pub contact: FiniteGapPlasticAdvance,
    /// Load request, actual residual, tested bounds, and numerical work.
    pub load: FiniteGapPlasticLoadReceipt,
}

/// Geometry continuation at a maintained load, with separate actuator work.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLoadRecessionAdvance {
    /// Accepted-to-candidate recession with normal and parameter work separated.
    pub continuation: FiniteGapPlasticContinuation,
    /// Load-control evidence for the same physical candidate.
    pub load: FiniteGapPlasticLoadReceipt,
}

/// Numerical evidence for one bounded force-control operation.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLoadReceipt {
    /// Content identity of the request, bounds, and complete physical response.
    pub identity: ContentHash,
    /// Original immutable physical model.
    pub model_identity: ContentHash,
    /// Accepted history shared by every trial.
    pub predecessor_state_identity: ContentHash,
    /// Selected candidate history.
    pub candidate_state_identity: ContentHash,
    /// Original native normal/plastic receipt identity.
    pub response_identity: ContentHash,
    /// Requested nonnegative compressive force [N].
    pub target_force_n: f64,
    /// Integrated force of the returned pressure field [N].
    pub achieved_force_n: f64,
    /// Signed achieved minus requested force [N].
    pub force_residual_n: f64,
    /// Final tested search bounds, or the explicitly retained/unloaded point [m].
    /// These are numerical search metadata, not an independent rigorous enclosure.
    pub approach_bracket_m: [f64; 2],
    /// Native resultants at the reported tested bounds [N].
    pub force_bracket_n: [f64; 2],
    /// Native equilibrium trials; exact history reuse requires zero trials.
    pub equilibrium_evaluations: usize,
    /// Complete numerical controls consumed by this operation.
    pub controls: FiniteGapPlasticLoadControls,
    /// Position-selection convention actually used.
    pub disposition: FiniteGapPlasticLoadDisposition,
}

/// Load-control refusal without publishing or adopting any trial history.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGapPlasticLoadError {
    /// Invalid force, search controls, or recession input.
    InvalidInput { field: &'static str },
    /// Original contact, ownership, cancellation, or finite-domain refusal.
    Contact(FiniteGapPlasticError),
    /// The requested force exceeds even the sum of all admitted interior caps.
    /// This upper bound is necessary, not a claim that it can be attained.
    TargetExceedsCapacity {
        target_force_n: f64,
        capacity_upper_bound_n: f64,
    },
    /// The tested approach interval does not realize the requested force.
    TargetOutsideApproachBounds {
        target_force_n: f64,
        approach_bracket_m: [f64; 2],
        force_bracket_n: [f64; 2],
    },
    /// The exact requested resultant lies on a genuine all-bound pressure plateau.
    NonUniqueLoad {
        target_force_n: f64,
        plateau_force_n: f64,
        /// Dual-feasible approach interval. The upper endpoint can be unbounded.
        plateau_approach_m: [f64; 2],
    },
    /// Work/resolution limits or an unresolved numerical response prevent balance.
    DidNotConverge {
        equilibrium_evaluations: usize,
        approach_bracket_m: [f64; 2],
        best_force_residual_n: Option<f64>,
    },
    /// Recession-at-load requires the accepted state already to carry that load.
    PredecessorLoadMismatch {
        target_force_n: f64,
        predecessor_force_n: f64,
        force_tolerance_n: f64,
    },
}

impl From<FiniteGapPlasticError> for FiniteGapPlasticLoadError {
    fn from(value: FiniteGapPlasticError) -> Self {
        Self::Contact(value)
    }
}

impl fmt::Display for FiniteGapPlasticLoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Contact(error) => write!(f, "{error}"),
            Self::InvalidInput { field } => write!(f, "invalid finite-gap load-control {field}"),
            Self::TargetExceedsCapacity {
                target_force_n,
                capacity_upper_bound_n,
            } => write!(
                f,
                "requested force {target_force_n} N exceeds interior pressure-cap bound {capacity_upper_bound_n} N"
            ),
            Self::TargetOutsideApproachBounds {
                target_force_n,
                approach_bracket_m,
                force_bracket_n,
            } => write!(
                f,
                "requested force {target_force_n} N is outside tested forces {force_bracket_n:?} N at approaches {approach_bracket_m:?} m"
            ),
            Self::NonUniqueLoad {
                target_force_n,
                plateau_approach_m,
                ..
            } => write!(
                f,
                "force {target_force_n} N does not select a unique approach on plastic plateau {plateau_approach_m:?} m"
            ),
            Self::DidNotConverge {
                equilibrium_evaluations,
                approach_bracket_m,
                best_force_residual_n,
            } => write!(
                f,
                "load control unresolved after {equilibrium_evaluations} equilibria in {approach_bracket_m:?} m; best force residual {best_force_residual_n:?} N"
            ),
            Self::PredecessorLoadMismatch {
                target_force_n,
                predecessor_force_n,
                force_tolerance_n,
            } => write!(
                f,
                "maintained force {target_force_n} N differs from accepted force {predecessor_force_n} N beyond tolerance {force_tolerance_n} N"
            ),
        }
    }
}

impl core::error::Error for FiniteGapPlasticLoadError {}

struct LoadTrial {
    approach_m: f64,
    force_n: f64,
    contact: FiniteGapPlasticAdvance,
}

struct LoadProgress {
    target_force_n: f64,
    force_tolerance_n: f64,
    controls: FiniteGapPlasticLoadControls,
    evaluations: usize,
    bracket_m: [f64; 2],
    best_residual_n: Option<f64>,
    zero_pressure: Vec<f64>,
}

impl LoadProgress {
    fn unresolved(&self) -> FiniteGapPlasticLoadError {
        FiniteGapPlasticLoadError::DidNotConverge {
            equilibrium_evaluations: self.evaluations,
            approach_bracket_m: self.bracket_m,
            best_force_residual_n: self.best_residual_n,
        }
    }
}

impl FiniteGapPlasticModel {
    /// Find an admitted approach for a prescribed nonnegative normal resultant.
    ///
    /// All trials start from `before`'s plastic/recession/work history. No trial
    /// is accepted until the force criterion and uniqueness checks succeed.
    /// Changing force does not specify its path, so this method does not invent
    /// a normal-work integral for the finite load change.
    pub fn advance_at_load(
        &self,
        before: &FiniteGapPlasticState,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLoadAdvance, FiniteGapPlasticLoadError> {
        self.load_at_recession(before, &before.recession_m, target_force_n, controls, cx)
    }

    /// Apply recession while maintaining the predecessor's compressive load.
    ///
    /// The accepted resultant must already match `target_force_n` within the
    /// specified tolerance. All trials combine the NEW recession with ORIGINAL
    /// accepted plastic depth and work. The returned normal work is the constant
    /// controlled force times approach change; parameter work subtracts it from
    /// elastic-storage change plus plastic work. Geometry-path refinement and
    /// force tolerances still limit this endpoint bookkeeping.
    pub fn recede_at_load(
        &self,
        before: &FiniteGapPlasticState,
        recession_increment_m: &[f64],
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLoadRecessionAdvance, FiniteGapPlasticLoadError> {
        checkpoint(cx)?;
        self.check_state_owner(before)?;
        let tolerance = validate_controls(target_force_n, controls)?;
        let predecessor_force_n = self.accepted_force(before)?;
        if (predecessor_force_n - target_force_n).abs() > tolerance {
            return Err(FiniteGapPlasticLoadError::PredecessorLoadMismatch {
                target_force_n,
                predecessor_force_n,
                force_tolerance_n: tolerance,
            });
        }
        if recession_increment_m.len() != self.gap_m.len()
            || recession_increment_m
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid_load("finite nonnegative recession field"));
        }
        let mut recession = before.recession_m.clone();
        let mut actual_increment = vec![0.0; recession.len()];
        for (i, requested) in recession_increment_m.iter().enumerate() {
            recession[i] += requested;
            if !recession[i].is_finite() {
                return Err(invalid_load("finite accumulated recession"));
            }
            actual_increment[i] = recession[i] - before.recession_m[i];
            if *requested > 0.0 && actual_increment[i] == 0.0 {
                return Err(invalid_load("representable recession increment"));
            }
        }
        let result = self.load_at_recession(before, &recession, target_force_n, controls, cx)?;
        let continuation = self.finish_load_recession(
            before,
            result.contact,
            actual_increment,
            target_force_n,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(FiniteGapPlasticLoadRecessionAdvance {
            continuation,
            load: result.load,
        })
    }

    fn load_at_recession(
        &self,
        before: &FiniteGapPlasticState,
        recession_m: &[f64],
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLoadAdvance, FiniteGapPlasticLoadError> {
        checkpoint(cx)?;
        self.check_state_owner(before)?;
        let force_tolerance_n = validate_controls(target_force_n, controls)?;
        let capacity_upper_bound_n = self.interior_force_capacity_bound();
        if target_force_n > capacity_upper_bound_n {
            return Err(FiniteGapPlasticLoadError::TargetExceedsCapacity {
                target_force_n,
                capacity_upper_bound_n,
            });
        }
        let mut progress = LoadProgress {
            target_force_n,
            force_tolerance_n,
            controls,
            evaluations: 0,
            bracket_m: [controls.minimum_approach_m, controls.maximum_approach_m],
            best_residual_n: None,
            zero_pressure: vec![0.0; self.gap_m.len()],
        };
        let previous_force = self.accepted_force(before)?;
        if recession_m == before.recession_m.as_slice()
            && before.approach_m >= controls.minimum_approach_m
            && before.approach_m <= controls.maximum_approach_m
            && (previous_force - target_force_n).abs() <= force_tolerance_n
            && ((target_force_n > 0.0 && previous_force > 0.0)
                || (target_force_n == 0.0 && before.pressure_pa.iter().all(|p| *p == 0.0)))
        {
            let (_, opening) = self.stored_energy_and_opening(before, cx)?;
            let contact = self.unchanged_advance(before, opening, cx)?;
            progress.bracket_m = [before.approach_m; 2];
            return self.finish_load(
                before,
                LoadTrial {
                    approach_m: before.approach_m,
                    force_n: previous_force,
                    contact,
                },
                &progress,
                [previous_force; 2],
                FiniteGapPlasticLoadDisposition::RetainedExistingEquilibrium,
                cx,
            );
        }
        let mut lower = self.load_trial(
            before,
            recession_m,
            controls.minimum_approach_m,
            &mut progress,
            cx,
        )?;
        if target_force_n == 0.0 {
            if lower.contact.candidate.pressure_pa.iter().any(|p| *p > 0.0) {
                return Err(FiniteGapPlasticLoadError::TargetOutsideApproachBounds {
                    target_force_n,
                    approach_bracket_m: [lower.approach_m; 2],
                    force_bracket_n: [lower.force_n; 2],
                });
            }
            progress.bracket_m = [lower.approach_m; 2];
            let forces = [lower.force_n; 2];
            return self.finish_load(
                before,
                lower,
                &progress,
                forces,
                FiniteGapPlasticLoadDisposition::UnloadedAtLowerApproach,
                cx,
            );
        }
        if controls.minimum_approach_m == controls.maximum_approach_m {
            if self.load_trial_balanced(before, recession_m, &lower, &progress)? {
                let forces = [lower.force_n; 2];
                return self.finish_load(
                    before,
                    lower,
                    &progress,
                    forces,
                    FiniteGapPlasticLoadDisposition::ForceBalanced,
                    cx,
                );
            }
            return Err(progress.unresolved());
        }
        let mut upper = self.load_trial(
            before,
            recession_m,
            controls.maximum_approach_m,
            &mut progress,
            cx,
        )?;
        if upper.force_n + force_tolerance_n < lower.force_n {
            return Err(progress.unresolved());
        }
        // A repeated force with free pressure variables is not proof of a
        // physical plastic plateau. Do not manufacture an approach inside it.
        if lower.force_n == upper.force_n
            && lower.force_n > 0.0
            && (lower.force_n - target_force_n).abs() <= force_tolerance_n
            && upper.approach_m - lower.approach_m > controls.approach_tolerance_m
        {
            self.refuse_physical_plateau(before, recession_m, &lower, &progress)?;
            self.refuse_physical_plateau(before, recession_m, &upper, &progress)?;
            return Err(progress.unresolved());
        }
        let forces = [lower.force_n, upper.force_n];
        if self.load_trial_balanced(before, recession_m, &lower, &progress)? {
            return self.finish_load(
                before,
                lower,
                &progress,
                forces,
                FiniteGapPlasticLoadDisposition::ForceBalanced,
                cx,
            );
        }
        if self.load_trial_balanced(before, recession_m, &upper, &progress)? {
            return self.finish_load(
                before,
                upper,
                &progress,
                forces,
                FiniteGapPlasticLoadDisposition::ForceBalanced,
                cx,
            );
        }
        if target_force_n < lower.force_n || target_force_n > upper.force_n {
            return Err(FiniteGapPlasticLoadError::TargetOutsideApproachBounds {
                target_force_n,
                approach_bracket_m: progress.bracket_m,
                force_bracket_n: forces,
            });
        }
        let mut bisect_next = false;
        loop {
            checkpoint(cx)?;
            let width = upper.approach_m - lower.approach_m;
            if width <= controls.approach_tolerance_m || upper.force_n <= lower.force_n {
                return Err(progress.unresolved());
            }
            let midpoint = f64::midpoint(lower.approach_m, upper.approach_m);
            let fraction = (target_force_n - lower.force_n) / (upper.force_n - lower.force_n);
            let secant = fraction.mul_add(width, lower.approach_m);
            let approach = if !bisect_next
                && secant.is_finite()
                && secant > lower.approach_m
                && secant < upper.approach_m
            {
                secant
            } else {
                midpoint
            };
            if approach <= lower.approach_m || approach >= upper.approach_m {
                return Err(progress.unresolved());
            }
            let trial = self.load_trial(before, recession_m, approach, &mut progress, cx)?;
            if trial.force_n < lower.force_n - force_tolerance_n
                || trial.force_n > upper.force_n + force_tolerance_n
            {
                return Err(progress.unresolved());
            }
            let repeated_free_response = ((trial.force_n == lower.force_n
                && trial.approach_m - lower.approach_m > controls.approach_tolerance_m)
                || (trial.force_n == upper.force_n
                    && upper.approach_m - trial.approach_m > controls.approach_tolerance_m))
                && (trial.force_n - target_force_n).abs() <= force_tolerance_n;
            if repeated_free_response {
                self.refuse_physical_plateau(before, recession_m, &trial, &progress)?;
                return Err(progress.unresolved());
            }
            if self.load_trial_balanced(before, recession_m, &trial, &progress)? {
                return self.finish_load(
                    before,
                    trial,
                    &progress,
                    [lower.force_n, upper.force_n],
                    FiniteGapPlasticLoadDisposition::ForceBalanced,
                    cx,
                );
            }
            if trial.force_n < target_force_n {
                lower = trial;
            } else {
                upper = trial;
            }
            progress.bracket_m = [lower.approach_m, upper.approach_m];
            // At least every other stalled secant step becomes a bisection;
            // interior secants remain fast on the piecewise-linear contact law.
            bisect_next = upper.approach_m - lower.approach_m > 0.75 * width;
        }
    }

    fn load_trial(
        &self,
        before: &FiniteGapPlasticState,
        recession_m: &[f64],
        approach_m: f64,
        progress: &mut LoadProgress,
        cx: &Cx<'_>,
    ) -> Result<LoadTrial, FiniteGapPlasticLoadError> {
        if progress.evaluations >= progress.controls.maximum_equilibrium_evaluations {
            return Err(progress.unresolved());
        }
        progress.evaluations += 1;
        // A cold feasible pressure guess avoids accepting the OLD pressure
        // merely because its closure error is below the native tolerance.
        // Otherwise a small load/geometry change can create a false force
        // plateau. Plastic/recession/work history still comes only from before.
        let contact = self.advance_with_history(
            before,
            approach_m,
            recession_m,
            &progress.zero_pressure,
            cx,
        )?;
        let force_n = contact.receipt.normal.normal_force_n;
        let residual = (force_n - progress.target_force_n).abs();
        progress.best_residual_n = Some(
            progress
                .best_residual_n
                .map_or(residual, |old| old.min(residual)),
        );
        Ok(LoadTrial {
            approach_m,
            force_n,
            contact,
        })
    }

    fn load_trial_balanced(
        &self,
        before: &FiniteGapPlasticState,
        recession_m: &[f64],
        trial: &LoadTrial,
        progress: &LoadProgress,
    ) -> Result<bool, FiniteGapPlasticLoadError> {
        if trial.force_n <= 0.0
            || (trial.force_n - progress.target_force_n).abs() > progress.force_tolerance_n
        {
            return Ok(false);
        }
        if self.has_free_pressure(trial) {
            let capped_force = trial
                .contact
                .candidate
                .pressure_pa
                .iter()
                .zip(&self.config.hardness_pa)
                .filter_map(|(pressure, cap)| (*pressure == *cap).then_some(*pressure))
                .sum::<f64>()
                * (self.config.grid.cell_width_m * self.config.grid.cell_depth_m);
            // A loose force tolerance must not accept the far end of a plastic
            // plateau merely because a newly touching cell carries a tiny
            // positive load. Continue searching/refuse unresolved; free cells
            // alone do not prove nonuniqueness at this approximate response.
            return Ok(capped_force == 0.0 || capped_force != progress.target_force_n);
        }
        // Near a cap, a merely close force must not hide an unattainable target.
        if trial.force_n != progress.target_force_n {
            return Ok(false);
        }
        let interval = self.refuse_physical_plateau(before, recession_m, trial, progress)?;
        Ok(interval.is_some_and(|bounds| {
            trial.approach_m >= bounds[0] - progress.controls.approach_tolerance_m
                && trial.approach_m <= bounds[1] + progress.controls.approach_tolerance_m
        }))
    }

    fn has_free_pressure(&self, trial: &LoadTrial) -> bool {
        trial
            .contact
            .candidate
            .pressure_pa
            .iter()
            .zip(&self.config.hardness_pa)
            .any(|(pressure, cap)| *pressure > 0.0 && *pressure < *cap)
    }

    fn refuse_physical_plateau(
        &self,
        before: &FiniteGapPlasticState,
        recession_m: &[f64],
        trial: &LoadTrial,
        progress: &LoadProgress,
    ) -> Result<Option<[f64; 2]>, FiniteGapPlasticLoadError> {
        if trial.force_n != progress.target_force_n
            || trial.force_n == 0.0
            || self.has_free_pressure(trial)
        {
            return Ok(None);
        }
        let mut lower = f64::NEG_INFINITY;
        let mut upper = f64::INFINITY;
        for (i, pressure) in trial.contact.candidate.pressure_pa.iter().enumerate() {
            // The dual interval uses ORIGINAL plastic depth, not the candidate
            // dents that would vary arbitrarily across the force plateau.
            let onset = self.gap_m[i]
                + before.plastic_depth_m[i]
                + recession_m[i]
                + trial.contact.receipt.normal.elastic_opening_m[i];
            if *pressure == self.config.hardness_pa[i] {
                lower = lower.max(onset);
            } else if *pressure == 0.0 {
                upper = upper.min(onset);
            }
        }
        if !lower.is_finite() || upper.is_nan() {
            return Err(progress.unresolved());
        }
        if upper - lower > progress.controls.approach_tolerance_m {
            return Err(FiniteGapPlasticLoadError::NonUniqueLoad {
                target_force_n: progress.target_force_n,
                plateau_force_n: trial.force_n,
                plateau_approach_m: [lower, upper],
            });
        }
        Ok(Some([lower, upper]))
    }

    fn accepted_force(
        &self,
        state: &FiniteGapPlasticState,
    ) -> Result<f64, FiniteGapPlasticLoadError> {
        let area = self.config.grid.cell_width_m * self.config.grid.cell_depth_m;
        let force = state.pressure_pa.iter().sum::<f64>() * area;
        if !force.is_finite() {
            return Err(invalid_load("finite accepted normal force"));
        }
        Ok(force)
    }

    fn interior_force_capacity_bound(&self) -> f64 {
        let grid = self.config.grid;
        let ring = self.config.boundary_clearance_cells;
        let area = grid.cell_width_m * grid.cell_depth_m;
        let mut bound = 0.0;
        for (i, hardness) in self.config.hardness_pa.iter().enumerate() {
            let x = i % grid.cells_x;
            let y = i / grid.cells_x;
            if x >= ring && y >= ring && x < grid.cells_x - ring && y < grid.cells_y - ring {
                // Round each positive contribution and addition upward so this
                // necessary capacity bound cannot exclude a force by roundoff.
                bound = (bound + (hardness * area).next_up()).next_up();
            }
        }
        bound
    }

    fn finish_load(
        &self,
        before: &FiniteGapPlasticState,
        trial: LoadTrial,
        progress: &LoadProgress,
        force_bracket_n: [f64; 2],
        disposition: FiniteGapPlasticLoadDisposition,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLoadAdvance, FiniteGapPlasticLoadError> {
        let force_residual_n = trial.force_n - progress.target_force_n;
        let mut bytes = Vec::with_capacity(256);
        for identity in [
            self.identity,
            before.identity,
            trial.contact.candidate.identity,
            trial.contact.receipt.identity,
        ] {
            bytes.extend_from_slice(identity.as_bytes());
        }
        for value in [
            progress.target_force_n,
            trial.force_n,
            force_residual_n,
            progress.bracket_m[0],
            progress.bracket_m[1],
            force_bracket_n[0],
            force_bracket_n[1],
            progress.controls.minimum_approach_m,
            progress.controls.maximum_approach_m,
            progress.controls.absolute_force_tolerance_n,
            progress.controls.relative_force_tolerance,
            progress.controls.approach_tolerance_m,
        ] {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        bytes.extend_from_slice(&(progress.evaluations as u64).to_le_bytes());
        bytes.extend_from_slice(
            &(progress.controls.maximum_equilibrium_evaluations as u64).to_le_bytes(),
        );
        bytes.push(disposition as u8);
        let identity = hash_domain(
            "org.frankensim.fs-contact.finite-gap-plastic-load.v1",
            &bytes,
        );
        checkpoint(cx)?;
        Ok(FiniteGapPlasticLoadAdvance {
            load: FiniteGapPlasticLoadReceipt {
                identity,
                model_identity: self.identity,
                predecessor_state_identity: before.identity,
                candidate_state_identity: trial.contact.candidate.identity,
                response_identity: trial.contact.receipt.identity,
                target_force_n: progress.target_force_n,
                achieved_force_n: trial.force_n,
                force_residual_n,
                approach_bracket_m: progress.bracket_m,
                force_bracket_n,
                equilibrium_evaluations: progress.evaluations,
                controls: progress.controls,
                disposition,
            },
            contact: trial.contact,
        })
    }
}

fn validate_controls(
    force_n: f64,
    controls: FiniteGapPlasticLoadControls,
) -> Result<f64, FiniteGapPlasticLoadError> {
    if !force_n.is_finite() || force_n < 0.0 {
        return Err(invalid_load("nonnegative finite target force"));
    }
    if !controls.minimum_approach_m.is_finite()
        || !controls.maximum_approach_m.is_finite()
        || controls.minimum_approach_m > controls.maximum_approach_m
        || !(controls.maximum_approach_m - controls.minimum_approach_m).is_finite()
    {
        return Err(invalid_load("finite ordered approach bounds"));
    }
    if !controls.absolute_force_tolerance_n.is_finite()
        || controls.absolute_force_tolerance_n <= 0.0
        || !controls.relative_force_tolerance.is_finite()
        || controls.relative_force_tolerance < 0.0
        || controls.relative_force_tolerance >= 1.0
        || !controls.approach_tolerance_m.is_finite()
        || controls.approach_tolerance_m <= 0.0
        || controls.maximum_equilibrium_evaluations == 0
    {
        return Err(invalid_load(
            "positive finite load tolerances and work bound",
        ));
    }
    let tolerance =
        controls.absolute_force_tolerance_n + controls.relative_force_tolerance * force_n;
    if !tolerance.is_finite() {
        return Err(invalid_load("finite combined force tolerance"));
    }
    Ok(tolerance)
}

fn invalid_load(field: &'static str) -> FiniteGapPlasticLoadError {
    FiniteGapPlasticLoadError::InvalidInput { field }
}
