//! A finite normal loading spring in series with hardness-limited contact.
//!
//! The equilibrium is `F_contact(delta) = K * (position - delta)`, with finite
//! positive `K`. Its residual increases strictly with approach even on a fully
//! plastic contact-force plateau, so the loading device determines the missing
//! displacement. A detached equilibrium has zero spring compression and lies
//! exactly at the prescribed drive position. Coordinates retain the original
//! gap datum; the loading spring does not change the indentation law.
//!
//! Every pressure trial starts from the original accepted irreversible history.
//! Recession at a fixed drive position exchanges work with the spring; its
//! normal work is minus the spring-storage change, separate from geometry work.
//! A change of drive position returns an endpoint equilibrium, not an inferred
//! actuator-work integral, inertia, damping, or a dynamic loading-machine model.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::Cx;

use super::{
    FiniteGapPlasticAdvance, FiniteGapPlasticContinuation, FiniteGapPlasticError,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError, FiniteGapPlasticModel,
    FiniteGapPlasticState, checkpoint, load::validate_controls,
};

/// Linear normal loading device, with position in the original approach chart.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticDrive {
    /// Finite positive series stiffness [N/m]. No infinite-stiffness sentinel.
    pub stiffness_n_per_m: f64,
    /// Prescribed unloaded position relative to the original gap minimum [m].
    pub position_m: f64,
}

/// How the loading-device equilibrium was obtained.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum FiniteGapPlasticDriveDisposition {
    /// Contact and spring force balance within the reported tolerance.
    ForceBalanced = 0,
    /// The unchanged accepted contact already satisfies this drive boundary.
    RetainedExistingEquilibrium = 1,
}

/// Contact history selected by the finite loading-device stiffness.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticDriveAdvance {
    /// Original native contact response, without accepted search-trial history.
    pub contact: FiniteGapPlasticAdvance,
    /// Loading-device equilibrium and spring storage for the same candidate.
    pub drive: FiniteGapPlasticDriveReceipt,
}

/// Recession while the original loading spring and its position remain fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticDriveRecessionAdvance {
    /// Contact continuation; its normal work comes from the loading spring.
    pub continuation: FiniteGapPlasticContinuation,
    /// Spring force, storage, and balance at the continued equilibrium.
    pub drive: FiniteGapPlasticDriveReceipt,
}

/// Explicit series-spring boundary, balance, storage, and numerical work.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticDriveReceipt {
    /// Content identity of the boundary, controls, and original native response.
    pub identity: ContentHash,
    /// Immutable contact model used by every trial.
    pub model_identity: ContentHash,
    /// Original accepted history used by every trial.
    pub predecessor_state_identity: ContentHash,
    /// Selected physical candidate.
    pub candidate_state_identity: ContentHash,
    /// Original normal/plastic response identity.
    pub response_identity: ContentHash,
    /// Finite-stiffness loading boundary used for this equilibrium.
    pub boundary: FiniteGapPlasticDrive,
    /// Selected rigid approach in the original gap chart [m].
    pub approach_m: f64,
    /// Integrated force of the native candidate pressure field [N].
    pub achieved_force_n: f64,
    /// Compressive spring force `K * (position - approach)` [N].
    pub drive_force_n: f64,
    /// Signed contact force minus spring force [N].
    pub force_residual_n: f64,
    /// Accepted residual threshold: absolute plus relative times the larger
    /// magnitude of the two endpoint forces [N].
    pub force_tolerance_n: f64,
    /// Reversible storage in the loading spring only [J].
    pub spring_energy_j: f64,
    /// Final tested approach interval; an exact retained state reports a point.
    /// The upper bound is limited to the unloaded drive position. These are
    /// numerical search bounds, not a rigorous physical position enclosure.
    pub approach_bracket_m: [f64; 2],
    /// Signed contact-minus-spring residuals at the reported tested bounds [N].
    pub balance_residual_bracket_n: [f64; 2],
    /// Number of native pressure equilibria attempted, including rejected trials.
    pub equilibrium_evaluations: usize,
    /// Original numerical controls. The relative force tolerance uses the
    /// larger contact/spring force rather than a prescribed constant load.
    pub controls: FiniteGapPlasticLoadControls,
    /// Actual position-selection convention.
    pub disposition: FiniteGapPlasticDriveDisposition,
}

/// Refusal without adopting any pressure-search or recession history.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGapPlasticDriveError {
    /// Invalid loading spring, search controls, recession, or finite arithmetic.
    InvalidInput { field: &'static str },
    /// Original contact ownership, geometry, convergence, or cancellation refusal.
    Contact(FiniteGapPlasticError),
    /// The admitted approach bounds do not contain a spring/contact balance.
    TargetOutsideApproachBounds {
        approach_bracket_m: [f64; 2],
        balance_residual_bracket_n: [f64; 2],
    },
    /// Numerical resolution or the work budget expired before force balance.
    DidNotConverge {
        equilibrium_evaluations: usize,
        approach_bracket_m: [f64; 2],
        best_force_residual_n: Option<f64>,
    },
    /// Fixed-drive recession requires an already equilibrated predecessor.
    PredecessorDriveMismatch {
        contact_force_n: f64,
        drive_force_n: f64,
        force_tolerance_n: f64,
    },
}

impl From<FiniteGapPlasticError> for FiniteGapPlasticDriveError {
    fn from(error: FiniteGapPlasticError) -> Self {
        Self::Contact(error)
    }
}

impl fmt::Display for FiniteGapPlasticDriveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Contact(error) => write!(f, "{error}"),
            Self::InvalidInput { field } => write!(f, "invalid finite-gap normal drive {field}"),
            Self::TargetOutsideApproachBounds {
                approach_bracket_m,
                balance_residual_bracket_n,
            } => write!(
                f,
                "normal drive has no admitted balance in {approach_bracket_m:?} m; contact-minus-spring residuals {balance_residual_bracket_n:?} N"
            ),
            Self::DidNotConverge {
                equilibrium_evaluations,
                approach_bracket_m,
                best_force_residual_n,
            } => write!(
                f,
                "normal drive unresolved after {equilibrium_evaluations} equilibria in {approach_bracket_m:?} m; best force residual {best_force_residual_n:?} N"
            ),
            Self::PredecessorDriveMismatch {
                contact_force_n,
                drive_force_n,
                force_tolerance_n,
            } => write!(
                f,
                "accepted force {contact_force_n} N differs from fixed-drive force {drive_force_n} N beyond tolerance {force_tolerance_n} N"
            ),
        }
    }
}

impl core::error::Error for FiniteGapPlasticDriveError {}

struct DriveTrial {
    contact: FiniteGapPlasticAdvance,
    drive_force_n: f64,
    residual_n: f64,
    tolerance_n: f64,
}

impl DriveTrial {
    fn balanced(&self) -> bool {
        balanced(
            self.contact.receipt.normal.normal_force_n,
            self.drive_force_n,
            self.residual_n,
            self.tolerance_n,
        )
    }

    fn approach_m(&self) -> f64 {
        self.contact.candidate.approach_m
    }
}

struct DriveProgress {
    boundary: FiniteGapPlasticDrive,
    controls: FiniteGapPlasticLoadControls,
    evaluations: usize,
    bracket_m: [f64; 2],
    best_residual_n: Option<f64>,
    zero_pressure: Vec<f64>,
}

impl DriveProgress {
    fn unresolved(&self) -> FiniteGapPlasticDriveError {
        FiniteGapPlasticDriveError::DidNotConverge {
            equilibrium_evaluations: self.evaluations,
            approach_bracket_m: self.bracket_m,
            best_force_residual_n: self.best_residual_n,
        }
    }

    fn outside(&self, residuals: [f64; 2]) -> FiniteGapPlasticDriveError {
        FiniteGapPlasticDriveError::TargetOutsideApproachBounds {
            approach_bracket_m: self.bracket_m,
            balance_residual_bracket_n: residuals,
        }
    }
}

impl FiniteGapPlasticModel {
    /// Equilibrate contact against a finite-stiffness normal loading device.
    ///
    /// The contact-force plateau is allowed: positive loading stiffness makes
    /// the combined residual strictly increasing. Each trial uses original
    /// plastic depth, recession, and accumulated work. No actuator-work path
    /// is inferred when the prescribed loading-device position changes.
    pub fn advance_with_drive(
        &self,
        before: &FiniteGapPlasticState,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticDriveAdvance, FiniteGapPlasticDriveError> {
        self.drive_at_recession(before, &before.recession_m, drive, controls, cx)
    }

    /// Remove material while holding this loading spring and its position fixed.
    ///
    /// The predecessor must already balance the requested drive. The work
    /// delivered to contact is `mean(spring force) * delta(approach)`, exactly
    /// minus the spring-energy change for a linear spring at fixed position.
    /// It is an internal transfer when contact and spring share the accounting
    /// boundary; the held actuator does zero work. Geometry parameter work is
    /// still separate, and neither account is automatically assigned to heat.
    pub fn recede_with_drive(
        &self,
        before: &FiniteGapPlasticState,
        recession_increment_m: &[f64],
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticDriveRecessionAdvance, FiniteGapPlasticDriveError> {
        checkpoint(cx)?;
        self.check_state_owner(before)?;
        validate_drive(drive, controls)?;
        let previous_force = self.drive_accepted_force(before, cx)?;
        let previous_spring_force = spring_force(drive, before.approach_m)?;
        let tolerance = force_tolerance(previous_force, previous_spring_force, controls)?;
        if !balanced(
            previous_force,
            previous_spring_force,
            previous_force - previous_spring_force,
            tolerance,
        ) {
            return Err(FiniteGapPlasticDriveError::PredecessorDriveMismatch {
                contact_force_n: previous_force,
                drive_force_n: previous_spring_force,
                force_tolerance_n: tolerance,
            });
        }
        if recession_increment_m.len() != self.gap_m.len() {
            return Err(invalid_drive("recession field shape"));
        }
        let mut recession = before.recession_m.clone();
        let mut actual_increment = vec![0.0; recession.len()];
        for (i, requested) in recession_increment_m.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            if !requested.is_finite() || *requested < 0.0 {
                return Err(invalid_drive("finite nonnegative recession field"));
            }
            recession[i] += requested;
            if !recession[i].is_finite() {
                return Err(invalid_drive("finite accumulated recession"));
            }
            actual_increment[i] = recession[i] - before.recession_m[i];
            if *requested > 0.0 && actual_increment[i] == 0.0 {
                return Err(invalid_drive("representable recession increment"));
            }
        }
        let next = self.drive_at_recession(before, &recession, drive, controls, cx)?;
        let mean_spring_force = 0.5 * previous_spring_force + 0.5 * next.drive.drive_force_n;
        let normal_work_j =
            mean_spring_force * (next.contact.candidate.approach_m - before.approach_m);
        if !normal_work_j.is_finite() {
            return Err(invalid_drive("finite fixed-drive spring work"));
        }
        let continuation = self.finish_boundary_recession(
            before,
            next.contact,
            actual_increment,
            normal_work_j,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(FiniteGapPlasticDriveRecessionAdvance {
            continuation,
            drive: next.drive,
        })
    }

    fn drive_at_recession(
        &self,
        before: &FiniteGapPlasticState,
        recession_m: &[f64],
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticDriveAdvance, FiniteGapPlasticDriveError> {
        checkpoint(cx)?;
        self.check_state_owner(before)?;
        validate_drive(drive, controls)?;
        // Nonadhesive equilibrium cannot stretch the spring into tension. If
        // the entire admitted interval is above its unloaded position, test
        // the lower endpoint and report that one-sided refusal.
        let upper_approach = controls
            .maximum_approach_m
            .min(drive.position_m)
            .max(controls.minimum_approach_m);
        let mut progress = DriveProgress {
            boundary: drive,
            controls,
            evaluations: 0,
            bracket_m: [controls.minimum_approach_m, upper_approach],
            best_residual_n: None,
            zero_pressure: vec![0.0; self.gap_m.len()],
        };
        if recession_m == before.recession_m.as_slice()
            && before.approach_m >= progress.bracket_m[0]
            && before.approach_m <= progress.bracket_m[1]
        {
            let force = self.drive_accepted_force(before, cx)?;
            let drive_force_n = spring_force(drive, before.approach_m)?;
            let tolerance_n = force_tolerance(force, drive_force_n, controls)?;
            let residual_n = force - drive_force_n;
            if balanced(force, drive_force_n, residual_n, tolerance_n) {
                let (_, opening) = self.stored_energy_and_opening(before, cx)?;
                let contact = self.unchanged_advance(before, opening, cx)?;
                progress.bracket_m = [before.approach_m; 2];
                return self.finish_drive(
                    before,
                    DriveTrial {
                        contact,
                        drive_force_n,
                        residual_n,
                        tolerance_n,
                    },
                    &progress,
                    [residual_n; 2],
                    FiniteGapPlasticDriveDisposition::RetainedExistingEquilibrium,
                    cx,
                );
            }
        }
        let mut lower = self.drive_trial(
            before,
            recession_m,
            progress.bracket_m[0],
            &mut progress,
            cx,
        )?;
        if lower.balanced() {
            progress.bracket_m = [lower.approach_m(); 2];
            let residuals = [lower.residual_n; 2];
            return self.finish_drive(
                before,
                lower,
                &progress,
                residuals,
                FiniteGapPlasticDriveDisposition::ForceBalanced,
                cx,
            );
        }
        if lower.residual_n >= 0.0 || upper_approach == controls.minimum_approach_m {
            progress.bracket_m = [lower.approach_m(); 2];
            return Err(progress.outside([lower.residual_n; 2]));
        }
        let mut upper = self.drive_trial(before, recession_m, upper_approach, &mut progress, cx)?;
        if upper.balanced() {
            let residuals = [lower.residual_n, upper.residual_n];
            return self.finish_drive(
                before,
                upper,
                &progress,
                residuals,
                FiniteGapPlasticDriveDisposition::ForceBalanced,
                cx,
            );
        }
        if upper.residual_n < 0.0 {
            return Err(progress.outside([lower.residual_n, upper.residual_n]));
        }
        let mut bisect_next = false;
        loop {
            checkpoint(cx)?;
            let width = upper.approach_m() - lower.approach_m();
            if width <= controls.approach_tolerance_m {
                return Err(progress.unresolved());
            }
            // Scale opposite-sign residuals before forming a secant fraction;
            // their unscaled difference could overflow despite finite endpoints.
            let scale = (-lower.residual_n).max(upper.residual_n);
            let left = -lower.residual_n / scale;
            let right = upper.residual_n / scale;
            let fraction = left / (left + right);
            let secant = fraction.mul_add(width, lower.approach_m());
            let approach = if !bisect_next
                && secant.is_finite()
                && secant > lower.approach_m()
                && secant < upper.approach_m()
            {
                secant
            } else {
                f64::midpoint(lower.approach_m(), upper.approach_m())
            };
            if approach <= lower.approach_m() || approach >= upper.approach_m() {
                return Err(progress.unresolved());
            }
            let trial = self.drive_trial(before, recession_m, approach, &mut progress, cx)?;
            if (trial.residual_n < lower.residual_n
                && lower.residual_n - trial.residual_n > lower.tolerance_n.max(trial.tolerance_n))
                || (trial.residual_n > upper.residual_n
                    && trial.residual_n - upper.residual_n
                        > upper.tolerance_n.max(trial.tolerance_n))
            {
                return Err(progress.unresolved());
            }
            if trial.balanced() {
                let residuals = [lower.residual_n, upper.residual_n];
                return self.finish_drive(
                    before,
                    trial,
                    &progress,
                    residuals,
                    FiniteGapPlasticDriveDisposition::ForceBalanced,
                    cx,
                );
            }
            if trial.residual_n < 0.0 {
                lower = trial;
            } else {
                upper = trial;
            }
            progress.bracket_m = [lower.approach_m(), upper.approach_m()];
            bisect_next = upper.approach_m() - lower.approach_m() > 0.75 * width;
        }
    }

    fn drive_trial(
        &self,
        before: &FiniteGapPlasticState,
        recession_m: &[f64],
        approach_m: f64,
        progress: &mut DriveProgress,
        cx: &Cx<'_>,
    ) -> Result<DriveTrial, FiniteGapPlasticDriveError> {
        if progress.evaluations >= progress.controls.maximum_equilibrium_evaluations {
            return Err(progress.unresolved());
        }
        progress.evaluations += 1;
        // Reuse a cold feasible pressure vector: a warm guess within closure
        // tolerance can otherwise freeze force across a small approach change.
        let contact = self.advance_with_history(
            before,
            approach_m,
            recession_m,
            &progress.zero_pressure,
            cx,
        )?;
        let force = contact.receipt.normal.normal_force_n;
        let drive_force_n = spring_force(progress.boundary, approach_m)?;
        let residual_n = force - drive_force_n;
        if !residual_n.is_finite() {
            return Err(invalid_drive("finite contact-minus-spring residual"));
        }
        let tolerance_n = force_tolerance(force, drive_force_n, progress.controls)?;
        progress.best_residual_n = Some(
            progress
                .best_residual_n
                .map_or(residual_n.abs(), |old| old.min(residual_n.abs())),
        );
        Ok(DriveTrial {
            contact,
            drive_force_n,
            residual_n,
            tolerance_n,
        })
    }

    fn drive_accepted_force(
        &self,
        before: &FiniteGapPlasticState,
        cx: &Cx<'_>,
    ) -> Result<f64, FiniteGapPlasticDriveError> {
        let mut pressure_sum = 0.0;
        for pressures in before.pressure_pa.chunks(64) {
            checkpoint(cx)?;
            for pressure in pressures {
                pressure_sum += pressure;
            }
        }
        let force = pressure_sum * (self.config.grid.cell_width_m * self.config.grid.cell_depth_m);
        if !force.is_finite() {
            return Err(invalid_drive("finite accepted normal force"));
        }
        Ok(force)
    }

    fn finish_drive(
        &self,
        before: &FiniteGapPlasticState,
        trial: DriveTrial,
        progress: &DriveProgress,
        residuals: [f64; 2],
        disposition: FiniteGapPlasticDriveDisposition,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticDriveAdvance, FiniteGapPlasticDriveError> {
        let approach_m = trial.approach_m();
        let force = trial.contact.receipt.normal.normal_force_n;
        let spring_energy_j =
            (0.5 * trial.drive_force_n) * (progress.boundary.position_m - approach_m);
        if !spring_energy_j.is_finite() || spring_energy_j < 0.0 {
            return Err(invalid_drive("finite nonnegative loading-spring energy"));
        }
        let mut bytes = Vec::with_capacity(320);
        for identity in [
            self.identity,
            before.identity,
            trial.contact.candidate.identity,
            trial.contact.receipt.identity,
        ] {
            bytes.extend_from_slice(identity.as_bytes());
        }
        for value in [
            progress.boundary.stiffness_n_per_m,
            progress.boundary.position_m,
            approach_m,
            force,
            trial.drive_force_n,
            trial.residual_n,
            trial.tolerance_n,
            spring_energy_j,
            progress.bracket_m[0],
            progress.bracket_m[1],
            residuals[0],
            residuals[1],
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
            "org.frankensim.fs-contact.finite-gap-plastic-drive.v1",
            &bytes,
        );
        checkpoint(cx)?;
        Ok(FiniteGapPlasticDriveAdvance {
            drive: FiniteGapPlasticDriveReceipt {
                identity,
                model_identity: self.identity,
                predecessor_state_identity: before.identity,
                candidate_state_identity: trial.contact.candidate.identity,
                response_identity: trial.contact.receipt.identity,
                boundary: progress.boundary,
                approach_m,
                achieved_force_n: force,
                drive_force_n: trial.drive_force_n,
                force_residual_n: trial.residual_n,
                force_tolerance_n: trial.tolerance_n,
                spring_energy_j,
                approach_bracket_m: progress.bracket_m,
                balance_residual_bracket_n: residuals,
                equilibrium_evaluations: progress.evaluations,
                controls: progress.controls,
                disposition,
            },
            contact: trial.contact,
        })
    }
}

fn validate_drive(
    drive: FiniteGapPlasticDrive,
    controls: FiniteGapPlasticLoadControls,
) -> Result<(), FiniteGapPlasticDriveError> {
    validate_controls(0.0, controls).map_err(|error| match error {
        FiniteGapPlasticLoadError::InvalidInput { field } => invalid_drive(field),
        _ => invalid_drive("finite force-search controls"),
    })?;
    if !drive.stiffness_n_per_m.is_finite()
        || drive.stiffness_n_per_m <= 0.0
        || !drive.position_m.is_finite()
    {
        return Err(invalid_drive(
            "finite positive stiffness and finite position",
        ));
    }
    Ok(())
}

fn spring_force(
    drive: FiniteGapPlasticDrive,
    approach_m: f64,
) -> Result<f64, FiniteGapPlasticDriveError> {
    let force = drive.stiffness_n_per_m * (drive.position_m - approach_m);
    if !force.is_finite() {
        return Err(invalid_drive("finite loading-spring force"));
    }
    Ok(force)
}

fn force_tolerance(
    contact_force: f64,
    drive_force: f64,
    controls: FiniteGapPlasticLoadControls,
) -> Result<f64, FiniteGapPlasticDriveError> {
    let tolerance = controls.absolute_force_tolerance_n
        + controls.relative_force_tolerance * contact_force.abs().max(drive_force.abs());
    if !tolerance.is_finite() {
        return Err(invalid_drive("finite combined force tolerance"));
    }
    Ok(tolerance)
}

fn balanced(contact_force: f64, drive_force: f64, residual: f64, tolerance: f64) -> bool {
    drive_force >= 0.0 && residual.abs() <= tolerance
        // A detached spring/contact balance has a unique zero-compression
        // position; a small force tolerance must not invent a slack interval.
        && (contact_force > 0.0 || drive_force == 0.0)
}

fn invalid_drive(field: &'static str) -> FiniteGapPlasticDriveError {
    FiniteGapPlasticDriveError::InvalidInput { field }
}
