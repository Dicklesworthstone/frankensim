//! Original attached Maugis potential for native mechanical evolution.
//!
//! This owner explicitly selects the stable contacting branch. An opt-in
//! extension follows its continuous zero-solid-area junction into the original
//! deformed noncontact adhesive disk, when that branch exists. Neither mode
//! crosses its terminal fold or turns quasistatic release into an inertial
//! impulse. The retained original source supplies all forces and absolute energy.

use super::continuation::{
    MaugisContinuationBudget, MaugisContinuationError, MaugisEquilibrium, MaugisNoncontactLimits,
};
use super::{MaugisError, MaugisWorkEstimate, SourceMaugisSphereContact};
use core::f64::consts::PI;
use fs_blake3::ContentHash;

/// Explicit attached-branch choice with bounded original-force integration.
#[derive(Clone, Debug, PartialEq)]
pub struct MaugisAttachedPotential {
    source: SourceMaugisSphereContact,
    budget: MaugisContinuationBudget,
    noncontact: Option<MaugisNoncontactLimits>,
}

/// Work and the conjugate mean force on one straight approach interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaugisMechanicalPath {
    /// Mean signed normal force, N; attraction is negative.
    pub mean_force_n: f64,
    /// Original force integral, J, and a numerical quadrature error estimate.
    /// Each integrand evaluation uses at most two bounded original inversions.
    pub work: MaugisWorkEstimate,
    /// Original endpoint-potential difference minus integrated work, J.
    pub energy_residual_j: f64,
}

impl MaugisAttachedPotential {
    /// Select attachment explicitly. This does not infer previous snap history.
    pub fn new(
        source: SourceMaugisSphereContact,
        budget: MaugisContinuationBudget,
    ) -> Result<Self, MaugisContinuationError> {
        if !budget.work.absolute_tolerance_j.is_finite()
            || budget.work.absolute_tolerance_j <= 0.0
            || !budget.work.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&budget.work.relative_tolerance)
            || budget.work.maximum_evaluations < 5
            || !budget.energy_absolute_tolerance_j.is_finite()
            || budget.energy_absolute_tolerance_j < 0.0
            || !budget.energy_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&budget.energy_relative_tolerance)
        {
            return Err(MaugisError::InvalidInput {
                field: "attached mechanical work budget",
            }
            .into());
        }
        Ok(Self {
            source,
            budget,
            noncontact: None,
        })
    }

    /// Explicitly allow loss and recovery of SOLID contact along the original
    /// locally stable connected branch. Across its zero-area junction, force
    /// and absolute energy are continuous; no snap, impulse or loss is added.
    /// The deformed adhesive disk retains attraction with zero friction area.
    ///
    /// Sources without a stable noncontact interval refuse. The lower disk
    /// fold still refuses further separation: a dynamical transition to the
    /// detached phase needs its own physical internal-motion/loss model.
    pub fn with_noncontact_continuation(mut self) -> Result<Self, MaugisContinuationError> {
        self.noncontact = Some(self.source.noncontact_limits()?.ok_or(
            MaugisError::InvalidInput {
                field: "source has no stable continuous noncontact branch",
            },
        )?);
        Ok(self)
    }

    /// Explicit branch policy and cached physical bounds. None retains the
    /// original contact-only admission; the original source identity is shared
    /// because values on the common contacting domain are unchanged.
    #[must_use]
    pub const fn noncontact_continuation(&self) -> Option<MaugisNoncontactLimits> {
        self.noncontact
    }
    /// Original immutable finite-zone law, including its material receipts.
    #[must_use]
    pub const fn source(&self) -> &SourceMaugisSphereContact {
        &self.source
    }
    /// Original model identity; this adapter never promotes source authority.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.source.identity()
    }
    /// Explicit work and independent energy admission policy.
    #[must_use]
    pub const fn budget(&self) -> MaugisContinuationBudget {
        self.budget
    }
    /// Source-scale displacement [m] for bounded native numerical probes.
    #[must_use]
    pub const fn indentation_scale_m(&self) -> f64 {
        self.source.indentation_scale_m
    }
    /// Original energy and force on the explicitly admitted connected branch.
    /// Its terminal endpoint is finite; further separation refuses. Selecting
    /// the continuous noncontact extension never selects the detached phase.
    pub fn evaluate(
        &self,
        indentation_m: f64,
    ) -> Result<MaugisEquilibrium, MaugisContinuationError> {
        if let Some(limits) = self.noncontact
            && indentation_m < limits.first_contact_indentation_m
        {
            self.source.evaluate_noncontact(indentation_m)
        } else {
            self.source.evaluate_contact_energy(indentation_m)
        }
    }
    /// Integrate the actual force along the commanded straight approach path.
    /// Symmetric, centered quadrature preserves a resolved mean effort even
    /// when endpoint-energy subtraction would lose all significant digits.
    /// Numerical work estimates are not certified interval bounds.
    pub fn path(
        &self,
        from_m: f64,
        to_m: f64,
    ) -> Result<MaugisMechanicalPath, MaugisContinuationError> {
        if let Some(limits) = self.noncontact
            && from_m.min(to_m) < limits.first_contact_indentation_m
        {
            return self.connected_path(from_m, to_m, limits.first_contact_indentation_m);
        }
        self.contact_path(from_m, to_m)
    }

    fn contact_path(
        &self,
        from_m: f64,
        to_m: f64,
    ) -> Result<MaugisMechanicalPath, MaugisContinuationError> {
        let from = self.evaluate(from_m)?;
        if from_m == to_m {
            return Ok(MaugisMechanicalPath {
                mean_force_n: from.normal_force_n(),
                work: MaugisWorkEstimate {
                    work_j: 0.0,
                    estimated_absolute_error_j: 0.0,
                    evaluations: 0,
                },
                energy_residual_j: 0.0,
            });
        }
        let to = self.evaluate(to_m)?;
        let lo = from_m.min(to_m);
        let hi = from_m.max(to_m);
        let center = lo.midpoint(hi);
        let half = 0.5 * (hi - lo);
        let delta = to_m - from_m;
        let central_force = self.source.evaluate_attached(center)?.normal_force_n;
        let force_scale = self.source.force_scale_n;
        let work_scale = force_scale * delta.abs();
        if !work_scale.is_finite() || work_scale <= 0.0 {
            return Err(MaugisError::Unrepresentable {
                field: "attached mechanical work scale",
            }
            .into());
        }
        let absolute = self.budget.work.absolute_tolerance_j / work_scale;
        let relative = self.budget.work.relative_tolerance;
        let (correction, error, evaluations) = super::work::integrate_interval(
            0.0,
            1.0,
            absolute + relative * (central_force / force_scale).abs(),
            relative,
            self.budget.work.maximum_evaluations,
            |u| {
                if u == 0.0 {
                    return Ok(0.0);
                }
                let (left, right) = if u == 1.0 {
                    (from.normal_force_n(), to.normal_force_n())
                } else {
                    (
                        self.source
                            .evaluate_attached(center - u * half)?
                            .normal_force_n,
                        self.source
                            .evaluate_attached(center + u * half)?
                            .normal_force_n,
                    )
                };
                Ok(0.5
                    * ((left - central_force) / force_scale
                        + (right - central_force) / force_scale))
            },
        )?;
        let mean_force_n = central_force + correction * force_scale;
        self.finish_path(
            from,
            to,
            delta,
            mean_force_n,
            (error * work_scale, evaluations),
        )
    }

    fn connected_path(
        &self,
        from_m: f64,
        to_m: f64,
        touch_m: f64,
    ) -> Result<MaugisMechanicalPath, MaugisContinuationError> {
        // Admit both original endpoints before integrating. A failed trial
        // cannot extrapolate past the noncontact fold or switch to detachment.
        let from = self.evaluate(from_m)?;
        let to = self.evaluate(to_m)?;
        let delta = to_m - from_m;
        if delta == 0.0 {
            return self.finish_path(from, to, 0.0, from.normal_force_n(), (0.0, 0));
        }
        let lo = from_m.min(to_m);
        let hi = from_m.max(to_m);
        let disk_end = hi.min(touch_m);
        let disk_a = self.source.evaluate_noncontact(lo)?;
        let disk_b = self.source.evaluate_noncontact(disk_end)?;
        let mean_disk = self.noncontact_mean_force(disk_a, disk_b)?;
        if hi <= touch_m {
            return self.finish_path(from, to, delta, mean_disk, (0.0, 0));
        }
        // Split at the genuine continuous zero-area junction, so adaptive
        // quadrature only sees the original contacting branch. The disk part
        // is exact and consumes no quadrature evaluations or hidden budget.
        let contact = self.contact_path(touch_m, hi)?;
        let span = hi - lo;
        let mean_force_n =
            mean_disk * ((touch_m - lo) / span) + contact.mean_force_n * ((hi - touch_m) / span);
        self.finish_path(
            from,
            to,
            delta,
            mean_force_n,
            (
                contact.work.estimated_absolute_error_j,
                contact.work.evaluations,
            ),
        )
    }

    // On the original disk, d(delta)/dc = d0*(c-c_fold) and
    // F/F0 = -pi*lambda*c^2/2. Integrate in x=c-c_fold, then cancel
    // the common interval length before division. Every remaining polynomial
    // term is nonnegative, including arbitrarily short paths near the fold.
    fn noncontact_mean_force(
        &self,
        a: MaugisEquilibrium,
        b: MaugisEquilibrium,
    ) -> Result<f64, MaugisContinuationError> {
        let fold = (8.0 / (3.0 * PI)) * self.source.lambda;
        let fold_radius = fold * self.source.radius_scale_m;
        let x = (a.cohesive_radius_m() - fold_radius) / self.source.radius_scale_m;
        let y = (b.cohesive_radius_m() - fold_radius) / self.source.radius_scale_m;
        let sum = x + y;
        let second_over_first = if sum == 0.0 {
            0.0
        } else {
            (2.0 / 3.0) * (x * x + x * y + y * y) / sum
        };
        let shape = fold * fold + 2.0 * fold * second_over_first + 0.5 * (x * x + y * y);
        let force = -((PI / 2.0) * self.source.lambda * self.source.force_scale_n) * shape;
        if !force.is_finite() || force >= 0.0 || x < 0.0 || y < 0.0 {
            return Err(MaugisError::Unrepresentable {
                field: "noncontact mechanical mean effort",
            }
            .into());
        }
        Ok(force)
    }

    fn finish_path(
        &self,
        from: MaugisEquilibrium,
        to: MaugisEquilibrium,
        delta: f64,
        mean_force_n: f64,
        (estimated_absolute_error_j, evaluations): (f64, usize),
    ) -> Result<MaugisMechanicalPath, MaugisContinuationError> {
        let work_j = mean_force_n * delta;
        if delta != 0.0 && mean_force_n != 0.0 && work_j == 0.0 {
            return Err(MaugisError::Unrepresentable {
                field: "mechanical path work increment",
            }
            .into());
        }
        let energy_residual_j = (to.potential_j() - from.potential_j()) - work_j;
        let tolerance = self.budget.energy_absolute_tolerance_j
            + self.budget.energy_relative_tolerance
                * from
                    .potential_j()
                    .abs()
                    .max(to.potential_j().abs())
                    .max(work_j.abs());
        if [
            mean_force_n,
            work_j,
            energy_residual_j,
            tolerance,
            estimated_absolute_error_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || energy_residual_j.abs() > tolerance
        {
            return Err(MaugisContinuationError::EnergyClosure {
                residual_j: energy_residual_j,
                tolerance_j: tolerance,
            });
        }
        Ok(MaugisMechanicalPath {
            mean_force_n,
            work: MaugisWorkEstimate {
                work_j,
                estimated_absolute_error_j,
                evaluations,
            },
            energy_residual_j,
        })
    }
    /// Complete original energy remainder for the supplied mechanical midpoint.
    /// The supplied midpoint can differ from the scalar endpoint average because
    /// it comes through an affine multi-body map. Its actual force is retained.
    pub fn midpoint_energy_remainder(
        &self,
        from_m: f64,
        to_m: f64,
        midpoint_m: f64,
    ) -> Result<f64, MaugisContinuationError> {
        let midpoint = self.evaluate(midpoint_m)?;
        let path = self.path(from_m, to_m)?;
        let value = (path.mean_force_n - midpoint.normal_force_n()) * (to_m - from_m);
        if !value.is_finite() {
            return Err(MaugisError::Unrepresentable {
                field: "attached mechanical energy remainder",
            }
            .into());
        }
        Ok(value)
    }
}
