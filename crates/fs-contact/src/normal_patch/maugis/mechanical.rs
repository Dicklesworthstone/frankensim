//! Original attached Maugis potential for native mechanical evolution.
//!
//! This owner explicitly selects the stable contacting branch. It does not
//! continue through its fold, select a noncontact state, or turn quasistatic
//! release into an inertial impulse. Those changes need an additional physical
//! owner. The retained original source supplies all forces and absolute energy.

use super::continuation::{MaugisContinuationBudget, MaugisContinuationError, MaugisEquilibrium};
use super::{MaugisError, MaugisWorkEstimate, SourceMaugisSphereContact};
use fs_blake3::ContentHash;

/// Explicit attached-branch choice with bounded original-force integration.
#[derive(Clone, Debug, PartialEq)]
pub struct MaugisAttachedPotential {
    source: SourceMaugisSphereContact,
    budget: MaugisContinuationBudget,
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
        Ok(Self { source, budget })
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
    /// Original attached energy and force, including the finite branch endpoint.
    /// A request beyond that endpoint refuses rather than switching phase.
    pub fn evaluate(
        &self,
        indentation_m: f64,
    ) -> Result<MaugisEquilibrium, MaugisContinuationError> {
        self.source.evaluate_contact_energy(indentation_m)
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
        let work_j = mean_force_n * delta;
        let energy_residual_j = (to.potential_j() - from.potential_j()) - work_j;
        let tolerance = self.budget.energy_absolute_tolerance_j
            + self.budget.energy_relative_tolerance
                * from
                    .potential_j()
                    .abs()
                    .max(to.potential_j().abs())
                    .max(work_j.abs());
        let estimated_absolute_error_j = error * work_scale;
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
