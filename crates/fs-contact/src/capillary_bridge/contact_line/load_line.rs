//! Finite mechanical compliance coupled to the original dynamic wetting law.
//!
//! The scalar mechanical equation is h=h_free+S*F_mean. Every force probe
//! advances the SAME immutable contact-line state over the SAME physical dt;
//! nonlinear iteration never becomes additional physical relaxation. The
//! original source supplies surface work, line friction and its implicit force
//! derivative. Positive local compatibility slope and propagated numerical
//! error are admitted independently of the constitutive solve.
//!
//! The outer evaluation ceiling counts all original `advance` calls, including
//! rejected probes. Each call retains the source's own bounded geometry and
//! contact-line budgets, so total work is bounded by their product. A native
//! receiving owner reserves one additional call for its actual-motion audit.
//! This is local continuation, not a global uniqueness or trajectory proof.

use super::{ContactLineError, ContactLineState, ContactLineStep, SourceContactLineBridge};
use crate::capillary_bridge::BridgeError;
use crate::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_exec::CancelGate;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineLoadLineConfig {
    pub mechanical: BridgeLoadLineConfig,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ContactLineLoadLineError {
    Invalid(&'static str),
    Source(ContactLineError),
    Cancelled,
    Budget,
    GapIncrement,
    ForceAccuracy {
        estimated_error_n: f64,
        tolerance_n: f64,
    },
    CompatibilityAccuracy {
        estimated_error_m: f64,
        tolerance_m: f64,
    },
    CompatibilitySlope {
        slope: f64,
        minimum: f64,
    },
    NotConverged {
        residual_m: f64,
        tolerance_m: f64,
        iterations: usize,
    },
}
impl core::fmt::Display for ContactLineLoadLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "dynamic contact-line mechanical coupling refused: {self:?}"
        )
    }
}
impl std::error::Error for ContactLineLoadLineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Source(error) => Some(error),
            _ => None,
        }
    }
}
impl From<ContactLineError> for ContactLineLoadLineError {
    fn from(error: ContactLineError) -> Self {
        Self::Source(error)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineLoadLineStep {
    pub bridge: ContactLineStep,
    pub compatibility_residual_m: f64,
    pub compatibility_tolerance_m: f64,
    pub estimated_compatibility_error_m: f64,
    pub compatibility_slope: f64,
    pub iterations: usize,
    /// Constitutive transition calls; each has its original inner work ceiling.
    pub evaluations: usize,
}

/// Allocation-free continuation borrowing the original source and all history.
#[derive(Debug)]
pub struct ContactLineBridgeLoadLine<'a> {
    source: &'a SourceContactLineBridge,
    compliance_m_per_n: f64,
    config: ContactLineLoadLineConfig,
}
impl<'a> ContactLineBridgeLoadLine<'a> {
    pub fn new(
        source: &'a SourceContactLineBridge,
        compliance_m_per_n: f64,
        config: ContactLineLoadLineConfig,
    ) -> Result<Self, ContactLineLoadLineError> {
        let m = config.mechanical;
        if !compliance_m_per_n.is_finite()
            || compliance_m_per_n < 0.0
            || !(1..=128).contains(&m.max_iterations)
            || m.max_evaluations == 0
            || !(1..=32).contains(&m.max_backtracks)
            || !positive(m.gap_absolute_tolerance_m)
            || !positive(m.gap_relative_tolerance)
            || m.gap_relative_tolerance >= 1.0
            || !positive(m.maximum_gap_increment_m)
            || !positive(m.minimum_compatibility_slope)
            || m.minimum_compatibility_slope > 1.0
            || !positive(config.force_absolute_tolerance_n)
            || !config.force_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&config.force_relative_tolerance)
        {
            return Err(ContactLineLoadLineError::Invalid(
                "finite nonnegative compliance and positive mechanical/accuracy budgets required",
            ));
        }
        Ok(Self {
            source,
            compliance_m_per_n,
            config,
        })
    }
    pub const fn source(&self) -> &'a SourceContactLineBridge {
        self.source
    }
    pub const fn compliance_m_per_n(&self) -> f64 {
        self.compliance_m_per_n
    }
    pub const fn config(&self) -> ContactLineLoadLineConfig {
        self.config
    }

    /// An independent original constitutive proposal on final physical motion.
    pub fn advance_actual(
        &self,
        previous: ContactLineState,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<ContactLineStep, ContactLineLoadLineError> {
        poll(gate)?;
        self.check_increment(previous, gap_m)?;
        let step = self.source.advance(previous, gap_m, dt_s, gate)?;
        self.admit_force_error(&step)?;
        poll(gate)?;
        Ok(step)
    }

    pub fn trial(
        &self,
        previous: ContactLineState,
        free_gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<ContactLineLoadLineStep, ContactLineLoadLineError> {
        poll(gate)?;
        if !free_gap_m.is_finite() {
            return Err(ContactLineLoadLineError::Invalid(
                "finite free gap required",
            ));
        }
        let initial = if self.compliance_m_per_n == 0.0 {
            free_gap_m
        } else {
            previous.gap_m()
        };
        let mut evaluations = 0;
        let mut current =
            self.evaluate(previous, initial, free_gap_m, dt_s, &mut evaluations, gate)?;
        for iteration in 0..=self.config.mechanical.max_iterations {
            let residual = current.compatibility_residual_m.abs();
            if residual + current.estimated_compatibility_error_m
                <= current.compatibility_tolerance_m
            {
                poll(gate)?;
                current.iterations = iteration;
                current.evaluations = evaluations;
                return Ok(current);
            }
            if iteration == self.config.mechanical.max_iterations {
                return Err(ContactLineLoadLineError::NotConverged {
                    residual_m: residual,
                    tolerance_m: current.compatibility_tolerance_m,
                    iterations: iteration,
                });
            }
            let correction =
                checked(-current.compatibility_residual_m / current.compatibility_slope)?;
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..self.config.mechanical.max_backtracks {
                poll(gate)?;
                let candidate = current.bridge.state().gap_m() + fraction * correction;
                match self.evaluate(
                    previous,
                    candidate,
                    free_gap_m,
                    dt_s,
                    &mut evaluations,
                    gate,
                ) {
                    Ok(next)
                        if next.compatibility_residual_m.abs() < residual
                            || next.compatibility_residual_m.abs()
                                + next.estimated_compatibility_error_m
                                <= next.compatibility_tolerance_m =>
                    {
                        accepted = Some(next);
                        break;
                    }
                    Ok(_) => {}
                    Err(error) if outside_trial_domain(&error) => {}
                    Err(error) => return Err(error),
                }
                fraction *= 0.5;
            }
            current = accepted.ok_or(ContactLineLoadLineError::NotConverged {
                residual_m: residual,
                tolerance_m: current.compatibility_tolerance_m,
                iterations: iteration + 1,
            })?;
        }
        unreachable!("bounded contact-line continuation returns at its iteration ceiling")
    }

    fn check_increment(
        &self,
        previous: ContactLineState,
        gap_m: f64,
    ) -> Result<(), ContactLineLoadLineError> {
        let dh = gap_m - previous.gap_m();
        if !gap_m.is_finite()
            || !dh.is_finite()
            || dh.abs() > self.config.mechanical.maximum_gap_increment_m
        {
            return Err(ContactLineLoadLineError::GapIncrement);
        }
        Ok(())
    }
    fn admit_force_error(&self, step: &ContactLineStep) -> Result<f64, ContactLineLoadLineError> {
        let error = checked(step.estimated_force_error_n)?;
        if error < 0.0 {
            return Err(ContactLineLoadLineError::Invalid(
                "negative force error estimate",
            ));
        }
        let tolerance = checked(
            self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance * step.mean_normal_force_n.abs(),
        )?;
        if error > tolerance {
            return Err(ContactLineLoadLineError::ForceAccuracy {
                estimated_error_n: error,
                tolerance_n: tolerance,
            });
        }
        Ok(error)
    }
    fn evaluate(
        &self,
        previous: ContactLineState,
        gap_m: f64,
        free_gap_m: f64,
        dt_s: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<ContactLineLoadLineStep, ContactLineLoadLineError> {
        poll(gate)?;
        if *evaluations >= self.config.mechanical.max_evaluations {
            return Err(ContactLineLoadLineError::Budget);
        }
        *evaluations += 1;
        let bridge = self.advance_actual(previous, gap_m, dt_s, gate)?;
        let response = checked(self.compliance_m_per_n * bridge.mean_normal_force_n)?;
        let slope =
            checked(1.0 - self.compliance_m_per_n * bridge.mean_force_gap_derivative_n_per_m)?;
        if slope < self.config.mechanical.minimum_compatibility_slope {
            return Err(ContactLineLoadLineError::CompatibilitySlope {
                slope,
                minimum: self.config.mechanical.minimum_compatibility_slope,
            });
        }
        let residual = checked((bridge.state().gap_m() - free_gap_m) - response)?;
        let scale = (gap_m - previous.gap_m())
            .abs()
            .max((free_gap_m - previous.gap_m()).abs())
            .max(response.abs());
        let tolerance = checked(
            self.config.mechanical.gap_absolute_tolerance_m
                + self.config.mechanical.gap_relative_tolerance * scale,
        )?;
        let error = checked(self.compliance_m_per_n * self.admit_force_error(&bridge)?)?;
        if error > tolerance {
            return Err(ContactLineLoadLineError::CompatibilityAccuracy {
                estimated_error_m: error,
                tolerance_m: tolerance,
            });
        }
        Ok(ContactLineLoadLineStep {
            bridge,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
            estimated_compatibility_error_m: error,
            compatibility_slope: slope,
            iterations: 0,
            evaluations: 0,
        })
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn checked(x: f64) -> Result<f64, ContactLineLoadLineError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(ContactLineLoadLineError::Invalid(
            "nonfinite mechanical coupling value",
        ))
    }
}
fn poll(gate: &CancelGate) -> Result<(), ContactLineLoadLineError> {
    if gate.is_requested() {
        Err(ContactLineLoadLineError::Cancelled)
    } else {
        Ok(())
    }
}
fn outside_trial_domain(error: &ContactLineLoadLineError) -> bool {
    matches!(
        error,
        ContactLineLoadLineError::GapIncrement
            | ContactLineLoadLineError::CompatibilitySlope { .. }
            | ContactLineLoadLineError::Source(ContactLineError::Domain(_))
            | ContactLineLoadLineError::Source(ContactLineError::Original(
                BridgeError::Invalid { .. } | BridgeError::Limit { .. }
            ))
    )
}
