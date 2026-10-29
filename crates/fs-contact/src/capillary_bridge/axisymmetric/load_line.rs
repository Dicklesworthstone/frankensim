//! The actual axisymmetric liquid bridge against a mechanical load line.
//!
//! The mechanical owner supplies its free gap and complete projected held-force
//! compliance S. This owner solves h=h_free+S*F_mean with the original source's
//! energy-secant force and its actual converged profile gap. Every proposal
//! retains the original source, branch, stability, quadrature and work gates.
//! There is no liquid time integrator: the receiving mechanics must establish
//! that the original quasistatic bridge assumptions apply at its physical rate.
//!
//! A finite-step force uncertainty is estimated from the original energy error
//! divided by the actual gap increment. At exact rest, the original endpoint
//! force and half endpoint tangent are used with their quadrature estimate.
//! A distinct but unresolved physical increment refuses; it is never replaced
//! by rest, a force projection, or an invented rupture. The error estimates and
//! positive local compatibility slope are not global uniqueness, interval, or
//! trajectory-accuracy certificates. Mechanical timestep refinement can address
//! an attractive instability; tighter source accuracy may be needed for small
//! energy differences, including the accepted initial state.

use super::{
    AxisymmetricError, AxisymmetricState, AxisymmetricStep, AxisymmetricWorkBudget,
    SourceAxisymmetricCapillaryBridge,
};
use crate::capillary_bridge::BridgeError;
use crate::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_exec::CancelGate;

/// Separate outer mechanical and original equilibrium/work accuracy budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricLoadLineConfig {
    pub mechanical: BridgeLoadLineConfig,
    pub work: AxisymmetricWorkBudget,
    /// Admission for the estimated error of the work-conjugate mean force.
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AxisymmetricLoadLineError {
    Invalid(&'static str),
    Bridge(AxisymmetricError),
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
impl core::fmt::Display for AxisymmetricLoadLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "axisymmetric capillary load line refused: {self:?}")
    }
}
impl std::error::Error for AxisymmetricLoadLineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Bridge(error) => Some(error),
            _ => None,
        }
    }
}
impl From<AxisymmetricError> for AxisymmetricLoadLineError {
    fn from(error: AxisymmetricError) -> Self {
        Self::Bridge(error)
    }
}

/// Immutable physical successor and independently inspectable acceptance data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricLoadLineStep {
    pub bridge: AxisymmetricStep,
    /// Uses the actual converged gap, not the requested equilibrium coordinate.
    pub compatibility_residual_m: f64,
    pub compatibility_tolerance_m: f64,
    /// Numerical force and profile-gap estimates propagated into compatibility.
    pub estimated_compatibility_error_m: f64,
    pub compatibility_slope: f64,
    pub estimated_force_error_n: f64,
    pub force_error_tolerance_n: f64,
    pub iterations: usize,
    /// All original transition proposals, including failed backtracking probes.
    /// Each has its own explicit equilibrium root/quadrature ceilings.
    pub evaluations: usize,
}

/// Allocation-free local continuation borrowing the complete original source.
#[derive(Debug)]
pub struct AxisymmetricBridgeLoadLine<'a> {
    source: &'a SourceAxisymmetricCapillaryBridge,
    compliance_m_per_n: f64,
    config: AxisymmetricLoadLineConfig,
}
impl<'a> AxisymmetricBridgeLoadLine<'a> {
    pub fn new(
        source: &'a SourceAxisymmetricCapillaryBridge,
        compliance_m_per_n: f64,
        config: AxisymmetricLoadLineConfig,
    ) -> Result<Self, AxisymmetricLoadLineError> {
        let mechanical = config.mechanical;
        config.work.equilibrium.validate()?;
        if !compliance_m_per_n.is_finite()
            || compliance_m_per_n < 0.0
            || !(1..=128).contains(&mechanical.max_iterations)
            || mechanical.max_evaluations == 0
            || !(1..=32).contains(&mechanical.max_backtracks)
            || !positive(mechanical.gap_absolute_tolerance_m)
            || !positive(mechanical.gap_relative_tolerance)
            || mechanical.gap_relative_tolerance >= 1.0
            || !positive(mechanical.maximum_gap_increment_m)
            || !positive(mechanical.minimum_compatibility_slope)
            || mechanical.minimum_compatibility_slope > 1.0
            || !positive(config.work.absolute_tolerance_j)
            || !fraction(config.work.relative_tolerance)
            || !positive(config.force_absolute_tolerance_n)
            || !fraction(config.force_relative_tolerance)
        {
            return Err(AxisymmetricLoadLineError::Invalid(
                "nonnegative finite compliance and explicit positive mechanical/source/force budgets required",
            ));
        }
        Ok(Self {
            source,
            compliance_m_per_n,
            config,
        })
    }
    #[must_use]
    pub const fn source(&self) -> &'a SourceAxisymmetricCapillaryBridge {
        self.source
    }
    #[must_use]
    pub const fn compliance_m_per_n(&self) -> f64 {
        self.compliance_m_per_n
    }
    #[must_use]
    pub const fn config(&self) -> AxisymmetricLoadLineConfig {
        self.config
    }

    /// One original constitutive transition on actual staged mechanical motion.
    /// A native owner reserves one evaluation for this final independent check.
    /// This applies the same source, physical-increment and force-error gates as
    /// root proposals, without performing a second equilibrium solve.
    pub fn advance_actual(
        &self,
        previous: AxisymmetricState,
        requested_gap_m: f64,
        gate: &CancelGate,
    ) -> Result<AxisymmetricStep, AxisymmetricLoadLineError> {
        poll(gate)?;
        self.check_increment(previous, requested_gap_m)?;
        let bridge =
            self.source
                .continue_to_gap(previous, requested_gap_m, self.config.work, gate)?;
        self.check_increment(previous, bridge.state().profile().gap_m)?;
        self.admit_force_error(previous, &bridge)?;
        poll(gate)?;
        Ok(bridge)
    }

    /// Pure observer; does not solve a shape or mutate either source state.
    /// The value is an estimate, not an interval-certified force enclosure.
    pub fn estimated_force_error_n(
        &self,
        previous: AxisymmetricState,
        bridge: &AxisymmetricStep,
    ) -> Result<f64, AxisymmetricLoadLineError> {
        self.source.require_state(previous)?;
        let next = bridge.state();
        self.source.require_state(next)?;
        let dh = next.profile().gap_m - previous.profile().gap_m;
        if dh == 0.0 {
            // F=-pi*gamma*a*(2+q), a proportional to the inverse cube root of
            // the positive volume quadrature. q is the represented state.
            checked(
                bridge.mean_normal_force_n.abs()
                    * (next.report().quadrature_relative_error_estimate / 3.0
                        + 32.0 * f64::EPSILON),
            )
        } else {
            checked(
                bridge.estimated_work_error_j / dh.abs()
                    + 16.0 * f64::EPSILON * bridge.mean_normal_force_n.abs(),
            )
        }
    }

    /// Start from the accepted physical gap, retaining the connected branch.
    /// Zero compliance is an explicit prescribed-motion fixture. A proposal is
    /// publishable only after the receiving mechanics audits its actual work.
    pub fn trial(
        &self,
        previous: AxisymmetricState,
        free_gap_m: f64,
        gate: &CancelGate,
    ) -> Result<AxisymmetricLoadLineStep, AxisymmetricLoadLineError> {
        poll(gate)?;
        if !free_gap_m.is_finite() {
            return Err(AxisymmetricLoadLineError::Invalid(
                "finite free gap required",
            ));
        }
        let initial = if self.compliance_m_per_n == 0.0 {
            free_gap_m
        } else {
            previous.profile().gap_m
        };
        let mut evaluations = 0;
        let mut current = self.evaluate(previous, initial, free_gap_m, &mut evaluations, gate)?;
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
                return Err(AxisymmetricLoadLineError::NotConverged {
                    residual_m: residual,
                    tolerance_m: current.compatibility_tolerance_m,
                    iterations: iteration,
                });
            }
            let correction =
                checked(-current.compatibility_residual_m / current.compatibility_slope)?;
            let actual_gap = current.bridge.state().profile().gap_m;
            let mut scale = 1.0;
            let mut accepted = None;
            for _ in 0..self.config.mechanical.max_backtracks {
                poll(gate)?;
                let candidate = actual_gap + scale * correction;
                match self.evaluate(previous, candidate, free_gap_m, &mut evaluations, gate) {
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
                scale *= 0.5;
            }
            current = accepted.ok_or(AxisymmetricLoadLineError::NotConverged {
                residual_m: residual,
                tolerance_m: current.compatibility_tolerance_m,
                iterations: iteration + 1,
            })?;
        }
        unreachable!("bounded axisymmetric continuation returns at its iteration ceiling")
    }

    fn check_increment(
        &self,
        previous: AxisymmetricState,
        gap: f64,
    ) -> Result<(), AxisymmetricLoadLineError> {
        let increment = gap - previous.profile().gap_m;
        if !gap.is_finite()
            || !increment.is_finite()
            || increment.abs() > self.config.mechanical.maximum_gap_increment_m
        {
            return Err(AxisymmetricLoadLineError::GapIncrement);
        }
        Ok(())
    }
    fn admit_force_error(
        &self,
        previous: AxisymmetricState,
        bridge: &AxisymmetricStep,
    ) -> Result<(f64, f64), AxisymmetricLoadLineError> {
        let error = self.estimated_force_error_n(previous, bridge)?;
        let tolerance = checked(
            self.config.force_absolute_tolerance_n
                + self.config.force_relative_tolerance * bridge.mean_normal_force_n.abs(),
        )?;
        if error > tolerance {
            return Err(AxisymmetricLoadLineError::ForceAccuracy {
                estimated_error_n: error,
                tolerance_n: tolerance,
            });
        }
        Ok((error, tolerance))
    }
    fn evaluate(
        &self,
        previous: AxisymmetricState,
        requested_gap: f64,
        free: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<AxisymmetricLoadLineStep, AxisymmetricLoadLineError> {
        poll(gate)?;
        if *evaluations >= self.config.mechanical.max_evaluations {
            return Err(AxisymmetricLoadLineError::Budget);
        }
        *evaluations += 1;
        let bridge = self.advance_actual(previous, requested_gap, gate)?;
        let state = bridge.state();
        let gap = state.profile().gap_m;
        let increment = gap - previous.profile().gap_m;
        let response = checked(self.compliance_m_per_n * bridge.mean_normal_force_n)?;
        let slope =
            checked(1.0 - self.compliance_m_per_n * bridge.mean_force_gap_derivative_n_per_m)?;
        if slope < self.config.mechanical.minimum_compatibility_slope {
            return Err(AxisymmetricLoadLineError::CompatibilitySlope {
                slope,
                minimum: self.config.mechanical.minimum_compatibility_slope,
            });
        }
        let residual = checked((gap - free) - response)?;
        let scale = increment
            .abs()
            .max((free - previous.profile().gap_m).abs())
            .max(response.abs());
        let tolerance = checked(
            self.config.mechanical.gap_absolute_tolerance_m
                + self.config.mechanical.gap_relative_tolerance * scale,
        )?;
        let (force_error, force_tolerance) = self.admit_force_error(previous, &bridge)?;
        // h=a*H, with relative estimates from positive height/volume integrals.
        let profile_gap_error = checked(
            gap * (4.0 / 3.0 * state.report().quadrature_relative_error_estimate
                + 32.0 * f64::EPSILON),
        )?;
        let compatibility_error =
            checked(profile_gap_error + self.compliance_m_per_n * force_error)?;
        if compatibility_error > tolerance {
            return Err(AxisymmetricLoadLineError::CompatibilityAccuracy {
                estimated_error_m: compatibility_error,
                tolerance_m: tolerance,
            });
        }
        Ok(AxisymmetricLoadLineStep {
            bridge,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
            estimated_compatibility_error_m: compatibility_error,
            compatibility_slope: slope,
            estimated_force_error_n: force_error,
            force_error_tolerance_n: force_tolerance,
            iterations: 0,
            evaluations: 0,
        })
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn fraction(value: f64) -> bool {
    value.is_finite() && (0.0..1.0).contains(&value)
}
fn checked(value: f64) -> Result<f64, AxisymmetricLoadLineError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(AxisymmetricLoadLineError::Invalid(
            "nonfinite derived load-line value",
        ))
    }
}
fn poll(gate: &CancelGate) -> Result<(), AxisymmetricLoadLineError> {
    if gate.is_requested() {
        Err(AxisymmetricLoadLineError::Cancelled)
    } else {
        Ok(())
    }
}
fn outside_trial_domain(error: &AxisymmetricLoadLineError) -> bool {
    matches!(
        error,
        AxisymmetricLoadLineError::GapIncrement
            | AxisymmetricLoadLineError::CompatibilitySlope { .. }
            | AxisymmetricLoadLineError::Bridge(
                AxisymmetricError::GapOutsideBranch { .. }
                    | AxisymmetricError::StabilityUnestablished { .. }
                    | AxisymmetricError::TurningBranch
                    | AxisymmetricError::Physical(
                        BridgeError::Invalid { .. } | BridgeError::Limit { .. }
                    )
            )
    )
}
