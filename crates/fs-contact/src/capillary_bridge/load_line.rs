//! Work-conjugate capillary adhesion on an original held-force mechanical map.
//!
//! The mechanical owner supplies its free endpoint gap and full projected
//! compliance S. The coupled equation is h=h_free+S*F_mean(h), with the SAME
//! analytic path force used for body work, including rupture and coalescence.
//! Original liquid/source/topology admission is never bypassed by root trials.
//! This bounded local continuation is not a global uniqueness or temporal
//! accuracy certificate. A vanishing/negative compatibility slope refuses the
//! sample for physical timestep refinement; no false detached or clipped state
//! is substituted for an unresolved attractive instability.

use super::{BridgeError, BridgeState, BridgeStep, SourcePlanarCapillaryBridge};
use fs_exec::CancelGate;

/// Numerical work and physical increment admission for one coupled sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BridgeLoadLineConfig {
    pub max_iterations: usize,
    /// Every original bridge proposal, including rejected line-search probes.
    pub max_evaluations: usize,
    pub max_backtracks: usize,
    pub gap_absolute_tolerance_m: f64,
    pub gap_relative_tolerance: f64,
    pub maximum_gap_increment_m: f64,
    /// Lower bound on 1-S*dF_mean/dh, in (0,1].
    pub minimum_compatibility_slope: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BridgeLoadLineError {
    Invalid(&'static str),
    Bridge(BridgeError),
    Cancelled,
    Budget,
    GapIncrement,
    /// The selected discrete attractive branch needs timestep refinement.
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
impl core::fmt::Display for BridgeLoadLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "capillary load line refused: {self:?}")
    }
}
impl std::error::Error for BridgeLoadLineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Bridge(error) => Some(error),
            _ => None,
        }
    }
}
impl From<BridgeError> for BridgeLoadLineError {
    fn from(error: BridgeError) -> Self {
        Self::Bridge(error)
    }
}

/// Immutable bridge successor and the independently measurable compatibility.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BridgeLoadLineStep {
    pub bridge: BridgeStep,
    pub compatibility_residual_m: f64,
    pub compatibility_tolerance_m: f64,
    pub compatibility_slope: f64,
    pub iterations: usize,
    pub evaluations: usize,
}

/// Allocation-free scalar solve borrowing the complete original source owner.
#[derive(Debug)]
pub struct PlanarBridgeLoadLine<'a> {
    source: &'a SourcePlanarCapillaryBridge,
    compliance_m_per_n: f64,
    config: BridgeLoadLineConfig,
}
impl<'a> PlanarBridgeLoadLine<'a> {
    pub fn new(
        source: &'a SourcePlanarCapillaryBridge,
        compliance_m_per_n: f64,
        config: BridgeLoadLineConfig,
    ) -> Result<Self, BridgeLoadLineError> {
        if !compliance_m_per_n.is_finite()
            || compliance_m_per_n < 0.0
            || !(1..=128).contains(&config.max_iterations)
            || config.max_evaluations == 0
            || !(1..=32).contains(&config.max_backtracks)
            || !positive(config.gap_absolute_tolerance_m)
            || !positive(config.gap_relative_tolerance)
            || config.gap_relative_tolerance >= 1.0
            || !positive(config.maximum_gap_increment_m)
            || !positive(config.minimum_compatibility_slope)
            || config.minimum_compatibility_slope > 1.0
        {
            return Err(BridgeLoadLineError::Invalid(
                "finite nonnegative compliance and explicit positive work/compatibility budgets required",
            ));
        }
        Ok(Self {
            source,
            compliance_m_per_n,
            config,
        })
    }
    #[must_use]
    pub const fn source(&self) -> &'a SourcePlanarCapillaryBridge {
        self.source
    }
    #[must_use]
    pub const fn compliance_m_per_n(&self) -> f64 {
        self.compliance_m_per_n
    }
    #[must_use]
    pub const fn config(&self) -> BridgeLoadLineConfig {
        self.config
    }

    /// The original accepted state is immutable. Only a complete native
    /// mechanics and energy audit authorizes publication of the returned state.
    /// The free gap may be outside the source range; all evaluated candidates
    /// must satisfy the original source and continuum/rate/energy limits.
    pub fn trial(
        &self,
        previous: BridgeState,
        free_gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<BridgeLoadLineStep, BridgeLoadLineError> {
        poll(gate)?;
        if !free_gap_m.is_finite() || !positive(dt_s) {
            return Err(BridgeLoadLineError::Invalid(
                "finite free gap and positive timestep required",
            ));
        }
        // A zero-compliance fixture is an explicit prescribed trajectory.
        // Otherwise start from the accepted physical state, which avoids
        // assuming an unconstrained predictor satisfies the liquid's domain.
        let initial = if self.compliance_m_per_n == 0.0 {
            free_gap_m
        } else {
            previous.gap_m()
        };
        let mut evaluations = 0;
        let mut current =
            self.evaluate(previous, initial, free_gap_m, dt_s, &mut evaluations, gate)?;
        for iteration in 0..=self.config.max_iterations {
            let residual = current.compatibility_residual_m.abs();
            if residual <= current.compatibility_tolerance_m {
                poll(gate)?;
                current.iterations = iteration;
                current.evaluations = evaluations;
                return Ok(current);
            }
            if iteration == self.config.max_iterations {
                return Err(BridgeLoadLineError::NotConverged {
                    residual_m: residual,
                    tolerance_m: current.compatibility_tolerance_m,
                    iterations: iteration,
                });
            }
            let correction =
                checked(-current.compatibility_residual_m / current.compatibility_slope)?;
            let old_gap = current.bridge.state().gap_m();
            let mut scale = 1.0;
            let mut accepted = None;
            for _ in 0..self.config.max_backtracks {
                poll(gate)?;
                let candidate = old_gap + scale * correction;
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
            current = accepted.ok_or(BridgeLoadLineError::NotConverged {
                residual_m: residual,
                tolerance_m: current.compatibility_tolerance_m,
                iterations: iteration + 1,
            })?;
        }
        unreachable!("bounded capillary continuation returns at the iteration ceiling")
    }

    fn evaluate(
        &self,
        previous: BridgeState,
        gap: f64,
        free: f64,
        dt: f64,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<BridgeLoadLineStep, BridgeLoadLineError> {
        poll(gate)?;
        let increment = gap - previous.gap_m();
        if !gap.is_finite()
            || !increment.is_finite()
            || increment.abs() > self.config.maximum_gap_increment_m
        {
            return Err(BridgeLoadLineError::GapIncrement);
        }
        if *evaluations >= self.config.max_evaluations {
            return Err(BridgeLoadLineError::Budget);
        }
        *evaluations += 1;
        let bridge = self.source.advance(previous, gap, dt, gate)?;
        let response = checked(self.compliance_m_per_n * bridge.mean_normal_force_n)?;
        let slope =
            checked(1.0 - self.compliance_m_per_n * bridge.mean_force_gap_derivative_n_per_m)?;
        if slope < self.config.minimum_compatibility_slope {
            return Err(BridgeLoadLineError::CompatibilitySlope {
                slope,
                minimum: self.config.minimum_compatibility_slope,
            });
        }
        let residual = checked((gap - free) - response)?;
        let scale = increment
            .abs()
            .max((free - previous.gap_m()).abs())
            .max(response.abs());
        let tolerance = checked(
            self.config.gap_absolute_tolerance_m + self.config.gap_relative_tolerance * scale,
        )?;
        Ok(BridgeLoadLineStep {
            bridge,
            compatibility_residual_m: residual,
            compatibility_tolerance_m: tolerance,
            compatibility_slope: slope,
            iterations: 0,
            evaluations: 0,
        })
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn checked(value: f64) -> Result<f64, BridgeLoadLineError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BridgeLoadLineError::Invalid(
            "nonfinite derived load-line value",
        ))
    }
}
fn poll(gate: &CancelGate) -> Result<(), BridgeLoadLineError> {
    if gate.is_requested() {
        Err(BridgeLoadLineError::Cancelled)
    } else {
        Ok(())
    }
}
fn outside_trial_domain(error: &BridgeLoadLineError) -> bool {
    matches!(
        error,
        BridgeLoadLineError::GapIncrement
            | BridgeLoadLineError::CompatibilitySlope { .. }
            | BridgeLoadLineError::Bridge(
                BridgeError::Invalid { .. }
                    | BridgeError::Limit { .. }
                    | BridgeError::NoConnectedProfile
                    | BridgeError::OverlappingCaps
            )
    )
}
