//! First-order instrument response driven by the original compiled sensor row.
//!
//! The retained reading obeys `tau * reading' = equilibrium - reading`, where
//! `equilibrium = H * field + offset` includes the declared mounting model.
//! Held and affine reconstructions are explicit: these are exact integrations
//! of those inputs, not claims about an unresolved physical trajectory between
//! accepted solver states. No measurement noise is drawn and no measurement is
//! assimilated; an uncertain initial instrument memory remains the caller's
//! inference problem.

use core::fmt;

use super::{CompiledSensorOperator, SensorError};

/// One accepted instrument memory, retaining its original compiled identity.
///
/// A step returns a new private candidate. Discarding it leaves the accepted
/// reading, target and clock untouched, so solver rejection and cancellation
/// cannot partially advance instrumentation.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorState<'a> {
    operator: &'a CompiledSensorOperator,
    time_s: f64,
    reading: f64,
    equilibrium_reading: f64,
}

/// Refusal before publication of an instrument state.
#[derive(Debug, Clone, PartialEq)]
pub enum SensorTransientError {
    /// The original field/operator evaluation refused.
    Sensor(SensorError),
    /// The clock, initial condition or response arithmetic is invalid.
    Invalid(&'static str),
    /// Cooperative cancellation observed at a bounded checkpoint.
    Cancelled,
}

impl fmt::Display for SensorTransientError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Sensor(error) => error.fmt(f),
            Self::Invalid(reason) => f.write_str(reason),
            Self::Cancelled => f.write_str("sensor evolution cancelled before publication"),
        }
    }
}

impl core::error::Error for SensorTransientError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Sensor(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SensorError> for SensorTransientError {
    fn from(error: SensorError) -> Self {
        Self::Sensor(error)
    }
}

impl<'a> SensorState<'a> {
    /// Bind an explicit initial reading to the same field and simulation clock.
    ///
    /// First-order instruments may start out of equilibrium. An instantaneous
    /// instrument must start at its compiled steady reading; an arbitrary
    /// retained memory would contradict that declared model.
    pub fn new(
        operator: &'a CompiledSensorOperator,
        time_s: f64,
        initial_field: &[f64],
        initial_reading: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SensorTransientError> {
        poll(cancelled)?;
        if !time_s.is_finite() || !initial_reading.is_finite() {
            return Err(SensorTransientError::Invalid(
                "initial sensor time and reading must be finite",
            ));
        }
        let equilibrium_reading =
            operator.predict_with_checkpoint(initial_field, &mut || poll(cancelled))?;
        if operator.dynamics().is_declared_instantaneous() && initial_reading != equilibrium_reading
        {
            return Err(SensorTransientError::Invalid(
                "instantaneous sensor must start at its equilibrium reading",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            operator,
            time_s,
            reading: initial_reading,
            equilibrium_reading,
        })
    }

    /// Original compiled support, mounting model, dynamics and identity.
    #[must_use]
    pub const fn operator(&self) -> &'a CompiledSensorOperator {
        self.operator
    }

    /// Exact accepted simulation time in seconds.
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Actual instrument reading, including its retained response lag.
    #[must_use]
    pub const fn reading(&self) -> f64 {
        self.reading
    }

    /// Compiled steady reading of the last accepted field.
    #[must_use]
    pub const fn equilibrium_reading(&self) -> f64 {
        self.equilibrium_reading
    }

    /// Evolve with the supplied endpoint field held over the whole interval.
    ///
    /// A jump from the retained equilibrium target is explicit in this method.
    /// `endpoint_time_s` must be finite, strictly later than the accepted clock,
    /// and separated by a finite duration. The exact supplied clock is retained.
    pub fn step_held_to_candidate(
        &self,
        endpoint_time_s: f64,
        endpoint_field: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SensorTransientError> {
        self.step(
            endpoint_time_s,
            endpoint_time_s - self.time_s,
            endpoint_field,
            false,
            cancelled,
        )
    }

    /// Evolve with an affine target between the accepted and endpoint fields.
    ///
    /// The starting target comes from this state's accepted field, so it cannot
    /// be replaced accidentally by a stale or unrelated starting snapshot. This
    /// is exact for an affine `H * field + offset`, not necessarily for the
    /// trajectory of a nonlinear thermal or fluid solver between endpoints.
    pub fn step_affine_to_candidate(
        &self,
        endpoint_time_s: f64,
        endpoint_field: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SensorTransientError> {
        self.step_affine_interval_candidate(
            endpoint_time_s,
            endpoint_time_s - self.time_s,
            endpoint_field,
            cancelled,
        )
    }

    /// Evolve over the physical owner's authoritative affine-input duration.
    ///
    /// Some solvers integrate a chosen `duration_s` and round the sum when
    /// publishing the endpoint clock. Others clip a step to an exact event
    /// clock and subtract to obtain its duration. Admit either convention
    /// exactly, preserving both the native evolution duration and endpoint.
    pub fn step_affine_interval_candidate(
        &self,
        endpoint_time_s: f64,
        duration_s: f64,
        endpoint_field: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SensorTransientError> {
        self.step(endpoint_time_s, duration_s, endpoint_field, true, cancelled)
    }

    fn step(
        &self,
        endpoint_time_s: f64,
        duration_s: f64,
        endpoint_field: &[f64],
        affine: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, SensorTransientError> {
        poll(cancelled)?;
        if !endpoint_time_s.is_finite()
            || endpoint_time_s <= self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || (self.time_s + duration_s != endpoint_time_s
                && endpoint_time_s - self.time_s != duration_s)
        {
            return Err(SensorTransientError::Invalid(
                "sensor endpoint and finite positive duration must agree with the accepted clock",
            ));
        }
        let target = self
            .operator
            .predict_with_checkpoint(endpoint_field, &mut || poll(cancelled))?;
        let reading = match self.operator.dynamics().time_constant_s() {
            None => target,
            Some(tau) => {
                let x = duration_s / tau;
                if x == 0.0 {
                    return Err(SensorTransientError::Invalid(
                        "sensor duration divided by time constant is not representable",
                    ));
                }
                let response = -fs_math::det::expm1(-x);
                // Keep exp independently: 1-response can already round to
                // zero while a large retained reading still has a finite tail.
                let memory_weight = fs_math::det::exp(-x);
                let start = if affine {
                    self.equilibrium_reading
                } else {
                    target
                };
                // b = 1 - (1-exp(-x))/x; its small-x series avoids cancelling
                // almost equal numbers for fast accepted solver substeps.
                let (start_weight, endpoint_weight) = if !affine {
                    (0.0, response)
                } else if x < 1.0e-3 {
                    let endpoint = x
                        * (0.5
                            + x * (-1.0 / 6.0 + x * (1.0 / 24.0 + x * (-1.0 / 120.0 + x / 720.0))));
                    (response - endpoint, endpoint)
                } else {
                    // This expression also retains the O(1/x) starting-target
                    // contribution when the endpoint weight rounds to one.
                    (response / x - memory_weight, 1.0 - response / x)
                };
                let reading = if self.reading == start && start == target {
                    target
                } else {
                    memory_weight.mul_add(
                        self.reading,
                        start_weight.mul_add(start, endpoint_weight * target),
                    )
                };
                if !reading.is_finite() {
                    return Err(SensorTransientError::Invalid(
                        "sensor response arithmetic exceeded its finite range",
                    ));
                }
                reading
            }
        };
        poll(cancelled)?;
        Ok(Self {
            operator: self.operator,
            time_s: endpoint_time_s,
            reading,
            equilibrium_reading: target,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SensorTransientError> {
    if cancelled() {
        Err(SensorTransientError::Cancelled)
    } else {
        Ok(())
    }
}
