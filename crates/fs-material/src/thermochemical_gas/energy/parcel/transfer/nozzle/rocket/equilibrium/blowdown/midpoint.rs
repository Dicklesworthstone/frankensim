//! Bounded implicit midpoint for the original finite-reservoir outflow owners.
//!
//! Solve dm = dt * rate(m_old - dm/2). The caller supplies its actual
//! thermodynamic/nozzle evaluation; this scalar owner never reconstructs gas,
//! changes composition or accounts for energy. Publication remains with the
//! original reservoir after its complete endpoint admission.

/// Actual rate and optional nonnegative derivative with respect to tank mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OutflowRate {
    pub mass_flow_kg_s: f64,
    /// d(mdot)/d(m_remaining) [1/s]. Used only to propose a safeguarded Newton
    /// step; the actual mass-rate residual is the sole convergence criterion.
    /// `None` uses bisection, including an unrepresentable derivative near rest.
    pub derivative_per_s: Option<f64>,
}

/// Numerical and physical inventory bounds, shared by all outflow consumers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OutflowMidpointSolve {
    pub duration_s: f64,
    /// Fraction of the OLD mass, in (0,0.5]. A larger requested step refuses.
    pub maximum_outgoing_fraction: f64,
    pub mass_tolerance_kg: f64,
    /// Between one and 128 actual midpoint rate evaluations.
    pub maximum_iterations: usize,
}

/// Accepted scalar mass quadrature; this alone is not a physical gas endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OutflowMidpoint {
    pub remaining_mass_kg: f64,
    /// Computed from the actual representable old and remaining masses.
    pub expelled_mass_kg: f64,
    pub midpoint_mass_kg: f64,
    pub midpoint_rate: OutflowRate,
    pub mass_residual_kg: f64,
    pub iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum OutflowMidpointError<E> {
    Invalid(&'static str),
    Source(E),
    /// Largest duration at the authored outgoing-fraction boundary for this
    /// midpoint rule. This is a retry bound, not a temporal accuracy estimate.
    InventoryLimit { maximum_duration_s: f64 },
    NoConvergence { iterations: usize, mass_residual_kg: f64 },
    Cancelled,
}

/// Resolve an actual implicit midpoint mass equation without clipping a step.
///
/// The optional Newton tangent need not be exact: bracketed bisection remains
/// available and every accepted point closes its actual source rate. A zero
/// initial rate returns the identical inventory without a trial. Finite
/// positive flow that cannot change the represented mass refuses.
pub fn solve<E>(
    initial_mass_kg: f64,
    initial_rate_kg_s: f64,
    controls: OutflowMidpointSolve,
    mut rate_at_mass: impl FnMut(f64, &mut dyn FnMut() -> bool) -> Result<OutflowRate, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<OutflowMidpoint, OutflowMidpointError<E>> {
    use OutflowMidpointError as Error;
    if cancelled() {
        return Err(Error::Cancelled);
    }
    if !positive(initial_mass_kg)
        || !initial_rate_kg_s.is_finite()
        || initial_rate_kg_s < 0.0
        || !positive(controls.duration_s)
        || !positive(controls.maximum_outgoing_fraction)
        || controls.maximum_outgoing_fraction > 0.5
        || !positive(controls.mass_tolerance_kg)
        || !(1..=128).contains(&controls.maximum_iterations)
    {
        return Err(Error::Invalid("invalid finite outflow inventory, time or work bounds"));
    }
    if initial_rate_kg_s == 0.0 {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        return Ok(OutflowMidpoint {
            remaining_mass_kg: initial_mass_kg,
            expelled_mass_kg: 0.0,
            midpoint_mass_kg: initial_mass_kg,
            midpoint_rate: OutflowRate { mass_flow_kg_s: 0.0, derivative_per_s: Some(0.0) },
            mass_residual_kg: 0.0,
            iterations: 0,
        });
    }
    let limit = initial_mass_kg * controls.maximum_outgoing_fraction;
    let prediction = controls.duration_s * initial_rate_kg_s;
    if !positive(limit) || initial_mass_kg - limit == initial_mass_kg || !positive(prediction) {
        return Err(Error::Invalid("unrepresentable outgoing mass allowance or rate increment"));
    }
    let (mut low, mut high) = (0.0, limit);
    let mut dm = prediction.min(limit);
    let mut last_residual = f64::INFINITY;
    let mut iterations = 0;
    for iteration in 1..=controls.maximum_iterations {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        iterations = iteration;
        let remaining = initial_mass_kg - dm;
        let actual = initial_mass_kg - remaining;
        let midpoint = initial_mass_kg - 0.5 * actual;
        if !positive(remaining) || !positive(actual)
            || midpoint == initial_mass_kg || midpoint == remaining {
            return Err(Error::Invalid("distinct representable remaining and midpoint masses required"));
        }
        let rate = rate_at_mass(midpoint, &mut cancelled).map_err(Error::Source)?;
        if !rate.mass_flow_kg_s.is_finite() || rate.mass_flow_kg_s < 0.0
            || rate.derivative_per_s.is_some_and(|s| !s.is_finite() || s < 0.0) {
            return Err(Error::Invalid("finite nonnegative outflow rate and optional tangent required"));
        }
        let requested = controls.duration_s * rate.mass_flow_kg_s;
        if !requested.is_finite() || (rate.mass_flow_kg_s > 0.0 && requested == 0.0) {
            return Err(Error::Invalid("unrepresentable midpoint rate increment"));
        }
        last_residual = actual - requested;
        if last_residual.abs() <= controls.mass_tolerance_kg && actual <= limit {
            if cancelled() {
                return Err(Error::Cancelled);
            }
            return Ok(OutflowMidpoint {
                remaining_mass_kg: remaining,
                expelled_mass_kg: actual,
                midpoint_mass_kg: midpoint,
                midpoint_rate: rate,
                mass_residual_kg: last_residual,
                iterations,
            });
        }
        if last_residual < 0.0 {
            if dm == limit {
                let maximum_duration_s = actual / rate.mass_flow_kg_s;
                if !positive(maximum_duration_s) {
                    return Err(Error::Invalid("unrepresentable inventory retry duration"));
                }
                return Err(Error::InventoryLimit { maximum_duration_s });
            }
            low = dm;
        } else {
            high = dm;
        }
        let newton = rate.derivative_per_s.and_then(|derivative| {
            let slope = 1.0 + 0.5 * controls.duration_s * derivative;
            let next = actual - last_residual / slope;
            (slope.is_finite() && slope >= 1.0 && next.is_finite() && next > low && next < high)
                .then_some(next)
        });
        let next = newton.unwrap_or(low + 0.5 * (high - low));
        if next == dm || next == low || next == high {
            break;
        }
        dm = next;
    }
    Err(Error::NoConvergence { iterations, mass_residual_kg: last_residual })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
mod tests;
