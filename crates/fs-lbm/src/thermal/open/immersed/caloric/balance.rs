//! Bounded scalar balance used by the original native caloric body. This
//! numerical component owns neither source data, physical state nor time.

#[derive(Clone, Debug, PartialEq)]
pub enum CaloricBalanceError<S> {
    Material(S),
    /// No root fits the admitted source enthalpy span; no endpoint clamping.
    OutsideEnthalpyRange {
        bound_j_kg: f64,
        residual: f64,
    },
    NonConvergence {
        iterations: usize,
        residual: f64,
    },
    Unrepresentable(&'static str),
}

pub(super) fn solve_balance<S, E>(
    old_h: f64,
    bounds: [f64; 2],
    mass_per_energy: f64,
    source: f64,
    conductance: f64,
    trace: f64,
    tolerance: f64,
    maximum_iterations: usize,
    evaluate: &mut impl FnMut(f64) -> Result<(f64, f64), CaloricBalanceError<S>>,
    poll: &mut impl FnMut() -> Result<(), E>,
    map: &impl Fn(CaloricBalanceError<S>) -> E,
) -> Result<(f64, f64), E> {
    let mut lo = bounds[0];
    let mut hi = bounds[1];
    let mut h = old_h;
    let mut residual = f64::INFINITY;
    for iteration in 0..maximum_iterations {
        poll()?;
        let (temperature, slope) = evaluate(h).map_err(map)?;
        let storage = mass_per_energy * (h - old_h);
        let heat = conductance * (temperature - trace);
        residual = (storage + heat) - source;
        let derivative = mass_per_energy + conductance * slope;
        if !temperature.is_finite()
            || !slope.is_finite()
            || slope < 0.0
            || !storage.is_finite()
            || !heat.is_finite()
            || !residual.is_finite()
            || !derivative.is_finite()
            || derivative <= 0.0
        {
            return Err(map(CaloricBalanceError::Unrepresentable(
                "caloric endpoint balance",
            )));
        }
        if residual.abs() <= tolerance {
            return Ok((h, temperature));
        }
        if residual < 0.0 {
            lo = h;
        } else {
            hi = h;
        }
        // Only the required source-domain edge is queried. The bulk fluid's
        // temperature need not be a valid material temperature for this body.
        if iteration == 0 {
            let bound = if residual < 0.0 { hi } else { lo };
            poll()?;
            let (t, _) = evaluate(bound).map_err(map)?;
            let r = mass_per_energy * (bound - old_h) + conductance * (t - trace) - source;
            if !r.is_finite() {
                return Err(map(CaloricBalanceError::Unrepresentable(
                    "caloric bracket balance",
                )));
            }
            if r.abs() <= tolerance {
                return Ok((bound, t));
            }
            if (residual < 0.0 && r < 0.0) || (residual > 0.0 && r > 0.0) {
                return Err(map(CaloricBalanceError::OutsideEnthalpyRange {
                    bound_j_kg: bound,
                    residual: r,
                }));
            }
        }
        let newton = h - residual / derivative;
        // Newton stays inside the admitted monotone bracket. A rejected
        // proposal uses a midpoint without forming an overflowing span.
        let middle = 0.5 * lo + 0.5 * hi;
        let next = if newton.is_finite() && newton > lo && newton < hi {
            newton
        } else {
            middle
        };
        if next == h || next <= lo || next >= hi {
            return Err(map(CaloricBalanceError::NonConvergence {
                iterations: iteration + 1,
                residual,
            }));
        }
        h = next;
    }
    Err(map(CaloricBalanceError::NonConvergence {
        iterations: maximum_iterations,
        residual,
    }))
}

#[cfg(test)]
mod tests;
