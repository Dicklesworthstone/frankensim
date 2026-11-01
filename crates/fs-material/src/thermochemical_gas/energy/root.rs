//! Shared bounded inversion for energy and entropy on an increasing gas branch.
//! The constitutive owner supplies the residual, its exact thermal tangent and
//! EOS pressure. This module does not supply a heat-capacity or gas law.

#[derive(Clone, Copy)]
pub(in crate::thermochemical_gas) struct Sample {
    pub residual: f64,
    pub tangent: f64,
    pub pressure_pa: f64,
}
pub(in crate::thermochemical_gas) struct Solution {
    pub temperature_k: f64,
    pub sample: Sample,
    pub iterations: usize,
}
pub(in crate::thermochemical_gas) enum Error<E> {
    Evaluation(E),
    Invalid(&'static str),
    Cancelled,
    OutsideBracket { lower: f64, upper: f64 },
    NotConverged { iterations: usize, residual: f64 },
}

// Callers validate their units, gas-temperature envelope and iteration cap.
// A small bracket is never a substitute for satisfying the physical residual.
pub(in crate::thermochemical_gas) fn solve<E>(
    lo: f64, hi: f64, tolerance: f64, maximum_iterations: usize,
    evaluate: impl FnMut(f64) -> Result<Sample, E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Solution, Error<E>> {
    solve_with_floor(lo, hi, tolerance, maximum_iterations, evaluate, cancelled, |_, _| 0.0)
}

/// [`solve`] with a representation floor. A temperature moves in quanta
/// (cv times the spacing of floats at T), so a residual below half a quantum
/// cannot be improved by ANY representable temperature. Where `solve` gives up
/// because the bracket cannot shrink further (the stagnation exit), this takes
/// the better of the two final bracket ends and accepts it when its residual is
/// within `floor(temperature, sample)`. Everything `solve` accepts is returned
/// bit for bit (the floor is consulted nowhere else), every refusal `solve`
/// raises with a larger final residual is raised unchanged, and a zero floor
/// is `solve`.
pub(in crate::thermochemical_gas) fn solve_with_floor<E>(
    mut lo: f64, mut hi: f64, tolerance: f64, maximum_iterations: usize,
    mut evaluate: impl FnMut(f64) -> Result<Sample, E>,
    cancelled: &mut impl FnMut() -> bool,
    floor: impl Fn(f64, &Sample) -> f64,
) -> Result<Solution, Error<E>> {
    let mut sample = |t| {
        if cancelled() { return Err(Error::Cancelled); }
        let value = evaluate(t).map_err(Error::Evaluation)?;
        if !value.residual.is_finite() || !value.tangent.is_finite() || value.tangent <= 0.0
            || !value.pressure_pa.is_finite() || value.pressure_pa <= 0.0
        {
            return Err(Error::Invalid("finite residual, positive tangent and representable EOS pressure required"));
        }
        Ok(value)
    };
    let lower = sample(lo)?;
    let upper = sample(hi)?;
    let (flo, fhi) = (lower.residual, upper.residual);
    if flo > fhi { return Err(Error::OutsideBracket { lower: flo, upper: fhi }); }
    for (temperature_k, value) in [(lo, lower), (hi, upper)] {
        if value.residual.abs() <= tolerance {
            return Ok(Solution { temperature_k, sample: value, iterations: 0 });
        }
    }
    if flo > 0.0 || fhi < 0.0 {
        return Err(Error::OutsideBracket { lower: flo, upper: fhi });
    }
    let (mut lo_sample, mut hi_sample) = (lower, upper);
    let mut t = lo + 0.5 * (hi - lo);
    let mut residual = flo;
    for iteration in 1..=maximum_iterations {
        let value = sample(t)?;
        residual = value.residual;
        if residual.abs() <= tolerance {
            return Ok(Solution { temperature_k: t, sample: value, iterations: iteration });
        }
        if residual < 0.0 { lo = t; lo_sample = value; } else { hi = t; hi_sample = value; }
        let proposed = t - residual / value.tangent;
        let margin = 0.1 * (hi - lo);
        let next = if proposed.is_finite() && proposed > lo + margin && proposed < hi - margin {
            proposed
        } else { lo + 0.5 * (hi - lo) };
        if next <= lo || next >= hi || next == t {
            // The bracket cannot shrink: its ends are adjacent representable
            // temperatures (or the iterate sits on one). The better end is the
            // closest representable endpoint; accept it only within the floor.
            let (end, end_sample) = if lo_sample.residual.abs() <= hi_sample.residual.abs() {
                (lo, lo_sample)
            } else {
                (hi, hi_sample)
            };
            if end_sample.residual.abs() <= floor(end, &end_sample) {
                return Ok(Solution { temperature_k: end, sample: end_sample, iterations: iteration });
            }
            return Err(Error::NotConverged { iterations: iteration, residual });
        }
        t = next;
    }
    Err(Error::NotConverged { iterations: maximum_iterations, residual })
}
