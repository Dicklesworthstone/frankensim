//! Bounded model-order selection over the existing vector-fit owner.
//!
//! On a 4*N+1 lattice: 4k are training samples, 4k+2 select the order,
//! and odd samples are a final audit. Audit values never select poles, order,
//! normalization or a retry. Failure of that last audit is terminal.
use super::{Error, C64, DiscreteStateSpace, fit_observer,
    MAX_RELATIVE_ERROR, MAX_RMS_ERROR};

/// Cold work bound, not a pole count silently chosen from instrument identity.
pub const MAX_ORDER: usize = 32;
pub const MAX_INTERVALS: usize = 128;

#[derive(Debug, Clone, Copy)]
pub(super) struct FitReport {
    pub order: usize,
    pub attempts: usize,
    pub selection_maximum: f64,
    pub selection_rms: f64,
    pub audit_maximum: f64,
    pub audit_rms: f64,
}

/// Realized filter and independent frequency-sample audit. This is not a
/// continuous-band error bound, a passivity certificate or a measured fit.
pub(super) fn fit(omega: &[f64], values: &[C64], dt: f64, max_order: usize)
    -> Result<(DiscreteStateSpace, FitReport), Error>
{
    admit(omega, values, dt, max_order)?;
    let (filter, mut report, scale) = select(omega, values, dt, max_order)?;
    let (maximum, rms) = audit(&filter, omega, values, scale)?;
    report.audit_maximum = maximum;
    report.audit_rms = rms;
    Ok((filter, report))
}

fn admit(omega: &[f64], values: &[C64], dt: f64, max_order: usize) -> Result<(), Error> {
    if omega.len() != values.len() || omega.len() < 17 || omega.len() > 4*MAX_INTERVALS+1
        || (omega.len()-1)%4 != 0 || !dt.is_finite() || dt <= 0.0
        || !(2..=MAX_ORDER).contains(&max_order) || max_order%2 != 0
        // Conservative overdetermination: at least 2*order+1 complex training samples.
        || (omega.len()-1)/4 < 2*max_order
        || omega.iter().enumerate().any(|(i,&w)| !w.is_finite() || w <= 0.0
            || w*dt >= core::f64::consts::PI || i > 0 && w <= omega[i-1])
        || values.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
    {
        return Err("observer order selection requires a bounded 4*N+1 grid, finite samples, an even order in 2..=32 and at least 2*order+1 training frequencies".into());
    }
    Ok(())
}

// The old fitter's even/odd split becomes training/selection, not final audit.
// Reusing it retains conjugation, prewarping, proper stable realization, finite
// scaling and the original error thresholds. No second fitting law is introduced.
fn select(omega: &[f64], values: &[C64], dt: f64, max_order: usize)
    -> Result<(DiscreteStateSpace, FitReport, f64), Error>
{
    let frequencies: Vec<_> = omega.iter().step_by(2).copied().collect();
    let response: Vec<_> = values.iter().step_by(2).copied().collect();
    let scale = response.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale == 0.0 && response.iter().any(|v|v.re != 0.0 || v.im != 0.0) {
        return Err("observer training/selection scale overflow or underflow".into());
    }
    let mut last = String::new();
    for (attempt, order) in (2..=max_order).step_by(2).enumerate() {
        match fit_observer(&frequencies, &response, dt, order) {
            Ok((filter, maximum, rms)) => return Ok((filter, FitReport {
                order: if scale == 0.0 { 0 } else { order }, attempts: attempt+1,
                selection_maximum: maximum, selection_rms: rms,
                audit_maximum: 0.0, audit_rms: 0.0,
            }, scale)),
            Err(error) => last = error.to_string(),
        }
    }
    Err(format!("observer order budget exhausted at {max_order}; last candidate: {last}").into())
}

fn audit(filter: &DiscreteStateSpace, omega: &[f64], values: &[C64], scale: f64)
    -> Result<(f64, f64), Error>
{
    let mut maximum = 0.0_f64;
    let mut squares = 0.0;
    let mut count = 0;
    for i in (1..omega.len()).step_by(2) {
        if scale == 0.0 && (values[i].re != 0.0 || values[i].im != 0.0) {
            return Err("observer independent audit sees nonzero data against zero training data".into());
        }
        let defect = (filter.eval(omega[i])?.conj()-values[i]).abs();
        let error = if scale == 0.0 {
            if defect == 0.0 { 0.0 } else { f64::INFINITY }
        } else { defect/scale };
        if !error.is_finite() { return Err("observer independent audit has nonfinite error or a nonzero response against zero training data".into()); }
        maximum = maximum.max(error); squares += error*error; count += 1;
    }
    let rms = (squares/count as f64).sqrt();
    if maximum > MAX_RELATIVE_ERROR || rms > MAX_RMS_ERROR {
        return Err(format!("observer independent audit refused: max={maximum}, rms={rms}; limits={MAX_RELATIVE_ERROR}/{MAX_RMS_ERROR}; audit samples were not reused for order selection").into());
    }
    Ok((maximum, rms))
}

#[cfg(test)]
mod tests;
