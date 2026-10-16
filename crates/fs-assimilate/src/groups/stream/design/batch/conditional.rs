//! Square-root Joseph conditioning of feature rows with a shared noise column.
use super::super::{dot, finite, norm, product, reserve, square};
use crate::{AssimError, WorkProgress};

pub(super) struct Conditional {
    pub residual: Vec<f64>,
    pub reduction: f64,
    pub innovation_variance: f64,
}

pub(super) fn condition(
    target: &[f64],
    sensor: &[f64],
    noise_variance: f64,
    progress: &mut WorkProgress<'_, '_>,
) -> Result<Conditional, AssimError> {
    debug_assert_eq!(target.len(), sensor.len());
    let mut measurement = reserve(sensor.len() + 1)?;
    measurement.extend_from_slice(sensor);
    let noise_root = noise_variance.sqrt();
    measurement.push(noise_root);
    let innovation_root = norm(&measurement, progress)?;
    let innovation_variance = square(innovation_root)?;
    let target_root = norm(target, progress)?;
    let mut residual = reserve(target.len() + 1)?;
    if target_root == 0.0 {
        residual.resize(target.len() + 1, 0.0);
        return Ok(Conditional { residual, reduction: 0.0, innovation_variance });
    }
    let mut tn = reserve(target.len())?;
    let mut sn = reserve(sensor.len())?;
    for (t, s) in target.iter().zip(sensor) {
        tn.push(*t / target_root);
        sn.push(*s / innovation_root);
        progress.scalar("batch probe normalization", 2)?;
    }
    let correlation = dot(&tn, &sn, progress)?;
    if correlation.abs() > 1.0 { return Err(AssimError::CovarianceCertificationUnresolved); }
    let reduction = square(product(target_root, correlation)?)?;
    if correlation == 0.0 {
        // Preserve the original bits when this observation has no influence.
        residual.extend_from_slice(target);
        residual.push(0.0);
    } else {
        for (t, s) in tn.iter().zip(&sn) {
            residual.push(product(target_root, (-correlation).mul_add(*s, *t))?);
            progress.scalar("batch probe Joseph action", 3)?;
        }
        // The minus sign is essential for cross-covariance between TWO rows.
        residual.push(product(product(target_root, -correlation)?, noise_root / innovation_root)?);
    }
    let posterior_root = norm(&residual, progress)?;
    if posterior_root > target_root || reduction > square(target_root)? {
        return Err(AssimError::CovarianceCertificationUnresolved);
    }
    finite(posterior_root)?;
    Ok(Conditional { residual, reduction, innovation_variance })
}
