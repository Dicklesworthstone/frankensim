//! Shared native reference weights and physical-occupancy integration.
use super::{Sdf3Elasticity, failure};
use crate::{EvaluationStop, SolveControl};

pub(super) fn mass_weights(op: &impl Sdf3Elasticity) -> Result<&[f64], EvaluationStop> {
    let weights = op
        .reference_cell_masses_kg()
        .ok_or_else(|| failure("sdf3-material-mass-unavailable"))?;
    if weights.len() != op.cells() || weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
        return Err(failure("sdf3-material-mass-weights"));
    }
    Ok(weights)
}

// Both kg and normalized volume are linear in PHYSICAL occupancy. The native
// owner supplies the weights; SIMP stiffness and its ersatz floor do not enter.
pub(super) fn weighted_occupancy(
    weights: &[f64],
    projected: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    if weights.len() != projected.len() || weights.is_empty() {
        return Err(failure("sdf3-inventory-shape"));
    }
    let mut value = 0.0;
    for (i, (&weight, &rho)) in weights.iter().zip(projected).enumerate() {
        if i % 1024 == 0 {
            control.checkpoint("sdf3-inventory")?;
        }
        if !weight.is_finite() || weight <= 0.0 || !rho.is_finite() || !(0.0..=1.0).contains(&rho) {
            return Err(failure("sdf3-inventory-input"));
        }
        let term = weight * rho;
        if !term.is_finite() || (rho != 0.0 && term == 0.0) {
            return Err(failure("sdf3-inventory-range"));
        }
        value += term;
        if !value.is_finite() {
            return Err(failure("sdf3-inventory-range"));
        }
    }
    Ok(value)
}

pub(super) fn volume_fraction(
    op: &impl Sdf3Elasticity,
    projected: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    let volumes = op.volumes();
    if volumes.len() != op.cells() || volumes.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(failure("sdf3-inventory-volume-weights"));
    }
    let total: f64 = volumes.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(failure("sdf3-inventory-volume-weights"));
    }
    let weights: Vec<_> = volumes.iter().map(|v| v / total).collect();
    weighted_occupancy(&weights, projected, control)
}
