//! Equal/opposite half-cell flux through one ideal-contact material interface.
use super::*;

/// kg water / (s * unit activity), including the physical area and dry density.
pub(super) fn half_conductance(
    layer: &MoistureLayer,
    diffusion: &KirchhoffMoistureDiffusivity,
    a: f64,
    b: f64,
) -> Result<f64, MoistureError> {
    let left = layer.sorption.moisture(a)?;
    let right = layer.sorption.moisture(b)?;
    let mean = diffusion.mean_diffusivity(left, right)?;
    if mean == 0.0 {
        return Ok(0.0);
    }
    let capacity = layer.sorption.mean_capacity(a, b)?;
    let dx = layer.spacing();
    let value = 2.0 * layer.slab.dry_mass_per_cell_kg / dx / dx * mean * capacity;
    if !positive(value) {
        return Err(MoistureError::Invalid(
            "layered half-cell conductance overflow or underflow",
        ));
    }
    Ok(value)
}

fn series(left: f64, right: f64) -> Result<(f64, f64), MoistureError> {
    if left == 0.0 {
        return Ok((0.0, if right == 0.0 { 0.5 } else { 1.0 }));
    }
    if right == 0.0 {
        return Ok((0.0, 0.0));
    }
    let coupling = if left <= right {
        let ratio = left / right;
        (left / (1.0 + ratio), 1.0 / (1.0 + ratio))
    } else {
        let ratio = right / left;
        (right / (1.0 + ratio), ratio / (1.0 + ratio))
    };
    if !positive(coupling.0) {
        return Err(MoistureError::Invalid(
            "positive layered interface conductance underflow",
        ));
    }
    Ok(coupling)
}

fn close(a: f64, b: f64) -> bool {
    a == b || (a - b).abs() / a.abs().max(b.abs()) <= 128.0 * f64::EPSILON
}

/// A common interface activity, never a common dry-basis moisture content.
/// The physical residual is monotone because each source D and dM/da is
/// nonnegative. Bisection is limited to 128 steps, including degenerate D=0.
pub(super) fn join(
    layers: [&MoistureLayer; 2],
    profiles: [&KirchhoffMoistureDiffusivity; 2],
    activity: [f64; 2],
    cx: &Cx<'_>,
) -> Result<(f64, f64), MoistureError> {
    let evaluate = |fraction: f64| -> Result<(f64, f64), MoistureError> {
        let a = fraction.mul_add(activity[1] - activity[0], activity[0]);
        let left = half_conductance(layers[0], profiles[0], activity[0], a)?;
        let right = half_conductance(layers[1], profiles[1], a, activity[1])?;
        series(left, right)
    };
    if activity[0] == activity[1] {
        return Ok((evaluate(0.5)?.0, activity[0]));
    }
    let mut bracket: [f64; 2] = [0.0, 1.0];
    let mut iterations = 0;
    for iteration in 1..=128 {
        poll(cx)?;
        iterations = iteration;
        let fraction = bracket[0] + 0.5 * (bracket[1] - bracket[0]);
        let coupling = evaluate(fraction)?;
        let checked = evaluate(coupling.1)?;
        if close(coupling.0, checked.0) && close(coupling.1, checked.1) {
            let surface = checked.1.mul_add(activity[1] - activity[0], activity[0]);
            return Ok((checked.0, surface));
        }
        if fraction == bracket[0] || fraction == bracket[1] {
            break;
        }
        if fraction > coupling.1 {
            bracket[1] = fraction;
        } else {
            bracket[0] = fraction;
        }
    }
    Err(MoistureError::LayeredInterfaceNotConverged {
        iterations,
        activity_bracket: bracket.map(|f| f.mul_add(activity[1] - activity[0], activity[0])),
    })
}
