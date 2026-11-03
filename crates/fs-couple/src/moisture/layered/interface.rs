//! Equal/opposite Kirchhoff fluxes through material and transfer resistances.
use super::*;

/// kg water / (s * unit activity), including the physical area and dry density.
pub(super) fn half_conductance(
    layer: &MoistureLayer,
    diffusion: &KirchhoffMoistureDiffusivity,
    a: f64,
    b: f64,
    pullback: f64,
) -> Result<f64, MoistureError> {
    let left = layer.sorption.moisture(a)?;
    let right = layer.sorption.moisture(b)?;
    let mean = diffusion.mean_diffusivity(left, right)?;
    if mean == 0.0 {
        return Ok(0.0);
    }
    let capacity = layer.sorption.mean_capacity(a, b)?;
    let dx = layer.spacing();
    let value = 2.0 * layer.slab.dry_mass_per_cell_kg / dx / dx * mean * capacity * pullback;
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
    pullback: [f64; 2],
    cx: &Cx<'_>,
) -> Result<(f64, f64), MoistureError> {
    let evaluate = |fraction: f64| -> Result<(f64, f64), MoistureError> {
        let a = fraction.mul_add(activity[1] - activity[0], activity[0]);
        let left = half_conductance(layers[0], profiles[0], activity[0], a, pullback[0])?;
        let right = half_conductance(layers[1], profiles[1], a, activity[1], pullback[1])?;
        series(left, right)
    };
    locate(activity, evaluate, cx)
}

fn locate(
    activity: [f64; 2],
    evaluate: impl Fn(f64) -> Result<(f64, f64), MoistureError>,
    cx: &Cx<'_>,
) -> Result<(f64, f64), MoistureError> {
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

/// Eliminate a surface law against the original nonlinear half-cell integral.
pub(super) fn surface(
    layer: &MoistureLayer,
    profile: &KirchhoffMoistureDiffusivity,
    activity: [f64; 2],
    conductance: Option<f64>,
    pullback: f64,
    cx: &Cx<'_>,
) -> Result<(f64, f64), MoistureError> {
    match conductance {
        None => Ok((
            half_conductance(layer, profile, activity[0], activity[1], pullback)?,
            activity[1],
        )),
        Some(0.0) => Ok((0.0, activity[0])),
        Some(h) => locate(
            activity,
            |fraction| {
                let a = fraction.mul_add(activity[1] - activity[0], activity[0]);
                series(
                    half_conductance(layer, profile, activity[0], a, pullback)?,
                    h,
                )
            },
            cx,
        ),
    }
}

/// Stable three-resistance series law and left/right half-cell activity drops.
/// Small drops stay separate even when the two represented surface activities
/// coincide to floating precision at a very high transfer conductance.
fn triple(left: f64, middle: f64, right: f64) -> Result<(f64, [f64; 2]), MoistureError> {
    // A trial may lie on a zero-D plateau even when the complete cell-to-cell
    // path conducts. Its limiting series map must move the surface bracket
    // through the plateau, not form 0/0 or reject the underlying material.
    if left == 0.0 {
        return Ok((0.0, if right == 0.0 { [0.5, 0.5] } else { [1.0, 0.0] }));
    }
    if right == 0.0 {
        return Ok((0.0, [0.0, 1.0]));
    }
    let scale = left.min(middle).min(right);
    let drops = [scale / left, scale / middle, scale / right];
    let sum = drops[0] + drops[1] + drops[2];
    let g = scale / sum;
    if !positive(g) {
        return Err(MoistureError::Invalid(
            "finite activity-transfer conductance underflow",
        ));
    }
    Ok((g, [drops[0] / sum, drops[2] / sum]))
}

/// A finite material join has two activities and a single water flux. Eliminate
/// the transfer jump analytically, leaving one monotone scalar root for the
/// left surface. Its analytic differential conductances accelerate a bounded
/// bracket; no nested roots, dense mesh system or diffusion splitting is used.
pub(super) fn transfer_join(
    layers: [&MoistureLayer; 2],
    profiles: [&KirchhoffMoistureDiffusivity; 2],
    activity: [f64; 2],
    conductance: Option<f64>,
    pullback: [f64; 2],
    cx: &Cx<'_>,
) -> Result<(f64, [f64; 2]), MoistureError> {
    let Some(h) = conductance else {
        let (g, a) = join(layers, profiles, activity, pullback, cx)?;
        return Ok((g, [a; 2]));
    };
    if h == 0.0 {
        return Ok((0.0, activity));
    }
    let delta = activity[1] - activity[0];
    let gl = half_conductance(
        layers[0],
        profiles[0],
        activity[0],
        activity[1],
        pullback[0],
    )?;
    let gr = half_conductance(
        layers[1],
        profiles[1],
        activity[0],
        activity[1],
        pullback[1],
    )?;
    // A completely impermeable half path carries no water. Its interface
    // activities follow the transmitting neighbor and the finite transfer law.
    if gl == 0.0 || gr == 0.0 {
        let a = if gr == 0.0 { activity[0] } else { activity[1] };
        return Ok((0.0, [a; 2]));
    }
    let initial = triple(gl, h, gr)?;
    if delta == 0.0 {
        return Ok((initial.0, activity));
    }
    let surfaces = |drops: [f64; 2]| {
        [
            drops[0].mul_add(delta, activity[0]),
            (-drops[1]).mul_add(delta, activity[1]),
        ]
    };
    let evaluate = |drops: [f64; 2]| -> Result<(f64, [f64; 2]), MoistureError> {
        let a = surfaces(drops);
        triple(
            half_conductance(layers[0], profiles[0], activity[0], a[0], pullback[0])?,
            h,
            half_conductance(layers[1], profiles[1], a[1], activity[1], pullback[1])?,
        )
    };
    // The affine case needs only one consistency check and retains a stable
    // series calculation even when one physical resistance is extremely small.
    let mut coupling = initial;
    let mut bracket = [0.0, 1.0];
    let mut fraction = coupling.1[0];
    let mut iterations = 0;
    for iteration in 1..=128 {
        iterations = iteration;
        poll(cx)?;
        let checked = evaluate(coupling.1)?;
        if close(coupling.0, checked.0)
            && close(coupling.1[0], checked.1[0])
            && close(coupling.1[1], checked.1[1])
        {
            return Ok((checked.0, surfaces(checked.1)));
        }
        let left_surface = fraction.mul_add(delta, activity[0]);
        let left = half_conductance(
            layers[0],
            profiles[0],
            activity[0],
            left_surface,
            pullback[0],
        )?;
        let flux = left * fraction;
        let right_fraction = fraction + flux / h;
        if !right_fraction.is_finite() || right_fraction >= 1.0 {
            bracket[1] = fraction;
            fraction = bracket[0] + 0.5 * (bracket[1] - bracket[0]);
            coupling = evaluate([fraction, 1.0 - fraction])?;
            continue;
        }
        let right_surface = right_fraction.mul_add(delta, activity[0]);
        let right = half_conductance(
            layers[1],
            profiles[1],
            right_surface,
            activity[1],
            pullback[1],
        )?;
        if left == 0.0 && right == 0.0 {
            // Both actual half paths are impermeable at this common surface,
            // so the transfer jump and all three physical fluxes are zero.
            return Ok((0.0, [left_surface, right_surface]));
        }
        coupling = triple(left, h, right)?;
        let residual = flux - right * (1.0 - right_fraction);
        if residual < 0.0 {
            bracket[0] = fraction;
        } else {
            bracket[1] = fraction;
        }
        let dl = half_conductance(
            layers[0],
            profiles[0],
            left_surface,
            left_surface,
            pullback[0],
        )?;
        let dr = half_conductance(
            layers[1],
            profiles[1],
            right_surface,
            right_surface,
            pullback[1],
        )?;
        let slope = dl + dr + (dl / h) * dr;
        let trial = fraction - residual / slope;
        let next =
            if slope.is_finite() && positive(slope) && trial > bracket[0] && trial < bracket[1] {
                trial
            } else {
                bracket[0] + 0.5 * (bracket[1] - bracket[0])
            };
        if next == fraction {
            break;
        }
        fraction = next;
    }
    Err(MoistureError::LayeredInterfaceNotConverged {
        iterations,
        activity_bracket: bracket.map(|f| f.mul_add(delta, activity[0])),
    })
}
