//! Natural Navier/free-slip traces for the compatible Cartesian stress.
//!
//! At a wall, eliminate its tangential fluid trace from viscous cell
//! dissipation plus mu/b times the squared slip. The centered off-diagonal
//! derivatives satisfy G_tn + alpha_n G_nt = beta_tn. Only each unordered
//! derivative pair is coupled, including corners and singleton dimensions.
//! Complement weights are retained independently: 1-alpha loses finite
//! information for either extremely small or extremely large slip lengths.

use super::{Viscosity, no_slip_traction};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::compressible) enum Slip {
    /// Natural zero tangential traction, not an even tangential ghost.
    Free,
    /// Finite nonnegative slip length; zero is exactly prescribed velocity.
    Length(f64),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::compressible) enum GradientSample<const D: usize> {
    /// An actual neighbor, periodic neighbor, or the center at a zero-flux end.
    Cell([f64; D]),
    Wall { velocity: [f64; D], slip: Slip },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::compressible) struct WallStress<const D: usize> {
    /// Global tau[a, axis], oriented along the positive coordinate axis.
    pub traction: [f64; D],
    /// Physical interfacial friction heat per area, excluding wall work.
    pub slip_dissipation: f64,
}

impl Slip {
    // Return phi=d/(d+b) and alpha=b/(d+b) independently. Neither fraction is
    // recovered by subtracting the other from one.
    fn weights(self, distance: f64) -> Result<(f64, f64), &'static str> {
        if !distance.is_finite() || distance <= 0.0 {
            return Err("positive finite compatible wall distance required");
        }
        match self {
            Self::Free => Ok((0.0, 1.0)),
            Self::Length(length) => {
                if !length.is_finite() || length < 0.0 {
                    return Err("finite nonnegative compatible slip length required");
                }
                if length == 0.0 {
                    return Ok((1.0, 0.0));
                }
                let (phi, alpha) = if length <= distance {
                    let ratio = length / distance;
                    (1.0 / (1.0 + ratio), ratio / (1.0 + ratio))
                } else {
                    let ratio = distance / length;
                    (ratio / (1.0 + ratio), 1.0 / (1.0 + ratio))
                };
                if phi == 0.0 || alpha == 0.0 {
                    return Err("finite Navier geometric weight is unrepresentable");
                }
                Ok((phi, alpha))
            }
        }
    }
}

// The series length need not fit in f64. A second ordering also avoids an
// intermediate overflowing quotient when division by the factor <=2 rescues it.
fn divide_length(value: f64, distance: f64, length: f64) -> f64 {
    let scale = distance.max(length);
    let denominator = distance / scale + length / scale;
    let scaled = value / scale;
    if scaled.is_finite() {
        scaled / denominator
    } else {
        (value / denominator) / scale
    }
}

fn product_over_length(a: f64, b: f64, distance: f64, length: f64) -> f64 {
    let (larger, smaller) = (a.max(b.abs()), a.min(b.abs()));
    let quotient = divide_length(larger, distance, length);
    let value = if quotient.is_finite() {
        quotient * smaller
    } else {
        let other = divide_length(smaller, distance, length);
        if other.is_finite() {
            other * larger
        } else {
            divide_length(larger * smaller, distance, length)
        }
    };
    value.copysign(b)
}

fn average(a: f64, b: f64) -> f64 {
    let sum = a + b;
    if sum.is_finite() { 0.5 * sum } else { 0.5 * a + 0.5 * b }
}

fn numerator(x: f64, y: f64, alpha: f64, delta: f64) -> f64 {
    if delta < 0.5 {
        delta.mul_add(y, x - y)
    } else {
        (-alpha).mul_add(y, x)
    }
}

// One uncoupled equation fixes its derivative exactly. Evaluate the other
// equation without overflowing x-y when its weighted difference still fits.
fn uncoupled_difference(x: f64, y: f64, alpha: f64, delta: f64) -> f64 {
    let direct = numerator(x, y, alpha, delta);
    if direct.is_finite() {
        direct
    } else {
        let scale = x.abs().max(y.abs());
        numerator(x / scale, y / scale, alpha, delta) * scale
    }
}

// Compute numerator*scale/determinant without destroying a small normalized
// RHS or overflowing a quotient that a small scale would bring back into range.
fn rescale(numerator: f64, scale: f64, determinant: f64) -> f64 {
    let quotient = numerator / determinant;
    if quotient.is_finite() {
        quotient * scale
    } else {
        let other = scale / determinant;
        if other.is_finite() {
            numerator * other
        } else {
            (numerator * scale) / determinant
        }
    }
}

/// Source-independent, stack-only physical velocity-gradient reconstruction.
///
/// Return indexing is [velocity component][derivative]. The callback receives
/// axis and side (0 lower, 1 upper). Reflecting sides must supply the actual
/// physical wall; a clamped open side is Cell(center), not Free. Both Cell
/// samples retain the original centered arithmetic. Wall axes use relative
/// increments before elimination, avoiding subtraction of almost equal ghosts.
///
/// A pair of singleton axes with pure free slip on both ends has an
/// undetermined rigid-rotation mode. Only a zero provisional RHS is admitted
/// there, and its minimum-rotation representative is zero. Finite lengths never
/// silently become that null mode. No implicit friction floor is introduced.
pub(in crate::compressible) fn velocity_gradient<const D: usize>(
    spacing: [f64; D],
    velocity: [f64; D],
    mut sample: impl FnMut(usize, usize) -> GradientSample<D>,
) -> Result<[[f64; D]; D], &'static str> {
    if (D != 2 && D != 3) || !velocity.iter().all(|v| v.is_finite()) {
        return Err("compatible velocity-gradient dimensions or state");
    }
    let mut gradient = [[0.0; D]; D];
    let mut alpha = [0.0; D];
    let mut delta = [1.0; D];
    let mut pure_free = [false; D];
    for axis in 0..D {
        let h = spacing[axis];
        if !h.is_finite() || h <= 0.0 {
            return Err("positive finite compatible grid spacing required");
        }
        let lower = sample(axis, 0);
        let upper = sample(axis, 1);
        if let (GradientSample::Cell(lo), GradientSample::Cell(hi)) = (lower, upper) {
            for a in 0..D {
                gradient[a][axis] = (hi[a] - lo[a]) / (2.0 * h);
            }
            continue;
        }
        let distance = 0.5 * h;
        let mut derivative = [[0.0; D]; 2];
        let mut phi = [1.0; 2];
        let mut response = [0.0; 2];
        pure_free[axis] = matches!(lower, GradientSample::Wall { slip: Slip::Free, .. })
            && matches!(upper, GradientSample::Wall { slip: Slip::Free, .. });
        for (side, value) in [lower, upper].into_iter().enumerate() {
            let sign = if side == 0 { -1.0 } else { 1.0 };
            match value {
                GradientSample::Cell(other) => {
                    for a in 0..D {
                        derivative[side][a] = sign * (other[a] - velocity[a]) / h;
                    }
                }
                GradientSample::Wall { velocity: wall, slip } => {
                    if wall[axis] != 0.0 || !wall.iter().all(|v| v.is_finite()) {
                        return Err("compatible wall requires finite tangential velocity");
                    }
                    (phi[side], response[side]) = slip.weights(distance)?;
                    for a in 0..D {
                        let relative = wall[a] - velocity[a];
                        derivative[side][a] = if a == axis {
                            sign * relative / distance
                        } else {
                            match slip {
                                Slip::Free => 0.0,
                                Slip::Length(length) => {
                                    let value = divide_length(relative, distance, length);
                                    if value == 0.0 && relative != 0.0 {
                                        return Err("relative Navier wall gradient is unrepresentable");
                                    }
                                    sign * value
                                }
                            }
                        };
                    }
                }
            }
        }
        alpha[axis] = 0.5 * (response[0] + response[1]);
        if alpha[axis] == 0.0 && (response[0] > 0.0 || response[1] > 0.0) {
            return Err("finite Navier gradient response is unrepresentable");
        }
        delta[axis] = 0.5 * (phi[0] + phi[1]);
        if delta[axis] == 0.0 && !pure_free[axis] {
            return Err("finite Navier gradient complement is unrepresentable");
        }
        for a in 0..D {
            gradient[a][axis] = average(derivative[0][a], derivative[1][a]);
            if gradient[a][axis] == 0.0 && derivative[0][a] + derivative[1][a] != 0.0 {
                return Err("relative Navier centered gradient is unrepresentable");
            }
        }
    }
    if !gradient.iter().flatten().all(|v| v.is_finite()) {
        return Err("nonfinite provisional compatible velocity gradient");
    }
    for n in 0..D {
        for t in n + 1..D {
            let (x, y) = (gradient[t][n], gradient[n][t]);
            let (gn, gt) = if alpha[n] == 0.0 {
                (x, if alpha[t] == 0.0 { y } else { uncoupled_difference(y, x, alpha[t], delta[t]) })
            } else if alpha[t] == 0.0 {
                (uncoupled_difference(x, y, alpha[n], delta[n]), y)
            } else if delta[n] == 0.0 && delta[t] == 0.0 {
                if !pure_free[n] || !pure_free[t] || x != 0.0 || y != 0.0 {
                    return Err("inconsistent natural-slip rotation null mode");
                }
                (0.0, 0.0)
            } else {
                // This positive form retains tiny complements without
                // subtracting a product rounded to one.
                let determinant = delta[n] + alpha[n] * delta[t];
                if !determinant.is_finite() || determinant <= 0.0 {
                    return Err("compatible Navier gradient system is unrepresentable");
                }
                let scale = x.abs().max(y.abs());
                if scale == 0.0 {
                    (0.0, 0.0)
                } else {
                    let (xn, yn) = (x / scale, y / scale);
                    (
                        rescale(numerator(xn, yn, alpha[n], delta[n]), scale, determinant),
                        rescale(numerator(yn, xn, alpha[t], delta[t]), scale, determinant),
                    )
                }
            };
            gradient[t][n] = gn;
            gradient[n][t] = gt;
        }
    }
    if !gradient.iter().flatten().all(|v| v.is_finite()) {
        return Err("nonfinite compatible velocity gradient");
    }
    Ok(gradient)
}

/// Positive-axis wall traction and physical Navier interface dissipation.
///
/// Gradient indexing is component, derivative and must come from the same
/// natural trace elimination as interior stress. Cell shear/bulk viscosities
/// own mechanics; the separate original wall-temperature query owns conduction.
/// For outward sign s, finite b gives t_out=mu*(Uw-u+s*d*C_tu_n)/(d+b)
/// and trace-Uw=-b/(d+b)*(Uw-u+s*d*C_tu_n). Their negative product is
/// interfacial heat, computed without division by mu. A zero shear coefficient
/// therefore has exactly zero tangential traction and interface dissipation.
///
/// Positive-axis conservative energy flux is -Uw dot traction
/// + s*(1-chi)*slip_dissipation + original conductive heat. The retained fluid
/// heat is chi times interface dissipation; no second volume source is added.
#[allow(clippy::too_many_arguments)]
pub(in crate::compressible) fn wall_stress<const D: usize>(
    axis: usize,
    outward_sign: f64,
    wall_distance: f64,
    velocity: [f64; D],
    wall_velocity: [f64; D],
    slip: Slip,
    cell: Viscosity,
    gradient: impl Fn(usize, usize) -> f64,
) -> Result<WallStress<D>, &'static str> {
    if (D != 2 && D != 3) || axis >= D || !matches!(outward_sign, -1.0 | 1.0)
        || !cell.valid() || !velocity.iter().chain(&wall_velocity).all(|v| v.is_finite())
        || wall_velocity[axis] != 0.0
    {
        return Err("compatible wall stress dimensions, state or coefficients");
    }
    let (_, response) = slip.weights(wall_distance)?;
    if matches!(slip, Slip::Length(length) if length == 0.0) {
        let traction = no_slip_traction(
            axis, outward_sign, wall_distance, velocity, wall_velocity, cell, gradient,
        );
        return if traction.iter().all(|v| v.is_finite()) {
            Ok(WallStress { traction, slip_dissipation: 0.0 })
        } else {
            Err("nonfinite compatible no-slip traction")
        };
    }
    let mut result = WallStress { traction: [0.0; D], slip_dissipation: 0.0 };
    let divergence: f64 = (0..D).filter(|&a| a != axis).map(|a| gradient(a, a)).sum();
    result.traction[axis] = ((4.0 / 3.0) * cell.shear + cell.bulk)
        * (-outward_sign * velocity[axis] / wall_distance)
        + (cell.bulk - (2.0 / 3.0) * cell.shear) * divergence;
    if let Slip::Length(length) = slip {
        if cell.shear != 0.0 {
            for a in 0..D {
                if a == axis { continue; }
                let relative = (outward_sign * wall_distance)
                    .mul_add(gradient(axis, a), wall_velocity[a] - velocity[a]);
                if !relative.is_finite() {
                    return Err("nonfinite compatible Navier relative velocity");
                }
                let outward = product_over_length(cell.shear, relative, wall_distance, length);
                if outward == 0.0 && relative != 0.0 {
                    return Err("physical Navier traction is unrepresentable");
                }
                let slip_velocity = -response * relative;
                if slip_velocity == 0.0 && relative != 0.0 && outward != 0.0 {
                    return Err("physical Navier slip trace is unrepresentable");
                }
                result.traction[a] = outward_sign * outward;
                result.slip_dissipation -= outward * slip_velocity;
            }
        }
    }
    if !result.traction.iter().all(|v| v.is_finite())
        || !result.slip_dissipation.is_finite() || result.slip_dissipation < 0.0
    {
        return Err("nonfinite compatible Navier wall transport");
    }
    Ok(result)
}
