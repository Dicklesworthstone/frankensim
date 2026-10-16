//! Scalar-gap condensation of an arbitrary linear mechanical midpoint response.
//! Pressure and swept volume use the same signed projection. Gap/contact use
//! their own projection: modes that hardly move the slit can still sweep fluid.
use super::dynamic::ApertureDrive;
use super::paired::PairedApertureLimits;
use crate::acoustic_realize::AcousticRealizeError;
use crate::unilateral_contact::distributed::ApertureContactStep;
use fs_dcontact::Obstacle;

pub(super) const MAX_MODES: usize = 32;
#[derive(Clone, Copy)]
pub(super) struct UniformTrial {
    pub velocity: [f64; MAX_MODES],
    pub dp: f64,
    pub residual: f64,
    pub scale: f64,
}
fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Condense v = free + mobility*(gap_shape*contact - pressure_area*dp).
/// The caller owns midpoint mechanics, material memories and final work checks.
/// Generalizes the paired-aperture root algorithm with explicit signed projections.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve(
    free: &[f64], mobility: &[f64], area: &[f64], gap_shape: &[f64],
    gap: f64, dt: f64, width: f64, rho: f64, z: f64, drive: ApertureDrive,
    limits: PairedApertureLimits, obstacle: &Obstacle,
) -> Result<UniformTrial, AcousticRealizeError> {
    let n = free.len();
    if n == 0 || n > MAX_MODES || [mobility.len(), area.len(), gap_shape.len()]
        .iter().any(|&len| len != n) {
        return Err(invalid("uniform aperture needs matching bounded modal projections"));
    }
    let mut pressure_response = [0.0; MAX_MODES];
    let mut contact_response = [0.0; MAX_MODES];
    for i in 0..n {
        pressure_response[i] = mobility[i]*area[i];
        contact_response[i] = mobility[i]*gap_shape[i];
    }
    solve_response(free, &pressure_response[..n], &contact_response[..n], area,
        gap_shape, gap, dt, width, rho, z, drive, limits, obstacle)
}

/// Same scalar contact/flow solve with an already condensed FULL mechanical
/// response: v = free + contact_response*F - pressure_response*dp. A symmetric
/// passive mechanical owner provides both columns from the SAME factorization.
/// Off-diagonal damping therefore does not need another Bernoulli/contact law.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_response(
    free: &[f64], pressure_response: &[f64], contact_response: &[f64],
    area: &[f64], gap_shape: &[f64], gap: f64, dt: f64, width: f64, rho: f64,
    z: f64, drive: ApertureDrive, limits: PairedApertureLimits, obstacle: &Obstacle,
) -> Result<UniformTrial, AcousticRealizeError> {
    let n = free.len();
    if n == 0 || n > MAX_MODES
        || [pressure_response.len(), contact_response.len(), area.len(), gap_shape.len()]
            .iter().any(|&len| len != n)
        || free.iter().chain(pressure_response).chain(contact_response).chain(area)
            .chain(gap_shape).any(|x| !x.is_finite()) {
        return Err(invalid("uniform aperture needs matching finite mechanical response columns"));
    }
    let contact = ApertureContactStep::new(obstacle, gap)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let pressure_mobility = (0..n).map(|i| pressure_response[i]*area[i]).sum::<f64>();
    let denominator = 1.0 + z*pressure_mobility;
    if !pressure_mobility.is_finite() || pressure_mobility < 0.0
        || !denominator.is_finite() || denominator <= 0.0 {
        return Err(invalid("uniform aperture swept-flow response overflowed"));
    }
    let mut evaluations = 0;
    let mut evaluate = |velocity: f64| -> Result<UniformTrial, AcousticRealizeError> {
        evaluations += 1;
        if evaluations > limits.max_evaluations {
            return Err(invalid("uniform aperture nonlinear work budget exhausted"));
        }
        let next_gap = gap + dt*velocity;
        let force = contact.response(next_gap, velocity)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?.force;
        let mut velocities = [0.0; MAX_MODES];
        let mut moving = drive.body_flow_m3_s;
        for i in 0..n {
            velocities[i] = free[i] + contact_response[i]*force;
            moving -= area[i]*velocities[i];
        }
        let source = (drive.upstream_pressure_pa - 2.0*drive.incoming_pressure_pa - z*moving)/denominator;
        let opening = f64::midpoint(gap, next_gap).max(0.0);
        let dp = super::moving::fixed_opening_pressure_drop(width, opening, rho,
            z/denominator, source)?;
        let mut residual = velocity;
        let mut scale = velocity.abs();
        for i in 0..n {
            velocities[i] -= pressure_response[i]*dp;
            residual -= gap_shape[i]*velocities[i];
            scale += (gap_shape[i]*velocities[i]).abs() + (gap_shape[i]*free[i]).abs();
        }
        if ![velocity, next_gap, dp, residual, scale].iter().all(|x| x.is_finite())
            || velocities[..n].iter().any(|x| !x.is_finite()) {
            return Err(invalid("uniform aperture trial left the finite physical domain"));
        }
        Ok(UniformTrial { velocity: velocities, dp, residual, scale })
    };
    let resolved = |t: UniformTrial| t.residual.abs() <= limits.velocity_tolerance_m_s
        + limits.relative_tolerance*t.scale;
    let initial_velocity = free.iter().zip(gap_shape).map(|(v,g)| v*g).sum();
    let initial = evaluate(initial_velocity)?;
    if resolved(initial) { return Ok(initial); }
    let mut span = free.iter().zip(gap_shape).map(|(v,g)| (v*g).abs()).sum::<f64>().max(0.01);
    let (mut lo, mut hi) = (-span, span);
    let (mut left, mut right) = (evaluate(lo)?, evaluate(hi)?);
    while left.residual > 0.0 || right.residual < 0.0 {
        span *= 2.0;
        if !span.is_finite() { return Err(invalid("uniform aperture root has no finite bracket")); }
        if left.residual > 0.0 { lo = -span; left = evaluate(lo)?; }
        if right.residual < 0.0 { hi = span; right = evaluate(hi)?; }
    }
    loop {
        if resolved(left) { return Ok(left); }
        if resolved(right) { return Ok(right); }
        let mid = f64::midpoint(lo, hi);
        if mid == lo || mid == hi { return Err(invalid("uniform aperture has unresolved gap momentum")); }
        let next = evaluate(mid)?;
        if resolved(next) { return Ok(next); }
        if next.residual < 0.0 { lo = mid; left = next; }
        else { hi = mid; right = next; }
    }
}
