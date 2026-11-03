//! Exact planar two-vortex flow for the shared algebraic-core kernel.
//!
//! Splitting the vortex Hamiltonian into pair terms retains its linear and
//! angular impulse invariants. Each pair keeps its separation and therefore
//! rotates with constant angular speed `(Gamma_i + Gamma_j)/(2 pi (d²+a²))`.
//! Opposite equal strengths instead translate as a dipole. The displacement
//! form used here is continuous through that zero-total-circulation limit.
//!
//! See J. Frank, *Symplectic and energy conserving discretizations of
//! Hamiltonian PDEs*, pp. 11–12,
//! <https://webspace.science.uu.nl/~zegel101/NUMpde2025/Lecture8Extra2.pdf>.
//! The same radial pair solution applies to the declared constant shared
//! algebraic core; it does not describe viscous core evolution or boundaries.

use crate::{VortexParticle, binary_parts};
use fs_math::det;

/// Refusal of a finite, exact two-vortex subflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairFlowError {
    /// Positions, strengths and time must be finite; the common core must
    /// additionally be nonnegative.
    InvalidInput,
    /// The separation, rotation angle or final position exceeds the admitted
    /// floating-point domain. Neither input particle is mutated.
    Unrepresentable,
}

impl core::fmt::Display for PairFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => "vortex pair needs finite inputs and a nonnegative common core",
            Self::Unrepresentable => "vortex pair flow exceeds the finite numerical domain",
        })
    }
}

impl core::error::Error for PairFlowError {}

/// Advance an isolated pair exactly for the crate's regularized Biot–Savart
/// model, up to floating-point roundoff, without allocation.
///
/// Positive circulation induces counterclockwise velocity. Both particles use
/// the same `core` in metres; circulations are in m²/s and `dt` is seconds.
/// Negative time is supported. Zero time or exact coincidence is an identity,
/// consistently with [`crate::induced_velocity`]. Circulation is unchanged.
///
/// The flow preserves `sum Gamma*x`, `sum Gamma*y`, and `sum Gamma*|x|²`.
/// These are vortex impulse moments, not the momentum inventory of an
/// arbitrary finite fluid control volume. Applying the pair map repeatedly
/// within a larger system requires a declared ordering and time integrator;
/// [`crate::conservative::simulate_with_cx`] supplies a symmetric composition.
///
/// # Errors
/// Refuses invalid scalars or nonrepresentable separation/angle/output. All
/// fields are validated even for zero time. This fixed-work leaf has no heap
/// allocation; the multi-particle entry owns execution and memory admission.
pub fn advance_pair(
    particles: [VortexParticle; 2],
    dt: f64,
    core: f64,
) -> Result<[VortexParticle; 2], PairFlowError> {
    if !dt.is_finite()
        || !core.is_finite()
        || core < 0.0
        || particles.iter().any(|p| {
            !p.circulation.is_finite() || p.pos.iter().any(|coordinate| !coordinate.is_finite())
        })
    {
        return Err(PairFlowError::InvalidInput);
    }
    advance_validated_pair(particles, dt, core, 1.0)
}

pub(super) fn advance_validated_pair(
    mut particles: [VortexParticle; 2],
    dt: f64,
    core: f64,
    time_fraction: f64,
) -> Result<[VortexParticle; 2], PairFlowError> {
    if dt == 0.0 || particles.iter().all(|p| p.circulation == 0.0) {
        return Ok(particles);
    }
    let r = [
        particles[0].pos[0] - particles[1].pos[0],
        particles[0].pos[1] - particles[1].pos[1],
    ];
    if r.iter().any(|x| !x.is_finite()) {
        return Err(PairFlowError::Unrepresentable);
    }
    if r == [0.0; 2] {
        return Ok(particles);
    }
    let scale = r[0].abs().max(r[1].abs()).max(core);
    let normalized_r = [r[0] / scale, r[1] / scale];
    let normalized_core = core / scale;
    // Applying the fixed composition fraction in the denominator keeps a
    // subnormal dt from rounding to zero before circulation can scale it.
    let denominator = core::f64::consts::TAU / time_fraction
        * (normalized_r[0] * normalized_r[0]
            + normalized_r[1] * normalized_r[1]
            + normalized_core * normalized_core);
    let gamma = [particles[0].circulation, particles[1].circulation];
    let total = gamma[0] + gamma[1];
    let angle = if total.is_finite() {
        scaled_flow(total, dt, 1.0, scale, denominator, 2)
    } else {
        // Like-signed finite strengths may overflow only when summed. Halve
        // first and restore the factor inside the scaled product.
        scaled_flow(
            0.5 * gamma[0] + 0.5 * gamma[1],
            dt,
            2.0,
            scale,
            denominator,
            2,
        )
    };
    if !angle.is_finite() {
        return Err(PairFlowError::Unrepresentable);
    }
    let (sinc, cosc) = rotation_functions(angle);
    let shape = [
        -cosc * normalized_r[0] - sinc * normalized_r[1],
        sinc * normalized_r[0] - cosc * normalized_r[1],
    ];
    for axis in 0..2 {
        particles[0].pos[axis] += scaled_flow(gamma[1], dt, shape[axis], scale, denominator, 1);
        particles[1].pos[axis] -= scaled_flow(gamma[0], dt, shape[axis], scale, denominator, 1);
    }
    if particles
        .iter()
        .any(|p| p.pos.iter().any(|x| !x.is_finite()))
    {
        return Err(PairFlowError::Unrepresentable);
    }
    Ok(particles)
}

// sinc(theta), (1-cos(theta))/theta. Taylor evaluation avoids cancellation
// near a dipole, without imposing a dimensional circulation threshold.
fn rotation_functions(angle: f64) -> (f64, f64) {
    if angle.abs() < 0.001 {
        let square = angle * angle;
        (
            1.0 + square * (-1.0 / 6.0 + square * (1.0 / 120.0 - square / 5040.0)),
            angle * (0.5 + square * (-1.0 / 24.0 + square * (1.0 / 720.0 - square / 40320.0))),
        )
    } else {
        let (sin, cos) = det::sin_cos(angle);
        (sin / angle, (1.0 - cos) / angle)
    }
}

// gamma * dt * component / (scale^power * denominator). Only exceptional
// ranges need exponent arithmetic; ordinary pairs use the direct expression.
fn scaled_flow(
    gamma: f64,
    dt: f64,
    component: f64,
    scale: f64,
    denominator: f64,
    power: i32,
) -> f64 {
    if gamma == 0.0 || dt == 0.0 || component == 0.0 {
        return 0.0;
    }
    let divisor = if power == 2 {
        scale * scale * denominator
    } else {
        scale * denominator
    };
    let rate = gamma / divisor;
    let duration = rate * dt;
    let ordinary = duration * component;
    if divisor.is_normal() && rate.is_normal() && duration.is_normal() && ordinary.is_normal() {
        return ordinary;
    }
    let sign = if gamma.is_sign_negative() ^ dt.is_sign_negative() ^ component.is_sign_negative() {
        -1.0
    } else {
        1.0
    };
    let (g, ge) = binary_parts(gamma.abs());
    let (t, te) = binary_parts(dt.abs());
    let (c, ce) = binary_parts(component.abs());
    let (s, se) = binary_parts(scale);
    let divisor_mantissa = if power == 2 { s * s } else { s };
    let (mantissa, adjustment) = binary_parts(g * t * c / (divisor_mantissa * denominator));
    let exponent = ge + te + ce - power * se + adjustment;
    let magnitude = if exponent > 1023 {
        f64::INFINITY
    } else if exponent < -1075 {
        0.0
    } else if exponent >= -1022 {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else {
        let shifted = f64::from_bits(((exponent + 1023 + 54) as u64) << 52);
        (mantissa * shifted) * (1.0 / 18_014_398_509_481_984.0)
    };
    magnitude.copysign(sign)
}
