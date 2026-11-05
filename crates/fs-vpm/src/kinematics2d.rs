//! Analytic derivatives of the original 2-D algebraic-core Biot–Savart field.
//!
//! These are Eulerian samples of the represented velocity, with
//! `omega_z = Gamma a² / (pi (r²+a²)²)` per source. They are distinct from
//! the quadrature vorticity `Gamma/area` used by particle-strength exchange.
//! No Gaussian core, inferred particle area or finite-difference scale enters.

use super::{VelocityKernel, VortexParticle, binary_parts, scale_binary};
use core::f64::consts::PI;

/// One exact-kernel sample, up to floating-point arithmetic and source summation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ParticleKinematics2 {
    /// Induced velocity [m/s], including each signed circulation.
    pub velocity: [f64; 2],
    /// `gradient[i][j] = partial velocity[i] / partial position[j]` [1/s].
    pub gradient: [[f64; 2]; 2],
    /// Analytic curl `dv/dx - du/dy` [1/s], positive counterclockwise.
    /// Evaluated directly from the SAME kernel, before subtracting nearly
    /// equal strain terms; a thin core's small far-field vorticity is retained
    /// even when the rounded Jacobian entries have identical bit patterns.
    pub vorticity: f64,
}

/// Velocity, full spatial Jacobian and vorticity at one point, in one source pass.
///
/// The velocity reuses [`super::induced_velocity`]'s exact kernel and source
/// accumulation order. The analytic derivative has zero divergence. At a
/// positive-core particle center velocity is zero but vorticity is finite;
/// a nonzero point vortex at a zero-core coincident probe is singular and
/// returns `None`. Away from zero-core point vortices the ordinary curl is
/// zero; their distributional circulation is not an ordinary sampled scalar.
/// A zero-strength source contributes no field, including at coincidence.
///
/// Inputs must be finite and `core` nonnegative. `None` also denotes an
/// unrepresentable coordinate difference, contribution or accumulated field.
/// Binary scaling retains representable derivatives when a squared radius,
/// fourth power or intermediate quotient alone overflows or underflows.
///
/// O(N) work, O(1) stack storage and no allocation or source mutation. Like
/// the existing per-probe velocity and exact 3-D kinematics helpers, this
/// function has no `Cx` or cancellation receipt; callers bound bulk queries.
#[must_use]
pub fn induced_kinematics(
    particles: &[VortexParticle],
    point: [f64; 2],
    core: f64,
) -> Option<ParticleKinematics2> {
    if !core.is_finite() || core < 0.0 || point.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let kernel = VelocityKernel::new(core);
    let mut result = ParticleKinematics2::default();
    for source in particles {
        if !source.circulation.is_finite() || source.pos.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let r = [point[0] - source.pos[0], point[1] - source.pos[1]];
        if r.iter().any(|v| !v.is_finite()) {
            return None;
        }
        if source.circulation == 0.0 {
            continue;
        }
        let scale = r[0].abs().max(r[1].abs()).max(core);
        if scale == 0.0 {
            return None;
        }
        let x = r[0] / scale;
        let y = r[1] / scale;
        let a = core / scale;
        let radius2 = x * x + y * y + a * a;
        let normalization = PI * radius2 * radius2;

        // With q=r²+a² and C=Gamma/(2 pi), the Jacobian is
        // [ 2Cxy/q², C(y²-x²-a²)/q² ;
        //   C(y²-x²+a²)/q², -2Cxy/q² ].
        // Form the trace-free strain and rotation separately. Factoring
        // y²-x² as (|y|-|x|)(|y|+|x|) avoids cancellation of rounded squares;
        // the sum is normalized before addition, so it cannot overflow.
        let strain_xx = product(source.circulation, r[0], r[1], scale, 4, normalization);
        let strain_xy = product(
            source.circulation,
            r[1].abs() - r[0].abs(),
            y.abs() + x.abs(),
            scale,
            3,
            2.0 * normalization,
        );
        let (curl_m, curl_e) =
            product_parts(source.circulation, core, core, scale, 4, normalization);
        let curl = scale_binary(curl_m, curl_e);
        // Restore each output from the unrounded exponent form. Halving a
        // rounded subnormal curl would introduce a second rounding in spin.
        let spin = scale_binary(curl_m, curl_e - 1);
        let velocity = kernel.contribution(r, source.circulation);
        result.velocity[0] += velocity[0];
        result.velocity[1] += velocity[1];
        result.gradient[0][0] += strain_xx;
        result.gradient[0][1] += strain_xy - spin;
        result.gradient[1][0] += strain_xy + spin;
        result.gradient[1][1] = -result.gradient[0][0];
        result.vorticity += curl;
        if !result
            .velocity
            .iter()
            .chain(result.gradient.iter().flatten())
            .chain(core::iter::once(&result.vorticity))
            .all(|v| v.is_finite())
        {
            return None;
        }
    }
    Some(result)
}

// Gamma*a*b / (scale^power * normalization), keeping every original factor
// until exponent assembly. Squaring a normalized tiny core first would erase
// a representable curl when a large circulation subsequently restores it.
fn product(circulation: f64, a: f64, b: f64, scale: f64, power: i32, normalization: f64) -> f64 {
    let (mantissa, exponent) = product_parts(circulation, a, b, scale, power, normalization);
    scale_binary(mantissa, exponent)
}

fn product_parts(
    circulation: f64,
    a: f64,
    b: f64,
    scale: f64,
    power: i32,
    normalization: f64,
) -> (f64, i32) {
    let negative = circulation.is_sign_negative() ^ a.is_sign_negative() ^ b.is_sign_negative();
    let sign = if negative { -1.0 } else { 1.0 };
    if circulation == 0.0 || a == 0.0 || b == 0.0 {
        return (0.0_f64.copysign(sign), 0);
    }
    let (gm, ge) = binary_parts(circulation.abs());
    let (am, ae) = binary_parts(a.abs());
    let (bm, be) = binary_parts(b.abs());
    let (sm, se) = binary_parts(scale);
    let scale2 = sm * sm;
    let scale_power = if power == 3 {
        scale2 * sm
    } else {
        scale2 * scale2
    };
    let fraction = (gm * am * bm) / (scale_power * normalization);
    let (mantissa, correction) = binary_parts(fraction);
    (
        mantissa.copysign(sign),
        ge + ae + be - power * se + correction,
    )
}
