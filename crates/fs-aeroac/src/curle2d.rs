//! 2D frequency-domain Curle radiation: compact dipole sources over
//! the OUTGOING cylindrical Green's function.
//!
//! With the workspace `e^{-i omega t}` convention the free-space 2D
//! Helmholtz Green's function (`(nabla^2 + k^2) G = -delta`) is
//! `G(r) = (i/4) H0^(1)(kr)`, and a compact rigid body with unsteady
//! surface force spectrum `F_hat` (force per unit span, the Curle
//! surface-pressure dipole at low Mach) radiates
//!
//! `p_hat(x) = -F_hat_i d/dx_i G(|x - y|)
//!           = (i k / 4) H1^(1)(kr) (r_hat . F_hat)`
//!
//! since `dH0/dr = -k H1`. This is the 2D LINE-source form —
//! cylindrical `1/sqrt(r)` spreading, `cos(theta)` directivity — and
//! the trap the bead's polish round names (naive 3D `1/r` formulas
//! over 2D data give wrong spectral slopes) is structurally excluded
//! by working in the frequency domain with the Hankel kernel.
//! Everything here is a SHAPE/SCALING authority
//! ([`crate::SCOPE_STATEMENT`]); absolute SPL additionally needs a
//! 2D-to-3D span correction.

use crate::bessel::hankel1_outgoing;
use crate::{AeroacError, SCOPE_STATEMENT};
use fs_math::c64::C64;
use fs_math::det;

/// One radiated-pressure evaluation.
#[derive(Debug, Clone)]
pub struct DipoleField {
    /// Complex pressure amplitude at the observer [Pa per unit-span
    /// force], `e^{-i omega t}` convention.
    pub pressure: C64,
    /// Source-observer distance [m].
    pub radius: f64,
    /// The honest-scope statement, embedded in every output (the
    /// marketing-mutation guard asserts it).
    pub scope: &'static str,
}

/// Radiated pressure of a compact dipole with force-per-unit-span
/// spectrum `force` [N/m] at wavenumber `k` [1/m], observer at `obs`,
/// source at `src` (both [m]).
///
/// # Errors
/// [`AeroacError::NonFinite`] on bad inputs;
/// [`AeroacError::InvalidParameter`] for `k <= 0` or coincident
/// observer/source (the field is singular at the source).
pub fn dipole_pressure(
    force: [C64; 2],
    k: f64,
    obs: [f64; 2],
    src: [f64; 2],
) -> Result<DipoleField, AeroacError> {
    if !k.is_finite()
        || obs.iter().chain(&src).any(|v| !v.is_finite())
        || force.iter().any(|c| !c.re.is_finite() || !c.im.is_finite())
    {
        return Err(AeroacError::NonFinite {
            what: "dipole inputs",
        });
    }
    if k <= 0.0 {
        return Err(AeroacError::InvalidParameter {
            what: "wavenumber must be positive",
        });
    }
    let dx = obs[0] - src[0];
    let dy = obs[1] - src[1];
    let r = det::sqrt(dx * dx + dy * dy);
    if r == 0.0 {
        return Err(AeroacError::InvalidParameter {
            what: "observer coincides with the source (singular field)",
        });
    }
    let rhat = [dx / r, dy / r];
    let radial_force = force[0].scale(rhat[0]) + force[1].scale(rhat[1]);
    // (i k / 4) H1(kr) (rhat . F)
    let h1 = hankel1_outgoing(k * r);
    let ik4 = C64::new(0.0, k / 4.0);
    Ok(DipoleField {
        pressure: ik4 * h1 * radial_force,
        radius: r,
        scope: SCOPE_STATEMENT,
    })
}
